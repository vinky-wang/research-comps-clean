# Fitting models -------------------------------------------------------------
# Self-contained shared-field model fitting; no old project scripts are loaded
# [ORIGINAL] Original model choices; [EXTENSION] corrected projectors/separate noise;

# Corrected reproduction of the Melding, PSgeo, and PSmelding models ------------
#
# This file provides functions for fitting three spatial models with INLA-SPDE:
#
#   Melding   : Gaussian point measurements + Gaussian areal measurements
#   PSgeo     : Gaussian point measurements + point-location process
#   PSmelding : Both measurement types + point-location process
#
# All three models use a common latent Matérn spatial field. The preferential-
# sampling models additionally use that field to explain sampling locations
#
# Important implementation choices:
#   * Preserve the original constrained alpha = 2 Matérn SPDE and INLA priors
#   * Estimate the observation-noise precision(s), rather than fixing them
#   * Use separate Gaussian likelihoods for point and areal measurements
#   * Correct the areal projector and point-process integration weights
#   * Save the complete INLA return, likelihood mapping, and fit diagnostics
#
# This file only defines functions. Sourcing it does not attach packages, change
# global options, or fit a model. configure_runtime() checks packages once per run
#


# Method references (shared-field simulation model):
# Zhong et al. (2025), Sections 2.3, 3.3-3.4, equations (3), (4), (10), and
# Figure 3: https://doi.org/10.1093/jrsssa/qnae058
# Simpson et al. (2016), equations (3)-(4), Section 5 and Figure 1:
# https://doi.org/10.1093/biomet/asv064
# Moraga et al. (2017), Section 2.2, p. 30:
# https://doi.org/10.1016/j.spasta.2017.04.006
# These numerical methods do not determine the priors or the simulation design

# 1. Project the spatial field onto areal observations -------------------------

# Construct the matrix that maps mesh-vertex values to areal averages
#
# [PAPER] Moraga Section 2.2 and Zhong Section 3.4/Figure 3
# For area B_j, the approximation used by Moraga et al. (2017) is
#
#     mean_field(B_j) ~= (1 / H_j) * sum_{k in V_j} field(vertex_k),
#
# where V_j is the set of mesh vertices intersecting B_j and H_j = |V_j|
# Thus row j of the projector contains 1/H_j at those mesh vertices and zero
# elsewhere. This is an equal-vertex-weight approximation to a block average,
# NOT an exact area integral; its accuracy depends on mesh resolution
# This approximates |B_j|^(-1) integral_B_j phi(s) ds, consistent with Zhong
# equation (3). Do not multiply by |B_j|: the response is an average, not a total
# In particular, these 1/H_j weights are distinct from the LGCP dual-cell areas
#
# Area membership is calculated independently for each polygon. In particular,
# a vertex lying on a shared boundary can contribute to multiple adjacent areas
#
# Returns:
#   A       : sparse (number of areas) x (number of mesh vertices) projector
#   members : list of mesh-vertex indices belonging to each area
#   counts  : number of vertices used in each areal average
# [EXTENSION] Each area independently retains vertices on shared polygon boundaries
build_areal_projector <- function(mesh, boundaryregion, dearea) {
  # Put the mesh vertices into an sf object in the same coordinate system as
  # the study region, so polygon/point intersections are geometrically valid
  mesh_vertices <- sf::st_as_sf(
    data.frame(x = mesh$loc[, 1], y = mesh$loc[, 2]),
    coords = c("x", "y"),
    crs = sf::st_crs(boundaryregion)
  )

  # Each list element contains *all* vertices intersecting one areal polygon
  # st_make_valid() handles invalid polygon geometries before intersection
  area_members <- sf::st_intersects(
    sf::st_make_valid(dearea),
    mesh_vertices
  )
  vertex_counts <- lengths(area_members)

  if (any(vertex_counts == 0L)) {
    stop("An areal unit contains no mesh vertices; refine the mesh.")
  }

  # Sparse matrix entries: row = area, column = vertex, value = 1/H_j
  # The resulting rows represent averages, so each must sum to one
  area_projection <- Matrix::sparseMatrix(
    i = rep(seq_along(area_members), vertex_counts),
    j = unlist(area_members, use.names = FALSE),
    x = rep(1 / vertex_counts, vertex_counts),
    dims = c(nrow(dearea), mesh$n)
  )

  # Positive vertex_counts and the 1/H_j construction ensure unit row sums
  # The fitted bundle also saves these row sums for inspection

  list(
    A = area_projection,
    members = area_members,
    counts = vertex_counts
  )
}


# 2. Construct integration weights for the point-location process -------------

# [PAPER] Simpson equations (3)-(4), Section 5; Zhong Section 3.3/equation (10)
# This is the mesh-based rule of Simpson et al., NOT classical Simpson's 1/3 rule
# Approximate the spatial integral in the log-Gaussian Cox process likelihood:
#
#   log L = sum_{observed points s_i} log(lambda(s_i))
#           - integral_D lambda(s) ds  + constant,
#
# using one quadrature weight w_k per mesh vertex. The approximation is
#
#   integral_D lambda(s) ds ~= sum_k w_k * lambda(vertex_k)
#
# Each mesh triangle is divided into three quadrilaterals, one associated with
# each triangle vertex. A quadrilateral is bounded by that vertex, the midpoints
# of its two incident edges, and the triangle centroid. After clipping all these
# pieces to the study region D, their areas are summed by associated vertex
# These are vertex-centered (barycentric dual-cell) quadrature weights
# Each piece is exactly the part of its triangle where that vertex's linear
# basis function is maximal. Integration nodes are the MESH VERTICES; triangle
# centroids define the cell boundaries, not the intensity evaluation locations
# Thus w_k = |V_k intersect D|, unnormalised, and sum(w_k) = |D|
#
# Vertices outside D can receive zero weight. They are excluded from the Poisson
# integration rows, while the original mesh indexing is retained throughout
#
# This is a self-contained adaptation of the bounded numerical construction in
# R/simulation.R; that script is deliberately not sourced here
#
# Returns:
#   weights     : integration weight for every mesh vertex, including zeros
#   active      : indices of vertices with positive integration weight
#   A           : sparse selector matrix for the active mesh vertices
#   domain_area : geometric area of the observation window D
# [EXTENSION] Clip dual cells to the observation window and omit zero-exposure integration rows
build_point_process_quadrature <- function(mesh, boundaryregion) {
  mesh_triangles <- mesh$graph$tv
  number_of_pieces <- 3L * nrow(mesh_triangles)

  dual_pieces <- vector("list", number_of_pieces)
  piece_vertex_id <- integer(number_of_pieces)

  # Build three pieces per triangle. The pieces collectively cover each
  # triangle; grouping pieces with the same vertex produces its dual cell
  for (triangle_index in seq_len(nrow(mesh_triangles))) {
    triangle_vertices <- mesh$loc[
      mesh_triangles[triangle_index, ],
      1:2,
      drop = FALSE
    ]
    triangle_centre <- colMeans(triangle_vertices)

    for (corner_index in 1:3) {
      piece_index <- 3L * (triangle_index - 1L) + corner_index
      next_corner <- corner_index %% 3L + 1L
      previous_corner <- (corner_index + 1L) %% 3L + 1L

      # Repeat the first coordinate to close the sf polygon ring
      piece_coordinates <- rbind(
        triangle_vertices[corner_index, ],
        (triangle_vertices[corner_index, ] +
            triangle_vertices[next_corner, ]) / 2,
        triangle_centre,
        (triangle_vertices[corner_index, ] +
            triangle_vertices[previous_corner, ]) / 2,
        triangle_vertices[corner_index, ]
      )

      dual_pieces[[piece_index]] <- sf::st_polygon(
        list(piece_coordinates)
      )
      piece_vertex_id[piece_index] <- mesh_triangles[
        triangle_index, corner_index
      ]
    }
  }

  # Union the observation-window polygons before clipping, so overlapping
  # boundary components do not cause the same quadrature area to be counted
  # more than once
  observation_domain <- sf::st_sf(
    geometry = sf::st_union(
      sf::st_geometry(sf::st_make_valid(boundaryregion))
    )
  )

  dual_geometry <- sf::st_sf(
    vertex_id = piece_vertex_id,
    geometry = sf::st_sfc(
      dual_pieces,
      crs = sf::st_crs(boundaryregion)
    )
  )

  clipped_pieces <- suppressWarnings(
    sf::st_intersection(dual_geometry, observation_domain)
  )
  clipped_areas <- as.numeric(sf::st_area(clipped_pieces))

  if (any(!is.finite(clipped_areas)) || any(clipped_areas < 0)) {
    stop("Dual-cell clipping produced invalid areas.")
  }

  # Sum the positive areas assigned to each vertex. Repeated row/column
  # indices in sparseMatrix() are combined, which aggregates multiple pieces
  # from neighbouring triangles into a single vertex-level weight
  has_positive_area <- clipped_areas > 0
  integration_weights <- as.numeric(
    Matrix::sparseMatrix(
      i = clipped_pieces$vertex_id[has_positive_area],
      j = rep(1L, sum(has_positive_area)),
      x = clipped_areas[has_positive_area],
      dims = c(mesh$n, 1L)
    )
  )

  observation_domain_area <- as.numeric(
    sf::st_area(observation_domain)
  )

  # The dual pieces should form a partition of D, up to numerical tolerance
  # Checking their total weight helps detect missing or duplicated area
  if (length(observation_domain_area) != 1L ||
      !is.finite(observation_domain_area) ||
      observation_domain_area <= 0 ||
      abs(sum(integration_weights) - observation_domain_area) > 1e-8) {
    stop("Dual-cell integration weights do not cover the observation domain.")
  }

  active_vertices <- which(integration_weights > 0)

  list(
    weights = integration_weights,
    active = active_vertices,
    A = Matrix::Diagonal(mesh$n)[active_vertices, , drop = FALSE],
    domain_area = observation_domain_area
  )
}


# 3. Fit a single Melding, PSgeo, or PSmelding model -----------------------------

# Arguments:
#   model          : "Melding", "PSgeo", or "PSmelding"
#   depoint        : sf point observations with numeric column `value`
#   dearea         : sf areal observations with numeric column `value`;
#                    unused by PSgeo
#   dppoint        : sf locations where the response surface is predicted
#   boundaryregion : sf observation-window polygons in a projected, planar CRS
#   mesh           : INLA-SPDE triangulation covering the point observations
#                    and prediction locations
#   config         : settings controlling noise estimation, runtime, and seed
#   on_inla_return : optional callback to checkpoint the *unmodified* INLA fit
#
# The returned list retains the full INLA fit, prediction summaries, design
# matrices, diagnostics, and configuration
fit_spatial_model <- function(
  model,
  depoint,
  dearea,
  dppoint,
  boundaryregion,
  mesh,
  config = list(),
  on_inla_return = NULL
) {
  # 3.1. Validate configuration and spatial inputs -----------------------------

  if (length(model) != 1L ||
      !model %in% c("Melding", "PSgeo", "PSmelding")) {
    stop("model must be one of Melding, PSgeo, PSmelding.")
  }
  if (!is.list(config)) {
    stop("config must be a list.")
  }
  if (!is.null(on_inla_return) && !is.function(on_inla_return)) {
    stop("on_inla_return must be NULL or a function.")
  }

  # Reuse shared defaults; these extra fields describe this individual fit
  default_config <- c(
    default_simulation_config$fitting,
    list(
      tau2_s = default_simulation_config$generating$noise$tau2_s[1L],
      tau2_B = default_simulation_config$generating$noise$tau2_B[1L],
      seed = NULL
    )
  )
  config <- utils::modifyList(default_config, config, keep.null = TRUE)

  # Generating variances are recorded as metadata; they do not set the fitted
  # precisions. Their validity is checked once in validate_config()

  number_of_threads <- config$inla_threads
  if (length(number_of_threads) != 1L ||
      !is.finite(number_of_threads) ||
      number_of_threads < 1 ||
      number_of_threads != floor(number_of_threads)) {
    stop("config$inla_threads must be a positive integer.")
  }

  # configure_runtime() checks packages before the study starts. Explicit
  # package::function calls also fail clearly if this helper is used alone

  # These flags determine which data types and likelihood components are used
  # PSgeo omits areal data; Melding does not model the point-location process
  use_areal_data <- model != "PSgeo"
  use_point_pattern <- model != "Melding"

  if (is.null(mesh) || !inherits(dppoint, "sf") || !nrow(dppoint)) {
    stop("Supply a mesh and sf prediction points.")
  }
  if (!inherits(depoint, "sf") || !nrow(depoint) ||
      !"value" %in% names(depoint) || any(!is.finite(depoint$value))) {
    stop("Point observations must be sf points with finite value entries.")
  }
  if (use_areal_data &&
      (!inherits(dearea, "sf") || !nrow(dearea) ||
        !"value" %in% names(dearea) || any(!is.finite(dearea$value)))) {
    stop("Areal observations must be sf polygons with finite value entries.")
  }
  if (!inherits(boundaryregion, "sf") ||
      isTRUE(sf::st_is_longlat(boundaryregion))) {
    stop("The boundary must be sf geometry in planar coordinates.")
  }
  if (!isTRUE(sf::st_crs(depoint) == sf::st_crs(boundaryregion)) ||
      !isTRUE(sf::st_crs(dppoint) == sf::st_crs(boundaryregion)) ||
      (use_areal_data &&
        !isTRUE(sf::st_crs(dearea) == sf::st_crs(boundaryregion)))) {
    stop("Point, area, prediction and boundary coordinates must share a CRS.")
  }

  # 3.2. Protect the caller's random-number stream -----------------------------

  # An optional seed makes a fit reproducible without permanently changing the
  # caller's .Random.seed. on.exit() restores even an initially absent seed
  if (!is.null(config$seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous_seed <- if (had_seed) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }

    on.exit(
      {
        if (had_seed) {
          assign(".Random.seed", previous_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv,
          inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )

    set.seed(config$seed)
  }

  # Record wall-clock, user CPU, and system CPU time as named quantities
  summarize_runtime <- function(time_difference) {
    c(
      elapsed = unname(time_difference[["elapsed"]]),
      user = unname(
        time_difference[["user.self"]] + time_difference[["user.child"]]
      ),
      system = unname(
        time_difference[["sys.self"]] + time_difference[["sys.child"]]
      )
    )
  }
  fit_started <- proc.time()

  # 3.3. Set up the Matérn SPDE and spatial projectors --------------------------

  # Retain the repository's constrained Matérn field (alpha = 2). No PC priors
  # are introduced
  spde <- INLA::inla.spde2.matern(
    mesh = mesh,
    alpha = 2,
    constr = TRUE
  )

  # The preferential models use 'i' for the spatial field so that 'j' can be
  # its INLA copy in the point-pattern likelihood. Melding uses the name 's'
  spatial_effect_name <- if (use_point_pattern) {
    "i"
  } else {
    "s"
  }

  # Standard SPDE point projectors interpolate the latent mesh weights using
  # piecewise-linear (barycentric) basis functions at each requested location
  A_point <- INLA::inla.spde.make.A(
    mesh = mesh,
    loc = sf::st_coordinates(depoint)[, 1:2, drop = FALSE]
  )
  A_prediction <- INLA::inla.spde.make.A(
    mesh = mesh,
    loc = sf::st_coordinates(dppoint)[, 1:2, drop = FALSE]
  )

  if (any(abs(Matrix::rowSums(A_point) - 1) > 1e-7) ||
      any(abs(Matrix::rowSums(A_prediction) - 1) > 1e-7)) {
    stop("The mesh must cover all observation and prediction locations.")
  }

  areal_projector <- if (use_areal_data) {
    build_areal_projector(mesh, boundaryregion, dearea)
  } else {
    NULL
  }
  point_pattern_quadrature <- if (use_point_pattern) {
    build_point_process_quadrature(mesh, boundaryregion)
  } else {
    NULL
  }

  # 3.4. Specify likelihood families and their response columns ----------------

  # Point and areal responses have separate estimated Gaussian precisions
  # PSgeo uses only the point-response precision
  # The sampling-location process, when present, receives a Poisson family
  number_of_gaussian_families <- if (use_areal_data) {
    2L
  } else {
    1L
  }
  likelihood_families <- rep("gaussian", number_of_gaussian_families)
  if (use_point_pattern) {
    likelihood_families <- c(likelihood_families, "poisson")
  }

  # Column indices refer to INLA's multiple-likelihood response matrix
  point_column <- 1L
  area_column <- if (use_areal_data) {
    number_of_gaussian_families
  } else {
    NA_integer_
  }
  pattern_column <- if (use_point_pattern) {
    length(likelihood_families)
  } else {
    NA_integer_
  }

  # A row contributes to exactly one likelihood: the observed response goes
  # into that family's column, and the other columns contain NA. For a single
  # family, INLA accepts a response vector instead of a one-column matrix
  make_response <- function(values, column) {
    if (length(likelihood_families) == 1L) {
      return(values)
    }

    response <- matrix(
      NA_real_,
      nrow = length(values),
      ncol = length(likelihood_families)
    )
    response[, column] <- values
    response
  }

  # Reusable stack constructor for Gaussian measurement and prediction rows
  # 'b0' is the response intercept; the named spatial effect is represented by
  # mesh indices and linked to the observations through their projector A
  make_gaussian_stack <- function(tag, values, projector, column) {
    INLA::inla.stack(
      tag = tag,
      data = list(
        y = make_response(values, column),
        e = rep(NA_real_, nrow(projector))
      ),
      A = list(1, projector),
      effects = list(
        data.frame(b0 = rep(1, nrow(projector))),
        stats::setNames(list(seq_len(mesh$n)), spatial_effect_name)
      )
    )
  }

  # 3.5. Assemble point, point-pattern, areal, and prediction stacks ------------

  # The tags are part of the existing interface; the 'pred1' tag is used later
  # to recover rows in the fitted-value summaries corresponding to dppoint
  component_stacks <- list(
    make_gaussian_stack("est1", depoint$value, A_point, point_column)
  )
  A_pattern <- NULL

  if (use_point_pattern) {
    # A_pattern evaluates the copied field at (a) integration mesh vertices
    # and (b) the *observed* sampling locations, in that order
    A_pattern <- rbind(point_pattern_quadrature$A, A_point)
    number_of_integration_rows <- length(point_pattern_quadrature$active)

    # [PAPER] Simpson equation (4), as used by Zhong equation (10):
    #   log L = sum_i eta(s_i) - sum_k w_k exp(eta(v_k)) + constant,
    # with eta = alpha + gamma*phi. The event term uses only observed point
    # locations; areal centroids and prediction sites are not point events
    # INLA uses the kernel y*eta - E*exp(eta) for this E=0 pseudo-data device
    # It is not an ordinary dpois(1, lambda=0) observation
    # Standard Poisson pseudo-data representation of the point-process term:
    #   * Quadrature rows: count = 0, exposure = integration weight w_k
    #   * Observed locations: count = 1, exposure = 0
    # Together these represent the negative integral and positive point-sum
    # terms of the log-Gaussian Cox process likelihood, respectively
    pattern_response <- c(
      rep(0, number_of_integration_rows),
      rep(1, nrow(depoint))
    )
    pattern_exposure <- c(
      point_pattern_quadrature$weights[point_pattern_quadrature$active],
      rep(0, nrow(depoint))
    )

    component_stacks <- c(
      component_stacks,
      list(
        INLA::inla.stack(
          tag = "est1pp",
          data = list(
            y = make_response(pattern_response, pattern_column),
            e = pattern_exposure
          ),
          A = list(1, A_pattern),
          effects = list(
            data.frame(b0pp = rep(1, nrow(A_pattern))),
            j = seq_len(mesh$n)
          )
        )
      )
    )
  }

  if (use_areal_data) {
    component_stacks <- c(
      component_stacks,
      list(
        make_gaussian_stack(
          "est2", dearea$value, areal_projector$A, area_column
        )
      )
    )
  }

  # Missing responses request fitted latent-surface values at prediction sites;
  # they do not add observations to the fitted likelihood
  component_stacks <- c(
    component_stacks,
    list(
      make_gaussian_stack(
        "pred1", rep(NA_real_, nrow(dppoint)), A_prediction, point_column
      )
    )
  )

  combined_stack <- do.call(INLA::inla.stack, component_stacks)
  stack_data <- INLA::inla.stack.data(combined_stack)

  # 3.6. Specify the latent model and INLA arguments ---------------------------

  # Response surface:              b0 + phi(s)
  # Preferential point intensity:  b0pp + gamma * phi(s)
  #
  # f(j, copy = "i", fixed = FALSE) uses the same underlying spatial field
  # as the response, with an estimated multiplicative copy coefficient gamma
  # [PAPER] Zhong equations (3)-(4): both Gaussian supports share b0 and the
  # SAME response field. Only the point-process likelihood uses its gamma copy
  # No separate areal field or estimated areal scaling coefficient is added in
  # this shared-field simulation model (the paper's application is different)
  # There are no extra spatial fields or covariate terms in this reproduction
  # The estimated copy coefficient gamma uses INLA's default prior
  model_formula <- if (use_point_pattern) {
    y ~ 0 + b0 + b0pp + f(i, model = spde) + f(j, copy = "i", fixed = FALSE)
  } else {
    y ~ 0 + b0 + f(s, model = spde)
  }

  # Store the SPDE object in the formula's environment. This makes the model
  # formula self-contained for the INLA call and its saved configuration
  formula_environment <- new.env(parent = asNamespace("INLA"))
  formula_environment$spde <- spde
  environment(model_formula) <- formula_environment

  # Apply the same variance adjustment to every model for comparable numerics
  inla_controls <- list(control.vb = list(strategy = "variance"))

  posterior_quantiles <- c(
    0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5,
    0.75, 0.9, 0.95, 0.975, 0.99, 0.995
  )

  # Keep the full latent configurations and marginal predictor distributions
  # Those quantities are required for later joint posterior simulations/scores,
  # whereas summary.fitted.values alone only provides marginal summaries
  inla_arguments <- list(
    formula = model_formula,
    family = likelihood_families,
    data = stack_data,
    quantiles = posterior_quantiles,
    control.compute = list(
      config = TRUE,
      return.marginals = TRUE,
      return.marginals.predictor = TRUE
    ),
    control.predictor = list(
      compute = TRUE,
      link = 1,
      A = INLA::inla.stack.A(combined_stack)
    ),
    control.inla = inla_controls,
    num.threads = paste(number_of_threads, number_of_threads, sep = ":"),
    verbose = isTRUE(config$verbose),
    safe = FALSE,  # Make one fitting attempt, rather than an automatic retry
    .parent.frame = formula_environment
  )

  # Only Poisson rows use these exposure values. Gaussian stacks carry NA
  if (use_point_pattern) {
    inla_arguments$E <- stack_data$e
  }

  # 3.7. Document the expected likelihood-to-precision mapping -----------------

  # Store this mapping *before* fitting, rather than assuming that the first
  # hyperparameter returned by INLA always belongs to the point observations
  # The first Gaussian precision belongs to points; when areas are included,
  # the second Gaussian precision belongs to areas
  expected_precision_names <- c(
    "Precision for the Gaussian observations",
    if (number_of_gaussian_families == 2L) {
      "Precision for the Gaussian observations[2]"
    }
  )

  likelihood_map <- data.frame(
    component = c(
      "point",
      if (use_areal_data) {
        "area"
      },
      if (use_point_pattern) {
        "point_pattern"
      }
    ),
    family = c(
      "gaussian",
      if (use_areal_data) {
        "gaussian"
      },
      if (use_point_pattern) {
        "poisson"
      }
    ),
    column = c(
      point_column,
      if (use_areal_data) {
        area_column
      },
      if (use_point_pattern) {
        pattern_column
      }
    ),
    precision_name = c(
      expected_precision_names[point_column],
      if (use_areal_data) {
        expected_precision_names[area_column]
      },
      if (use_point_pattern) {
        NA_character_
      }
    ),
    fixed = FALSE,
    fixed_variance = NA_real_,
    stringsAsFactors = FALSE
  )

  # 3.8. Fit once ------------------------------------------------------------

  # Collect warnings without printing them during a large simulation study
  # If INLA fails, retain the original condition inside a classified error
  fit_warnings <- character()
  inla_started <- proc.time()

  fit <- withCallingHandlers(
    tryCatch(
      do.call(INLA::inla, inla_arguments),
      error = function(error_condition) {
        if (inherits(error_condition, "interrupt")) {
          stop(error_condition)
        }

        stop(structure(
          list(
            message = conditionMessage(error_condition),
            call = conditionCall(error_condition),
            parent = error_condition,
            origin = "INLA",
            warnings = unique(fit_warnings),
            diagnostics = list(
              stage = "INLA",
              runtime = summarize_runtime(proc.time() - fit_started),
              inla_runtime = summarize_runtime(proc.time() - inla_started)
            )
          ),
          class = unique(c(
            "psmelding_inla_error",
            class(error_condition)
          ))
        ))
      }
    ),
    warning = function(warning_condition) {
      fit_warnings <<- c(
        fit_warnings,
        conditionMessage(warning_condition)
      )
      invokeRestart("muffleWarning")
    }
  )

  inla_runtime <- summarize_runtime(proc.time() - inla_started)
  runtime <- summarize_runtime(proc.time() - fit_started)

  # 3.9. Checkpoint the raw fit before inspecting its output --------------------

  # The callback sees the complete, unmodified INLA return and enough context
  # to reconstruct the fit. It runs immediately after INLA returns
  # A failed checkpoint is reported explicitly to avoid silently losing a fit
  if (!is.null(on_inla_return)) {
    tryCatch(
      on_inla_return(
        fit = fit,
        context = list(
          model = model,
          inputs = list(
            depoint = depoint,
            dearea = if (use_areal_data) {
              dearea
            } else {
              NULL
            },
            dppoint = dppoint,
            boundaryregion = boundaryregion,
            mesh = mesh
          ),
          stack = combined_stack,
          spde = spde,
          spde_field_name = spatial_effect_name,
          projection_matrices = list(
            Ae1 = A_point,
            Ae2 = if (use_areal_data) {
              areal_projector$A
            } else {
              NULL
            },
            Ap1 = A_prediction,
            A.pp = A_pattern
          ),
          integration_weights = if (use_point_pattern) {
            point_pattern_quadrature$weights
          } else {
            NULL
          },
          inla_arguments = inla_arguments,
          config = config,
          warnings = unique(fit_warnings),
          runtime = runtime,
          inla_runtime = inla_runtime,
          likelihood_map_expected = likelihood_map
        )
      ),
      error = function(error_condition) {
        stop(structure(
          list(
            message = paste(
              "Could not preserve the original INLA return:",
              conditionMessage(error_condition)
            ),
            call = NULL,
            parent = error_condition
          ),
          class = c("psmelding_checkpoint_error", "error", "condition")
        ))
      }
    )
  }

  # 3.10. Validate the return and locate prediction rows -----------------------

  # INLA's stack may reorder rows when combining stacks. Use its stored tag
  # indices instead of assuming prediction rows occur at the end of the fit
  prediction_index <- INLA::inla.stack.index(
    combined_stack,
    tag = "pred1"
  )$data

  # Missing joint configurations affect joint scoring, but should not discard
  # available marginal predictions; record the issue as a diagnostic instead
  has_joint_configurations <- !is.null(fit$misc$configs) &&
    length(fit$misc$configs) > 0L

  if (length(prediction_index) != nrow(dppoint)) {
    stop("Prediction index does not align with dppoint.")
  }

  returned_precision_names <- grep(
    "^Precision for the Gaussian observations",
    rownames(fit$summary.hyperpar),
    value = TRUE
  )
  likelihood_mapping_valid <- setequal(
    returned_precision_names,
    expected_precision_names
  )

  # Do not invent or substitute precision labels when INLA returns unexpected
  # names. Preserve both lists and flag the mismatch for downstream scoring

  # 3.11. Extract posterior predictions at dppoint ------------------------------

  # These are posterior summaries of the fitted latent response surface at
  # prediction points, not draws of new observations with measurement noise
  prediction_summary <- fit$summary.fitted.values[
    prediction_index,
    ,
    drop = FALSE
  ]

  predictions <- dppoint
  predictions$pred_mean <- prediction_summary[, "mean"]
  predictions$pred_med <- prediction_summary[, "0.5quant"]
  predictions$sd <- prediction_summary[, "sd"]
  predictions$pred_ll <- prediction_summary[, "0.025quant"]
  predictions$pred_ul <- prediction_summary[, "0.975quant"]

  # Retain the complete fit and the objects needed for scoring and diagnosis
  list(
    model = model,
    fit = fit,
    predictions = predictions,
    prediction_index = prediction_index,
    prediction_coordinates = sf::st_coordinates(dppoint)[
      , 1:2, drop = FALSE
    ],
    stack = combined_stack,
    spde = spde,
    spde_field_name = spatial_effect_name,
    projection_matrices = list(
      Ae1 = A_point,
      Ae2 = if (use_areal_data) {
        areal_projector$A
      } else {
        NULL
      },
      Ap1 = A_prediction,
      A.pp = A_pattern
    ),
    integration_weights = if (use_point_pattern) {
      point_pattern_quadrature$weights
    } else {
      NULL
    },
    formula = model_formula,
    inla_arguments = inla_arguments,
    likelihood_map = likelihood_map,
    runtime = runtime,
    inla_runtime = inla_runtime,
    warnings = unique(fit_warnings),

    diagnostics = list(
      has_joint_configurations = has_joint_configurations,
      likelihood_mapping_valid = likelihood_mapping_valid,
      expected_precision_names = expected_precision_names,
      returned_precision_names = returned_precision_names,
      has_predictor_marginals = length(fit$marginals.linear.predictor) > 0L,
      finite_prediction_summaries = all(is.finite(as.matrix(prediction_summary))),
      mode_status = fit$mode$mode.status,
      inla_cpu_used = fit$cpu.used,
      area_vertex_counts = if (use_areal_data) {
        areal_projector$counts
      } else {
        NULL
      },
      area_projector_row_sums = if (use_areal_data) {
        Matrix::rowSums(areal_projector$A)
      } else {
        NULL
      },
      integration_weight_sum = if (use_point_pattern) {
        sum(point_pattern_quadrature$weights)
      } else {
        NULL
      },
      integration_domain_area = if (use_point_pattern) {
        point_pattern_quadrature$domain_area
      } else {
        NULL
      },
      active_integration_vertices = if (use_point_pattern) {
        point_pattern_quadrature$active
      } else {
        NULL
      }
    ),

    config = config
  )
}


# 4. Public study interface ---------------------------------------------------
# build_mesh() creates the triangulation once per study run. fit_model() adapts
# a generated dataset and study config to the single-model fitting function
# Keep mesh controls separate from the grid used to generate the true surface
build_mesh <- function(config) {
  generating <- config$generating
  boundary <- rbind(
    c(generating$x_limits[1], generating$y_limits[1]),
    c(generating$x_limits[2], generating$y_limits[1]),
    c(generating$x_limits[2], generating$y_limits[2]),
    c(generating$x_limits[1], generating$y_limits[2]),
    c(generating$x_limits[1], generating$y_limits[1])
  )
  INLA::inla.mesh.2d(
    loc.domain = boundary,
    offset = config$mesh$offset,
    max.edge = config$mesh$max_edge,
    cutoff = config$mesh$cutoff
  )
}

fit_model <- function(model, dataset, mesh, config, on_inla_return = NULL) {
  fitting <- config$fitting
  fitting$tau2_s <- dataset$scenario$tau2_s
  fitting$tau2_B <- dataset$scenario$tau2_B
  job <- dataset$scenario
  fitting$seed <- sampling_seed(config, job, "fit", model)
  fit_spatial_model(
    model = model,
    depoint = dataset$points,
    dearea = dataset$areas,
    dppoint = dataset$prediction_points,
    boundaryregion = dataset$boundary,
    mesh = mesh,
    config = fitting,
    on_inla_return = on_inla_return
  )
}
