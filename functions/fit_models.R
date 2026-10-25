# Create spatial mesh
create_mesh <- function(max_edge = c(0.05, 0.25)) {
  # Boundary of D = [0,1] x [0,1]
  boundary = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))
  
  # Construct triangular mesh
  # Allow max.edge to vary for mesh sensitivity
  # Fix offset and cutoff specifications
  mesh = INLA::inla.mesh.2d(loc.domain = boundary,
                            max.edge = max_edge,
                            offset = c(0.1, 0.35),
                            cutoff = 0.01)
  
  return(mesh)
}



# Create SPDE model
create_spde <- function(mesh, prior_range = NULL, prior_sigma = NULL) {
  # Require both PC priors or neither
  if (is.null(prior_range) != is.null(prior_sigma)) {
    stop("Provide both prior_range and prior_sigma, or neither.")
  }
  
  if (is.null(prior_range)) {
    # Matern SPDE with integrate-to-zero constraint
    spde = INLA::inla.spde2.matern(mesh = mesh,
                                   alpha = 2,
                                   constr = TRUE)
  } else {
    # Matern SPDE with penalized complexity priors
    # Allow for prior sensitivity
    spde = INLA::inla.spde2.pcmatern(mesh = mesh,
                                     alpha = 2,
                                     prior.range = prior_range,
                                     prior.sigma = prior_sigma,
                                     constr = TRUE)
  }
  
  return(spde)
}



# Create projection matrices
create_projections <- function(data, mesh, areal_method = c("vertices", "grid")) {
  # Allow previous method name
  if (identical(areal_method, "mesh_vertices")) {
    areal_method = "vertices"
  }
  areal_method = match.arg(areal_method)
  
  # Observed point locations
  points = data$points[, c("x", "y")]
  
  # Prediction grid
  grid = data$prediction[, c("x", "y")]
  
  # Point observation projection matrix
  A_points = INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(points))
  
  # Areal observation projection matrix
  J = nrow(data$areas)
  side = as.integer(sqrt(J))
  
  if (side^2 != J) {
    stop("The number of areas must be a perfect square.")
  }
  
  if (areal_method == "grid") {
    # Assign prediction grid cells to areas
    ix = pmin(side, floor(grid$x * side) + 1L)
    iy = pmin(side, floor(grid$y * side) + 1L)
    
    area_id = (iy - 1L) * side + ix
    
    # Average over prediction grid cells
    A_areas = INLA::inla.spde.make.A(mesh = mesh,
                                     loc = as.matrix(grid),
                                     block = area_id,
                                     n.block = J,
                                     block.rescale = "count")
  } else {
    # Identify mesh vertices within each area
    area_locations = lapply(seq_len(J), function(j) {
      area = data$areas[j, ]
      
      inside = (mesh$loc[, 1] >= area$xmin & mesh$loc[, 1] <= area$xmax &
                  mesh$loc[, 2] >= area$ymin & mesh$loc[, 2] <= area$ymax)
      
      if (!any(inside)) {
        stop("Area ", j, " contains no mesh vertices.")
      }
      
      mesh$loc[inside, 1:2, drop = FALSE]
    })
    
    # Assign selected mesh vertices to areas
    area_id = rep(seq_len(J),
                  times = vapply(area_locations, nrow, integer(1)))
    
    # Average over mesh vertices
    A_areas = INLA::inla.spde.make.A(mesh = mesh,
                                     loc = do.call(rbind, area_locations),
                                     block = area_id,
                                     n.block = J,
                                     block.rescale = "sum")
  }
  
  # Prediction projection matrix
  A_prediction = INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(grid))
  
  return(list(points = A_points,
              areas = A_areas,
              prediction = A_prediction))
}



# Create Gaussian observation stack
create_gaussian_stack <- function(data, projections, spde,
                                  include_areas = TRUE, preferential = FALSE,
                                  separate_noise = FALSE) {
  # Format response for one, two, or three likelihoods
  make_response <- function(values, tag) {
    n = length(values)
    
    if (separate_noise && include_areas) {
      if (tag == "points") {
        y = cbind(values, rep(NA_real_, n))
      } else if (tag == "areas") {
        y = cbind(rep(NA_real_, n), values)
      } else {
        y = matrix(NA_real_, nrow = n, ncol = 2L)
      }
    } else {
      y = matrix(values, ncol = 1L)
    }
    
    if (preferential) {
      y = cbind(y, rep(NA_real_, n))
    }
    
    if (ncol(y) == 1L) {
      return(as.numeric(y[, 1]))
    }
    
    return(y)
  }
  
  # Construct one Gaussian data stack
  make_stack <- function(values, A, tag) {
    n = length(values)
    
    stack = INLA::inla.stack(
      tag = tag,
      data = list(
        y = make_response(values, tag),
        exposure = rep(NA_real_, n)
      ),
      A = list(1, A),
      effects = list(
        data.frame(
          b0 = rep(1, n),
          b0pp = rep(0, n)
        ),
        i = seq_len(spde$n.spde)
      )
    )
    
    return(stack)
  }
  
  # Point observations
  stacks = list(
    make_stack(values = data$points$value,
               A = projections$points,
               tag = "points")
  )
  
  # Areal observations
  if (include_areas) {
    stacks[[length(stacks) + 1L]] = make_stack(
      values = data$areas$value,
      A = projections$areas,
      tag = "areas"
    )
  }
  
  # Prediction locations
  stacks[[length(stacks) + 1L]] = make_stack(
    values = rep(NA_real_, nrow(data$prediction)),
    A = projections$prediction,
    tag = "prediction"
  )
  
  # Combine stacks
  stack = do.call(INLA::inla.stack, stacks)
  
  return(stack)
}



# Create mesh-dual integration weights
create_mesh_dual_weights <- function(mesh) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("The sf package is required for mesh-dual integration.")
  }
  
  if (mesh$manifold != "R2") {
    stop("Mesh-dual integration requires a two-dimensional mesh.")
  }
  
  # Mesh vertices and triangles
  vertices = mesh$loc[, 1:2, drop = FALSE]
  triangles = mesh$graph$tv
  
  # Divide each triangle into three dual-cell pieces
  polygons = vector("list", 3L * nrow(triangles))
  vertex_id = integer(length(polygons))
  index = 0L
  
  for (t in seq_len(nrow(triangles))) {
    nodes = triangles[t, ]
    xy = vertices[nodes, , drop = FALSE]
    centre = colMeans(xy)
    
    for (k in seq_len(3L)) {
      next_k = c(2L, 3L, 1L)[k]
      prev_k = c(3L, 1L, 2L)[k]
      
      polygon = rbind(
        xy[k, ],
        (xy[k, ] + xy[next_k, ]) / 2,
        centre,
        (xy[k, ] + xy[prev_k, ]) / 2,
        xy[k, ]
      )
      
      index = index + 1L
      polygons[[index]] = sf::st_polygon(list(polygon))
      vertex_id[index] = nodes[k]
    }
  }
  
  # Construct dual-cell pieces
  dual = sf::st_sf(
    vertex_id = vertex_id,
    geometry = sf::st_sfc(polygons)
  )
  
  # Study domain D = [0,1] x [0,1]
  boundary = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))
  
  domain = sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(boundary))
    )
  )
  
  # Intersect dual-cell pieces with the study domain
  clipped = suppressWarnings(
    sf::st_intersection(dual, domain)
  )
  
  # Sum clipped areas for each mesh vertex
  weights = numeric(mesh$n)
  
  if (nrow(clipped) > 0L) {
    area_sum = tapply(as.numeric(sf::st_area(clipped)),
                      clipped$vertex_id, sum)
    
    weights[as.integer(names(area_sum))] = as.numeric(area_sum)
  }
  
  # Check integration weights
  if (any(!is.finite(weights)) ||
      any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-6) {
    stop("Mesh-dual integration weights are invalid.")
  }
  
  return(weights)
}



# Create preferential sampling stack
create_point_process_stack <- function(
    data, projections, spde, mesh = NULL,
    separate_noise = FALSE,
    point_process_method = c("mesh_dual", "grid_counts"),
    integration_weights = NULL
) {
  point_process_method = match.arg(point_process_method)
  
  # Number of Gaussian likelihood columns
  n_gaussian = if (separate_noise) 2L else 1L
  
  if (point_process_method == "grid_counts") {
    # Prediction grid
    grid = data$prediction[, c("x", "y")]
    
    # Number of grid cells along each axis
    x_values = sort(unique(grid$x))
    y_values = sort(unique(grid$y))
    
    nx = length(x_values)
    ny = length(y_values)
    
    if (nx * ny != nrow(grid)) {
      stop("The prediction grid must be rectangular and complete.")
    }
    
    # Identify the grid cell containing each sampled point
    ix = pmin(nx, floor(data$points$x * nx) + 1L)
    iy = pmin(ny, floor(data$points$y * ny) + 1L)
    
    point_id = (iy - 1L) * nx + ix
    
    # Grid-cell IDs in the original row order
    grid_ix = match(grid$x, x_values)
    grid_iy = match(grid$y, y_values)
    
    grid_id = (grid_iy - 1L) * nx + grid_ix
    selected = match(point_id, grid_id)
    
    if (anyNA(selected)) {
      stop("Some sampled points could not be assigned to grid cells.")
    }
    
    # Number of sampled points in each grid cell
    y_pp = tabulate(selected, nbins = nrow(grid))
    
    # Cell areas for D = [0,1] x [0,1]
    exposure = rep(1 / nrow(grid), nrow(grid))
    
    # Projection matrix
    A_pp = projections$prediction
    
  } else {
    if (is.null(mesh)) {
      stop("A mesh is required for mesh-dual integration.")
    }
    
    # Compute integration weights unless supplied
    if (is.null(integration_weights)) {
      integration_weights = create_mesh_dual_weights(mesh)
    }
    
    if (spde$n.spde != mesh$n ||
        length(integration_weights) != mesh$n ||
        any(!is.finite(integration_weights)) ||
        any(integration_weights < 0) ||
        abs(sum(integration_weights) - 1) > 1e-6) {
      stop("Invalid mesh or mesh-dual integration weights.")
    }
    
    # Integration points at mesh vertices
    A_integration = Matrix::Diagonal(x = rep(1, mesh$n))
    
    # Combine integration points and observed point locations
    A_pp = rbind(A_integration, projections$points)
    
    # Zero counts at mesh vertices, one at each observed point
    y_pp = c(rep(0, mesh$n), rep(1, nrow(data$points)))
    
    # Integration weights and zero exposure for observed events
    exposure = c(integration_weights, rep(0, nrow(data$points)))
  }
  
  # Construct response matrix
  n = length(y_pp)
  
  y = matrix(NA_real_, nrow = n, ncol = n_gaussian + 1L)
  y[, n_gaussian + 1L] = y_pp
  
  # Poisson observation stack
  stack = INLA::inla.stack(
    tag = "sampling",
    data = list(
      y = y,
      exposure = exposure
    ),
    A = list(1, A_pp),
    effects = list(
      data.frame(
        b0 = rep(0, n),
        b0pp = rep(1, n)
      ),
      j = seq_len(spde$n.spde)
    )
  )
  
  return(stack)
}



# Fit one model
fit_model <- function(
    data,
    model = c("geo", "ps_geo",
              "melding", "ps_melding"),
    spde,
    projections,
    mesh = NULL,
    prior_gamma = c(1, 10),
    prior_noise = NULL,
    prior_noise_B = prior_noise,
    separate_noise = FALSE,
    point_process_method = c("mesh_dual", "grid_counts"),
    integration_weights = NULL,
    fixed_precision = 0.001,
    compute_criteria = FALSE
) {
  model = match.arg(model)
  point_process_method = match.arg(point_process_method)
  
  # Determine model components
  preferential = model %in% c("ps_geo", "ps_melding")
  include_areas = model %in% c("melding", "ps_melding")
  
  use_separate_noise = separate_noise && include_areas
  
  # Gaussian likelihood stack
  stack = create_gaussian_stack(
    data = data,
    projections = projections,
    spde = spde,
    include_areas = include_areas,
    preferential = preferential,
    separate_noise = use_separate_noise
  )
  
  if (preferential) {
    # Add preferential sampling likelihood
    sampling_stack = create_point_process_stack(
      data = data,
      projections = projections,
      spde = spde,
      mesh = mesh,
      separate_noise = use_separate_noise,
      point_process_method = point_process_method,
      integration_weights = integration_weights
    )
    
    stack = INLA::inla.stack(stack, sampling_stack)
    
    # Shared spatial field with estimated preferentiality
    # Normal prior parameters are mean and precision
    formula = y ~ 0 + b0 + b0pp +
      f(i, model = spde) +
      f(j, copy = "i", hyper = list(
        beta = list(
          fixed = FALSE,
          prior = "normal",
          param = prior_gamma
        )
      ))
    
  } else {
    # Non-preferential spatial model
    formula = y ~ 0 + b0 + f(i, model = spde)
  }
  
  # Specify likelihood families
  n_gaussian = if (use_separate_noise) 2L else 1L
  family = rep("gaussian", n_gaussian)
  
  if (preferential) {
    family = c(family, "poisson")
  }
  
  if (length(family) == 1L) {
    family = family[[1L]]
  }
  
  # Observation-error precision priors
  make_gaussian_prior <- function(prior) {
    if (is.null(prior)) {
      return(list())
    }
    
    return(list(
      hyper = list(
        prec = list(prior = "pc.prec", param = prior)
      )
    ))
  }
  
  family_prior = list(make_gaussian_prior(prior_noise))
  
  if (use_separate_noise) {
    family_prior = c(family_prior,
                     list(make_gaussian_prior(prior_noise_B)))
  }
  
  if (preferential) {
    family_prior = c(family_prior, list(list()))
  }
  
  if (length(family_prior) == 1L) {
    family_prior = family_prior[[1L]]
  }
  
  # Prepare INLA data
  stack_data = INLA::inla.stack.data(stack)
  
  # Fit model
  result = INLA::inla(
    formula = formula,
    family = family,
    data = stack_data,
    E = if (preferential) stack_data$exposure else NULL,
    control.family = family_prior,
    control.fixed = list(mean = 0, prec = fixed_precision),
    control.predictor = list(
      A = INLA::inla.stack.A(stack),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(
      dic = compute_criteria,
      waic = compute_criteria
    )
  )
  
  # Extract posterior predictions
  index = INLA::inla.stack.index(stack, tag = "prediction")$data
  summary_pred = result$summary.linear.predictor[index, ]
  
  prediction = data.frame(
    x = data$prediction$x,
    y = data$prediction$y,
    mean = summary_pred$mean,
    median = summary_pred[["0.5quant"]],
    sd = summary_pred$sd,
    q025 = summary_pred[["0.025quant"]],
    q975 = summary_pred[["0.975quant"]]
  )
  
  return(list(model = model,
              fit = result,
              prediction = prediction))
}



# Fit all selected models
fit_models <- function(
    data,
    models = c("geo", "ps_geo", "melding", "ps_melding"),
    mesh = NULL,
    max_edge = c(0.05, 0.25),
    prior_range = NULL,
    prior_sigma = NULL,
    prior_gamma = c(1, 10),
    prior_noise = NULL,
    prior_noise_B = prior_noise,
    separate_noise = FALSE,
    point_process_method = c("mesh_dual", "grid_counts"),
    areal_method = c("vertices", "grid"),
    fixed_precision = 0.001,
    compute_criteria = FALSE
) {
  models = match.arg(
    models,
    choices = c("geo", "ps_geo", "melding", "ps_melding"),
    several.ok = TRUE
  )
  
  point_process_method = match.arg(point_process_method)
  
  # Allow previous areal-method name
  if (identical(areal_method, "mesh_vertices")) {
    areal_method = "vertices"
  }
  areal_method = match.arg(areal_method)
  
  # Create mesh unless supplied
  if (is.null(mesh)) {
    mesh = create_mesh(max_edge = max_edge)
  }
  
  # Create SPDE model
  spde = create_spde(mesh = mesh,
                     prior_range = prior_range,
                     prior_sigma = prior_sigma)
  
  # Create projection matrices
  projections = create_projections(data = data,
                                   mesh = mesh,
                                   areal_method = areal_method)
  
  # Compute mesh-dual weights when required
  integration_weights = NULL
  
  if (point_process_method == "mesh_dual" &&
      any(models %in% c("ps_geo", "ps_melding"))) {
    integration_weights = create_mesh_dual_weights(mesh)
  }
  
  # Fit selected models
  fits = lapply(models, function(model) {
    fit_model(
      data = data,
      model = model,
      spde = spde,
      projections = projections,
      mesh = mesh,
      prior_gamma = prior_gamma,
      prior_noise = prior_noise,
      prior_noise_B = prior_noise_B,
      separate_noise = separate_noise,
      point_process_method = point_process_method,
      integration_weights = integration_weights,
      fixed_precision = fixed_precision,
      compute_criteria = compute_criteria
    )
  })
  
  names(fits) = models
  
  return(list(
    models = fits,
    mesh = mesh,
    spde = spde,
    point_process_method = point_process_method,
    areal_method = areal_method,
    integration_weights = integration_weights
  ))
}



# Fit models with logging 
# Fit each model separately so one failure does not stop the others
fit_models_logged <- function(data, models, mesh, replication, seed, ...) {
  attempts <- lapply(models, function(model) {
    message("  Fitting ", model)
    
    warnings <- character()
    messages <- character()
    error_message <- NA_character_
    
    # Fit model and capture R warnings, messages, and errors
    start <- Sys.time()
    
    result <- tryCatch(
      withCallingHandlers(
        fit_models(
          data = data,
          models = model,
          mesh = mesh,
          ...
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          messages <<- c(messages, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      ),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )
    
    elapsed_sec <- as.numeric(
      difftime(Sys.time(), start, units = "secs")
    )
    
    # Check numerical validity
    predictions_finite <- FALSE
    intervals_valid <- FALSE
    hyperparameters_finite <- FALSE
    gamma_present <- FALSE
    inla_sec <- NA_real_
    
    if (!is.null(result)) {
      fit <- result$models[[model]]$fit
      pred <- result$models[[model]]$prediction
      required <- c("mean", "median", "sd", "q025", "q975")
      
      predictions_finite <- !is.null(pred) &&
        nrow(pred) == nrow(data$prediction) &&
        all(required %in% names(pred)) &&
        all(is.finite(as.matrix(pred[, required, drop = FALSE])))
      
      intervals_valid <- predictions_finite &&
        all(pred$sd >= 0) &&
        all(pred$q025 <= pred$median) &&
        all(pred$median <= pred$q975)
      
      hyperparameters_finite <- !is.null(fit$summary.hyperpar) &&
        all(is.finite(as.matrix(fit$summary.hyperpar)))
      
      gamma_present <- "Beta for j" %in% rownames(fit$summary.hyperpar)
      
      if ("Total" %in% names(fit$cpu.used)) {
        inla_sec <- as.numeric(fit$cpu.used[["Total"]])
      }
    }
    
    # Preferential sampling models must have a gamma estimate
    gamma_valid <- !(model %in% c("ps_melding", "ps_geo")) || gamma_present
    
    # Classify fitting outcome
    status <- if (is.null(result)) {
      "error"
    } else if (!predictions_finite || !intervals_valid ||
               !hyperparameters_finite || !gamma_valid) {
      "invalid"
    } else if (length(warnings) > 0L) {
      "warning"
    } else {
      "ok"
    }
    
    # Record one log row per attempted model
    log <- data.frame(
      replication = replication,
      seed = seed,
      scenario = data$scenario,
      model = model,
      status = status,
      included = status %in% c("ok", "warning"),
      elapsed_sec = elapsed_sec,
      inla_sec = inla_sec,
      mesh_vertices = mesh$n,
      n_warnings = length(warnings),
      n_messages = length(messages),
      predictions_finite = predictions_finite,
      intervals_valid = intervals_valid,
      hyperparameters_finite = hyperparameters_finite,
      gamma_present = gamma_present,
      error_message = error_message,
      warning_text = paste(unique(warnings), collapse = " | "),
      message_text = paste(messages, collapse = " | ")
    )
    
    return(list(result = result, log = log))
  })
  
  names(attempts) <- models
  
  # Retain successful and warning-bearing fits for evaluation
  included <- vapply(attempts, function(x) x$log$included, logical(1))
  
  fits <- lapply(models[included], function(model) {
    attempts[[model]]$result$models[[model]]
  })
  
  names(fits) <- models[included]
  
  # Retain the SPDE specification for parameter recovery
  spde <- NULL
  if (any(included)) {
    first_model <- models[which(included)[1]]
    spde <- attempts[[first_model]]$result$spde
  }
  
  # Combine valid models for prediction evaluation
  results <- list(models = fits, mesh = mesh, spde = spde)
  
  fit_log <- do.call(rbind, lapply(attempts, function(x) x$log))
  rownames(fit_log) <- NULL
  
  return(list(results = results, log = fit_log))
}
