# ==================================================
# PSmelding Model Configuration
# ==================================================

# 1. Default simulation configuration ------------------------------------------

default_simulation_config <- list(
  package_library = Sys.getenv(
    "RESEARCH_COMPS_R_LIB",
    unset = "/Users/vinky/research-comps/.r-lib/4.2-arm64"
  ),
  models = c("Melding", "PSgeo", "PSmelding"),
  n_replicates = 100L,  
  base_seed = 2026L,  

  # Specifications for data generating parameters ------------------------------
  generating = list(
    # D = [0,1] X [0,1]
    x_limits = c(0, 1),
    y_limits = c(0, 1),
    grid_n = 50L,

    # Z(s) = mu + phi(s), where phi(s) has a Matern covariance:
    #   sigma2 : marginal variance of phi(s)
    #   nu     : Matern smoothness (fixed at 1 in this experiment)
    #   rho    : practical spatial correlation range
    mu = 0,
    sigma2 = 1,
    nu = 1,
    rho = 0.2,

    # Constant-preferentiality sampling model:
    #   log(lambda(s)) = alpha + gamma * phi(s)
    # when locations are sampled conditional on a fixed point count I, the
    # intensity is normalized to sampling probabilities
    # The common factor exp(alpha) then cancels
    # It does not cancel when the total Poisson count is itself random, 
    # so retain alpha as part of the model specification
    alpha = 0.05,

    # I = number of observed point responses
    # J = total number of areal units
    I = 100L,
    J = 4L,

    # Preferential sampling parameter
    gamma = 1,

    # NULL: use the constant-gamma scenarios above
    # A non-NULL will be used for extension
    sampling = NULL,

    # Independent Gaussian measurement-error variances for each data source:
    #   epsilon_s,i ~ N(0, tau2_s) at point locations
    #   epsilon_B,j ~ N(0, tau2_B) for areal averages
    noise = data.frame(
      tau2_s = 0.1,
      tau2_B = 0.1
    )
  ),

  # Specifications for INLA-SPDE -----------------------------------------------
  fitting = list(
    # Default priors
    #   inla.spde2.matern(mesh, alpha = 2, constr = TRUE)
    #    intercept: N(0, 1000)
    #    Gaussian observation precision: Gamma(1, 5e-5)
    #    estimated copy coefficient: N(1, 0.1)

    inla_threads = 1L,
    verbose = FALSE
  ),
  # offset   : mesh extension beyond the interior boundary
  # max_edge : maximum triangle-edge lengths
  # cutoff   : minimum spacing used to avoid overly close mesh vertices
  mesh = list(
    offset = c(0.1, 0.35), 
    max_edge = c(0.05, 0.25), 
    cutoff = 0.01 
  ),

  # Number of joint posterior draws used for latent-surface CRPS 
  scoring = list(joint_draws = 1000L)
)


# 2. Complete simulation configuration  ----------------------------------------

# Merge custom with default simulation configuration 
merge_config <- function(config, overrides, path = "config") {
  if (!is.list(overrides) || is.data.frame(overrides)) {
    stop(path, " overrides must be a named list.")
  }
  if (!length(overrides)) {
    return(config)
  }
  if (is.null(names(overrides)) || anyNA(names(overrides)) ||
      any(!nzchar(names(overrides))) || anyDuplicated(names(overrides))) {
    stop(path, " overrides must have unique, nonempty names.")
  }
  unknown <- setdiff(names(overrides), names(config))
  if (length(unknown)) {
    stop("Unknown settings in ", path, ": ", paste(unknown, collapse = ", "))
  }
  for (name in names(overrides)) {
    value <- overrides[[name]]
    if (is.list(config[[name]]) && !is.data.frame(config[[name]])) {
      value <- merge_config(config[[name]], value, paste(path, name, sep = "$"))
    }
    config[name] <- list(value)
  }
  config
}

# Build a single complete simulation configuration
study_config <- function(settings, section, project_dir) {
  experiments <- settings$experiments
  settings$experiments <- NULL
  config <- merge_config(default_simulation_config, settings)

  # Apply experiment overrides only for studies that define them
  if (!is.null(experiments)) {
    if (!is.list(experiments) || length(section) != 1L ||
        !section %in% names(experiments)) {
      stop("Unknown experiment: ", section)
    }
    config <- merge_config(config, experiments[[section]])
  }

  # Complete internal noise labels and project paths after merging settings
  noise <- config$generating$noise
  if (is.data.frame(noise) && !"noise_id" %in% names(noise)) {
    noise$noise_id <- if (nrow(noise) == 1L) {
      "equal"
    } else {
      paste0("noise_", seq_len(nrow(noise)))
    }
    config$generating$noise <- noise
  }
  config$project_dir <- normalizePath(project_dir, mustWork = TRUE)

  config
}


# 3. Validate a study configuration ------------------------------------------
validate_config <- function(config) {
  
  # Helper functions ---------------------------------------------------------
  positive <- function(x) {
    # One or more finite, strictly positive numeric values
    is.numeric(x) && length(x) > 0L &&
      all(is.finite(x)) && all(x > 0)
  }
  
  integer_scalar <- function(x) {
    # Exactly one positive integer-valued number 
    positive(x) && length(x) == 1L && x == floor(x)
  }
  
  integer_vector <- function(x) {
    # One or more distinct positive integer-valued numbers
    # anyDuplicated(x) returns 0 when x contains no repeated values
    positive(x) && all(x == floor(x)) && !anyDuplicated(x)
  }
  
  generating <- config$generating
  
  # Check replications, random seed, and data sizes 
  stopifnot(
    integer_scalar(config$n_replicates),
    integer_scalar(config$base_seed),
    config$base_seed < 2147483647,
    integer_vector(generating$I),
    integer_vector(generating$J),
    integer_scalar(generating$grid_n),
    generating$grid_n >= 2L
  )
  
  # Check GRF 
  # sigma2 = marginal variance of phi(s), which must be positive
  # rho    = practical spatial range; several distinct values are allowed
  # nu     = Matern smoothness, fixed at 1 in this simulation study
  stopifnot(
    length(generating$mu) == 1L,
    is.finite(generating$mu),
    length(generating$alpha) == 1L,
    is.finite(generating$alpha),
    positive(generating$sigma2),
    length(generating$sigma2) == 1L,
    positive(generating$rho),
    !anyDuplicated(generating$rho),
    identical(as.numeric(generating$nu), 1)
  )
  
  # Check spatial domain
  # rectangle defined by x_limits and y_limits
  # each interval must have two finite endpoints, with lower < upper
  # Areal units form a sqrt(J) by sqrt(J) partition
  for (limits in list(generating$x_limits, generating$y_limits)) {
    stopifnot(
      length(limits) == 2L,
      all(is.finite(limits)),
      diff(limits) > 0
    )
  }
  areas_per_axis <- sqrt(generating$J)
  
  if (any(areas_per_axis != floor(areas_per_axis)) ||
      any(generating$grid_n %% areas_per_axis != 0)) {
    stop("J must be a perfect square whose side divides grid_n.")
  }
  
  # Check preferential sampling mechanism
  # sampling_scenarios() creates one row for each sampling scenario:
  # constant:
  #   log(lambda(s)) = alpha + gamma * phi(s)
  #
  # piecewise [EXTENSION]:
  #   log(lambda(s)) = alpha + gamma(s) * phi(s),
  #   where gamma(s) equals gamma_L or gamma_R on opposite sides of split_x
  #
  # independent_field [EXTENSION]:
  #   log(lambda(s)) = alpha + gamma * phi(s) + delta * U(s),
  #   where U(s) is an additional GRF independent of phi(s)
  sampling <- sampling_scenarios(config)
  stopifnot(
    nrow(sampling) > 0L,
    !anyNA(sampling$sampling_id),
    !anyDuplicated(sampling$sampling_id),
    all(grepl("^[A-Za-z0-9_-]+$", sampling$sampling_id)),
    all(sampling$sampling_type %in% c(
      "constant", "piecewise", "independent_field"
    ))
  )
  
  for (row in seq_len(nrow(sampling))) {
    mechanism <- sampling[row, ]
    
    if (mechanism$sampling_type == "piecewise") {
      # Both regional preferentiality coefficients must be finite
      # The dividing coordinate must lie strictly within the x-domain,
      # so that it separates the domain into two nonempty subregions
      stopifnot(
        all(is.finite(c(
          mechanism$gamma_L,
          mechanism$gamma_R,
          mechanism$split_x
        ))),
        mechanism$split_x > generating$x_limits[1],
        mechanism$split_x < generating$x_limits[2]
      )
    } else {
      # Constant and independent-field mechanisms use a single gamma
      stopifnot(is.finite(mechanism$gamma))
    }
    
    if (mechanism$sampling_type == "independent_field") {
      # delta >= 0 controls the contribution of the additional field U(s)
      # rho_U and sigma2_U define its spatial range and marginal variance
      stopifnot(
        is.finite(mechanism$delta),
        mechanism$delta >= 0,
        positive(mechanism$rho_U),
        positive(mechanism$sigma2_U)
      )
    } else {
      # Other mechanisms must not accidentally include the extra field
      stopifnot(is.finite(mechanism$delta), mechanism$delta == 0)
    }
  }
  
  # Check point and areal measurement-error scenarios --------------------------
  # Each row specifies ONE pair of observation-error variances
  #   tau2_s = variance for point measurements
  #   tau2_B = variance for areal measurements
  noise <- generating$noise
  
  stopifnot(
    is.data.frame(noise),
    nrow(noise) > 0L,
    all(c("noise_id", "tau2_s", "tau2_B") %in% names(noise)),
    is.character(noise$noise_id),
    !anyNA(noise$noise_id),
    !anyDuplicated(noise$noise_id),
    all(grepl("^[A-Za-z0-9_-]+$", noise$noise_id)),
    positive(noise$tau2_s),
    positive(noise$tau2_B)
  )
  
  # Check models and computational settings ------------------------------------
  stopifnot(
    is.character(config$models),
    length(config$models) > 0L,
    !anyDuplicated(config$models),
    all(config$models %in% c("Melding", "PSgeo", "PSmelding"))
  )
  
  # INLA settings and posterior draws
  stopifnot(
    integer_scalar(config$fitting$inla_threads),
    integer_scalar(config$scoring$joint_draws),
    config$scoring$joint_draws >= 2L
  )
  
  # offset and max_edge each require an inner and outer value
  # cutoff is a single positive minimum-separation setting for mesh vertices
  stopifnot(
    length(config$mesh$offset) == 2L,
    positive(config$mesh$offset),
    length(config$mesh$max_edge) == 2L,
    positive(config$mesh$max_edge),
    length(config$mesh$cutoff) == 1L,
    positive(config$mesh$cutoff)
  )
  
  # Return the original configuration without printing it to the console
  invisible(config)
}


# 4. Prepare the R/INLA environment ------------------------------------------
configure_runtime <- function(config) {
  if (dir.exists(config$package_library)) {
    .libPaths(c(config$package_library, .libPaths()))
  }

  # Enforce the R version for which the original local environment was built
  if (getRversion() < "4.2" || getRversion() >= "4.3") {
    stop("Use the existing R 4.2 launcher; this environment was built for R 4.2.")
  }

  if (Sys.info()[["sysname"]] == "Darwin") {
    # Help the compiled INLA components find R and existing GCC libraries
    # Preserve any other DYLD_LIBRARY_PATH entries already configured
    Sys.setenv(INLA_RHOME = R.home())

    gcc_dir <- "/opt/homebrew/opt/gcc/lib/gcc/current"
    existing_paths <- strsplit(
      Sys.getenv("DYLD_LIBRARY_PATH"), ":", fixed = TRUE
    )[[1]]

    Sys.setenv(
      DYLD_LIBRARY_PATH = paste(
        unique(c(gcc_dir, existing_paths[nzchar(existing_paths)])),
        collapse = ":"
      )
    )
  }

  required_packages <- c(
    "INLA", "RandomFields", "raster", "sp", "sf", "Matrix",
    "spatstat.geom", "spatstat.random"
  )

  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages)) {
    stop("Missing packages: ", paste(missing_packages, collapse = ", "))
  }

  inla_version <- as.character(utils::packageVersion("INLA"))

  if (utils::compareVersion(inla_version, "23.04.24") != 0L) {
    stop(
      "This baseline preserves INLA 23.04.24 defaults; ",
      "review priors before changing versions."
    )
  }

  invisible(required_packages)
}



# 5. Define sampling scenarios and assemble the simulation plan ----------------

# Build table of sampling mechanisms and parameters  ---------------------------
sampling_scenarios <- function(config) {
  
  required_columns <- c(
    "sampling_id", 
    "sampling_type", 
    "gamma", 
    "gamma_L", 
    "gamma_R", 
    "split_x",
    "delta", 
    "rho_U", 
    "sigma2_U"
  )
  
  sampling <- config$generating$sampling
  
  # If there is no explicit extension table, create one constant-gamma
  # scenario for each value in config$generating$gamma
  if (is.null(sampling)) {
    gamma_values <- config$generating$gamma
    
    sampling <- data.frame(
      sampling_id = paste0("gamma_", gsub("-", "minus", gamma_values)),
      sampling_type = "constant",
      gamma = gamma_values,
      gamma_L = NA_real_,
      gamma_R = NA_real_,
      split_x = NA_real_,
      delta = 0,
      rho_U = NA_real_,
      sigma2_U = NA_real_
    )
    
    # Safe file naming
    sampling$sampling_id <- gsub("[.]", "p", sampling$sampling_id)
  }
  
  # Check that the table has the expected columns and numeric parameter types
  stopifnot(
    is.data.frame(sampling),
    all(required_columns %in% names(sampling)),
    is.character(sampling$sampling_id),
    is.character(sampling$sampling_type),
    all(vapply(
      sampling[required_columns[-c(1, 2)]],
      is.numeric,
      logical(1)
    ))
  )
  
  # Fix the column order for compatibility with the simulation runner
  sampling[required_columns]
}


# Generate seeds while preserving paired comparisons ---------------------------
sampling_seed <- function(config, scenario, component, model = NULL) {
  # Equivalent constant mechanisms share seeds, including the extension baselines
  mechanism <- scenario$sampling_type
  gamma <- scenario$gamma
  if (mechanism == "piecewise" && scenario$gamma_L == scenario$gamma_R) {
    mechanism <- "constant"
    gamma <- scenario$gamma_L
  }
  if (mechanism == "independent_field" && scenario$delta == 0) {
    mechanism <- "constant"
  }

  key <- list(component, scenario$rho, scenario$replication, scenario$I)
  if (component != "locations") {
    # Labels should not distinguish the same equal-noise baseline across studies
    noise_id <- if (scenario$tau2_s == .1 && scenario$tau2_B == .1) {
      "equal"
    } else {
      scenario$noise_id
    }
    key <- c(key, list(if (model == "PSgeo") "no_areas" else scenario$J,
      gamma, noise_id, model))
  } else {
    key <- c(key, list(gamma))
  }

  if (mechanism != "constant") {
    sampling_key <- paste(mechanism, scenario$gamma_L, scenario$gamma_R,
      scenario$split_x, scenario$delta, scenario$rho_U, scenario$sigma2_U,
      sep = "|")
    key <- c(key, list(sampling_key))
  }
  do.call(simulation_seed, c(list(config), key))
}


# Expand the configured experimental factors into simulation cases

scenario_plan <- function(config) {
  validate_config(config)
  
  generating <- config$generating
  sampling <- sampling_scenarios(config)
  
  # Cross the independent experimental factors ---------------------------------
  # expand.grid() returns ALL combinations of:
  #   - replication number;
  #   - Matérn spatial range rho;
  #   - number of point observations I;
  #   - number of areal observations J;
  #   - one ROW of the sampling-mechanism table;
  #   - one noise scenario, identified by noise_id.
  # Preserve the factor ordering: it determines the order of case_id values
  scenarios <- expand.grid(
    replication = seq_len(config$n_replicates),
    rho = generating$rho,
    I = generating$I,
    J = generating$J,
    sampling_row = seq_len(nrow(sampling)),
    noise_id = generating$noise$noise_id,
    stringsAsFactors = FALSE
  )
  
  # Expand sampling-row indices into actual mechanism parameters ---------------
  # The temporary sampling_row column selects the correct row of the table
  scenarios <- cbind(
    scenarios,
    sampling[scenarios$sampling_row, , drop = FALSE]
  )
  scenarios$sampling_row <- NULL
  rownames(scenarios) <- NULL
  
  # Label each scenario with its experiment or use "development" by default
  scenarios$experiment <- if (is.null(config$section)) {
    "development"
  } else {
    config$section
  }
  
  # Record theoretical preferential-dependence strength ------------------------
  # For the independent-field mechanism:
  #   log(lambda(s)) = alpha + gamma * phi(s) + delta * U(s),
  # where phi and U are independent, zero-mean Gaussian fields with
  # variances sigma2 and sigma2_U, respectively
  # At any fixed location s:
  #   Corr(log(lambda(s)), phi(s))
  #     = gamma * sqrt(sigma2) /
  #       sqrt(gamma^2 * sigma2 + delta^2 * sigma2_U).
  scenarios$log_intensity_phi_cor <- NA_real_
  independent_rows <- scenarios$sampling_type == "independent_field"
  
  intensity_sd <- sqrt(
    scenarios$gamma[independent_rows]^2 * generating$sigma2 +
      scenarios$delta[independent_rows]^2 *
      scenarios$sigma2_U[independent_rows]
  )
  
  scenarios$log_intensity_phi_cor[independent_rows] <- ifelse(
    intensity_sd > 0,
    scenarios$gamma[independent_rows] * sqrt(generating$sigma2) / intensity_sd,
    NA_real_
  )
  
  # Attach the paired point and areal observation-noise variances --------------
  # match() selects the ONE row associated with each noise_id
  # This preserves the intended pair (tau2_s, tau2_B) for every simulation scenario
  noise_rows <- match(scenarios$noise_id, generating$noise$noise_id)
  scenarios$tau2_s <- generating$noise$tau2_s[noise_rows]
  scenarios$tau2_B <- generating$noise$tau2_B[noise_rows]
  scenarios$case_id <- sprintf("case_%04d", seq_len(nrow(scenarios)))
  
  scenarios
}


# 6. Reproducible random seeds ------------------------------------------------

# The same base seed and labels give the same seed, regardless of run order
# Different labels can occasionally produce the same seed
simulation_seed <- function(config, ...) {
  hash <- as.double(config$base_seed)
  key <- paste(..., collapse = "|")

  # Convert labels to a seed within R's integer range
  for (character_code in utf8ToInt(key)) {
    hash <- (hash * 131 + character_code) %% 2147483646
  }

  as.integer(hash + 1)
}


# Run code with a chosen seed, then restore the previous random-number state
with_simulation_seed <- function(seed, code) {
  # Save the random-number state, if one exists
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  previous_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }

  # Restore the original state or remove the temporary one, even after an error
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })

  set.seed(seed)

  # Evaluate code after setting the seed
  force(code)
}
