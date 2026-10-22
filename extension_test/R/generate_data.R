# Data generation ------------------------------------------------------------
# [ORIGINAL] Adapted from PSmelding-master/functions/samplegenerator.r:
# latt_generation(), pcoxsample(), datagenerator(), and areasample()
# No original scripts are sourced. Sourcing this file does not generate data

# Read in order: sampling formula -> latent fields -> one observed dataset
# Inputs come from config.R; outputs feed fit_model() and evaluate_fit()
# Raster values, coordinates and truth vectors must keep the same cell order

# 1. Define the sampling intensity --------------------------------------------
# sampling_config describes ONE scenario: gamma, gamma_L/R, alpha and delta are scalars
# Expand experiment vectors through scenario_plan() before evaluating them
# These functions are deterministic: they neither draw fields nor sample sites
gamma_surface <- function(s1, sampling_config) {
  if (!is.numeric(s1) || any(!is.finite(s1))) {
    stop("s1 must contain finite coordinates.")
  }
  mechanism <- sampling_config$sampling_type
  if (!is.character(mechanism) || length(mechanism) != 1L || is.na(mechanism)) {
    stop("Supply one sampling_type.")
  }
  scalar <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop(name, " must be one finite number; resolve the scenario first.")
    }
    value
  }
  spatial <- function() {
    split_x <- if (is.null(sampling_config$split_x)) {
      .5
    } else {
      scalar(sampling_config$split_x, "split_x")
    }
    ifelse(
      s1 < split_x,
      scalar(sampling_config$gamma_L, "gamma_L"),
      scalar(sampling_config$gamma_R, "gamma_R")
    )
  }
  switch(
    mechanism,
    constant = rep(scalar(sampling_config$gamma, "gamma"), length(s1)),
    piecewise = spatial(),  # Backward-compatible name for the same mechanism
    independent_field = rep(scalar(sampling_config$gamma, "gamma"), length(s1)),
    stop("Unknown sampling_type: ", mechanism)
  )
}

log_intensity <- function(coords, phi, U = NULL, sampling_config) {
  coords <- as.matrix(coords)
  if (!is.numeric(coords) || ncol(coords) != 2L || any(!is.finite(coords))) {
    stop("coords must have two finite numeric columns, in the same order as phi and U.")
  }
  if (!is.numeric(phi) || length(phi) != nrow(coords) || any(!is.finite(phi))) {
    stop("phi must contain one finite field value per coordinate.")
  }
  if (!is.numeric(sampling_config$alpha) ||
      length(sampling_config$alpha) != 1L || !is.finite(sampling_config$alpha) ||
      !is.numeric(sampling_config$delta) ||
      length(sampling_config$delta) != 1L || !is.finite(sampling_config$delta)) {
    stop("alpha and delta must each be one finite number.")
  }
  if (is.null(U)) {
    if (sampling_config$delta != 0) {
      stop("Supply U when delta is nonzero.")
    }
    U <- rep(0, nrow(coords))
  }
  if (!is.numeric(U) || length(U) != nrow(coords) || any(!is.finite(U))) {
    stop("U must contain one finite field value per coordinate.")
  }
  gamma_s <- gamma_surface(coords[, 1], sampling_config)
  # phi is the zero-mean response field (Z-mu). U is generated independently
  # upstream; this function applies the intensity formula without transforming fields
  value <- sampling_config$alpha + (gamma_s * phi + sampling_config$delta * U)
  if (any(!is.finite(value))) {
    stop("Nonfinite log intensity.")
  }
  value
}

# 2. Generate one response surface Z(s) = mu + phi(s) ---------------------------
# Evaluate the Matérn field at cell centres. The seed depends on the replication
# and range, allowing the same field to be reused across sampling/noise settings
generate_surface <- function(config, replication, component = "surface") {
  generating <- config$generating
  cell_centres <- function(limits) {
    limits[1] + (seq_len(generating$grid_n) - .5) *
      diff(limits) / generating$grid_n
  }
  coordinates <- expand.grid(x = cell_centres(generating$x_limits), y = cell_centres(generating$y_limits))
  seed <- simulation_seed(config, component, generating$rho, replication)
  process <- with_simulation_seed(
    seed,
    {
      # [ORIGINAL] For this nu=1 parameterization, RandomFields scale = rho/2
      covariance <- RandomFields::RMmatern(
        nu = generating$nu,
        var = generating$sigma2,
        scale = generating$rho / 2
      )
      RandomFields::RFsimulate(covariance + RandomFields::RMtrend(mean = generating$mu), x = coordinates)
    }
  )
  xyz <- data.frame(x = process@coords[, 1], y = process@coords[, 2], z = process$variable1)
  surface <- raster::rasterFromXYZ(xyz)
  # [EXTENSION] The simulation uses a planar unit system, not longitude/latitude
  raster::crs(surface) <- NA_character_
  if (any(!is.finite(raster::values(surface)))) {
    stop("Nonfinite generated field.")
  }
  # No realized-field centering or rescaling: mu and sigma2 are ensemble values
  list(raster = surface, seed = seed)
}

# [EXTENSION] Independent clustering fields use a separate random seed and the same
# Matérn convention as phi. No centering, rescaling, or sample orthogonalization:
# independence is a generating assumption, not zero empirical correlation
# 3. Generate the optional independent clustering field U(s) ------------------
# Cache one field per (rho_U, sigma2_U) pair. Varying delta reweights that same
# realization, so comparisons isolate the effect of unexplained clustering
generate_sampling_fields <- function(config, replication) {
  sampling <- sampling_scenarios(config)
  settings <- unique(sampling[sampling$sampling_type == "independent_field",
    c("rho_U", "sigma2_U"), drop = FALSE])
  fields <- list()
  for (row in seq_len(nrow(settings))) {
    current <- config
    current$generating$mu <- 0
    current$generating$rho <- settings$rho_U[row]
    current$generating$sigma2 <- settings$sigma2_U[row]
    key <- paste(settings$rho_U[row], settings$sigma2_U[row], sep = "|")
    fields[[key]] <- generate_surface(current, replication, component = "clustering_U")
  }
  fields
}

# The returned arrays use the response raster's cell order. Only sampling
# locations depend on gamma(s) or U(s); responses remain mu + phi + noise
# 4. Evaluate the intensity and record its generating truth -------------------
# A piecewise mechanism has a scalar gamma truth only when both sides agree
# Otherwise the fitted constant gamma has no single generating value to score
sampling_truth <- function(config, scenario, field) {
  coordinates <- raster::xyFromCell(field$raster, seq_len(raster::ncell(field$raster)))
  phi <- as.numeric(raster::values(field$raster)) - config$generating$mu
  sampling_config <- as.list(scenario)
  sampling_config$alpha <- config$generating$alpha
  gamma_s <- gamma_surface(coordinates[, 1], sampling_config)
  gamma_truth <- scenario$gamma
  U <- rep(0, length(phi))
  U_seed <- NA_integer_
  if (scenario$sampling_type == "piecewise") {
    gamma_truth <- if (scenario$gamma_L == scenario$gamma_R) {
      scenario$gamma_L
    } else {
      NA_real_
    }
  } else if (scenario$sampling_type == "independent_field") {
    key <- paste(scenario$rho_U, scenario$sigma2_U, sep = "|")
    independent_fields <- field$independent_fields
    if (is.null(independent_fields)) {
      independent_fields <- generate_sampling_fields(config, scenario$replication)
    }
    independent <- independent_fields[[key]]
    stopifnot(
      !is.null(independent),
      identical(
        raster::xyFromCell(independent$raster, seq_len(raster::ncell(independent$raster))),
        coordinates
      )
    )
    U <- as.numeric(raster::values(independent$raster))
    U_seed <- independent$seed
  }
  full_log_intensity <- log_intensity(coordinates, phi, U, sampling_config)
  # [ORIGINAL] Conditional on I, alpha cancels. Evaluate with alpha=0 for weights
  # instead of subtracting alpha afterward, preserving the existing numerics
  sampling_config$alpha <- 0
  log_weight <- log_intensity(coordinates, phi, U, sampling_config)
  list(
    gamma_s = gamma_s,
    gamma = gamma_truth,
    gamma_L = scenario$gamma_L,
    gamma_R = scenario$gamma_R,
    split_x = scenario$split_x,
    phi = phi,
    U = U,
    U_seed = U_seed,
    log_weight = log_weight,
    log_intensity = full_log_intensity,
    # This is the generating intensity before conditioning on the fixed point count
    # Keep log intensity as well, since exp() can overflow for extreme parameters
    intensity = exp(full_log_intensity),
    mechanism = scenario$sampling_type,
    gamma_note = if (is.na(gamma_truth)) {
      "No single true gamma: fitted constant coefficient approximates spatially varying preferentiality."
    } else
    if (scenario$sampling_type == "independent_field" && scenario$delta > 0) {
      "True response-field coefficient; fitted intensity omits the independent U field."
    } else {
      ""
    }
  )
}

# 5. Generate locations, observations, and prediction targets ------------------
# One row of scenario_plan() supplies scalar settings for this dataset
# The returned list keeps observed data separate from latent scoring truth
generate_dataset <- function(config, scenario, field = NULL) {
  generating <- config$generating
  generating$rho <- scenario$rho
  generating$I <- scenario$I
  generating$J <- scenario$J
  config$generating <- generating
  if (is.null(field)) {
    field <- generate_surface(config, scenario$replication)
  }
  surface <- field$raster
  # Separate seeds pair comparisons deliberately: point noise does not depend
  # on J; areal noise does not depend on I. Noise variances scale the same
  # standard-normal draws when their settings change
  seeds <- list(
    surface = field$seed,
    locations = sampling_seed(config, scenario, "locations"),
    point_noise = simulation_seed(config, "point-noise", scenario$rho, scenario$replication, scenario$I),
    area_noise = simulation_seed(config, "area-noise", scenario$rho, scenario$replication, scenario$J)
  )

  # [ORIGINAL] Condition on I; the constant alpha cancels on normalization
  # [EXTENSION] log(lambda)=alpha+gamma(s)*phi(s)+delta*U(s); phi=Z-mu
  # [EXTENSION] Subtracting the maximum log weight avoids overflow; fields are unchanged
  sampling <- sampling_truth(config, scenario, field)
  seeds$clustering_U <- sampling$U_seed
  log_weight <- sampling$log_weight
  weight <- exp(log_weight - max(log_weight))
  # Raster cells run from the top row downward; spatstat images use increasing
  # y coordinates. Reverse the rows and y coordinates together to avoid a
  # vertically reflected sampling density
  weight_matrix <- matrix(weight, nrow = nrow(surface), ncol = ncol(surface), byrow = TRUE)
  density <- spatstat.geom::im(
    weight_matrix[nrow(surface):1L, , drop = FALSE],
    xcol = raster::xFromCol(surface, seq_len(ncol(surface))),
    yrow = rev(raster::yFromRow(surface, seq_len(nrow(surface))))
  )
  window <- spatstat.geom::owin(xrange = generating$x_limits, yrange = generating$y_limits)
  sites <- with_simulation_seed(
    seeds$locations,
    spatstat.random::rpoint(n = generating$I, f = density, win = window)
  )
  locations <- cbind(x = sites$x, y = sites$y)
  # Point observations use the raster value at each sampled location, plus
  # independent Gaussian measurement error with variance tau2_s
  point_truth <- raster::extract(surface, locations)
  point_values <- point_truth + with_simulation_seed(seeds$point_noise,
    stats::rnorm(generating$I, sd = sqrt(scenario$tau2_s)))
  points <- sf::st_as_sf(
    data.frame(x = locations[, 1], y = locations[, 2],
      truth = point_truth, value = point_values),
    coords = c("x", "y")
  )

  # Equal-sized raster cells align exactly with these square areal blocks
  # Their arithmetic average is the area-weighted average on this grid
  blocks <- raster::aggregate(surface, fact = generating$grid_n / sqrt(generating$J), fun = mean)
  names(blocks) <- "truth"
  areas <- sf::st_as_sf(raster::rasterToPolygons(blocks))
  areas <- sf::st_set_crs(areas, sf::NA_crs_)
  areas$value <- areas$truth + with_simulation_seed(seeds$area_noise,
    stats::rnorm(generating$J, sd = sqrt(scenario$tau2_B)))

  coordinates <- raster::xyFromCell(surface, seq_len(raster::ncell(surface)))
  colnames(coordinates) <- c("x", "y")
  prediction_points <- sf::st_as_sf(as.data.frame(coordinates), coords = c("x", "y"))
  Z <- as.numeric(raster::values(surface))
  # [EXTENSION] Score predictions against latent truth Z without drawing
  # additional held-out responses; observed point/area noise remains unchanged
  outline <- rbind(
    c(generating$x_limits[1], generating$y_limits[1]),
    c(generating$x_limits[2], generating$y_limits[1]),
    c(generating$x_limits[2], generating$y_limits[2]),
    c(generating$x_limits[1], generating$y_limits[2]),
    c(generating$x_limits[1], generating$y_limits[1])
  )
  boundary <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(outline))))
  # Check generated values after the random draws and spatial conversions;
  # configuration validation cannot detect failures in those operations
  stopifnot(all(is.finite(points$value)), all(is.finite(areas$value)),
    all(is.finite(Z)))
  list(
    points = points,
    areas = areas,
    prediction_points = prediction_points,
    boundary = boundary,
    field = field,
    scenario = scenario,
    seeds = seeds,
    truth = list(
      Z = Z,
      coordinates = coordinates,
      sampling = sampling,
      parameters = list(
        mu = generating$mu,
        sigma2 = generating$sigma2,
        nu = generating$nu,
        rho = generating$rho,
        theta = 8 * generating$sigma2 / generating$rho^2,
        alpha = NA_real_,
        gamma = sampling$gamma,
        tau2_s = scenario$tau2_s,
        tau2_B = scenario$tau2_B
      )
    )
  )
}
