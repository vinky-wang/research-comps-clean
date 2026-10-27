# Generate latent field
generate_field <- function(rho) {
  # Simulate Matern Gaussian random field
  field = spatstat.random::rGRFmatern(
    W = spatstat.geom::owin(xrange = c(0, 1), yrange = c(0, 1)),
    mu = 0,
    var = 1,
    scale = rho / 2,
    nu = 1,
    dimyx = c(50, 50)
  )
  
  # Convert to data frame
  latent = as.data.frame(field)
  names(latent) = c("x", "y", "phi")
  latent$Z = latent$phi
  
  return(latent)
}



# Generate point observations
generate_points <- function(field, I, tau2_s, log_lambda) {
  # Fixed grid-cell width, D = [0,1]^2 and 50 x 50 grid
  h = 1 / 50
  
  # Compute sampling probabilities conditional on exactly I points
  if (length(log_lambda) != nrow(field) ||
      any(!is.finite(log_lambda))) {
    stop("Invalid log-intensity values.")
  }
  
  prob = exp(log_lambda - max(log_lambda))
  prob = prob / sum(prob)
  
  # Sample point locations
  selected = sample.int(n = nrow(field),
                        size = I,
                        replace = TRUE,
                        prob = prob)
  
  # Add uniform jitter within selected grid cells
  # Use each selected cell's field value as the approximate point truth
  points = data.frame(
    x = field$x[selected] + runif(I, -h / 2, h / 2),
    y = field$y[selected] + runif(I, -h / 2, h / 2),
    true = field$Z[selected]
  )
  
  # Generate point observations
  points$value = points$true + rnorm(I, sd = sqrt(tau2_s))
  
  return(list(points = points, sampling_probability = prob))
}



# Generate areal observations
generate_areas <- function(field, J, tau2_B) {
  # Number of areas along each axis
  if (length(J) != 1L || !is.finite(J) ||
      J < 1 || J != floor(J)) {
    stop("J must be a positive integer.")
  }
  
  side = as.integer(sqrt(J))
  
  if (side^2 != J || 50L %% side != 0L) {
    stop("J must be a perfect square and sqrt(J) must divide 50.")
  }
  
  # Assign grid cells to areas
  ix = pmin(side, floor(field$x * side) + 1L)
  iy = pmin(side, floor(field$y * side) + 1L)
  area_id = (iy - 1L) * side + ix
  
  # Construct area boundaries
  area_grid = expand.grid(col = seq_len(side), row = seq_len(side))
  
  areas = data.frame(
    area_id = seq_len(J),
    xmin = (area_grid$col - 1) / side,
    xmax = area_grid$col / side,
    ymin = (area_grid$row - 1) / side,
    ymax = area_grid$row / side
  )
  
  # Compute true areal averages
  areas$true = vapply(seq_len(J),
                      function(j) mean(field$Z[area_id == j]),
                      numeric(1))
  
  # Generate areal observations
  areas$value = areas$true + rnorm(J, sd = sqrt(tau2_B))
  
  return(areas)
}



# Generate simulated data
generate_data <- function(
    rho = 0.2,
    gamma = 1,
    I = 100L,
    J = 25L,
    tau2_s = 0.1,
    tau2_B = 0.1,
    seed = 123L,
    scenario = c("constant", "spatially_varying", "independent_field"),
    delta = 0.5,
    beta_psi = 0.5,
    rho_psi = 0.1
) {
  scenario = match.arg(scenario)
  
  set.seed(seed)
  seeds = sample.int(.Machine$integer.max, 4L)
  
  # Generate latent field
  set.seed(seeds[1])
  field = generate_field(rho)
  
  # Specify preferential sampling mechanism
  if (scenario == "constant") {
    log_lambda = gamma * field$phi
    
  } else if (scenario == "spatially_varying") {
    # Divide domain into five horizontal regions
    region = pmin(5L, floor(5 * field$y) + 1L)
    
    # Region-specific preferentiality coefficients
    gamma_s = gamma + delta * (region - 3) / 2
    
    log_lambda = gamma_s * field$phi
    
    # Store regional information
    field$region = region
    field$gamma_s = gamma_s
    
  } else if (scenario == "independent_field") {
    # Generate an additional independent Gaussian field
    set.seed(seeds[2])
    field2 = generate_field(rho_psi)
    
    stopifnot(identical(field[, c("x", "y")],
                        field2[, c("x", "y")]))
    
    field$psi = field2$phi
    
    log_lambda = gamma * field$phi + beta_psi * field$psi
  }
  
  # Generate point observations
  set.seed(seeds[3])
  point_data = generate_points(field = field,
                               I = I,
                               tau2_s = tau2_s,
                               log_lambda = log_lambda)
  
  # Generate areal observations
  set.seed(seeds[4])
  areas = generate_areas(field = field, J = J, tau2_B = tau2_B)
  
  # Store sampling information
  field$log_lambda = log_lambda
  field$sampling_probability = point_data$sampling_probability
  
  # Return simulated data
  return(list(
    latent = field,
    points = point_data$points,
    areas = areas,
    prediction = data.frame(
      x = field$x,
      y = field$y,
      true = field$Z
    ),
    scenario = scenario,
    seed = seed
  ))
}