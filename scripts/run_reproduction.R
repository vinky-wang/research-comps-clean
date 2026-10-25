source("functions/generate_data.R")
source("functions/create_summary.R")
source("functions/fit_models.R")


# Generate simulated data ------------------------------------------------------

# Preferential sampling scenarios
scenarios <- c(
  "Constant gamma" = "constant",
  "Spatially varying gamma" = "spatially_varying",
  "Independent latent field" = "independent_field"
)

# Generate data for each scenario
simulations <- lapply(scenarios, function(scenario) {
  generate_data(
    rho = 0.2,
    gamma = 1,
    I = 100,
    J = 25,
    tau2_s = 0.1,
    tau2_B = 0.1,
    seed = 123,
    scenario = scenario,
    delta = 1,
    beta_psi = 1,
    rho_psi = 0.2
  )
})

# Extract simulation outputs
field <- simulations[[1]]$latent
log_lambda <- lapply(simulations, function(d) d$latent$log_lambda)
points <- lapply(simulations, function(d) d$points)

# Select one simulated dataset
dataset <- simulations[["Constant gamma"]]


# Plot simulated data ----------------------------------------------------------

# Preferential sampling scenarios
fig1 <- plot_sampling_scenarios(
  field = field,
  scenarios = log_lambda,
  points = points,
  limits = c(-5, 5)
)

fig1

# Simulated surface, areal observations, and point observations
fig2 <- plot_simulated_data(dataset)

fig2


# Fit models across 10 replications --------------------------------------------

# Keep the same mesh and model specifications across replications
model_names <- c("ps_melding", "melding", "ps_geo")
mesh <- create_mesh()

n_replications <- 2L
simulation_seeds <- 123L + seq_len(n_replications) - 1L

# Set output directories
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path("results", "reproduction")
fits_dir <- file.path(output_dir, "fits")

dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)

# Generate independent datasets and fit all four models to each
replication_results <- lapply(seq_len(n_replications), function(replication) {
  message("Replication ", replication, " of ", n_replications)
  
  # Generate a dataset under the constant-gamma scenario
  simulated <- generate_data(
    rho = 0.2,
    gamma = 1,
    I = 100,
    J = 25,
    tau2_s = 0.1,
    tau2_B = 0.1,
    seed = simulation_seeds[replication],
    scenario = "constant"
  )
  
  # Fit and log all four models on the same dataset
  logged <- fit_models_logged(
    data = simulated,
    models = model_names,
    mesh = mesh,
    replication = replication,
    seed = simulation_seeds[replication],
    point_process_method = "mesh_dual",
    areal_method = "vertices",
    separate_noise = FALSE,
    prior_gamma = c(1, 10)
  )
  
  # Save simulated data, fitted models, and fitting diagnostics
  saveRDS(
    list(
      data = simulated,
      results = logged$results,
      log = logged$log
    ),
    file = file.path(
      fits_dir,
      sprintf("rep_%03d.rds", replication)
    )
  )
  
  # Evaluate only numerically valid fits
  metrics <- NULL
  
  if (length(logged$results$models) > 0L) {
    metrics <- evaluate_replication(
      results = logged$results,
      truth = simulated$prediction$true,
      tau2_s = 0.1,
      gamma_true = 1
    )
    
    metrics$replication <- replication
  }
  
  return(list(metrics = metrics, log = logged$log))
})

# Combine replication-level metrics
replication_metrics <- do.call(rbind, lapply(replication_results, function(x) {
  x$metrics
}))

# Combine fitting logs, including errors and invalid fits
fit_log <- do.call(rbind, lapply(replication_results, function(x) {
  x$log
}))

rownames(fit_log) <- NULL

if (is.null(replication_metrics)) {
  stop("No valid model fits were available for evaluation. Inspect fit_log.")
}

rownames(replication_metrics) <- NULL

# Inspect replication-level metrics and fitting diagnostics
head(replication_metrics)
head(fit_log)

# Summarize fitting outcomes
fit_status <- summarize_fit_log(fit_log, models = model_names)

fit_status

# Inspect fits that require investigation
subset(fit_log, status != "ok")


# Summarize simulation results -------------------------------------------------

summary_results <- summarize_replications(
  metrics = replication_metrics,
  gamma_true = 1
)

# Main predictive performance
summary_results$prediction_n
summary_results$prediction_mean
summary_results$prediction_sd
summary_results$prediction_mcse

# Zhong's original Wasserstein calculation
summary_results$prediction_reproduction
summary_results$prediction_reproduction_mcse

# Posterior uncertainty diagnostics
summary_results$prediction_diagnostics
summary_results$prediction_diagnostics_mcse

# Preferentiality parameter recovery
summary_results$gamma
summary_results$gamma_mcse


# Record simulation settings ---------------------------------------------------

run_info <- list(
  r_version = R.version.string,
  inla_version = INLA::inla.version("version"),
  session_info = utils::sessionInfo(),
  models = model_names,
  n_replications = n_replications,
  simulation_seeds = simulation_seeds,
  generating = list(
    rho = 0.2,
    gamma = 1,
    I = 100,
    J = 25,
    tau2_s = 0.1,
    tau2_B = 0.1,
    scenario = "constant"
  ),
  fits_dir = fits_dir,
  mesh_vertices = mesh$n,
  point_process_method = "mesh_dual",
  areal_method = "vertices",
  separate_noise = FALSE,
  prior_gamma = c(1, 10)
)

run_info


# Save simulation results ------------------------------------------------------

# Use the same run identifier as the fitted model directory
output_file <- file.path(output_dir, "summary.rds")

saveRDS(
  list(
    replication_metrics = replication_metrics,
    fit_log = fit_log,
    fit_status = fit_status,
    summary_results = summary_results,
    run_info = run_info
  ),
  file = output_file
)

message("Saved pilot results to: ", output_file)
