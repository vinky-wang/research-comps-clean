source("functions/generate_data.R")
source("functions/fit_models.R")


# Set up extension scenarios ---------------------------------------------------

# Preferential sampling mechanisms (constant-gamma baseline omitted)
scenarios <- c("spatially_varying", "independent_field")

# Areal measurement noise (baseline tau2_B = 0.1 omitted)
tau2_B_values <- c(0.025, 0.05, 0.2, 0.4)

# Preferential sampling strength (baseline gamma = 1 omitted)
gamma_values <- c(0, 0.5, 1.5, 2)

# Keep the same mesh, models, and fitting specifications
model_names <- c("ps_melding", "melding", "ps_geo")
noise_models <- c("ps_melding", "melding")
mesh <- create_mesh()

n_replications <- 100L
replications <- seq_len(2L)
simulation_seeds <- 123L + seq_len(n_replications) - 1L

output_dir <- file.path("results", "extensions")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# Record simulation settings before fitting -----------------------------------

run_info <- list(
  r_version = R.version.string,
  inla_version = INLA::inla.version("version"),
  session_info = utils::sessionInfo(),
  scenarios = scenarios,
  models = model_names,
  n_replications = n_replications,
  simulation_seeds = simulation_seeds,
  generating = list(
    rho = 0.2,
    sigma2 = 1,
    gamma = 1,
    I = 100,
    J = 25,
    tau2_s = 0.1,
    tau2_B = 0.1,
    delta = 1,
    beta_psi = 1,
    rho_psi = 0.2
  ),
  measurement_noise = list(
    scenario = "constant",
    tau2_B_values = tau2_B_values,
    ps_geo_reused = TRUE,
    ps_geo_source_tau2_B = tau2_B_values[1]
  ),
  gamma_strength = list(
    scenario = "constant",
    gamma_values = gamma_values
  ),
  mesh_vertices = mesh$n,
  point_process_method = "mesh_dual",
  areal_method = "vertices",
  separate_noise = TRUE,
  prior_gamma = c(1, 10)
)

saveRDS(run_info, file = file.path(output_dir, "run_info.rds"))

run_info$assigned_replications <- replications
saveRDS(run_info, file = file.path(output_dir, "run_info_laptop1.rds"))


# Vary preferential sampling mechanisms ---------------------------------------

for (scenario in scenarios) {
  message("Scenario: ", scenario)
  
  scenario_dir <- file.path(output_dir, scenario)
  
  for (replication in replications) {
    message("  Replication ", replication, " of ", n_replications)
    
    simulated <- generate_data(
      rho = 0.2,
      gamma = 1,
      I = 100,
      J = 25,
      tau2_s = 0.1,
      tau2_B = 0.1,
      seed = simulation_seeds[replication],
      scenario = scenario,
      delta = 1,
      beta_psi = 1,
      rho_psi = 0.2
    )
    
    logged <- fit_models_logged(
      data = simulated,
      models = model_names,
      mesh = mesh,
      replication = replication,
      seed = simulation_seeds[replication],
      point_process_method = "mesh_dual",
      areal_method = "vertices",
      separate_noise = TRUE,
      prior_gamma = c(1, 10)
    )
    
    # Save fits and log before starting the next replication
    save_replication(
      data = simulated,
      results = logged$results,
      log = logged$log,
      output_dir = scenario_dir,
      replication = replication
    )
    
    rm(simulated, logged)
  }
}


# Vary areal measurement noise ------------------------------------------------

noise_output_dir <- file.path(output_dir, "measurement_noise")

# Fit PS geo once at the first noise level, then reuse it at the others
for (replication in replications) {
  message("Measurement-noise replication ", replication, " of ", n_replications)
  
  seed <- simulation_seeds[replication]
  reference <- NULL
  logged_geo <- NULL
  
  for (k in seq_along(tau2_B_values)) {
    tau2_B <- tau2_B_values[k]
    message("  Areal noise variance: ", tau2_B)
    
    simulated <- generate_data(
      rho = 0.2,
      gamma = 1,
      I = 100,
      J = 25,
      tau2_s = 0.1,
      tau2_B = tau2_B,
      seed = seed,
      scenario = "constant"
    )
    
    if (k == 1L) {
      reference <- simulated
    } else {
      # Changing areal noise must leave the point data and truth unchanged
      stopifnot(
        isTRUE(all.equal(simulated$points, reference$points)),
        isTRUE(all.equal(simulated$prediction, reference$prediction))
      )
    }
    
    models_to_fit <- if (k == 1L) model_names else noise_models
    
    logged <- fit_models_logged(
      data = simulated,
      models = models_to_fit,
      mesh = mesh,
      replication = replication,
      seed = seed,
      point_process_method = "mesh_dual",
      areal_method = "vertices",
      separate_noise = TRUE,
      prior_gamma = c(1, 10)
    )
    
    if (k == 1L) {
      # Retain the point-only fit for subsequent noise levels
      logged_geo <- list(
        models = logged$results$models["ps_geo"],
        spde = logged$results$spde,
        log = logged$log[logged$log$model == "ps_geo", , drop = FALSE]
      )
      
      results <- logged$results
      fit_log <- logged$log
      fit_log$reused <- FALSE
    } else {
      results <- list(
        models = c(logged$results$models, logged_geo$models),
        mesh = mesh,
        spde = logged$results$spde
      )
      
      if (is.null(results$spde)) {
        results$spde <- logged_geo$spde
      }
      
      fit_log <- rbind(logged$log, logged_geo$log)
      fit_log$reused <- fit_log$model == "ps_geo"
    }
    
    fit_log$tau2_B <- tau2_B
    rownames(fit_log) <- NULL
    
    condition_dir <- file.path(noise_output_dir, paste0("tau2_B_", tau2_B))
    
    save_replication(
      data = simulated,
      results = results,
      log = fit_log,
      output_dir = condition_dir,
      replication = replication
    )
    
    rm(simulated, logged, results, fit_log)
  }
  
  rm(reference, logged_geo)
}


# Vary preferential sampling strength ----------------------------------------

gamma_output_dir <- file.path(output_dir, "gamma_strength")

for (replication in replications) {
  message("Gamma-strength replication ", replication, " of ", n_replications)
  
  seed <- simulation_seeds[replication]
  
  # Reference field and areal observations are shared across gamma values
  reference <- generate_data(
    rho = 0.2,
    gamma = 1,
    I = 100,
    J = 25,
    tau2_s = 0.1,
    tau2_B = 0.1,
    seed = seed,
    scenario = "constant"
  )
  
  for (gamma_value in gamma_values) {
    message("  Preferentiality parameter: ", gamma_value)
    
    simulated <- generate_data(
      rho = 0.2,
      gamma = gamma_value,
      I = 100,
      J = 25,
      tau2_s = 0.1,
      tau2_B = 0.1,
      seed = seed,
      scenario = "constant"
    )
    
    stopifnot(
      isTRUE(all.equal(simulated$prediction, reference$prediction)),
      isTRUE(all.equal(simulated$areas, reference$areas))
    )
    
    logged <- fit_models_logged(
      data = simulated,
      models = model_names,
      mesh = mesh,
      replication = replication,
      seed = seed,
      point_process_method = "mesh_dual",
      areal_method = "vertices",
      separate_noise = TRUE,
      prior_gamma = c(1, 10)
    )
    
    fit_log <- logged$log
    fit_log$gamma_true <- gamma_value
    
    condition_dir <- file.path(gamma_output_dir, paste0("gamma_", gamma_value))
    
    save_replication(
      data = simulated,
      results = logged$results,
      log = fit_log,
      output_dir = condition_dir,
      replication = replication
    )
    
    rm(simulated, logged, fit_log)
  }
  
  rm(reference)
}

message("Finished fitting. Run scripts/summarize_results.R for summaries.")
