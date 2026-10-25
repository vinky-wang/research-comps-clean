source("functions/generate_data.R")
source("functions/fit_models.R")


# Set up reproduction scenarios ------------------------------------------------

# Number of observations
I_values <- c(100, 250)
J_values <- c(4, 25, 100)

# Spatial range
rho_values <- c(0.1, 0.2, 0.4)

# Preferential sampling strength
gamma_values <- c(0, 1)

# All simulation configurations
scenarios <- expand.grid(
  rho = rho_values,
  gamma = gamma_values,
  I = I_values,
  J = J_values
)

# Keep the same mesh and model specifications across replications
model_names <- c("ps_melding", "melding", "ps_geo")
mesh <- create_mesh()

n_replications <- 2L
simulation_seeds <- 123L + seq_len(n_replications) - 1L

# Set output directory
# Keep the full study separate from existing pilot results
output_dir <- file.path("results", "reproduction", "main")
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
    sigma2 = 1,
    nu = 1,
    tau2_s = 0.1,
    tau2_B = 0.1,
    scenario = "constant"
  ),
  mesh_vertices = mesh$n,
  point_process_method = "mesh_dual",
  areal_method = "vertices",
  separate_noise = TRUE,
  prior_gamma = c(1, 10)
)

saveRDS(
  run_info,
  file = file.path(output_dir, "run_info.rds")
)


# Fit and save models across scenarios ----------------------------------------

for (k in seq_len(nrow(scenarios))) {
  scenario <- scenarios[k, ]
  
  rho <- scenario$rho
  gamma <- scenario$gamma
  I <- scenario$I
  J <- scenario$J
  
  message(
    "Scenario ", k, " of ", nrow(scenarios),
    ": rho = ", rho,
    ", gamma = ", gamma,
    ", I = ", I,
    ", J = ", J
  )
  
  # Separate directory for each configuration
  scenario_dir <- file.path(
    output_dir,
    paste0("rho_", rho),
    paste0("gamma_", gamma),
    paste0("I_", I, "_J_", J)
  )
  
  for (replication in seq_len(n_replications)) {
    message("  Replication ", replication, " of ", n_replications)
    
    # Generate simulated data
    simulated <- generate_data(
      rho = rho,
      gamma = gamma,
      I = I,
      J = J,
      tau2_s = 0.1,
      tau2_B = 0.1,
      seed = simulation_seeds[replication],
      scenario = "constant"
    )
    
    # Fit and log all three models
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
    
    # Record the generating configuration in the fitting log
    logged$log$rho <- rho
    logged$log$gamma_true <- gamma
    logged$log$I <- I
    logged$log$J <- J
    
    # Save fitted models and log immediately after each replication
    save_replication(
      data = simulated,
      results = logged$results,
      log = logged$log,
      output_dir = scenario_dir,
      replication = replication
    )
    
    # Release large fitted objects before the next replication
    rm(simulated, logged)
  }
}

message("Finished fitting. Run scripts/summarize_results.R for summaries.")