source("functions/generate_data.R")
source("functions/fit_models.R")
source("functions/create_summary.R")


# Prior sensitivity ------------------------------------------------------------

# Common fitting specifications
mesh = create_mesh(max_edge = c(0.05, 0.25))
areal_method = "vertices"


## Areal-noise ---------------------------------------------------------------

# Simulation settings
prior_replications = seq_len(10L)
prior_tau2_B_values = c(0.025, 0.05, 0.1, 0.2, 0.4)
prior_models = c("melding", "ps_melding")

# PC priors: P(sigma_B > U) = 0.05
prior_specs = list(
  default = NULL,
  pc_U0.5 = c(0.5, 0.05),
  pc_U1 = c(1, 0.05)
)

prior_output_dir = file.path("results", "prerun", "tau2_B_prior_sensitivity")
dir.create(prior_output_dir, recursive = TRUE, showWarnings = FALSE)

# Record experiment settings
saveRDS(
  list(
    replications = prior_replications,
    tau2_B_values = prior_tau2_B_values,
    models = prior_models,
    priors = prior_specs,
    mesh_vertices = mesh$n,
    areal_method = areal_method,
    prior_noise = NULL,
    prior_gamma = c(1, 10)
  ),
  file = file.path(prior_output_dir, "run_info.rds")
)


# Fit models 

for (replication in prior_replications) {
  seed = 123L + replication - 1L
  
  for (tau2_B in prior_tau2_B_values) {
    
    # Generate data
    simulated = generate_data(
      rho = 0.2,
      gamma = 1,
      I = 100L,
      J = 25L,
      tau2_s = 0.1,
      tau2_B = tau2_B,
      seed = seed,
      scenario = "constant"
    )
    
    A_B = create_projections(simulated, mesh, areal_method)$areas
    
    for (prior_name in names(prior_specs)) {
      
      condition_dir = file.path(
        prior_output_dir,
        paste0("tau2_B_", tau2_B),
        prior_name
      )
      dir.create(condition_dir, recursive = TRUE, showWarnings = FALSE)
      
      fit_path = file.path(
        condition_dir,
        sprintf("rep_%03d.rds", replication)
      )
      
      if (file.exists(fit_path)) next
      
      message(
        "Replication ", replication,
        ", tau2_B = ", tau2_B,
        ", prior = ", prior_name
      )
      
      logged = fit_models_logged(
        data = simulated,
        models = prior_models,
        mesh = mesh,
        replication = replication,
        seed = seed,
        point_process_method = "mesh_dual",
        areal_method = areal_method,
        separate_noise = TRUE,
        prior_gamma = c(1, 10),
        prior_noise = NULL,
        prior_noise_B = prior_specs[[prior_name]]
      )
      
      diagnostics = do.call(rbind, lapply(prior_models, function(model) {
        summarise_fit(
          data = simulated,
          logged = logged,
          model = model,
          mesh = mesh,
          A_B = A_B,
          replication = replication,
          tau2_B = tau2_B,
          mesh_name = "standard",
          areal_method = areal_method,
          prior_name = prior_name
        )
      }))
      
      saveRDS(
        list(
          data = simulated,
          results = logged$results,
          log = logged$log,
          diagnostics = diagnostics
        ),
        file = fit_path
      )
      
      rm(logged, diagnostics)
    }
  }
}


# Summarize results 

fit_files = list.files(
  prior_output_dir,
  pattern = "^rep_[0-9]+\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

prior_diagnostics = do.call(rbind, lapply(fit_files, function(path) {
  readRDS(path)$diagnostics
}))
rownames(prior_diagnostics) = NULL

# Absolute log-error in areal-noise variance
prior_diagnostics$abs_log_tau2_B_error = ifelse(
  is.finite(prior_diagnostics$tau2_B_median) &
    prior_diagnostics$tau2_B_median > 0,
  abs(log(
    prior_diagnostics$tau2_B_median /
      prior_diagnostics$tau2_B_true
  )),
  NA_real_
)

write.csv(
  prior_diagnostics,
  file.path(prior_output_dir, "fit_diagnostics.csv"),
  row.names = FALSE
)

# Summarize by true variance, prior, and model
groups = split(
  prior_diagnostics,
  interaction(
    prior_diagnostics$tau2_B_true,
    prior_diagnostics$areal_prior,
    prior_diagnostics$model,
    drop = TRUE
  )
)

prior_summary = do.call(rbind, lapply(groups, function(d) {
  data.frame(
    tau2_B_true = d$tau2_B_true[1],
    prior = d$areal_prior[1],
    model = d$model[1],
    n_attempted = nrow(d),
    n_included = sum(d$included),
    n_near_zero = sum(d$near_zero_tau2_B, na.rm = TRUE),
    median_tau2_B = safe_median(d$tau2_B_median),
    median_abs_log_error = safe_median(d$abs_log_tau2_B_error),
    tau2_B_coverage = mean(d$tau2_B_covered, na.rm = TRUE),
    median_prediction_mse = safe_median(d$prediction_mse),
    median_prediction_coverage = safe_median(d$prediction_coverage)
  )
}))
rownames(prior_summary) = NULL

write.csv(
  prior_summary,
  file.path(prior_output_dir, "comparison_summary.csv"),
  row.names = FALSE
)

print(prior_summary)
message("Areal-noise prior sensitivity saved to: ", prior_output_dir)


baseline = subset(prior_diagnostics, tau2_B_true == 0.1)

# Point-error variance
print(aggregate(
  tau2_s_median ~ model + areal_prior,
  data = baseline,
  FUN = median
))

# Preferentiality coefficient
print(aggregate(
  gamma_median ~ areal_prior,
  data = subset(baseline, model == "ps_melding"),
  FUN = median
))


## Varying preferentiality ----------------------------------------------------

# Simulation settings
gamma_replications = seq_len(10L)
gamma_values = c(0, 0.5, 1, 1.5, 2)

# Normal priors: mean and precision
gamma_priors = list(
  original = c(1, 10),
  weak_center1 = c(1, 1),
  weak_center0 = c(0, 1)
)

# Fixed areal-noise prior
gamma_prior_noise_B = c(1, 0.05)

gamma_output_dir = file.path(
  "results", "prerun", "gamma_prior_sensitivity"
)
dir.create(gamma_output_dir, recursive = TRUE, showWarnings = FALSE)

# Record experiment settings
saveRDS(
  list(
    replications = gamma_replications,
    gamma_values = gamma_values,
    priors = gamma_priors,
    models = c("melding", "ps_melding"),
    mesh_vertices = mesh$n,
    areal_method = areal_method,
    prior_noise = NULL,
    prior_noise_B = gamma_prior_noise_B
  ),
  file = file.path(gamma_output_dir, "run_info.rds")
)


# Extract preferentiality diagnostics 

summarise_gamma_fit = function(data, logged, model, replication,
                               gamma_true, gamma_prior, A_B) {
  
  diagnostic = summarise_fit(
    data = data,
    logged = logged,
    model = model,
    mesh = mesh,
    A_B = A_B,
    replication = replication,
    tau2_B = 0.1,
    mesh_name = "standard",
    areal_method = areal_method,
    prior_name = "pc_U1"
  )
  
  diagnostic$gamma_true = gamma_true
  diagnostic$gamma_prior = gamma_prior
  diagnostic$gamma_sd = NA_real_
  diagnostic$gamma_q025 = NA_real_
  diagnostic$gamma_q975 = NA_real_
  diagnostic$gamma_covered = NA
  
  if (model == "ps_melding") {
    fitted = logged$results$models[[model]]
    
    if (!is.null(fitted)) {
      hp = fitted$fit$summary.hyperpar
      
      if ("Beta for j" %in% rownames(hp) &&
          all(c("sd", "0.025quant", "0.975quant") %in% colnames(hp))) {
        
        diagnostic$gamma_sd = hp["Beta for j", "sd"]
        diagnostic$gamma_q025 = hp["Beta for j", "0.025quant"]
        diagnostic$gamma_q975 = hp["Beta for j", "0.975quant"]
        
        diagnostic$gamma_covered =
          diagnostic$gamma_q025 <= gamma_true &&
          gamma_true <= diagnostic$gamma_q975
      }
    }
  }
  
  return(diagnostic)
}


# Fit models 

for (replication in gamma_replications) {
  seed = 123L + replication - 1L
  
  for (gamma_true in gamma_values) {
    
    # Generate data
    simulated = generate_data(
      rho = 0.2,
      gamma = gamma_true,
      I = 100L,
      J = 25L,
      tau2_s = 0.1,
      tau2_B = 0.1,
      seed = seed,
      scenario = "constant"
    )
    
    A_B = create_projections(simulated, mesh, areal_method)$areas
    
    condition_dir = file.path(
      gamma_output_dir,
      paste0("gamma_", gamma_true)
    )
    
    # Fit Melding once per dataset
    melding_path = file.path(
      condition_dir,
      "melding",
      sprintf("rep_%03d.rds", replication)
    )
    
    if (!file.exists(melding_path)) {
      dir.create(
        dirname(melding_path),
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      message(
        "Replication ", replication,
        ", gamma = ", gamma_true,
        ", model = melding"
      )
      
      logged = fit_models_logged(
        data = simulated,
        models = "melding",
        mesh = mesh,
        replication = replication,
        seed = seed,
        point_process_method = "mesh_dual",
        areal_method = areal_method,
        separate_noise = TRUE,
        prior_noise = NULL,
        prior_noise_B = gamma_prior_noise_B
      )
      
      diagnostics = summarise_gamma_fit(
        data = simulated,
        logged = logged,
        model = "melding",
        replication = replication,
        gamma_true = gamma_true,
        gamma_prior = "not_applicable",
        A_B = A_B
      )
      
      saveRDS(
        list(
          data = simulated,
          results = logged$results,
          log = logged$log,
          diagnostics = diagnostics
        ),
        file = melding_path
      )
      
      rm(logged, diagnostics)
    }
    
    # Fit PS Melding under each gamma prior
    for (prior_name in names(gamma_priors)) {
      
      fit_path = file.path(
        condition_dir,
        "ps_melding",
        prior_name,
        sprintf("rep_%03d.rds", replication)
      )
      
      dir.create(
        dirname(fit_path),
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      if (file.exists(fit_path)) next
      
      message(
        "Replication ", replication,
        ", gamma = ", gamma_true,
        ", prior = ", prior_name
      )
      
      logged = fit_models_logged(
        data = simulated,
        models = "ps_melding",
        mesh = mesh,
        replication = replication,
        seed = seed,
        point_process_method = "mesh_dual",
        areal_method = areal_method,
        separate_noise = TRUE,
        prior_gamma = gamma_priors[[prior_name]],
        prior_noise = NULL,
        prior_noise_B = gamma_prior_noise_B
      )
      
      diagnostics = summarise_gamma_fit(
        data = simulated,
        logged = logged,
        model = "ps_melding",
        replication = replication,
        gamma_true = gamma_true,
        gamma_prior = prior_name,
        A_B = A_B
      )
      
      saveRDS(
        list(
          data = simulated,
          results = logged$results,
          log = logged$log,
          diagnostics = diagnostics
        ),
        file = fit_path
      )
      
      rm(logged, diagnostics)
    }
  }
}


# Summarize results 

fit_files = list.files(
  gamma_output_dir,
  pattern = "^rep_[0-9]+\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

gamma_diagnostics = do.call(rbind, lapply(fit_files, function(path) {
  readRDS(path)$diagnostics
}))
rownames(gamma_diagnostics) = NULL

# Compare prediction MSE against Melding on the same dataset
melding_reference = subset(
  gamma_diagnostics,
  model == "melding",
  select = c("replication", "gamma_true", "prediction_mse")
)

names(melding_reference)[3] = "melding_prediction_mse"

gamma_diagnostics = merge(
  gamma_diagnostics,
  melding_reference,
  by = c("replication", "gamma_true"),
  all.x = TRUE
)

gamma_diagnostics$delta_prediction_mse =
  gamma_diagnostics$prediction_mse -
  gamma_diagnostics$melding_prediction_mse

write.csv(
  gamma_diagnostics,
  file.path(gamma_output_dir, "fit_diagnostics.csv"),
  row.names = FALSE
)

# Summarize by true gamma, prior, and model
groups = split(
  gamma_diagnostics,
  interaction(
    gamma_diagnostics$gamma_true,
    gamma_diagnostics$gamma_prior,
    gamma_diagnostics$model,
    drop = TRUE
  )
)

gamma_summary = do.call(rbind, lapply(groups, function(d) {
  data.frame(
    gamma_true = d$gamma_true[1],
    prior = d$gamma_prior[1],
    model = d$model[1],
    n_attempted = nrow(d),
    n_included = sum(d$included),
    median_gamma = safe_median(d$gamma_median),
    gamma_rmse = if (all(is.na(d$gamma_median))) NA_real_ else
      sqrt(mean((d$gamma_median - d$gamma_true)^2, na.rm = TRUE)),
    median_gamma_sd = safe_median(d$gamma_sd),
    gamma_coverage = if (all(is.na(d$gamma_covered))) NA_real_ else
      mean(d$gamma_covered, na.rm = TRUE),
    median_tau2_s = safe_median(d$tau2_s_median),
    median_tau2_B = safe_median(d$tau2_B_median),
    median_prediction_mse = safe_median(d$prediction_mse),
    median_prediction_coverage = safe_median(d$prediction_coverage),
    median_delta_mse = safe_median(d$delta_prediction_mse)
  )
}))
rownames(gamma_summary) = NULL

write.csv(
  gamma_summary,
  file.path(gamma_output_dir, "comparison_summary.csv"),
  row.names = FALSE
)

print(gamma_summary)
message("Gamma prior sensitivity saved to: ", gamma_output_dir)




# Mesh sensitivity ------------------------------------------------------------

# Simulation settings
mesh_replications = seq_len(10L)
mesh_gamma_values = c(1, 2)
mesh_models = c("melding", "ps_melding")

# Mesh resolutions
mesh_specs = list(
  standard = c(0.05, 0.25),
  finer = c(0.035, 0.175),
  finest = c(0.025, 0.125)
)

meshes = lapply(mesh_specs, create_mesh)

# Keep other fitting specifications fixed
mesh_areal_method = "vertices"
mesh_prior_gamma = c(1, 1)
mesh_prior_noise_B = c(1, 0.05)

mesh_output_dir = file.path("results", "prerun", "mesh_sensitivity")
dir.create(mesh_output_dir, recursive = TRUE, showWarnings = FALSE)

# Record experiment settings
saveRDS(
  list(
    replications = mesh_replications,
    gamma_values = mesh_gamma_values,
    models = mesh_models,
    mesh_specs = mesh_specs,
    mesh_vertices = vapply(meshes, function(x) x$n, numeric(1)),
    areal_method = mesh_areal_method,
    prior_gamma = mesh_prior_gamma,
    prior_noise = NULL,
    prior_noise_B = mesh_prior_noise_B
  ),
  file = file.path(mesh_output_dir, "run_info.rds")
)


# Check areal projections ------------------------------------------------------

# Use one dataset to compare projection accuracy
reference = generate_data(
  rho = 0.2,
  gamma = 1,
  I = 100L,
  J = 25L,
  tau2_s = 0.1,
  tau2_B = 0.1,
  seed = 123L,
  scenario = "constant"
)

areas = reference$areas

# Exact area averages of x + y and x^2 + y^2
linear_true = (
  areas$xmin + areas$xmax +
    areas$ymin + areas$ymax
) / 2

quadratic_true = (
  areas$xmin^2 + areas$xmin * areas$xmax + areas$xmax^2 +
    areas$ymin^2 + areas$ymin * areas$ymax + areas$ymax^2
) / 3

projection_checks = do.call(rbind, lapply(names(meshes), function(mesh_name) {
  current_mesh = meshes[[mesh_name]]
  A_B = create_projections(
    reference, current_mesh, mesh_areal_method
  )$areas
  
  xy = current_mesh$loc[, 1:2, drop = FALSE]
  
  linear_estimate = as.numeric(A_B %*% (xy[, 1] + xy[, 2]))
  quadratic_estimate = as.numeric(
    A_B %*% (xy[, 1]^2 + xy[, 2]^2)
  )
  
  data.frame(
    mesh = mesh_name,
    mesh_vertices = current_mesh$n,
    max_row_sum_error = max(abs(Matrix::rowSums(A_B) - 1)),
    linear_rmse = sqrt(mean((linear_estimate - linear_true)^2)),
    quadratic_rmse = sqrt(mean((quadratic_estimate - quadratic_true)^2))
  )
}))

write.csv(
  projection_checks,
  file.path(mesh_output_dir, "projection_checks.csv"),
  row.names = FALSE
)

print(projection_checks)


# Fit models -------------------------------------------------------------------

for (replication in mesh_replications) {
  seed = 123L + replication - 1L
  
  for (gamma_true in mesh_gamma_values) {
    
    # Generate data once for each replication and gamma
    simulated = generate_data(
      rho = 0.2,
      gamma = gamma_true,
      I = 100L,
      J = 25L,
      tau2_s = 0.1,
      tau2_B = 0.1,
      seed = seed,
      scenario = "constant"
    )
    
    for (mesh_name in names(meshes)) {
      current_mesh = meshes[[mesh_name]]
      
      fit_path = file.path(
        mesh_output_dir,
        paste0("gamma_", gamma_true),
        mesh_name,
        sprintf("rep_%03d.rds", replication)
      )
      
      dir.create(
        dirname(fit_path),
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      if (file.exists(fit_path)) next
      
      message(
        "Replication ", replication,
        ", gamma = ", gamma_true,
        ", mesh = ", mesh_name
      )
      
      A_B = create_projections(
        simulated, current_mesh, mesh_areal_method
      )$areas
      
      logged = fit_models_logged(
        data = simulated,
        models = mesh_models,
        mesh = current_mesh,
        replication = replication,
        seed = seed,
        point_process_method = "mesh_dual",
        areal_method = mesh_areal_method,
        separate_noise = TRUE,
        prior_gamma = mesh_prior_gamma,
        prior_noise = NULL,
        prior_noise_B = mesh_prior_noise_B
      )
      
      diagnostics = do.call(rbind, lapply(mesh_models, function(model) {
        
        diagnostic = summarise_fit(
          data = simulated,
          logged = logged,
          model = model,
          mesh = current_mesh,
          A_B = A_B,
          replication = replication,
          tau2_B = 0.1,
          mesh_name = mesh_name,
          areal_method = mesh_areal_method,
          prior_name = "pc_U1"
        )
        
        diagnostic$gamma_true = gamma_true
        diagnostic$gamma_q025 = NA_real_
        diagnostic$gamma_q975 = NA_real_
        diagnostic$gamma_covered = NA
        
        if (model == "ps_melding") {
          fitted = logged$results$models[[model]]
          
          if (!is.null(fitted)) {
            hp = fitted$fit$summary.hyperpar
            
            if ("Beta for j" %in% rownames(hp) &&
                all(c("0.025quant", "0.975quant") %in% colnames(hp))) {
              
              diagnostic$gamma_q025 = hp["Beta for j", "0.025quant"]
              diagnostic$gamma_q975 = hp["Beta for j", "0.975quant"]
              
              diagnostic$gamma_covered =
                diagnostic$gamma_q025 <= gamma_true &&
                gamma_true <= diagnostic$gamma_q975
            }
          }
        }
        
        return(diagnostic)
      }))
      
      saveRDS(
        list(
          data = simulated,
          results = logged$results,
          log = logged$log,
          diagnostics = diagnostics
        ),
        file = fit_path
      )
      
      rm(logged, diagnostics)
    }
  }
}


# Summarize results ------------------------------------------------------------

fit_files = list.files(
  mesh_output_dir,
  pattern = "^rep_[0-9]+\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

mesh_diagnostics = do.call(rbind, lapply(fit_files, function(path) {
  readRDS(path)$diagnostics
}))
rownames(mesh_diagnostics) = NULL

# Compare each mesh against the standard mesh
standard_reference = subset(
  mesh_diagnostics,
  mesh == "standard",
  select = c(
    "replication", "gamma_true", "model",
    "prediction_mse", "gamma_median"
  )
)

names(standard_reference)[4:5] = c(
  "standard_prediction_mse",
  "standard_gamma_median"
)

mesh_diagnostics = merge(
  mesh_diagnostics,
  standard_reference,
  by = c("replication", "gamma_true", "model"),
  all.x = TRUE
)

mesh_diagnostics$delta_prediction_mse =
  mesh_diagnostics$prediction_mse -
  mesh_diagnostics$standard_prediction_mse

mesh_diagnostics$delta_gamma_abs_error =
  abs(mesh_diagnostics$gamma_median - mesh_diagnostics$gamma_true) -
  abs(mesh_diagnostics$standard_gamma_median - mesh_diagnostics$gamma_true)

write.csv(
  mesh_diagnostics,
  file.path(mesh_output_dir, "fit_diagnostics.csv"),
  row.names = FALSE
)

# Summarize by true gamma, mesh, and model
groups = split(
  mesh_diagnostics,
  interaction(
    mesh_diagnostics$gamma_true,
    mesh_diagnostics$mesh,
    mesh_diagnostics$model,
    drop = TRUE
  )
)

mesh_summary = do.call(rbind, lapply(groups, function(d) {
  data.frame(
    gamma_true = d$gamma_true[1],
    mesh = d$mesh[1],
    mesh_vertices = d$mesh_vertices[1],
    model = d$model[1],
    n_attempted = nrow(d),
    n_included = sum(d$included),
    median_gamma = safe_median(d$gamma_median),
    gamma_rmse = if (all(is.na(d$gamma_median))) NA_real_ else
      sqrt(mean((d$gamma_median - d$gamma_true)^2, na.rm = TRUE)),
    gamma_coverage = if (all(is.na(d$gamma_covered))) NA_real_ else
      mean(d$gamma_covered, na.rm = TRUE),
    median_tau2_s = safe_median(d$tau2_s_median),
    median_tau2_B = safe_median(d$tau2_B_median),
    tau2_B_coverage = mean(d$tau2_B_covered, na.rm = TRUE),
    median_prediction_mse = safe_median(d$prediction_mse),
    median_prediction_coverage = safe_median(d$prediction_coverage),
    median_delta_mse = safe_median(d$delta_prediction_mse),
    median_delta_gamma_error = safe_median(d$delta_gamma_abs_error),
    median_elapsed_sec = safe_median(d$elapsed_sec)
  )
}))
rownames(mesh_summary) = NULL

write.csv(
  mesh_summary,
  file.path(mesh_output_dir, "comparison_summary.csv"),
  row.names = FALSE
)

print(mesh_summary)
message("Mesh sensitivity results saved to: ", mesh_output_dir)



# Areal-method sensitivity -----------------------------------------------------

# Simulation settings
areal_replications = mesh_replications
areal_gamma_values = mesh_gamma_values
areal_models = mesh_models
areal_methods = c("vertices", "grid")

# Use the same meshes and priors
areal_meshes = meshes
areal_prior_gamma = mesh_prior_gamma
areal_prior_noise_B = mesh_prior_noise_B

areal_output_dir = file.path(
  "results", "prerun", "areal_method_sensitivity"
)
dir.create(areal_output_dir, recursive = TRUE, showWarnings = FALSE)

# Record experiment settings
saveRDS(
  list(
    replications = areal_replications,
    gamma_values = areal_gamma_values,
    models = areal_models,
    methods = areal_methods,
    mesh_specs = mesh_specs,
    mesh_vertices = vapply(areal_meshes, function(x) x$n, numeric(1)),
    prior_gamma = areal_prior_gamma,
    prior_noise = NULL,
    prior_noise_B = areal_prior_noise_B,
    vertices_results_dir = mesh_output_dir
  ),
  file = file.path(areal_output_dir, "run_info.rds")
)


# Check areal projections ------------------------------------------------------

# Compare both methods using the same reference dataset
areal_projection_checks = do.call(rbind, lapply(
  names(areal_meshes), function(mesh_name) {
    
    current_mesh = areal_meshes[[mesh_name]]
    xy = current_mesh$loc[, 1:2, drop = FALSE]
    
    do.call(rbind, lapply(areal_methods, function(method) {
      
      A_B = create_projections(
        reference, current_mesh, method
      )$areas
      
      linear_estimate = as.numeric(
        A_B %*% (xy[, 1] + xy[, 2])
      )
      
      quadratic_estimate = as.numeric(
        A_B %*% (xy[, 1]^2 + xy[, 2]^2)
      )
      
      data.frame(
        mesh = mesh_name,
        mesh_vertices = current_mesh$n,
        areal_method = method,
        max_row_sum_error = max(
          abs(Matrix::rowSums(A_B) - 1)
        ),
        linear_rmse = sqrt(mean(
          (linear_estimate - linear_true)^2
        )),
        quadratic_rmse = sqrt(mean(
          (quadratic_estimate - quadratic_true)^2
        ))
      )
    }))
  }
))
rownames(areal_projection_checks) = NULL

write.csv(
  areal_projection_checks,
  file.path(areal_output_dir, "projection_checks.csv"),
  row.names = FALSE
)

print(areal_projection_checks)


# Fit models -------------------------------------------------------------------

# The vertices fits were saved during mesh sensitivity.
# Only fit models using grid-based areal integration.

for (replication in areal_replications) {
  seed = 123L + replication - 1L
  
  for (gamma_true in areal_gamma_values) {
    
    # Generate the same data used in mesh sensitivity
    simulated = generate_data(
      rho = 0.2,
      gamma = gamma_true,
      I = 100L,
      J = 25L,
      tau2_s = 0.1,
      tau2_B = 0.1,
      seed = seed,
      scenario = "constant"
    )
    
    for (mesh_name in names(areal_meshes)) {
      current_mesh = areal_meshes[[mesh_name]]
      
      fit_path = file.path(
        areal_output_dir,
        paste0("gamma_", gamma_true),
        mesh_name,
        "grid",
        sprintf("rep_%03d.rds", replication)
      )
      
      dir.create(
        dirname(fit_path),
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      if (file.exists(fit_path)) next
      
      message(
        "Replication ", replication,
        ", gamma = ", gamma_true,
        ", mesh = ", mesh_name,
        ", areal method = grid"
      )
      
      # Grid-based areal projection
      A_B = create_projections(
        simulated, current_mesh, "grid"
      )$areas
      
      logged = fit_models_logged(
        data = simulated,
        models = areal_models,
        mesh = current_mesh,
        replication = replication,
        seed = seed,
        point_process_method = "mesh_dual",
        areal_method = "grid",
        separate_noise = TRUE,
        prior_gamma = areal_prior_gamma,
        prior_noise = NULL,
        prior_noise_B = areal_prior_noise_B
      )
      
      diagnostics = do.call(rbind, lapply(areal_models, function(model) {
        
        diagnostic = summarise_fit(
          data = simulated,
          logged = logged,
          model = model,
          mesh = current_mesh,
          A_B = A_B,
          replication = replication,
          tau2_B = 0.1,
          mesh_name = mesh_name,
          areal_method = "grid",
          prior_name = "pc_U1"
        )
        
        diagnostic$gamma_true = gamma_true
        diagnostic$gamma_q025 = NA_real_
        diagnostic$gamma_q975 = NA_real_
        diagnostic$gamma_covered = NA
        
        if (model == "ps_melding") {
          fitted = logged$results$models[[model]]
          
          if (!is.null(fitted)) {
            hp = fitted$fit$summary.hyperpar
            
            if ("Beta for j" %in% rownames(hp) &&
                all(c("0.025quant", "0.975quant") %in% colnames(hp))) {
              
              diagnostic$gamma_q025 = hp["Beta for j", "0.025quant"]
              diagnostic$gamma_q975 = hp["Beta for j", "0.975quant"]
              
              diagnostic$gamma_covered =
                diagnostic$gamma_q025 <= gamma_true &&
                gamma_true <= diagnostic$gamma_q975
            }
          }
        }
        
        return(diagnostic)
      }))
      
      saveRDS(
        list(
          data = simulated,
          results = logged$results,
          log = logged$log,
          diagnostics = diagnostics
        ),
        file = fit_path
      )
      
      rm(logged, diagnostics)
    }
  }
}


# Summarize results ------------------------------------------------------------

# Existing vertices fits
vertex_files = list.files(
  mesh_output_dir,
  pattern = "^rep_[0-9]+\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

# New grid fits
grid_files = list.files(
  areal_output_dir,
  pattern = "^rep_[0-9]+\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(vertex_files) == 0L || length(grid_files) == 0L) {
  stop("Both vertices and grid results are required for comparison.")
}

# Combine diagnostics from both methods
areal_diagnostics = do.call(rbind, lapply(
  c(vertex_files, grid_files),
  function(path) readRDS(path)$diagnostics
))
rownames(areal_diagnostics) = NULL

# Retain the settings considered in this experiment
areal_diagnostics = subset(
  areal_diagnostics,
  gamma_true %in% areal_gamma_values &
    mesh %in% names(areal_meshes) &
    areal_method %in% areal_methods
)

# Reference: vertices method on the same mesh and dataset
vertex_reference = subset(
  areal_diagnostics,
  areal_method == "vertices",
  select = c(
    "replication", "gamma_true", "mesh", "model",
    "prediction_mse", "gamma_median", "tau2_B_median"
  )
)

names(vertex_reference)[5:7] = c(
  "vertex_prediction_mse",
  "vertex_gamma_median",
  "vertex_tau2_B_median"
)

areal_diagnostics = merge(
  areal_diagnostics,
  vertex_reference,
  by = c("replication", "gamma_true", "mesh", "model"),
  all.x = TRUE
)

# Paired differences relative to vertices
# Negative values indicate improvement
areal_diagnostics$delta_prediction_mse =
  areal_diagnostics$prediction_mse -
  areal_diagnostics$vertex_prediction_mse

areal_diagnostics$delta_gamma_abs_error =
  abs(areal_diagnostics$gamma_median - areal_diagnostics$gamma_true) -
  abs(areal_diagnostics$vertex_gamma_median - areal_diagnostics$gamma_true)

areal_diagnostics$delta_tau2_B_abs_error =
  abs(areal_diagnostics$tau2_B_median - areal_diagnostics$tau2_B_true) -
  abs(areal_diagnostics$vertex_tau2_B_median - areal_diagnostics$tau2_B_true)

write.csv(
  areal_diagnostics,
  file.path(areal_output_dir, "fit_diagnostics.csv"),
  row.names = FALSE
)

# Summarize by true gamma, mesh, areal method, and model
groups = split(
  areal_diagnostics,
  interaction(
    areal_diagnostics$gamma_true,
    areal_diagnostics$mesh,
    areal_diagnostics$areal_method,
    areal_diagnostics$model,
    drop = TRUE
  )
)

areal_summary = do.call(rbind, lapply(groups, function(d) {
  data.frame(
    gamma_true = d$gamma_true[1],
    mesh = d$mesh[1],
    mesh_vertices = d$mesh_vertices[1],
    areal_method = d$areal_method[1],
    model = d$model[1],
    n_attempted = nrow(d),
    n_included = sum(d$included),
    
    median_gamma = safe_median(d$gamma_median),
    gamma_rmse = if (all(is.na(d$gamma_median))) NA_real_ else
      sqrt(mean((d$gamma_median - d$gamma_true)^2, na.rm = TRUE)),
    gamma_coverage = if (all(is.na(d$gamma_covered))) NA_real_ else
      mean(d$gamma_covered, na.rm = TRUE),
    
    median_tau2_s = safe_median(d$tau2_s_median),
    median_tau2_B = safe_median(d$tau2_B_median),
    tau2_B_coverage = mean(d$tau2_B_covered, na.rm = TRUE),
    
    median_prediction_mse = safe_median(d$prediction_mse),
    median_prediction_coverage = safe_median(d$prediction_coverage),
    
    median_delta_mse = safe_median(d$delta_prediction_mse),
    median_delta_gamma_error = safe_median(d$delta_gamma_abs_error),
    median_delta_tau2_B_error = safe_median(d$delta_tau2_B_abs_error),
    
    median_elapsed_sec = safe_median(d$elapsed_sec)
  )
}))
rownames(areal_summary) = NULL

write.csv(
  areal_summary,
  file.path(areal_output_dir, "comparison_summary.csv"),
  row.names = FALSE
)

print(areal_summary)
message("Areal-method sensitivity results saved to: ", areal_output_dir)




# Point-process likelihood -----------------------------------------------------

## Oracle preferentiality recovery --------------------------------------------

# Simulation settings
pp_replications = seq_len(100L)
pp_gamma_values = c(0, 0.5, 1, 1.5, 2)
pp_gamma_interval = c(-5, 7)

pp_output_dir = file.path(
  "results", "prerun", "point_process_likelihood"
)
dir.create(pp_output_dir, recursive = TRUE, showWarnings = FALSE)

# Record experiment settings
saveRDS(
  list(
    replications = pp_replications,
    gamma_values = pp_gamma_values,
    I = 100L,
    rho = 0.2,
    scenario = "constant",
    gamma_interval = pp_gamma_interval,
    likelihood = "conditional_raster",
    latent_field_known = TRUE
  ),
  file = file.path(pp_output_dir, "run_info.rds")
)


# Oracle likelihood ------------------------------------------------------------

fit_oracle_gamma = function(phi_grid, phi_points, gamma_true, interval) {
  
  I = length(phi_points)
  
  # Negative log-likelihood, conditional on I and the true field
  negative_loglik = function(gamma) {
    eta = gamma * phi_grid
    eta_max = max(eta)
    
    # Numerically stable log(mean(exp(eta)))
    log_integral = eta_max + log(mean(exp(eta - eta_max)))
    
    I * log_integral - gamma * sum(phi_points)
  }
  
  # Maximum likelihood estimate
  fitted = optimize(
    negative_loglik,
    interval = interval
  )
  
  gamma_mle = fitted$minimum
  
  # Profile-likelihood 95% confidence interval
  cutoff = fitted$objective + qchisq(0.95, df = 1) / 2
  
  profile_equation = function(gamma) {
    negative_loglik(gamma) - cutoff
  }
  
  gamma_q025 = if (
    gamma_mle > interval[1] &&
    profile_equation(interval[1]) > 0
  ) {
    uniroot(
      profile_equation,
      interval = c(interval[1], gamma_mle)
    )$root
  } else {
    NA_real_
  }
  
  gamma_q975 = if (
    gamma_mle < interval[2] &&
    profile_equation(interval[2]) > 0
  ) {
    uniroot(
      profile_equation,
      interval = c(gamma_mle, interval[2])
    )$root
  } else {
    NA_real_
  }
  
  # Observed information
  eta = gamma_mle * phi_grid
  weights = exp(eta - max(eta))
  weights = weights / sum(weights)
  
  weighted_mean = sum(weights * phi_grid)
  
  information = I * sum(
    weights * (phi_grid - weighted_mean)^2
  )
  
  gamma_se = if (is.finite(information) && information > 0) {
    1 / sqrt(information)
  } else {
    NA_real_
  }
  
  # Score evaluated at the true gamma
  eta_true = gamma_true * phi_grid
  probabilities = exp(eta_true - max(eta_true))
  probabilities = probabilities / sum(probabilities)
  
  score_true = sum(phi_points) -
    I * sum(probabilities * phi_grid)
  
  # Record diagnostics
  data.frame(
    gamma_mle = gamma_mle,
    gamma_se = gamma_se,
    gamma_q025 = gamma_q025,
    gamma_q975 = gamma_q975,
    gamma_covered = if (
      is.finite(gamma_q025) && is.finite(gamma_q975)
    ) {
      gamma_q025 <= gamma_true && gamma_true <= gamma_q975
    } else {
      NA
    },
    score_at_truth = score_true,
    at_boundary = (
      gamma_mle <= interval[1] + 1e-3 ||
        gamma_mle >= interval[2] - 1e-3
    )
  )
}


# Estimate gamma under the oracle likelihood ----------------------------------

for (replication in pp_replications) {
  seed = 123L + replication - 1L
  
  for (gamma_true in pp_gamma_values) {
    
    fit_path = file.path(
      pp_output_dir,
      "oracle",
      paste0("gamma_", gamma_true),
      sprintf("rep_%03d.rds", replication)
    )
    
    dir.create(
      dirname(fit_path),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    if (file.exists(fit_path)) next
    
    # Generate data using the original sampling mechanism
    simulated = generate_data(
      rho = 0.2,
      gamma = gamma_true,
      I = 100L,
      J = 25L,
      tau2_s = 0.1,
      tau2_B = 0.1,
      seed = seed,
      scenario = "constant"
    )
    
    # True field values on the raster and at selected point cells
    phi_grid = simulated$latent$phi
    phi_points = simulated$points$true
    
    # Estimate gamma with the true field known
    diagnostic = fit_oracle_gamma(
      phi_grid = phi_grid,
      phi_points = phi_points,
      gamma_true = gamma_true,
      interval = pp_gamma_interval
    )
    
    diagnostic$replication = replication
    diagnostic$seed = seed
    diagnostic$gamma_true = gamma_true
    
    saveRDS(diagnostic, file = fit_path)
  }
  
  if (replication %% 10L == 0L) {
    message("Completed oracle replication ", replication)
  }
}


# Summarize oracle results -----------------------------------------------------

oracle_files = list.files(
  file.path(pp_output_dir, "oracle"),
  pattern = "^rep_[0-9]+\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

oracle_diagnostics = do.call(rbind, lapply(
  oracle_files,
  readRDS
))
rownames(oracle_diagnostics) = NULL

write.csv(
  oracle_diagnostics,
  file.path(pp_output_dir, "oracle_diagnostics.csv"),
  row.names = FALSE
)

# Summarize by true preferentiality
groups = split(
  oracle_diagnostics,
  oracle_diagnostics$gamma_true
)

oracle_summary = do.call(rbind, lapply(groups, function(d) {
  data.frame(
    gamma_true = d$gamma_true[1],
    n_replications = nrow(d),
    median_gamma = safe_median(d$gamma_mle),
    gamma_bias = mean(d$gamma_mle - d$gamma_true),
    gamma_rmse = sqrt(mean(
      (d$gamma_mle - d$gamma_true)^2
    )),
    median_gamma_se = safe_median(d$gamma_se),
    gamma_coverage = mean(d$gamma_covered, na.rm = TRUE),
    mean_score_at_truth = mean(d$score_at_truth),
    n_boundary = sum(d$at_boundary)
  )
}))
rownames(oracle_summary) = NULL

write.csv(
  oracle_summary,
  file.path(pp_output_dir, "oracle_summary.csv"),
  row.names = FALSE
)

print(oracle_summary)
message("Oracle likelihood results saved to: ", pp_output_dir)


## Compare oracle and fitted PS Melding ---------------------------------------

# Compare the same replications used in the previous experiments
# Use the finest mesh and grid-based areal integration

fitted_reference = subset(
  areal_diagnostics,
  model == "ps_melding" &
    mesh == "finest" &
    areal_method == "grid",
  select = c(
    "replication", "gamma_true", "gamma_median"
  )
)

oracle_comparison = merge(
  oracle_diagnostics,
  fitted_reference,
  by = c("replication", "gamma_true")
)

# Compare the oracle MLE with the fitted posterior median
oracle_comparison_summary = aggregate(
  cbind(gamma_mle, gamma_median) ~ gamma_true,
  data = oracle_comparison,
  FUN = median
)

names(oracle_comparison_summary)[2:3] = c(
  "oracle_mle",
  "ps_melding_median"
)

write.csv(
  oracle_comparison,
  file.path(pp_output_dir, "oracle_vs_fitted.csv"),
  row.names = FALSE
)

print(oracle_comparison_summary)


# Point-process likelihood -----------------------------------------------------

## Oracle preferentiality recovery --------------------------------------------

# Simulation settings
pp_replications = seq_len(100L)
pp_gamma_values = c(0, 0.5, 1, 1.5, 2)
pp_gamma_interval = c(-5, 7)

pp_output_dir = file.path(
  "results", "prerun", "point_process_likelihood"
)
dir.create(pp_output_dir, recursive = TRUE, showWarnings = FALSE)

# Record experiment settings
saveRDS(
  list(
    replications = pp_replications,
    gamma_values = pp_gamma_values,
    I = 100L,
    rho = 0.2,
    scenario = "constant",
    gamma_interval = pp_gamma_interval,
    likelihood = "conditional_raster",
    latent_field_known = TRUE
  ),
  file = file.path(pp_output_dir, "run_info.rds")
)


# Oracle likelihood ------------------------------------------------------------

fit_oracle_gamma = function(phi_grid, phi_points, gamma_true, interval) {
  
  I = length(phi_points)
  
  # Negative log-likelihood, conditional on I and the true field
  negative_loglik = function(gamma) {
    eta = gamma * phi_grid
    eta_max = max(eta)
    
    # Numerically stable log(mean(exp(eta)))
    log_integral = eta_max + log(mean(exp(eta - eta_max)))
    
    I * log_integral - gamma * sum(phi_points)
  }
  
  # Maximum likelihood estimate
  fitted = optimize(
    negative_loglik,
    interval = interval
  )
  
  gamma_mle = fitted$minimum
  
  # Profile-likelihood 95% confidence interval
  cutoff = fitted$objective + qchisq(0.95, df = 1) / 2
  
  profile_equation = function(gamma) {
    negative_loglik(gamma) - cutoff
  }
  
  gamma_q025 = if (
    gamma_mle > interval[1] &&
    profile_equation(interval[1]) > 0
  ) {
    uniroot(
      profile_equation,
      interval = c(interval[1], gamma_mle)
    )$root
  } else {
    NA_real_
  }
  
  gamma_q975 = if (
    gamma_mle < interval[2] &&
    profile_equation(interval[2]) > 0
  ) {
    uniroot(
      profile_equation,
      interval = c(gamma_mle, interval[2])
    )$root
  } else {
    NA_real_
  }
  
  # Observed information
  eta = gamma_mle * phi_grid
  weights = exp(eta - max(eta))
  weights = weights / sum(weights)
  
  weighted_mean = sum(weights * phi_grid)
  
  information = I * sum(
    weights * (phi_grid - weighted_mean)^2
  )
  
  gamma_se = if (is.finite(information) && information > 0) {
    1 / sqrt(information)
  } else {
    NA_real_
  }
  
  # Score evaluated at the true gamma
  eta_true = gamma_true * phi_grid
  probabilities = exp(eta_true - max(eta_true))
  probabilities = probabilities / sum(probabilities)
  
  score_true = sum(phi_points) -
    I * sum(probabilities * phi_grid)
  
  # Record diagnostics
  data.frame(
    gamma_mle = gamma_mle,
    gamma_se = gamma_se,
    gamma_q025 = gamma_q025,
    gamma_q975 = gamma_q975,
    gamma_covered = if (
      is.finite(gamma_q025) && is.finite(gamma_q975)
    ) {
      gamma_q025 <= gamma_true && gamma_true <= gamma_q975
    } else {
      NA
    },
    score_at_truth = score_true,
    at_boundary = (
      gamma_mle <= interval[1] + 1e-3 ||
        gamma_mle >= interval[2] - 1e-3
    )
  )
}


# Estimate gamma under the oracle likelihood ----------------------------------

for (replication in pp_replications) {
  seed = 123L + replication - 1L
  
  for (gamma_true in pp_gamma_values) {
    
    fit_path = file.path(
      pp_output_dir,
      "oracle",
      paste0("gamma_", gamma_true),
      sprintf("rep_%03d.rds", replication)
    )
    
    dir.create(
      dirname(fit_path),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    if (file.exists(fit_path)) next
    
    # Generate data using the original sampling mechanism
    simulated = generate_data(
      rho = 0.2,
      gamma = gamma_true,
      I = 100L,
      J = 25L,
      tau2_s = 0.1,
      tau2_B = 0.1,
      seed = seed,
      scenario = "constant"
    )
    
    # True field values on the raster and at selected point cells
    phi_grid = simulated$latent$phi
    phi_points = simulated$points$true
    
    # Estimate gamma with the true field known
    diagnostic = fit_oracle_gamma(
      phi_grid = phi_grid,
      phi_points = phi_points,
      gamma_true = gamma_true,
      interval = pp_gamma_interval
    )
    
    diagnostic$replication = replication
    diagnostic$seed = seed
    diagnostic$gamma_true = gamma_true
    
    saveRDS(diagnostic, file = fit_path)
  }
  
  if (replication %% 10L == 0L) {
    message("Completed oracle replication ", replication)
  }
}


# Summarize oracle results -----------------------------------------------------

oracle_files = list.files(
  file.path(pp_output_dir, "oracle"),
  pattern = "^rep_[0-9]+\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

oracle_diagnostics = do.call(rbind, lapply(
  oracle_files,
  readRDS
))
rownames(oracle_diagnostics) = NULL

write.csv(
  oracle_diagnostics,
  file.path(pp_output_dir, "oracle_diagnostics.csv"),
  row.names = FALSE
)

# Summarize by true preferentiality
groups = split(
  oracle_diagnostics,
  oracle_diagnostics$gamma_true
)

oracle_summary = do.call(rbind, lapply(groups, function(d) {
  data.frame(
    gamma_true = d$gamma_true[1],
    n_replications = nrow(d),
    median_gamma = safe_median(d$gamma_mle),
    gamma_bias = mean(d$gamma_mle - d$gamma_true),
    gamma_rmse = sqrt(mean(
      (d$gamma_mle - d$gamma_true)^2
    )),
    median_gamma_se = safe_median(d$gamma_se),
    gamma_coverage = mean(d$gamma_covered, na.rm = TRUE),
    mean_score_at_truth = mean(d$score_at_truth),
    n_boundary = sum(d$at_boundary)
  )
}))
rownames(oracle_summary) = NULL

write.csv(
  oracle_summary,
  file.path(pp_output_dir, "oracle_summary.csv"),
  row.names = FALSE
)

print(oracle_summary)
message("Oracle likelihood results saved to: ", pp_output_dir)


## Compare oracle and fitted PS Melding ---------------------------------------

# Compare the same replications used in the previous experiments
# Use the finest mesh and grid-based areal integration

fitted_reference = subset(
  areal_diagnostics,
  model == "ps_melding" &
    mesh == "finest" &
    areal_method == "grid",
  select = c(
    "replication", "gamma_true", "gamma_median"
  )
)

oracle_comparison = merge(
  oracle_diagnostics,
  fitted_reference,
  by = c("replication", "gamma_true")
)

# Compare the oracle MLE with the fitted posterior median
oracle_comparison_summary = aggregate(
  cbind(gamma_mle, gamma_median) ~ gamma_true,
  data = oracle_comparison,
  FUN = median
)

names(oracle_comparison_summary)[2:3] = c(
  "oracle_mle",
  "ps_melding_median"
)

write.csv(
  oracle_comparison,
  file.path(pp_output_dir, "oracle_vs_fitted.csv"),
  row.names = FALSE
)

print(oracle_comparison_summary)


## Point-process approximation sensitivity ------------------------------------

# Simulation settings
pp_fit_replications = seq_len(10L)
pp_fit_gamma_values = c(1, 2)
pp_fit_models = "ps_melding"
pp_fit_methods = c("mesh_dual", "grid_counts")

# Keep meshes, areal method, and priors fixed
pp_fit_meshes = meshes
pp_fit_areal_method = "grid"
pp_fit_prior_gamma = c(1, 1)
pp_fit_prior_noise_B = c(1, 0.05)

pp_fit_output_dir = file.path(
  pp_output_dir, "method_sensitivity"
)
dir.create(pp_fit_output_dir, recursive = TRUE, showWarnings = FALSE)

# Record experiment settings
saveRDS(
  list(
    replications = pp_fit_replications,
    gamma_values = pp_fit_gamma_values,
    models = pp_fit_models,
    methods = pp_fit_methods,
    mesh_specs = mesh_specs,
    mesh_vertices = vapply(pp_fit_meshes, function(x) x$n, numeric(1)),
    areal_method = pp_fit_areal_method,
    prior_gamma = pp_fit_prior_gamma,
    prior_noise = NULL,
    prior_noise_B = pp_fit_prior_noise_B,
    mesh_dual_results_dir = areal_output_dir
  ),
  file = file.path(pp_fit_output_dir, "run_info.rds")
)


# Fit models -------------------------------------------------------------------

# Reuse mesh_dual fits from areal-method sensitivity.
# Only fit the grid_counts point-process approximation.

for (replication in pp_fit_replications) {
  seed = 123L + replication - 1L
  
  for (gamma_true in pp_fit_gamma_values) {
    
    # Generate the same data used in previous experiments
    simulated = generate_data(
      rho = 0.2,
      gamma = gamma_true,
      I = 100L,
      J = 25L,
      tau2_s = 0.1,
      tau2_B = 0.1,
      seed = seed,
      scenario = "constant"
    )
    
    for (mesh_name in names(pp_fit_meshes)) {
      current_mesh = pp_fit_meshes[[mesh_name]]
      
      fit_path = file.path(
        pp_fit_output_dir,
        paste0("gamma_", gamma_true),
        mesh_name,
        "grid_counts",
        sprintf("rep_%03d.rds", replication)
      )
      
      dir.create(
        dirname(fit_path),
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      if (file.exists(fit_path)) next
      
      message(
        "Replication ", replication,
        ", gamma = ", gamma_true,
        ", mesh = ", mesh_name,
        ", point-process method = grid_counts"
      )
      
      # Grid-based areal projection
      A_B = create_projections(
        simulated, current_mesh, pp_fit_areal_method
      )$areas
      
      # Fit PS Melding using grid-cell counts
      logged = fit_models_logged(
        data = simulated,
        models = pp_fit_models,
        mesh = current_mesh,
        replication = replication,
        seed = seed,
        point_process_method = "grid_counts",
        areal_method = pp_fit_areal_method,
        separate_noise = TRUE,
        prior_gamma = pp_fit_prior_gamma,
        prior_noise = NULL,
        prior_noise_B = pp_fit_prior_noise_B
      )
      
      # Summarize model fit
      diagnostic = summarise_fit(
        data = simulated,
        logged = logged,
        model = pp_fit_models,
        mesh = current_mesh,
        A_B = A_B,
        replication = replication,
        tau2_B = 0.1,
        mesh_name = mesh_name,
        areal_method = pp_fit_areal_method,
        prior_name = "pc_U1"
      )
      
      diagnostic$gamma_true = gamma_true
      diagnostic$gamma_q025 = NA_real_
      diagnostic$gamma_q975 = NA_real_
      diagnostic$gamma_covered = NA
      
      # Preferentiality credible interval
      fitted = logged$results$models[[pp_fit_models]]
      
      if (!is.null(fitted)) {
        hp = fitted$fit$summary.hyperpar
        
        if ("Beta for j" %in% rownames(hp) &&
            all(c("0.025quant", "0.975quant") %in% colnames(hp))) {
          
          diagnostic$gamma_q025 = hp["Beta for j", "0.025quant"]
          diagnostic$gamma_q975 = hp["Beta for j", "0.975quant"]
          
          diagnostic$gamma_covered =
            diagnostic$gamma_q025 <= gamma_true &&
            gamma_true <= diagnostic$gamma_q975
        }
      }
      
      saveRDS(
        list(
          data = simulated,
          results = logged$results,
          log = logged$log,
          diagnostics = diagnostic
        ),
        file = fit_path
      )
      
      rm(logged, diagnostic)
    }
  }
}


# Summarize results ------------------------------------------------------------

# Existing mesh_dual fits from areal-method sensitivity
mesh_dual_files = list.files(
  areal_output_dir,
  pattern = "^rep_[0-9]+\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

# New grid_counts fits
grid_counts_files = list.files(
  pp_fit_output_dir,
  pattern = "^rep_[0-9]+\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(mesh_dual_files) == 0L ||
    length(grid_counts_files) == 0L) {
  stop("Both mesh_dual and grid_counts results are required.")
}

# Extract existing PS Melding results
mesh_dual_diagnostics = do.call(rbind, lapply(
  mesh_dual_files,
  function(path) {
    
    d = readRDS(path)$diagnostics
    d = subset(d, model == "ps_melding")
    d$point_process_method = "mesh_dual"
    
    return(d)
  }
))

# Extract new grid_counts results
grid_counts_diagnostics = do.call(rbind, lapply(
  grid_counts_files,
  function(path) {
    
    d = readRDS(path)$diagnostics
    d$point_process_method = "grid_counts"
    
    return(d)
  }
))

# Combine both point-process methods
pp_fit_diagnostics = rbind(
  mesh_dual_diagnostics,
  grid_counts_diagnostics
)
rownames(pp_fit_diagnostics) = NULL

# Retain the settings considered in this experiment
pp_fit_diagnostics = subset(
  pp_fit_diagnostics,
  gamma_true %in% pp_fit_gamma_values &
    mesh %in% names(pp_fit_meshes) &
    areal_method == pp_fit_areal_method &
    model == pp_fit_models &
    point_process_method %in% pp_fit_methods
)

# Reference: mesh_dual on the same mesh and dataset
mesh_dual_reference = subset(
  pp_fit_diagnostics,
  point_process_method == "mesh_dual",
  select = c(
    "replication", "gamma_true", "mesh", "model",
    "prediction_mse", "gamma_median", "tau2_B_median"
  )
)

names(mesh_dual_reference)[5:7] = c(
  "mesh_dual_prediction_mse",
  "mesh_dual_gamma_median",
  "mesh_dual_tau2_B_median"
)

pp_fit_diagnostics = merge(
  pp_fit_diagnostics,
  mesh_dual_reference,
  by = c("replication", "gamma_true", "mesh", "model"),
  all.x = TRUE
)

# Paired differences relative to mesh_dual
# Negative values indicate improvement
pp_fit_diagnostics$delta_prediction_mse =
  pp_fit_diagnostics$prediction_mse -
  pp_fit_diagnostics$mesh_dual_prediction_mse

pp_fit_diagnostics$delta_gamma_abs_error =
  abs(pp_fit_diagnostics$gamma_median - pp_fit_diagnostics$gamma_true) -
  abs(pp_fit_diagnostics$mesh_dual_gamma_median - pp_fit_diagnostics$gamma_true)

pp_fit_diagnostics$delta_tau2_B_abs_error =
  abs(pp_fit_diagnostics$tau2_B_median - pp_fit_diagnostics$tau2_B_true) -
  abs(pp_fit_diagnostics$mesh_dual_tau2_B_median - pp_fit_diagnostics$tau2_B_true)

write.csv(
  pp_fit_diagnostics,
  file.path(pp_fit_output_dir, "fit_diagnostics.csv"),
  row.names = FALSE
)

# Summarize by true gamma, mesh, and point-process method
groups = split(
  pp_fit_diagnostics,
  interaction(
    pp_fit_diagnostics$gamma_true,
    pp_fit_diagnostics$mesh,
    pp_fit_diagnostics$point_process_method,
    drop = TRUE
  )
)

pp_fit_summary = do.call(rbind, lapply(groups, function(d) {
  data.frame(
    gamma_true = d$gamma_true[1],
    mesh = d$mesh[1],
    mesh_vertices = d$mesh_vertices[1],
    point_process_method = d$point_process_method[1],
    model = d$model[1],
    
    n_attempted = nrow(d),
    n_included = sum(d$included),
    
    median_gamma = safe_median(d$gamma_median),
    gamma_rmse = if (all(is.na(d$gamma_median))) NA_real_ else
      sqrt(mean((d$gamma_median - d$gamma_true)^2, na.rm = TRUE)),
    gamma_coverage = if (all(is.na(d$gamma_covered))) NA_real_ else
      mean(d$gamma_covered, na.rm = TRUE),
    
    median_tau2_s = safe_median(d$tau2_s_median),
    median_tau2_B = safe_median(d$tau2_B_median),
    tau2_B_coverage = mean(d$tau2_B_covered, na.rm = TRUE),
    
    median_prediction_mse = safe_median(d$prediction_mse),
    median_prediction_coverage = safe_median(d$prediction_coverage),
    
    median_delta_mse = safe_median(d$delta_prediction_mse),
    median_delta_gamma_error = safe_median(d$delta_gamma_abs_error),
    median_delta_tau2_B_error = safe_median(d$delta_tau2_B_abs_error),
    
    median_elapsed_sec = safe_median(d$elapsed_sec)
  )
}))
rownames(pp_fit_summary) = NULL

write.csv(
  pp_fit_summary,
  file.path(pp_fit_output_dir, "comparison_summary.csv"),
  row.names = FALSE
)

print(pp_fit_summary)
message("Point-process sensitivity results saved to: ", pp_fit_output_dir)