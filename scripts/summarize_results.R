source("functions/create_summary.R")

# Summarize the completed first-laptop run
# Change to seq_len(100L) after combining both laptops
replications <- seq_len(50L)
summarize_reproduction <- FALSE


# Read simulation settings ----------------------------------------------------

default_settings <- list(
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
)

read_run_info <- function(output_dir, defaults) {
  file <- file.path(output_dir, "run_info.rds")
  
  if (file.exists(file)) {
    return(readRDS(file))
  }
  
  warning("run_info.rds is missing from ", output_dir,
          "; using the original simulation settings. Verify these values.")
  
  return(defaults)
}

make_parameter_truth <- function(settings,
                                 gamma = settings$gamma,
                                 tau2_B = settings$tau2_B) {
  sigma2 <- if (is.null(settings$sigma2)) 1 else settings$sigma2
  
  return(c(
    gamma = gamma,
    rho = settings$rho,
    sigma2 = sigma2,
    tau2_s = settings$tau2_s,
    tau2_B = tau2_B
  ))
}


# Load and summarize saved replications ---------------------------------------

summarize_saved_run <- function(
    output_dir,
    gamma_true = 1,
    tau2_s = 0.1,
    parameter_truth = c(gamma = 1, rho = 0.2, sigma2 = 1,
                        tau2_s = 0.1, tau2_B = 0.1),
    replications = NULL,
    n_replications = NULL,
    models = c("ps_melding", "melding", "ps_geo")
) {
  fits_dir <- file.path(output_dir, "fits")
  log_dir <- file.path(output_dir, "log")
  
  fit_files <- sort(list.files(
    fits_dir,
    pattern = "^rep_[0-9]+\\.rds$",
    full.names = TRUE
  ))
  
  log_files <- sort(list.files(
    log_dir,
    pattern = "^rep_[0-9]+\\.rds$",
    full.names = TRUE
  ))
  
  # Restrict summaries to the requested replication IDs
  if (is.null(replications) && !is.null(n_replications)) {
    replications <- seq_len(n_replications)
  }
  
  if (!is.null(replications)) {
    expected <- sprintf("rep_%03d.rds", replications)
    fit_files <- file.path(fits_dir, intersect(expected, basename(fit_files)))
    log_files <- file.path(log_dir, intersect(expected, basename(log_files)))
    
    if (length(fit_files) < length(replications)) {
      message("Found ", length(fit_files), " of ", length(replications),
              " saved fits in: ", output_dir)
    }
  }
  
  if (length(fit_files) == 0L && length(log_files) == 0L) {
    message("No saved replications found in: ", output_dir)
    return(NULL)
  }
  
  # Read the separate fitting logs, including unsuccessful fits
  logs <- lapply(log_files, function(file) {
    tryCatch(
      readRDS(file),
      error = function(e) {
        warning("Could not load log ", file, ": ", conditionMessage(e))
        NULL
      }
    )
  })
  
  valid_logs <- !vapply(logs, is.null, logical(1))
  log_files <- log_files[valid_logs]
  logs <- logs[valid_logs]
  
  # Older fit files may have only an embedded log
  missing_logs <- setdiff(basename(fit_files), basename(log_files))
  
  if (length(missing_logs) > 0L) {
    message("Reading embedded logs for ", length(missing_logs),
            " replication(s) in: ", output_dir)
    
    logs <- c(logs, lapply(missing_logs, function(filename) {
      file <- file.path(fits_dir, filename)
      
      tryCatch(
        readRDS(file)$log,
        error = function(e) {
          warning("Could not recover embedded log from ", file, ": ",
                  conditionMessage(e))
          NULL
        }
      )
    }))
  }
  
  fit_log <- do.call(rbind, logs)
  
  if (!is.null(fit_log)) {
    rownames(fit_log) <- NULL
  }
  
  fit_status <- if (!is.null(fit_log)) {
    summarize_fit_log(fit_log, models = models)
  } else {
    NULL
  }
  
  # Evaluate each saved fit independently to avoid retaining all INLA objects
  evaluated <- lapply(fit_files, function(file) {
    saved <- tryCatch(
      readRDS(file),
      error = function(e) {
        warning("Could not load ", file, ": ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(saved) || length(saved$results$models) == 0L) {
      return(list(metrics = NULL, parameter_metrics = NULL))
    }
    
    replication <- saved$log$replication[1]
    
    metrics <- tryCatch(
      evaluate_replication(
        results = saved$results,
        truth = saved$data$prediction$true,
        tau2_s = tau2_s,
        gamma_true = gamma_true
      ),
      error = function(e) {
        warning("Prediction evaluation skipped for ", file, ": ",
                conditionMessage(e))
        NULL
      }
    )
    
    if (!is.null(metrics)) {
      metrics$replication <- replication
    }
    
    parameter_metrics <- tryCatch(
      evaluate_parameter_replication(
        results = saved$results,
        truth = parameter_truth,
        separate_noise = TRUE
      ),
      error = function(e) {
        warning("Parameter evaluation skipped for ", file, ": ",
                conditionMessage(e))
        NULL
      }
    )
    
    if (!is.null(parameter_metrics)) {
      parameter_metrics$replication <- replication
    }
    
    return(list(
      metrics = metrics,
      parameter_metrics = parameter_metrics
    ))
  })
  
  # Combine replication-level results
  replication_metrics <- do.call(rbind, lapply(evaluated, function(x) {
    x$metrics
  }))
  
  parameter_metrics <- do.call(rbind, lapply(evaluated, function(x) {
    x$parameter_metrics
  }))
  
  summary_results <- NULL
  parameter_summary <- NULL
  
  if (!is.null(replication_metrics)) {
    rownames(replication_metrics) <- NULL
    summary_results <- summarize_replications(
      metrics = replication_metrics,
      gamma_true = gamma_true
    )
  }
  
  if (!is.null(parameter_metrics)) {
    rownames(parameter_metrics) <- NULL
    parameter_summary <- summarize_parameter_recovery(parameter_metrics)
  }
  
  return(list(
    replication_metrics = replication_metrics,
    parameter_metrics = parameter_metrics,
    fit_log = fit_log,
    fit_status = fit_status,
    summary_results = summary_results,
    parameter_summary = parameter_summary
  ))
}


# Combine model summaries across extension settings ---------------------------

combine_summaries <- function(
    results, values, column, metric, component = "summary_results"
) {
  comparison <- lapply(seq_along(results), function(k) {
    result <- results[[k]]
    
    if (is.null(result) || is.null(result[[component]])) {
      return(NULL)
    }
    
    summary <- result[[component]][[metric]]
    
    if (is.null(summary)) {
      return(NULL)
    }
    
    comparison <- data.frame(setting = values[k], summary)
    names(comparison)[1] <- column
    
    return(comparison)
  })
  
  comparison <- do.call(rbind, comparison)
  
  if (!is.null(comparison)) {
    rownames(comparison) <- NULL
  }
  
  return(comparison)
}


# Summarize baseline reproduction ---------------------------------------------

reproduction <- NULL

if (summarize_reproduction) {
  reproduction_dir <- file.path("results", "reproduction")
  reproduction_info <- read_run_info(
    reproduction_dir,
    defaults = list(generating = default_settings)
  )
  
  reproduction_settings <- reproduction_info$generating
  
  reproduction <- summarize_saved_run(
    output_dir = reproduction_dir,
    gamma_true = reproduction_settings$gamma,
    tau2_s = reproduction_settings$tau2_s,
    parameter_truth = make_parameter_truth(reproduction_settings),
    n_replications = reproduction_info$n_replications
  )
  
  if (!is.null(reproduction)) {
    print(reproduction$fit_status)
    print(reproduction$summary_results$prediction_mean)
    print(reproduction$summary_results$prediction_mcse)
    print(reproduction$summary_results$gamma)
    print(reproduction$summary_results$gamma_mcse)
    
    if (!is.null(reproduction$parameter_summary)) {
      print(reproduction$parameter_summary$parameter_mean)
      print(reproduction$parameter_summary$parameter_mcse)
    }
  }
}


# Summarize preferential sampling mechanisms ---------------------------------

extensions_dir <- file.path("results", "extensions")
extension_info <- read_run_info(
  extensions_dir,
  defaults = list(
    generating = default_settings,
    scenarios = c("spatially_varying", "independent_field"),
    measurement_noise = list(tau2_B_values = c(0.025, 0.05, 0.2, 0.4)),
    gamma_strength = list(gamma_values = c(0, 0.5, 1.5, 2))
  )
)

extension_settings <- extension_info$generating
scenarios <- extension_info$scenarios

extension_results <- lapply(scenarios, function(scenario) {
  # A single true gamma does not exist under spatially varying gamma(s)
  gamma_true <- if (scenario == "spatially_varying") {
    NA_real_
  } else {
    extension_settings$gamma
  }
  
  parameter_truth <- make_parameter_truth(extension_settings, gamma = gamma_true)
  
  if (scenario == "spatially_varying") {
    parameter_truth <- c(parameter_truth, delta = extension_settings$delta)
  }
  
  if (scenario == "independent_field") {
    parameter_truth <- c(
      parameter_truth,
      beta_psi = extension_settings$beta_psi,
      rho_psi = extension_settings$rho_psi
    )
  }
  
  summarize_saved_run(
    output_dir = file.path(extensions_dir, scenario),
    gamma_true = gamma_true,
    tau2_s = extension_settings$tau2_s,
    parameter_truth = parameter_truth,
    replications = replications
  )
})

names(extension_results) <- scenarios

prediction_mean <- combine_summaries(
  extension_results, scenarios, "scenario", "prediction_mean"
)

prediction_n <- combine_summaries(
  extension_results, scenarios, "scenario", "prediction_n"
)

prediction_mcse <- combine_summaries(
  extension_results, scenarios, "scenario", "prediction_mcse"
)

scenario_parameter_mean <- combine_summaries(
  extension_results, scenarios, "scenario", "parameter_mean",
  component = "parameter_summary"
)

scenario_parameter_mcse <- combine_summaries(
  extension_results, scenarios, "scenario", "parameter_mcse",
  component = "parameter_summary"
)

print(prediction_n)
print(prediction_mean)
print(prediction_mcse)
print(scenario_parameter_mean)
print(scenario_parameter_mcse)


# Summarize varying areal measurement noise ----------------------------------

tau2_B_values <- extension_info$measurement_noise$tau2_B_values

noise_results <- lapply(tau2_B_values, function(tau2_B) {
  summarize_saved_run(
    output_dir = file.path(
      extensions_dir, "measurement_noise", paste0("tau2_B_", tau2_B)
    ),
    gamma_true = extension_settings$gamma,
    tau2_s = extension_settings$tau2_s,
    parameter_truth = make_parameter_truth(extension_settings, tau2_B = tau2_B),
    replications = replications
  )
})

names(noise_results) <- paste0("tau2_B_", tau2_B_values)

noise_prediction_mean <- combine_summaries(
  noise_results, tau2_B_values, "tau2_B", "prediction_mean"
)

noise_prediction_n <- combine_summaries(
  noise_results, tau2_B_values, "tau2_B", "prediction_n"
)

noise_prediction_mcse <- combine_summaries(
  noise_results, tau2_B_values, "tau2_B", "prediction_mcse"
)

noise_parameter_mean <- combine_summaries(
  noise_results, tau2_B_values, "tau2_B", "parameter_mean",
  component = "parameter_summary"
)

noise_parameter_mcse <- combine_summaries(
  noise_results, tau2_B_values, "tau2_B", "parameter_mcse",
  component = "parameter_summary"
)

print(noise_prediction_n)
print(noise_prediction_mean)
print(noise_prediction_mcse)
if (!is.null(noise_parameter_mean)) {
  print(subset(noise_parameter_mean, parameter %in% c("tau2_s", "tau2_B")))
  print(subset(noise_parameter_mcse, parameter %in% c("tau2_s", "tau2_B")))
}


# Summarize varying preferential sampling strength ----------------------------

gamma_values <- extension_info$gamma_strength$gamma_values

gamma_results <- lapply(gamma_values, function(gamma_value) {
  summarize_saved_run(
    output_dir = file.path(
      extensions_dir, "gamma_strength", paste0("gamma_", gamma_value)
    ),
    gamma_true = gamma_value,
    tau2_s = extension_settings$tau2_s,
    parameter_truth = make_parameter_truth(
      extension_settings, gamma = gamma_value
    ),
    replications = replications
  )
})

names(gamma_results) <- paste0("gamma_", gamma_values)

gamma_prediction_mean <- combine_summaries(
  gamma_results, gamma_values, "gamma_true", "prediction_mean"
)

gamma_prediction_n <- combine_summaries(
  gamma_results, gamma_values, "gamma_true", "prediction_n"
)

gamma_prediction_mcse <- combine_summaries(
  gamma_results, gamma_values, "gamma_true", "prediction_mcse"
)

gamma_recovery <- combine_summaries(
  gamma_results, gamma_values, "gamma_true", "gamma"
)

gamma_recovery_mcse <- combine_summaries(
  gamma_results, gamma_values, "gamma_true", "gamma_mcse"
)

gamma_parameter_mean <- combine_summaries(
  gamma_results, gamma_values, "gamma_true", "parameter_mean",
  component = "parameter_summary"
)

gamma_parameter_mcse <- combine_summaries(
  gamma_results, gamma_values, "gamma_true", "parameter_mcse",
  component = "parameter_summary"
)

print(gamma_prediction_n)
print(gamma_prediction_mean)
print(gamma_prediction_mcse)
print(gamma_recovery)
print(gamma_recovery_mcse)
if (!is.null(gamma_parameter_mean)) {
  print(subset(gamma_parameter_mean, parameter == "gamma"))
  print(subset(gamma_parameter_mcse, parameter == "gamma"))
}


# Inspect parameter recovery and fitting diagnostics --------------------------

# Fitting outcomes for each extension setting
print(lapply(extension_results, function(x) x$fit_status))
print(lapply(noise_results, function(x) x$fit_status))
print(lapply(gamma_results, function(x) x$fit_status))

# Examine warnings, errors, and numerically invalid fits
problem_fits <- function(x) {
  if (is.null(x) || is.null(x$fit_log)) {
    return(NULL)
  }
  
  return(subset(x$fit_log, status != "ok" | n_warnings > 0L))
}

print(lapply(extension_results, function(x) utils::head(problem_fits(x), 10L)))
print(lapply(noise_results, function(x) utils::head(problem_fits(x), 10L)))
print(lapply(gamma_results, function(x) utils::head(problem_fits(x), 10L)))

# Examples:
# reproduction$parameter_summary$parameter_mean
# reproduction$summary_results$prediction_reproduction
# extension_results$spatially_varying$replication_metrics
# noise_results$tau2_B_0.025$fit_status
# gamma_results$gamma_2$parameter_summary$parameter_mcse

