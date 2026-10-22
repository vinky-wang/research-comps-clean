# Summarize retained prior-sensitivity checkpoints without refitting
# Usage: Rscript report_prior_sensitivity.R [study_root] [report_directory]
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else file.path(getwd(), "prior_sensitivity")
output <- if (length(args) >= 2L) args[2L] else file.path(root, "reports")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
plan <- read.csv(file.path(root, "plan.csv"), stringsAsFactors = FALSE)
settings <- readRDS(file.path(root, "config.rds"))
screens <- settings$screens
records <- setNames(lapply(plan$job_id, function(id) {
  path <- file.path(root, "checkpoints", paste0(id, ".rds"))
  if (file.exists(path)) readRDS(path) else NULL
}), plan$job_id)
fields <- c("posterior_mean", "median", "lower", "upper")
all_fields <- c(fields, "posterior_sd")
write_table <- function(x, file) write.csv(x, file.path(output, file), row.names = FALSE)
bind_rows <- function(x) {
  x <- Filter(function(a) !is.null(a) && nrow(a) > 0L, x)
  if (length(x)) do.call(rbind, x) else data.frame()
}
finite_range <- function(x) if (length(x) && all(is.finite(x))) diff(range(x)) else NA_real_
finite_max <- function(x) if (length(x) && all(is.finite(x))) max(x) else NA_real_
within_limit <- function(x, limit) length(x) == 1L && is.finite(x) && x <= limit
diagnostic <- function(r, name) {
  value <- r$diagnostics[[name]]
  if (is.null(value) || length(value) != 1L) NA_real_ else as.numeric(value)
}
boundary_probability <- function(r) {
  if (!is.data.frame(r$boundary)) return(NA_real_)
  b <- r$boundary[r$boundary$component == "area" & abs(r$boundary$threshold - .001) < 1e-12, ]
  if (nrow(b) != 1L) NA_real_ else b$probability_below
}
parameters_for <- function(model) c("mu", "tau2_s", "tau2_B", "rho", "sigma2", "theta",
                                   if (model == "PSmelding") "gamma")
valid_parameters <- function(r, wanted) {
  p <- r$parameters
  if (!is.data.frame(p) || !all(c("parameter", all_fields) %in% names(p))) return(FALSE)
  p <- p[p$parameter %in% wanted, , drop = FALSE]
  if (nrow(p) != length(wanted) || anyDuplicated(p$parameter) ||
      !all(wanted %in% p$parameter) || !all(is.finite(as.matrix(p[all_fields])))) return(FALSE)
  if ("summary_valid" %in% names(p) && !all(p$summary_valid %in% TRUE)) return(FALSE)
  if ("coherence_warning" %in% names(p) && !all(p$coherence_warning %in% FALSE)) return(FALSE)
  positive <- p$parameter %in% c("tau2_s", "tau2_B", "rho", "sigma2", "theta")
  all(p$posterior_sd > 0, p$lower <= p$median, p$median <= p$upper,
      p$lower[positive] > 0, p$posterior_mean[positive] > 0)
}
parameter_quality_flags <- function(r, wanted) {
  p <- r$parameters
  if (!is.data.frame(p) || !"parameter" %in% names(p)) return(data.frame(
    summary_valid_flag_available = FALSE, coherence_flag_available = FALSE,
    declared_invalid_parameters = NA_integer_, parameter_coherence_warnings = NA_integer_,
    unknown_parameter_quality_flags = NA_integer_))
  p <- p[p$parameter %in% wanted, , drop = FALSE]
  has_valid <- "summary_valid" %in% names(p)
  has_coherence <- "coherence_warning" %in% names(p)
  data.frame(summary_valid_flag_available = has_valid, coherence_flag_available = has_coherence,
    declared_invalid_parameters = if (has_valid) sum(p$summary_valid %in% FALSE) else 0L,
    parameter_coherence_warnings = if (has_coherence) sum(p$coherence_warning %in% TRUE) else 0L,
    unknown_parameter_quality_flags = if (has_valid || has_coherence)
      sum(is.na(p$summary_valid)) + sum(is.na(p$coherence_warning)) else 0L)
}

# Retain every planned attempt, including incomplete or failed attempts
job_rows <- lapply(seq_len(nrow(plan)), function(i) {
  r <- records[[plan$job_id[i]]]
  values <- if (isTRUE(r$score_success) && is.data.frame(r$metrics)) r$metrics[1L, ] else
    data.frame(MSE = NA_real_, CRPS_Z_joint = NA_real_, CRPS_MCSE_approx = NA_real_, CRPS_Z_gaussian = NA_real_)
  data.frame(plan[i, ], recorded = !is.null(r), fit_success = isTRUE(r$fit_success),
    score_success = isTRUE(r$score_success),
    parameters_valid = valid_parameters(r, parameters_for(plan$model[i])),
    parameter_quality_flags(r, parameters_for(plan$model[i])),
    mode_status = diagnostic(r, "mode_status"),
    negative_eigenvalues = diagnostic(r, "max_negative_eigenvalues"),
    log_posterior_mode = diagnostic(r, "log_posterior_mode"),
    integration_excess = diagnostic(r, "max_integration_log_relative_density"),
    max_stdev_correction = diagnostic(r, "max_stdev_correction"),
    areal_probability_below_0p001 = boundary_probability(r), values,
    fit_seconds = if (is.null(r$fit_seconds)) NA_real_ else r$fit_seconds,
    warnings = paste(r$warnings, collapse = " | "),
    error = if (is.null(r$error)) "" else r$error, stringsAsFactors = FALSE)
})
jobs <- bind_rows(job_rows)
write_table(jobs, "all_attempts.csv")
group_plan <- unique(plan[c("dataset", "model", "prior_id")])
group_plan <- group_plan[order(match(group_plan$dataset, c("baseline", "clustering_delta2")),
  match(group_plan$model, c("Melding", "PSmelding")), group_plan$prior_id), ]
parameter_spreads <- list()
pairwise_parameters <- list()
pairwise_predictions <- list()
group_rows <- vector("list", nrow(group_plan))
for (g in seq_len(nrow(group_plan))) {
  key <- group_plan[g, ]
  selected <- jobs$dataset == key$dataset & jobs$model == key$model & jobs$prior_id == key$prior_id
  j <- jobs[selected, ]
  r <- records[j$job_id]
  wanted <- parameters_for(key$model)
  p <- bind_rows(lapply(seq_along(r), function(i) {
    if (!is.data.frame(r[[i]]$parameters)) return(NULL)
    z <- r[[i]]$parameters
    z <- z[z$parameter %in% wanted, , drop = FALSE]
    data.frame(start_area = j$start_area[i], z)
  }))
  local_pairs <- list()
  for (parameter in wanted) {
    z <- if (nrow(p)) p[p$parameter == parameter, , drop = FALSE] else data.frame()
    complete <- nrow(z) == 3L && !anyDuplicated(z$start_area) &&
      all(is.finite(as.matrix(z[all_fields])))
    row <- data.frame(key, parameter = parameter, n_summaries = nrow(z))
    for (field in all_fields) row[[paste0(field, "_absolute_spread")]] <-
      if (complete) finite_range(z[[field]]) else NA_real_
    parameter_spreads[[length(parameter_spreads) + 1L]] <- row
    if (nrow(z) >= 2L) {
      combinations <- combn(seq_len(nrow(z)), 2L)
      for (k in seq_len(ncol(combinations))) {
        a <- z[combinations[1L, k], ]; b <- z[combinations[2L, k], ]
        pooled_sd <- sqrt((a$posterior_sd^2 + b$posterior_sd^2) / 2)
        for (field in fields) {
          shift <- abs(a[[field]] - b[[field]])
          scaled <- if (is.finite(pooled_sd) && pooled_sd > 0) shift / pooled_sd else NA_real_
          local_pairs[[length(local_pairs) + 1L]] <- data.frame(key, parameter = parameter,
            statistic = field, start_1 = a$start_area, start_2 = b$start_area,
            value_1 = a[[field]], value_2 = b[[field]], absolute_shift = shift,
            pooled_posterior_sd = pooled_sd, shift_pooled_sd = scaled,
            screen_pass = within_limit(scaled, screens$parameter_shift_pooled_sd))
        }
      }
    }
  }
  pp <- bind_rows(local_pairs)
  pairwise_parameters[[g]] <- pp
  pred_pairs <- list()
  combinations <- combn(seq_along(r), 2L)
  for (k in seq_len(ncol(combinations))) {
    a <- combinations[1L, k]; b <- combinations[2L, k]
    va <- r[[a]]$prediction_mean; vb <- r[[b]]$prediction_mean
    rms <- if (length(va) == 2500L && length(vb) == 2500L && all(is.finite(c(va, vb))))
      sqrt(mean((va - vb)^2)) else NA_real_
    pred_pairs[[k]] <- data.frame(key, start_1 = j$start_area[a], start_2 = j$start_area[b],
      prediction_mean_RMS = rms, screen_pass = within_limit(rms, screens$prediction_RMS))
  }
  pr <- bind_rows(pred_pairs)
  pairwise_predictions[[g]] <- pr
  mode_spread <- finite_range(j$log_posterior_mode)
  parameter_shift <- if (nrow(pp) == length(wanted) * 3L * length(fields))
    finite_max(pp$shift_pooled_sd) else NA_real_
  statistic_shift <- function(statistics) {
    rows <- pp[pp$statistic %in% statistics, , drop = FALSE]
    if (nrow(rows) == length(wanted) * 3L * length(statistics))
      finite_max(rows$shift_pooled_sd) else NA_real_
  }
  mean_shift <- statistic_shift("posterior_mean")
  median_shift <- statistic_shift("median")
  endpoint_shift <- statistic_shift(c("lower", "upper"))
  prediction_rms <- finite_max(pr$prediction_mean_RMS)
  boundary_spread <- finite_range(j$areal_probability_below_0p001)
  component <- c(all_attempts_recorded = nrow(j) == 3L && all(j$recorded),
    all_fitted_and_scored = all(j$fit_success & j$score_success),
    required_parameters_valid = all(j$parameters_valid),
    mode_status_pass = all(is.finite(j$mode_status) & j$mode_status == 0),
    Hessian_pass = all(is.finite(j$negative_eigenvalues) & j$negative_eigenvalues == 0),
    integration_pass = all(is.finite(j$integration_excess) & j$integration_excess <= screens$positive_integration_log_density),
    log_mode_spread_pass = within_limit(mode_spread, screens$log_mode_spread),
    parameter_stability_pass = within_limit(parameter_shift, screens$parameter_shift_pooled_sd),
    prediction_stability_pass = within_limit(prediction_rms, screens$prediction_RMS),
    boundary_stability_pass = within_limit(boundary_spread, screens$boundary_probability_spread))
  worst <- if (nrow(pp) && any(is.finite(pp$shift_pooled_sd))) pp[which.max(pp$shift_pooled_sd), ] else NULL
  group_rows[[g]] <- data.frame(key, planned = nrow(j), recorded = sum(j$recorded),
    fitted = sum(j$fit_success), scored = sum(j$score_success),
    declared_invalid_parameter_summaries = sum(j$declared_invalid_parameters),
    parameter_coherence_warnings = sum(j$parameter_coherence_warnings),
    unknown_parameter_quality_flags = sum(j$unknown_parameter_quality_flags),
    all_components_pass = all(component), as.list(component),
    log_mode_spread = mode_spread, maximum_parameter_shift_pooled_sd = parameter_shift,
    maximum_parameter_mean_shift_pooled_sd = mean_shift,
    maximum_parameter_median_shift_pooled_sd = median_shift,
    maximum_interval_endpoint_shift_pooled_sd = endpoint_shift,
    mean_only_stability_pass = within_limit(mean_shift, screens$parameter_shift_pooled_sd),
    worst_parameter = if (is.null(worst)) NA_character_ else worst$parameter,
    worst_statistic = if (is.null(worst)) NA_character_ else worst$statistic,
    max_prediction_mean_RMS = prediction_rms, areal_boundary_probability_spread = boundary_spread,
    MSE_start_spread = finite_range(j$MSE), CRPS_start_spread = finite_range(j$CRPS_Z_joint),
    max_CRPS_MCSE_approx = finite_max(j$CRPS_MCSE_approx),
    failed_components = paste(names(component)[!component], collapse = " | "), stringsAsFactors = FALSE)
}
groups <- bind_rows(group_rows)
write_table(groups, "stability_by_prior.csv")
write_table(bind_rows(parameter_spreads), "parameter_absolute_spreads.csv")
write_table(bind_rows(pairwise_parameters), "parameter_start_comparisons.csv")
write_table(bind_rows(pairwise_predictions), "prediction_start_comparisons.csv")

