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

# Compare priors at matched initialization; starts are not simulation replicates
contrasts <- lapply(which(jobs$prior_id != "A"), function(i) {
  b <- jobs[i, ]
  a <- jobs[jobs$dataset == b$dataset & jobs$model == b$model & jobs$start_area == b$start_area & jobs$prior_id == "A", ]
  stopifnot(nrow(a) == 1L)
  ra <- records[[a$job_id]]; rb <- records[[b$job_id]]
  valid <- isTRUE(a$score_success) && isTRUE(b$score_success)
  mse_difference <- if (valid) b$MSE - a$MSE else NA_real_
  crps_difference <- if (valid) b$CRPS_Z_joint - a$CRPS_Z_joint else NA_real_
  relative <- function(change, reference) if (is.finite(reference) && reference > 0) change / reference else NA_real_
  mse_relative <- relative(mse_difference, a$MSE)
  crps_relative <- relative(crps_difference, a$CRPS_Z_joint)
  ba <- ra$batch_crps; bb <- rb$batch_crps
  mcse <- if (length(ba) == 10L && length(bb) == 10L && all(is.finite(c(ba, bb)))) sd(bb - ba) / sqrt(10) else NA_real_
  a_stable <- groups$all_components_pass[groups$dataset == b$dataset & groups$model == b$model & groups$prior_id == "A"]
  data.frame(dataset = b$dataset, model = b$model, prior_id = b$prior_id, start_area = b$start_area,
    reference_job = a$job_id, alternative_job = b$job_id, both_scored = valid,
    reference_A_stability_pass = a_stable, MSE_A = a$MSE, MSE_alternative = b$MSE,
    MSE_difference = mse_difference, MSE_relative_change = mse_relative,
    MSE_relative_material = if (is.finite(mse_relative)) abs(mse_relative) > screens$predictive_relative_change else NA,
    CRPS_A = a$CRPS_Z_joint, CRPS_alternative = b$CRPS_Z_joint,
    CRPS_difference = crps_difference, CRPS_relative_change = crps_relative,
    CRPS_relative_material = if (is.finite(crps_relative)) abs(crps_relative) > screens$predictive_relative_change else NA,
    CRPS_absolute_material = if (is.finite(crps_difference)) abs(crps_difference) > screens$crps_absolute_change else NA,
    CRPS_difference_MCSE_approx = mcse,
    MCSE_method = "SD of 10 matched batch CRPS differences / sqrt(10); approximate block proxy",
    stringsAsFactors = FALSE)
})
contrasts <- bind_rows(contrasts)
write_table(contrasts, "matched_prior_predictive_contrasts.csv")
parameter_contrasts <- bind_rows(lapply(seq_len(nrow(contrasts)), function(i) {
  c <- contrasts[i, ]
  a <- records[[c$reference_job]]$parameters
  b <- records[[c$alternative_job]]$parameters
  if (!is.data.frame(a) || !is.data.frame(b)) return(NULL)
  wanted <- parameters_for(c$model)
  bind_rows(lapply(wanted, function(parameter) {
    pa <- a[a$parameter == parameter, , drop = FALSE]
    pb <- b[b$parameter == parameter, , drop = FALSE]
    if (nrow(pa) != 1L || nrow(pb) != 1L) return(NULL)
    row <- data.frame(c[c("dataset", "model", "prior_id", "start_area", "reference_job", "alternative_job")],
      parameter = parameter, stringsAsFactors = FALSE)
    for (field in all_fields) {
      row[[paste0(field, "_A")]] <- pa[[field]]
      row[[paste0(field, "_alternative")]] <- pb[[field]]
      row[[paste0(field, "_difference")]] <- pb[[field]] - pa[[field]]
    }
    row
  }))
}))
write_table(parameter_contrasts, "matched_prior_parameter_contrasts.csv")
contrast_groups <- unique(contrasts[c("dataset", "model", "prior_id")])
contrast_ranges <- bind_rows(lapply(seq_len(nrow(contrast_groups)), function(i) {
  key <- contrast_groups[i, ]
  z <- contrasts[contrasts$dataset == key$dataset & contrasts$model == key$model & contrasts$prior_id == key$prior_id, ]
  row <- data.frame(key, n_scored_pairs = sum(z$both_scored))
  for (field in c("MSE_difference", "MSE_relative_change", "CRPS_difference", "CRPS_relative_change", "CRPS_difference_MCSE_approx")) {
    valid <- all(is.finite(z[[field]]))
    row[[paste0(field, "_min")]] <- if (valid) min(z[[field]]) else NA_real_
    row[[paste0(field, "_max")]] <- if (valid) max(z[[field]]) else NA_real_
  }
  row
}))
write_table(contrast_ranges, "predictive_contrast_ranges.csv")

# Draw point summaries for every initialization rather than selecting one fit
colors <- c(A = "#475569", B = "#0072B2", C = "#D55E00")
panel_plan <- unique(plan[c("dataset", "model")])
panel_plan <- panel_plan[order(match(panel_plan$dataset, c("baseline", "clustering_delta2")),
  match(panel_plan$model, c("Melding", "PSmelding"))), ]
label_panel <- function(key) paste(if (key$dataset == "baseline") "Baseline" else "Independent clustering: delta = 2", key$model, sep = " | ")
job_x <- function(j) match(j$prior_id, c("A", "B", "C")) + c(-.18, 0, .18)[match(j$start_area, c(.001, .1, .4))]
plot_title <- if (all(jobs$recorded)) "All 36 attempts recorded" else paste(sum(jobs$recorded), "of 36 attempts recorded: preliminary")
parameter_rows <- bind_rows(lapply(seq_len(nrow(plan)), function(i) {
  r <- records[[plan$job_id[i]]]
  if (is.data.frame(r$parameters)) data.frame(plan[i, ], r$parameters, row.names = NULL) else NULL
}))
if (nrow(parameter_rows)) {
  write_table(parameter_rows, "all_parameter_summaries.csv")
  png(file.path(output, "areal_variance_by_prior.png"), width = 1800, height = 1500, res = 160)
  par(mfrow = c(2, 2), mar = c(4, 5, 3, 1), oma = c(4, 0, 4, 0), las = 1)
  for (i in seq_len(nrow(panel_plan))) {
    key <- panel_plan[i, ]
    z <- parameter_rows[parameter_rows$dataset == key$dataset & parameter_rows$model == key$model & parameter_rows$parameter == "tau2_B", ]
    y <- unlist(z[c("lower", "upper", "posterior_mean", "median")])
    y <- y[is.finite(y) & y > 0]
    lim <- if (length(y)) range(c(y, .001, .1)) else c(1e-5, 1)
    plot(NA, xlim = c(.5, 3.5), ylim = lim, log = "y", xaxt = "n", xlab = "Observation-noise prior", ylab = "Areal variance (log scale)", main = label_panel(key))
    axis(1, at = 1:3, labels = c("A: original", "B: PC 0.316", "C: PC 0.5"))
    abline(h = .1, lty = 2, col = "#6B7280"); abline(h = .001, lty = 3, col = "#B91C1C")
    if (nrow(z)) {
      x <- job_x(z); color <- colors[z$prior_id]
      segments(x, z$lower, x, z$upper, col = color, lwd = 1.8)
      segments(x - .025, z$lower, x + .025, z$lower, col = color)
      segments(x - .025, z$upper, x + .025, z$upper, col = color)
      points(x, z$posterior_mean, pch = 16, col = color, cex = 1)
      points(x, z$median, pch = 4, col = color, cex = 1.2, lwd = 1.7)
    }
  }
  mtext(paste("Areal variance posterior summaries |", plot_title), outer = TRUE, side = 3, line = 1.5, cex = 1.2, font = 2)
  mtext("Circle: mean; cross: median; line: 95% credible interval; dashed: truth 0.1; dotted: boundary 0.001", outer = TRUE, side = 1, line = 1, cex = .9)
  mtext("Within each prior, left to right: initial areal variance 0.001, 0.1, 0.4", outer = TRUE, side = 1, line = 2.5, cex = .9)
  dev.off()
}

png(file.path(output, "latent_predictive_scores.png"), width = 2100, height = 1700, res = 160)
par(mfrow = c(4, 2), mar = c(3.5, 5, 2.6, 1), oma = c(4, 0, 3, 0), las = 1)
for (i in seq_len(nrow(panel_plan))) {
  key <- panel_plan[i, ]; z <- jobs[jobs$dataset == key$dataset & jobs$model == key$model, ]
  for (metric in c("MSE", "CRPS_Z_joint")) {
    y <- z[[metric]]
    error <- if (metric == "CRPS_Z_joint") 2 * z$CRPS_MCSE_approx else rep(0, nrow(z))
    candidates <- c(y - error, y + error); candidates <- candidates[is.finite(candidates)]
    lim <- if (length(candidates)) range(candidates) else c(0, 1)
    if (diff(lim) == 0) lim <- lim + c(-1, 1) * max(.001, abs(lim[1]) * .02)
    lim <- lim + c(-1, 1) * diff(lim) * .12
    plot(NA, xlim = c(.5, 3.5), ylim = lim, xaxt = "n", xlab = "", ylab = if (metric == "MSE") "Latent MSE" else "Latent joint CRPS", main = label_panel(key))
    axis(1, at = 1:3, labels = c("A: original", "B: PC 0.316", "C: PC 0.5"))
    x <- job_x(z)
    if (metric == "CRPS_Z_joint") segments(x, y - error, x, y + error, col = colors[z$prior_id], lwd = 1.5)
    points(x, y, pch = c(16, 17, 15)[match(z$start_area, c(.001, .1, .4))], col = colors[z$prior_id], cex = 1.1)
  }
}
mtext(paste("Latent-surface prediction |", plot_title), outer = TRUE, side = 3, line = 1, cex = 1.3, font = 2)
mtext("Within each prior, left to right: starts 0.001 (circle), 0.1 (triangle), 0.4 (square)", outer = TRUE, side = 1, line = 1, cex = .95)
mtext("CRPS bars: +/- 2 approximate batch MCSE; diagnostic only, not a credible interval or across-replicate uncertainty", outer = TRUE, side = 1, line = 2.5, cex = .9)
dev.off()

columns <- c("all_fitted_and_scored", "required_parameters_valid", "mode_status_pass", "Hessian_pass", "integration_pass", "log_mode_spread_pass", "parameter_stability_pass", "prediction_stability_pass", "boundary_stability_pass")
labels <- c("Fit / score", "Parameters\nvalid", "Mode = 0", "Hessian", "Integration", "Log-mode\nspread", "Parameter\nshift", "Prediction\nRMS", "Boundary\nprobability")
png(file.path(output, "initialization_stability.png"), width = 2100, height = 1500, res = 160)
par(mar = c(8, 17, 5, 1), xpd = NA)
plot(NA, xlim = c(.5, length(columns) + .5), ylim = c(.5, nrow(groups) + .5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
for (i in seq_len(nrow(groups))) for (j in seq_along(columns)) {
  pass <- isTRUE(groups[[columns[j]]][i]); y <- nrow(groups) + 1L - i
  rect(j - .46, y - .43, j + .46, y + .43, col = if (pass) "#CBE9DF" else "#F5C9C2", border = "white")
  text(j, y, if (pass) "PASS" else "REVIEW", cex = .8, col = if (pass) "#065F46" else "#9A3412")
}
axis(1, at = seq_along(columns), labels = labels, tick = FALSE, cex.axis = .9)
axis(2, at = rev(seq_len(nrow(groups))), labels = paste(ifelse(groups$dataset == "baseline", "Baseline", "Clustering delta 2"), groups$model, paste0("prior ", groups$prior_id), sep = " | "), tick = FALSE, las = 1, cex.axis = .85)
title(main = paste("Initialization stability screens |", plot_title), line = 2)
mtext("Each row includes all three starts; REVIEW also includes missing or invalid evidence", side = 1, line = 5, cex = .95)
mtext("Thresholds are prospective practical screens, not a proof of numerical convergence", side = 1, line = 6.5, cex = .95)
dev.off()

status <- c(paste("Planned attempts:", nrow(plan)), paste("Recorded attempts:", sum(jobs$recorded)),
  paste("Successful fits:", sum(jobs$fit_success)), paste("Successfully scored:", sum(jobs$score_success)),
  paste("Groups passing every prospective screen:", sum(groups$all_components_pass), "of", nrow(groups)),
  "", "All priors and starts are retained; no fit was selected to replace an original fit",
  "Do not compare log-posterior-mode values across priors; only within-prior spreads are used",
  "Matched prior contrasts retain initialization dependence; starts are not independent simulation replicates",
  "CRPS MCSE and paired contrast MCSE are approximate 10-batch proxies for 1000-draw Monte Carlo error",
  "The fitted clustering model omits U(s); its gamma is not interpreted as a correctly specified recovery target",
  "The group pass flag requires all 3 attempted fits, valid summaries, and every listed component screen",
  "Existing summary_valid and coherence_warning flags are respected; flagged numeric values remain in the comparisons",
  "A passing group is not a convergence guarantee; no production-prior recommendation is automated")
writeLines(status, file.path(output, "comparison_status.txt"))
cat(paste(status[1:5], collapse = "\n"), "\n")
