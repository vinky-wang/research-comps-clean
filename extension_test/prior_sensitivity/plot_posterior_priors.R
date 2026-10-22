# Plot saved posterior marginals without fitting or selecting among attempts
# Arguments: study directory, figure directory, common areal start
args <- commandArgs(trailingOnly = TRUE)
study_dir <- if (length(args)) args[[1]] else
  "/Users/vinky/research-comps-clean/extension_test/prior_sensitivity"
figure_dir <- if (length(args) >= 2L) args[[2]] else file.path(study_dir, "reports")
common_start <- if (length(args) >= 3L) as.numeric(args[[3]]) else 0.1
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(utils::compareVersion(as.character(utils::packageVersion("INLA")), "23.04.24") == 0L)

plan <- utils::read.csv(file.path(study_dir, "plan.csv"), stringsAsFactors = FALSE)
log_v <- seq(-6, 1, length.out = 1401)
v <- 10^log_v
curves <- list()
availability <- list()

for (i in seq_len(nrow(plan))) {
  job <- plan[i, ]
  checkpoint <- file.path(study_dir, "checkpoints", paste0(job$job_id, ".rds"))
  status <- "not yet recorded"
  if (file.exists(checkpoint)) {
    record <- readRDS(checkpoint)
    status <- if (isTRUE(record$fit_success)) "fit saved" else "fit failed"
    if (isTRUE(record$fit_success)) {
      bundle <- readRDS(record$fit_file)
      for (component in c("point", "area")) {
        precision_name <- bundle$likelihood_map$precision_name[
          bundle$likelihood_map$component == component
        ]
        stopifnot(length(precision_name) == 1L)
        marginal <- bundle$fit$internal.marginals.hyperpar[[
          sub("^Precision", "Log precision", precision_name)
        ]]
        stopifnot(is.matrix(marginal), ncol(marginal) == 2L)
        # theta = log(q), x = log10(v) = -theta / log(10)
        # Density of x includes |d theta/d x| = log(10)
        density <- log(10) * INLA::inla.dmarginal(-log(10) * log_v, marginal)
        # The monotone transformation gives P(V <= v) = 1 - F_theta(-log(v))
        cdf <- pmin(1, pmax(0, 1 - INLA::inla.pmarginal(-log(v), marginal)))
        curves[[length(curves) + 1L]] <- data.frame(
          job, component = component, log10_variance = log_v,
          log10_variance_density = density, variance_cdf = cdf,
          mode_status = record$diagnostics$mode_status, row.names = NULL
        )
      }
      rm(bundle)
      gc(verbose = FALSE)
    }
  }
  availability[[i]] <- data.frame(job, plot_status = status)
}

availability <- do.call(rbind, availability)
utils::write.csv(availability, file.path(figure_dir, "posterior_plot_availability.csv"), row.names = FALSE)
if (!length(curves)) stop("No completed saved fits available to plot")
curves <- do.call(rbind, curves)
utils::write.csv(curves, file.path(figure_dir, "posterior_plot_curves.csv"), row.names = FALSE)

prior_ids <- c("A", "B", "C")
prior_colors <- c(A = "#3B4B63", B = "#007F86", C = "#B05B00")
prior_labels <- c("A: Gamma(1, 0.00005)", "B: PC(0.316, 0.5)", "C: PC(0.5, 0.5)")
starts <- c(0.001, 0.1, 0.4)
start_colors <- c("#97366B", "#147A8C", "#AD6719")
models <- c("Melding", "PSmelding")
datasets <- c(baseline = "Baseline", clustering_delta2 = "Independent clustering: delta = 2")

begin_figure <- function(path, ncol, height = 1450) {
  grDevices::png(path, width = if (ncol == 2) 1900 else 2350, height = height, res = 160)
  layout(rbind(matrix(seq_len(2L * ncol), 2L, byrow = TRUE),
               rep(2L * ncol + 1L, ncol)), heights = c(0.44, 0.44, 0.12))
  par(mar = c(4.1, 4.6, 2.8, 1), oma = c(0, 0, 3.7, 0),
      mgp = c(2.8, 0.8, 0), las = 1, cex = 0.95)
}

finish_figure <- function(title, subtitle, labels, colors, lty) {
  mtext(title, side = 3, outer = TRUE, line = 2, font = 2, cex = 1.2)
  mtext(subtitle, side = 3, outer = TRUE, line = 0.6, cex = 0.9)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = labels, col = colors, lty = lty, lwd = 2.6,
         horiz = TRUE, bty = "n", cex = 0.95)
  text(0.5, 0.13, "All completed attempts retained, including numerically flagged fits",
       cex = 0.85)
  grDevices::dev.off()
}

for (dataset in names(datasets)) {
  common <- curves[curves$dataset == dataset & curves$start_area == common_start, ]
  begin_figure(file.path(figure_dir, paste0("posterior_noise_density_", dataset, ".png")), 2)
  for (model in models) {
    for (component in c("point", "area")) {
      panel <- common[common$model == model & common$component == component, ]
      upper <- if (nrow(panel)) max(panel$log10_variance_density, na.rm = TRUE) * 1.08 else 1
      plot(NA, xlim = range(log_v), ylim = c(0, max(upper, 0.01)),
           xlab = expression(log[10](tau^2)),
           ylab = expression(paste("Density of ", log[10](tau^2))),
           main = paste(model, "-", if (component == "point") "Point variance" else "Areal variance"))
      abline(v = c(-3, -1), lty = c(2, 3), col = "#AAAAAA")
      for (prior in prior_ids) {
        d <- panel[panel$prior_id == prior, ]
        if (nrow(d)) lines(d$log10_variance, d$log10_variance_density,
                          col = prior_colors[prior], lwd = 2.6)
      }
      missing <- setdiff(prior_ids, unique(panel$prior_id))
      if (length(missing)) mtext(paste("Pending/missing:", paste(missing, collapse = ", ")),
                                side = 3, line = 0, cex = 0.7, col = "#777777")
    }
  }
  finish_figure(
    paste(datasets[[dataset]], "- Observation-noise posteriors"),
    paste0("Common areal variance start = ", common_start,
           "; density includes log-scale Jacobian; dashed v = 0.001, dotted truth v = 0.1"),
    prior_labels, prior_colors, rep(1, 3)
  )

  area <- curves[curves$dataset == dataset & curves$component == "area", ]
  begin_figure(file.path(figure_dir, paste0("posterior_areal_cdf_starts_", dataset, ".png")), 3)
  for (model in models) {
    for (prior in prior_ids) {
      panel <- area[area$model == model & area$prior_id == prior, ]
      plot(NA, xlim = range(log_v), ylim = c(0, 1),
           xlab = expression(log[10](v)), ylab = expression(Pr(tau[B]^2 <= v)),
           main = paste(model, "- Prior", prior))
      abline(v = c(-3, -1), lty = c(2, 3), col = "#AAAAAA")
      abline(h = 0.5, col = "#EEEEEE")
      for (j in seq_along(starts)) {
        d <- panel[panel$start_area == starts[j], ]
        if (nrow(d)) lines(d$log10_variance, d$variance_cdf,
                          col = start_colors[j], lty = j, lwd = 2.5)
      }
      missing <- setdiff(starts, unique(panel$start_area))
      if (length(missing)) mtext(paste("Pending/missing:", paste(missing, collapse = ", ")),
                                side = 3, line = 0, cex = 0.7, col = "#777777")
    }
  }
  finish_figure(
    paste(datasets[[dataset]], "- Areal variance across initializations"),
    "Exact monotone transformation of saved log-precision marginal; dashed v = 0.001, dotted truth v = 0.1",
    paste("Areal variance start =", starts), start_colors, seq_along(starts)
  )
}

writeLines(c(
  "Posterior plots use completed saved attempts only; no model is refitted",
  "No attempt is removed based on mode status, Hessian diagnostics or posterior values",
  paste("Common-start prior comparison uses areal variance start", common_start),
  "Every point variance start is 0.1",
  "Noise density panels show x = log10(v), not a variance density with log-scaled axes",
  "With theta = log(q), f_x(x) = log(10) * f_theta(-log(10)*x)",
  "CDF panels use P(V <= v) = 1 - F_theta(-log(v)), requiring no density Jacobian",
  "Both transformations use INLA interpolation of the saved marginal, not transformed posterior means",
  "Curves show variances from 10^-6 through 10^1; mass outside this window is not renormalized",
  "Numerically invalid or poorly approximated saved marginals remain numerically questionable after transformation",
  "Missing attempts are labelled and listed in posterior_plot_availability.csv",
  "All plotted values are exported in posterior_plot_curves.csv"
), file.path(figure_dir, "posterior_plot_notes.txt"))
cat("Completed fit curves:", sum(availability$plot_status == "fit saved"), "/", nrow(plan), "\n")
cat("Figures saved to:", normalizePath(figure_dir), "\n")
