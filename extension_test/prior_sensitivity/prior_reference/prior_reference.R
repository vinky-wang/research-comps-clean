# Observation-noise priors on precision q = 1 / v, with v = tau^2
# Run with the existing R 4.2 launcher and INLA 23.04.24
# Optional first argument selects the output directory

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else file.path(
  "/Users/vinky/research-comps-clean/extension_test/prior_sensitivity",
  "prior_reference"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expected_inla <- "23.04.24"
installed_inla <- as.character(utils::packageVersion("INLA"))
stopifnot(utils::compareVersion(installed_inla, expected_inla) == 0L)

priors <- data.frame(
  prior_id = c("A", "B", "C"),
  label = c("A: Gamma(1, 0.00005) precision", "B: PC(0.316, 0.5)", "C: PC(0.5, 0.5)"),
  inla_prior = c("loggamma", "pc.prec", "pc.prec"),
  parameter_1 = c(1, 0.316, 0.5),
  parameter_2 = c(0.00005, 0.5, 0.5),
  stringsAsFactors = FALSE
)

variance_cdf <- function(v, prior) {
  if (prior$prior_id == "A") {
    return(stats::pgamma(1 / v, shape = 1, rate = 0.00005, lower.tail = FALSE))
  }
  lambda <- -log(prior$parameter_2) / prior$parameter_1
  -expm1(-lambda * sqrt(v))
}

variance_quantile <- function(p, prior) {
  if (prior$prior_id == "A") {
    return(1 / stats::qgamma(1 - p, shape = 1, rate = 0.00005))
  }
  lambda <- -log(prior$parameter_2) / prior$parameter_1
  (stats::qexp(p, rate = lambda))^2
}

log10_variance_density <- function(x, prior) {
  v <- 10^x
  if (prior$prior_id == "A") {
    log_density_v <- stats::dgamma(1 / v, shape = 1, rate = 0.00005, log = TRUE) -
      2 * log(v)
  } else {
    lambda <- -log(prior$parameter_2) / prior$parameter_1
    log_density_v <- log(lambda) - log(2) - log(v) / 2 - lambda * sqrt(v)
  }
  # The Jacobian makes this the density of log10(v), not the density of v
  exp(log_density_v + log(log(10)) + log(v))
}

# Check the PC formulas against the installed INLA implementation
test_variance <- 10^seq(-7, 2, length.out = 101)
test_probability <- c(0.025, 0.5, 0.975)
checks <- do.call(rbind, lapply(2:3, function(i) {
  p <- priors[i, ]
  analytic_cdf <- variance_cdf(test_variance, p)
  inla_cdf <- 1 - INLA::inla.pc.pprec(
    1 / test_variance, u = p$parameter_1, alpha = p$parameter_2
  )
  analytic_quantile <- variance_quantile(test_probability, p)
  inla_quantile <- 1 / INLA::inla.pc.qprec(
    1 - test_probability, u = p$parameter_1, alpha = p$parameter_2
  )
  analytic_log_density <- log10_variance_density(log10(test_variance), p)
  inla_log_density <- INLA::inla.pc.dprec(
    1 / test_variance, u = p$parameter_1, alpha = p$parameter_2
  ) * log(10) / test_variance
  data.frame(
    prior_id = p$prior_id,
    max_cdf_error = max(abs(analytic_cdf - inla_cdf)),
    max_quantile_error = max(abs(analytic_quantile - inla_quantile)),
    max_log10_density_error = max(abs(analytic_log_density - inla_log_density)),
    check_passed = isTRUE(all.equal(analytic_cdf, inla_cdf)) &&
      isTRUE(all.equal(analytic_quantile, inla_quantile)) &&
      isTRUE(all.equal(analytic_log_density, inla_log_density))
  )
}))
stopifnot(all(checks$check_passed))

summary <- do.call(rbind, lapply(seq_len(nrow(priors)), function(i) {
  p <- priors[i, ]
  lambda <- if (p$prior_id == "A") NA_real_ else -log(p$parameter_2) / p$parameter_1
  quantiles <- variance_quantile(c(0.025, 0.5, 0.975), p)
  cbind(p, data.frame(
    sd_exponential_rate = lambda,
    variance_mean = if (p$prior_id == "A") Inf else 2 / lambda^2,
    variance_median = quantiles[2],
    variance_q025 = quantiles[1],
    variance_q975 = quantiles[3],
    probability_variance_lt_0p001 = variance_cdf(0.001, p),
    probability_variance_lt_0p01 = variance_cdf(0.01, p),
    probability_variance_lt_0p1 = variance_cdf(0.1, p),
    probability_variance_gt_1 = 1 - variance_cdf(1, p)
  ))
}))
utils::write.csv(summary, file.path(output_dir, "prior_summary.csv"), row.names = FALSE)
utils::write.csv(checks, file.path(output_dir, "inla_formula_checks.csv"), row.names = FALSE)

x <- seq(-7, 2, length.out = 1801)
curves <- do.call(rbind, lapply(seq_len(nrow(priors)), function(i) {
  p <- priors[i, ]
  data.frame(
    prior_id = p$prior_id,
    log10_variance = x,
    variance = 10^x,
    log10_variance_density = log10_variance_density(x, p),
    variance_cdf = variance_cdf(10^x, p)
  )
}))
utils::write.csv(curves, file.path(output_dir, "prior_curves.csv"), row.names = FALSE)

colors <- c("#3B4B63", "#007F86", "#B05B00")
draw_figure <- function() {
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(0.9, 0.1))
  par(mar = c(5.1, 5.1, 3.5, 1), oma = c(0, 0, 2.2, 0),
      mgp = c(3.1, 0.8, 0), las = 1, cex = 0.95)
  plot(NA, xlim = range(x), ylim = c(0, 0.95), xlab = expression(log[10](v)),
       ylab = expression(paste("Density of ", log[10](v))),
       main = "Prior density on log-variance scale")
  abline(v = log10(c(0.001, 0.1)), lty = c(2, 3), col = "#999999")
  for (i in seq_len(nrow(priors))) {
    lines(x, log10_variance_density(x, priors[i, ]), col = colors[i], lwd = 2.7)
  }
  legend("topright", legend = c("v = 0.001", "v = 0.1"), lty = c(2, 3),
         col = "#777777", bty = "n", cex = 0.85)
  plot(NA, xlim = range(x), ylim = c(0, 1), xlab = expression(log[10](v)),
       ylab = expression(Pr(tau^2 <= v)), main = "Prior cumulative probability")
  abline(v = log10(c(0.001, 0.1)), lty = c(2, 3), col = "#999999")
  abline(h = 0.5, col = "#DDDDDD")
  for (i in seq_len(nrow(priors))) {
    lines(x, variance_cdf(10^x, priors[i, ]), col = colors[i], lwd = 2.7)
    points(-3, variance_cdf(0.001, priors[i, ]), col = colors[i], pch = 16)
  }
  mtext("Observation variance priors: v = tau^2 = 1/q", outer = TRUE,
        side = 3, line = 0.3, cex = 1.25, font = 2)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = priors$label, col = colors, lwd = 2.7,
         horiz = TRUE, bty = "n", cex = 0.94)
}
grDevices::png(file.path(output_dir, "observation_variance_priors.png"),
              width = 2000, height = 1050, res = 160)
draw_figure()
grDevices::dev.off()
grDevices::pdf(file.path(output_dir, "observation_variance_priors.pdf"),
              width = 12.5, height = 6.5625)
draw_figure()
grDevices::dev.off()

writeLines(c(
  "Observation-noise prior reference",
  "",
  paste("Verified against installed INLA", installed_inla, "using inla.pc.pprec, inla.pc.qprec and inla.pc.dprec"),
  "This is an analytic prior comparison, without model fitting or data conditioning",
  "Each prior is applied independently to point and areal observation precisions",
  "",
  "Let q = 1/v and v = tau^2 be the observation variance",
  "A: loggamma c(1, 0.00005) means q ~ Gamma(shape = 1, rate = 0.00005)",
  "   v has inverse-Gamma(shape = 1, scale = 0.00005) density",
  "   F(v) = exp(-0.00005/v), so the variance mean is infinite",
  "B/C: pc.prec c(u, alpha) means P(tau > u) = alpha",
  "   tau ~ Exponential(lambda), lambda = -log(alpha)/u",
  "   F(v) = 1 - exp(-lambda*sqrt(v)); E(v) = 2/lambda^2",
  "   With alpha = 0.5, median(v) = u^2",
  "   B has exact median 0.316^2 = 0.099856; C has exact median 0.25",
  "",
  sprintf("P(v < 0.001): A = %.6f, B = %.6f, C = %.6f",
          summary$probability_variance_lt_0p001[1],
          summary$probability_variance_lt_0p001[2],
          summary$probability_variance_lt_0p001[3]),
  "PC precision priors retain a base model at zero observation SD",
  "Their variance densities also diverge at zero, but their integrated probability",
  "below 0.001 is far smaller than under A; density height is not probability mass",
  "A concentrates much more prior probability near very small variance while retaining",
  "a heavy upper tail with no finite mean; this is distinct from a posterior diagnosis",
  "",
  "Figure left: density of x = log10(v), including the Jacobian log(10)*v",
  "It is not a variance-density curve displayed on a logarithmic axis",
  "Figure right: variance CDF evaluated at v = 10^x",
  "Full numeric summaries and formula-check errors are in the accompanying CSV files"
), file.path(output_dir, "prior_reference.txt"))

print(summary)
print(checks)
cat("Saved prior reference to:", normalizePath(output_dir), "\n")
