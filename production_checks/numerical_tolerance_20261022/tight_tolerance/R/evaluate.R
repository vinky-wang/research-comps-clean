# Evaluation helpers ----------------------------------------------------------
# This file contains its own numerical helpers; no old project files are sourced
# Checked transformations are retained from psmelding_numerics.R, with scoring
# functions adapted from replicate_psmelding_corrected_metrics.R
# Public entry points: evaluate_fit() and summarize_metrics() below; numerical helpers follow
# [EXTENSION] CRPS targets the latent response surface Z(s) only
# CRPS_Z_gaussian uses a Gaussian approximation to each predictive marginal
# CRPS_Z_joint uses joint posterior draws, retaining hyperparameter uncertainty
# No future noisy responses are generated or scored

# 1. Evaluate one saved model -------------------------------------------------
# Order: align predictions -> parameter summaries -> marginal scores -> joint
# scores -> review flags. A score is invalid only if its own required inputs
# fail; unrelated warnings remain visible without deleting usable scores
evaluate_fit <- function(bundle, truth, config, seed = 1L) {
  alignment <- check_prediction_alignment(bundle, truth)
  if (!alignment$coordinates_valid || !alignment$truth_valid) {
    stop("Predictions and truth must have identical finite coordinates and ordering.")
  }
  parameters <- score_parameters(bundle, truth,
    list(nu = config$generating$nu, score_conditioned_alpha = FALSE))
  # [EXTENSION] Preserve fitted gamma summaries under misspecification, but do not
  # invent a scalar truth for unequal left/right preferential coefficients
  if (!is.null(truth$sampling) && nzchar(truth$sampling$gamma_note)) {
    selected <- parameters$parameter %in% c("gamma", "native:Beta for j")
    parameters$note[selected] <- truth$sampling$gamma_note
  }
  native <- check_native_parameter_summaries(bundle$fit)
  prediction <- bundle$predictions
  Z <- truth$Z
  n <- length(Z)
  means <- prediction$pred_mean
  medians <- prediction$pred_med
  sds <- prediction$sd
  rows <- pointwise <- list()
  # Evaluate each score independently. Store its failure reason beside NA
  # instead of turning a partial scoring failure into a missing replicate
  record <- function(name, expression) {
    reason <- ""
    values <- tryCatch(
      force(expression),
      error = function(error) {
        reason <<- conditionMessage(error)
        rep(NA_real_, n)
      }
    )
    valid <- is.numeric(values) && length(values) == n && all(is.finite(values))
    if (!valid && !nzchar(reason)) {
      reason <- "Nonfinite pointwise score or invalid required inputs."
    }
    pointwise[[name]] <<- if (valid) values else rep(NA_real_, n)
    rows[[name]] <<- data.frame(
      region = "all",
      n_locations = n,
      metric = name,
      value = if (valid) mean(values) else NA_real_,
      valid = valid,
      issue = reason
    )
  }
  record("MSE", {
    stopifnot(alignment$mean_valid)
    (means - Z)^2
  })
  record("MAE", {
    stopifnot(alignment$median_valid)
    abs(medians - Z)
  })
  # [ORIGINAL] [EXTENSION] Keep the paper's squared W2 definition, using standard deviations
  record("WD_paper", {
    stopifnot(alignment$mean_valid, alignment$sd_valid)
    (means - Z)^2 + (sds - sqrt(truth$parameters$tau2_s))^2
  })
  record("CRPS_Z_gaussian", {
    stopifnot(alignment$mean_valid, alignment$sd_valid)
    gaussian_crps(means, sds, Z)
  })
  lower <- prediction$pred_ll
  upper <- prediction$pred_ul
  interval_valid <- length(lower) == n && length(upper) == n &&
    all(is.finite(lower)) && all(is.finite(upper)) && all(lower <= upper)
  record("Z_coverage95", {
    stopifnot(interval_valid)
    as.numeric(lower <= Z & Z <= upper)
  })
  record("Z_IS95", {
    stopifnot(interval_valid)
    interval_score(lower, upper, Z)
  })

  # [EXTENSION] Draw the latent surface with hyperparameter uncertainty
  # Measurement noise is estimated in fitting, but no extra Y draws are needed
  drawing_error <- NULL
  record("CRPS_Z_joint", {
    Z_draws <- tryCatch(
      draw_latent_predictions(bundle, config$scoring$joint_draws, seed)$Z,
      error = function(error) {
        drawing_error <<- conditionMessage(error)
        stop(error)
      }
    )
    if (nrow(Z_draws) != n) {
      stop("Latent posterior draws and truth must have the same number of locations.")
    }
    sample_crps(Z_draws, Z)
  })
  metrics <- do.call(rbind, rows)
  coordinates <- bundle$prediction_coordinates
  regions <- rep("all", n)
  if (identical(truth$sampling$mechanism, "piecewise")) {
    regions <- ifelse(coordinates[, 1] < truth$sampling$split_x, "left", "right")
    # Reuse the pointwise scores and joint draws for both halves of the domain
    for (region in c("left", "right")) {
      selected <- regions == region
      regional <- do.call(rbind, lapply(names(pointwise), function(name) {
        values <- pointwise[[name]][selected]
        valid <- length(values) > 0L && all(is.finite(values))
        data.frame(region = region, n_locations = sum(selected), metric = name,
          value = if (valid) mean(values) else NA_real_, valid = valid,
          issue = if (valid) "" else "No valid pointwise scores for this region")
      }))
      metrics <- rbind(metrics, regional)
    }
  }
  pointwise <- data.frame(
    x = coordinates[, 1], y = coordinates[, 2], region = regions,
    truth = Z, mean_error = means - Z, median_error = medians - Z,
    as.data.frame(pointwise)
  )
  rownames(metrics) <- NULL
  mode_status <- bundle$fit$mode$mode.status
  if (length(mode_status) != 1L) {
    mode_status <- NA_real_
  }
  hessian_warning <- any(grepl("ensure spd for Hessian", bundle$fit$logfile))
  needs_review <- !isTRUE(mode_status == 0) || hessian_warning ||
    !alignment$valid || any(!parameters$summary_valid | parameters$coherence_warning) ||
    any(!native$summary_valid | native$coherence_warning) || any(!metrics$valid)
  metrics$needs_review <- needs_review
  # Review flags do not exclude otherwise valid quantities from primary summaries
  list(
    metrics = metrics,
    pointwise = pointwise,
    parameters = parameters,
    diagnostics = list(
      mode_status = mode_status,
      needs_review = needs_review,
      prediction = alignment,
      hessian_adjustment_reported = hessian_warning,
      joint_drawing_error = drawing_error,
      noise_marginals = attr(parameters, "noise_diagnostics")
    ),
    metadata = list(
      seed = seed,
      joint_draws = config$scoring$joint_draws,
      scoring = "Joint-draw CRPS uses an unbiased pairwise-distance estimator; finite-draw scores can be negative.",
      alpha = "Not assessed against a generating value under fixed-count sampling."
    )
  )
}

# 2. Summarize replications ---------------------------------------------------
# MCSE = sample SD / sqrt(number of valid replicate scores)
# The 2.5%/97.5% columns describe the spread of scores across replications;
# they are not posterior credible intervals or confidence limits for the mean
summarize_metrics <- function(metrics, n_requested) {
  keys <- c(
    "experiment",
    "rho",
    "I",
    "J",
    "sampling_id",
    "sampling_type",
    "gamma",
    "gamma_L",
    "gamma_R",
    "split_x",
    "delta",
    "rho_U",
    "sigma2_U",
    "log_intensity_phi_cor",
    "noise_id",
    "tau2_s",
    "tau2_B",
    "model",
    "region",
    "metric"
  )
  # NA marks inapplicable parameters (e.g. a scalar gamma in a piecewise case)
  # paste keeps those groups; interaction() would silently discard NA rows
  groups <- split(seq_len(nrow(metrics)), do.call(paste, c(metrics[keys], sep = "|")))
  do.call(
    rbind,
    lapply(
      groups,
      function(index) {
        x <- metrics[index, , drop = FALSE]
        values <- x$value[x$valid & is.finite(x$value)]
        count <- length(values)
        cbind(
          x[1L, keys, drop = FALSE],
          n_requested = n_requested,
          n_scored = nrow(x),
          n_valid = count,
          n_invalid = nrow(x) - count,
          n_flagged = sum(x$needs_review),
          mean = if (count) {
            mean(values)
          } else {
            NA_real_
          },
          mcse = if (count > 1L) {
            stats::sd(values) / sqrt(count)
          } else {
            NA_real_
          },
          lower025 = if (count) {
            unname(stats::quantile(values, .025))
          } else {
            NA_real_
          },
          upper975 = if (count) {
            unname(stats::quantile(values, .975))
          } else {
            NA_real_
          }
        )
      }
    )
  )
}


# 3. Convert Gaussian precision posteriors to noise variances -----------------
# Keep the numerical checks in this section: inverse moments can be sensitive
# to tiny low-precision tails. Finite INLA summary values alone do not establish
# that E[1/q] or Var(1/q) has been computed reliably
# Work on eta = log(precision), so tau^2 = exp(-eta). All variance moments
# and quantiles below use one normalized density on the saved eta support
# These helpers never refit a model, draw samples, or change a fitted object

# With eta = log(q), measurement-error variance is exp(-eta), not exp(eta)
variance_from_log_precision <- function(eta) {
  if (!is.numeric(eta) || !length(eta) || any(!is.finite(eta))) {
    stop("Log precision must contain finite numeric values.")
  }
  value <- exp(-eta)
  if (any(!is.finite(value)) || any(value <= 0)) {
    stop("Log precision cannot be represented as a positive finite variance.")
  }
  value
}

# Represent an unusable marginal explicitly while retaining diagnostic details
invalid_noise_summary <- function(reason, diagnostics = list(), probs = c(.025, .5, .975)) {
  list(
    summary = setNames(
      rep(NA_real_, 5L),
      c("posterior_mean", "posterior_sd", "lower", "median", "upper")
    ),
    quantiles = rep(NA_real_, length(probs)),
    mode = NA_real_,
    valid = FALSE,
    reason = reason,
    diagnostics = diagnostics
  )
}

# Validate the saved density before interpolation; zero endpoints require
# separate tail bounds, while zero densities inside its support are rejected
validate_noise_marginal <- function(marginal, name, allow_endpoint_zeros = FALSE) {
  if (is.list(marginal) && all(c("x", "y") %in% names(marginal))) {
    marginal <- cbind(marginal$x, marginal$y)
  }
  if (!is.matrix(marginal) || ncol(marginal) != 2L || nrow(marginal) < 3L ||
      any(!is.finite(marginal)) || any(diff(marginal[, 1L]) <= 0) ||
      any(marginal[, 2L] < 0) || !any(marginal[, 2L] > 0)) {
    stop("Invalid saved ", name, " marginal: require ordered finite coordinates and nonnegative density.")
  }
  positive <- which(marginal[, 2L] > 0)
  if (length(positive) < 3L || any(marginal[seq.int(min(positive), max(positive)), 2L] <= 0)) {
    stop("Invalid saved ", name, " marginal: interior zero density or fewer than three positive knots.")
  }
  if (!allow_endpoint_zeros && any(marginal[, 2L] == 0)) {
    stop("Zero endpoint density requires explicit tail validation.")
  }
  marginal
}

# Bound the mass and inverse moments contributed by saved zero-endpoint
# intervals. These bounds concern the saved support, not unknown external tails
bound_zero_density_tails <- function(marginal) {
  positive <- which(marginal[, 2L] > 0)
  first <- min(positive)
  last <- max(positive)
  bounds <- list()
  if (first > 1L) {
    bounds$left <- list(
      interval = marginal[c(1L, first), 1L],
      endpoint_density = marginal[first, 2L],
      log_integral_upper_bound = log(marginal[first, 1L] - marginal[1L, 1L]) +
        log(marginal[first, 2L]) - (0:2) * marginal[1L, 1L]
    )
  }
  if (last < nrow(marginal)) {
    bounds$right <- list(
      interval = marginal[c(last, nrow(marginal)), 1L],
      endpoint_density = marginal[last, 2L],
      log_integral_upper_bound = log(marginal[nrow(marginal), 1L] - marginal[last, 1L]) +
        log(marginal[last, 2L]) - (0:2) * marginal[last, 1L]
    )
  }
  list(
    positive = positive,
    removed_indices = which(marginal[, 2L] == 0),
    tails = bounds,
    relative_tolerance = 1e-12,
    assumption = paste(
      "On an interval between a saved zero-density endpoint and the adjacent positive knot,",
      "use a nonnegative monotone boundary connection bounded above by that positive density.",
      "The bound applies to this explicit boundary representation inside saved support,",
      "not to unknown posterior mass outside the saved marginal. No positive density floor is added."
    )
  )
}

interpolate_log_precision_density <- function(marginal, factor = 60L) {
  # Match INLA's inla.smarginal() log-density interpolation: adapt the
  # interpolation axis to the saved knot spacing, then map it back to eta
  x <- marginal[, 1L]
  y <- marginal[, 2L]
  dx <- length(x) * diff(x) / mean(diff(x))
  axis <- c(0, cumsum(sqrt(dx)))
  axis <- min(x) + diff(range(x)) * axis / max(axis)
  from_axis <- stats::splinefun(axis, x, method = "hyman")
  log_density <- stats::splinefun(axis, log(y), method = "fmm")
  smooth_axis <- seq(min(axis), max(axis), length.out = factor * length(x))
  stats::splinefun(from_axis(smooth_axis), log_density(smooth_axis), method = "fmm")
}

integrate_noise_moments <- function(axis, log_density, powers = 0:2) {
  log_y <- log_density(axis)
  if (any(!is.finite(log_y))) {
    stop("Nonfinite interpolated log density.")
  }
  log_integrals <- vapply(
    powers,
    function(k) {
      values <- log_y - k * axis
      shift <- max(values)
      scaled <- exp(values - shift)
      integral <- sum(diff(axis) * (head(scaled, -1L) + tail(scaled, -1L)) / 2)
      if (!is.finite(integral) || integral <= 0) {
        stop("Log-density integral has no positive finite mass.")
      }
      shift + log(integral)
    },
    numeric(1)
  )
  scaled <- exp(log_y - max(log_y))
  mass <- diff(axis) * (head(scaled, -1L) + tail(scaled, -1L)) / 2
  cdf <- c(0, cumsum(mass)) / sum(mass)
  list(log_integrals = log_integrals, axis = axis, cdf = cdf)
}

noise_variance_quantiles <- function(integral, probs) {
  # The transformation is decreasing; invert at 1-p before transforming
  keep <- !duplicated(integral$cdf)
  eta <- stats::approx(
    integral$cdf[keep],
    integral$axis[keep],
    xout = 1 - probs,
    method = "linear",
    rule = 1
  )$y
  eta[probs == 0] <- max(integral$axis)
  eta[probs == 1] <- min(integral$axis)
  variance_from_log_precision(eta)
}

# Summarize the finite mixture of saved INLA configurations as a diagnostic
# comparison with the marginal-density calculation; it does not replace it
summarize_noise_configurations <- function(fit, precision_name) {
  diagnostics <- list(
    precision_name = precision_name,
    weight_method = "normalized exp(configuration log.posterior), as in INLA posterior sampler"
  )
  tryCatch(
    {
      internal_name <- sub("^Precision", "Log precision", precision_name)
      diagnostics$internal_name <- internal_name
      configs <- fit$misc$configs
      if (is.null(configs$config) || !length(configs$config) ||
          length(configs$config) != configs$nconfig) {
        stop("Missing saved joint configurations.")
      }
      transformations <- fit$misc$from.theta
      index <- match(internal_name, names(transformations))
      if (is.na(index)) {
        stop("Cannot map the precision to its saved internal transformation.")
      }
      eta <- vapply(
        configs$config,
        function(config) {
          if (!identical(names(config$theta)[index], internal_name)) {
            stop("Configuration theta names do not match the precision transformation.")
          }
          as.numeric(config$theta[index])
        },
        numeric(1)
      )
      precision <- vapply(
        eta,
        function(value) {
          as.numeric(do.call(transformations[[index]], list(value)))
        },
        numeric(1)
      )
      if (any(!is.finite(precision)) || any(precision <= 0) ||
          any(abs(log(precision) - eta) > 1e-8 * pmax(1, abs(eta)))) {
        stop("Saved from.theta is inconsistent with the log-precision mapping.")
      }
      variance <- variance_from_log_precision(eta)
      log_weights <- vapply(configs$config, function(config) config$log.posterior, numeric(1))
      if (any(is.na(log_weights)) || any(log_weights == Inf) || !any(is.finite(log_weights))) {
        stop("Invalid joint configuration weights.")
      }
      weights <- exp(log_weights - max(log_weights))
      weights <- weights / sum(weights)
      mean <- sum(weights * variance)
      sd <- sqrt(sum(weights * (variance - mean)^2))
      ordering <- order(variance)
      cdf <- cumsum(weights[ordering])
      q <- vapply(c(.025, .5, .975), function(p) variance[ordering][which(cdf >= p)[1L]], numeric(1))
      summary <- c(posterior_mean = mean, posterior_sd = sd,
        lower = q[1L], median = q[2L], upper = q[3L])
      if (any(!is.finite(summary)) || mean <= 0 || any(q <= 0)) {
        stop("Invalid configuration-weighted variance summary.")
      }
      diagnostics$config_count <- length(weights)
      diagnostics$effective_configurations <- 1 / sum(weights^2)
      diagnostics$internal_support <- range(eta)
      diagnostics$variance_support <- range(variance)
      list(
        summary = summary,
        valid = TRUE,
        reason = NA_character_,
        diagnostics = diagnostics,
        values = variance,
        weights = weights
      )
    },
    error = function(error) invalid_noise_summary(conditionMessage(error), diagnostics)
  )
}

# Main noise-summary routine: map q to eta, integrate exp(-eta) and exp(-2 eta),
# invert the CDF for quantiles, and check resolution/tail sensitivity
# E[1/q] is generally different from 1/E[q]
summarize_noise_variance <- function(fit, precision_name, probs = c(.025, .5, .975)) {
  diagnostics <- list(
    precision_name = precision_name,
    method = "normalized saved log-precision density; CDF inversion and inverse moments"
  )
  tryCatch(
    {
      if (length(precision_name) != 1L || is.na(precision_name) ||
          !grepl("^Precision for the Gaussian", precision_name)) {
        stop("Expected one named Gaussian precision marginal.")
      }
      if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
          any(probs < 0 | probs > 1)) {
        stop("Quantile probabilities must lie in [0, 1].")
      }
      # 3A. Read the saved precision density and express it on eta = log(q)
    # The Jacobian is q: f_eta(eta) = q * f_q(q)
    raw <- validate_noise_marginal(
      fit$marginals.hyperpar[[precision_name]],
      "precision",
      allow_endpoint_zeros = TRUE
    )
      if (any(raw[, 1L] <= 0)) {
        stop("Saved precision support must be positive.")
      }
      internal_name <- sub("^Precision", "Log precision", precision_name)
      diagnostics$internal_name <- internal_name
      internal <- fit$internal.marginals.hyperpar[[internal_name]]
      mapped <- cbind(log(raw[, 1L]), raw[, 1L] * raw[, 2L])
      if (is.null(internal)) {
        internal <- validate_noise_marginal(mapped, "log precision", allow_endpoint_zeros = TRUE)
        diagnostics$representation <- "log x and x*density derived exactly from saved precision knots"
        diagnostics$jacobian_max_relative_error <- 0
      } else {
        internal <- validate_noise_marginal(internal, "log precision", allow_endpoint_zeros = TRUE)
        if (!identical(dim(mapped), dim(internal))) {
          stop("Saved internal and external precision marginals have different knots.")
        }
        coordinate_error <- max(abs(mapped[, 1L] - internal[, 1L]) / pmax(1, abs(internal[, 1L])))
        if (!identical(mapped[, 2L] == 0, internal[, 2L] == 0)) {
          stop("Saved internal and external marginals have different zero-density knots.")
        }
        positive <- internal[, 2L] > 0
        density_error <- max(abs(mapped[positive, 2L] / internal[positive, 2L] - 1))
        diagnostics$jacobian_max_relative_error <- max(coordinate_error, density_error)
        if (diagnostics$jacobian_max_relative_error > 1e-8) {
          stop("Saved internal marginal fails the external-precision Jacobian check.")
        }
        diagnostics$representation <- "saved internal log-precision marginal"
      }
      diagnostics$original_internal_support <- range(internal[, 1L])
      # 3B. Bound zero-density endpoint intervals before trimming them. This
    # prevents negligible probability tails from hiding large inverse moments
    endpoint_bounds <- bound_zero_density_tails(internal)
      diagnostics$endpoint_zero_handling <- endpoint_bounds
      internal <- internal[endpoint_bounds$positive, , drop = FALSE]
      bounds <- range(internal[, 1L])
      diagnostics$internal_support <- bounds
      diagnostics$variance_support <- sort(variance_from_log_precision(bounds))
      density <- interpolate_log_precision_density(internal)
      # 3C. Integrate normalization, first and second inverse moments on two
    # grids. Compare both moments and quantiles before accepting the result
    probabilities <- c(.025, .5, .975, probs)
      coarse <- integrate_noise_moments(
        sort(unique(c(internal[, 1L], seq(bounds[1L], bounds[2L], length.out = 8193L)))),
        density
      )
      fine <- integrate_noise_moments(
        sort(unique(c(internal[, 1L], seq(bounds[1L], bounds[2L], length.out = 16385L)))),
        density
      )
      coarse_moments <- exp(coarse$log_integrals[-1L] - coarse$log_integrals[1L])
      moments <- exp(fine$log_integrals[-1L] - fine$log_integrals[1L])
      q <- noise_variance_quantiles(fine, probabilities)
      diagnostics$grid_relative_error <- max(
        abs(coarse_moments / moments - 1),
        abs(noise_variance_quantiles(coarse, probabilities) / q - 1)
      )
      if (!is.finite(diagnostics$grid_relative_error) || diagnostics$grid_relative_error > 2e-4) {
        stop("Log-precision moments or quantiles failed grid convergence.")
      }

      # Independent adaptive integration checks the same density, interval by
      # interval so narrow peaks cannot be skipped by a whole-support integral
      adaptive_log <- vapply(
        0:2,
        function(k) {
          shift <- max(density(fine$axis) - k * fine$axis)
          values <- vapply(
            seq_len(nrow(internal) - 1L),
            function(j) {
              stats::integrate(
                function(eta) exp(density(eta) - k * eta - shift),
                internal[j, 1L],
                internal[j + 1L, 1L],
                rel.tol = 1e-8,
                abs.tol = 1e-12,
                subdivisions = 200L,
                stop.on.error = TRUE
              )$value
            },
            numeric(1)
          )
          shift + log(sum(values))
        },
        numeric(1)
      )
      adaptive_moments <- exp(adaptive_log[-1L] - adaptive_log[1L])
      if (length(endpoint_bounds$removed_indices)) {
        # Each inverse moment needs its own tail check. A small probability mass
        # alone is insufficient when exp(-eta) becomes large in a left tail
        tail_logs <- lapply(endpoint_bounds$tails, function(x) x$log_integral_upper_bound)
        combined <- vapply(
          seq_along(adaptive_log),
          function(k) {
            values <- vapply(tail_logs, function(x) x[k], numeric(1))
            max(values) + log(sum(exp(values - max(values))))
          },
          numeric(1)
        )
        relative_bounds <- exp(combined - adaptive_log)
        names(relative_bounds) <- c("mass", "first_inverse_moment", "second_inverse_moment")
        diagnostics$endpoint_zero_handling$relative_integral_upper_bounds <- relative_bounds
        if (any(!is.finite(relative_bounds)) ||
            any(relative_bounds > endpoint_bounds$relative_tolerance)) {
          stop("Zero endpoint tails have material or unbounded inverse-moment mass; cannot trim safely.")
        }
        diagnostics$endpoint_zero_handling$method <- "Drop matching contiguous zero endpoints after negligible mass and inverse-moment tail bounds."
      } else {
        diagnostics$endpoint_zero_handling$method <- "No zero endpoint density."
        diagnostics$endpoint_zero_handling$relative_integral_upper_bounds <- c(
          mass = 0, first_inverse_moment = 0, second_inverse_moment = 0
        )
      }
      diagnostics$adaptive_relative_error <- max(abs(expm1(fine$log_integrals - adaptive_log)))
      diagnostics$log_normalizing_integral <- fine$log_integrals[1L]
      diagnostics$adaptive_moments <- setNames(adaptive_moments, c("mean", "second_moment"))
      if (!is.finite(diagnostics$adaptive_relative_error) || diagnostics$adaptive_relative_error > 1e-5) {
        stop("Log-precision integrals failed independent adaptive validation.")
      }
      # Check that the interpolation itself is resolved, separately from the
      # quadrature and CDF grids used to integrate that interpolation
      denser <- integrate_noise_moments(fine$axis, interpolate_log_precision_density(internal, factor = 120L))
      diagnostics$interpolation_relative_error <- max(
        abs(expm1(fine$log_integrals - denser$log_integrals)),
        abs(q / noise_variance_quantiles(denser, probabilities) - 1)
      )
      if (!is.finite(diagnostics$interpolation_relative_error) ||
          diagnostics$interpolation_relative_error > 2e-4) {
        stop("Log-precision density failed interpolation convergence.")
      }
      if (any(!is.finite(moments)) || any(moments <= 0)) {
        stop("Variance moments are not positive and finite.")
      }
      # 3D. Convert raw moments to a standard deviation. Only floating-point
    # roundoff may be clipped; a materially negative variance is an error
    posterior_variance <- moments[2L] - moments[1L]^2
      roundoff_tolerance <- 100 * .Machine$double.eps * max(moments[2L], moments[1L]^2)
      if (!is.finite(posterior_variance) || posterior_variance < -roundoff_tolerance) {
        stop("Variance second moment is inconsistent with its mean.")
      }
      diagnostics$variance_roundoff_clipped <- posterior_variance < 0
      posterior_variance <- max(0, posterior_variance)
      if (any(!is.finite(q)) || any(q <= 0) ||
          any(diff(q[order(probabilities)]) < 0) ||
          any(q < diagnostics$variance_support[1L] | q > diagnostics$variance_support[2L])) {
        stop("Variance quantiles failed support or ordering validation.")
      }
      summary <- c(
        posterior_mean = moments[1L],
        posterior_sd = sqrt(posterior_variance),
        lower = q[1L],
        median = q[2L],
        upper = q[3L]
      )

      # A variance density has the extra Jacobian exp(eta), because
      # f_tau2(exp(-eta)) = f_eta(eta) * exp(eta)
      objective <- function(eta) density(eta) + eta
      peak <- which.max(objective(fine$axis))
      search <- fine$axis[c(max(1L, peak - 1L), min(length(fine$axis), peak + 1L))]
      candidates <- c(bounds, stats::optimize(objective, search, maximum = TRUE)$maximum)
      mode <- variance_from_log_precision(candidates[which.max(objective(candidates))])
      diagnostics$boundary_density_relative_to_peak <- exp(density(bounds) - max(density(fine$axis)))
      diagnostics$native_precision_summary <- if (precision_name %in% rownames(fit$summary.hyperpar)) {
        fit$summary.hyperpar[precision_name, , drop = FALSE]
      } else {
        NULL
      }
      # 3E. Compare with the saved joint configurations for diagnosis. Keep
    # this comparison separate from the marginal summary used for reporting
    configurations <- summarize_noise_configurations(fit, precision_name)
      diagnostics$configuration_valid <- configurations$valid
      diagnostics$configuration_reason <- configurations$reason
      diagnostics$configuration_summary <- configurations$summary
      diagnostics$configuration_diagnostics <- configurations$diagnostics
      if (isTRUE(configurations$valid)) {
        means <- c(summary["posterior_mean"], configurations$summary["posterior_mean"])
        sds <- c(summary["posterior_sd"], configurations$summary["posterior_sd"])
        diagnostics$marginal_joint_mean_ratio <- max(means) / min(means)
        denominator <- sqrt(sum(sds^2))
        diagnostics$marginal_joint_standardized_difference <- if (denominator > 0) {
          abs(diff(means)) / denominator
        } else if (diff(means) == 0) {
          0
        } else {
          Inf
        }
      } else {
        diagnostics$marginal_joint_mean_ratio <- NA_real_
        diagnostics$marginal_joint_standardized_difference <- NA_real_
      }
      list(
        summary = summary,
        quantiles = unname(q[-(1:3)]),
        mode = unname(mode),
        valid = TRUE,
        reason = NA_character_,
        diagnostics = diagnostics
      )
    },
    error = function(error) invalid_noise_summary(conditionMessage(error), diagnostics, probs)
  )
}

# 4. Check saved summaries and prediction alignment ---------------------------
# These checks flag invalid supports/quantiles and inconsistent moments
# They return flags rather than discarding all outputs from the fit
check_parameter_summaries <- function(parameters) {
  required <- c("parameter", "posterior_mean", "posterior_sd", "lower", "median", "upper")
  if (!is.data.frame(parameters) || !all(required %in% names(parameters))) {
    stop("Parameter quality requires parameter, mean, SD and three quantile columns.")
  }
  out <- parameters
  n <- nrow(out)
  positive <- out$parameter %in% c("tau2_s", "tau2_B", "tau2_shared", "rho", "theta", "sigma2") |
    grepl("^(native:)?(Precision|Range|Stdev) for ", out$parameter)
  positive[is.na(positive)] <- FALSE
  mn <- out$posterior_mean
  sd <- out$posterior_sd
  lo <- out$lower
  md <- out$median
  hi <- out$upper
  out$positive_support <- positive
  out$sd_valid <- is.finite(sd) & sd >= 0
  out$mean_valid <- is.finite(mn) & out$sd_valid & (!positive | mn > 0)
  out$median_valid <- is.finite(md) & (!positive | md > 0)
  out$interval_valid <- is.finite(lo) & is.finite(hi) & lo <= hi &
    (!positive | (lo > 0 & hi > 0))
  quantiles_valid <- out$interval_valid & out$median_valid & lo <= md & md <= hi
  out$summary_valid <- out$mean_valid & quantiles_valid
  out$truth_available <- if ("truth" %in% names(out)) {
    is.finite(out$truth)
  } else {
    rep(FALSE, n)
  }
  reasons <- rep("", n)
  add_issue <- function(condition, label) {
    hit <- which(condition %in% TRUE)
    reasons[hit] <<- ifelse(nzchar(reasons[hit]), paste(reasons[hit], label, sep = "; "), label)
  }
  add_issue(!is.finite(mn), "nonfinite_mean")
  add_issue(!is.finite(sd), "nonfinite_sd")
  add_issue(is.finite(sd) & sd < 0, "negative_sd")
  add_issue(positive & is.finite(mn) & mn <= 0, "nonpositive_mean")
  add_issue(!is.finite(lo) | !is.finite(md) | !is.finite(hi), "nonfinite_quantiles")
  add_issue(is.finite(lo) & is.finite(md) & is.finite(hi) & (lo > md | md > hi),
    "unordered_quantiles")
  add_issue(
    positive & ((is.finite(lo) & lo <= 0) | (is.finite(md) & md <= 0) |
        (is.finite(hi) & hi <= 0)),
    "nonpositive_quantiles"
  )

  # [EXTENSION] Necessary moment/quantile bounds, allowing 1% numerical tolerance
  # Coherence warnings retain means in the primary analysis. A missing truth
  # target is intentional for some native parameters and is not a fit defect
  eligible <- out$summary_valid
  median_warning <- eligible & abs(mn - md) > 1.01 * sd
  lower_warning <- eligible & pmax(mn - lo, 0) > 1.01 * sqrt(39) * sd
  upper_warning <- eligible & pmax(hi - mn, 0) > 1.01 * sqrt(39) * sd
  positive_bound <- .475 * lo + .475 * md + .025 * hi
  positive_warning <- eligible & positive & mn * 1.01 < positive_bound
  add_issue(median_warning, "mean_median_sd_incoherence")
  add_issue(lower_warning, "lower_tail_moment_incoherence")
  add_issue(upper_warning, "upper_tail_moment_incoherence")
  add_issue(positive_warning, "positive_quantile_mean_incoherence")
  out$coherence_warning <- median_warning | lower_warning | upper_warning | positive_warning
  out$issue <- reasons
  out
}

# Apply the same diagnostics to the native fixed-effect/hyperparameter table
# as to the interpretable parameters reported by score_parameters()
check_native_parameter_summaries <- function(fit) {
  rows <- lapply(
    c("summary.fixed", "summary.hyperpar"),
    function(name) {
      x <- fit[[name]]
      if (is.null(x) || !nrow(x)) {
        return(NULL)
      }
      column <- function(name) if (name %in% colnames(x)) {
        as.numeric(x[, name])
      } else
      {
        rep(NA_real_, nrow(x))
      }
      data.frame(
        parameter = paste0("native:", rownames(x)),
        scale = "native_INLA",
        posterior_mean = column("mean"),
        posterior_sd = column("sd"),
        lower = column("0.025quant"),
        median = column("0.5quant"),
        upper = column("0.975quant"),
        truth = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(check_parameter_summaries(data.frame(
      parameter = character(),
      posterior_mean = numeric(),
      posterior_sd = numeric(),
      lower = numeric(),
      median = numeric(),
      upper = numeric()
    )))
  }
  check_parameter_summaries(do.call(rbind, rows))
}

# Compare coordinates as well as row counts. Otherwise, a reordered raster
# can produce plausible-looking scores against the wrong truth locations
check_prediction_alignment <- function(bundle, truth) {
  z <- if (is.list(truth)) {
    truth$Z
  } else {
    truth
  }
  actual <- bundle$prediction_coordinates
  expected <- if (is.list(truth)) {
    truth$coordinates
  } else {
    NULL
  }
  issue <- character()
  numeric_coordinates <- function(x) {
    is.matrix(x) && is.numeric(x) && ncol(x) == 2L && all(is.finite(x))
  }
  coordinates_valid <- numeric_coordinates(actual) && numeric_coordinates(expected) &&
    identical(dim(actual), dim(expected)) && nrow(actual) == length(z) &&
    isTRUE(all.equal(unname(actual), unname(expected), tolerance = 1e-10))
  if (!coordinates_valid) {
    issue <- c(issue, "prediction_coordinates_missing_or_misaligned")
  }
  pred <- bundle$predictions
  geometry_checked <- inherits(pred, "sf")
  if (geometry_checked) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      coordinates_valid <- FALSE
      issue <- c(issue, "sf_unavailable_for_geometry_check")
    } else {
      geometry <- tryCatch(sf::st_coordinates(pred)[, 1:2, drop = FALSE], error = function(e) NULL)
      geometry_valid <- numeric_coordinates(geometry) && numeric_coordinates(actual) &&
        identical(dim(geometry), dim(actual)) &&
        isTRUE(all.equal(unname(geometry), unname(actual), tolerance = 1e-10))
      coordinates_valid <- coordinates_valid && geometry_valid
      if (!geometry_valid) {
        issue <- c(issue, "prediction_geometry_misaligned")
      }
    }
  }
  values <- if (inherits(pred, "sf") && requireNamespace("sf", quietly = TRUE)) {
    sf::st_drop_geometry(pred)
  } else {
    pred
  }
  truth_valid <- is.numeric(z) && length(z) > 0L && all(is.finite(z))
  rows_valid <- is.data.frame(values) && nrow(values) == length(z)
  valid_column <- function(name, nonnegative = FALSE) {
    x <- values[[name]]
    rows_valid && is.numeric(x) && length(x) == length(z) && all(is.finite(x)) &&
      (!nonnegative || all(x >= 0))
  }
  mean_valid <- valid_column("pred_mean")
  median_valid <- valid_column("pred_med")
  sd_valid <- valid_column("sd", nonnegative = TRUE)
  summaries_valid <- mean_valid && median_valid && sd_valid
  if (!truth_valid) {
    issue <- c(issue, "invalid_prediction_truth")
  }
  if (!rows_valid) {
    issue <- c(issue, "prediction_row_count_misaligned")
  }
  if (!mean_valid) {
    issue <- c(issue, "invalid_prediction_mean")
  }
  if (!median_valid) {
    issue <- c(issue, "invalid_prediction_median")
  }
  if (!sd_valid) {
    issue <- c(issue, "invalid_prediction_sd")
  }
  list(
    valid = coordinates_valid && truth_valid && summaries_valid,
    coordinates_valid = coordinates_valid,
    truth_valid = truth_valid,
    rows_valid = rows_valid,
    mean_valid = mean_valid,
    median_valid = median_valid,
    sd_valid = sd_valid,
    summaries_valid = summaries_valid,
    geometry_checked = geometry_checked,
    issue = paste(unique(issue), collapse = "; ")
  )
}

# 5. Basic scoring formulas and posterior-draw seeds --------------------------
# Normalize the posterior-sampling seed to the supported integer range
posterior_draw_seed <- function(seed) {
  if (length(seed) != 1L || !is.finite(seed) || seed < 1 || seed != floor(seed)) {
    stop("seed must be a positive integer.")
  }
  as.integer((as.double(seed) - 1) %% 2147483646 + 1)
}

# E|X| for X ~ Normal(mean, sd^2), including the sd=0 point-mass limit
# Used in the closed-form Gaussian approximation to latent-surface CRPS
expected_absolute_normal <- function(mean, sd) {
  if (any(!is.finite(mean)) || any(!is.finite(sd)) || any(sd < 0)) {
    stop("Finite normal means and nonnegative standard deviations are required.")
  }
  n <- max(length(mean), length(sd))
  if (!length(mean) || !length(sd) ||
      !all(c(length(mean), length(sd)) %in% c(1L, n))) {
    stop("Incompatible lengths.")
  }
  mean <- rep_len(mean, n)
  sd <- rep_len(sd, n)
  ans <- abs(mean)
  ok <- sd > 0
  z <- mean[ok] / sd[ok]
  ans[ok] <- 2 * sd[ok] * stats::dnorm(z) + mean[ok] * (2 * stats::pnorm(z) - 1)
  ans
}

# Gaussian CRPS = E|X-y| - 0.5 E|X-Xprime|, evaluated in closed form
gaussian_crps <- function(mean, sd, observed) {
  expected_absolute_normal(mean - observed, sd) - sd / sqrt(pi)
}

# Central (1-alpha) interval width plus a penalty when truth lies outside it
interval_score <- function(lower, upper, truth, alpha = .05) {
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie in (0,1).")
  }
  if (any(lower > upper, na.rm = TRUE)) {
    stop("Interval lower bound exceeds upper bound.")
  }
  upper - lower + (2 / alpha) * pmax(lower - truth, 0) +
    (2 / alpha) * pmax(truth - upper, 0)
}

# 6. Match parameter summaries to their observation likelihoods ---------------
# Use the saved mapping, never hyperparameter position: shared and separate
# noise models return different numbers of Gaussian precisions
get_gaussian_likelihoods <- function(bundle) {
  map <- bundle$likelihood_map
  required <- c("component", "family", "column", "precision_name", "fixed", "fixed_variance")
  if (!is.data.frame(map) || !all(required %in% names(map))) {
    stop("Saved bundle must include the complete likelihood_map.")
  }
  if (anyNA(map$component) || anyDuplicated(map$component)) {
    stop("Likelihood components must be distinct and nonmissing.")
  }
  if (sum(map$component == "point") != 1L) {
    stop("The likelihood_map must identify exactly one point likelihood.")
  }
  noise <- map[map$component %in% c("point", "area"), , drop = FALSE]
  if (anyNA(noise$family) || any(noise$family != "gaussian")) {
    stop("Point and area likelihoods must be Gaussian.")
  }
  if (!is.logical(noise$fixed) || anyNA(noise$fixed)) {
    stop("Likelihood fixed flags must be nonmissing logical values.")
  }
  if (!is.numeric(noise$column) || any(!is.finite(noise$column)) ||
      any(noise$column < 1) || any(noise$column != floor(noise$column))) {
    stop("Likelihood columns must be positive integers.")
  }
  estimated <- !noise$fixed
  if (anyNA(noise$precision_name[estimated]) ||
      any(!nzchar(noise$precision_name[estimated]))) {
    stop("Each estimated Gaussian likelihood needs its saved precision name.")
  }
  if (any(noise$fixed) &&
      (any(!is.finite(noise$fixed_variance[noise$fixed])) ||
        any(noise$fixed_variance[noise$fixed] <= 0))) {
    stop("Fixed Gaussian variances must be positive and finite.")
  }
  if (!identical(bundle$model, "PSgeo") && !"area" %in% noise$component) {
    stop("A melding fit must identify its area likelihood.")
  }
  if (identical(bundle$model, "PSgeo") && "area" %in% noise$component) {
    stop("PSgeo must not contain an area likelihood.")
  }
  noise[match(c("point", "area"), noise$component, nomatch = 0L), , drop = FALSE]
}

# 7. Score interpretable model parameters -------------------------------------
# Resolve truths without assigning a false target to the conditioned sampling
# intercept or to a constant gamma in a spatially heterogeneous mechanism
resolve_scoring_truth <- function(truth, config) {
  z <- if (is.list(truth)) {
    truth$Z
  } else {
    truth
  }
  if (!is.numeric(z) || !length(z) || any(!is.finite(z))) {
    stop("truth must supply a finite numeric latent surface Z.")
  }
  values <- config
  if (is.list(truth) && !is.null(truth$parameters)) {
    values <- utils::modifyList(values, as.list(truth$parameters))
  }
  scalar <- function(nm, default = NA_real_) {
    x <- values[[nm]]
    if (is.null(x)) {
      default
    } else if (length(x) == 1L && is.numeric(x)) {
      as.numeric(x)
    } else
    {
      stop("Parameter truth '", nm, "' must be a scalar.")
    }
  }
  pars <- setNames(
    vapply(
      c("mu", "rho", "sigma2", "theta", "gamma", "alpha", "tau2_s", "tau2_B"),
      scalar,
      numeric(1)
    ),
    c("mu", "rho", "sigma2", "theta", "gamma", "alpha", "tau2_s", "tau2_B")
  )
  if (!is.finite(pars[["theta"]]) && is.finite(pars[["rho"]])) {
    pars[["theta"]] <- 8 * scalar("sigma2", 1) / pars[["rho"]]^2
  }
  if (!is.finite(pars[["tau2_s"]]) || pars[["tau2_s"]] <= 0) {
    stop("Supply a positive true point-noise variance tau2_s.")
  }

  # [EXTENSION] A fixed point count cancels the intensity intercept
  # There is no identifiable true b0pp = .05 for this sampling design
  if (!isTRUE(config$score_conditioned_alpha)) {
    pars[["alpha"]] <- NA_real_
  }
  list(Z = z, parameters = pars)
}

# Transform the entire marginal, not just its mean. Reverse quantiles for
# decreasing transformations so lower/median/upper keep their meaning
summarize_transformed_marginal <- function(
  marginal, transform = identity, decreasing = FALSE
) {
  if (is.null(marginal) || !nrow(marginal)) {
    stop("Required INLA marginal was not saved.")
  }
  qs <- INLA::inla.qmarginal(c(.025, .5, .975), marginal)
  qs <- transform(qs)
  if (decreasing) {
    qs <- rev(qs)
  }
  mn <- INLA::inla.emarginal(transform, marginal)
  second <- INLA::inla.emarginal(function(x) transform(x)^2, marginal)
  c(
    posterior_mean = mn,
    posterior_sd = sqrt(max(0, second - mn^2)),
    lower = qs[1L],
    median = qs[2L],
    upper = qs[3L]
  )
}

# Report native summaries plus mu, gamma, range, microergodic parameter, and
# measurement-error variances. Missing truths remain NA for coverage scoring
score_parameters <- function(bundle, truth, config = list()) {
  fit <- bundle$fit
  p <- resolve_scoring_truth(truth, config)$parameters
  rows <- list()
  noise_diagnostics <- list()
  add <- function(
    name, value, actual = NA_real_, scale = "derived", source = name, note = ""
  ) {
    value <- tryCatch(
      value,
      error = function(error) {
        out <- setNames(rep(NA_real_, 5L), c("posterior_mean", "posterior_sd", "lower", "median", "upper"))
        attr(out, "summary_error") <- conditionMessage(error)
        out
      }
    )
    if (!is.null(attr(value, "summary_error"))) {
      note <- paste(note, "Summary failed:", attr(value, "summary_error"))
    }
    r <- data.frame(
      parameter = name,
      scale = scale,
      source = source,
      posterior_mean = unname(value[["posterior_mean"]]),
      posterior_sd = unname(value[["posterior_sd"]]),
      lower = unname(value[["lower"]]),
      median = unname(value[["median"]]),
      upper = unname(value[["upper"]]),
      truth = actual,
      note = note,
      stringsAsFactors = FALSE
    )
    r$mean_error <- r$posterior_mean - actual
    r$median_error <- r$median - actual
    interval_ok <- check_parameter_summaries(r)$interval_valid
    r$covered <- if (is.finite(actual) && interval_ok) {
      r$lower <= actual && actual <= r$upper
    } else {
      NA
    }

    # [EXTENSION] Score against the true parameter using the paper's unweighted IS
    r$IS_paper <- if (is.finite(actual) && interval_ok) {
      interval_score(r$lower, r$upper, actual)
    } else {
      NA_real_
    }

    rows[[length(rows) + 1L]] <<- r
  }
  summary_value <- function(row) {
    c(
      posterior_mean = unname(row[["mean"]]),
      posterior_sd = unname(row[["sd"]]),
      lower = unname(row[["0.025quant"]]),
      median = unname(row[["0.5quant"]]),
      upper = unname(row[["0.975quant"]])
    )
  }
  field <- if (!is.null(bundle$spde_field_name)) {
    bundle$spde_field_name
  } else
  if (identical(bundle$model, "Melding")) {
    "s"
  } else {
    "i"
  }
  theta1 <- paste("Theta1 for", field)
  theta2 <- paste("Theta2 for", field)
  native_truth <- c(
    b0 = p[["mu"]],
    b0pp = p[["alpha"]],
    "Beta for j" = p[["gamma"]]
  )
  native_truth[theta1] <- if (is.finite(p[["theta"]]) && p[["theta"]] > 0) {
    -.5 * log(4 * pi * p[["theta"]])
  } else {
    NA_real_
  }
  native_truth[theta2] <- if (is.finite(p[["rho"]]) && p[["rho"]] > 0) {
    log(sqrt(8) / p[["rho"]])
  } else {
    NA_real_
  }
  noise_map <- get_gaussian_likelihoods(bundle)
  noise_truth <- c(point = p[["tau2_s"]], area = p[["tau2_B"]])
  precision_names <- unique(noise_map$precision_name[!noise_map$fixed])
  for (name in precision_names) {
    actual <- noise_truth[noise_map$component[!noise_map$fixed & noise_map$precision_name == name]]
    # Explicitly shared precision has no single true value if the noises differ
    native_truth[name] <- if (all(is.finite(actual)) && length(unique(actual)) == 1L) {
      1 / actual[1L]
    } else {
      NA_real_
    }
  }
  for (block in c("summary.fixed", "summary.hyperpar")) {
    tab <- fit[[block]]
    if (is.null(tab)) {
      next
    }
    for (nm in rownames(tab)) {
      actual <- native_truth[nm]
      if (!length(actual)) {
        actual <- NA_real_
      }
      note <- if (nm == "b0pp" && !isTRUE(config$score_conditioned_alpha)) {
        "No true intensity intercept under fixed-count sampling; not scored against alpha=.05."
      } else {
        ""
      }

      add(
        paste0("native:", nm),
        summary_value(tab[nm, , drop = FALSE]),
        unname(actual),
        "native_INLA",
        nm,
        note
      )
    }
  }
  aliases <- c(mu = "b0", alpha = "b0pp", gamma = "Beta for j")
  for (name in names(aliases)) {
    alias <- aliases[[name]]
    tab <- if (alias %in% rownames(fit$summary.fixed)) {
      fit$summary.fixed
    } else {
      fit$summary.hyperpar
    }
    if (!is.null(tab) && alias %in% rownames(tab)) {
      add(
        name,
        summary_value(tab[alias, , drop = FALSE]),
        p[[name]],
        "interpretable",
        alias,
        if (name == "alpha" && !isTRUE(config$score_conditioned_alpha)) {
          "Conditioned point count: alpha has no identifiable true value."
        } else {
          ""
        }
      )
    }
  }
  marginals <- fit$marginals.hyperpar
  validate_spde_transform <- function() {
    if (!is.null(config$nu) && config$nu != 1) {
      stop("Microergodic/range transformation requires the repository nu=1 SPDE.")
    }
    spde <- bundle[["spde", exact = TRUE]]
    B0 <- spde$param.inla$B0
    B1 <- spde$param.inla$B1
    if (is.null(B0) || is.null(B1) || ncol(B0) != 3L || ncol(B1) != 3L ||
        any(abs(sweep(B0, 2L, c(0, 1, 0), "-")) > 1e-12) ||
        any(abs(sweep(B1, 2L, c(0, 0, 2), "-")) > 1e-12)) {
      stop("Saved SPDE does not use the original log(tau), log(kappa) parameterization.")
    }
    invisible(TRUE)
  }
  if (theta1 %in% rownames(fit$summary.hyperpar) || !is.null(marginals[[theta1]])) {
    # [ORIGINAL] The default SPDE has alpha = 2, d = 2, and Theta1 = log(tau)
    # [EXTENSION] Compute E[theta] from its marginal, not from E[Theta1]
    add(
      "theta",
      {
        validate_spde_transform()
        summarize_transformed_marginal(
          marginals[[theta1]],
          function(x) exp(-2 * x) / (4 * pi),
          decreasing = TRUE
        )
      },
      p[["theta"]],
      source = theta1,
      note = "nu=1,d=2: theta=1/(4*pi*exp(2*Theta1)); marginal transformed."
    )
  }
  if (theta2 %in% rownames(fit$summary.hyperpar) || !is.null(marginals[[theta2]])) {
    add(
      "rho",
      {
        validate_spde_transform()
        summarize_transformed_marginal(
          marginals[[theta2]],
          function(x) sqrt(8) * exp(-x),
          decreasing = TRUE
        )
      },
      p[["rho"]],
      source = theta2,
      note = "nu=1: rho=sqrt(8)*exp(-Theta2); marginal transformed."
    )
  }
  # INLA's log-variance marginal is a linear combination of the joint SPDE
  # hyperparameters, retaining dependence between log(tau) and log(kappa)
  add(
    "sigma2",
    {
      validate_spde_transform()
      spatial <- INLA::inla.spde2.result(fit, field, bundle$spde, do.transform = FALSE)
      summarize_transformed_marginal(spatial$marginals.log.variance.nominal[[1L]], exp)
    },
    p[["sigma2"]],
    source = "SPDE log nominal variance marginal",
    note = "Transform the joint-derived log-variance marginal, not posterior means; nominal covariance variance."
  )
  for (j in seq_len(nrow(noise_map))) {
    component <- noise_map$component[j]
    parameter <- if (component == "point") {
      "tau2_s"
    } else {
      "tau2_B"
    }
    actual <- noise_truth[[component]]
    if (noise_map$fixed[j]) {
      variance <- noise_map$fixed_variance[j]
      add(
        parameter,
        c(
          posterior_mean = variance,
          posterior_sd = 0,
          lower = variance,
          median = variance,
          upper = variance
        ),
        actual,
        source = paste(component, "fixed Gaussian variance"),
        note = "Fixed in fitting; degenerate interval."
      )
    } else {
      name <- noise_map$precision_name[j]
      components <- noise_map$component[!noise_map$fixed & noise_map$precision_name == name]
      if (length(unique(noise_truth[components])) > 1L) {
        actual <- NA_real_
      }
      noise <- summarize_noise_variance(fit, name)
      noise_diagnostics[[component]] <- noise
      add(
        parameter,
        noise$summary,
        actual,
        source = name,
        note = paste(
          component,
          "likelihood: E[1/precision], not 1/E[precision]; mapping saved in likelihood_map."
        )
      )
    }
  }
  ans <- do.call(rbind, rows)
  rownames(ans) <- NULL
  ans <- check_parameter_summaries(ans)
  attr(ans, "noise_diagnostics") <- noise_diagnostics
  ans
}

# 8. Draw the latent response surface from the joint posterior ----------------
# Rows are prediction locations; columns are joint posterior draws. The saved
# stack indices locate Z in APredictor, which already includes the projector A
draw_latent_predictions <- function(bundle, n = 100L, seed = 1L) {
  if (length(n) != 1L || !is.finite(n) || n < 2 || n != floor(n)) {
    stop("At least two joint draws are required.")
  }
  fit <- bundle$fit
  if (is.null(fit$misc$configs)) {
    stop("Saved fit lacks control.compute$config=TRUE configurations.")
  }
  idx <- as.integer(bundle$prediction_index)
  if (!length(idx) || anyNA(idx) || any(idx < 1L) || anyDuplicated(idx)) {
    stop("Invalid saved prediction_index.")
  }
  contents <- fit$misc$configs$contents

  # [EXTENSION] Stack data rows index APredictor after applying the A matrix
  # Predictor indexes the smaller unprojected latent linear predictor
  tag <- if ("APredictor" %in% contents$tag) {
    "APredictor"
  } else {
    "Predictor"
  }
  count <- contents$length[match(tag, contents$tag)]
  if (length(count) != 1L || !is.finite(count) || any(idx > count)) {
    stop("Prediction indices do not align with saved ", tag, " configurations.")
  }
  selection <- setNames(list(idx), tag)
  samples <- with_simulation_seed(
    posterior_draw_seed(seed),
    INLA::inla.posterior.sample(
      n = as.integer(n),
      result = fit,
      selection = selection,
      intern = TRUE,
      add.names = TRUE,
      seed = posterior_draw_seed(seed),
      num.threads = "1:1",
      parallel.configs = FALSE
    )
  )
  wanted <- paste0(tag, ":", idx)
  Z <- vapply(
    samples,
    function(s) {
      index <- match(wanted, rownames(s$latent))
      if (anyNA(index)) {
        stop("A joint draw is missing required named predictor rows.")
      }
      as.numeric(s$latent[index, 1L])
    },
    numeric(length(idx))
  )
  Z <- matrix(Z, nrow = length(idx), ncol = as.integer(n))
  if (any(!is.finite(Z))) {
    stop("Nonfinite joint latent draws.")
  }
  list(
    Z = Z,
    metadata = list(
      seed = posterior_draw_seed(seed),
      n = as.integer(n),
      predictor_tag = tag,
      prediction_index = idx,
      sampling = "Latent response surface from INLA joint configurations; no observation-noise draws.",
      num_threads = "1:1",
      intern = TRUE
    )
  )
}

# 9. Compute CRPS from draws --------------------------------------------------
# Estimate E|X-y| - 0.5 E|X-Xprime| using distinct posterior-draw pairs
# Here y is the known latent truth Z(s), with no added measurement noise
# This unbiased finite-draw estimator can be negative
sample_crps <- function(draws, observed) {
  # sum_{i < j}|x_i - x_j| = sum_i(2i - M - 1)x_(i)
  # Sorting avoids an M x M distance matrix at every spatial location
  draws <- as.matrix(draws)
  n <- ncol(draws)
  if (n < 2L || nrow(draws) != length(observed)) {
    stop("Invalid sample-CRPS inputs.")
  }
  sorted <- t(apply(draws, 1L, sort))
  half_pair <- as.numeric(sorted %*% (2 * seq_len(n) - n - 1)) / (n * (n - 1))
  rowMeans(abs(draws - observed)) - half_pair
}
