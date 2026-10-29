# General utilities -----------------------------------------------------------


safe_mean = function(x) {
  x = x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}


safe_median = function(x) {
  x = x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}


safe_mcse = function(x) {
  x = x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}


bind_rows = function(x) {
  x = Filter(function(d) is.data.frame(d) && nrow(d) > 0L, x)
  if (!length(x)) return(data.frame())
  output = do.call(rbind, x)
  rownames(output) = NULL
  output
}


group_summary_rows = function(data, group_vars) {
  stopifnot(all(group_vars %in% names(data)))
  grouping = do.call(
    interaction, c(data[group_vars], list(drop = TRUE, lex.order = TRUE))
  )
  split(data, grouping, drop = TRUE)
}


interval_score = function(lower, upper, truth, alpha = 0.05) {
  (upper - lower) +
    (2 / alpha) * pmax(lower - truth, 0) +
    (2 / alpha) * pmax(truth - upper, 0)
}


# Prerun extensions -----------------------------------------------------------


summarise_fit = function(data, logged, model, mesh, A_B, replication,
                         tau2_B, mesh_name, areal_method, prior_name) {

  log_row = logged$log[logged$log$model == model, , drop = FALSE]
  fitted = logged$results$models[[model]]

  # Initialize diagnostics
  row = data.frame(
    replication = replication,
    seed = 123L + replication - 1L,
    tau2_B_true = tau2_B,
    mesh = mesh_name,
    mesh_vertices = mesh$n,
    areal_method = areal_method,
    areal_prior = prior_name,
    model = model,
    status = log_row$status,
    included = log_row$included,
    n_warnings = log_row$n_warnings,
    elapsed_sec = log_row$elapsed_sec,
    error_message = log_row$error_message,
    point_noise_realized = mean((data$points$value - data$points$true)^2),
    areal_noise_realized = mean((data$areas$value - data$areas$true)^2),
    n_gaussian_precisions = NA_integer_,
    tau2_s_median = NA_real_,
    tau2_B_median = NA_real_,
    tau2_B_q025 = NA_real_,
    tau2_B_q975 = NA_real_,
    tau2_B_covered = NA,
    near_zero_tau2_B = NA,
    gamma_median = NA_real_,
    prediction_mse = NA_real_,
    prediction_coverage = NA_real_,
    areal_truth_rmse = NA_real_,
    hyperparameter_names = NA_character_
  )

  if (is.null(fitted)) return(row)

  fit = fitted$fit
  hp = fit$summary.hyperpar

  row$hyperparameter_names = paste(rownames(hp), collapse = " | ")

  # Extract observation-error variances
  gaussian_rows = grep(
    "^Precision for the Gaussian",
    rownames(hp),
    ignore.case = TRUE
  )

  row$n_gaussian_precisions = length(gaussian_rows)

  if (length(gaussian_rows) == 2L &&
      all(c("0.025quant", "0.5quant", "0.975quant") %in% colnames(hp))) {

    point_precision = hp[gaussian_rows[1], "0.5quant"]
    area_precision = hp[gaussian_rows[2], ]

    row$tau2_s_median = 1 / point_precision
    row$tau2_B_median = 1 / area_precision[["0.5quant"]]

    # Reciprocal transformation reverses quantile order
    row$tau2_B_q025 = 1 / area_precision[["0.975quant"]]
    row$tau2_B_q975 = 1 / area_precision[["0.025quant"]]

    row$tau2_B_covered = row$tau2_B_q025 <= tau2_B &&
      tau2_B <= row$tau2_B_q975

    row$near_zero_tau2_B = row$tau2_B_median < 0.001
  }

  # Extract preferentiality estimate
  gamma_row = which(rownames(hp) == "Beta for j")

  if (length(gamma_row) == 1L) {
    row$gamma_median = hp[gamma_row, "0.5quant"]
  }

  # Evaluate spatial predictions
  pred = fitted$prediction

  if (!is.null(pred) && nrow(pred) == nrow(data$prediction)) {
    row$prediction_mse = mean(
      (pred$mean - data$prediction$true)^2
    )

    row$prediction_coverage = mean(
      pred$q025 <= data$prediction$true &
        data$prediction$true <= pred$q975
    )
  }

  # Evaluate estimated areal signal
  spatial_mean = fit$summary.random$i$mean

  if (!is.null(spatial_mean) &&
      length(spatial_mean) == ncol(A_B) &&
      "b0" %in% rownames(fit$summary.fixed)) {

    areal_estimate = fit$summary.fixed["b0", "mean"] +
      as.numeric(A_B %*% spatial_mean)

    row$areal_truth_rmse = sqrt(
      mean((areal_estimate - data$areas$true)^2)
    )
  }

  row
}


summarize_fit_log = function(fit_log, models = unique(fit_log$model)) {
  comparison = lapply(models, function(name) {
    x = fit_log[fit_log$model == name, , drop = FALSE]

    # Reused fits are not additional fitting attempts
    reused = if ("reused" %in% names(x)) {
      !is.na(x$reused) & x$reused
    } else rep(FALSE, nrow(x))

    attempted = x[!reused, , drop = FALSE]

    data.frame(
      model = name,
      attempted = nrow(attempted),
      included = sum(attempted$included),
      warnings = sum(attempted$status == "warning"),
      invalid = sum(attempted$status == "invalid"),
      errors = sum(attempted$status == "error"),
      reused = sum(reused),
      mean_elapsed_sec = if (nrow(attempted)) {
        mean(attempted$elapsed_sec)
      } else NA_real_
    )
  })

  comparison = do.call(rbind, comparison)
  rownames(comparison) = NULL

  comparison
}


summarize_marginal = function(marginal) {
  if (is.null(marginal)) stop("Posterior marginal is unavailable.")

  posterior_mean = INLA::inla.emarginal(function(x) x, marginal)
  posterior_second = INLA::inla.emarginal(function(x) x^2, marginal)
  quantiles = INLA::inla.qmarginal(c(0.025, 0.975), marginal)

  c(
    mean = posterior_mean,
    sd = sqrt(max(0, posterior_second - posterior_mean^2)),
    q025 = quantiles[1],
    q975 = quantiles[2]
  )
}


evaluate = function(pred, truth) {
  c(
    rmse = sqrt(mean((pred$mean - truth)^2)),
    coverage = mean(truth >= pred$q025 & truth <= pred$q975),
    mean_sd = mean(pred$sd)
  )
}


compare_predictions = function(pred1, pred2) {
  stopifnot(isTRUE(all.equal(pred1[, c("x", "y")], pred2[, c("x", "y")])))

  difference = pred1$mean - pred2$mean

  c(
    mean_difference = mean(difference),
    mean_abs_difference = mean(abs(difference)),
    max_abs_difference = max(abs(difference)),
    correlation = cor(pred1$mean, pred2$mean)
  )
}


evaluate_parameter_replication = function(results, truth, separate_noise = TRUE) {
  comparison = lapply(names(results$models), function(name) {
    fit = results$models[[name]]$fit
    hyperpar = fit$marginals.hyperpar

    if (is.null(hyperpar)) {
      stop("Posterior hyperparameter marginals are missing for ", name, ".")
    }

    # Collect posterior marginals on the parameter scales of interest
    marginals = list()

    if ("Beta for j" %in% names(hyperpar)) {
      marginals$gamma = hyperpar[["Beta for j"]]
    }

    # Transform Gaussian observation precisions into variances
    precision_names = grep(
      "^Precision for .*Gaussian",
      names(hyperpar),
      value = TRUE,
      ignore.case = TRUE
    )

    if (length(precision_names) < 1L) {
      stop("Gaussian observation precision is missing for ", name, ".")
    }

    marginals$tau2_s = INLA::inla.tmarginal(
      function(x) 1 / x,
      hyperpar[[precision_names[1]]]
    )

    if (separate_noise && name %in% c("melding", "ps_melding")) {
      if (length(precision_names) < 2L) {
        stop("Areal observation precision is missing for ", name, ".")
      }

      marginals$tau2_B = INLA::inla.tmarginal(
        function(x) 1 / x,
        hyperpar[[precision_names[2]]]
      )
    }

    # Transform SPDE parameters into nominal range and variance
    if (!is.null(results$spde)) {
      spde_result = INLA::inla.spde2.result(
        inla = fit,
        name = "i",
        spde = results$spde,
        do.transform = TRUE
      )

      marginals$rho = spde_result$marginals.range.nominal[[1]]
      marginals$sigma2 = spde_result$marginals.variance.nominal[[1]]
    }

    # Only evaluate parameters that exist in the fitted model
    parameters = intersect(names(marginals), names(truth))
    parameters = parameters[is.finite(truth[parameters])]

    do.call(rbind, lapply(parameters, function(parameter) {
      posterior = summarize_marginal(marginals[[parameter]])
      true = as.numeric(truth[[parameter]])

      data.frame(
        model = name,
        parameter = parameter,
        true = true,
        estimate = posterior[["mean"]],
        posterior_sd = posterior[["sd"]],
        q025 = posterior[["q025"]],
        q975 = posterior[["q975"]],
        bias = posterior[["mean"]] - true,
        squared_error = (posterior[["mean"]] - true)^2,
        coverage = as.numeric(true >= posterior[["q025"]] &
                                true <= posterior[["q975"]]),
        interval_score = interval_score(
          posterior[["q025"]], posterior[["q975"]], true
        )
      )
    }))
  })

  comparison = do.call(rbind, comparison)

  if (!is.null(comparison)) {
    rownames(comparison) = NULL
  }

  comparison
}


summarize_ps = function(results, truth) {
  comparison = lapply(results, function(result) {
    fit = result$models$ps_melding$fit
    pred = result$models$ps_melding$prediction

    gamma = fit$summary.hyperpar["Beta for j", ]

    c(
      gamma = gamma[["mean"]],
      gamma_sd = gamma[["sd"]],
      gamma_q025 = gamma[["0.025quant"]],
      gamma_q975 = gamma[["0.975quant"]],
      evaluate(pred, truth)
    )
  })

  do.call(rbind, comparison)
}


summarize_replications = function(metrics, gamma_true = 1) {
  # Monte Carlo standard error across independent replications
  mcse = function(x) {
    if (length(x) < 2L) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  }

  # Average predictive performance
  metric_names = c("mse", "mae", "wd2", "bias", "coverage", "interval_score")

  prediction_n = stats::aggregate(
    metrics[, "mse", drop = FALSE],
    by = list(model = metrics$model),
    FUN = length
  )

  names(prediction_n)[2] = "n_replications"

  prediction_mean = stats::aggregate(
    metrics[, metric_names, drop = FALSE],
    by = list(model = metrics$model),
    FUN = mean
  )

  # Between-replication standard deviations
  prediction_sd = stats::aggregate(
    metrics[, metric_names, drop = FALSE],
    by = list(model = metrics$model),
    FUN = sd
  )

  # Monte Carlo standard errors of average predictive performance
  prediction_mcse = stats::aggregate(
    metrics[, metric_names, drop = FALSE],
    by = list(model = metrics$model),
    FUN = mcse
  )

  # Zhong's original Wasserstein calculation for reproduction
  prediction_reproduction = stats::aggregate(
    metrics[, "wd_zhong", drop = FALSE],
    by = list(model = metrics$model),
    FUN = mean
  )

  prediction_reproduction_mcse = stats::aggregate(
    metrics[, "wd_zhong", drop = FALSE],
    by = list(model = metrics$model),
    FUN = mcse
  )

  # Additional uncertainty diagnostics
  prediction_diagnostics = stats::aggregate(
    metrics[, c("mean_interval_width", "mean_sd"), drop = FALSE],
    by = list(model = metrics$model),
    FUN = mean
  )

  prediction_diagnostics_mcse = stats::aggregate(
    metrics[, c("mean_interval_width", "mean_sd"), drop = FALSE],
    by = list(model = metrics$model),
    FUN = mcse
  )

  # Recovery of gamma for preferential sampling models
  ps_metrics = metrics[!is.na(metrics$gamma_mean), , drop = FALSE]
  gamma_summary = NULL
  gamma_mcse = NULL

  # A constant gamma has no true value under spatially varying gamma(s)
  if (is.finite(gamma_true) && nrow(ps_metrics) > 0L) {
    gamma_summary = do.call(rbind, lapply(
      split(ps_metrics, ps_metrics$model),
      function(x) {
        data.frame(
          model = x$model[1],
          gamma_mean = mean(x$gamma_mean),
          gamma_bias = mean(x$gamma_mean - gamma_true),
          gamma_rmse = sqrt(mean((x$gamma_mean - gamma_true)^2)),
          gamma_coverage = mean(x$gamma_coverage),
          gamma_interval_score = mean(x$gamma_interval_score),
          mean_gamma_sd = mean(x$gamma_sd)
        )
      }
    ))

    gamma_mcse = do.call(rbind, lapply(
      split(ps_metrics, ps_metrics$model),
      function(x) {
        gamma_sq_error = (x$gamma_mean - gamma_true)^2
        gamma_rmse = sqrt(mean(gamma_sq_error))

        # Delta-method MCSE for gamma RMSE
        gamma_rmse_mcse = if (gamma_rmse > 0) {
          mcse(gamma_sq_error) / (2 * gamma_rmse)
        } else {
          NA_real_
        }

        data.frame(
          model = x$model[1],
          gamma_mean = mcse(x$gamma_mean),
          gamma_bias = mcse(x$gamma_mean - gamma_true),
          gamma_rmse = gamma_rmse_mcse,
          gamma_coverage = mcse(x$gamma_coverage),
          gamma_interval_score = mcse(x$gamma_interval_score),
          mean_gamma_sd = mcse(x$gamma_sd)
        )
      }
    ))
  }

  if (!is.null(gamma_summary)) {
    rownames(gamma_summary) = NULL
  }

  if (!is.null(gamma_mcse)) {
    rownames(gamma_mcse) = NULL
  }

  list(
    prediction_n = prediction_n,
    prediction_mean = prediction_mean,
    prediction_sd = prediction_sd,
    prediction_mcse = prediction_mcse,
    prediction_reproduction = prediction_reproduction,
    prediction_reproduction_mcse = prediction_reproduction_mcse,
    prediction_diagnostics = prediction_diagnostics,
    prediction_diagnostics_mcse = prediction_diagnostics_mcse,
    gamma = gamma_summary,
    gamma_mcse = gamma_mcse
  )
}


summarize_parameter_recovery = function(metrics) {
  if (is.null(metrics) || nrow(metrics) == 0L) {
    return(list(parameter_mean = NULL, parameter_mcse = NULL))
  }

  mcse = function(x) {
    if (length(x) < 2L) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  }

  groups = split(
    metrics,
    list(metrics$model, metrics$parameter),
    drop = TRUE
  )

  parameter_mean = do.call(rbind, lapply(groups, function(x) {
    data.frame(
      model = x$model[1],
      parameter = x$parameter[1],
      true = x$true[1],
      n_replications = nrow(x),
      estimate = mean(x$estimate),
      bias = mean(x$bias),
      rmse = sqrt(mean(x$squared_error)),
      coverage = mean(x$coverage),
      interval_score = mean(x$interval_score),
      mean_posterior_sd = mean(x$posterior_sd)
    )
  }))

  parameter_mcse = do.call(rbind, lapply(groups, function(x) {
    rmse = sqrt(mean(x$squared_error))

    data.frame(
      model = x$model[1],
      parameter = x$parameter[1],
      estimate = mcse(x$estimate),
      bias = mcse(x$bias),
      rmse = if (rmse > 0) mcse(x$squared_error) / (2 * rmse) else NA_real_,
      coverage = mcse(x$coverage),
      interval_score = mcse(x$interval_score),
      mean_posterior_sd = mcse(x$posterior_sd)
    )
  }))

  rownames(parameter_mean) = NULL
  rownames(parameter_mcse) = NULL

  list(
    parameter_mean = parameter_mean,
    parameter_mcse = parameter_mcse
  )
}


# Summarize results -----------------------------------------------------------


extract_parameters = function(model, separate_noise = TRUE) {
  output = c(
    gamma_median = NA_real_,
    tau2_median = NA_real_,
    tau2_s_median = NA_real_,
    tau2_B_median = NA_real_,
    tau2_B_q025 = NA_real_,
    tau2_B_q975 = NA_real_
  )

  hp = model$fit$summary.hyperpar
  needed = c("0.025quant", "0.5quant", "0.975quant")
  if (is.null(hp) || !all(needed %in% colnames(hp))) return(output)

  if ("Beta for j" %in% rownames(hp)) {
    output["gamma_median"] = hp["Beta for j", "0.5quant"]
  }

  precision_rows = grep(
    "^Precision for .*Gaussian", rownames(hp),
    value = TRUE, ignore.case = TRUE
  )

  if (length(precision_rows) >= 1L) {
    variance_median = 1 / hp[precision_rows[1], "0.5quant"]
    output["tau2_median"] = variance_median
    output["tau2_s_median"] = variance_median
  }

  if (separate_noise && length(precision_rows) >= 2L) {
    output["tau2_B_median"] = 1 / hp[precision_rows[2], "0.5quant"]
    # Inverting a precision reverses the quantile endpoints
    output["tau2_B_q025"] = 1 / hp[precision_rows[2], "0.975quant"]
    output["tau2_B_q975"] = 1 / hp[precision_rows[2], "0.025quant"]
  }

  output
}


evaluate_replication = function(results, truth, tau2_s = 0.1,
                                gamma_true = 1, mu_true = 0) {
  comparison = lapply(names(results$models), function(name) {
    model = results$models[[name]]
    pred = model$prediction

    if (nrow(pred) != length(truth)) {
      stop("Prediction and truth lengths do not match.")
    }

    if (is.null(pred$median)) {
      stop("Posterior median is missing. Refit using the updated fit_models.R.")
    }

    # Predictive accuracy and uncertainty
    error = pred$mean - truth

    metrics = data.frame(
      model = name,
      mse = mean(error^2),
      mae = mean(abs(pred$median - truth)),
      bias = mean(error),
      # Match Zhong's released evaluation code
      wd_zhong = mean(error^2 + pred$sd^2 + tau2_s -
                        2 * tau2_s * pred$sd^2),
      # Squared Gaussian Wasserstein distance as defined in the paper
      wd2 = mean(error^2 + (pred$sd - sqrt(tau2_s))^2),
      coverage = mean(truth >= pred$q025 & truth <= pred$q975),
      interval_score = mean(interval_score(pred$q025, pred$q975, truth)),
      mean_sd = mean(pred$sd),
      mean_interval_width = mean(pred$q975 - pred$q025),
      gamma_mean = NA_real_,
      gamma_sd = NA_real_,
      gamma_coverage = NA_real_,
      gamma_interval_score = NA_real_,
      mu_mean = NA_real_,
      mu_true = mu_true,
      mu_coverage = NA_real_,
      mu_interval_score = NA_real_
    )

    # Preferentiality parameter, when present
    hyperpar = model$fit$summary.hyperpar

    if ("Beta for j" %in% rownames(hyperpar)) {
      gamma = hyperpar["Beta for j", ]

      metrics$gamma_mean = gamma[["mean"]]
      metrics$gamma_sd = gamma[["sd"]]
      metrics$gamma_coverage = as.numeric(
        gamma_true >= gamma[["0.025quant"]] &
          gamma_true <= gamma[["0.975quant"]]
      )
      metrics$gamma_interval_score = interval_score(
        gamma[["0.025quant"]], gamma[["0.975quant"]], gamma_true
      )
    }

    # Population mean: use the posterior interval for the intercept
    fixed = model$fit$summary.fixed
    if (!is.null(fixed) && "b0" %in% rownames(fixed) &&
        all(c("mean", "0.025quant", "0.975quant") %in% colnames(fixed))) {
      mu = fixed["b0", ]
      metrics$mu_mean = mu[["mean"]]
      metrics$mu_coverage = as.numeric(
        mu_true >= mu[["0.025quant"]] &
          mu_true <= mu[["0.975quant"]]
      )
      metrics$mu_interval_score = interval_score(
        mu[["0.025quant"]], mu[["0.975quant"]], mu_true
      )
    }

    metrics
  })

  do.call(rbind, comparison)
}


summarize_model_metrics = function(metrics, group_vars) {
  if (!nrow(metrics)) return(data.frame())
  groups = group_summary_rows(metrics, c(group_vars, "model"))

  bind_rows(lapply(groups, function(d) {
    data.frame(
      d[1, group_vars, drop = FALSE],
      model = d$model[1],
      n_included = nrow(d),
      mean_mse = safe_mean(d$mse),
      mcse_mse = safe_mcse(d$mse),
      mean_mae = safe_mean(d$mae),
      mcse_mae = safe_mcse(d$mae),
      mean_bias = safe_mean(d$bias),
      mean_wd2 = safe_mean(d$wd2),
      mcse_wd2 = safe_mcse(d$wd2),
      mean_wd_zhong = safe_mean(d$wd_zhong),
      mcse_wd_zhong = safe_mcse(d$wd_zhong),
      mean_coverage = safe_mean(d$coverage),
      mean_interval_score = safe_mean(d$interval_score),
      mcse_interval_score = safe_mcse(d$interval_score),
      mean_interval_width = safe_mean(d$mean_interval_width),
      median_gamma = safe_median(d$gamma_median),
      gamma_bias = safe_mean(d$gamma_median - d$gamma_true),
      gamma_rmse = sqrt(safe_mean((d$gamma_median - d$gamma_true)^2)),
      gamma_coverage = safe_mean(d$gamma_coverage),
      mean_gamma_interval_score = safe_mean(d$gamma_interval_score),
      mean_mu = safe_mean(d$mu_mean),
      mu_bias = safe_mean(d$mu_mean - d$mu_true),
      mu_coverage = safe_mean(d$mu_coverage),
      mean_mu_interval_score = safe_mean(d$mu_interval_score),
      median_tau2 = safe_median(d$tau2_median),
      median_tau2_s = safe_median(d$tau2_s_median),
      median_tau2_B = safe_median(d$tau2_B_median),
      tau2_B_bias = safe_mean(d$tau2_B_median - d$tau2_B_true),
      tau2_B_coverage = safe_mean(d$tau2_B_coverage),
      mean_tau2_B_interval_score = safe_mean(d$tau2_B_interval_score),
      row.names = NULL
    )
  }))
}


summarize_fit_times = function(fit_log, group_vars) {
  if (!nrow(fit_log)) return(data.frame())
  groups = group_summary_rows(fit_log, c(group_vars, "model"))

  bind_rows(lapply(groups, function(d) {
    reused = d$reused %in% TRUE
    attempted = d[!reused, , drop = FALSE]
    data.frame(
      d[1, group_vars, drop = FALSE],
      model = d$model[1],
      # Reused rows retain the original fit duration. Report per-fit time,
      # while attempt/failure counts continue to exclude reuse.
      median_elapsed_sec = safe_median(d$elapsed_sec),
      n_attempted = nrow(attempted),
      n_reused = sum(reused),
      n_failed = sum(attempted$status %in% c("invalid", "error")),
      row.names = NULL
    )
  }))
}


summarize_paired_metrics = function(metrics, group_vars) {
  if (!nrow(metrics)) {
    return(list(paired_metrics = data.frame(), paired_summary = data.frame()))
  }

  measure_names = c("mse", "mae", "wd2", "wd_zhong", "interval_score")
  id_vars = c(group_vars, "replication")

  paired_metrics = reshape(
    metrics[c(id_vars, "model", measure_names)],
    idvar = id_vars, timevar = "model", direction = "wide"
  )

  for (measure in measure_names) {
    for (name in c("ps_melding", "melding", "ps_geo")) {
      column = paste(measure, name, sep = ".")
      if (!column %in% names(paired_metrics)) {
        paired_metrics[[column]] = NA_real_
      }
    }
  }

  for (measure in measure_names) {
    paired_metrics[[paste0("delta_", measure)]] =
      paired_metrics[[paste0(measure, ".melding")]] -
      paired_metrics[[paste0(measure, ".ps_melding")]]
  }

  paired_metrics$delta_mse_vs_geo =
    paired_metrics$mse.ps_geo - paired_metrics$mse.ps_melding
  paired_metrics$relative_mse_reduction = ifelse(
    is.finite(paired_metrics$mse.melding) & paired_metrics$mse.melding > 0,
    paired_metrics$delta_mse / paired_metrics$mse.melding,
    NA_real_
  )

  groups = group_summary_rows(paired_metrics, group_vars)
  paired_summary = bind_rows(lapply(groups, function(d) {
    mean_delta = safe_mean(d$delta_mse)
    mcse_delta = safe_mcse(d$delta_mse)
    valid = is.finite(d$delta_mse)

    data.frame(
      d[1, group_vars, drop = FALSE],
      n_paired = sum(valid),
      mean_delta_mse = mean_delta,
      mcse_delta_mse = mcse_delta,
      ci_low = mean_delta - 1.96 * mcse_delta,
      ci_high = mean_delta + 1.96 * mcse_delta,
      proportion_ps_better = if (any(valid)) {
        mean(d$delta_mse[valid] > 0)
      } else {
        NA_real_
      },
      mean_relative_mse_reduction = safe_mean(d$relative_mse_reduction),
      mean_delta_mae = safe_mean(d$delta_mae),
      mean_delta_wd2 = safe_mean(d$delta_wd2),
      mean_delta_wd_zhong = safe_mean(d$delta_wd_zhong),
      mean_delta_interval_score = safe_mean(d$delta_interval_score),
      n_paired_vs_geo = sum(is.finite(d$delta_mse_vs_geo)),
      mean_delta_mse_vs_geo = safe_mean(d$delta_mse_vs_geo),
      row.names = NULL
    )
  }))

  list(paired_metrics = paired_metrics, paired_summary = paired_summary)
}


summarize_pairwise_gains = function(paired_metrics, group_vars) {
  if (!nrow(paired_metrics)) return(data.frame())
  stopifnot(all(c(group_vars, "replication") %in% names(paired_metrics)),
            !anyDuplicated(paired_metrics[c(group_vars, "replication")]))
  measures = c("mse", "mae", "wd_zhong", "wd2", "interval_score")
  groups = group_summary_rows(paired_metrics, group_vars)

  bind_rows(lapply(groups, function(d) {
    bind_rows(lapply(c("melding", "ps_geo"), function(comparator) {
      bind_rows(lapply(measures, function(measure) {
        columns = paste(measure, c(comparator, "ps_melding"), sep = ".")
        if (!all(columns %in% names(d))) return(data.frame())
        baseline = d[[columns[1]]]
        adjusted = d[[columns[2]]]
        valid = is.finite(baseline) & is.finite(adjusted)
        delta = baseline[valid] - adjusted[valid]
        mean_gain = safe_mean(delta)
        mcse = safe_mcse(delta)
        relative = if (measure == "mse") {
          positive = valid & baseline > 0
          safe_mean((baseline[positive] - adjusted[positive]) / baseline[positive])
        } else NA_real_

        data.frame(
          d[1, group_vars, drop = FALSE],
          comparator = comparator, metric = measure,
          n_paired = length(delta), mean_gain = mean_gain, mcse_gain = mcse,
          ci_low = mean_gain - 1.96 * mcse,
          ci_high = mean_gain + 1.96 * mcse,
          proportion_ps_better = if (length(delta)) mean(delta > 0) else NA_real_,
          mean_relative_mse_reduction = relative,
          row.names = NULL
        )
      }))
    }))
  }))
}


format_pairwise_gain_table = function(data, group_vars = character()) {
  if (!nrow(data)) return(data.frame())
  metric_labels = c(
    mse = "MSE", mae = "MAE", wd_zhong = "WD (Zhong code)",
    wd2 = "W2^2 (corrected)", interval_score = "Interval score"
  )
  ordering = c(data[group_vars], list(
    match(data$comparator, c("melding", "ps_geo")),
    match(data$metric, names(metric_labels))
  ))
  data = data[do.call(order, ordering), , drop = FALSE]
  table = data.frame(
    Comparison = unname(c(
      melding = "Melding - PSmelding", ps_geo = "PSgeo - PSmelding"
    )[data$comparator]),
    Metric = unname(metric_labels[data$metric]),
    Pairs = data$n_paired,
    "Mean gain (95% MC CI)" = format_mc_interval(data$mean_gain, data$mcse_gain),
    "PSmelding wins (%)" = round(100 * data$proportion_ps_better, 1),
    check.names = FALSE
  )
  if (length(group_vars)) {
    table = cbind(data[group_vars], table)
    # All-condition tables are displayed separately by comparison, for MSE only.
    if (length(unique(data$comparator)) == 1L) table$Comparison = NULL
    if (length(unique(data$metric)) == 1L) table$Metric = NULL
  }
  rownames(table) = NULL
  table
}


summarize_completion = function(completion, group_vars) {
  if (!nrow(completion)) return(data.frame())
  groups = group_summary_rows(completion, group_vars)

  bind_rows(lapply(groups, function(d) {
    data.frame(
      d[1, group_vars, drop = FALSE],
      n_target = max(d$n_target),
      n_assigned = sum(d$n_assigned),
      n_readable = sum(d$n_readable),
      n_complete = sum(d$n_complete),
      n_issues = sum(d$n_issues),
      row.names = NULL
    )
  }))
}


record_issue = function(condition, replication, model, issue) {
  data.frame(
    condition, replication = replication, model = model, issue = issue,
    row.names = NULL, check.names = FALSE
  )
}


add_condition = function(data, condition) {
  for (name in names(condition)) {
    data[[name]] = condition[[name]][1]
  }
  data
}


read_saved_condition = function(
    directory, condition, replications, model_names,
    tau2_s_true, gamma_true, tau2_B_true, mu_true = 0,
    separate_noise = TRUE,
    n_target = length(replications)
) {
  metric_list = list()
  log_list = list()
  issue_list = list()
  n_readable = 0L
  n_complete = 0L

  for (replication in replications) {
    filename = sprintf("rep_%03d.rds", replication)
    fit_path = file.path(directory, "fits", filename)
    log_path = file.path(directory, "log", filename)

    saved = if (file.exists(fit_path)) {
      tryCatch(readRDS(fit_path), error = function(e) e)
    } else {
      NULL
    }

    if (is.null(saved)) {
      issue_list[[length(issue_list) + 1L]] = record_issue(
        condition, replication, NA_character_, "missing fit file"
      )
    } else if (inherits(saved, "error")) {
      issue_list[[length(issue_list) + 1L]] = record_issue(
        condition, replication, NA_character_,
        paste("unreadable fit:", conditionMessage(saved))
      )
    }

    # Read the independent log when available; older fits can embed the log
    fit_log = if (file.exists(log_path)) {
      tryCatch(readRDS(log_path), error = function(e) NULL)
    } else {
      NULL
    }

    external_log_ok = is.data.frame(fit_log) && nrow(fit_log) > 0L &&
      "model" %in% names(fit_log)
    if (!external_log_ok) {
      issue_list[[length(issue_list) + 1L]] = record_issue(
        condition, replication, NA_character_,
        "missing or unreadable separate log file"
      )
    }

    if (!external_log_ok && is.list(saved) &&
        is.data.frame(saved$log)) {
      fit_log = saved$log
    }

    if (is.data.frame(fit_log) && nrow(fit_log) > 0L &&
        "model" %in% names(fit_log)) {
      if (!"reused" %in% names(fit_log)) fit_log$reused = FALSE
      if (!"elapsed_sec" %in% names(fit_log)) fit_log$elapsed_sec = NA_real_
      if (!"status" %in% names(fit_log)) fit_log$status = NA_character_
      if (!"included" %in% names(fit_log)) fit_log$included = NA
      fit_log$replication = replication
      fit_log = add_condition(fit_log, condition)
      log_list[[length(log_list) + 1L]] = fit_log[c(
        names(condition), "replication", "model", "status",
        "included", "elapsed_sec", "reused"
      )]
    }

    if (!is.list(saved) || inherits(saved, "error")) next
    if (is.null(saved$data$prediction$true) ||
        !is.list(saved$results$models)) {
      issue_list[[length(issue_list) + 1L]] = record_issue(
        condition, replication, NA_character_,
        "missing prediction truth or model results"
      )
      next
    }

    n_readable = n_readable + 1L
    n_valid_models = 0L

    for (model_name in model_names) {
      model = saved$results$models[[model_name]]

      if (!is.list(model) || is.null(model$fit) ||
          is.null(model$prediction)) {
        issue_list[[length(issue_list) + 1L]] = record_issue(
          condition, replication, model_name, "missing or failed model"
        )
        next
      }

      one_model = list(models = setNames(list(model), model_name))
      metrics = tryCatch(
        evaluate_replication(
          results = one_model,
          truth = saved$data$prediction$true,
          tau2_s = tau2_s_true,
          gamma_true = gamma_true,
          mu_true = mu_true
        ),
        error = function(e) e
      )

      if (inherits(metrics, "error")) {
        issue_list[[length(issue_list) + 1L]] = record_issue(
          condition, replication, model_name,
          paste("evaluation error:", conditionMessage(metrics))
        )
        next
      }

      parameters = extract_parameters(model, separate_noise = separate_noise)
      for (name in names(parameters)) {
        metrics[[name]] = unname(parameters[name])
      }

      metrics$gamma_true = gamma_true
      metrics$mu_true = mu_true
      metrics$tau2_B_true = tau2_B_true
      metrics$tau2_B_coverage = if (separate_noise &&
                                    is.finite(tau2_B_true) &&
                                    is.finite(metrics$tau2_B_q025) &&
                                    is.finite(metrics$tau2_B_q975)) {
        as.numeric(tau2_B_true >= metrics$tau2_B_q025 &&
                     tau2_B_true <= metrics$tau2_B_q975)
      } else {
        NA_real_
      }
      metrics$tau2_B_interval_score = if (separate_noise &&
                                          is.finite(tau2_B_true) &&
                                          is.finite(metrics$tau2_B_q025) &&
                                          is.finite(metrics$tau2_B_q975)) {
        interval_score(metrics$tau2_B_q025, metrics$tau2_B_q975,
                       tau2_B_true)
      } else {
        NA_real_
      }

      metrics$replication = replication
      metrics = add_condition(metrics, condition)
      metric_list[[length(metric_list) + 1L]] = metrics
      n_valid_models = n_valid_models + 1L
    }

    if (n_valid_models == length(model_names) && external_log_ok) {
      n_complete = n_complete + 1L
    }
    rm(saved)
  }

  list(
    metrics = bind_rows(metric_list),
    log = bind_rows(log_list),
    issues = bind_rows(issue_list),
    completion = data.frame(
      condition,
      n_target = n_target,
      n_assigned = length(replications),
      n_readable = n_readable,
      n_complete = n_complete,
      n_issues = length(issue_list),
      row.names = NULL, check.names = FALSE
    )
  )
}


format_summary_table = function(data, columns, labels, order_by = NULL) {
  if (!is.data.frame(data) || !nrow(data)) return(data.frame())
  stopifnot(length(columns) == length(labels), all(columns %in% names(data)))

  output = data[columns]
  names(output) = labels
  if ("Model" %in% names(output)) {
    output$Model = unname(
      summary_model_labels()[as.character(output$Model)]
    )
  }
  if (!is.null(order_by)) {
    output = output[do.call(order, output[order_by]), , drop = FALSE]
  }
  rownames(output) = NULL
  output
}


format_fitting_time_table = function(summary) {
  return(format_summary_table(
    summary, c("model", "median_elapsed_sec"),
    c("Model", "Median seconds")
  ))
}


make_reproduction_tables = function(
    info, completion, selected_summary = data.frame()
) {
  design = info$scenarios
  tables = list(
    design = data.frame(
      Setting = c(
        "Conditions", "Target replications per condition",
        "Spatial ranges", "Preferentiality strengths", "Point counts",
        "Areal counts", "Mesh max.edge", "Mesh vertices",
        "Point-process integration", "Areal projection",
        "Separate point/areal noise", "Prior for gamma"
      ),
      Value = c(
        nrow(design), info$n_replications,
        paste(unique(design$rho), collapse = ", "),
        paste(unique(design$gamma), collapse = ", "),
        paste(unique(design$I), collapse = ", "),
        paste(unique(design$J), collapse = ", "),
        paste(info$mesh_max_edge, collapse = ", "),
        info$mesh_vertices,
        info$point_process_method,
        info$areal_method,
        isTRUE(info$separate_noise),
        paste(info$prior_gamma, collapse = ", ")
      )
    ),
    completion = completion,
    incomplete_conditions = subset(completion, n_complete < n_target)
  )

  if (nrow(selected_summary)) {
    tables$accuracy = format_summary_table(
      selected_summary,
      c("model", "n_included", "mean_mse", "mcse_mse",
        "mean_mae", "mean_wd2"),
      c("Model", "n", "MSE", "MCSE", "MAE", "W2^2")
    )
    tables$predictive_uncertainty = format_summary_table(
      selected_summary,
      c("model", "mean_coverage", "mean_interval_score",
        "mean_interval_width"),
      c("Model", "CP", "Inscore", "Mean width")
    )
    tables$mu_parameters = format_summary_table(
      selected_summary,
      c("model", "mean_mu", "mu_bias", "mu_coverage",
        "mean_mu_interval_score"),
      c("Model", "Mean mu", "Bias", "95% CP", "Inscore")
    )
    tables$gamma_parameters = format_summary_table(
      selected_summary,
      c("model", "median_gamma", "gamma_bias", "gamma_rmse",
        "gamma_coverage", "mean_gamma_interval_score", "median_tau2"),
      c("Model", "Median gamma", "Gamma bias", "Gamma RMSE",
        "95% CP", "Inscore", "Shared tau2")
    )
    if ("median_elapsed_sec" %in% names(selected_summary)) {
      tables$timing = format_fitting_time_table(selected_summary)
    }
  }

  tables
}


make_extension_tables = function(
    info, completion, summary = data.frame(),
    paired_summary = data.frame()
) {
  tables = list(
    design = data.frame(
      Setting = c(
        "Target replications", "Baseline rho", "Point count",
        "Areal count", "Baseline gamma", "Point noise variance",
        "Areal noise variance", "Mesh max.edge",
        "Point-process integration", "Areal projection",
        "Separate point/areal noise", "Prior for gamma"
      ),
      Value = c(
        info$n_replications,
        info$generating$rho,
        info$generating$I,
        info$generating$J,
        info$generating$gamma,
        info$generating$tau2_s,
        info$generating$tau2_B,
        paste(info$mesh_max_edge, collapse = ", "),
        info$point_process_method,
        info$areal_method,
        isTRUE(info$separate_noise),
        paste(info$prior_gamma, collapse = ", ")
      )
    ),
    completion = completion
  )

  if (!nrow(summary)) return(tables)

  tables$parameter_uncertainty = format_summary_table(
    summary,
    c("experiment", "condition", "model", "mu_coverage",
      "mean_mu_interval_score", "gamma_coverage",
      "mean_gamma_interval_score"),
    c("Experiment", "Condition", "Model", "Mu CP", "Mu Inscore",
      "Gamma CP", "Gamma Inscore"),
    order_by = c("Experiment", "Condition", "Model")
  )

  gamma_summary = subset(summary, experiment == "gamma_strength")
  gamma_paired = subset(paired_summary, experiment == "gamma_strength")
  if (nrow(gamma_summary)) {
    tables$gamma_accuracy = format_summary_table(
      gamma_summary,
      c("gamma_true", "model", "n_included", "mean_mse",
        "mcse_mse", "mean_wd2", "mean_coverage", "mean_interval_score"),
      c("True gamma", "Model", "n", "MSE", "MCSE",
        "W2^2", "95% CP", "Inscore"),
      order_by = c("True gamma", "Model")
    )
    tables$gamma_recovery = format_summary_table(
      subset(gamma_summary, model %in% c("ps_melding", "ps_geo")),
      c("gamma_true", "model", "median_gamma", "gamma_coverage",
        "mean_gamma_interval_score"),
      c("True gamma", "Model", "Median gamma", "95% CP", "Inscore")
    )
  }
  if (nrow(gamma_paired)) {
    tables$gamma_paired = format_summary_table(
      gamma_paired,
      c("gamma_true", "n_paired", "mean_delta_mse",
        "ci_low", "ci_high", "proportion_ps_better"),
      c("True gamma", "n pairs", "MSE gain", "CI lower",
        "CI upper", "P(PS better)")
    )
  }

  noise_summary = subset(summary, experiment == "measurement_noise")
  noise_paired = subset(paired_summary, experiment == "measurement_noise")
  if (nrow(noise_summary)) {
    tables$noise_accuracy = format_summary_table(
      noise_summary,
      c("tau2_B_true", "model", "n_included", "mean_mse",
        "mcse_mse", "mean_wd2", "mean_coverage", "mean_interval_score"),
      c("True areal variance", "Model", "n", "MSE",
        "MCSE", "W2^2", "95% CP", "Inscore"),
      order_by = c("True areal variance", "Model")
    )
    tables$noise_recovery = format_summary_table(
      subset(noise_summary, model %in% c("ps_melding", "melding")),
      c("tau2_B_true", "model", "median_tau2_B",
        "tau2_B_bias", "tau2_B_coverage", "mean_tau2_B_interval_score"),
      c("True areal variance", "Model", "Median fitted variance",
        "Bias", "95% CP", "Inscore")
    )
  }
  if (nrow(noise_paired)) {
    tables$noise_paired = format_summary_table(
      noise_paired,
      c("tau2_B_true", "n_paired", "mean_delta_mse",
        "ci_low", "ci_high", "mean_relative_mse_reduction"),
      c("True areal variance", "n pairs", "MSE gain",
        "CI lower", "CI upper", "Relative gain"),
      order_by = "True areal variance"
    )
  }

  mechanism_summary = subset(summary, experiment == "mechanism")
  mechanism_paired = subset(paired_summary, experiment == "mechanism")
  if (nrow(mechanism_summary)) {
    tables$mechanism_accuracy = format_summary_table(
      mechanism_summary,
      c("condition", "model", "n_included", "mean_mse",
        "mcse_mse", "mean_wd2", "mean_coverage", "mean_interval_score"),
      c("Mechanism", "Model", "n", "MSE",
        "MCSE", "W2^2", "95% CP", "Inscore")
    )
  }
  if (nrow(mechanism_paired)) {
    tables$mechanism_paired = format_summary_table(
      mechanism_paired,
      c("condition", "n_paired", "mean_delta_mse", "ci_low",
        "ci_high", "mean_delta_interval_score"),
      c("Mechanism", "n pairs", "MSE gain",
        "CI lower", "CI upper", "Interval-score gain")
    )
  }
  tables
}


# Format a mean and its 95 percent Monte Carlo confidence interval
format_mc_interval = function(estimate, mcse) {
  result = rep(NA_character_, length(estimate))
  valid_mean = is.finite(estimate)
  valid_interval = valid_mean & is.finite(mcse)

  result[valid_mean] = sprintf("%.3f", estimate[valid_mean])
  result[valid_interval] = sprintf(
    "%.3f (%.3f, %.3f)",
    estimate[valid_interval],
    estimate[valid_interval] - 1.96 * mcse[valid_interval],
    estimate[valid_interval] + 1.96 * mcse[valid_interval]
  )

  result
}


# Plotting --------------------------------------------------------------------


summary_model_labels = function() {
  c(
    ps_melding = "PSmelding",
    melding = "Melding",
    ps_geo = "PSgeo"
  )
}


summary_model_colours = function() {
  c(
    ps_melding = "#2E7ED1",
    melding = "#F86022",
    ps_geo = "#666666"
  )
}


spatial_theme = function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "grey50",
        fill = NA,
        linewidth = 0.4
      ),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 11
      )
    )
}


spatial_fill_scale = function(limits) {
  ggplot2::scale_fill_stepsn(
    colours = c(
      "#2E7ED1", # Blue
      "#80BCE3", # Light blue
      "#C6DCE1", # Pale blue
      "#F7ECA0", # Pale yellow
      "#FDC94B", # Yellow
      "#F86022", # Orange
      "#C21B1F"  # Red
    ),
    limits = limits,
    oob = scales::squish,
    show.limits = TRUE,
    name = "Value"
  )
}


plot_raster_map = function(data, value_col, title, limits) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = x,
      y = y,
      fill = .data[[value_col]]
    )
  ) +
    ggplot2::geom_raster() +
    spatial_fill_scale(limits) +
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::labs(title = title) +
    spatial_theme()
}


plot_point_map = function(data, title, limits) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = x, y = y, fill = value)
  ) +
    ggplot2::geom_point(
      shape = 21,
      colour = "black",
      size = 1.5,
      stroke = 0.3
    ) +
    spatial_fill_scale(limits) +
    ggplot2::guides(fill = "none") +
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::labs(title = title) +
    spatial_theme()
}


plot_sampling_scenarios = function(field, scenarios, points, limits = NULL) {
  # Check inputs
  if (!is.list(scenarios) ||
      !is.list(points) ||
      length(scenarios) == 0L ||
      is.null(names(scenarios)) ||
      !identical(names(scenarios), names(points))) {
    stop("scenarios and points must be named lists with matching names.")
  }

  if (!all(vapply(
    scenarios,
    function(x) {
      is.numeric(x) &&
        length(x) == nrow(field) &&
        all(is.finite(x))
    },
    logical(1)
  ))) {
    stop("Each scenario must contain log-intensity values aligned with field.")
  }

  # Common colour limits
  if (is.null(limits)) {

    max_value = max(1, ceiling(max(
      abs(unlist(scenarios)),
      abs(unlist(lapply(points, function(p) p$value)))
    )))

    limits = c(-max_value, max_value)
  }

  # Organize log-intensity surfaces
  intensity_data = do.call(rbind, lapply(names(scenarios), function(name) {

    data.frame(
      x = field$x,
      y = field$y,
      value = scenarios[[name]],
      scenario = name,
      panel = "Log-intensity"
    )
  }))

  # Organize point observations
  point_data = do.call(rbind, lapply(names(points), function(name) {

    data.frame(
      x = points[[name]]$x,
      y = points[[name]]$y,
      value = points[[name]]$value,
      scenario = name,
      panel = "Point observations"
    )
  }))

  # Fix row and column order
  intensity_data$scenario = factor(
    intensity_data$scenario,
    levels = names(scenarios)
  )

  point_data$scenario = factor(
    point_data$scenario,
    levels = names(scenarios)
  )

  intensity_data$panel = factor(
    intensity_data$panel,
    levels = c("Log-intensity", "Point observations")
  )

  point_data$panel = factor(
    point_data$panel,
    levels = c("Log-intensity", "Point observations")
  )

  # Plot
  ggplot2::ggplot() +

    # Log-intensity surfaces
    ggplot2::geom_raster(
      data = intensity_data,
      ggplot2::aes(x, y, fill = value)
    ) +

    # Point observations
    ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(x, y, fill = value),
      shape = 21,
      colour = "black",
      size = 1.5,
      stroke = 0.3
    ) +

    # Shared colour scale
    spatial_fill_scale(limits) +

    # Spatial domain
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +

    # Facet rows and columns
    ggplot2::facet_grid(
      panel ~ scenario,
      switch = "y"
    ) +

    # Plotting theme
    spatial_theme() +

    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(
        size = 11,
        face = "plain"
      ),
      strip.text.y.left = ggplot2::element_text(
        angle = 90,
        size = 10,
        face = "plain"
      ),
      strip.placement = "outside",
      panel.spacing = grid::unit(0.7, "lines"),
      legend.position = "none"
    )
}


plot_areal_map = function(
    areas,
    value_col = "value",
    title = "Areal observations",
    limits = NULL
) {
  # Colour limits
  if (is.null(limits)) {
    max_value = max(1, ceiling(max(abs(areas[[value_col]]))))
    limits = c(-max_value, max_value)
  }

  ggplot2::ggplot(
    areas,
    ggplot2::aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = .data[[value_col]]
    )
  ) +
    ggplot2::geom_rect(
      colour = "grey50",
      linewidth = 0.4
    ) +
    spatial_fill_scale(limits) +
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::labs(title = title) +
    spatial_theme() +
    ggplot2::theme(legend.position = "none")
}


plot_simulated_data = function(data, limits = NULL) {
  field = data$latent
  areas = data$areas
  points = data$points

  # Common colour limits
  if (is.null(limits)) {

    max_value = max(1, ceiling(max(
      abs(field$Z),
      abs(areas$value),
      abs(points$value)
    )))

    limits = c(-max_value, max_value)
  }

  # Latent surface
  p1 = plot_raster_map(
    data = field,
    value_col = "Z",
    title = "Simulated surface",
    limits = limits
  ) +
    ggplot2::theme(legend.position = "none")

  # Areal observations
  p2 = plot_areal_map(
    areas = areas,
    title = "Areal observations",
    limits = limits
  )

  # Point observations
  p3 = plot_point_map(
    data = points,
    title = "Point observations",
    limits = limits
  )

  # Combine into one row
  patchwork::wrap_plots(p1, p2, p3, nrow = 1)
}


zhong_prediction_reference = function(J_values = c(4, 25, 100)) {
  # Zhong et al. (2025), Supplementary Material, Table 1, preferential sampling:
  # I = 100, gamma = 1. The rho = 0.2 rows also appear in main-paper Table 2.
  # These are published means, not medians or replication-level observations.
  reference = data.frame(
    rho = rep(c(0.1, 0.2, 0.4), 7),
    gamma = 1,
    I = 100L,
    J = c(rep(c(4, 25, 100), each = 3), rep(0, 3),
          rep(c(4, 25, 100), each = 3)),
    model = rep(c("melding", "ps_geo", "ps_melding"), c(9, 3, 9)),
    mean_mse = c(
      1.05, 0.47, 0.15, 0.59, 0.26, 0.09, 0.38, 0.14, 0.04,
      0.78, 0.40, 0.16,
      0.68, 0.35, 0.14, 0.55, 0.26, 0.09, 0.38, 0.14, 0.04
    ),
    mean_mae = c(
      0.80, 0.52, 0.29, 0.60, 0.39, 0.23, 0.48, 0.29, 0.16,
      0.69, 0.48, 0.29,
      0.65, 0.45, 0.28, 0.58, 0.39, 0.23, 0.48, 0.29, 0.16
    ),
    # Published WD is compared with the released-code score, never corrected wd2.
    mean_wd_zhong = c(
      1.60, 0.84, 0.31, 1.26, 0.56, 0.19, 0.83, 0.28, 0.08,
      1.27, 0.71, 0.31,
      1.22, 0.66, 0.27, 1.05, 0.50, 0.18, 0.71, 0.27, 0.08
    )
  )
  reference$published_J = reference$J
  # PS Geo uses points only; repeat its J = 0 reference across displayed panels.
  point_only = subset(reference, model == "ps_geo")
  bind_rows(c(
    list(subset(reference, model != "ps_geo" & J %in% J_values)),
    lapply(unique(J_values), function(J) {
      d = point_only
      d$J = J
      d
    })
  ))
}


plot_reproduction_performance = function(
    metrics, gamma = 1, I = 100,
    rho_values = c(0.1, 0.2, 0.4),
    J_values = c(4, 25, 100),
    plot_metrics = c("MSE" = "mse", "MAE" = "mae", "W2²" = "wd2"),
    benchmark_reference = NULL
) {
  # Use replication-level scores for the reproduced preferential sampling cases
  if (!nrow(metrics)) return(NULL)
  stopifnot(length(plot_metrics) > 0L, !is.null(names(plot_metrics)),
            all(plot_metrics %in% c("mse", "mae", "wd2", "wd_zhong")))
  required = c("rho", "gamma", "I", "J", "model", unname(plot_metrics))
  stopifnot(all(required %in% names(metrics)))

  figure_data = metrics[
    metrics$gamma == gamma & metrics$I == I &
      metrics$rho %in% rho_values & metrics$J %in% J_values &
      metrics$model %in% names(summary_model_labels()),
    , drop = FALSE
  ]
  if (!nrow(figure_data)) return(NULL)

  # Use the requested measures so MSE and MAE/WD can be displayed separately.
  plot_data = bind_rows(lapply(names(plot_metrics), function(label) {
    data.frame(
      rho = figure_data$rho,
      J = figure_data$J,
      model = figure_data$model,
      metric = label,
      score = figure_data[[plot_metrics[[label]]]]
    )
  }))
  plot_data = plot_data[is.finite(plot_data$score), , drop = FALSE]
  if (!nrow(plot_data)) return(NULL)

  plot_data$rho = factor(plot_data$rho, levels = rho_values)
  plot_data$J = factor(
    plot_data$J, levels = J_values,
    labels = paste0(J_values, " areas\n", I, " points")
  )
  plot_data$model = factor(
    plot_data$model, levels = names(summary_model_labels())
  )
  plot_data$metric = factor(plot_data$metric, levels = names(plot_metrics))
  # Fixed model offsets align boxes and reference markers even if a fit is missing.
  model_offsets = setNames(
    (seq_along(summary_model_labels()) - 2) * 0.8 / 3,
    names(summary_model_labels())
  )
  plot_data$x_position = as.numeric(plot_data$rho) +
    unname(model_offsets[as.character(plot_data$model)])

  plot = ggplot2::ggplot(
    plot_data, ggplot2::aes(
      x = x_position, y = score, fill = model,
      group = interaction(rho, model)
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.7 / 3,
      colour = "grey25", linewidth = 0.3,
      outlier.colour = "grey40", outlier.size = 0.7
    ) +
    ggplot2::facet_grid(metric ~ J, scales = "free_y", drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(rho_values), labels = rho_values,
      limits = c(0.5, length(rho_values) + 0.5)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        ps_melding = "#1B9E77", melding = "#D95F02", ps_geo = "#7570B3"
      ),
      breaks = names(summary_model_labels()),
      labels = summary_model_labels()
    ) +
    ggplot2::labs(
      x = expression(rho),
      y = if (identical(unname(plot_metrics), "mse")) "MSE" else "Score",
      fill = "Model"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "grey45", fill = NA, linewidth = 0.4
      ),
      strip.background = ggplot2::element_rect(fill = "#EDF1F3"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (!is.null(benchmark_reference)) {
    # Match both metric and full scenario key. wd2 has no published benchmark.
    reference = bind_rows(lapply(names(plot_metrics), function(label) {
      measure = plot_metrics[[label]]
      reference_column = paste0("mean_", measure)
      if (measure == "wd2" || !reference_column %in% names(benchmark_reference)) {
        return(data.frame())
      }
      d = merge(
        benchmark_reference,
        unique(figure_data[
          is.finite(figure_data[[measure]]), c("rho", "gamma", "I", "J", "model")
        ]),
        by = c("rho", "gamma", "I", "J", "model")
      )
      if (!nrow(d)) return(data.frame())
      d$metric = label
      d$reference_score = d[[reference_column]]
      d[is.finite(d$reference_score), , drop = FALSE]
    }))
    if (nrow(reference)) {
      reference$rho = factor(reference$rho, levels = rho_values)
      reference$J = factor(
        reference$J, levels = J_values,
        labels = paste0(J_values, " areas\n", I, " points")
      )
      reference$model = factor(reference$model, levels = names(model_offsets))
      reference$metric = factor(
        reference$metric, levels = names(plot_metrics)
      )
      reference$x_position = as.numeric(reference$rho) +
        unname(model_offsets[as.character(reference$model)])
      means = stats::aggregate(
        score ~ rho + J + model + metric + x_position,
        data = plot_data[plot_data$metric %in% reference$metric, , drop = FALSE],
        FUN = mean
      )
      plot = plot +
        ggplot2::geom_point(
          data = means, ggplot2::aes(shape = "Reproduction mean"),
          fill = "white", colour = "black", size = 2.2, stroke = 0.6
        ) +
        ggplot2::geom_point(
          data = reference,
          ggplot2::aes(y = reference_score, shape = "Zhong mean"),
          fill = "black", colour = "white", size = 3, stroke = 0.4
        ) +
        ggplot2::scale_shape_manual(
          name = NULL,
          values = c("Reproduction mean" = 21, "Zhong mean" = 23)
        ) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(order = 1),
          shape = ggplot2::guide_legend(
            order = 2, override.aes = list(
              fill = c("white", "black"), colour = c("black", "white")
            )
          )
        ) +
        ggplot2::labs(
          subtitle = "Boxes show reproduction replications; markers compare means",
          caption = paste(
            "Zhong et al. (2025), Supplementary Table 1 (preferential sampling).",
            "PS Geo's point-only reference is repeated across areal counts.", sep = "\n"
          )
        )
    }
  }
  plot
}


plot_reproduction_mse_gain = function(paired_summary) {
  if (!nrow(paired_summary)) return(NULL)
  if (nrow(paired_summary) == 1L) {
    return(
      ggplot2::ggplot(
        paired_summary,
        ggplot2::aes(x = "Selected condition", y = mean_delta_mse)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::geom_pointrange(
          ggplot2::aes(ymin = ci_low, ymax = ci_high),
          colour = "#2E7ED1"
        ) +
        ggplot2::labs(
          x = NULL, y = "MSE(Melding) - MSE(PSmelding)"
        ) +
        ggplot2::theme_minimal()
    )
  }

  ggplot2::ggplot(
    paired_summary,
    ggplot2::aes(
      x = factor(J), y = mean_delta_mse,
      colour = factor(rho), group = factor(rho)
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0.12
    ) +
    ggplot2::facet_grid(gamma ~ I, labeller = ggplot2::label_both) +
    ggplot2::labs(
      x = "Number of areal observations (J)",
      y = "MSE(Melding) - MSE(PSmelding)", colour = expression(rho)
    ) +
    ggplot2::theme_minimal()
}


plot_gamma_prediction_mse = function(summary) {
  if (!nrow(summary)) return(NULL)
  ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x = gamma_true, y = mean_mse, colour = model, group = model
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_mse - 1.96 * mcse_mse,
        ymax = mean_mse + 1.96 * mcse_mse
      ),
      width = 0.025, na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = summary_model_colours(), labels = summary_model_labels()
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(summary$gamma_true))) +
    ggplot2::labs(
      x = expression(gamma), y = "Mean prediction MSE", colour = "Model"
    ) +
    ggplot2::theme_minimal()
}


plot_noise_prediction_mse = function(summary) {
  if (!nrow(summary)) return(NULL)
  ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x = tau2_B_true, y = mean_mse, colour = model, group = model
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_mse - 1.96 * mcse_mse,
        ymax = mean_mse + 1.96 * mcse_mse
      ),
      width = 0.01, na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_colour_manual(
      values = summary_model_colours(), labels = summary_model_labels()
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(summary$tau2_B_true))) +
    ggplot2::labs(
      x = expression(tau[B]^2), y = "Mean prediction MSE", colour = "Model",
      subtitle = "Extension runs; error bars are 95% Monte Carlo intervals"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
}


plot_gamma_parameter_recovery = function(summary) {
  if (!nrow(summary)) return(NULL)
  recovery = subset(
    summary,
    model %in% c("ps_melding", "ps_geo") &
      is.finite(median_gamma) & is.finite(gamma_true)
  )
  if (!nrow(recovery)) return(NULL)

  ggplot2::ggplot(
    recovery,
    ggplot2::aes(
      x = gamma_true, y = median_gamma, colour = model, group = model
    )
  ) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(
      values = summary_model_colours(), labels = summary_model_labels()
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(recovery$gamma_true))) +
    ggplot2::labs(
      x = "True gamma", y = "Median posterior gamma", colour = "Model"
    ) +
    ggplot2::theme_minimal()
}


plot_gamma_mse_gain = function(paired_summary) {
  if (!nrow(paired_summary)) return(NULL)
  ggplot2::ggplot(
    paired_summary,
    ggplot2::aes(x = gamma_true, y = mean_delta_mse)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_line(colour = "#2E7ED1") +
    ggplot2::geom_point(colour = "#2E7ED1", size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0.025
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(paired_summary$gamma_true))) +
    ggplot2::labs(
      x = expression(gamma),
      y = "MSE(Melding) - MSE(PSmelding)"
    ) +
    ggplot2::theme_minimal()
}


plot_noise_mse_gain = function(paired_summary) {
  if (!nrow(paired_summary)) return(NULL)
  ggplot2::ggplot(
    paired_summary,
    ggplot2::aes(x = tau2_B_true, y = mean_delta_mse)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_line(colour = "#F86022") +
    ggplot2::geom_point(colour = "#F86022", size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0.01
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(paired_summary$tau2_B_true))) +
    ggplot2::labs(
      x = expression(tau[B]^2),
      y = "MSE(Melding) - MSE(PSmelding)"
    ) +
    ggplot2::theme_minimal()
}


plot_mechanism_predictive_performance = function(metrics) {
  if (!nrow(metrics)) return(NULL)
  conditions = c("constant", "independent_field", "spatially_varying")
  measures = c("mse", "coverage", "interval_score")
  stopifnot(all(c("experiment", "condition", "model", measures) %in% names(metrics)))
  metrics = subset(metrics, experiment == "mechanism" & condition %in% conditions &
                     model %in% names(summary_model_labels()))
  if (!nrow(metrics)) return(NULL)
  groups = group_summary_rows(metrics, c("condition", "model"))
  summary = bind_rows(lapply(groups, function(d) {
    bind_rows(lapply(measures, function(measure) {
      values = d[[measure]]
      values = values[is.finite(values)]
      if (!length(values)) return(data.frame())
      if (measure == "coverage") values = 100 * values
      estimate = safe_mean(values)
      mcse = safe_mcse(values)
      data.frame(
        condition = d$condition[1], model = d$model[1], metric = measure,
        n = length(values), estimate = estimate,
        low = estimate - 1.96 * mcse, high = estimate + 1.96 * mcse
      )
    }))
  }))
  if (!nrow(summary)) return(NULL)
  summary$Condition = factor(
    summary$condition, levels = conditions,
    labels = c("Constant (reference)", "Independent field", "Spatially varying")
  )
  summary$Measure = factor(
    summary$metric, levels = measures,
    labels = c("Prediction MSE\nLower is better", "Coverage (%)\nTarget: 95%",
               "Interval score\nLower is better")
  )
  summary$Model = factor(
    summary$model, levels = rev(names(summary_model_labels())),
    labels = rev(unname(summary_model_labels()))
  )
  nominal = data.frame(
    Measure = factor("Coverage (%)\nTarget: 95%", levels = levels(summary$Measure)),
    target = 95
  )
  counts = range(summary$n)
  count_note = if (counts[1] == counts[2]) {
    sprintf("%d replications per model and mechanism.", counts[1])
  } else {
    sprintf("%d-%d usable replications per model, mechanism and metric.", counts[1], counts[2])
  }
  ggplot2::ggplot(summary, ggplot2::aes(x = estimate, y = Model, colour = model)) +
    ggplot2::geom_vline(
      data = nominal, ggplot2::aes(xintercept = target),
      linetype = "dashed", colour = "grey50"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = low, xmax = high), orientation = "y",
      width = 0.14, linewidth = 0.65
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_grid(Condition ~ Measure, scales = "free_x") +
    ggplot2::scale_colour_manual(values = summary_model_colours(), guide = "none") +
    ggplot2::labs(
      x = NULL, y = NULL,
      subtitle = "Means and 95% Monte Carlo intervals across replications",
      caption = paste(count_note,
        "Coverage is for nominal 95% spatial predictive intervals. Dashed line: 95% target.",
        sep = "\n")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#EDF1F3", colour = NA),
      panel.spacing = grid::unit(1.1, "lines")
    )
}


extract_mesh_edges = function(mesh) {
  triangles = mesh$graph$tv
  edge_ids = rbind(
    triangles[, c(1, 2)],
    triangles[, c(2, 3)],
    triangles[, c(3, 1)]
  )
  data.frame(
    x = mesh$loc[edge_ids[, 1], 1],
    y = mesh$loc[edge_ids[, 1], 2],
    xend = mesh$loc[edge_ids[, 2], 1],
    yend = mesh$loc[edge_ids[, 2], 2]
  )
}


plot_reproduction_mesh = function(mesh, data) {
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = extract_mesh_edges(mesh),
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = "grey80", linewidth = 0.2
    ) +
    ggplot2::geom_rect(
      data = data$areas,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA, colour = "black", linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = data$points,
      ggplot2::aes(x = x, y = y),
      colour = "red", size = 1.2
    ) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "black", fill = NA, linewidth = 0.7
      )
    )
}


plot_extension_mesh = function(mesh, areas, simulations) {
  scenario_points = bind_rows(lapply(names(simulations), function(name) {
    data.frame(
      x = simulations[[name]]$points$x,
      y = simulations[[name]]$points$y,
      scenario = name
    )
  }))
  scenario_points$scenario = factor(
    scenario_points$scenario, levels = names(simulations)
  )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = extract_mesh_edges(mesh),
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = "grey80", linewidth = 0.2
    ) +
    ggplot2::geom_rect(
      data = areas,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA, colour = "#F86022", linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = scenario_points,
      ggplot2::aes(x = x, y = y),
      colour = "#2E7ED1", size = 1.1, alpha = 0.85
    ) +
    ggplot2::facet_wrap(~ scenario, nrow = 1) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "grey50", fill = NA, linewidth = 0.5
      )
    )
}


read_summary_csv = function(directory, name) {
  path = file.path(directory, paste0(name, ".csv"))

  if (!file.exists(path)) {
    return(data.frame())
  }

  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


# Summary output helpers -------------------------------------------------------

save_tables = function(tables, directory) {
  table_dir = file.path(directory, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  available = names(tables)[vapply(tables, function(x) {
    is.data.frame(x) && nrow(x) > 0L
  }, logical(1))]

  for (name in available) {
    utils::write.csv(
      tables[[name]], file.path(table_dir, paste0(name, ".csv")),
      row.names = FALSE, na = ""
    )
  }
  available
}


save_plot = function(plot, directory, name, width = 7, height = 4) {
  if (is.null(plot)) return(FALSE)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = file.path(directory, paste0(name, ".png")),
    plot = plot, width = width, height = height, dpi = 300
  )
  TRUE
}


save_report_table = function(report_index, summary_dir, section, name, data) {
  if (!nrow(data)) return(report_index)
  destination = file.path(summary_dir, section, "tables")
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    data, file.path(destination, paste0(name, ".csv")), row.names = FALSE
  )
  report_index[[section]]$tables = unique(c(
    report_index[[section]]$tables, name
  ))
  report_index
}


save_report_plot = function(
    report_index, summary_dir, section, name, plot, width, height
) {
  destination = file.path(summary_dir, section)
  if (save_plot(plot, destination, name, width = width, height = height)) {
    report_index[[section]]$figures = unique(c(
      report_index[[section]]$figures, name
    ))
  }
  report_index
}


save_metrics = function(
    directory, completion, issues, log,
    metrics = data.frame(), summary = data.frame(),
    paired_metrics = data.frame(), paired_summary = data.frame()
) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    completion, file.path(directory, "completion.csv"), row.names = FALSE
  )
  utils::write.csv(
    issues, file.path(directory, "issues.csv"), row.names = FALSE
  )
  if (nrow(log)) {
    utils::write.csv(
      log, file.path(directory, "fit_log.csv"), row.names = FALSE
    )
  }
  if (nrow(metrics)) {
    utils::write.csv(
      metrics, file.path(directory, "replication_metrics.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      summary, file.path(directory, "model_summary.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      paired_metrics, file.path(directory, "paired_replications.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      paired_summary, file.path(directory, "paired_summary.csv"),
      row.names = FALSE
    )
  }
}


empty_index = function() {
  list(
    found = FALSE,
    ready = FALSE,
    focus = NULL,
    tables = character(),
    figures = character(),
    findings = character(),
    n_complete = 0L,
    n_assigned = 0L,
    n_conditions = 0L,
    n_fully_complete = 0L,
    n_issues = 0L
  )
}


# Parameter recovery from saved posterior quantiles ----------------------------

extract_recovery_quantiles = function(results, separate_noise = TRUE) {
  quantile_columns = c("0.025quant", "0.5quant", "0.975quant")
  blc_names = rownames(results$spde$param.inla$BLC)
  bind_rows(lapply(names(results$models), function(model) {
    fit = results$models[[model]]$fit
    if (is.null(fit)) return(data.frame())
    hp = fit$summary.hyperpar
    precision = grep("^Precision for .*Gaussian", rownames(hp), ignore.case = TRUE)
    blc = fit$summary.spde2.blc$i
    get_quantiles = function(table, row) {
      if (is.null(table) || length(row) != 1L || is.na(row) ||
          !all(quantile_columns %in% colnames(table))) return(rep(NA_real_, 3))
      as.numeric(table[row, quantile_columns])
    }
    spatial_quantiles = function(name) {
      if (is.null(blc) || nrow(blc) != length(blc_names)) return(rep(NA_real_, 3))
      get_quantiles(blc, match(name, blc_names))
    }
    point = get_quantiles(hp, precision[1])
    values = list()
    values[[if (separate_noise) "tau2_s" else "tau2"]] = 1 / rev(point)
    if (separate_noise) {
      values$tau2_B = if (model == "ps_geo") rep(NA_real_, 3) else {
        1 / rev(get_quantiles(hp, precision[2]))
      }
    }
    # Monotone transforms preserve quantiles; decreasing transforms reverse them.
    values$rho = exp(spatial_quantiles("range.nominal.1"))
    values$sigma2 = exp(spatial_quantiles("variance.nominal.1"))
    # For this stationary 2D Matern model with nu = 1,
    # sigma^2 * kappa^2 = 1 / (4*pi*tau_spde^2).
    # Check the saved joint linear-combination definition before using log(tau).
    B = results$spde$param.inla$BLC
    theta_rows = c("variance.nominal.1", "kappa.1", "tau.1")
    theta_ok = all(theta_rows %in% rownames(B)) && isTRUE(all.equal(
      as.numeric(B[theta_rows[1], ] + 2 * B[theta_rows[2], ]),
      c(-log(4 * pi), rep(0, ncol(B) - 1)) - 2 * B[theta_rows[3], ],
      check.attributes = FALSE
    ))
    values$theta = if (theta_ok) {
      exp(-log(4 * pi) - 2 * rev(spatial_quantiles("tau.1")))
    } else rep(NA_real_, 3)
    bind_rows(lapply(names(values), function(parameter) {
      q = values[[parameter]]
      valid = all(is.finite(q)) && all(q > 0) && all(diff(q) >= 0)
      data.frame(
        model = model, parameter = parameter,
        estimate = if (valid) q[2] else NA_real_,
        q025 = if (valid) q[1] else NA_real_,
        q975 = if (valid) q[3] else NA_real_,
        status = if (parameter == "tau2_B" && model == "ps_geo") {
          "not estimated"
        } else if (valid) "ok" else "posterior quantiles unavailable"
      )
    }))
  }))
}


read_recovery_cache = function(job, separate_noise) {
  source_info = file.info(job$path)
  fingerprint = c(size = source_info$size, mtime = as.numeric(source_info$mtime))
  cached = if (file.exists(job$cache)) {
    tryCatch(readRDS(job$cache), error = function(e) NULL)
  } else NULL
  if (!is.null(cached) && identical(cached$version, 1L) &&
      identical(cached$fingerprint, fingerprint) &&
      identical(cached$separate_noise, separate_noise)) return(cached)
  extracted = tryCatch({
    saved = readRDS(job$path)
    extract_recovery_quantiles(saved$results, separate_noise)
  }, error = function(e) conditionMessage(e))
  cached = list(
    version = 1L, fingerprint = fingerprint, separate_noise = separate_noise,
    data = if (is.data.frame(extracted)) extracted else data.frame(),
    error = if (is.character(extracted)) extracted else NA_character_
  )
  dir.create(dirname(job$cache), recursive = TRUE, showWarnings = FALSE)
  temporary = paste0(job$cache, ".", Sys.getpid(), ".tmp")
  saveRDS(cached, temporary)
  if (!file.rename(temporary, job$cache)) stop("Cannot save parameter cache: ", job$cache)
  cached
}


refresh_parameter_recovery = function(results_dir, summary_dir, section, workers = 2L) {
  info = readRDS(file.path(results_dir, "run_info.rds"))
  metrics = read_summary_csv(summary_dir, "replication_metrics")
  group_vars = if (section == "reproduction") c("rho", "gamma", "I", "J") else {
    c("experiment", "condition")
  }
  if (!nrow(metrics)) return(invisible(NULL))
  keys = c(group_vars, "replication")
  runs = unique(metrics[keys])
  relative = if (section == "reproduction") {
    file.path(paste0("rho_", runs$rho), paste0("gamma_", runs$gamma),
              paste0("I_", runs$I, "_J_", runs$J))
  } else {
    ifelse(runs$experiment == "mechanism", runs$condition,
           file.path(runs$experiment, runs$condition))
  }
  jobs = lapply(seq_len(nrow(runs)), function(i) list(
    path = file.path(results_dir, relative[i], "fits", sprintf("rep_%03d.rds", runs$replication[i])),
    cache = file.path(summary_dir, "parameter_cache", relative[i], sprintf("rep_%03d.rds", runs$replication[i]))
  ))
  message("Extracting parameter recovery for ", section, ": ", length(jobs), " saved fits (cached files are reused)")
  batches = split(seq_along(jobs), rep(seq_len(workers), length.out = length(jobs)))
  extract_batch = function(indices) {
    lapply(seq_along(indices), function(k) {
      i = indices[k]
      x = read_recovery_cache(jobs[[i]], isTRUE(info$separate_noise))
      if (k %% 25L == 0L || k == length(indices)) {
        message(section, " worker ", Sys.getpid(), ": ", k, "/", length(indices))
      }
      list(index = i, result = x)
    })
  }
  results = if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(batches, extract_batch, mc.cores = workers)
  } else lapply(batches, extract_batch)
  if (any(vapply(results, function(x) is.null(x) || inherits(x, "try-error"), logical(1)))) {
    stop("Parameter extraction worker failed; completed per-file caches are retained")
  }
  results = unlist(results, recursive = FALSE)
  issues = bind_rows(lapply(results, function(x) {
    if (is.na(x$result$error)) return(data.frame())
    cbind(runs[x$index, , drop = FALSE], issue = x$result$error)
  }))
  raw = bind_rows(lapply(results, function(x) {
    if (!nrow(x$result$data)) return(data.frame())
    cbind(runs[rep(x$index, nrow(x$result$data)), , drop = FALSE], x$result$data)
  }))
  # Use the same included model/replication set as the predictive summaries.
  raw = merge(raw, unique(metrics[c(keys, "model")]), by = c(keys, "model"), sort = FALSE)
  rho_true = if (section == "reproduction") raw$rho else rep(info$generating$rho, nrow(raw))
  areal_true = rep(info$generating$tau2_B, nrow(raw))
  if (section == "extensions") {
    truth = unique(metrics[c(keys, "tau2_B_true")])
    raw = merge(raw, truth, by = keys, sort = FALSE)
    rho_true = rep(info$generating$rho, nrow(raw))
    areal_true = raw$tau2_B_true
  }
  # The simulator uses scale=rho/2 and nu=1: kappa=sqrt(8)/rho.
  true_values = list(
    tau2 = rep(info$generating$tau2_s, nrow(raw)),
    tau2_s = rep(info$generating$tau2_s, nrow(raw)), tau2_B = areal_true,
    rho = rho_true, sigma2 = rep(info$generating$sigma2, nrow(raw)),
    theta = info$generating$sigma2 * 8 / rho_true^2
  )
  if (!isTRUE(info$separate_noise) &&
      !isTRUE(all.equal(info$generating$tau2_s, info$generating$tau2_B))) {
    stop("Shared variance has no unique truth when generating noise variances differ")
  }
  raw$true = vapply(seq_len(nrow(raw)), function(i) true_values[[raw$parameter[i]]][i], numeric(1))
  raw$bias = raw$estimate - raw$true
  raw$squared_error = raw$bias^2
  raw$coverage = ifelse(is.finite(raw$q025) & is.finite(raw$q975),
                        as.numeric(raw$q025 <= raw$true & raw$true <= raw$q975), NA_real_)
  raw$interval_score = interval_score(raw$q025, raw$q975, raw$true)
  raw = raw[do.call(order, raw[c(group_vars, "model", "parameter", "replication")]), ]
  summary = summarize_parameter_recovery(raw, group_vars)
  utils::write.csv(raw, file.path(summary_dir, "parameter_replications.csv"), row.names = FALSE, na = "")
  utils::write.csv(summary, file.path(summary_dir, "parameter_recovery.csv"), row.names = FALSE, na = "")
  if (!nrow(issues)) issues = data.frame(issue = character())
  utils::write.csv(issues, file.path(summary_dir, "parameter_issues.csv"), row.names = FALSE, na = "")
  invisible(summary)
}


summarize_parameter_recovery = function(data, group_vars) {
  bind_rows(lapply(group_summary_rows(data, c(group_vars, "model", "parameter")), function(d) {
    stopifnot(length(unique(d$true)) == 1L)
    valid = is.finite(d$estimate) & is.finite(d$q025) & is.finite(d$q975)
    cbind(d[1, c(group_vars, "model", "parameter"), drop = FALSE], data.frame(
      true = d$true[1], n_available = nrow(d), n = sum(valid),
      median_estimate = safe_median(d$estimate[valid]),
      bias = safe_mean(d$bias[valid]), rmse = sqrt(safe_mean(d$squared_error[valid])),
      coverage = safe_mean(d$coverage[valid]),
      mean_interval_score = safe_mean(d$interval_score[valid]),
      status = if (all(d$status == "not estimated")) "not estimated" else {
        if (all(valid)) "ok" else "incomplete posterior summaries"
      }
    ))
  }))
}


format_parameter_recovery = function(data, condition_columns = character()) {
  if (!nrow(data)) return(data.frame())
  labels = c(tau2 = "Shared tau^2", tau2_s = "Point tau_s^2", tau2_B = "Areal tau_B^2",
             rho = "Range rho", sigma2 = "Spatial sigma^2", theta = "Microergodic theta")
  data = data[do.call(order, c(data[condition_columns], list(
    match(data$parameter, names(labels)), data$true,
    match(data$model, names(summary_model_labels()))
  ))), ]
  out = data.frame(
    Parameter = unname(labels[data$parameter]), Truth = data$true,
    Model = unname(summary_model_labels()[data$model]), n = data$n,
    `Median estimate` = data$median_estimate, Bias = data$bias, RMSE = data$rmse,
    `95% coverage` = data$coverage, `Interval score` = data$mean_interval_score,
    Status = ifelse(data$status == "ok", "", data$status), check.names = FALSE
  )
  if (all(data$status == "ok")) out$Status = NULL
  if (length(condition_columns)) {
    conditions = data[condition_columns]
    column_labels = c(rho = "True rho", gamma = "True gamma", I = "I", J = "J",
                      condition = "Condition", experiment = "Experiment")
    names(conditions) = unname(column_labels[condition_columns])
    if ("Condition" %in% names(conditions)) {
      mechanism_labels = c(constant = "Constant", independent_field = "Independent field",
                           spatially_varying = "Spatially varying")
      matched = conditions$Condition %in% names(mechanism_labels)
      conditions$Condition[matched] = unname(mechanism_labels[conditions$Condition[matched]])
      conditions$Condition = sub("^tau2_B_", "tau_B^2 = ", conditions$Condition)
      conditions$Condition = sub("^gamma_", "gamma = ", conditions$Condition)
    }
    out = cbind(conditions, out)
  }
  out
}


plot_noise_parameter_recovery = function(data) {
  d = subset(data, experiment == "measurement_noise" & parameter == "tau2_B" &
               is.finite(estimate) & model %in% c("melding", "ps_melding"))
  if (!nrow(d)) return(NULL)
  ggplot2::ggplot(d, ggplot2::aes(
    x = true, y = estimate, group = factor(true), colour = model
  )) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_boxplot(width = 0.16, outlier.size = 0.6) +
    ggplot2::facet_wrap(~ model, labeller = ggplot2::as_labeller(summary_model_labels())) +
    ggplot2::scale_x_log10(breaks = sort(unique(d$true))) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = summary_model_colours(), labels = summary_model_labels()) +
    ggplot2::labs(x = expression("True areal variance " * tau[B]^2),
      y = "Posterior median areal variance",
      subtitle = "Across-replication distributions; dashed line: estimated = true (both axes logarithmic)") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "none")
}


# PDF tables wrap within the text width and can continue over page breaks.
# Very wide exports are split into column panels, retaining condition/model keys.
print_report_table = function(data, digits = 3) {
  if (!knitr::is_latex_output()) {
    print(knitr::kable(data, digits = digits, row.names = FALSE))
    return(invisible(data))
  }
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("PDF table layout requires the installed kableExtra package")
  }
  keys = intersect(c("Experiment", "Condition", "True rho", "True gamma", "I", "J",
                     "experiment", "condition", "rho", "gamma", "Model", "model",
                     "Parameter", "parameter"), names(data))
  # Leave room for at least four result columns per panel.
  keys = head(keys, 7L)
  columns = setdiff(names(data), keys)
  panels = if (ncol(data) <= 11L) list(names(data)) else {
    lapply(split(columns, ceiling(seq_along(columns) / (11L - length(keys)))),
           function(x) c(keys, x))
  }
  for (i in seq_along(panels)) {
    d = data[panels[[i]]]
    if (length(panels) > 1L) {
      cat(sprintf("\n\\par\\noindent\\textit{Column panel %d of %d}\\par\n", i, length(panels)))
    }
    # Allocate width by the longest word, allowing long paths/identifiers to wrap.
    weight = vapply(seq_along(d), function(j) {
      values = if (is.numeric(d[[j]])) format(round(d[[j]], digits), trim = TRUE) else as.character(d[[j]])
      words = unlist(strsplit(c(names(d)[j], values), "[[:space:]_/]+"))
      words = words[!is.na(words)]
      max(4, min(24, max(nchar(words), 0)))
    }, numeric(1))
    # Available width is \linewidth less two 3pt column paddings per column.
    widths = sprintf("\\dimexpr %.5f\\linewidth - %.5fpt\\relax",
                      weight / sum(weight), 6 * ncol(d) * weight / sum(weight))
    tab = knitr::kable(d, format = "latex", booktabs = TRUE, longtable = TRUE,
                      digits = digits, row.names = FALSE, linesep = "")
    tab = kableExtra::kable_styling(tab, latex_options = "repeat_header",
                                  repeat_header_text = "")
    tab = as.character(tab)
    alignment = ifelse(vapply(d, is.numeric, logical(1)), "raggedleft", "raggedright")
    spec = paste0(">{\\", alignment, "\\arraybackslash}p{", widths, "}", collapse = "")
    opening = regexpr("\\\\begin\\{longtable\\}(?:\\[[^]]*\\])?\\{[^\n]*\\}", tab, perl = TRUE)
    if (opening[1] < 0L) stop("Cannot find the LaTeX longtable column specification")
    regmatches(tab, opening) = paste0("\\begin{longtable}{", spec, "}")
    # Escaped underscores and slashes otherwise form unbreakable file/column names.
    tab = gsub("\\_", "\\_\\allowbreak{}", as.character(tab), fixed = TRUE)
    tab = gsub("/", "/\\allowbreak{}", tab, fixed = TRUE)
    cat("\n\\begingroup\n\\fontsize{8.5}{10.5}\\selectfont\n",
        "\\setlength{\\tabcolsep}{3pt}\n\\renewcommand{\\arraystretch}{1.12}\n",
        tab, "\n\\endgroup\n\n", sep = "")
  }
  invisible(data)
}
