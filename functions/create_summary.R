# Helper functions ------------------------------------------------------------

safe_median = function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}


# Extract model diagnostics
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
  
  return(row)
}






# Plotting theme

spatial_theme <- function() {
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


# Colour scale
spatial_fill_scale <- function(limits) {
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

# Plot raster
plot_raster_map <- function(data, value_col, title, limits) {
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


# Plot points
plot_point_map <- function(data, title, limits) {
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


# Plot simulation scenarios
plot_sampling_scenarios <- function(field, scenarios, points, limits = NULL) {
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
    
    max_value <- max(1, ceiling(max(
      abs(unlist(scenarios)),
      abs(unlist(lapply(points, function(p) p$value)))
    )))
    
    limits <- c(-max_value, max_value)
  }
  
  # Organize log-intensity surfaces
  intensity_data <- do.call(rbind, lapply(names(scenarios), function(name) {
    
    data.frame(
      x = field$x,
      y = field$y,
      value = scenarios[[name]],
      scenario = name,
      panel = "Log-intensity"
    )
  }))
  
  # Organize point observations
  point_data <- do.call(rbind, lapply(names(points), function(name) {
    
    data.frame(
      x = points[[name]]$x,
      y = points[[name]]$y,
      value = points[[name]]$value,
      scenario = name,
      panel = "Point observations"
    )
  }))
  
  # Fix row and column order
  intensity_data$scenario <- factor(
    intensity_data$scenario,
    levels = names(scenarios)
  )
  
  point_data$scenario <- factor(
    point_data$scenario,
    levels = names(scenarios)
  )
  
  intensity_data$panel <- factor(
    intensity_data$panel,
    levels = c("Log-intensity", "Point observations")
  )
  
  point_data$panel <- factor(
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


# Plot areal observations
plot_areal_map <- function(
    areas,
    value_col = "value",
    title = "Areal observations",
    limits = NULL
) {
  # Colour limits
  if (is.null(limits)) {
    max_value <- max(1, ceiling(max(abs(areas[[value_col]]))))
    limits <- c(-max_value, max_value)
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

# Plot simulated data
plot_simulated_data <- function(data, limits = NULL) {
  field <- data$latent
  areas <- data$areas
  points <- data$points
  
  # Common colour limits
  if (is.null(limits)) {
    
    max_value <- max(1, ceiling(max(
      abs(field$Z),
      abs(areas$value),
      abs(points$value)
    )))
    
    limits <- c(-max_value, max_value)
  }
  
  # Latent surface
  p1 <- plot_raster_map(
    data = field,
    value_col = "Z",
    title = "Simulated surface",
    limits = limits
  ) +
    ggplot2::theme(legend.position = "none")
  
  # Areal observations
  p2 <- plot_areal_map(
    areas = areas,
    title = "Areal observations",
    limits = limits
  )
  
  # Point observations
  p3 <- plot_point_map(
    data = points,
    title = "Point observations",
    limits = limits
  )
  
  # Combine into one row
  patchwork::wrap_plots(p1, p2, p3, nrow = 1)
}


# Evaluation helpers -----------------------------------------------------------

# Prediction accuracy
evaluate <- function(pred, truth) {
  return(c(
    rmse = sqrt(mean((pred$mean - truth)^2)),
    coverage = mean(truth >= pred$q025 & truth <= pred$q975),
    mean_sd = mean(pred$sd)
  ))
}

# Summarize PS melding results
summarize_ps <- function(results, truth) {
  comparison <- lapply(results, function(result) {
    fit <- result$models$ps_melding$fit
    pred <- result$models$ps_melding$prediction
    
    gamma <- fit$summary.hyperpar["Beta for j", ]
    
    c(
      gamma = gamma[["mean"]],
      gamma_sd = gamma[["sd"]],
      gamma_q025 = gamma[["0.025quant"]],
      gamma_q975 = gamma[["0.975quant"]],
      evaluate(pred, truth)
    )
  })
  
  return(do.call(rbind, comparison))
}

# Compare posterior mean predictions
compare_predictions <- function(pred1, pred2) {
  stopifnot(
    isTRUE(all.equal(
      pred1[, c("x", "y")],
      pred2[, c("x", "y")]
    ))
  )
  
  difference <- pred1$mean - pred2$mean
  
  return(c(
    mean_difference = mean(difference),
    mean_abs_difference = mean(abs(difference)),
    max_abs_difference = max(abs(difference)),
    correlation = cor(pred1$mean, pred2$mean)
  ))
}



# Interval score for a prediction or credible interval
interval_score <- function(lower, upper, truth, alpha = 0.05) {
  return((upper - lower) +
           (2 / alpha) * pmax(lower - truth, 0) +
           (2 / alpha) * pmax(truth - upper, 0))
}


# Evaluate all models for one simulation replication
evaluate_replication <- function(results, truth, tau2_s = 0.1, gamma_true = 1) {
  comparison <- lapply(names(results$models), function(name) {
    model <- results$models[[name]]
    pred <- model$prediction
    
    if (nrow(pred) != length(truth)) {
      stop("Prediction and truth lengths do not match.")
    }
    
    if (is.null(pred$median)) {
      stop("Posterior median is missing. Refit using the updated fit_models.R.")
    }
    
    # Predictive accuracy and uncertainty
    error <- pred$mean - truth
    
    metrics <- data.frame(
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
      gamma_interval_score = NA_real_
    )
    
    # Preferentiality parameter, when present
    hyperpar <- model$fit$summary.hyperpar
    
    if ("Beta for j" %in% rownames(hyperpar)) {
      gamma <- hyperpar["Beta for j", ]
      
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
    
    return(metrics)
  })
  
  return(do.call(rbind, comparison))
}



# Summarize model performance across replications
summarize_replications <- function(metrics, gamma_true = 1) {
  # Monte Carlo standard error across independent replications
  mcse <- function(x) {
    if (length(x) < 2L) {
      return(NA_real_)
    }
    
    return(stats::sd(x) / sqrt(length(x)))
  }
  
  # Average predictive performance
  metric_names <- c("mse", "mae", "wd2", "bias", "coverage", "interval_score")
  
  prediction_n <- stats::aggregate(
    metrics[, "mse", drop = FALSE],
    by = list(model = metrics$model),
    FUN = length
  )
  
  names(prediction_n)[2] <- "n_replications"
  
  prediction_mean <- stats::aggregate(
    metrics[, metric_names, drop = FALSE],
    by = list(model = metrics$model),
    FUN = mean
  )
  
  # Between-replication standard deviations
  prediction_sd <- stats::aggregate(
    metrics[, metric_names, drop = FALSE],
    by = list(model = metrics$model),
    FUN = sd
  )
  
  # Monte Carlo standard errors of average predictive performance
  prediction_mcse <- stats::aggregate(
    metrics[, metric_names, drop = FALSE],
    by = list(model = metrics$model),
    FUN = mcse
  )
  
  # Zhong's original Wasserstein calculation for reproduction
  prediction_reproduction <- stats::aggregate(
    metrics[, "wd_zhong", drop = FALSE],
    by = list(model = metrics$model),
    FUN = mean
  )
  
  prediction_reproduction_mcse <- stats::aggregate(
    metrics[, "wd_zhong", drop = FALSE],
    by = list(model = metrics$model),
    FUN = mcse
  )
  
  # Additional uncertainty diagnostics
  prediction_diagnostics <- stats::aggregate(
    metrics[, c("mean_interval_width", "mean_sd"), drop = FALSE],
    by = list(model = metrics$model),
    FUN = mean
  )
  
  prediction_diagnostics_mcse <- stats::aggregate(
    metrics[, c("mean_interval_width", "mean_sd"), drop = FALSE],
    by = list(model = metrics$model),
    FUN = mcse
  )
  
  # Recovery of gamma for preferential sampling models
  ps_metrics <- metrics[!is.na(metrics$gamma_mean), , drop = FALSE]
  gamma_summary <- NULL
  gamma_mcse <- NULL
  
  # A constant gamma has no true value under spatially varying gamma(s)
  if (is.finite(gamma_true) && nrow(ps_metrics) > 0L) {
    gamma_summary <- do.call(rbind, lapply(
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
    
    gamma_mcse <- do.call(rbind, lapply(
      split(ps_metrics, ps_metrics$model),
      function(x) {
        gamma_sq_error <- (x$gamma_mean - gamma_true)^2
        gamma_rmse <- sqrt(mean(gamma_sq_error))
        
        # Delta-method MCSE for gamma RMSE
        gamma_rmse_mcse <- if (gamma_rmse > 0) {
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
    rownames(gamma_summary) <- NULL
  }
  
  if (!is.null(gamma_mcse)) {
    rownames(gamma_mcse) <- NULL
  }
  
  return(list(
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
  ))
}


# Summarize fitting outcomes ---------------------------------------------------

summarize_fit_log <- function(fit_log, models = unique(fit_log$model)) {
  comparison <- lapply(models, function(name) {
    x <- fit_log[fit_log$model == name, , drop = FALSE]
    
    # Reused fits are not additional fitting attempts
    reused <- if ("reused" %in% names(x)) {
      !is.na(x$reused) & x$reused
    } else {
      rep(FALSE, nrow(x))
    }
    
    attempted <- x[!reused, , drop = FALSE]
    
    data.frame(
      model = name,
      attempted = nrow(attempted),
      included = sum(attempted$included),
      warnings = sum(attempted$status == "warning"),
      invalid = sum(attempted$status == "invalid"),
      errors = sum(attempted$status == "error"),
      reused = sum(reused),
      mean_elapsed_sec = if (nrow(attempted) > 0L) {
        mean(attempted$elapsed_sec)
      } else {
        NA_real_
      }
    )
  })
  
  comparison <- do.call(rbind, comparison)
  rownames(comparison) <- NULL
  
  return(comparison)
}


# Summarize a posterior marginal ----------------------------------------------

summarize_marginal <- function(marginal) {
  if (is.null(marginal)) {
    stop("Posterior marginal is unavailable.")
  }
  
  posterior_mean <- INLA::inla.emarginal(function(x) x, marginal)
  posterior_second <- INLA::inla.emarginal(function(x) x^2, marginal)
  quantiles <- INLA::inla.qmarginal(c(0.025, 0.975), marginal)
  
  return(c(
    mean = posterior_mean,
    sd = sqrt(max(0, posterior_second - posterior_mean^2)),
    q025 = quantiles[1],
    q975 = quantiles[2]
  ))
}


# Evaluate parameter recovery for one simulation replication ------------------

evaluate_parameter_replication <- function(results, truth, separate_noise = TRUE) {
  comparison <- lapply(names(results$models), function(name) {
    fit <- results$models[[name]]$fit
    hyperpar <- fit$marginals.hyperpar
    
    if (is.null(hyperpar)) {
      stop("Posterior hyperparameter marginals are missing for ", name, ".")
    }
    
    # Collect posterior marginals on the parameter scales of interest
    marginals <- list()
    
    if ("Beta for j" %in% names(hyperpar)) {
      marginals$gamma <- hyperpar[["Beta for j"]]
    }
    
    # Transform Gaussian observation precisions into variances
    precision_names <- grep(
      "^Precision for .*Gaussian",
      names(hyperpar),
      value = TRUE,
      ignore.case = TRUE
    )
    
    if (length(precision_names) < 1L) {
      stop("Gaussian observation precision is missing for ", name, ".")
    }
    
    marginals$tau2_s <- INLA::inla.tmarginal(
      function(x) 1 / x,
      hyperpar[[precision_names[1]]]
    )
    
    if (separate_noise && name %in% c("melding", "ps_melding")) {
      if (length(precision_names) < 2L) {
        stop("Areal observation precision is missing for ", name, ".")
      }
      
      marginals$tau2_B <- INLA::inla.tmarginal(
        function(x) 1 / x,
        hyperpar[[precision_names[2]]]
      )
    }
    
    # Transform SPDE parameters into nominal range and variance
    if (!is.null(results$spde)) {
      spde_result <- INLA::inla.spde2.result(
        inla = fit,
        name = "i",
        spde = results$spde,
        do.transform = TRUE
      )
      
      marginals$rho <- spde_result$marginals.range.nominal[[1]]
      marginals$sigma2 <- spde_result$marginals.variance.nominal[[1]]
    }
    
    # Only evaluate parameters that exist in the fitted model
    parameters <- intersect(names(marginals), names(truth))
    parameters <- parameters[is.finite(truth[parameters])]
    
    do.call(rbind, lapply(parameters, function(parameter) {
      posterior <- summarize_marginal(marginals[[parameter]])
      true <- as.numeric(truth[[parameter]])
      
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
  
  comparison <- do.call(rbind, comparison)
  
  if (!is.null(comparison)) {
    rownames(comparison) <- NULL
  }
  
  return(comparison)
}


# Summarize parameter recovery across replications -----------------------------

summarize_parameter_recovery <- function(metrics) {
  if (is.null(metrics) || nrow(metrics) == 0L) {
    return(list(parameter_mean = NULL, parameter_mcse = NULL))
  }
  
  mcse <- function(x) {
    if (length(x) < 2L) {
      return(NA_real_)
    }
    
    return(stats::sd(x) / sqrt(length(x)))
  }
  
  groups <- split(
    metrics,
    list(metrics$model, metrics$parameter),
    drop = TRUE
  )
  
  parameter_mean <- do.call(rbind, lapply(groups, function(x) {
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
  
  parameter_mcse <- do.call(rbind, lapply(groups, function(x) {
    rmse <- sqrt(mean(x$squared_error))
    
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
  
  rownames(parameter_mean) <- NULL
  rownames(parameter_mcse) <- NULL
  
  return(list(
    parameter_mean = parameter_mean,
    parameter_mcse = parameter_mcse
  ))
}
