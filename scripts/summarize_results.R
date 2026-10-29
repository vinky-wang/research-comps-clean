source("functions/create_summary.R")


# Settings --------------------------------------------------------------------

reproduction_dir = "results/reproduction"
extension_dir = "results/extensions"
summary_dir = "results/summary"
update_tables = FALSE
# Set TRUE to refresh parameter recovery from saved fits (unchanged files are cached).
update_parameter_tables = update_tables
parameter_workers = 2L
create_plots = FALSE
update_reproduction_tables = TRUE
update_report_tables = TRUE
create_report_plots = TRUE

# Select the individual figures to regenerate
make_plots = list(
  reproduction_performance = TRUE,
  reproduction_mean_mse_comparison = TRUE,
  reproduction_paired_gain = TRUE,
  reproduction_simulated_data = TRUE,
  reproduction_mesh = TRUE,
  gamma_prediction_mse = TRUE,
  gamma_recovery = TRUE,
  gamma_paired_gain = TRUE,
  noise_prediction_mse = TRUE,
  noise_paired_gain = TRUE,
  mechanism_predictive_performance = TRUE,
  extension_sampling_scenarios = TRUE,
  extension_mesh = TRUE
)

reproduction_summary_dir = file.path(summary_dir, "reproduction")
extension_summary_dir = file.path(summary_dir, "extensions")
dir.create(reproduction_summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extension_summary_dir, recursive = TRUE, showWarnings = FALSE)


# Summarize results -----------------------------------------------------------

index_path = file.path(summary_dir, "report_index.rds")

if (update_tables) {
  report_index = list(
    generated_at = Sys.time(),
    reproduction = empty_index(),
    extensions = empty_index()
  )

  # Keep references to saved figures when only some plots are regenerated
  if (file.exists(index_path)) {
    previous_index = readRDS(index_path)
    for (section in c("reproduction", "extensions")) {
      figures = previous_index[[section]]$figures
      if (length(figures)) {
        existing = file.exists(file.path(
          summary_dir, section, paste0(figures, ".png")
        ))
        report_index[[section]]$figures = figures[existing]
      }
    }
  }
} else {
  if (!file.exists(index_path)) {
    stop("No existing summary index: run with update_tables = TRUE first")
  }
  report_index = readRDS(index_path)
}


# Reproduction ----------------------------------------------------------------

reproduction_info_path = file.path(reproduction_dir, "run_info.rds")
if (!file.exists(reproduction_info_path)) {
  # Support the earlier completed single-condition run
  single_case_dir = "results/reproduction/rho_0.2_gamma_1_I_100_J_25"
  if (file.exists(file.path(single_case_dir, "run_info.rds"))) {
    reproduction_info_path = file.path(single_case_dir, "run_info.rds")
  }
}

if (update_tables && file.exists(reproduction_info_path)) {
  message("Summarizing reproduction...")
  info = readRDS(reproduction_info_path)
  design = info$scenarios
  stopifnot(
    is.data.frame(design),
    all(c("rho", "gamma", "I", "J") %in% names(design))
  )
  replications = if (!is.null(info$assigned_replications)) {
    info$assigned_replications
  } else {
    seq_len(info$n_replications)
  }
  results_dir = dirname(reproduction_info_path)

  runs = lapply(seq_len(nrow(design)), function(i) {
    condition = design[i, , drop = FALSE]
    directory = file.path(
      results_dir,
      paste0("rho_", condition$rho),
      paste0("gamma_", condition$gamma),
      paste0("I_", condition$I, "_J_", condition$J)
    )
    message("  Reading ", basename(directory))
    read_saved_condition(
      directory = directory,
      condition = condition,
      replications = replications,
      model_names = info$models,
      tau2_s_true = info$generating$tau2_s,
      gamma_true = condition$gamma,
      tau2_B_true = info$generating$tau2_B,
      mu_true = info$generating$mu,
      separate_noise = isTRUE(info$separate_noise),
      n_target = info$n_replications
    )
  })

  metrics = bind_rows(lapply(runs, function(x) x$metrics))
  log = bind_rows(lapply(runs, function(x) x$log))
  issues = bind_rows(lapply(runs, function(x) x$issues))
  completion = summarize_completion(
    bind_rows(lapply(runs, function(x) x$completion)),
    group_vars = c("rho", "gamma", "I", "J")
  )

  ready = nrow(metrics) > 0L
  summary = data.frame()
  paired_metrics = data.frame()
  paired_summary = data.frame()
  selected_summary = data.frame()
  selected_paired = data.frame()
  focus = NULL

  if (ready) {
    summary = summarize_model_metrics(metrics, c("rho", "gamma", "I", "J"))
    timing = summarize_fit_times(log, c("rho", "gamma", "I", "J"))
    if (nrow(timing)) {
      summary = merge(
        summary, timing,
        by = c("rho", "gamma", "I", "J", "model"),
        all.x = TRUE, sort = FALSE
      )
    }
    paired = summarize_paired_metrics(metrics, c("rho", "gamma", "I", "J"))
    paired_metrics = paired$paired_metrics
    paired_summary = paired$paired_summary

    available = unique(summary[c("rho", "gamma", "I", "J")])
    focus = subset(available, rho == 0.2 & gamma == 1 & I == 100 & J == 25)
    if (!nrow(focus)) focus = available[1, , drop = FALSE]
    focus = focus[1, , drop = FALSE]
    selected_summary = subset(
      summary,
      rho == focus$rho & gamma == focus$gamma &
        I == focus$I & J == focus$J
    )
    selected_paired = subset(
      paired_summary,
      rho == focus$rho & gamma == focus$gamma &
        I == focus$I & J == focus$J
    )
  }

  save_metrics(
    reproduction_summary_dir, completion, issues, log,
    metrics, summary, paired_metrics, paired_summary
  )
  tables = make_reproduction_tables(
    info, completion, selected_summary
  )
  if (nrow(issues)) tables$issues_preview = utils::head(issues, 10)

  report_index$reproduction$found = TRUE
  report_index$reproduction$ready = ready
  report_index$reproduction$focus = if (!is.null(focus)) {
    as.list(focus[1, , drop = FALSE])
  } else {
    NULL
  }
  report_index$reproduction$tables = save_tables(
    tables, reproduction_summary_dir
  )
  report_index$reproduction$n_complete = sum(completion$n_complete)
  report_index$reproduction$n_assigned = sum(completion$n_assigned)
  report_index$reproduction$n_conditions = nrow(completion)
  report_index$reproduction$n_fully_complete = sum(
    completion$n_complete == completion$n_target
  )
  report_index$reproduction$n_issues = nrow(issues)

  if (ready && nrow(selected_paired)) {
    d = selected_paired[1, ]
    report_index$reproduction$findings = sprintf(
      paste0(
        "For the detailed reproduction condition (rho = %s, gamma = %s, ",
        "I = %s, J = %s), the mean paired MSE gain is %.3f ",
        "(95%% Monte Carlo interval [%.3f, %.3f]; %d paired replications)."
      ),
      focus$rho, focus$gamma, focus$I, focus$J,
      d$mean_delta_mse, d$ci_low, d$ci_high, d$n_paired
    )
  }
} else if (update_tables) {
  message("Reproduction metadata not found: ", reproduction_info_path)
}


# Extensions ------------------------------------------------------------------

extension_info_path = file.path(extension_dir, "run_info.rds")
if (update_tables && file.exists(extension_info_path)) {
  message("Summarizing extensions...")
  info = readRDS(extension_info_path)
  replications = if (!is.null(info$assigned_replications)) {
    info$assigned_replications
  } else {
    seq_len(info$n_replications)
  }

  scenarios = as.character(info$scenarios)
  tau2_B_values = info$measurement_noise$tau2_B_values
  gamma_values = info$gamma_strength$gamma_values
  design = data.frame(
    experiment = c(
      rep("mechanism", length(scenarios)),
      rep("measurement_noise", length(tau2_B_values)),
      rep("gamma_strength", length(gamma_values))
    ),
    condition = c(
      scenarios,
      paste0("tau2_B_", tau2_B_values),
      paste0("gamma_", gamma_values)
    ),
    gamma_true = c(
      ifelse(scenarios == "spatially_varying", NA_real_,
             info$generating$gamma),
      rep(info$generating$gamma, length(tau2_B_values)),
      gamma_values
    ),
    tau2_B_true = c(
      rep(info$generating$tau2_B, length(scenarios)),
      tau2_B_values,
      rep(info$generating$tau2_B, length(gamma_values))
    )
  )

  runs = lapply(seq_len(nrow(design)), function(i) {
    condition = design[i, , drop = FALSE]
    if (condition$experiment == "mechanism") {
      directory = file.path(extension_dir, condition$condition)
    } else {
      directory = file.path(
        extension_dir, condition$experiment, condition$condition
      )
    }
    message("  Reading ", condition$experiment, " / ", condition$condition)
    read_saved_condition(
      directory = directory,
      condition = condition,
      replications = replications,
      model_names = info$models,
      tau2_s_true = info$generating$tau2_s,
      gamma_true = condition$gamma_true,
      tau2_B_true = condition$tau2_B_true,
      mu_true = 0,
      separate_noise = isTRUE(info$separate_noise),
      n_target = info$n_replications
    )
  })

  metrics = bind_rows(lapply(runs, function(x) x$metrics))
  log = bind_rows(lapply(runs, function(x) x$log))
  issues = bind_rows(lapply(runs, function(x) x$issues))
  completion = summarize_completion(
    bind_rows(lapply(runs, function(x) x$completion)),
    group_vars = c("experiment", "condition")
  )

  ready = nrow(metrics) > 0L
  summary = data.frame()
  paired_metrics = data.frame()
  paired_summary = data.frame()

  if (ready) {
    summary = summarize_model_metrics(metrics, c("experiment", "condition"))
    timing = summarize_fit_times(log, c("experiment", "condition"))
    if (nrow(timing)) {
      summary = merge(
        summary, timing,
        by = c("experiment", "condition", "model"),
        all.x = TRUE, sort = FALSE
      )
    }
    paired = summarize_paired_metrics(metrics, c("experiment", "condition"))
    paired_metrics = paired$paired_metrics
    paired_summary = paired$paired_summary

    condition_truth = unique(design[c(
      "experiment", "condition", "gamma_true", "tau2_B_true"
    )])
    summary = merge(
      summary, condition_truth,
      by = c("experiment", "condition"), all.x = TRUE, sort = FALSE
    )
    paired_summary = merge(
      paired_summary, condition_truth,
      by = c("experiment", "condition"), all.x = TRUE, sort = FALSE
    )
  }

  save_metrics(
    extension_summary_dir, completion, issues, log,
    metrics, summary, paired_metrics, paired_summary
  )
  tables = make_extension_tables(info, completion, summary, paired_summary)
  if (nrow(issues)) tables$issues_preview = utils::head(issues, 10)

  report_index$extensions$found = TRUE
  report_index$extensions$ready = ready
  report_index$extensions$tables = save_tables(tables, extension_summary_dir)
  report_index$extensions$n_complete = sum(completion$n_complete)
  report_index$extensions$n_assigned = sum(completion$n_assigned)
  report_index$extensions$n_conditions = nrow(completion)
  report_index$extensions$n_fully_complete = sum(
    completion$n_complete == completion$n_target
  )
  report_index$extensions$n_issues = nrow(issues)

  if (ready) {
    for (experiment_name in c(
      "gamma_strength", "measurement_noise", "mechanism"
    )) {
      d = subset(
        paired_summary,
        experiment == experiment_name & is.finite(mean_delta_mse)
      )
      if (!nrow(d)) next
      experiment_label = switch(
        experiment_name,
        gamma_strength = "preferentiality-strength",
        measurement_noise = "areal-noise",
        mechanism = "sampling-mechanism"
      )
      report_index$extensions$findings = c(
        report_index$extensions$findings,
        sprintf(
          "Across the %s extensions, mean paired MSE gains range from %.3f to %.3f.",
          experiment_label, min(d$mean_delta_mse), max(d$mean_delta_mse)
        )
      )
    }
  }
} else if (update_tables) {
  message("Extension metadata not found: ", extension_info_path)
}


# Save summary tables before optional figure creation -------------------------

if (update_tables) {
  saveRDS(report_index, index_path)
  message("Summaries saved to ", summary_dir)
}


# Recover range and variance parameters using the same included replications.
if (update_parameter_tables) {
  for (section in c("reproduction", "extensions")) {
    if (!isTRUE(report_index[[section]]$ready)) next
    refresh_parameter_recovery(
      if (section == "reproduction") reproduction_dir else extension_dir,
      file.path(summary_dir, section), section, workers = parameter_workers
    )
  }
}


# Refresh reproduction fitting times from cached logs -------------------------

if ((update_tables || update_reproduction_tables || update_report_tables) &&
    isTRUE(report_index$reproduction$ready)) {
  reproduction_summary = read_summary_csv(reproduction_summary_dir, "model_summary")
  reproduction_log = read_summary_csv(reproduction_summary_dir, "fit_log")
  if (nrow(reproduction_summary) && nrow(reproduction_log)) {
    timing_keys = c("rho", "gamma", "I", "J", "model")
    timing = summarize_fit_times(reproduction_log, c("rho", "gamma", "I", "J"))
    timing_columns = setdiff(names(timing), timing_keys)
    reproduction_summary = merge(
      reproduction_summary[setdiff(names(reproduction_summary), timing_columns)],
      timing, by = timing_keys, all.x = TRUE, sort = FALSE
    )
    utils::write.csv(
      reproduction_summary,
      file.path(reproduction_summary_dir, "model_summary.csv"), row.names = FALSE
    )
    focus = report_index$reproduction$focus
    if (!is.null(focus)) {
      selected_timing = merge(
        reproduction_summary, as.data.frame(focus),
        by = c("rho", "gamma", "I", "J")
      )
      report_index = save_report_table(
        report_index, summary_dir, "reproduction", "timing",
        format_fitting_time_table(selected_timing)
      )
    }
    saveRDS(report_index, index_path)
  }
}


# Reproduction accuracy and Zhong comparison tables --------------------------

# Refresh these tables from the saved model summaries without reading fits
if (update_reproduction_tables &&
    isTRUE(report_index$reproduction$ready)) {
  reproduction_summary = read_summary_csv(
    reproduction_summary_dir, "model_summary"
  )

  required_columns = c(
    "rho", "gamma", "I", "J", "model", "n_included",
    "mean_mse", "mcse_mse", "mean_mae", "mcse_mae",
    "mean_wd2", "mcse_wd2"
  )

  if (nrow(reproduction_summary) &&
      all(required_columns %in% names(reproduction_summary))) {

    # Order scenarios, then models within each scenario
    model_order = c("ps_melding", "melding", "ps_geo")
    reproduction_summary = reproduction_summary[order(
      reproduction_summary$rho,
      reproduction_summary$gamma,
      reproduction_summary$I,
      reproduction_summary$J,
      match(reproduction_summary$model, model_order)
    ), , drop = FALSE]

    # Full reproduction accuracy table in the layout of the reference table
    accuracy = data.frame(
      rho = reproduction_summary$rho,
      gamma = reproduction_summary$gamma,
      I = reproduction_summary$I,
      J = reproduction_summary$J,
      Model = unname(summary_model_labels()[reproduction_summary$model]),
      n = reproduction_summary$n_included,
      MSE = format_mc_interval(
        reproduction_summary$mean_mse, reproduction_summary$mcse_mse
      ),
      MAE = format_mc_interval(
        reproduction_summary$mean_mae, reproduction_summary$mcse_mae
      ),
      "W2^2" = format_mc_interval(
        reproduction_summary$mean_wd2, reproduction_summary$mcse_wd2
      ),
      check.names = FALSE
    )

    # Published Table 2 MSE for rho = 0.2, gamma = 1, I = 100
    # PS Geo is a point-only model, so the J = 0 reference repeats across J
    zhong_mse = subset(zhong_prediction_reference(), rho == 0.2)
    names(zhong_mse)[names(zhong_mse) == "mean_mse"] = "zhong_mean_mse"

    scenario_columns = c("rho", "gamma", "I", "J", "model")
    reproduction_key = do.call(paste, c(
      reproduction_summary[scenario_columns], sep = "_"
    ))
    zhong_key = do.call(paste, c(zhong_mse[scenario_columns], sep = "_"))
    published_mse = zhong_mse$zhong_mean_mse[
      match(reproduction_key, zhong_key)
    ]

    # Keep numeric columns for future plots and calculations
    mse_by_J_comparison = data.frame(
      rho = reproduction_summary$rho,
      gamma = reproduction_summary$gamma,
      I = reproduction_summary$I,
      J = reproduction_summary$J,
      model = reproduction_summary$model,
      n = reproduction_summary$n_included,
      reproduction_mean_mse = reproduction_summary$mean_mse,
      reproduction_mcse = reproduction_summary$mcse_mse,
      reproduction_ci_lower = reproduction_summary$mean_mse -
        1.96 * reproduction_summary$mcse_mse,
      reproduction_ci_upper = reproduction_summary$mean_mse +
        1.96 * reproduction_summary$mcse_mse,
      zhong_mean_mse = published_mse,
      difference_from_zhong = reproduction_summary$mean_mse - published_mse
    )

    table_dir = file.path(reproduction_summary_dir, "tables")
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

    utils::write.csv(
      accuracy, file.path(table_dir, "accuracy.csv"), row.names = FALSE
    )
    utils::write.csv(
      mse_by_J_comparison,
      file.path(table_dir, "mse_by_J_comparison.csv"),
      row.names = FALSE
    )

    report_index$reproduction$tables = unique(c(
      report_index$reproduction$tables,
      "accuracy", "mse_by_J_comparison"
    ))
    saveRDS(report_index, index_path)

    message("Updated all-condition reproduction accuracy and MSE tables")
  } else {
    warning("Missing required reproduction model_summary columns for tables")
  }
}


# Reproduction paired comparisons from cached replication-level scores --------

if ((update_tables || update_reproduction_tables || update_report_tables) &&
    isTRUE(report_index$reproduction$ready)) {
  paired_replications = read_summary_csv(
    reproduction_summary_dir, "paired_replications"
  )
  condition_columns = c("rho", "gamma", "I", "J")
  pairwise_gains = summarize_pairwise_gains(
    paired_replications, condition_columns
  )
  if (nrow(pairwise_gains)) {
    utils::write.csv(
      pairwise_gains,
      file.path(reproduction_summary_dir, "pairwise_gain_summary.csv"),
      row.names = FALSE
    )
    focus = report_index$reproduction$focus
    if (!is.null(focus)) {
      focused_gains = merge(
        pairwise_gains, as.data.frame(focus), by = condition_columns
      )
      report_index = save_report_table(
        report_index, summary_dir, "reproduction", "paired_comparisons",
        format_pairwise_gain_table(focused_gains)
      )
    }
    for (baseline in c("melding", "ps_geo")) {
      mse_gains = subset(pairwise_gains, comparator == baseline & metric == "mse")
      report_index = save_report_table(
        report_index, summary_dir, "reproduction",
        paste0("paired_mse_", baseline),
        format_pairwise_gain_table(mse_gains, condition_columns)
      )
    }
    report_index$reproduction$tables = setdiff(
      report_index$reproduction$tables,
      c("paired_gain", "paired_other", "paired_vs_geo", "all_conditions")
    )
    saveRDS(report_index, index_path)
    message("Updated reproduction paired summaries for Melding and PSgeo")
  }
}


# Reproduction plots ----------------------------------------------------------

# Replication-level accuracy figures are also part of the main report.
if ((create_plots || create_report_plots) &&
    isTRUE(make_plots$reproduction_performance) &&
    isTRUE(report_index$reproduction$ready)) {
  reproduction_metrics = read_summary_csv(
    reproduction_summary_dir, "replication_metrics"
  )
  p_reproduction_mse = plot_reproduction_performance(
    reproduction_metrics, gamma = 1, I = 100,
    plot_metrics = c("MSE" = "mse"),
    benchmark_reference = zhong_prediction_reference()
  )
  p_reproduction_mae_wd = plot_reproduction_performance(
    reproduction_metrics, gamma = 1, I = 100,
    plot_metrics = c("MAE" = "mae", "WD (Zhong code)" = "wd_zhong"),
    benchmark_reference = zhong_prediction_reference()
  )
  report_index = save_report_plot(
    report_index, summary_dir, "reproduction", "prediction_mse",
    p_reproduction_mse, width = 10, height = 4.8
  )
  report_index = save_report_plot(
    report_index, summary_dir, "reproduction", "prediction_mae_wd",
    p_reproduction_mae_wd, width = 10, height = 6.5
  )
  if (!is.null(p_reproduction_mse) && !is.null(p_reproduction_mae_wd)) {
    report_index$reproduction$figures = setdiff(
      report_index$reproduction$figures,
      c("prediction_performance", "prediction_mae_w2")
    )
  }
}


if (create_plots && isTRUE(report_index$reproduction$ready)) {
  reproduction_metrics = read_summary_csv(
    reproduction_summary_dir, "replication_metrics"
  )
  reproduction_summary = read_summary_csv(
    reproduction_summary_dir, "model_summary"
  )
  reproduction_paired = read_summary_csv(
    reproduction_summary_dir, "paired_summary"
  )

  # Mean MSE by J compared with Zhong et al Table 2
  if (isTRUE(make_plots$reproduction_mean_mse_comparison) &&
      nrow(reproduction_summary)) {
    reproduction_mse = subset(
      reproduction_summary,
      rho == 0.2 & gamma == 1 & I == 100 &
        J %in% c(4, 25, 100) &
        model %in% c("ps_melding", "melding", "ps_geo") &
        is.finite(mean_mse)
    )

    if (nrow(reproduction_mse)) {
      if (any(reproduction_mse$n_included != 100L)) {
        warning("Some plotted reproduction means use fewer than 100 replications")
      }

      reproduction_mse = reproduction_mse[
        , c("J", "model", "mean_mse"), drop = FALSE
      ]
      reproduction_mse$source = "Reproduction"

      # Published preferential-sampling means for rho = 0.2 and I = 100
      # PS Geo has no areal observations so its J = 0 value is repeated
      zhong_mse = subset(
        zhong_prediction_reference(), rho == 0.2,
        select = c("J", "model", "mean_mse")
      )
      zhong_mse$source = "Zhong et al. (2025)"

      mse_comparison = rbind(reproduction_mse, zhong_mse)
      mse_comparison$J = factor(mse_comparison$J, levels = c(4, 25, 100))
      mse_comparison$model = factor(
        mse_comparison$model, levels = names(summary_model_labels())
      )
      mse_comparison$source = factor(
        mse_comparison$source,
        levels = c("Reproduction", "Zhong et al. (2025)")
      )

      reproduction_line = subset(
        mse_comparison, source == "Reproduction"
      )
      zhong_line = subset(
        mse_comparison, source == "Zhong et al. (2025)"
      )

      p_reproduction_mse_comparison = ggplot2::ggplot() +
        ggplot2::geom_line(
          data = reproduction_line,
          ggplot2::aes(
            x = J, y = mean_mse, colour = model, group = model
          ),
          linewidth = 0.9
        ) +
        ggplot2::geom_point(
          data = reproduction_line,
          ggplot2::aes(
            x = J, y = mean_mse, colour = model
          ),
          shape = 16, size = 2.6
        ) +
        ggplot2::geom_line(
          data = zhong_line,
          ggplot2::aes(
            x = J, y = mean_mse, colour = model, group = model
          ),
          linetype = "dashed", linewidth = 0.7,
          alpha = 0.85, show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data = zhong_line,
          ggplot2::aes(
            x = J, y = mean_mse, colour = model
          ),
          shape = 8, size = 3.1, stroke = 0.85,
          show.legend = FALSE
        ) +
        ggplot2::scale_colour_manual(
          values = summary_model_colours(),
          breaks = names(summary_model_labels()),
          labels = summary_model_labels()
        ) +
        ggplot2::labs(
          subtitle = paste(
            "Solid lines and circles: reproduction",
            "| Dashed lines and stars: Zhong et al. (2025)"
          ),
          x = "Number of areal observations (J)",
          y = "Mean MSE",
          colour = NULL
        ) +
        ggplot2::guides(
          colour = ggplot2::guide_legend(
            nrow = 1,
            override.aes = list(linetype = "solid", shape = 16)
          )
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_line(
            colour = "grey90", linewidth = 0.3
          ),
          legend.position = "top",
          legend.justification = "left",
          plot.subtitle = ggplot2::element_text(
            size = 10, colour = "grey35",
            margin = ggplot2::margin(b = 8)
          )
        )

      if (save_plot(
        p_reproduction_mse_comparison, reproduction_summary_dir,
        "mse_by_J_comparison", width = 8, height = 5
      )) {
        report_index$reproduction$figures = unique(c(
          report_index$reproduction$figures, "mse_by_J_comparison"
        ))
      }
    }
  }

  # Paired MSE improvement from PS Melding over Melding
  if (isTRUE(make_plots$reproduction_paired_gain)) {
    p_reproduction_gain = plot_reproduction_mse_gain(reproduction_paired)
    if (save_plot(
      p_reproduction_gain, reproduction_summary_dir,
      "paired_mse_gain", width = 8, height = 5
    )) {
      report_index$reproduction$figures = unique(c(
        report_index$reproduction$figures, "paired_mse_gain"
      ))
    }
  }

  # Example latent surface, areal observations and SPDE mesh
  make_reproduction_examples = isTRUE(make_plots$reproduction_simulated_data) ||
    isTRUE(make_plots$reproduction_mesh)

  if (make_reproduction_examples && !is.null(report_index$reproduction$focus)) {
    info_path = file.path(reproduction_dir, "run_info.rds")
    focus = report_index$reproduction$focus
    if (file.exists(info_path)) {
      info = readRDS(info_path)
      replications = if (!is.null(info$assigned_replications)) {
        info$assigned_replications
      } else {
        seq_len(info$n_replications)
      }
      focus_dir = file.path(
        reproduction_dir,
        paste0("rho_", focus$rho),
        paste0("gamma_", focus$gamma),
        paste0("I_", focus$I, "_J_", focus$J)
      )
      example_paths = file.path(
        focus_dir, "fits", sprintf("rep_%03d.rds", replications)
      )
      example_paths = example_paths[file.exists(example_paths)]

      if (length(example_paths)) {
        saved = tryCatch(readRDS(example_paths[1]), error = function(e) NULL)
        if (is.list(saved) && !is.null(saved$data$latent)) {
          dataset = saved$data

          if (isTRUE(make_plots$reproduction_simulated_data)) {
            p_reproduction_data = plot_simulated_data(dataset)
            if (save_plot(
              p_reproduction_data, reproduction_summary_dir,
              "simulated_data", width = 12, height = 4
            )) {
              report_index$reproduction$figures = unique(c(
                report_index$reproduction$figures, "simulated_data"
              ))
            }
          }

          if (isTRUE(make_plots$reproduction_mesh)) {
            boundary = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))
            mesh = INLA::inla.mesh.2d(
              loc.domain = boundary,
              max.edge = info$mesh_max_edge,
              offset = info$mesh_offset,
              cutoff = info$mesh_cutoff
            )
            p_reproduction_mesh = plot_reproduction_mesh(mesh, dataset)
            if (save_plot(
              p_reproduction_mesh, reproduction_summary_dir,
              "mesh_areas_points", width = 6, height = 6
            )) {
              report_index$reproduction$figures = unique(c(
                report_index$reproduction$figures, "mesh_areas_points"
              ))
            }
          }
        }
      }
    }
  }
}


# Refresh extension display tables from the current cached extension results.
if (update_report_tables && isTRUE(report_index$extensions$ready)) {
  tables = make_extension_tables(
    readRDS(file.path(extension_dir, "run_info.rds")),
    read_summary_csv(extension_summary_dir, "completion"),
    read_summary_csv(extension_summary_dir, "model_summary"),
    read_summary_csv(extension_summary_dir, "paired_summary")
  )
  report_index$extensions$tables = unique(c(
    report_index$extensions$tables, save_tables(tables, extension_summary_dir)
  ))
}


# Extension plots -------------------------------------------------------------

# Predictive accuracy and calibration with the constant-mechanism reference.
if ((create_plots || create_report_plots) &&
    isTRUE(make_plots$mechanism_predictive_performance) &&
    isTRUE(report_index$extensions$ready)) {
  extension_metrics = read_summary_csv(extension_summary_dir, "replication_metrics")
  p_mechanism_performance = plot_mechanism_predictive_performance(extension_metrics)
  report_index = save_report_plot(
    report_index, summary_dir, "extensions", "mechanism_predictive_performance",
    p_mechanism_performance, width = 10, height = 7.2
  )
  if (!is.null(p_mechanism_performance)) {
    report_index$extensions$figures = setdiff(
      report_index$extensions$figures, "mechanism_prediction_mse"
    )
  }
}


if ((create_plots || create_report_plots) && isTRUE(report_index$extensions$ready)) {
  extension_summary = read_summary_csv(
    extension_summary_dir, "model_summary"
  )
  extension_paired = read_summary_csv(
    extension_summary_dir, "paired_summary"
  )

  if (nrow(extension_summary)) {
    gamma_summary = subset(extension_summary, experiment == "gamma_strength")

    # Mean prediction MSE across preferentiality strengths
    if (isTRUE(make_plots$gamma_prediction_mse)) {
      p_gamma_mse = plot_gamma_prediction_mse(gamma_summary)
      if (save_plot(p_gamma_mse, extension_summary_dir, "gamma_prediction_mse")) {
        report_index$extensions$figures = unique(c(
          report_index$extensions$figures, "gamma_prediction_mse"
        ))
      }
    }

    # Mean prediction MSE across areal-noise levels, including the baseline.
    if (isTRUE(make_plots$noise_prediction_mse)) {
      noise_summary = subset(extension_summary, experiment == "measurement_noise")
      report_index = save_report_plot(
        report_index, summary_dir, "extensions", "noise_prediction_mse",
        plot_noise_prediction_mse(noise_summary), width = 9, height = 4.7
      )
    }

    # Recovery of the preferentiality parameter
    if (isTRUE(make_plots$gamma_recovery)) {
      p_gamma_recovery = plot_gamma_parameter_recovery(gamma_summary)
      if (save_plot(p_gamma_recovery, extension_summary_dir, "gamma_recovery")) {
        report_index$extensions$figures = unique(c(
          report_index$extensions$figures, "gamma_recovery"
        ))
      }
    }


  }

  if (nrow(extension_paired)) {
    gamma_paired = subset(extension_paired, experiment == "gamma_strength")
    noise_paired = subset(extension_paired, experiment == "measurement_noise")

    # Paired MSE gain across preferentiality strengths
    if (isTRUE(make_plots$gamma_paired_gain)) {
      p_gamma_gain = plot_gamma_mse_gain(gamma_paired)
      if (save_plot(p_gamma_gain, extension_summary_dir, "gamma_paired_gain")) {
        report_index$extensions$figures = unique(c(
          report_index$extensions$figures, "gamma_paired_gain"
        ))
      }
    }

    # Paired MSE gain across areal measurement-error variances
    if (isTRUE(make_plots$noise_paired_gain)) {
      p_noise_gain = plot_noise_mse_gain(noise_paired)
      if (save_plot(p_noise_gain, extension_summary_dir, "noise_paired_gain")) {
        report_index$extensions$figures = unique(c(
          report_index$extensions$figures, "noise_paired_gain"
        ))
      }
    }
  }
}


# Extension illustrations -----------------------------------------------------

make_extension_examples = isTRUE(make_plots$extension_sampling_scenarios) ||
  isTRUE(make_plots$extension_mesh)

if (create_plots && isTRUE(report_index$extensions$found) &&
    make_extension_examples) {
  source("functions/generate_data.R")
  source("functions/fit_models.R")

  scenarios = c(
    "Constant gamma" = "constant",
    "Spatially varying gamma" = "spatially_varying",
    "Independent latent field" = "independent_field"
  )
  simulations = lapply(scenarios, function(scenario) {
    generate_data(
      rho = 0.2, gamma = 1, I = 100, J = 25,
      tau2_s = 0.1, tau2_B = 0.1, seed = 123,
      scenario = scenario, delta = 1, beta_psi = 1, rho_psi = 0.2
    )
  })

  # Point observations across the three sampling scenarios
  if (isTRUE(make_plots$extension_sampling_scenarios)) {
    field = simulations[[1]]$latent
    log_lambda = lapply(simulations, function(d) d$latent$log_lambda)
    points = lapply(simulations, function(d) d$points)
    p_sampling_scenarios = plot_sampling_scenarios(
      field = field, scenarios = log_lambda, points = points,
      limits = c(-5, 5)
    )
    if (save_plot(
      p_sampling_scenarios, extension_summary_dir,
      "sampling_scenarios", width = 12, height = 4
    )) {
      report_index$extensions$figures = unique(c(
        report_index$extensions$figures, "sampling_scenarios"
      ))
    }
  }

  # Mesh, areal units and point observations across sampling scenarios
  if (isTRUE(make_plots$extension_mesh)) {
    extension_info_path = file.path(extension_dir, "run_info.rds")
    if (file.exists(extension_info_path)) {
      extension_info = readRDS(extension_info_path)
      mesh = create_mesh(max_edge = extension_info$mesh_max_edge)
      p_extension_mesh = plot_extension_mesh(
        mesh, areas = simulations[["Constant gamma"]]$areas,
        simulations = simulations
      )
      if (save_plot(
        p_extension_mesh, extension_summary_dir,
        "mesh_areas_points", width = 12, height = 4
      )) {
        report_index$extensions$figures = unique(c(
          report_index$extensions$figures, "mesh_areas_points"
        ))
      }
    }
  }
}


# Save report index -----------------------------------------------------------

if (create_plots) {
  saveRDS(report_index, index_path)
  message("Plots saved to ", summary_dir)
}

# Comprehensive and main-report results ---------------------------------------

# These summaries use the cached CSV files and do not load fitted INLA objects
if (update_report_tables || create_report_plots) {

  # Preserve all available measures and actual per-model replication counts
  model_columns = c(
    "rho", "gamma", "I", "J", "experiment", "condition",
    "gamma_true", "tau2_B_true", "model", "n_included",
    "mean_mse", "mcse_mse", "mse_ci_low", "mse_ci_high",
    "mean_mae", "mcse_mae", "mean_wd2", "mcse_wd2",
    "mean_coverage", "mean_interval_width", "mean_interval_score",
    "mean_mu", "mu_bias", "mu_coverage", "mean_mu_interval_score",
    "median_gamma", "gamma_bias", "gamma_coverage",
    "mean_gamma_interval_score", "median_tau2", "median_tau2_s",
    "median_tau2_B", "tau2_B_bias", "tau2_B_coverage",
    "median_elapsed_sec", "n_attempted", "n_reused", "n_failed"
  )
  paired_columns = c(
    "rho", "gamma", "I", "J", "experiment", "condition",
    "gamma_true", "tau2_B_true", "n_paired", "mean_delta_mse",
    "mcse_delta_mse", "ci_low", "ci_high", "proportion_ps_better",
    "mean_relative_mse_reduction", "mean_delta_mae", "mean_delta_wd2",
    "mean_delta_interval_score", "n_paired_vs_geo",
    "mean_delta_mse_vs_geo"
  )

  report_models = list()
  report_pairs = list()
  for (section in c("reproduction", "extensions")) {
    if (!isTRUE(report_index[[section]]$ready)) next
    destination = file.path(summary_dir, section)
    model_data = read_summary_csv(destination, "model_summary")
    paired_data = read_summary_csv(destination, "paired_summary")

    if (nrow(model_data) && all(c(
      "mean_mse", "mcse_mse", "model", "n_included"
    ) %in% names(model_data))) {
      model_data$mse_ci_low = model_data$mean_mse -
        1.96 * model_data$mcse_mse
      model_data$mse_ci_high = model_data$mean_mse +
        1.96 * model_data$mcse_mse
      report_models[[section]] = model_data
      if (update_report_tables) {
        ordered_columns = c(
          intersect(model_columns, names(model_data)),
          setdiff(names(model_data), model_columns)
        )
        report_index = save_report_table(
          report_index, summary_dir,
          section, "report_all_model_results",
          model_data[ordered_columns]
        )
      }
    }

    if (nrow(paired_data) && all(c(
      "n_paired", "mean_delta_mse", "ci_low", "ci_high"
    ) %in% names(paired_data))) {
      report_pairs[[section]] = paired_data
      if (update_report_tables) {
        ordered_columns = c(
          intersect(paired_columns, names(paired_data)),
          setdiff(names(paired_data), paired_columns)
        )
        report_index = save_report_table(
          report_index, summary_dir,
          section, "report_all_paired_results",
          paired_data[ordered_columns]
        )
      }
    }
  }

  model_labels = summary_model_labels()
  model_colours = summary_model_colours()
  main_theme = ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )

  # Reproduction: summary of every spatial range and preferentiality setting
  reproduction_pairs = report_pairs$reproduction
  if (!is.null(reproduction_pairs) && nrow(reproduction_pairs) &&
      all(c("rho", "gamma", "I", "J") %in% names(reproduction_pairs))) {
    valid_pairs = reproduction_pairs[
      is.finite(reproduction_pairs$mean_delta_mse) &
        reproduction_pairs$n_paired > 0, , drop = FALSE
    ]
    grouped = split(valid_pairs, interaction(
      valid_pairs$rho, valid_pairs$gamma, valid_pairs$I, drop = TRUE
    ))
    reproduction_overview = bind_rows(lapply(grouped, function(d) {
      data.frame(
        rho = d$rho[1], gamma = d$gamma[1], I = d$I[1],
        n_J = length(unique(d$J)),
        n_positive_gain = sum(d$mean_delta_mse > 0),
        n_mc_ci_above_zero = sum(d$ci_low > 0, na.rm = TRUE),
        min_gain = min(d$mean_delta_mse),
        max_gain = max(d$mean_delta_mse),
        min_n_paired = min(d$n_paired),
        max_n_paired = max(d$n_paired)
      )
    }))
    if (update_report_tables) {
      report_index = save_report_table(
        report_index, summary_dir,
        "reproduction", "report_condition_overview", reproduction_overview
      )
    }

    if (create_report_plots && nrow(valid_pairs)) {
      gain_data = valid_pairs
      gain_data$J = factor(
        gain_data$J, levels = sort(unique(gain_data$J))
      )
      gain_data$Gamma = factor(
        paste0("gamma = ", gain_data$gamma)
      )
      gain_data$rho_group = factor(
        gain_data$rho, levels = sort(unique(gain_data$rho))
      )
      rho_colours = setNames(
        grDevices::hcl.colors(nlevels(gain_data$rho_group), "Dark 3"),
        levels(gain_data$rho_group)
      )
      gain_dodge = ggplot2::position_dodge(width = 0.55)

      p_reproduction_all_gains = ggplot2::ggplot(
        gain_data, ggplot2::aes(
          x = J, y = mean_delta_mse,
          colour = rho_group, group = rho_group
        )
      ) +
        ggplot2::geom_hline(
          yintercept = 0, colour = "grey55", linetype = "dashed"
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = ci_low, ymax = ci_high),
          position = gain_dodge, width = 0.14, linewidth = 0.5
        ) +
        ggplot2::geom_point(position = gain_dodge, size = 2.5) +
        ggplot2::facet_grid(I ~ Gamma) +
        ggplot2::scale_colour_manual(values = rho_colours) +
        ggplot2::labs(
          x = "Number of areal observations (J)",
          y = "Paired MSE gain: Melding - PS Melding",
          colour = "Spatial range (rho)",
          subtitle = "All available reproduction conditions; bars are 95% Monte Carlo intervals"
        ) + main_theme
      report_index = save_report_plot(
        report_index, summary_dir,
        "reproduction", "report_all_paired_gains",
        p_reproduction_all_gains, 10, 5.2
      )
    }
  }

  # Only the exact published Table 2 conditions have Zhong MSE references
  zhong_path = file.path(
    reproduction_summary_dir, "tables", "mse_by_J_comparison.csv"
  )
  if (isTRUE(report_index$reproduction$ready) && file.exists(zhong_path)) {
    zhong = utils::read.csv(zhong_path, stringsAsFactors = FALSE)
    needed = c(
      "rho", "gamma", "I", "J", "model", "n",
      "reproduction_mean_mse", "reproduction_ci_lower",
      "reproduction_ci_upper", "zhong_mean_mse"
    )
    if (all(needed %in% names(zhong))) {
      zhong = subset(
        zhong,
        rho == 0.2 & gamma == 1 & I == 100 &
          J %in% c(4, 25, 100) & is.finite(zhong_mean_mse) &
          is.finite(reproduction_mean_mse)
      )
    } else {
      zhong = data.frame()
    }

    if (nrow(zhong)) {
      zhong$delta = zhong$reproduction_mean_mse - zhong$zhong_mean_mse
      zhong$delta_low = zhong$reproduction_ci_lower - zhong$zhong_mean_mse
      zhong$delta_high = zhong$reproduction_ci_upper - zhong$zhong_mean_mse
      zhong$percent_difference = 100 * zhong$delta / zhong$zhong_mean_mse
      zhong$Model = factor(
        unname(model_labels[zhong$model]), levels = unname(model_labels)
      )
      zhong = zhong[order(zhong$Model, zhong$J), , drop = FALSE]

      if (update_report_tables) {
        zhong_table = data.frame(
          Model = as.character(zhong$Model), J = zhong$J, n = zhong$n,
          `Reproduction MSE (95% MC CI)` = sprintf(
            "%.3f (%.3f, %.3f)", zhong$reproduction_mean_mse,
            zhong$reproduction_ci_lower, zhong$reproduction_ci_upper
          ),
          `Zhong MSE` = zhong$zhong_mean_mse,
          `Difference` = zhong$delta,
          `Difference (%)` = zhong$percent_difference,
          check.names = FALSE
        )
        report_index = save_report_table(
          report_index, summary_dir,
          "reproduction", "report_zhong_comparison", zhong_table
        )
      }

      if (create_report_plots) {
        zhong$J = factor(zhong$J, levels = c(4, 25, 100))
        p_agreement = ggplot2::ggplot(
          zhong, ggplot2::aes(x = J, colour = Model, group = Model)
        ) +
          ggplot2::geom_line(
            ggplot2::aes(y = reproduction_mean_mse), linewidth = 0.75
          ) +
          ggplot2::geom_errorbar(
            ggplot2::aes(
              ymin = reproduction_ci_lower, ymax = reproduction_ci_upper
            ), width = 0.1, linewidth = 0.4
          ) +
          ggplot2::geom_point(
            ggplot2::aes(y = reproduction_mean_mse), size = 2.4
          ) +
          ggplot2::geom_line(
            ggplot2::aes(y = zhong_mean_mse), linetype = "dashed",
            linewidth = 0.6
          ) +
          ggplot2::geom_point(
            ggplot2::aes(y = zhong_mean_mse), shape = 8, size = 3
          ) +
          ggplot2::facet_wrap(~ Model, nrow = 1) +
          ggplot2::scale_colour_manual(
            values = setNames(model_colours, unname(model_labels))
          ) +
          ggplot2::labs(
            x = "Number of areal observations (J)",
            y = "Mean MSE",
            subtitle = paste(
              "Solid circles with 95% Monte Carlo intervals: reproduction",
              "| Dashed lines and stars: Zhong et al."
            )
          ) + main_theme + ggplot2::theme(legend.position = "none")
        report_index = save_report_plot(
          report_index, summary_dir,
          "reproduction", "report_mse_agreement", p_agreement, 11, 4.3
        )

        p_difference = ggplot2::ggplot(
          zhong, ggplot2::aes(x = J, y = delta, colour = Model)
        ) +
          ggplot2::geom_hline(
            yintercept = 0, linetype = "dashed", colour = "grey55"
          ) +
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = delta_low, ymax = delta_high),
            width = 0.12, linewidth = 0.55
          ) +
          ggplot2::geom_point(size = 2.6) +
          ggplot2::facet_wrap(~ Model, nrow = 1) +
          ggplot2::scale_colour_manual(
            values = setNames(model_colours, unname(model_labels))
          ) +
          ggplot2::labs(
            x = "Number of areal observations (J)",
            y = "Reproduction mean MSE - Zhong mean MSE",
            subtitle = paste(
              "Zero denotes agreement",
              "| Intervals reflect reproduction Monte Carlo error only"
            )
          ) + main_theme + ggplot2::theme(legend.position = "none")
        report_index = save_report_plot(
          report_index, summary_dir,
          "reproduction", "report_mse_difference", p_difference, 11, 4.3
        )
      }
    }
  }

  # Reproduction parameter-recovery table for the central published setting
  reproduction_models = report_models$reproduction
  if (update_report_tables && !is.null(reproduction_models) &&
      all(c("rho", "gamma", "I", "J") %in% names(reproduction_models))) {
    parameters = subset(
      reproduction_models, rho == 0.2 & I == 100 & J == 25 &
        gamma %in% c(0, 1)
    )
    parameter_columns = c(
      "gamma", "model", "n_included", "mu_bias", "mu_coverage",
      "mean_mu_interval_score", "median_gamma", "gamma_bias",
      "gamma_coverage", "mean_gamma_interval_score", "median_tau2"
    )
    report_index = save_report_table(
      report_index, summary_dir,
      "reproduction", "report_parameter_recovery",
      parameters[intersect(parameter_columns, names(parameters))]
    )
  }

  # Extensions: condition-level gain patterns across the three experiments
  extension_pairs = report_pairs$extensions
  if (!is.null(extension_pairs) && nrow(extension_pairs) &&
      all(c("experiment", "condition") %in% names(extension_pairs))) {
    valid_pairs = extension_pairs[
      is.finite(extension_pairs$mean_delta_mse) &
        extension_pairs$n_paired > 0, , drop = FALSE
    ]
    grouped = split(valid_pairs, valid_pairs$experiment)
    extension_overview = bind_rows(lapply(grouped, function(d) {
      data.frame(
        experiment = d$experiment[1], n_conditions = nrow(d),
        n_positive_gain = sum(d$mean_delta_mse > 0),
        n_mc_ci_above_zero = sum(d$ci_low > 0, na.rm = TRUE),
        min_gain = min(d$mean_delta_mse),
        max_gain = max(d$mean_delta_mse),
        min_n_paired = min(d$n_paired),
        max_n_paired = max(d$n_paired)
      )
    }))
    if (update_report_tables) {
      report_index = save_report_table(
        report_index, summary_dir,
        "extensions", "report_experiment_overview", extension_overview
      )
    }

    if (create_report_plots && nrow(valid_pairs)) {
      gain_data = valid_pairs
      gain_data$Experiment = factor(
        gain_data$experiment,
        levels = c("gamma_strength", "measurement_noise", "mechanism"),
        labels = c("Preferentiality strength", "Areal noise", "Mechanism")
      )
      gain_data = gain_data[order(
        gain_data$Experiment, gain_data$gamma_true,
        gain_data$tau2_B_true, gain_data$condition
      ), , drop = FALSE]
      gain_data$Condition = factor(
        gain_data$condition,
        levels = rev(unique(gain_data$condition))
      )

      p_extension_gains = ggplot2::ggplot(
        gain_data,
        ggplot2::aes(x = mean_delta_mse, y = Condition)
      ) +
        ggplot2::geom_vline(
          xintercept = 0, linetype = "dashed", colour = "grey55"
        ) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = ci_low, xend = ci_high,
            y = Condition, yend = Condition
          ), linewidth = 0.6, colour = "#2E7ED1"
        ) +
        ggplot2::geom_point(size = 2.2, colour = "#2E7ED1") +
        ggplot2::facet_grid(
          Experiment ~ ., scales = "free_y", space = "free_y"
        ) +
        ggplot2::labs(
          x = "Paired MSE gain: Melding - PS Melding",
          y = NULL,
          subtitle = "Positive values favor PS Melding; bars are 95% Monte Carlo intervals"
        ) + main_theme
      report_index = save_report_plot(
        report_index, summary_dir,
        "extensions", "report_all_paired_gains", p_extension_gains,
        10, min(12, max(6, 0.4 * nrow(gain_data) + 2))
      )
    }
  }

  # Extensions: mechanism accuracy, predictive coverage and gamma recovery
  extension_models = report_models$extensions
  if (!is.null(extension_models) &&
      all(c("experiment", "condition") %in% names(extension_models))) {
    mechanism = subset(
      extension_models,
      experiment == "mechanism" & is.finite(mean_mse)
    )
    gamma_models = subset(
      extension_models,
      experiment == "gamma_strength" &
        model %in% c("ps_melding", "ps_geo")
    )
    if (update_report_tables) {
      mechanism_columns = c(
        "condition", "model", "n_included", "mean_mse",
        "mse_ci_low", "mse_ci_high", "mean_mae", "mean_wd2",
        "mean_coverage", "mean_interval_score"
      )
      report_index = save_report_table(
        report_index, summary_dir,
        "extensions", "report_mechanism_accuracy",
        mechanism[intersect(mechanism_columns, names(mechanism))]
      )
      gamma_columns = c(
        "gamma_true", "model", "n_included", "median_gamma",
        "gamma_bias", "gamma_coverage", "mean_gamma_interval_score"
      )
      report_index = save_report_table(
        report_index, summary_dir,
        "extensions", "report_gamma_parameter_recovery",
        gamma_models[intersect(gamma_columns, names(gamma_models))]
      )
    }

    if (create_report_plots && nrow(mechanism)) {
      mechanism$Scenario = factor(
        mechanism$condition,
        levels = c("constant", "spatially_varying", "independent_field"),
        labels = c("Constant", "Spatially varying", "Independent field")
      )
      mechanism$Model = factor(
        unname(model_labels[mechanism$model]), levels = unname(model_labels)
      )
      dodge = ggplot2::position_dodge(width = 0.55)

      p_mechanism = ggplot2::ggplot(
        mechanism,
        ggplot2::aes(x = Scenario, y = mean_mse, colour = Model)
      ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = mse_ci_low, ymax = mse_ci_high),
          width = 0.12, position = dodge, linewidth = 0.5
        ) +
        ggplot2::geom_point(size = 2.8, position = dodge) +
        ggplot2::scale_colour_manual(
          values = setNames(model_colours, unname(model_labels))
        ) +
        ggplot2::labs(
          x = "Generating sampling mechanism", y = "Mean prediction MSE",
          colour = "Model",
          subtitle = "Error bars denote 95% Monte Carlo intervals"
        ) + main_theme
      report_index = save_report_plot(
        report_index, summary_dir,
        "extensions", "report_mechanism_mse", p_mechanism, 9, 4.8
      )

      if ("mean_coverage" %in% names(mechanism)) {
        coverage_data = mechanism[is.finite(mechanism$mean_coverage), ]
        if (nrow(coverage_data)) {
          p_coverage = ggplot2::ggplot(
            coverage_data,
            ggplot2::aes(x = Scenario, y = mean_coverage, colour = Model)
          ) +
            ggplot2::geom_hline(
              yintercept = 0.95, linetype = "dashed", colour = "grey55"
            ) +
            ggplot2::geom_point(size = 2.8, position = dodge) +
            ggplot2::scale_colour_manual(
              values = setNames(model_colours, unname(model_labels))
            ) +
            ggplot2::coord_cartesian(ylim = c(0, 1)) +
            ggplot2::labs(
              x = "Generating sampling mechanism",
              y = "Mean 95% predictive interval coverage",
              colour = "Model",
              subtitle = "Dashed line indicates nominal 95% coverage"
            ) + main_theme
          report_index = save_report_plot(
            report_index, summary_dir,
            "extensions", "report_mechanism_coverage", p_coverage, 9, 4.8
          )
        }
      }
    }
  }

  # Range and variance recovery: compact report tables and complete exports.
  for (section in c("reproduction", "extensions")) {
    destination = file.path(summary_dir, section)
    recovery = read_summary_csv(destination, "parameter_recovery")
    if (!nrow(recovery)) next
    groups = if (section == "reproduction") c("rho", "gamma", "I", "J") else {
      c("experiment", "condition")
    }
    if (update_report_tables) {
      report_index = save_report_table(
        report_index, summary_dir, section, "report_all_parameter_recovery",
        format_parameter_recovery(recovery, groups)
      )
      spatial = subset(recovery, parameter %in% c("sigma2", "theta"))
      spatial = if (section == "reproduction") {
        subset(spatial, gamma == 1 & I == 100 & J == 25)
      } else subset(spatial, experiment == "mechanism")
      report_index = save_report_table(
        report_index, summary_dir, section, "report_spatial_diagnostics",
        format_parameter_recovery(spatial, if (section == "reproduction") "rho" else "condition")
      )
      if (section == "reproduction") {
        focus_recovery = subset(recovery, rho == 0.2 & I == 100 & J == 25 &
                                  parameter %in% c("tau2", "rho"))
        report_index = save_report_table(
          report_index, summary_dir, section, "report_noise_range_recovery",
          format_parameter_recovery(focus_recovery, "gamma")
        )
        range_recovery = subset(recovery, gamma == 1 & I == 100 & J == 25 & parameter == "rho")
        report_index = save_report_table(
          report_index, summary_dir, section, "report_range_recovery",
          format_parameter_recovery(range_recovery)
        )
      } else {
        noise_recovery = subset(recovery, experiment == "measurement_noise" &
                                  parameter %in% c("tau2_s", "tau2_B"))
        report_index = save_report_table(
          report_index, summary_dir, section, "report_noise_recovery",
          format_parameter_recovery(noise_recovery, "condition")
        )
        baseline_recovery = subset(recovery, experiment == "mechanism" &
                                     condition == "constant" &
                                     parameter %in% c("tau2_s", "tau2_B", "rho"))
        report_index = save_report_table(
          report_index, summary_dir, section, "report_noise_range_recovery",
          format_parameter_recovery(baseline_recovery)
        )
        mechanism_range = subset(recovery, experiment == "mechanism" & parameter == "rho")
        report_index = save_report_table(
          report_index, summary_dir, section, "report_range_recovery",
          format_parameter_recovery(mechanism_range, "condition")
        )
      }
    }
    if (section == "extensions" && create_report_plots) {
      recovery_replications = read_summary_csv(destination, "parameter_replications")
      if (nrow(recovery_replications)) {
        report_index = save_report_plot(
          report_index, summary_dir, section, "noise_parameter_recovery",
          plot_noise_parameter_recovery(recovery_replications), width = 9, height = 5
        )
      }
    }
  }

  # Retire mixed-study displays now that the extension baselines are fitted.
  report_index$extensions$tables = setdiff(
    report_index$extensions$tables,
    c("report_gamma_with_reproduction", "report_noise_with_reproduction",
      "report_gamma_paired_with_reproduction", "report_noise_paired_with_reproduction")
  )
  report_index$extensions$figures = setdiff(
    report_index$extensions$figures,
    c("report_gamma_reference_mse", "report_gamma_reference_gain",
      "report_noise_reference_mse", "report_noise_reference_gain",
      "report_gamma_reference_recovery")
  )

  # Make the added tables and figures visible to the R Markdown report
  saveRDS(report_index, index_path)
  message("Comprehensive tables and main-report figures updated from cached summaries")
}

message("Knit scripts/summarize_results.Rmd to create the PDF")
