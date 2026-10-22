#!/usr/bin/env Rscript
# Shared simulation engine ---------------------------------------------------
# Use run_reproduction.R or run_extension.R for study-specific designs
# This file defines the shared functions; launch a study through either entry script

.runner_directory <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  path <- if (length(paths)) {
    tail(paths, 1)[[1]]
  } else if (length(arg)) {
    sub("^--file=", "", arg[1])
  } else {
    stop("Cannot locate this runner.")
  }
  dirname(normalizePath(path, mustWork = TRUE))
})

# 1. Load the four modules into one private environment -----------------------
# Functions share this environment, so changing settings does not require
# sourcing helpers into the user's global workspace
load_simulation <- function(directory = .runner_directory) {
  env <- new.env(parent = globalenv())
  for (name in c("config.R", "generate_data.R", "fit_models.R", "evaluate.R")) {
    sys.source(file.path(directory, name), envir = env)
  }
  env
}

# 2. Save complete objects without overwriting an existing result --------------
# Write to a temporary file beside the destination, then rename it. This avoids
# exposing a half-written RDS file if serialization fails
save_new_rds <- function(object, path) {
  if (file.exists(path)) {
    stop("Refusing to replace: ", path)
  }
  temporary <- tempfile(".saving-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not save: ", path)
  }
  invisible(path)
}

# PSgeo does not use J or areal measurement error. Reuse only when the full
# location mechanism, point variance, and underlying field are identical
# 3. Identify reusable PSgeo fits and count actual fitting attempts -----------
point_data_key <- function(plan) {
  columns <- c(
    "rho",
    "replication",
    "I",
    "sampling_type",
    "gamma",
    "gamma_L",
    "gamma_R",
    "split_x",
    "delta",
    "rho_U",
    "sigma2_U",
    "tau2_s"
  )
  do.call(paste, c(plan[columns], sep = "|"))
}

planned_fit_count <- function(config, plan) {
  # Melding/PSmelding use every dataset; PSgeo uses each distinct point dataset
  melding_fits <- nrow(plan) * sum(config$models != "PSgeo")
  psgeo_fits <- if ("PSgeo" %in% config$models) {
    length(unique(point_data_key(plan)))
  } else {
    0L
  }
  melding_fits + psgeo_fits
}

# 4. Run the study: field -> dataset -> model -> scores ------------------------
# Each run gets its own directory. Preserve failures and successful fits so a
# numerical problem in one model does not erase work from other scenarios
run_simulation <- function(config, engine = load_simulation(),
  shared = new.env(parent = emptyenv())) {
  # scenario_plan() validates the entire configuration before runtime setup
  plan <- engine$scenario_plan(config)
  packages <- engine$configure_runtime(config)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  run_dir <- file.path(config$output_dir,
    paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  if (dir.exists(run_dir) || !dir.create(run_dir)) {
    stop("Cannot create a fresh run directory.")
  }
  for (directory in c("simulation_data", "predictions", "inla_fits", "diagnostics")) {
    dir.create(file.path(run_dir, directory))
  }
  # Record settings and provenance once per run
  sources <- file.path(
    .runner_directory,
    c("config.R", "generate_data.R", "fit_models.R", "evaluate.R", "run_simulation.R")
  )
  if (!is.null(config$entry_script)) {
    sources <- c(sources, config$entry_script)
  }
  run_config <- list(
    config = config,
    plan = plan,
    source_md5 = tools::md5sum(sources),
    planned_physical_fits = planned_fit_count(config, plan),
    session = utils::sessionInfo(),
    package_versions = vapply(packages,
      function(p) as.character(utils::packageVersion(p)), character(1)),
    started_at = Sys.time()
  )
  writeLines(c(
    capture.output(print(run_config$session)),
    "", "Simulation package versions:",
    paste(names(run_config$package_versions), run_config$package_versions)
  ), file.path(run_dir, "session_info.txt"))
  message("Saving this run to ", run_dir)

  # Accumulate one status row per requested model comparison. Metrics and
  # parameter scores remain separate tables because they have different units
  status <- metrics <- parameters <- list()
  collect_score <- function(score, scenario, model, source_case, reused) {
    labels <- data.frame(scenario, model = model, source_case = source_case, fit_reused = reused)
    metrics[[length(metrics) + 1L]] <<- cbind(labels[rep(1L, nrow(score$metrics)), ], score$metrics)
    parameters[[length(parameters) + 1L]] <<- cbind(
      labels[rep(1L, nrow(score$parameters)), ],
      score$parameters
    )
  }
  # Reuse the same mesh and response fields across experiments run together
  mesh_settings <- list(config$mesh, config$generating$x_limits, config$generating$y_limits)
  if (!identical(shared$mesh_settings, mesh_settings)) {
    shared$mesh <- engine$build_mesh(config)
    shared$mesh_settings <- mesh_settings
  }
  mesh <- shared$mesh
  # Record resolved priors and the actual mesh even if every fitting attempt fails
  prior_defaults <- INLA::inla.models()
  fixed_defaults <- INLA::inla.set.control.fixed.default()
  run_config$fitting_setup <- list(
    mesh = mesh,
    mesh_settings = config$mesh,
    spde = INLA::inla.spde2.matern(mesh, alpha = 2, constr = TRUE),
    spde_settings = list(alpha = 2, constr = TRUE, prior.range = NULL, prior.sigma = NULL),
    priors = list(
      # b0 and b0pp are explicit regressors, so the ordinary fixed-effect prior applies
      response_intercept = list(distribution = "normal", mean = fixed_defaults$mean,
        precision = fixed_defaults$prec),
      sampling_intercept = list(distribution = "normal", mean = fixed_defaults$mean,
        precision = fixed_defaults$prec),
      gaussian_precision = prior_defaults$likelihood$gaussian$hyper$theta1[
        c("name", "prior", "param", "fixed")
      ],
      gaussian_precision_note = "Separate point and areal precisions; PSgeo has points only",
      preferential_coefficient = list(
        prior = prior_defaults$latent$copy$hyper$theta$prior,
        param = prior_defaults$latent$copy$hyper$theta$param,
        fixed = FALSE,
        models = c("PSgeo", "PSmelding")
      )
    ),
    numerical_controls = list(control.vb = list(strategy = "variance"))
  )
  save_new_rds(run_config, file.path(run_dir, "config.rds"))
  field_settings <- list(config$base_seed, config$generating[c(
    "x_limits", "y_limits", "grid_n", "mu", "sigma2", "nu"
  )])
  if (!identical(shared$field_settings, field_settings)) {
    shared$fields <- list()
    shared$field_settings <- field_settings
  }
  # [EXTENSION] Share the field across gamma/noise/I/J comparisons at a given range
  # generate_dataset() uses seeds that keep points independent of J and areas
  # independent of I. PSgeo reuse includes failed attempts: no automatic retry
  fields <- unique(plan[c("rho", "replication")])
  for (field_row in seq_len(nrow(fields))) {
    rho <- fields$rho[field_row]
    replication <- fields$replication[field_row]
    field_config <- config
    field_config$generating$rho <- rho
    field_key <- paste(rho, replication, sep = "|")
    if (is.null(shared$fields[[field_key]])) {
      shared$fields[[field_key]] <- engine$generate_surface(field_config, replication)
    }
    field <- shared$fields[[field_key]]
    field$independent_fields <- engine$generate_sampling_fields(field_config, replication)
    psgeo_cache <- list()
    for (row in which(plan$replication == replication & plan$rho == rho)) {
      scenario <- plan[row, , drop = FALSE]
      dataset <- engine$generate_dataset(config, scenario, field)
      data_file <- file.path("simulation_data", paste0(scenario$case_id, ".rds"))
      save_new_rds(dataset, file.path(run_dir, data_file))
      for (model in config$models) {
        fit_id <- paste(scenario$case_id, model, sep = "_")
        folder <- file.path(run_dir, "diagnostics", fit_id)
        dir.create(folder)
        fit_file <- file.path("inla_fits", paste0(fit_id, ".rds"))
        raw_file <- file.path("inla_fits", paste0(fit_id, "_raw.rds"))
        prediction_file <- file.path("predictions", paste0(fit_id, ".rds"))
        prediction_scores_file <- file.path("predictions", paste0(fit_id, "_scores.rds"))
        point_key <- point_data_key(scenario)
        cached <- if (model == "PSgeo") {
          psgeo_cache[[point_key]]
        } else {
          NULL
        }
        # Reused PSgeo comparisons point to the first attempt, including a
        # failed attempt. This avoids changing the experiment through retries
        if (!is.null(cached)) {
          save_new_rds(
            list(
              source_case = cached$case_id,
              source_directory = file.path("..", paste0(cached$case_id, "_PSgeo")),
              reason = "Identical point data, prediction targets, and prior; PSgeo does not use areas."
            ),
            file.path(folder, "reuse.rds")
          )
          entry <- cached$status
          entry[names(scenario)] <- scenario
          entry$fit_reused <- TRUE
          entry$simulation_file <- data_file
          if (!is.null(cached$score)) {
            collect_score(cached$score, scenario, model, cached$case_id, TRUE)
            save_new_rds(cached$score, file.path(folder, "scores.rds"))
          }
          status[[length(status) + 1L]] <- entry
          next
        }
        stage <- "fitting"
        score <- bundle <- failure <- NULL
        attempt_warnings <- character()
        attempt_started <- proc.time()
        message(scenario$case_id, ": ", model)
        fit_seed <- engine$sampling_seed(config, scenario, "fit", model)
        score_seed <- engine$sampling_seed(config, scenario, "scores", model)
        save_new_rds(
          list(model = model, scenario = scenario,
            seeds = list(data = dataset$seeds, fitting = fit_seed, scoring = score_seed),
            started_at = Sys.time()),
          file.path(folder, "attempt.rds")
        )
        # Preserve the full INLA return before processing predictions or scores
        # Remove the raw checkpoint only after the complete fit bundle is saved
        entry <- withCallingHandlers(tryCatch(
          {
            bundle <- engine$fit_model(
              model,
              dataset,
              mesh,
              config,
              on_inla_return = function(fit, context) {
                save_new_rds(list(fit = fit, context = context), file.path(run_dir, raw_file))
              }
            )
            save_new_rds(bundle, file.path(run_dir, fit_file))
            unlink(file.path(run_dir, raw_file)) # Complete fit is now durable
            stage <- "predictions"
            save_new_rds(
              list(scenario = scenario, model = model, predictions = bundle$predictions),
              file.path(run_dir, prediction_file)
            )
            stage <- "scoring"
            score <- engine$evaluate_fit(bundle, dataset$truth, config, score_seed)
            save_new_rds(score, file.path(folder, "scores.rds"))
            save_new_rds(
              list(scenario = scenario, model = model, pointwise = score$pointwise),
              file.path(run_dir, prediction_scores_file)
            )
            collect_score(score, scenario, model, scenario$case_id, FALSE)
            data.frame(
              scenario,
              model = model,
              source_case = scenario$case_id,
              fit_reused = FALSE,
              fit_saved = TRUE,
              scored = TRUE,
              mode_status = score$diagnostics$mode_status,
              needs_review = score$diagnostics$needs_review,
              message = ""
            )
          },
          error = function(error) {
            failure <<- list(
              stage = stage,
              message = conditionMessage(error),
              classes = class(error),
              warnings = unique(c(attempt_warnings, error$warnings)),
              runtime = (proc.time() - attempt_started)[c("user.self", "sys.self", "elapsed")],
              fitting_diagnostics = error$diagnostics,
              occurred_at = Sys.time()
            )
            save_new_rds(failure, file.path(folder, "error.rds"))
            data.frame(
              scenario,
              model = model,
              source_case = scenario$case_id,
              fit_reused = FALSE,
              fit_saved = file.exists(file.path(run_dir, fit_file)) ||
                file.exists(file.path(run_dir, raw_file)),
              scored = FALSE,
              mode_status = NA_real_,
              needs_review = TRUE,
              message = conditionMessage(error)
            )
          }
        ), warning = function(warning) {
          attempt_warnings <<- c(attempt_warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        })
        elapsed <- unname((proc.time() - attempt_started)["elapsed"])
        # Keep diagnostics for successful fits and failed fitting/scoring attempts
        save_new_rds(
          list(
            model = model,
            scenario = scenario,
            fit_saved = entry$fit_saved,
            scored = entry$scored,
            needs_review = entry$needs_review,
            warnings = unique(c(attempt_warnings, bundle$warnings, failure$warnings)),
            elapsed_seconds = elapsed,
            fitting = if (!is.null(bundle)) bundle$diagnostics else failure$fitting_diagnostics,
            fitting_runtime = bundle$inla_runtime,
            hessian_adjustment_reported = if (!is.null(bundle)) {
              any(grepl("ensure spd for Hessian", bundle$fit$logfile))
            } else {
              NA
            },
            scoring = score$diagnostics,
            error = failure
          ),
          file.path(folder, "diagnostics.rds")
        )
        entry$elapsed_seconds <- elapsed
        entry$warnings <- paste(unique(c(attempt_warnings, bundle$warnings, failure$warnings)),
          collapse = " | ")
        entry$hessian_adjustment_reported <- if (!is.null(bundle)) {
          any(grepl("ensure spd for Hessian", bundle$fit$logfile))
        } else {
          NA
        }
        entry$failed_stage <- if (is.null(failure)) "" else failure$stage
        if (!is.null(bundle$diagnostics$mode_status)) {
          entry$mode_status <- bundle$diagnostics$mode_status
        }
        # Paths are relative to this run directory; reused fits retain source paths
        entry$simulation_file <- data_file
        entry$prediction_file <- if (file.exists(file.path(run_dir, prediction_file))) prediction_file else NA_character_
        entry$prediction_scores_file <- if (file.exists(file.path(run_dir, prediction_scores_file))) {
          prediction_scores_file
        } else {
          NA_character_
        }
        entry$fit_file <- if (file.exists(file.path(run_dir, fit_file))) {
          fit_file
        } else if (file.exists(file.path(run_dir, raw_file))) {
          raw_file
        } else {
          NA_character_
        }
        entry$diagnostics_directory <- file.path("diagnostics", fit_id)
        status[[length(status) + 1L]] <- entry
        if (model == "PSgeo") {
          psgeo_cache[[point_key]] <- list(
            case_id = scenario$case_id,
            status = entry,
            score = if (isTRUE(entry$scored)) {
              score
            } else {
              NULL
            }
          )
        }
      }
      utils::write.csv(do.call(rbind, status), file.path(run_dir, "fit_diagnostics.csv"), row.names = FALSE)
    }
  }
  # 4E. Export comparison tables. summarize_metrics() reports valid-score
  # counts as well as averages, so missing results remain visible
  status <- do.call(rbind, status)
  utils::write.csv(status, file.path(run_dir, "fit_diagnostics.csv"), row.names = FALSE)
  if (length(metrics)) {
    metrics <- do.call(rbind, metrics)
    utils::write.csv(metrics, file.path(run_dir, "predictive_scores.csv"), row.names = FALSE)
    utils::write.csv(
      engine$summarize_metrics(metrics, config$n_replicates),
      file.path(run_dir, "predictive_score_summary.csv"),
      row.names = FALSE
    )
  }
  if (length(parameters)) {
    utils::write.csv(do.call(rbind, parameters),
      file.path(run_dir, "parameter_estimates.csv"), row.names = FALSE)
  }
  message(
    sum(status$scored),
    "/",
    nrow(status),
    " comparison rows scored; ",
    sum(!status$fit_reused),
    " physical fitting attempts. See fit_diagnostics.csv."
  )
  invisible(list(directory = run_dir, status = status))
}

# 5. Shared command-line handling --------------------------------------------
# The default is a preview. Only --run starts fitting; --test reduces
# numerical resolution for a workflow check and saves under a separate folder
study_main <- function(args, config_factory, entry_script,
  sections = c("main", "appendix"), default_section = "main",
  test_all_conditions = FALSE) {
  allowed <- args %in% c("--run", "--test", "--help") |
    args %in% paste0("--section=", c(sections, "all")) |
    grepl("^--replicates=[1-9][0-9]*$", args)
  section_args <- grep("^--section=", args, value = TRUE)
  replication_args <- grep("^--replicates=", args, value = TRUE)
  usage <- paste(
    basename(entry_script),
    "[--run]",
    paste0("[--section=", paste(c(sections, "all"), collapse = "|"), "]"),
    "[--replicates=N] [--test]"
  )
  if ("--help" %in% args) {
    cat(usage, "\n")
    return(invisible(NULL))
  }
  if (any(!allowed) || length(section_args) > 1L || length(replication_args) > 1L) {
    stop(usage)
  }
  section <- if (length(section_args)) {
    sub("^--section=", "", section_args)
  } else {
    default_section
  }
  sections <- if (section == "all") {
    sections
  } else {
    section
  }
  engine <- load_simulation()
  results <- list()
  shared <- new.env(parent = emptyenv())
  for (selected in sections) {
    config <- config_factory(selected, dirname(dirname(entry_script)))
    config$entry_script <- entry_script
    if (length(replication_args)) {
      config$n_replicates <- as.numeric(sub("^--replicates=", "", replication_args))
    }
    if ("--test" %in% args) {
      config$n_replicates <- 1L
      config$generating$rho <- config$generating$rho[1L]
      config$generating$I <- 30L
      config$generating$J <- 4L
      config$generating$grid_n <- 10L
      if (!test_all_conditions) {
        config$generating$gamma <- if (1 %in% config$generating$gamma) {
          1
        } else {
          config$generating$gamma[1L]
        }
        config$generating$noise <- config$generating$noise[1L, , drop = FALSE]
        if (!is.null(config$generating$sampling)) {
          config$generating$sampling <- config$generating$sampling[1L, , drop = FALSE]
        }
      }
      config$mesh$max_edge <- c(.25, .5)
      config$scoring$joint_draws <- 50L
      config$output_dir <- file.path(config$output_dir, "test")
    }
    plan <- engine$scenario_plan(config)
    cat("\n", config$study, " / ", selected, "\n", sep = "")
    print(
      unique(plan[c(
        "rho",
        "I",
        "J",
        "sampling_id",
        "gamma",
        "gamma_L",
        "gamma_R",
        "delta",
        "rho_U",
        "sigma2_U",
        "tau2_s",
        "tau2_B",
        "log_intensity_phi_cor"
      )]),
      row.names = FALSE
    )
    cat(
      "Replications per scenario:",
      config$n_replicates,
      "| Physical fits:",
      planned_fit_count(config, plan),
      "\n"
    )
    cat("Output:", config$output_dir, "\n")
    if ("--run" %in% args) {
      results[[selected]] <- run_simulation(config, engine, shared)
    }
    else {
      results[[selected]] <- config
    }
  }
  if (!"--run" %in% args) {
    cat("\nPreview only. Add --run to start fitting.\n")
  }
  else if (any(vapply(results, function(x) any(!x$status$scored), logical(1)))) {
    stop("Some jobs failed; inspect fit_diagnostics.csv and the preserved errors in the printed output directories.")
  }
  invisible(results)
}
