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
  shared = new.env(parent = emptyenv()), run_dir = NULL) {
  # scenario_plan() validates the entire configuration before runtime setup
  plan <- selected_scenario_plan(config, engine)
  packages <- engine$configure_runtime(config)
  versions <- runner_package_versions(packages)
  source_hashes <- runner_source_hashes()
  identity <- run_identity(config, source_hashes, versions)
  if (is.null(run_dir)) {
    run_dir <- file.path(config$output_dir,
      paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
    if (dir.exists(run_dir)) stop("Run directory already exists: ", run_dir)
  }
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  release_lock <- acquire_run_lock(run_dir)
  on.exit(release_lock(), add = TRUE)
  for (directory in c("simulation_data", "predictions", "inla_fits", "diagnostics")) {
    dir.create(file.path(run_dir, directory), showWarnings = FALSE)
  }
  config_file <- file.path(run_dir, "config.rds")
  resuming <- file.exists(config_file)
  if (resuming) {
    run_config <- readRDS(config_file)
    if (!identical(run_config$identity, identity) || !identical(run_config$plan, plan)) {
      stop("Resume settings, source code or runtime differ from this batch; use a new run ID")
    }
    if (file.exists(file.path(run_dir, "batch_manifest.rds"))) {
      status <- verify_completed_run(run_dir)
      message("Verified completed batch: ", run_dir)
      return(invisible(list(directory = run_dir, status = status)))
    }
  } else {
    run_config <- list(
      config = config, plan = plan, identity = identity,
      source_sha256 = source_hashes,
      planned_physical_fits = planned_fit_count(config, plan),
      session = utils::sessionInfo(), package_versions = versions,
      started_at = Sys.time()
    )
    writeLines(c(capture.output(print(run_config$session)), "",
      "Simulation package versions:", paste(names(versions), versions)),
      file.path(run_dir, "session_info.txt"))
  }
  message(if (resuming) "Resuming " else "Saving this run to ", run_dir)

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
  if (resuming) {
    shared$mesh <- run_config$fitting_setup$mesh
    shared$mesh_settings <- mesh_settings
  } else if (!identical(shared$mesh_settings, mesh_settings)) {
    shared$mesh <- engine$build_mesh(config)
    shared$mesh_settings <- mesh_settings
  }
  mesh <- shared$mesh
  # Record resolved priors and the actual mesh even if every fitting attempt fails
  if (!resuming) {
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
        gaussian_precision = c(
          list(name = prior_defaults$likelihood$gaussian$hyper$theta1$name),
          config$fitting$noise_prior,
          list(fixed = FALSE)
        ),
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
    save_new_rds(run_config, config_file)
  }
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
      data_file <- file.path("simulation_data", paste0(scenario$case_id, ".rds"))
      data_manifest <- file.path(run_dir, "simulation_data",
        paste0(scenario$case_id, "_manifest.rds"))
      if (file.exists(data_manifest)) {
        verify_inventory(run_dir, readRDS(data_manifest))
      }
      dataset <- if (file.exists(file.path(run_dir, data_file))) {
        readRDS(file.path(run_dir, data_file))
      } else {
        value <- engine$generate_dataset(config, scenario, field)
        save_new_rds(value, file.path(run_dir, data_file))
        value
      }
      if (!file.exists(data_manifest)) {
        save_new_rds(inventory_files(run_dir, data_file), data_manifest)
      }
      for (model in config$models) {
        fit_id <- paste(scenario$case_id, model, sep = "_")
        folder <- file.path(run_dir, "diagnostics", fit_id)
        dir.create(folder, showWarnings = FALSE)
        fit_file <- file.path("inla_fits", paste0(fit_id, ".rds"))
        raw_file <- file.path("inla_fits", paste0(fit_id, "_raw.rds"))
        prediction_file <- file.path("predictions", paste0(fit_id, ".rds"))
        prediction_scores_file <- file.path("predictions", paste0(fit_id, "_scores.rds"))
        point_key <- point_data_key(scenario)
        completed <- read_job_checkpoint(run_dir, folder)
        if (!is.null(completed)) {
          entry <- completed$status
          status[[length(status) + 1L]] <- entry
          if (!is.null(completed$score)) {
            collect_score(completed$score, scenario, model, entry$source_case, entry$fit_reused)
          }
          if (model == "PSgeo" && !isTRUE(entry$fit_reused)) {
            psgeo_cache[[point_key]] <- list(case_id = scenario$case_id,
              status = entry, score = completed$score)
          }
          next
        }
        cached <- if (model == "PSgeo") {
          psgeo_cache[[point_key]]
        } else {
          NULL
        }
        # Reused PSgeo comparisons point to the first attempt, including a
        # failed attempt. This avoids changing the experiment through retries
        if (!is.null(cached)) {
          if (!file.exists(file.path(folder, "reuse.rds"))) save_new_rds(
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
            if (!file.exists(file.path(folder, "scores.rds"))) {
              save_new_rds(cached$score, file.path(folder, "scores.rds"))
            }
          }
          finish_job_checkpoint(entry, run_dir, folder)
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
        previous_attempt <- file.exists(file.path(folder, "attempt.rds"))
        if (!previous_attempt) save_new_rds(
          list(model = model, scenario = scenario,
            seeds = list(data = dataset$seeds, fitting = fit_seed, scoring = score_seed),
            started_at = Sys.time()),
          file.path(folder, "attempt.rds")
        )
        # Preserve the full INLA return before processing predictions or scores
        # Remove the raw checkpoint only after the complete fit bundle is saved
        entry <- withCallingHandlers(tryCatch(
          {
            if (file.exists(file.path(folder, "error.rds"))) {
              failure <- readRDS(file.path(folder, "error.rds"))
              stop(failure$message)
            }
            if (file.exists(file.path(run_dir, fit_file))) {
              bundle <- readRDS(file.path(run_dir, fit_file))
            } else {
              if (previous_attempt) {
                stage <- "interrupted"
                stop("Interrupted fitting attempt retained for review; no complete fit is available for scoring and no automatic refit was attempted")
              }
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
            }
            stage <- "predictions"
            if (!file.exists(file.path(run_dir, prediction_file))) save_new_rds(
              list(scenario = scenario, model = model, predictions = bundle$predictions),
              file.path(run_dir, prediction_file)
            )
            stage <- "scoring"
            score <- if (file.exists(file.path(folder, "scores.rds"))) {
              readRDS(file.path(folder, "scores.rds"))
            } else {
              value <- engine$evaluate_fit(bundle, dataset$truth, config, score_seed)
              save_new_rds(value, file.path(folder, "scores.rds"))
              value
            }
            if (!file.exists(file.path(run_dir, prediction_scores_file))) save_new_rds(
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
            if (file.exists(file.path(folder, "error.rds"))) {
              failure <<- readRDS(file.path(folder, "error.rds"))
            } else {
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
            }
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
        if (!file.exists(file.path(folder, "diagnostics.rds"))) save_new_rds(
          list(
            model = model,
            scenario = scenario,
            fit_saved = entry$fit_saved,
            scored = entry$scored,
            needs_review = entry$needs_review,
            warnings = unique(c(attempt_warnings, bundle$warnings, failure$warnings)),
            elapsed_seconds = elapsed,
            elapsed_scope = if (previous_attempt) {
              "Resume session only; original INLA runtime is retained separately when available"
            } else "Complete attempt in this session",
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
        # Recover original timing and warnings if interruption followed diagnostics
        saved_diagnostics <- readRDS(file.path(folder, "diagnostics.rds"))
        entry$elapsed_seconds <- saved_diagnostics$elapsed_seconds
        entry$elapsed_scope <- saved_diagnostics$elapsed_scope
        entry$warnings <- paste(saved_diagnostics$warnings, collapse = " | ")
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
        finish_job_checkpoint(entry, run_dir, folder)
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
      atomic_write_csv(do.call(rbind, status), file.path(run_dir, "fit_diagnostics.csv"))
    }
  }
  # 4E. Export comparison tables. summarize_metrics() reports valid-score
  # counts as well as averages, so missing results remain visible
  status <- do.call(rbind, status)
  atomic_write_csv(status, file.path(run_dir, "fit_diagnostics.csv"))
  if (length(metrics)) {
    metrics <- do.call(rbind, metrics)
    atomic_write_csv(metrics, file.path(run_dir, "predictive_scores.csv"))
    atomic_write_csv(
      engine$summarize_metrics(metrics, length(unique(plan$replication))),
      file.path(run_dir, "predictive_score_summary.csv")
    )
  }
  if (length(parameters)) {
    atomic_write_csv(do.call(rbind, parameters),
      file.path(run_dir, "parameter_estimates.csv"))
  }
  finalize_run_manifest(run_dir, config, status)
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


# 5. Batch checkpoints and verification --------------------------------------
# Write derived tables atomically; original fitting artifacts use save_new_rds
atomic_write_csv <- function(object, path) {
  temporary <- tempfile(".saving-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop("Could not save: ", path)
  invisible(path)
}

atomic_save_rds <- function(object, path) {
  temporary <- tempfile(".saving-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop("Could not save: ", path)
  invisible(path)
}

# Assign case IDs from the full study before selecting a laptop or batch
selected_scenario_plan <- function(config, engine) {
  plan <- engine$scenario_plan(config)
  ids <- config$replicate_ids
  if (is.null(ids)) return(plan)
  if (!is.numeric(ids) || !length(ids) || anyNA(ids) ||
      any(!is.finite(ids)) || any(ids != floor(ids)) || anyDuplicated(ids) ||
      any(ids < 1L | ids > config$n_replicates)) {
    stop("replicate_ids must be distinct integers between 1 and n_replicates.")
  }
  plan[plan$replication %in% ids, , drop = FALSE]
}

runner_source_hashes <- function() {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package digest is required for SHA-256 verification.")
  }
  names <- c("config.R", "generate_data.R", "fit_models.R", "evaluate.R",
    "run_simulation.R", "run_reproduction.R", "run_extension.R")
  paths <- file.path(.runner_directory, names)
  if (!all(file.exists(paths))) stop("Cannot find all seven simulation scripts.")
  stats::setNames(vapply(paths, function(path) {
    digest::digest(file = path, algo = "sha256")
  }, character(1)), names)
}

# Include installed dependencies without depending on which namespaces are loaded
runner_package_versions <- function(packages) {
  installed <- utils::installed.packages()
  installed <- installed[!duplicated(installed[, "Package"]), , drop = FALSE]
  packages <- unique(c(packages, "digest"))
  if (!all(packages %in% installed[, "Package"])) {
    stop("A required package is missing from the configured R libraries.")
  }
  dependencies <- tools::package_dependencies(packages, db = installed,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)
  packages <- sort(unique(c(packages, unlist(dependencies, use.names = FALSE))))
  packages <- setdiff(packages, "R")
  if (!all(packages %in% installed[, "Package"])) {
    stop("An installed simulation package has a missing dependency.")
  }
  stats::setNames(as.character(installed[match(packages, installed[, "Package"]), "Version"]), packages)
}

# Local paths may differ across laptops; scientific and numerical settings may not
run_identity <- function(config, source_hashes, versions) {
  config[c("project_dir", "output_dir", "package_library", "entry_script")] <- NULL
  list(config = config, source_sha256 = source_hashes,
    package_versions = versions,
    R = list(version = as.character(getRversion()), platform = R.version$platform),
    RNGkind = RNGkind())
}

# One writer owns a batch; an abandoned lock is recoverable only on its own host
acquire_run_lock <- function(run_dir) {
  lock <- file.path(run_dir, ".lock")
  owner_file <- file.path(lock, "owner.rds")
  hostname <- unname(Sys.info()["nodename"])
  read_owner <- function() {
    tryCatch(readRDS(owner_file), error = function(error) NULL)
  }
  dead_local_owner <- function(owner) {
    if (is.null(owner) || !identical(owner$hostname, hostname) ||
        length(owner$pid) != 1L || is.na(owner$pid) || owner$pid < 1L) {
      return(FALSE)
    }
    alive <- tryCatch(tools::pskill(owner$pid, signal = 0L),
      error = function(error) NA)
    identical(alive, FALSE)
  }
  if (!dir.create(lock, showWarnings = FALSE)) {
    if (!dead_local_owner(read_owner())) {
      stop("Batch is locked by a live, remote, or unverified owner: ", lock)
    }
    # Claim stale-lock recovery and recheck the owner before moving the lock
    # A second process cannot reclaim a newly acquired lock using old owner data
    reclaim <- file.path(lock, ".reclaim")
    if (!dir.create(reclaim, showWarnings = FALSE)) {
      stop("Another process is recovering this lock: ", lock)
    }
    if (!dead_local_owner(read_owner())) {
      unlink(reclaim, recursive = TRUE)
      stop("The lock owner changed during recovery: ", lock)
    }
    history <- file.path(run_dir, "lock_history")
    dir.create(history, showWarnings = FALSE)
    archived_lock <- tempfile("stale-", tmpdir = history)
    if (!file.rename(lock, archived_lock)) {
      unlink(reclaim, recursive = TRUE)
      stop("Could not preserve the abandoned lock: ", lock)
    }
    if (!dir.create(lock, showWarnings = FALSE)) {
      stop("Another process acquired the batch lock: ", lock)
    }
  }
  token <- basename(tempfile("owner-"))
  owner <- list(hostname = hostname, pid = Sys.getpid(), token = token,
    acquired_at = Sys.time())
  tryCatch(save_new_rds(owner, owner_file), error = function(error) {
    unlink(lock, recursive = TRUE)
    stop(error)
  })
  function() {
    current <- read_owner()
    if (!is.null(current) && identical(current$token, token)) {
      unlink(lock, recursive = TRUE)
    }
    invisible(NULL)
  }
}

inventory_files <- function(run_dir, files) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package digest is required for SHA-256 verification.")
  }
  files <- sort(unique(as.character(files)))
  if (anyNA(files) || any(!nzchar(files)) ||
      any(grepl("^/|(^|/)\\.\\.(/|$)", files))) {
    stop("Inventoried files must have safe paths relative to the batch.")
  }
  root <- normalizePath(run_dir, mustWork = TRUE)
  paths <- file.path(root, files)
  info <- file.info(paths)
  if (anyNA(info$size) || any(info$isdir)) {
    stop("A required output file is missing: ", paste(files[is.na(info$size) | info$isdir], collapse = ", "))
  }
  resolved <- vapply(paths, normalizePath, character(1), mustWork = TRUE)
  if (any(!startsWith(resolved, paste0(root, "/")))) {
    stop("An output file points outside its batch directory.")
  }
  data.frame(path = files, bytes = as.numeric(info$size),
    sha256 = vapply(paths, function(path) {
      digest::digest(file = path, algo = "sha256")
    }, character(1)), row.names = NULL, stringsAsFactors = FALSE)
}

verify_inventory <- function(run_dir, inventory) {
  if (!is.data.frame(inventory) ||
      !all(c("path", "bytes", "sha256") %in% names(inventory)) ||
      anyDuplicated(inventory$path)) {
    stop("Invalid output inventory.")
  }
  current <- inventory_files(run_dir, inventory$path)
  expected <- inventory[match(current$path, inventory$path), , drop = FALSE]
  changed <- is.na(expected$bytes) | is.na(expected$sha256) |
    current$bytes != expected$bytes | current$sha256 != expected$sha256
  if (any(changed)) {
    stop("Output verification failed: ", paste(current$path[changed], collapse = ", "))
  }
  invisible(TRUE)
}

# The terminal status is written last by the runner, after the fit attempt ends
# A recorded failure is terminal too; resuming never selectively retries failures
finish_job_checkpoint <- function(entry, run_dir, folder) {
  status_file <- file.path(folder, "status.rds")
  if (!file.exists(status_file)) {
    save_new_rds(entry, status_file)
  } else if (!identical(readRDS(status_file), entry)) {
    stop("Refusing to change a completed fit status: ", folder)
  }
  columns <- c("simulation_file", "fit_file", "prediction_file", "prediction_scores_file")
  if (isTRUE(entry$scored) && any(vapply(columns, function(name) {
    value <- entry[[name]]
    length(value) != 1L || is.na(value) || !nzchar(value)
  }, logical(1)))) {
    stop("A scored fit is missing a required artifact reference: ", folder)
  }
  files <- unlist(entry[intersect(columns, names(entry))], use.names = FALSE)
  files <- files[!is.na(files) & nzchar(files)]
  root <- normalizePath(run_dir, mustWork = TRUE)
  normalized_folder <- normalizePath(folder, mustWork = TRUE)
  if (!startsWith(normalized_folder, paste0(root, "/"))) {
    stop("The job checkpoint folder must be inside its batch directory.")
  }
  relative_folder <- substring(normalized_folder, nchar(root) + 2L)
  own_files <- list.files(folder, recursive = TRUE, all.files = FALSE)
  own_files <- own_files[own_files != "checkpoint_manifest.rds"]
  files <- c(files, file.path(relative_folder, own_files))
  if (isTRUE(entry$scored)) {
    diagnostics_folder <- entry$diagnostics_directory
    if (is.null(diagnostics_folder) || is.na(diagnostics_folder) || !nzchar(diagnostics_folder)) {
      diagnostics_folder <- relative_folder
    }
    files <- c(files, file.path(relative_folder, "scores.rds"),
      file.path(diagnostics_folder, "diagnostics.rds"))
  }
  inventory <- inventory_files(run_dir, files)
  manifest <- file.path(folder, "checkpoint_manifest.rds")
  if (!file.exists(manifest)) {
    save_new_rds(inventory, manifest)
  } else {
    verify_inventory(run_dir, readRDS(manifest))
  }
  invisible(inventory)
}

read_job_checkpoint <- function(run_dir, folder) {
  status_file <- file.path(folder, "status.rds")
  if (!file.exists(status_file)) return(NULL)
  entry <- readRDS(status_file)
  manifest <- file.path(folder, "checkpoint_manifest.rds")
  if (!file.exists(manifest)) {
    # Recover an interruption between saving the terminal status and its hashes
    finish_job_checkpoint(entry, run_dir, folder)
  }
  verify_inventory(run_dir, readRDS(manifest))
  list(status = entry, score = if (isTRUE(entry$scored)) {
    readRDS(file.path(folder, "scores.rds"))
  } else NULL)
}

finalize_run_manifest <- function(run_dir, config, status) {
  if (!is.data.frame(status) || !nrow(status) || anyNA(status$scored)) {
    stop("A completed batch must contain terminal comparison statuses.")
  }
  atomic_save_rds(status, file.path(run_dir, "status.rds"))
  files <- list.files(run_dir, recursive = TRUE, all.files = FALSE)
  files <- files[!files %in% c("file_manifest.csv", "batch_manifest.rds") &
    !startsWith(files, "lock_history/")]
  inventory <- inventory_files(run_dir, files)
  manifest_file <- file.path(run_dir, "file_manifest.csv")
  atomic_write_csv(inventory, manifest_file)
  verify_inventory(run_dir, inventory)
  manifest <- list(
    status = if (all(status$scored)) "complete" else "complete_with_failures",
    experiment = config$section,
    replicate_ids = sort(unique(status$replication)),
    comparisons = nrow(status),
    physical_attempts = sum(!status$fit_reused),
    scored = sum(status$scored),
    failures = sum(!status$scored),
    needs_review = sum(status$needs_review, na.rm = TRUE),
    verification = "Local files checked against SHA-256; Drive upload not verified",
    file_manifest_sha256 = digest::digest(file = manifest_file, algo = "sha256"),
    verified_at = Sys.time()
  )
  atomic_save_rds(manifest, file.path(run_dir, "batch_manifest.rds"))
  invisible(manifest)
}

verify_completed_run <- function(run_dir) {
  manifest <- readRDS(file.path(run_dir, "batch_manifest.rds"))
  if (!manifest$status %in% c("complete", "complete_with_failures")) {
    stop("This batch has no valid completion record.")
  }
  manifest_file <- file.path(run_dir, "file_manifest.csv")
  if (!identical(digest::digest(file = manifest_file, algo = "sha256"),
      manifest$file_manifest_sha256)) {
    stop("The batch file manifest has changed.")
  }
  inventory <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
  files <- list.files(run_dir, recursive = TRUE, all.files = FALSE)
  files <- files[!files %in% c("file_manifest.csv", "batch_manifest.rds") &
    !startsWith(files, "lock_history/")]
  if (!setequal(files, inventory$path)) {
    stop("The completed batch has missing or unexpected files.")
  }
  verify_inventory(run_dir, inventory)
  status <- readRDS(file.path(run_dir, "status.rds"))
  metadata <- readRDS(file.path(run_dir, "config.rds"))
  plan <- metadata$plan
  models <- metadata$config$models
  if (!is.data.frame(status) || !is.data.frame(plan) || !nrow(plan) || !length(models) ||
      !all(c("experiment", "case_id", "model", "replication", "scored", "fit_reused", "needs_review") %in% names(status)) ||
      anyNA(status$scored) || anyNA(status$fit_reused)) {
    stop("The batch is missing its plan or terminal job records.")
  }
  expected <- expand.grid(row = seq_len(nrow(plan)), model = models, stringsAsFactors = FALSE)
  expected_keys <- paste(plan$experiment[expected$row], plan$case_id[expected$row], expected$model, sep = "|")
  actual_keys <- paste(status$experiment, status$case_id, status$model, sep = "|")
  if (anyDuplicated(actual_keys) || !setequal(expected_keys, actual_keys)) {
    stop("Completed batch job records do not match the saved plan.")
  }
  actual_counts <- list(comparisons = nrow(status),
    physical_attempts = sum(!status$fit_reused), scored = sum(status$scored),
    failures = sum(!status$scored), needs_review = sum(status$needs_review, na.rm = TRUE))
  expected_status <- if (all(status$scored)) "complete" else "complete_with_failures"
  if (!identical(manifest[names(actual_counts)], actual_counts) ||
      !identical(manifest$status, expected_status) ||
      !identical(manifest$experiment, metadata$config$section) ||
      !identical(manifest$replicate_ids, sort(unique(status$replication)))) {
    stop("Batch manifest counts or labels disagree with the saved job records.")
  }
  status
}

# 6. Shared command-line handling --------------------------------------------
# Named runs resume automatically after checking their settings and artifacts
# --test uses a smaller design and saves under a separate test folder
study_main <- function(args, config_factory, entry_script,
  sections = c("main", "appendix"), default_section = "main",
  test_all_conditions = FALSE) {
  flags <- c("--run", "--test", "--help", "--verify-only")
  options <- c("section", "replicates", "replicate-start", "replicate-end",
    "batch-size", "run-id", "laptop", "output-root")
  keys <- sub("=.*$", "", sub("^--", "", args))
  usage <- paste(basename(entry_script), "[--run | --verify-only]",
    paste0("[--section=", paste(c(sections, "all"), collapse = "|"), "]"),
    "[--replicates=N] [--replicate-start=N --replicate-end=N]",
    "[--run-id=ID --laptop=laptop_1|laptop_2 --batch-size=5]",
    "[--output-root=PATH] [--test]")
  if ("--help" %in% args) {
    cat(usage, "\nNamed runs resume automatically; --replicates sets the full study size\n")
    return(invisible(NULL))
  }
  allowed <- args %in% flags | (grepl("^--[^=]+=.+$", args) & keys %in% options)
  if (any(!allowed) || anyDuplicated(keys) ||
      all(c("--run", "--verify-only") %in% args)) stop(usage)
  option <- function(name, default = NULL) {
    value <- args[keys == name]
    if (length(value)) sub("^[^=]+=", "", value) else default
  }
  integer_option <- function(name, default = NULL) {
    value <- option(name)
    if (is.null(value)) return(default)
    if (!grepl("^[1-9][0-9]*$", value) ||
        !is.finite(as.numeric(value)) || as.numeric(value) > .Machine$integer.max) {
      stop("--", name, " must be a positive integer")
    }
    as.integer(value)
  }
  section <- option("section", default_section)
  if (!section %in% c(sections, "all")) stop(usage)
  sections <- if (section == "all") sections else section
  run_id <- option("run-id")
  laptop <- option("laptop")
  batch_size <- integer_option("batch-size", 5L)
  start <- integer_option("replicate-start")
  end <- integer_option("replicate-end")
  if (xor(is.null(start), is.null(end))) stop("Specify both replicate-start and replicate-end")
  if (!is.null(start) && start > end) stop("replicate-start exceeds replicate-end")
  if (!is.null(run_id) && !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id)) {
    stop("Use letters, digits, underscores or hyphens in run-id")
  }
  if (!is.null(run_id) && (is.null(laptop) || !laptop %in% c("laptop_1", "laptop_2"))) {
    stop("Named runs require --laptop=laptop_1 or --laptop=laptop_2")
  }
  if (is.null(run_id) && any(c("batch-size", "laptop", "output-root", "verify-only") %in% keys)) {
    stop("Batch, output-root and verification options require --run-id")
  }
  if (!is.null(run_id) && !"--test" %in% args && is.null(start)) {
    stop("Named laptop runs require an explicit replicate-start and replicate-end")
  }
  if ("--test" %in% args && (!is.null(start) || "replicates" %in% keys)) {
    stop("--test already selects one small replicate; omit replicate-count/range options")
  }
  engine <- load_simulation()
  results <- list()
  shared <- new.env(parent = emptyenv())
  for (selected in sections) {
    config <- config_factory(selected, dirname(dirname(entry_script)))
    config$entry_script <- entry_script
    config$n_replicates <- integer_option("replicates", config$n_replicates)
    if ("--test" %in% args) {
      config$n_replicates <- 1L
      config$generating$rho <- config$generating$rho[1L]
      config$generating$I <- 30L
      config$generating$J <- 4L
      config$generating$grid_n <- 10L
      if (!test_all_conditions) {
        config$generating$gamma <- if (1 %in% config$generating$gamma) 1 else config$generating$gamma[1L]
        config$generating$noise <- config$generating$noise[1L, , drop = FALSE]
        if (!is.null(config$generating$sampling)) {
          config$generating$sampling <- config$generating$sampling[1L, , drop = FALSE]
        }
      }
      config$mesh$max_edge <- c(.25, .5)
      config$scoring$joint_draws <- 50L
      config$output_dir <- file.path(config$output_dir, "test")
    }
    if (!is.null(end) && end > config$n_replicates) {
      stop("replicate_ids exceed the full study size")
    }
    ids <- if (is.null(start)) seq_len(config$n_replicates) else seq.int(start, end)
    config$replicate_ids <- ids
    plan <- selected_scenario_plan(config, engine)
    cat("\n", config$study, " / ", selected, "\n", sep = "")
    print(unique(plan[c("rho", "I", "J", "sampling_id", "gamma", "gamma_L", "gamma_R",
      "delta", "rho_U", "sigma2_U", "tau2_s", "tau2_B", "log_intensity_phi_cor")]), row.names = FALSE)
    cat("Global replicates:", min(ids), "-", max(ids),
      "| Full study size:", config$n_replicates,
      "| Physical fits:", planned_fit_count(config, plan), "\n")
    batches <- if (is.null(run_id)) list(ids) else split(ids, (seq_along(ids) - 1L) %/% batch_size)
    for (batch in batches) {
      batch_config <- config
      batch_config$replicate_ids <- unname(batch)
      run_dir <- NULL
      result_key <- selected
      if (!is.null(run_id)) {
        output_root <- option("output-root", file.path(config$project_dir, "results",
          sub("^run_", "", tools::file_path_sans_ext(basename(entry_script)))))
        if ("--test" %in% args) output_root <- file.path(output_root, "test")
        batch_id <- sprintf("batch_%03d_%03d", min(batch), max(batch))
        run_dir <- file.path(path.expand(output_root), run_id, laptop, batch_id, selected)
        batch_config$batch <- list(run_id = run_id, laptop = laptop, batch_id = batch_id)
        batch_config$output_dir <- run_dir
        result_key <- paste(selected, batch_id, sep = "/")
      }
      cat("Output:", if (is.null(run_dir)) batch_config$output_dir else run_dir, "\n")
      if ("--run" %in% args) {
        results[[result_key]] <- run_simulation(batch_config, engine, shared, run_dir)
      } else if ("--verify-only" %in% args) {
        saved <- readRDS(file.path(run_dir, "config.rds"))
        packages <- engine$configure_runtime(batch_config)
        versions <- runner_package_versions(packages)
        expected <- run_identity(batch_config, runner_source_hashes(), versions)
        if (!identical(saved$identity, expected) ||
            !identical(saved$plan, selected_scenario_plan(batch_config, engine))) {
          stop("Verification settings, source code or runtime differ: ", run_dir)
        }
        release_lock <- acquire_run_lock(run_dir)
        status <- tryCatch(verify_completed_run(run_dir), finally = release_lock())
        results[[result_key]] <- list(directory = run_dir, status = status)
      } else {
        results[[result_key]] <- batch_config
      }
      if (!is.null(run_dir) && any(c("--run", "--verify-only") %in% args)) {
        # One record per laptop/batch/experiment avoids concurrent CSV appends
        manifest <- readRDS(file.path(run_dir, "batch_manifest.rds"))
        records_dir <- file.path(path.expand(output_root), run_id, "manifests")
        dir.create(records_dir, recursive = TRUE, showWarnings = FALSE)
        record <- data.frame(run_id = run_id, laptop = laptop, batch_id = batch_id,
          experiment = selected, replicate_start = min(batch), replicate_end = max(batch),
          comparisons = manifest$comparisons, scored = manifest$scored,
          failures = manifest$failures, needs_review = manifest$needs_review,
          status = manifest$status,
          batch_directory = file.path(laptop, batch_id, selected),
          batch_manifest_sha256 = digest::digest(file = file.path(run_dir, "batch_manifest.rds"), algo = "sha256"),
          verification = manifest$verification)
        atomic_write_csv(record, file.path(records_dir,
          paste0(laptop, "_", batch_id, "_", selected, ".csv")))
      }
    }
  }
  if (!any(c("--run", "--verify-only") %in% args)) {
    cat("\nPreview only. Add --run to start fitting\n")
  } else if (any(vapply(results, function(x) any(!x$status$scored), logical(1)))) {
    stop("All available attempts were retained; some jobs failed or were interrupted. Inspect fit_diagnostics.csv; resume does not retry failures")
  }
  invisible(results)
}
