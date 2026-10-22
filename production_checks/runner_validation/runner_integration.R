work <- normalizePath(getwd(), mustWork = TRUE)
setwd(work)
source("R/run_extension.R")
engine <- load_simulation()
out <- file.path(work, "production_checks", "runner_validation", paste0("outputs_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out, recursive = TRUE)
checks <- character()
check <- function(value, label) {
  stopifnot(isTRUE(value))
  checks <<- c(checks, label)
  cat("PASS:", label, "\n")
}
expect_error <- function(code, pattern, label) {
  message <- tryCatch({ force(code); NA_character_ }, error = conditionMessage)
  check(!is.na(message) && grepl(pattern, message), label)
}

# Partition the existing complete 100-replicate plans, keeping all case IDs
configs <- lapply(names(extension_settings$experiments), extension_config, project_dir = work)
names(configs) <- names(extension_settings$experiments)
counts <- c(0L, 0L)
for (name in names(configs)) {
  config <- configs[[name]]
  full <- engine$scenario_plan(config)
  halves <- lapply(list(1:50, 51:100), function(ids) {
    config$replicate_ids <- ids
    plan <- selected_scenario_plan(config, engine)
    check(identical(plan, full[full$replication %in% ids, , drop = FALSE]),
      paste(name, "full-plan IDs retained", min(ids)))
    for (component in c("locations", "fit", "scores")) {
      selected_seeds <- vapply(seq_len(nrow(plan)), function(i)
        engine$sampling_seed(config, plan[i, , drop = FALSE], component, "PSmelding"), integer(1))
      original <- full[match(plan$case_id, full$case_id), , drop = FALSE]
      full_seeds <- vapply(seq_len(nrow(original)), function(i)
        engine$sampling_seed(config, original[i, , drop = FALSE], component, "PSmelding"), integer(1))
      check(identical(selected_seeds, full_seeds), paste(name, component, "seeds retained", min(ids)))
    }
    plan
  })
  counts <- counts + vapply(halves, function(plan) planned_fit_count(config, plan), numeric(1))
  check(!length(intersect(halves[[1]]$case_id, halves[[2]]$case_id)), paste(name, "halves disjoint"))
}
check(identical(as.integer(counts), c(1450L, 1450L)), "1,450 physical fits per laptop")

preview <- NULL
capture.output(preview <- study_main(c("--section=all", "--replicate-start=51", "--replicate-end=100",
  "--run-id=integration", "--laptop=laptop_2", "--batch-size=5", paste0("--output-root=", out)),
  extension_config, normalizePath("R/run_extension.R"), sections = names(configs), default_section = "all"))
check(length(preview) == 40L, "named CLI creates 10 planned batches per experiment")
check(all(vapply(preview, function(x) x$n_replicates == 100L && length(x$replicate_ids) == 5L &&
  min(x$replicate_ids) >= 51L && max(x$replicate_ids) <= 100L, logical(1))), "CLI keeps full study size and global ranges")
check(!dir.exists(file.path(out, "integration")), "preview writes no run artifacts")
expect_error(study_main(c("--run-id=bad", "--laptop=laptop_1", "--replicate-start=101", "--replicate-end=102"),
  extension_config, normalizePath("R/run_extension.R"), sections = names(configs), default_section = "strength"), "replicate_ids", "out-of-range CLI rejected")
expect_error(study_main(c("--run-id=bad", "--laptop=laptop_1"), extension_config,
  normalizePath("R/run_extension.R"), sections = names(configs), default_section = "strength"),
  "explicit replicate-start", "named production run requires explicit global range")
expect_error(study_main(c("--run-id=bad", "--laptop=laptop_1", "--replicate-start=1", "--replicate-end=2147483647"),
  extension_config, normalizePath("R/run_extension.R"), sections = names(configs), default_section = "strength"),
  "replicate_ids", "huge out-of-range end rejected before sequence allocation")
capture.output(test_preview <- study_main(c("--test", "--run-id=small", "--laptop=laptop_1"),
  extension_config, normalizePath("R/run_extension.R"), sections = names(configs), default_section = "strength"))
check(length(test_preview) == 1L && test_preview[[1]]$n_replicates == 1L &&
  identical(test_preview[[1]]$replicate_ids, 1L), "named test preview safely selects its own one-replicate range")

# Real production mesh and prior initialization, fake model operations
config <- configs$strength
config$generating$gamma <- 1
config$n_replicates <- 2L
config$replicate_ids <- 2L
config$models <- "Melding"
config$output_dir <- out
packages <- engine$configure_runtime(config)
mesh <- engine$build_mesh(config)
check(mesh$n > 1000L, "real production mesh initialized")
state <- new.env(parent = emptyenv())
state$fits <- state$scores <- state$datasets <- 0L
state$mode <- "success"
fake <- load_simulation()
fake$build_mesh <- function(config) mesh
fake$generate_surface <- function(config, replication) list(replication = replication)
fake$generate_sampling_fields <- function(config, replication) list()
fake$generate_dataset <- function(config, scenario, field) {
  state$datasets <- state$datasets + 1L
  list(scenario = scenario, truth = c(0, 1), seeds = list(surface = fake$simulation_seed(config, "surface", scenario$rho, scenario$replication)))
}
fake$fit_model <- function(model, dataset, mesh, config, on_inla_return) {
  state$fits <- state$fits + 1L
  if (state$mode == "fit_failure") stop("Deliberate test fit failure")
  if (state$mode == "saved_warning") warning("Preserve this original fitting warning")
  fit <- list(logfile = character())
  on_inla_return(fit, list(model = model))
  if (state$mode == "interrupt_raw") stop(structure(list(message = "Deliberate interruption after raw return", call = NULL), class = c("interrupt", "condition")))
  list(fit = fit, predictions = data.frame(mean = c(0, 1)),
    diagnostics = list(mode_status = 0), warnings = character(), inla_runtime = .1)
}
fake$evaluate_fit <- function(bundle, truth, config, score_seed) {
  state$scores <- state$scores + 1L
  if (state$mode == "interrupt_score") stop(structure(list(message = "Deliberate post-fit interruption", call = NULL), class = c("interrupt", "condition")))
  list(metrics = data.frame(region = "all", metric = "MSE", value = 0, valid = TRUE, needs_review = FALSE),
    parameters = data.frame(parameter = "mu", mean = 0), pointwise = data.frame(error = c(0, 0)),
    diagnostics = list(mode_status = 0, needs_review = FALSE))
}
file_hashes <- function(path) {
  files <- sort(list.files(path, recursive = TRUE, full.names = TRUE))
  setNames(vapply(files, function(f) digest::digest(file = f, algo = "sha256"), character(1)), substring(files, nchar(path) + 2L))
}

success_dir <- file.path(out, "success")
first <- run_simulation(config, fake, run_dir = success_dir)
check(nrow(first$status) == 1L && all(first$status$scored) && state$fits == 1L, "successful job checkpoints and finalizes")
check(identical(first$status$elapsed_scope, "Complete attempt in this session"), "uninterrupted job elapsed scope is explicit")
frozen <- readRDS(file.path(success_dir, "config.rds"))
check(identical(frozen$fitting_setup$priors$gaussian_precision$prior, "pc.prec") &&
  identical(frozen$fitting_setup$priors$gaussian_precision$param, c(.316, .5)), "actual prior B frozen")
before <- file_hashes(success_dir)
run_simulation(config, fake, run_dir = success_dir)
check(identical(before, file_hashes(success_dir)) && state$fits == 1L && state$scores == 1L, "completed resume leaves every artifact byte-identical")
changed <- config
changed$base_seed <- changed$base_seed + 1L
expect_error(run_simulation(changed, fake, run_dir = success_dir), "differ", "changed configuration rejected")
source_fixture <- file.path(out, "changed_source")
dir.create(source_fixture)
file.copy(list.files(file.path(work, "R"), full.names = TRUE), source_fixture)
cat("\n# test-only source mutation\n", file = file.path(source_fixture, "config.R"), append = TRUE)
original_directory <- .runner_directory
.runner_directory <- source_fixture
expect_error(run_simulation(config, fake, run_dir = success_dir), "differ", "changed source rejected")
.runner_directory <- original_directory
release <- acquire_run_lock(success_dir)
expect_error(run_simulation(config, fake, run_dir = success_dir), "locked", "second writer rejected")
release()

failure_dir <- file.path(out, "failure")
state$mode <- "fit_failure"
failed <- run_simulation(config, fake, run_dir = failure_dir)
check(!failed$status$scored && failed$status$failed_stage == "fitting", "fit failure retained as terminal status")
before <- file_hashes(failure_dir)
fit_count <- state$fits
state$mode <- "success"
run_simulation(config, fake, run_dir = failure_dir)
check(state$fits == fit_count && identical(before, file_hashes(failure_dir)), "failed attempt never silently retried or replaced")

interrupt_dir <- file.path(out, "interrupted_after_fit")
state$mode <- "interrupt_score"
interrupted <- tryCatch({run_simulation(config, fake, run_dir = interrupt_dir); FALSE}, interrupt = function(c) TRUE)
check(interrupted && !dir.exists(file.path(interrupt_dir, ".lock")) &&
  file.exists(file.path(interrupt_dir, "inla_fits", "case_0002_Melding.rds")), "interruption preserves complete fit and releases lock")
fit_count <- state$fits
fit_hash <- digest::digest(file = file.path(interrupt_dir, "inla_fits", "case_0002_Melding.rds"), algo = "sha256")
state$mode <- "success"
resumed <- run_simulation(config, fake, run_dir = interrupt_dir)
check(all(resumed$status$scored) && state$fits == fit_count &&
  identical(fit_hash, digest::digest(file = file.path(interrupt_dir, "inla_fits", "case_0002_Melding.rds"), algo = "sha256")), "post-fit interruption resumes scoring without refitting")
check(grepl("Resume session only", resumed$status$elapsed_scope), "resumed job timing is labelled as partial session")

corrupt_data_dir <- file.path(out, "corrupt_incomplete_dataset")
state$mode <- "interrupt_score"
interrupted <- tryCatch({run_simulation(config, fake, run_dir = corrupt_data_dir); FALSE}, interrupt = function(c) TRUE)
check(interrupted && file.exists(file.path(corrupt_data_dir, "simulation_data", "case_0002_manifest.rds")),
  "dataset hashes exist before an interrupted job completes")
cat("corrupt dataset", file = file.path(corrupt_data_dir, "simulation_data", "case_0002.rds"), append = TRUE)
state$mode <- "success"
fit_count <- state$fits
score_count <- state$scores
expect_error(run_simulation(config, fake, run_dir = corrupt_data_dir), "verification failed", "corrupted dataset in incomplete job is rejected before loading")
check(state$fits == fit_count && state$scores == score_count && !dir.exists(file.path(corrupt_data_dir, ".lock")),
  "incomplete dataset corruption causes no fit or scoring and releases the lock")

diagnostics_dir <- file.path(out, "interrupted_diagnostics")
original_finish <- finish_job_checkpoint
finish_job_checkpoint <- function(entry, run_dir, folder) {
  if (basename(run_dir) == "interrupted_diagnostics") stop(structure(list(message = "Interrupt after saved diagnostics", call = NULL), class = c("interrupt", "condition")))
  original_finish(entry, run_dir, folder)
}
state$mode <- "saved_warning"
interrupted <- tryCatch({run_simulation(config, fake, run_dir = diagnostics_dir); FALSE}, interrupt = function(c) TRUE)
finish_job_checkpoint <- original_finish
diagnostics_file <- file.path(diagnostics_dir, "diagnostics", "case_0002_Melding", "diagnostics.rds")
check(interrupted && file.exists(diagnostics_file), "diagnostics survive interruption before terminal checkpoint")
original_diagnostics <- readRDS(diagnostics_file)
original_diagnostics_hash <- digest::digest(file = diagnostics_file, algo = "sha256")
fit_count <- state$fits
score_count <- state$scores
state$mode <- "success"
diagnostics_resume <- run_simulation(config, fake, run_dir = diagnostics_dir)
check(state$fits == fit_count && state$scores == score_count &&
  identical(original_diagnostics_hash, digest::digest(file = diagnostics_file, algo = "sha256")) &&
  identical(diagnostics_resume$status$elapsed_seconds, original_diagnostics$elapsed_seconds) &&
  identical(diagnostics_resume$status$elapsed_scope, original_diagnostics$elapsed_scope) &&
  grepl("Preserve this original fitting warning", diagnostics_resume$status$warnings),
  "resuming after diagnostics preserves original timing, scope and warnings without fitting or scoring")

raw_dir <- file.path(out, "interrupted_raw_fit")
state$mode <- "interrupt_raw"
interrupted <- tryCatch({run_simulation(config, fake, run_dir = raw_dir); FALSE}, interrupt = function(c) TRUE)
raw_file <- file.path(raw_dir, "inla_fits", "case_0002_Melding_raw.rds")
check(interrupted && file.exists(raw_file), "raw INLA return survives interrupted fit postprocessing")
fit_count <- state$fits
raw_hash <- digest::digest(file = raw_file, algo = "sha256")
state$mode <- "success"
raw_resume <- run_simulation(config, fake, run_dir = raw_dir)
check(state$fits == fit_count && !raw_resume$status$scored && raw_resume$status$failed_stage == "interrupted" &&
  identical(raw_hash, digest::digest(file = raw_file, algo = "sha256")), "incomplete fit is retained for review without automatic retry")

lock_dir <- file.path(out, "lock_cases")
dir.create(file.path(lock_dir, ".lock"), recursive = TRUE)
saveRDS(list(hostname = "different-test-host", pid = 2147483647L), file.path(lock_dir, ".lock", "owner.rds"))
expect_error(acquire_run_lock(lock_dir), "locked", "unverified remote lock cannot be reclaimed")
saveRDS(list(hostname = unname(Sys.info()["nodename"]), pid = 2147483647L), file.path(lock_dir, ".lock", "owner.rds"))
release <- acquire_run_lock(lock_dir)
check(dir.exists(file.path(lock_dir, "lock_history")) && length(list.files(file.path(lock_dir, "lock_history"), recursive = TRUE)) >= 1L,
  "dead local lock is archived before recovery")
release()

psgeo <- config
psgeo$models <- "PSgeo"
psgeo$generating$J <- c(4L, 25L)
fit_count <- state$fits
psgeo_dir <- file.path(out, "psgeo")
reused <- run_simulation(psgeo, fake, run_dir = psgeo_dir)
check(nrow(reused$status) == 2L && sum(reused$status$fit_reused) == 1L && state$fits == fit_count + 1L && all(reused$status$scored), "PSgeo is fitted once across areal counts")
before <- file_hashes(psgeo_dir)
run_simulation(psgeo, fake, run_dir = psgeo_dir)
check(identical(before, file_hashes(psgeo_dir)), "PSgeo reused checkpoint verifies unchanged")

# Execute the shared CLI through its real batch loop, substituting only the engine
original_loader <- load_simulation
load_simulation <- function(directory = .runner_directory) fake
factory <- function(section, project_dir) {
  value <- config
  value$section <- section
  value$project_dir <- project_dir
  value
}
cli_args <- c("--section=strength", "--replicate-start=1", "--replicate-end=2", "--batch-size=1",
  "--run-id=cli_run", "--laptop=laptop_1", paste0("--output-root=", out))
fit_count <- state$fits
capture.output(cli_runs <- study_main(c(cli_args, "--run"), factory, normalizePath("R/run_extension.R"), sections = "strength", default_section = "strength"))
check(length(cli_runs) == 2L && state$fits == fit_count + 2L &&
  all(vapply(cli_runs, function(x) all(x$status$scored), logical(1))), "named CLI runs both requested global batches")
cli_root <- file.path(out, "cli_run")
records_files <- list.files(file.path(cli_root, "manifests"), full.names = TRUE)
records <- do.call(rbind, lapply(records_files, read.csv, stringsAsFactors = FALSE))
check(length(records_files) == 2L && nrow(records) == 2L &&
  all(records$run_id == "cli_run" & records$laptop == "laptop_1" & records$experiment == "strength" &
    records$status == "complete" & records$comparisons == 1L & records$scored == 1L & records$failures == 0L) &&
  identical(sort(records$replicate_start), 1:2) && identical(records$replicate_start, records$replicate_end),
  "central manifests contain one complete record for each laptop batch experiment")
check(all(grepl("Local files checked", records$verification) & grepl("Drive upload not verified", records$verification)) &&
  all(vapply(seq_len(nrow(records)), function(i) identical(records$batch_manifest_sha256[i],
    digest::digest(file = file.path(cli_root, records$batch_directory[i], "batch_manifest.rds"), algo = "sha256")), logical(1))),
  "central records bind saved batch manifests and claim only local verification")
before <- file_hashes(cli_root)
capture.output(cli_verified <- study_main(c(cli_args, "--verify-only"), factory, normalizePath("R/run_extension.R"), sections = "strength", default_section = "strength"))
check(identical(before, file_hashes(cli_root)) && length(cli_verified) == 2L && state$fits == fit_count + 2L,
  "CLI verify-only checks completed batches without changing files or fitting")
capture.output(cli_resumed <- study_main(c(cli_args, "--run"), factory, normalizePath("R/run_extension.R"), sections = "strength", default_section = "strength"))
check(identical(before, file_hashes(cli_root)) && length(cli_resumed) == 2L && state$fits == fit_count + 2L,
  "repeated named CLI command resumes both batches unchanged")
summary <- read.csv(file.path(cli_runs[[1]]$directory, "predictive_score_summary.csv"))
check(all(summary$n_requested == 1L), "batch summaries count selected replicates only")
load_simulation <- original_loader

corrupt <- file.path(success_dir, first$status$prediction_file)
cat("corruption", file = corrupt, append = TRUE)
expect_error(run_simulation(config, fake, run_dir = success_dir), "verification failed", "corrupted saved prediction rejected")
check(!dir.exists(file.path(success_dir, ".lock")), "verification failure releases writer lock")

writeLines(c(paste(length(checks), "checks passed"), checks, paste("Artifacts:", out)), file.path(work, "production_checks", "runner_validation", "results.txt"))
cat("\nALL", length(checks), "CHECKS PASSED\nArtifacts:", out, "\n")
