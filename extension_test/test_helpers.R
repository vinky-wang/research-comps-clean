# Isolated test driver; the copied R modules are an unchanged source snapshot
source('R/run_extension.R', local = TRUE)
test_root <- normalizePath(getwd())
engine <- load_simulation()
sections <- names(extension_settings$experiments)
configs <- setNames(lapply(sections, function(section) {
  cfg <- extension_config(section, test_root)
  cfg$n_replicates <- 1L
  cfg
}), sections)
shared <- new.env(parent = emptyenv())
for (directory in c('results', 'checkpoints', 'reports')) dir.create(directory, showWarnings = FALSE)

verify_snapshot <- function() {
  expected <- readRDS('source_snapshot.rds')
  actual <- tools::md5sum(names(expected))
  if (!identical(unname(expected), unname(actual))) stop('The source snapshot changed; use a new test folder')
}

job_plan <- function() {
  do.call(rbind, lapply(sections, function(section) {
    plan <- engine$scenario_plan(configs[[section]])
    do.call(rbind, lapply(seq_len(nrow(plan)), function(i) {
      do.call(rbind, lapply(configs[[section]]$models, function(model) {
        data.frame(plan[i,], model = model,
          job_id = paste(section, plan$case_id[i], model, sep = '_'))
      }))
    }))
  }))
}

# Each physical fit has its own durable completion record
# Restarting skips completed attempts, including failures; no selective retries
run_job <- function(job, cfg, phase) {
  id <- paste(phase, job$job_id, sep = '_')
  checkpoint <- file.path('checkpoints', paste0(id, '.rds'))
  if (file.exists(checkpoint)) {
    message('Already recorded: ', id)
    return(readRDS(checkpoint))
  }
  cfg$models <- job$model
  cfg$output_dir <- file.path(test_root, 'results', phase, job$job_id)
  cfg$entry_script <- file.path(test_root, 'R/run_extension.R')
  # Keep the original case ID and replicate ID while executing one job
  scenario <- engine$scenario_plan(cfg)
  scenario <- scenario[scenario$case_id == job$case_id, , drop = FALSE]
  stopifnot(nrow(scenario) == 1L)
  local_engine <- load_simulation()
  local_engine$scenario_plan <- function(config) scenario
  previous <- list.files(cfg$output_dir, pattern = paste0('^', job$case_id, '_', job$model, '\\.rds$'),
    recursive = TRUE, full.names = TRUE)
  previous <- previous[grepl('/inla_fits/', previous, fixed = TRUE)]
  recovered <- length(previous) > 0L
  if (recovered) {
    # An interrupted scoring/output stage can reuse the complete saved fit
    prior_fit <- tail(sort(previous), 1L)
    local_engine$fit_model <- function(...) readRDS(prior_fit)
  } else if (dir.exists(cfg$output_dir) && length(list.files(cfg$output_dir, recursive = TRUE))) {
    result <- list(job = job, phase = phase, directory = NA_character_, status = NULL,
      error = 'Interrupted attempt without a complete fit; inspect retained raw checkpoints before explicitly retrying',
      recovered_fit = FALSE, elapsed_seconds = 0)
    save_new_rds(result, checkpoint)
    return(result)
  }
  disk <- system2('df', c('-k', shQuote(test_root)), stdout = TRUE)
  available <- as.numeric(strsplit(trimws(tail(disk, 1L)), '[[:space:]]+')[[1]][4])
  if (is.finite(available) && available < 5 * 1024^2) {
    stop('Less than 5 GiB free; preserved outputs must be moved before continuing')
  }
  started <- proc.time()
  result <- tryCatch({
    run <- run_simulation(cfg, local_engine, shared)
    list(job = job, phase = phase, directory = run$directory, status = run$status,
      error = NULL, recovered_fit = recovered)
  }, error = function(error) list(job = job, phase = phase, directory = NA_character_,
    status = NULL, error = conditionMessage(error), recovered_fit = recovered))
  result$elapsed_seconds <- unname((proc.time() - started)['elapsed'])
  save_new_rds(result, checkpoint)
  message('Finished ', id, ' in ', round(result$elapsed_seconds, 1), ' seconds')
  gc()
  result
}

run_pilot <- function() {
  plan <- job_plan()
  stopifnot(nrow(plan) == 29L)
  if (!file.exists('pilot_plan.csv')) write.csv(plan, 'pilot_plan.csv', row.names = FALSE)
  for (i in seq_len(nrow(plan))) {
    run_job(plan[i,], configs[[plan$experiment[i]]], 'pilot')
  }
}

selected_jobs <- function() {
  plan <- job_plan()
  plan[(plan$experiment == 'strength' & plan$gamma == 1 & !is.na(plan$gamma)) |
    (plan$experiment == 'spatial_gamma' & plan$sampling_id == 'sign_reversal') |
    (plan$experiment == 'clustering' & plan$delta == 2), , drop = FALSE]
}

read_job <- function(phase, id) {
  path <- file.path('checkpoints', paste0(phase, '_', id, '.rds'))
  if (file.exists(path)) readRDS(path) else NULL
}

job_files <- function(record) {
  if (is.null(record) || is.null(record$status) || !isTRUE(record$status$fit_saved)) return(NULL)
  x <- record$status
  if (is.na(x$fit_file) || grepl('_raw.rds$', x$fit_file)) return(NULL)
  list(fit = file.path(record$directory, x$fit_file),
    data = file.path(record$directory, x$simulation_file),
    scores = file.path(record$directory, x$diagnostics_directory, 'scores.rds'))
}

run_stability <- function() {
  chosen <- selected_jobs()
  stopifnot(nrow(chosen) == 5L, !anyNA(chosen$job_id))
  for (i in seq_len(nrow(chosen))) {
    job <- chosen[i,]
    base <- read_job('pilot', job$job_id)
    files <- job_files(base)
    if (is.null(files)) {
      message('No complete pilot fit for stability check: ', job$job_id)
      next
    }
    cfg <- configs[[job$experiment]]
    cfg$mesh$max_edge <- cfg$mesh$max_edge * .75
    cfg$mesh$cutoff <- cfg$mesh$cutoff * .75
    run_job(job, cfg, 'finer_mesh')
    checkpoint <- file.path('checkpoints', paste0('crps_', job$job_id, '.rds'))
    if (file.exists(checkpoint)) next
    started <- proc.time()
    result <- tryCatch({
      bundle <- readRDS(files$fit)
      data <- readRDS(files$data)
      old <- readRDS(files$scores)
      draws <- engine$draw_latent_predictions(bundle, n = 5000L, seed = old$metadata$seed)
      values <- engine$sample_crps(draws$Z, data$truth$Z)
      list(job = job, error = NULL, n_draws = 5000L, seed = old$metadata$seed,
        pointwise_crps = values, mean_crps = mean(values),
        pilot_crps = old$metrics$value[old$metrics$metric == 'CRPS_Z_joint' & old$metrics$region == 'all'])
    }, error = function(error) list(job = job, error = conditionMessage(error)))
    result$elapsed_seconds <- unname((proc.time() - started)['elapsed'])
    save_new_rds(result, checkpoint)
    message('Saved high-draw CRPS check: ', job$job_id)
    gc()
  }
}

# Keep failed fits visible and distinguish completion from numerical validity
write_reports <- function() {
  paths <- list.files('checkpoints', pattern = '^(pilot|finer_mesh)_.*\\.rds$', full.names = TRUE)
  records <- lapply(paths, readRDS)
  if (!length(records)) return(invisible(NULL))
  rows <- lapply(records, function(record) {
    status <- record$status
    files <- job_files(record)
    size <- if (!is.na(record$directory)) {
      sum(file.info(list.files(record$directory, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE)
    } else 0
    metrics <- params <- NULL
    if (!is.null(files) && file.exists(files$scores)) {
      score <- readRDS(files$scores)
      metrics <- score$metrics
      params <- score$parameters
    }
    data.frame(phase = record$phase, job_id = record$job$job_id,
      experiment = record$job$experiment, model = record$job$model,
      fit_saved = !is.null(status) && isTRUE(status$fit_saved),
      scored = !is.null(status) && isTRUE(status$scored),
      needs_review = is.null(status) || isTRUE(status$needs_review),
      mode_status = if (is.null(status)) NA_real_ else status$mode_status,
      hessian_adjustment = if (is.null(status)) NA else status$hessian_adjustment_reported,
      warnings = if (is.null(status)) '' else status$warnings,
      invalid_scores = if (is.null(metrics)) NA_integer_ else sum(!metrics$valid),
      invalid_parameter_summaries = if (is.null(params)) NA_integer_ else sum(!params$summary_valid),
      parameter_coherence_flags = if (is.null(params)) NA_integer_ else sum(params$coherence_warning),
      elapsed_seconds = record$elapsed_seconds, bytes = size,
      recovered_fit = record$recovered_fit,
      error = if (is.null(record$error)) {if (is.null(status)) '' else status$message} else record$error,
      directory = record$directory)
  })
  status <- do.call(rbind, rows)
  write.csv(status, 'reports/job_diagnostics.csv', row.names = FALSE)
  pilot <- status[status$phase == 'pilot',]
  estimate <- data.frame(completed_pilot_jobs = nrow(pilot), expected_pilot_jobs = 29,
    pilot_seconds = sum(pilot$elapsed_seconds), pilot_bytes = sum(pilot$bytes),
    projected_full_hours = sum(pilot$elapsed_seconds) * 100 / 3600,
    projected_full_GiB = sum(pilot$bytes) * 100 / 1024^3,
    note = 'Linear extrapolation from one replicate; includes per-job test overhead; not an upper bound')
  write.csv(estimate, 'reports/runtime_storage_estimate.csv', row.names = FALSE)
  # Compare paired predictions and summaries on the same generated datasets
  comparisons <- list(); parameter_changes <- list()
  for (id in selected_jobs()$job_id) {
    b <- job_files(read_job('pilot', id)); f <- job_files(read_job('finer_mesh', id))
    if (is.null(b) || is.null(f) || !file.exists(b$scores) || !file.exists(f$scores)) next
    base <- readRDS(b$fit); fine <- readRDS(f$fit)
    bs <- readRDS(b$scores); fs <- readRDS(f$scores)
    stopifnot(identical(readRDS(b$data)$points, readRDS(f$data)$points))
    rms <- sqrt(mean((base$predictions$pred_mean - fine$predictions$pred_mean)^2))
    get_crps <- function(s) s$metrics$value[s$metrics$metric == 'CRPS_Z_joint' & s$metrics$region == 'all']
    change <- get_crps(fs) - get_crps(bs)
    comparisons[[id]] <- data.frame(job_id = id, base_vertices = base$spde$n.spde,
      finer_vertices = fine$spde$n.spde, prediction_mean_rms_change = rms,
      joint_crps_change = change, prediction_screen_pass = rms <= .05,
      crps_screen_pass = is.finite(change) && abs(change) <= .01)
    wanted <- c('mu','gamma','tau2_s','tau2_B','rho','sigma2','theta')
    a <- bs$parameters[bs$parameters$parameter %in% wanted,]
    z <- fs$parameters[match(a$parameter, fs$parameters$parameter),]
    shift <- (z$posterior_mean-a$posterior_mean)/a$posterior_sd
    parameter_changes[[id]] <- data.frame(job_id = id, parameter = a$parameter,
      base_mean = a$posterior_mean, finer_mean = z$posterior_mean,
      shift_in_base_posterior_sd = shift, screen_pass = is.finite(shift) & abs(shift) <= .25)
  }
  if(length(comparisons)) write.csv(do.call(rbind, comparisons), 'reports/mesh_comparison.csv', row.names = FALSE)
  if(length(parameter_changes)) write.csv(do.call(rbind, parameter_changes), 'reports/mesh_parameter_changes.csv', row.names = FALSE)
  crps <- lapply(list.files('checkpoints', pattern = '^crps_.*\\.rds$', full.names = TRUE), readRDS)
  if(length(crps)) write.csv(do.call(rbind,lapply(crps,function(x) {
    ok <- is.null(x$error)
    data.frame(job_id = x$job$job_id, n_draws = if(ok) x$n_draws else NA,
      pilot_crps = if(ok) x$pilot_crps else NA, high_draw_crps = if(ok) x$mean_crps else NA,
      change = if(ok) x$mean_crps-x$pilot_crps else NA,
      screen_pass = ok && is.finite(x$mean_crps-x$pilot_crps) && abs(x$mean_crps-x$pilot_crps)<=.01,
      elapsed_seconds = x$elapsed_seconds, error = if(ok) '' else x$error)
  })), 'reports/crps_stability.csv',row.names=FALSE)
  # Paired PSmelding-minus-Melding differences; negative favours PSmelding
  paired <- list()
  plan <- job_plan()
  for(section in c('strength','areal_noise','clustering')) for(case in unique(plan$case_id[plan$experiment==section])) {
    b <- job_files(read_job('pilot',paste(section,case,'Melding',sep='_')))
    p <- job_files(read_job('pilot',paste(section,case,'PSmelding',sep='_')))
    if(is.null(b)||is.null(p)||!file.exists(b$scores)||!file.exists(p$scores)) next
    bm <- readRDS(b$scores)$metrics; pm <- readRDS(p$scores)$metrics
    bidx <- match(paste(pm$region,pm$metric),paste(bm$region,bm$metric))
    paired[[paste(section,case)]] <- data.frame(experiment=section,case_id=case,
      replication=1,region=pm$region,metric=pm$metric,
      Melding=bm$value[bidx],PSmelding=pm$value,difference=pm$value-bm$value[bidx],
      valid_pair=pm$valid & bm$valid[bidx],needs_review=pm$needs_review | bm$needs_review[bidx])
  }
  if(length(paired)) write.csv(do.call(rbind,paired),'reports/paired_model_differences.csv',row.names=FALSE)
  writeLines(c('EXTENSION TEST STATUS',
    paste('Pilot jobs recorded:',nrow(pilot),'of 29'),
    paste('Pilot jobs scored:',sum(pilot$scored)),
    paste('Pilot review flags:',sum(pilot$needs_review)),
    paste('Finer mesh jobs recorded:',sum(status$phase=='finer_mesh'),'of 5'),
    paste('High-draw CRPS checks recorded:',length(crps),'of 5'),
    'See the CSV reports for validity, mesh and CRPS screens, and resource estimates',
    'One replicate is a workflow/numerical pilot, not an estimate of simulation-study performance',
    'Review flags are retained; successful process exit is not a convergence certificate'),
    'reports/status.txt')
}
