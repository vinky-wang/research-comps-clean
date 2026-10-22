# Controlled diagnostic fits only; original fitted objects are never replaced
# Run from extension_test using the repository's R 4.2 launcher
source('test_helpers.R')
verify_snapshot()
engine$configure_runtime(configs[[1]])
root <- file.path(test_root, 'investigation')
for (dir in c('optimizer/fits', 'optimizer/checkpoints', 'optimizer/logs', 'reports')) {
  dir.create(file.path(root, dir), recursive = TRUE, showWarnings = FALSE)
}
Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
  VECLIB_MAXIMUM_THREADS='1')
ids <- c('strength_case_0003_Melding', 'clustering_case_0003_Melding',
         'strength_case_0003_PSmelding')
plan <- expand.grid(job_id=ids, phase=c('pilot','finer_mesh'),
  setting=c('tight_saved','tight_area_0p001','tight_area_0p1'), stringsAsFactors=FALSE)
plan <- rbind(data.frame(job_id=ids[1], phase='pilot',
  setting=c('default_repeat_1','default_repeat_2')), plan,
  data.frame(job_id=ids[1], phase='pilot', setting='tight_saved_repeat'))
plan$run_id <- paste(plan$phase, plan$job_id, plan$setting, sep='__')
plan$tolerance <- ifelse(grepl('^default',plan$setting), .005, .0001)
plan$initial_area_variance <- ifelse(plan$setting=='tight_area_0p001', .001,
  ifelse(plan$setting=='tight_area_0p1', .1, NA_real_))
write.csv(plan, file.path(root,'optimizer/plan.csv'), row.names=FALSE)
saveRDS(list(control.inla=INLA::inla.set.control.inla.default(),
  control.mode=INLA::inla.set.control.mode.default(),
  env=Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                  'VECLIB_MAXIMUM_THREADS'))), file.path(root,'optimizer/runtime_controls.rds'))
writeLines(capture.output(sessionInfo()),file.path(root,'optimizer/session_info.txt'))

# Hash every original fitted object and source before and after these runs
original_files <- c(list.files('R',full.names=TRUE), 'test_helpers.R',
  list.files('results',pattern='\\.rds$',recursive=TRUE,full.names=TRUE))
integrity_file <- file.path(root,'optimizer/original_md5.rds')
if (!file.exists(integrity_file)) saveRDS(tools::md5sum(original_files),integrity_file)

save_record <- function(record, path) {
  if (file.exists(path)) stop('Refusing to replace a recorded attempt: ',path)
  saveRDS(record,path)
}

one_run <- function(job) {
  checkpoint <- file.path(root,'optimizer/checkpoints',paste0(job$run_id,'.rds'))
  if (file.exists(checkpoint)) return(readRDS(checkpoint))
  files <- job_files(read_job(job$phase,job$job_id))
  b <- readRDS(files$fit); d <- readRDS(files$data); old <- readRDS(files$scores)
  args <- b$inla_arguments
  args$num.threads <- '1:1'
  # Leave approximation strategy, priors, integration rule and VB correction unchanged
  if (!grepl('^default',job$setting)) {
    initial_theta <- b$fit$mode$theta
    if (is.finite(job$initial_area_variance)) {
      area_name <- b$likelihood_map$precision_name[b$likelihood_map$component=='area']
      internal_name <- sub('^Precision','Log precision',area_name)
      stopifnot(internal_name %in% names(initial_theta))
      initial_theta[internal_name] <- -log(job$initial_area_variance)
    }
    args$control.mode <- list(theta=initial_theta, x=b$fit$mode$x,
                              restart=TRUE, fixed=FALSE)
    args$control.inla$tolerance <- .0001
  }
  output <- file.path(root,'optimizer/fits',paste0(job$run_id,'.rds'))
  if (file.exists(output)) stop('Unrecorded fit already exists; preserve and inspect: ',output)
  record <- list(job=job, source=files, seed=b$config$seed,
    control.mode=args$control.mode, control.inla=args$control.inla,
    original_mode=b$fit$mode, fit_file=output)
  warnings <- character(); start <- proc.time()
  message('Starting ',job$run_id)
  fit <- tryCatch(withCallingHandlers(
    engine$with_simulation_seed(b$config$seed, do.call(INLA::inla,args)),
    warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')}
  ),error=function(e)e)
  record$elapsed_seconds <- unname((proc.time()-start)['elapsed'])
  record$warnings <- unique(warnings)
  record$success <- !inherits(fit,'error')
  if (!record$success) {
    record$error <- conditionMessage(fit)
    save_record(record,checkpoint)
    message('Recorded failure: ',record$error)
    return(record)
  }
  # Save the entire new return before extracting summaries
  b$fit <- fit; b$inla_arguments <- args; b$warnings <- unique(warnings)
  b$runtime <- c(elapsed=record$elapsed_seconds)
  b$inla_runtime <- b$runtime
  b$diagnostics$mode_status <- fit$mode$mode.status
  b$diagnostics$inla_cpu_used <- fit$cpu.used
  p <- fit$summary.fitted.values[b$prediction_index,,drop=FALSE]
  for (mapping in list(c('pred_mean','mean'),c('pred_med','0.5quant'),c('sd','sd'),
                      c('pred_ll','0.025quant'),c('pred_ul','0.975quant'))) {
    b$predictions[[mapping[1]]] <- p[,mapping[2]]
  }
  b$investigation <- list(run=job,source=files,seed=b$config$seed)
  saveRDS(b,output)
  writeLines(fit$logfile,file.path(root,'optimizer/logs',paste0(job$run_id,'.txt')))
  record$parameters <- tryCatch(engine$score_parameters(b,d$truth,list(nu=1)),
                                error=function(e)conditionMessage(e))
  record$mode <- fit$mode
  record$prediction_mean <- b$predictions$pred_mean
  record$prediction_rms_change <- sqrt(mean((record$prediction_mean-
    readRDS(files$fit)$predictions$pred_mean)^2))
  record$log_posterior_change <- fit$mode$log.posterior.mode-record$original_mode$log.posterior.mode
  neg <- grep('number of negative eigenvalues=',fit$logfile,value=TRUE)
  neg <- as.integer(sub('.*number of negative eigenvalues=\\[([0-9]+)\\].*','\\1',neg))
  record$max_negative_eigenvalues <- if(length(neg)) max(neg,na.rm=TRUE) else NA_integer_
  record$covariance_eigenvalues <- eigen(fit$misc$cov.intern,symmetric=TRUE,only.values=TRUE)$values
  record$covariance <- fit$misc$cov.intern
  record$nfunc <- fit$misc$nfunc
  record$original_parameters <- old$parameters
  record$boundary <- lapply(c('point','area'),function(component) {
    nm <- b$likelihood_map$precision_name[b$likelihood_map$component==component]
    marginal <- fit$internal.marginals.hyperpar[[sub('^Precision','Log precision',nm)]]
    thresholds <- c(.0001,.001,.01,.025,.1)
    data.frame(component=component,threshold=thresholds,
      probability_below=1-INLA::inla.pmarginal(-log(thresholds),marginal))
  })
  save_record(record,checkpoint)
  message('Completed ',job$run_id,': mode=',fit$mode$mode.status,
    ', negative Hessian eigenvalues=',record$max_negative_eigenvalues,
    ', seconds=',round(record$elapsed_seconds,1))
  record
}

write_reports <- function(records) {
  rows <- lapply(records,function(r) {
    data.frame(run_id=r$job$run_id,phase=r$job$phase,job_id=r$job$job_id,
      setting=r$job$setting,success=r$success,
      mode_status=if(r$success) r$mode$mode.status else NA,
      max_negative_eigenvalues=if(r$success)r$max_negative_eigenvalues else NA,
      log_posterior_mode=if(r$success)r$mode$log.posterior.mode else NA,
      log_posterior_change=if(r$success)r$log_posterior_change else NA,
      prediction_rms_change=if(r$success)r$prediction_rms_change else NA,
      elapsed_seconds=r$elapsed_seconds,
      warnings=paste(r$warnings,collapse=' | '),
      error=if(is.null(r$error)) '' else r$error)
  })
  write.csv(do.call(rbind,rows),file.path(root,'reports/optimizer_runs.csv'),row.names=FALSE)
  params <- lapply(records,function(r) if(isTRUE(r$success) && is.data.frame(r$parameters)) {
    x <- r$parameters
    x$run_id <- r$job$run_id
    x$original_mean <- r$original_parameters$posterior_mean[match(x$parameter,r$original_parameters$parameter)]
    x$original_sd <- r$original_parameters$posterior_sd[match(x$parameter,r$original_parameters$parameter)]
    x$shift_in_original_sd <- (x$posterior_mean-x$original_mean)/x$original_sd
    x
  })
  write.csv(do.call(rbind,params),file.path(root,'reports/optimizer_parameters.csv'),row.names=FALSE)
  boundary <- lapply(records,function(r) if(isTRUE(r$success)) {
    x <- do.call(rbind,r$boundary); x$run_id <- r$job$run_id; x
  })
  write.csv(do.call(rbind,boundary),file.path(root,'reports/optimizer_boundary.csv'),row.names=FALSE)
}

records <- list()
for (i in seq_len(nrow(plan))) {
  records[[i]] <- one_run(plan[i,,drop=FALSE])
  write_reports(records)
  gc()
}
expected <- readRDS(integrity_file)
unchanged <- identical(unname(expected),unname(tools::md5sum(names(expected))))
writeLines(paste('All original fitted objects, data, results and snapshot R scripts unchanged:',unchanged),
  file.path(root,'reports/original_integrity.txt'))
stopifnot(unchanged)
cat('Completed ',length(records),' recorded attempts\n')
