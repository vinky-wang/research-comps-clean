# Run from extension_test with the existing R 4.2 launcher
# Each attempt is retained separately; no original result is replaced
source('test_helpers.R')
verify_snapshot()
engine$configure_runtime(configs[[1]])
root <- file.path(test_root, 'prior_sensitivity')
for (dir in c('fits','checkpoints','logs','predictions','posterior_draws','simulation_data','reports')) {
  dir.create(file.path(root,dir),recursive=TRUE,showWarnings=FALSE)
}
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
priors <- list(A=list(prior='loggamma',param=c(1,5e-5)),
               B=list(prior='pc.prec',param=c(.316,.5)),
               C=list(prior='pc.prec',param=c(.5,.5)))
plan <- expand.grid(dataset=c('baseline','clustering_delta2'),model=c('Melding','PSmelding'),
                    prior_id=names(priors),start_area=c(.001,.1,.4),stringsAsFactors=FALSE)
plan$source_id <- paste(ifelse(plan$dataset=='baseline','strength_case_0003','clustering_case_0003'),
                        plan$model,sep='_')
plan$job_id <- paste(plan$dataset,plan$model,plan$prior_id,
                     paste0('area_',gsub('[.]','p',plan$start_area)),sep='__')
stopifnot(nrow(plan)==36L)
write.csv(plan,file.path(root,'plan.csv'),row.names=FALSE)
settings <- list(priors=priors,plan=plan,phase='finer_mesh',point_start=.1,
  hyperparameter_tolerance=1e-4,h=.005,threads='1:1',joint_draws=1000L,
  crps_batches=10L,prior_applies_to=c('point','area'),
  screens=list(log_mode_spread=.01,parameter_shift_pooled_sd=.1,
    prediction_RMS=.01,boundary_probability_spread=.01,
    positive_integration_log_density=.001,predictive_relative_change=.05,
    crps_absolute_change=.01),
  scope='One saved replicate, two datasets, both models, finer mesh only')
saveRDS(settings,file.path(root,'config.rds'))
writeLines(capture.output(sessionInfo()),file.path(root,'session_info.txt'))
saveRDS(list(control.inla=INLA::inla.set.control.inla.default(),
  control.mode=INLA::inla.set.control.mode.default()),file.path(root,'runtime_defaults.rds'))

# Preserve generating truth once per dataset and verify common geometry
reference_mesh <- NULL
for(dataset in unique(plan$dataset)) {
  rows <- plan[plan$dataset==dataset,]
  for(model in unique(rows$model)) {
    id <- rows$source_id[rows$model==model][1]
    files <- job_files(read_job('finer_mesh',id))
    b <- readRDS(files$fit); d <- readRDS(files$data)
    if(is.null(reference_mesh)) reference_mesh <- b$spde$mesh
    stopifnot(identical(reference_mesh$loc,b$spde$mesh$loc),
      identical(reference_mesh$graph$tv,b$spde$mesh$graph$tv))
    path <- file.path(root,'simulation_data',paste0(dataset,'.rds'))
    if(file.exists(path)) {
      previous <- readRDS(path)
      stopifnot(identical(previous$points,d$points),identical(previous$areas,d$areas),
                identical(previous$truth,d$truth))
    } else saveRDS(d,path)
    rm(b,d);gc()
  }
}
saveRDS(reference_mesh,file.path(root,'mesh.rds'));rm(reference_mesh);gc()
original_files <- c(list.files('R',full.names=TRUE),'test_helpers.R',
  list.files('results',pattern='\\.rds$',recursive=TRUE,full.names=TRUE),
  list.files('investigation/optimizer/fits',pattern='\\.rds$',full.names=TRUE))
integrity_file <- file.path(root,'original_md5.rds')
if(!file.exists(integrity_file))saveRDS(tools::md5sum(original_files),integrity_file)

atomic_save <- function(object,path) {
  temporary <- paste0(path,'.tmp')
  saveRDS(object,temporary)
  if(!file.rename(temporary,path)) stop('Could not complete output: ',path)
}
finite_max <- function(x) if(any(is.finite(x)))max(x[is.finite(x)]) else NA_real_

fit_diagnostics <- function(fit) {
  log <- fit$logfile
  neg <- grep('number of negative eigenvalues=',log,value=TRUE,fixed=TRUE)
  neg <- as.numeric(sub('.*number of negative eigenvalues=\\[([0-9]+)\\].*','\\1',neg))
  rel <- grep('config .*log\\(rel.dens\\)=',log,value=TRUE)
  rel <- suppressWarnings(as.numeric(sub('.*log\\(rel.dens\\)= *([-0-9.]+).*','\\1',rel)))
  strategy <- grep('Optimizing pred-strategy',log,value=TRUE,fixed=TRUE)
  list(mode_status=fit$mode$mode.status,max_negative_eigenvalues=finite_max(neg),
    log_posterior_mode=fit$mode$log.posterior.mode,
    max_integration_log_relative_density=finite_max(rel),
    integration_point_exceeds_mode=any(rel>.001,na.rm=TRUE),
    max_stdev_correction=finite_max(c(fit$misc$stdev.corr.negative,fit$misc$stdev.corr.positive)),
    nfunc=fit$misc$nfunc,predictor_strategy=if(length(strategy))
      sub('.*choose\\[([^]]+)\\].*','\\1',tail(strategy,1)) else NA_character_)
}

run_one <- function(job) {
  checkpoint <- file.path(root,'checkpoints',paste0(job$job_id,'.rds'))
  if(file.exists(checkpoint)) {message('Already recorded: ',job$job_id);return(readRDS(checkpoint))}
  disk <- system2('df',c('-k',shQuote(root)),stdout=TRUE)
  available <- as.numeric(strsplit(trimws(tail(disk,1)),'[[:space:]]+')[[1]][4])
  if(is.finite(available)&&available<5*1024^2)stop('Less than 5 GiB free; preserve existing results before continuing')
  files <- job_files(read_job('finer_mesh',job$source_id))
  b <- readRDS(files$fit)
  d <- readRDS(file.path(root,'simulation_data',paste0(job$dataset,'.rds')))
  args <- b$inla_arguments
  controls <- if(is.null(args$control.family))replicate(length(args$family),list(),simplify=FALSE) else args$control.family
  stopifnot(length(controls)==length(args$family))
  initial <- b$fit$mode$theta
  for(component in c('point','area')) {
    map <- b$likelihood_map[b$likelihood_map$component==component,]
    stopifnot(nrow(map)==1L,map$family=='gaussian',!map$fixed)
    variance_start <- if(component=='point') .1 else job$start_area
    internal <- sub('^Precision','Log precision',map$precision_name)
    stopifnot(internal %in% names(initial))
    initial[internal] <- -log(variance_start)
    prior <- priors[[job$prior_id]]
    prior$initial <- -log(variance_start);prior$fixed <- FALSE
    controls[[map$column]]$hyper <- list(prec=prior)
  }
  args$control.family <- controls
  args$control.mode <- list(theta=initial,x=b$fit$mode$x,restart=TRUE,fixed=FALSE)
  args$control.inla$tolerance <- 1e-4
  args$num.threads <- '1:1'
  old_score <- readRDS(files$scores)
  record <- list(job=job,source=files,fit_seed=b$config$seed,
    score_seed=old_score$metadata$seed,control.family=controls,
    control.mode=args$control.mode,control.inla=args$control.inla,
    fit_success=FALSE,score_success=FALSE,error='',warnings=character())
  fit_path <- file.path(root,'fits',paste0(job$job_id,'.rds'))
  record$fit_file <- fit_path
  started <- proc.time()
  message('Starting ',job$job_id)
  if(file.exists(fit_path)) {
    # Recover scoring after an interruption without fitting a second time
    b <- readRDS(fit_path);fit <- b$fit;record$recovered_fit <- TRUE
    record$fit_seconds <- b$prior_sensitivity$fit_seconds
    record$warnings <- b$warnings
  } else {
    fit_start <- proc.time();warnings <- character()
    fit <- tryCatch(withCallingHandlers(
      engine$with_simulation_seed(b$config$seed,do.call(INLA::inla,args)),
      warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
      error=function(e)e)
    record$fit_seconds <- unname((proc.time()-fit_start)['elapsed'])
    record$warnings <- unique(warnings);record$recovered_fit <- FALSE
    if(inherits(fit,'error')) {
      record$error <- conditionMessage(fit)
      record$elapsed_seconds <- unname((proc.time()-started)['elapsed'])
      atomic_save(record,checkpoint);message('Recorded fit failure: ',record$error);return(record)
    }
    b$fit <- fit;b$inla_arguments <- args;b$warnings <- record$warnings
    b$runtime <- b$inla_runtime <- c(elapsed=record$fit_seconds)
    b$diagnostics$mode_status <- fit$mode$mode.status
    b$diagnostics$inla_cpu_used <- fit$cpu.used
    tab <- fit$summary.fitted.values[b$prediction_index,,drop=FALSE]
    for(mapping in list(c('pred_mean','mean'),c('pred_med','0.5quant'),c('sd','sd'),
                        c('pred_ll','0.025quant'),c('pred_ul','0.975quant')))b$predictions[[mapping[1]]] <- tab[,mapping[2]]
    b$prior_sensitivity <- list(job=job,source=files,fit_seconds=record$fit_seconds,
      fit_seed=record$fit_seed,score_seed=record$score_seed)
    atomic_save(b,fit_path)
  }
  record$fit_success <- TRUE
  record$mode <- fit$mode
  record$diagnostics <- fit_diagnostics(fit)
  record$prediction_mean <- b$predictions$pred_mean
  record$prediction_sd <- b$predictions$sd
  writeLines(fit$logfile,file.path(root,'logs',paste0(job$job_id,'.txt')))
  atomic_save(b$predictions,file.path(root,'predictions',paste0(job$job_id,'.rds')))
  record$parameters <- tryCatch(engine$score_parameters(b,d$truth,list(nu=1)),error=function(e)conditionMessage(e))
  record$boundary <- tryCatch(do.call(rbind,lapply(c('point','area'),function(component){
    nm <- b$likelihood_map$precision_name[b$likelihood_map$component==component]
    marginal <- fit$internal.marginals.hyperpar[[sub('^Precision','Log precision',nm)]]
    thresholds <- c(.001,.01,.025,.1)
    data.frame(component=component,threshold=thresholds,
      probability_below=pmin(1,pmax(0,1-INLA::inla.pmarginal(-log(thresholds),marginal))))
  })),error=function(e)conditionMessage(e))
  score <- tryCatch({
    alignment <- engine$check_prediction_alignment(b,d$truth)
    stopifnot(alignment$coordinates_valid,alignment$truth_valid,alignment$mean_valid)
    draw_path <- file.path(root,'posterior_draws',paste0(job$job_id,'.rds'))
    draws <- if(file.exists(draw_path))readRDS(draw_path) else {
      z <- engine$draw_latent_predictions(b,1000L,record$score_seed)
      atomic_save(z,draw_path);z
    }
    crps <- engine$sample_crps(draws$Z,d$truth$Z)
    mse <- (b$predictions$pred_mean-d$truth$Z)^2
    permutation <- engine$with_simulation_seed(record$score_seed,sample.int(ncol(draws$Z)))
    batches <- split(permutation,rep(seq_len(10),each=100))
    batch_crps <- vapply(batches,function(index)
      mean(engine$sample_crps(draws$Z[,index,drop=FALSE],d$truth$Z)),numeric(1))
    stopifnot(all(is.finite(mse)),all(is.finite(crps)))
    list(metrics=data.frame(MSE=mean(mse),CRPS_Z_joint=mean(crps),
      CRPS_MCSE_approx=sd(batch_crps)/sqrt(length(batch_crps)),
      CRPS_Z_gaussian=mean(engine$gaussian_crps(b$predictions$pred_mean,b$predictions$sd,d$truth$Z))),
      batch_crps=batch_crps,
      pointwise=data.frame(d$truth$coordinates,truth=d$truth$Z,MSE=mse,CRPS_Z_joint=crps),
      draw_file=draw_path)
  },error=function(e)e)
  if(inherits(score,'error')) record$error <- paste('Scoring:',conditionMessage(score)) else {
    record$score_success <- TRUE;record$metrics <- score$metrics;record$batch_crps <- score$batch_crps
    record$pointwise_file <- file.path(root,'predictions',paste0(job$job_id,'_scores.rds'))
    atomic_save(score$pointwise,record$pointwise_file)
    record$draw_file <- score$draw_file
  }
  record$elapsed_seconds <- unname((proc.time()-started)['elapsed'])
  atomic_save(record,checkpoint)
  message('Completed ',job$job_id,': mode=',record$diagnostics$mode_status,
    ', negative Hessian=',record$diagnostics$max_negative_eigenvalues,
    ', scored=',record$score_success,', elapsed=',round(record$elapsed_seconds,1),'s')
  record
}

write_reports <- function(records) {
  status <- lapply(records,function(r){
    diag <- r$diagnostics
    data.frame(r$job,fit_success=r$fit_success,score_success=r$score_success,
      fit_seconds=r$fit_seconds,elapsed_seconds=r$elapsed_seconds,
      mode_status=if(is.null(diag))NA else diag$mode_status,
      negative_eigenvalues=if(is.null(diag))NA else diag$max_negative_eigenvalues,
      log_posterior_mode=if(is.null(diag))NA else diag$log_posterior_mode,
      integration_excess=if(is.null(diag))NA else diag$max_integration_log_relative_density,
      max_stdev_correction=if(is.null(diag))NA else diag$max_stdev_correction,
      predictor_strategy=if(is.null(diag))NA else diag$predictor_strategy,
      warnings=paste(r$warnings,collapse=' | '),error=r$error)
  })
  write.csv(do.call(rbind,status),file.path(root,'reports/fit_diagnostics.csv'),row.names=FALSE)
  params <- lapply(records,function(r)if(is.data.frame(r$parameters))cbind(r$job,r$parameters))
  write.csv(do.call(rbind,params),file.path(root,'reports/parameter_estimates.csv'),row.names=FALSE)
  boundary <- lapply(records,function(r)if(is.data.frame(r$boundary))cbind(r$job,r$boundary))
  write.csv(do.call(rbind,boundary),file.path(root,'reports/boundary_probabilities.csv'),row.names=FALSE)
  metrics <- lapply(records,function(r)if(isTRUE(r$score_success))cbind(r$job,r$metrics))
  write.csv(do.call(rbind,metrics),file.path(root,'reports/predictive_scores.csv'),row.names=FALSE)
}

records <- list()
for(i in seq_len(nrow(plan))) {
  records[[i]] <- run_one(plan[i,,drop=FALSE]);write_reports(records);gc()
}
expected <- readRDS(integrity_file)
unchanged <- identical(unname(expected),unname(tools::md5sum(names(expected))))
writeLines(paste('All original fits, simulation data, results and snapshot source unchanged:',unchanged),
           file.path(root,'reports/original_integrity.txt'))
stopifnot(unchanged)
cat('Completed all ',length(records),' planned attempts\n')
