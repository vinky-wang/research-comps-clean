# Six controlled checks using the frozen current main engine and production mesh
root <- '/tmp/psmelding-production-numerics-20261022'
setwd(root)
source('R/run_extension.R')
engine <- load_simulation()
for (directory in c('fits','checkpoints','logs','predictions','posterior_draws','simulation_data','reports'))
  dir.create(file.path(root,directory),recursive=TRUE,showWarnings=FALSE)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
configs <- list(baseline=extension_config('strength',root),
                clustering_delta2=extension_config('clustering',root))
engine$configure_runtime(configs[[1]])
mesh <- engine$build_mesh(configs[[1]])
stopifnot(mesh$n==1676L)
plan <- expand.grid(dataset=c('baseline','clustering_delta2'),start_area=c(.001,.1,.4),stringsAsFactors=FALSE)
plan$model <- ifelse(plan$dataset=='baseline','PSmelding','Melding')
plan$prior_id <- 'B'
plan$source_id <- paste(ifelse(plan$dataset=='baseline','strength_case_0003','clustering_case_0003'),plan$model,sep='_')
plan$job_id <- paste(plan$dataset,plan$model,'B',paste0('area_',gsub('[.]','p',plan$start_area)),sep='__')
write.csv(plan,'plan.csv',row.names=FALSE)
settings <- list(priors=list(B=configs[[1]]$fitting$noise_prior),plan=plan,
  phase='production_mesh',mesh_vertices=mesh$n,configs=configs,
  point_start=.1,threads='1:1',joint_draws=1000L,crps_batches=10L,
  numerical_settings='Unchanged current main control.inla, with controlled starts only',
  screens=list(log_mode_spread=.01,parameter_shift_pooled_sd=.1,prediction_RMS=.01,
    boundary_probability_spread=.01,positive_integration_log_density=.001),
  scope='Two saved datasets from one replicate, six fits; diagnostic screening only')
saveRDS(settings,'config.rds');saveRDS(mesh,'mesh.rds')
writeLines(capture.output(sessionInfo()),'session_info.txt')
saveRDS(INLA::inla.set.control.inla.default(),'runtime_defaults.rds')

original_root <- '/Users/vinky/research-comps-clean/extension_test'
reference_files <- list()
for (dataset in unique(plan$dataset)) {
  id <- plan$source_id[plan$dataset==dataset][1]
  record <- readRDS(file.path(original_root,'checkpoints',paste0('pilot_',id,'.rds')))
  files <- list(fit=file.path(record$directory,record$status$fit_file),
    data=file.path(record$directory,record$status$simulation_file),
    scores=file.path(record$directory,record$status$diagnostics_directory,'scores.rds'))
  reference_files[[dataset]] <- files
  b <- readRDS(files$fit);d <- readRDS(files$data)
  stopifnot(identical(mesh$loc,b$spde$mesh$loc),identical(mesh$graph$tv,b$spde$mesh$graph$tv))
  saveRDS(d,file.path('simulation_data',paste0(dataset,'.rds')))
}
original_files <- unlist(reference_files,use.names=FALSE)
original_md5 <- tools::md5sum(original_files)
saveRDS(list(files=reference_files,md5=original_md5),'input_provenance.rds')
saveRDS(tools::md5sum(list.files('R',full.names=TRUE)),'source_snapshot.rds')

# Inject only controlled initial values into the unchanged fitting function
# The likelihood matrices and main numerical controls must match the saved pilot
fit_body <- as.list(body(engine$fit_spatial_model))
fit_position <- which(vapply(fit_body,function(x)
  is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('fit')),logical(1)))
stopifnot(length(fit_position)==1L)
injection <- quote({
  stopifnot(identical(inla_arguments$control.predictor$A,controlled_reference$control.predictor$A),
    identical(inla_arguments$data,controlled_reference$data),
    identical(inla_arguments$control.inla,controlled_reference$control.inla),
    identical(inla_arguments$num.threads,'1:1'))
  inla_arguments$control.mode <- controlled_mode
  for (component in c('point','area')) {
    map <- likelihood_map[likelihood_map$component==component,]
    stopifnot(nrow(map)==1L,map$family=='gaussian',!map$fixed)
    start <- if(component=='point') .1 else controlled_start_area
    inla_arguments$control.family[[map$column]]$hyper$prec$initial <- -log(start)
  }
})
body(engine$fit_spatial_model) <- as.call(append(fit_body,list(injection),after=fit_position-1L))
finite_max <- function(x) if(any(is.finite(x)))max(x[is.finite(x)]) else NA_real_
fit_diagnostics <- function(fit) {
  log <- fit$logfile
  neg <- grep('number of negative eigenvalues=',log,value=TRUE,fixed=TRUE)
  neg <- as.numeric(sub('.*number of negative eigenvalues=\\[([0-9]+)\\].*','\\1',neg))
  rel <- grep('config .*log\\(rel.dens\\)=',log,value=TRUE)
  rel <- suppressWarnings(as.numeric(sub('.*log\\(rel.dens\\)= *([-0-9.]+).*','\\1',rel)))
  list(mode_status=fit$mode$mode.status,max_negative_eigenvalues=finite_max(neg),
    log_posterior_mode=fit$mode$log.posterior.mode,
    max_integration_log_relative_density=finite_max(rel),
    integration_point_exceeds_mode=any(rel>.001,na.rm=TRUE),
    max_stdev_correction=finite_max(c(fit$misc$stdev.corr.negative,fit$misc$stdev.corr.positive)),
    nfunc=fit$misc$nfunc)
}
atomic_save <- function(object,path) {
  temporary <- paste0(path,'.tmp');saveRDS(object,temporary)
  if(!file.rename(temporary,path))stop('Could not save ',path)
}
run_one <- function(job) {
  checkpoint <- file.path(root,'checkpoints',paste0(job$job_id,'.rds'))
  if(file.exists(checkpoint))return(readRDS(checkpoint))
  files <- reference_files[[job$dataset]]
  source <- readRDS(files$fit)
  data <- readRDS(file.path(root,'simulation_data',paste0(job$dataset,'.rds')))
  cfg <- configs[[job$dataset]]
  initial <- source$fit$mode$theta
  for(component in c('point','area')) {
    map <- source$likelihood_map[source$likelihood_map$component==component,]
    internal <- sub('^Precision','Log precision',map$precision_name)
    initial[internal] <- -log(if(component=='point') .1 else job$start_area)
  }
  engine$controlled_mode <- list(theta=initial,x=source$fit$mode$x,restart=TRUE,fixed=FALSE)
  engine$controlled_start_area <- job$start_area
  engine$controlled_reference <- source$inla_arguments
  fit_seed <- engine$sampling_seed(cfg,data$scenario,'fit',job$model)
  score_seed <- engine$sampling_seed(cfg,data$scenario,'scores',job$model)
  stopifnot(identical(fit_seed,source$config$seed))
  record <- list(job=job,source=files,fit_seed=fit_seed,score_seed=score_seed,
    control.mode=engine$controlled_mode,control.inla=source$inla_arguments$control.inla,
    control.family=list(prior=cfg$fitting$noise_prior),fit_success=FALSE,score_success=FALSE,
    error='',warnings=character())
  fit_path <- file.path(root,'fits',paste0(job$job_id,'.rds'))
  record$fit_file <- fit_path
  started <- proc.time();message('Starting ',job$job_id)
  if(file.exists(fit_path)) {
    b <- readRDS(fit_path)
  } else {
    b <- tryCatch(engine$fit_model(job$model,data,mesh,cfg,
      on_inla_return=function(fit,context)atomic_save(list(fit=fit,context=context),
        file.path(root,'fits',paste0(job$job_id,'_raw.rds')))),error=function(e)e)
    if(inherits(b,'error')) {
      record$error<-conditionMessage(b);record$elapsed_seconds<-unname((proc.time()-started)['elapsed'])
      record$fit_seconds<-record$elapsed_seconds;atomic_save(record,checkpoint)
      message('Fit failed: ',record$error);return(record)
    }
    atomic_save(b,fit_path)
  }
  record$fit_seconds <- b$inla_runtime[['elapsed']]
  record$warnings <- b$warnings
  record$fit_success <- TRUE;record$mode<-b$fit$mode
  record$diagnostics<-fit_diagnostics(b$fit)
  record$prediction_mean<-b$predictions$pred_mean;record$prediction_sd<-b$predictions$sd
  writeLines(b$fit$logfile,file.path(root,'logs',paste0(job$job_id,'.txt')))
  atomic_save(b$predictions,file.path(root,'predictions',paste0(job$job_id,'.rds')))
  record$parameters <- tryCatch(engine$score_parameters(b,data$truth,list(nu=1)),error=function(e)conditionMessage(e))
  record$boundary <- tryCatch(do.call(rbind,lapply(c('point','area'),function(component){
    nm<-b$likelihood_map$precision_name[b$likelihood_map$component==component]
    marginal<-b$fit$internal.marginals.hyperpar[[sub('^Precision','Log precision',nm)]]
    thresholds<-c(.001,.01,.025,.1)
    data.frame(component=component,threshold=thresholds,
      probability_below=pmin(1,pmax(0,1-INLA::inla.pmarginal(-log(thresholds),marginal))))
  })),error=function(e)conditionMessage(e))
  score <- tryCatch({
    alignment<-engine$check_prediction_alignment(b,data$truth)
    stopifnot(alignment$coordinates_valid,alignment$truth_valid,alignment$mean_valid)
    draw_path<-file.path(root,'posterior_draws',paste0(job$job_id,'.rds'))
    draws<-if(file.exists(draw_path))readRDS(draw_path) else {
      z<-engine$draw_latent_predictions(b,1000L,score_seed);atomic_save(z,draw_path);z
    }
    crps<-engine$sample_crps(draws$Z,data$truth$Z)
    mse<-(b$predictions$pred_mean-data$truth$Z)^2
    permutation<-engine$with_simulation_seed(score_seed,sample.int(ncol(draws$Z)))
    batches<-split(permutation,rep(seq_len(10),each=100))
    batch_crps<-vapply(batches,function(index)
      mean(engine$sample_crps(draws$Z[,index,drop=FALSE],data$truth$Z)),numeric(1))
    stopifnot(all(is.finite(mse)),all(is.finite(crps)))
    list(metrics=data.frame(MSE=mean(mse),CRPS_Z_joint=mean(crps),
      CRPS_MCSE_approx=sd(batch_crps)/sqrt(length(batch_crps)),
      CRPS_Z_gaussian=mean(engine$gaussian_crps(b$predictions$pred_mean,b$predictions$sd,data$truth$Z))),
      batch_crps=batch_crps,pointwise=data.frame(data$truth$coordinates,truth=data$truth$Z,MSE=mse,CRPS_Z_joint=crps),
      draw_file=draw_path)
  },error=function(e)e)
  if(inherits(score,'error'))record$error<-paste('Scoring:',conditionMessage(score)) else {
    record$score_success<-TRUE;record$metrics<-score$metrics;record$batch_crps<-score$batch_crps
    record$pointwise_file<-file.path(root,'predictions',paste0(job$job_id,'_scores.rds'))
    atomic_save(score$pointwise,record$pointwise_file);record$draw_file<-score$draw_file
  }
  record$elapsed_seconds<-unname((proc.time()-started)['elapsed'])
  atomic_save(record,checkpoint)
  message('Completed ',job$job_id,': mode=',record$diagnostics$mode_status,
    ', negative Hessian=',record$diagnostics$max_negative_eigenvalues,
    ', scored=',record$score_success,', elapsed=',round(record$elapsed_seconds,1),'s')
  record
}
if('--prepare-only' %in% commandArgs(trailingOnly=TRUE)) {
  cat('Prepared six-fit production-mesh plan; no fits run\n')
} else {
  records<-list()
  for(i in seq_len(nrow(plan))) {records[[i]]<-run_one(plan[i,,drop=FALSE]);gc()}
  saveRDS(records,'records.rds')
  unchanged<-identical(unname(original_md5),unname(tools::md5sum(names(original_md5))))
  writeLines(paste('Original fit, data and score files unchanged:',unchanged),'reports/original_integrity.txt')
  stopifnot(unchanged)
  cat('Completed ',length(records),' controlled attempts\n')
}
