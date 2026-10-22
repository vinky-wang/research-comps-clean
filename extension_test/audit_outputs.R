# Review saved outputs without refitting or changing the source snapshot
source('test_helpers.R')
verify_snapshot()
paths <- list.files('checkpoints', pattern = '^(pilot|finer_mesh)_.*\\.rds$', full.names = TRUE)
records <- lapply(paths, readRDS)
rows <- list(); hessian_rows <- list(); excerpts <- character(); reference <- NULL
for (record in records) {
  files <- job_files(record)
  if (is.null(files)) next
  bundle <- readRDS(files$fit)
  data <- readRDS(files$data)
  if (is.null(reference)) reference <- list(prior_mean=bundle$spde$param.inla$theta.mu,
    prior_precision=bundle$spde$param.inla$theta.Q, controls=bundle$inla_arguments$control.inla)
  score <- if(file.exists(files$scores)) readRDS(files$scores) else NULL
  gamma <- if(is.null(score)) NULL else score$parameters[score$parameters$parameter=='gamma',]
  heterogeneous <- data$scenario$sampling_type=='piecewise' && data$scenario$gamma_L != data$scenario$gamma_R
  rows[[length(rows)+1L]] <- data.frame(phase=record$phase,job_id=record$job$job_id,
    n_point=nrow(data$points),n_area=nrow(data$areas),n_prediction=nrow(bundle$predictions),
    latent_surface_retained=length(data$truth$Z)==2500,
    predictor_marginals_retained=length(bundle$fit$marginals.linear.predictor)>0,
    joint_configurations_retained=!is.null(bundle$fit$misc$configs),
    same_prior_mean=identical(bundle$spde$param.inla$theta.mu,reference$prior_mean),
    same_prior_precision=identical(bundle$spde$param.inla$theta.Q,reference$prior_precision),
    same_numerical_controls=identical(bundle$inla_arguments$control.inla,reference$controls),
    correct_heterogeneous_gamma_truth=if(heterogeneous && record$job$model=='PSmelding')
      !is.null(gamma) && nrow(gamma)==1 && is.na(gamma$truth) && is.na(gamma$covered) else TRUE,
    mesh_vertices=bundle$spde$n.spde)
  log <- bundle$fit$logfile
  h <- grep('ensure spd for Hessian|negative eigenvalues|mode seems',log,value=TRUE)
  values <- grep('number of negative eigenvalues=',log,value=TRUE)
  values <- as.integer(sub('.*number of negative eigenvalues=\\[([0-9]+)\\].*','\\1',values))
  hessian_rows[[length(hessian_rows)+1L]] <- data.frame(phase=record$phase,job_id=record$job$job_id,
    mode_status=bundle$fit$mode$mode.status,spd_check_logged=any(grepl('ensure spd for Hessian',log)),
    negative_eigenvalue_checks=length(values),
    max_negative_eigenvalues=if(length(values)) max(values,na.rm=TRUE) else NA_integer_,
    final_negative_eigenvalues=if(length(values)) tail(values,1) else NA_integer_,
    mode_seems_fine=any(grepl('mode seems fine',log,fixed=TRUE)))
  excerpts <- c(excerpts,paste(record$phase,record$job$job_id),h,'')
  rm(bundle,data); gc()
}
if(length(rows)) write.csv(do.call(rbind,rows),'reports/output_audit.csv',row.names=FALSE)
if(length(hessian_rows)) write.csv(do.call(rbind,hessian_rows),'reports/hessian_review.csv',row.names=FALSE)
writeLines(excerpts,'reports/hessian_excerpts.txt')
ids <- c('strength_case_0003_Melding','areal_noise_case_0003_Melding',
 'spatial_gamma_case_0001_PSmelding','clustering_case_0001_Melding')
reference <- NULL; baseline <- list()
for(id in ids) {
  files <- job_files(read_job('pilot',id))
  if(is.null(files)) next
  data <- readRDS(files$data)
  if(is.null(reference)) reference <- data
  baseline[[id]] <- data.frame(job_id=id,
    same_Z=identical(data$truth$Z,reference$truth$Z),
    same_points=identical(data$points,reference$points),
    same_areas=identical(data$areas,reference$areas),
    same_location_seed=identical(data$seeds$locations,reference$seeds$locations))
}
if(length(baseline)) write.csv(do.call(rbind,baseline),'reports/shared_baseline_audit.csv',row.names=FALSE)
original <- file.path('..','R',basename(list.files('R',full.names=TRUE)))
unchanged <- identical(unname(tools::md5sum(original)),unname(tools::md5sum(list.files('R',full.names=TRUE))))
writeLines(paste('Original main R scripts unchanged:',unchanged),'reports/source_integrity.txt')
cat('Saved output, prior, baseline and Hessian audits\n')
# Post-hoc repeatability diagnostic prompted by different baseline mode flags
# These comparisons use already-computed fits; they add no selected retries
repeatability <- list()
for(model in c('Melding','PSmelding')) {
  base_id <- paste('strength_case_0003',model,sep='_')
  bf <- job_files(read_job('pilot',base_id))
  if(is.null(bf) || !file.exists(bf$scores)) next
  b <- readRDS(bf$fit); bs <- readRDS(bf$scores)
  targets <- c(paste('areal_noise_case_0003',model,sep='_'),
    paste('clustering_case_0001',model,sep='_'),
    if(model=='PSmelding') 'spatial_gamma_case_0001_PSmelding')
  for(id in targets) {
    ff <- job_files(read_job('pilot',id))
    if(is.null(ff) || !file.exists(ff$scores)) next
    f <- readRDS(ff$fit); fs <- readRDS(ff$scores)
    wanted <- c('mu','gamma','tau2_s','tau2_B','rho','sigma2','theta')
    p <- bs$parameters[bs$parameters$parameter %in% wanted,]
    q <- fs$parameters[match(p$parameter,fs$parameters$parameter),]
    shifts <- abs((q$posterior_mean-p$posterior_mean)/p$posterior_sd)
    repeatability[[id]] <- data.frame(model=model,reference_job=base_id,comparison_job=id,
      same_R_seed=identical(b$config$seed,f$config$seed),
      reference_mode=b$fit$mode$mode.status,comparison_mode=f$fit$mode$mode.status,
      prediction_mean_rms_change=sqrt(mean((b$predictions$pred_mean-f$predictions$pred_mean)^2)),
      max_parameter_mean_shift_in_reference_sd=max(shifts,na.rm=TRUE))
    rm(f);gc()
  }
  rm(b);gc()
}
if(length(repeatability)) write.csv(do.call(rbind,repeatability),
 'reports/baseline_fit_repeatability.csv',row.names=FALSE)
