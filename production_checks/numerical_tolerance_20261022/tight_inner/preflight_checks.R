# Read-only validation of all six inputs and starts; never fits.
source('/tmp/psmelding-tight-inner-20261022/run_checks.R')
stopifnot('--prepare-only' %in% commandArgs(trailingOnly=TRUE))
preflight <- list()
for(i in seq_len(nrow(plan))) {
  job <- plan[i,,drop=FALSE]
  files <- reference_files[[job$dataset]]
  source <- readRDS(files$fit)
  prior_record <- readRDS(file.path(prior_root,'checkpoints',paste0(job$job_id,'.rds')))
  prior_fit <- readRDS(file.path(prior_root,'fits',paste0(job$job_id,'.rds')))
  d <- readRDS(file.path(root,'simulation_data',paste0(job$dataset,'.rds')))
  initial <- source$fit$mode$theta
  for(component in c('point','area')) {
    map <- source$likelihood_map[source$likelihood_map$component==component,]
    initial[sub('^Precision','Log precision',map$precision_name)] <- -log(if(component=='point') .1 else job$start_area)
  }
  mode <- list(theta=initial,x=source$fit$mode$x,restart=TRUE,fixed=FALSE)
  fit_seed <- engine$sampling_seed(configs[[job$dataset]],d$scenario,'fit',job$model)
  score_seed <- engine$sampling_seed(configs[[job$dataset]],d$scenario,'scores',job$model)
  stopifnot(identical(mode,prior_record$control.mode),
    identical(d,readRDS(file.path(prior_root,'simulation_data',paste0(job$dataset,'.rds')))),
    identical(mesh$loc,prior_fit$spde$mesh$loc),identical(mesh$graph$tv,prior_fit$spde$mesh$graph$tv),
    identical(source$inla_arguments$control.inla,prior_fit$inla_arguments$control.inla),
    identical(source$inla_arguments$control.inla,list(control.vb=list(strategy='variance'))),
    identical(fit_seed,prior_record$fit_seed),identical(score_seed,prior_record$score_seed),
    identical(configs[[job$dataset]]$fitting$noise_prior,prior_record$control.family$prior))
  tight_fit <- readRDS(file.path('/tmp/psmelding-tight-tolerance-20261022','fits',paste0(job$job_id,'.rds')))
  tight_record <- readRDS(file.path('/tmp/psmelding-tight-tolerance-20261022','checkpoints',paste0(job$job_id,'.rds')))
  stopifnot(identical(tight_fit$inla_arguments$control.inla,list(control.vb=list(strategy='variance'),tolerance=1e-4)),
    identical(tight_fit$inla_arguments$control.compute,source$inla_arguments$control.compute),
    identical(tight_record$control.mode,mode),
    identical(tight_record$fit_seed,fit_seed),identical(tight_record$score_seed,score_seed),
    identical(tight_fit$inla_arguments$data,prior_fit$inla_arguments$data),
    identical(tight_fit$inla_arguments$control.predictor,prior_fit$inla_arguments$control.predictor),
    identical(tight_fit$spde$mesh$loc,mesh$loc),identical(tight_fit$spde$mesh$graph$tv,mesh$graph$tv))
  fixed_fit <- readRDS(file.path('/tmp/psmelding-fixed-backend-20261022','fits',paste0(job$job_id,'.rds')))
  fixed_record <- readRDS(file.path('/tmp/psmelding-fixed-backend-20261022','checkpoints',paste0(job$job_id,'.rds')))
  stopifnot(identical(fixed_fit$inla_arguments$control.inla,list(control.vb=list(strategy='variance'),tolerance=1e-4)),
    identical(fixed_fit$inla_arguments$control.compute,c(source$inla_arguments$control.compute,list(internal.opt=FALSE))),
    identical(fixed_record$control.mode,mode),
    identical(fixed_record$fit_seed,fit_seed),identical(fixed_record$score_seed,score_seed),
    identical(fixed_fit$inla_arguments$data,prior_fit$inla_arguments$data),
    identical(fixed_fit$inla_arguments$control.predictor,prior_fit$inla_arguments$control.predictor),
    identical(fixed_fit$spde$mesh$loc,mesh$loc),identical(fixed_fit$spde$mesh$graph$tv,mesh$graph$tv))
  preflight[[i]] <- list(job=job,fit_seed=fit_seed,score_seed=score_seed,control.mode=mode,
    reference_control.inla=source$inla_arguments$control.inla,
    control.inla=c(source$inla_arguments$control.inla,list(tolerance=1e-4,tolerance.step=1e-5)),
    control.compute=c(source$inla_arguments$control.compute,list(internal.opt=FALSE)),
    inputs_and_initial_values_identical=TRUE)
}
saveRDS(preflight,'preflight_checks.rds')
writeLines(capture.output(lapply(preflight,function(x) {
  x$control.mode$x <- paste(length(x$control.mode$x),'original latent initial values; exact match verified')
  x
})),'preflight_checks.txt')
cat('All six inputs/starts/seeds match all retained phases; controls match fixed-backend phase except tolerance.step=1e-5.\n')
