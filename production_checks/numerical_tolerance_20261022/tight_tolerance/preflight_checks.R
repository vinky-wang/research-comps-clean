# Read-only validation of all six inputs and starts; never fits.
source('/tmp/psmelding-tight-tolerance-20261022/run_checks.R')
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
  preflight[[i]] <- list(job=job,fit_seed=fit_seed,score_seed=score_seed,control.mode=mode,
    reference_control.inla=source$inla_arguments$control.inla,
    control.inla=c(source$inla_arguments$control.inla,list(tolerance=1e-4)),
    inputs_and_initial_values_identical=TRUE)
}
saveRDS(preflight,'preflight_checks.rds')
writeLines(capture.output(lapply(preflight,function(x) {
  x$control.mode$x <- paste(length(x$control.mode$x),'original latent initial values; exact match verified')
  x
})),'preflight_checks.txt')
cat('All six data/mesh/starts/seeds/prior/default controls exactly match the prior diagnostic.\n')
