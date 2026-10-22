# Read-only comparisons with the retained fixed-backend diagnostic.
root <- '/tmp/psmelding-tight-inner-20261022'
prior <- '/tmp/psmelding-fixed-backend-20261022'
setwd(root)
records <- readRDS('records.rds')
settings_rows <- list()
controls <- list()
object_hash <- function(x) digest::digest(x,algo='sha256')
for(i in seq_along(records)) {
  r <- records[[i]]
  id <- r$job$job_id
  new <- readRDS(file.path(root,'fits',paste0(id,'.rds')))
  old <- readRDS(file.path(prior,'fits',paste0(id,'.rds')))
  old_record <- readRDS(file.path(prior,'checkpoints',paste0(id,'.rds')))
  new_args <- new$inla_arguments
  old_args <- old$inla_arguments
  new_outer_only <- new_args$control.inla
  new_outer_only$tolerance.step <- NULL
  compare_names <- c('family','data','control.predictor','control.compute','control.fixed',
    'control.family','control.mode','num.threads','safe','verbose')
  same <- vapply(compare_names,function(nm)identical(new_args[[nm]],old_args[[nm]]),logical(1))
  stopifnot(all(same),identical(new_outer_only,old_args$control.inla),
    identical(new_args$control.inla$tolerance.step,1e-5),
    identical(new_args$control.compute$internal.opt,FALSE),
    identical(new_args$control.inla$tolerance,1e-4),
    identical(new$spde$mesh$loc,old$spde$mesh$loc),
    identical(new$spde$mesh$graph$tv,old$spde$mesh$graph$tv),
    identical(r$fit_seed,old_record$fit_seed),identical(r$score_seed,old_record$score_seed))
  components <- list(control.inla=new_args$control.inla,reference_control.inla=old_args$control.inla,
    reference_control.compute=old_args$control.compute,
    control.mode=new_args$control.mode,control.family=new_args$control.family,
    control.compute=new_args$control.compute,control.predictor=new_args$control.predictor,
    family=new_args$family,data=new_args$data,
    mesh=list(loc=new$spde$mesh$loc,tv=new$spde$mesh$graph$tv),
    seeds=list(fit=r$fit_seed,score=r$score_seed))
  controls[[id]] <- components
  settings_rows[[i]] <- data.frame(r$job,all_other_checked_arguments_identical=all(same),
    only_inla_control_change='tolerance.step: default 0.001 -> explicit 0.00001',
    fit_seed=r$fit_seed,score_seed=r$score_seed,
    component=names(components),sha256=vapply(components,object_hash,character(1)),row.names=NULL)
}
saveRDS(controls,'exact_fitting_controls.rds')
write.csv(do.call(rbind,settings_rows),'reports/settings_hashes.csv',row.names=FALSE)
old_groups <- read.csv(file.path(prior,'reports/stability_by_prior.csv'))
new_groups <- read.csv('reports/stability_by_prior.csv')
summary <- merge(old_groups,new_groups,by=c('dataset','model','prior_id'),suffixes=c('_outer_only_fixed_backend','_tight_inner_fixed_backend'))
write.csv(summary,'reports/outer_vs_inner_stability.csv',row.names=FALSE)
old_jobs <- read.csv(file.path(prior,'reports/all_attempts.csv'))
new_jobs <- read.csv('reports/all_attempts.csv')
jobs <- merge(old_jobs,new_jobs,by=c('dataset','model','prior_id','start_area','source_id','job_id'),
  suffixes=c('_outer_only_fixed_backend','_tight_inner_fixed_backend'))
write.csv(jobs,'reports/outer_vs_inner_attempts.csv',row.names=FALSE)
cat('All six retained fits match their fixed-backend references for every checked argument; only tolerance.step differs.\n')
print(new_groups[,c('dataset','all_components_pass','maximum_parameter_mean_shift_pooled_sd',
  'maximum_interval_endpoint_shift_pooled_sd','max_prediction_mean_RMS','log_mode_spread','failed_components')],row.names=FALSE)

reference_roots <- c(default_0p005='/Users/vinky/research-comps-clean/production_checks/numerical_stability',
  tight_0p0001='/tmp/psmelding-tight-tolerance-20261022',fixed_backend=prior,tight_inner=root)
four_phases <- do.call(rbind,lapply(names(reference_roots),function(phase) {
  x <- read.csv(file.path(reference_roots[[phase]],'reports/stability_by_prior.csv'))
  data.frame(phase=phase,x)
}))
write.csv(four_phases,'reports/four_phase_stability.csv',row.names=FALSE)
