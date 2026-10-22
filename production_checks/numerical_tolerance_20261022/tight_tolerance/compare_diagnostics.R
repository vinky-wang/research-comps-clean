# Read-only comparisons with the retained default-tolerance diagnostic.
root <- '/tmp/psmelding-tight-tolerance-20261022'
prior <- '/Users/vinky/research-comps-clean/production_checks/numerical_stability'
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
  new_untightened <- new_args$control.inla
  new_untightened$tolerance <- NULL
  compare_names <- c('family','data','control.predictor','control.compute','control.fixed',
    'control.family','control.mode','num.threads','safe','verbose')
  same <- vapply(compare_names,function(nm)identical(new_args[[nm]],old_args[[nm]]),logical(1))
  stopifnot(all(same),identical(new_untightened,old_args$control.inla),
    identical(new_args$control.inla$tolerance,1e-4),
    identical(new$spde$mesh$loc,old$spde$mesh$loc),
    identical(new$spde$mesh$graph$tv,old$spde$mesh$graph$tv),
    identical(r$fit_seed,old_record$fit_seed),identical(r$score_seed,old_record$score_seed))
  components <- list(control.inla=new_args$control.inla,reference_control.inla=old_args$control.inla,
    control.mode=new_args$control.mode,control.family=new_args$control.family,
    control.compute=new_args$control.compute,control.predictor=new_args$control.predictor,
    family=new_args$family,data=new_args$data,
    mesh=list(loc=new$spde$mesh$loc,tv=new$spde$mesh$graph$tv),
    seeds=list(fit=r$fit_seed,score=r$score_seed))
  controls[[id]] <- components
  settings_rows[[i]] <- data.frame(r$job,all_other_checked_arguments_identical=all(same),
    only_inla_control_change='tolerance: default 0.005 -> explicit 0.0001',
    fit_seed=r$fit_seed,score_seed=r$score_seed,
    component=names(components),sha256=vapply(components,object_hash,character(1)),row.names=NULL)
}
saveRDS(controls,'exact_fitting_controls.rds')
write.csv(do.call(rbind,settings_rows),'reports/settings_hashes.csv',row.names=FALSE)
old_groups <- read.csv(file.path(prior,'reports/stability_by_prior.csv'))
new_groups <- read.csv('reports/stability_by_prior.csv')
summary <- merge(old_groups,new_groups,by=c('dataset','model','prior_id'),suffixes=c('_default_0p005','_tight_0p0001'))
write.csv(summary,'reports/default_vs_tight_stability.csv',row.names=FALSE)
old_jobs <- read.csv(file.path(prior,'reports/all_attempts.csv'))
new_jobs <- read.csv('reports/all_attempts.csv')
jobs <- merge(old_jobs,new_jobs,by=c('dataset','model','prior_id','start_area','source_id','job_id'),
  suffixes=c('_default_0p005','_tight_0p0001'))
write.csv(jobs,'reports/default_vs_tight_attempts.csv',row.names=FALSE)
cat('All six retained fits match their default-tolerance references for every checked argument; only tolerance differs.\n')
print(new_groups[,c('dataset','all_components_pass','maximum_parameter_mean_shift_pooled_sd',
  'maximum_interval_endpoint_shift_pooled_sd','max_prediction_mean_RMS','log_mode_spread','failed_components')],row.names=FALSE)
