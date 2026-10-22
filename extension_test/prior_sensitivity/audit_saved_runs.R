# Read-only audit of saved fits; outputs are separate audit reports
# Run after the requested prior/model combinations have completed
args <- commandArgs(trailingOnly=TRUE)
setwd(if(length(args)) args[1] else getwd())
source('test_helpers.R')
engine$configure_runtime(configs[[1]])
root <- file.path(test_root,'prior_sensitivity')
out <- if(length(args)>1) args[2] else file.path(root,'reports','audit')
dir.create(out,showWarnings=FALSE,recursive=TRUE)
settings <- readRDS(file.path(root,'config.rds'))
plan <- read.csv(file.path(root,'plan.csv'),stringsAsFactors=FALSE)
files <- list.files(file.path(root,'checkpoints'),pattern='\\.rds$',full.names=TRUE)
records <- lapply(files,readRDS)
names(records) <- vapply(records,function(r)r$job$job_id,character(1))
stopifnot(length(records)>0)

# Cache only the original starting values and fixed model components
source_ids <- unique(vapply(records,function(r)r$job$source_id,character(1)))
sources <- setNames(lapply(source_ids,function(id) {
  f <- job_files(read_job('finer_mesh',id))
  b <- readRDS(f$fit);sc <- readRDS(f$scores)
  list(theta=b$fit$mode$theta,x=b$fit$mode$x,fit_seed=b$config$seed,
    score_seed=sc$metadata$seed,likelihood_map=b$likelihood_map,
    family=b$inla_arguments$family,control.inla=b$inla_arguments$control.inla,
    mesh_loc=b$spde$mesh$loc,mesh_tv=b$spde$mesh$graph$tv,
    spde_prior_mean=b$spde$param.inla$theta.mu,
    spde_prior_precision=b$spde$param.inla$theta.Q,
    spde_B0=b$spde$param.inla$B0,spde_B1=b$spde$param.inla$B1,
    projection_matrices=b$projection_matrices,
    integration_weights=b$integration_weights,
    formula=deparse(b$inla_arguments$formula),
    control.fixed=b$inla_arguments$control.fixed)
}),source_ids)
same_num <- function(a,b,tol=1e-12) isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE))
controls_ok <- function(controls,job,s) {
  expected <- settings$priors[[job$prior_id]]
  if(length(controls)!=length(s$family))return(FALSE)
  checks <- vapply(c('point','area'),function(component) {
    map <- s$likelihood_map[s$likelihood_map$component==component,]
    p <- controls[[map$column]]$hyper$prec
    initial <- -log(if(component=='point')settings$point_start else job$start_area)
    identical(p$prior,expected$prior)&&same_num(p$param,expected$param)&&
      same_num(p$initial,initial)&&identical(p$fixed,FALSE)
  },logical(1))
  other <- which(s$family!='gaussian')
  all(checks)&&all(vapply(controls[other],function(x)length(x)==0L,logical(1)))
}
checkpoint_rows <- lapply(records,function(r) {
  s <- sources[[r$job$source_id]]
  initial <- s$theta
  for(component in c('point','area')) {
    map <- s$likelihood_map[s$likelihood_map$component==component,]
    nm <- sub('^Precision','Log precision',map$precision_name)
    initial[nm] <- -log(if(component=='point')settings$point_start else r$job$start_area)
  }
  controls <- s$control.inla;controls$tolerance <- settings$hyperparameter_tolerance
  data.frame(job_id=r$job$job_id,source_id=r$job$source_id,
    prior_and_noise_starts_correct=controls_ok(r$control.family,r$job,s),
    complete_theta_start_correct=identical(r$control.mode$theta,initial),
    latent_start_unchanged=identical(r$control.mode$x,s$x),
    restart_true=identical(r$control.mode$restart,TRUE),
    hyperparameters_not_fixed=identical(r$control.mode$fixed,FALSE),
    fit_seed_unchanged=identical(r$fit_seed,s$fit_seed),
    score_seed_unchanged=identical(r$score_seed,s$score_seed),
    numerical_controls_correct=identical(r$control.inla,controls))
})
checkpoint_rows <- do.call(rbind,checkpoint_rows)
write.csv(checkpoint_rows,file.path(out,'all_checkpoint_settings.csv'),row.names=FALSE)

# Read no more than one full fitted bundle per prior/model pair
successful <- Filter(function(r)isTRUE(r$fit_success)&&file.exists(r$fit_file),records)
keys <- vapply(successful,function(r)paste(r$job$prior_id,r$job$model,sep=':'),character(1))
selected <- successful[!duplicated(keys)]
section_block <- function(log,column) {
  starts <- grep('parse section=.*name=\\[INLA.Data',log)
  target <- grep(paste0('name=[INLA.Data',column,']'),log,fixed=TRUE)
  if(length(target)!=1L)return(character())
  stop <- grep('parse section=',log,fixed=TRUE)
  following <- stop[stop>target]
  log[seq.int(target,if(length(following))min(following)-1L else length(log))]
}
log_prior_ok <- function(log,column,expected,initial) {
  x <- section_block(log,column)
  p <- grep('PRIOR0->name=',x,value=TRUE,fixed=TRUE)
  p <- sub('.*PRIOR0->name=\\[([^]]+)\\].*','\\1',p)
  z <- grep('PRIOR0->PARAMETERS0=',x,value=TRUE,fixed=TRUE)
  z <- sub('.*PRIOR0->PARAMETERS0=\\[([^]]+)\\].*','\\1',z)
  z <- suppressWarnings(as.numeric(strsplit(trimws(z),'[,[:space:]]+')[[1]]))
  i <- grep('initialise log_precision[',x,value=TRUE,fixed=TRUE)
  i <- suppressWarnings(as.numeric(sub('.*initialise log_precision\\[([^]]+)\\].*','\\1',i)))
  clean <- function(a)gsub('[^a-z0-9]','',tolower(a))
  length(p)==1L&&identical(clean(p),clean(expected$prior))&&
    same_num(z,expected$param)&&length(i)==1L&&abs(i-initial)<1e-5
}
full_rows <- lapply(selected,function(r) {
  s <- sources[[r$job$source_id]];b <- readRDS(r$fit_file)
  expected <- settings$priors[[r$job$prior_id]]
  logchecks <- vapply(c('point','area'),function(component) {
    m <- s$likelihood_map[s$likelihood_map$component==component,]
    initial <- -log(if(component=='point')settings$point_start else r$job$start_area)
    log_prior_ok(b$fit$logfile,m$column,expected,initial)
  },logical(1))
  row <- data.frame(job_id=r$job$job_id,prior_id=r$job$prior_id,model=r$job$model,
    mesh_vertices=b$spde$mesh$n,
    actual_prior_arguments_correct=controls_ok(b$inla_arguments$control.family,r$job,s),
    actual_prior_logs_correct=all(logchecks),
    family_order_unchanged=identical(b$inla_arguments$family,s$family),
    mesh_unchanged=identical(b$spde$mesh$loc,s$mesh_loc)&&identical(b$spde$mesh$graph$tv,s$mesh_tv),
    spatial_prior_mean_unchanged=identical(b$spde$param.inla$theta.mu,s$spde_prior_mean),
    spatial_prior_precision_unchanged=identical(b$spde$param.inla$theta.Q,s$spde_prior_precision),
    spatial_parameterization_unchanged=identical(b$spde$param.inla$B0,s$spde_B0)&&identical(b$spde$param.inla$B1,s$spde_B1),
    projectors_unchanged=identical(b$projection_matrices,s$projection_matrices),
    integration_weights_unchanged=identical(b$integration_weights,s$integration_weights),
    model_formula_unchanged=identical(deparse(b$inla_arguments$formula),s$formula),
    fixed_effect_controls_unchanged=identical(b$inla_arguments$control.fixed,s$control.fixed),
    single_thread_requested=identical(b$inla_arguments$num.threads,'1:1'),
    gradient_h_005_logged=any(grepl('Gradient is computed using Central difference with step-length 0.005000',b$fit$logfile,fixed=TRUE)),
    epsx_0001_logged=any(grepl('epsx = 0.0001',b$fit$logfile,fixed=TRUE)),
    epsg_0001_logged=any(grepl('epsg = 0.0001',b$fit$logfile,fixed=TRUE)))
  rm(b);gc();row
})
full_rows <- do.call(rbind,full_rows)
write.csv(full_rows,file.path(out,'sampled_full_fit_settings.csv'),row.names=FALSE)
logical_ok <- function(d)all(unlist(d[vapply(d,is.logical,logical(1))]),na.rm=FALSE)
summary <- c(paste('Saved checkpoint records audited:',nrow(checkpoint_rows)),
 paste('Checkpoint setting checks all pass:',logical_ok(checkpoint_rows)),
 paste('Full fit bundles audited:',nrow(full_rows)),
 paste('Full-fit actual argument/log/fixed-model checks all pass:',logical_ok(full_rows)),
 'At most one full fitted bundle read per prior/model pair',
 'This audit does not assess statistical quality of posterior estimates',
 'No fitted object, source, or study output was modified')
writeLines(summary,file.path(out,'audit_summary.txt'));cat(paste(summary,collapse='\n'),'\n')
if(!logical_ok(checkpoint_rows)||!logical_ok(full_rows))stop('One or more audit checks failed; inspect individual CSV flags')
