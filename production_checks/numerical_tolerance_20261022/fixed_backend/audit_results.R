# Inspect retained fits only; this script does not call INLA fitting
root <- '/tmp/psmelding-fixed-backend-20261022'
setwd(root)
records <- readRDS('records.rds')
parameters <- do.call(rbind,lapply(records,function(r)data.frame(r$job,r$parameters,row.names=NULL)))
write.csv(parameters,'reports/parameter_estimates.csv',row.names=FALSE)
boundary <- do.call(rbind,lapply(records,function(r)data.frame(r$job,r$boundary,row.names=NULL)))
write.csv(boundary,'reports/boundary_probabilities.csv',row.names=FALSE)
small <- list();audit <- list();support <- list()
for(i in seq_along(records)) {
  r <- records[[i]];b <- readRDS(r$fit_file);args<-b$inla_arguments
  gaussian <- which(args$family=='gaussian')
  stopifnot(length(gaussian)==2L,
    all(vapply(args$control.family[gaussian],function(x)
      identical(x$hyper$prec$prior,'pc.prec')&&identical(x$hyper$prec$param,c(.316,.5))&&
      identical(x$hyper$prec$fixed,FALSE),logical(1))),
    identical(args$control.inla,list(control.vb=list(strategy='variance'),tolerance=1e-4)),
    identical(args$control.compute,list(config=TRUE,return.marginals=TRUE,return.marginals.predictor=TRUE,internal.opt=FALSE)),
    identical(args$num.threads,'1:1'),b$spde$mesh$n==1676L,
    identical(args$control.mode,r$control.mode),args$safe==FALSE)
  strategy <- grep('Optimizing pred-strategy',b$fit$logfile,value=TRUE,fixed=TRUE)
  covariance <- b$fit$misc$cov.intern
  small[[r$job$job_id]]<-list(theta=b$fit$mode$theta,covariance=covariance,
    stdev_negative=b$fit$misc$stdev.corr.negative,stdev_positive=b$fit$misc$stdev.corr.positive,
    nconfig=b$fit$misc$configs$nconfig)
  audit[[i]]<-data.frame(r$job,mesh_vertices=b$spde$mesh$n,threads=args$num.threads,
    noise_prior=args$control.family[[gaussian[1]]]$hyper$prec$prior,
    prior_u=.316,prior_alpha=.5,
    point_initial_variance=exp(-args$control.family[[gaussian[1]]]$hyper$prec$initial),
    area_initial_variance=exp(-args$control.family[[gaussian[2]]]$hyper$prec$initial),
    configured_tolerance='explicit 0.0001',internal_opt=args$control.compute$internal.opt,integration='auto (CCD at these dimensions)',
    predictor_strategy=if(length(strategy))sub('.*choose\\[([^]]+)\\].*','\\1',tail(strategy,1)) else NA,
    safe=args$safe,nconfig=b$fit$misc$configs$nconfig,nfunc=b$fit$misc$nfunc)
  for(component in c('point','area')) {
    map<-b$likelihood_map[b$likelihood_map$component==component,]
    m<-b$fit$internal.marginals.hyperpar[[sub('^Precision','Log precision',map$precision_name)]]
    support[[length(support)+1L]]<-data.frame(r$job,component=component,
      minimum_log_precision=min(m[,1]),maximum_log_precision=max(m[,1]),
      threshold_log_precision=-log(.001),
      threshold_outside_tabulated_support=-log(.001)>max(m[,1]))
  }
}
write.csv(do.call(rbind,audit),'reports/fit_settings_audit.csv',row.names=FALSE)
write.csv(do.call(rbind,support),'reports/boundary_support_audit.csv',row.names=FALSE)
saveRDS(small,'reports/optimizer_details.rds')
pairs<-list()
for(dataset in unique(vapply(records,function(r)r$job$dataset,character(1)))) {
  ids<-vapply(Filter(function(r)r$job$dataset==dataset,records),function(r)r$job$job_id,character(1))
  combinations<-combn(ids,2L)
  for(k in seq_len(ncol(combinations))) {
    x<-small[[combinations[1,k]]];y<-small[[combinations[2,k]]]
    pairs[[length(pairs)+1L]]<-data.frame(dataset=dataset,job1=combinations[1,k],job2=combinations[2,k],
      max_internal_mode_difference=max(abs(x$theta-y$theta)),
      relative_covariance_Frobenius_difference=sqrt(sum((x$covariance-y$covariance)^2))/
        sqrt(mean(c(sum(x$covariance^2),sum(y$covariance^2)))),
      max_negative_stdev_correction_difference=max(abs(x$stdev_negative-y$stdev_negative)),
      max_positive_stdev_correction_difference=max(abs(x$stdev_positive-y$stdev_positive)))
  }
}
write.csv(do.call(rbind,pairs),'reports/optimizer_start_comparisons.csv',row.names=FALSE)
cat('Validated all six retained fit settings; no fits run\n')
