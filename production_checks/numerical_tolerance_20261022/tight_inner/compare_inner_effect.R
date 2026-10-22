# Read retained records and compact optimizer summaries; no fitting.
root <- '/tmp/psmelding-tight-inner-20261022'
prior <- '/tmp/psmelding-fixed-backend-20261022'
a <- readRDS(file.path(root,'records.rds'))
b <- readRDS(file.path(prior,'records.rds'))
pa <- readRDS(file.path(root,'reports/optimizer_details.rds'))
pb <- readRDS(file.path(prior,'reports/optimizer_details.rds'))
fields <- c('posterior_mean','posterior_sd','median','lower','upper')
rows <- lapply(seq_along(a),function(i) {
  x <- a[[i]];y <- b[[i]];id <- x$job$job_id
  stopifnot(identical(id,y$job$job_id))
  data.frame(x$job,
    parameter_summaries_identical=identical(x$parameters[,fields],y$parameters[,fields]),
    parameter_max_absolute_difference=max(abs(as.matrix(x$parameters[,fields])-as.matrix(y$parameters[,fields]))),
    mode_theta_identical=identical(x$mode$theta,y$mode$theta),
    mode_theta_max_absolute_difference=max(abs(x$mode$theta-y$mode$theta)),
    mode_x_identical=identical(x$mode$x,y$mode$x),
    mode_x_max_absolute_difference=max(abs(x$mode$x-y$mode$x)),
    covariance_identical=identical(pa[[id]]$covariance,pb[[id]]$covariance),
    covariance_max_absolute_difference=max(abs(pa[[id]]$covariance-pb[[id]]$covariance)),
    prediction_mean_identical=identical(x$prediction_mean,y$prediction_mean),
    prediction_mean_max_absolute_difference=max(abs(x$prediction_mean-y$prediction_mean)),
    prediction_sd_identical=identical(x$prediction_sd,y$prediction_sd),
    log_mode_difference=x$diagnostics$log_posterior_mode-y$diagnostics$log_posterior_mode)
})
write.csv(do.call(rbind,rows),file.path(root,'reports/inner_tolerance_effect.csv'),row.names=FALSE)
print(do.call(rbind,rows),row.names=FALSE)
