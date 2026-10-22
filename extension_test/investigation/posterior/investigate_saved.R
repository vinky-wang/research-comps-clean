# Read retained fits and record posterior/optimizer diagnostics without refitting
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else getwd()
setwd(normalizePath(root, mustWork = TRUE))
source('test_helpers.R')
engine$configure_runtime(configs[[1]])
out <- if (length(args) > 1L) args[2] else file.path(getwd(), 'investigation', 'posterior')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ids <- c('strength_case_0003_Melding', 'clustering_case_0003_Melding', 'strength_case_0003_PSmelding')
targets <- expand.grid(phase = c('pilot', 'finer_mesh'), job_id = ids, stringsAsFactors = FALSE)
targets <- rbind(targets, data.frame(phase = 'pilot', job_id = 'areal_noise_case_0003_Melding'))
noise_rows <- list(); param_rows <- list(); diagnostic_rows <- list(); theta_rows <- list(); covariance_rows <- list(); retained_marginals <- list(); excerpts <- character()
trap <- function(x,y) sum(diff(x)*(head(y,-1)+tail(y,-1))/2)
for (i in seq_len(nrow(targets))) {
  phase <- targets$phase[i]; id <- targets$job_id[i]; key <- paste(phase,id,sep='/')
  record <- read_job(phase,id); files <- job_files(record)
  b <- readRDS(files$fit); f <- b$fit; s <- readRDS(files$scores)
  pars <- s$parameters
  param_rows[[key]] <- cbind(data.frame(phase=phase,job_id=id),pars)
  logs <- f$logfile
  negatives <- grep('number of negative eigenvalues=', logs, value=TRUE)
  negatives <- as.integer(sub('.*number of negative eigenvalues=\\[([0-9]+)\\].*','\\1',negatives))
  cv <- f$misc$cov.intern
  eig <- eigen(cv,symmetric=TRUE,only.values=TRUE)$values
  config_log_density <- vapply(f$misc$configs$config, function(x) x$log.posterior.orig, numeric(1))
  config_mode_index <- which.min(vapply(f$misc$configs$config, function(x) max(abs(x$theta - f$mode$theta)), numeric(1)))
  diagnostic_rows[[key]] <- data.frame(phase=phase, job_id=id, mode_status=f$mode$mode.status,
    nfunc=f$misc$nfunc, max_configuration_log_density_above_mode=max(config_log_density)-config_log_density[config_mode_index], log_posterior_mode=f$mode$log.posterior.mode, mesh_vertices=b$spde$n.spde,
    negative_eigenvalue_checks=length(negatives), max_negative_eigenvalues=if(length(negatives))max(negatives) else NA,
    final_negative_eigenvalues=if(length(negatives))tail(negatives,1) else NA,
    saved_covariance_smallest_eigenvalue=min(eig), saved_covariance_largest_eigenvalue=max(eig),
    saved_covariance_condition=max(eig)/min(eig), max_abs_correlation=max(abs(f$misc$cor.intern[row(f$misc$cor.intern)!=col(f$misc$cor.intern)])),
    asymmetry_correction_min=min(c(f$misc$stdev.corr.negative,f$misc$stdev.corr.positive)),
    asymmetry_correction_max=max(c(f$misc$stdev.corr.negative,f$misc$stdev.corr.positive)),
    fit_ok=if(length(f$ok))isTRUE(f$ok) else NA, warnings=paste(b$warnings,collapse='; '),
    fit_seed=b$config$seed, num_threads=b$inla_arguments$num.threads)
  theta_rows[[key]] <- data.frame(phase=phase,job_id=id, parameter=names(f$mode$theta), mode_theta=unname(f$mode$theta),
    native_marginal_mean=f$internal.summary.hyperpar[names(f$mode$theta),'mean'])
  cvnames <- f$misc$theta.tags
  for (j in seq_len(nrow(cv))) for(k in seq_len(ncol(cv))) {
    covariance_rows[[paste(key,j,k)]] <- data.frame(phase=phase,job_id=id,parameter_1=cvnames[j],parameter_2=cvnames[k],covariance=cv[j,k],correlation=f$misc$cor.intern[j,k])
  }
  excerpt <- grep('Hessian|eigenvalue|mode seems|mode.status|Optim|optim|gradient|Newton|negative|correction|\"threads\"|Smart optimize', logs, value=TRUE, ignore.case=TRUE)
  excerpts <- c(excerpts,paste('FIT',key),excerpt,'')
  retained_marginals[[key]] <- list(native=f$marginals.hyperpar, internal=f$internal.marginals.hyperpar, variance=list(),spatial=list())
  noise_map <- engine$get_gaussian_likelihoods(b)
  for (j in seq_len(nrow(noise_map))) {
    nm <- noise_map$precision_name[j]; parameter <- if(noise_map$component[j]=='area')'tau2_B' else 'tau2_s'
    summary <- engine$summarize_noise_variance(f,nm)
    stopifnot(summary$valid)
    eta <- f$internal.marginals.hyperpar[[sub('^Precision','Log precision',nm)]]
    eta <- eta[eta[,2]>0,,drop=FALSE]
    ld <- engine$interpolate_log_precision_density(eta)
    axis <- sort(unique(c(eta[,1],seq(min(eta[,1]),max(eta[,1]),length.out=32769L))))
    integ <- engine$integrate_noise_moments(axis,ld)
    cdf <- approxfun(integ$axis,integ$cdf,rule=2)
    threshold <- c(.001,.01,.025)
    prob <- 1-cdf(-log(threshold))
    # Density of log variance is f_eta(-log variance), while variance density
    # includes a 1/variance Jacobian; these modes need not coincide
    log_y <- ld(axis); dens <- exp(log_y-integ$log_integrals[1])
    v <- exp(-axis); ord <- order(v)
    retained_marginals[[key]]$variance[[parameter]] <- data.frame(variance=v[ord],log_variance=-axis[ord],density_log_variance=dens[ord],density_variance=(dens/v)[ord])
    native <- pars[pars$parameter==parameter,]
    noise_rows[[paste(key,parameter)]] <- data.frame(phase=phase,job_id=id,parameter=parameter,
      mean=native$posterior_mean,sd=native$posterior_sd,q025=native$lower,median=native$median,q975=native$upper,
      variance_density_mode=summary$mode,variance_at_log_variance_density_mode=exp(-axis[which.max(log_y)]),
      variance_at_joint_hyperparameter_mode=exp(-f$mode$theta[sub('^Precision','Log precision',nm)]),
      prob_lt_0p001=prob[1],prob_lt_0p01=prob[2],prob_lt_0p025=prob[3],
      saved_variance_support_min=min(v),saved_variance_support_max=max(v),
      log_density_endpoint_rel_peak_low_variance=exp(tail(log_y,1)-max(log_y)),
      log_density_endpoint_rel_peak_high_variance=exp(head(log_y,1)-max(log_y)),
      saved_log_density_integral=exp(integ$log_integrals[1]),
      endpoint_trim_relative_mass_bound=unname(summary$diagnostics$endpoint_zero_handling$relative_integral_upper_bounds['mass']))
  }
  spatial <- INLA::inla.spde2.result(f,b$spde_field_name,b$spde,do.transform=FALSE)
  retained_marginals[[key]]$spatial$sigma2 <- INLA::inla.tmarginal(exp, spatial$marginals.log.variance.nominal[[1]])
  retained_marginals[[key]]$spatial$rho <- INLA::inla.tmarginal(function(x)sqrt(8)*exp(-x),f$marginals.hyperpar[[paste('Theta2 for',b$spde_field_name)]])
  retained_marginals[[key]]$spatial$gamma <- f$marginals.hyperpar[['Beta for j']]
  rm(b,f,s);gc()
}
write.csv(do.call(rbind,noise_rows),file.path(out,'noise_boundary.csv'),row.names=FALSE)
write.csv(do.call(rbind,param_rows),file.path(out,'parameter_summaries.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostic_rows),file.path(out,'optimizer_diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,theta_rows),file.path(out,'joint_mode_hyperparameters.csv'),row.names=FALSE)
write.csv(do.call(rbind,covariance_rows),file.path(out,'internal_covariance.csv'),row.names=FALSE)
writeLines(excerpts,file.path(out,'optimizer_log_excerpts.txt'))
saveRDS(retained_marginals,file.path(out,'selected_marginals.rds'))
# Present mass on log10 variance, with the true generating variance in black
png(file.path(out,'noise_posterior_distributions.png'),width=1800,height=1500,res=170)
par(mfrow=c(3,2),mar=c(4.1,4.5,3.4,1),oma=c(0,0,2,0))
labels <- c('Baseline Melding','Delta=2 Melding','Baseline PSmelding')
for(i in seq_along(ids))for(parameter in c('tau2_s','tau2_B')) {
  a<-retained_marginals[[paste('pilot',ids[i],sep='/')]]$variance[[parameter]]
  z<-retained_marginals[[paste('finer_mesh',ids[i],sep='/')]]$variance[[parameter]]
  yl<-max(c(a$density_log_variance,z$density_log_variance))*log(10)
  plot(log10(a$variance),a$density_log_variance*log(10),type='l',col='#2166ac',lwd=2,xlim=if(parameter=='tau2_s')c(-2,-.3) else c(-6,.5),ylim=c(0,yl),
    xlab=expression(log[10](variance)),ylab=expression('Density of '*log[10](variance)),main=paste(labels[i],if(parameter=='tau2_s')'point error' else 'areal error'))
  lines(log10(z$variance),z$density_log_variance*log(10),col='#b2182b',lwd=2,lty=2)
  abline(v=log10(.1),lty=3)
  legend('topright',c('Pilot mesh','Finer mesh','Generating variance 0.1'),col=c('#2166ac','#b2182b','black'),lty=c(1,2,3),lwd=c(2,2,1),bty='n',cex=.65)
}
mtext('Density is shown on log variance; peaks are not modes of variance density',outer=TRUE,cex=.85)
dev.off()
png(file.path(out,'spatial_hyperparameter_distributions.png'),width=1800,height=1500,res=170)
par(mfrow=c(3,2),mar=c(4.1,4.5,3.4,1))
for(i in seq_along(ids))for(parameter in c('rho','sigma2')) {
 a<-retained_marginals[[paste('pilot',ids[i],sep='/')]]$spatial[[parameter]]
 z<-retained_marginals[[paste('finer_mesh',ids[i],sep='/')]]$spatial[[parameter]]
 xr<-range(c(INLA::inla.qmarginal(c(.001,.999),a),INLA::inla.qmarginal(c(.001,.999),z)))
 plot(a,type='l',col='#2166ac',lwd=2,xlim=xr,ylim=c(0,max(a[,2],z[,2])),xlab=parameter,ylab='Posterior density',main=paste(labels[i],parameter))
 lines(z,col='#b2182b',lwd=2,lty=2)
 abline(v=if(parameter=='rho').2 else 1,lty=3)
 legend('topright',c('Pilot mesh','Finer mesh','Generating value'),col=c('#2166ac','#b2182b','black'),lty=c(1,2,3),lwd=c(2,2,1),bty='n',cex=.7)
}
dev.off()
writeLines(c('Posterior boundary probabilities use the retained finite-support log-precision marginal, normalized and integrated on a dense grid',
'They do not establish behavior beyond saved marginal support or validate the INLA approximation',
'No continuous precision/variance posterior puts point mass at zero variance',
'The variance mode, variance at log-variance mode, and variance at joint hyperparameter mode are reported separately',
'The saved cov.intern is post-Hessian-correction output; its positive eigenvalues do not negate original negative-eigenvalue diagnostics',
'All calculations are diagnostic reanalysis of saved fits; no fit was replaced or refitted'),file.path(out,'methods.txt'))


# Compare identical saved baseline data through the mode and curvature stages
b<-readRDS(job_files(read_job('pilot','strength_case_0003_Melding'))$fit)$fit
f<-readRDS(job_files(read_job('pilot','areal_noise_case_0003_Melding'))$fit)$fit
sink(file.path(out,'baseline_saved_comparison.txt'))
cat('Joint hyperparameter modes\n');print(cbind(strength=b$mode$theta,areal_noise=f$mode$theta,difference=f$mode$theta-b$mode$theta),digits=12)
cat('\nLog posterior mode difference:',f$mode$log.posterior.mode-b$mode$log.posterior.mode,'\n')
cat('\nMaximum absolute saved internal covariance difference:',max(abs(b$misc$cov.intern-f$misc$cov.intern)),'\n')
cat('Maximum absolute saved correlation difference:',max(abs(b$misc$cor.intern-f$misc$cor.intern)),'\n')
cat('\nCovariance eigenvalues\n');print(cbind(strength=b$misc$cov.intern.eigenvalues,areal_noise=f$misc$cov.intern.eigenvalues),digits=12)
cat('\nAbsolute overlaps of eigenvectors (rows strength; columns areal_noise)\n');print(abs(t(b$misc$cov.intern.eigenvectors)%*%f$misc$cov.intern.eigenvectors),digits=12)
cat('\nDirectional standard deviation correction factors\n');print(cbind(strength_neg=b$misc$stdev.corr.negative,strength_pos=b$misc$stdev.corr.positive,areal_neg=f$misc$stdev.corr.negative,areal_pos=f$misc$stdev.corr.positive),digits=12)
for(x in list(strength=b,areal_noise=f)) {
 log<-x$logfile; idx<-grep('Eigenvectors of the Hessian',log);cat('\nRaw Hessian output excerpt\n');cat(log[seq(idx-4,idx+25)],sep='\n')
}
sink()

# Confirm the installed likelihood prior and record its variance-scale implications
hp <- INLA::inla.models()$likelihood$gaussian$hyper$theta1
stopifnot(hp$prior == 'loggamma', identical(as.numeric(hp$param), c(1,5e-5)))
prior <- data.frame(parameter='noise_variance', precision_prior='Gamma(shape=1, rate=0.00005)',
 variance_prior='InverseGamma(shape=1, scale=0.00005)', variance_density_mode=5e-5/2,
 median=5e-5/log(2), q025=5e-5/(-log(.025)), q975=5e-5/(-log(.975)),
 prob_lt_0p001=exp(-5e-5/.001), prob_lt_0p01=exp(-5e-5/.01), prob_lt_0p025=exp(-5e-5/.025))
write.csv(prior,file.path(out,'noise_prior_reference.csv'),row.names=FALSE)
prior_evidence <- c(paste('Installed INLA:',utils::packageVersion('INLA')),
 capture.output(str(hp)), '', 'Saved baseline Melding log excerpts:',
 grep('PRIOR0->PARAMETERS0=|PRIOR0->name=|PRIOR0->from_theta=',b$logfile,value=TRUE,fixed=FALSE))
writeLines(prior_evidence,file.path(out,'noise_prior_evidence.txt'))
cat('Saved posterior investigation complete\n')
