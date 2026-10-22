# Read-only optimizer comparison; outputs are separate from fitted objects
args <- commandArgs(trailingOnly=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
root <- if(length(args)) args[1] else normalizePath(file.path(dirname(script_file), ".."))
out <- if(length(args)>1) args[2] else file.path(root,"reports")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
settings <- new.env()
sys.source(file.path(dirname(root),"R/config.R"),envir=settings)
.libPaths(c(settings$default_simulation_config$package_library,.libPaths()))
records <- lapply(list.files(file.path(root,"optimizer/checkpoints"),
                  pattern="\\.rds$",full.names=TRUE),readRDS)
records <- Filter(function(r)isTRUE(r$success),records)
wanted <- c("mu","gamma","rho","sigma2","theta","tau2_s","tau2_B")
finite_max <- function(x) if(any(is.finite(x)))max(x[is.finite(x)]) else NA_real_
pairs <- directions <- originals <- list()
groups <- unique(vapply(records,function(r)paste(r$job$phase,r$job$job_id,sep="__"),""))
for(group in groups) {
  rr <- Filter(function(r) paste(r$job$phase,r$job$job_id,sep="__")==group,records)
  b <- readRDS(rr[[1]]$source$fit)
  original <- list(job=rr[[1]]$job,mode=b$fit$mode,
      parameters=rr[[1]]$original_parameters,prediction_mean=b$predictions$pred_mean)
  original$job$setting <- "original"
  original$fit_file <- rr[[1]]$source$fit
  originals[[group]] <- original
  compare <- c(list(original),rr)
  if(length(compare)<2) next
  for(ij in combn(seq_along(compare),2,simplify=FALSE)) {
    a<-compare[[ij[1]]];z<-compare[[ij[2]]]
    p<-a$parameters[a$parameters$parameter %in% wanted,]
    q<-z$parameters[match(p$parameter,z$parameters$parameter),]
    osd<-original$parameters$posterior_sd[match(p$parameter,original$parameters$parameter)]
    shifts<-(q$posterior_mean-p$posterior_mean)/osd
    pairs[[length(pairs)+1L]]<-data.frame(
      phase=a$job$phase,job_id=a$job$job_id,reference=a$job$setting,comparison=z$job$setting,
      reference_mode_status=a$mode$mode.status,comparison_mode_status=z$mode$mode.status,
      mode_theta_max_abs=finite_max(abs(z$mode$theta-a$mode$theta)),
      log_posterior_difference=z$mode$log.posterior.mode-a$mode$log.posterior.mode,
      prediction_mean_RMS=sqrt(mean((z$prediction_mean-a$prediction_mean)^2)),
      parameter_mean_max_abs_shift_in_original_sd=finite_max(abs(shifts)),
      max_shift_parameter=p$parameter[which.max(abs(shifts))],
      area_mean_a=p$posterior_mean[p$parameter=="tau2_B"],
      area_mean_b=q$posterior_mean[p$parameter=="tau2_B"],
      area_median_a=p$median[p$parameter=="tau2_B"],
      area_median_b=q$median[p$parameter=="tau2_B"])
  }
  rm(b);gc()
}
chosen_strategy <- function(log,pattern) {
  row <- grep(pattern,log,value=TRUE,fixed=TRUE)
  if(!length(row)) return(NA_character_)
  sub(".*choose\\[([^]]+)\\].*","\\1",tail(row,1))
}
for(r in c(originals,records)) {
  b<-readRDS(r$fit_file);f<-b$fit
  rel<-grep("config .*log\\(rel.dens\\)=",f$logfile,value=TRUE)
  rel<-suppressWarnings(as.numeric(sub(".*log\\(rel.dens\\)= *([-0-9.]+).*","\\1",rel)))
  values<-f$misc$cov.intern.eigenvalues
  n<-length(values)
  neg<-grep("number of negative eigenvalues=",f$logfile,value=TRUE,fixed=TRUE)
  neg<-suppressWarnings(as.numeric(sub(".*number of negative eigenvalues=\\[([0-9]+)\\].*","\\1",neg)))
  directions[[length(directions)+1L]]<-data.frame(
    phase=r$job$phase,job_id=r$job$job_id,setting=r$job$setting,
    mode_status=f$mode$mode.status,negative_eigenvalues=finite_max(neg),
    nfunc=f$misc$nfunc,
    predictor_strategy=chosen_strategy(f$logfile,"Optimizing pred-strategy"),
    Qx_strategy=chosen_strategy(f$logfile,"Optimizing Qx-strategy"),
    log_posterior_mode=f$mode$log.posterior.mode,
    max_integration_log_relative_density=finite_max(rel),
    integration_point_exceeds_mode=any(rel>1e-3,na.rm=TRUE),
    max_stdev_correction=finite_max(c(f$misc$stdev.corr.negative,f$misc$stdev.corr.positive)),
    covariance_min_eigenvalue=min(values),
    covariance_max_eigenvalue=max(values),
    eigenvalue_duplicate_min_relative_gap=if(n>1)min(diff(sort(values))/max(values))else NA)
  rm(b,f);gc()
}
if(length(pairs))write.csv(do.call(rbind,pairs),file.path(out,"optimizer_pairwise.csv"),row.names=FALSE)
if(length(directions))write.csv(do.call(rbind,directions),file.path(out,"optimizer_direction_diagnostics.csv"),row.names=FALSE)
cat("Compared ",length(records)," successful controlled fits\n")

