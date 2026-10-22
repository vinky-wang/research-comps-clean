# Read-only boundary/marginal audit after all planned prior runs have completed
# No refitting, posterior draws, or modifications to study outputs
setwd('/Users/vinky/research-comps-clean/extension_test')
source('test_helpers.R')
engine$configure_runtime(configs[[1]])
root <- file.path(test_root,'prior_sensitivity')
out <- file.path(root,'reports','boundary_audit')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
plan <- read.csv(file.path(root,'plan.csv'),stringsAsFactors=FALSE)
paths <- list.files(file.path(root,'checkpoints'),pattern='\\.rds$',full.names=TRUE)
records <- lapply(paths,readRDS)
ids <- vapply(records,function(r)r$job$job_id,character(1))
stopifnot(nrow(plan)==36L,length(records)==36L,!anyDuplicated(ids),setequal(ids,plan$job_id))
thresholds <- c(.001,.01,.025,.1)
value <- function(x,default=NA_real_) if(is.null(x)||length(x)!=1L)default else unname(x)
flag <- function(x) isTRUE(x)

inspect_component <- function(b,r,component) {
  name <- b$likelihood_map$precision_name[b$likelihood_map$component==component]
  eta_name <- sub('^Precision','Log precision',name)
  parameter <- if(component=='point')'tau2_s' else 'tau2_B'
  p <- if(is.data.frame(r$parameters))r$parameters[r$parameters$parameter==parameter,,drop=FALSE] else NULL
  nd <- attr(r$parameters,'noise_diagnostics')[[component]]
  diag <- nd$diagnostics
  marginal <- b$fit$internal.marginals.hyperpar[[eta_name]]
  check <- tryCatch(engine$validate_noise_marginal(marginal,eta_name,allow_endpoint_zeros=TRUE),error=function(e)e)
  valid <- !inherits(check,'error')
  base <- data.frame(job_id=r$job$job_id,dataset=r$job$dataset,model=r$job$model,
    prior_id=r$job$prior_id,start_area=r$job$start_area,component=component,
    precision_name=name,internal_marginal_valid=valid,
    marginal_error=if(valid)'' else conditionMessage(check),
    parameter_summary_valid=if(is.null(p))NA else value(p$summary_valid,NA),
    parameter_coherence_warning=if(is.null(p))NA else value(p$coherence_warning,NA),
    noise_summary_valid=if(is.null(nd))NA else flag(nd$valid),
    noise_summary_reason=if(is.null(nd))'No saved noise_diagnostics attribute' else value(nd$reason,''),
    mean=if(is.null(p))NA else value(p$posterior_mean),
    median=if(is.null(p))NA else value(p$median),
    lower=if(is.null(p))NA else value(p$lower),
    upper=if(is.null(p))NA else value(p$upper),
    grid_relative_error=value(diag$grid_relative_error),
    adaptive_relative_error=value(diag$adaptive_relative_error),
    interpolation_relative_error=value(diag$interpolation_relative_error),
    jacobian_max_relative_error=value(diag$jacobian_max_relative_error),
    marginal_joint_mean_ratio=value(diag$marginal_joint_mean_ratio),
    marginal_joint_standardized_difference=value(diag$marginal_joint_standardized_difference),
    joint_configuration_valid=if(is.null(diag$configuration_valid))NA else flag(diag$configuration_valid),
    effective_joint_configurations=value(diag$configuration_diagnostics$effective_configurations),
    saved_log_normalizing_integral=value(diag$log_normalizing_integral))
  base$eta_min <- base$eta_max <- base$variance_support_min <- base$variance_support_max <- NA_real_
  base$left_endpoint_relative_density <- base$right_endpoint_relative_density <- NA_real_
  base$zero_endpoint_knots <- NA_integer_
  smooth <- NULL
  if(valid) {
    marginal <- check
    bounds <- range(marginal[,1L])
    base$eta_min <- bounds[1L];base$eta_max <- bounds[2L]
    vs <- sort(exp(-bounds));base$variance_support_min <- vs[1L];base$variance_support_max <- vs[2L]
    rel <- marginal[c(1L,nrow(marginal)),2L]/max(marginal[,2L])
    base$left_endpoint_relative_density <- rel[1L];base$right_endpoint_relative_density <- rel[2L]
    base$zero_endpoint_knots <- sum(marginal[,2L]==0)
    smooth <- tryCatch({
      # Use the same normalized interpolated density as the noise-parameter summaries
      # Endpoint-zero intervals are excluded here only for a diagnostic comparison
      # Their saved moment-bound validity is separately reported by noise_summary_valid
      positive <- marginal[marginal[,2L]>0,,drop=FALSE]
      bounds <- range(positive[,1L])
      density <- engine$interpolate_log_precision_density(positive)
      integrals <- lapply(c(8193L,16385L),function(n) {
        axis <- sort(unique(c(positive[,1L],seq(bounds[1L],bounds[2L],length.out=n),
          -log(thresholds)[-log(thresholds)>bounds[1L]&-log(thresholds)<bounds[2L]])))
        engine$integrate_noise_moments(axis,density,powers=0)
      })
      probability <- function(z) 1-stats::approx(z$axis,z$cdf,xout=-log(thresholds),rule=2)$y
      list(coarse=probability(integrals[[1]]),fine=probability(integrals[[2]]),
        grid_change=max(abs(probability(integrals[[1]])-probability(integrals[[2]]))))
    },error=function(e)e)
  }
  base$boundary_cdf_grid_change <- if(is.list(smooth)&&!inherits(smooth,'error'))smooth$grid_change else NA_real_
  base$boundary_cdf_error <- if(inherits(smooth,'error'))conditionMessage(smooth) else ''
  boundary <- lapply(seq_along(thresholds),function(i) {
    t <- thresholds[i];eta <- -log(t)
    stored <- if(is.data.frame(r$boundary)) r$boundary$probability_below[r$boundary$component==component&r$boundary$threshold==t] else NULL
    raw <- if(valid)tryCatch(1-INLA::inla.pmarginal(eta,marginal),error=function(e)NA_real_) else NA_real_
    grid <- if(is.list(smooth)&&!inherits(smooth,'error'))smooth$fine[i] else NA_real_
    location <- if(!valid)'invalid_marginal' else if(eta<base$eta_min)'below_saved_log_precision_support' else
      if(eta>base$eta_max)'above_saved_log_precision_support' else 'inside_saved_log_precision_support'
    outside <- valid&&location!='inside_saved_log_precision_support'
    issue <- character()
    if(!valid)issue <- c(issue,'invalid_saved_marginal')
    if(!flag(base$noise_summary_valid))issue <- c(issue,'noise_summary_not_valid')
    if(!flag(base$parameter_summary_valid))issue <- c(issue,'parameter_summary_not_valid')
    if(flag(base$parameter_coherence_warning))issue <- c(issue,'parameter_coherence_warning')
    if(outside)issue <- c(issue,'threshold_outside_saved_support_probability_not_exact_tail_statement')
    if(is.finite(raw)&&(raw<0||raw>1))issue <- c(issue,'unclamped_probability_outside_unit_interval')
    if(!is.finite(raw))issue <- c(issue,'nonfinite_boundary_probability')
    if(is.finite(raw)&&is.finite(grid)&&abs(raw-grid)>.001)issue <- c(issue,'cdf_interpolation_difference_over_0p001')
    if(valid&&max(base$left_endpoint_relative_density,base$right_endpoint_relative_density)>.001)
      issue <- c(issue,'endpoint_density_over_0p001_of_peak')
    data.frame(job_id=r$job$job_id,dataset=r$job$dataset,model=r$job$model,
      prior_id=r$job$prior_id,start_area=r$job$start_area,component=component,
      threshold=t,log_precision_threshold=eta,threshold_location=location,
      threshold_outside_saved_support=outside,
      stored_probability=value(stored),unclamped_inla_probability=raw,
      normalized_saved_support_probability=grid,
      inla_vs_normalized_support_difference=if(is.finite(raw)&&is.finite(grid))raw-grid else NA_real_,
      stored_vs_recomputed_difference=if(length(stored)==1L&&is.finite(raw))stored-raw else NA_real_,
      caution=length(issue)>0L,caution_reasons=paste(issue,collapse=' | '))
  })
  list(marginal=base,boundary=do.call(rbind,boundary))
}
all_marginals <- all_boundaries <- failed <- list()
for(i in seq_along(records)) {
  r <- records[[i]]
  if(!isTRUE(r$fit_success)||!file.exists(r$fit_file)) {
    failed[[length(failed)+1L]] <- data.frame(job_id=r$job$job_id,reason='No successful retained fit')
    next
  }
  b <- readRDS(r$fit_file)
  for(component in c('point','area')) {
    x <- inspect_component(b,r,component)
    all_marginals[[length(all_marginals)+1L]] <- x$marginal
    all_boundaries[[length(all_boundaries)+1L]] <- x$boundary
  }
  rm(b);gc()
}
marginals <- do.call(rbind,all_marginals);boundaries <- do.call(rbind,all_boundaries)
write.csv(marginals,file.path(out,'noise_marginal_audit.csv'),row.names=FALSE)
write.csv(boundaries,file.path(out,'boundary_probability_audit.csv'),row.names=FALSE)
if(length(failed))write.csv(do.call(rbind,failed),file.path(out,'unavailable_fits.csv'),row.names=FALSE)
area001 <- boundaries[boundaries$component=='area'&boundaries$threshold==.001,,drop=FALSE]
notes <- c(
 paste('Planned/completed checkpoint records:',nrow(plan),length(records)),
 paste('Retained fits inspected once:',nrow(marginals)/2),
 paste('Noise marginal rows:',nrow(marginals)),
 paste('Areal .001 threshold cautions:',sum(area001$caution)),
 paste('Areal .001 thresholds outside saved support:',sum(area001$threshold_outside_saved_support)),
 paste('Noise summaries invalid or unavailable:',sum(is.na(marginals$noise_summary_valid)|!marginals$noise_summary_valid)),
 'Boundary events use P(variance<t)=P(log precision>-log(t))',
 'CDF comparisons condition on the normalized interpolated SAVED marginal support',
 'No support check can bound unknown posterior mass outside the saved support',
 'A computed zero or one outside that support is not an exact posterior tail statement',
 'Small endpoint density is reassuring but is not a proof of negligible missing tail mass',
 'Caution flags are diagnostic prompts and do not automatically invalidate a fit',
 'CDF difference .001 and endpoint-to-peak density .001 are post-hoc numerical audit screens',
 'Joint-configuration/marginal differences are recorded, not used as automatic failure rules',
 'All original data, fitted objects, and study outputs were left unchanged')
writeLines(notes,file.path(out,'boundary_audit_summary.txt'));cat(paste(notes,collapse='\n'),'\n')
