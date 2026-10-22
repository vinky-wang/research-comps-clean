Fixed-backend diagnostic, 2026-10-22

Scope and exact settings
Exactly six controlled fits: baseline PSmelding and delta=2 Melding, each with
area-variance starts 0.001, 0.1, 0.4 and point-variance start 0.1. The original
1676-vertex mesh, saved datasets, prior B (separate Gaussian precisions,
pc.prec c(0.316,0.5)), latent/hyperparameter initial vectors, fit/scoring seeds,
1:1 threading, and 1000 joint latent draws are unchanged. All six ran
sequentially in one R 4.2 process. No full study was launched.

The only fitting-control change from the retained tighter-tolerance phase is
control.compute$internal.opt=FALSE. The full compute list is config=TRUE,
return.marginals=TRUE, return.marginals.predictor=TRUE, internal.opt=FALSE.
control.inla remains list(control.vb=list(strategy='variance'),tolerance=1e-4).
The original pilot/default phase remains the source of starting vectors;
preflight also checks equality to the immediately preceding tight phase.

Results
All six fits/scoring calls succeeded; optimizer statuses are 0, Hessian
checks have zero negative eigenvalues, and integration/summary validity pass.
Delta=2 Melding passes all stability screens. Baseline PSmelding still fails
parameter and prediction stability. The fixed-backend change does not solve
this observed baseline start sensitivity and is not adopted here.

Maximum within-group spread                         baseline       delta=2
Parameter mean shift / pooled posterior SD          0.200251       0.00995179
Interval endpoint shift / pooled posterior SD       0.645217       0.0318773
Prediction mean RMS                                 0.0108707      0.000968191
Log-posterior mode spread                            0.00196114     0.000371993

The worst baseline endpoint remains theta upper. Screens are unchanged:
0.1 pooled SD, 0.01 prediction RMS and 0.01 log-mode spread. Boundary
probability spread is zero, with tabulated marginal support recorded
separately. Predictor_strategy is NA because internal optimization is
explicitly disabled and no timing-selection log line exists; this is not a
failed fit. Results on these two datasets do not prove study-wide stability.

Evidence
reports/stability_by_prior.csv and all_attempts.csv: complete verdict/status.
reports/parameter_start_comparisons.csv and prediction_start_comparisons.csv:
all planned pairwise comparisons.
reports/fit_settings_audit.csv: actual compute controls and numerical settings.
reports/settings_hashes.csv and exact_fitting_controls.rds: exact controls,
seeds, likelihood data/matrices, mesh, initial vectors and SHA-256 digests.
reports/auto_vs_fixed_backend_stability.csv and corresponding attempts file:
comparison to the immediately preceding tighter-tolerance phase.
reports/three_phase_stability.csv: default, tighter outer tolerance, fixed
backend results together. preflight_checks.rds/txt/log validates inputs before
fitting; reports/reference_integrity.csv/txt verifies previous phases.
sha256_manifest.csv inventories local artifacts; the manifest excludes itself.

Path provenance and preservation
Original execution root:
  /tmp/psmelding-fixed-backend-20261022
Planned durable archive root:
  /Users/vinky/research-comps-clean/production_checks/numerical_tolerance_20261022/fixed_backend
Immediately preceding phase execution root:
  /tmp/psmelding-tight-tolerance-20261022
Its existing, verified durable archive:
  /Users/vinky/research-comps-clean/production_checks/numerical_tolerance_20261022/tight_tolerance
Default reference:
  /Users/vinky/research-comps-clean/production_checks/numerical_stability
Original pilot inputs:
  /Users/vinky/research-comps-clean/extension_test

Completed scripts, RDS files, recorded calls and metadata retain original
absolute /tmp paths. Copying does not rewrite them. If /tmp is removed, use
the durable root plus the original relative suffix: a record's fit_file can
be read as file.path(archive_root,'fits',basename(record$fit_file)); remap
checkpoints, predictions, posterior_draws and simulation_data similarly.
path_mapping.csv records these mappings. All raw and complete fit files are
present; an old provenance string does not indicate a missing artifact.

report_checks.R accepts the archive root and a fresh output directory as its
two arguments and reads checkpoints without fitting. Other audit scripts
contain original roots; audit a scratch copy or remap paths in memory rather
than rewriting preserved artifacts. Do not execute run_checks.R merely to
audit a copy. The parent handles durable copying and publishing; this local
manifest makes no claim of Drive upload.
