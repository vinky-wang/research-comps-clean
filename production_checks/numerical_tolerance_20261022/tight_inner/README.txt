Final inner-tolerance diagnostic, 2026-10-22

Scope and exact settings
Exactly six controlled fits: baseline PSmelding and delta=2 Melding, each with
area-variance starts 0.001, 0.1, 0.4 and point-variance start 0.1. The original
1676-vertex mesh, saved datasets, prior B (separate Gaussian precisions,
pc.prec c(0.316,0.5)), latent/hyperparameter initial vectors, fit/scoring seeds,
1:1 threading, and 1000 joint latent draws are unchanged. All six ran
sequentially in one R 4.2 process. No full study was launched.

The only fitting-control change from the retained fixed-backend phase is
control.inla$tolerance.step=1e-5. The full numerical controls are:
control.inla = list(control.vb=list(strategy='variance'),tolerance=1e-4,
                    tolerance.step=1e-5)
control.compute = list(config=TRUE,return.marginals=TRUE,
                       return.marginals.predictor=TRUE,internal.opt=FALSE)
All six native INLA logs confirm tolerance_step=1e-05 and internal.opt=[0].
The separately logged tolerance_func remains 0.0005, unchanged from phase 2.
Original pilot/default fits remain the source of initial vectors; preflight
also checks equality to both immediately preceding diagnostic phases.

Results and stopping point
All six fits/scoring calls succeeded, optimizer statuses are 0, Hessian checks
have zero negative eigenvalues, and integration/summary validity pass.
Delta=2 Melding passes all screens. Baseline PSmelding still fails parameter
and prediction stability. This is the last authorized bounded fitting phase;
no further fits or production-setting changes are made by this diagnostic.

Maximum within-group spread                         baseline       delta=2
Parameter mean shift / pooled posterior SD          0.200251       0.00995179
Interval endpoint shift / pooled posterior SD       0.645217       0.0318773
Prediction mean RMS                                 0.0108707      0.000968191
Log-posterior mode spread                            0.00196114     0.000371993

For every corresponding fit, the saved parameter summaries, theta/x modes,
internal covariance, prediction means/SDs and log-posterior mode are exactly
identical to the fixed-backend phase. Changing tolerance.step had no measured
effect on these retained deterministic outputs. reports/inner_tolerance_effect.csv
checks each fit explicitly. The failed baseline endpoints cannot be explained
by CRPS sampling noise; parameter summaries and prediction means are deterministic.

The worst baseline endpoint is theta upper. Screens are unchanged: 0.1 pooled
SD, 0.01 prediction RMS and 0.01 log-mode spread. Boundary-probability spread
is zero; tabulated support is recorded separately. Predictor_strategy is NA
because internal optimization is disabled and no timing-selection line exists.
These two-dataset diagnostics do not establish stability across the full study.

Evidence
reports/stability_by_prior.csv and all_attempts.csv: full verdict/status.
reports/parameter_start_comparisons.csv and prediction_start_comparisons.csv:
all planned start comparisons. reports/fit_settings_audit.csv: exact controls
and native-log confirmation of inner tolerance/internal optimization.
reports/settings_hashes.csv and exact_fitting_controls.rds: exact controls,
initial values, seeds, likelihood data/matrices, mesh and SHA-256 digests.
reports/outer_vs_inner_stability.csv and corresponding attempts file:
comparison against the preceding fixed-backend phase.
reports/four_phase_stability.csv: original default, tighter outer tolerance,
fixed backend, and tighter inner tolerance together.
reports/inner_tolerance_effect.csv: exact equality checks versus phase 2.
preflight_checks.rds/txt/log validates the six inputs before fitting.
reports/reference_integrity.csv/txt verifies all previous diagnostic files.
sha256_manifest.csv inventories local artifacts and excludes itself.

Path provenance and preservation
Original execution root:
  /tmp/psmelding-tight-inner-20261022
Planned durable archive root:
  /Users/vinky/research-comps-clean/production_checks/numerical_tolerance_20261022/tight_inner
Previous phase execution roots:
  /tmp/psmelding-fixed-backend-20261022
  /tmp/psmelding-tight-tolerance-20261022
Their existing, verified durable archives:
  /Users/vinky/research-comps-clean/production_checks/numerical_tolerance_20261022/fixed_backend
  /Users/vinky/research-comps-clean/production_checks/numerical_tolerance_20261022/tight_tolerance
Default reference:
  /Users/vinky/research-comps-clean/production_checks/numerical_stability
Original pilot inputs:
  /Users/vinky/research-comps-clean/extension_test

Completed scripts, RDS files, recorded calls and metadata retain original
absolute /tmp paths. Copying does not rewrite them. If /tmp is removed, use
the durable root plus the original relative suffix: read a record's fit_file
as file.path(archive_root,'fits',basename(record$fit_file)); remap checkpoints,
predictions, posterior_draws and simulation_data similarly. path_mapping.csv
records the mappings. All raw and complete fit files are present; an original
provenance string does not indicate a missing artifact.

report_checks.R accepts archive root and fresh output directory arguments,
and reads checkpoints without fitting. Other audit scripts contain original
roots; use a scratch copy or remap paths in memory rather than altering the
preserved files. Do not execute run_checks.R merely to audit a copy.
The parent task handles durable copying and publishing. This local manifest
does not claim that a Drive upload has occurred.
