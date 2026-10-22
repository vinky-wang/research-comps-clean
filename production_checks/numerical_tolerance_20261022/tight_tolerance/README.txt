Tighter outer-tolerance diagnostic, 2026-10-22

Scope
Exactly six controlled fits: baseline PSmelding and delta=2 Melding, each with
area-variance starts 0.001, 0.1, 0.4 and point-variance start 0.1. The original
1676-vertex mesh, two saved datasets, prior B (separate Gaussian precisions,
pc.prec c(0.316,0.5)), latent/hyperparameter initial vectors, fit/scoring seeds,
1:1 threading, and 1000 joint latent draws are unchanged. The sole fitting
control change from the previous diagnostic is control.inla$tolerance=1e-4.
All six ran sequentially in one R 4.2 process. No full study was launched.

Results
All six fits/scoring calls succeeded; all optimizer statuses are 0, all Hessian
checks have zero negative eigenvalues, and integration/summary validity pass.
Delta=2 Melding passes all stability screens. Baseline PSmelding still fails
parameter and prediction stability. These are diagnostics on two datasets,
not proof of stability over the planned study.

Maximum within-group spread                         baseline       delta=2
Parameter mean shift / pooled posterior SD          0.181759       0.0216263
Interval endpoint shift / pooled posterior SD       0.536625       0.0470447
Prediction mean RMS                                 0.0110239      0.00255696
Log-posterior mode spread                            0.00199360     0.00128950

The previous default-0.005 mean shifts were 0.194782 and 0.272128, and endpoint
shifts were 0.505484 and 0.629806. Baseline's worst endpoint is theta upper.
The prior screens remain unchanged: 0.1 pooled SD, 0.01 prediction RMS and
0.01 log-mode spread. Boundary probability spread is zero; boundary support
and probabilities are tabulated separately rather than treated as proof of
numerical accuracy. Predictor strategy is recorded, not assumed causal.

Evidence
reports/stability_by_prior.csv: verdict and all component screens.
reports/all_attempts.csv: every planned attempt, including status and scores.
reports/parameter_start_comparisons.csv: all parameter/statistic/start pairs.
reports/prediction_start_comparisons.csv: all prediction comparisons.
reports/fit_settings_audit.csv: actual controls and runtime strategies.
reports/settings_hashes.csv and exact_fitting_controls.rds: exact controls,
initial values, seeds, likelihood data/matrices, mesh and SHA-256 digests.
reports/default_vs_tight_stability.csv and default_vs_tight_attempts.csv:
comparison to the original default-tolerance diagnostic.
preflight_checks.rds/txt/log: validation before fitting.
source_comparison.csv: seven current main source hashes; fitting, data,
configuration, evaluation and extension modules match the earlier diagnostic.
Only run_simulation.R has newer runner infrastructure; this diagnostic calls
fit_model directly and uses the same load_simulation module-loading routine.
reports/original_integrity.txt: original pilot inputs remained unchanged.
reports/previous_diagnostic_integrity.csv/txt: previous diagnostic preserved.
sha256_manifest.csv: local artifact file integrity (manifest excludes itself).

Path provenance and preservation
Original execution root:
  /tmp/psmelding-tight-tolerance-20261022
Planned durable archive root:
  /Users/vinky/research-comps-clean/production_checks/numerical_tolerance_20261022/tight_tolerance
Default-tolerance reference root:
  /Users/vinky/research-comps-clean/production_checks/numerical_stability
Original saved pilot root:
  /Users/vinky/research-comps-clean/extension_test

Completed scripts, RDS files, recorded calls and metadata deliberately retain
original absolute /tmp paths. Copying this directory does not rewrite them.
If /tmp is removed, use the durable root plus the original relative suffix:
for example a record's fit_file is read as file.path(archive_root, 'fits',
basename(record$fit_file)); apply the same prefix substitution to checkpoints,
predictions, posterior_draws and simulation_data. path_mapping.csv records
this mapping. These are provenance strings, not missing fit artifacts.

report_checks.R accepts the archive root and a fresh report output directory
as its first two arguments and reads all checkpoints without fitting.
The other retained audit scripts contain the original root. For a future
read-only audit, work on a scratch copy or remap that root and record paths
in memory; do not rewrite preserved artifacts. Do not execute run_checks.R
merely to audit a copy: it is the original fitting script.

The parent task handles durable copying and publishing. This directory and
its manifest describe the local diagnostic only; they do not claim a Drive
upload has occurred.
