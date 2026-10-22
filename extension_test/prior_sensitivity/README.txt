OBSERVATION-NOISE PRIOR SENSITIVITY

Start with findings.txt for the completed 36-fit analysis
This folder contains a small diagnostic experiment, not a production study
Original source code, simulated data and fitted objects remain unchanged

DESIGN
Two saved datasets: baseline gamma=1 and independent clustering delta=2
Both Melding and PSmelding, using the same saved finer mesh
Three priors applied independently to point and areal precisions
Three areal-variance starts per dataset/model/prior:0.001,0.1,0.4
See analysis_plan.txt for controls, practical screens and interpretation limits

CONTENTS
config.rds / plan.csv: analysis settings and 36 planned fits
simulation_data/: retained data, true fields, generating parameters and seeds
mesh.rds: common saved triangulation
fits/: all new complete INLA bundles, without replacing originals
checkpoints/: durable results for each attempted fit and scoring operation
logs/: original INLA diagnostic output for each attempt
predictions/: pointwise latent-surface summaries and prediction scores
posterior_draws/: 1000 latent-surface draws per fit used for CRPS
reports/: parameter estimates, probability near zero, MSE/CRPS and diagnostics
prior_reference/: analytic prior calibration and plots checked against INLA
review/: independent code/syntax checks
session_info.txt / runtime_defaults.rds: installed environment and defaults

RUNNING
From extension_test:
/Users/vinky/research-comps/R/run_r42 --vanilla prior_sensitivity/run_prior_sensitivity.R

Completed attempts, including failures, are skipped on restart
A saved full fit can finish an interrupted scoring step without refitting
All original hashes are checked after completion
Do not remove checkpoints to silently replace unwanted results
The driver stops before a new fit if less than 5 GiB remains free

INTERPRETATION
Use newly controlled prior-A runs as the reference, matched by starting value
Do not compare raw log-posterior-mode values across different priors
Do not treat optimizer starts as independent simulation replicates
CRPS targets latent Z only; no future observation noise is added
Approximate CRPS Monte Carlo uncertainty is separate from simulation-study
uncertainty, which cannot be estimated from one underlying field realization

REPORTS WITHOUT REFITTING
After all planned attempts finish, run these from extension_test:
/Users/vinky/research-comps/R/run_r42 --vanilla prior_sensitivity/report_prior_sensitivity.R
/Users/vinky/research-comps/R/run_r42 --vanilla prior_sensitivity/plot_posterior_priors.R
/Users/vinky/research-comps/R/run_r42 --vanilla prior_sensitivity/audit_saved_runs.R
/Users/vinky/research-comps/R/run_r42 --vanilla prior_sensitivity/audit_boundaries.R

KEY REPORT FILES
reports/stability_by_prior.csv: each numerical screen, with failures kept visible
reports/parameter_estimates.csv: all posterior summaries and generating truths
reports/areal_variance_summary.csv: focused areal summaries, probabilities and scores
reports/boundary_probabilities.csv: INLA-marginal probabilities near zero
reports/matched_prior_predictive_contrasts.csv: B/C versus A at matched starts
reports/parameter_absolute_spreads.csv: variation over initializations
reports/predictive_scores.csv: MSE and latent CRPS for every attempt
reports/posterior_noise_density_*.png: compare priors at the common start 0.1
reports/posterior_areal_cdf_starts_*.png: compare starts within every prior
reports/audit/: settings, likelihood mapping and fixed-model checks
reports/boundary_audit/: numerical support and boundary-CDF checks

Zero/one boundary probabilities are numerical marginal approximations, not
proof that the exact posterior tail probability is zero/one
