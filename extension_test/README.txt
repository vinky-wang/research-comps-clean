EXTENSION TEST

This isolated folder contains an unchanged snapshot of the six extension R modules.
The main research-comps-clean/R scripts and existing results are not modified.

Run from this folder with the existing R 4.2 launcher:
/Users/vinky/research-comps/R/run_r42 --vanilla run_tests.R --phase=all

Available phases: pilot, stability, report, all (default)

Pilot: one replicate of each of the 16 extension conditions, 29 model fits
D=[0,1]^2, I=100, J=25, rho=.2, sigma2=1, nu=1, N=2500
Original production mesh and 1000 posterior draws are retained

Stability: baseline gamma=1 (both models), sign reversal (PSmelding), and
independent clustering delta=2 (both models), for five additional fits
Both mesh maximum-edge values and cutoff are multiplied by .75
Domain, data, priors, seeds, numerical controls and other settings are unchanged
Five pilot fits are rescored with 5000 posterior draws without refitting

Results are under results/pilot and results/finer_mesh
checkpoints contains one durable result per attempted job and CRPS check
reports contains diagnostic tables, paired comparisons, and resource estimates
All latent fields, predictions and full fits are retained

RESTARTS AND BATCHING
Re-running the same command skips jobs already recorded, including failures
An interrupted job with a complete saved fit can finish scoring without refitting
An interrupted job with only a raw checkpoint is preserved and flagged for inspection
No failed fit is automatically retried or silently replaced
Do not edit the snapshot after starting: verify_snapshot checks source hashes
Do not delete checkpoints to manufacture additional replicates

The main repository runner still has no automatic resume or replicate-offset option
--replicates=N means replicates 1 through N, not N additional independent draws
For production, use one R=100 run unless explicit disjoint replicate IDs are implemented
This test's job checkpoints demonstrate restartable execution without changing production

RESOURCE PROTECTION
The driver checks free disk before each model job and stops below 5 GiB
Runtime/storage projections are linear extrapolations, not guarantees
The pilot saves settings per physical job, so it has slightly more metadata overhead
than the main runner, which saves settings once per experiment
RESULTS AND AUDIT
Start with reports/findings.txt for the completed test findings and remaining issues
reports/restart_check.txt records the verified completed-job restart check
audit_outputs.R inspects saved outputs without refitting models
No full 100-replicate production study was started

FOLLOW-UP NUMERICAL INVESTIGATION
See investigation/findings.txt for the saved-posterior, areal-projection and
21-fit controlled optimizer investigation; all originals remain preserved

OBSERVATION-NOISE PRIOR SENSITIVITY
See prior_sensitivity/findings.txt for the completed 36-fit comparison of the
original and two PC observation-noise priors, with all starts and results retained
