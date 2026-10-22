TWO-LAPTOP EXTENSION RUNNER
Version: production_v2_20261022
Status: VALIDATION_PENDING — full simulations have not started

What changed
The seven main R scripts remain in R/. Only run_simulation.R changed
The full R=100 plan is built before selecting global replicate IDs, preserving
case IDs, seeds and common latent fields. Named runs use five replicates per
batch by default. Each experiment has its own output directory within a batch

Planned layout
results/extension/production_v2_20261022/
  laptop_1/batch_001_005/strength/ ... all four experiments ...
  laptop_1/batch_006_010/ ... through batch_046_050/
  laptop_2/batch_051_055/ ... through batch_096_100/
  manifests/laptop_1_batch_001_005_strength.csv ...

Each experiment directory retains simulation_data/, inla_fits/, predictions/,
diagnostics/, config.rds, session_info.txt and CSV summaries, together with
file_manifest.csv and batch_manifest.rds. Central records describe local
verification; they do not claim the upload to Drive has been verified

Preview the assignments now (no fitting)
Run from the project root using the configured R 4.2 environment
On this laptop its launcher is /Users/vinky/research-comps/R/run_r42
On the second laptop use its R 4.2 Rscript and matching package library

Laptop 1:
Rscript R/run_extension.R --run-id=production_v2_20261022 --laptop=laptop_1 --replicate-start=1 --replicate-end=50 --batch-size=5

Laptop 2:
Rscript R/run_extension.R --run-id=production_v2_20261022 --laptop=laptop_2 --replicate-start=51 --replicate-end=100 --batch-size=5

The default extension section is all four experiments. Each laptop has 1450
planned fits, 145 per five-replicate batch across the four experiments
Keep --replicates at its default 100; it sets the full study size, not a range
Named non-test laptop runs require an explicit start and end
--output-root=PATH can choose a different parent of the versioned run folder
Neither a laptop name nor its assigned range changes base_seed=2026

Starting and resuming AFTER numerical clearance
Add --run to the same command. Repeating that exact command resumes its named
directories. Keep run ID, ranges, batch size, source code and settings fixed
Completed jobs are verified and skipped. A complete saved fit can finish scoring
without refitting. Recorded numerical failures remain failures. If interruption
occurred before a complete fit bundle was saved, any raw return is retained and
the attempt is marked interrupted for review; it is not automatically refitted
If scoring resumes, elapsed_scope identifies resumed-session timing and the
original INLA runtime is retained separately when available

Source hashes, scientific settings, R/package versions and RNG kind must match
before resume. Checksums protect saved datasets and completed job outputs
A lock prevents two processes from writing the same experiment directory
Abandoned locks are recovered only on the original host after its process ends
Do not delete checkpoints or locks to force a rerun; preserve attempts

Verification
Replace --run with --verify-only in the same named command after it completes
This checks settings, source/runtime identity, every expected job, file coverage
and SHA256 hashes. A complete-with-failures batch remains visibly unsuccessful
After uploading, verify downloaded/copied bytes against the manifests before
marking the Drive copy verified. Do not let both laptops append to one shared CSV
Each laptop writes separate central records and separate output directories

Two-laptop reproducibility check
production_checks/check_laptop.R performs two baseline fits (Melding and
PSmelding) at full production resolution and uses only global replicate 1
It is separate from the production plan. All generated data/fits are retained

Laptop 1 (completed local check is saved under production_checks/laptop_1_check):
Rscript production_checks/check_laptop.R --run --output-dir=production_checks/laptop_1_check

Copy that check folder to laptop 2, then run:
Rscript production_checks/check_laptop.R --run --output-dir=production_checks/laptop_2_check --compare=/path/to/copied/laptop_1_check

The comparison checks source/settings/runtime, actual data/seeds, mesh and priors
exactly, and writes parameter and score comparisons for review. Matching inputs
do not imply bit-identical INLA posterior calculations on different computers
Record and review numerical differences; no automatic production launch occurs

Remaining numerical issue
Six controlled prior-B fits on the original 1676-vertex mesh all returned mode
status 0, but baseline PSmelding and delta=2 Melding both failed initialization-
stability screens. Prior B alone has not cleared the current numerical settings
See numerical_stability/findings.txt and its preserved six-fit audit
The earlier finer-mesh/tighter-tolerance results cannot clear this configuration
Keep the full launch pending numerical follow-up and the second-laptop check

Preservation
This version is a validation snapshot. Changes to model/numerical settings after
these checks need a new frozen snapshot and a new run ID. Original production_v1
Drive artifacts and earlier diagnostic fits remain unchanged
