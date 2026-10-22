EXTENSION NUMERICAL INVESTIGATION

Start with findings.txt for the completed conclusions and remaining issues
Original project R files, simulated data and fitted objects are preserved
All new optimizer attempts are retained; no preferred result replaces an original

CONTENTS
posterior/: saved hyperparameter distributions, boundary checks, prior reference,
            original optimizer diagnostics and figures
projection/: true areal averages versus equal-vertex and exact FE projections,
             area-level results, validation and figures
optimizer/: 21 planned controlled fits, starting values, logs and completion records
reports/: controlled-fit summaries, pairwise comparisons and integrity result
methods/: independent check of posterior summaries and optimizer diagnostics
optimizer_plan.txt: fixed investigation design recorded before these runs
run_optimizer.R: isolated controlled-fit driver
plot_optimizer.R: figure of areal-variance summaries across starts

REPRODUCE FROM extension_test
Use /Users/vinky/research-comps/R/run_r42 --vanilla followed by one of:
  investigation/posterior/investigate_saved.R
  investigation/projection/projection_diagnostics.R
  investigation/run_optimizer.R
  investigation/methods/compare_optimizer.R
  investigation/plot_optimizer.R

The optimizer driver skips completed records, including failures
Never delete records to silently replace unfavorable outcomes
The diagnostic scripts rewrite only their own reports/figures
No full simulation study is started by these scripts

INTERPRETATION
One replicate diagnoses workflow/numerical sensitivity; it cannot establish
parameter recovery, empirical coverage, or model performance across simulations
Exact FE integration here applies to reconstructed nodal truth, not the
unobserved continuous generating field between raster cell centres
Tightening hyperparameter tolerance does not tighten all INLA approximations
Near-zero variance is not an imposed numerical hard bound
Full fitted objects are retained so these questions can be revisited
