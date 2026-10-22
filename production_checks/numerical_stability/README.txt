PRODUCTION-MESH PRIOR-B NUMERICAL CHECK

Purpose
Check initialization sensitivity in the current extension fitting configuration
before any full simulation study is launched

Design
Two retained datasets from replicate 1: gamma=1 baseline (PSmelding), and
independent clustering delta=2 (Melding). The point and areal generating
variances are both .1, with I=100, J=25, N=2500 and rho=.2
No response surface or observations are regenerated

Six controlled fits: each dataset/model pair uses initial areal variance
.001, .1 and .4. Initial point variance is .1 in all fits. The other starting
hyperparameters and latent x are taken from that model's original same-mesh
pilot fit and held fixed across the three starts. INLA re-optimizes all of them
This checks the areal-variance initialization direction, not every possible
starting value or the untouched INLA default initialization

Current production settings
Observation precisions each independently use prior B, pc.prec c(.316,.5)
Separate point and areal observation variances are fitted
Mesh offset c(.1,.35), max.edge c(.05,.25), cutoff .01, with 1676 vertices
The regenerated mesh matches the original pilot coordinates and triangulation
The likelihood matrices/data and control.inla are asserted identical to the
original same-mesh pilot before fitting, apart from explicit prior/initial values
Current main control.inla uses variance VB; installed default tolerance=.005
Default integration/approximation strategies are auto, h=.005, safe=FALSE
No tighter tolerance, finer mesh, alternate integration or automatic retries
One INLA/BLAS/OpenMP thread, fixed fit and scoring seeds within each dataset

Source and outputs
R/ is the frozen main-engine snapshot taken for these checks
source_snapshot.rds contains MD5 hashes of those copied scripts
run_checks.R alters only the in-memory fitting function to insert initial values
It does not modify the copied source or the working project
config.rds, runtime_defaults.rds and session_info.txt record settings/environment
simulation_data/ retains both input datasets, including the expensive truth/data
fits/ retains complete model bundles and unmodified raw INLA returns
posterior_draws/ retains 1000 joint latent-surface draws per fit
predictions/ retains pointwise summaries and latent scores
checkpoints/ records each attempt, including failures, starts, seeds and diagnostics
reports/ contains the comparisons and input-integrity verification
Earlier results and source files are never overwritten

Screens
All six attempts must be recorded, fitted, scored and have valid summaries
Mode status=0, no negative Hessian eigenvalues, and no integration point more
than .001 log units above the reported mode
Within a dataset/model: log-mode spread <=.01; maximum pairwise posterior
mean/median/95% endpoint shift <=.1 pooled posterior SD; pointwise posterior
mean RMS difference <=.01; spread of P(areal variance<.001) <=.01
These are practical screens, not a proof of convergence or frequentist validity
CRPS Monte Carlo SE uses 10 blocks of 100 draws and is an approximate proxy
The three starts are diagnostic repeats, not independent simulation replicates

Commands
/Users/vinky/research-comps/R/run_r42 --vanilla /tmp/psmelding-production-numerics-20261022/run_checks.R
/Users/vinky/research-comps/R/run_r42 --vanilla /tmp/psmelding-production-numerics-20261022/report_checks.R /tmp/psmelding-production-numerics-20261022

This diagnostic folder does not authorize or start a full production run
