# PSMelding: Simulation Reproduction and Extensions

This repository contains code for reproducing and extending the selected work of [Zhong et al. (2025)](https://github.com/ruimanzhong/PSmelding) as part of the PhD research comprehensive examination in the Department of Statistical Sciences at the University of Toronto.

The study compares three models:

- **PS Melding:** combines point and areal observations and models preferential sampling locations.
- **Melding:** combines point and areal observations without modeling preferential sampling.
- **PS Geo:** models preferentially sampled point observations without areal data.

## Project structure

```text
functions/
  generate_data.R          # Simulate spatial fields, point data and areal data
  fit_models.R             # Build the mesh, projections and INLA models
  create_summary.R         # Evaluate fits and construct tables/figures

scripts/
  run_reproduction.R       # Reproduce Zhong et al.'s simulation settings
  run_extensions.R         # Additional simulation experiments
  prerun_extensions.R      # Optional mesh/prior sensitivity checks
  summarize_results.R      # Summarize completed fits and save results
```

Run all commands from the **project root**, which contains `functions/` and `scripts/`. The `results/` directories are created by the scripts.

## Requirements

- R and the **INLA** package.
- R packages used in the code: `spatstat.geom`, `spatstat.random`, `Matrix`, `sf`, `ggplot2`, `patchwork` and `scales`.

## Simulation design

The spatial domain is $D=[0,1]^2$. Data are generated from a Matérn Gaussian field (smoothness $\nu=1$), evaluated on a $50\times50$ grid. The default setting uses $\mu=0$, $\sigma^2=1$, $\rho=0.2$, $I=100$ point observations, $J=25$ areal observations, and observation-error variances $\tau_s^2=\tau_B^2=0.1$.

| Study | Conditions examined |
|---|---|
| Reproduction | Spatial range $\rho$, preferentiality $\gamma$, point count $I$, and areal count $J$ |
| Preferentiality strength | Vary $\gamma$ under the constant-preferentiality model |
| Areal measurement noise | Vary $\tau_B^2$ |
| Sampling mechanism | Constant $\gamma$, spatially varying $\gamma(s)$, and an additional independent latent field |

The scripts define the exact parameter grids and replication counts at the top of each file. The main simulation runs are configured for **100 replications per condition**, using seeds beginning at 123.

The reproduction uses the authors' mesh, likelihood approximation and prior choices. The extensions use the fitting specifications selected from the preliminary sensitivity checks.

## Run the analyses

1. **Run the reproduction:**

   ```bash
   Rscript scripts/run_reproduction.R
   ```

2. **Run the extensions:**

   ```bash
   Rscript scripts/run_extensions.R
   ```

3. **Update the summary tables and figures:**

   In `scripts/summarize_results.R`, set:

   ```r
   update_tables = TRUE
   create_plots = TRUE
   update_reproduction_tables = TRUE
   update_report_tables = TRUE
   create_report_plots = TRUE
   ```

   Then run:

   ```bash
   Rscript scripts/summarize_results.R
   ```


## Reference

Zhong et al. (2025), *PSmelding* ([paper](https://arxiv.org/abs/2309.03316); [original code](https://github.com/ruimanzhong/PSmelding)).
