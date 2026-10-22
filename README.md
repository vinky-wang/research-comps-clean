# research-comps-clean

R code for reproducing and extending Zhong et al. (2025), with separate reproduction and extension entry scripts

**Status: development snapshot, 22 October 2026. Full extension simulations have not started.** The two-laptop batching and output preservation are implemented. Numerical acceptance is still pending for baseline PSmelding, and laptop 2 has not yet been independently checked

## Clone and environment

This is a private repository; authenticate with a GitHub account that has access

```sh
git clone https://github.com/vinky-wang/research-comps-clean.git
cd research-comps-clean
```

The reference environment is **R 4.2.3, Apple silicon/arm64, INLA 23.04.24**. Package versions are recorded in [environment/package_versions.csv](environment/package_versions.csv) and the reference session in [environment/session_info.txt](environment/session_info.txt)

Cloning does not install R, the package library, or native libraries. Use the matching pinned environment on laptop 2 and set its actual package-library path. The fallback path in `R/config.R` belongs to laptop 1

```sh
export RESEARCH_COMPS_R_LIB="/absolute/path/to/.r-lib/4.2-arm64"
export R_LIBS_USER="$RESEARCH_COMPS_R_LIB"
RC_RSCRIPT="/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/bin/Rscript"
```

The launcher in the older `research-comps` repository is not required when using this explicit Rscript path. The configured macOS runtime also expects GCC libraries under `/opt/homebrew/opt/gcc/lib/gcc/current`

## Main scripts

| File | Purpose |
| --- | --- |
| `R/config.R` | Configuration, validation, scenario plans and seeds |
| `R/generate_data.R` | Latent fields, sampling locations, point and areal observations |
| `R/fit_models.R` | Mesh, SPDE representation and INLA model fitting |
| `R/evaluate.R` | Posterior summaries, parameter recovery and latent-surface scores |
| `R/run_simulation.R` | Shared execution, output saving, batching and resume verification |
| `R/run_reproduction.R` | Corrected reproduction settings |
| `R/run_extension.R` | Four extension experiments |

Reproduction retains the original observation-precision Gamma prior. Extensions use prior B, `pc.prec` with `c(0.316, 0.5)`, for each separately estimated observation precision. The model and numerical settings in the main scripts have not been replaced with the diagnostic trial settings

## Preview laptop assignments without fitting

Run from the project root. These commands omit `--run`, so they only preview the work

```sh
# Laptop 1: global replicates 1–50
"$RC_RSCRIPT" --vanilla R/run_extension.R \
  --run-id=production_v2_20261022 --laptop=laptop_1 \
  --replicate-start=1 --replicate-end=50 --batch-size=5

# Laptop 2: global replicates 51–100
"$RC_RSCRIPT" --vanilla R/run_extension.R \
  --run-id=production_v2_20261022 --laptop=laptop_2 \
  --replicate-start=51 --replicate-end=100 --batch-size=5
```

Each assignment covers all four experiments and 1,450 model fits. Keep the full study size at 100 replicates and the same base seed on both laptops. A new accepted numerical configuration will require a new frozen run version before production

## Laptop comparison and current numerical findings

The clone includes `production_checks/laptop_1_check/laptop_check.rds`, a small reference containing input identity, data, mesh and summary results. `production_checks/check_laptop.R --compare` needs only this record

When an independent baseline check is appropriate, use a fresh laptop-2 output directory that has never contained copied laptop-1 fits

```sh
"$RC_RSCRIPT" --vanilla production_checks/check_laptop.R --run \
  --output-dir=production_checks/laptop_2_independent_check \
  --compare=production_checks/laptop_1_check
```

This performs two baseline fits with the current main-code settings. It does not launch the study or automatically approve numerical agreement. A copied completed output directory can resume without fitting, so it must not be used for the independent check

The historical `laptop_preflight_check.R` and the full preflight in `laptop_2_apple_silicon.txt` also verify the original `fit_run` inventory. That large directory is excluded from Git and must be transferred separately to use those inventory checks. The compact reference alone is sufficient for the comparison command above

The [latest numerical findings](production_checks/numerical_tolerance_20261022/findings.txt) supersede the earlier numerical-status paragraph in `production_checks/README.txt`. Tighter tolerance cleared the tested delta=2 Melding group; baseline PSmelding still failed the initialization-stability screens. Do not interpret successful startup or optimizer status as production clearance

## Preserved work and generated files

`extension_test/` and `production_checks/` retain diagnostic scripts, archived code snapshots, notes and compact result tables. Historical scripts and manifests intentionally retain original laptop-1 or `/tmp` paths and references to large artifacts; consult their path mappings before rerunning them

Git excludes full INLA objects, latent surfaces, pointwise predictions, posterior draws, large generated run artifacts, package libraries and article PDFs. These files remain on laptop 1; no results were deleted. Existing manifests describe the full local diagnostic archives, not only the subset included in this clone

The latest compact reports are also available in the private [Google Drive analysis folder](https://drive.google.com/drive/folders/1h1h-0CJ1fiBJ7nxRKUbjYD7Ubwy7jBsC). The Drive reports archive does not include the full fitted objects
