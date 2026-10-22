stopifnot(
  as.character(getRversion()) == "4.2.3",
  R.version$platform == "aarch64-apple-darwin20",
  dir.exists(Sys.getenv("RESEARCH_COMPS_R_LIB"))
)
source("R/run_extension.R")
engine <- load_simulation()
config <- extension_config("strength")
config$generating$gamma <- 1
config$replicate_ids <- 1L
packages <- engine$configure_runtime(config)
reference_dir <- Sys.getenv("RC_REFERENCE")
reference <- readRDS(file.path(reference_dir, "laptop_check.rds"))
reference_status <- verify_completed_run(file.path(reference_dir, "fit_run"))
stopifnot(all(reference_status$scored))
identity <- run_identity(config, runner_source_hashes(),
                         runner_package_versions(packages))
stopifnot(identical(identity, reference$identity))
print(Sys.info())
print(sessionInfo())
print(sf::sf_extSoftVersion())
print(INLA::inla.getOption("inla.call"))
cat("Input identity and transferred fit inventory match; no fits performed.\n")
