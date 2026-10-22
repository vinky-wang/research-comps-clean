# Run from the project root with the same R 4.2 environment on both laptops
# This check uses one baseline replicate at full production resolution
source("R/run_extension.R")

args <- commandArgs(trailingOnly = TRUE)
allowed <- args == "--run" | grepl("^--(output-dir|compare)=.+$", args)
keys <- sub("=.*$", "", args)
if (any(!allowed) || anyDuplicated(keys)) {
  stop("Use --output-dir=PATH [--run] [--compare=OTHER_LAPTOP_CHECK]")
}
option <- function(name, default = NULL) {
  value <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(value)) sub("^[^=]+=", "", value) else default
}
output <- path.expand(option("output-dir", "production_checks/local_check"))
reference <- option("compare")
engine <- load_simulation()
config <- extension_config("strength")
config$generating$gamma <- 1
config$replicate_ids <- 1L
config$entry_script <- normalizePath("R/run_extension.R")
config$output_dir <- output
dir.create(output, recursive = TRUE, showWarnings = FALSE)
packages <- engine$configure_runtime(config)
versions <- runner_package_versions(packages)

if ("--run" %in% args) {
  result <- run_simulation(config, engine, run_dir = file.path(output, "fit_run"))
  if (any(!result$status$scored)) stop("The laptop check contains failed fits")
}
if (!file.exists(file.path(output, "fit_run", "batch_manifest.rds"))) {
  stop("Run this check with --run first; it performs two baseline fits")
}
fit_dir <- file.path(output, "fit_run")
verify_completed_run(fit_dir)
saved <- readRDS(file.path(fit_dir, "config.rds"))
current_identity <- run_identity(config, runner_source_hashes(), versions)
if (!identical(saved$identity, current_identity)) stop("The check's code or runtime changed")
dataset <- readRDS(file.path(fit_dir, "simulation_data", paste0(saved$plan$case_id, ".rds")))
mesh <- saved$fitting_setup$mesh
record <- list(
  identity = current_identity,
  data = list(seeds = dataset$seeds, Z = dataset$truth$Z,
    coordinates = dataset$truth$coordinates,
    point_coordinates = sf::st_coordinates(dataset$points),
    point_response = dataset$points$value, point_truth = dataset$points$truth,
    area_response = dataset$areas$value, area_truth = dataset$areas$truth),
  mesh = list(loc = mesh$loc, triangles = mesh$graph$tv),
  priors = saved$fitting_setup$priors,
  parameters = read.csv(file.path(fit_dir, "parameter_estimates.csv")),
  scores = read.csv(file.path(fit_dir, "predictive_scores.csv")),
  diagnostics = read.csv(file.path(fit_dir, "fit_diagnostics.csv")),
  session = sessionInfo()
)
record_file <- file.path(output, "laptop_check.rds")
if (!file.exists(record_file)) save_new_rds(record, record_file)
if (!is.null(reference)) {
  other <- readRDS(file.path(path.expand(reference), "laptop_check.rds"))
  checks <- data.frame(check = c("code/settings/runtime", "data/seeds", "mesh", "priors"),
    identical = c(identical(record$identity, other$identity),
      identical(record$data, other$data), identical(record$mesh, other$mesh),
      identical(record$priors, other$priors)))
  atomic_write_csv(checks, file.path(output, "laptop_comparison.csv"))
  parameters <- merge(record$parameters, other$parameters, by = c("model", "parameter"),
    suffixes = c("_local", "_reference"))
  scale <- sqrt((parameters$posterior_sd_local^2 + parameters$posterior_sd_reference^2) / 2)
  parameters$mean_difference_pooled_sd <- abs(parameters$posterior_mean_local -
    parameters$posterior_mean_reference) / scale
  atomic_write_csv(parameters, file.path(output, "parameter_comparison.csv"))
  atomic_write_csv(merge(record$scores, other$scores, by = c("model", "region", "metric"),
    suffixes = c("_local", "_reference")), file.path(output, "score_comparison.csv"))
  print(checks, row.names = FALSE)
  if (!all(checks$identical)) stop("Laptop inputs or environments differ; inspect laptop_comparison.csv")
  message("Inputs match; review parameter/score differences and diagnostics before production")
} else {
  message("Saved ", record_file, "; copy this check folder to laptop 2 for --compare")
}
