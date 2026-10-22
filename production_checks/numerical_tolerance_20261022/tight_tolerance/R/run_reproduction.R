# ==================================================
# Reproduce Simulations in Zhong et al (2025) 
# ==================================================

# Run from the research-comps-clean project root
source("R/run_simulation.R", local = TRUE)

reproduction_settings <- list(
  generating = list(
    I = c(100L, 250L),
    J = c(4, 25, 100),
    rho = c(.1, .2, .4),
    sigma2 = 1,
    gamma = c(0, 1),
    noise = data.frame(tau2_s = .1, tau2_B = .1) 
  ),
  # [ORIGINAL] Gamma(shape, rate) prior for observation precision q = 1/tau^2
  fitting = list(
    noise_prior = list(prior = "loggamma", param = c(1, 5e-5))
  ),
  models = c("Melding", "PSgeo", "PSmelding"),
  n_replicates = 100L,
  base_seed = 2026L
)

reproduction_config <- function(
  section = "reproduction",
  project_dir = getwd(),
  settings = reproduction_settings
) {
  engine <- load_simulation()
  config <- engine$study_config(settings, section, project_dir)
  config$study <- "Corrected Zhong reproduction"
  config$section <- section
  config$output_dir <- file.path(config$project_dir, "results", "reproduction")
  config
}

# Run from the command line
if (sys.nframe() == 0L) {
  study_main(
    commandArgs(trailingOnly = TRUE),
    reproduction_config,
    normalizePath("R/run_reproduction.R", mustWork = TRUE),
    sections = "reproduction",
    default_section = "reproduction"
  )
}
