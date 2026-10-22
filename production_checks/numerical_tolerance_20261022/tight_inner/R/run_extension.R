# ==================================================
# Extend Simulations in Zhong et al (2025)
# ==================================================

# Run from the research-comps-clean project root
source("R/run_simulation.R", local = TRUE)

# Shared settings for all four experiments
extension_settings <- list(
  generating = list(
    x_limits = c(0, 1),
    y_limits = c(0, 1),
    grid_n = 50L,  # N = 2500 simulation and prediction locations
    I = 100L,
    J = 25L,  # 5 x 5 areas
    rho = .2,
    sigma2 = 1,
    nu = 1,
    mu = 0,
    alpha = .05,  # Cancels when sampling a fixed number of locations
    gamma = 1,
    sampling = NULL,
    noise = data.frame(tau2_s = .1, tau2_B = .1)
  ),
  # [EXTENSION] Prior B for both separately estimated observation precisions
  # P(observation SD > .316) = .5, giving median variance .099856
  fitting = list(
    noise_prior = list(prior = "pc.prec", param = c(.316, .5))
  ),
  models = c("Melding", "PSmelding"),
  n_replicates = 100L,
  base_seed = 2026L,

  # Each experiment changes only the settings listed below
  experiments = list(
    # Preferential strength, including gamma=0 as a non-preferential reference
    strength = list(
      generating = list(gamma = c(0, .5, 1, 1.5, 2))
    ),

    # Areal measurement precision with fixed point variance and gamma=1
    areal_noise = list(
      generating = list(
        noise = data.frame(
          noise_id = c("area_0025", "area_005", "area_01", "area_02", "area_04"),
          tau2_s = .1,
          tau2_B = c(.025, .05, .1, .2, .4)
        )
      )
    ),

    # Spatially varying preferentiality, fitted with one constant gamma
    spatial_gamma = list(
      models = "PSmelding",
      generating = list(
        # x < split_x uses gamma_L; x >= split_x uses gamma_R
        # NA marks parameters unused by this sampling mechanism
        sampling = data.frame(
          sampling_id = c("constant_1", "mild_heterogeneity", "sign_reversal"),
          sampling_type = "piecewise",
          gamma = NA_real_,
          gamma_L = c(1, .5, -.5),
          gamma_R = c(1, 1.5, 2.5),
          split_x = .5,
          delta = 0,
          rho_U = NA_real_,
          sigma2_U = NA_real_
        )
      )
    ),

    # Independent spatial clustering, fitted without the additional U field
    clustering = list(
      generating = list(
        # Unit-variance independent U gives Corr(log intensity, phi)=1/sqrt(1+delta^2)
        sampling = data.frame(
          sampling_id = c("delta_0", "delta_1", "delta_2"),
          sampling_type = "independent_field",
          gamma = 1,
          gamma_L = NA_real_,
          gamma_R = NA_real_,
          split_x = NA_real_,
          delta = c(0, 1, 2),
          rho_U = .2,
          sigma2_U = 1
        )
      )
    )
  )
)

extension_config <- function(
  section = "strength",
  project_dir = getwd(),
  settings = extension_settings
) {
  engine <- load_simulation()
  config <- engine$study_config(settings, section, project_dir)
  config$study <- "Preferential sampling extensions"
  config$section <- section
  config$output_dir <- file.path(config$project_dir, "results", "extension", section)
  config
}

# Run from the command line
if (sys.nframe() == 0L) {
  study_main(
    commandArgs(trailingOnly = TRUE),
    extension_config,
    normalizePath("R/run_extension.R", mustWork = TRUE),
    sections = names(extension_settings$experiments),
    default_section = "all",
    test_all_conditions = TRUE
  )
}
