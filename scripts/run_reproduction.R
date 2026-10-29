source("functions/generate_data.R")
source("functions/fit_models.R")


# Set up reproduction scenarios -----------------------------------------------

rho_values = c(0.1, 0.2, 0.4)
gamma_values = c(0, 1)
I_values = c(100L, 250L)
J_values = c(4L, 25L, 100L)

scenarios = expand.grid(
  rho = rho_values,
  gamma = gamma_values,
  I = I_values,
  J = J_values
)


# Model specifications ---------------------------------------------------------

# Authors' mesh, likelihood, and prior specifications
model_names = c("ps_melding", "melding", "ps_geo")
areal_models = c("ps_melding", "melding")

mesh_max_edge = c(0.05, 0.25)
mesh_offset = c(0.1, 0.35)
mesh_cutoff = 0.01

boundary = rbind(
  c(0, 0),
  c(1, 0),
  c(1, 1),
  c(0, 1),
  c(0, 0)
)

mesh = INLA::inla.mesh.2d(
  loc.domain = boundary,
  max.edge = mesh_max_edge,
  offset = mesh_offset,
  cutoff = mesh_cutoff
)

point_process_method = "mesh_dual"
areal_method = "vertices"

separate_noise = FALSE
prior_gamma = c(1, 10)
prior_noise = NULL
prior_noise_B = NULL

# Authors' ordinary Matérn SPDE with a sum-to-zero constraint
create_spde = function(mesh, prior_range = NULL, prior_sigma = NULL) {
  stopifnot(is.null(prior_range), is.null(prior_sigma))

  INLA::inla.spde2.matern(
    mesh = mesh,
    alpha = 2,
    constr = TRUE
  )
}


# Replications -----------------------------------------------------------------

n_replications = 100L
replications = seq_len(n_replications)

simulation_seeds = 123L + replications - 1L

output_dir = file.path("results", "reproduction")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# Record simulation settings ---------------------------------------------------

run_info = list(
  r_version = R.version.string,
  inla_version = INLA::inla.version("version"),
  session_info = utils::sessionInfo(),

  scenarios = scenarios,
  models = model_names,

  n_replications = n_replications,
  assigned_replications = replications,
  simulation_seeds = simulation_seeds,

  generating = list(
    mu = 0,
    sigma2 = 1,
    nu = 1,
    alpha = 0.05,
    tau2_s = 0.1,
    tau2_B = 0.1,
    scenario = "constant"
  ),

  mesh_max_edge = mesh_max_edge,
  mesh_offset = mesh_offset,
  mesh_cutoff = mesh_cutoff,
  mesh_vertices = mesh$n,

  spde = list(
    alpha = 2,
    constr = TRUE,
    prior_range = NULL,
    prior_sigma = NULL
  ),

  point_process_method = point_process_method,
  areal_method = areal_method,
  separate_noise = separate_noise,

  prior_gamma = prior_gamma,
  prior_noise = prior_noise,
  prior_noise_B = prior_noise_B,

  ps_geo_reused_across_J = TRUE,
  ps_geo_reference_J = J_values[1]
)

run_info_path = file.path(output_dir, "run_info.rds")

# Prevent mixing results from incompatible simulation specifications
if (file.exists(run_info_path)) {

  previous_info = readRDS(run_info_path)

  settings = c(
    "scenarios",
    "models",
    "n_replications",
    "assigned_replications",
    "simulation_seeds",
    "generating",
    "mesh_max_edge",
    "mesh_offset",
    "mesh_cutoff",
    "mesh_vertices",
    "spde",
    "point_process_method",
    "areal_method",
    "separate_noise",
    "prior_gamma",
    "prior_noise",
    "prior_noise_B",
    "ps_geo_reused_across_J",
    "ps_geo_reference_J"
  )

  for (setting in settings) {
    if (!identical(previous_info[[setting]], run_info[[setting]])) {
      stop(
        "Incompatible reproduction settings: ",
        setting,
        ". Choose a new output directory before fitting."
      )
    }
  }

} else {
  saveRDS(run_info, file = run_info_path)
}


# Check previously saved replications -----------------------------------------

# A replication is complete only if both its fit and log are readable
is_complete = function(condition_dir, replication) {

  filename = sprintf("rep_%03d.rds", replication)

  fit_path = file.path(condition_dir, "fits", filename)
  log_path = file.path(condition_dir, "log", filename)

  if (!file.exists(fit_path) || !file.exists(log_path)) {
    return(FALSE)
  }

  saved = tryCatch(
    readRDS(fit_path),
    error = function(e) NULL
  )

  fit_log = tryCatch(
    readRDS(log_path),
    error = function(e) NULL
  )

  if (is.null(saved) ||
      !is.list(saved$results$models) ||
      !all(model_names %in% names(saved$results$models)) ||
      !is.data.frame(fit_log) ||
      !all(model_names %in% fit_log$model)) {
    return(FALSE)
  }

  all(vapply(
    saved$results$models[model_names],
    function(x) {
      is.list(x) &&
        !is.null(x$fit) &&
        !is.null(x$prediction)
    },
    logical(1)
  ))
}


# Fit and save models ----------------------------------------------------------

# Group by rho, gamma, and I so PS Geo can be reused across J
configurations = unique(scenarios[c("rho", "gamma", "I")])

for (k in seq_len(nrow(configurations))) {

  configuration = configurations[k, ]

  rho = configuration$rho
  gamma = configuration$gamma
  I = configuration$I

  message(
    "Configuration ", k, " of ", nrow(configurations),
    ": rho = ", rho,
    ", gamma = ", gamma,
    ", I = ", I
  )

  for (replication in replications) {

    message(
      "  Replication ", replication,
      " of ", n_replications
    )

    seed = simulation_seeds[replication]

    reference = NULL
    logged_geo = NULL

    for (j in seq_along(J_values)) {

      J = J_values[j]

      condition_dir = file.path(
        output_dir,
        paste0("rho_", rho),
        paste0("gamma_", gamma),
        paste0("I_", I, "_J_", J)
      )

      # Skip completed replications
      if (is_complete(condition_dir, replication)) {

        # Retrieve PS Geo from the first J for reuse
        if (j == 1L) {

          fit_path = file.path(
            condition_dir,
            "fits",
            sprintf("rep_%03d.rds", replication)
          )

          saved = readRDS(fit_path)
          reference = saved$data

          geo_log = saved$log[
            saved$log$model == "ps_geo",
            ,
            drop = FALSE
          ]

          geo_log$reused = NULL
          geo_log$rho = NULL
          geo_log$gamma_true = NULL
          geo_log$I = NULL
          geo_log$J = NULL

          logged_geo = list(
            models = saved$results$models["ps_geo"],
            log = geo_log
          )

          rm(saved, geo_log)
        }

        next
      }

      # Generate data
      simulated = generate_data(
        rho = rho,
        gamma = gamma,
        I = I,
        J = J,
        tau2_s = 0.1,
        tau2_B = 0.1,
        seed = seed,
        scenario = "constant"
      )

      # Verify the same point observations and latent field across J
      if (is.null(reference)) {

        reference = simulated

      } else {

        stopifnot(
          isTRUE(all.equal(
            simulated$points,
            reference$points
          )),
          isTRUE(all.equal(
            simulated$prediction,
            reference$prediction
          ))
        )
      }

      # Fit PS Geo only for the first J
      models_to_fit = if (j == 1L) {
        model_names
      } else {
        areal_models
      }

      logged = fit_models_logged(
        data = simulated,
        models = models_to_fit,
        mesh = mesh,
        replication = replication,
        seed = seed,
        point_process_method = point_process_method,
        areal_method = areal_method,
        separate_noise = separate_noise,
        prior_gamma = prior_gamma,
        prior_noise = prior_noise,
        prior_noise_B = prior_noise_B
      )

      if (j == 1L) {

        # Store PS Geo for reuse at larger J
        geo_log = logged$log[
          logged$log$model == "ps_geo",
          ,
          drop = FALSE
        ]

        logged_geo = list(
          models = logged$results$models["ps_geo"],
          log = geo_log
        )

        results = logged$results
        fit_log = logged$log

        fit_log$reused = FALSE

      } else {

        # Reuse the PS Geo fit from J = 4
        if (is.null(logged_geo) ||
            is.null(logged_geo$models[["ps_geo"]])) {
          stop(
            "Cannot reuse PS Geo at replication ",
            replication
          )
        }

        results = logged$results

        results$models = c(
          results$models,
          logged_geo$models
        )

        geo_log = logged_geo$log

        geo_log$reused = NULL
        geo_log$rho = NULL
        geo_log$gamma_true = NULL
        geo_log$I = NULL
        geo_log$J = NULL

        fit_log = rbind(
          logged$log,
          geo_log
        )

        fit_log$reused = fit_log$model == "ps_geo"
      }

      # Record simulation parameters
      fit_log$rho = rho
      fit_log$gamma_true = gamma
      fit_log$I = I
      fit_log$J = J

      rownames(fit_log) = NULL

      # Save replication
      save_replication(
        data = simulated,
        results = results,
        log = fit_log,
        output_dir = condition_dir,
        replication = replication
      )

      rm(
        simulated,
        logged,
        results,
        fit_log,
        geo_log
      )
    }

    rm(reference, logged_geo)
  }
}

message("Finished fitting. Run scripts/summarize_results.R for summaries.")