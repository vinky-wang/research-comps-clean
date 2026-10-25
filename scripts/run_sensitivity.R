source("functions/generate_data.R")
source("functions/create_summary.R")
source("functions/fit_models.R")


# Generate one simulated dataset ----------------------------------------------

dataset <- generate_data(
  rho = 0.2,
  gamma = 1,
  I = 100,
  J = 25,
  tau2_s = 0.1,
  tau2_B = 0.1,
  seed = 123,
  scenario = "constant"
)

truth <- dataset$prediction$true


# Check mesh and projection matrices ------------------------------------------

# Mesh sensitivity -------------------------------------------------------------

# Mesh resolutions
mesh_settings <- list(
  coarse = c(0.10, 0.25),
  medium = c(0.05, 0.25),
  fine = c(0.025, 0.25)
)

# Create meshes
meshes <- lapply(mesh_settings, function(max_edge) {
  create_mesh(
    max_edge = max_edge
  )
})

# Compare number of mesh vertices
sapply(meshes, function(mesh) mesh$n)

# Plot meshes
old_par <- par(mfrow = c(1, 3))

for (name in names(meshes)) {
  plot(meshes[[name]], main = name)
}

par(old_par)

# Select medium mesh for subsequent comparisons
mesh <- meshes$medium
plot(mesh)

# Create SPDE model
spde <- create_spde(mesh)

# Create projection matrices
A <- create_projections(
  data = dataset,
  mesh = mesh,
  areal_method = "vertices"
)

# Check dimensions
dim(A$points)
dim(A$areas)
dim(A$prediction)

# Check areal averaging weights
range(Matrix::rowSums(A$areas))


# Check likelihood stacks ------------------------------------------------------

# Gaussian observation stacks
gaussian_shared <- create_gaussian_stack(
  data = dataset,
  projections = A,
  spde = spde,
  include_areas = TRUE,
  preferential = TRUE,
  separate_noise = FALSE
)

gaussian_separate <- create_gaussian_stack(
  data = dataset,
  projections = A,
  spde = spde,
  include_areas = TRUE,
  preferential = TRUE,
  separate_noise = TRUE
)

# Mesh-dual integration weights
weights <- create_mesh_dual_weights(mesh)

c(
  total_weight = sum(weights),
  positive_weights = sum(weights > 0),
  mesh_vertices = mesh$n
)

# Check integration of constant and linear surfaces
c(
  integral_1 = sum(weights),
  integral_x = sum(weights * mesh$loc[, 1]),
  integral_y = sum(weights * mesh$loc[, 2])
)

# Preferential sampling stacks
sampling_grid <- create_point_process_stack(
  data = dataset,
  projections = A,
  spde = spde,
  mesh = mesh,
  point_process_method = "grid_counts"
)

sampling_dual <- create_point_process_stack(
  data = dataset,
  projections = A,
  spde = spde,
  mesh = mesh,
  point_process_method = "mesh_dual",
  integration_weights = weights
)

# Check likelihood dimensions
rbind(
  gaussian_shared = dim(INLA::inla.stack.data(gaussian_shared)$y),
  gaussian_separate = dim(INLA::inla.stack.data(gaussian_separate)$y),
  sampling_grid = dim(INLA::inla.stack.data(sampling_grid)$y),
  sampling_dual = dim(INLA::inla.stack.data(sampling_dual)$y)
)

# Check sampled points and integration weights
check_sampling_stack <- function(stack) {
  stack_data <- INLA::inla.stack.data(stack)
  poisson_col <- ncol(stack_data$y)
  
  return(c(
    points = sum(stack_data$y[, poisson_col], na.rm = TRUE),
    exposure = sum(stack_data$exposure)
  ))
}

rbind(
  grid_counts = check_sampling_stack(sampling_grid),
  mesh_dual = check_sampling_stack(sampling_dual)
)


# Fit all four models ----------------------------------------------------------

# Primary fitting specification:
# - Mesh-dual point-process likelihood
# - Mesh-vertex areal integration
# - Shared Gaussian observation precision
# - Constrained Matern SPDE

results_main <- fit_models(
  data = dataset,
  mesh = mesh,
  point_process_method = "mesh_dual",
  areal_method = "vertices",
  separate_noise = FALSE,
  prior_gamma = c(1, 10)
)

# Model summary
summary(results_main$models$ps_melding$fit)

# Posterior predictions
head(results_main$models$ps_melding$prediction)

# Compare prediction accuracy
model_metrics <- t(vapply(
  results_main$models,
  function(model) evaluate(model$prediction, truth),
  numeric(3)
))

# Baseline predictor: Z(s) = 0
rmse_baseline <- sqrt(mean(truth^2))

rbind(
  baseline = c(
    rmse = rmse_baseline,
    coverage = NA_real_,
    mean_sd = NA_real_
  ),
  model_metrics
)

# Posterior estimate of gamma
results_main$models$ps_melding$fit$
  summary.hyperpar["Beta for j", ]


# Prior sensitivity ------------------------------------------------------------

# Gamma prior: Normal(mean, precision)
# Default prior: mean = 1, precision = 10
# Compare changes in prior mean and precision separately

# Shift prior mean while holding precision fixed
results_gamma_mean <- fit_models(
  data = dataset,
  models = "ps_melding",
  mesh = mesh,
  point_process_method = "mesh_dual",
  areal_method = "vertices",
  separate_noise = FALSE,
  prior_gamma = c(0, 10)
)

# Reduce prior precision while holding mean fixed
results_gamma_precision <- fit_models(
  data = dataset,
  models = "ps_melding",
  mesh = mesh,
  point_process_method = "mesh_dual",
  areal_method = "vertices",
  separate_noise = FALSE,
  prior_gamma = c(1, 1)
)

# Shift prior mean and reduce precision
results_weak <- fit_models(
  data = dataset,
  models = "ps_melding",
  mesh = mesh,
  point_process_method = "mesh_dual",
  areal_method = "vertices",
  separate_noise = FALSE,
  prior_gamma = c(0, 1)
)

# Compare posterior estimates of gamma and predictions
summarize_ps(
  results = list(
    default = results_main,
    shifted_mean = results_gamma_mean,
    lower_precision = results_gamma_precision,
    weak = results_weak
  ),
  truth = truth
)


# Observation-noise sensitivity ------------------------------------------------

# Separate Gaussian precisions for point and areal observations
results_separate <- fit_models(
  data = dataset,
  models = "ps_melding",
  mesh = mesh,
  point_process_method = "mesh_dual",
  areal_method = "vertices",
  separate_noise = TRUE,
  prior_gamma = c(1, 10)
)

# Compare shared and separate Gaussian precisions
summarize_ps(
  results = list(
    shared = results_main,
    separate = results_separate
  ),
  truth = truth
)

# Examine estimated observation-error precisions
results_separate$models$ps_melding$fit$summary.hyperpar


# Point-process approximation sensitivity -------------------------------------

# Grid-count approximation
results_counts <- fit_models(
  data = dataset,
  models = "ps_melding",
  mesh = mesh,
  point_process_method = "grid_counts",
  areal_method = "vertices",
  separate_noise = FALSE,
  prior_gamma = c(1, 10)
)

# Compare point-process approximations
summarize_ps(
  results = list(
    grid_counts = results_counts,
    mesh_dual = results_main
  ),
  truth = truth
)

# Compare posterior mean predictions
pred_counts <- results_counts$models$ps_melding$prediction
pred_dual <- results_main$models$ps_melding$prediction

compare_predictions(pred_counts, pred_dual)


# Areal projection sensitivity -------------------------------------------------

# Grid-based areal projection
A_grid <- create_projections(
  data = dataset,
  mesh = mesh,
  areal_method = "grid"
)

# Check dimensions
dim(A_grid$areas)
dim(A$areas)

# Check averaging weights
rbind(
  grid = range(Matrix::rowSums(A_grid$areas)),
  mesh_vertices = range(Matrix::rowSums(A$areas))
)

# Check areal projection for a linear surface
area_centres <- data.frame(
  x = (dataset$areas$xmin + dataset$areas$xmax) / 2,
  y = (dataset$areas$ymin + dataset$areas$ymax) / 2
)

check_areal_projection <- function(A) {
  projected_x <- as.vector(A$areas %*% mesh$loc[, 1])
  projected_y <- as.vector(A$areas %*% mesh$loc[, 2])
  
  return(c(
    x = max(abs(projected_x - area_centres$x)),
    y = max(abs(projected_y - area_centres$y))
  ))
}

rbind(
  grid = check_areal_projection(A_grid),
  mesh_vertices = check_areal_projection(A)
)

# Fit using grid-based areal integration
results_grid <- fit_models(
  data = dataset,
  models = "ps_melding",
  mesh = mesh,
  point_process_method = "mesh_dual",
  areal_method = "grid",
  separate_noise = FALSE,
  prior_gamma = c(1, 10)
)

# Compare areal integration methods
summarize_ps(
  results = list(
    grid = results_grid,
    mesh_vertices = results_main
  ),
  truth = truth
)

# Compare posterior mean predictions
pred_grid <- results_grid$models$ps_melding$prediction
pred_vertices <- results_main$models$ps_melding$prediction

compare_predictions(pred_grid, pred_vertices)


# SPDE constraint sensitivity --------------------------------------------------

# Compare constr = FALSE and constr = TRUE
# Hold mesh, likelihoods, projections, and priors fixed

# Construct unconstrained Matern SPDE
spde_unconstrained <- INLA::inla.spde2.matern(
  mesh = mesh,
  alpha = 2,
  constr = FALSE
)

# Fit melding and PS melding with unconstrained SPDE
constraint_models <- c("melding", "ps_melding")

fits_unconstrained <- lapply(constraint_models, function(model) {
  fit_model(
    data = dataset,
    model = model,
    spde = spde_unconstrained,
    projections = A,
    mesh = mesh,
    prior_gamma = c(1, 10),
    separate_noise = FALSE,
    point_process_method = "mesh_dual",
    integration_weights = weights
  )
})

names(fits_unconstrained) <- constraint_models

results_unconstrained <- list(
  models = fits_unconstrained,
  mesh = mesh,
  spde = spde_unconstrained
)

# Compare prediction accuracy for both models
rbind(
  melding_unconstrained = evaluate(
    results_unconstrained$models$melding$prediction, truth
  ),
  melding_constrained = evaluate(
    results_main$models$melding$prediction, truth
  ),
  ps_melding_unconstrained = evaluate(
    results_unconstrained$models$ps_melding$prediction, truth
  ),
  ps_melding_constrained = evaluate(
    results_main$models$ps_melding$prediction, truth
  )
)

# Compare posterior estimates of gamma
summarize_ps(
  results = list(
    unconstrained = results_unconstrained,
    constrained = results_main
  ),
  truth = truth
)

# Compare posterior mean predictions
pred_unconstrained <- results_unconstrained$models$ps_melding$prediction
pred_constrained <- results_main$models$ps_melding$prediction

compare_predictions(pred_unconstrained, pred_constrained)

# Compare fixed effects
rbind(
  unconstrained = results_unconstrained$models$ps_melding$fit$
    summary.fixed["b0", ],
  constrained = results_main$models$ps_melding$fit$
    summary.fixed["b0", ]
)

# Compare estimated hyperparameters
results_unconstrained$models$ps_melding$fit$summary.hyperpar
results_main$models$ps_melding$fit$summary.hyperpar


# Spatial prior sensitivity ----------------------------------------------------

# Compare default Matern priors with penalized complexity priors
# Hold mesh, constraint, and likelihood specifications fixed

# PC prior:
# P(range < 0.2) = 0.5
# P(sigma > 1) = 0.5

results_pc <- fit_models(
  data = dataset,
  models = "ps_melding",
  mesh = mesh,
  prior_range = c(0.2, 0.5),
  prior_sigma = c(1, 0.5),
  point_process_method = "mesh_dual",
  areal_method = "vertices",
  separate_noise = FALSE,
  prior_gamma = c(1, 10)
)

# Compare default Matern and PC priors
summarize_ps(
  results = list(
    default = results_main,
    pc_prior = results_pc
  ),
  truth = truth
)

# Compare posterior mean predictions
pred_pc <- results_pc$models$ps_melding$prediction

compare_predictions(pred_constrained, pred_pc)


# Spatial prior and constraint interaction ------------------------------------

# Construct PC-prior SPDE without integrate-to-zero constraint
spde_pc_unconstrained <- INLA::inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(0.2, 0.5),
  prior.sigma = c(1, 0.5),
  constr = FALSE
)

# Fit PS melding with unconstrained PC-prior SPDE
results_pc_unconstrained <- list(
  models = list(
    ps_melding = fit_model(
      data = dataset,
      model = "ps_melding",
      spde = spde_pc_unconstrained,
      projections = A,
      mesh = mesh,
      prior_gamma = c(1, 10),
      separate_noise = FALSE,
      point_process_method = "mesh_dual",
      integration_weights = weights
    )
  )
)

# Compare all four spatial-prior specifications
summarize_ps(
  results = list(
    matern_unconstrained = results_unconstrained,
    matern_constrained = results_main,
    pc_unconstrained = results_pc_unconstrained,
    pc_constrained = results_pc
  ),
  truth = truth
)


# Fitted mesh sensitivity ------------------------------------------------------

# Compare PS melding under different mesh resolutions
# Hold point-process quadrature, areal integration, and priors fixed
# Use the same 50 x 50 grid-count approximation for all meshes

# Fit coarse mesh
results_coarse <- fit_models(
  data = dataset,
  models = "ps_melding",
  mesh = meshes$coarse,
  point_process_method = "grid_counts",
  areal_method = "grid",
  separate_noise = FALSE,
  prior_gamma = c(1, 10)
)

# Fit medium mesh
results_medium <- fit_models(
  data = dataset,
  models = "ps_melding",
  mesh = meshes$medium,
  point_process_method = "grid_counts",
  areal_method = "grid",
  separate_noise = FALSE,
  prior_gamma = c(1, 10)
)

# Fit fine mesh
results_fine <- fit_models(
  data = dataset,
  models = "ps_melding",
  mesh = meshes$fine,
  point_process_method = "grid_counts",
  areal_method = "grid",
  separate_noise = FALSE,
  prior_gamma = c(1, 10)
)

# Collect fitted models
results_mesh <- list(
  coarse = results_coarse,
  medium = results_medium,
  fine = results_fine
)

# Compare mesh vertices and posterior estimates
mesh_comparison <- summarize_ps(
  results = results_mesh,
  truth = truth
)

mesh_comparison <- data.frame(
  mesh_vertices = sapply(meshes, function(mesh) mesh$n),
  mesh_comparison
)

mesh_comparison

# Compare posterior mean predictions
pred_coarse <- results_mesh$coarse$models$ps_melding$prediction
pred_medium <- results_mesh$medium$models$ps_melding$prediction
pred_fine <- results_mesh$fine$models$ps_melding$prediction

rbind(
  coarse_vs_medium = compare_predictions(pred_coarse, pred_medium),
  medium_vs_fine = compare_predictions(pred_medium, pred_fine),
  coarse_vs_fine = compare_predictions(pred_coarse, pred_fine)
)


# INLA approximation sensitivity ----------------------------------------------

# Compare default INLA inference with variance correction
# Use the non-preferential melding model
# Hold mesh, SPDE, projections, and priors fixed

# Construct Gaussian stack without preferential sampling
melding_stack <- create_gaussian_stack(
  data = dataset,
  projections = A,
  spde = spde,
  include_areas = TRUE,
  preferential = FALSE,
  separate_noise = FALSE
)

# Fit melding with INLA variance correction
fit_melding_vb <- INLA::inla(
  formula = y ~ 0 + b0 + f(i, model = spde),
  family = "gaussian",
  data = INLA::inla.stack.data(melding_stack),
  control.fixed = list(
    mean = 0,
    prec = 0.001
  ),
  control.predictor = list(
    A = INLA::inla.stack.A(melding_stack),
    compute = TRUE,
    link = 1
  ),
  control.inla = list(
    control.vb = list(strategy = "variance")
  ),
  control.compute = list(
    dic = FALSE,
    waic = FALSE
  )
)

# Extract posterior predictions
index_vb <- INLA::inla.stack.index(
  melding_stack,
  tag = "prediction"
)$data

summary_vb <- fit_melding_vb$summary.linear.predictor[index_vb, ]

pred_vb <- data.frame(
  x = dataset$prediction$x,
  y = dataset$prediction$y,
  mean = summary_vb$mean,
  sd = summary_vb$sd,
  q025 = summary_vb[["0.025quant"]],
  q975 = summary_vb[["0.975quant"]]
)

# Compare prediction accuracy
pred_default <- results_main$models$melding$prediction

rbind(
  default = evaluate(pred_default, truth),
  variance_correction = evaluate(pred_vb, truth)
)

# Compare posterior mean predictions
compare_predictions(pred_default, pred_vb)

# Compare fixed effects
rbind(
  default = results_main$models$melding$fit$
    summary.fixed["b0", ],
  variance_correction = fit_melding_vb$
    summary.fixed["b0", ]
)

# Compare hyperparameters
results_main$models$melding$fit$summary.hyperpar
fit_melding_vb$summary.hyperpar
