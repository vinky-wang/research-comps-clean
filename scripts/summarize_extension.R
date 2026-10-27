
source("functions/generate_data.R")
source("functions/create_summary.R")


# Example  ---------------------------------------------------------------------

# Preferential sampling scenarios
scenarios <- c(
  "Constant gamma" = "constant",
  "Spatially varying gamma" = "spatially_varying",
  "Independent latent field" = "independent_field"
)

# Generate data for each scenario
simulations <- lapply(scenarios, function(scenario) {
  generate_data(
    rho = 0.2,
    gamma = 1,
    I = 100,
    J = 25,
    tau2_s = 0.1,
    tau2_B = 0.1,
    seed = 123,
    scenario = scenario,
    delta = 1,
    beta_psi = 1,
    rho_psi = 0.2
  )
})

# Extract simulation outputs
field <- simulations[[1]]$latent
log_lambda <- lapply(simulations, function(d) d$latent$log_lambda)
points <- lapply(simulations, function(d) d$points)

# Select one simulated dataset
dataset <- simulations[["Constant gamma"]]
truth <- dataset$prediction$true


# Plot of extended scenarios
plot_extended <- plot_simulated_data(dataset)

plot_extended <- plot_sampling_scenarios(
  field = field,
  scenarios = log_lambda,
  points = points,
  limits = c(-5, 5)
)




# Compare mesh, areal units, and points across scenarios -----------------------
mesh <- create_mesh()

# Extract triangle edges from the INLA mesh
triangles <- mesh$graph$tv
edge_ids <- rbind(
  triangles[, c(1, 2)],
  triangles[, c(2, 3)],
  triangles[, c(3, 1)]
)

mesh_edges <- data.frame(
  x = mesh$loc[edge_ids[, 1], 1],
  y = mesh$loc[edge_ids[, 1], 2],
  xend = mesh$loc[edge_ids[, 2], 1],
  yend = mesh$loc[edge_ids[, 2], 2]
)

# Combine points, keeping the scenario names for faceting
scenario_points <- do.call(rbind, lapply(names(simulations), function(scenario) {
  data.frame(
    x = simulations[[scenario]]$points$x,
    y = simulations[[scenario]]$points$y,
    scenario = scenario
  )
}))
scenario_points$scenario <- factor(
  scenario_points$scenario,
  levels = names(simulations)
)

plot_mesh_areas_points <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = mesh_edges,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey80", linewidth = 0.2
  ) +
  ggplot2::geom_rect(
    data = dataset$areas,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA, colour = "black", linewidth = 0.7
  ) +
  ggplot2::geom_point(
    data = scenario_points,
    ggplot2::aes(x = x, y = y),
    colour = "red", size = 1.2
  ) +
  ggplot2::facet_wrap(~ scenario, nrow = 1) +
  ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  ggplot2::labs(x = "", y = "") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      colour = "black", fill = NA, linewidth = 0.7
    )
  )

print(plot_mesh_areas_points)
