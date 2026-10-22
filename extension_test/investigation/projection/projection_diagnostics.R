# Saved truth and mesh diagnostic only; no model objects are changed
# Run from extension_test, or pass its path as the first argument
# An optional second argument selects a separate output directory
args <- commandArgs(trailingOnly = TRUE)
test_root <- normalizePath(if (length(args)) args[1] else getwd(), mustWork = TRUE)
stopifnot(file.exists(file.path(test_root, 'test_helpers.R')))
setwd(test_root)
source('test_helpers.R')
engine$configure_runtime(configs[[1]])
out <- if (length(args) >= 2L) args[2] else file.path(test_root, 'investigation', 'projection')
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# A linear function integrates exactly as area times its value at the centroid
# Clip each mesh triangle to each areal polygon, then integrate barycentric bases
exact_fe_projector <- function(mesh, areas) {
  tv <- mesh$graph$tv
  triangles <- sf::st_sf(triangle_id = seq_len(nrow(tv)), geometry = sf::st_sfc(
    lapply(seq_len(nrow(tv)), function(i) {
      p <- mesh$loc[tv[i, ], 1:2, drop = FALSE]
      sf::st_polygon(list(rbind(p, p[1, ])))
    }), crs = sf::st_crs(areas)))
  target <- sf::st_sf(area_id = seq_len(nrow(areas)), geometry = sf::st_geometry(areas))
  fragments <- suppressWarnings(sf::st_intersection(triangles, target))
  fragment_area <- as.numeric(sf::st_area(fragments))
  keep <- fragment_area > 1e-15
  fragments <- fragments[keep, ]
  fragment_area <- fragment_area[keep]
  centroid <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(fragments))))
  area <- as.numeric(sf::st_area(areas))
  barycentric <- t(vapply(seq_len(nrow(fragments)), function(i) {
    vertex <- tv[fragments$triangle_id[i], ]
    triangle <- mesh$loc[vertex, 1:2, drop = FALSE]
    solve(rbind(1, t(triangle)), c(1, centroid[i, ]))
  }, numeric(3)))
  if (any(barycentric < -1e-9) || any(barycentric > 1 + 1e-9)) {
    stop('Invalid clipped-fragment barycentric weights')
  }
  weights <- barycentric * (fragment_area / area[fragments$area_id])
  A <- Matrix::sparseMatrix(i = rep(fragments$area_id, each = 3),
    j = as.vector(t(tv[fragments$triangle_id, , drop = FALSE])),
    x = as.vector(t(weights)), dims = c(nrow(areas), mesh$n))
  list(A = A, fragments = fragments, area = area,
    covered_area = as.numeric(rowsum(fragment_area, fragments$area_id)))
}

# Values exist at cell centres only; these are diagnostic reconstructions
# All exterior coordinates are clamped to the nearest available centre
nodal_truth <- function(r, xy, method) {
  values <- matrix(as.numeric(raster::values(r)), nrow = nrow(r), ncol = ncol(r), byrow = TRUE)
  x <- raster::xFromCol(r, seq_len(ncol(r)))
  y <- raster::yFromRow(r, seq_len(nrow(r)))
  u <- pmin(pmax((xy[, 1] - min(x)) / raster::res(r)[1], 0), ncol(r) - 1)
  v <- pmin(pmax((max(y) - xy[, 2]) / raster::res(r)[2], 0), nrow(r) - 1)
  if (method == 'nearest_cell') return(values[cbind(floor(v + .5) + 1L, floor(u + .5) + 1L)])
  col <- pmin(floor(u) + 1L, ncol(r) - 1L)
  row <- pmin(floor(v) + 1L, nrow(r) - 1L)
  dx <- u - (col - 1L)
  dy <- v - (row - 1L)
  (1 - dx) * (1 - dy) * values[cbind(row, col)] +
    dx * (1 - dy) * values[cbind(row, col + 1L)] +
    (1 - dx) * dy * values[cbind(row + 1L, col)] +
    dx * dy * values[cbind(row + 1L, col + 1L)]
}

id <- 'strength_case_0003_Melding'
data <- readRDS(job_files(read_job('pilot', id))$data)
areas <- data$areas
centres <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(areas))))
raster_points <- sf::st_as_sf(as.data.frame(data$truth$coordinates), coords = c('x', 'y'))
cells <- sf::st_intersects(areas, raster_points)
true_means <- vapply(cells, function(index) mean(data$truth$Z[index]), numeric(1))
stopifnot(max(abs(true_means - areas$truth)) < 1e-12, all(lengths(cells) == 100L))

all_rows <- list(); meshes <- list(); audit <- list(); bundle_records <- list()
for (phase in c('pilot', 'finer_mesh')) {
  paths <- job_files(read_job(phase, id))
  bundle <- readRDS(paths$fit)
  check_data <- readRDS(paths$data)
  stopifnot(identical(data$truth$Z, check_data$truth$Z), identical(areas$truth, check_data$areas$truth))
  mesh <- bundle$spde$mesh
  equal <- bundle$projection_matrices$Ae2
  fe <- exact_fe_projector(mesh, areas)
  xy <- mesh$loc[, 1:2]
  area_centres_exact <- cbind(as.numeric(fe$A %*% xy[, 1]), as.numeric(fe$A %*% xy[, 2]))
  stopifnot(max(abs(Matrix::rowSums(fe$A) - 1)) < 1e-10,
    max(abs(area_centres_exact - centres)) < 1e-10,
    max(abs(fe$covered_area - fe$area)) < 1e-10)
  member_count <- Matrix::colSums(equal > 0)
  active_fe <- Matrix::colSums(fe$A) > 0
  outside <- xy[, 1] < -1e-10 | xy[, 1] > 1 + 1e-10 | xy[, 2] < -1e-10 | xy[, 2] > 1 + 1e-10
  audit[[phase]] <- data.frame(phase = phase, vertices = mesh$n, triangles = nrow(mesh$graph$tv),
    cells_per_area = min(lengths(cells)), max_abs_truth_reconstruction_error = max(abs(true_means - areas$truth)),
    equal_min_vertex_count = min(Matrix::rowSums(equal > 0)), equal_max_vertex_count = max(Matrix::rowSums(equal > 0)),
    equal_max_row_sum_error = max(abs(Matrix::rowSums(equal) - 1)),
    shared_boundary_vertices = sum(member_count > 1), shared_boundary_memberships = sum(pmax(member_count - 1, 0)),
    exact_max_row_sum_error = max(abs(Matrix::rowSums(fe$A) - 1)),
    exact_max_linear_function_error = max(abs(area_centres_exact - centres)),
    exact_max_coverage_error = max(abs(fe$covered_area - fe$area)),
    exact_nonzero_count = length(fe$A@x), outside_vertices_with_exact_weight = sum(active_fe & outside),
    outside_total_exact_weight = sum(fe$A[, outside, drop = FALSE]))
  for (method in c('nearest_cell', 'bilinear')) {
    z <- nodal_truth(data$field$raster, xy, method)
    equal_average <- as.numeric(equal %*% z)
    fe_average <- as.numeric(fe$A %*% z)
    all_rows[[paste(phase, method)]] <- data.frame(phase = phase, interpolation = method,
      area_id = seq_len(nrow(areas)), x = centres[, 1], y = centres[, 2], area = fe$area,
      vertex_count = Matrix::rowSums(equal > 0), true_average = true_means,
      equal_vertex_average = equal_average, exact_fe_average = fe_average,
      equal_minus_truth = equal_average - true_means, fe_minus_truth = fe_average - true_means,
      equal_minus_fe = equal_average - fe_average, observed_noise = areas$value - true_means)
  }
  meshes[[phase]] <- list(mesh = mesh, equal_A = equal, exact_fe_A = fe$A,
    area_coverage = fe$covered_area, fragments = fe$fragments)
  bundle_records[[phase]] <- list(file = paths$fit, data = paths$data)
  rm(bundle); gc()
}
rows <- do.call(rbind, all_rows)
summary <- do.call(rbind, lapply(split(rows, list(rows$phase, rows$interpolation), drop = TRUE), function(r) {
  do.call(rbind, lapply(c('equal_minus_truth', 'fe_minus_truth', 'equal_minus_fe'), function(error) {
    value <- r[[error]]
    data.frame(phase = r$phase[1], interpolation = r$interpolation[1], contrast = error,
      mean_error = mean(value), RMSE = sqrt(mean(value^2)), MAE = mean(abs(value)),
      max_abs_error = max(abs(value)), mean_squared_error = mean(value^2),
      rms_relative_to_true_areal_noise_sd = sqrt(mean(value^2)) / sqrt(.1))
  }))
}))
mesh_changes <- do.call(rbind, lapply(c('nearest_cell', 'bilinear'), function(method) {
  a <- rows[rows$phase == 'pilot' & rows$interpolation == method, ]
  b <- rows[rows$phase == 'finer_mesh' & rows$interpolation == method, ]
  stopifnot(identical(a$area_id, b$area_id))
  data.frame(interpolation = method, area_id = a$area_id, x = a$x, y = a$y,
    equal_finer_minus_base = b$equal_vertex_average - a$equal_vertex_average,
    exact_fe_finer_minus_base = b$exact_fe_average - a$exact_fe_average,
    quadrature_finer_minus_base = b$equal_minus_fe - a$equal_minus_fe)
}))
shared <- list()
for (phase in c('pilot', 'finer_mesh')) for (job in c('strength_case_0003_PSmelding', 'clustering_case_0003_Melding', 'clustering_case_0003_PSmelding')) {
  paths <- job_files(read_job(phase, job))
  d <- readRDS(paths$data)
  b <- readRDS(paths$fit)
  shared[[paste(phase, job)]] <- data.frame(phase = phase, job_id = job,
    same_truth_surface = identical(data$truth$Z, d$truth$Z), same_true_areal_means = identical(areas$truth, d$areas$truth),
    same_observed_areal_values = identical(areas$value, d$areas$value),
    same_mesh_coordinates = identical(meshes[[phase]]$mesh$loc, b$spde$mesh$loc),
    same_equal_vertex_projector = identical(meshes[[phase]]$equal_A, b$projection_matrices$Ae2))
  rm(b, d); gc()
}
write.csv(rows, file.path(out, 'projection_area_comparison.csv'), row.names = FALSE)
write.csv(summary, file.path(out, 'projection_summary.csv'), row.names = FALSE)
write.csv(do.call(rbind, audit), file.path(out, 'projection_audit.csv'), row.names = FALSE)
write.csv(mesh_changes, file.path(out, 'projection_mesh_changes.csv'), row.names = FALSE)
write.csv(do.call(rbind, shared), file.path(out, 'projection_shared_truth_check.csv'), row.names = FALSE)
saveRDS(list(areas = areas, raster = data$field$raster, comparisons = rows, mesh_changes = mesh_changes,
  meshes = meshes, audits = audit, source = bundle_records), file.path(out, 'projection_diagnostics.rds'))

png(file.path(out, 'areal_average_scatter.png'), width = 1800, height = 1500, res = 180)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
for (method in c('nearest_cell', 'bilinear')) for (phase in c('pilot', 'finer_mesh')) {
  r <- rows[rows$phase == phase & rows$interpolation == method, ]
  limits <- range(r$true_average, r$equal_vertex_average, r$exact_fe_average)
  plot(r$true_average, r$equal_vertex_average, pch = 16, col = '#0072B2', xlim = limits, ylim = limits,
    xlab = 'Stored true areal average', ylab = 'Average of interpolated nodal truth',
    main = paste(phase, '/', method))
  points(r$true_average, r$exact_fe_average, pch = 1, col = '#D55E00')
  abline(0, 1, lty = 2, col = 'gray50')
  legend('topleft', c('Equal vertex A_B', 'Exact FE integration'), pch = c(16, 1),
    col = c('#0072B2', '#D55E00'), bty = 'n', cex = .8)
}
dev.off()
png(file.path(out, 'areal_projection_error_maps.png'), width = 2100, height = 1500, res = 180)
par(mfrow = c(2, 3), mar = c(2.5, 2.5, 3, 1), oma = c(0, 0, 2, 0))
r <- rows[rows$interpolation == 'bilinear', ]
limit <- max(abs(unlist(r[c('equal_minus_truth', 'fe_minus_truth', 'equal_minus_fe')])))
breaks <- seq(-limit, limit, length.out = 102)
colours <- grDevices::colorRampPalette(c('#2166AC', 'white', '#B2182B'))(101)
for (phase in c('pilot', 'finer_mesh')) for (error in c('equal_minus_truth', 'fe_minus_truth', 'equal_minus_fe')) {
  current <- r[r$phase == phase, ]
  plot(sf::st_geometry(areas), col = colours[findInterval(current[[error]], breaks, all.inside = TRUE)],
    border = 'grey70', main = paste(phase, error, sep = '\n'), axes = TRUE)
  text(current$x, current$y, sprintf('%.3f', current[[error]]), cex = .65)
}
mtext(sprintf('Bilinear reconstruction; common colour scale +/- %.3f', limit), outer = TRUE, line = 0)
dev.off()
writeLines(c(
  'AREAL PROJECTION DIAGNOSTIC',
  'Purpose: compare the saved equal-vertex A_B with exact integration of the same piecewise-linear FE field',
  'The stored generating field is available only at 2500 cell centres',
  'Nearest-cell and bilinear reconstructions define diagnostic nodal values; neither recovers unobserved continuous-field truth',
  'Both reconstructions clamp coordinates to the nearest available cell centre at/outside the boundary',
  'Exact FE integration clips every mesh triangle to every overlapping area, integrates each linear basis by its clipped-polygon centroid and area, and divides by block area',
  'Each exact row sums to one, each block is covered, and linear x/y coordinate functions integrate to the polygon centroid',
  'All 25 original block averages match the mean of their 100 stored raster cells exactly within numerical tolerance',
  'equal_minus_fe isolates the equal-vertex quadrature difference conditional on the selected nodal values',
  'fe_minus_truth contains FE representation, raster interpolation and boundary extension effects',
  'equal_minus_truth combines both effects, including their interaction/cancellation',
  'This diagnostic does not refit models and cannot causally attribute the areal-noise variance jump to projection or optimizer instability',
  'Read projection_summary.csv, projection_audit.csv, projection_shared_truth_check.csv, and the pointwise projection_area_comparison.csv',
  'No saved fit, original code, or original output was modified'), file.path(out, 'README.txt'))
print(summary)
print(do.call(rbind, audit))
print(do.call(rbind, shared))
