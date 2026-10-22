AREAL PROJECTION DIAGNOSTIC
Purpose: compare the saved equal-vertex A_B with exact integration of the same piecewise-linear FE field
The stored generating field is available only at 2500 cell centres
Nearest-cell and bilinear reconstructions define diagnostic nodal values; neither recovers unobserved continuous-field truth
Both reconstructions clamp coordinates to the nearest available cell centre at/outside the boundary
Exact FE integration clips every mesh triangle to every overlapping area, integrates each linear basis by its clipped-polygon centroid and area, and divides by block area
Each exact row sums to one, each block is covered, and linear x/y coordinate functions integrate to the polygon centroid
All 25 original block averages match the mean of their 100 stored raster cells exactly within numerical tolerance
equal_minus_fe isolates the equal-vertex quadrature difference conditional on the selected nodal values
fe_minus_truth contains FE representation, raster interpolation and boundary extension effects
equal_minus_truth combines both effects, including their interaction/cancellation
This diagnostic does not refit models and cannot causally attribute the areal-noise variance jump to projection or optimizer instability
Read projection_summary.csv, projection_audit.csv, projection_shared_truth_check.csv, and the pointwise projection_area_comparison.csv
No saved fit, original code, or original output was modified
