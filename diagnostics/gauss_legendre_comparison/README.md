# Gauss–Legendre, direct_sum and pure spectral comparison

Residuals are vector L2 differences divided by the reference vector L2 norm over identical observation grids. The reference is higher-order composite Gauss–Legendre integration, checked against a second resolution; it is a numerical reference, not an exact solution. Independent kernel, geometry and unit checks live in `tests/test_gauss_legendre.py`.

Pure spectral runs set `auto_mode=0`, `edge_correction=0`, `hybrid_mode=0` and `reg_lambda=0`. Their source mesh uses horizontal refine factor 4 and 12 radial samples. Direct summation runs retain card radial count 4. Each timing is one warmed-build Python call including subprocess startup and file I/O; it is not a repeated benchmark.

The cases span magnetic surface/orbital/far fields, a concave magnetic polygon and gravity. Low harmonic degrees cannot resolve the narrow near-surface sources. Spectral gravity includes the degree-zero mass term from the signed sum of discretized source masses, with degrees 1 through lmax fitted to the remaining field. The direct_sum magnetic solver orients side-wall normals using polygon winding, including at concave edges. The notched footprint is independently checked against a sum of nonoverlapping rectangles in the tests.

Compiler: GNU Fortran (Homebrew GCC 15.2.0_1) 15.2.0. Platform: macOS-26.6.2-arm64-arm-64bit-Mach-O.

## Reference convergence

| Case | Unit | Peak magnitude | Reference change | Converged |
|---|---|---:|---:|---|
| magnetic_block_surface | nt | 1109.56 | 3.337e-14 | True |
| magnetic_block_orbital | nt | 117.27 | 1.039e-15 | True |
| magnetic_concave_polygon_orbital | nt | 88.8105 | 5.702e-16 | True |
| gravity_sloping_polygon_orbital | mgal | 15.2554 | 8.935e-16 | True |
| magnetic_block_far | nt | 0.000103364 | 6.335e-16 | True |
| gravity_block_far | mgal | 0.00232062 | 6.253e-16 | True |

## Resolution sweeps

| Case | Solver | Setting | Relative L2 | Vector RMSE | Time (s) |
|---|---|---|---:|---:|---:|
| magnetic_block_surface | gauss_legendre | order_4 | 1.544e+00 | 8.060e+02 | 0.0149 |
| magnetic_block_surface | gauss_legendre | order_8 | 2.426e+00 | 1.267e+03 | 0.0143 |
| magnetic_block_surface | gauss_legendre | order_16 | 4.106e-01 | 2.144e+02 | 0.0144 |
| magnetic_block_surface | gauss_legendre | order_32 | 4.996e-02 | 2.609e+01 | 0.0172 |
| magnetic_block_surface | gauss_legendre | order_64 | 6.499e-04 | 3.394e-01 | 0.0447 |
| magnetic_block_surface | gauss_legendre | order_8_panels_4 | 1.406e-03 | 7.343e-01 | 0.0176 |
| magnetic_block_surface | direct_sum | refine_1 | 1.742e-01 | 9.096e+01 | 0.0139 |
| magnetic_block_surface | direct_sum | refine_2 | 2.388e-02 | 1.247e+01 | 0.0144 |
| magnetic_block_surface | direct_sum | refine_4 | 7.489e-03 | 3.911e+00 | 0.0144 |
| magnetic_block_surface | spectral | degree_4 | 1.000e+00 | 5.223e+02 | 0.1291 |
| magnetic_block_surface | spectral | degree_8 | 1.000e+00 | 5.223e+02 | 0.1330 |
| magnetic_block_surface | spectral | degree_12 | 1.000e+00 | 5.223e+02 | 0.1638 |
| magnetic_block_surface | spectral | degree_24 | 9.940e-01 | 5.191e+02 | 1.4099 |
| magnetic_block_orbital | gauss_legendre | order_4 | 4.924e-03 | 3.429e-01 | 0.0146 |
| magnetic_block_orbital | gauss_legendre | order_8 | 4.465e-06 | 3.110e-04 | 0.0137 |
| magnetic_block_orbital | gauss_legendre | order_16 | 1.763e-12 | 1.228e-10 | 0.0138 |
| magnetic_block_orbital | gauss_legendre | order_32 | 1.445e-15 | 1.007e-13 | 0.0174 |
| magnetic_block_orbital | gauss_legendre | order_64 | 8.236e-16 | 5.736e-14 | 0.0443 |
| magnetic_block_orbital | gauss_legendre | order_8_panels_4 | 8.595e-15 | 5.987e-13 | 0.0179 |
| magnetic_block_orbital | direct_sum | refine_1 | 1.233e-03 | 8.586e-02 | 0.0136 |
| magnetic_block_orbital | direct_sum | refine_2 | 3.153e-04 | 2.196e-02 | 0.0135 |
| magnetic_block_orbital | direct_sum | refine_4 | 1.186e-04 | 8.260e-03 | 0.0140 |
| magnetic_block_orbital | spectral | degree_4 | 1.001e+00 | 6.975e+01 | 0.1254 |
| magnetic_block_orbital | spectral | degree_8 | 1.004e+00 | 6.991e+01 | 0.1292 |
| magnetic_block_orbital | spectral | degree_12 | 1.002e+00 | 6.982e+01 | 0.1627 |
| magnetic_block_orbital | spectral | degree_24 | 8.097e-01 | 5.640e+01 | 1.4164 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_4 | 3.383e-03 | 1.711e-01 | 0.0159 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_8 | 2.734e-06 | 1.383e-04 | 0.0148 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_16 | 1.247e-12 | 6.307e-11 | 0.0155 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_32 | 1.531e-15 | 7.746e-14 | 0.0256 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_64 | 7.480e-16 | 3.784e-14 | 0.1072 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_8_panels_4 | 5.510e-15 | 2.787e-13 | 0.0265 |
| magnetic_concave_polygon_orbital | direct_sum | refine_1 | 1.375e-03 | 6.953e-02 | 0.0145 |
| magnetic_concave_polygon_orbital | direct_sum | refine_2 | 3.366e-04 | 1.703e-02 | 0.0146 |
| magnetic_concave_polygon_orbital | direct_sum | refine_4 | 1.321e-04 | 6.683e-03 | 0.0150 |
| magnetic_concave_polygon_orbital | spectral | degree_4 | 1.001e+00 | 5.066e+01 | 0.0990 |
| magnetic_concave_polygon_orbital | spectral | degree_8 | 1.004e+00 | 5.079e+01 | 0.1023 |
| magnetic_concave_polygon_orbital | spectral | degree_12 | 1.003e+00 | 5.071e+01 | 0.1285 |
| magnetic_concave_polygon_orbital | spectral | degree_24 | 1.229e+00 | 6.215e+01 | 1.2944 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_4 | 1.618e-04 | 1.499e-03 | 0.0145 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_8 | 1.816e-08 | 1.683e-07 | 0.0140 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_16 | 2.486e-15 | 2.304e-14 | 0.0145 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_32 | 1.804e-15 | 1.672e-14 | 0.0205 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_64 | 1.095e-15 | 1.015e-14 | 0.0623 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_8_panels_4 | 2.808e-15 | 2.603e-14 | 0.0201 |
| gravity_sloping_polygon_orbital | direct_sum | refine_1 | 3.237e-03 | 3.000e-02 | 0.0137 |
| gravity_sloping_polygon_orbital | direct_sum | refine_2 | 8.207e-04 | 7.606e-03 | 0.0137 |
| gravity_sloping_polygon_orbital | direct_sum | refine_4 | 2.629e-04 | 2.436e-03 | 0.0148 |
| gravity_sloping_polygon_orbital | spectral | degree_4 | 9.929e-01 | 9.202e+00 | 0.0489 |
| gravity_sloping_polygon_orbital | spectral | degree_8 | 9.800e-01 | 9.083e+00 | 0.0510 |
| gravity_sloping_polygon_orbital | spectral | degree_12 | 9.661e-01 | 8.954e+00 | 0.0700 |
| gravity_sloping_polygon_orbital | spectral | degree_24 | 8.893e-01 | 8.242e+00 | 1.0767 |
| magnetic_block_far | gauss_legendre | order_4 | 2.024e-15 | 2.046e-19 | 0.0149 |
| magnetic_block_far | gauss_legendre | order_8 | 2.445e-15 | 2.472e-19 | 0.0139 |
| magnetic_block_far | gauss_legendre | order_16 | 1.408e-15 | 1.424e-19 | 0.0145 |
| magnetic_block_far | gauss_legendre | order_32 | 4.566e-17 | 4.617e-21 | 0.0179 |
| magnetic_block_far | gauss_legendre | order_64 | 6.903e-16 | 6.980e-20 | 0.0451 |
| magnetic_block_far | gauss_legendre | order_8_panels_4 | 2.516e-15 | 2.544e-19 | 0.0182 |
| magnetic_block_far | direct_sum | refine_1 | 1.019e-03 | 1.031e-07 | 0.0141 |
| magnetic_block_far | direct_sum | refine_2 | 8.847e-04 | 8.947e-08 | 0.0138 |
| magnetic_block_far | direct_sum | refine_4 | 8.231e-04 | 8.323e-08 | 0.0144 |
| magnetic_block_far | spectral | degree_4 | 4.137e-02 | 4.184e-06 | 0.1274 |
| magnetic_block_far | spectral | degree_8 | 5.076e-04 | 5.133e-08 | 0.1300 |
| magnetic_block_far | spectral | degree_12 | 4.312e-06 | 4.361e-10 | 0.1628 |
| magnetic_block_far | spectral | degree_24 | 9.844e-07 | 9.955e-11 | 1.4227 |
| gravity_block_far | gauss_legendre | order_4 | 2.079e-15 | 4.823e-18 | 0.0144 |
| gravity_block_far | gauss_legendre | order_8 | 2.459e-15 | 5.705e-18 | 0.0148 |
| gravity_block_far | gauss_legendre | order_16 | 1.407e-15 | 3.264e-18 | 0.0148 |
| gravity_block_far | gauss_legendre | order_32 | 8.292e-17 | 1.924e-19 | 0.0178 |
| gravity_block_far | gauss_legendre | order_64 | 6.953e-16 | 1.613e-18 | 0.0386 |
| gravity_block_far | gauss_legendre | order_8_panels_4 | 2.518e-15 | 5.841e-18 | 0.0174 |
| gravity_block_far | direct_sum | refine_1 | 7.437e-07 | 1.725e-09 | 0.0143 |
| gravity_block_far | direct_sum | refine_2 | 7.042e-07 | 1.634e-09 | 0.0147 |
| gravity_block_far | direct_sum | refine_4 | 7.656e-07 | 1.776e-09 | 0.0156 |
| gravity_block_far | spectral | degree_4 | 5.276e-03 | 1.224e-05 | 0.0920 |
| gravity_block_far | spectral | degree_8 | 3.727e-05 | 8.648e-08 | 0.0944 |
| gravity_block_far | spectral | degree_12 | 2.033e-07 | 4.716e-10 | 0.1201 |
| gravity_block_far | spectral | degree_24 | 7.128e-08 | 1.654e-10 | 1.2872 |

![Convergence and cost](convergence.png)

Exact models, options, per-component RMSE and solver diagnostics are saved in `results.json`; summed XYZ fields are saved in `fields.npz`. Residual units follow each case. `results.csv` provides the compact sweep table. Plot labels show the first and final node orders, direct_sum refine factors or spectral degrees. The separate diamond labelled 8×4 uses order 8 with four panels in each dimension.
