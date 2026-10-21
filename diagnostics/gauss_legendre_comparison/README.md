# Gauss–Legendre, direct_sum and pure spectral comparison

Residuals are vector L2 differences divided by the reference vector L2 norm over identical observation grids. The reference is higher-order composite Gauss–Legendre integration, checked against a second resolution. it is a numerical reference, not an exact solution. Independent kernel, geometry and unit checks live in `tests/test_gauss_legendre.py`.

Pure spectral runs use the native Gauss–Legendre/Fourier transform with extended-range Schmidt semi-normalization without the Condon-Shortley phase and set `auto_mode=0`, `edge_correction=0`, `hybrid_mode=0` and `reg_lambda=0`. Fit latitudes are Gauss–Legendre nodes and longitude counts round upward to powers of two. Their source mesh uses horizontal refine factor 4 and 12 radial samples. Point-source runs retain card radial count 4. Each timing is one warmed-build Python call including subprocess startup and file I/O. It is not a repeated benchmark.

The cases span magnetic surface/orbital/far fields, a concave magnetic polygon and gravity. Low harmonic degrees cannot resolve the narrow near-surface sources. Spectral gravity includes the degree-zero mass term from the signed sum of discretized source masses, with degrees 1 through lmax fitted to the remaining field. The direct_sum magnetic solver orients side-wall normals using polygon winding, including at concave edges. The notched footprint is independently checked against a sum of nonoverlapping rectangles in the tests.

Compiler: GNU Fortran (Homebrew GCC 16.2.0) 16.2.0. Platform: macOS-26.6.2-arm64-arm-64bit-Mach-O.

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
| magnetic_block_surface | gauss_legendre | order_4 | 1.544e+00 | 8.060e+02 | 0.0136 |
| magnetic_block_surface | gauss_legendre | order_8 | 2.426e+00 | 1.267e+03 | 0.0130 |
| magnetic_block_surface | gauss_legendre | order_16 | 4.106e-01 | 2.144e+02 | 0.0130 |
| magnetic_block_surface | gauss_legendre | order_32 | 4.996e-02 | 2.609e+01 | 0.0162 |
| magnetic_block_surface | gauss_legendre | order_64 | 6.499e-04 | 3.394e-01 | 0.0432 |
| magnetic_block_surface | gauss_legendre | order_8_panels_4 | 1.406e-03 | 7.343e-01 | 0.0166 |
| magnetic_block_surface | direct_sum | refine_1 | 1.742e-01 | 9.096e+01 | 0.0129 |
| magnetic_block_surface | direct_sum | refine_2 | 2.388e-02 | 1.247e+01 | 0.0126 |
| magnetic_block_surface | direct_sum | refine_4 | 7.489e-03 | 3.911e+00 | 0.0141 |
| magnetic_block_surface | spectral | degree_4 | 1.000e+00 | 5.222e+02 | 0.5033 |
| magnetic_block_surface | spectral | degree_8 | 1.000e+00 | 5.223e+02 | 0.1642 |
| magnetic_block_surface | spectral | degree_12 | 1.000e+00 | 5.223e+02 | 0.1771 |
| magnetic_block_surface | spectral | degree_24 | 1.005e+00 | 5.248e+02 | 0.6383 |
| magnetic_block_orbital | gauss_legendre | order_4 | 4.924e-03 | 3.429e-01 | 0.0141 |
| magnetic_block_orbital | gauss_legendre | order_8 | 4.465e-06 | 3.110e-04 | 0.0130 |
| magnetic_block_orbital | gauss_legendre | order_16 | 1.763e-12 | 1.228e-10 | 0.0131 |
| magnetic_block_orbital | gauss_legendre | order_32 | 1.445e-15 | 1.007e-13 | 0.0165 |
| magnetic_block_orbital | gauss_legendre | order_64 | 8.236e-16 | 5.736e-14 | 0.0429 |
| magnetic_block_orbital | gauss_legendre | order_8_panels_4 | 8.595e-15 | 5.987e-13 | 0.0167 |
| magnetic_block_orbital | direct_sum | refine_1 | 1.233e-03 | 8.586e-02 | 0.0128 |
| magnetic_block_orbital | direct_sum | refine_2 | 3.153e-04 | 2.196e-02 | 0.0131 |
| magnetic_block_orbital | direct_sum | refine_4 | 1.186e-04 | 8.260e-03 | 0.0137 |
| magnetic_block_orbital | spectral | degree_4 | 1.001e+00 | 6.972e+01 | 0.1661 |
| magnetic_block_orbital | spectral | degree_8 | 1.002e+00 | 6.981e+01 | 0.1640 |
| magnetic_block_orbital | spectral | degree_12 | 1.003e+00 | 6.983e+01 | 0.1790 |
| magnetic_block_orbital | spectral | degree_24 | 9.854e-01 | 6.863e+01 | 0.6451 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_4 | 3.383e-03 | 1.711e-01 | 0.0141 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_8 | 2.734e-06 | 1.383e-04 | 0.0130 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_16 | 1.247e-12 | 6.307e-11 | 0.0143 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_32 | 1.531e-15 | 7.746e-14 | 0.0242 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_64 | 7.480e-16 | 3.784e-14 | 0.1056 |
| magnetic_concave_polygon_orbital | gauss_legendre | order_8_panels_4 | 5.510e-15 | 2.787e-13 | 0.0250 |
| magnetic_concave_polygon_orbital | direct_sum | refine_1 | 1.375e-03 | 6.953e-02 | 0.0130 |
| magnetic_concave_polygon_orbital | direct_sum | refine_2 | 3.366e-04 | 1.703e-02 | 0.0129 |
| magnetic_concave_polygon_orbital | direct_sum | refine_4 | 1.321e-04 | 6.683e-03 | 0.0135 |
| magnetic_concave_polygon_orbital | spectral | degree_4 | 1.001e+00 | 5.063e+01 | 0.1257 |
| magnetic_concave_polygon_orbital | spectral | degree_8 | 1.002e+00 | 5.070e+01 | 0.1273 |
| magnetic_concave_polygon_orbital | spectral | degree_12 | 1.003e+00 | 5.072e+01 | 0.1365 |
| magnetic_concave_polygon_orbital | spectral | degree_24 | 9.810e-01 | 4.962e+01 | 0.4849 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_4 | 1.618e-04 | 1.499e-03 | 0.0138 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_8 | 1.816e-08 | 1.683e-07 | 0.0126 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_16 | 2.486e-15 | 2.304e-14 | 0.0134 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_32 | 1.804e-15 | 1.672e-14 | 0.0186 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_64 | 1.095e-15 | 1.015e-14 | 0.0621 |
| gravity_sloping_polygon_orbital | gauss_legendre | order_8_panels_4 | 2.808e-15 | 2.603e-14 | 0.0188 |
| gravity_sloping_polygon_orbital | direct_sum | refine_1 | 3.237e-03 | 3.000e-02 | 0.0128 |
| gravity_sloping_polygon_orbital | direct_sum | refine_2 | 8.207e-04 | 7.606e-03 | 0.0131 |
| gravity_sloping_polygon_orbital | direct_sum | refine_4 | 2.629e-04 | 2.436e-03 | 0.0131 |
| gravity_sloping_polygon_orbital | spectral | degree_4 | 9.933e-01 | 9.207e+00 | 0.0571 |
| gravity_sloping_polygon_orbital | spectral | degree_8 | 9.816e-01 | 9.098e+00 | 0.0578 |
| gravity_sloping_polygon_orbital | spectral | degree_12 | 9.647e-01 | 8.941e+00 | 0.0618 |
| gravity_sloping_polygon_orbital | spectral | degree_24 | 8.571e-01 | 7.944e+00 | 0.2046 |
| magnetic_block_far | gauss_legendre | order_4 | 2.024e-15 | 2.046e-19 | 0.0138 |
| magnetic_block_far | gauss_legendre | order_8 | 2.445e-15 | 2.472e-19 | 0.0127 |
| magnetic_block_far | gauss_legendre | order_16 | 1.408e-15 | 1.424e-19 | 0.0138 |
| magnetic_block_far | gauss_legendre | order_32 | 4.566e-17 | 4.617e-21 | 0.0168 |
| magnetic_block_far | gauss_legendre | order_64 | 6.903e-16 | 6.980e-20 | 0.0431 |
| magnetic_block_far | gauss_legendre | order_8_panels_4 | 2.516e-15 | 2.544e-19 | 0.0177 |
| magnetic_block_far | direct_sum | refine_1 | 1.019e-03 | 1.031e-07 | 0.0129 |
| magnetic_block_far | direct_sum | refine_2 | 8.847e-04 | 8.947e-08 | 0.0136 |
| magnetic_block_far | direct_sum | refine_4 | 8.231e-04 | 8.323e-08 | 0.0135 |
| magnetic_block_far | spectral | degree_4 | 4.137e-02 | 4.183e-06 | 0.1650 |
| magnetic_block_far | spectral | degree_8 | 5.075e-04 | 5.132e-08 | 0.1658 |
| magnetic_block_far | spectral | degree_12 | 4.311e-06 | 4.359e-10 | 0.1802 |
| magnetic_block_far | spectral | degree_24 | 2.894e-08 | 2.927e-12 | 0.6415 |
| gravity_block_far | gauss_legendre | order_4 | 2.079e-15 | 4.823e-18 | 0.0139 |
| gravity_block_far | gauss_legendre | order_8 | 2.459e-15 | 5.705e-18 | 0.0128 |
| gravity_block_far | gauss_legendre | order_16 | 1.407e-15 | 3.264e-18 | 0.0133 |
| gravity_block_far | gauss_legendre | order_32 | 8.292e-17 | 1.924e-19 | 0.0155 |
| gravity_block_far | gauss_legendre | order_64 | 6.953e-16 | 1.613e-18 | 0.0373 |
| gravity_block_far | gauss_legendre | order_8_panels_4 | 2.518e-15 | 5.841e-18 | 0.0159 |
| gravity_block_far | direct_sum | refine_1 | 7.437e-07 | 1.725e-09 | 0.0129 |
| gravity_block_far | direct_sum | refine_2 | 7.042e-07 | 1.634e-09 | 0.0134 |
| gravity_block_far | direct_sum | refine_4 | 7.656e-07 | 1.776e-09 | 0.0145 |
| gravity_block_far | spectral | degree_4 | 5.275e-03 | 1.224e-05 | 0.1160 |
| gravity_block_far | spectral | degree_8 | 3.727e-05 | 8.647e-08 | 0.1161 |
| gravity_block_far | spectral | degree_12 | 2.043e-07 | 4.740e-10 | 0.1238 |
| gravity_block_far | spectral | degree_24 | 1.921e-08 | 4.456e-11 | 0.4407 |

![Convergence and cost](convergence.png)

Exact models, options, per-component RMSE and solver diagnostics are saved in `results.json`. Summed XYZ fields are saved in `fields.npz`. Residual units follow each case. `results.csv` provides the compact sweep table. Plot labels show the first and final node orders, direct_sum refine factors or spectral degrees. The separate diamond labelled 8×4 uses order 8 with four panels in each dimension.
