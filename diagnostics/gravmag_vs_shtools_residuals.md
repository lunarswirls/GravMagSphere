# GravMag Spectral vs SHTOOLS Residuals

This table records pairwise prediction differences for the specified cases. Magnetic errors use nT and gravity errors use mGal. These are model-to-model discrepancies, not errors against an exact physical solution; see [residual methodology](external_solver_residuals_method.md).

- Residual definition: `gravmag_sphere_spectral - shtools_lsq`
- lmax: `24`

| Case | RMSE(Br) | RMSE(Btheta) | RMSE(Bphi) | RMSE(Btot) | Plot |
|---|---:|---:|---:|---:|---|
| gravmag_sphere_1body_grav_fixedlim_incna_decna | 8.062342e+00 | 2.216018e+00 | 5.200876e+00 | 8.079534e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_grav_fixedlim_incna_decna_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_grav_polygon_incna_decna | 8.062342e+00 | 2.216018e+00 | 5.200876e+00 | 8.079534e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_grav_polygon_incna_decna_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_mag_fixedlim_inc-45_dec120_orient | 4.157158e+00 | 3.928429e+00 | 4.160588e+00 | 5.844792e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_mag_fixedlim_inc-45_dec120_orient_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_mag_fixedlim_inc0_dec0_orient | 1.060396e+01 | 2.987727e+00 | 2.143517e+00 | 8.665554e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_mag_fixedlim_inc0_dec0_orient_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_mag_fixedlim_inc0_dec90_orient | 2.450778e+00 | 1.253089e+00 | 4.364119e+00 | 4.967284e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_mag_fixedlim_inc0_dec90_orient_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_mag_fixedlim_inc45_dec45_orient | 8.168921e+00 | 4.187384e+00 | 3.258126e+00 | 7.822949e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_mag_fixedlim_inc45_dec45_orient_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base | 5.785288e+00 | 5.944657e+00 | 1.324621e+00 | 7.005404e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_mag_fixedlim_inc90_dec0_orient | 3.647763e+00 | 6.924468e+00 | 1.185981e+00 | 5.730440e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_orient_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_mag_fixedlim_inc90_dec0_weak | 8.614731e-01 | 1.488693e+00 | 3.311515e-01 | 1.061427e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_weak_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge | 1.393902e+02 | 7.475404e+01 | 8.025794e+01 | 1.645855e+02 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge_weak | 2.865486e+01 | 2.206457e+01 | 1.495318e+01 | 3.602616e+01 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge_weak_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_1body_mag_polygon_inc90_dec0_base | 3.311341e+00 | 5.656927e+00 | 1.305149e+00 | 4.305883e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_1body_mag_polygon_inc90_dec0_base_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_3body_mag_polygon_incmix_decmix | 1.141458e+01 | 4.541497e+00 | 4.248007e+00 | 1.275963e+01 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_3body_mag_polygon_incmix_decmix_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_3body_mag_polygon_incmix_decmix_weak | 2.642591e+00 | 1.323186e+00 | 1.036070e+00 | 3.075874e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_3body_mag_polygon_incmix_decmix_weak_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_5body_mag_polygon_incmix_decmix | 5.510561e+00 | 2.579863e+00 | 3.894314e+00 | 6.634913e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_5body_mag_polygon_incmix_decmix_gravmag_minus_shtools_residuals.png |
| gravmag_sphere_5body_mag_polygon_incmix_decmix_weak | 1.322569e+00 | 6.944255e-01 | 9.874385e-01 | 1.627702e+00 | /Users/danywaller/code/GravMagPy/diagnostics/gravmag_vs_shtools_residuals/gravmag_sphere_5body_mag_polygon_incmix_decmix_weak_gravmag_minus_shtools_residuals.png |

## Methodology Differences and Residual Drivers

- `gravmag_sphere_spectral` in this comparison uses GravMagSphere's physical-source-to-SH pipeline with its own regularization controls
- `shtools_lsq` uses direct SH least-squares expansion from sampled field values

Residuals increase when:

- Polygon boundaries are sharp relative to resolved harmonic degree
- Solver-specific regularization suppresses high-degree content differently
- Multi-body anomalies compete for the same low-degree coefficients in a global fit

Large-complex polygon cases are worst because discontinuous geometry excites high-degree terms that truncated SH models cannot represent without ringing. Sad :(
