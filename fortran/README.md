# GravMagSphere Fortran backends

`orbital.f90` provides a C-interoperable double-precision magnetic volume kernel for arbitrary orbital positions. The GravMagPy frontend compiles and calls it through `gravmagpy.predict_field` and `gravmagpy.fit_sources`.

The same module exposes `dipole_field` for evaluating fitted equivalent-source moments at arbitrary positions through `predict_equivalent_field`.

It also exposes `monopole_kernel` for the regional Ravat-style L2resL1 pipeline.
The kernel maps equivalent charges in nT km² to XYZ fields in nT; Python handles
endpoint frame rotations, difference covariance, robust IRLS and surface-Br
regularization. See [the method documentation](../docs/ravat.md).

`gravmag_sphere_direct_sum.f90` evaluates extended blocks and polygonal bodies by direct summation. Its magnetic branch uses surface-charge quadrature, placing equivalent charges at surface-patch centers and triangle centroids. Its gravity branch uses midpoint-volume quadrature, placing each cell's mass at its midpoint. The Python target is `direct_sum`; `direct` remains a selector alias.

The Python model interface (`block`, `polygon`, `equivalent_source_grid`, `sphere_model`, and `run_sphere_model`) generates cards internally and runs `gravmag_sphere_direct_sum`, `gravmag_sphere_spectral` or `gravmag_sphere_gauss_legendre`. `fit_equivalent_sources` calls `gravmag_sphere_dipole_grid_fit` using normalized observation CSVs. Lower-level `build_fortran`, `run_fortran`, and `run_grid_model` expose compilation, positional arguments, and card execution. The [shell build workflow](../examples/README.md) lives in `examples`.

The five executable programs share `gravmag_paths.f90`. Omitted or empty output arguments create files in an `output/` folder beside the input; explicit arguments select a destination. Compile this utility before each executable when building manually. Python and shell builders include it automatically. See [output path conventions](../docs/output_paths.md) for filenames and positional argument examples.

`gravmag_sphere_gauss_legendre.f90` dispatches three options: 1 is the original nested volume method, 2 is the recent Schmidt harmonic method, and 3 is the fully normalized geodetic harmonic method. Compile `gravmag_paths.f90`, `gravmag_sphere_harmonics.f90`, `gravmag_sphere_transform.f90`, `gravmag_sphere_normalized.f90`, then `gravmag_sphere_gauss_legendre.f90`. Enable OpenMP for the native transforms. Python and shell builders include these dependencies. No external Fourier library is needed.

Option 1 emits one XYZ block per body. Options 2/3 emit the collective XYZ field once. Magnetic units are nT and gravity units are mGal. Default normalization is 1, so existing positional calls retain the original volume calculation.

```bash
bash examples/run_gravmag_sphere_gauss_legendre.sh 1737.4 sources.in "" 8 16 16 2
```

Arguments are radius, input, optional output, radial order, latitude order, longitude order and subdivisions, followed by optional normalization `1|2|3`. Zero orders use card 3; subdivisions defaults to 1. Omitted output becomes `output/<input-stem>_quadrature.txt`. `GRAVMAG_DIAGNOSTICS=1` reports body timing, node counts, effective orders, subdivisions and integrated volume. The general `run_input_to_xyz.sh gauss_legendre` and shell batch launcher also support this solver: their existing `source_nr/source_nlat/source_nlon` options specify node orders, and `refine_factor` specifies subdivisions (default 2 in those launchers). Python option names and method limits are documented in [sphere.md](../docs/sphere.md).

For options 2/3, use `0 0 0 1` in the four volume-control positions, put normalization at position 8, then append the original harmonic controls (`lmax`, `refine_factor`, `ntheta_fit`, `nphi_fit`, `reg_lambda`, `reg_power`, `source_nlat`, `source_nlon`, `source_nr`, `auto_mode`, `joint_strength`, `edge_correction`, `hybrid_mode`, `hybrid_band_deg`, `complex_vertex_threshold`, `hybrid_transition_deg`). The volume positions must keep these placeholders in normalized modes. Harmonic controls can be omitted to use their existing defaults.

```bash
# Schmidt, degree 24, pure harmonic evaluation:
bash examples/run_gravmag_sphere_gauss_legendre.sh 1737.4 sources.in "" 0 0 0 1 2 24 2 72 144 0.2 4 0 0 0 0 1 0 0
# Fully normalized: same controls, normalization 3:
bash examples/run_gravmag_sphere_gauss_legendre.sh 1737.4 sources.in "" 0 0 0 1 3 24 2 72 144 0.2 4 0 0 0 0 1 0 0
```

`GRAVMAG_SPECTRAL_MAX_MIB` retains its name for compatibility and limits the workspace of normalized Gauss–Legendre transforms. Historical Schmidt validation is retained under `diagnostics/native_high_degree/`; reruns now select Gauss–Legendre option 2.
