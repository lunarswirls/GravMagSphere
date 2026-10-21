# GravMagSphere Fortran backends

`orbital.f90` provides a C-interoperable double-precision magnetic volume kernel for arbitrary orbital positions. The GravMagPy frontend compiles and calls it through `gravmagpy.predict_field` and `gravmagpy.fit_sources`.

The same module exposes `dipole_field` for evaluating fitted equivalent-source moments at arbitrary positions through `predict_equivalent_field`.

`gravmag_sphere_direct_sum.f90` evaluates extended blocks and polygonal bodies by direct summation. Its magnetic branch uses surface-charge quadrature, placing equivalent charges at surface-patch centers and triangle centroids. Its gravity branch uses midpoint-volume quadrature, placing each cell's mass at its midpoint. The Python target is `direct_sum`; `direct` remains a selector alias.

The Python model interface (`block`, `polygon`, `equivalent_source_grid`, `sphere_model`, and `run_sphere_model`) generates cards internally and runs `gravmag_sphere_direct_sum`, `gravmag_sphere_spectral` or `gravmag_sphere_gauss_legendre`. `fit_equivalent_sources` calls `gravmag_sphere_dipole_grid_fit` using normalized observation CSVs. Lower-level `build_fortran`, `run_fortran`, and `run_grid_model` expose compilation, positional arguments, and card execution. The [shell build workflow](../examples/README.md) lives in `examples`.

The five executable programs share `gravmag_paths.f90`. Omitted or empty output arguments create files in an `output/` folder beside the input; explicit arguments select a destination. Compile this utility before each executable when building manually. Python and shell builders include it automatically. See [output path conventions](../docs/output_paths.md) for filenames and positional argument examples.

The spectral executable additionally requires `gravmag_sphere_harmonics.f90` followed by `gravmag_sphere_transform.f90` before `gravmag_sphere_spectral.f90`. These provide Fukushima-style extended-range Schmidt recurrences and the native Gauss–Legendre/Fourier projection. Enable the compiler's OpenMP option for parallel execution. The Python and shell builders include these modules and flags. No external Fourier library is needed.

`gravmag_sphere_gauss_legendre.f90` is the `gauss_legendre` target. It restores the original nested Gauss–Legendre volume method with double precision, generated nodes/weights, optional composite panels and polygon boundary integration. It uses the current magnetic/gravity cards and emits one XYZ block per body, in nT or mGal. It is distinct from the spectral `gravmag_sphere_spectral.f90` program.

```bash
bash examples/run_gravmag_sphere_gauss_legendre.sh 1737.4 sources.in "" 8 16 16 2
```

Arguments are radius, input, optional output, radial order, latitude order, longitude order and subdivisions. Zero orders use card 3; subdivisions defaults to 1. Omitted output becomes `output/<input-stem>_quadrature.txt`. `GRAVMAG_DIAGNOSTICS=1` reports body timing, node counts, effective orders, subdivisions and integrated volume. The general `run_input_to_xyz.sh gauss_legendre` and shell batch launcher also support this solver: their existing `source_nr/source_nlat/source_nlon` options specify node orders, and `refine_factor` specifies subdivisions (default 2 in those launchers). Python option names and method limits are documented in [sphere.md](../docs/sphere.md).
