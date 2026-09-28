# GravMagSphere diagnostics

This folder contains compiler benchmarks, direct_sum/spectral comparisons, equivalent-source comparisons, formatted inputs, and reports/figures. These workflows assess numerical behavior and modeling assumptions; the external comparison adapters are separate from the package's inversion interfaces.

## Contents

- `external_solver_inputs/`: sixteen formatted comparison cases, observations, metadata, and saved solver outputs
- `diagnostic_fortran.py`: compiler benchmarking and spectral-degree ringing sweeps
- `external_solver_compare.py`: external input generation and SciPy/SHTOOLS/Harmonica comparisons
- `hybrid_vertex_sweep.py`: polygon-complexity sweep
- `shtools_operating_map.py`: accuracy/runtime operating maps
- `complexlarge_direct_vs_spectral_analysis.py`: boundary-visibility analysis
- `plot_fixed_polygon_compare.py`: fixed-limit/polygon field comparison plot
- `run_comparison_tests.sh`: comparison-suite runner
- `gauss_legendre_compare.py`: node-order, panel, direct_sum-mesh and pure spectral degree sweeps
- `gauss_legendre_comparison/`: executed three-solver report, figure, inputs, fields and machine-readable metrics
- Markdown reports and residual/side-by-side plot directories: reference results for the recorded case settings

Relative paths for benchmark, comparison, and sweep options are resolved from the repository root. The fixed/polygon plotting script uses caller-relative input paths and defaults its image to `figs/fixed_vs_polygon_components.png` under the fixed input's owning case directory. The suite runner can be invoked from any directory.

## Environment and usage

The core diagnostics need NumPy, SciPy, and Matplotlib; SHTOOLS/Harmonica suites additionally need `pyshtools`/`harmonica`. The runner skips missing optional backends unless `--strict-deps` is supplied. Full benchmark/sweep runs rebuild executables and can take substantial time; `--help` is non-mutating.

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python diagnostics/external_solver_compare.py --help
bash diagnostics/run_comparison_tests.sh --help
```

## Workflows

### Gauss–Legendre restoration and three-solver comparison

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python -m unittest discover -s tests -p test_gauss_legendre.py -v
/Users/danywaller/code/venvs/gravmagpy/bin/python diagnostics/gauss_legendre_compare.py
bash diagnostics/run_comparison_tests.sh --python /Users/danywaller/code/venvs/gravmagpy/bin/python --tests quadrature --skip-build
```

The script accepts an optional output directory and a non-mutating `--help`; edit its `settings` dictionary for other sweep resolutions. It compiles through the package cache before timing runs. Default outputs are in [gauss_legendre_comparison](gauss_legendre_comparison/README.md): Markdown, CSV, JSON, compressed XYZ arrays, six card inputs and a convergence/cost figure. A failed reference-resolution check returns nonzero after saving the report.

Tests independently compare the restored tensor quadrature with the orbital backend, exact polygon volume, decomposed concave bodies, a point-mass limit, longitude wrapping and multiple source summation. Direct summation comparisons cover magnetic and gravity fields, including concave polygon winding and refinement. High-altitude magnetic and gravity cases check spectral-degree convergence; antipodal gravity sources isolate the degree-zero mass term and check its sign, distance scaling and preservation under regularization. The full sweep also reports near-surface resolution limits. Spectral auto selection, edge corrections and hybrid mode are explicitly disabled for this comparison.

### 1) Compiler runtime benchmark + hotspot report

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python diagnostics/diagnostic_fortran.py benchmark \
  --compilers gfortran,ifx,ifort,flang-new,nvfortran \
  --repeats 3
```

Outputs:

- `diagnostics/fortran_benchmark_report.json`
- `diagnostics/fortran_benchmark_report.md`

The benchmark report includes:

- wall time per compiler/case
- slowest solver stage (from internal Fortran timers)
- largest memory bucket estimate

### 2) Spherical harmonic `lmax` ringing sweep

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python diagnostics/diagnostic_fortran.py lmax-sweep \
  --compiler gfortran \
  --lmin 8 --lmax 96 --lstep 8
```

Outputs:

- `diagnostics/lmax_sweep.csv`
- `diagnostics/lmax_sweep.md`

The sweep compares spectral results to a direct_sum-solver baseline and reports
edge-focused ringing metrics (`edge_amp`, `edge_overshoot`) per `lmax`.

### 3) Per-run diagnostics from Fortran executables

Set `GRAVMAG_DIAGNOSTICS=1` to print stage timings and memory estimates:

```bash
GRAVMAG_DIAGNOSTICS=1 examples/gravmag_sphere_spectral 1737.4 examples/lunar_examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in diagnostics/diag_gauss_xyz.txt 24 2 72 144 0.2 4.0 0 0 0
GRAVMAG_DIAGNOSTICS=1 examples/gravmag_sphere_direct_sum 1737.4 examples/lunar_examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in diagnostics/diag_direct_xyz.txt 2 0 0 0
```

Diagnostics lines are emitted as `DIAG|<solver>|<category>|<key>=<value>`.

### 4) Run all solver comparison suites

`run_comparison_tests.sh` orchestrates all comparison workflows:

- external solver comparison (Fortran spectral vs SciPy SH LSQ vs SHTOOLS LSQ)
- GravMagSphere-vs-SHTOOLS residual suite
- GravMagSphere-vs-Harmonica XYZ residual suite
- complexlarge direct_sum-vs-spectral boundary analysis
- hybrid-vs-spectral vertex sweep

```bash
# run all comparison suites
bash diagnostics/run_comparison_tests.sh --python /Users/danywaller/code/venvs/gravmagpy/bin/python

# run selected suites only
bash diagnostics/run_comparison_tests.sh --python /Users/danywaller/code/venvs/gravmagpy/bin/python \
  --tests external,shtools,harmonica

# limit to specific lunar examples
bash diagnostics/run_comparison_tests.sh --python /Users/danywaller/code/venvs/gravmagpy/bin/python \
  --tests external,harmonica \
  --examples gravmag_sphere_1body_mag_polygon_inc90_dec0_base,gravmag_sphere_3body_mag_polygon_incmix_decmix
```

Notes:

- Set `--no-clean` to keep intermediate txt/csv artifacts
- Set `--strict-deps` to fail immediately if Python dependencies are missing

### 5) Fixed-limit versus polygon maps

After generating the corresponding BRTP tables:

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python diagnostics/plot_fixed_polygon_compare.py \
  --fixed examples/lunar_examples/output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_brtp.txt \
  --polygon examples/lunar_examples/output/gravmag_sphere_1body_mag_polygon_inc90_dec0_base_brtp.txt
```

Use `--out` to select a figure path explicitly. See [residual methodology](external_solver_residuals_method.md), [solver comparison](external_solver_comparison.md), and [coordinate conventions](gravmagsphere_inclination_declination_physical_audit.md) for interpretation. Benchmark timings and residual tables describe their recorded configurations, not universal accuracy or performance guarantees.
