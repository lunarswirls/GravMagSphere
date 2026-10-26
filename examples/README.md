# GravMagSphere Examples

Example Fortran + Python workflow for spherical gravity/magnetic forward modeling.

- Fortran solvers produce cartesian field components (`Fx Fy Fz`) on a lon/lat grid
- A Fortran converter rotates cartesian components (`Bx By Bz`) to spherical components (`Br Btheta Bphi`)
- Python scripts create 2D PNG maps and 3D Plotly HTML views

GravMagPy (`gravmagpy`) is the Python frontend for GravMagSphere. Reusable compilation/run functions are exposed by `gravmagpy.utils`, and the dictionary-based Python model interface is described in [the Python modeling guide](../docs/sphere.md).

## Directory layout

- `lunar_examples/`, `earth_crust_cases/`, and other case directories: input card files (`*.in`)
- `<input-directory>/output`: numeric outputs (`*.txt`, `*.csv`)
- `<input-directory>/figs`: figures (`*.png`, `*.html`)

Solver comparison scripts, formatted external inputs, and reports are in the top-level [diagnostics folder](../diagnostics/README.md). See [output path conventions](../docs/output_paths.md) for default filenames and explicit overrides.

## Build

```bash
cd /Users/danywaller/code/GravMagPy/examples
./build_gravmag_tools.sh
```

Build outputs:

- `./gravmag_sphere_direct_sum`
- `./gravmag_sphere_spectral`
- `./gravmag_sphere_gauss_legendre`
- `./gravmag_xyz_to_brtp`

## Compiled executables

### 1) Direct summation solver (Fortran)

`gravmag_sphere_direct_sum` computes XYZ fields by summing magnetic surface charges or gravity point masses

```bash
./gravmag_sphere_direct_sum <radius_km> <input.in> [output_xyz.txt] [refine_factor] [source_nlat] [source_nlon] [source_nr]
```

Example:

```bash
./gravmag_sphere_direct_sum 1737.4 \
  lunar_examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in \
  lunar_examples/output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_xyz.txt \
  2 0 0 0
```

### 2) Spectral solver (Fortran)

`gravmag_sphere_spectral` fits spherical harmonic coefficients and evaluates XYZ fields

```bash
./gravmag_sphere_spectral <radius_km> <input.in> [output_xyz.txt] \
  [lmax] [refine_factor] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] \
  [source_nlat] [source_nlon] [source_nr] [auto_mode] [joint_strength] [edge_correction] \
  [hybrid_mode] [hybrid_band_deg] [complex_vertex_threshold] [hybrid_transition_deg]
```

Defaults:

- `auto_mode=1` (pilot sweep auto-selects `lmax/reg_lambda/reg_power`)
- `joint_strength=1.0` (joint Br/Btheta/Bphi fit)
- `edge_correction=1` (local edge RBF correction)
- `hybrid_mode=1` (auto: enable hybrid for complex or multi-body inputs)
- `hybrid_band_deg=1.5`
- `complex_vertex_threshold=12`
- `hybrid_transition_deg=0.75`

Example:

```bash
./gravmag_sphere_spectral 1737.4 \
  lunar_examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in \
  lunar_examples/output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_spectral_xyz.txt \
  24 2 72 144 0.2 4.0 0 0 0 1 1.0 1 1 1.5 12 0.75
```

### Spectral Flags and Use Cases

`gravmag_sphere_spectral`, `run_input_to_xyz.sh`, and `run_gravmag_sphere_spectral.sh` expose:

- `auto_mode` (`0|1`):
  - `1` auto-selects `(lmax, reg_lambda, reg_power)` from a pilot sweep.
  - `0` uses user-provided seed values directly.
- `joint_strength` (`>=0`):
  - weight for joint vector fit (`Br/Btheta/Bphi`) in coefficient solve.
  - `0` reduces to radial-only fit.
- `edge_correction` (`0|1`) and `--no-local-correction`:
  - enable/disable local edge RBF correction.
- `hybrid_mode` (`0|1|2`):
  - `0`: spectral-only evaluation.
  - `1`: auto hybrid (enabled when input is complex or multi-body).
  - `2`: force hybrid for all inputs.
- `hybrid_band_deg`:
  - distance from polygon edge where direct_sum kernel result dominates.
- `hybrid_transition_deg`:
  - blend width between direct_sum and spectral predictions.
- `complex_vertex_threshold`:
  - polygon vertex threshold used by complexity classifier in hybrid auto mode.

Typical use cases:

```bash
# simple body, spectral-only evaluation
./run_input_to_xyz.sh spectral 1737.4 lunar_examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in \
  lunar_examples/output/simple_spectral_only_xyz.txt \
  --auto-mode 0 --lmax 24 --reg-lambda 0.2 --reg-power 4 \
  --edge-correction 0 --hybrid-mode 0

# complex polygon with automatic selection, local correction, and hybrid evaluation
./run_input_to_xyz.sh spectral 1737.4 lunar_examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in \
  lunar_examples/output/complex_auto_hybrid_xyz.txt

# force hybrid and widen the direct_sum edge band
./run_input_to_xyz.sh spectral 1737.4 lunar_examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in \
  lunar_examples/output/complex_force_hybrid_xyz.txt \
  --hybrid-mode 2 --hybrid-band-deg 2.0 --hybrid-transition-deg 1.0

# collective fit for bodies sharing the card-2 grid
./run_input_to_xyz.sh spectral 1737.4 lunar_examples/gravmag_sphere_3body_mag_polygon_incmix_decmix.in \
  lunar_examples/output/multibody_collective_xyz.txt --auto-mode 1
```

### 3) XYZ -> spherical converter (Fortran)

`gravmag_xyz_to_brtp` converts XYZ components to `Br Btheta Bphi`

```bash
./gravmag_xyz_to_brtp <input_xyz.txt> [output_brtp.txt]
```

Example:

```bash
./gravmag_xyz_to_brtp \
  lunar_examples/output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_xyz.txt
```

Omitting an output argument uses the input-relative default. To supply later positional controls with a default output, pass `""` in the output slot. The converter writes the example above to `lunar_examples/output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_brtp.txt`.

## Batch shell scripts

### Single-step wrappers

1. `run_input_to_xyz.sh` (solver step)

```bash
./run_input_to_xyz.sh <solver> <R_sphere_km> <input.in> [xyz_output] \
  [refine_factor] [lmax] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] \
  [source_nlat] [source_nlon] [source_nr] [auto_mode] [joint_strength] [edge_correction] \
  [hybrid_mode] [hybrid_band_deg] [complex_vertex_threshold] [hybrid_transition_deg] \
  [--lmax N] [--reg-lambda X] [--reg-power P] [--auto-mode 0|1] [--joint-strength X] \
  [--edge-correction 0|1] [--hybrid-mode 0|1|2] [--hybrid-band-deg X] \
  [--complex-vertex-threshold N] [--hybrid-transition-deg X] [--no-local-correction]
```

Examples:

```bash
./run_input_to_xyz.sh direct_sum 1737.4 \
  lunar_examples/gravmag_sphere_1body_grav_fixedlim_incna_decna.in

./run_input_to_xyz.sh spectral 1737.4 \
  lunar_examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in

./run_input_to_xyz.sh spectral 1737.4 \
  lunar_examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in \
  --lmax 30 --reg-lambda 0.1 --reg-power 4 --hybrid-mode 0 --no-local-correction
```

2. `run_xyz_to_brtp.sh` (conversion step)

```bash
./run_xyz_to_brtp.sh <input_xyz.txt> [output_brtp.txt]
```

### End-to-end wrappers

1. `run_gravmag_end_to_end.sh` (single case, direct_sum solver)

```bash
./run_gravmag_end_to_end.sh <R_sphere_km> <input.in> [xyz_out] [brtp_out] \
  [refine_factor] [source_nlat] [source_nlon] [source_nr]
```

2. `run_all_examples_end_to_end.sh` (all examples)

```bash
./run_all_examples_end_to_end.sh [solver] [R_sphere_km] [examples_dir] [output_dir] \
  [refine_factor] [lmax] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] \
  [source_nlat] [source_nlon] [source_nr]
```

Example:

```bash
./run_all_examples_end_to_end.sh spectral 1737.4 lunar_examples "" 2 24 72 144 0.2 4.0 0 0 0
```

### Solver-specific wrappers

```bash
./run_gravmag_sphere_direct_sum.sh <R_sphere_km> <input.in> [xyz_output] [refine] [source_nlat] [source_nlon] [source_nr]
./run_gravmag_sphere_spectral.sh <R_sphere_km> <input.in> [xyz_output] [lmax] [refine] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] [source_nlat] [source_nlon] [source_nr] [auto_mode] [joint_strength] [edge_correction] [hybrid_mode] [hybrid_band_deg] [complex_vertex_threshold] [hybrid_transition_deg] [--lmax N] [--reg-lambda X] [--reg-power P] [--auto-mode 0|1] [--joint-strength X] [--edge-correction 0|1] [--hybrid-mode 0|1|2] [--hybrid-band-deg X] [--complex-vertex-threshold N] [--hybrid-transition-deg X] [--no-local-correction]
```

## Visualization scripts

### 2D maps for all examples

`run_all_examples_brtp.py` can run solver + converter + plotting in one Python workflow

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python run_all_examples_brtp.py \
  --solver spectral \
  --rsphere-km 1737.4 \
  --examples-dir lunar_examples
```

Outputs:

- numeric tables in `lunar_examples/output`
- PNG figures in `lunar_examples/figs` (`*_brtp_2x2.png`)

The default example directory is `lunar_examples`. Batch directory overrides are relative to the launcher directory unless absolute paths are supplied.

### 3D Plotly for one source case

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python plot_single_case_3d_plotly.py \
  lunar_examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in \
  lunar_examples/output/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge_brtp.txt \
  --sparse-step 8 --max-vectors 450 --vector-scale 0.06
```

Omitting the HTML argument writes `<input-stem>_3d.html` under the input directory's `figs/`. An explicit HTML path is relative to the caller's working directory unless absolute.

### Fixed-limit vs polygon comparison plot

Use `diagnostics/plot_fixed_polygon_compare.py` with the two field tables. The [diagnostics guide](../diagnostics/README.md#5-fixed-limit-versus-polygon-maps) provides the repository-root command and output convention.

### Helper API for custom plotting

Use `gravmag_plot_helpers.py` from scripts or notebook:

```python
from pathlib import Path
from gravmag_plot_helpers import plot_example

fig = plot_example(
  "gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base",
  root_dir=Path("."),
  examples_dir="lunar_examples",
)
figs_dir = Path("lunar_examples/figs")
figs_dir.mkdir(parents=True, exist_ok=True)
fig.savefig(figs_dir / "custom_example.png", dpi=250)
```

## Notebook

Open:

- `gravmag_examples_visualization.ipynb`

The notebook is organized as:

1. setup/import cell
2. one markdown explanation + one plotting cell per example
3. each plot cell renders Br/Btheta/Bphi/Btot for that case

Run with Jupyter from this directory:

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python -m jupyter notebook gravmag_examples_visualization.ipynb
```

## Typical workflow

```bash
# 1) build
./build_gravmag_tools.sh

# 2) run all lunar_examples (numeric output)
./run_all_examples_end_to_end.sh spectral

# 3) generate figures
/Users/danywaller/code/venvs/gravmagpy/bin/python run_all_examples_brtp.py --solver spectral
```
