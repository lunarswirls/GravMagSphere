# Input-relative output paths

Numeric products use `output/` and figures use `figs/` under the input file's directory. The location of the launcher and the current working directory do not determine these defaults.

For example:

```text
examples/earth_crust_cases/
  bangui_equivalent_dipoles.in
  output/
    bangui_equivalent_dipoles.txt
  figs/
    bangui_equivalent_dipoles_bxyz.png
```

Directories are created when files are written. A derived input already inside an `output/` or `figs/` tree belongs to that tree's parent case directory: conversion does not create `output/output/`, and plotting a saved equivalent-source CSV uses sibling `figs/` rather than `output/figs/`. For multiple observation CSVs, the first input determines the default case directory; Python observation plots/saves use the first retained input in `source_file` metadata.

Explicit destinations override defaults. Multi-case diagnostic workflows use their configured diagnostic workspace destinations.

The Reiner Gamma case uses `reiner_gamma_test/output/` for tables and caches and singular `reiner_gamma_test/fig/` for its case-specific plots. Other case defaults use `figs/`.

## Fortran executables

| Program | Omitted output argument |
|---|---|
| `gravmag_sphere_direct_sum R input.in` | `input-directory/output/input.txt` |
| `gravmag_sphere_spectral R input.in` | `input-directory/output/input_spectral.txt` |
| `gravmag_sphere_gauss_legendre R input.in` | `input-directory/output/input_quadrature.txt` |
| `gravmag_xyz_to_brtp input_xyz.txt` | `case-directory/output/input_brtp.txt` |
| `gravmag_sphere_dipole_grid_fit R observations.csv` | `input-directory/output/observations_dipole_fit_predictions.csv` and `observations_dipole_fit_dipoles.csv` |

Pass an empty string (`""`) in an output argument's positional slot to request its default while supplying later solver options. Output paths are command-line arguments, not input-card fields. Shell wrappers use `_xyz`, `_spectral`, and `_brtp` filename suffixes according to the selected workflow.

`fortran/gravmag_paths.f90` handles defaults and directory creation for the five executables. Builders include it automatically. Custom compile commands must compile it before the program using it, with the same module include/output directory. Directory creation uses POSIX `mkdir`, with paths shell-quoted to preserve spaces and special characters. The orbital shared-library kernel returns arrays and does not create output directories.

## Python

```python
from gravmagpy import run_grid_model
from gravmagpy.plotting import plot_fields
from gravmagpy.utils import artifact_path

input_path = "examples/earth_crust_cases/bangui_equivalent_dipoles.in"
run = run_grid_model(input_path, radius_km=6371.2)
figure = plot_fields(input_path, run["output_path"])
chosen_figure = artifact_path(input_path, kind="figs", suffix="_comparison.png")
```

`run_grid_model` accepts an omitted output path. `plot_fields`, `plot_fit`, and `plot_equivalent_maps` accept an omitted image path. `save_fit` accepts an omitted directory when observation `source_file` metadata are available. In-memory observations/solutions without file provenance require an explicit save or figure destination. `run_sphere_model` and `fit_equivalent_sources` return arrays without automatically persisting files; supply `output_dir` to save their products.

The batch Python and Bash runners default to `<examples-dir>/output` and, for the Python plotting runner, `<examples-dir>/figs`. Explicit batch directory overrides are launcher-relative when supplied as relative paths; absolute paths are used as supplied. Single-case file arguments are relative to the caller's working directory.

Each single-case Bash launcher resolves its input-relative paths directly, accounting for `output/` and `figs/` ancestors to avoid nested output directories.

LPMAG wake examples save under the configured `wake_dir/output/<workflow>` and `wake_dir/figs/<workflow>`. `plot_equivalent_geotiffs.py` reads `wake_dir/output/equivalent_wake/dipoles.csv` and writes into `wake_dir/figs/equivalent_geotiffs`; set `source_path` to select a different saved solution. Synthetic and dictionary-generated examples have no external input file and use local `output/` and `figs/` directories. See the architecture's [LPMAG reference workflow](architecture.md#lpmag-reference-workflow) for input conventions.
