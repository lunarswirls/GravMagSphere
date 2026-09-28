# Earth crust equivalent-source cases

These inputs provide three terrestrial magnetic benchmarks for a sphere radius of `6371.2 km`.

The GravMag Sphere `.in` format describes finite magnetized bodies rather than zero-volume point dipoles. Each fixed-limit body in these files represents one equivalent-source cell. With the Gauss solver, setting the three source controls to `1 1 1` reduces every cell to one volume dipole at its midpoint.

## Cases

### `bangui_equivalent_dipoles.in`

Default numeric output for the direct_sum executable is `output/bangui_equivalent_dipoles.txt` inside this case directory. Python figures default to `figs/` here. Explicit paths override these defaults; see [output path conventions](../../docs/output_paths.md).

- 163 half-degree source cells within 400 km of `4.37 deg N, 18.56 deg E`
- source layer from 3.0 to 7.5 km depth
- published disc-model magnetization of `10 A/m`
- local geophysical direction `I=+25 deg`, `D=-18 deg`, converted separately at every source to the solver's global Cartesian angle convention
- observations at 4 km altitude on a quarter-degree regional grid

This is a discretized version of the published 800 km diameter, 4.5 km thick Bangui disc interpretation. It is the most physically constrained case in this directory.

### `meyer_global_32400_dipoles.in`

- exactly 32,400 sources on a global `2 deg x 2 deg` grid
- two layers spanning 0–20 km and 20–40 km depth
- induced directions follow a centered axial-dipole field and are converted separately to global Cartesian angles
- ten reproducible proxy magnetization classes from 0.05 to 1.10 A/m
- observations at 400 km altitude on a 5-degree global grid

The grid, layer count and source count follow the published Meyer et al. geometry. Class placement and strengths are explicit synthetic proxies, not a transcription of the published thickness/susceptibility classes or individual moments. This is a solver stress case, not a reproduction of the published field solution.

### `mayhew_north_america_equivalent_layer.in`

- 1,218 approximately equal-area sources over North America
- about 220 km north-south spacing, with longitude spacing increased toward the pole
- 40 km magnetic layer and axial-dipole induced directions
- reproducible regional proxy magnetizations for the Canadian Shield, western thermal provinces, central continent and Appalachians
- observations at 400 km altitude on a 2-degree regional grid

The source spacing, equal-area design and layer thickness follow the published Mayhew configuration. Magnetization contrasts are documented proxies, not fitted or transcribed published source moments.

## Generate

Run the generator from the repository root using the selected virtual environment:

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python examples/earth_crust_cases/generate_earth_crust_cases.py
```

## Run

From `examples`, one-dipole-per-cell Gauss configurations are:

```bash
./gravmag_sphere_spectral 6371.2 earth_crust_cases/bangui_equivalent_dipoles.in "" 20 1 21 41 0.05 4 1 1 1 0
./gravmag_sphere_spectral 6371.2 earth_crust_cases/mayhew_north_america_equivalent_layer.in "" 12 1 15 31 0.05 4 1 1 1 0
./gravmag_sphere_spectral 6371.2 earth_crust_cases/meyer_global_32400_dipoles.in "" 5 1 7 13 0.05 4 1 1 1 0
```

The global case is intentionally much larger. Start with low harmonic degree and coarse fit controls, then increase them while monitoring memory.

## References

- Meyer, J., Hufen, J.-H., Siebert, M., and Hahn, A. (1983), *Investigations of the internal geomagnetic field by means of a global model of the Earth's crust*
- Mayhew, M. A. and Galliher, S. C. (1982), *An equivalent layer magnetization model for the United States derived from MAGSAT data*, doi:10.1029/GL009i004p00311
- Ravat, D. et al. (1993), *Improvement of equivalent source inversion technique with a more symmetric dipole distribution model*, doi:10.1016/0031-9201(93)90012-X
- Regan, R. D. and Marsh, B. D. (1982), *The Bangui magnetic anomaly: its geological origin*, doi:10.1029/JB087iB02p01107
