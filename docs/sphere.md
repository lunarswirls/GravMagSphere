# Python modeling with GravMag Sphere

GravMag Sphere Fortran programs live directly in `fortran/` and support both executable/input-card workflows and a named, dictionary-based Python interface. Python handles model construction, validation, file translation, and result assembly. Fortran performs the forward calculations and equivalent-dipole fit.

## Blocks, polygons, and equivalent volume grids

```python
import numpy as np
from gravmagpy import block, polygon, equivalent_source_grid, sphere_model, run_sphere_model

sources = [
    block(latitude_deg=(6, 7), longitude_deg=(-61, -60), depth_km=(5, 10),
          magnetization_a_m=(1, -0.5, 0.2), name="western source"),
    polygon(vertices_lat_lon=[(7, -60), (8, -59), (7, -58)], depth_km=(3, 8),
            magnetization_a_m=(-1, 0.2, 0.5), name="eastern source"),
]
model = sphere_model(sources, latitude_deg=np.linspace(5, 9, 21),
                     longitude_deg=np.linspace(-62, -57, 26), altitude_km=40, radius_km=1737.4)
result = run_sphere_model(model, solver="point_source", solver_options={"refine_factor": 2})
bx, by, bz, btot = (result[key] for key in ("bx_nt", "by_nt", "bz_nt", "btot_nt"))
```

`block` and `polygon` return editable dictionaries, and `sphere_model` copies and validates them. Validation is repeated before export/execution, so edits to a model are checked. Models can also be serialized with the standard `json` module. Unknown source or solver-option keys are rejected to catch misspellings.

Latitude/longitude are degrees, east-positive. Positions, source depths, and reference radius use km. Depth bounds are `(top, bottom)` below the reference surface, not below the spacecraft. Magnetization is a uniform **planet-fixed Cartesian vector in A/m**, not a local north/east/down vector. For gravity, supply `density_kg_m3` instead of magnetization. Magnetic and gravity sources require separate models because their output units differ.

Observation axes must be increasing, uniformly spaced one-dimensional arrays. Use an unwrapped longitude axis for an antimeridian-crossing grid, for example `[179, 180, 181]`; do not repeat the same meridian at both ends of a global grid. Individual sources must avoid the poles and span less than 180 degrees longitude. Split larger regions into cells. Polygon vertices are ordered boundary coordinates, optionally with the first vertex repeated to close the polygon; use simple, non-self-intersecting footprints.

`equivalent_source_grid` generates finite-volume cells from edges:

```python
magnetization = np.zeros((1, 2, 3, 3))
magnetization[..., 0] = [1, -1, 1]
sources = equivalent_source_grid(
    latitude_edges_deg=[6, 7, 8], longitude_edges_deg=[-62, -61, -60, -59],
    depth_edges_km=[5, 10], magnetization_a_m=magnetization,
    mesh={"radial": 4, "latitude": 12, "longitude": 12},
)
```

Property arrays broadcast to `(ndepth, nlat, nlon, 3)` for magnetization or `(ndepth, nlat, nlon)` for density. A single vector/scalar gives a uniform layer. Cells are returned in depth-major, then latitude, then longitude order. This function constructs a source model; it does not infer the properties from observations. A finite-volume magnetization in A/m is distinct from a point-dipole moment in A m².

## Solver controls and results

`run_sphere_model` generates input cards and executes the selected Fortran tool. No user-managed files are required. Pass `output_dir="output/my_model"` to retain `model.in` and `field.txt` for inspection or plotting.

Point-source options: `refine_factor`, `source_nlat`, `source_nlon`, and `source_nr`. The last three default to zero, meaning the per-source `mesh` settings are used. Default mesh is `radial=4, latitude=16, longitude=16`; point_source horizontal sampling is multiplied by `refine_factor`.

The point_source magnetic solver orients side-wall surface charges from the polygon winding in unwrapped longitude/latitude coordinates. This supports concave footprints in either vertex order, with optional repeated closing vertices. Horizontal and radial refinement still control integration accuracy.

Spectral options use the Fortran names: `lmax`, `refine_factor`, `ntheta_fit`, `nphi_fit`, `reg_lambda`, `reg_power`, `source_nlat`, `source_nlon`, `source_nr`, `auto_mode`, `joint_strength`, `edge_correction`, `hybrid_mode`, `hybrid_band_deg`, `complex_vertex_threshold`, and `hybrid_transition_deg`. Defaults match the executable. Set `auto_mode=0` to prevent its automatic parameter selection. Compiler/cache/timeout options can be passed through as keyword arguments.

Spectral gravity includes the degree-zero radial acceleration `-G * total_mass / radius**2`, using the signed sum of all discretized source masses. This coefficient is fixed before fitting degrees 1 through `lmax`, so regularization preserves the mass term. Source-mesh resolution controls the mass approximation. Magnetic fields retain zero degree-zero coefficient.

`solver="gauss_legendre"` uses nested longitude, latitude and radius Gauss–Legendre volume integration. It supports magnetic and gravity blocks and simple concave polygons in double precision through `gravmag_sphere_gauss_legendre`. The separate `gravmag_sphere_spectral` executable performs the spherical-harmonic fit.

```python
result = run_sphere_model(model, solver="gauss_legendre", solver_options={
    "radial_order": 8, "latitude_order": 16, "longitude_order": 16,
    "subdivisions": 2,
})
```

The three orders default to zero, inheriting each source's `mesh` radial/latitude/longitude counts (card 3). Explicit orders are integers 1–256. `subdivisions` defaults to 1 and divides each integration interval into that many equal panels in every dimension; it accepts integers 1–256. Cost grows approximately as the product of the three orders times `subdivisions**3`, multiplied by the observation count and polygon longitude slabs. Increase orders or split intervals until the fields converge, particularly for shallow sources and surface observations.

Polygon edges are straight in unwrapped longitude/latitude coordinates. Integration splits at vertex longitudes and pairs all latitude crossings, avoiding a masked bounding-box raster. Simple polygons may be clockwise or counterclockwise, closed or open, and cross the longitude seam. Holes and self-intersecting footprints are outside the supported geometry. Sources must avoid the poles and span less than 180° longitude. The observation sphere must lie strictly above every source top (`altitude_km + depth_top_km > 0`); surface observations are supported for buried sources.

This restores the numerical method from `gravmag_sphere_brtp.f` at historical Git revision `71d13a1`, using computed roots/weights, modern SI properties and XYZ output. It does not reproduce the old single-precision program or its obsolete card conventions bit for bit. Cards 4 and 6 remain compatibility metadata. See the [executed comparison report](../diagnostics/gauss_legendre_comparison/README.md) for convergence, timing and spectral limitations.

Returned `field` has shape `(nlat, nlon, 3)` in the input grid order; `total` has shape `(nlat, nlon)`. Magnetic aliases are `bx_nt`, `by_nt`, `bz_nt`, `btot_nt`; gravity aliases are `gx_mgal`, `gy_mgal`, `gz_mgal`, `gtot_mgal`. Point-source and Gauss–Legendre body contributions are summed as vectors before calculating total magnitude. Spectral output already represents the collective model and is not multiplied by the source count. `raw_output`, `stdout`, and `stderr` are available for diagnostics.

For magnetic plotting, retain outputs and call `gravmagpy.plotting.plot_fields(result['input_path'], result['output_path'])`. Figures default to the owning input directory's `figs/`; an explicit third `image_path` argument overrides it. `run_grid_model(input_path)` defaults numeric output to `output/<input-stem>.txt` for the point_source solver. See [output paths](output_paths.md).

## Input-card interoperability

```python
from gravmagpy import read_model_input, write_model_input, run_grid_model

model = read_model_input("sources.in", radius_km=1737.4)
model["sources"][0]["magnetization_a_m"] = [1, 0, 0]
write_model_input(model, "edited.in")
run_grid_model("edited.in", radius_km=1737.4)
```

Import/export preserves physical block/polygon geometry, titles, and compatibility values from cards 4 and 6. Magnetic amplitude/inclination/declination are converted to/from Cartesian magnetization; numerical equivalence rather than identical formatting is intended. Dictionary import requires a shared observation grid and the geometry/mesh constraints described above. Card workflows outside those dictionary constraints use `read_input`, `run_grid_model`, the executables, or shell runners directly.

## Fit an equivalent point-dipole model

```python
from gravmagpy import load_lpmag_csv, fit_equivalent_sources

observations = load_lpmag_csv(paths, latitude_range=(-27.5, -22.5),
                             longitude_range=(-107.5, -102.5), sigma_nt=1.0)
result = fit_equivalent_sources(
    observations, depth_km=10, spacing_deg=(1, 1), layers=1,
    padding_deg=(0.5, 0.5), regularization=0.0, output_dir="output/equivalent_fit",
)
```

This calls `gravmag_sphere_dipole_grid_fit`, fitting three Cartesian moments per grid point jointly to all observation altitudes. Grid coverage is determined from the observations plus latitude/longitude padding. `depth_km` locates the shallowest source layer below the reference surface; deeper layers use `layer_spacing_km`. Python translates this to the executable's depth-below-lowest-observation convention. `max_memory_mib` preflights the dense design and normal-equation arrays; coarse spacing is important for larger regions.

The returned `source_table` includes longitude, latitude, radius, surface depth, and `mx_am2`, `my_am2`, `mz_am2` moments in A m². `predictions` retains input observation order and, when available, file/row/time provenance. `predicted_nt`, `residual_nt`, `rmse_nt`, and `component_rmse_nt` are available as arrays/scalars. Optional output includes `dipoles.csv` and `predictions.csv`. Dipoles are not silently exported as finite-volume `.in` bodies: that conversion requires an explicit volume and magnetization assumption.

The Fortran fitter solves unweighted normal equations without background offsets. Nonuniform observation uncertainties are rejected. It does not use the nonlinear geometry optimizer's priors, background groups, robust losses, or bounds. Source depth/grid geometry remain fixed during this fit; this is an equivalent representation, not evidence of unique geological geometry.

### Damping convention

The Fortran diagonal penalty is

```text
(regularization + 1e-12) * max(mean(diag(transpose(g) @ g)), 1e-30)
```

Here `g` maps A m² moments to Tesla. The `1e-30` scale floor can dominate at orbital distances; consequently the executable's default `regularization=1e-4` can strongly suppress the modeled field. The Python wrapper defaults to **zero additional damping**, but the fixed `1e-12` term applies at that setting. This is not a promise of zero bias or exact recovery. Compare residuals and moment sizes across depths, spacings, and damping choices.

Runnable workflows are `examples/python_models.py` and `examples/lpmag/equivalent_wake.py`.

Use `predict_equivalent_field` or `equivalent_field_grid` to evaluate fitted moments at new positions/altitudes. [GeoTIFF context mapping](maps.md) displays spectral/context panels with field contours and source locations; see `examples/lpmag/plot_equivalent_geotiffs.py`.
