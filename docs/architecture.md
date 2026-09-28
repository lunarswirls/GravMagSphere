# GravMagSphere architecture

GravMagSphere is the complete Fortran magnetic and gravity solver and inversion toolkit. GravMagPy (`gravmagpy`) is its Python frontend for source construction, spacecraft-data ingestion, inversion, and visualization. Users work with dictionaries, arrays, and tables; adapters handle input cards and executable calls when needed.

The numerical foundation uses a configurable spherical reference body. The orbital-data adapter targets wake-sorted Lunar Prospector magnetic measurements. The development direction is a body-independent equivalent-source framework with interchangeable forward backends and inversion methods. Those extensions are described separately from the implemented capabilities below.

## Repository components

```text
src/gravmagpy/              reusable Python package
  data.py                  LPMAG ingestion, selection, and observation provenance
  geometry.py              coordinates, spherical blocks, volumes, and geographic grids
  forward.py               in-process Fortran magnetic-field evaluation
  inversion.py             bounded block-geometry and magnetization fitting
  sphere.py                dictionary-based grid models and solver adapters
  equivalent.py            fixed-grid dipole fitting and prediction
  cards.py                 input-card parsing, export, and grid execution
  plotting.py              grid-field and orbital-fit figures
  maps.py                  equivalent-field maps and georeferenced raster panels
  utils/                   compilation, execution, path, and coordinate utilities
fortran/                   numerical kernels, executable programs, and path utilities
examples/                  configured workflows, input models, and notebooks
diagnostics/               cross-solver comparisons, benchmarks, and reference products
tests/                     numerical, interface, georeferencing, and path checks
data/spice_kernels/         kernels and metakernel for coordinate workflows
docs/                      modeling, data, mapping, and execution guides
```

The package owns validation and public interfaces. Fortran owns the numerical kernels. Examples choose data paths, source hypotheses, and run settings; the package does not import or execute them. Diagnostic integrations are comparison tools, not automatically available production inversion backends.

## Data and model contracts

### Observations

The magnetic inversion interfaces consume an observation dictionary with:

- `xyz_km`: `(n, 3)` Cartesian positions, retaining each sample's actual altitude
- `field_nt`: `(n, 3)` magnetic vectors in the same body-fixed axes as the positions
- `sigma_nt`: positive component uncertainties, broadcastable to `(n, 3)`
- `radius_km`: the spherical reference radius
- `table`, `frame`, and `report`: sample provenance, coordinate labeling, and selection diagnostics

`load_lpmag_csv` normalizes body-fixed SEL columns and retains source filenames, row indices, and available timestamps. It supports configurable column mappings and geographic position columns, but it is an LPMAG/SEL adapter, not a general mission-frame transformation service. Field variability such as `Brms` is a quality-selection quantity, not automatically an independent measurement uncertainty. The [LPMAG reference workflow](#lpmag-reference-workflow) supplies the input template and practical constraints for the mission-neutral development plan.

Python modeling positions and depths use km; magnetic fields use nT; magnetization uses A/m. Equivalent dipole moments use A m² and are not interchangeable with magnetization. The dipole CSV stores radius in metres, with explicit conversion at the prediction boundary. Grid gravity models use density contrast in kg/m³ and return acceleration in mGal.

### Sources

There are three distinct source representations:

- Orbital inversion blocks describe center latitude/longitude, angular widths, top depth, thickness, and uniform body-fixed Cartesian magnetization
- Grid models contain block or polygon dictionaries, material properties, a reference radius, and a regular observation grid; `equivalent_source_grid` constructs finite-volume cells
- Equivalent dipole solutions contain source positions and fitted Cartesian moments, optionally organized into depth layers

These representations have separate adapters and validation. A point-dipole solution is not implicitly converted into a volume model: that needs an explicit volume and material assumption. `read_model_input` and `write_model_input` connect supported dictionary models to Fortran input cards. Lower-level card/executable interfaces support card-based workflows. See [source modeling](sphere.md).

### Predictions and results

Orbital predictions are ordered `(n, 3)` vectors. Grid predictions are `(nlat, nlon, 3)` arrays with coordinate axes. Contributions are summed as vectors before total-field magnitude is computed. Bx, By, and Bz are body-fixed Cartesian components, not local north/east/down components.

Nonlinear fit results include predictions, residuals, source parameters, background estimates, convergence status, active bounds, and Jacobian-rank diagnostics. Equivalent-grid results include source moments, observation-aligned predictions, residuals, and grid settings. These dictionaries are related but do not yet implement a single solver-independent result schema.

## Numerical backends

| Build target | Interface | Role and supported evaluation |
|---|---|---|
| `orbital` | C-interoperable shared library | Double-precision magnetic spherical-block integration and point-dipole prediction at arbitrary external Cartesian positions |
| `direct_sum` | `gravmag_sphere_direct_sum` executable | Input-card block/polygon magnetic or gravity fields on a regular grid; single-precision field calculations |
| `spectral` | `gravmag_sphere_spectral` executable | Grid fields using a spherical-harmonic representation, with configurable local and hybrid corrections |
| `gauss_legendre` | `gravmag_sphere_gauss_legendre` executable | Double-precision magnetic/gravity volume quadrature for blocks and simple polygons, with configurable node orders and composite panels |
| `xyz_to_brtp` | `gravmag_xyz_to_brtp` executable | Rotation of Cartesian field tables into spherical components; no source inversion |
| `dipole_grid` | `gravmag_sphere_dipole_grid_fit` executable | Linear fitting of vector dipole moments on a fixed subsurface grid from orbital magnetic samples |

The orbital volume kernel integrates magnetic dipole contributions with tensor Gauss–Legendre quadrature and the spherical volume Jacobian. The direct_sum magnetic grid solver uses surface-charge calculations. Their outputs are checked against independent references and against each other under suitable refinement; their discretizations and precision are not assumed identical.

`forward.py` loads the shared library with `ctypes`. Contiguous `float64` arrays cross the C interface, with Python `(n, 3)` storage corresponding to Fortran `(3, n)` arrays. Optimizer evaluations reuse the loaded library and evaluate spacecraft positions directly, without interpolating a fixed-altitude map or launching an executable for every residual calculation.

`sphere.py` validates grid models, generates cards, translates named options, and collects summed field arrays. `equivalent.py` adapts the fixed-grid fitter's depth convention and evaluates fitted moments at additional positions through the shared-library dipole kernel.

## Fitting and mapping workflows

### Block-geometry inversion

`fit_sources` uses bounded `scipy.optimize.least_squares` to fit parameters selected by the bounds dictionaries. Unselected parameters and source count remain fixed. Available parameters include block position, widths, depth, thickness, and magnetization. The objective supports component uncertainty weighting, parameter priors, robust losses, and optional constant background vectors per file or observation group.

Every residual evaluation uses the samples' individual positions. The quadrature rule retains a fixed node topology as block dimensions change. Bounds keep source volumes admissible throughout the search. Convergence and Jacobian diagnostics describe the chosen optimization problem; they do not establish a unique subsurface interpretation.

### Fixed-grid equivalent sources

`fit_equivalent_sources` holds source geometry fixed and solves for vector moments. The Fortran implementation uses dense, unweighted regularized normal equations, includes a fixed diagonal damping floor, and has no fitted background offsets. The Python adapter rejects nonuniform uncertainties instead of silently discarding them. Its regularization parameter is not interchangeable with spectral damping; see [the damping convention](sphere.md#damping-convention).

`predict_equivalent_field` evaluates fitted moments at new positions. `equivalent_field_grid` evaluates a north-up geographic grid at a chosen altitude. This is model-based field continuation, not interpolation of mixed-altitude measurements. Depth, spacing, and damping remain modeling choices whose effects need validation.

### Visualization

`plotting.py` renders component grids and observation-aligned fit comparisons. `maps.py` places predicted fields and source locations alongside GeoTIFF products using their coordinate reference systems, affine transforms, masks, and scale metadata. It uses a spherical planetary geographic CRS and checks the raster's reference-body radius. Rasterio and Matplotlib are optional dependencies. This geographic map interface is regional, not a general polar or irregular-body mapping system. See [mapping](maps.md).

## Build, execution, and storage

`build_fortran` compiles the requested target with a compatible Fortran compiler. Source contents, compiler identity/version, flags, and platform determine the cache key. Temporary compilation directories isolate module files; completed binaries are published into the cache. The default cache is under the system temporary directory, and `build_dir` selects an explicit location.

The checkout uses `fortran/`; installed wheels contain sources under `gravmagpy/_fortran`. Compilation occurs on use, not on package import. A compiler and its runtime must be available on the execution machine. NumPy, SciPy, and pandas are core dependencies; plotting and geospatial integrations are optional extras.

`run_fortran` passes arguments as a subprocess argument list and checks both return status and Fortran error messages. Executable programs use `gravmag_paths.f90` for default destinations and POSIX directory creation. Bash launchers resolve paths directly. Numeric products default to input-relative `output/` and figures to `figs/`; Reiner Gamma's case scripts use `fig/`. Explicit destinations override defaults. In-memory modeling calls do not implicitly persist their temporary cards or results. See [output paths](output_paths.md).

## Path forward: body-independent, multi-solver equivalent sources

The LPMAG workflow below is implemented and provides the reference case for further development. The six numbered sections describe proposed interfaces and capabilities, not an implemented unified API. A **forward backend** evaluates a physical source representation; an **inversion method** estimates its parameters. These choices should be independent where their capabilities permit, rather than treating every backend as a complete, interchangeable fitting workflow.

### LPMAG reference workflow

Wake-sorted Lunar Prospector magnetic CSVs provide a concrete template for ingestion, multi-altitude fitting, and model-based mapping. The configured examples use the editable `wake_dir`:

```text
/Users/danywaller/Projects/moon/lpmag_l1b_5s_avg/lp_mag_shadow_state
```

The measurement CSVs are external inputs. [The synthetic example](../examples/lpmag/synthetic_fit.py) generates a compatible dataset at 20, 40, and 80 km without requiring mission files. [The block-fit example](../examples/lpmag/fit_wake.py) and [the equivalent-grid example](../examples/lpmag/equivalent_wake.py) use the configured wake directory. Numeric products are written to `wake_dir/output/<workflow>` and figures to `wake_dir/figs/<workflow>`.

#### Input columns and coordinates

| Columns | Meaning and units |
|---|---|
| `X_SEL`, `Y_SEL`, `Z_SEL` | body-fixed Cartesian position, km |
| `Bx_SEL`, `By_SEL`, `Bz_SEL` | magnetic vector in the same axes, nT |
| `utc` | optional timestamp retained as provenance |
| `wake_t_km`, `wake_rperp_km` | optional wake-selection metadata, km |
| `Brms` | optional field variability for quality selection, nT |

SEL is body-fixed: +z points north, +x intersects longitude zero, and +y intersects 90 degrees east, following the [PDS lunar-coordinate description](https://pds-ppi.igpp.ucla.edu/data/lp-mag-calibrated/document/data-lunarcrds-desc.txt). Positions and vectors must use this same frame. SSE components require a time-dependent rotation before ingestion; renaming columns does not transform them.

For Cartesian inputs, altitude is `norm([X_SEL, Y_SEL, Z_SEL]) - radius_km`, using a default lunar radius of 1737.4 km. A CSV altitude column does not override that calculation. Geographic position columns can instead be selected explicitly:

```python
from gravmagpy import load_lpmag_csv

observations = load_lpmag_csv(
    "selected_wake.csv",
    columns={
        "lat_deg": "latitude_deg", "lon_deg": "longitude_deg", "altitude_km": "height_km",
        "bx_nt": "Bx", "by_nt": "By", "bz_nt": "Bz",
    },
)
```

This mapping still requires body-fixed vectors, degrees for latitude/longitude, km for height, and nT for field. Convert metre/Tesla inputs before ingestion. The loader does not infer a frame rotation or unit conversion from column names.

#### Selection, uncertainty, and provenance

With `wake_only=True`, rows must have finite positive `wake_t_km` and finite nonnegative `wake_rperp_km` when both columns are present. This consumes an upstream wake classification; it does not recompute the wake cone. Without both columns, the filename must contain `wake` to declare preselection, or the caller must explicitly use `wake_only=False`. A filename containing only `nightside` is insufficient.

Non-numeric or nonfinite required values are excluded, and selection counts are returned in `observations["report"]`. The normalized table preserves the source filename, zero-based source row, available timestamp, positions, fields, and uncertainties. Header-only files are allowed; an entirely empty selection raises an error. Duplicate input paths are rejected, but overlapping files can contain duplicate measurements, which the caller must resolve.

Latitude, longitude, and altitude ranges select a region. A longitude interval with its lower endpoint greater than its upper endpoint crosses the antimeridian. `max_brms_nt` applies an absolute variability cutoff. `Brms` is not automatically interpreted as a calibrated component uncertainty. Supply a scalar or three-component `sigma_nt`, or three `sigma_columns` in nT; `sigma_nt=1` gives unweighted component residuals in nT.

#### Fitting and interpretation

`fit_sources` minimizes `(predicted - observed) / sigma` for all three components at every sample's actual position. Source dictionaries use `lat_deg`/`lon_deg` for the center, `lat_width_deg`/`lon_width_deg` for full angular extents, `depth_top_km` and `thickness_km` for radial geometry, and `mx_a_m`/`my_a_m`/`mz_a_m` for uniform body-fixed magnetization. Bounds select free parameters; unspecified parameters stay fixed.

Optional priors follow the source list-of-dictionaries structure, with `(mean, standard_deviation)` pairs. Their standardized residuals are appended to the data residuals. `loss` and `f_scale` control the least-squares loss; a nonlinear loss applies to prior residuals as well as data residuals.

`fit_background=True` estimates one constant Cartesian vector per file, or per label in `observations["groups"]`. Such offsets can absorb broad crustal signals and do not represent a time-varying external-field model. The `background_only_rmse_nt` baseline helps separate improvement due to sources from improvement due to offsets alone.

Depth, thickness, lateral extent, and magnetization can trade off. Multiple altitudes constrain spatial decay but do not establish uniqueness. The real-data block example fixes thickness while fitting location, widths, top depth, and magnetization. Check alternative initial guesses and bounds, and increase `quadrature_order` to assess numerical stability. `success` indicates optimizer termination, not a proven global optimum; active bounds, singular values, and Jacobian rank are local diagnostics and can be affected by priors.

`save_fit` writes `sources.csv`, `predictions.csv`, and `fit.json`. Predictions retain observation order and measured, crustal, background, modeled, and residual components; residuals are modeled minus observed. Btot is the vector magnitude, not a fourth independent fit residual. `plot_fit` compares components and magnitudes along the sample sequence. `write_source_input` exports blocks on a configurable constant-altitude map grid for a separate forward run; the inversion uses the individual orbital altitudes.

The same observations can feed `fit_equivalent_sources`, whose fixed geometry, unweighted objective, and damping differ from the block optimizer; see [equivalent-source fitting](sphere.md#fit-an-equivalent-point-dipole-model). Saved moments support [GeoTIFF context maps](maps.md) and prediction at other altitudes. Generalizing this workflow requires explicit body/frame metadata, independent data-selection policies, and shared solver/result contracts as described below.

### 1. Explicit planetary and observation context

Introduce a validated body dictionary shared by ingestion, source construction, forward modeling, and mapping. It should identify the body-fixed frame, reference shape, size parameters, longitude convention, latitude definition, and depth/altitude convention. Named presets should be conveniences backed by recorded values; a body-independent entrypoint should not silently assume lunar radius or SEL coordinates.

Start with the supported spherical geometry. Extend shape adapters to ellipsoids, relief surfaces, and irregular meshes as corresponding source/evaluation backends become available. Changing a radius alone must not imply support for non-spherical geometry. Validate source burial and each backend's observation-domain restrictions against the selected reference surface.

Define a mission-neutral observation schema with positions, observable type, vector basis, units, timestamps, quality flags, uncertainty information, and provenance. Keep LPMAG selection rules in its adapter. Add mission/frame adapters that transform positions and vectors together; record any SPICE kernels or rotation models used. Magnetic vectors, scalar field magnitudes, and gravity observations need distinct observable definitions.

### 2. Common source and forward-operator contracts

Describe source geometry separately from coefficients, distinguishing dipole moments, volume magnetization, density contrast, and harmonic coefficients. Preserve units and normalization through every adapter. Use dictionaries for public configuration and arrays/tables for numerical payloads.

Add a backend registry with capability declarations: supported physics, source families, reference shapes, observables, precision, arbitrary-position evaluation, derivatives, and memory limits. Unsupported combinations should fail explicitly. A grid-only backend should not silently provide interpolated orbital predictions as though it evaluated the requested positions.

For models linear in their coefficients, provide forward and adjoint actions (`matvec` and `rmatvec`), with optional matrix assembly for small problems. Expose geometry derivatives separately. Begin with the Fortran dipole and volume kernels; admit spectral or diagnostic integrations after they satisfy the relevant contracts. Different source bases may share a prediction interface without sharing coefficient meanings.

### 3. Shared inversion and regularization layer

Make uncertainty weighting, background models, regularization, and termination diagnostics independent of the forward backend. For fixed geometry, a common objective can be expressed as:

```text
minimize over coefficients m and background parameters a:
    ||W (G m + H a - d)||² + lambda² ||L (m - m0)||²
```

Here `G` is the source operator, `H` the background design, `W` the uncertainty whitening operator, and `L` a declared damping or spatial-regularization operator. Record coefficient scaling and penalty normalization so regularization is interpretable across resolutions and bodies. Keep numerical stabilization distinct from scientific regularization.

Provide small-problem reference solves and matrix-free iterative methods for larger layers, with an explicit memory budget. Support component selection, nonuniform uncertainties, optional covariance models, and robust residual treatment where the method allows them. Scalar total-field magnitude is nonlinear in vector-source coefficients; it needs an appropriate nonlinear observable or an explicitly documented linearization, not the vector linear-fit assumption.

### 4. Geometry estimation above equivalent-source fitting

Build geometry search on the common coefficient-fitting layer. A first implementation can search depth, layer spacing, or source extent while solving linear coefficients at each trial geometry. This separates the geometric hypothesis from its best-fitting moments. Generalize block geometry and adaptive source meshes after their operators and derivatives are validated.

Expose source-count or resolution selection as model-selection choices rather than implicit changes during a fit. Distinguish prescribed remanent magnetization from induced responses, which need an explicit inducing-field model. Add background time dependence through declared nuisance models rather than treating all orbital signal as crustal structure.

### 5. Consistent results and reproducibility

Return a common solution dictionary containing body/frame metadata, source geometry, coefficient definitions, solver configuration, units, predictions, residuals, convergence information, and provenance. Prediction at another altitude should consume this saved solution without mission-specific assumptions.

Record backend/version and build fingerprints, input selection, reference-shape parameters, coordinate transformations, regularization, and any external-field treatment. Store large arrays separately from readable metadata. Keep physical evaluation independent from GeoTIFF rendering and output-directory choices.

### 6. Validation across bodies and solvers

Use independent point-source references, superposition, frame-rotation checks, unit-conversion checks, adjoint tests, and mesh/order convergence. Include multiple reference radii and observation-height/source-depth ratios; add polar, longitude-seam, and non-spherical cases as their capabilities are implemented.

Compare backends at the same positions and in the same frame, using consistent units and declared source/discretization assumptions. Do not equate independent component interpolation with a physically coupled magnetic source model. Assess inversion with synthetic recovery, held-out passes and altitude bands, noise/background sensitivity, and changes in depth, resolution, and regularization.

The first delivery milestone is a spherical-body dipole problem that accepts explicit planetary context and normalized observations, switches compatible fitting methods without changing ingestion or plotting, and predicts at held-out altitudes. Support for arbitrary planetary bodies then expands through explicit shape, frame, mission, and physics capabilities—not radius substitution alone. Equivalent sources remain a field representation; geological geometry requires additional constraints and uncertainty analysis.
