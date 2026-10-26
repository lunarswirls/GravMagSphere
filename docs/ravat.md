# Regional Ravat-style L2resL1 pipeline

This implementation adds equivalent magnetic monopoles to GravMagSphere's native
`orbital.f90` library and uses SciPy augmented least squares for the inversion.
It reproduces the published **method structure**, rather than claiming to recover
the authors' exact source locations, regularization schedules, or numerical results.
The existing vector-dipole ridge fitter is a different method and is not used.

References: Ravat, Purucker & Olsen (2020), JGR Planets,
[10.1029/2019JE006187](https://doi.org/10.1029/2019JE006187), sections 2, 3 and 5.1;
[data archive](https://uknowledge.uky.edu/ees_data/2/). Ravat et al. refer to
Olsen et al. (2017) for further details of the equivalent-monopole/IRLS formulation.

## What is implemented

A source amplitude q has units nT km². Its mathematical potential is q/|x-s|
and its field is q(x-s)/|x-s|³, with positions in km. These are equivalent
sources, not physical isolated magnetic charges. The native kernel returns XYZ
columns; each orbital endpoint is rotated into its own local r/theta/phi frame
before differencing. Theta points south, phi east, and r outward.

Let G map source amplitudes to finite along-track vector differences, A map them
to surface radial field Br, L whiten the endpoint-difference covariance, and
w_j be normalized positive surface quadrature weights. Stage 1 minimizes

    sum_i rho([L(G m2 - d)]_i) + alpha_l2² sum_j w_j (A m2)_j².

Stage 2 fixes m2, sets d_res = d - G m2, and minimizes

    sum_i rho([L(G m1 - d_res)]_i)
        + 2 alpha_l1² sum_j w_j sqrt((A m1)_j² + epsilon_nt²).

Thus the L1-like penalty is on **surface radial field**, not source amplitudes,
and in the second stage it applies to the **correction**, not the fixed baseline.
The factor 2 defines this implementation's alpha convention and gives model IRLS
weights 1/sqrt(Br²+epsilon²). This is a smooth approximation to L1; epsilon is
explicit, fixed per run, and should be included in a sensitivity study.

The data loss is Tukey biweight with c=4.5 on whitened residuals. Its IRLS weights
are (1-(z/c)²)² for |z|<c and zero otherwise. Every stage starts with an unweighted
solve, then updates robust weights. Backtracking prevents later accepted objective
increases; convergence requires stability of the objective, predicted data and
regularized surface field. Charge changes alone are not a reliable convergence
criterion when equivalent sources are nonunique. Histories report rank, effective
row count, objective, field/charge changes, and stalled steps. Tukey is nonconvex:
initialization and uncertainties matter, and a failure to retain data is explicit.

The returned total vector is B2+B1 and magnitude is **norm(B2+B1)**. Magnitudes
are never added as scalars. This avoids introducing the scalar/vector inconsistency
observed in the downloaded regional DAT, whose origin remains unconfirmed.

## Orbital input and selection

Use `load_lpmag_csv` output: planet-fixed XYZ positions in km, vector fields in nT,
positive independent XYZ-component uncertainties, and `table` metadata with `utc`,
`source_file`, and original `source_row`. The CLI CSV defaults are `X_SEL`, `Y_SEL`,
`Z_SEL`, `Bx_SEL`, `By_SEL`, `Bz_SEL`, and `utc`. UTC must be present and valid.
For other headers or per-component uncertainties, use the Python loader's column
mapping/sigma options. Supply optional per-sample `track_ids` to prevent pairing
across orbit/pass boundaries within one file.

The default pipeline forms **every valid consecutive five-second pair** (5±0.6 s),
never crossing file, track, original-row, or cadence gaps. It does not sort samples
or bridge rows removed by quality selection. Differences remain in nT, not nT/km
or nT/s. No polynomial detrending or arbitrary per-orbit baseline subtraction is
performed. Components at each endpoint use that endpoint's spherical basis.

Endpoint XYZ uncertainties are rotated with full 3×3 covariance. Consecutive
differences have cross covariance **minus the covariance of their shared endpoint**.
A block-tridiagonal Cholesky factorization and forward solve whiten complete
segments without allocating a full dense track covariance. Separate segments
restart the factorization. This correctly propagates assumed independent endpoint
measurement noise; it does not model serial instrumental/external errors. The
paper describes variance weights but not this exact full-covariance implementation.
The lower-level API also supports `pairing='disjoint'` for independent pairs and
the optional joint mean/difference data mode.

For real LP data, apply selection upstream using the paper's section 3.1 criteria:
low-altitude 1999 extended mission data, wake/tail selection with the stated solar
geometry, its polar exception beyond ±75°, and consistent body-fixed positions.
This pipeline does not download PDS data, derive Sun positions, perform instrument
calibration, or automatically recreate that selection. The CLI requires
`--selection-confirmed`, and applies a 12–48 km altitude filter. Its differences-only
fit has no vector-field observations or ad-hoc detrending. Generic loader wake
flags or a filename containing 'wake' do not establish the paper's selection.
Use a regional, padded subset of observations; the CLI does not automatically
crop files to the source bounds.

**The published surface DAT and COF files are model products, not orbital inputs.**
They can validate maps or support separately labelled synthetic experiments, but
cannot reproduce the original inversion's information content.

## Source geometry, tuning and limitations

`regional_monopole_grid` creates approximately equally spaced latitude rows with
longitude spacing adjusted by cos(latitude); it is **not** the paper's Leopardi
equal-area distribution. Sources lie radially below the 1738.1 × 1736.0 km lunar
ellipsoid. By default, regularization samples are directly above source locations,
with normalized solid-angle cell weights. This is a regional approximation to the
published surface average, not an exact integral over the whole Moon or exact
ellipsoid area. The API accepts independent surface locations/weights: rerun with
a denser quadrature to check sensitivity to unsampled features.

Section 5.1 mentions 7 km spacing / 5 km depth, while the Figure 5 caption specifies
8.7 km / 10 km. Both geometry choices are configurable; neither is asserted to be
the unique exact figure configuration. The CLI defaults to 15 km spacing / 10 km
depth for an inexpensive demonstration. Source bounds should extend beyond the
interpretation region to reduce edge effects. Sweep source geometry, damping and
epsilon, validate withheld complete passes, and inspect N–S artifacts and model
stability. No automatic tuning criterion is claimed to replace the paper's visual
choice of a stable downward-continued field.

Alphas are **not numerically transferable** from the paper: normalized quadrature,
source units, whitening, sample counts and smooth-L1 convention affect their scale;
stage-1 and stage-2 alphas also have different units. Real-data CLI runs require
both values explicitly. Model residuals alone cannot establish resolved surface
wavelengths or physically meaningful magnetization/depth.

This is a dense, memory-budgeted **regional** implementation. It is not a scalable
replica of the authors' 35,000–100,000-source global inversions or their 84-subset
merging. There is no global flux/zero-net-charge constraint, no global tiling,
no SH coefficient export, and no formal covariance calculation. Local equivalent
charges should not be interpreted or extrapolated as a global magnetic model.

## Run

From the repository root (requires gfortran, NumPy, SciPy, pandas; Matplotlib for plots):

```bash
.venv/bin/python examples/ravat_l2resl1/run.py --demo --plot
```

The demo uses stronger L1 damping to avoid large unsupported edge peaks; low data misfit alone is insufficient for downward continuation. This uses a reproducible **analytic dipole** truth at 30 km altitude, adds noise and
small pass offsets, and fits equivalent monopoles. It is a synthetic method check,
not a reconstruction of the published Reiner Gamma inversion. The supplied notebook
`examples/ravat_l2resl1_pipeline.ipynb` runs this same example and displays its stages.

For already selected real observations (damping values below are illustrative only):

```bash
.venv/bin/python examples/ravat_l2resl1/run.py \
  --observations /path/to/selected_regional_lp.csv --selection-confirmed \
  --bounds 295 0 310 15 --spacing-km 8.7 --depth-km 10 \
  --alpha-l2 0.2 --alpha-l1 2.0 --sigma-nt 0.1 \
  --iterations 200 --memory-mib 2048 --plot
```

The source bounds may need padding for a chosen interpretation window. At fine
source spacing, large files can exceed the dense memory budget; downselect to a
smaller padded region or coarsen the grid. Do not silently decimate differences to
hide resource limits or tune on withheld passes.

Outputs include reloadable `solution.npz`, `monopoles.csv`, endpoint pair indices
and predictions, the selected observation index, both iteration histories,
`fit.json`, `run.json` with input hashes/configuration, and ellipsoid surface maps
of baseline, correction and the vector sum. `success=False` means a stage did not
meet the convergence criteria, even if finite maps were saved. Check convergence
before using the output. To reload for new positions:

```python
import numpy as np
from gravmagpy import predict_ravat_field
with np.load('examples/ravat_l2resl1/output/solution.npz') as solution:
    prediction = predict_ravat_field(xyz_km, solution)
```

Validation: `PYTHONPATH=src .venv/bin/python -m unittest discover -s tests -p test_ravat.py -v`.

## MoonBroom filtered wake directory

The pipeline notebook now defaults to the user-selected directory:

    /Users/danywaller/Projects/moon/lpmag_l1b_5s_avg/lp_mag_nightside_filtering/filtered_wake

`load_moonbroom_wake` selects the `alt_lt_40km_brms_90pct` child, reads only dated
`ma*_nightside_wake_alt_lt_40km_brms_90pct.csv` files, and excludes aggregate,
35-km, 80-percentile and detrended variants to avoid duplicates or mixing products.
The time window is January 1 through July 28, 1999. Body-fixed Cartesian positions
are used to recompute location and radial altitude, and samples are restricted to
12–40 km above the 1737.4 km reference sphere and the configured geographic bounds.
The input products already carry their own altitude/Brms/wake cuts; those are
recorded as a MoonBroom selection, not asserted identical to the paper's selection.

The default Reiner Gamma window, 299–305°E and 4–11°N, retains 550 samples from
10 dated files and yields 443 valid consecutive difference pairs (244 if disjoint pairing is chosen). Sigma=1 nT is an explicit
assumed component uncertainty, not estimated from Brms. The default damping values
are an exploratory starting point, not a tuned published-model reproduction.
With a 1000-iteration cap, both stages converge for this dataset; the notebook
reports the actual histories and stores outputs under `output/filtered_wake`.

```bash
.venv/bin/python examples/ravat_l2resl1/run.py \
  --wake-dir /Users/danywaller/Projects/moon/lpmag_l1b_5s_avg/lp_mag_nightside_filtering/filtered_wake \
  --alpha-l2 0.2 --alpha-l1 2 --sigma-nt 1 --iterations 1000 --plot \
  --output examples/ravat_l2resl1/output/filtered_wake
```

Alternatively pass the child directory directly. This explicit wake-directory
route does not need `--selection-confirmed`: it records the actual known filter
policy in the run metadata instead of claiming a full paper-style selection.
Only retained input files are hashed; `input_report.csv` also lists the dates/files
with no regional samples. Input CSVs are read without modification.

## Joint vector and difference inversion

The optional `data_mode='vectors_and_differences'` remains available. The current
notebook, CLI and Python fit API use `data_mode='differences'` for the requested
gradient-based preprocessing. The CLI/notebook select every adjacent pair for that
mode and disjoint pairs for joint mode; the lower-level pairing API defaults to
disjoint for compatibility. CLI `--data-mode vectors` fits vector samples directly.
Joint outputs use `output/filtered_wake/vectors_and_differences`; the current
gradient-only outputs use `output/filtered_wake/differences`.

The joint mode fits **Br, Btheta and Bphi** at every orbital sample and their
along-track differences, using the same equivalent monopole amplitudes for every
component. This is an extension of the preferred differences-only configuration
in the paper. The L2/L1 model penalties remain on surface Br, as documented above;
adding vector data does not imply fitting three independent source models or
changing the model norm to a component-magnitude penalty.

Differences derived from endpoint measurements are not independent extra data.
Naively stacking B_left, B_right, and B_right-B_left would give a singular joint
covariance and would overcount data if diagonal errors were assumed. Instead,
each disjoint pair contributes the invertible six-coordinate basis

    mean = (B_left + B_right)/2
    delta = B_right - B_left

with full 6×6 covariance

    C_mean  = (C_left + C_right)/4
    C_delta = C_left + C_right
    C_mean,delta = (C_right - C_left)/2.

Each C is the full 3×3 spherical-frame covariance at its own endpoint. The complete
block matrix is Cholesky-whitened; the cross block is retained even when endpoint
uncertainties differ. Unpaired samples contribute their complete three-component
vectors with their own covariance. Pair means are component averages, not vectors
asserted to live at a new midpoint location; the design matrix uses the same
endpoint averaging operation.

For Gaussian least squares this preserves exactly the original vector information:
the whitened normal matrix, right-hand side and data quadratic norm agree with a
vectors-only formulation. The unit test verifies these identities with unequal,
anisotropic endpoint errors and unpaired data. With Tukey loss, component-wise
robust weights depend on the whitened basis: mean and difference channels can be
downweighted differently. This is explicit, rather than claiming that the derived
differences supply additional independent measurements. No per-pass external-field
background is currently estimated.

The selected 550 samples contribute 1650 independent whitened rows, including
244 mean/difference blocks and 62 unpaired vectors. `vectors.csv` records every
observed/predicted/residual component, `differences.csv` records the paired fields,
and `component_metrics.csv` gives separate vector and difference bias/RMSE for
Br, Btheta and Bphi. Surface-stage figures now display all three components and
magnitude, with the final magnitude still calculated after vector addition.

## Gradient preprocessing audit and remaining replication limits

`preprocessing.json` and the notebook show component conventions, pair selection,
finite-difference sign, cadence tolerance and uncertainty assumptions. The saved
`differences.csv` includes both UTCs, original rows, source files, elapsed seconds
and spatial separation for every pair, alongside observed and predicted component
differences. The notebook also plots the input components and their five-second
differences over one contiguous segment. All three differences enter the objective;
absolute-vector plots/metrics are diagnostic only in this mode.

Ravat et al. section 3.1 and Figure 1b specify five-second differences without the
polynomial/ad-hoc corrections used in the vector-data experiments. Section 3.1
also discusses fitted internal/external dipole removal following Purucker &
Nicholas (2010). Original fitted coefficients and windows are unavailable here:
the pipeline does not invent a regional polynomial/dipole replacement and records
this missing step explicitly. It also cannot restore discarded 40–48 km data,
Brms-rejected samples, or the paper's polar exception from prefiltered files.
The existing MoonBroom wake selection is not asserted equivalent to the published
wake/tail selection. Retaining those limits is deliberate provenance, not a claim
of exact reproduction of every preprocessing detail.
