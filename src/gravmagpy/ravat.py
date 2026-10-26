"""Regional Ravat-style L2 + residual-L1 equivalent-monopole inversion.

This reproduces the documented objective structure, not unpublished tuning or
source placement. Magnetic charges are mathematical equivalent sources, in
nT km^2. Differences are finite endpoint differences in local spherical frames,
not derivatives divided by distance/time. See docs/ravat.md for assumptions.
"""

import ctypes
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import lstsq

from .geometry import spherical_to_cartesian
from .data import load_lpmag_csv
from .utils import build_fortran


def load_moonbroom_wake(directory, *, bounds=(299, 4, 305, 11), sigma_nt=1.0,
                       start="1999-01-01", stop="1999-07-29"):
    """Load undetrended per-day filtered wake files; start inclusive, stop exclusive.

    Accept the filtered_wake root or its alt_lt_40km_brms_90pct child.
    Preserve file/row identities and reject cadence gaps later in pairing.
    Existing MoonBroom wake/Brms cuts are not asserted identical to Ravat's.
    """
    directory = Path(directory).expanduser().resolve()
    child = directory/"alt_lt_40km_brms_90pct"
    input_dir = child if child.is_dir() else directory
    pattern = "ma*_nightside_wake_alt_lt_40km_brms_90pct.csv"
    files = sorted(input_dir.glob(pattern))
    start_time, stop_time = pd.Timestamp(start, tz="UTC"), pd.Timestamp(stop, tz="UTC")
    if start_time >= stop_time:
        raise ValueError("start must precede exclusive stop")
    # These products use maYYMMDD names. Select intersecting days before IO.
    files = [p for p in files if start_time.floor('D') <=
             pd.to_datetime(p.name[2:8], format='%y%m%d', utc=True) < stop_time]
    if not files:
        raise FileNotFoundError(f"no dated per-day {pattern} files in {input_dir}")
    west, south, east, north = map(float, bounds)
    if east < west:
        east += 360
    low, high = (west+180)%360-180, (east+180)%360-180
    observations = load_lpmag_csv(files, wake_only=True, sigma_nt=sigma_nt,
        latitude_range=(south, north), longitude_range=(low, high), altitude_range_km=(12,40))
    dates = pd.to_datetime(observations['table'].utc, utc=True)
    mask = ((dates >= start_time) & (dates < stop_time)).to_numpy()
    for key in ('xyz_km', 'field_nt', 'sigma_nt'):
        observations[key] = observations[key][mask]
    observations['table'] = observations['table'].loc[mask].reset_index(drop=True)
    if not len(observations['table']):
        raise ValueError("no wake observations remain in the configured region/time interval")
    for report in observations['report']:
        report['retained_rows'] = int((observations['table'].source_file == report['path']).sum())
    observations['selection'] = dict(directory=str(input_dir), pattern=pattern,
        start_inclusive=start_time.isoformat(), stop_exclusive=stop_time.isoformat(),
        bounds=list(bounds), input_products="undetrended filtered wake: alt_lt_40km_brms_90pct",
        additional_radial_altitude_km=[12,40], radial_reference_km=1737.4,
        sigma_nt=sigma_nt, sigma_is_assumed=True,
        caveat="MoonBroom selection; equivalence to the paper's wake/tail/quality selection has not been established")
    return observations


def _xyz(values, name):
    values = np.ascontiguousarray(values, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] != 3 or not len(values) or not np.isfinite(values).all():
        raise ValueError(f"{name} must be finite nonempty (n,3) positions in km")
    if np.any(np.linalg.norm(values, axis=1) == 0):
        raise ValueError(f"{name} cannot include the origin")
    return values


def spherical_basis(xyz_km):
    """Rows of each rotation matrix are outward r, southward theta, eastward phi."""
    xyz = _xyz(xyz_km, "xyz_km")
    lon = np.arctan2(xyz[:, 1], xyz[:, 0])
    lat = np.arctan2(xyz[:, 2], np.hypot(xyz[:, 0], xyz[:, 1]))
    cl, sl, cp, sp = np.cos(lat), np.sin(lat), np.cos(lon), np.sin(lon)
    return np.stack((np.column_stack((cl*cp, cl*sp, sl)),
                     np.column_stack((sl*cp, sl*sp, -cl)),
                     np.column_stack((-sp, cp, np.zeros_like(sp)))), axis=1)


@lru_cache(maxsize=8)
def _backend(library):
    backend = ctypes.CDLL(str(library))
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    backend.monopole_kernel.argtypes = [ctypes.c_int, ctypes.c_int, array, array, array,
                                        ctypes.POINTER(ctypes.c_int)]
    backend.monopole_kernel.restype = None
    return backend


def monopole_kernel(xyz_km, source_xyz_km, *, library=None, max_memory_mib=512):
    """GravMagSphere native XYZ field kernel, shape (observations,3,sources)."""
    xyz, sources = _xyz(xyz_km, "xyz_km"), _xyz(source_xyz_km, "source_xyz_km")
    if not np.isfinite(max_memory_mib) or max_memory_mib <= 0:
        raise ValueError("max_memory_mib must be positive")
    if max(len(xyz), len(sources)) > np.iinfo(np.int32).max:
        raise ValueError("dimensions exceed the native integer range")
    if xyz.size*len(sources)*8 > max_memory_mib*1024**2:
        raise ValueError("monopole kernel exceeds memory limit; reduce region or batch predictions")
    kernel = np.empty((len(xyz), 3, len(sources)))
    status = ctypes.c_int()
    _backend(library or build_fortran("orbital")).monopole_kernel(
        len(xyz), len(sources), xyz, sources, kernel, ctypes.byref(status))
    if status.value or not np.isfinite(kernel).all():
        raise ValueError(f"monopole evaluation failed (status {status.value}); check source/observer separation")
    return kernel


def ellipsoid_points(lat_deg, lon_deg, *, a_km=1738.1, b_km=1736.0, altitude_km=0):
    """Geocentric latitude; altitude is radial distance above the ellipsoid."""
    if not np.isfinite([a_km, b_km, altitude_km]).all() or not a_km >= b_km > 0 or altitude_km < 0:
        raise ValueError("need a>=b>0 and nonnegative radial altitude")
    lat, lon = np.broadcast_arrays(np.asarray(lat_deg, float), np.asarray(lon_deg, float))
    if not np.isfinite([lat, lon]).all() or np.any(np.abs(lat) > 90):
        raise ValueError("invalid latitude/longitude")
    angle = np.deg2rad(lat)
    radius = a_km*b_km/np.sqrt((b_km*np.cos(angle))**2+(a_km*np.sin(angle))**2)
    return spherical_to_cartesian(lat, lon, radius+altitude_km)


def regional_monopole_grid(bounds, *, spacing_km=8.7, depth_km=10, a_km=1738.1, b_km=1736.0):
    """Approximately equal-spacing stagger-free latitude rows, not Leopardi tiling.

    bounds=(west,south,east,north), continuous longitude with span <=180 degrees.
    Source depth is radial below the ellipsoid. Returned surface rows have
    solid-angle cell weights (not exact ellipsoid-area weights).
    """
    bounds = np.asarray(bounds, float)
    if bounds.shape != (4,) or not np.isfinite(bounds).all():
        raise ValueError("bounds must be finite (west,south,east,north)")
    west, south, east, north = bounds
    if east < west:
        east += 360
    if not (0 < east-west <= 180 and -90 < south < north < 90):
        raise ValueError("regional bounds must avoid poles and span at most 180 degrees")
    if not np.isfinite([spacing_km, depth_km]).all() or spacing_km <= 0 or not 0 < depth_km < b_km:
        raise ValueError("need positive spacing and source depth less than polar radius")
    step = np.rad2deg(spacing_km/a_km)
    nlat = max(1, int(np.ceil((north-south)/step)))
    lat_edges = np.linspace(south, north, nlat+1)
    rows = []
    for low, high in zip(lat_edges[:-1], lat_edges[1:]):
        lat = (low+high)/2
        nlon = max(1, int(np.ceil((east-west)*np.cos(np.deg2rad(lat))/step)))
        weight = (np.sin(np.deg2rad(high))-np.sin(np.deg2rad(low)))*np.deg2rad((east-west)/nlon)
        for lon in west+(np.arange(nlon)+0.5)*(east-west)/nlon:
            rows.append((lat, lon, weight))
    table = pd.DataFrame(rows, columns=["lat_deg", "lon_deg", "area_weight"])
    surface = ellipsoid_points(table.lat_deg, table.lon_deg, a_km=a_km, b_km=b_km)
    radius = np.linalg.norm(surface, axis=1)
    sources = surface*((radius-depth_km)/radius)[:, None]
    table["radius_km"] = radius-depth_km
    return {"source_xyz_km": sources, "surface_xyz_km": surface,
            "area_weights": table.area_weight.to_numpy(), "table": table}


def along_track_differences(observations, *, cadence_s=5, tolerance_s=0.6, pairing="disjoint"):
    """Build five-second endpoint differences, rejecting time/file/row gaps.

    Requires load_lpmag_csv metadata: utc, source_file, source_row. Optional
    observations['track_ids'] supplies per-sample pass identities. Pairing is
    pairing='adjacent' uses every valid consecutive pair, as in a sampled track.
    Shared endpoint covariance is then handled in the differences-only solve.
    pairing='disjoint' retains independent pairs and supports joint mean/delta
    fitting. XYZ component errors are rotated in full. No detrending is applied.
    """
    if observations.get("frame", "planet_fixed") not in ("planet_fixed", "SEL"):
        raise ValueError("differences require body-fixed Cartesian positions and vectors")
    if pairing not in ('disjoint','adjacent'):
        raise ValueError("pairing must be disjoint or adjacent")
    xyz = _xyz(observations["xyz_km"], "xyz_km")
    field = np.asarray(observations["field_nt"], float)
    sigma = np.broadcast_to(np.asarray(observations.get("sigma_nt", 1), float), xyz.shape)
    if field.shape != xyz.shape or not np.isfinite(field).all() or not np.isfinite(sigma).all() or np.any(sigma <= 0):
        raise ValueError("field and positive XYZ uncertainties must match positions")
    if not np.isfinite([cadence_s, tolerance_s]).all() or cadence_s <= 0 or not 0 <= tolerance_s < cadence_s:
        raise ValueError("invalid difference cadence/tolerance")
    table = pd.DataFrame(observations["table"])
    if len(table) != len(xyz) or not {"utc", "source_file", "source_row"}.issubset(table):
        raise ValueError("time/file/original-row metadata are required")
    times = pd.to_datetime(table.utc, utc=True, errors="raise")
    if times.isna().any():
        raise ValueError("utc must be present for every sample")
    seconds = times.array.as_unit("ns").asi8/1e9
    original_rows = pd.to_numeric(table.source_row, errors="raise").to_numpy()
    if not np.isfinite(original_rows).all() or np.any(original_rows != np.floor(original_rows)):
        raise ValueError("source_row must contain original integer row identifiers")
    tracks = np.asarray(observations.get("track_ids", table.source_file)).astype(str)
    if tracks.shape != (len(xyz),):
        raise ValueError("track_ids must identify each sample")
    pairs, i = [], 0
    while i < len(xyz)-1:
        j = i+1
        dt = seconds[j]-seconds[i]
        if (tracks[i] == tracks[j] and table.source_file.iloc[i] == table.source_file.iloc[j]
                and original_rows[j]-original_rows[i] == 1 and abs(dt-cadence_s) <= tolerance_s):
            pairs.append((i, j))
            i += 2 if pairing == 'disjoint' else 1
        else:
            i += 1
    if not pairs:
        raise ValueError("no contiguous pairs match the requested cadence")
    pairs = np.asarray(pairs, dtype=int)
    basis = spherical_basis(xyz)
    vector = np.einsum('nci,ni->nc', basis, field)
    covariance = np.einsum('nci,ni,ndi->ncd', basis, sigma**2, basis)
    left, right = pairs.T
    diff_covariance = covariance[left]+covariance[right]
    pair_table = pd.DataFrame(dict(left_index=left,right_index=right,
        source_file=table.source_file.to_numpy()[left],
        left_source_row=original_rows[left],right_source_row=original_rows[right],
        left_utc=times.astype(str).to_numpy()[left],right_utc=times.astype(str).to_numpy()[right],
        dt_s=seconds[right]-seconds[left],
        separation_km=np.linalg.norm(xyz[right]-xyz[left],axis=1)))
    return {"xyz_km": xyz, "pairs": pairs,
            "vector_nt": vector, "vector_covariance_nt2": covariance,
            "difference_nt": vector[right]-vector[left],
            "covariance_nt2": diff_covariance, "pair_table": pair_table,
            "preprocessing": dict(pairing=pairing,cadence_s=cadence_s,tolerance_s=tolerance_s,
                components=['Br','Btheta','Bphi'], difference_sign='later minus earlier',
                units='nT per finite five-second interval; not divided by separation or time',
                detrending='none',ad_hoc_pass_rejection='none beyond existing input filters',
                uncertainties='assumed endpoint XYZ sigma, full rotation and shared-endpoint propagation',
                dipole_term_removal='not applied: original fitted coefficients/windows unavailable'),
            "whitener": np.linalg.inv(np.linalg.cholesky(diff_covariance))}


def _difference_data_system(spherical, differences):
    """Whiten consecutive differences using block-tridiagonal covariance.

    Cov(delta_i,delta_(i-1))=-C_shared. Block Cholesky and forward substitution
    avoid allocating the complete track covariance. Separate segments restart.
    """
    pairs = np.asarray(differences['pairs'])
    left,right = pairs.T
    raw_g = spherical[right]-spherical[left]
    raw_d = np.asarray(differences['difference_nt'])
    if len(np.unique(pairs)) == pairs.size:
        w = differences['whitener']
        return (np.einsum('nci,nis->ncs',w,raw_g).reshape(-1,spherical.shape[2]),
                np.einsum('nci,ni->nc',w,raw_d).ravel())
    cov = np.asarray(differences.get('vector_covariance_nt2'),float)
    if cov.shape != (len(spherical),3,3) or not np.isfinite(cov).all():
        raise ValueError('overlapping differences require original endpoint covariance')
    output_g,output_d = np.empty_like(raw_g),np.empty_like(raw_d)
    used = set()
    previous_cholesky = None
    for k,(a,b) in enumerate(pairs):
        continues = k > 0 and a == pairs[k-1,1]
        if b in used or (a in used and not continues):
            raise ValueError('overlapping pairs must form ordered, nonbranching tracks')
        used.update((a,b))
        diagonal = cov[a]+cov[b]
        if continues:
            subdiagonal = np.linalg.solve(previous_cholesky,(-cov[a]).T).T
            diagonal = diagonal-subdiagonal@subdiagonal.T
            rhs_g = raw_g[k]-subdiagonal@output_g[k-1]
            rhs_d = raw_d[k]-subdiagonal@output_d[k-1]
        else:
            rhs_g,rhs_d = raw_g[k],raw_d[k]
        previous_cholesky = np.linalg.cholesky(diagonal)
        output_g[k] = np.linalg.solve(previous_cholesky,rhs_g)
        output_d[k] = np.linalg.solve(previous_cholesky,rhs_d)
    return output_g.reshape(-1,spherical.shape[2]),output_d.ravel()


def _vector_data_system(spherical, differences, mode):
    """Whiten vectors, or pair means + differences, without duplicate information.

    [B_left,B_right,Delta B] has a singular covariance since Delta B is derived.
    The invertible [mean,Delta B] basis retains all six independent measurements.
    Full mean/difference cross covariance is retained for unequal endpoint errors.
    """
    vector = np.asarray(differences.get("vector_nt"), float)
    covariance = np.asarray(differences.get("vector_covariance_nt2"), float)
    n, _, ns = spherical.shape
    if vector.shape != (n,3) or covariance.shape != (n,3,3) or not np.isfinite(vector).all() or not np.isfinite(covariance).all():
        raise ValueError("vector modes require original spherical vectors and endpoint covariance")
    if not np.allclose(covariance, covariance.transpose(0,2,1), rtol=1e-12, atol=1e-12):
        raise ValueError("endpoint covariance must be symmetric positive definite")
    # Validate positive definiteness also for endpoints entering a paired block.
    np.linalg.cholesky(covariance)
    if mode == "vectors":
        whitener = np.linalg.inv(np.linalg.cholesky(covariance))
        return (np.einsum('nci,nis->ncs',whitener,spherical).reshape(-1,ns),
                np.einsum('nci,ni->nc',whitener,vector).ravel())
    left, right = np.asarray(differences['pairs']).T
    if len(np.unique(differences['pairs'])) != np.asarray(differences['pairs']).size:
        raise ValueError('joint mean/difference fitting requires disjoint pairs; choose pairing="disjoint"')
    if not np.allclose(vector[right]-vector[left], differences['difference_nt'], atol=1e-10, rtol=1e-10):
        raise ValueError("joint mode requires differences derived from the same vector measurements")
    total = covariance[left]+covariance[right]
    cross = (covariance[right]-covariance[left])/2
    joint_covariance = np.empty((len(left),6,6))
    joint_covariance[:,:3,:3],joint_covariance[:,3:,3:] = total/4,total
    joint_covariance[:,:3,3:],joint_covariance[:,3:,:3] = cross,cross.transpose(0,2,1)
    whitener = np.linalg.inv(np.linalg.cholesky(joint_covariance))
    raw_g = np.concatenate(((spherical[left]+spherical[right])/2,
                             spherical[right]-spherical[left]),axis=1)
    raw_d = np.concatenate(((vector[left]+vector[right])/2, vector[right]-vector[left]),axis=1)
    matrices = [np.einsum('nci,nis->ncs',whitener,raw_g).reshape(-1,ns)]
    targets = [np.einsum('nci,ni->nc',whitener,raw_d).ravel()]
    unpaired = np.setdiff1d(np.arange(n), np.asarray(differences['pairs']).ravel())
    if len(unpaired):
        single_w = np.linalg.inv(np.linalg.cholesky(covariance[unpaired]))
        matrices.append(np.einsum('nci,nis->ncs',single_w,spherical[unpaired]).reshape(-1,ns))
        targets.append(np.einsum('nci,ni->nc',single_w,vector[unpaired]).ravel())
    return np.vstack(matrices), np.concatenate(targets)


def _irls(g, data, radial, area, *, alpha, l1, epsilon_nt, max_iterations, tolerance, tukey_c):
    """Small/regional dense augmented least squares; no squared normal equations."""
    def objective(m):
        z = g@m-data
        az = np.abs(z)
        if tukey_c is None:
            misfit = np.sum(z*z)
        else:
            u = np.minimum(az/tukey_c, 1)
            misfit = np.sum(tukey_c**2/3*(1-(1-u*u)**3))
        br = radial@m
        penalty = np.sum(area*(2*np.sqrt(br*br+epsilon_nt**2) if l1 else br*br))
        return float(misfit+alpha**2*penalty)

    m = np.zeros(g.shape[1])
    history = []
    converged = False
    for iteration in range(max_iterations):
        residual = g@m-data
        weights = np.ones_like(data)
        # First solve is unweighted: otherwise a zero model can reject the signal.
        if tukey_c is not None and iteration:
            u = residual/tukey_c
            weights = np.where(np.abs(u) < 1, (1-u*u)**2, 0)
            if np.count_nonzero(weights) < 3:
                raise ValueError("Tukey weighting rejected nearly all data; check uncertainties/initial damping")
        model_weights = area.copy()
        if l1:
            model_weights /= np.sqrt((radial@m)**2+epsilon_nt**2)
        data_scale = np.sqrt(weights)
        penalty = alpha*np.sqrt(model_weights)[:, None]*radial
        augmented = np.vstack((g*data_scale[:, None], penalty))
        rhs = np.r_[data*data_scale, np.zeros(len(radial))]
        candidate, _, rank, _ = lstsq(augmented, rhs, lapack_driver="gelsy")
        value_before = objective(m)
        # Backtracking handles nonconvex robust loss without accepting an increase.
        step = 1.0
        full_candidate = candidate.copy()
        while iteration and objective(candidate) > value_before+1e-10*max(1, value_before) and step > 2**-20:
            step *= 0.5
            candidate = m+step*(full_candidate-m)
        if iteration and objective(candidate) > value_before+1e-10*max(1, value_before):
            candidate = m.copy()
            stalled = True
        else:
            stalled = False
        change = np.linalg.norm(candidate-m)/max(1, np.linalg.norm(m))
        field_change = max(np.linalg.norm(g@(candidate-m))/max(1, np.linalg.norm(g@m)),
                           np.linalg.norm(np.sqrt(area)*(radial@(candidate-m)))/
                           max(1, np.linalg.norm(np.sqrt(area)*(radial@m))))
        objective_change = abs(objective(candidate)-value_before)/max(1, abs(value_before))
        m = candidate
        history.append({"iteration": iteration+1, "objective": objective(m), "relative_change": float(change),
                        "relative_field_change": float(field_change), "relative_objective_change": float(objective_change),
                        "active_rows": int(np.count_nonzero(weights)), "rank": int(rank), "stalled": stalled})
        # Equivalent charges can be nonunique. Convergence is defined in data
        # and regularization-field space, together with objective stability.
        if iteration > 0 and objective_change < tolerance and field_change < tolerance and not stalled:
            converged = True
            break
    return m, pd.DataFrame(history), converged


def fit_ravat_l2resl1(differences, source_xyz_km, surface_xyz_km, *, alpha_l2, alpha_l1,
                     area_weights=None, epsilon_nt=0.1, max_iterations=40, tolerance=1e-5,
                     tukey_c=4.5, max_memory_mib=512, library=None, data_mode="differences"):
    """Robust L2 surface-Br baseline plus residual-fitted L1 surface-Br correction.

    Minimize robust whitened data misfit plus alpha_l2^2 mean(Br_baseline^2),
    then fit d-G*m_baseline with alpha_l1^2 * 2*mean(sqrt(Br_correction^2+eps^2)).
    The factor 2 yields IRLS weights 1/sqrt(Br^2+eps^2). Alpha values therefore
    have different units between stages and are not the paper's tuning values.
    No L1 penalty on source amplitudes; the correction's surface Br is penalized.
    data_mode='vectors_and_differences' fits pair means and differences with full
    covariance, plus unpaired vectors; this does not double-count measurements.
    For Gaussian loss it carries the same information as data_mode='vectors'.
    Tukey weights operate in the selected whitened basis and are not invariant
    under changing that basis. Surface-Br regularization is unchanged.
    """
    parameters = [alpha_l2, alpha_l1, epsilon_nt, tolerance, max_memory_mib]
    if not np.isfinite(parameters).all() or np.any(np.asarray(parameters) <= 0):
        raise ValueError("alphas, epsilon, tolerance and memory limit must be positive")
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations < 2:
        raise ValueError("max_iterations must be an integer >=2")
    if tukey_c is not None and (not np.isfinite(tukey_c) or tukey_c <= 0):
        raise ValueError("Tukey cutoff must be positive or None")
    if data_mode not in ("differences", "vectors", "vectors_and_differences"):
        raise ValueError("invalid data_mode")
    xyz = _xyz(differences["xyz_km"], "xyz_km")
    sources = _xyz(source_xyz_km, "source_xyz_km")
    surface = _xyz(surface_xyz_km, "surface_xyz_km")
    pairs = np.asarray(differences["pairs"])
    if pairs.ndim != 2 or pairs.shape[1] != 2 or not len(pairs) or not np.issubdtype(pairs.dtype, np.integer):
        raise ValueError("pairs must be a nonempty integer (n,2) array")
    if pairs.min() < 0 or pairs.max() >= len(xyz) or np.any(pairs[:,0] >= pairs[:,1]):
        raise ValueError("pair endpoints must be valid and ordered earlier/later")
    data = np.asarray(differences["difference_nt"], float)
    whiten = np.asarray(differences["whitener"], float)
    if data.shape != (len(pairs), 3) or whiten.shape != (len(pairs),3,3) or not np.isfinite(data).all() or not np.isfinite(whiten).all():
        raise ValueError("invalid differences or covariance whitener")
    if np.any(np.linalg.norm(sources, axis=1) >= np.linalg.norm(surface, axis=1).min()):
        raise ValueError("all sources must lie below all regularization surface radii")
    if np.any(np.linalg.norm(xyz, axis=1) <= np.linalg.norm(sources, axis=1).max()):
        raise ValueError("orbital samples must lie outside the source shell")
    # Conservative peak budget includes kernel rotations, differences, augmented
    # copies and factorization workspace. This implementation is intentionally regional.
    estimated = 8*(12*len(xyz)*len(sources)+12*len(surface)*len(sources)+
                   12*len(pairs)*len(sources)+8*len(sources)**2)/1024**2
    if estimated > max_memory_mib:
        raise ValueError(f"dense regional solve exceeds memory budget: approximately {estimated:.3f} MiB; reduce source count/region")
    area = np.ones(len(surface)) if area_weights is None else np.asarray(area_weights, float)
    if area.shape != (len(surface),) or not np.isfinite(area).all() or np.any(area <= 0):
        raise ValueError("surface area weights must be finite and positive")
    area = area/area.sum()
    library = library or build_fortran("orbital")
    kernel = monopole_kernel(xyz, sources, library=library, max_memory_mib=max_memory_mib)
    spherical = np.einsum('nci,nis->ncs', spherical_basis(xyz), kernel)
    left, right = pairs.T
    g_raw = spherical[right]-spherical[left]
    del kernel
    if data_mode == "differences":
        g, target = _difference_data_system(spherical,differences)
    else:
        g, target = _vector_data_system(spherical, differences, data_mode)
    surface_kernel = monopole_kernel(surface, sources, library=library, max_memory_mib=max_memory_mib)
    radial = np.einsum('ni,nis->ns', spherical_basis(surface)[:,0], surface_kernel)
    del surface_kernel
    options = dict(epsilon_nt=epsilon_nt, max_iterations=max_iterations, tolerance=tolerance, tukey_c=tukey_c)
    base, hist2, ok2 = _irls(g, target, radial, area, alpha=alpha_l2, l1=False, **options)
    correction, hist1, ok1 = _irls(g, target-g@base, radial, area, alpha=alpha_l1, l1=True, **options)
    charges = base+correction
    prediction = np.einsum('ncs,s->nc', g_raw, charges)
    result = {"source_xyz_km": sources, "baseline_charges_nt_km2": base,
            "correction_charges_nt_km2": correction, "charges_nt_km2": charges,
            "predicted_difference_nt": prediction, "residual_difference_nt": prediction-data,
            "baseline_difference_nt": np.einsum('ncs,s->nc', g_raw, base),
            "surface_br_baseline_nt": radial@base, "surface_br_correction_nt": radial@correction,
            "history_l2": hist2, "history_l1": hist1, "success": bool(ok2 and ok1),
            "converged_l2": ok2, "converged_l1": ok1,
            "alpha_l2": float(alpha_l2), "alpha_l1": float(alpha_l1), "epsilon_nt": float(epsilon_nt),
            "tukey_c": tukey_c, "estimated_memory_mib": float(estimated),
            "data_mode": data_mode, "independent_data_rows": int(len(target)),
            "predicted_vector_nt": np.einsum('ncs,s->nc',spherical,charges)}
    if 'vector_nt' in differences:
        result['residual_vector_nt'] = result['predicted_vector_nt']-differences['vector_nt']
    return result


def predict_ravat_field(xyz_km, solution, *, part="total", batch_size=256, library=None):
    """Predict XYZ and local spherical vectors; magnitude after vector addition."""
    xyz = _xyz(xyz_km, "xyz_km")
    key = {"total": "charges_nt_km2", "baseline": "baseline_charges_nt_km2",
           "correction": "correction_charges_nt_km2"}.get(part)
    if key is None:
        raise ValueError("part must be total, baseline, or correction")
    if not isinstance(batch_size, (int, np.integer)) or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    sources = _xyz(solution["source_xyz_km"], "source_xyz_km")
    charges = np.asarray(solution[key], float)
    if charges.shape != (len(sources),) or not np.isfinite(charges).all():
        raise ValueError("invalid equivalent charges")
    library = library or build_fortran("orbital")
    field = np.empty_like(xyz)
    for start in range(0, len(xyz), batch_size):
        sl = slice(start, start+batch_size)
        field[sl] = np.einsum('ncs,s->nc', monopole_kernel(xyz[sl], sources, library=library), charges)
    spherical = np.einsum('nci,ni->nc', spherical_basis(xyz), field)
    return {"xyz_nt": field, "rtp_nt": spherical, "total_nt": np.linalg.norm(field, axis=1)}


def save_ravat_fit(solution, differences, directory):
    """Persist reloadable equivalent charges, endpoint diagnostics and IRLS history."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(solution["source_xyz_km"], columns=["x_km", "y_km", "z_km"])
    for name in ("baseline_charges_nt_km2", "correction_charges_nt_km2", "charges_nt_km2"):
        table[name] = solution[name]
    table.to_csv(directory/"monopoles.csv", index=False)
    np.savez_compressed(directory/"solution.npz", **{key: solution[key] for key in
        ("source_xyz_km", "baseline_charges_nt_km2", "correction_charges_nt_km2", "charges_nt_km2")})
    pairs = differences.get('pair_table',pd.DataFrame(differences["pairs"], columns=["left_index", "right_index"])).copy()
    for j, component in enumerate(("r", "theta", "phi")):
        pairs[f"observed_delta_b{component}_nt"] = differences["difference_nt"][:,j]
        pairs[f"predicted_delta_b{component}_nt"] = solution["predicted_difference_nt"][:,j]
    pairs.to_csv(directory/"differences.csv", index=False)
    (directory/'preprocessing.json').write_text(json.dumps(differences.get('preprocessing',{}),indent=2)+'\n')
    metrics = []
    if 'vector_nt' in differences:
        vectors = pd.DataFrame(differences['xyz_km'],columns=['x_km','y_km','z_km'])
        for j, component in enumerate(('r','theta','phi')):
            vectors[f'observed_b{component}_nt'] = differences['vector_nt'][:,j]
            vectors[f'predicted_b{component}_nt'] = solution['predicted_vector_nt'][:,j]
            vectors[f'residual_b{component}_nt'] = solution['residual_vector_nt'][:,j]
        vectors.to_csv(directory/'vectors.csv',index=False)
    for kind, residual_key in [('vector','residual_vector_nt'),('difference','residual_difference_nt')]:
        if residual_key not in solution:
            continue
        for j, component in enumerate(('Br','Btheta','Bphi')):
            residual = solution[residual_key][:,j]
            metrics.append(dict(data=kind, component=component,bias_nt=float(np.mean(residual)),
                                rmse_nt=float(np.sqrt(np.mean(residual**2)))))
    pd.DataFrame(metrics).to_csv(directory/'component_metrics.csv',index=False)
    for suffix in ("l2", "l1"):
        solution[f"history_{suffix}"].to_csv(directory/f"history_{suffix}.csv", index=False)
    names = ("success", "converged_l2", "converged_l1", "alpha_l2", "alpha_l1", "epsilon_nt", "tukey_c", "estimated_memory_mib",
             "data_mode", "independent_data_rows")
    (directory/"fit.json").write_text(json.dumps({key: solution[key] for key in names}, indent=2)+"\n")
