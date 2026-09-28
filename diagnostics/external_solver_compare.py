#!/usr/bin/env python3
"""
Format GravMag Sphere lunar_examples for external solvers and compare predictions

Subcommands:
  - repos:    write a markdown list of candidate external solver repositories
  - format:   run point_source baseline and export solver-friendly CSV/JSON files
  - compare:  compare baseline against available solver backends
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import warnings
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.special import lpmv


repo_root = Path(__file__).resolve().parent.parent
example_root = repo_root / "examples"
COMPONENTS = ("br", "btheta", "bphi", "btot")
XYZ_COMPONENTS = ("bx", "by", "bz", "btot")
METRIC_STATS = ("rmse", "mae", "max_abs", "r2")
RE_GAUSS_META = re.compile(
    r"lmax=\s*(?P<lmax>\d+)\s+reg_lambda=\s*(?P<reg_lambda>[-+0-9.eE]+)\s+reg_power=\s*(?P<reg_power>[-+0-9.eE]+)"
)


def run_cmd(cmd: Sequence[str], cwd: Path, env: Dict[str, str] | None = None) -> str:
    proc = subprocess.run(
        list(cmd),
        cwd=str(cwd),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if proc.returncode != 0:
        joined = " ".join(cmd)
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {joined}\n"
            f"--- begin output ---\n{proc.stdout}\n--- end output ---"
        )
    return proc.stdout


def parse_list_ints(raw: str) -> List[int]:
    out: List[int] = []
    for part in raw.split(","):
        s = part.strip()
        if not s:
            continue
        out.append(int(s))
    return out


def parse_list_floats(raw: str) -> List[float]:
    out: List[float] = []
    for part in raw.split(","):
        s = part.strip()
        if not s:
            continue
        out.append(float(s))
    return out


def noncomment_lines(path: Path) -> List[str]:
    lines: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#") or s.startswith("!"):
                continue
            lines.append(s)
    return lines


def read_example_metadata(path: Path, rsphere_km: float) -> Dict[str, object]:
    lines = noncomment_lines(path)
    if len(lines) < 6:
        raise RuntimeError(f"Could not parse first body cards from {path}")
    title = lines[0]
    c2 = lines[1].split()
    c5 = lines[4].split()
    if len(c2) < 7:
        raise RuntimeError(f"Card 2 parse failed for {path}")
    if len(c5) < 5:
        raise RuntimeError(f"Card 5 parse failed for {path}")
    elvo_km = float(c2[4])
    ifield = int(float(c5[0]))
    return {
        "title": title,
        "mode": "magnetic" if ifield == 2 else "gravity",
        "lat0_deg": float(c2[0]),
        "lon0_deg": float(c2[1]),
        "dlat_deg": float(c2[2]),
        "dlon_deg": float(c2[3]),
        "elvo_km": elvo_km,
        "nlat": int(float(c2[5])),
        "nlon": int(float(c2[6])),
        "radius_m": (rsphere_km + elvo_km) * 1000.0,
    }


def parse_gauss_header_params(xyz_path: Path) -> Dict[str, float]:
    params: Dict[str, float] = {}
    if not xyz_path.exists():
        return params
    with xyz_path.open("r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s.startswith("#"):
                break
            m = RE_GAUSS_META.search(s)
            if not m:
                continue
            params["lmax_used"] = float(m.group("lmax"))
            params["reg_lambda_used"] = float(m.group("reg_lambda"))
            params["reg_power_used"] = float(m.group("reg_power"))
            break
    return params


def load_brtp_aggregate(path: Path) -> Dict[str, np.ndarray]:
    arr = np.loadtxt(path, comments="#")
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] < 7:
        raise RuntimeError(f"Unexpected BRTP format in {path}")

    lon = np.round(arr[:, 1], 6)
    lat = np.round(arr[:, 2], 6)
    br = arr[:, 3]
    bt = arr[:, 4]
    bp = arr[:, 5]

    accum: Dict[Tuple[float, float], List[float]] = {}
    for lo, la, bri, bti, bpi in zip(lon, lat, br, bt, bp):
        key = (float(lo), float(la))
        if key not in accum:
            accum[key] = [0.0, 0.0, 0.0]
        accum[key][0] += float(bri)
        accum[key][1] += float(bti)
        accum[key][2] += float(bpi)

    keys = sorted(accum.keys())
    lon_out = np.array([k[0] for k in keys], dtype=float)
    lat_out = np.array([k[1] for k in keys], dtype=float)
    br_out = np.array([accum[k][0] for k in keys], dtype=float)
    bt_out = np.array([accum[k][1] for k in keys], dtype=float)
    bp_out = np.array([accum[k][2] for k in keys], dtype=float)
    btot_out = np.sqrt(br_out * br_out + bt_out * bt_out + bp_out * bp_out)
    return {
        "lon_deg": lon_out,
        "lat_deg": lat_out,
        "br": br_out,
        "btheta": bt_out,
        "bphi": bp_out,
        "btot": btot_out,
    }


def load_xyz_aggregate(path: Path) -> Dict[str, np.ndarray]:
    arr = np.loadtxt(path, comments="#")
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] < 7:
        raise RuntimeError(f"Unexpected XYZ format in {path}")

    lon = np.round(arr[:, 1], 6)
    lat = np.round(arr[:, 2], 6)
    bx = arr[:, 3]
    by = arr[:, 4]
    bz = arr[:, 5]

    accum: Dict[Tuple[float, float], List[float]] = {}
    for lo, la, bxi, byi, bzi in zip(lon, lat, bx, by, bz):
        key = (float(lo), float(la))
        if key not in accum:
            accum[key] = [0.0, 0.0, 0.0]
        accum[key][0] += float(bxi)
        accum[key][1] += float(byi)
        accum[key][2] += float(bzi)

    keys = sorted(accum.keys())
    lon_out = np.array([k[0] for k in keys], dtype=float)
    lat_out = np.array([k[1] for k in keys], dtype=float)
    bx_out = np.array([accum[k][0] for k in keys], dtype=float)
    by_out = np.array([accum[k][1] for k in keys], dtype=float)
    bz_out = np.array([accum[k][2] for k in keys], dtype=float)
    btot_out = np.sqrt(bx_out * bx_out + by_out * by_out + bz_out * bz_out)
    return {
        "lon_deg": lon_out,
        "lat_deg": lat_out,
        "bx": bx_out,
        "by": by_out,
        "bz": bz_out,
        "btot": btot_out,
    }


def load_observations_csv(path: Path) -> Dict[str, np.ndarray]:
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return {
        "lon_deg": arr[:, 0],
        "lat_deg": arr[:, 1],
        "radius_m": arr[:, 2],
        "br": arr[:, 3],
        "btheta": arr[:, 4],
        "bphi": arr[:, 5],
        "btot": arr[:, 6],
    }


def load_observations_xyz_csv(path: Path) -> Dict[str, np.ndarray]:
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return {
        "lon_deg": arr[:, 0],
        "lat_deg": arr[:, 1],
        "radius_m": arr[:, 2],
        "bx": arr[:, 3],
        "by": arr[:, 4],
        "bz": arr[:, 5],
        "btot": arr[:, 6],
    }


def write_observations_csv(path: Path, data: Dict[str, np.ndarray], radius_m: float) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["longitude_deg", "latitude_deg", "radius_m", "Br", "Btheta", "Bphi", "Btot"])
        for i in range(data["lon_deg"].size):
            w.writerow(
                [
                    f"{data['lon_deg'][i]:.6f}",
                    f"{data['lat_deg'][i]:.6f}",
                    f"{radius_m:.3f}",
                    f"{data['br'][i]:.9e}",
                    f"{data['btheta'][i]:.9e}",
                    f"{data['bphi'][i]:.9e}",
                    f"{data['btot'][i]:.9e}",
                ]
            )


def write_observations_xyz_csv(path: Path, data: Dict[str, np.ndarray], radius_m: float) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["longitude_deg", "latitude_deg", "radius_m", "Bx", "By", "Bz", "Btot"])
        for i in range(data["lon_deg"].size):
            w.writerow(
                [
                    f"{data['lon_deg'][i]:.6f}",
                    f"{data['lat_deg'][i]:.6f}",
                    f"{radius_m:.3f}",
                    f"{data['bx'][i]:.9e}",
                    f"{data['by'][i]:.9e}",
                    f"{data['bz'][i]:.9e}",
                    f"{data['btot'][i]:.9e}",
                ]
            )


def write_single_component_csv(
    path: Path,
    lon_deg: np.ndarray,
    lat_deg: np.ndarray,
    value: np.ndarray,
    radius_m: float,
    value_name: str,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["longitude_deg", "latitude_deg", "radius_m", value_name])
        for lo, la, v in zip(lon_deg, lat_deg, value):
            w.writerow([f"{lo:.6f}", f"{la:.6f}", f"{radius_m:.3f}", f"{v:.9e}"])


def discover_examples(examples_dir: Path, selected: Sequence[str]) -> List[Path]:
    all_examples = sorted(examples_dir.glob("*.in"))
    if not selected:
        return all_examples
    lookup = {p.name: p for p in all_examples}
    out: List[Path] = []
    for s in selected:
        p = lookup.get(s)
        if p is None:
            raise FileNotFoundError(f"Example not found in {examples_dir}: {s}")
        out.append(p)
    return out


def format_examples(
    *,
    examples_dir: Path,
    out_root: Path,
    rsphere_km: float,
    refine_factor: int,
    selected: Sequence[str],
) -> None:
    out_root.mkdir(parents=True, exist_ok=True)
    run_cmd(["./build_gravmag_tools.sh"], cwd=example_root)
    examples = discover_examples(examples_dir, selected)

    for example in examples:
        stem = example.stem
        case_dir = out_root / stem
        case_dir.mkdir(parents=True, exist_ok=True)
        xyz_out = case_dir / f"{stem}_direct_xyz.txt"
        brtp_out = case_dir / f"{stem}_direct_brtp.txt"

        env = os.environ.copy()
        env["SKIP_BUILD"] = "1"
        run_cmd(
            [
                "./run_input_to_xyz.sh",
                "point_source",
                f"{rsphere_km}",
                str(example),
                str(xyz_out),
                "--refine-factor",
                f"{refine_factor}",
            ],
            cwd=example_root,
            env=env,
        )
        run_cmd(["./run_xyz_to_brtp.sh", str(xyz_out), str(brtp_out)], cwd=example_root, env=env)

        meta = read_example_metadata(example, rsphere_km)
        agg = load_brtp_aggregate(brtp_out)
        agg_xyz = load_xyz_aggregate(xyz_out)
        radius_m = float(meta["radius_m"])

        shutil.copy2(example, case_dir / example.name)
        write_observations_csv(case_dir / "observations_spherical.csv", agg, radius_m)
        write_observations_xyz_csv(case_dir / "observations_xyz.csv", agg_xyz, radius_m)
        write_observations_xyz_csv(case_dir / "harmonica_xyz_input.csv", agg_xyz, radius_m)
        write_single_component_csv(
            case_dir / "shtools_lsq_input.csv",
            agg["lon_deg"],
            agg["lat_deg"],
            agg["br"],
            radius_m,
            "Br",
        )
        write_single_component_csv(
            case_dir / "scipy_sh_lsq_input.csv",
            agg["lon_deg"],
            agg["lat_deg"],
            agg["br"],
            radius_m,
            "Br",
        )
        meta_out = dict(meta)
        meta_out.update(
            {
                "example_path": str(example),
                "baseline_xyz": str(xyz_out),
                "baseline_brtp": str(brtp_out),
                "refine_factor": refine_factor,
                "rsphere_km": rsphere_km,
            }
        )
        with (case_dir / "metadata.json").open("w", encoding="utf-8") as f:
            json.dump(meta_out, f, indent=2)


def align_component_maps_generic(
    ref: Dict[str, np.ndarray],
    pred: Dict[str, np.ndarray],
    components: Sequence[str],
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    ref_map: Dict[Tuple[float, float], Dict[str, float]] = {}
    pred_map: Dict[Tuple[float, float], Dict[str, float]] = {}

    for i in range(ref["lon_deg"].size):
        key = (float(ref["lon_deg"][i]), float(ref["lat_deg"][i]))
        ref_map[key] = {c: float(ref[c][i]) for c in components}
    for i in range(pred["lon_deg"].size):
        key = (float(pred["lon_deg"][i]), float(pred["lat_deg"][i]))
        pred_map[key] = {c: float(pred[c][i]) for c in components}

    keys = sorted(set(ref_map.keys()) & set(pred_map.keys()))
    if not keys:
        raise RuntimeError("No overlapping lon/lat coordinates for comparison.")

    lon = np.array([k[0] for k in keys], dtype=float)
    lat = np.array([k[1] for k in keys], dtype=float)
    ref_aligned = {c: np.array([ref_map[k][c] for k in keys], dtype=float) for c in components}
    pred_aligned = {c: np.array([pred_map[k][c] for k in keys], dtype=float) for c in components}
    return lon, lat, ref_aligned, pred_aligned


def align_component_maps(
    ref: Dict[str, np.ndarray],
    pred: Dict[str, np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    return align_component_maps_generic(ref, pred, COMPONENTS)


def align_component_maps_xyz(
    ref: Dict[str, np.ndarray],
    pred: Dict[str, np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    return align_component_maps_generic(ref, pred, XYZ_COMPONENTS)


def metric_values(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    residual = y_pred - y_true
    rmse = float(np.sqrt(np.mean(residual * residual)))
    mae = float(np.mean(np.abs(residual)))
    max_abs = float(np.max(np.abs(residual)))
    ss_res = float(np.sum(residual * residual))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    r2 = float("nan") if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
    return {"rmse": rmse, "mae": mae, "max_abs": max_abs, "r2": r2}


def all_component_metrics_generic(
    ref: Dict[str, np.ndarray],
    pred: Dict[str, np.ndarray],
    components: Sequence[str],
) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for c in components:
        m = metric_values(ref[c], pred[c])
        for stat in METRIC_STATS:
            out[f"{stat}_{c}"] = m[stat]
    return out


def all_component_metrics(ref: Dict[str, np.ndarray], pred: Dict[str, np.ndarray]) -> Dict[str, float]:
    return all_component_metrics_generic(ref, pred, COMPONENTS)


def all_component_metrics_xyz(ref: Dict[str, np.ndarray], pred: Dict[str, np.ndarray]) -> Dict[str, float]:
    return all_component_metrics_generic(ref, pred, XYZ_COMPONENTS)


def design_matrix_real_sh(lat_deg: np.ndarray, lon_deg: np.ndarray, lmax: int) -> np.ndarray:
    lat_rad = np.deg2rad(lat_deg)
    lon_rad = np.deg2rad(lon_deg)
    x = np.sin(lat_rad)  # cos(colatitude)
    cols: List[np.ndarray] = []
    for l in range(0, lmax + 1):
        for m in range(0, l + 1):
            plm = lpmv(m, l, x)
            cols.append(plm * np.cos(m * lon_rad))
            if m > 0:
                cols.append(plm * np.sin(m * lon_rad))
    return np.column_stack(cols)


def solve_scipy_sh_lsq_component(
    lat_deg: np.ndarray,
    lon_deg: np.ndarray,
    data: np.ndarray,
    lmax: int,
    ridge_lambda: float,
) -> np.ndarray:
    a = design_matrix_real_sh(lat_deg, lon_deg, lmax)
    if ridge_lambda > 0.0:
        ata = a.T @ a + ridge_lambda * np.eye(a.shape[1], dtype=float)
        atd = a.T @ data
        coef = np.linalg.solve(ata, atd)
    else:
        coef, *_ = np.linalg.lstsq(a, data, rcond=None)
    return a @ coef


def solve_shtools_lsq_component(
    lat_deg: np.ndarray,
    lon_deg: np.ndarray,
    data: np.ndarray,
    lmax: int,
) -> np.ndarray:
    import pyshtools as pysh  # type: ignore

    cilm, _chi2 = pysh.expand.SHExpandLSQ(data, lat_deg, lon_deg, lmax)
    g = pysh.expand.LSQ_G(lat_deg, lon_deg, lmax)
    coeffs: List[float] = []
    for l in range(0, lmax + 1):
        for m in range(0, l + 1):
            coeffs.append(float(cilm[0, l, m]))
        for m in range(1, l + 1):
            coeffs.append(float(cilm[1, l, m]))
    coeff_vec = np.array(coeffs, dtype=float)
    return g @ coeff_vec


def predict_scipy_components(
    obs: Dict[str, np.ndarray],
    lmax: int,
    ridge_lambda: float,
) -> Dict[str, np.ndarray]:
    out: Dict[str, np.ndarray] = {
        "lon_deg": obs["lon_deg"].copy(),
        "lat_deg": obs["lat_deg"].copy(),
    }
    out["br"] = solve_scipy_sh_lsq_component(obs["lat_deg"], obs["lon_deg"], obs["br"], lmax, ridge_lambda)
    out["btheta"] = solve_scipy_sh_lsq_component(
        obs["lat_deg"], obs["lon_deg"], obs["btheta"], lmax, ridge_lambda
    )
    out["bphi"] = solve_scipy_sh_lsq_component(obs["lat_deg"], obs["lon_deg"], obs["bphi"], lmax, ridge_lambda)
    out["btot"] = np.sqrt(out["br"] * out["br"] + out["btheta"] * out["btheta"] + out["bphi"] * out["bphi"])
    return out


def predict_shtools_components(obs: Dict[str, np.ndarray], lmax: int) -> Dict[str, np.ndarray]:
    out: Dict[str, np.ndarray] = {
        "lon_deg": obs["lon_deg"].copy(),
        "lat_deg": obs["lat_deg"].copy(),
    }
    out["br"] = solve_shtools_lsq_component(obs["lat_deg"], obs["lon_deg"], obs["br"], lmax)
    out["btheta"] = solve_shtools_lsq_component(obs["lat_deg"], obs["lon_deg"], obs["btheta"], lmax)
    out["bphi"] = solve_shtools_lsq_component(obs["lat_deg"], obs["lon_deg"], obs["bphi"], lmax)
    out["btot"] = np.sqrt(out["br"] * out["br"] + out["btheta"] * out["btheta"] + out["bphi"] * out["bphi"])
    return out


def run_fortran_spectral_xyz_output(
    *,
    example_path: Path,
    case_dir: Path,
    rsphere_km: float,
    refine_factor: int,
    ntheta_fit: int,
    nphi_fit: int,
    auto_mode: bool,
    lmax_seed: int,
    reg_lambda_seed: float,
    reg_power_seed: float,
) -> Tuple[Path, Dict[str, float]]:
    xyz_out = case_dir / f"{example_path.stem}_spectral_xyz.txt"
    env = os.environ.copy()
    env["SKIP_BUILD"] = "1"
    cmd = [
        "./run_input_to_xyz.sh",
        "spectral",
        f"{rsphere_km}",
        str(example_path),
        str(xyz_out),
    ]
    if auto_mode:
        # Keep spectral tuning in auto mode; only override sampling if requested.
        if refine_factor != 2:
            cmd.extend(["--refine-factor", f"{refine_factor}"])
        if ntheta_fit != 72:
            cmd.extend(["--ntheta-fit", f"{ntheta_fit}"])
        if nphi_fit != 144:
            cmd.extend(["--nphi-fit", f"{nphi_fit}"])
    else:
        # use the manual seed settings
        cmd.extend(
            [
                f"{refine_factor}",
                f"{lmax_seed}",
                f"{ntheta_fit}",
                f"{nphi_fit}",
                f"{reg_lambda_seed}",
                f"{reg_power_seed}",
                "0",
                "0",
                "0",
                "0",
            ]
        )

    run_cmd(cmd, cwd=example_root, env=env)
    hdr = parse_gauss_header_params(xyz_out)
    if "lmax_used" not in hdr:
        hdr["lmax_used"] = float(lmax_seed)
    if "reg_lambda_used" not in hdr:
        hdr["reg_lambda_used"] = float(reg_lambda_seed)
    if "reg_power_used" not in hdr:
        hdr["reg_power_used"] = float(reg_power_seed)
    return xyz_out, hdr


def run_fortran_spectral_prediction(
    *,
    example_path: Path,
    case_dir: Path,
    rsphere_km: float,
    refine_factor: int,
    ntheta_fit: int,
    nphi_fit: int,
    auto_mode: bool,
    lmax_seed: int,
    reg_lambda_seed: float,
    reg_power_seed: float,
) -> Tuple[Dict[str, np.ndarray], Dict[str, float]]:
    xyz_out, hdr = run_fortran_spectral_xyz_output(
        example_path=example_path,
        case_dir=case_dir,
        rsphere_km=rsphere_km,
        refine_factor=refine_factor,
        ntheta_fit=ntheta_fit,
        nphi_fit=nphi_fit,
        auto_mode=auto_mode,
        lmax_seed=lmax_seed,
        reg_lambda_seed=reg_lambda_seed,
        reg_power_seed=reg_power_seed,
    )
    brtp_out = case_dir / f"{example_path.stem}_spectral_brtp.txt"
    env = os.environ.copy()
    env["SKIP_BUILD"] = "1"
    run_cmd(["./run_xyz_to_brtp.sh", str(xyz_out), str(brtp_out)], cwd=example_root, env=env)
    return load_brtp_aggregate(brtp_out), hdr


def run_fortran_spectral_prediction_xyz(
    *,
    example_path: Path,
    case_dir: Path,
    rsphere_km: float,
    refine_factor: int,
    ntheta_fit: int,
    nphi_fit: int,
    auto_mode: bool,
    lmax_seed: int,
    reg_lambda_seed: float,
    reg_power_seed: float,
) -> Tuple[Dict[str, np.ndarray], Dict[str, float]]:
    xyz_out, hdr = run_fortran_spectral_xyz_output(
        example_path=example_path,
        case_dir=case_dir,
        rsphere_km=rsphere_km,
        refine_factor=refine_factor,
        ntheta_fit=ntheta_fit,
        nphi_fit=nphi_fit,
        auto_mode=auto_mode,
        lmax_seed=lmax_seed,
        reg_lambda_seed=reg_lambda_seed,
        reg_power_seed=reg_power_seed,
    )
    return load_xyz_aggregate(xyz_out), hdr


def local_xy_from_lonlat(
    lon_deg: np.ndarray,
    lat_deg: np.ndarray,
    lon0_deg: float,
    lat0_deg: float,
    radius_m: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    lon_delta_rad = np.deg2rad(lon_deg - lon0_deg)
    lat_delta_rad = np.deg2rad(lat_deg - lat0_deg)
    lat0_rad = np.deg2rad(lat0_deg)
    x = radius_m * np.cos(lat0_rad) * lon_delta_rad
    y = radius_m * lat_delta_rad
    z = np.zeros_like(x)
    return x, y, z


def estimate_mean_grid_spacing_m(
    lon_deg: np.ndarray,
    lat_deg: np.ndarray,
    lat0_deg: float,
    radius_m: float,
) -> float:
    ulon = np.unique(np.round(lon_deg, 6))
    ulat = np.unique(np.round(lat_deg, 6))
    dlon = float(np.median(np.abs(np.diff(ulon)))) if ulon.size > 1 else 0.1
    dlat = float(np.median(np.abs(np.diff(ulat)))) if ulat.size > 1 else 0.1
    lat0_rad = np.deg2rad(lat0_deg)
    dx = abs(radius_m * np.cos(lat0_rad) * np.deg2rad(dlon))
    dy = abs(radius_m * np.deg2rad(dlat))
    spacing = max(1.0, float(np.mean([dx, dy])))
    return spacing


def predict_harmonica_xyz_components(
    obs_xyz: Dict[str, np.ndarray],
    *,
    lon0_deg: float,
    lat0_deg: float,
    radius_m: float,
    damping: float,
    depth_m: float,
    block_size_m: float,
) -> Dict[str, np.ndarray]:
    import harmonica as hm  # type: ignore

    easting, northing, upward = local_xy_from_lonlat(
        obs_xyz["lon_deg"], obs_xyz["lat_deg"], lon0_deg, lat0_deg, radius_m
    )
    coords = (easting, northing, upward)
    out: Dict[str, np.ndarray] = {
        "lon_deg": obs_xyz["lon_deg"].copy(),
        "lat_deg": obs_xyz["lat_deg"].copy(),
    }
    damping_value: float | None = None if damping <= 0.0 else float(damping)
    for comp in ("bx", "by", "bz"):
        eqs = hm.EquivalentSources(depth=float(depth_m), damping=damping_value, block_size=float(block_size_m))
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=FutureWarning,
                message=r"The provided callable <function median.*",
            )
            eqs.fit(coords, obs_xyz[comp])
        out[comp] = np.asarray(eqs.predict(coords), dtype=float)
    out["btot"] = np.sqrt(out["bx"] * out["bx"] + out["by"] * out["by"] + out["bz"] * out["bz"])
    return out


def format_harmonica_params_label(damping_used: float, depth_m_used: float) -> str:
    return f"damping={damping_used:g}, depth_m={depth_m_used:g}"


def choose_best_harmonica_prediction(
    obs_xyz: Dict[str, np.ndarray],
    *,
    lon0_deg: float,
    lat0_deg: float,
    radius_m: float,
    damping_grid: Sequence[float],
    depth_factor_grid: Sequence[float],
    block_factor: float,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, float], Dict[str, float]]:
    spacing_m = estimate_mean_grid_spacing_m(obs_xyz["lon_deg"], obs_xyz["lat_deg"], lat0_deg, radius_m)
    block_size_m = max(1.0, float(block_factor) * spacing_m)
    best_err = float("inf")
    best_tuple: (
        Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, float], Dict[str, float]]
        | None
    ) = None

    for depth_factor in depth_factor_grid:
        depth_m = max(1.0, float(depth_factor) * spacing_m)
        for damping in damping_grid:
            pred = predict_harmonica_xyz_components(
                obs_xyz,
                lon0_deg=lon0_deg,
                lat0_deg=lat0_deg,
                radius_m=radius_m,
                damping=float(damping),
                depth_m=depth_m,
                block_size_m=block_size_m,
            )
            lon, lat, ref_aligned, pred_aligned = align_component_maps_xyz(obs_xyz, pred)
            m = all_component_metrics_xyz(ref_aligned, pred_aligned)
            err = float(m["rmse_btot"])
            if err < best_err:
                best_err = err
                params = {
                    "harmonica_damping_used": float(damping),
                    "harmonica_depth_factor_used": float(depth_factor),
                    "harmonica_depth_m_used": float(depth_m),
                    "harmonica_block_size_m_used": float(block_size_m),
                }
                best_tuple = (lon, lat, ref_aligned, pred_aligned, m, params)
    if best_tuple is None:
        raise RuntimeError("No successful Harmonica equivalent-source predictions in optimization grid.")
    return best_tuple


def to_grid(lon: np.ndarray, lat: np.ndarray, values: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    lon_r = np.round(lon, 6)
    lat_r = np.round(lat, 6)
    ulon = np.unique(lon_r)
    ulat = np.unique(lat_r)
    ulon.sort()
    ulat.sort()
    i_lon = {x: i for i, x in enumerate(ulon)}
    i_lat = {y: i for i, y in enumerate(ulat)}
    z = np.full((ulat.size, ulon.size), np.nan, dtype=float)
    for x, y, v in zip(lon_r, lat_r, values):
        z[i_lat[y], i_lon[x]] = v
    if np.isnan(z).any():
        raise RuntimeError("Residual grid has missing cells.")
    return ulon, ulat, z


def plot_residuals(
    *,
    lon: np.ndarray,
    lat: np.ndarray,
    ref: Dict[str, np.ndarray],
    pred: Dict[str, np.ndarray],
    solver: str,
    case: str,
    params_label: str,
    out_path: Path,
) -> None:
    if "MPLCONFIGDIR" not in os.environ:
        mpl_dir = repo_root / "diagnostics" / ".mplconfig"
        mpl_dir.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_dir)

    import matplotlib

    if os.environ.get("DISPLAY", "") == "" and os.environ.get("MPLBACKEND", "") == "":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    comp_labels = [("br", "Br"), ("btheta", "Btheta"), ("bphi", "Bphi"), ("btot", "Btot")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    axes = axes.ravel()

    for ax, (key, label) in zip(axes, comp_labels):
        residual = pred[key] - ref[key]
        ulon, ulat, z = to_grid(lon, lat, residual)
        vmax = float(np.nanmax(np.abs(z)))
        if vmax <= 0.0:
            vmax = 1.0e-12
        im = ax.pcolormesh(ulon, ulat, z, shading="nearest", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(f"{label} residual")
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        ax.set_aspect("equal", adjustable="box")
        cbar = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.12, shrink=0.85)
        cbar.set_label("Pred - baseline")

    fig.suptitle(f"{case}: {solver} residuals (Br/Btheta/Bphi/Btot) [{params_label}]")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_side_by_side_solutions(
    *,
    lon: np.ndarray,
    lat: np.ndarray,
    gravmag: Dict[str, np.ndarray],
    shtools: Dict[str, np.ndarray],
    case: str,
    gravmag_params_label: str,
    shtools_params_label: str,
    out_path: Path,
) -> None:
    if "MPLCONFIGDIR" not in os.environ:
        mpl_dir = repo_root / "diagnostics" / ".mplconfig"
        mpl_dir.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_dir)

    import matplotlib

    if os.environ.get("DISPLAY", "") == "" and os.environ.get("MPLBACKEND", "") == "":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    comp_labels = [("br", "Br"), ("btheta", "Btheta"), ("bphi", "Bphi"), ("btot", "Btot")]
    fig, axes = plt.subplots(4, 2, figsize=(12, 14), constrained_layout=True)

    for i, (key, label) in enumerate(comp_labels):
        ax_g = axes[i, 0]
        ax_s = axes[i, 1]
        ulon_g, ulat_g, z_g = to_grid(lon, lat, gravmag[key])
        ulon_s, ulat_s, z_s = to_grid(lon, lat, shtools[key])

        vmin = float(min(np.nanmin(z_g), np.nanmin(z_s)))
        vmax = float(max(np.nanmax(z_g), np.nanmax(z_s)))
        if vmax <= vmin:
            scale = max(1.0e-12, abs(vmax))
            vmin = -scale
            vmax = scale

        im_g = ax_g.pcolormesh(ulon_g, ulat_g, z_g, shading="nearest", cmap="RdBu_r", vmin=vmin, vmax=vmax)
        ax_s.pcolormesh(ulon_s, ulat_s, z_s, shading="nearest", cmap="RdBu_r", vmin=vmin, vmax=vmax)

        ax_g.set_ylabel("Latitude [deg]")
        ax_g.set_xlabel("Longitude [deg]")
        ax_s.set_xlabel("Longitude [deg]")
        ax_g.set_aspect("equal", adjustable="box")
        ax_s.set_aspect("equal", adjustable="box")
        ax_g.set_title(f"{label} GravMagSphere")
        ax_s.set_title(f"{label} SHTOOLS")

        cbar = fig.colorbar(im_g, ax=[ax_g, ax_s], orientation="horizontal", pad=0.1, shrink=0.95)
        cbar.set_label(label)

    fig.suptitle(
        f"{case}: GravMagSphere vs SHTOOLS fields "
        f"[gravmag: {gravmag_params_label} | shtools: {shtools_params_label}]"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_residuals_xyz(
    *,
    lon: np.ndarray,
    lat: np.ndarray,
    ref: Dict[str, np.ndarray],
    pred: Dict[str, np.ndarray],
    solver: str,
    case: str,
    params_label: str,
    out_path: Path,
    residual_label: str,
) -> None:
    if "MPLCONFIGDIR" not in os.environ:
        mpl_dir = repo_root / "diagnostics" / ".mplconfig"
        mpl_dir.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_dir)

    import matplotlib

    if os.environ.get("DISPLAY", "") == "" and os.environ.get("MPLBACKEND", "") == "":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    comp_labels = [("bx", "Bx"), ("by", "By"), ("bz", "Bz"), ("btot", "Btot")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    axes = axes.ravel()

    for ax, (key, label) in zip(axes, comp_labels):
        residual = pred[key] - ref[key]
        ulon, ulat, z = to_grid(lon, lat, residual)
        vmax = float(np.nanmax(np.abs(z)))
        if vmax <= 0.0:
            vmax = 1.0e-12
        im = ax.pcolormesh(ulon, ulat, z, shading="nearest", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(f"{label} residual")
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        ax.set_aspect("equal", adjustable="box")
        cbar = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.12, shrink=0.85)
        cbar.set_label(residual_label)

    fig.suptitle(f"{case}: {solver} residuals (Bx/By/Bz/Btot) [{params_label}]")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_side_by_side_solutions_xyz(
    *,
    lon: np.ndarray,
    lat: np.ndarray,
    left: Dict[str, np.ndarray],
    right: Dict[str, np.ndarray],
    case: str,
    left_name: str,
    right_name: str,
    left_params_label: str,
    right_params_label: str,
    out_path: Path,
) -> None:
    if "MPLCONFIGDIR" not in os.environ:
        mpl_dir = repo_root / "diagnostics" / ".mplconfig"
        mpl_dir.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_dir)

    import matplotlib

    if os.environ.get("DISPLAY", "") == "" and os.environ.get("MPLBACKEND", "") == "":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    comp_labels = [("bx", "Bx"), ("by", "By"), ("bz", "Bz"), ("btot", "Btot")]
    fig, axes = plt.subplots(4, 2, figsize=(12, 14), constrained_layout=True)

    for i, (key, label) in enumerate(comp_labels):
        ax_l = axes[i, 0]
        ax_r = axes[i, 1]
        ulon_l, ulat_l, z_l = to_grid(lon, lat, left[key])
        ulon_r, ulat_r, z_r = to_grid(lon, lat, right[key])

        vmin = float(min(np.nanmin(z_l), np.nanmin(z_r)))
        vmax = float(max(np.nanmax(z_l), np.nanmax(z_r)))
        if vmax <= vmin:
            scale = max(1.0e-12, abs(vmax))
            vmin = -scale
            vmax = scale

        im_l = ax_l.pcolormesh(ulon_l, ulat_l, z_l, shading="nearest", cmap="RdBu_r", vmin=vmin, vmax=vmax)
        ax_r.pcolormesh(ulon_r, ulat_r, z_r, shading="nearest", cmap="RdBu_r", vmin=vmin, vmax=vmax)

        ax_l.set_ylabel("Latitude [deg]")
        ax_l.set_xlabel("Longitude [deg]")
        ax_r.set_xlabel("Longitude [deg]")
        ax_l.set_aspect("equal", adjustable="box")
        ax_r.set_aspect("equal", adjustable="box")
        ax_l.set_title(f"{label} {left_name}")
        ax_r.set_title(f"{label} {right_name}")

        cbar = fig.colorbar(im_l, ax=[ax_l, ax_r], orientation="horizontal", pad=0.1, shrink=0.95)
        cbar.set_label(label)

    fig.suptitle(
        f"{case}: {left_name} vs {right_name} fields "
        f"[{left_name}: {left_params_label} | {right_name}: {right_params_label}]"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def format_params_label(lmax_used: float, reg_lambda_used: float, reg_power_used: float) -> str:
    return f"lmax={int(round(lmax_used))}, reg_lambda={reg_lambda_used:g}, reg_power={reg_power_used:g}"


def choose_best_scipy_prediction(
    obs: Dict[str, np.ndarray],
    *,
    lmax_grid: Sequence[int],
    ridge_grid: Sequence[float],
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, float], Dict[str, float]]:
    best_err = float("inf")
    best_tuple: Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, float], Dict[str, float]] | None = None
    for lm in lmax_grid:
        for ridge in ridge_grid:
            pred = predict_scipy_components(obs, lmax=int(lm), ridge_lambda=float(ridge))
            lon, lat, ref_aligned, pred_aligned = align_component_maps(obs, pred)
            m = all_component_metrics(ref_aligned, pred_aligned)
            err = float(m["rmse_btot"])
            if err < best_err:
                best_err = err
                params = {
                    "lmax_used": float(lm),
                    "reg_lambda_used": float(ridge),
                    "reg_power_used": 0.0,
                }
                best_tuple = (lon, lat, ref_aligned, pred_aligned, m, params)
    if best_tuple is None:
        raise RuntimeError("No successful SciPy SH LSQ predictions in optimization grid.")
    return best_tuple


def choose_best_shtools_prediction(
    obs: Dict[str, np.ndarray],
    *,
    lmax_grid: Sequence[int],
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, float], Dict[str, float]]:
    best_err = float("inf")
    best_tuple: Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, float], Dict[str, float]] | None = None
    for lm in lmax_grid:
        pred = predict_shtools_components(obs, lmax=int(lm))
        lon, lat, ref_aligned, pred_aligned = align_component_maps(obs, pred)
        m = all_component_metrics(ref_aligned, pred_aligned)
        err = float(m["rmse_btot"])
        if err < best_err:
            best_err = err
            params = {
                "lmax_used": float(lm),
                "reg_lambda_used": 0.0,
                "reg_power_used": 0.0,
            }
            best_tuple = (lon, lat, ref_aligned, pred_aligned, m, params)
    if best_tuple is None:
        raise RuntimeError("No successful SHTOOLS LSQ predictions in optimization grid.")
    return best_tuple


def compare_cases(
    *,
    formatted_root: Path,
    report_md: Path,
    report_csv: Path,
    residual_dir: Path,
    side_by_side_dir: Path,
    lmax: int,
    ridge_lambda: float,
    scipy_lmax_grid: Sequence[int],
    scipy_ridge_grid: Sequence[float],
    shtools_lmax_grid: Sequence[int],
    fortran_auto_mode: bool,
    rsphere_km: float,
    refine_factor: int,
    ntheta_fit: int,
    nphi_fit: int,
    reg_lambda: float,
    reg_power: float,
    selected: Sequence[str],
) -> None:
    case_dirs = sorted([p for p in formatted_root.iterdir() if p.is_dir()])
    if selected:
        selected_set = set(selected)
        case_dirs = [p for p in case_dirs if p.name in selected_set]

    rows: List[Dict[str, object]] = []

    for case_dir in case_dirs:
        fortran_pred_raw: Dict[str, np.ndarray] | None = None
        fortran_used: Dict[str, float] | None = None
        obs_csv = case_dir / "observations_spherical.csv"
        meta_json = case_dir / "metadata.json"
        if not obs_csv.exists() or not meta_json.exists():
            continue

        with meta_json.open("r", encoding="utf-8") as f:
            meta = json.load(f)
        example_path = Path(meta["example_path"])
        obs = load_observations_csv(obs_csv)

        # 1) Fortran spectral backend
        try:
            pred, used = run_fortran_spectral_prediction(
                example_path=example_path,
                case_dir=case_dir,
                rsphere_km=rsphere_km,
                refine_factor=refine_factor,
                ntheta_fit=ntheta_fit,
                nphi_fit=nphi_fit,
                auto_mode=fortran_auto_mode,
                lmax_seed=lmax,
                reg_lambda_seed=reg_lambda,
                reg_power_seed=reg_power,
            )
            lon, lat, ref_aligned, pred_aligned = align_component_maps(obs, pred)
            m = all_component_metrics(ref_aligned, pred_aligned)
            plot_path = residual_dir / f"{case_dir.name}_fortran_spectral_residuals.png"
            plot_residuals(
                lon=lon,
                lat=lat,
                ref=ref_aligned,
                pred=pred_aligned,
                solver="fortran_spectral",
                case=case_dir.name,
                params_label=format_params_label(
                    used["lmax_used"], used["reg_lambda_used"], used["reg_power_used"]
                ),
                out_path=plot_path,
            )
            rows.append(
                {
                    "case": case_dir.name,
                    "solver": "fortran_spectral",
                    "status": "ok",
                    "plot": str(plot_path),
                    "side_by_side_plot": "",
                    **used,
                    **m,
                }
            )
            fortran_pred_raw = pred
            fortran_used = used
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "case": case_dir.name,
                    "solver": "fortran_spectral",
                    "status": f"error: {exc}",
                    "side_by_side_plot": "",
                }
            )

        # 2) SciPy SH LSQ backend (optimized over lmax/ridge grid)
        try:
            lon, lat, ref_aligned, pred_aligned, m, used = choose_best_scipy_prediction(
                obs, lmax_grid=scipy_lmax_grid, ridge_grid=scipy_ridge_grid
            )
            plot_path = residual_dir / f"{case_dir.name}_scipy_sh_lsq_residuals.png"
            plot_residuals(
                lon=lon,
                lat=lat,
                ref=ref_aligned,
                pred=pred_aligned,
                solver="scipy_sh_lsq",
                case=case_dir.name,
                params_label=format_params_label(
                    used["lmax_used"], used["reg_lambda_used"], used["reg_power_used"]
                ),
                out_path=plot_path,
            )
            rows.append(
                {
                    "case": case_dir.name,
                    "solver": "scipy_sh_lsq",
                    "status": "ok",
                    "plot": str(plot_path),
                    "side_by_side_plot": "",
                    **used,
                    **m,
                }
            )
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "case": case_dir.name,
                    "solver": "scipy_sh_lsq",
                    "status": f"error: {exc}",
                    "side_by_side_plot": "",
                }
            )

        # 3) SHTOOLS LSQ backend (optimized over lmax grid)
        try:
            lon, lat, ref_aligned, pred_aligned, m, used = choose_best_shtools_prediction(
                obs, lmax_grid=shtools_lmax_grid
            )
            plot_path = residual_dir / f"{case_dir.name}_shtools_lsq_residuals.png"
            plot_residuals(
                lon=lon,
                lat=lat,
                ref=ref_aligned,
                pred=pred_aligned,
                solver="shtools_lsq",
                case=case_dir.name,
                params_label=format_params_label(
                    used["lmax_used"], used["reg_lambda_used"], used["reg_power_used"]
                ),
                out_path=plot_path,
            )

            side_plot = ""
            if fortran_pred_raw is not None and fortran_used is not None:
                shtools_raw = predict_shtools_components(obs, lmax=int(round(used["lmax_used"])))
                lon_ss, lat_ss, gravmag_aligned, shtools_aligned = align_component_maps(fortran_pred_raw, shtools_raw)
                side_plot_path = side_by_side_dir / f"{case_dir.name}_gravmagsphere_vs_shtools.png"
                plot_side_by_side_solutions(
                    lon=lon_ss,
                    lat=lat_ss,
                    gravmag=gravmag_aligned,
                    shtools=shtools_aligned,
                    case=case_dir.name,
                    gravmag_params_label=format_params_label(
                        fortran_used["lmax_used"],
                        fortran_used["reg_lambda_used"],
                        fortran_used["reg_power_used"],
                    ),
                    shtools_params_label=format_params_label(
                        used["lmax_used"], used["reg_lambda_used"], used["reg_power_used"]
                    ),
                    out_path=side_plot_path,
                )
                side_plot = str(side_plot_path)

            rows.append(
                {
                    "case": case_dir.name,
                    "solver": "shtools_lsq",
                    "status": "ok",
                    "plot": str(plot_path),
                    "side_by_side_plot": side_plot,
                    **used,
                    **m,
                }
            )
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "case": case_dir.name,
                    "solver": "shtools_lsq",
                    "status": f"error: {exc}",
                    "side_by_side_plot": "",
                }
            )

    report_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["case", "solver", "status", "plot", "side_by_side_plot", "lmax_used", "reg_lambda_used", "reg_power_used"] + [
        f"{s}_{c}" for c in COMPONENTS for s in METRIC_STATS
    ]
    with report_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    report_md.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# External Solver Comparison",
        "",
        f"- fortran mode: `{'auto' if fortran_auto_mode else 'manual-seeded'}`",
        f"- fortran seed (manual mode): lmax `{lmax}`, reg_lambda `{reg_lambda}`, reg_power `{reg_power}`",
        f"- scipy optimization grid: lmax `{','.join(str(v) for v in scipy_lmax_grid)}`, ridge `{','.join(f'{v:g}' for v in scipy_ridge_grid)}`",
        f"- shtools optimization grid: lmax `{','.join(str(v) for v in shtools_lmax_grid)}`",
        f"- residual plots: `{residual_dir}`",
        f"- GravMagSphere vs SHTOOLS side-by-side plots: `{side_by_side_dir}`",
        "",
        "| Case | Solver | Status | lmax_used | reg_lambda_used | reg_power_used | RMSE(Br) | RMSE(Btheta) | RMSE(Bphi) | RMSE(Btot) | R2(Btot) | Residual plot | Side-by-side plot |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        if str(row["status"]) != "ok":
            lines.append(f"| {row['case']} | {row['solver']} | {row['status']} | - | - | - | - | - | - | - | - | - | - |")
        else:
            lines.append(
                "| {case} | {solver} | {status} | {lmax_used:.0f} | {reg_lambda_used:.6g} | {reg_power_used:.6g} | "
                "{rmse_br:.6e} | {rmse_btheta:.6e} | {rmse_bphi:.6e} | {rmse_btot:.6e} | {r2_btot:.6f} | {plot} | {side_by_side_plot} |".format(
                    **row
                )
            )
    lines.append("")
    lines.append("Residuals are `prediction - point_source-baseline` on the shared lon/lat grid.")
    lines.append("SciPy and SHTOOLS rows use per-case optimized settings from the configured search grids.")
    report_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def ensure_observations_xyz_csv(case_dir: Path, meta: Dict[str, object]) -> Path:
    obs_xyz_csv = case_dir / "observations_xyz.csv"
    if obs_xyz_csv.exists():
        return obs_xyz_csv

    baseline_xyz: Path | None = None
    baseline_raw = str(meta.get("baseline_xyz", "")).strip()
    if baseline_raw:
        candidate = Path(baseline_raw)
        baseline_xyz = candidate if candidate.is_absolute() else repo_root / candidate
    if baseline_xyz is None or not baseline_xyz.exists():
        fallback = case_dir / f"{case_dir.name}_direct_xyz.txt"
        if fallback.exists():
            baseline_xyz = fallback
    if baseline_xyz is None or not baseline_xyz.exists():
        raise RuntimeError(f"Could not locate baseline xyz output for {case_dir}")

    radius_m = float(meta.get("radius_m", 0.0))
    if radius_m <= 0.0:
        raise RuntimeError(f"Missing or invalid radius_m in metadata for {case_dir}")

    agg_xyz = load_xyz_aggregate(baseline_xyz)
    write_observations_xyz_csv(obs_xyz_csv, agg_xyz, radius_m)
    write_observations_xyz_csv(case_dir / "harmonica_xyz_input.csv", agg_xyz, radius_m)
    return obs_xyz_csv


def compare_cases_harmonica_xyz(
    *,
    formatted_root: Path,
    report_md: Path,
    report_csv: Path,
    residual_dir: Path,
    side_by_side_dir: Path,
    fortran_auto_mode: bool,
    rsphere_km: float,
    refine_factor: int,
    ntheta_fit: int,
    nphi_fit: int,
    lmax: int,
    reg_lambda: float,
    reg_power: float,
    harmonica_damping_grid: Sequence[float],
    harmonica_depth_factor_grid: Sequence[float],
    harmonica_block_factor: float,
    selected: Sequence[str],
) -> None:
    case_dirs = sorted([p for p in formatted_root.iterdir() if p.is_dir()])
    if selected:
        selected_set = set(selected)
        case_dirs = [p for p in case_dirs if p.name in selected_set]

    rows: List[Dict[str, object]] = []
    for case_dir in case_dirs:
        meta_json = case_dir / "metadata.json"
        if not meta_json.exists():
            continue
        with meta_json.open("r", encoding="utf-8") as f:
            meta = json.load(f)

        obs_xyz_csv = ensure_observations_xyz_csv(case_dir, meta)
        obs_xyz = load_observations_xyz_csv(obs_xyz_csv)
        example_path = Path(str(meta["example_path"]))
        lon0_deg = float(meta.get("lon0_deg", float(np.mean(obs_xyz["lon_deg"]))))
        lat0_deg = float(meta.get("lat0_deg", float(np.mean(obs_xyz["lat_deg"]))))
        radius_m = float(meta.get("radius_m", (rsphere_km * 1000.0)))

        fortran_pred_raw: Dict[str, np.ndarray] | None = None
        harmonica_pred_raw: Dict[str, np.ndarray] | None = None
        fortran_used: Dict[str, float] | None = None
        harmonica_used: Dict[str, float] | None = None

        # 1) GravMagSphere spectral vs point_source baseline (XYZ)
        try:
            pred_fortran, used_fortran = run_fortran_spectral_prediction_xyz(
                example_path=example_path,
                case_dir=case_dir,
                rsphere_km=rsphere_km,
                refine_factor=refine_factor,
                ntheta_fit=ntheta_fit,
                nphi_fit=nphi_fit,
                auto_mode=fortran_auto_mode,
                lmax_seed=lmax,
                reg_lambda_seed=reg_lambda,
                reg_power_seed=reg_power,
            )
            lon, lat, ref_aligned, pred_aligned = align_component_maps_xyz(obs_xyz, pred_fortran)
            m = all_component_metrics_xyz(ref_aligned, pred_aligned)
            plot_path = residual_dir / f"{case_dir.name}_fortran_spectral_xyz_residuals.png"
            plot_residuals_xyz(
                lon=lon,
                lat=lat,
                ref=ref_aligned,
                pred=pred_aligned,
                solver="fortran_spectral_xyz",
                case=case_dir.name,
                params_label=format_params_label(
                    used_fortran["lmax_used"], used_fortran["reg_lambda_used"], used_fortran["reg_power_used"]
                ),
                out_path=plot_path,
                residual_label="Pred - baseline",
            )
            row = {
                "case": case_dir.name,
                "solver": "fortran_spectral_vs_direct_xyz",
                "status": "ok",
                "plot": str(plot_path),
                "side_by_side_plot": "",
                "lmax_used": used_fortran["lmax_used"],
                "reg_lambda_used": used_fortran["reg_lambda_used"],
                "reg_power_used": used_fortran["reg_power_used"],
                "harmonica_damping_used": float("nan"),
                "harmonica_depth_factor_used": float("nan"),
                "harmonica_depth_m_used": float("nan"),
                "harmonica_block_size_m_used": float("nan"),
            }
            row.update(m)
            rows.append(row)
            fortran_pred_raw = pred_fortran
            fortran_used = used_fortran
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "case": case_dir.name,
                    "solver": "fortran_spectral_vs_direct_xyz",
                    "status": f"error: {exc}",
                    "plot": "",
                    "side_by_side_plot": "",
                    "lmax_used": float("nan"),
                    "reg_lambda_used": float("nan"),
                    "reg_power_used": float("nan"),
                    "harmonica_damping_used": float("nan"),
                    "harmonica_depth_factor_used": float("nan"),
                    "harmonica_depth_m_used": float("nan"),
                    "harmonica_block_size_m_used": float("nan"),
                }
            )

        # 2) Harmonica equivalent sources vs point_source baseline (XYZ)
        try:
            lon, lat, ref_aligned, pred_aligned, m, used_h = choose_best_harmonica_prediction(
                obs_xyz,
                lon0_deg=lon0_deg,
                lat0_deg=lat0_deg,
                radius_m=radius_m,
                damping_grid=harmonica_damping_grid,
                depth_factor_grid=harmonica_depth_factor_grid,
                block_factor=harmonica_block_factor,
            )
            plot_path = residual_dir / f"{case_dir.name}_harmonica_eqs_xyz_residuals.png"
            plot_residuals_xyz(
                lon=lon,
                lat=lat,
                ref=ref_aligned,
                pred=pred_aligned,
                solver="harmonica_eqs_xyz",
                case=case_dir.name,
                params_label=format_harmonica_params_label(
                    used_h["harmonica_damping_used"], used_h["harmonica_depth_m_used"]
                ),
                out_path=plot_path,
                residual_label="Pred - baseline",
            )
            row = {
                "case": case_dir.name,
                "solver": "harmonica_eqs_vs_direct_xyz",
                "status": "ok",
                "plot": str(plot_path),
                "side_by_side_plot": "",
                "lmax_used": float("nan"),
                "reg_lambda_used": float("nan"),
                "reg_power_used": float("nan"),
                "harmonica_damping_used": used_h["harmonica_damping_used"],
                "harmonica_depth_factor_used": used_h["harmonica_depth_factor_used"],
                "harmonica_depth_m_used": used_h["harmonica_depth_m_used"],
                "harmonica_block_size_m_used": used_h["harmonica_block_size_m_used"],
            }
            row.update(m)
            rows.append(row)
            harmonica_pred_raw = {
                "lon_deg": lon.copy(),
                "lat_deg": lat.copy(),
                "bx": pred_aligned["bx"].copy(),
                "by": pred_aligned["by"].copy(),
                "bz": pred_aligned["bz"].copy(),
                "btot": pred_aligned["btot"].copy(),
            }
            harmonica_used = used_h
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "case": case_dir.name,
                    "solver": "harmonica_eqs_vs_direct_xyz",
                    "status": f"error: {exc}",
                    "plot": "",
                    "side_by_side_plot": "",
                    "lmax_used": float("nan"),
                    "reg_lambda_used": float("nan"),
                    "reg_power_used": float("nan"),
                    "harmonica_damping_used": float("nan"),
                    "harmonica_depth_factor_used": float("nan"),
                    "harmonica_depth_m_used": float("nan"),
                    "harmonica_block_size_m_used": float("nan"),
                }
            )

        # 3) GravMagSphere spectral - Harmonica residuals + side-by-side fields
        if fortran_pred_raw is not None and harmonica_pred_raw is not None and fortran_used is not None and harmonica_used is not None:
            lon_gh, lat_gh, fortran_aligned, harmonica_aligned = align_component_maps_xyz(
                fortran_pred_raw, harmonica_pred_raw
            )
            m_gh = all_component_metrics_xyz(harmonica_aligned, fortran_aligned)
            plot_path = residual_dir / f"{case_dir.name}_gravmagsphere_minus_harmonica_xyz_residuals.png"
            plot_residuals_xyz(
                lon=lon_gh,
                lat=lat_gh,
                ref=harmonica_aligned,
                pred=fortran_aligned,
                solver="gravmagsphere_minus_harmonica_xyz",
                case=case_dir.name,
                params_label=(
                    f"gravmag[{format_params_label(fortran_used['lmax_used'], fortran_used['reg_lambda_used'], fortran_used['reg_power_used'])}]"
                    f"; harmonica[{format_harmonica_params_label(harmonica_used['harmonica_damping_used'], harmonica_used['harmonica_depth_m_used'])}]"
                ),
                out_path=plot_path,
                residual_label="GravMag - Harmonica",
            )
            side_plot_path = side_by_side_dir / f"{case_dir.name}_gravmagsphere_vs_harmonica_xyz.png"
            plot_side_by_side_solutions_xyz(
                lon=lon_gh,
                lat=lat_gh,
                left=fortran_aligned,
                right=harmonica_aligned,
                case=case_dir.name,
                left_name="GravMagSphere spectral",
                right_name="Harmonica",
                left_params_label=format_params_label(
                    fortran_used["lmax_used"], fortran_used["reg_lambda_used"], fortran_used["reg_power_used"]
                ),
                right_params_label=format_harmonica_params_label(
                    harmonica_used["harmonica_damping_used"], harmonica_used["harmonica_depth_m_used"]
                ),
                out_path=side_plot_path,
            )
            row = {
                "case": case_dir.name,
                "solver": "fortran_minus_harmonica_xyz",
                "status": "ok",
                "plot": str(plot_path),
                "side_by_side_plot": str(side_plot_path),
                "lmax_used": fortran_used["lmax_used"],
                "reg_lambda_used": fortran_used["reg_lambda_used"],
                "reg_power_used": fortran_used["reg_power_used"],
                "harmonica_damping_used": harmonica_used["harmonica_damping_used"],
                "harmonica_depth_factor_used": harmonica_used["harmonica_depth_factor_used"],
                "harmonica_depth_m_used": harmonica_used["harmonica_depth_m_used"],
                "harmonica_block_size_m_used": harmonica_used["harmonica_block_size_m_used"],
            }
            row.update(m_gh)
            rows.append(row)
        else:
            rows.append(
                {
                    "case": case_dir.name,
                    "solver": "fortran_minus_harmonica_xyz",
                    "status": "error: missing successful fortran/harmonica predictions",
                    "plot": "",
                    "side_by_side_plot": "",
                    "lmax_used": float("nan"),
                    "reg_lambda_used": float("nan"),
                    "reg_power_used": float("nan"),
                    "harmonica_damping_used": float("nan"),
                    "harmonica_depth_factor_used": float("nan"),
                    "harmonica_depth_m_used": float("nan"),
                    "harmonica_block_size_m_used": float("nan"),
                }
            )

    fieldnames = [
        "case",
        "solver",
        "status",
        "plot",
        "side_by_side_plot",
        "lmax_used",
        "reg_lambda_used",
        "reg_power_used",
        "harmonica_damping_used",
        "harmonica_depth_factor_used",
        "harmonica_depth_m_used",
        "harmonica_block_size_m_used",
    ] + [f"{s}_{c}" for c in XYZ_COMPONENTS for s in METRIC_STATS]

    report_csv.parent.mkdir(parents=True, exist_ok=True)
    with report_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    report_md.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# GravMagSphere Spectral vs Harmonica (XYZ)",
        "",
        f"- fortran mode: `{'auto' if fortran_auto_mode else 'manual-seeded'}`",
        f"- fortran seed (manual mode): lmax `{lmax}`, reg_lambda `{reg_lambda}`, reg_power `{reg_power}`",
        f"- harmonica damping grid: `{','.join(f'{v:g}' for v in harmonica_damping_grid)}`",
        f"- harmonica depth-factor grid (x mean grid spacing): `{','.join(f'{v:g}' for v in harmonica_depth_factor_grid)}`",
        f"- harmonica block factor (x mean grid spacing): `{harmonica_block_factor:g}`",
        f"- residual plots: `{residual_dir}`",
        f"- GravMagSphere vs Harmonica side-by-side plots: `{side_by_side_dir}`",
        "",
        "| Case | Solver | Status | lmax_used | reg_lambda_used | reg_power_used | damping | depth_m | block_m | RMSE(Bx) | RMSE(By) | RMSE(Bz) | RMSE(Btot) | R2(Btot) | Residual plot | Side-by-side plot |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        if str(row["status"]) != "ok":
            lines.append(
                f"| {row['case']} | {row['solver']} | {row['status']} | - | - | - | - | - | - | - | - | - | - | - | - |"
            )
        else:
            lines.append(
                "| {case} | {solver} | {status} | {lmax_used:.0f} | {reg_lambda_used:.6g} | {reg_power_used:.6g} | "
                "{harmonica_damping_used:.6g} | {harmonica_depth_m_used:.6g} | {harmonica_block_size_m_used:.6g} | "
                "{rmse_bx:.6e} | {rmse_by:.6e} | {rmse_bz:.6e} | {rmse_btot:.6e} | {r2_btot:.6f} | {plot} | {side_by_side_plot} |".format(
                    **row
                )
            )
    lines.append("")
    lines.append("Residual definitions:")
    lines.append("- `fortran_spectral_vs_direct_xyz`: `spectral - direct_baseline`.")
    lines.append("- `harmonica_eqs_vs_direct_xyz`: `harmonica - direct_baseline`.")
    lines.append("- `fortran_minus_harmonica_xyz`: `spectral - harmonica`.")
    report_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_repo_catalog(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = """# External Solver Repository Candidates

These repositories are the closest external counterparts for spherical harmonic
or equivalent-source workflows compatible with GravMag Sphere outputs.

## 1) SHTOOLS / pyshtools
- GitHub: https://github.com/SHTOOLS/SHTOOLS
- Docs: https://shtools.github.io/SHTOOLS/
- Relevant APIs:
  - https://shtools.github.io/SHTOOLS/pyshexpandlsq.html
  - https://shtools.github.io/SHTOOLS/pylsq_g.html

## 2) Harmonica (Fatiando a Terra)
- GitHub: https://github.com/fatiando/harmonica
- Docs: https://www.fatiando.org/harmonica/latest/

## 3) SimPEG (potential fields)
- GitHub: https://github.com/simpeg/simpeg
- Docs: https://docs.simpeg.xyz/latest/content/api/simpeg.potential_fields.html
- Note: SimPEG is a broader forward/inversion framework rather than a point_source
  drop-in equivalent-source interpolator.
"""
    path.write_text(text, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="External solver formatting/comparison for GravMag Sphere lunar_examples."
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_repos = sub.add_parser("repos", help="Write markdown file with candidate external solver repositories.")
    p_repos.add_argument(
        "--out",
        default="diagnostics/external_solver_repos.md",
        help="Output markdown path, relative to repository root.",
    )

    p_fmt = sub.add_parser("format", help="Format lunar_examples for external solver pipelines.")
    p_fmt.add_argument("--examples-dir", "--lunar_examples-dir", default="examples/lunar_examples")
    p_fmt.add_argument("--out-root", default="diagnostics/external_solver_inputs")
    p_fmt.add_argument("--rsphere-km", type=float, default=1737.4)
    p_fmt.add_argument("--refine-factor", type=int, default=2)
    p_fmt.add_argument("--examples", "--lunar_examples", nargs="*", default=[], help="Optional list of example filenames.")

    p_cmp = sub.add_parser("compare", help="Run comparisons across available solver backends.")
    p_cmp.add_argument("--formatted-root", default="diagnostics/external_solver_inputs")
    p_cmp.add_argument("--report-md", default="diagnostics/external_solver_comparison.md")
    p_cmp.add_argument("--report-csv", default="diagnostics/external_solver_comparison.csv")
    p_cmp.add_argument(
        "--residual-dir",
        default="diagnostics/external_solver_residuals",
        help="Directory for residual PNG plots.",
    )
    p_cmp.add_argument(
        "--side-by-side-dir",
        default="diagnostics/gravmag_vs_shtools_side_by_side",
        help="Directory for GravMagSphere vs SHTOOLS side-by-side solution PNG plots.",
    )
    p_cmp.add_argument("--examples", "--lunar_examples", nargs="*", default=[], help="Optional list of case directory names.")
    p_cmp.add_argument("--lmax", type=int, default=24, help="Fortran seed lmax (manual mode fallback).")
    p_cmp.add_argument(
        "--ridge-lambda",
        type=float,
        default=1.0e-8,
        help="Single-value scipy ridge parameter (used when scipy-ridge-grid is omitted).",
    )
    p_cmp.add_argument("--scipy-lmax-grid", default="12,18,24,30,36")
    p_cmp.add_argument("--scipy-ridge-grid", default="0,1e-10,1e-8,1e-6")
    p_cmp.add_argument("--shtools-lmax-grid", default="12,18,24,30,36")
    p_cmp.add_argument("--fortran-auto-mode", type=int, default=1, help="1=auto mode, 0=manual seeded mode.")
    p_cmp.add_argument("--rsphere-km", type=float, default=1737.4)
    p_cmp.add_argument("--refine-factor", type=int, default=2)
    p_cmp.add_argument("--ntheta-fit", type=int, default=72)
    p_cmp.add_argument("--nphi-fit", type=int, default=144)
    p_cmp.add_argument("--reg-lambda", type=float, default=0.2)
    p_cmp.add_argument("--reg-power", type=float, default=4.0)

    p_cmp_h = sub.add_parser(
        "compare-harmonica-xyz",
        help="Run GravMagSphere spectral vs Harmonica equivalent-source comparisons in XYZ components.",
    )
    p_cmp_h.add_argument("--formatted-root", default="diagnostics/external_solver_inputs")
    p_cmp_h.add_argument("--report-md", default="diagnostics/gravmag_vs_harmonica_xyz.md")
    p_cmp_h.add_argument("--report-csv", default="diagnostics/gravmag_vs_harmonica_xyz.csv")
    p_cmp_h.add_argument(
        "--residual-dir",
        default="diagnostics/gravmag_vs_harmonica_xyz_residuals",
        help="Directory for residual PNG plots.",
    )
    p_cmp_h.add_argument(
        "--side-by-side-dir",
        default="diagnostics/gravmag_vs_harmonica_xyz_side_by_side",
        help="Directory for GravMagSphere vs Harmonica side-by-side solution PNG plots.",
    )
    p_cmp_h.add_argument("--examples", "--lunar_examples", nargs="*", default=[], help="Optional list of case directory names.")
    p_cmp_h.add_argument("--fortran-auto-mode", type=int, default=1, help="1=auto mode, 0=manual seeded mode.")
    p_cmp_h.add_argument("--rsphere-km", type=float, default=1737.4)
    p_cmp_h.add_argument("--refine-factor", type=int, default=2)
    p_cmp_h.add_argument("--ntheta-fit", type=int, default=72)
    p_cmp_h.add_argument("--nphi-fit", type=int, default=144)
    p_cmp_h.add_argument("--lmax", type=int, default=24, help="Fortran seed lmax (manual mode fallback).")
    p_cmp_h.add_argument("--reg-lambda", type=float, default=0.2)
    p_cmp_h.add_argument("--reg-power", type=float, default=4.0)
    p_cmp_h.add_argument("--harmonica-damping-grid", default="0,1e-8,1e-6,1e-4")
    p_cmp_h.add_argument("--harmonica-depth-factor-grid", default="1,2,4")
    p_cmp_h.add_argument("--harmonica-block-factor", type=float, default=4.0)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd == "repos":
        out = Path(args.out)
        if not out.is_absolute():
            out = repo_root / out
        write_repo_catalog(out)
    elif args.cmd == "format":
        examples_dir = Path(args.examples_dir)
        if not examples_dir.is_absolute():
            examples_dir = repo_root / examples_dir
        out_root = Path(args.out_root)
        if not out_root.is_absolute():
            out_root = repo_root / out_root
        format_examples(
            examples_dir=examples_dir,
            out_root=out_root,
            rsphere_km=float(args.rsphere_km),
            refine_factor=int(args.refine_factor),
            selected=args.examples,
        )
    elif args.cmd == "compare":
        formatted_root = Path(args.formatted_root)
        if not formatted_root.is_absolute():
            formatted_root = repo_root / formatted_root
        report_md = Path(args.report_md)
        if not report_md.is_absolute():
            report_md = repo_root / report_md
        report_csv = Path(args.report_csv)
        if not report_csv.is_absolute():
            report_csv = repo_root / report_csv
        residual_dir = Path(args.residual_dir)
        if not residual_dir.is_absolute():
            residual_dir = repo_root / residual_dir
        side_by_side_dir = Path(args.side_by_side_dir)
        if not side_by_side_dir.is_absolute():
            side_by_side_dir = repo_root / side_by_side_dir
        scipy_lmax_grid = parse_list_ints(str(args.scipy_lmax_grid))
        shtools_lmax_grid = parse_list_ints(str(args.shtools_lmax_grid))
        if not scipy_lmax_grid:
            scipy_lmax_grid = [int(args.lmax)]
        if not shtools_lmax_grid:
            shtools_lmax_grid = [int(args.lmax)]
        if str(args.scipy_ridge_grid).strip():
            scipy_ridge_grid = parse_list_floats(str(args.scipy_ridge_grid))
        else:
            scipy_ridge_grid = [float(args.ridge_lambda)]
        if not scipy_ridge_grid:
            scipy_ridge_grid = [float(args.ridge_lambda)]
        run_cmd(["./build_gravmag_tools.sh"], cwd=example_root)
        compare_cases(
            formatted_root=formatted_root,
            report_md=report_md,
            report_csv=report_csv,
            residual_dir=residual_dir,
            side_by_side_dir=side_by_side_dir,
            lmax=int(args.lmax),
            ridge_lambda=float(args.ridge_lambda),
            scipy_lmax_grid=scipy_lmax_grid,
            scipy_ridge_grid=scipy_ridge_grid,
            shtools_lmax_grid=shtools_lmax_grid,
            fortran_auto_mode=(int(args.fortran_auto_mode) == 1),
            rsphere_km=float(args.rsphere_km),
            refine_factor=int(args.refine_factor),
            ntheta_fit=int(args.ntheta_fit),
            nphi_fit=int(args.nphi_fit),
            reg_lambda=float(args.reg_lambda),
            reg_power=float(args.reg_power),
            selected=args.examples,
        )
    elif args.cmd == "compare-harmonica-xyz":
        formatted_root = Path(args.formatted_root)
        if not formatted_root.is_absolute():
            formatted_root = repo_root / formatted_root
        report_md = Path(args.report_md)
        if not report_md.is_absolute():
            report_md = repo_root / report_md
        report_csv = Path(args.report_csv)
        if not report_csv.is_absolute():
            report_csv = repo_root / report_csv
        residual_dir = Path(args.residual_dir)
        if not residual_dir.is_absolute():
            residual_dir = repo_root / residual_dir
        side_by_side_dir = Path(args.side_by_side_dir)
        if not side_by_side_dir.is_absolute():
            side_by_side_dir = repo_root / side_by_side_dir
        damping_grid = parse_list_floats(str(args.harmonica_damping_grid))
        depth_factor_grid = parse_list_floats(str(args.harmonica_depth_factor_grid))
        if not damping_grid:
            damping_grid = [0.0]
        if not depth_factor_grid:
            depth_factor_grid = [1.0]
        run_cmd(["./build_gravmag_tools.sh"], cwd=example_root)
        compare_cases_harmonica_xyz(
            formatted_root=formatted_root,
            report_md=report_md,
            report_csv=report_csv,
            residual_dir=residual_dir,
            side_by_side_dir=side_by_side_dir,
            fortran_auto_mode=(int(args.fortran_auto_mode) == 1),
            rsphere_km=float(args.rsphere_km),
            refine_factor=int(args.refine_factor),
            ntheta_fit=int(args.ntheta_fit),
            nphi_fit=int(args.nphi_fit),
            lmax=int(args.lmax),
            reg_lambda=float(args.reg_lambda),
            reg_power=float(args.reg_power),
            harmonica_damping_grid=damping_grid,
            harmonica_depth_factor_grid=depth_factor_grid,
            harmonica_block_factor=float(args.harmonica_block_factor),
            selected=args.examples,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
