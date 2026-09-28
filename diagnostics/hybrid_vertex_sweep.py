#!/usr/bin/env python3
"""
Run a synthetic-geometry sweep to determine when the hybrid GravMagSphere
solver should be preferred over spectral-only evaluation.

For each generated polygon case (increasing vertex count), the script runs:
1) direct_sum solver (reference)
2) spectral solver with hybrid disabled
3) spectral solver with hybrid in auto mode

Outputs:
- diagnostics/hybrid_vertex_sweep/hybrid_vertex_sweep_detailed.csv
- diagnostics/hybrid_vertex_sweep/hybrid_vertex_sweep_by_vertex.csv
- diagnostics/hybrid_vertex_sweep/hybrid_vertex_sweep.md
- diagnostics/hybrid_vertex_sweep/hybrid_vertex_sweep.png
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np


repo_root = Path(__file__).resolve().parent.parent
example_root = repo_root / "examples"
DEFAULT_OUTDIR = repo_root / "diagnostics" / "hybrid_vertex_sweep"
GEN_EXAMPLES_DIR = repo_root / "diagnostics" / "_hybrid_vertex_sweep_examples"
WORK_DIR = repo_root / "diagnostics" / "_hybrid_vertex_sweep_workspace"
COMPONENTS = ("bx", "by", "bz", "btot")
RE_GAUSS = re.compile(
    r"lmax=\s*(?P<lmax>\d+)\s+reg_lambda=\s*(?P<reg_lambda>[-+0-9.eE]+)\s+reg_power=\s*(?P<reg_power>[-+0-9.eE]+)"
)
RE_HYBRID = re.compile(
    r"hybrid_mode=\s*(?P<mode>\d+)\s+input_complex=\s*(?P<input_complex>[TF])\s+"
    r"hybrid_active=\s*(?P<active>[01])\s+hybrid_band_deg=\s*(?P<band>[-+0-9.eE]+)\s+"
    r"hybrid_transition_deg=\s*(?P<transition>[-+0-9.eE]+)\s+hybrid_cells=\s*(?P<cells>\d+)"
)
RE_DIAG = re.compile(r"^DIAG\|gauss\|time\|([^=]+)=(.+)$")


@dataclass
class SolverRun:
    wall_s: float
    stdout: str
    xyz_path: Path
    data: Dict[str, np.ndarray]
    lmax_used: float
    reg_lambda_used: float
    reg_power_used: float
    input_complex: bool
    hybrid_active: bool
    hybrid_mode: int
    hybrid_cells: int
    diag_total_s: float
    diag_hybrid_s: float


def parse_list_ints(raw: str) -> List[int]:
    out: List[int] = []
    for part in raw.split(","):
        s = part.strip()
        if not s:
            continue
        out.append(int(s))
    return out


def parse_list_str(raw: str) -> List[str]:
    out: List[str] = []
    for part in raw.split(","):
        s = part.strip()
        if not s:
            continue
        out.append(s)
    return out


def run_cmd(cmd: Sequence[str], cwd: Path, env: Dict[str, str]) -> Tuple[float, str]:
    t0 = time.perf_counter()
    proc = subprocess.run(
        list(cmd),
        cwd=str(cwd),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    wall_s = time.perf_counter() - t0
    if proc.returncode != 0:
        joined = " ".join(cmd)
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {joined}\n"
            f"--- begin output ---\n{proc.stdout}\n--- end output ---"
        )
    return wall_s, proc.stdout


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


def align_xyz_maps(
    ref: Dict[str, np.ndarray],
    pred: Dict[str, np.ndarray],
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    ref_map: Dict[Tuple[float, float], Dict[str, float]] = {}
    pred_map: Dict[Tuple[float, float], Dict[str, float]] = {}
    for i in range(ref["lon_deg"].size):
        key = (float(ref["lon_deg"][i]), float(ref["lat_deg"][i]))
        ref_map[key] = {c: float(ref[c][i]) for c in COMPONENTS}
    for i in range(pred["lon_deg"].size):
        key = (float(pred["lon_deg"][i]), float(pred["lat_deg"][i]))
        pred_map[key] = {c: float(pred[c][i]) for c in COMPONENTS}

    keys = sorted(set(ref_map.keys()) & set(pred_map.keys()))
    if not keys:
        raise RuntimeError("No overlapping lon/lat samples.")

    ref_out: Dict[str, np.ndarray] = {}
    pred_out: Dict[str, np.ndarray] = {}
    for c in COMPONENTS:
        ref_out[c] = np.array([ref_map[k][c] for k in keys], dtype=float)
        pred_out[c] = np.array([pred_map[k][c] for k in keys], dtype=float)
    return ref_out, pred_out


def metric_values(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    residual = y_pred - y_true
    rmse = float(np.sqrt(np.mean(residual * residual)))
    mae = float(np.mean(np.abs(residual)))
    max_abs = float(np.max(np.abs(residual)))
    ss_res = float(np.sum(residual * residual))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    r2 = float("nan") if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
    return {"rmse": rmse, "mae": mae, "max_abs": max_abs, "r2": r2}


def component_metrics(
    ref: Dict[str, np.ndarray],
    pred: Dict[str, np.ndarray],
    suffix: str,
) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for c in COMPONENTS:
        m = metric_values(ref[c], pred[c])
        out[f"rmse_{c}_{suffix}"] = m["rmse"]
        out[f"mae_{c}_{suffix}"] = m["mae"]
        out[f"max_abs_{c}_{suffix}"] = m["max_abs"]
        out[f"r2_{c}_{suffix}"] = m["r2"]
    return out


def parse_header_and_diag(xyz_path: Path, stdout: str) -> Dict[str, float | int | bool]:
    out: Dict[str, float | int | bool] = {
        "lmax_used": float("nan"),
        "reg_lambda_used": float("nan"),
        "reg_power_used": float("nan"),
        "input_complex": False,
        "hybrid_active": False,
        "hybrid_mode": -1,
        "hybrid_cells": 0,
        "diag_total_s": float("nan"),
        "diag_hybrid_s": 0.0,
    }

    if xyz_path.exists():
        with xyz_path.open("r", encoding="utf-8") as f:
            for raw in f:
                s = raw.strip()
                if not s.startswith("#"):
                    break
                m = RE_GAUSS.search(s)
                if m:
                    out["lmax_used"] = float(m.group("lmax"))
                    out["reg_lambda_used"] = float(m.group("reg_lambda"))
                    out["reg_power_used"] = float(m.group("reg_power"))
                h = RE_HYBRID.search(s)
                if h:
                    out["hybrid_mode"] = int(h.group("mode"))
                    out["input_complex"] = (h.group("input_complex") == "T")
                    out["hybrid_active"] = (int(h.group("active")) == 1)
                    out["hybrid_cells"] = int(h.group("cells"))

    for line in stdout.splitlines():
        d = RE_DIAG.match(line.strip())
        if not d:
            continue
        k, v = d.groups()
        try:
            fv = float(v.strip())
        except ValueError:
            continue
        if k == "body_total_s":
            out["diag_total_s"] = fv
        elif k == "hybrid_direct_eval_s":
            out["diag_hybrid_s"] = fv

    return out


def generate_polygon(
    n_vertices: int,
    shape: str,
    center_lon_deg: float,
    center_lat_deg: float,
    radius_deg: float,
    waviness: float,
    waviness_harmonic: int,
) -> Tuple[np.ndarray, np.ndarray]:
    ang = np.linspace(0.0, 2.0 * math.pi, n_vertices, endpoint=False)
    if shape == "regular":
        r = np.full_like(ang, radius_deg)
    elif shape == "wavy":
        harm = max(3, waviness_harmonic)
        r = radius_deg * (1.0 + waviness * np.cos(float(harm) * ang))
    else:
        raise ValueError(f"Unsupported shape '{shape}'")

    lon = center_lon_deg + r * np.cos(ang)
    lat = center_lat_deg + r * np.sin(ang)

    lon = np.concatenate([lon, lon[:1]])
    lat = np.concatenate([lat, lat[:1]])
    return lat, lon


def write_example_input(
    out_path: Path,
    *,
    title: str,
    lat0_deg: float,
    lon0_deg: float,
    dlat_deg: float,
    dlon_deg: float,
    nlat: int,
    nlon: int,
    nr: int,
    ntheta: int,
    nphi: int,
    mag_amp_apm: float,
    inc_deg: float,
    dec_deg: float,
    depth_top_km: float,
    depth_bot_km: float,
    poly_lat: np.ndarray,
    poly_lon: np.ndarray,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write(f"{title}\n")
        f.write(f"{lat0_deg:.6f} {lon0_deg:.6f} {dlat_deg:.6f} {dlon_deg:.6f} 0.0 {nlat:d} {nlon:d}\n")
        f.write(f"{nr:d} {ntheta:d} {nphi:d} 0\n")
        f.write("0.1 1.50 -2.0 2.0 0.080\n")
        f.write(f"2 1 {mag_amp_apm:.6f} {inc_deg:.6f} {dec_deg:.6f}\n")
        f.write("0 1 1 0.0 0.0\n")
        f.write(f"{poly_lat.size:d} {depth_top_km:.6f} {depth_bot_km:.6f}\n")
        for la, lo in zip(poly_lat, poly_lon):
            f.write(f"{la: .6f} {lo: .6f}\n")


def build_once(env: Dict[str, str]) -> None:
    run_cmd(["./build_gravmag_tools.sh"], cwd=example_root, env=env)


def run_direct(
    *,
    example_path: Path,
    output_xyz: Path,
    rsphere_km: float,
    refine_factor: int,
    source_nlat: int,
    source_nlon: int,
    source_nr: int,
    env: Dict[str, str],
) -> SolverRun:
    cmd = [
        "./run_input_to_xyz.sh",
        "direct_sum",
        f"{rsphere_km}",
        str(example_path),
        str(output_xyz),
        "--refine-factor",
        f"{refine_factor}",
        "--source-nlat",
        f"{source_nlat}",
        "--source-nlon",
        f"{source_nlon}",
        "--source-nr",
        f"{source_nr}",
    ]
    wall_s, stdout = run_cmd(cmd, cwd=example_root, env=env)
    data = load_xyz_aggregate(output_xyz)
    return SolverRun(
        wall_s=wall_s,
        stdout=stdout,
        xyz_path=output_xyz,
        data=data,
        lmax_used=float("nan"),
        reg_lambda_used=float("nan"),
        reg_power_used=float("nan"),
        input_complex=False,
        hybrid_active=False,
        hybrid_mode=-1,
        hybrid_cells=0,
        diag_total_s=float("nan"),
        diag_hybrid_s=0.0,
    )


def run_spectral(
    *,
    example_path: Path,
    output_xyz: Path,
    rsphere_km: float,
    refine_factor: int,
    lmax_seed: int,
    ntheta_fit: int,
    nphi_fit: int,
    reg_lambda_seed: float,
    reg_power_seed: float,
    source_nlat: int,
    source_nlon: int,
    source_nr: int,
    auto_mode: int,
    joint_strength: float,
    edge_correction: int,
    hybrid_mode: int,
    hybrid_band_deg: float,
    complex_vertex_threshold: int,
    hybrid_transition_deg: float,
    env: Dict[str, str],
) -> SolverRun:
    cmd = [
        "./run_input_to_xyz.sh",
        "spectral",
        f"{rsphere_km}",
        str(example_path),
        str(output_xyz),
        "--refine-factor",
        f"{refine_factor}",
        "--lmax",
        f"{lmax_seed}",
        "--ntheta-fit",
        f"{ntheta_fit}",
        "--nphi-fit",
        f"{nphi_fit}",
        "--reg-lambda",
        f"{reg_lambda_seed}",
        "--reg-power",
        f"{reg_power_seed}",
        "--source-nlat",
        f"{source_nlat}",
        "--source-nlon",
        f"{source_nlon}",
        "--source-nr",
        f"{source_nr}",
        "--auto-mode",
        f"{auto_mode}",
        "--joint-strength",
        f"{joint_strength}",
        "--edge-correction",
        f"{edge_correction}",
        "--hybrid-mode",
        f"{hybrid_mode}",
        "--hybrid-band-deg",
        f"{hybrid_band_deg}",
        "--complex-vertex-threshold",
        f"{complex_vertex_threshold}",
        "--hybrid-transition-deg",
        f"{hybrid_transition_deg}",
    ]
    wall_s, stdout = run_cmd(cmd, cwd=example_root, env=env)
    data = load_xyz_aggregate(output_xyz)
    meta = parse_header_and_diag(output_xyz, stdout)
    return SolverRun(
        wall_s=wall_s,
        stdout=stdout,
        xyz_path=output_xyz,
        data=data,
        lmax_used=float(meta["lmax_used"]),
        reg_lambda_used=float(meta["reg_lambda_used"]),
        reg_power_used=float(meta["reg_power_used"]),
        input_complex=bool(meta["input_complex"]),
        hybrid_active=bool(meta["hybrid_active"]),
        hybrid_mode=int(meta["hybrid_mode"]),
        hybrid_cells=int(meta["hybrid_cells"]),
        diag_total_s=float(meta["diag_total_s"]),
        diag_hybrid_s=float(meta["diag_hybrid_s"]),
    )


def ensure_plot(
    out_png: Path,
    summary_rows: List[Dict[str, float | int | str]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_shape: Dict[str, Dict[str, List[Tuple[float, float]]]] = {}
    for row in summary_rows:
        shape = str(row["shape"])
        if shape not in by_shape:
            by_shape[shape] = {
                "spec_rmse": [],
                "hyb_rmse": [],
                "spec_wall": [],
                "hyb_wall": [],
                "impr": [],
            }
        v = float(row["n_vertices"])
        by_shape[shape]["spec_rmse"].append((v, float(row["rmse_btot_spec"])))
        by_shape[shape]["hyb_rmse"].append((v, float(row["rmse_btot_hybrid"])))
        by_shape[shape]["spec_wall"].append((v, float(row["wall_s_spec"])))
        by_shape[shape]["hyb_wall"].append((v, float(row["wall_s_hybrid"])))
        by_shape[shape]["impr"].append((v, 100.0 * float(row["rmse_improvement_frac"])))

    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    for shape, vals in sorted(by_shape.items()):
        for key in vals:
            vals[key].sort(key=lambda t: t[0])
        x = [v for v, _ in vals["spec_rmse"]]
        axs[0].plot(x, [y for _, y in vals["spec_rmse"]], "o-", label=f"{shape} spectral")
        axs[0].plot(x, [y for _, y in vals["hyb_rmse"]], "s--", label=f"{shape} hybrid")
        axs[1].plot(x, [y for _, y in vals["spec_wall"]], "o-", label=f"{shape} spectral")
        axs[1].plot(x, [y for _, y in vals["hyb_wall"]], "s--", label=f"{shape} hybrid")
        axs[2].plot(x, [y for _, y in vals["impr"]], "o-", label=shape)

    axs[0].set_title("Btot RMSE vs vertex count")
    axs[0].set_xlabel("vertices")
    axs[0].set_ylabel("RMSE (nT)")
    axs[0].grid(True, alpha=0.3)

    axs[1].set_title("Wall time vs vertex count")
    axs[1].set_xlabel("vertices")
    axs[1].set_ylabel("wall time (s)")
    axs[1].grid(True, alpha=0.3)

    axs[2].set_title("Hybrid RMSE improvement (%)")
    axs[2].set_xlabel("vertices")
    axs[2].set_ylabel("100 * (RMSEspec-RMSEhyb)/RMSEspec")
    axs[2].axhline(0.0, color="k", lw=1.0, alpha=0.4)
    axs[2].grid(True, alpha=0.3)

    handles, labels = axs[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, ncol=2, loc="upper center")
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.93])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def choose_threshold(
    by_vertex_rows: List[Dict[str, float | int]],
    *,
    min_improvement_frac: float,
    max_runtime_ratio: float,
    min_fraction_better: float,
    runtime_penalty: float,
) -> int:
    if not by_vertex_rows:
        return 12

    sorted_rows = sorted(by_vertex_rows, key=lambda r: int(r["n_vertices"]))
    viable = [
        int(r["n_vertices"])
        for r in sorted_rows
        if float(r["mean_rmse_improvement_frac"]) >= min_improvement_frac
        and float(r["mean_runtime_ratio_hybrid_over_spec"]) <= max_runtime_ratio
        and float(r["fraction_cases_hybrid_better"]) >= min_fraction_better
    ]
    if viable:
        return min(viable)

    best_score = -1.0e30
    best_v = int(sorted_rows[0]["n_vertices"])
    for r in sorted_rows:
        imp = float(r["mean_rmse_improvement_frac"])
        run_ratio = float(r["mean_runtime_ratio_hybrid_over_spec"])
        score = imp - runtime_penalty * max(0.0, run_ratio - 1.0)
        if score > best_score:
            best_score = score
            best_v = int(r["n_vertices"])
    return best_v


def main() -> None:
    parser = argparse.ArgumentParser(description="Hybrid solver vertex sweep.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTDIR))
    parser.add_argument("--shapes", default="regular,wavy", help="Comma-separated: regular,wavy")
    parser.add_argument("--vertices", default="4,6,8,10,12,14,16,20,24,30,36,48")
    parser.add_argument("--rsphere-km", type=float, default=1737.4)
    parser.add_argument("--refine-factor", type=int, default=2)
    parser.add_argument("--lmax", type=int, default=18)
    parser.add_argument("--ntheta-fit", type=int, default=48)
    parser.add_argument("--nphi-fit", type=int, default=96)
    parser.add_argument("--reg-lambda", type=float, default=0.2)
    parser.add_argument("--reg-power", type=float, default=4.0)
    parser.add_argument("--source-nlat", type=int, default=12)
    parser.add_argument("--source-nlon", type=int, default=12)
    parser.add_argument("--source-nr", type=int, default=2)
    parser.add_argument("--auto-mode", type=int, default=0, choices=(0, 1))
    parser.add_argument("--joint-strength", type=float, default=1.0)
    parser.add_argument("--edge-correction", type=int, default=1, choices=(0, 1))
    parser.add_argument("--hybrid-band-deg", type=float, default=1.5)
    parser.add_argument("--hybrid-transition-deg", type=float, default=0.75)
    parser.add_argument("--complex-vertex-threshold", type=int, default=12)
    parser.add_argument("--radius-deg", type=float, default=2.0)
    parser.add_argument("--waviness", type=float, default=0.30)
    parser.add_argument("--waviness-harmonic", type=int, default=5)
    parser.add_argument("--obs-lat0", type=float, default=-6.0)
    parser.add_argument("--obs-lon0", type=float, default=-6.0)
    parser.add_argument("--obs-dlat", type=float, default=0.2)
    parser.add_argument("--obs-dlon", type=float, default=0.2)
    parser.add_argument("--obs-nlat", type=int, default=61)
    parser.add_argument("--obs-nlon", type=int, default=61)
    parser.add_argument("--inc-deg", type=float, default=30.0)
    parser.add_argument("--dec-deg", type=float, default=210.0)
    parser.add_argument("--mag-amp", type=float, default=0.25)
    parser.add_argument("--depth-top-km", type=float, default=0.5)
    parser.add_argument("--depth-bot-km", type=float, default=8.0)
    parser.add_argument("--min-improvement-frac", type=float, default=0.10)
    parser.add_argument("--max-runtime-ratio", type=float, default=2.5)
    parser.add_argument("--min-fraction-better", type=float, default=0.5)
    parser.add_argument("--runtime-penalty", type=float, default=0.1)
    parser.add_argument("--skip-build", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    if not out_dir.is_absolute():
        out_dir = repo_root / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    GEN_EXAMPLES_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["SKIP_BUILD"] = "1"
    env["GRAVMAG_DIAGNOSTICS"] = "1"
    if "MPLCONFIGDIR" not in env:
        mpl = repo_root / "diagnostics" / ".mplconfig"
        mpl.mkdir(parents=True, exist_ok=True)
        env["MPLCONFIGDIR"] = str(mpl)

    if not args.skip_build:
        env_build = env.copy()
        env_build["SKIP_BUILD"] = "0"
        build_once(env_build)

    shapes = parse_list_str(args.shapes)
    vertices = sorted(set(parse_list_ints(args.vertices)))
    if not shapes:
        raise RuntimeError("No shapes configured.")
    if not vertices:
        raise RuntimeError("No vertex counts configured.")

    detailed_rows: List[Dict[str, float | int | str]] = []
    summary_rows: List[Dict[str, float | int | str]] = []
    by_vertex_rows: List[Dict[str, float | int]] = []

    for shape in shapes:
        for n_vertices in vertices:
            case_id = f"{shape}_v{n_vertices:03d}"
            case_workspace = WORK_DIR / case_id
            case_workspace.mkdir(parents=True, exist_ok=True)
            example_path = GEN_EXAMPLES_DIR / f"{case_id}.in"

            poly_lat, poly_lon = generate_polygon(
                n_vertices=n_vertices,
                shape=shape,
                center_lon_deg=0.0,
                center_lat_deg=0.0,
                radius_deg=float(args.radius_deg),
                waviness=float(args.waviness),
                waviness_harmonic=int(args.waviness_harmonic),
            )
            write_example_input(
                example_path,
                title=f"Hybrid sweep case: {shape}, vertices={n_vertices}",
                lat0_deg=float(args.obs_lat0),
                lon0_deg=float(args.obs_lon0),
                dlat_deg=float(args.obs_dlat),
                dlon_deg=float(args.obs_dlon),
                nlat=int(args.obs_nlat),
                nlon=int(args.obs_nlon),
                nr=int(args.source_nr),
                ntheta=int(args.source_nlat),
                nphi=int(args.source_nlon),
                mag_amp_apm=float(args.mag_amp),
                inc_deg=float(args.inc_deg),
                dec_deg=float(args.dec_deg),
                depth_top_km=float(args.depth_top_km),
                depth_bot_km=float(args.depth_bot_km),
                poly_lat=poly_lat,
                poly_lon=poly_lon,
            )

            direct_xyz = case_workspace / "direct_xyz.txt"
            spec_xyz = case_workspace / "spectral_nohybrid_xyz.txt"
            hyb_xyz = case_workspace / "spectral_hybrid_auto_xyz.txt"

            direct_sum = run_direct(
                example_path=example_path,
                output_xyz=direct_xyz,
                rsphere_km=float(args.rsphere_km),
                refine_factor=int(args.refine_factor),
                source_nlat=int(args.source_nlat),
                source_nlon=int(args.source_nlon),
                source_nr=int(args.source_nr),
                env=env,
            )
            spec = run_spectral(
                example_path=example_path,
                output_xyz=spec_xyz,
                rsphere_km=float(args.rsphere_km),
                refine_factor=int(args.refine_factor),
                lmax_seed=int(args.lmax),
                ntheta_fit=int(args.ntheta_fit),
                nphi_fit=int(args.nphi_fit),
                reg_lambda_seed=float(args.reg_lambda),
                reg_power_seed=float(args.reg_power),
                source_nlat=int(args.source_nlat),
                source_nlon=int(args.source_nlon),
                source_nr=int(args.source_nr),
                auto_mode=int(args.auto_mode),
                joint_strength=float(args.joint_strength),
                edge_correction=int(args.edge_correction),
                hybrid_mode=0,
                hybrid_band_deg=float(args.hybrid_band_deg),
                complex_vertex_threshold=int(args.complex_vertex_threshold),
                hybrid_transition_deg=float(args.hybrid_transition_deg),
                env=env,
            )
            hybrid = run_spectral(
                example_path=example_path,
                output_xyz=hyb_xyz,
                rsphere_km=float(args.rsphere_km),
                refine_factor=int(args.refine_factor),
                lmax_seed=int(args.lmax),
                ntheta_fit=int(args.ntheta_fit),
                nphi_fit=int(args.nphi_fit),
                reg_lambda_seed=float(args.reg_lambda),
                reg_power_seed=float(args.reg_power),
                source_nlat=int(args.source_nlat),
                source_nlon=int(args.source_nlon),
                source_nr=int(args.source_nr),
                auto_mode=int(args.auto_mode),
                joint_strength=float(args.joint_strength),
                edge_correction=int(args.edge_correction),
                hybrid_mode=1,
                hybrid_band_deg=float(args.hybrid_band_deg),
                complex_vertex_threshold=int(args.complex_vertex_threshold),
                hybrid_transition_deg=float(args.hybrid_transition_deg),
                env=env,
            )

            ref_spec, pred_spec = align_xyz_maps(direct_sum.data, spec.data)
            ref_hyb, pred_hyb = align_xyz_maps(direct_sum.data, hybrid.data)

            m_spec = component_metrics(ref_spec, pred_spec, "spec")
            m_hyb = component_metrics(ref_hyb, pred_hyb, "hybrid")

            rmse_spec = float(m_spec["rmse_btot_spec"])
            rmse_hyb = float(m_hyb["rmse_btot_hybrid"])
            rmse_improvement_frac = 0.0
            if rmse_spec > 0.0:
                rmse_improvement_frac = (rmse_spec - rmse_hyb) / rmse_spec

            runtime_ratio = hybrid.wall_s / max(spec.wall_s, 1.0e-12)
            hybrid_better = (
                rmse_improvement_frac >= float(args.min_improvement_frac)
                and runtime_ratio <= float(args.max_runtime_ratio)
            )

            detailed_rows.append(
                {
                    "case": case_id,
                    "shape": shape,
                    "n_vertices": n_vertices,
                    "solver": "spectral_nohybrid",
                    "wall_s": spec.wall_s,
                    "diag_total_s": spec.diag_total_s,
                    "diag_hybrid_s": spec.diag_hybrid_s,
                    "lmax_used": spec.lmax_used,
                    "reg_lambda_used": spec.reg_lambda_used,
                    "reg_power_used": spec.reg_power_used,
                    "input_complex": int(spec.input_complex),
                    "hybrid_active": int(spec.hybrid_active),
                    "hybrid_mode": spec.hybrid_mode,
                    "hybrid_cells": spec.hybrid_cells,
                    **m_spec,
                }
            )
            detailed_rows.append(
                {
                    "case": case_id,
                    "shape": shape,
                    "n_vertices": n_vertices,
                    "solver": "spectral_hybrid_auto",
                    "wall_s": hybrid.wall_s,
                    "diag_total_s": hybrid.diag_total_s,
                    "diag_hybrid_s": hybrid.diag_hybrid_s,
                    "lmax_used": hybrid.lmax_used,
                    "reg_lambda_used": hybrid.reg_lambda_used,
                    "reg_power_used": hybrid.reg_power_used,
                    "input_complex": int(hybrid.input_complex),
                    "hybrid_active": int(hybrid.hybrid_active),
                    "hybrid_mode": hybrid.hybrid_mode,
                    "hybrid_cells": hybrid.hybrid_cells,
                    **m_hyb,
                }
            )
            summary_rows.append(
                {
                    "case": case_id,
                    "shape": shape,
                    "n_vertices": n_vertices,
                    "input_complex": int(hybrid.input_complex),
                    "hybrid_active": int(hybrid.hybrid_active),
                    "rmse_btot_spec": rmse_spec,
                    "rmse_btot_hybrid": rmse_hyb,
                    "rmse_improvement_frac": rmse_improvement_frac,
                    "runtime_ratio_hybrid_over_spec": runtime_ratio,
                    "hybrid_better": int(hybrid_better),
                    "wall_s_spec": spec.wall_s,
                    "wall_s_hybrid": hybrid.wall_s,
                    "lmax_used": hybrid.lmax_used,
                    "reg_lambda_used": hybrid.reg_lambda_used,
                    "reg_power_used": hybrid.reg_power_used,
                }
            )

            print(
                f"[hybrid-vertex-sweep] {case_id}: "
                f"RMSEspec={rmse_spec:.4e}, RMSEhyb={rmse_hyb:.4e}, "
                f"impr={100.0 * rmse_improvement_frac:+.2f}%, "
                f"runtime_ratio={runtime_ratio:.2f}, "
                f"input_complex={int(hybrid.input_complex)}, hybrid_active={int(hybrid.hybrid_active)}"
            )

    vertices_sorted = sorted(set(int(r["n_vertices"]) for r in summary_rows))
    for n_vertices in vertices_sorted:
        rows_v = [r for r in summary_rows if int(r["n_vertices"]) == n_vertices]
        mean_impr = float(np.mean([float(r["rmse_improvement_frac"]) for r in rows_v]))
        mean_runtime = float(np.mean([float(r["runtime_ratio_hybrid_over_spec"]) for r in rows_v]))
        frac_better = float(np.mean([float(r["hybrid_better"]) for r in rows_v]))
        frac_complex = float(np.mean([float(r["input_complex"]) for r in rows_v]))
        by_vertex_rows.append(
            {
                "n_vertices": n_vertices,
                "mean_rmse_improvement_frac": mean_impr,
                "mean_runtime_ratio_hybrid_over_spec": mean_runtime,
                "fraction_cases_hybrid_better": frac_better,
                "fraction_cases_classified_complex": frac_complex,
            }
        )

    recommended_threshold = choose_threshold(
        by_vertex_rows,
        min_improvement_frac=float(args.min_improvement_frac),
        max_runtime_ratio=float(args.max_runtime_ratio),
        min_fraction_better=float(args.min_fraction_better),
        runtime_penalty=float(args.runtime_penalty),
    )

    detailed_csv = out_dir / "hybrid_vertex_sweep_detailed.csv"
    by_vertex_csv = out_dir / "hybrid_vertex_sweep_by_vertex.csv"
    summary_md = out_dir / "hybrid_vertex_sweep.md"
    summary_png = out_dir / "hybrid_vertex_sweep.png"

    detailed_fields: List[str] = [
        "case",
        "shape",
        "n_vertices",
        "solver",
        "wall_s",
        "diag_total_s",
        "diag_hybrid_s",
        "lmax_used",
        "reg_lambda_used",
        "reg_power_used",
        "input_complex",
        "hybrid_active",
        "hybrid_mode",
        "hybrid_cells",
    ]
    for c in COMPONENTS:
        for stat in ("rmse", "mae", "max_abs", "r2"):
            detailed_fields.append(f"{stat}_{c}_spec")
            detailed_fields.append(f"{stat}_{c}_hybrid")

    with detailed_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=detailed_fields)
        w.writeheader()
        for r in detailed_rows:
            row = {k: r.get(k, "") for k in detailed_fields}
            w.writerow(row)

    with by_vertex_csv.open("w", newline="", encoding="utf-8") as f:
        fields = [
            "n_vertices",
            "mean_rmse_improvement_frac",
            "mean_runtime_ratio_hybrid_over_spec",
            "fraction_cases_hybrid_better",
            "fraction_cases_classified_complex",
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in sorted(by_vertex_rows, key=lambda x: int(x["n_vertices"])):
            w.writerow(r)

    ensure_plot(summary_png, summary_rows)

    lines: List[str] = []
    lines.append("# Hybrid Vertex Sweep")
    lines.append("")
    lines.append("## Setup")
    lines.append(f"- shapes: `{','.join(shapes)}`")
    lines.append(f"- vertices: `{','.join(str(v) for v in vertices_sorted)}`")
    lines.append(
        f"- spectral seeds: lmax `{args.lmax}`, reg_lambda `{args.reg_lambda}`, reg_power `{args.reg_power}`, auto_mode `{args.auto_mode}`"
    )
    lines.append(
        f"- hybrid params: band `{args.hybrid_band_deg}` deg, transition `{args.hybrid_transition_deg}` deg, "
        f"complex_vertex_threshold(test setting) `{args.complex_vertex_threshold}`"
    )
    lines.append(
        f"- recommendation criteria: improvement >= `{args.min_improvement_frac:.3f}`, "
        f"runtime_ratio <= `{args.max_runtime_ratio:.3f}`, fraction_better >= `{args.min_fraction_better:.3f}`"
    )
    lines.append("")
    lines.append("## Recommendation")
    lines.append(f"- recommended `complex_vertex_threshold`: **{recommended_threshold}**")
    lines.append("")
    lines.append("## By Vertex (mean across shapes)")
    lines.append("")
    lines.append(
        "| n_vertices | mean RMSE improvement frac | mean runtime ratio (hyb/spec) | fraction hybrid better | fraction classified complex |"
    )
    lines.append("|---:|---:|---:|---:|---:|")
    for r in sorted(by_vertex_rows, key=lambda x: int(x["n_vertices"])):
        lines.append(
            "| {n_vertices} | {mean_rmse_improvement_frac:.4f} | {mean_runtime_ratio_hybrid_over_spec:.4f} | "
            "{fraction_cases_hybrid_better:.4f} | {fraction_cases_classified_complex:.4f} |".format(**r)
        )
    lines.append("")
    lines.append("## Per Case")
    lines.append("")
    lines.append(
        "| case | input_complex | hybrid_active | RMSEspec(Btot) | RMSEhyb(Btot) | improvement frac | runtime ratio |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for r in sorted(summary_rows, key=lambda x: (str(x["shape"]), int(x["n_vertices"]))):
        lines.append(
            "| {case} | {input_complex} | {hybrid_active} | {rmse_btot_spec:.4e} | {rmse_btot_hybrid:.4e} | "
            "{rmse_improvement_frac:.4f} | {runtime_ratio_hybrid_over_spec:.4f} |".format(**r)
        )
    lines.append("")
    lines.append("## Files")
    lines.append(f"- detailed csv: `{detailed_csv}`")
    lines.append(f"- by-vertex csv: `{by_vertex_csv}`")
    lines.append(f"- plot: `{summary_png}`")

    summary_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[hybrid-vertex-sweep] wrote: {detailed_csv}")
    print(f"[hybrid-vertex-sweep] wrote: {by_vertex_csv}")
    print(f"[hybrid-vertex-sweep] wrote: {summary_png}")
    print(f"[hybrid-vertex-sweep] wrote: {summary_md}")
    print(f"[hybrid-vertex-sweep] recommended complex_vertex_threshold={recommended_threshold}")


if __name__ == "__main__":
    main()
