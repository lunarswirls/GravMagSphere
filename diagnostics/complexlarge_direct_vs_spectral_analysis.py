#!/usr/bin/env python3
"""Analyze complexlarge direct_sum vs spectral boundary visibility in XYZ components."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


COMPONENTS = ("bx", "by", "bz", "btot")


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

    accum: Dict[Tuple[float, float], np.ndarray] = {}
    for lo, la, bxi, byi, bzi in zip(lon, lat, bx, by, bz):
        key = (float(lo), float(la))
        if key not in accum:
            accum[key] = np.zeros(3, dtype=float)
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


def align_maps(
    a: Dict[str, np.ndarray], b: Dict[str, np.ndarray]
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    a_idx = {(float(a["lon_deg"][i]), float(a["lat_deg"][i])): i for i in range(a["lon_deg"].size)}
    b_idx = {(float(b["lon_deg"][i]), float(b["lat_deg"][i])): i for i in range(b["lon_deg"].size)}
    keys = sorted(set(a_idx.keys()) & set(b_idx.keys()))
    if not keys:
        raise RuntimeError("No overlapping lon/lat cells for alignment.")

    lon = np.array([k[0] for k in keys], dtype=float)
    lat = np.array([k[1] for k in keys], dtype=float)
    aa = {c: np.array([a[c][a_idx[k]] for k in keys], dtype=float) for c in COMPONENTS}
    bb = {c: np.array([b[c][b_idx[k]] for k in keys], dtype=float) for c in COMPONENTS}
    return lon, lat, aa, bb


def to_grid(lon: np.ndarray, lat: np.ndarray, val: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ulon = np.unique(np.round(lon, 6))
    ulat = np.unique(np.round(lat, 6))
    ulon.sort()
    ulat.sort()
    ilon = {x: i for i, x in enumerate(ulon)}
    ilat = {y: i for i, y in enumerate(ulat)}
    z = np.full((ulat.size, ulon.size), np.nan, dtype=float)
    for lo, la, vv in zip(lon, lat, val):
        z[ilat[float(np.round(la, 6))], ilon[float(np.round(lo, 6))]] = float(vv)
    if np.isnan(z).any():
        raise RuntimeError("Grid conversion produced missing cells.")
    return ulon, ulat, z


def metric_summary(ref: Dict[str, np.ndarray], pred: Dict[str, np.ndarray]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for c in COMPONENTS:
        r = pred[c] - ref[c]
        out[f"rmse_{c}"] = float(np.sqrt(np.mean(r * r)))
        out[f"mae_{c}"] = float(np.mean(np.abs(r)))
        out[f"max_abs_{c}"] = float(np.max(np.abs(r)))
    return out


def robust_limits(z1: np.ndarray, z2: np.ndarray, percentile: float = 99.5) -> Tuple[float, float]:
    vals = np.concatenate([z1.ravel(), z2.ravel()])
    lo = float(np.nanpercentile(vals, 100.0 - percentile))
    hi = float(np.nanpercentile(vals, percentile))
    if hi <= lo:
        s = max(abs(hi), 1.0e-12)
        return -s, s
    return lo, hi


def plot_side_by_side_fields(
    *,
    lon: np.ndarray,
    lat: np.ndarray,
    direct_sum: Dict[str, np.ndarray],
    spectral: Dict[str, np.ndarray],
    out_path: Path,
    title: str,
) -> None:
    import matplotlib

    if os.environ.get("DISPLAY", "") == "" and os.environ.get("MPLBACKEND", "") == "":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [("bx", "Bx"), ("by", "By"), ("bz", "Bz"), ("btot", "Btot")]
    fig, axes = plt.subplots(4, 2, figsize=(12, 14), constrained_layout=True)
    for i, (k, lbl) in enumerate(labels):
        ax_l = axes[i, 0]
        ax_r = axes[i, 1]
        ulon, ulat, z_d = to_grid(lon, lat, direct_sum[k])
        _, _, z_s = to_grid(lon, lat, spectral[k])
        vmin, vmax = robust_limits(z_d, z_s)
        im = ax_l.pcolormesh(ulon, ulat, z_d, shading="nearest", cmap="RdBu_r", vmin=vmin, vmax=vmax)
        ax_r.pcolormesh(ulon, ulat, z_s, shading="nearest", cmap="RdBu_r", vmin=vmin, vmax=vmax)
        ax_l.set_title(f"{lbl} direct_sum")
        ax_r.set_title(f"{lbl} spectral")
        ax_l.set_xlabel("Longitude [deg]")
        ax_r.set_xlabel("Longitude [deg]")
        ax_l.set_ylabel("Latitude [deg]")
        ax_l.set_aspect("equal", adjustable="box")
        ax_r.set_aspect("equal", adjustable="box")
        cbar = fig.colorbar(im, ax=[ax_l, ax_r], orientation="horizontal", pad=0.1, shrink=0.95)
        cbar.set_label(lbl)
    fig.suptitle(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def plot_residuals(
    *,
    lon: np.ndarray,
    lat: np.ndarray,
    ref: Dict[str, np.ndarray],
    pred: Dict[str, np.ndarray],
    out_path: Path,
    title: str,
    grad_mask: np.ndarray | None = None,
) -> None:
    import matplotlib

    if os.environ.get("DISPLAY", "") == "" and os.environ.get("MPLBACKEND", "") == "":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [("bx", "Bx"), ("by", "By"), ("bz", "Bz"), ("btot", "Btot")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    axes = axes.ravel()
    for ax, (k, lbl) in zip(axes, labels):
        ulon, ulat, z = to_grid(lon, lat, pred[k] - ref[k])
        vmax = float(np.nanpercentile(np.abs(z), 99.5))
        vmax = max(vmax, 1.0e-12)
        im = ax.pcolormesh(ulon, ulat, z, shading="nearest", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(f"{lbl} residual")
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        ax.set_aspect("equal", adjustable="box")
        if grad_mask is not None and k == "btot":
            ax.contour(ulon, ulat, grad_mask.astype(float), levels=[0.5], colors="k", linewidths=0.6)
        cbar = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.12, shrink=0.85)
        cbar.set_label("Pred - direct_sum")
    fig.suptitle(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def compute_gradient_mask(lon: np.ndarray, lat: np.ndarray, btot: np.ndarray, percentile: float = 90.0) -> np.ndarray:
    ulon, ulat, z = to_grid(lon, lat, btot)
    dlon = float(np.median(np.diff(ulon))) if ulon.size > 1 else 1.0
    dlat = float(np.median(np.diff(ulat))) if ulat.size > 1 else 1.0
    gy, gx = np.gradient(z, dlat, dlon)
    grad = np.sqrt(gx * gx + gy * gy)
    threshold = float(np.nanpercentile(grad, percentile))
    return grad >= threshold


def summarize_boundary_contrast(
    lon: np.ndarray,
    lat: np.ndarray,
    direct_sum: Dict[str, np.ndarray],
    spectral: Dict[str, np.ndarray],
    spectral_nolocal: Dict[str, np.ndarray],
) -> Dict[str, float]:
    mask = compute_gradient_mask(lon, lat, direct_sum["btot"], percentile=90.0)
    _, _, r_s = to_grid(lon, lat, spectral["btot"] - direct_sum["btot"])
    _, _, r_n = to_grid(lon, lat, spectral_nolocal["btot"] - direct_sum["btot"])
    out: Dict[str, float] = {}
    for name, r in (("spectral", r_s), ("spectral_nolocal", r_n)):
        near = float(np.nanmean(np.abs(r[mask])))
        far = float(np.nanmean(np.abs(r[~mask])))
        out[f"{name}_mean_abs_residual_near_edge"] = near
        out[f"{name}_mean_abs_residual_far_edge"] = far
        out[f"{name}_near_to_far_ratio"] = near / far if far > 0.0 else float("nan")
    return out


def write_summary_md(
    *,
    out_path: Path,
    run_dir: Path,
    metrics_spec: Dict[str, float],
    metrics_spec_nolocal: Dict[str, float],
    metrics_diff_local: Dict[str, float],
    edge_summary: Dict[str, float],
    side_by_side_path: Path,
    residual_local_path: Path,
    residual_nolocal_path: Path,
) -> None:
    lines = [
        "# Complexlarge Direct summation vs Spectral Boundary Analysis",
        "",
        f"- run directory: `{run_dir}`",
        f"- side-by-side fields: `{side_by_side_path}`",
        f"- residuals (spectral - direct_sum): `{residual_local_path}`",
        f"- residuals (spectral no-local-correction - direct_sum): `{residual_nolocal_path}`",
        "",
        "## Metrics: Spectral - Direct summation",
        "",
        "| Component | RMSE | MAE | MaxAbs |",
        "|---|---:|---:|---:|",
    ]
    for c in COMPONENTS:
        lines.append(
            f"| {c} | {metrics_spec[f'rmse_{c}']:.6e} | {metrics_spec[f'mae_{c}']:.6e} | {metrics_spec[f'max_abs_{c}']:.6e} |"
        )
    lines += [
        "",
        "## Metrics: Spectral(no local correction) - Direct summation",
        "",
        "| Component | RMSE | MAE | MaxAbs |",
        "|---|---:|---:|---:|",
    ]
    for c in COMPONENTS:
        lines.append(
            f"| {c} | {metrics_spec_nolocal[f'rmse_{c}']:.6e} | {metrics_spec_nolocal[f'mae_{c}']:.6e} | {metrics_spec_nolocal[f'max_abs_{c}']:.6e} |"
        )
    lines += [
        "",
        "## Metrics: Spectral(local) - Spectral(no local correction)",
        "",
        "| Component | RMSE | MAE | MaxAbs |",
        "|---|---:|---:|---:|",
    ]
    for c in COMPONENTS:
        lines.append(
            f"| {c} | {metrics_diff_local[f'rmse_{c}']:.6e} | {metrics_diff_local[f'mae_{c}']:.6e} | {metrics_diff_local[f'max_abs_{c}']:.6e} |"
        )
    lines += [
        "",
        "## Edge-Region Contrast (using top-10% direct_sum Btot gradient as boundary proxy)",
        "",
        "| Model | MeanAbsResidual Near Edge | MeanAbsResidual Far From Edge | Near/Far Ratio |",
        "|---|---:|---:|---:|",
        (
            f"| spectral | {edge_summary['spectral_mean_abs_residual_near_edge']:.6e} | "
            f"{edge_summary['spectral_mean_abs_residual_far_edge']:.6e} | "
            f"{edge_summary['spectral_near_to_far_ratio']:.3f} |"
        ),
        (
            f"| spectral_nolocal | {edge_summary['spectral_nolocal_mean_abs_residual_near_edge']:.6e} | "
            f"{edge_summary['spectral_nolocal_mean_abs_residual_far_edge']:.6e} | "
            f"{edge_summary['spectral_nolocal_near_to_far_ratio']:.3f} |"
        ),
        "",
        "## Interpretation",
        "",
        (
            "The visible boundary is primarily a high-gradient edge effect from approximating a sharp polygonal source with a finite-order spectral model. "
            "Residuals concentrate near edge-like zones much more strongly than away from edges."
        ),
        (
            "Local edge correction changes the solution significantly near boundaries (see local-vs-nolocal metrics), "
            "which can make boundary contrast more visible even when broad ringing is reduced."
        ),
    ]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Complexlarge direct_sum vs spectral boundary analysis.")
    parser.add_argument(
        "--run-dir",
        default="diagnostics/complexlarge_direct_vs_spectral",
        help="Directory containing direct_sum/spectral xyz outputs (relative to repository root unless absolute).",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    run_dir = Path(args.run_dir)
    if not run_dir.is_absolute():
        run_dir = root / run_dir

    direct_path = run_dir / "complexlarge_direct_xyz.txt"
    spectral_path = run_dir / "complexlarge_spectral_xyz.txt"
    spectral_nolocal_path = run_dir / "complexlarge_spectral_nolocal_xyz.txt"

    direct_raw = load_xyz_aggregate(direct_path)
    spectral_raw = load_xyz_aggregate(spectral_path)
    spectral_nolocal_raw = load_xyz_aggregate(spectral_nolocal_path)

    lon, lat, direct_sum, spectral = align_maps(direct_raw, spectral_raw)
    _, _, _, spectral_nolocal = align_maps(direct_raw, spectral_nolocal_raw)
    _, _, spectral_local, spectral_nolocal_aligned = align_maps(spectral_raw, spectral_nolocal_raw)

    metrics_spec = metric_summary(direct_sum, spectral)
    metrics_spec_nolocal = metric_summary(direct_sum, spectral_nolocal)
    metrics_diff_local = metric_summary(spectral_local, spectral_nolocal_aligned)
    edge_summary = summarize_boundary_contrast(lon, lat, direct_sum, spectral, spectral_nolocal)

    grad_mask = compute_gradient_mask(lon, lat, direct_sum["btot"], percentile=90.0)

    side_by_side_path = run_dir / "complexlarge_direct_vs_spectral_side_by_side.png"
    residual_local_path = run_dir / "complexlarge_spectral_minus_direct_residuals.png"
    residual_nolocal_path = run_dir / "complexlarge_spectral_nolocal_minus_direct_residuals.png"
    summary_md = run_dir / "complexlarge_direct_vs_spectral_findings.md"

    plot_side_by_side_fields(
        lon=lon,
        lat=lat,
        direct_sum=direct_sum,
        spectral=spectral,
        out_path=side_by_side_path,
        title="Complexlarge: GravMagSphere direct_sum vs spectral (Bx/By/Bz/Btot)",
    )
    plot_residuals(
        lon=lon,
        lat=lat,
        ref=direct_sum,
        pred=spectral,
        out_path=residual_local_path,
        title="Complexlarge: spectral - direct_sum residuals (with local edge correction)",
        grad_mask=grad_mask,
    )
    plot_residuals(
        lon=lon,
        lat=lat,
        ref=direct_sum,
        pred=spectral_nolocal,
        out_path=residual_nolocal_path,
        title="Complexlarge: spectral(no local correction) - direct_sum residuals",
        grad_mask=grad_mask,
    )
    write_summary_md(
        out_path=summary_md,
        run_dir=run_dir,
        metrics_spec=metrics_spec,
        metrics_spec_nolocal=metrics_spec_nolocal,
        metrics_diff_local=metrics_diff_local,
        edge_summary=edge_summary,
        side_by_side_path=side_by_side_path,
        residual_local_path=residual_local_path,
        residual_nolocal_path=residual_nolocal_path,
    )

    print(f"Wrote: {side_by_side_path}")
    print(f"Wrote: {residual_local_path}")
    print(f"Wrote: {residual_nolocal_path}")
    print(f"Wrote: {summary_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
