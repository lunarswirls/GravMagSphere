#!/usr/bin/env python3
"""
Sweep GravMag spectral solver parameters against SHTOOLS target
Builds operating map for runtime/accuracy tradeoffs

THIS IS A WORK IN PROGRESS!!!!!!!
Do not touch unless you know what you are doing!
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import time
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


repo_root = Path(__file__).resolve().parent.parent
example_root = repo_root / "examples"
COMPONENTS = ("br", "btheta", "bphi", "btot")
DIAG_RE = re.compile(r"^DIAG\|gauss\|time\|([^=]+)=(.+)$")


def parse_list_floats(raw: str) -> List[float]:
    out: List[float] = []
    for part in raw.split(","):
        s = part.strip()
        if not s:
            continue
        out.append(float(s))
    return out


def parse_list_ints(raw: str) -> List[int]:
    out: List[int] = []
    for part in raw.split(","):
        s = part.strip()
        if not s:
            continue
        out.append(int(s))
    return out


def run_cmd(cmd: Sequence[str], cwd: Path, env: Dict[str, str] | None = None) -> Tuple[float, str]:
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


def parse_gauss_diag_times(output: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for line in output.splitlines():
        m = DIAG_RE.match(line.strip())
        if not m:
            continue
        key, raw = m.groups()
        try:
            out[key] = float(raw.strip())
        except ValueError:
            continue
    return out


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


def align_components(
    ref: Dict[str, np.ndarray],
    pred: Dict[str, np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    ref_map: Dict[Tuple[float, float], Dict[str, float]] = {}
    pred_map: Dict[Tuple[float, float], Dict[str, float]] = {}
    for i in range(ref["lon_deg"].size):
        k = (float(ref["lon_deg"][i]), float(ref["lat_deg"][i]))
        ref_map[k] = {c: float(ref[c][i]) for c in COMPONENTS}
    for i in range(pred["lon_deg"].size):
        k = (float(pred["lon_deg"][i]), float(pred["lat_deg"][i]))
        pred_map[k] = {c: float(pred[c][i]) for c in COMPONENTS}

    keys = sorted(set(ref_map.keys()) & set(pred_map.keys()))
    if not keys:
        raise RuntimeError("No overlapping lon/lat samples.")

    lon = np.array([k[0] for k in keys], dtype=float)
    lat = np.array([k[1] for k in keys], dtype=float)
    ref_a = {c: np.array([ref_map[k][c] for k in keys], dtype=float) for c in COMPONENTS}
    pred_a = {c: np.array([pred_map[k][c] for k in keys], dtype=float) for c in COMPONENTS}
    return lon, lat, ref_a, pred_a


def metric_values(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    r = y_pred - y_true
    rmse = float(np.sqrt(np.mean(r * r)))
    mae = float(np.mean(np.abs(r)))
    max_abs = float(np.max(np.abs(r)))
    ss_res = float(np.sum(r * r))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    r2 = float("nan") if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
    return {"rmse": rmse, "mae": mae, "max_abs": max_abs, "r2": r2}


def component_metrics(ref: Dict[str, np.ndarray], pred: Dict[str, np.ndarray], suffix: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for c in COMPONENTS:
        m = metric_values(ref[c], pred[c])
        out[f"rmse_{c}_{suffix}"] = m["rmse"]
        out[f"mae_{c}_{suffix}"] = m["mae"]
        out[f"max_abs_{c}_{suffix}"] = m["max_abs"]
        out[f"r2_{c}_{suffix}"] = m["r2"]
    return out


def predict_shtools_components(obs: Dict[str, np.ndarray], lmax: int) -> Dict[str, np.ndarray]:
    import pyshtools as pysh  # type: ignore

    lat = obs["lat_deg"]
    lon = obs["lon_deg"]
    out: Dict[str, np.ndarray] = {"lon_deg": lon.copy(), "lat_deg": lat.copy()}

    def solve_comp(data: np.ndarray) -> np.ndarray:
        cilm, _chi2 = pysh.expand.SHExpandLSQ(data, lat, lon, lmax)
        g = pysh.expand.LSQ_G(lat, lon, lmax)
        coeffs: List[float] = []
        for l in range(0, lmax + 1):
            for m in range(0, l + 1):
                coeffs.append(float(cilm[0, l, m]))
            for m in range(1, l + 1):
                coeffs.append(float(cilm[1, l, m]))
        return g @ np.array(coeffs, dtype=float)

    out["br"] = solve_comp(obs["br"])
    out["btheta"] = solve_comp(obs["btheta"])
    out["bphi"] = solve_comp(obs["bphi"])
    out["btot"] = np.sqrt(out["br"] * out["br"] + out["btheta"] * out["btheta"] + out["bphi"] * out["bphi"])
    return out


def run_baseline_direct(
    *,
    workspace: Path,
    example: Path,
    rsphere_km: float,
    refine_factor: int,
) -> Dict[str, np.ndarray]:
    xyz = workspace / "baseline_direct_xyz.txt"
    brtp = workspace / "baseline_direct_brtp.txt"
    env = os.environ.copy()
    env["SKIP_BUILD"] = "1"
    run_cmd(
        [
            "./run_input_to_xyz.sh",
            "direct_sum",
            f"{rsphere_km}",
            str(example),
            str(xyz),
            f"{refine_factor}",
            "24",
            "72",
            "144",
            "0.2",
            "4.0",
            "0",
            "0",
            "0",
        ],
        cwd=example_root,
        env=env,
    )
    run_cmd(["./run_xyz_to_brtp.sh", str(xyz), str(brtp)], cwd=example_root, env=env)
    return load_brtp_aggregate(brtp)


def run_fortran_spectral(
    *,
    workspace: Path,
    example: Path,
    rsphere_km: float,
    refine_factor: int,
    lmax: int,
    ntheta_fit: int,
    nphi_fit: int,
    reg_lambda: float,
    reg_power: float,
) -> Tuple[Dict[str, np.ndarray], float, Dict[str, float]]:
    xyz = workspace / f"spec_l{lmax}_lam{reg_lambda:g}_pow{reg_power:g}_xyz.txt"
    brtp = workspace / f"spec_l{lmax}_lam{reg_lambda:g}_pow{reg_power:g}_brtp.txt"
    env = os.environ.copy()
    env["SKIP_BUILD"] = "1"
    env["GRAVMAG_DIAGNOSTICS"] = "1"
    wall_s, out = run_cmd(
        [
            "./run_input_to_xyz.sh",
            "spectral",
            f"{rsphere_km}",
            str(example),
            str(xyz),
            f"{refine_factor}",
            f"{lmax}",
            f"{ntheta_fit}",
            f"{nphi_fit}",
            f"{reg_lambda}",
            f"{reg_power}",
            "0",
            "0",
            "0",
        ],
        cwd=example_root,
        env=env,
    )
    run_cmd(["./run_xyz_to_brtp.sh", str(xyz), str(brtp)], cwd=example_root, env=env)
    pred = load_brtp_aggregate(brtp)
    diag = parse_gauss_diag_times(out)
    return pred, wall_s, diag


def pareto_front(rows: List[Dict[str, object]], err_key: str, time_key: str = "wall_s") -> List[Dict[str, object]]:
    pts = sorted(rows, key=lambda r: (float(r[time_key]), float(r[err_key])))
    out: List[Dict[str, object]] = []
    best_err = float("inf")
    for r in pts:
        e = float(r[err_key])
        if e < best_err:
            out.append(r)
            best_err = e
    return out


def choose_for_runtime(rows: List[Dict[str, object]], err_key: str, budget: float) -> Dict[str, object] | None:
    candidates = [r for r in rows if float(r["wall_s"]) <= budget]
    if not candidates:
        return None
    return sorted(candidates, key=lambda r: (float(r[err_key]), -int(r["lmax"]), float(r["wall_s"])))[0]


def choose_for_accuracy(rows: List[Dict[str, object]], err_key: str, target: float) -> Dict[str, object] | None:
    candidates = [r for r in rows if float(r[err_key]) <= target]
    if not candidates:
        return None
    return sorted(candidates, key=lambda r: (float(r["wall_s"]), -int(r["lmax"]), float(r[err_key])))[0]


def make_plots(
    rows: List[Dict[str, object]],
    *,
    err_key: str,
    plot_dir: Path,
    lmax_values: List[int],
    reg_lambdas: List[float],
    reg_powers: List[float],
) -> None:
    if "MPLCONFIGDIR" not in os.environ:
        mpl_dir = repo_root / "diagnostics" / ".mplconfig"
        mpl_dir.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_dir)

    import matplotlib

    if os.environ.get("DISPLAY", "") == "" and os.environ.get("MPLBACKEND", "") == "":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_dir.mkdir(parents=True, exist_ok=True)

    x = np.array([float(r["wall_s"]) for r in rows], dtype=float)
    y = np.array([float(r[err_key]) for r in rows], dtype=float)
    c = np.array([int(r["lmax"]) for r in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    sc = ax.scatter(x, y, c=c, cmap="viridis", s=28)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("lmax")
    ax.set_xlabel("Runtime [s]")
    ax.set_ylabel(err_key)
    ax.set_title("Runtime vs accuracy sweep")
    fig.savefig(plot_dir / "runtime_vs_accuracy.png", dpi=200)
    plt.close(fig)

    for rp in reg_powers:
        mat = np.full((len(reg_lambdas), len(lmax_values)), np.nan, dtype=float)
        for i, lam in enumerate(reg_lambdas):
            for j, lm in enumerate(lmax_values):
                hits = [r for r in rows if int(r["lmax"]) == lm and float(r["reg_lambda"]) == lam and float(r["reg_power"]) == rp]
                if hits:
                    mat[i, j] = float(hits[0][err_key])
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        im = ax.imshow(mat, origin="lower", aspect="auto", cmap="magma")
        ax.set_xticks(np.arange(len(lmax_values)))
        ax.set_xticklabels([str(v) for v in lmax_values], rotation=45, ha="right")
        ax.set_yticks(np.arange(len(reg_lambdas)))
        ax.set_yticklabels([f"{v:g}" for v in reg_lambdas])
        ax.set_xlabel("lmax")
        ax.set_ylabel("reg_lambda")
        ax.set_title(f"{err_key} heatmap (reg_power={rp:g})")
        cb = fig.colorbar(im, ax=ax)
        cb.set_label(err_key)
        fig.savefig(plot_dir / f"heatmap_{err_key}_pow{rp:g}.png", dpi=200)
        plt.close(fig)


def sweep(args: argparse.Namespace) -> None:
    example = Path(args.example)
    if not example.is_absolute():
        example = repo_root / example
    if not example.exists():
        raise FileNotFoundError(f"Example input not found: {example}")

    lmax_values = list(range(args.lmax_min, args.lmax_max + 1, args.lmax_step))
    reg_lambdas = parse_list_floats(args.reg_lambdas)
    reg_powers = parse_list_floats(args.reg_powers)
    runtime_budgets = parse_list_floats(args.runtime_budgets)
    accuracy_targets = parse_list_floats(args.accuracy_targets)
    err_component = args.error_component.lower()
    if err_component not in COMPONENTS:
        raise ValueError(f"error_component must be one of {COMPONENTS}, got: {err_component}")
    err_key = f"rmse_{err_component}_shtools"

    workspace = Path(args.workspace)
    if not workspace.is_absolute():
        workspace = repo_root / workspace
    workspace = workspace / example.stem
    workspace.mkdir(parents=True, exist_ok=True)

    run_cmd(["./build_gravmag_tools.sh"], cwd=example_root)
    baseline = run_baseline_direct(
        workspace=workspace,
        example=example,
        rsphere_km=args.rsphere_km,
        refine_factor=args.refine_factor,
    )

    shtools_cache: Dict[int, Dict[str, np.ndarray]] = {}
    for lm in lmax_values:
        shtools_cache[lm] = predict_shtools_components(baseline, lm)

    rows: List[Dict[str, object]] = []
    total = len(lmax_values) * len(reg_lambdas) * len(reg_powers)
    done = 0
    for lm in lmax_values:
        for lam in reg_lambdas:
            for rp in reg_powers:
                done += 1
                row: Dict[str, object] = {
                    "case": example.stem,
                    "lmax": lm,
                    "reg_lambda": lam,
                    "reg_power": rp,
                }
                try:
                    pred, wall_s, diag = run_fortran_spectral(
                        workspace=workspace,
                        example=example,
                        rsphere_km=args.rsphere_km,
                        refine_factor=args.refine_factor,
                        lmax=lm,
                        ntheta_fit=args.ntheta_fit,
                        nphi_fit=args.nphi_fit,
                        reg_lambda=lam,
                        reg_power=rp,
                    )
                    _lon, _lat, ref_d, pred_d = align_components(baseline, pred)
                    _lon2, _lat2, ref_s, pred_s = align_components(shtools_cache[lm], pred)
                    row.update(component_metrics(ref_d, pred_d, "direct_sum"))
                    row.update(component_metrics(ref_s, pred_s, "shtools"))
                    row["wall_s"] = wall_s
                    for k, v in diag.items():
                        row[f"diag_{k}"] = v
                    row["status"] = "ok"
                except Exception as exc:  # noqa: BLE001
                    row["status"] = f"error: {exc}"
                rows.append(row)
                print(f"[{done}/{total}] lmax={lm} lambda={lam:g} power={rp:g} -> {row['status']}")

    ok_rows = [r for r in rows if r.get("status") == "ok"]
    if not ok_rows:
        raise RuntimeError("No successful runs in sweep.")

    frontier = pareto_front(ok_rows, err_key=err_key)
    runtime_recs = [(b, choose_for_runtime(ok_rows, err_key, b)) for b in runtime_budgets]
    accuracy_recs = [(a, choose_for_accuracy(ok_rows, err_key, a)) for a in accuracy_targets]

    out_csv = Path(args.output_csv)
    if not out_csv.is_absolute():
        out_csv = repo_root / out_csv
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = sorted(set().union(*(r.keys() for r in rows)))
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    out_md = Path(args.output_md)
    if not out_md.is_absolute():
        out_md = repo_root / out_md
    out_md.parent.mkdir(parents=True, exist_ok=True)

    plot_dir = Path(args.plot_dir)
    if not plot_dir.is_absolute():
        plot_dir = repo_root / plot_dir
    make_plots(
        ok_rows,
        err_key=err_key,
        plot_dir=plot_dir,
        lmax_values=lmax_values,
        reg_lambdas=reg_lambdas,
        reg_powers=reg_powers,
    )

    def fmt_rec(r: Dict[str, object] | None) -> str:
        if r is None:
            return "no feasible configuration"
        return (
            f"lmax={int(r['lmax'])}, reg_lambda={float(r['reg_lambda']):g}, reg_power={float(r['reg_power']):g}, "
            f"wall_s={float(r['wall_s']):.3f}, {err_key}={float(r[err_key]):.6e}"
        )

    lines: List[str] = []
    lines.append("# SHTOOLS Operating Map")
    lines.append("")
    lines.append(f"- Example: `{example}`")
    lines.append(f"- Sweep size: `{total}` runs, successful `{len(ok_rows)}`")
    lines.append(f"- Target error metric: `{err_key}`")
    lines.append("")
    lines.append("## Runtime Budget Recommendations")
    lines.append("")
    lines.append("| Budget (s) | Recommended (lmax, reg_lambda, reg_power) |")
    lines.append("|---:|---|")
    for budget, rec in runtime_recs:
        lines.append(f"| {budget:.3f} | {fmt_rec(rec)} |")
    lines.append("")
    lines.append("## Accuracy Target Recommendations")
    lines.append("")
    lines.append("| Target RMSE | Recommended (lmax, reg_lambda, reg_power) |")
    lines.append("|---:|---|")
    for target, rec in accuracy_recs:
        lines.append(f"| {target:.6e} | {fmt_rec(rec)} |")
    lines.append("")
    lines.append("## Pareto Frontier (runtime vs target error)")
    lines.append("")
    lines.append("| lmax | reg_lambda | reg_power | wall_s | target_error |")
    lines.append("|---:|---:|---:|---:|---:|")
    for r in frontier:
        lines.append(
            f"| {int(r['lmax'])} | {float(r['reg_lambda']):g} | {float(r['reg_power']):g} | "
            f"{float(r['wall_s']):.6f} | {float(r[err_key]):.6e} |"
        )
    lines.append("")
    lines.append("## Best Configurations By Target Error")
    lines.append("")
    best_sorted = sorted(ok_rows, key=lambda r: float(r[err_key]))[:20]
    lines.append("| rank | lmax | reg_lambda | reg_power | wall_s | target_error |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for i, r in enumerate(best_sorted, 1):
        lines.append(
            f"| {i} | {int(r['lmax'])} | {float(r['reg_lambda']):g} | {float(r['reg_power']):g} | "
            f"{float(r['wall_s']):.6f} | {float(r[err_key]):.6e} |"
        )
    lines.append("")
    lines.append(f"Plots written to `{plot_dir}`.")

    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    out_json = out_md.with_suffix(".json")
    summary = {
        "example": str(example),
        "err_key": err_key,
        "runtime_recommendations": [
            {"budget_s": b, "recommendation": rec} for b, rec in runtime_recs
        ],
        "accuracy_recommendations": [
            {"target_rmse": a, "recommendation": rec} for a, rec in accuracy_recs
        ],
        "pareto_frontier": frontier,
        "output_csv": str(out_csv),
        "output_md": str(out_md),
        "plot_dir": str(plot_dir),
    }
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Wrote sweep CSV: {out_csv}")
    print(f"Wrote operating map markdown: {out_md}")
    print(f"Wrote operating map json: {out_json}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Build operating map vs SHTOOLS target for GravMag spectral solver.")
    p.add_argument(
        "--example",
        default="examples/lunar_examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in",
        help="Example .in file to sweep.",
    )
    p.add_argument("--rsphere-km", type=float, default=1737.4)
    p.add_argument("--refine-factor", type=int, default=2)
    p.add_argument("--ntheta-fit", type=int, default=72)
    p.add_argument("--nphi-fit", type=int, default=144)
    p.add_argument("--lmax-min", type=int, default=12)
    p.add_argument("--lmax-max", type=int, default=48)
    p.add_argument("--lmax-step", type=int, default=6)
    p.add_argument("--reg-lambdas", default="0,0.02,0.05,0.1,0.2,0.5,1,2")
    p.add_argument("--reg-powers", default="2,4,6,8")
    p.add_argument("--error-component", default="btot", help="One of: br,btheta,bphi,btot")
    p.add_argument("--runtime-budgets", default="0.5,1,2,4,8,16,32")
    p.add_argument("--accuracy-targets", default="1,2,5,10,20,50,100")
    p.add_argument("--workspace", default="diagnostics/operating_map_workspace")
    p.add_argument("--output-csv", default="diagnostics/shtools_operating_map_grid.csv")
    p.add_argument("--output-md", default="diagnostics/shtools_operating_map.md")
    p.add_argument("--plot-dir", default="diagnostics/shtools_operating_map_plots")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    sweep(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
