#!/usr/bin/env python3
"""
Fortran diagnostics for GravMag Sphere:
- compiler runtime benchmarks with stage-level timings/memory estimates
- spherical harmonic lmax sweep with edge-ringing diagnostics
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - explicit runtime dependency message
    raise SystemExit("This script requires numpy. Install it in your Python environment.") from exc


OPENMP_FLAGS: Dict[str, List[str]] = {
    "gfortran": ["-fopenmp"],
    "ifx": ["-qopenmp"],
    "ifort": ["-qopenmp"],
    "flang-new": ["-fopenmp"],
    "nvfortran": ["-mp"],
}

SUPPORTED_COMPILERS = tuple(OPENMP_FLAGS.keys())
DIAG_PATTERN = re.compile(r"^DIAG\|([^|]+)\|([^|]+)\|([^=]+)=(.+)$")


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    solver: str


def solver_diag_key(solver: str) -> str:
    return solver


def run_cmd(
    cmd: Sequence[str],
    cwd: Path,
    env: Dict[str, str] | None = None,
    check: bool = True,
) -> Tuple[float, str]:
    start = time.perf_counter()
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    elapsed = time.perf_counter() - start
    if check and proc.returncode != 0:
        command = shlex.join(cmd)
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {command}\n"
            f"--- begin output ---\n{proc.stdout}\n--- end output ---"
        )
    return elapsed, proc.stdout


def detect_compilers(requested: Sequence[str]) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    for compiler in requested:
        path = shutil.which(compiler)
        if path:
            out.append((compiler, path))
    return out


def parse_diag_lines(output: str) -> Dict[str, Dict[str, Dict[str, float]]]:
    parsed: Dict[str, Dict[str, Dict[str, float]]] = {}
    for line in output.splitlines():
        m = DIAG_PATTERN.match(line.strip())
        if not m:
            continue
        solver, category, key, value_raw = m.groups()
        parsed.setdefault(solver, {}).setdefault(category, {})
        value_raw = value_raw.strip()
        try:
            if re.fullmatch(r"[-+]?\d+", value_raw):
                value = float(int(value_raw))
            else:
                value = float(value_raw)
        except ValueError:
            continue
        parsed[solver][category][key] = value
    return parsed


def build_binaries(
    root: Path,
    compiler: str,
    compiler_path: str,
    opt_flags: Sequence[str],
    extra_flags: Sequence[str],
) -> Dict[str, object]:
    build_root = root / "build" / "diagnostics" / compiler
    mod_dir = build_root / "mod"
    bin_dir = build_root / "bin"
    mod_dir.mkdir(parents=True, exist_ok=True)
    bin_dir.mkdir(parents=True, exist_ok=True)

    source_dir = root / "fortran"
    state = source_dir / "gravmag_sphere_state.f90"
    subs = source_dir / "gravmag_sphere_subs.f90"
    physics = source_dir / "gravmag_sphere_physics.f90"
    direct_sum = source_dir / "gravmag_sphere_direct_sum.f90"
    spectral = source_dir / "gravmag_sphere_spectral.f90"
    conv = source_dir / "gravmag_xyz_to_brtp.f90"
    paths = source_dir / "gravmag_paths.f90"

    omp_flags = OPENMP_FLAGS.get(compiler, [])
    common = [compiler_path, "-std=f2008", *opt_flags, *extra_flags]

    compile_times: Dict[str, float] = {}

    cmd_direct_sum = [
        *common,
        *omp_flags,
        f"-J{mod_dir}",
        f"-I{mod_dir}",
        str(paths),
        str(state),
        str(subs),
        str(physics),
        str(direct_sum),
        "-o",
        str(bin_dir / "gravmag_sphere_direct_sum"),
    ]
    compile_times["gravmag_sphere_direct_sum"], _ = run_cmd(cmd_direct_sum, cwd=root)

    cmd_spectral = [
        *common,
        f"-J{mod_dir}", f"-I{mod_dir}", str(paths),
        str(spectral),
        "-o",
        str(bin_dir / "gravmag_sphere_spectral"),
    ]
    compile_times["gravmag_sphere_spectral"], _ = run_cmd(cmd_spectral, cwd=root)

    cmd_conv = [
        *common,
        f"-J{mod_dir}", f"-I{mod_dir}", str(paths),
        str(conv),
        "-o",
        str(bin_dir / "gravmag_xyz_to_brtp"),
    ]
    compile_times["gravmag_xyz_to_brtp"], _ = run_cmd(cmd_conv, cwd=root)

    return {
        "build_root": str(build_root),
        "binaries": {
            "direct_sum": str(bin_dir / "gravmag_sphere_direct_sum"),
            "spectral": str(bin_dir / "gravmag_sphere_spectral"),
            "convert": str(bin_dir / "gravmag_xyz_to_brtp"),
        },
        "compile_times_s": compile_times,
    }


def noncomment_lines(path: Path) -> List[str]:
    lines: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s or s.startswith("#") or s.startswith("!"):
                continue
            lines.append(s)
    return lines


def parse_first_body_geometry(path: Path) -> Dict[str, object]:
    lines = noncomment_lines(path)
    if len(lines) < 7:
        raise ValueError(f"Could not parse cards in {path}")

    c2 = lines[1].split()
    c3 = lines[2].split()
    nblim = int(float(c3[3]))
    dlat = abs(float(c2[2]))
    dlon = abs(float(c2[3]))

    if nblim == 1:
        card7 = lines[6].split()
        lat_max, lat_min = float(card7[0]), float(card7[1])
        lon_max, lon_min = float(card7[2]), float(card7[3])
        poly_lat = np.array([lat_min, lat_max, lat_max, lat_min, lat_min], dtype=float)
        poly_lon = np.array([lon_min, lon_min, lon_max, lon_max, lon_min], dtype=float)
    else:
        head = lines[6].split()
        npts = int(float(head[0]))
        poly_lat = np.zeros(npts, dtype=float)
        poly_lon = np.zeros(npts, dtype=float)
        for i in range(npts):
            p1, p2 = (float(v) for v in lines[7 + i].split()[:2])
            if abs(p1) <= 90.0 and abs(p2) > 90.0:
                poly_lat[i], poly_lon[i] = p1, p2
            elif abs(p1) > 90.0 and abs(p2) <= 90.0:
                poly_lon[i], poly_lat[i] = p1, p2
            else:
                poly_lat[i], poly_lon[i] = p1, p2

    return {"poly_lat": poly_lat, "poly_lon": poly_lon, "dlat": dlat, "dlon": dlon}


def wrap180(values: np.ndarray) -> np.ndarray:
    return np.mod(values + 180.0, 360.0) - 180.0


def unwrap_near(values: np.ndarray, ref: float) -> np.ndarray:
    out = wrap180(values)
    out = np.where(out - ref >= 180.0, out - 360.0, out)
    out = np.where(out - ref < -180.0, out + 360.0, out)
    return out


def point_edge_distance_deg(
    lon: np.ndarray,
    lat: np.ndarray,
    poly_lon: np.ndarray,
    poly_lat: np.ndarray,
) -> np.ndarray:
    ref = float(wrap180(np.array([np.mean(poly_lon)]))[0])
    lon_u = unwrap_near(lon, ref)
    poly_lon_u = unwrap_near(poly_lon, ref)

    lat0_rad = math.radians(float(np.mean(poly_lat)))
    scale = math.cos(lat0_rad)

    x = lon_u * scale
    y = lat
    px = poly_lon_u * scale
    py = poly_lat
    if px[0] != px[-1] or py[0] != py[-1]:
        px = np.append(px, px[0])
        py = np.append(py, py[0])

    dist = np.full(x.shape[0], np.inf, dtype=float)
    for i in range(px.size - 1):
        x1, y1 = px[i], py[i]
        x2, y2 = px[i + 1], py[i + 1]
        dx = x2 - x1
        dy = y2 - y1
        denom = dx * dx + dy * dy
        if denom <= 0.0:
            d = np.hypot(x - x1, y - y1)
        else:
            t = ((x - x1) * dx + (y - y1) * dy) / denom
            t = np.clip(t, 0.0, 1.0)
            qx = x1 + t * dx
            qy = y1 + t * dy
            d = np.hypot(x - qx, y - qy)
        dist = np.minimum(dist, d)
    return dist


def load_brtp_total(path: Path) -> Dict[Tuple[float, float], float]:
    arr = np.loadtxt(path, comments="#")
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] < 4:
        raise ValueError(f"Unexpected BRTP format in {path}")

    lon = np.round(arr[:, 1], 6)
    lat = np.round(arr[:, 2], 6)
    br = arr[:, 3]
    out: Dict[Tuple[float, float], float] = {}
    for lo, la, brv in zip(lon, lat, br):
        key = (float(lo), float(la))
        out[key] = out.get(key, 0.0) + float(brv)
    return out


def align_fields(
    baseline: Dict[Tuple[float, float], float],
    candidate: Dict[Tuple[float, float], float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    keys = sorted(set(baseline.keys()) & set(candidate.keys()))
    if not keys:
        raise ValueError("No overlapping lon/lat points between baseline and candidate fields.")
    lon = np.array([k[0] for k in keys], dtype=float)
    lat = np.array([k[1] for k in keys], dtype=float)
    residual = np.array([candidate[k] - baseline[k] for k in keys], dtype=float)
    baseline_vals = np.array([baseline[k] for k in keys], dtype=float)
    return lon, lat, residual, baseline_vals


def rms(values: np.ndarray) -> float:
    if values.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(values * values)))


def benchmark_subcommand(args: argparse.Namespace) -> Dict[str, object]:
    root = Path(__file__).resolve().parent.parent
    diagnostics_dir = root / "diagnostics"
    diagnostics_dir.mkdir(parents=True, exist_ok=True)

    requested = [s.strip() for s in args.compilers.split(",") if s.strip()]
    compilers = detect_compilers(requested)
    if not compilers:
        raise RuntimeError(f"No requested compilers found on PATH: {requested}")

    opt_flags = shlex.split(args.opt_flags)
    extra_flags = shlex.split(args.extra_flags) if args.extra_flags else []

    benchmark_cases = [
        BenchmarkCase(name="direct_polygon", solver="direct_sum"),
        BenchmarkCase(name="spectral_polygon", solver="spectral"),
    ]

    output_root = root / "diagnostics" / "runs"
    output_root.mkdir(parents=True, exist_ok=True)
    example = Path(args.example)
    if not example.is_absolute():
        example = root / example
    if not example.exists():
        raise FileNotFoundError(f"Example input not found: {example}")

    rsphere = str(args.rsphere_km)

    report: Dict[str, object] = {
        "timestamp": int(time.time()),
        "root": str(root),
        "example": str(example),
        "compilers": [],
    }

    for compiler, compiler_path in compilers:
        build_info = build_binaries(root, compiler, compiler_path, opt_flags, extra_flags)
        binaries = build_info["binaries"]  # type: ignore[assignment]
        compiler_entry: Dict[str, object] = {
            "compiler": compiler,
            "compiler_path": compiler_path,
            "build": build_info,
            "cases": [],
        }
        env = os.environ.copy()
        env["GRAVMAG_DIAGNOSTICS"] = "1"
        env["OMP_NUM_THREADS"] = str(args.omp_threads)

        for case in benchmark_cases:
            run_records: List[Dict[str, object]] = []
            for rep in range(args.repeats):
                run_dir = output_root / compiler / case.name / f"rep{rep + 1}"
                run_dir.mkdir(parents=True, exist_ok=True)
                xyz_out = run_dir / f"{case.name}.xyz.txt"

                if case.solver == "direct_sum":
                    cmd = [
                        binaries["direct_sum"],  # type: ignore[index]
                        rsphere,
                        str(example),
                        str(xyz_out),
                        str(args.refine_factor),
                        "0",
                        "0",
                        "0",
                    ]
                else:
                    cmd = [
                        binaries["spectral"],  # type: ignore[index]
                        rsphere,
                        str(example),
                        str(xyz_out),
                        str(args.lmax),
                        str(args.refine_factor),
                        str(args.ntheta_fit),
                        str(args.nphi_fit),
                        str(args.reg_lambda),
                        str(args.reg_power),
                        "0",
                        "0",
                        "0",
                    ]

                elapsed, output = run_cmd(cmd, cwd=root, env=env, check=True)
                diag = parse_diag_lines(output)
                run_records.append(
                    {
                        "rep": rep + 1,
                        "elapsed_s": elapsed,
                        "diag": diag,
                    }
                )

            wall_times = [float(r["elapsed_s"]) for r in run_records]
            summary: Dict[str, object] = {
                "name": case.name,
                "solver": case.solver,
                "runs": run_records,
                "elapsed_s": {
                    "min": min(wall_times),
                    "median": statistics.median(wall_times),
                    "max": max(wall_times),
                },
            }

            # Aggregate stage timing + memory summaries from DIAG lines
            stage_totals: Dict[str, List[float]] = {}
            mem_totals: Dict[str, List[float]] = {}
            for run in run_records:
                diag = run["diag"]  # type: ignore[assignment]
                solver_diag = (
                    diag.get(solver_diag_key(case.solver), {}) if isinstance(diag, dict) else {}
                )
                for key, value in solver_diag.get("time", {}).items():
                    stage_totals.setdefault(key, []).append(float(value))
                for key, value in solver_diag.get("memory", {}).items():
                    mem_totals.setdefault(key, []).append(float(value))

            if stage_totals:
                stage_means = {k: statistics.mean(v) for k, v in stage_totals.items() if v}
                summary["time_stage_mean_s"] = stage_means
                phase_only = {k: v for k, v in stage_means.items() if k != "body_total_s"}
                if phase_only:
                    slowest = max(phase_only.items(), key=lambda kv: kv[1])
                    summary["slowest_stage"] = {"name": slowest[0], "mean_s": slowest[1]}
            if mem_totals:
                mem_means = {k: statistics.mean(v) for k, v in mem_totals.items() if v}
                summary["memory_stage_mean_mib"] = mem_means
                peak = max(mem_means.items(), key=lambda kv: kv[1])
                summary["largest_memory_bucket"] = {"name": peak[0], "mean_mib": peak[1]}

            compiler_entry["cases"].append(summary)  # type: ignore[arg-type]

        report["compilers"].append(compiler_entry)  # type: ignore[arg-type]

    output_json = Path(args.output_json)
    if not output_json.is_absolute():
        output_json = root / output_json
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")

    output_md = Path(args.output_md)
    if not output_md.is_absolute():
        output_md = root / output_md
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text(render_benchmark_markdown(report), encoding="utf-8")

    print(f"Wrote benchmark report JSON: {output_json}")
    print(f"Wrote benchmark report markdown: {output_md}")
    return report


def render_benchmark_markdown(report: Dict[str, object]) -> str:
    lines: List[str] = []
    lines.append("# Fortran Compiler Benchmark")
    lines.append("")
    lines.append(f"- Example: `{report['example']}`")
    lines.append("")
    lines.append("| Compiler | Case | Median wall (s) | Slowest stage | Largest memory bucket |")
    lines.append("|---|---:|---:|---|---|")

    for compiler in report.get("compilers", []):
        comp_name = compiler.get("compiler", "unknown")
        for case in compiler.get("cases", []):
            elapsed = case.get("elapsed_s", {}).get("median", float("nan"))
            slowest = case.get("slowest_stage", {})
            slowest_txt = (
                f"{slowest.get('name', 'n/a')} ({slowest.get('mean_s', float('nan')):.3f}s)"
                if slowest
                else "n/a"
            )
            mem = case.get("largest_memory_bucket", {})
            mem_txt = (
                f"{mem.get('name', 'n/a')} ({mem.get('mean_mib', float('nan')):.2f} MiB)"
                if mem
                else "n/a"
            )
            lines.append(
                f"| {comp_name} | {case.get('name', 'n/a')} | {float(elapsed):.3f} | "
                f"{slowest_txt} | {mem_txt} |"
            )
    lines.append("")
    lines.append(
        "Notes: stage and memory values come from solver DIAG instrumentation (`GRAVMAG_DIAGNOSTICS=1`)."
    )
    return "\n".join(lines) + "\n"


def lmax_sweep_subcommand(args: argparse.Namespace) -> Dict[str, object]:
    root = Path(__file__).resolve().parent.parent
    diagnostics_dir = root / "diagnostics"
    diagnostics_dir.mkdir(parents=True, exist_ok=True)

    compiler = args.compiler
    compiler_path = shutil.which(compiler)
    if not compiler_path:
        raise RuntimeError(f"Compiler not found on PATH: {compiler}")

    opt_flags = shlex.split(args.opt_flags)
    extra_flags = shlex.split(args.extra_flags) if args.extra_flags else []
    build_info = build_binaries(root, compiler, compiler_path, opt_flags, extra_flags)
    binaries = build_info["binaries"]  # type: ignore[assignment]

    example = Path(args.example)
    if not example.is_absolute():
        example = root / example
    if not example.exists():
        raise FileNotFoundError(f"Example input not found: {example}")

    geom = parse_first_body_geometry(example)
    poly_lon = geom["poly_lon"]  # type: ignore[assignment]
    poly_lat = geom["poly_lat"]  # type: ignore[assignment]
    spacing_deg = max(float(geom["dlat"]), float(geom["dlon"]))
    edge_band = args.edge_band_mult * spacing_deg
    far_band = args.far_band_mult * spacing_deg

    run_root = root / "diagnostics" / "lmax_sweep_runs" / compiler
    run_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["GRAVMAG_DIAGNOSTICS"] = "1"
    env["OMP_NUM_THREADS"] = str(args.omp_threads)

    direct_xyz = run_root / "baseline_direct.xyz.txt"
    direct_brtp = run_root / "baseline_direct.brtp.txt"

    rsphere = str(args.rsphere_km)
    direct_cmd = [
        binaries["direct_sum"],  # type: ignore[index]
        rsphere,
        str(example),
        str(direct_xyz),
        str(args.refine_factor),
        "0",
        "0",
        "0",
    ]
    run_cmd(direct_cmd, cwd=root, env=env, check=True)
    run_cmd([binaries["convert"], str(direct_xyz), str(direct_brtp)], cwd=root, env=env, check=True)  # type: ignore[index]
    baseline = load_brtp_total(direct_brtp)

    lvals = list(range(args.lmin, args.lmax + 1, args.lstep))
    rows: List[Dict[str, object]] = []
    pass_lmax: List[int] = []

    signal_scale = np.percentile(np.abs(np.array(list(baseline.values()), dtype=float)), 95)
    signal_scale = float(max(signal_scale, 1.0e-12))

    for lmax in lvals:
        xyz_out = run_root / f"spectral_l{lmax}.xyz.txt"
        brtp_out = run_root / f"spectral_l{lmax}.brtp.txt"

        spectral_cmd = [
            binaries["spectral"],  # type: ignore[index]
            rsphere,
            str(example),
            str(xyz_out),
            str(lmax),
            str(args.refine_factor),
            str(args.ntheta_fit),
            str(args.nphi_fit),
            str(args.reg_lambda),
            str(args.reg_power),
            "0",
            "0",
            "0",
        ]
        elapsed, output = run_cmd(spectral_cmd, cwd=root, env=env, check=True)
        diag = parse_diag_lines(output)
        run_cmd([binaries["convert"], str(xyz_out), str(brtp_out)], cwd=root, env=env, check=True)  # type: ignore[index]
        candidate = load_brtp_total(brtp_out)

        lon, lat, residual, baseline_vals = align_fields(baseline, candidate)
        edge_dist = point_edge_distance_deg(lon, lat, poly_lon, poly_lat)
        edge_mask = edge_dist <= edge_band
        far_mask = edge_dist >= far_band
        if np.count_nonzero(far_mask) < 8:
            far_mask = ~edge_mask

        edge_rms = rms(residual[edge_mask])
        far_rms = rms(residual[far_mask])
        global_rmse = rms(residual)
        edge_amp = edge_rms / max(far_rms, 1.0e-12)
        edge_overshoot = float(np.max(np.abs(residual[edge_mask])) / signal_scale) if np.any(edge_mask) else float("inf")

        ringing_ok = bool(
            edge_amp <= args.edge_amp_threshold and edge_overshoot <= args.overshoot_threshold
        )
        if ringing_ok:
            pass_lmax.append(lmax)

        row = {
            "lmax": lmax,
            "wall_s": elapsed,
            "global_rmse": global_rmse,
            "edge_rms": edge_rms,
            "far_rms": far_rms,
            "edge_amp": edge_amp,
            "edge_overshoot": edge_overshoot,
            "ringing_ok": ringing_ok,
            "diag": diag.get(solver_diag_key("spectral"), {}),
        }
        rows.append(row)

    recommended = max(pass_lmax) if pass_lmax else None

    output_csv = Path(args.output_csv)
    if not output_csv.is_absolute():
        output_csv = root / output_csv
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "lmax",
                "wall_s",
                "global_rmse",
                "edge_rms",
                "far_rms",
                "edge_amp",
                "edge_overshoot",
                "ringing_ok",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row["lmax"],
                    f"{row['wall_s']:.6f}",
                    f"{row['global_rmse']:.6e}",
                    f"{row['edge_rms']:.6e}",
                    f"{row['far_rms']:.6e}",
                    f"{row['edge_amp']:.6f}",
                    f"{row['edge_overshoot']:.6f}",
                    int(bool(row["ringing_ok"])),
                ]
            )

    output_md = Path(args.output_md)
    if not output_md.is_absolute():
        output_md = root / output_md
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text(
        render_lmax_markdown(
            compiler=compiler,
            example=example,
            rows=rows,
            recommended=recommended,
            edge_band=edge_band,
            far_band=far_band,
            edge_amp_threshold=args.edge_amp_threshold,
            overshoot_threshold=args.overshoot_threshold,
        ),
        encoding="utf-8",
    )

    print(f"Wrote lmax sweep CSV: {output_csv}")
    print(f"Wrote lmax sweep markdown: {output_md}")
    if recommended is None:
        print("No lmax satisfied current ringing thresholds.")
    else:
        print(f"Recommended maximum lmax without edge-ringing (by current thresholds): {recommended}")

    return {
        "compiler": compiler,
        "example": str(example),
        "recommended_lmax": recommended,
        "rows": rows,
    }


def render_lmax_markdown(
    *,
    compiler: str,
    example: Path,
    rows: Sequence[Dict[str, object]],
    recommended: int | None,
    edge_band: float,
    far_band: float,
    edge_amp_threshold: float,
    overshoot_threshold: float,
) -> str:
    lines: List[str] = []
    lines.append("# Spherical Harmonic lmax Sweep")
    lines.append("")
    lines.append(f"- Compiler: `{compiler}`")
    lines.append(f"- Example: `{example}`")
    lines.append(f"- Edge band: `{edge_band:.3f}` deg")
    lines.append(f"- Far band: `{far_band:.3f}` deg")
    lines.append(f"- Ringing thresholds: edge_amp <= `{edge_amp_threshold}`, edge_overshoot <= `{overshoot_threshold}`")
    lines.append("")
    lines.append("| lmax | Wall (s) | RMSE | Edge amp | Edge overshoot | Ringing OK |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        lines.append(
            f"| {int(row['lmax'])} | {float(row['wall_s']):.3f} | {float(row['global_rmse']):.3e} | "
            f"{float(row['edge_amp']):.3f} | {float(row['edge_overshoot']):.3f} | "
            f"{'yes' if row['ringing_ok'] else 'no'} |"
        )
    lines.append("")
    if recommended is None:
        lines.append("Recommended max lmax: none (thresholds too strict for this setup).")
    else:
        lines.append(f"Recommended max lmax: **{recommended}**")
    lines.append(
        "Caution: this is empirical for the selected case, solver regularization, mesh, and thresholds."
    )
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Fortran diagnostics for GravMag Sphere.")
    sub = parser.add_subparsers(dest="command", required=True)

    bench = sub.add_parser("benchmark", help="Benchmark available compilers and report hotspots.")
    bench.add_argument("--compilers", default="gfortran,ifx,ifort,flang-new,nvfortran")
    bench.add_argument("--repeats", type=int, default=3)
    bench.add_argument(
        "--example",
        default="examples/lunar_examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in",
    )
    bench.add_argument("--rsphere-km", type=float, default=1737.4)
    bench.add_argument("--refine-factor", type=int, default=2)
    bench.add_argument("--lmax", type=int, default=24)
    bench.add_argument("--ntheta-fit", type=int, default=72)
    bench.add_argument("--nphi-fit", type=int, default=144)
    bench.add_argument("--reg-lambda", type=float, default=0.2)
    bench.add_argument("--reg-power", type=float, default=4.0)
    bench.add_argument("--omp-threads", type=int, default=1)
    bench.add_argument("--opt-flags", default="-O2")
    bench.add_argument("--extra-flags", default="")
    bench.add_argument("--output-json", default="diagnostics/fortran_benchmark_report.json")
    bench.add_argument("--output-md", default="diagnostics/fortran_benchmark_report.md")

    sweep = sub.add_parser("lmax-sweep", help="Sweep spectral lmax and score edge ringing.")
    sweep.add_argument("--compiler", default="gfortran")
    sweep.add_argument(
        "--example",
        default="examples/lunar_examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in",
    )
    sweep.add_argument("--rsphere-km", type=float, default=1737.4)
    sweep.add_argument("--refine-factor", type=int, default=2)
    sweep.add_argument("--lmin", type=int, default=8)
    sweep.add_argument("--lmax", type=int, default=96)
    sweep.add_argument("--lstep", type=int, default=8)
    sweep.add_argument("--ntheta-fit", type=int, default=72)
    sweep.add_argument("--nphi-fit", type=int, default=144)
    sweep.add_argument("--reg-lambda", type=float, default=0.2)
    sweep.add_argument("--reg-power", type=float, default=4.0)
    sweep.add_argument("--edge-band-mult", type=float, default=1.5)
    sweep.add_argument("--far-band-mult", type=float, default=4.0)
    sweep.add_argument("--edge-amp-threshold", type=float, default=2.0)
    sweep.add_argument("--overshoot-threshold", type=float, default=0.15)
    sweep.add_argument("--omp-threads", type=int, default=1)
    sweep.add_argument("--opt-flags", default="-O2")
    sweep.add_argument("--extra-flags", default="")
    sweep.add_argument("--output-csv", default="diagnostics/lmax_sweep.csv")
    sweep.add_argument("--output-md", default="diagnostics/lmax_sweep.md")

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "benchmark":
        benchmark_subcommand(args)
    elif args.command == "lmax-sweep":
        lmax_sweep_subcommand(args)
    else:  # pragma: no cover
        parser.error(f"Unknown command: {args.command}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
