#!/usr/bin/env python3
"""
Run every example input through the Fortran pipeline and generate 2x2 plots :)

This script is the Python-side "batch driver":
1) compile required Fortran tools,
2) run either spectral or point_source solver to XYZ output,
3) convert XYZ -> Br/Btheta/Bphi/Btot,
4) plot summed multi-body response for each case!
"""

import argparse
import re
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def run_cmd(cmd, cwd):
    """Run a command and raise with stdout/stderr if it fails."""
    p = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )
    return p.stdout


def load_brtp(path: Path):
    """Load BRTP table with columns: body_id, lon, lat, Br, Btheta, Bphi, Btot."""
    arr = np.loadtxt(path, comments="#")
    if arr.ndim != 2 or arr.shape[1] < 7:
        raise ValueError(f"Unexpected BRTP output format in {path}")
    out = {
        "body_id": arr[:, 0].astype(int),
        "lon": arr[:, 1],
        "lat": arr[:, 2],
        "Br": arr[:, 3],
        "Btheta": arr[:, 4],
        "Bphi": arr[:, 5],
        "Btot": arr[:, 6],
    }
    return out


def aggregate_total(data):
    """Sum vector components over body_id at each unique lon/lat location."""
    lon_r = np.round(data["lon"], 6)
    lat_r = np.round(data["lat"], 6)
    pts = np.column_stack((lon_r, lat_r))
    uniq, inv = np.unique(pts, axis=0, return_inverse=True)

    br = np.zeros(uniq.shape[0], dtype=float)
    bt = np.zeros(uniq.shape[0], dtype=float)
    bp = np.zeros(uniq.shape[0], dtype=float)
    np.add.at(br, inv, data["Br"])
    np.add.at(bt, inv, data["Btheta"])
    np.add.at(bp, inv, data["Bphi"])
    btot = np.sqrt(br * br + bt * bt + bp * bp)

    return {
        "lon": uniq[:, 0],
        "lat": uniq[:, 1],
        "Br": br,
        "Btheta": bt,
        "Bphi": bp,
        "Btot": btot,
    }


def to_grid(lon, lat, values):
    """Convert scattered lon/lat/value rows into a dense regular grid."""
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
        missing = int(np.isnan(z).sum())
        raise RuntimeError(f"Grid has {missing} missing cells after gridding.")
    return ulon, ulat, z


def parse_unit_from_header(path: Path):
    """Infer display units from header tags in the BRTP text file."""
    unit = ""
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if "_nT" in s:
                    unit = "nT"
                elif "_mGal" in s or "_MGAL" in s:
                    unit = "mGal"
                continue
            break
    return unit


def parse_unit_from_fname(path: Path):
    """Infer display units from file name"""
    if "_mag_" in str(path):
        unit = "nT"
        return unit

    elif "_grav_" in str(path):
        unit = "mGal"
        return unit
    else:
        raise ValueError("Filename must contain either mag or grav file!")


def plot_brtp(data, out_png: Path, title: str, unit: str):
    """Create and save Br/Btheta/Bphi/Btot 2x2 panel figure."""
    keys = ["Br", "Btheta", "Bphi", "Btot"]
    labels = ["Br", "Btheta", "Bphi", "Btot"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    axes = axes.ravel()

    for ax, key, label in zip(axes, keys, labels):
        ulon, ulat, z = to_grid(data["lon"], data["lat"], data[key])
        if key == "Btot":
            vmin, vmax = 0.0, float(np.nanmax(z))
            cmap = "viridis"
        else:
            m = float(np.nanmax(np.abs(z)))
            vmin, vmax = -m, m
            cmap = "RdBu_r"

        im = ax.pcolormesh(ulon, ulat, z, shading="nearest", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(label)
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        ax.set_aspect("equal", adjustable="box")
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.12, shrink=0.85)
        cb.set_label(unit if unit else "units")

    fig.suptitle(title)
    fig.savefig(out_png, dpi=250)
    plt.close(fig)


def should_skip_example(path: Path):
    """Skip non-runnable documentation example files by name"""
    name = path.name.lower()
    # Skip purely commented documentation lunar_examples
    if "comment" in name:
        return True
    return False


def pretty_case_name(stem: str) -> str:
    """Convert machine-style case names into a nice plot title"""
    s = stem
    if s.startswith("gravmag_sphere_"):
        s = s[len("gravmag_sphere_") :]

    s = re.sub(r"(\d+)body", r"\1 body", s)
    s = s.replace("fixedlim", "fixed limits")
    s = s.replace("complexlarge", "large complex polygon")
    s = s.replace("incmix", "mixed inc")
    s = s.replace("decmix", "mixed dec")
    s = re.sub(r"inc(-?\d+)", r"inc \1°", s)
    s = re.sub(r"dec(-?\d+)", r"dec \1°", s)

    token_map = {
        "mag": "magnetic",
        "grav": "gravity",
        "orient": "orientation",
        "base": "baseline",
        "weak": "weak",
        "polygon": "polygon",
        "fixed": "fixed",
        "incna": "inc n/a",
        "decna": "dec n/a",
    }

    toks = []
    for tok in s.replace("_", " ").split():
        toks.append(token_map.get(tok, tok))

    out = " ".join(toks)
    out = re.sub(r"\s+", " ", out).strip()
    if out:
        out = out[0].upper() + out[1:]
    return out


def first_body_ifield(path: Path):
    """Read Card 5 ifield from the first body block in an input file"""
    lines = []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#") or s.startswith("!"):
                continue
            lines.append(s)
            if len(lines) >= 5:
                break
    if len(lines) < 5:
        raise RuntimeError(f"Could not parse first 5 cards in {path}")
    toks = lines[4].split()
    if not toks:
        raise RuntimeError(f"Could not parse Card 5 in {path}")
    return int(float(toks[0]))


def main():
    """CLI entrypoint for end-to-end solver + conversion + plotting."""
    ap = argparse.ArgumentParser(
        description="Run all example inputs, convert XYZ->Br/Btheta/Bphi, and plot Br/Btheta/Bphi/Btot."
    )
    ap.add_argument("--solver", choices=["spectral", "point_source", "gauss_legendre", "direct", "gauss"], default="spectral")
    ap.add_argument("--rsphere-km", type=float, default=1737.4)
    ap.add_argument("--refine-factor", type=int, default=2)
    ap.add_argument("--lmax", type=int, default=24)
    ap.add_argument("--ntheta-fit", type=int, default=72)
    ap.add_argument("--nphi-fit", type=int, default=144)
    ap.add_argument("--reg-lambda", type=float, default=0.2)
    ap.add_argument("--reg-power", type=float, default=4.0)
    ap.add_argument("--source-nlat", type=int, default=0)
    ap.add_argument("--source-nlon", type=int, default=0)
    ap.add_argument("--source-nr", type=int, default=0)
    ap.add_argument("--examples-dir", "--lunar_examples-dir", default="lunar_examples")
    ap.add_argument("--output-dir", default=None)
    ap.add_argument("--figs-dir", default=None)
    args = ap.parse_args()
    args.solver = {"direct": "point_source", "gauss": "spectral"}.get(args.solver, args.solver)

    cwd = Path(__file__).resolve().parent
    examples_dir = cwd / args.examples_dir
    output_dir = cwd / args.output_dir if args.output_dir is not None else examples_dir / "output"
    figs_dir = cwd / args.figs_dir if args.figs_dir is not None else examples_dir / "figs"
    output_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    # compile from the canonical fortran source directory
    run_cmd(["bash", str(cwd / "build_gravmag_tools.sh")], cwd)

    cases = sorted(examples_dir.glob("*.in"))
    cases = [c for c in cases if not should_skip_example(c)]
    if not cases:
        raise RuntimeError("No runnable example .in files found.")

    print(f"Found {len(cases)} runnable example cases.")

    for ex in cases:
        stem = ex.stem
        out_xyz = output_dir / f"{stem}_xyz.txt"
        out_brtp = output_dir / f"{stem}_brtp.txt"
        out_png = figs_dir / f"{stem}_brtp_2x2.png"

        # Keep this parse for early input-format validation (Card 5 exists)
        _ = first_body_ifield(ex)
        use_spectral = args.solver == "spectral"
        solver_label = args.solver

        print(f"[RUN] {ex.name} ({solver_label})")
        if use_spectral:
            run_cmd(
                [
                    "./gravmag_sphere_spectral",
                    f"{args.rsphere_km}",
                    str(ex),
                    str(out_xyz),
                    f"{args.lmax}",
                    f"{args.refine_factor}",
                    f"{args.ntheta_fit}",
                    f"{args.nphi_fit}",
                    f"{args.reg_lambda}",
                    f"{args.reg_power}",
                    f"{args.source_nlat}",
                    f"{args.source_nlon}",
                    f"{args.source_nr}",
                ],
                cwd,
            )
        elif args.solver == "gauss_legendre":
            # preserve the shared launcher controls: source counts are orders, refinement is panels
            run_cmd(["./gravmag_sphere_gauss_legendre", str(args.rsphere_km), str(ex), str(out_xyz),
                     str(args.source_nr), str(args.source_nlat), str(args.source_nlon), str(args.refine_factor)], cwd)
        else:
            run_cmd(
                [
                    "./gravmag_sphere_point_source",
                    f"{args.rsphere_km}",
                    str(ex),
                    str(out_xyz),
                    f"{args.refine_factor}",
                    f"{args.source_nlat}",
                    f"{args.source_nlon}",
                    f"{args.source_nr}",
                ],
                cwd,
            )

        run_cmd(["./gravmag_xyz_to_brtp", str(out_xyz), str(out_brtp)], cwd)

        data = load_brtp(out_brtp)
        total = aggregate_total(data)
        unit = parse_unit_from_header(out_brtp) or parse_unit_from_fname(out_brtp)
        human = pretty_case_name(stem)
        title = f"{human} [{solver_label}]: Br / Btheta / Bphi / Btot (summed over bodies)"
        plot_brtp(total, out_png, title, unit)
        print(f"[OK ] wrote {out_xyz}, {out_brtp}, {out_png}")

    print("All example cases completed.")


if __name__ == "__main__":
    main()
