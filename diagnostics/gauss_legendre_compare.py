"""compare quadrature, direct_sum and pure spectral fields at matched grid points"""

import csv
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter

import numpy as np

from gravmagpy import block, polygon, sphere_model, run_sphere_model, write_model_input, build_fortran


# edit these settings to extend the source, altitude or resolution sweep
settings = {
    "radius_km": 1737.4,
    "latitude_deg": np.linspace(5.5, 8.5, 7).tolist(),
    "longitude_deg": np.linspace(-60.5, -57.5, 7).tolist(),
    "orders": [4, 8, 16, 32, 64],
    "refine_factors": [1, 2, 4],
    "degrees": [4, 8, 12, 24],
    "reference_options": {"radial_order": 16, "latitude_order": 48, "longitude_order": 48, "subdivisions": 4},
    "check_options": {"radial_order": 12, "latitude_order": 32, "longitude_order": 32, "subdivisions": 4},
    "reference_tolerance": 1e-4,
}


def compare_fields(field, reference):
    """use vector differences so opposing component residuals cannot cancel"""
    difference = field-reference
    return {"relative_l2": float(np.linalg.norm(difference)/np.linalg.norm(reference)),
            "component_rmse": np.sqrt(np.mean(difference**2, axis=(0, 1))).tolist(),
            "vector_rmse": float(np.sqrt(np.mean(np.sum(difference**2, axis=-1)))),
            "max_vector_difference": float(np.max(np.linalg.norm(difference, axis=-1)))}


if "--help" in sys.argv or "-h" in sys.argv:
    print("usage: python diagnostics/gauss_legendre_compare.py [output_directory]\n"
          "runs magnetic and gravity block/polygon comparisons at surface, orbital and far-field altitudes\n"
          "edit the settings dictionary to change the sweep; default output: diagnostics/gauss_legendre_comparison\n"
          "all numeric residuals use summed xyz vectors; pure spectral disables auto, edge and hybrid corrections")
    raise SystemExit(0)
if len(sys.argv) > 2:
    raise SystemExit("supply at most one output directory")

root = Path(__file__).resolve().parents[1]
destination = Path(sys.argv[1]).expanduser().resolve() if len(sys.argv) == 2 else root/"diagnostics/gauss_legendre_comparison"
destination.mkdir(parents=True, exist_ok=True)
os.environ["GRAVMAG_DIAGNOSTICS"] = "1"
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# compile before timing so cache misses do not distort solver comparisons
binaries = {solver: build_fortran(solver) for solver in ("gauss_legendre", "direct_sum", "spectral")}
compiler = subprocess.run([os.environ.get("FC", "gfortran"), "--version"],
                          check=True, capture_output=True, text=True).stdout.splitlines()[0]
cases = []
for name, kind, altitude, geometry in (
        ("magnetic_block_surface", "magnetic", 0, "block"),
        ("magnetic_block_orbital", "magnetic", 30, "block"),
        ("magnetic_concave_polygon_orbital", "magnetic", 30, "polygon"),
        ("gravity_sloping_polygon_orbital", "gravity", 30, "triangle"),
        ("magnetic_block_far", "magnetic", 5000, "block"),
        ("gravity_block_far", "gravity", 5000, "block")):
    material = {"magnetization_a_m": [1, -2, 3]} if kind == "magnetic" else {"density_kg_m3": 400}
    if geometry == "block":
        source = block([6, 8], [-60, -58], [2, 8], **material)
    elif geometry == "polygon":
        source = polygon([[6, -60], [6, -58], [6.5, -58], [6.5, -59],
                          [7.5, -59], [7.5, -58], [8, -58], [8, -60]], [2, 8], **material)
    else:
        source = polygon([[6, -60], [8, -59.5], [7, -58]], [2, 8], **material)
    model = sphere_model([source], radius_km=settings["radius_km"], altitude_km=altitude,
                         latitude_deg=settings["latitude_deg"], longitude_deg=settings["longitude_deg"])
    cases.append({"name": name, "model": model})

rows = []
references = []
arrays = {}
for case in cases:
    name, model = case["name"], case["model"]
    print(f"comparing {name}", flush=True)
    write_model_input(model, destination/f"{name}.in")
    reference = run_sphere_model(model, solver="gauss_legendre", solver_options=settings["reference_options"])
    check = run_sphere_model(model, solver="gauss_legendre", solver_options=settings["check_options"])
    metrics = compare_fields(check["field"], reference["field"])
    references.append({"case": name, "unit": reference["unit"], "peak": float(np.max(reference["total"])),
                       "check_relative_l2": metrics["relative_l2"],
                       "converged": metrics["relative_l2"] < settings["reference_tolerance"]})
    arrays[f"{name}__reference"] = reference["field"]
    arrays[f"{name}__reference_check"] = check["field"]
    runs = [("gauss_legendre", f"order_{order}", {"radial_order": order, "latitude_order": order,
              "longitude_order": order}) for order in settings["orders"]]
    runs += [("gauss_legendre", "order_8_panels_4", {"radial_order": 8, "latitude_order": 8,
               "longitude_order": 8, "subdivisions": 4})]
    runs += [("direct_sum", f"refine_{refine}", {"refine_factor": refine}) for refine in settings["refine_factors"]]
    runs += [("spectral", f"degree_{degree}", {"lmax": degree, "ntheta_fit": max(24, 2*degree+2),
              "nphi_fit": max(48, 4*degree+4), "auto_mode": 0, "edge_correction": 0, "hybrid_mode": 0,
              "reg_lambda": 0, "refine_factor": 4, "source_nr": 12}) for degree in settings["degrees"]]
    for solver, label, options in runs:
        start = perf_counter()
        result = run_sphere_model(model, solver=solver, solver_options=options)
        elapsed = perf_counter()-start
        metrics = compare_fields(result["field"], reference["field"])
        rows.append({"case": name, "solver": solver, "setting": label, "unit": result["unit"],
                     "elapsed_s": elapsed, **metrics, "options": options,
                     "diagnostics": result["stdout"].splitlines(), "stderr": result["stderr"]})
        arrays[f"{name}__{solver}__{label}"] = result["field"]

np.savez_compressed(destination/"fields.npz", **arrays)
report = {"settings": settings, "platform": platform.platform(), "python": platform.python_version(),
          "numpy": np.__version__, "compiler": compiler, "harmonic_normalization": "schmidt", "condon_shortley_phase": False,
          "binaries": {key: str(value) for key, value in binaries.items()},
          "cases": cases, "references": references, "runs": rows}
(destination/"results.json").write_text(json.dumps(report, indent=2)+"\n")
columns = ["case", "solver", "setting", "unit", "relative_l2", "vector_rmse", "max_vector_difference", "elapsed_s"]
with (destination/"results.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

lines = ["# Gauss–Legendre, direct_sum and pure spectral comparison", "",
         "Residuals are vector L2 differences divided by the reference vector L2 norm over identical observation grids. "
         "The reference is higher-order composite Gauss–Legendre integration, checked against a second resolution; "
         "it is a numerical reference, not an exact solution. Independent kernel, geometry and unit checks live in "
         "`tests/test_gauss_legendre.py`.", "",
         "Pure spectral runs use Schmidt semi-normalization without the Condon-Shortley phase and set `auto_mode=0`, `edge_correction=0`, `hybrid_mode=0` and `reg_lambda=0`. "
         "Their source mesh uses horizontal refine factor 4 and 12 radial samples. Direct summation runs retain card radial count 4. "
         "Each timing is one warmed-build Python call including subprocess startup and file I/O; it is not a repeated benchmark.", "",
         "The cases span magnetic surface/orbital/far fields, a concave magnetic polygon and gravity. "
         "Low harmonic degrees cannot resolve the narrow near-surface sources. Spectral gravity includes the degree-zero "
         "mass term from the signed sum of discretized source masses, with degrees 1 through lmax fitted to the remaining field. "
         "The direct_sum magnetic solver orients side-wall normals using polygon winding, including at concave edges. "
         "The notched footprint is independently checked against a sum of nonoverlapping rectangles in the tests.", "",
         f"Compiler: {compiler}. Platform: {platform.platform()}.", "",
         "## Reference convergence", "", "| Case | Unit | Peak magnitude | Reference change | Converged |",
         "|---|---|---:|---:|---|"]
for item in references:
    lines.append(f"| {item['case']} | {item['unit']} | {item['peak']:.6g} | {item['check_relative_l2']:.3e} | {item['converged']} |")
lines += ["", "## Resolution sweeps", "", "| Case | Solver | Setting | Relative L2 | Vector RMSE | Time (s) |",
          "|---|---|---|---:|---:|---:|"]
for row in rows:
    lines.append(f"| {row['case']} | {row['solver']} | {row['setting']} | {row['relative_l2']:.3e} | "
                 f"{row['vector_rmse']:.3e} | {row['elapsed_s']:.4f} |")
lines += ["", "![Convergence and cost](convergence.png)", "", "Exact models, options, per-component RMSE and solver diagnostics "
          "are saved in `results.json`; summed XYZ fields are saved in `fields.npz`. Residual units follow each case. "
          "`results.csv` provides the compact sweep table. Plot labels show the first and final node orders, direct_sum refine "
          "factors or spectral degrees. The separate diamond labelled 8×4 uses order 8 with four panels in each dimension.", ""]
(destination/"README.md").write_text("\n".join(lines))

figure, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
for axis, case in zip(axes.flat, cases):
    for solver, marker in (("gauss_legendre", "o"), ("direct_sum", "s"), ("spectral", "^")):
        selected = [row for row in rows if row["case"] == case["name"] and row["solver"] == solver
                    and "panels" not in row["setting"]]
        curves = axis.loglog([row["elapsed_s"] for row in selected], [max(row["relative_l2"], 1e-15) for row in selected],
                             marker=marker, label=solver)
        for index, row in enumerate((selected[0], selected[-1])):
            label = row["setting"].split("_")[-1]
            axis.annotate(label, (row["elapsed_s"], max(row["relative_l2"], 1e-15)),
                          xytext=((-6, 7) if index == 0 else (5, -9)), ha="right" if index == 0 else "left",
                          textcoords="offset points", fontsize=7)
        if solver == "gauss_legendre":
            composite = next(row for row in rows if row["case"] == case["name"] and "panels" in row["setting"])
            point = (composite["elapsed_s"], max(composite["relative_l2"], 1e-15))
            axis.plot(*point, marker="D", markerfacecolor="none", color=curves[0].get_color())
            axis.annotate("8×4", point, xytext=(4, 10), textcoords="offset points", fontsize=7)
    axis.set(title=case["name"].replace("_", " "), xlabel="elapsed time (s)", ylabel="relative vector L2 difference")
    axis.grid(True, which="both", alpha=0.2)
    axis.margins(x=0.12, y=0.15)
axes.flat[0].legend(fontsize=8)
figure.savefig(destination/"convergence.png", dpi=180)
plt.close(figure)
print(f"saved comparison to {destination}")
if not all(item["converged"] for item in references):
    raise SystemExit("reference convergence check failed; increase reference and check resolution before interpreting residuals")
