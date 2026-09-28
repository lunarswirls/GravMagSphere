#!/usr/bin/env python3
import csv
import math
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent

# Table 4 values extracted from Hemingway and Garrick-Bethell (2012),
# "Magnetic field direction and lunar swirl morphology: Insights from Airy and Reiner Gamma"
# Columns are:
#   latitude_deg_N, longitude_deg_E
# Table-wide constants:
#   depth = 5 km
#   magnetic moment = 1.82e11 Am2
#   local inclination = +2 deg
#   local declination = -8 deg
TABLE4_ROWS = [
    (7.649, 301.460),
    (7.673, 301.418),
    (7.693, 301.371),
    (7.710, 301.320),
    (7.723, 301.264),
    (7.731, 301.205),
    (7.735, 301.144),
    (7.735, 301.080),
    (7.731, 301.015),
    (7.723, 300.950),
    (7.710, 300.885),
    (7.693, 300.821),
    (7.673, 300.759),
    (7.649, 300.700),
    (7.622, 300.645),
    (7.592, 300.593),
    (7.559, 300.546),
    (7.525, 300.505),
    (7.488, 300.469),
    (7.451, 300.439),
    (7.315, 301.525),
    (7.278, 301.495),
    (7.241, 301.459),
    (7.207, 301.418),
    (7.174, 301.371),
    (7.144, 301.319),
    (7.117, 301.264),
    (7.093, 301.205),
    (7.073, 301.143),
    (7.056, 301.079),
    (7.043, 301.014),
    (7.035, 300.949),
    (7.031, 300.884),
    (7.031, 300.820),
    (7.035, 300.759),
    (7.043, 300.700),
    (7.056, 300.644),
    (7.073, 300.593),
    (7.093, 300.546),
    (7.117, 300.504),
    (7.542, 301.774),
    (7.538, 301.819),
    (7.543, 301.878),
    (7.555, 301.943),
    (7.573, 302.001),
    (7.593, 302.043),
    (7.612, 302.062),
    (7.628, 302.056),
    (7.636, 302.024),
    (7.637, 301.973),
    (7.629, 301.910),
    (7.615, 301.847),
    (7.596, 301.794),
    (7.576, 301.759),
    (7.557, 301.749),
]

TABLE_HEADER = ["Latitude", "Longitude", "Depth", "Magnetic Moment", "Inclination", "Declination"]
TABLE_UNITS = ["deg N", "deg E", "km", "Am2", "deg", "deg"]

# Hardcoded run configuration for the Reiner Gamma Table 4 direct_sum-solver case.
RSPHERE_KM = 1737.4
ELVO_KM = 18.0
SURFACE_ELVO_KM = 0.0
OBS_SPACING_DEG = 0.02
LAT_MIN_DEG = 5.5
LAT_MAX_DEG = 9.0
LON_MIN_E_DEG = 299.5
LON_MAX_E_DEG = 303.5
THICKNESS_KM = 1.0
DEPTH_KM = 5.0
MAGNETIC_MOMENT_AM2 = 1.82e11
INCLINATION_LOCAL_DEG = 2.0
DECLINATION_LOCAL_DEG = -8.0


@dataclass
class DipoleRow:
    row_id: int
    latitude_deg_n: float
    longitude_deg_e: float
    longitude_deg_wrapped: float
    depth_km: float
    magnetic_moment_am2: float
    inclination_local_deg: float
    declination_local_deg: float
    inclination_global_deg: float
    declination_global_deg: float
    volume_m3: float
    magnetization_apm: float
    implied_moment_am2: float


def wrap180(lon_deg: float) -> float:
    return ((lon_deg + 180.0) % 360.0) - 180.0


def local_geophysical_to_global_angles(
    lat_deg: float,
    lon_deg_wrapped: float,
    inc_local_deg: float,
    dec_local_deg: float,
) -> tuple[float, float]:
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg_wrapped)
    inc = math.radians(inc_local_deg)
    dec = math.radians(dec_local_deg)

    north = math.cos(inc) * math.cos(dec)
    east = math.cos(inc) * math.sin(dec)
    down = math.sin(inc)

    e = (-math.sin(lon), math.cos(lon), 0.0)
    n = (-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat))
    u = (math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat))

    mx = north * n[0] + east * e[0] - down * u[0]
    my = north * n[1] + east * e[1] - down * u[1]
    mz = north * n[2] + east * e[2] - down * u[2]

    dec_global_deg = math.degrees(math.atan2(my, mx))
    inc_global_deg = math.degrees(math.atan2(mz, math.hypot(mx, my)))
    return inc_global_deg, dec_global_deg


def median_nearest_neighbor_deg(rows: list[tuple[float, float]]) -> float:
    pts = [(lat, wrap180(lon_e)) for lat, lon_e in rows]
    mean_lat = sum(lat for lat, _ in pts) / len(pts)
    scale = math.cos(math.radians(mean_lat))
    nearest = []
    for i, (lat_i, lon_i) in enumerate(pts):
        best = float("inf")
        for j, (lat_j, lon_j) in enumerate(pts):
            if i == j:
                continue
            dlat = lat_j - lat_i
            dlon = (lon_j - lon_i) * scale
            dist = math.hypot(dlat, dlon)
            if dist < best:
                best = dist
        nearest.append(best)
    return float(np.median(nearest))


def source_volume_m3(
    rsphere_km: float,
    depth_km: float,
    lat_deg: float,
    box_width_deg: float,
    thickness_km: float,
) -> float:
    r_m = (rsphere_km - depth_km) * 1000.0
    dlat_rad = math.radians(box_width_deg)
    dlon_rad = math.radians(box_width_deg)
    area_m2 = (r_m * dlat_rad) * (r_m * math.cos(math.radians(lat_deg)) * dlon_rad)
    return area_m2 * thickness_km * 1000.0


def polygon_area_m2(
    vertices_lat_lon: list[tuple[float, float]],
    rsphere_km: float,
    depth_km: float,
) -> float:
    lat0 = float(np.mean([lat for lat, _ in vertices_lat_lon]))
    lon0 = float(np.mean([lon for _, lon in vertices_lat_lon]))
    r_m = (rsphere_km - depth_km) * 1000.0
    x = []
    y = []
    for lat_deg, lon_deg in vertices_lat_lon:
        dlon_deg = wrap180(lon_deg - lon0)
        x.append(r_m * math.cos(math.radians(lat0)) * math.radians(dlon_deg))
        y.append(r_m * math.radians(lat_deg - lat0))
    x.append(x[0])
    y.append(y[0])
    area = 0.0
    for i in range(len(vertices_lat_lon)):
        area += x[i] * y[i + 1] - x[i + 1] * y[i]
    return 0.5 * abs(area)


def build_rows(
    table_rows: list[tuple[float, float]],
    rsphere_km: float,
    box_width_deg: float,
    thickness_km: float,
    depth_km: float,
    magnetic_moment_am2: float,
    inclination_local_deg: float,
    declination_local_deg: float,
    global_inc_deg_uniform: float,
    global_dec_deg_uniform: float,
    magnetization_apm_uniform: float,
) -> list[DipoleRow]:
    out = []
    for idx, (lat_deg, lon_deg_e) in enumerate(table_rows, start=1):
        lon_wrapped = wrap180(lon_deg_e)
        volume_m3 = source_volume_m3(rsphere_km, depth_km, lat_deg, box_width_deg, thickness_km)
        implied_moment_am2 = magnetization_apm_uniform * volume_m3
        out.append(
            DipoleRow(
                row_id=idx,
                latitude_deg_n=lat_deg,
                longitude_deg_e=lon_deg_e,
                longitude_deg_wrapped=lon_wrapped,
                depth_km=depth_km,
                magnetic_moment_am2=magnetic_moment_am2,
                inclination_local_deg=inclination_local_deg,
                declination_local_deg=declination_local_deg,
                inclination_global_deg=global_inc_deg_uniform,
                declination_global_deg=global_dec_deg_uniform,
                volume_m3=volume_m3,
                magnetization_apm=magnetization_apm_uniform,
                implied_moment_am2=implied_moment_am2,
            )
        )
    return out


def write_rows_csv(path: Path, rows: list[DipoleRow]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "row_id",
                "latitude_deg_n",
                "longitude_deg_e",
                "longitude_deg_wrapped",
                "depth_km",
                "magnetic_moment_am2",
                "inclination_local_deg",
                "declination_local_deg",
                "inclination_global_deg",
                "declination_global_deg",
                "volume_m3",
                "magnetization_apm",
                "implied_moment_am2",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.row_id,
                    f"{row.latitude_deg_n:.6f}",
                    f"{row.longitude_deg_e:.6f}",
                    f"{row.longitude_deg_wrapped:.6f}",
                    f"{row.depth_km:.3f}",
                    f"{row.magnetic_moment_am2:.6e}",
                    f"{row.inclination_local_deg:.6f}",
                    f"{row.declination_local_deg:.6f}",
                    f"{row.inclination_global_deg:.6f}",
                    f"{row.declination_global_deg:.6f}",
                    f"{row.volume_m3:.6e}",
                    f"{row.magnetization_apm:.6e}",
                    f"{row.implied_moment_am2:.6e}",
                ]
            )


def evaluation_grid_meta(
    obs_spacing_deg: float,
    lat_min_deg: float,
    lat_max_deg: float,
    lon_min_e_deg: float,
    lon_max_e_deg: float,
) -> dict[str, float]:
    lon0_wrapped = wrap180(lon_min_e_deg)
    lon1_wrapped = wrap180(lon_max_e_deg)
    if lon1_wrapped <= lon0_wrapped:
        raise ValueError("Evaluation longitude window crosses the wrap boundary; use a non-wrapping range.")

    nlat = int(math.floor((lat_max_deg - lat_min_deg) / obs_spacing_deg + 0.5)) + 1
    nlon = int(math.floor((lon1_wrapped - lon0_wrapped) / obs_spacing_deg + 0.5)) + 1
    return {
        "lat_min_deg": lat_min_deg,
        "lat_max_deg": lat_max_deg,
        "lon_min_e_deg": lon_min_e_deg,
        "lon_max_e_deg": lon_max_e_deg,
        "lon_min_wrapped_deg": lon0_wrapped,
        "lon_max_wrapped_deg": lon1_wrapped,
        "nlat": nlat,
        "nlon": nlon,
        "obs_spacing_deg": obs_spacing_deg,
    }


def write_input_file(
    path: Path,
    rows: list[DipoleRow],
    box_width_deg: float,
    thickness_km: float,
    elvo_km: float,
    obs_spacing_deg: float,
    lat_min_deg: float,
    lat_max_deg: float,
    lon_min_e_deg: float,
    lon_max_e_deg: float,
) -> dict[str, float]:
    grid_meta = evaluation_grid_meta(
        obs_spacing_deg,
        lat_min_deg,
        lat_max_deg,
        lon_min_e_deg,
        lon_max_e_deg,
    )
    half = 0.5 * box_width_deg
    half_thickness_km = 0.5 * thickness_km

    with path.open("w") as f:
        for row in rows:
            title = (
                f"Reiner Gamma Table 4 row {row.row_id} "
                f"(lat={row.latitude_deg_n:.3f}, lonE={row.longitude_deg_e:.3f}, "
                f"depth={row.depth_km:.2f} km)"
            )
            f.write(title + "\n")
            f.write(
                f"{grid_meta['lat_min_deg']:.6f}  {grid_meta['lon_min_wrapped_deg']:.6f}  "
                f"{obs_spacing_deg:.6f}  {obs_spacing_deg:.6f}  {elvo_km:.3f}  "
                f"{int(grid_meta['nlat']):d}  {int(grid_meta['nlon']):d}\n"
            )
            f.write("1 4 4 1\n")
            f.write("0.1  1.50  -2.0  2.0  0.080\n")
            f.write(
                f"2 1 {row.magnetization_apm:.8e} "
                f"{row.inclination_global_deg:.6f} {row.declination_global_deg:.6f}\n"
            )
            f.write("0 1 1 0.0 0.0\n")
            f.write(
                f"{row.latitude_deg_n + half:.6f}  {row.latitude_deg_n - half:.6f}  "
                f"{row.longitude_deg_wrapped + half:.6f}  {row.longitude_deg_wrapped - half:.6f}  "
                f"{row.depth_km - half_thickness_km:.6f}  {row.depth_km + half_thickness_km:.6f}\n\n"
            )

    return grid_meta


def write_polygon_input_file(
    path: Path,
    rows: list[DipoleRow],
    polygon_magnetization_apm: float,
    polygon_inc_deg: float,
    polygon_dec_deg: float,
    thickness_km: float,
    elvo_km: float,
    obs_spacing_deg: float,
    lat_min_deg: float,
    lat_max_deg: float,
    lon_min_e_deg: float,
    lon_max_e_deg: float,
) -> dict[str, float]:
    grid_meta = evaluation_grid_meta(
        obs_spacing_deg,
        lat_min_deg,
        lat_max_deg,
        lon_min_e_deg,
        lon_max_e_deg,
    )
    half_thickness_km = 0.5 * thickness_km
    npts = len(rows) + 1

    with path.open("w") as f:
        f.write("Reiner Gamma Table 4 connected polygon approximation\n")
        f.write(
            f"{grid_meta['lat_min_deg']:.6f}  {grid_meta['lon_min_wrapped_deg']:.6f}  "
            f"{obs_spacing_deg:.6f}  {obs_spacing_deg:.6f}  {elvo_km:.3f}  "
            f"{int(grid_meta['nlat']):d}  {int(grid_meta['nlon']):d}\n"
        )
        f.write("1 16 16 0\n")
        f.write("0.1  1.50  -2.0  2.0  0.080\n")
        f.write(f"2 1 {polygon_magnetization_apm:.8e} {polygon_inc_deg:.6f} {polygon_dec_deg:.6f}\n")
        f.write("0 1 1 0.0 0.0\n")
        f.write(f"{npts:d}  {rows[0].depth_km - half_thickness_km:.6f}  {rows[0].depth_km + half_thickness_km:.6f}\n")
        for row in rows:
            f.write(f"{row.latitude_deg_n:.6f}  {row.longitude_deg_wrapped:.6f}\n")
        f.write(f"{rows[0].latitude_deg_n:.6f}  {rows[0].longitude_deg_wrapped:.6f}\n")

    return grid_meta


def run(cmd: list[str], cwd: Path, env: dict[str, str]) -> None:
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def load_brtp_sum(path: Path) -> dict[str, np.ndarray]:
    arr = np.loadtxt(path, comments="#")
    if arr.ndim == 1:
        arr = arr[None, :]
    coords = np.column_stack([arr[:, 1], arr[:, 2]])
    uniq, inv = np.unique(coords, axis=0, return_inverse=True)

    br = np.zeros(len(uniq))
    bt = np.zeros(len(uniq))
    bp = np.zeros(len(uniq))
    np.add.at(br, inv, arr[:, 3])
    np.add.at(bt, inv, arr[:, 4])
    np.add.at(bp, inv, arr[:, 5])
    btot = np.sqrt(br * br + bt * bt + bp * bp)

    lon_u = np.unique(uniq[:, 0])
    lat_u = np.unique(uniq[:, 1])

    def reshape_field(values: np.ndarray) -> np.ndarray:
        grid = np.full((len(lat_u), len(lon_u)), np.nan)
        lon_to_ix = {v: i for i, v in enumerate(lon_u)}
        lat_to_iy = {v: i for i, v in enumerate(lat_u)}
        for (lo, la), val in zip(uniq, values):
            grid[lat_to_iy[la], lon_to_ix[lo]] = val
        return grid

    return {
        "lon": lon_u,
        "lat": lat_u,
        "Br": reshape_field(br),
        "Btheta": reshape_field(bt),
        "Bphi": reshape_field(bp),
        "Btot": reshape_field(btot),
    }


def read_grid_meta_from_input_file(path: Path) -> dict[str, float]:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"Input file {path} does not contain a valid grid header")
    vals = lines[1].split()
    if len(vals) < 7:
        raise ValueError(f"Grid header in {path} is incomplete")

    lat_min_deg = float(vals[0])
    lon_min_wrapped_deg = float(vals[1])
    dlat_deg = float(vals[2])
    dlon_deg = float(vals[3])
    nlat = int(float(vals[5]))
    nlon = int(float(vals[6]))

    return {
        "lat_min_deg": lat_min_deg,
        "lat_max_deg": lat_min_deg + dlat_deg * (nlat - 1),
        "lon_min_wrapped_deg": lon_min_wrapped_deg,
        "lon_max_wrapped_deg": lon_min_wrapped_deg + dlon_deg * (nlon - 1),
    }


def make_field_plot(path: Path, fields: dict[str, np.ndarray], subtitle: str) -> None:
    keys = ["Br", "Btheta", "Bphi", "Btot"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    lon = fields["lon"]
    lat = fields["lat"]

    for ax, key in zip(axes.flat, keys):
        data = fields[key]
        if key == "Btot":
            vmax = float(np.nanmax(data))
            im = ax.imshow(
                data,
                origin="lower",
                aspect="auto",
                extent=[lon.min(), lon.max(), lat.min(), lat.max()],
                cmap="viridis",
                vmin=0.0,
                vmax=vmax,
            )
        else:
            vmax = float(np.nanmax(np.abs(data)))
            im = ax.imshow(
                data,
                origin="lower",
                aspect="auto",
                extent=[lon.min(), lon.max(), lat.min(), lat.max()],
                cmap="RdBu_r",
                vmin=-vmax,
                vmax=vmax,
            )
        ax.set_title(key)
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        fig.colorbar(im, ax=ax, shrink=0.9, label="nT")

    fig.suptitle(f"Reiner Gamma Table 4 direct_sum-solver approximation\n{subtitle}")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def make_source_plot(path: Path, rows: list[DipoleRow], box_width_deg: float, input_path: Path) -> None:
    lon = np.array([r.longitude_deg_wrapped for r in rows])
    lat = np.array([r.latitude_deg_n for r in rows])
    mag = np.array([r.magnetization_apm for r in rows])
    grid_meta = read_grid_meta_from_input_file(input_path)

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    sc = ax.scatter(lon, lat, c=mag, cmap="magma", s=70, edgecolor="black", linewidth=0.3)
    ax.set_title(
        "Reiner Gamma Table 4 source centers\n"
        f"finite-box direct_sum approximation ({box_width_deg:.4f} deg x {box_width_deg:.4f} deg x 1.0 km)"
    )
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    ax.set_xlim(grid_meta["lon_min_wrapped_deg"], grid_meta["lon_max_wrapped_deg"])
    ax.set_ylim(grid_meta["lat_min_deg"], grid_meta["lat_max_deg"])
    fig.colorbar(sc, ax=ax, label="Magnetization [A/m]")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def make_btot_comparison_plot(
    path: Path,
    discrete_fields: dict[str, np.ndarray],
    polygon_fields: dict[str, np.ndarray],
) -> None:
    discrete_btot = discrete_fields["Btot"]
    polygon_btot = polygon_fields["Btot"]
    diff = polygon_btot - discrete_btot
    lon = discrete_fields["lon"]
    lat = discrete_fields["lat"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    vmax = float(max(np.nanmax(discrete_btot), np.nanmax(polygon_btot)))
    diff_abs = float(np.nanmax(np.abs(diff)))

    im0 = axes[0].imshow(
        discrete_btot,
        origin="lower",
        aspect="auto",
        extent=[lon.min(), lon.max(), lat.min(), lat.max()],
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
    )
    axes[0].set_title("Discrete Btot")
    axes[0].set_xlabel("Longitude [deg]")
    axes[0].set_ylabel("Latitude [deg]")
    fig.colorbar(im0, ax=axes[0], shrink=0.9, label="nT")

    im1 = axes[1].imshow(
        polygon_btot,
        origin="lower",
        aspect="auto",
        extent=[lon.min(), lon.max(), lat.min(), lat.max()],
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
    )
    axes[1].set_title("Polygon Btot")
    axes[1].set_xlabel("Longitude [deg]")
    axes[1].set_ylabel("Latitude [deg]")
    fig.colorbar(im1, ax=axes[1], shrink=0.9, label="nT")

    im2 = axes[2].imshow(
        diff,
        origin="lower",
        aspect="auto",
        extent=[lon.min(), lon.max(), lat.min(), lat.max()],
        cmap="RdBu_r",
        vmin=-diff_abs,
        vmax=diff_abs,
    )
    axes[2].set_title("Polygon - Discrete Btot")
    axes[2].set_xlabel("Longitude [deg]")
    axes[2].set_ylabel("Latitude [deg]")
    fig.colorbar(im2, ax=axes[2], shrink=0.9, label="nT")

    fig.suptitle("Reiner Gamma Table 4 comparison: connected polygon vs discrete sources")
    fig.savefig(path, dpi=220)
    plt.close(fig)


HERE.mkdir(parents=True, exist_ok=True)
(HERE / ".mplconfig").mkdir(exist_ok=True)

box_width_deg = 0.35 * median_nearest_neighbor_deg(TABLE4_ROWS)
mean_lat_deg = float(np.mean([lat for lat, _ in TABLE4_ROWS]))
mean_lon_wrapped_deg = float(np.mean([wrap180(lon_e) for _, lon_e in TABLE4_ROWS]))
uniform_inc_global_deg, uniform_dec_global_deg = local_geophysical_to_global_angles(
    mean_lat_deg,
    mean_lon_wrapped_deg,
    INCLINATION_LOCAL_DEG,
    DECLINATION_LOCAL_DEG,
)
ref_volume_m3 = source_volume_m3(
    RSPHERE_KM,
    DEPTH_KM,
    mean_lat_deg,
    box_width_deg,
    THICKNESS_KM,
)
uniform_magnetization_apm = MAGNETIC_MOMENT_AM2 / ref_volume_m3
rows = build_rows(
    TABLE4_ROWS,
    rsphere_km=RSPHERE_KM,
    box_width_deg=box_width_deg,
    thickness_km=THICKNESS_KM,
    depth_km=DEPTH_KM,
    magnetic_moment_am2=MAGNETIC_MOMENT_AM2,
    inclination_local_deg=INCLINATION_LOCAL_DEG,
    declination_local_deg=DECLINATION_LOCAL_DEG,
    global_inc_deg_uniform=uniform_inc_global_deg,
    global_dec_deg_uniform=uniform_dec_global_deg,
    magnetization_apm_uniform=uniform_magnetization_apm,
)

rows_csv = HERE / "reiner_gamma_table4_rows.csv"
input_path = HERE / "reiner_gamma_table4_direct.in"
surface_input_path = HERE / "reiner_gamma_table4_surface_direct.in"
xyz_path = HERE / "output" / "reiner_gamma_table4_direct_xyz.txt"
brtp_path = HERE / "output" / "reiner_gamma_table4_direct_brtp.txt"
field_plot = HERE / "figs" / "reiner_gamma_table4_direct_brtp_2x2.png"
source_plot = HERE / "figs" / "reiner_gamma_table4_source_layout.png"
xyz_path.parent.mkdir(parents=True, exist_ok=True)
field_plot.parent.mkdir(parents=True, exist_ok=True)

write_rows_csv(rows_csv, rows)
write_input_file(
    input_path,
    rows,
    box_width_deg=box_width_deg,
    thickness_km=THICKNESS_KM,
    elvo_km=ELVO_KM,
    obs_spacing_deg=OBS_SPACING_DEG,
    lat_min_deg=LAT_MIN_DEG,
    lat_max_deg=LAT_MAX_DEG,
    lon_min_e_deg=LON_MIN_E_DEG,
    lon_max_e_deg=LON_MAX_E_DEG,
)
write_input_file(
    surface_input_path,
    rows,
    box_width_deg=box_width_deg,
    thickness_km=THICKNESS_KM,
    elvo_km=SURFACE_ELVO_KM,
    obs_spacing_deg=OBS_SPACING_DEG,
    lat_min_deg=LAT_MIN_DEG,
    lat_max_deg=LAT_MAX_DEG,
    lon_min_e_deg=LON_MIN_E_DEG,
    lon_max_e_deg=LON_MAX_E_DEG,
)

env = os.environ.copy()

run([str(ROOT / "build_gravmag_tools.sh")], cwd=ROOT, env=env)
run(
    [
        str(ROOT / "gravmag_sphere_direct_sum"),
        f"{RSPHERE_KM}",
        str(input_path),
        str(xyz_path),
        "1",
        "0",
        "0",
        "0",
    ],
    cwd=ROOT,
    env=env,
)
run([str(ROOT / "gravmag_xyz_to_brtp"), str(xyz_path), str(brtp_path)], cwd=ROOT, env=env)

fields = load_brtp_sum(brtp_path)
subtitle = (
    f"55 Table-4 dipoles approximated as {box_width_deg:.4f} deg square x "
    f"{THICKNESS_KM:.1f} km boxes, altitude {ELVO_KM:.1f} km, "
    f"window {LAT_MIN_DEG:.1f}-{LAT_MAX_DEG:.1f} N / "
    f"{LON_MIN_E_DEG:.1f}-{LON_MAX_E_DEG:.1f} E"
)
make_field_plot(field_plot, fields, subtitle)
make_source_plot(source_plot, rows, box_width_deg, input_path)

print(f"Wrote: {rows_csv}")
print(f"Wrote: {input_path}")
print(f"Wrote: {surface_input_path}")
print(f"Wrote: {xyz_path}")
print(f"Wrote: {brtp_path}")
print(f"Wrote: {field_plot}")
print(f"Wrote: {source_plot}")
