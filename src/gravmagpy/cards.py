"""read supported fortran input cards and export fitted spherical blocks"""

from pathlib import Path

import numpy as np

from .geometry import source_array
from .utils import artifact_path, run_fortran
from .utils.fortran import _solver_name


def read_input(path):
    """read every magnetic or gravity body in a fortran .in card file"""
    lines = [line.strip() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith(("#", "!"))]
    bodies, index = [], 0
    while index < len(lines):
        try:
            title = lines[index]
            grid = list(map(float, lines[index + 1].split()))
            mesh = list(map(int, lines[index + 2].split()))
            integration = list(map(float, lines[index + 3].split()))
            property_values = list(map(float, lines[index + 4].split()))
            output = list(map(float, lines[index + 5].split()))
            geometry = list(map(float, lines[index + 6].split()))
            if len(grid) != 7 or len(mesh) != 4 or any(len(card) != 5 for card in (integration, property_values, output)):
                raise ValueError("invalid card length")
            if mesh[3] == 1:
                lat_max, lat_min, lon_max, lon_min, top, bottom = geometry
                polygon = [[lat_min, lon_min], [lat_min, lon_max], [lat_max, lon_max],
                           [lat_max, lon_min], [lat_min, lon_min]]
                index += 7
            elif mesh[3] == 0:
                count, top, bottom = geometry
                if count != int(count) or count < 3:
                    raise ValueError("invalid polygon vertex count")
                count = int(count)
                polygon = [list(map(float, lines[index + 7 + vertex].split())) for vertex in range(count)]
                if any(len(vertex) != 2 for vertex in polygon):
                    raise ValueError("polygon rows need latitude longitude")
                polygon = [vertex[::-1] if abs(vertex[0]) > 90 and abs(vertex[1]) <= 90 else vertex
                           for vertex in polygon]
                index += 7 + count
            else:
                raise ValueError("nblim must be 0 or 1")
            if bottom <= top or top < 0:
                raise ValueError("invalid source depths")
            bodies.append({"title": title, "grid": grid, "mesh": mesh, "property": property_values,
                           "integration": integration, "output": output, "geometry": geometry,
                           "polygon_lat_lon": polygon, "depth_top_km": top, "depth_bottom_km": bottom})
        except (IndexError, ValueError) as error:
            raise ValueError(f"invalid input block near noncomment line {index + 1}: {error}") from error
    if not bodies:
        raise ValueError("input file contains no bodies")
    return bodies


def write_source_input(sources, path, *, radius_km=1737.4, grid=None, mesh=(4, 16, 16)):
    """export fitted spherical blocks as fixed-limit .in bodies

    grid is (lat0, lon0, dlat, dlon, observation_altitude_km, nlat, nlon)
    the orbital fit never uses this constant-altitude grid; it is for mapping
    """
    rows = source_array(sources, radius_km)
    if grid is None:
        lat_low = np.min(rows[:, 0] - rows[:, 2]/2) - 1
        lon_low = np.min(rows[:, 1] - rows[:, 3]/2) - 1
        lat_high = np.max(rows[:, 0] + rows[:, 2]/2) + 1
        lon_high = np.max(rows[:, 1] + rows[:, 3]/2) + 1
        grid = (max(-89, lat_low), lon_low, (min(89, lat_high)-max(-89, lat_low))/60,
                (lon_high-lon_low)/60, 30.0, 61, 61)
    if len(grid) != 7 or not np.isfinite(grid).all() or any(int(value) != value or value < 1 for value in grid[5:]):
        raise ValueError("grid must contain five finite values and two positive integer counts")
    if len(mesh) != 3 or any(int(value) != value or value < 1 for value in mesh):
        raise ValueError("mesh must contain three positive integer counts")
    blocks = []
    for index, row in enumerate(rows):
        amplitude = np.linalg.norm(row[6:9])
        inclination = np.rad2deg(np.arctan2(row[8], np.hypot(row[6], row[7])))
        declination = np.rad2deg(np.arctan2(row[7], row[6])) % 360
        blocks.append("\n".join([
            f"fitted source {index + 1}",
            " ".join(f"{value:.12g}" for value in grid[:5]) + f" {int(grid[5])} {int(grid[6])}",
            " ".join(str(int(value)) for value in mesh) + " 1",
            "0.1 1.5 -2 2 0.08",
            f"2 1 {amplitude:.12g} {inclination:.12g} {declination:.12g}",
            "0 1 1 0 0",
            " ".join(f"{value:.12g}" for value in [row[0]+row[2]/2, row[0]-row[2]/2,
                     row[1]+row[3]/2, row[1]-row[3]/2, row[4], row[4]+row[5]]),
        ]))
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n\n".join(blocks) + "\n", encoding="utf-8")
    return path


def run_grid_model(input_path, output_path=None, *, radius_km=1737.4, solver="point_source", options=(), **build_options):
    """run an existing point_source, spectral, or gauss_legendre .in model"""
    solver = _solver_name(solver)
    suffixes = {"point_source": ".txt", "spectral": "_spectral.txt", "gauss_legendre": "_quadrature.txt"}
    if solver not in suffixes:
        raise ValueError("solver must be point_source, spectral, or gauss_legendre")
    input_path = Path(input_path).expanduser().resolve()
    output_path = (artifact_path(input_path, suffix=suffixes[solver])
                   if output_path is None else Path(output_path).expanduser().resolve())
    read_input(input_path)
    if input_path == output_path:
        raise ValueError("input and output paths must differ")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result = run_fortran(solver, [radius_km, input_path, output_path, *options], **build_options)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"solver did not write {output_path}")
    return {"output_path": output_path, "stdout": result.stdout, "stderr": result.stderr}
