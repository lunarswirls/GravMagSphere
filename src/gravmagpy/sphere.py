"""dictionary-based models for the supported gravmag sphere tools"""

from copy import deepcopy
from pathlib import Path
import tempfile

import numpy as np

from .cards import read_input, run_grid_model
from .utils.fortran import _solver_name


def _values(values, name, shape=None):
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all() or (shape is not None and values.shape != shape):
        raise ValueError(f"{name} must contain finite values with shape {shape}")
    return values


def _interval(values, name):
    values = _values(values, name, (2,))
    if values[0] >= values[1]:
        raise ValueError(f"{name} must be an increasing pair")
    return values.tolist()


def _material(magnetization_a_m, density_kg_m3):
    if (magnetization_a_m is None) == (density_kg_m3 is None):
        raise ValueError("supply exactly one of magnetization_a_m or density_kg_m3")
    if magnetization_a_m is not None:
        return {"magnetization_a_m": _values(magnetization_a_m, "magnetization_a_m", (3,)).tolist()}
    return {"density_kg_m3": float(_values(density_kg_m3, "density_kg_m3", ()))}


def block(latitude_deg, longitude_deg, depth_km, *, magnetization_a_m=None, density_kg_m3=None,
          name="block", mesh=None):
    """define a spherical block with increasing bounds and xyz magnetization or density"""
    source = {"name": name, "geometry": "block", "latitude_deg": _interval(latitude_deg, "latitude_deg"),
              "longitude_deg": _interval(longitude_deg, "longitude_deg"),
              "depth_km": _interval(depth_km, "depth_km"),
              **_material(magnetization_a_m, density_kg_m3)}
    if mesh is not None:
        source["mesh"] = dict(mesh)
    return source


def polygon(vertices_lat_lon, depth_km, *, magnetization_a_m=None, density_kg_m3=None,
            name="polygon", mesh=None):
    """define a polygon footprint extruded between top and bottom depths"""
    vertices = _values(vertices_lat_lon, "vertices_lat_lon")
    if vertices.ndim != 2 or vertices.shape[1] != 2 or len(vertices) < 3:
        raise ValueError("vertices_lat_lon must contain at least three latitude/longitude pairs")
    source = {"name": name, "geometry": "polygon", "vertices_lat_lon": vertices.tolist(),
              "depth_km": _interval(depth_km, "depth_km"),
              **_material(magnetization_a_m, density_kg_m3)}
    if mesh is not None:
        source["mesh"] = dict(mesh)
    return source


def equivalent_source_grid(latitude_edges_deg, longitude_edges_deg, depth_edges_km, *,
                           magnetization_a_m=None, density_kg_m3=None, mesh=None):
    """generate volume cells ordered by depth, latitude, then longitude

    properties broadcast to (ndepth, nlat, nlon, 3) for xyz magnetization or
    (ndepth, nlat, nlon) for density; a single vector/scalar makes a uniform layer
    these are finite volumes, not point dipoles with magnetic moments
    """
    edges = [_values(values, name) for values, name in zip(
        (depth_edges_km, latitude_edges_deg, longitude_edges_deg), ("depth edges", "latitude edges", "longitude edges"))]
    if any(values.ndim != 1 or len(values) < 2 or np.any(np.diff(values) <= 0) for values in edges):
        raise ValueError("cell edges must be strictly increasing one-dimensional arrays")
    shape = tuple(len(values)-1 for values in edges)
    if (magnetization_a_m is None) == (density_kg_m3 is None):
        raise ValueError("supply exactly one of magnetization_a_m or density_kg_m3")
    magnetic = magnetization_a_m is not None
    name = "magnetization_a_m" if magnetic else "density_kg_m3"
    values = _values(magnetization_a_m if magnetic else density_kg_m3, name)
    values = np.broadcast_to(values, shape + (3,) if magnetic else shape)
    sources = []
    for layer, lat, lon in np.ndindex(shape):
        sources.append(block(edges[1][lat:lat+2], edges[2][lon:lon+2], edges[0][layer:layer+2],
                             **{name: values[layer, lat, lon]}, mesh=mesh, name=f"cell {layer} {lat} {lon}"))
    return sources


def _axis(values, name):
    values = _values(values, name)
    if values.ndim != 1 or len(values) < 1:
        raise ValueError(f"{name} must be a nonempty one-dimensional array")
    if len(values) > 1:
        spacing = np.diff(values)
        if spacing[0] <= 0 or not np.allclose(spacing, spacing[0], rtol=1e-9, atol=1e-10):
            raise ValueError(f"{name} must be uniformly spaced and increasing")
        if spacing[0] < 1e-5:
            raise ValueError("grid spacing must be at least 1e-5 degrees for the fortran text output")
    return values


def sphere_model(sources, *, latitude_deg, longitude_deg, altitude_km=30.0, radius_km=1737.4):
    """create a validated model dictionary with a common regular observation grid"""
    model = {"sources": deepcopy(list(sources)), "radius_km": float(radius_km),
             "grid": {"latitude_deg": np.asarray(latitude_deg).tolist(),
                      "longitude_deg": np.asarray(longitude_deg).tolist(), "altitude_km": float(altitude_km)}}
    _model_cards(model)
    return model


def _model_cards(model):
    if set(model) != {"sources", "grid", "radius_km"}:
        raise ValueError("model requires sources, grid, and radius_km")
    radius = float(_values(model["radius_km"], "radius_km", ()))
    grid = model["grid"]
    if set(grid) != {"latitude_deg", "longitude_deg", "altitude_km"}:
        raise ValueError("grid requires latitude_deg, longitude_deg, and altitude_km")
    lat, lon = [_axis(grid[key], key) for key in ("latitude_deg", "longitude_deg")]
    altitude = float(_values(grid["altitude_km"], "altitude_km", ()))
    if radius <= 0 or altitude < 0 or np.any(np.abs(lat) > 90) or lon[-1]-lon[0] >= 360:
        raise ValueError("invalid radius, altitude, latitude, or duplicate global longitude endpoint")
    if not model["sources"]:
        raise ValueError("at least one source is required")
    grid_card = [lat[0], lon[0], lat[1]-lat[0] if len(lat) > 1 else 1.0,
                 lon[1]-lon[0] if len(lon) > 1 else 1.0, altitude, len(lat), len(lon)]
    cards, fields = [], set()
    for source in model["sources"]:
        allowed = {"name", "geometry", "latitude_deg", "longitude_deg", "vertices_lat_lon", "depth_km",
                   "magnetization_a_m", "density_kg_m3", "mesh", "card4", "card6"}
        if set(source) - allowed:
            raise ValueError(f"unknown source keys: {sorted(set(source)-allowed)}")
        name = source.get("name", "source")
        if not isinstance(name, str) or not name.strip() or len(name.encode()) > 255 or "\n" in name or "\r" in name or name.lstrip().startswith(("#", "!")):
            raise ValueError("source name must fit on a single noncomment fortran title line")
        top, bottom = _interval(source["depth_km"], "depth_km")
        if top < 0 or bottom >= radius:
            raise ValueError("source depths must satisfy 0 <= top < bottom < radius_km")
        material = _material(source.get("magnetization_a_m"), source.get("density_kg_m3"))
        if "magnetization_a_m" in material:
            mx, my, mz = material["magnetization_a_m"]
            property_card = [2, 1, np.linalg.norm([mx, my, mz]), np.rad2deg(np.arctan2(mz, np.hypot(mx, my))),
                             np.rad2deg(np.arctan2(my, mx))]
            fields.add("magnetic")
        else:
            property_card = [1, 0, material["density_kg_m3"], 0, 0]
            fields.add("gravity")
        vertices = []
        if source["geometry"] == "block":
            low_lat, high_lat = _interval(source["latitude_deg"], "latitude_deg")
            low_lon, high_lon = _interval(source["longitude_deg"], "longitude_deg")
            if low_lat <= -90 or high_lat >= 90 or high_lon-low_lon >= 180:
                raise ValueError("blocks must avoid poles and span less than 180 degrees longitude; split larger regions")
            if "vertices_lat_lon" in source:
                raise ValueError("blocks cannot also declare polygon vertices")
            geometry_card = [high_lat, low_lat, high_lon, low_lon, top, bottom]
            flag = 1
        elif source["geometry"] == "polygon":
            if "latitude_deg" in source or "longitude_deg" in source:
                raise ValueError("polygons use vertices_lat_lon, not block bounds")
            vertices = _values(source["vertices_lat_lon"], "vertices_lat_lon")
            if vertices.ndim != 2 or vertices.shape[1] != 2 or len(vertices) < 3:
                raise ValueError("polygons require at least three latitude/longitude pairs")
            vertices = vertices.copy()
            if np.array_equal(vertices[0], vertices[-1]):
                vertices = vertices[:-1]
            vertices[:, 1] = np.rad2deg(np.unwrap(np.deg2rad(vertices[:, 1])))
            if len(vertices) < 3 or np.any(np.abs(vertices[:, 0]) >= 90) or np.ptp(vertices[:, 1]) >= 180:
                raise ValueError("polygons must avoid poles and span less than 180 degrees longitude")
            area = np.sum(vertices[:, 0]*np.roll(vertices[:, 1], 1) - vertices[:, 1]*np.roll(vertices[:, 0], 1))
            if abs(area) < 1e-12:
                raise ValueError("polygon has zero signed area")
            geometry_card, flag = [len(vertices), top, bottom], 0
        else:
            raise ValueError("geometry must be block or polygon")
        mesh = {"radial": 4, "latitude": 16, "longitude": 16}
        if set(source.get("mesh", {})) - set(mesh):
            raise ValueError("mesh accepts radial, latitude, and longitude counts")
        mesh.update(source.get("mesh", {}))
        if any(not np.isfinite(value) or int(value) != value or value < 1 for value in mesh.values()):
            raise ValueError("mesh counts must be positive integers")
        card4 = _values(source.get("card4", [0.1, 1.5, -2, 2, 0.08]), "card4", (5,))
        card6 = _values(source.get("card6", [0, 1, 1, 0, 0]), "card6", (5,))
        if np.any(card6[:3] != card6[:3].astype(int)):
            raise ValueError("the first three card6 values must be integers")
        # integer-valued cards must use integer tokens for list-directed fortran reads
        lines = [name, " ".join(f"{value:.12g}" for value in grid_card),
                 " ".join(str(int(value)) for value in [*mesh.values(), flag])]
        lines += [" ".join(f"{value:.12g}" for value in values)
                  for values in (card4, property_card, card6, geometry_card, *vertices)]
        cards.append("\n".join(lines))
    if len(fields) != 1:
        raise ValueError("use separate models for magnetic and gravity sources; their units cannot be summed")
    return "\n\n".join(cards) + "\n", fields.pop()


def write_model_input(model, path):
    """export a python model as backwards-compatible fortran input cards"""
    content, field = _model_cards(model)
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def read_model_input(path, *, radius_km=1737.4):
    """import supported common-grid block/polygon cards into editable dictionaries

    titles and the compatibility fields on cards 4 and 6 are preserved
    use read_input/run_grid_model for files outside the modern model constraints
    """
    bodies = read_input(path)
    grid = bodies[0]["grid"]
    sources = []
    for body in bodies:
        if body["grid"] != grid:
            raise ValueError("model import requires the same observation grid for every body")
        mode, remanence, amplitude, inc, dec = body["property"]
        if mode == 2 and remanence == 1:
            inc, dec = np.deg2rad([inc, dec])
            material = {"magnetization_a_m": amplitude*np.array([np.cos(inc)*np.cos(dec), np.cos(inc)*np.sin(dec), np.sin(inc)])}
        elif mode == 1:
            material = {"density_kg_m3": amplitude}
        else:
            raise ValueError("supported fields are gravity and uniform remanent magnetization")
        common = {**material, "name": body["title"], "mesh": dict(zip(("radial", "latitude", "longitude"), body["mesh"][:3]))}
        depth = [body["depth_top_km"], body["depth_bottom_km"]]
        if body["mesh"][3] == 1:
            high_lat, low_lat, high_lon, low_lon = body["geometry"][:4]
            source = block([low_lat, high_lat], [low_lon, high_lon], depth, **common)
        else:
            source = polygon(body["polygon_lat_lon"], depth, **common)
        source.update(card4=body["integration"], card6=body["output"])
        sources.append(source)
    lat0, lon0, dlat, dlon, altitude, nlat, nlon = grid
    if int(nlat) != nlat or int(nlon) != nlon:
        raise ValueError("observation counts must be integers")
    return sphere_model(sources, latitude_deg=lat0 + dlat*np.arange(int(nlat)),
                        longitude_deg=lon0 + dlon*np.arange(int(nlon)), altitude_km=altitude, radius_km=radius_km)


def _solver_arguments(solver, options):
    direct_sum = {"refine_factor": 2, "source_nlat": 0, "source_nlon": 0, "source_nr": 0}
    spectral = {"lmax": 24, "refine_factor": 2, "ntheta_fit": 72, "nphi_fit": 144,
                "reg_lambda": 0.2, "reg_power": 4.0, "source_nlat": 0, "source_nlon": 0, "source_nr": 0,
                "auto_mode": 1, "joint_strength": 1.0, "edge_correction": 1, "hybrid_mode": 1,
                "hybrid_band_deg": 1.5, "complex_vertex_threshold": 12, "hybrid_transition_deg": 0.75}
    gauss_legendre = {"radial_order": 0, "latitude_order": 0, "longitude_order": 0, "subdivisions": 1}
    solvers = {"direct_sum": direct_sum, "spectral": spectral, "gauss_legendre": gauss_legendre}
    if solver not in solvers:
        raise ValueError("solver must be direct_sum, spectral, or gauss_legendre")
    defaults = solvers[solver]
    if set(options) - set(defaults):
        raise ValueError(f"unknown {solver} options: {sorted(set(options)-set(defaults))}")
    values = {**defaults, **options}
    for key, value in values.items():
        if not np.isscalar(value) or not np.isfinite(value) or value < 0:
            raise ValueError(f"{key} must be a finite nonnegative scalar")
        if isinstance(defaults[key], int):
            if int(value) != value:
                raise ValueError(f"{key} must be an integer")
            values[key] = int(value)
        if key in ("lmax", "refine_factor", "ntheta_fit", "nphi_fit", "complex_vertex_threshold", "subdivisions") and value < 1:
            raise ValueError(f"{key} must be positive")
        if solver == "gauss_legendre" and value > 256:
            raise ValueError(f"{key} must be at most 256")
        if key in ("auto_mode", "edge_correction") and value not in (0, 1):
            raise ValueError(f"{key} must be 0 or 1")
        if key == "hybrid_mode" and value not in (0, 1, 2):
            raise ValueError("hybrid_mode must be 0, 1, or 2")
    if solver == "spectral" and values["lmax"] + 6*values["auto_mode"] > 46339:
        raise ValueError("spectral degree including the automatic search must be at most 46339")
    return tuple(values.values())


def run_sphere_model(model, *, solver="direct_sum", solver_options=None, output_dir=None, **build_options):
    """evaluate a python model through gravmag sphere and return summed grid arrays

    output_dir optionally retains model.in and field.txt; otherwise cards and
    executable output are temporary implementation details
    field has shape (nlat, nlon, 3), total is the magnitude of the summed vector
    """
    solver = _solver_name(solver)
    content, field_type = _model_cards(model)
    arguments = _solver_arguments(solver, solver_options or {})
    with tempfile.TemporaryDirectory(prefix="gravmagpy-sphere-") as temporary:
        directory = Path(output_dir).expanduser().resolve() if output_dir is not None else Path(temporary)
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "model.in"
        path.write_text(content, encoding="utf-8")
        run = run_grid_model(path, directory / "field.txt", radius_km=model["radius_km"], solver=solver,
                             options=arguments, **build_options)
        table = np.loadtxt(run["output_path"], ndmin=2)
        lat = np.asarray(model["grid"]["latitude_deg"])
        lon = np.asarray(model["grid"]["longitude_deg"])
        count = len(lat)*len(lon)
        blocks = 1 if solver == "spectral" else len(model["sources"])
        if table.shape != (count*blocks, 7) or not np.isfinite(table).all():
            raise RuntimeError("unexpected gravmag sphere output dimensions or nonfinite fields")
        # all solvers emit latitude-major grids; only spectral combines bodies on disk
        field = table[:, 3:6].reshape(blocks, len(lat), len(lon), 3).sum(axis=0)
        unit = "nt" if field_type == "magnetic" else "mgal"
        prefix = "b" if field_type == "magnetic" else "g"
        result = {"field": field, "total": np.linalg.norm(field, axis=-1), "unit": unit,
                  "latitude_deg": lat.copy(), "longitude_deg": lon.copy(), "raw_output": table,
                  "model": deepcopy(model), "solver": solver, "stdout": run["stdout"], "stderr": run["stderr"]}
        result.update({f"{prefix}{axis}_{unit}": field[..., index] for index, axis in enumerate(("x", "y", "z"))})
        result[f"{prefix}tot_{unit}"] = result["total"]
        if solver == "spectral":
            result.update(normalization="schmidt", condon_shortley_phase=False,
                          spectral_algorithm="native_gauss_legendre_fourier",
                          harmonic_recurrence="fukushima_extended_range_radix_2_512")
        if output_dir is not None:
            result.update(input_path=path, output_path=run["output_path"])
    return result
