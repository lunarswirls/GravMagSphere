"""numerical, data-contract, and executable integration checks"""

from copy import deepcopy
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from gravmagpy import (
    block_volume, build_fortran, fit_sources, load_lpmag_csv, predict_field, read_input,
    run_grid_model, save_fit, spherical_to_cartesian, write_source_input,
)


def source():
    return {"lat_deg": 7.5, "lon_deg": -59.0, "lat_width_deg": 1.2, "lon_width_deg": 1.5,
            "depth_top_km": 8.0, "thickness_km": 6.0, "mx_a_m": 1.5, "my_a_m": -0.7, "mz_a_m": 0.9}


class orbital_tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = build_fortran("orbital")

    def test_far_field_matches_independent_point_dipole(self):
        body = dict(source(), lat_deg=0, lon_deg=0, lat_width_deg=0.001, lon_width_deg=0.001, thickness_km=0.01)
        xyz = np.array([[1850, 25, -30], [1850, -25, 30]], dtype=float)
        center = np.array([1737.4 - body["depth_top_km"] - body["thickness_km"] / 2, 0, 0])
        moment = np.array([body[key] for key in ("mx_a_m", "my_a_m", "mz_a_m")]) * block_volume([body])[0]
        delta = (xyz - center) * 1000
        distance = np.linalg.norm(delta, axis=1)
        expected = 100 * (3 * (delta @ moment)[:, None] * delta / distance[:, None]**5 - moment / distance[:, None]**3)
        actual = predict_field(xyz, [body], library=self.library)
        np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=1e-12)

    def test_quadrature_convergence_and_superposition(self):
        xyz = spherical_to_cartesian([7, 8, 9], [-60, -58, -59], [1760, 1770, 1790])
        body = source()
        low = predict_field(xyz, [body], quadrature_order=8, library=self.library)
        high = predict_field(xyz, [body], quadrature_order=12, library=self.library)
        np.testing.assert_allclose(low, high, rtol=1e-4, atol=1e-5)
        reverse = dict(body, **{key: -body[key] for key in ("mx_a_m", "my_a_m", "mz_a_m")})
        cancelled = predict_field(xyz, [body, reverse], quadrature_order=8, library=self.library)
        np.testing.assert_allclose(cancelled, 0, atol=1e-12)

    def test_csv_sel_altitude_wake_and_uncertainties(self):
        xyz = spherical_to_cartesian([7, 8, 9], [-59, -58, -57], [1757.4, 1777.4, 1817.4])
        table = pd.DataFrame(xyz, columns=["X_SEL", "Y_SEL", "Z_SEL"])
        table[["Bx_SEL", "By_SEL", "Bz_SEL"]] = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        table["altitude_km"] = 9999
        table["wake_t_km"] = [200, np.nan, 600]
        table["wake_rperp_km"] = [100, np.nan, 100]
        table["Brms"] = [0.1, 0.1, 100]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample_nightside.csv"
            table.to_csv(path, index=False)
            data = load_lpmag_csv(path, sigma_nt=[0.1, 0.2, 0.3])
            np.testing.assert_allclose(data["table"].altitude_km, [20, 80])
            np.testing.assert_allclose(data["sigma_nt"], [[0.1, 0.2, 0.3]] * 2)
            self.assertEqual(data["table"].source_row.tolist(), [0, 2])
            filtered = load_lpmag_csv(path, max_brms_nt=1)
            self.assertEqual(len(filtered["table"]), 1)
            with self.assertRaises(ValueError):
                load_lpmag_csv([path, path])
            table.drop(columns=["wake_t_km", "wake_rperp_km"]).to_csv(path, index=False)
            with self.assertRaises(ValueError):
                load_lpmag_csv(path)

    def test_spherical_csv_rejects_negative_altitude(self):
        mapping = {"lat_deg": "lat", "lon_deg": "lon", "altitude_km": "alt"}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "selected_wake.csv"
            pd.DataFrame({"lat": [0], "lon": [0], "alt": [-4000],
                          "Bx_SEL": [1], "By_SEL": [2], "Bz_SEL": [3]}).to_csv(path, index=False)
            with self.assertRaises(ValueError):
                load_lpmag_csv(path, columns=mapping)

    def test_multialtitude_geometry_recovery(self):
        lat, lon = np.meshgrid(np.linspace(5.5, 9.5, 9), np.linspace(-61, -57, 9))
        xyz = spherical_to_cartesian(np.tile(lat.ravel(), 3), np.tile(lon.ravel(), 3), np.repeat([1757.4, 1777.4, 1817.4], 81))
        body = source()
        field = predict_field(xyz, [body], quadrature_order=10, library=self.library)
        data = {"xyz_km": xyz, "field_nt": field, "radius_km": 1737.4, "sigma_nt": 0.5}
        initial = [dict(body, lat_deg=7.7, lon_deg=-58.8, depth_top_km=12, lat_width_deg=1.0, lon_width_deg=1.8)]
        original = deepcopy(initial)
        bounds = [{"lat_deg": (6, 9), "lon_deg": (-60.5, -57.5), "depth_top_km": (1, 25),
                   "lat_width_deg": (0.5, 2), "lon_width_deg": (0.5, 2.5)}]
        fitted = fit_sources(data, initial, bounds, quadrature_order=10, library=self.library)
        self.assertTrue(fitted["success"])
        self.assertLess(fitted["rmse_nt"], 1e-7)
        for key in bounds[0]:
            self.assertAlmostEqual(fitted["sources"][0][key], body[key], places=5)
        self.assertEqual(initial, original)
        self.assertEqual(fitted["jacobian_rank"], 5)
        with tempfile.TemporaryDirectory() as directory:
            save_fit(fitted, data, directory)
            self.assertTrue((Path(directory) / "fit.json").is_file())
            self.assertEqual(len(pd.read_csv(Path(directory) / "predictions.csv")), len(xyz))

    def test_magnetization_and_background_recovery(self):
        lat, lon = np.meshgrid(np.linspace(5, 10, 7), np.linspace(-62, -56, 7))
        xyz = spherical_to_cartesian(np.tile(lat.ravel(), 2), np.tile(lon.ravel(), 2), np.repeat([1760, 1800], 49))
        body = source()
        offsets = np.repeat([[1.0, -2.0, 0.5], [-1.0, 0.5, 2.0]], 49, axis=0)
        observed = predict_field(xyz, [body], library=self.library) + offsets
        data = {"xyz_km": xyz, "field_nt": observed, "radius_km": 1737.4, "groups": np.repeat(["a", "b"], 49)}
        initial = [dict(body, mx_a_m=0.8, my_a_m=-0.3, mz_a_m=0.4)]
        bounds = [{key: (-5, 5) for key in ("mx_a_m", "my_a_m", "mz_a_m")}]
        fitted = fit_sources(data, initial, bounds, fit_background=True, library=self.library)
        self.assertTrue(fitted["success"])
        np.testing.assert_allclose(fitted["background_nt"], offsets, atol=1e-7)
        self.assertLess(fitted["rmse_nt"], 1e-7)

    def test_independent_legacy_surface_charge_comparison(self):
        body = source()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.in"
            output = Path(directory) / "direct_sum.txt"
            write_source_input([body], path, grid=(6.5, -60, 0.5, 0.5, 30, 5, 5), mesh=(8, 80, 80))
            parsed = read_input(path)
            self.assertEqual(len(parsed), 1)
            self.assertEqual(parsed[0]["depth_bottom_km"], 14)
            run_grid_model(path, output, options=(1,))
            values = np.loadtxt(output)
            xyz = spherical_to_cartesian(values[:, 2], values[:, 1], 1767.4)
            predicted = predict_field(xyz, [body], quadrature_order=12, library=self.library)
            relative_error = np.linalg.norm(predicted-values[:, 3:6]) / np.linalg.norm(predicted)
            self.assertLess(relative_error, 0.01)

    def test_invalid_geometry_and_fortran_errors(self):
        with self.assertRaises(ValueError):
            predict_field([[1000, 0, 0]], [source()], library=self.library)
        with self.assertRaises(ValueError):
            predict_field([[1800, 0, 0]], [dict(source(), thickness_km=-1)], library=self.library)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "broken.in"
            path.write_text("broken\n")
            with self.assertRaises(ValueError):
                read_input(path)
        self.assertEqual(build_fortran("orbital"), self.library)
