"""modern gravmag sphere model and equivalent-source integration checks"""

from copy import deepcopy
from pathlib import Path
import tempfile
import unittest
import os
from unittest.mock import patch

import numpy as np

from gravmagpy import (block, polygon, equivalent_source_grid, sphere_model, run_sphere_model,
                       write_model_input, read_model_input, read_input, run_grid_model,
                       fit_equivalent_sources, spherical_to_cartesian, build_fortran)
from gravmagpy.utils.fortran import fortran_source_dir


def small_model(sources):
    return sphere_model(sources, latitude_deg=[6, 7, 8], longitude_deg=[-61, -60, -59, -58], altitude_km=40)


class sphere_tests(unittest.TestCase):
    def test_fortran_sources_are_peers(self):
        directory = fortran_source_dir()
        for solver in ("point_source", "spectral", "gauss_legendre"):
            self.assertTrue((directory / f"gravmag_sphere_{solver}.f90").is_file())
        self.assertTrue((directory / "gravmag_sphere_dipole_grid_fit.f90").is_file())
        self.assertFalse((directory / "legacy").exists())

    def test_legacy_solver_selectors_preserve_collective_results(self):
        source = block([6, 7], [-60, -59], [2, 8], magnetization_a_m=[1, 0, 0])
        model = small_model([source, source])
        for solver, alias, options in (
                ("point_source", "direct", {"refine_factor": 1}),
                ("spectral", "gauss", {"lmax": 3, "ntheta_fit": 8, "nphi_fit": 16,
                                      "auto_mode": 0, "edge_correction": 0, "hybrid_mode": 0})):
            with self.subTest(solver=solver), patch.dict(os.environ, {"GRAVMAG_DIAGNOSTICS": "1"}):
                self.assertEqual(build_fortran(solver), build_fortran(alias))
                expected = run_sphere_model(model, solver=solver, solver_options=options)
                actual = run_sphere_model(model, solver=alias, solver_options=options)
                np.testing.assert_array_equal(actual["field"], expected["field"])
                self.assertIn(f"DIAG|{solver}|", actual["stdout"])

    def test_blocks_and_polygons_round_trip(self):
        sources = [block([6, 7], [-60, -59], [2, 8], magnetization_a_m=[1, -2, 3]),
                   polygon([[7, -61], [8, -60], [7, -59]], [3, 9], magnetization_a_m=[-1, 2, 1])]
        sources[0].update(card4=[1, 2, 3, 4, 5], card6=[0, 1, 1, 10, 20])
        model = small_model(sources)
        original = deepcopy(model)
        with tempfile.TemporaryDirectory() as directory:
            path = write_model_input(model, Path(directory) / "sources.in")
            restored = read_model_input(path)
            self.assertEqual(restored["grid"], model["grid"])
            self.assertEqual(restored["sources"][0]["card4"], sources[0]["card4"])
            self.assertEqual(restored["sources"][0]["card6"], sources[0]["card6"])
            for expected, actual in zip(sources, restored["sources"]):
                self.assertEqual(expected["geometry"], actual["geometry"])
                np.testing.assert_allclose(actual["magnetization_a_m"], expected["magnetization_a_m"], atol=1e-10)
            run = run_sphere_model(restored, solver_options={"refine_factor": 1})
            self.assertEqual(run["field"].shape, (3, 4, 3))
            self.assertTrue(np.isfinite(run["field"]).all())
        self.assertEqual(model, original)

    def test_existing_input_card_numerical_compatibility(self):
        root = Path(__file__).resolve().parents[1]
        path = root / "examples/lunar_examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in"
        model = read_model_input(path)
        result = run_sphere_model(model, solver_options={"refine_factor": 1})
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "original.txt"
            run_grid_model(path, output, options=(1,))
            np.testing.assert_allclose(result["raw_output"], np.loadtxt(output), rtol=1e-6, atol=1e-8)

    def test_equivalent_volume_grid_and_vector_sum(self):
        magnetization = np.zeros((2, 2, 3, 3))
        magnetization[..., 0] = np.arange(12).reshape(2, 2, 3)
        cells = equivalent_source_grid([5, 6, 7], [-62, -61, -60, -59], [1, 3, 5], magnetization_a_m=magnetization)
        self.assertEqual(len(cells), 12)
        self.assertEqual(cells[-1]["depth_km"], [3, 5])
        self.assertEqual(cells[-1]["magnetization_a_m"], [11, 0, 0])
        source = block([6, 7], [-60, -59], [2, 8], magnetization_a_m=[1, 0, 0])
        reversed_source = dict(source, magnetization_a_m=[-1, 0, 0])
        run = run_sphere_model(small_model([source, reversed_source]), solver_options={"refine_factor": 1})
        np.testing.assert_allclose(run["field"], 0, atol=1e-6)
        np.testing.assert_allclose(run["btot_nt"], np.linalg.norm(run["field"], axis=-1))

    def test_gravity_and_named_spectral_options(self):
        gravity = block([6, 7], [-60, -59], [2, 8], density_kg_m3=400)
        run = run_sphere_model(small_model([gravity]), solver_options={"refine_factor": 1})
        self.assertEqual(run["unit"], "mgal")
        self.assertGreater(np.max(run["gtot_mgal"]), 0)
        magnetic = block([6, 7], [-60, -59], [2, 8], magnetization_a_m=[1, 0, 0],
                         mesh={"radial": 2, "latitude": 4, "longitude": 4})
        run = run_sphere_model(small_model([magnetic]), solver="spectral", solver_options={
            "lmax": 3, "ntheta_fit": 8, "nphi_fit": 16, "refine_factor": 1,
            "auto_mode": 0, "edge_correction": 0, "hybrid_mode": 0})
        self.assertEqual(run["bx_nt"].shape, (3, 4))
        self.assertTrue(np.isfinite(run["field"]).all())

    def test_invalid_modern_models(self):
        magnetic = block([6, 7], [-60, -59], [2, 8], magnetization_a_m=[1, 0, 0])
        with self.assertRaises(ValueError):
            block([6, 7], [-60, -59], [2, 8])
        with self.assertRaises(ValueError):
            sphere_model([magnetic], latitude_deg=[1, 2, 4], longitude_deg=[1, 2])
        with self.assertRaises(ValueError):
            small_model([magnetic, block([6, 7], [-60, -59], [2, 8], density_kg_m3=1)])
        with self.assertRaises(ValueError):
            small_model([dict(magnetic, depth_km=[2, 2000])])
        with self.assertRaises(ValueError):
            run_sphere_model(small_model([magnetic]), solver_options={"unknown": 1})
        with self.assertRaises(ValueError):
            small_model([dict(magnetic, name="bad\ncard")])

    def test_equivalent_dipole_depth_and_prediction(self):
        lat, lon = np.meshgrid(np.linspace(6.1, 7.9, 4), np.linspace(-60.9, -59.1, 4), indexing="ij")
        xyz = spherical_to_cartesian(np.tile(lat.ravel(), 2), np.tile(lon.ravel(), 2), np.repeat([1767.4, 1797.4], 16))
        source_xyz = spherical_to_cartesian(7, -60, 1727.4)
        moment = np.array([1e12, -2e12, 3e12])
        delta = (xyz-source_xyz)*1000
        distance = np.linalg.norm(delta, axis=1)
        observed = 100*(3*(delta @ moment)[:, None]*delta/distance[:, None]**5 - moment/distance[:, None]**3)
        observations = {"xyz_km": xyz, "field_nt": observed, "radius_km": 1737.4, "sigma_nt": 1.0}
        with tempfile.TemporaryDirectory() as directory:
            result = fit_equivalent_sources(observations, depth_km=10, spacing_deg=(1, 1),
                                            regularization=0.0, output_dir=directory)
            self.assertTrue((Path(directory) / "dipoles.csv").is_file())
            self.assertEqual(len(result["predictions"]), len(xyz))
            np.testing.assert_allclose(result["source_table"]["depth_km"], 10, atol=1e-10)
            self.assertAlmostEqual(result["depth_below_min_observation_km"], 40)
            self.assertLess(result["rmse_nt"], np.sqrt(np.mean(observed**2))*0.1)
            # independently reconstruct the fitted moments to verify units and row ordering
            sources = result["source_table"]
            positions = spherical_to_cartesian(sources.lat_deg, sources.lon_deg, sources.radius_m/1000)
            prediction = np.zeros_like(observed)
            for position, fitted_moment in zip(positions, result["moments_a_m2"]):
                delta = (xyz-position)*1000
                distance = np.linalg.norm(delta, axis=1)
                prediction += 100*(3*(delta @ fitted_moment)[:, None]*delta/distance[:, None]**5 - fitted_moment/distance[:, None]**3)
            np.testing.assert_allclose(result["predicted_nt"], prediction, rtol=1e-10, atol=1e-10)
        with self.assertRaises(ValueError):
            fit_equivalent_sources(dict(observations, sigma_nt=[1, 2, 3]))
        with self.assertRaises(ValueError):
            fit_equivalent_sources(observations, spacing_deg=(0.001, 0.001), max_memory_mib=1)
