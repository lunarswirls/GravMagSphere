"""quadrature accuracy, geometry, units and three-solver convergence diagnostics"""

from copy import deepcopy
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from gravmagpy import (block, polygon, sphere_model, run_sphere_model, predict_field,
                       spherical_to_cartesian, write_model_input, run_grid_model, build_fortran)


def make_model(sources, altitude=30, radius=1737.4):
    return sphere_model(sources, latitude_deg=[6, 7, 8], longitude_deg=[-60, -59, -58],
                        altitude_km=altitude, radius_km=radius)


def quadrature(model, order=12, subdivisions=1):
    return run_sphere_model(model, solver="gauss_legendre", solver_options={
        "radial_order": order, "latitude_order": order, "longitude_order": order,
        "subdivisions": subdivisions})


def relative_difference(field, reference):
    return np.linalg.norm(field-reference)/np.linalg.norm(reference)


class gauss_legendre_tests(unittest.TestCase):
    def test_original_tensor_rule_against_orbital_backend(self):
        source = block([6, 8], [-60, -58], [2, 8], magnetization_a_m=[1, -2, 3])
        orbital = {"lat_deg": 7, "lon_deg": -59, "lat_width_deg": 2, "lon_width_deg": 2,
                   "depth_top_km": 2, "thickness_km": 6, "mx_a_m": 1, "my_a_m": -2, "mz_a_m": 3}
        lat, lon = np.meshgrid([6, 7, 8], [-60, -59, -58], indexing="ij")
        for radius in (1737.4, 6371.0):
            for order in (1, 3, 8):
                with self.subTest(radius=radius, order=order):
                    model = make_model([source], radius=radius)
                    reference = predict_field(spherical_to_cartesian(lat.ravel(), lon.ravel(), radius+30),
                                              [orbital], radius_km=radius, quadrature_order=order)
                    np.testing.assert_allclose(quadrature(model, order)["field"].reshape(-1, 3), reference,
                                               rtol=2e-10, atol=2e-9)

    def test_node_order_and_subdivision_convergence(self):
        source = block([6, 8], [-60, -58], [2, 8], magnetization_a_m=[1, -2, 3])
        model = make_model([source], altitude=5)
        reference = quadrature(model, 40, 2)["field"]
        coarse = relative_difference(quadrature(model, 4)["field"], reference)
        medium = relative_difference(quadrature(model, 12)["field"], reference)
        fine = relative_difference(quadrature(model, 32)["field"], reference)
        split = relative_difference(quadrature(model, 4, 4)["field"], reference)
        self.assertLess(medium, coarse/10)
        self.assertLess(fine, medium/100)
        self.assertLess(fine, 1e-5)
        self.assertLess(split, coarse/100)

    def test_concave_polygon_orientation_closure_and_superposition(self):
        vertices = [[6, -60], [6, -58], [7, -58], [7, -59], [8, -59], [8, -60]]
        material = {"magnetization_a_m": [1, -2, 3]}
        pieces = [block([6, 7], [-60, -58], [2, 8], **material),
                  block([7, 8], [-60, -59], [2, 8], **material)]
        reference = quadrature(make_model(pieces), 24)
        self.assertEqual(reference["raw_output"].shape, (18, 7))
        for outline in (vertices, vertices[::-1], vertices + vertices[:1]):
            result = quadrature(make_model([polygon(outline, [2, 8], **material)]), 24)
            np.testing.assert_allclose(result["field"], reference["field"], rtol=2e-9, atol=2e-8)
            np.testing.assert_allclose(result["total"], np.linalg.norm(result["field"], axis=-1))
        # a meridian through the notch intersects two disjoint latitude intervals
        notched = [[6, -60], [6, -58], [6.5, -58], [6.5, -59],
                   [7.5, -59], [7.5, -58], [8, -58], [8, -60]]
        notched_pieces = [block([6, 8], [-60, -59], [2, 8], **material),
                          block([6, 6.5], [-59, -58], [2, 8], **material),
                          block([7.5, 8], [-59, -58], [2, 8], **material)]
        np.testing.assert_allclose(quadrature(make_model([polygon(notched, [2, 8], **material)]), 24)["field"],
                                   quadrature(make_model(notched_pieces), 24)["field"], rtol=2e-9, atol=2e-8)
        cancellation = quadrature(make_model([pieces[0], dict(pieces[0], magnetization_a_m=[-1, 2, -3])]))
        np.testing.assert_allclose(cancellation["field"], 0, atol=1e-9)

    def test_sloping_polygon_volume_and_longitude_seam(self):
        vertices = np.array([[6, 179.3], [8, 179.8], [7, 180.7]])
        wrapped = vertices.copy()
        wrapped[wrapped[:, 1] > 180, 1] -= 360
        model = sphere_model([polygon(wrapped, [2, 8], density_kg_m3=400)],
                             latitude_deg=[6, 7, 8], longitude_deg=[179.5, 180, 180.5])
        with patch.dict(os.environ, {"GRAVMAG_DIAGNOSTICS": "1"}):
            result = quadrature(model, 24)
        lat = np.deg2rad(vertices[:, 0])
        lon = np.deg2rad(vertices[:, 1])
        dlat = np.roll(lat, -1)-lat
        dlon = np.roll(lon, -1)-lon
        area = abs(np.sum(dlon*np.sin(lat+dlat/2)*np.sinc(dlat/(2*np.pi))))
        volume = area*((1737.4-2)**3-(1737.4-8)**3)*1e9/3
        reported = float(next(line.split("=")[1] for line in result["stdout"].splitlines() if "volume_m3=" in line))
        self.assertAlmostEqual(reported/volume, 1, places=10)
        unwrapped_model = deepcopy(model)
        unwrapped_model["sources"][0]["vertices_lat_lon"] = vertices.tolist()
        np.testing.assert_allclose(quadrature(unwrapped_model, 24)["field"], result["field"], rtol=1e-11, atol=1e-11)
        self.assertTrue(np.all(result["raw_output"][:, 1] < 180))

    def test_gravity_point_mass_sign_units_and_negative_density(self):
        source = block([6.9999, 7.0001], [-59.0001, -58.9999], [4.999, 5.001], density_kg_m3=400)
        model = make_model([source], altitude=100)
        result = quadrature(model)
        lat, lon = np.meshgrid([6, 7, 8], [-60, -59, -58], indexing="ij")
        observation = spherical_to_cartesian(lat, lon, 1837.4)*1000
        center = spherical_to_cartesian(7, -59, 1732.4)*1000
        delta = center-observation
        volume = (np.deg2rad(0.0002)*(np.sin(np.deg2rad(7.0001))-np.sin(np.deg2rad(6.9999))) *
                  ((1737.4-4.999)**3-(1737.4-5.001)**3)*1e9/3)
        expected = 6.67430e-6*400*volume*delta/np.linalg.norm(delta, axis=-1)[..., None]**3
        self.assertEqual(result["unit"], "mgal")
        np.testing.assert_allclose(result["field"], expected, rtol=2e-8, atol=1e-14)
        reverse = quadrature(make_model([dict(source, density_kg_m3=-400)], altitude=100))
        np.testing.assert_allclose(reverse["field"], -result["field"], rtol=1e-12)
        with tempfile.TemporaryDirectory() as directory:
            path = write_model_input(model, Path(directory)/"case.in")
            output = run_grid_model(path, solver="gauss_legendre")["output_path"]
            converted = Path(directory)/"spherical.txt"
            subprocess.run([str(build_fortran("xyz_to_brtp")), str(output), str(converted)],
                           check=True, capture_output=True)
            self.assertIn("Br_mGal", converted.read_text())
            self.assertIn("# solver=gauss_legendre", converted.read_text())

    def test_direct_sum_and_pure_spectral_far_field_convergence(self):
        source = block([6, 8], [-60, -58], [2, 8], magnetization_a_m=[1, -2, 3])
        model = make_model([source], altitude=5000)
        reference = quadrature(model)["field"]
        direct_sum = run_sphere_model(model, solver="direct_sum", solver_options={"refine_factor": 2})["field"]
        self.assertLess(relative_difference(direct_sum, reference), 0.005)
        residuals = []
        for degree in (4, 8, 12):
            result = run_sphere_model(model, solver="spectral", solver_options={
                "lmax": degree, "ntheta_fit": max(24, 2*degree+2), "nphi_fit": max(48, 4*degree+4),
                "auto_mode": 0, "edge_correction": 0, "hybrid_mode": 0, "reg_lambda": 0})
            residuals.append(relative_difference(result["field"], reference))
        self.assertLess(residuals[1], residuals[0]/10)
        self.assertLess(residuals[2], residuals[1]/10)
        self.assertLess(residuals[2], 1e-4)

    def test_direct_sum_gravity_and_polygon_orbital_agreement(self):
        for material in ({"magnetization_a_m": [1, -2, 3]}, {"density_kg_m3": 400}):
            source = polygon([[6, -60], [8, -59.5], [7, -58]], [2, 8], **material)
            model = make_model([source])
            reference = quadrature(model, 24)["field"]
            direct_sum = run_sphere_model(model, solver="direct_sum", solver_options={
                "source_nlat": 96, "source_nlon": 96, "source_nr": 12})["field"]
            self.assertLess(relative_difference(direct_sum, reference), 0.01)

    def test_direct_sum_concave_polygon_winding_closure_and_convergence(self):
        vertices = [[6, -60], [6, -58], [6.5, -58], [6.5, -59],
                    [7.5, -59], [7.5, -58], [8, -58], [8, -60]]
        material = {"magnetization_a_m": [1, -2, 3]}
        pieces = [block([6, 8], [-60, -59], [2, 8], **material),
                  block([6, 6.5], [-59, -58], [2, 8], **material),
                  block([7.5, 8], [-59, -58], [2, 8], **material)]
        reference = quadrature(make_model(pieces), 24)["field"]
        for outline in (vertices, vertices[::-1], vertices + vertices[:1], vertices[3:] + vertices[:3]):
            with self.subTest(outline=outline):
                model = make_model([polygon(outline, [2, 8], **material)])
                residuals = []
                for refine in (1, 4):
                    result = run_sphere_model(model, solver="direct_sum", solver_options={
                        "refine_factor": refine, "source_nr": 12})
                    residuals.append(relative_difference(result["field"], reference))
                self.assertLess(residuals[1], 0.002)
                self.assertLess(residuals[1], residuals[0]/2)

    def test_spectral_gravity_far_field_convergence(self):
        source = block([6, 8], [-60, -58], [2, 8], density_kg_m3=400)
        model = make_model([source], altitude=5000)
        reference = quadrature(model, 16)["field"]
        residuals = []
        for degree in (4, 8, 12):
            result = run_sphere_model(model, solver="spectral", solver_options={
                "lmax": degree, "ntheta_fit": max(24, 2*degree+2), "nphi_fit": max(48, 4*degree+4),
                "auto_mode": 0, "edge_correction": 0, "hybrid_mode": 0, "reg_lambda": 0,
                "refine_factor": 4, "source_nr": 12})
            self.assertEqual(result["unit"], "mgal")
            residuals.append(relative_difference(result["field"], reference))
        self.assertLess(residuals[1], residuals[0]/10)
        self.assertLess(residuals[2], residuals[1]/10)
        self.assertLess(residuals[2], 1e-5)

    def test_spectral_gravity_monopole_sign_mass_radius_and_regularization(self):
        # antipodal sources have no degree-one field, leaving only the mass term at lmax=1
        volume = (np.deg2rad(2)*(np.sin(np.deg2rad(8))-np.sin(np.deg2rad(6))) *
                  ((1737.4-2)**3-(1737.4-8)**3)*1e9/3)
        lat, lon = np.meshgrid([6, 7, 8], [-60, -59, -58], indexing="ij")
        radial = spherical_to_cartesian(lat, lon, 1)
        for density in (400, -400):
            sources = [block([6, 8], [-60, -58], [2, 8], density_kg_m3=density),
                       block([-8, -6], [120, 122], [2, 8], density_kg_m3=density)]
            for altitude in (5000, 10000):
                expected = -6.67430e-6*(2*density*volume)/((1737.4+altitude)*1000)**2*radial
                for regularization in (0, 1000):
                    with self.subTest(density=density, altitude=altitude, regularization=regularization):
                        result = run_sphere_model(make_model(sources, altitude=altitude), solver="spectral",
                                                  solver_options={"lmax": 1, "ntheta_fit": 24, "nphi_fit": 48,
                                                  "auto_mode": 0, "edge_correction": 0, "hybrid_mode": 0,
                                                  "reg_lambda": regularization, "source_nr": 12})
                        np.testing.assert_allclose(result["field"], expected, rtol=2e-6, atol=1e-12)

    def test_cards_defaults_overrides_paths_and_invalid_controls(self):
        source = block([6, 8], [-60, -58], [2, 8], magnetization_a_m=[1, -2, 3],
                       mesh={"radial": 3, "latitude": 3, "longitude": 3})
        model = make_model([source])
        with tempfile.TemporaryDirectory() as directory:
            path = write_model_input(model, Path(directory)/"case with spaces.in")
            output = run_grid_model(path, solver="gauss_legendre")["output_path"]
            self.assertEqual(output, path.parent/"output/case with spaces_quadrature.txt")
            np.testing.assert_allclose(np.loadtxt(output), quadrature(model, 3)["raw_output"], rtol=1e-12)
            for options in ({"subdivisions": 0}, {"latitude_order": 1.5}, {"radial_order": -1},
                            {"longitude_order": 257}, {"subdivisions": np.nan}, {"refine_factor": 2}):
                with self.assertRaises(ValueError):
                    run_sphere_model(model, solver="gauss_legendre", solver_options=options)
            # the executable rejects invalid controls independently of python validation
            run = subprocess.run([str(build_fortran("gauss_legendre")), "1737.4", str(path),
                                  str(Path(directory)/"bad.txt"), "0", "0", "0", "0"], capture_output=True)
            self.assertNotEqual(run.returncode, 0)
        with self.assertRaisesRegex(RuntimeError, "above the source top"):
            quadrature(make_model([dict(source, depth_km=[0, 8])], altitude=0))
