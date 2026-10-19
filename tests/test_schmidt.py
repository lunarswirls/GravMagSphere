"""independent schmidt functions, vector fields and fitted-source checks against shtools"""

import ctypes
import importlib.util
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from gravmagpy import block, sphere_model, run_sphere_model, spherical_to_cartesian


def coefficient_order(degree):
    return [(part, l, m) for l in range(1, degree+1) for m in range(l+1)
            for part in range(1 if m == 0 else 2)]


def spherical_axes(latitude, longitude):
    theta = np.deg2rad(90-np.asarray(latitude))
    phi = np.deg2rad(longitude)
    return np.stack((np.stack((np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)), axis=-1),
                     np.stack((np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)), axis=-1),
                     np.stack((-np.sin(phi), np.cos(phi), np.zeros_like(phi)), axis=-1)), axis=-2)


@unittest.skipUnless(importlib.util.find_spec("pyshtools"), "requires optional pyshtools validation dependency")
class schmidt_tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory(prefix="gravmagpy-schmidt-")
        # the numerical reference does not need the optional gmt plotting backend
        with patch.dict(sys.modules, {"pygmt": None}):
            import pyshtools

        cls.shtools = pyshtools
        cls.records = []
        root = Path(__file__).resolve().parents[1]
        binary = Path(cls.directory.name)/"schmidt_probe.so"
        subprocess.run([os.environ.get("FC", "gfortran"), "-std=f2008", "-O2", "-shared", "-fPIC",
                        f"-J{cls.directory.name}", f"-I{cls.directory.name}",
                        str(root/"fortran/gravmag_sphere_harmonics.f90"),
                        str(root/"tests/fixtures/schmidt_probe.f90"), "-o", str(binary)],
                       check=True, capture_output=True, text=True)
        cls.library = ctypes.CDLL(str(binary))
        cls.probe = cls.library.schmidt_probe
        array = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
        cls.probe.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, array, array, array, array]
        cls.probe.restype = None

    @classmethod
    def tearDownClass(cls):
        destination = os.environ.get("GRAVMAG_SCHMIDT_REPORT")
        if destination:
            path = Path(destination)
            path.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(cls.records).assign(shtools_version=cls.shtools.__version__).to_csv(path, index=False)
        cls.directory.cleanup()

    def evaluate(self, degree, theta, phi=0.0):
        arrays = [np.empty((degree+1, degree+1), order="F") for _ in range(3)]
        basis = np.empty((3, degree*(degree+2)), order="F")
        self.probe(degree, theta, phi, *arrays, basis)
        return *arrays, basis

    def compare(self, name, actual, reference, tolerance, **metadata):
        maximum = float(np.max(np.abs(actual-reference)))
        scaled = maximum/max(1.0, float(np.max(np.abs(reference))))
        self.records.append({"check": name, "max_absolute_difference": maximum,
                             "scaled_max_difference": scaled, "tolerance": tolerance, **metadata})
        self.assertTrue(np.isfinite(actual).all())
        self.assertLess(scaled, tolerance, f"{name}: {metadata}")

    def test_functions_and_derivatives_against_shtools(self):
        for degree in (1, 8, 24, 64, 128, 256):
            indices = [(l, m) for l in range(degree+1) for m in range(l+1)]
            for z in (-0.999999, -0.7, 0.0, 0.2, 0.999999):
                theta = np.arccos(z)
                p, dp, quotient, _ = self.evaluate(degree, theta)
                reference, derivative = self.shtools.legendre.PlmSchmidt_d1(degree, z, csphase=1, cnorm=0)
                selected = tuple(np.asarray(indices).T)
                metadata = {"degree": degree, "cos_colatitude": z}
                self.compare("legendre", p[selected], reference, 5e-10, **metadata)
                self.compare("colatitude_derivative", dp[selected], -np.sin(theta)*derivative, 5e-10, **metadata)
                nonzonal = np.array([m > 0 for _, m in indices])
                self.compare("azimuthal_quotient", quotient[selected][nonzonal],
                             reference[nonzonal]/np.sin(theta), 5e-10, **metadata)

    def test_field_synthesis_against_shtools(self):
        random = np.random.default_rng(281019)
        for degree in (1, 8, 24, 64):
            cilm = np.zeros((2, degree+1, degree+1))
            order = coefficient_order(degree)
            for part, l, m in order:
                cilm[part, l, m] = random.normal()/(l+1)**2
            coefficients = np.array([cilm[item] for item in order])
            for latitude, longitude in ((-89.9, -179), (-65, 73), (0, -35), (23, 147), (89.9, 10)):
                _, _, _, basis = self.evaluate(degree, np.deg2rad(90-latitude), np.deg2rad(longitude))
                reference = self.shtools.gravmag.MakeMagGridPoint(cilm, 1.0, 1.0, latitude, longitude)
                self.compare("vector_synthesis", basis@coefficients, reference, 5e-10,
                             degree=degree, latitude_deg=latitude, longitude_deg=longitude)

    def test_exact_and_near_poles_against_dipole(self):
        # degree-one coefficients correspond to a centered dipole in the cartesian basis
        moment = np.array([0.7, -1.3, 2.1])
        coefficients = moment[[2, 0, 1]]
        for theta in (0.0, 1e-12, 1e-8, np.pi-1e-8, np.pi-1e-12, np.pi):
            for longitude in (-150, 0, 67):
                latitude = 90-np.rad2deg(theta)
                _, _, _, basis = self.evaluate(1, theta, np.deg2rad(longitude))
                axes = spherical_axes(latitude, longitude)
                actual = axes.T@(basis@coefficients)
                reference = 3*np.dot(moment, axes[0])*axes[0]-moment
                self.compare("polar_dipole", actual, reference, 5e-14,
                             colatitude_rad=theta, longitude_deg=longitude)

    def test_public_gravity_and_magnetic_fits_against_shtools(self):
        for degree in (3, 24):
            with self.subTest(degree=degree):
                self.check_source_fit(degree)

    def test_automatic_selection_replays_with_corrections(self):
        for kind in ("gravity", "magnetic"):
            material = {"density_kg_m3": 400} if kind == "gravity" else {"magnetization_a_m": [1, -2, 3]}
            source = block([6, 8], [-60, -58], [2, 8], **material)
            model = sphere_model([source], latitude_deg=[6, 7, 8], longitude_deg=[-60, -59, -58],
                                 altitude_km=5000)
            options = {"lmax": 6, "ntheta_fit": 24, "nphi_fit": 48, "refine_factor": 1,
                       "source_nlat": 2, "source_nlon": 2, "source_nr": 2,
                       "auto_mode": 1, "edge_correction": 1, "hybrid_mode": 2}
            automatic = run_sphere_model(model, solver="spectral", solver_options=options)
            selected = re.search(r'Auto-selected spectral controls: lmax=\s*(\d+)\s+reg_lambda=\s*([\d.E+-]+)'
                                 r'\s+reg_power=\s*([\d.E+-]+)', automatic["stdout"])
            self.assertIsNotNone(selected)
            options.update(lmax=int(selected[1]), reg_lambda=float(selected[2]), reg_power=float(selected[3]),
                           auto_mode=0)
            explicit = run_sphere_model(model, solver="spectral", solver_options=options)
            np.testing.assert_array_equal(automatic["field"], explicit["field"])
            self.records.append({"check": "automatic_selection_replay", "field": kind,
                                 "degree": options["lmax"], "damping": options["reg_lambda"],
                                 "max_absolute_difference": 0.0, "tolerance": 0.0})

    def check_source_fit(self, degree):
        radius = 1737.4
        altitude = 5000.0
        observation_radius = (radius+altitude)*1000
        center = spherical_to_cartesian(7, -59, radius-5)*1000
        volume = ((radius-5)*1000)**2*np.cos(np.deg2rad(7))*np.deg2rad(2)**2*6000
        order = coefficient_order(degree)
        ntheta = max(12, 2*degree+2)
        nphi = max(24, 4*degree+4)
        theta, phi = np.meshgrid((np.arange(ntheta)+0.5)*np.pi/ntheta,
                                -np.pi+np.arange(nphi)*2*np.pi/nphi, indexing="ij")
        latitude = 90-np.rad2deg(theta.ravel())
        longitude = np.rad2deg(phi.ravel())
        axes = spherical_axes(latitude, longitude)
        displacement = observation_radius*axes[:, 0]-center
        distance = np.linalg.norm(displacement, axis=1)
        # use shtools functions and derivatives to assemble the independent fit operator
        design = np.zeros((ntheta, nphi, 3, len(order)))
        for row in range(ntheta):
            angle = theta[row, 0]
            p, derivative = self.shtools.legendre.PlmSchmidt_d1(degree, np.cos(angle), csphase=1, cnorm=0)
            for column, (part, l, m) in enumerate(order):
                packed = l*(l+1)//2+m
                cosine = np.cos(m*phi[row])
                sine = np.sin(m*phi[row])
                trig = cosine if part == 0 else sine
                design[row, :, 0, column] = (l+1)*p[packed]*trig
                design[row, :, 1, column] = np.sin(angle)*derivative[packed]*trig
                design[row, :, 2, column] = m*p[packed]/np.sin(angle)*(sine if part == 0 else -cosine)
        design = design.reshape(-1, len(order))
        weights = np.repeat(np.sin(theta.ravel()), 3)
        weighted = design*np.sqrt(weights[:, None])
        diagonal = np.sum(weighted**2, axis=0)
        degree_weights = np.array([(l/degree)**4 for _, l, _ in order])
        output_latitude, output_longitude = np.meshgrid([-60, 0, 60], [-120, 0, 120], indexing="ij")
        for kind in ("gravity", "magnetic"):
            if kind == "gravity":
                source = block([6, 8], [-60, -58], [2, 8], density_kg_m3=400)
                cartesian = -6.67430e-11*400*volume*displacement/distance[:, None]**3
                monopole = -6.67430e-11*400*volume/observation_radius**2
                unit_scale = 1e5
            else:
                source = block([6, 8], [-60, -58], [2, 8], magnetization_a_m=[1, -2, 3])
                moment = volume*np.array([1, -2, 3])
                cartesian = 1e-7*(3*(displacement@moment)[:, None]*displacement/distance[:, None]**5
                                  - moment/distance[:, None]**3)
                monopole = 0.0
                unit_scale = 1e9
            samples = np.einsum("nij,nj->ni", axes, cartesian)
            samples[:, 0] -= monopole
            model = sphere_model([source], latitude_deg=[-60, 0, 60], longitude_deg=[-120, 0, 120],
                                 altitude_km=altitude, radius_km=radius)
            for damping in (0.0, 0.2):
                penalty = damping*degree_weights*diagonal + np.finfo(float).eps*np.mean(diagonal)
                coefficients = np.linalg.lstsq(np.vstack((weighted, np.diag(np.sqrt(penalty)))),
                                               np.r_[np.sqrt(weights)*samples.ravel(), np.zeros(len(order))],
                                               rcond=None)[0]
                cilm = np.zeros((2, degree+1, degree+1))
                for item, value in zip(order, coefficients):
                    cilm[item] = value
                reference = []
                for lat, lon in zip(output_latitude.ravel(), output_longitude.ravel()):
                    vector = self.shtools.gravmag.MakeMagGridPoint(cilm, 1, 1, lat, lon)
                    vector[0] += monopole
                    reference.append(spherical_axes(lat, lon).T@vector*unit_scale)
                output = Path(self.directory.name)/f"{kind}_{damping}"
                result = run_sphere_model(model, solver="spectral", output_dir=output, solver_options={
                    "lmax": degree, "ntheta_fit": ntheta, "nphi_fit": nphi, "refine_factor": 1,
                    "source_nlat": 1, "source_nlon": 1, "source_nr": 1, "auto_mode": 0,
                    "edge_correction": 0, "hybrid_mode": 0, "reg_lambda": damping})
                reference = np.asarray(reference)
                actual = result["field"].reshape(-1, 3)
                relative = np.linalg.norm(actual-reference)/np.linalg.norm(reference)
                self.records.append({"check": "source_fit", "field": kind, "damping": damping, "degree": degree,
                                     "relative_vector_difference": relative, "tolerance": 3e-8})
                self.assertLess(relative, 3e-8)
                self.assertEqual(result["normalization"], "schmidt")
                self.assertFalse(result["condon_shortley_phase"])
                header = result["output_path"].read_text()
                self.assertIn("normalization=schmidt_semi_normalized condon_shortley_phase=excluded", header)
