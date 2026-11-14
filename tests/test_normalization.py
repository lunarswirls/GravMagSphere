"""Historical-method preservation and geodetic 4pi normalization regressions."""

import json
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np

from gravmagpy import build_fortran, run_sphere_model, write_model_input
from gravmagpy.sphere import _solver_arguments
from native_support import native_library, jacobi_reference
from test_schmidt import coefficient_order, independent_vector


class normalization_tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.baselines = json.loads((Path(__file__).parent/'fixtures/solver_method_baselines.json').read_text())
        cls.library = native_library()

    def test_preserve_original_spectral_volume_and_recent_schmidt_methods(self):
        # These expected fields come from separately compiled historical commits.
        for case in self.baselines['cases']:
            solver = 'spectral' if case['method'] == 'spectral' else 'gauss_legendre'
            options = dict(case['options'])
            if case['method'] == 'schmidt':
                options['normalization'] = 2
            with self.subTest(kind=case['kind'], method=case['method'], case=case['case']):
                result = run_sphere_model(case['model'], solver=solver, solver_options=options)
                np.testing.assert_allclose(result['field'], case['field'], rtol=2e-12, atol=1e-14)
                expected = 'schmidt' if case['method'] == 'schmidt' else 'unnormalized'
                self.assertEqual(result['normalization'], expected)

    def test_full_point_basis_against_scipy_and_multiprecision(self):
        rng = np.random.default_rng(14112026)
        for degree in (1, 8, 24):
            coefficients = rng.normal(size=degree*(degree+2))
            schmidt = coefficients*np.array([np.sqrt(2*l+1) for _, l, _ in coefficient_order(degree)])
            for theta, phi in ((1e-7, .3), (.7, -.9), (2.1, 1.2), (np.pi-1e-7, -.4)):
                actual = np.empty(3)
                self.library.normalized_point_probe(degree, theta, phi, 3, coefficients, actual)
                expected = independent_vector(degree, schmidt, theta, phi)
                # SciPy cos(theta) loses accuracy extremely near a pole; use moderate angles here.
                if 1e-5 < theta < np.pi-1e-5:
                    np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=1e-9)
        degree = 2800
        coefficients = np.zeros(degree*(degree+2))
        modes = [(degree, 0, .4, 0.), (degree, 1, -.7, .3),
                 (degree, 1204, 1., -.6), (degree, degree, .2, .9)]
        for l, m, g, h in modes:
            index = l*l+max(0, 2*m-1)-1
            coefficients[index] = g
            if m:
                coefficients[index+1] = h
        for theta, phi in ((0., .3), (1e-7, -.4), (np.pi/6, .9), (np.pi/2, 1.9), (np.pi, -.7)):
            expected = np.zeros(3)
            for l, m, g, h in modes:
                p, derivative, quotient = jacobi_reference(l, m, max(theta, 1e-100))
                cosine, sine = np.cos(m*phi), np.sin(m*phi)
                amplitude = g*cosine+h*sine
                expected += np.sqrt(2*l+1)*np.array([(l+1)*p*amplitude, -derivative*amplitude,
                                                    m*quotient*(g*sine-h*cosine)])
            actual = np.empty(3)
            self.library.normalized_point_probe(degree, theta, phi, 3, coefficients, actual)
            np.testing.assert_allclose(actual, expected, rtol=2e-9, atol=2e-8)

    def test_full_bandwidth_coefficient_roundtrip(self):
        rng = np.random.default_rng(15112026)
        for degree in (8, 64, 128):
            coefficients = rng.normal(size=degree*(degree+2))
            for joint in (0., .4, 1.):
                recovered = np.empty_like(coefficients)
                statistics = np.zeros(3)
                self.library.normalized_roundtrip_probe(degree, degree+1, 2*degree+2, joint, 3,
                                                       coefficients, recovered, statistics)
                self.assertLess(statistics[1], 2e-10)
                self.assertLess(statistics[2], 2e-10)

    def test_full_and_schmidt_fields_with_damping_corrections_and_collective_bodies(self):
        for case in self.baselines['cases']:
            if case['method'] != 'schmidt':
                continue
            with self.subTest(kind=case['kind'], case=case['case']), tempfile.TemporaryDirectory() as directory:
                result = run_sphere_model(case['model'], solver='gauss_legendre', output_dir=directory,
                                         solver_options=dict(case['options'], normalization=3))
                np.testing.assert_allclose(result['field'], case['field'], rtol=2e-11, atol=1e-13)
                self.assertEqual(result['normalization'], 'fully_normalized')
                self.assertFalse(result['condon_shortley_phase'])
                self.assertEqual(result['raw_output'].shape, (9, 7))
                header = result['output_path'].read_text()
                self.assertIn('normalization=fully_normalized_4pi', header)
                self.assertIn('# solver=gauss_legendre', header)

    def test_option_names_validation_and_executable_dispatch(self):
        for number, name in ((1, 'unnormalized'), (2, 'schmidt'), (3, 'fully_normalized'), (3, '4pi')):
            self.assertEqual(_solver_arguments('gauss_legendre', {'normalization': number}),
                             _solver_arguments('gauss_legendre', {'normalization': name}))
        for value in (0, 4, 2.5, np.nan, np.inf, None, [], 'unknown'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                _solver_arguments('gauss_legendre', {'normalization': value})
        for options in ({'normalization': 2, 'radial_order': 4}, {'normalization': 1, 'lmax': 4},
                        {'normalization': 3, 'lmax': 46334}, {'normalization': 3, 'edge_correction': 2}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                _solver_arguments('gauss_legendre', options)
        # The harmonic path must not inherit the volume rule's 256-order limit.
        self.assertEqual(_solver_arguments('gauss_legendre', {'normalization': 3, 'lmax': 600})[5], 600)
        case = self.baselines['cases'][0]
        with tempfile.TemporaryDirectory() as directory:
            path = write_model_input(case['model'], Path(directory)/'case.in')
            binary = str(build_fortran('gauss_legendre'))
            for value in ('0', '4', 'nan', 'unknown'):
                run = subprocess.run([binary, '1737.4', str(path), str(Path(directory)/'bad.txt'),
                                      '0', '0', '0', '1', value], capture_output=True)
                self.assertNotEqual(run.returncode, 0)
