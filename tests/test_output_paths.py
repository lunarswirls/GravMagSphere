"""input-relative numeric and figure destinations across python, bash, and fortran"""

import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

from gravmagpy import block, run_grid_model, sphere_model, write_model_input
from gravmagpy.utils import artifact_path, build_fortran, input_directory


class output_path_tests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="gravmagpy-paths-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.case_dir = self.root / "case with spaces ' $literal"
        self.case_dir.mkdir()
        self.cwd = self.root / "unrelated working directory"
        self.cwd.mkdir()
        self.repo = Path(__file__).resolve().parents[1]
        self.model = sphere_model([block([6, 7], [-60, -59], [2, 8], magnetization_a_m=[1, 2, 3])],
                                  latitude_deg=[6, 7], longitude_deg=[-60, -59], altitude_km=40)
        self.input_path = write_model_input(self.model, self.case_dir / "small.in")

    def run_command(self, command):
        environment = dict(os.environ, MPLCONFIGDIR=os.environ.get("MPLCONFIGDIR", str(self.root / "matplotlib")),
                           PYTHONPATH=str(self.repo / "src"))
        run = subprocess.run(list(map(str, command)), cwd=self.cwd, env=environment,
                             capture_output=True, text=True, timeout=60)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertNotIn("ERROR:", run.stdout + run.stderr)
        return run

    def make_csv(self):
        directory = self.root / "csv inputs"
        directory.mkdir(exist_ok=True)
        path = directory / "observations.csv"
        pd.DataFrame({"latitude_deg": [6, 7], "longitude_deg": [-60, -59], "altitude_km": [30, 40],
                      "bx_nt": [1, 2], "by_nt": [2, 3], "bz_nt": [3, 4]}).to_csv(path, index=False)
        return path

    def make_launchers(self):
        directory = self.root / "launchers"
        directory.mkdir()
        scripts = self.repo / "examples"
        for name in ("run_gravmag_sphere_point_source.sh", "run_gravmag_sphere_spectral.sh", "run_gravmag_sphere_gauss_legendre.sh",
                     "run_input_to_xyz.sh", "run_xyz_to_brtp.sh", "run_gravmag_end_to_end.sh",
                     "run_all_examples_end_to_end.sh", "run_all_examples_brtp.py"):
            shutil.copy2(scripts / name, directory / name)
        shutil.copy2(self.repo / "tests/fixtures/skip_build.sh", directory / "build_gravmag_tools.sh")
        (directory / "build_gravmag_tools.sh").chmod(0o755)
        for name, target in (("gravmag_sphere_point_source", "point_source"), ("gravmag_sphere_spectral", "spectral"),
                             ("gravmag_sphere_gauss_legendre", "gauss_legendre"),
                             ("gravmag_xyz_to_brtp", "xyz_to_brtp")):
            (directory / name).symlink_to(build_fortran(target))
        dipole_dir = directory / "dipole_fit_test"
        dipole_dir.mkdir()
        shutil.copy2(scripts / "dipole_fit_test/run_dipole_grid_fit.sh", dipole_dir)
        (dipole_dir / "gravmag_sphere_dipole_grid_fit").symlink_to(build_fortran("dipole_grid"))
        return directory

    def test_artifact_paths_and_derived_inputs(self):
        self.assertEqual(artifact_path(self.input_path), self.case_dir / "output/small.txt")
        self.assertEqual(artifact_path(self.input_path, kind="figs", suffix=".png"), self.case_dir / "figs/small.png")
        self.assertEqual(input_directory(self.case_dir / "output/run/dipoles.csv"), self.case_dir)
        self.assertEqual(input_directory(self.case_dir / "figs/old.png"), self.case_dir)
        self.assertFalse((self.case_dir / "output").exists())
        with self.assertRaises(ValueError):
            artifact_path(self.input_path, name="../escape.txt")

    def test_fortran_point_source_spectral_conversion_and_explicit_paths(self):
        point_source = build_fortran("point_source")
        self.run_command([point_source, 1737.4, self.input_path])
        output = self.case_dir / "output/small.txt"
        self.assertTrue(output.is_file())
        explicit = self.root / "custom numeric" / "chosen.txt"
        self.run_command([point_source, 1737.4, self.input_path, explicit])
        np.testing.assert_allclose(np.loadtxt(output), np.loadtxt(explicit))
        self.run_command([build_fortran("xyz_to_brtp"), output])
        self.assertTrue((self.case_dir / "output/small_brtp.txt").is_file())
        spectral = build_fortran("spectral")
        self.run_command([spectral, 1737.4, self.input_path, "", 3, 1, 8, 16, 0.2, 4, 0, 0, 0, 0, 0, 0, 0])
        self.assertTrue((self.case_dir / "output/small_spectral.txt").is_file())
        self.assertFalse((self.case_dir / "output/output").exists())
        self.assertFalse((self.cwd / "output").exists())

    def test_fortran_csv_defaults_use_first_input(self):
        path = self.make_csv()
        other = self.root / "second.csv"
        shutil.copy2(path, other)
        executable = build_fortran("dipole_grid")
        self.run_command([executable, 1737.4, f"{path},{other}"])
        for suffix in ("_dipole_fit_predictions.csv", "_dipole_fit_dipoles.csv"):
            self.assertTrue(artifact_path(path, suffix=suffix).is_file())
        self.assertFalse((self.root / "output").exists())

    def test_quadrature_fortran_and_shell_launchers(self):
        scripts = self.make_launchers()
        self.run_command([build_fortran("gauss_legendre"), 1737.4, self.input_path])
        expected = np.loadtxt(self.case_dir / "output/small_quadrature.txt")
        explicit = self.root / "quadrature output" / "chosen.txt"
        self.run_command(["bash", scripts / "run_gravmag_sphere_gauss_legendre.sh", 1737.4, self.input_path, explicit])
        np.testing.assert_allclose(np.loadtxt(explicit), expected)
        self.run_command(["bash", scripts / "run_input_to_xyz.sh", "gauss_legendre", 1737.4, self.input_path,
                          explicit, "--refine-factor", 1])
        np.testing.assert_allclose(np.loadtxt(explicit), expected)
        self.run_command(["bash", scripts / "run_all_examples_end_to_end.sh", "gauss_legendre", 1737.4, self.case_dir])
        self.assertTrue((self.case_dir / "output/small_brtp.txt").is_file())

    def test_shell_single_case_launchers_and_overrides(self):
        scripts = self.make_launchers()
        output = self.case_dir / "output"
        for name in ("run_gravmag_sphere_point_source.sh", "run_gravmag_end_to_end.sh"):
            self.run_command(["bash", scripts / name, 1737.4, self.input_path])
        self.run_command(["bash", scripts / "run_input_to_xyz.sh", "point_source", 1737.4, self.input_path])
        self.run_command(["bash", scripts / "run_gravmag_sphere_spectral.sh", 1737.4, self.input_path,
                          "", 3, 1, 8, 16, 0.2, 4, 0, 0, 0, 0, 0, 0, 0])
        self.run_command(["bash", scripts / "run_xyz_to_brtp.sh", output / "small_xyz.txt"])
        for name in ("small_xyz.txt", "small_brtp.txt", "small_spectral.txt"):
            self.assertTrue((output / name).is_file())
        explicit = self.root / "explicit" / "chosen.txt"
        self.run_command(["bash", scripts / "run_input_to_xyz.sh", "point_source", 1737.4, self.input_path, explicit])
        self.assertTrue(explicit.is_file())
        path = self.make_csv()
        self.run_command(["bash", scripts / "dipole_fit_test/run_dipole_grid_fit.sh", 1737.4, path])
        self.assertTrue(artifact_path(path, suffix="_dipole_fit_dipoles.csv").is_file())
        # either optional csv destination can be supplied independently
        self.run_command(["bash", scripts / "dipole_fit_test/run_dipole_grid_fit.sh", 1737.4, path, "", explicit])
        self.assertFalse((scripts / "output").exists())
        nested = output / "nested run"
        nested.mkdir()
        shutil.copy2(output / "small_xyz.txt", nested / "nested_xyz.txt")
        self.run_command(["bash", scripts / "run_xyz_to_brtp.sh", nested / "nested_xyz.txt"])
        self.assertTrue((output / "nested_brtp.txt").is_file())
        self.assertFalse((nested / "output").exists())

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib is optional")
    def test_python_grid_and_figure_defaults(self):
        from gravmagpy.plotting import plot_fields
        output = run_grid_model(self.input_path, options=(1,))["output_path"]
        self.assertEqual(output, self.case_dir / "output/small.txt")
        figure = plot_fields(self.input_path, output)
        self.assertEqual(figure, self.case_dir / "figs/small_bxyz.png")
        explicit = self.root / "chosen.png"
        self.assertEqual(plot_fields(self.input_path, output, explicit), explicit)
        run = self.run_command([sys.executable, self.repo / "examples/seafloor_spreading_cases/plot_bxyz.py",
                                self.input_path, output])
        self.assertIn(str(figure), run.stdout)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib is optional")
    def test_batch_launchers_use_input_directory(self):
        scripts = self.make_launchers()
        self.run_command(["bash", scripts / "run_all_examples_end_to_end.sh", "point_source", 1737.4, self.case_dir])
        self.run_command([sys.executable, scripts / "run_all_examples_brtp.py", "--solver", "point_source",
                          "--examples-dir", self.case_dir, "--refine-factor", 1])
        self.assertTrue((self.case_dir / "output/small_xyz.txt").is_file())
        self.assertTrue((self.case_dir / "figs/small_brtp_2x2.png").is_file())
        self.assertFalse((scripts / "output").exists())
        self.assertFalse((scripts / "figs").exists())
        self.run_command([sys.executable, scripts / "run_all_examples_brtp.py", "--solver", "point_source",
                          "--examples-dir", self.case_dir, "--output-dir", "custom_output", "--figs-dir", "custom_figs"])
        self.assertTrue((scripts / "custom_output/small_xyz.txt").is_file())
        self.assertTrue((scripts / "custom_figs/small_brtp_2x2.png").is_file())
        self.run_command([sys.executable, scripts / "run_all_examples_brtp.py", "--solver", "gauss_legendre",
                          "--examples-dir", self.case_dir, "--refine-factor", 1])
        np.testing.assert_allclose(np.loadtxt(self.case_dir / "output/small_xyz.txt"),
                                   np.loadtxt(run_grid_model(self.input_path, solver="gauss_legendre")["output_path"]))

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib is optional")
    def test_observation_save_and_plot_defaults(self):
        from gravmagpy import save_fit
        from gravmagpy.plotting import plot_fit
        path = self.make_csv()
        field = np.zeros((2, 3))
        observations = {"table": pd.DataFrame({"source_file": [str(path), str(path)]}),
                        "xyz_km": np.array([[1767.4, 0, 0], [1768.4, 0, 0]]), "field_nt": field}
        result = {"sources": [], "radius_km": 1737.4, "success": True, "rmse_nt": 0.0,
                  **{f"{key}_nt": field for key in ("predicted", "crustal", "background", "residual")}}
        directory = save_fit(result, observations)
        self.assertEqual(directory, path.parent / "output")
        self.assertTrue((directory / "predictions.csv").is_file())
        self.assertEqual(plot_fit(observations, result), path.parent / "figs/observations_fit.png")
        with self.assertRaises(ValueError):
            save_fit(result, {key: value for key, value in observations.items() if key != "table"})

    @unittest.skipUnless(importlib.util.find_spec("rasterio") and importlib.util.find_spec("matplotlib"),
                         "map dependencies are optional")
    def test_saved_equivalent_map_uses_sibling_figs(self):
        from gravmagpy.maps import plot_equivalent_maps
        source_path = self.case_dir / "output/run/dipoles.csv"
        source_path.parent.mkdir(parents=True)
        pd.DataFrame({"lat_deg": [0], "lon_deg": [0], "radius_m": [1727400],
                      "mx_am2": [1e10], "my_am2": [0], "mz_am2": [0]}).to_csv(source_path, index=False)
        figure = plot_equivalent_maps(source_path, [], bounds=(-2, -2, 2, 2), shape=(4, 4), dpi=60)
        self.assertEqual(figure, self.case_dir / "figs/dipoles_context_30km.png")
        self.assertTrue(figure.is_file())
