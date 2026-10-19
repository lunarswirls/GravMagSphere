"""build cached fortran targets and run gravmag sphere command-line programs"""

import hashlib
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile


def fortran_source_dir():
    """locate sources in a checkout or an installed wheel"""
    package = Path(__file__).resolve().parents[1]
    candidates = [package / "_fortran", package.parents[1] / "fortran"]
    for candidate in candidates:
        if (candidate / "orbital.f90").is_file():
            return candidate
    raise FileNotFoundError("Fortran sources are missing; reinstall gravmagpy from its source distribution or wheel")


def _solver_name(solver):
    """resolve legacy selectors to the method names used by the backends"""
    return {"direct": "direct_sum", "gauss": "spectral"}.get(solver, solver)


def build_fortran(target="orbital", *, compiler=None, build_dir=None, force=False):
    """compile a target once per source/compiler/flag fingerprint and return its path

    targets are orbital (shared library), direct_sum, spectral, gauss_legendre, xyz_to_brtp, and
    dipole_grid (executables); an installed gfortran compiler is required
    """
    target = _solver_name(target)
    targets = {
        "orbital": ["orbital.f90"],
        "direct_sum": [f"gravmag_sphere_{part}.f90" for part in ("state", "subs", "physics", "direct_sum")],
        "spectral": ["gravmag_sphere_harmonics.f90", "gravmag_sphere_spectral.f90"],
        "gauss_legendre": ["gravmag_sphere_gauss_legendre.f90"],
        "xyz_to_brtp": ["gravmag_xyz_to_brtp.f90"],
        "dipole_grid": ["gravmag_sphere_dipole_grid_fit.f90"],
    }
    if target not in targets:
        raise ValueError(f"unknown target {target!r}; choose from {tuple(targets)}")
    executable = shutil.which(str(compiler or os.environ.get("FC", "gfortran")))
    if executable is None:
        raise RuntimeError("gfortran was not found; install it or pass compiler='/path/to/gfortran'")
    version = subprocess.run([executable, "--version"], check=True, capture_output=True, text=True).stdout
    source_dir = fortran_source_dir()
    sources = [source_dir / name for name in targets[target]]
    if target != "orbital":
        sources.insert(0, source_dir / "gravmag_paths.f90")
    flags = ["-std=f2008", "-O3", "-ffree-line-length-none"]
    suffix = ".exe" if sys.platform == "win32" else ""
    if target == "orbital":
        flags += ["-shared", "-fPIC"]
        suffix = ".dylib" if sys.platform == "darwin" else ".dll" if sys.platform == "win32" else ".so"
    fingerprint = hashlib.sha256((executable + version + repr(flags) + platform.platform()).encode())
    for source in sources:
        fingerprint.update(source.read_bytes())
    root = Path(build_dir).expanduser().resolve() if build_dir else Path(tempfile.gettempdir()) / "gravmagpy-build"
    directory = root / fingerprint.hexdigest()[:20]
    directory.mkdir(parents=True, exist_ok=True)
    output = directory / (target + suffix)
    if output.is_file() and not force:
        return output
    # private module/output paths prevent concurrent builds from sharing intermediate files
    with tempfile.TemporaryDirectory(prefix="compile-", dir=directory) as temporary:
        temporary = Path(temporary)
        binary = temporary / output.name
        command = [executable, *flags, f"-J{temporary}", f"-I{temporary}", *map(str, sources), "-o", str(binary)]
        result = subprocess.run(command, cwd=temporary, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(f"Fortran build failed for {target}:\n{result.stdout}\n{result.stderr}")
        os.replace(binary, output)
    return output


def run_fortran(target, arguments, *, compiler=None, build_dir=None, cwd=None, timeout=300):
    """compile and run a gravmag sphere target with positional arguments

    return subprocess.CompletedProcess; treat fortran stop 'error' messages as
    failures even when the fortran runtime returns exit status zero
    """
    target = _solver_name(target)
    if target == "orbital":
        raise ValueError("orbital is a shared library; use gravmagpy.predict_field")
    executable = build_fortran(target, compiler=compiler, build_dir=build_dir)
    result = subprocess.run(
        [str(executable), *map(str, arguments)], cwd=cwd, timeout=timeout, capture_output=True, text=True
    )
    if result.returncode or "error" in (result.stdout + result.stderr).lower():
        raise RuntimeError(f"Fortran {target} failed:\n{result.stdout}\n{result.stderr}")
    return result
