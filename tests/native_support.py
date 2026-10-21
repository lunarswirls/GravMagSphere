"""independent numerical references and native validation bindings"""

import ctypes
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile

import mpmath as mp
import numpy as np
from scipy.special import gammaln, lpmv


def native_library():
    root = Path(__file__).resolve().parents[1]
    sources = [root/name for name in ('fortran/gravmag_sphere_harmonics.f90',
               'fortran/gravmag_sphere_transform.f90', 'tests/fixtures/native_probe.f90')]
    flags = ['-std=f2008', '-O3', '-shared', '-fPIC', '-fopenmp', '-ffree-line-length-none']
    fingerprint = hashlib.sha256(repr(flags).encode()+b''.join(path.read_bytes() for path in sources)).hexdigest()[:20]
    directory = Path(tempfile.gettempdir())/'gravmag-native-validation'/fingerprint
    directory.mkdir(parents=True, exist_ok=True)
    binary = directory/'native_probe.so'
    if not binary.exists():
        subprocess.run([os.environ.get('FC','gfortran'), *flags, f'-J{directory}', f'-I{directory}',
                        *map(str,sources), '-o', str(binary)],check=True,capture_output=True,text=True)
    library = ctypes.CDLL(str(binary))
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
    complex_array = np.ctypeslib.ndpointer(dtype=np.complex128, flags='C_CONTIGUOUS')
    library.native_roundtrip_probe.argtypes = [ctypes.c_int]*3+[ctypes.c_double,array,array,array]
    library.native_selected_probe.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_double,array]
    library.native_fft_probe.argtypes = [ctypes.c_int,complex_array,ctypes.c_int]
    library.native_nodes_probe.argtypes = [ctypes.c_int,array,array]
    library.native_point_probe.argtypes = [ctypes.c_int,ctypes.c_double,ctypes.c_double,array,array]
    for name in ('native_roundtrip_probe','native_selected_probe','native_fft_probe','native_nodes_probe','native_point_probe'):
        getattr(library,name).restype = None
    return library


def jacobi_reference(l,m,theta,precision=80):
    """evaluate an independent jacobi-polynomial identity in arbitrary precision"""
    with mp.workdps(precision):
        angle = mp.mpf(float(theta))
        sine, cosine = mp.sin(angle), mp.cos(angle)
        n = l-m
        amplitude = mp.sqrt((2-int(m==0))*mp.factorial(l-m)/mp.factorial(l+m))
        amplitude *= mp.factorial(l+m)/(2**m*mp.factorial(l))
        # parity avoids a poorly conditioned hypergeometric evaluation near the south pole
        jacobi = (-1 if cosine<0 and n%2 else 1)*mp.jacobi(n,m,m,abs(cosine))
        value = amplitude*sine**m*jacobi
        derivative = mp.mpf(0)
        if m:
            derivative = amplitude*m*sine**(m-1)*cosine*jacobi
        if n:
            shifted = (-1 if cosine<0 and (n-1)%2 else 1)*mp.jacobi(n-1,m+1,m+1,abs(cosine))
            derivative -= amplitude*sine**(m+1)*mp.mpf(l+m+1)/2*shifted
        quotient = value/sine if m else mp.mpf(0)
        return np.array([float(value),float(derivative),float(quotient)])


def scipy_functions(degree,theta):
    """low-degree reference using scipy's phase-bearing unnormalized functions"""
    values,derivatives = [],[]
    sine,cosine = np.sin(theta),np.cos(theta)
    for l in range(degree+1):
        for m in range(l+1):
            factor = (-1)**m*np.sqrt((2-int(m==0))*np.exp(gammaln(l-m+1)-gammaln(l+m+1)))
            value = lpmv(m,l,cosine)
            previous = lpmv(m,l-1,cosine) if l>m else 0.
            values.append(factor*value)
            derivatives.append(factor*(l*cosine*value-(l+m)*previous)/sine)
    return np.array(values),np.array(derivatives)
