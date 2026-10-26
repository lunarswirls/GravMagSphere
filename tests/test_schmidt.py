"""validate native schmidt recurrences and vector projections without shtools"""

import os
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.special import roots_legendre
from gravmagpy import block, sphere_model, run_sphere_model, spherical_to_cartesian
from native_support import native_library, jacobi_reference, scipy_functions


def coefficient_order(degree):
    return [(part,l,m) for l in range(1,degree+1) for m in range(l+1)
            for part in range(1 if m==0 else 2)]


def spherical_axes(latitude,longitude):
    theta=np.deg2rad(90-np.asarray(latitude))
    phi=np.deg2rad(longitude)
    return np.stack((np.stack((np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)),axis=-1),
                     np.stack((np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)),axis=-1),
                     np.stack((-np.sin(phi),np.cos(phi),np.zeros_like(phi)),axis=-1)),axis=-2)


def independent_vector(degree,coefficients,theta,phi):
    p,dp=scipy_functions(degree,theta)
    field=np.zeros(3)
    for coefficient,(part,l,m) in zip(coefficients,coefficient_order(degree)):
        index=l*(l+1)//2+m
        trig=np.cos(m*phi) if part==0 else np.sin(m*phi)
        azimuth=np.sin(m*phi) if part==0 else -np.cos(m*phi)
        field += coefficient*np.array([(l+1)*p[index]*trig,-dp[index]*trig,m*p[index]/np.sin(theta)*azimuth])
    return field


class schmidt_tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory=tempfile.TemporaryDirectory(prefix='gravmag-native-test-')
        cls.library=native_library()
        cls.records=[]

    @classmethod
    def tearDownClass(cls):
        destination=os.environ.get('GRAVMAG_SCHMIDT_REPORT')
        if destination:
            path=Path(destination)
            path.parent.mkdir(parents=True,exist_ok=True)
            pd.DataFrame(cls.records).to_csv(path,index=False)
        cls.directory.cleanup()

    def compare(self,name,actual,reference,tolerance,**metadata):
        maximum=float(np.max(np.abs(actual-reference)))
        scaled=maximum/max(1.,float(np.max(np.abs(reference))))
        self.records.append({'check':name,'max_absolute_difference':maximum,
                             'scaled_max_difference':scaled,'tolerance':tolerance,**metadata})
        self.assertTrue(np.isfinite(actual).all())
        self.assertLess(scaled,tolerance,f'{name}: {metadata}')

    def test_scaled_functions_against_arbitrary_precision(self):
        for degree in (1,8,24,64,128,256):
            for m in sorted(set((0,1,degree//2,degree))):
                for theta in (1e-7,.3,1.2,np.pi-1e-7):
                    actual=np.zeros(3)
                    self.library.native_selected_probe(degree,m,theta,actual)
                    reference=jacobi_reference(degree,m,theta)
                    for index,name in enumerate(('legendre','colatitude_derivative','azimuthal_quotient')):
                        self.compare(name,actual[index:index+1],reference[index:index+1],2e-10,degree=degree,order=m,theta=theta)

    def test_native_fft_against_numpy(self):
        random=np.random.default_rng(271019)
        for count in (2,16,128,4096,32768,65536):
            original=random.normal(size=count)+1j*random.normal(size=count)
            values=original.copy()
            self.library.native_fft_probe(count,values,0)
            self.compare('fft_forward',values,np.fft.fft(original)/count,2e-13,count=count)
            self.library.native_fft_probe(count,values,1)
            self.compare('fft_inverse',values,original,2e-13,count=count)

    def test_underflowed_sectoral_seed_recovers_at_high_degree(self):
        degree,order,theta=2800,1204,np.pi/6
        actual=np.zeros(3)
        self.library.native_selected_probe(degree,order,theta,actual)
        reference=jacobi_reference(degree,order,theta)
        self.assertLess(order*np.log(np.sin(theta)),np.log(np.nextafter(0.,1.)))
        self.assertGreater(abs(reference[0]),.01)
        self.compare('underflow_recovery',actual,reference,2e-10,degree=degree,order=order)

    def test_extreme_polar_angles_keep_scaled_seeds_finite(self):
        for theta in (1e-100,1e-300):
            for order in (0,1,2,64,256):
                actual=np.zeros(3)
                self.library.native_selected_probe(256,order,theta,actual)
                reference=jacobi_reference(256,order,theta)
                self.compare('extreme_polar_angle',actual,reference,2e-10,degree=256,order=order,theta=theta)
                np.testing.assert_allclose(actual,reference,rtol=2e-10,atol=1e-320)

    def test_gauss_nodes_against_scipy(self):
        for count in (8,64,513):
            theta=np.zeros(count)
            weights=np.zeros(count)
            self.library.native_nodes_probe(count,theta,weights)
            nodes,reference=roots_legendre(count)
            self.compare('gauss_nodes',np.cos(theta),nodes[::-1],5e-14,count=count)
            self.compare('gauss_weights',weights,reference[::-1],5e-13,count=count)
            self.assertAlmostEqual(weights.sum(),2.,places=13)

    def test_vector_synthesis_against_independent_references(self):
        random=np.random.default_rng(281019)
        for degree in (1,8,24):
            coefficients=random.normal(size=degree*(degree+2))
            for theta,phi in ((.03,-2.3),(.7,.9),(1.5,0.),(2.4,1.9),(np.pi-.03,-.5)):
                actual=np.zeros(3)
                self.library.native_point_probe(degree,theta,phi,coefficients,actual)
                self.compare('vector_synthesis',actual,independent_vector(degree,coefficients,theta,phi),2e-10,degree=degree)
        for degree in (600,1200,2800):
            coefficients=np.zeros(degree*(degree+2))
            modes=[(degree,0,.4,0.),(degree,1,-.7,.3),
                   (degree,int(.43*degree),1.,-.6),(degree,degree,.2,.9)]
            for l,m,g,h in modes:
                index=l*l+max(0,2*m-1)-1
                coefficients[index]=g
                if m:
                    coefficients[index+1]=h
            for theta,phi in ((1e-7,-2.3),(np.pi/6,.9),(np.pi/2,1.9)):
                reference=np.zeros(3)
                for l,m,g,h in modes:
                    p,derivative,quotient=jacobi_reference(l,m,theta)
                    cosine,sine=np.cos(m*phi),np.sin(m*phi)
                    amplitude=g*cosine+h*sine
                    reference+=np.array([(l+1)*p*amplitude,-derivative*amplitude,
                                         m*quotient*(g*sine-h*cosine)])
                actual=np.zeros(3)
                self.library.native_point_probe(degree,theta,phi,coefficients,actual)
                self.compare('point_synthesis_high_degree',actual,reference,2e-10,degree=degree,theta=theta)

    def test_vector_transform_roundtrip(self):
        random=np.random.default_rng(301019)
        for degree in (8,64,128):
            coefficients=random.normal(size=degree*(degree+2))
            for joint in (0.,.3,1.,5.):
                recovered=np.empty_like(coefficients)
                statistics=np.zeros(3)
                nphi=2**int(np.ceil(np.log2(2*degree+1)))
                self.library.native_roundtrip_probe(degree,degree+1,nphi,joint,coefficients,recovered,statistics)
                self.compare('coefficient_roundtrip',recovered,coefficients,2e-11,degree=degree,joint_strength=joint)

    def test_exact_and_near_poles_against_dipole(self):
        moment=np.array([.7,-1.3,2.1])
        coefficients=moment[[2,0,1]].copy()
        for theta in (0.,1e-12,1e-8,np.pi-1e-8,np.pi-1e-12,np.pi):
            for longitude in (-150,0,67):
                axes=spherical_axes(90-np.rad2deg(theta),longitude)
                field=np.zeros(3)
                self.library.native_point_probe(1,theta,np.deg2rad(longitude),coefficients,field)
                reference=3*np.dot(moment,axes[0])*axes[0]-moment
                self.compare('polar_dipole',axes.T@field,reference,5e-14,theta=theta,longitude=longitude)

    def test_public_gravity_and_magnetic_fits_against_scipy_svd(self):
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

    def test_workspace_preflight_and_degree_bounds(self):
        source=block([6,8],[-60,-58],[2,8],density_kg_m3=400)
        model=sphere_model([source],latitude_deg=[6,7],longitude_deg=[-60,-59],altitude_km=5000)
        with self.assertRaisesRegex(ValueError,'46339'):
            run_sphere_model(model,solver='spectral',solver_options={'lmax':46334})
        with patch.dict(os.environ,{'GRAVMAG_SPECTRAL_MAX_MIB':'.001'}):
            with self.assertRaisesRegex(RuntimeError,'estimated spectral workspace'):
                run_sphere_model(model,solver='spectral',solver_options={'lmax':24})
        with patch.dict(os.environ,{'GRAVMAG_SPECTRAL_MAX_MIB':'nan'}):
            with self.assertRaisesRegex(RuntimeError,'invalid spectral memory limit'):
                run_sphere_model(model,solver='spectral',solver_options={'lmax':24})

    def check_source_fit(self, degree):
        radius = 1737.4
        altitude = 5000.0
        observation_radius = (radius+altitude)*1000
        center = spherical_to_cartesian(7, -59, radius-5)*1000
        volume = ((radius-5)*1000)**2*np.cos(np.deg2rad(7))*np.deg2rad(2)**2*6000
        order = coefficient_order(degree)
        ntheta = max(12, 2*degree+2)
        nphi = max(24, 4*degree+4)
        nphi = 2**int(np.ceil(np.log2(nphi)))
        nodes,quadrature_weights = roots_legendre(ntheta)
        theta,phi = np.meshgrid(np.arccos(nodes[::-1]),np.arange(nphi)*2*np.pi/nphi,indexing="ij")
        quadrature_weights = quadrature_weights[::-1]
        latitude = 90-np.rad2deg(theta.ravel())
        longitude = np.rad2deg(phi.ravel())
        axes = spherical_axes(latitude, longitude)
        displacement = observation_radius*axes[:, 0]-center
        distance = np.linalg.norm(displacement, axis=1)
        # use scipy functions to assemble the independent weighted design matrix
        design = np.zeros((ntheta, nphi, 3, len(order)))
        for row in range(ntheta):
            angle = theta[row, 0]
            p, derivative = scipy_functions(degree,angle)
            for column, (part, l, m) in enumerate(order):
                packed = l*(l+1)//2+m
                cosine = np.cos(m*phi[row])
                sine = np.sin(m*phi[row])
                trig = cosine if part == 0 else sine
                design[row, :, 0, column] = (l+1)*p[packed]*trig
                design[row, :, 1, column] = -derivative[packed]*trig
                design[row, :, 2, column] = m*p[packed]/np.sin(angle)*(sine if part == 0 else -cosine)
        design = design.reshape(-1, len(order))
        weights = np.repeat(np.repeat(quadrature_weights,nphi),3)
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
                    vector = independent_vector(degree,coefficients,np.deg2rad(90-lat),np.deg2rad(lon))
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
                                     "relative_vector_difference": relative, "tolerance": 2e-10})
                self.assertLess(relative, 2e-10)
                self.assertEqual(result["normalization"], "schmidt")
                self.assertFalse(result["condon_shortley_phase"])
                header = result["output_path"].read_text()
                self.assertIn("normalization=schmidt_semi_normalized condon_shortley_phase=excluded", header)
