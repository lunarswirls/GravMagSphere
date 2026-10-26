"""Independent physical, covariance, and staged inversion checks."""
import tempfile
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

from gravmagpy import build_fortran, spherical_to_cartesian
from gravmagpy.ravat import (monopole_kernel, spherical_basis, along_track_differences,
    fit_ravat_l2resl1, predict_ravat_field, save_ravat_fit, regional_monopole_grid, _irls,
    load_moonbroom_wake, _vector_data_system, _difference_data_system)


class RavatTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.library = build_fortran('orbital')

    def test_native_kernel_matches_negative_potential_gradient(self):
        xyz = np.array([[1760.,10.,20.],[1750.,-30.,4.]])
        sources = np.array([[1720.,0.,0.],[1710.,-15.,10.]])
        actual = monopole_kernel(xyz,sources,library=self.library)
        eps = 1e-3
        expected = np.empty_like(actual)
        for j in range(3):
            offset = np.eye(3)[j]*eps
            plus = 1/np.linalg.norm(xyz[:,None,:]+offset-sources[None,:,:],axis=2)
            minus = 1/np.linalg.norm(xyz[:,None,:]-offset-sources[None,:,:],axis=2)
            expected[:,j,:] = -(plus-minus)/(2*eps)
        np.testing.assert_allclose(actual,expected,rtol=1e-8,atol=1e-12)
        with self.assertRaises(ValueError):
            monopole_kernel(sources,sources,library=self.library)

    def observations(self):
        lat = np.linspace(5.,10.,60)
        lon = 302+np.sin(np.arange(60))*0.5
        xyz = spherical_to_cartesian(lat,lon,1767.4)
        return dict(xyz_km=xyz, field_nt=np.zeros_like(xyz),sigma_nt=[0.1,0.2,0.3],
                    table=pd.DataFrame({'utc':pd.date_range('1999-01-01',periods=60,freq='5s',tz='UTC'),
                                        'source_file':'one','source_row':np.arange(60)}))

    def test_pairing_frames_covariance_and_gaps(self):
        obs = self.observations()
        # A constant spherical field cancels; constant Cartesian fields generally don't.
        rotation = spherical_basis(obs['xyz_km'])
        obs['field_nt'] = np.einsum('nci,c->ni',rotation,[5.,-2.,1.])
        obs['table'].loc[2,'source_row'] = 99
        obs['table'].loc[4,'source_file'] = 'two'
        differences = along_track_differences(obs)
        np.testing.assert_allclose(differences['difference_nt'],0,atol=1e-14)
        self.assertNotIn(2,differences['pairs'].ravel())
        self.assertNotIn(4,differences['pairs'].ravel())
        w = differences['whitener'];cov = differences['covariance_nt2']
        identity = np.einsum('nij,njk,nlk->nil',w,cov,w)
        np.testing.assert_allclose(identity,np.broadcast_to(np.eye(3),identity.shape),atol=1e-14)
        self.assertGreater(np.max(np.abs(cov[:,0,1])),0.001)
        self.assertEqual(len(np.unique(differences['pairs'])),differences['pairs'].size)

    def test_staged_fit_recovers_heldout_field_and_combines_vectors(self):
        obs = self.observations()
        sources = spherical_to_cartesian([6.5,8.5],[301.8,302.2],[1727.4,1727.4])
        charges = np.array([15000.,-13000.])
        # Independent Coulomb field to generate truth; do not call native kernel.
        delta = obs['xyz_km'][:,None,:]-sources[None,:,:]
        obs['field_nt'] = np.sum(charges[None,:,None]*delta/np.linalg.norm(delta,axis=2)[:,:,None]**3,axis=1)
        differences = along_track_differences(obs)
        surface = spherical_to_cartesian([6,7,8,9],[301.5,302,302.5,302],1737.4)
        fit = fit_ravat_l2resl1(differences,sources,surface,alpha_l2=0.001,alpha_l1=0.0001,
                               tukey_c=None,max_iterations=50,library=self.library)
        self.assertTrue(fit['success'])
        heldout = spherical_to_cartesian([6.2,7.7,9.1],[301.9,302.4,301.8],1757.4)
        delta = heldout[:,None,:]-sources[None,:,:]
        truth = np.sum(charges[None,:,None]*delta/np.linalg.norm(delta,axis=2)[:,:,None]**3,axis=1)
        total = predict_ravat_field(heldout,fit,library=self.library)
        np.testing.assert_allclose(total['xyz_nt'],truth,rtol=1e-4,atol=1e-5)
        baseline = predict_ravat_field(heldout,fit,part='baseline',library=self.library)
        correction = predict_ravat_field(heldout,fit,part='correction',library=self.library)
        np.testing.assert_allclose(total['xyz_nt'],baseline['xyz_nt']+correction['xyz_nt'],atol=1e-13)
        np.testing.assert_allclose(total['total_nt'],np.linalg.norm(total['rtp_nt'],axis=1),atol=1e-13)
        self.assertLess(np.linalg.norm(fit['residual_difference_nt']),np.linalg.norm(fit['baseline_difference_nt']-differences['difference_nt']))
        for stage in ['l1','l2']:
            self.assertTrue(np.all(np.diff(fit['history_'+stage].objective) <= 1e-7))
        with tempfile.TemporaryDirectory() as directory:
            save_ravat_fit(fit,differences,directory)
            with np.load(Path(directory)/'solution.npz') as reloaded:
                np.testing.assert_allclose(predict_ravat_field(heldout,reloaded,library=self.library)['xyz_nt'],total['xyz_nt'])
        with self.assertRaisesRegex(ValueError,'memory'):
            fit_ravat_l2resl1(differences,sources,surface,alpha_l2=1,alpha_l1=1,max_memory_mib=1e-5)

    def test_grid_source_depth_and_weights(self):
        grid = regional_monopole_grid((299,4,305,11),spacing_km=30,depth_km=10)
        np.testing.assert_allclose(np.linalg.norm(grid['surface_xyz_km'],axis=1)-
                                  np.linalg.norm(grid['source_xyz_km'],axis=1),10,atol=1e-12)
        expected = np.deg2rad(6)*(np.sin(np.deg2rad(11))-np.sin(np.deg2rad(4)))
        self.assertAlmostEqual(grid['area_weights'].sum(),expected)

    def test_l1_surface_penalty_matches_independent_optimizer(self):
        from scipy.optimize import minimize
        rng = np.random.default_rng(12)
        g, radial = rng.normal(size=(20,3)), rng.normal(size=(8,3))
        data, area = rng.normal(size=20), np.full(8,1/8)
        alpha, eps = 0.7, 0.1
        def objective(m):
            return np.sum((g@m-data)**2)+2*alpha**2*np.sum(area*np.sqrt((radial@m)**2+eps**2))
        reference = minimize(objective,np.zeros(3),method='BFGS',tol=1e-10)
        actual, _, success = _irls(g,data,radial,area,alpha=alpha,l1=True,epsilon_nt=eps,
                                   max_iterations=200,tolerance=1e-9,tukey_c=None)
        self.assertTrue(success)
        np.testing.assert_allclose(actual,reference.x,atol=1e-7)

    def test_tukey_rejects_gross_outlier(self):
        g, radial = np.ones((101,1)), np.ones((1,1))
        data = np.r_[np.ones(100),100.]
        actual, history, success = _irls(g,data,radial,np.ones(1),alpha=0.01,l1=False,
            epsilon_nt=0.1,max_iterations=50,tolerance=1e-8,tukey_c=4.5)
        self.assertTrue(success)
        self.assertAlmostEqual(actual[0],1.,places=5)
        self.assertEqual(history.active_rows.iloc[-1],100)

    def test_filtered_wake_directory_uses_dated_individual_products(self):
        xyz = spherical_to_cartesian([5,6,7],[302,302,302],1767.4)
        frame = pd.DataFrame(xyz,columns=['X_SEL','Y_SEL','Z_SEL'])
        frame[['Bx_SEL','By_SEL','Bz_SEL']] = 1.
        frame['wake_t_km'],frame['wake_rperp_km'] = 1000.,900.
        frame['utc'] = pd.date_range('1999-01-01',periods=3,freq='5s',tz='UTC')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            child = root/'alt_lt_40km_brms_90pct'
            child.mkdir()
            for day in ['990101','980101']:
                frame.to_csv(child/f'ma{day}_nightside_wake_alt_lt_40km_brms_90pct.csv',index=False)
            # Aggregate/detrended variants must not be combined with day files.
            (child/'combined_alt_lt_40km_brms_90pct.csv').write_text('not input\n')
            result = load_moonbroom_wake(root,bounds=(299,4,305,11),start='1999-01-01',stop='1999-01-02')
            self.assertEqual(len(result['xyz_km']),3)
            self.assertEqual(len(result['report']),1)
            self.assertEqual(result['table'].source_row.tolist(),[0,1,2])
            self.assertTrue(result['selection']['sigma_is_assumed'])
            self.assertEqual(len(along_track_differences(result)['pairs']),1)

    def test_joint_covariance_retains_vector_information_without_duplication(self):
        obs = self.observations()
        rng = np.random.default_rng(812)
        obs['field_nt'] = rng.normal(size=obs['xyz_km'].shape)
        # Unequal, anisotropic endpoint uncertainties imply nonzero mean/delta covariance.
        obs['sigma_nt'] = rng.uniform(0.1,1.,size=obs['xyz_km'].shape)
        obs['table'].loc[2,'source_row'] = 99  # Include unpaired samples.
        differences = along_track_differences(obs)
        spherical = rng.normal(size=(len(obs['xyz_km']),3,4))
        gv,dv = _vector_data_system(spherical,differences,'vectors')
        gj,dj = _vector_data_system(spherical,differences,'vectors_and_differences')
        self.assertEqual(len(dj),3*len(obs['xyz_km']))
        np.testing.assert_allclose(gj.T@gj,gv.T@gv,rtol=1e-12,atol=1e-10)
        np.testing.assert_allclose(gj.T@dj,gv.T@dv,rtol=1e-12,atol=1e-10)
        self.assertAlmostEqual(float(dj@dj),float(dv@dv),places=10)
        # Every observed component participates in the data system.
        for j in range(3):
            perturbed = dict(differences)
            perturbed['vector_nt'] = differences['vector_nt'].copy()
            perturbed['vector_nt'][:,j] += 1
            left,right = differences['pairs'].T
            perturbed['difference_nt'] = perturbed['vector_nt'][right]-perturbed['vector_nt'][left]
            _, new_target = _vector_data_system(spherical,perturbed,'vectors_and_differences')
            self.assertGreater(np.linalg.norm(new_target-dj),1)

    def test_all_adjacent_gradient_whitening_matches_full_covariance(self):
        from scipy.linalg import block_diag
        obs = self.observations()
        rng = np.random.default_rng(814)
        obs['field_nt'] = rng.normal(size=obs['xyz_km'].shape)
        obs['sigma_nt'] = rng.uniform(0.1,1.,size=obs['xyz_km'].shape)
        obs['table'].loc[2,'source_row'] = 99
        d = along_track_differences(obs,pairing='adjacent')
        self.assertEqual(len(d['pairs']),57)
        self.assertNotIn(2,d['pairs'].ravel())
        design = rng.normal(size=(len(obs['xyz_km']),3,4))
        g,target = _difference_data_system(design,d)
        difference = np.zeros((3*len(d['pairs']),3*len(obs['xyz_km'])))
        for k,(left,right) in enumerate(d['pairs']):
            difference[3*k:3*k+3,3*left:3*left+3] = -np.eye(3)
            difference[3*k:3*k+3,3*right:3*right+3] = np.eye(3)
        covariance = difference@block_diag(*d['vector_covariance_nt2'])@difference.T
        cholesky = np.linalg.cholesky(covariance)
        independent_g = np.linalg.solve(cholesky,difference@design.reshape(-1,4))
        independent_d = np.linalg.solve(cholesky,difference@d['vector_nt'].ravel())
        np.testing.assert_allclose(g,independent_g,rtol=1e-10,atol=1e-10)
        np.testing.assert_allclose(target,independent_d,rtol=1e-10,atol=1e-10)
        self.assertEqual(d['preprocessing']['detrending'],'none')
        # A component baseline cancels under finite differences, without a fit.
        altered = dict(obs)
        altered['field_nt'] = obs['field_nt']+np.einsum('nci,c->ni',
            spherical_basis(obs['xyz_km']),[100.,-50.,25.])
        shifted = along_track_differences(altered,pairing='adjacent')
        np.testing.assert_allclose(shifted['difference_nt'],d['difference_nt'],atol=1e-12)


if __name__ == '__main__':
    unittest.main()
