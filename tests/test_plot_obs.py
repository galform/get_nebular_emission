#python -m unittest tests/test_plot_obs.py 

import unittest
import pytest
import numpy as np
from numpy.testing import assert_allclose
import os

import gne.gne_const as c
import gne.gne_plot_obs as obs

ex_root = 'output/iz61/ivol'
ex_end = 'ex.hdf5'
required_file = ex_root+'0/'+ex_end
omega0 = 0.3
omegab = 0.02
lambda0 = 0.7
h0 = 0.67
metadata = {'redshift': 1.03,
    'omega0': omega0,'omegab': omegab,
    'lambda0': lambda0, 'h0': h0,}

@pytest.mark.skipif(not os.path.exists(required_file), 
                   reason=f"Test data file not found: {required_file}")

class TestPredict(unittest.TestCase):
    def test_get_obs_bpt(self):
        z=0.1
        valx,valy,obsdata=obs.get_obs_bpt(z,'NII')
        self.assertTrue(obsdata)
        self.assertAlmostEqual(valy[0],
                               np.log10(1.05924e-15/3.26841e-15),places=5)
        self.assertAlmostEqual(valx[0],
                               np.log10(3.95472e-15/1.00695e-14),places=5)
        valx,valy,obsdata=obs.get_obs_bpt(z,'SII')
        self.assertTrue(obsdata)
        self.assertAlmostEqual(valx[0],
                               np.log10(3.08815e-15/1.00695e-14),places=5)
        
        z=1.5
        valx,valy,obsdata=obs.get_obs_bpt(z,'NII')
        self.assertTrue(obsdata)
        self.assertAlmostEqual(valx[0],-0.5654081,places=5)
        self.assertAlmostEqual(valy[0],0.3271674,places=5)        
        valx,valy,obsdata=obs.get_obs_bpt(z,'SII')
        self.assertTrue(obsdata)
        self.assertAlmostEqual(valx[0],-0.4765437,places=5)
        self.assertAlmostEqual(valy[0],0.5962023,places=5)        

        z=2
        valx,valy,obsdata=obs.get_obs_bpt(z,'NII')
        self.assertFalse(obsdata)

    def test_get_pozzetti(self):
        xx = np.array([0.5, 1, 2, 3, 5])
        x  = np.log10(xx)-16
        h3 = h0*h0*h0
        y  = -3.2321765

        metatest = {'redshift': 2.5,
                    'omega0': omega0,'omegab': omegab,
                    'lambda0': lambda0, 'h0': h0,}
        valx, valy, pp = obs.get_pozzetti(metadata=metatest)
        self.assertTrue(pp)
        np.testing.assert_allclose(valx,x, atol=0.001)
        self.assertAlmostEqual(valy[0],y,places=5)
            
        for z in [0, 0.2, 2.6]:
            metatest = {'redshift': z,
                        'omega0': omega0,'omegab': omegab,
                        'lambda0': lambda0, 'h0': h0,}
            valx, valy, pp = obs.get_pozzetti(metadata=metatest)
            self.assertFalse(pp)
            self.assertAlmostEqual(valx,-999,places=5) 
            self.assertAlmostEqual(valy,-999,places=5) 
        
        om = str(omega0)
        ob = str(omegab)
        l0 = str(lambda0)
        hh = str(h0)
        cosmo_str = ('_m'+om+'_b'+ob+'_l'+l0+'_h'+hh).replace('.','p')
        nomtabMpc = 'pozzettiMpc'+cosmo_str+'.txt'
        pozzettiMpc = os.path.join('output',nomtabMpc)
        pp = os.path.isfile(pozzettiMpc)
        self.assertTrue(pp, f"File {pozzettiMpc} does not exist")

if __name__ == '__main__':
    unittest.main()
