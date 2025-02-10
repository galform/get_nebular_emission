#python3 -m unittest tests/test_plots.py 

from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose
import os

import src.gne_const as c
import src.gne_plots as plt

ex_root = 'output/iz61/GP20_31p25kpc_z0_example_vol'

class TestPredict(TestCase):
    def test_get_obs_bpt(self):
        z=0.1
        valx,valy,obsdata=plt.get_obs_bpt(z,'NII')
        self.assertTrue(obsdata)
        self.assertAlmostEqual(valy[0],
                               np.log10(1.05924e-15/3.26841e-15),places=5)
        self.assertAlmostEqual(valx[0],
                               np.log10(3.95472e-15/1.00695e-14),places=5)
        valx,valy,obsdata=plt.get_obs_bpt(z,'SII')
        self.assertTrue(obsdata)
        self.assertAlmostEqual(valx[0],
                               np.log10(3.08815e-15/1.00695e-14),places=5)
        
        z=1.5
        valx,valy,obsdata=plt.get_obs_bpt(z,'NII')
        self.assertTrue(obsdata)
        self.assertAlmostEqual(valx[0],-0.5654081,places=5)
        self.assertAlmostEqual(valy[0],0.3271674,places=5)        
        valx,valy,obsdata=plt.get_obs_bpt(z,'SII')
        self.assertTrue(obsdata)
        self.assertAlmostEqual(valx[0],-0.4765437,places=5)
        self.assertAlmostEqual(valy[0],0.5962023,places=5)        

        z=2
        valx,valy,obsdata=plt.get_obs_bpt(z,'NII')
        self.assertFalse(obsdata)

        
    def test_plot_bpts(self):
        outplot = 'output/iz61/plots/bpts.pdf'
        self.assertEqual(plt.plot_bpts(ex_root),outplot)
        self.addCleanup(os.remove,outplot)
        
if __name__ == '__main__':
    unittest.main()
