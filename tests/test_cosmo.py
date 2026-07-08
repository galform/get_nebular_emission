#python -m unittest tests/test_cosmo.py 

import unittest
import numpy as np
from numpy.testing import assert_allclose
import os

import gne.gne_const as c
import gne.gne_cosmology as cosmo

verbose = False

# Values from Reyes-Peraza+25
redshift = 1.321
omega0 = 0.3089
lambda0 = 1-omega0
omegab = 0.0486
h0 = 0.6774

class TestPredict(unittest.TestCase):
    def test_ndeg2nV(self):
        cosmo.set_cosmology(omega0=omega0, omegab=omegab,
                            lambda0=lambda0,h0=h0)
        # Table 1 in Reyes-Peraza+25
        z1 = 0.9; z2 = 1.8
        ndeg = [4377,2268,580.7]
        nMpc = [1.299*1e-3,6.731*1e-4,1.723*1e-4]

        # Conversion
        for ii in range(len(ndeg)):
            nV=cosmo.ndeg2nV(ndeg[ii],z1,z2,verbose=verbose)
            self.assertAlmostEqual(nV,nMpc[ii],places=4)

if __name__ == '__main__':
    unittest.main()
