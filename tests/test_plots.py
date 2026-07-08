#python -m unittest tests/test_plots.py 

import unittest
import pytest
import numpy as np
from numpy.testing import assert_allclose
import os

import gne.gne_const as c
import gne.gne_plots as plt

ex_root = 'output/iz61/ivol'
ex_end = 'ex.hdf5'
required_file = ex_root+'0/'+ex_end

@pytest.mark.skipif(not os.path.exists(required_file), 
                   reason=f"Test data file not found: {required_file}")

class TestPredict(unittest.TestCase):
    def test_contour2Dsigma(self):
        levels,colors = plt.contour2Dsigma()
        assert_allclose(levels[1:],c.sigma_2Dprobs,rtol=0.01)
        self.assertEqual(colors[0][-1],1.0)

        nl=3
        levels,colors = plt.contour2Dsigma(n_levels=nl)
        assert_allclose(levels[1:],c.sigma_2Dprobs[0:nl],rtol=0.01)

                
if __name__ == '__main__':
    unittest.main()
