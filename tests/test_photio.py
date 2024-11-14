#python3 -m unittest tests/test_photio.py 

from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose

import src.gne_const as c
import src.gne_photio as ph

class TestPredict(TestCase):
    def test_get_limits(self):
        lims = np.zeros(2);
        lims[0]= 10; lims[1]=10000
        assert_allclose(ph.get_limits(propname ='nH',photmod='gutkin16'),
                        lims,rtol=0.01)

    
        

if __name__ == '__main__':
    unittest.main()
