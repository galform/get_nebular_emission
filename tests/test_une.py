#python3 -m unittest tests/test_Lagn.py 

from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose

import src.gne_const as c
import src.gne_une as une

class TestPredict(TestCase):
    def test_Ledd(self):
        mbh0 = 1e-38; val0 = 1.26
        self.assertAlmostEqual(agn.get_Ledd(mbh0),val0,2)
    
        nlen = 3
        mbh,val = [np.zeros(nlen) for i in range(2)]
        mbh.fill(mbh0); val.fill(val0)
        assert_allclose(agn.get_Ledd(mbh),val,rtol=0.01)

    
        

if __name__ == '__main__':
    unittest.main()
