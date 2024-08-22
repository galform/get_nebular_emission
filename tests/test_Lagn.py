#python3 -m unittest tests/test_Lagn.py 

from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose

import src.gne_Lagn as agn

class TestPredict(TestCase):
    def test_Ledd(self):
        mbh0 = 1e-38; val0 = 1.26
        self.assertAlmostEqual(agn.get_Ledd(mbh0),val0,2)
    
        nlen = 2
        mbh,val = [np.zeros(3) for i in range(nlen)]
        mbh.fill(mbh0); val.fill(val0)
        assert_allclose(agn.get_Ledd(mbh),val,rtol=0.01)

    
    def test_acc_rate_edd(self):
        mbh0 = 1e8; val0 = 2.2
        self.assertAlmostEqual(agn.acc_rate_edd(mbh0),val0,1)
    
        nlen = 2
        mbh,val = [np.zeros(3) for i in range(nlen)]
        mbh.fill(mbh0); val.fill(val0)
        assert_allclose(agn.acc_rate_edd(mbh),val,rtol=0.1)
        

if __name__ == '__main__':
    unittest.main()
