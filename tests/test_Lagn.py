#python3 -m unittest tests/test_Lagn.py 

from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose

import src.gne_const as c
import src.gne_Lagn as agn

class TestPredict(TestCase):
    def test_Ledd(self):
        mbh0 = 1e-38; val0 = 1.26
        self.assertAlmostEqual(agn.get_Ledd(mbh0),val0,2)
    
        nlen = 3
        mbh,val = [np.zeros(nlen) for i in range(2)]
        mbh.fill(mbh0); val.fill(val0)
        assert_allclose(agn.get_Ledd(mbh),val,rtol=0.01)

    
    def test_acc_rate_edd(self):
        mbh0 = 1e8; val0 = 2.2
        self.assertAlmostEqual(agn.acc_rate_edd(mbh0),val0,1)

        nlen = 3
        mbh,val = [np.zeros(nlen) for i in range(2)]
        mbh.fill(mbh0); val.fill(val0)
        assert_allclose(agn.acc_rate_edd(mbh),val,rtol=0.1)

        
    def test_acc_rate_radio(self):
        mhot0=1e19; mbh0=1; k0=1; ke0=1; val0=1.0
        self.assertAlmostEqual(agn.acc_rate_radio(mhot0,mbh0,kagn=k0,
                                                  kagn_exp=ke0),val0,2)

        nlen = 3
        mhot,mbh,val = [np.zeros(nlen) for i in range(3)]
        mhot.fill(mbh0); mbh.fill(mbh0); val.fill(val0)
        assert_allclose(agn.acc_rate_radio(mhot0,mbh0,kagn=k0,
                                           kagn_exp=ke0),val0,rtol=0.01)

        
    def test_Lagn_M16(self):
        mdot0 = c.yr_to_s/(c.c*c.c*1e7*c.Msun)
        val0 = c.e_r_agn*(1-c.e_f_agn)
        self.assertAlmostEqual(agn.get_Lagn_M16(mdot0),val0,3)

        nlen = 3
        y,val = [np.zeros(nlen) for i in range(2)]
        y.fill(mdot0); val.fill(val0)
        assert_allclose(agn.get_Lagn_M16(y),val,rtol=0.001)

        
    def test_Lagn_H14(self):
        # fedd>0.1
        mdot0 = c.yr_to_s/(c.c*c.c*1e7*c.Msun); mbh0 = 1e-40; val0 = 0.085
        self.assertAlmostEqual(agn.get_Lagn_H14(mdot0,mbh0),val0,3)

        nlen = 3
        mdot,mbh,val = [np.zeros(nlen) for i in range(3)]
        mdot.fill(mdot0); mbh.fill(mbh0); val.fill(val0)
        assert_allclose(agn.get_Lagn_H14(mdot,mbh),val,rtol=0.001)
        
        # fedd<0.1
        mdot0 = c.yr_to_s/(c.c*c.c*1e7*c.Msun); mbh0 = 1e-38; val0 = 0.079
        self.assertAlmostEqual(agn.get_Lagn_H14(mdot0,mbh0),val0,2)
        
        mdot.fill(mdot0); mbh.fill(mbh0); val.fill(val0)
        assert_allclose(agn.get_Lagn_H14(mdot,mbh),val,rtol=0.01)
        

if __name__ == '__main__':
    unittest.main()
