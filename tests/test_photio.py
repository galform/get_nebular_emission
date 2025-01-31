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
    
    def test_Zgrid(self):
        grid_str = ['0001','002','014','030']
        nout = len(grid_str)
        zout = np.array([0.0001,0.002,0.014,0.030])
        lzout = np.array([np.log10(zout[i]) for i in range(nout)])

        nz,zgrid,lzgrid = ph.get_Zgrid(grid_str)
        self.assertEqual(nz,nout)
        np.testing.assert_allclose(zgrid, zout)
        np.testing.assert_allclose(lzgrid, lzout)


    def test_interp_u_z(self):
        grid = np.random.rand(5, 10, 3)
        grid = np.array([[[1. 2.]]]) ##here
        print('grid: ',grid,type(grid))
        u = np.array([1.0, 2.0, c.notnum, 4.0, 5.0])
        ud = np.array([0.3, 0.7, 0.5, 0.2, 0.8])
        iu = np.array([0, 9, 3, 1, 3]) 
        zd = np.array([0.4, 0.6, 0.3, 0.8, 0.5])
        iz = np.array([1, 3, 2, 0, 4]) 
    
        emline = ph.interp_u_z(grid, u, ud, iu, zd, iz)
        # Check that the output has the correct shape
        self.assertEqual(emline.shape, (3, 5))    
        # Check that the output values are correct
        expected_emline = np.array([[grid[1, 2, 0]*(1-ud[0])*(1-zd[0]),
                                     grid[3, 4, 0]*ud[1]*(1-zd[1]),
                                     grid[2, 3, 0]*(1-ud[2])*zd[2],
                                     grid[0, 1, 0]*ud[3]*(1-zd[3]),
                                     grid[4, 3, 0]*ud[4]*zd[4]],
                                    [grid[1, 2, 1]*(1-ud[0])*(1-zd[0]),
                                     grid[3, 4, 1]*ud[1]*(1-zd[1]),
                                     grid[2, 3, 1]*(1-ud[2])*zd[2],
                                     grid[0, 1, 1]*ud[3]*(1-zd[3]),
                                     grid[4, 3, 1]*ud[4]*zd[4]],
                                    [grid[1, 2, 2]*(1-ud[0])*(1-zd[0]),
                                     grid[3, 4, 2]*ud[1]*(1-zd[1]),
                                     grid[2, 3, 2]*(1-ud[2])*zd[2],
                                     grid[0, 1, 2]*ud[3]*(1-zd[3]),
                                     grid[4, 3, 2]*ud[4]*zd[4]]
                                    ])
        np.testing.assert_allclose(emline, expected_emline)

        
if __name__ == '__main__':
    unittest.main()
