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


    def test_get_lines_gutkin16(self):
        uu = np.array([[c.notnum,-3.7,-4.,-6.,-1.7,-1.7]])
        zz = np.array([[0.02,c.notnum,0.015,0.015,0.0001,2e-5]])
        nh = np.array([[1.,1.,200.,200.,5000.,5000.]])

        el = ph.get_lines_gutkin16(uu,nh,zz,verbose=True)
        # Check that the output has the correct shape
        self.assertEqual(el.shape, (1, 18, 6))    
        # Check that the boundaries are handled as expected
        np.testing.assert_allclose(el[0,:,0],el[0,:,1])
        np.testing.assert_allclose(el[0,:,2],el[0,:,3])
        np.testing.assert_allclose(el[0,:,4],el[0,:,5])

        
    def test_get_lines_feltre16(self):
        uu = np.array([[c.notnum,-3.7,-5.,-6.,-1.7,-1.7]])
        zz = np.array([[0.02,c.notnum,0.015,0.015,1e-4,0.]])
        nh = np.array([[1.,1.,200.,200.,2000.,2000.]])

        el = ph.get_lines_feltre16(uu,nh,zz,verbose=True)
        # Check that the output has the correct shape
        self.assertEqual(el.shape, (1, 20, 6))    
        # Check that the boundaries are handled as expected
        np.testing.assert_allclose(el[0,:,0],el[0,:,1])
        np.testing.assert_allclose(el[0,:,2],el[0,:,3])
        np.testing.assert_allclose(el[0,:,4],el[0,:,5])

        
    def test_get_lines(self):
        uu = np.array([[-3.7,-1.,0],[-3.7,-1.7,-1.7]])
        zz = np.array([[0.02,0.01,0.01],[0.02,0.03,0.07]])
        nh = np.array([[10,200,200],[5,5000,5000]])   
        el = ph.get_lines(uu,nh,zz, photmod='gutkin16', verbose=True)
        # Check that the output has the correct shape
        self.assertEqual(el.shape, (2, 18, 3))    
        # Check that the boundaries are handled as expected
        np.testing.assert_allclose(el[0,:,0],el[1,:,0])
        np.testing.assert_allclose(el[0,:,1],el[0,:,2])
        np.testing.assert_allclose(el[1,:,1],el[1,:,2])

        uu = np.array([[-3.7,-1.,0],[-3.7,-1.7,-1.7]])
        zz = np.array([[0.02,0.01,0.01],[0.02,0.07,1.]])
        nh = np.array([[100,200,200],[5,2000,2000]])   
        el = ph.get_lines(uu,nh,zz, photmod='feltre16', verbose=True)
        # Check that the output has the correct shape
        self.assertEqual(el.shape, (2, 20, 3))    
        # Check that the boundaries are handled as expected
        np.testing.assert_allclose(el[0,:,0],el[1,:,0])
        np.testing.assert_allclose(el[0,:,1],el[0,:,2])
        np.testing.assert_allclose(el[1,:,1],el[1,:,2])

        
if __name__ == '__main__':
    unittest.main()
