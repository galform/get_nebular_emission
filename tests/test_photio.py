#python -m unittest tests/test_photio.py 

import unittest
import shutil
import os
import numpy as np
from numpy.testing import assert_allclose

import gne.gne_const as c
import gne.gne_io as io
import gne.gne_photio as ph

class TestPredict(unittest.TestCase):
    decimals = 8
    atol = 1e-7

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures before running tests in the class"""
        cls.opath = 'data/example_data/iz61/'
        cls.root = 'GP20_31p25kpc_z0_example_vol0'
        cls.txtfile = cls.opath + cls.root + '.txt'
        cls.hf5file = cls.opath + cls.root + '.hdf5'
        cls.expath = 'output/test/'
        cls.zz = 0.
        cls.filenom = io.generate_header(cls.hf5file, cls.zz, 100,
                                         0.7, 0.4, 0.3, 0.6,
                                         1e8, 100, 1e6,
                                         outpath=cls.expath, verbose=False)
        cls.noms = ['photmod_sfr','xid_sfr','co_sfr','imf_cut_sfr','nH_sfr_cm3',
                    'photmod_NLR','xid_NLR','alpha_NLR','nH_NLR_cm3']
        cls.nval = ['gutkin16',0.3,1,100,c.nH_sfr_cm3,
                    'feltre16',0.5,-1.7,c.nH_NLR_cm3]
        cls.num = io.add2header(cls.filenom, cls.noms, cls.nval, verbose=False)

    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests in the class are finished"""
        if os.path.exists(cls.expath):
            try:
                shutil.rmtree(cls.expath)
            except (OSError, PermissionError) as e:
                print(f"Warning: Could not remove test directory {cls.expath}: {e}")
    
    
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
    
        el = ph.get_lines_gutkin16(uu,zz,self.filenom,lnH=nh,verbose=False)
        # Check that the output has the correct shape
        self.assertEqual(el.shape, (1, 18, 6))    
        # Check that the boundaries are handled as expected
        np.testing.assert_allclose(el[0,:,0],el[0,:,1])
        np.testing.assert_allclose(el[0,:,2],el[0,:,3])
        np.testing.assert_allclose(el[0,:,4],el[0,:,5])

        el = ph.get_lines_gutkin16(uu,zz,self.filenom,verbose=False)
        self.assertEqual(el.shape, (1, 18, 6))    
        
      
    def test_get_lines_feltre16(self):
        uu = np.array([[c.notnum,-3.7,-5.,-6.,-1.7,-1.7]])
        zz = np.array([[0.02,c.notnum,0.015,0.015,1e-4,0.]])
        nh = np.array([[1.,1.,200.,200.,2000.,2000.]])
    
        el = ph.get_lines_feltre16(uu,zz,self.filenom,lnH=nh,verbose=False)
        # Check that the output has the correct shape
        self.assertEqual(el.shape, (1, 20, 6))    
        # Check that the boundaries are handled as expected
        np.testing.assert_allclose(el[0,:,0],el[0,:,1])
        np.testing.assert_allclose(el[0,:,2],el[0,:,3])
        np.testing.assert_allclose(el[0,:,4],el[0,:,5])

        el = ph.get_lines_feltre16(uu,zz,self.filenom,verbose=False)
        self.assertEqual(el.shape, (1, 20, 6))    
        
        
    def test_get_lines(self):
        uu = np.array([[-3.7,-1.,0],[-3.7,-1.7,-1.7]])
        zz = np.array([[0.02,0.01,0.01],[0.02,0.03,0.07]])
        nh = np.array([[10,200,200],[5,5000,5000]])   
        el = ph.get_lines(uu,zz,self.filenom,lnH=nh,verbose=False)
        # Check that the output has the correct shape
        self.assertEqual(el.shape, (2, 18, 3))    
        # Check that the boundaries are handled as expected
        np.testing.assert_allclose(el[0,:,0],el[1,:,0])
        np.testing.assert_allclose(el[0,:,1],el[0,:,2])
        np.testing.assert_allclose(el[1,:,1],el[1,:,2])
    
        uu = np.array([[-3.7,-1.,0],[-3.7,-1.7,-1.7]])
        zz = np.array([[0.02,0.01,0.01],[0.02,0.07,1.]])
        nh = np.array([[100,200,200],[5,2000,2000]])
        el = ph.get_lines(uu,zz,self.filenom,lnH=nh,
                          origin='NLR',photmod='feltre16',verbose=False)
        # Check that the output has the correct shape
        self.assertEqual(el.shape, (2, 20, 3))    
        # Check that the boundaries are handled as expected
        np.testing.assert_allclose(el[0,:,0],el[1,:,0])
        np.testing.assert_allclose(el[0,:,1],el[0,:,2])
        np.testing.assert_allclose(el[1,:,1],el[1,:,2])

        
if __name__ == '__main__':
    unittest.main()
