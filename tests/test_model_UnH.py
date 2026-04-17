#python -m unittest tests/test_model_UnH.py 

import shutil
import unittest
import os
import numpy as np
import gne.gne_io as io
import gne.gne_const as c
import gne.gne_model_UnH as unh


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
        cls.noms = ['model_U_NLR','model_spec_NLR','alpha_NLR',
                    'T_NLR_K','epsilon_NLR','nH_NLR_cm3']
        cls.nval = ['panuzzo03', c.model_spec_agn, c.alpha_NLR_feltre16, 5000, 1, 1]
        cls.num = io.add2header(cls.filenom, cls.noms, cls.nval, verbose=False)
   
    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests in the class are finished"""
        if os.path.exists(cls.expath):
            try:
                shutil.rmtree(cls.expath)
            except (OSError, PermissionError) as e:
                print(f"Warning: Could not remove test directory {cls.expath}: {e}")
    
    def test_get_alphaB(self):
        self.assertAlmostEqual(unh.get_alphaB(5000.),
                               4.54e-13,self.decimals)
        self.assertAlmostEqual(unh.get_alphaB(1000.),
                               4.54e-13,self.decimals)
        val = (1.43e-13+2.59e-13)/2.
        self.assertAlmostEqual(unh.get_alphaB(15000.),
                               val,self.decimals)
        self.assertAlmostEqual(unh.get_alphaB(30000.),
                               1.43e-13,self.decimals)  

    def test_get_Q_agn(self):
        Lagn = np.array([0.])
        alpha =c.alpha_NLR_feltre16
        expected = Lagn
        with self.assertRaises(SystemExit):
            unh.get_Q_agn(Lagn,alpha,model_spec="invalid",
                          verbose=False)
        vals = unh.get_Q_agn(Lagn,alpha)
        np.testing.assert_allclose(vals,expected,atol=self.atol)

        
    def test_get_U_panuzzo03(self):
        Q = np.array([1])
        alphaB = 4.54e-13
        uu = (3/(4*c.c_cm))*np.power(3*alphaB*alphaB/(4*np.pi),1/3)
        expected = np.log10(uu) - np.log10(3)
        vals = unh.get_U_panuzzo03(Q,self.filenom)
        np.testing.assert_allclose(vals,expected,atol=self.atol)
        
        vals = unh.get_U_panuzzo03(Q,self.filenom,epsilon=1,nH=1)
        np.testing.assert_allclose(vals,expected,atol=self.atol)

    def test_get_UnH_agn(self):
        Lagn = np.array([0.])
        val1, val2 = unh.get_UnH_agn(Lagn,None,None,self.filenom,
                                     verbose=False)
        self.assertEqual(val1[0], -999.)
        self.assertEqual(val2, None)

        Lagn = np.array([0.,1.])
        val1, val2 = unh.get_UnH_agn(Lagn,None,None,self.filenom,
                                     verbose=False)
        expected = np.array([[-999.,-17.63925603]])
        np.testing.assert_allclose(val1,expected,atol=self.atol)

    
if __name__ == '__main__':
    unittest.main()
