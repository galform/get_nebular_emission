#python3 -m unittest tests/test_io.py 

import shutil
import unittest
from numpy.testing import assert_allclose

import src.gne_io as io

txtfile = 'data/example_data/iz61/GP20_31p25kpc_z0_example_vol0.txt'
hf5file = 'data/example_data/iz61/GP20_31p25kpc_z0_example_vol0.hdf5'

class TestStringMethods(unittest.TestCase):

    def test_outroot(self):
        expath = 'output/a'
        self.assertEqual(io.get_outroot('a/example.txt',100),
                         'output/iz100/example')
        self.assertEqual(io.get_outroot('example.txt',32,outpath=expath),
                         expath+'/iz32/example')
        shutil.rmtree(expath)

    def test_plotpath(self):
        expath = 'plots/'
        self.assertEqual(io.get_plotpath('root'),expath)
        shutil.rmtree(expath)
        expath = 'output/a/plots/'
        self.assertEqual(io.get_plotpath('output/a/root'),expath)
        self.assertEqual(io.get_plotpath('output/a/'),expath)
        shutil.rmtree(expath)
        
    #def test_outnom(self):
        #self.assertEqual(io.get_outnom('a/example.txt',100),
        #                 'output/iz100/example.hdf5')
        #self.assertEqual(io.get_outnom('example.txt',39,ftype='plots'),
        #                 'output/iz39/plots/bpt_example.pdf')

    def test_get_data_agnnH(self):
        incols = [6, 11]
        sel=[0]
        assert_allclose(io.get_data_agnnH(txtfile,'rscale',incols,
                                           selection=sel,inputformat='txt'),
                         [[6.69049152e+08],[3.02573503e-03]],rtol=0.01)  
        assert_allclose(io.get_data_agnnH(txtfile,'reff',incols,
                                           selection=sel,inputformat='txt'),
                         [[6.69049152e+08],[1.80317940e-03]],rtol=0.01)  
        assert_allclose(io.get_data_agnnH(txtfile,'r',incols,
                                           selection=sel,inputformat='txt'),
                         [[6.69049152e+08],[9.01589699e-04]],rtol=0.01)  

        incols = ['data/mgas_disk','data/rhm_disk',
                  'data/mgas_bulge','data/rhm_bulge']
        assert_allclose(io.get_data_agnnH(hf5file,'rscale',incols,
                                           selection=sel),
                         [[6.69049152e+08],[3.02573503e-03],
                          [0.00000000e+00],[0.00000000e+00]],rtol=0.01)
        
if __name__ == '__main__':
    unittest.main()
