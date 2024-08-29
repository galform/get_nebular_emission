#python3 -m unittest tests/test_io.py 

import shutil
import unittest
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
        incols = [True, 6, 11, 19, 12]
        self.assertEqual(io.get_data_agnnH(txtfile,incols,inputformat='txt'),
                         [6, 11, 19, 12])

        incols = [True, 'data/mgas_disk','data/rhm_disk',
          'data/mgas_bulge','data/rhm_bulge']
        self.assertEqual(io.get_data_agnnH(hf5file,incols),
                         ['data/mgas_disk','data/rhm_disk',
          'data/mgas_bulge','data/rhm_bulge'])
        
if __name__ == '__main__':
    unittest.main()
