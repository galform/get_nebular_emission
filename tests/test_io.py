#python3 -m unittest tests/test_io.py 

import shutil
import unittest
import src.gne_io as io

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

#    def test_isupper(self):
#        self.assertTrue('FOO'.isupper())
#        self.assertFalse('Foo'.isupper())
#
#    def test_split(self):
#        s = 'hello world'
#        self.assertEqual(s.split(), ['hello', 'world'])
#        # check that s.split fails when the separator is not a string
#        with self.assertRaises(TypeError):
#            s.split(2)

if __name__ == '__main__':
    unittest.main()
