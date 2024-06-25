#python3 -m unittest tests/test_io.py 

import src.gne_io as io

import unittest

class TestStringMethods(unittest.TestCase):

    def test_eq(self):
        self.assertEqual(io.get_outnom('a/example.txt',100),
                         'output/iz100/example.hdf5')
        self.assertEqual(io.get_outnom('example.txt',39,ftype='plots'),
                         'output/iz39/plots/bpt_example.pdf')

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
