#python3 -m unittest tests/test_stats.py 

import unittest
import numpy as np

import src.gne_stats as st

class TestPredict(unittest.TestCase):
    def test_locate_interval(self):
        edges = np.array([0,1,2,3,4])
        self.assertEqual(st.locate_interval(0.5,edges),0)
        vals = st.locate_interval(np.array([0,-0.2,3.5,4,5]),edges)
        i = 0
        for vexp in [0,-1,3,4,4]:
            self.assertEqual(vals[i],vexp); i += 1

        edges = np.array([-3,-2.,-1.])
        self.assertEqual(st.locate_interval(-2.5,edges),0)
        vals = st.locate_interval(np.array([-3,-3.2,-1,0]),edges)
        i = 0
        for vexp in [0,-1,2,2]:
            self.assertEqual(vals[i],vexp); i += 1

        
    def test_interpl_weights(self):
        edges = np.array([0,1,2])
        # Test floats
        xd, ix = st.interpl_weights(0.5, edges)
        self.assertAlmostEqual(xd, 0.5, places=7)
        self.assertEqual(ix, 0)
        # Test arrays
        x = np.array([0.5, 0, 2.0, 1., 2.2, 1.2])
        xd, ix = st.interpl_weights(x, edges)
        expected_xd = [0.5, 0., 1., 0., 1., 0.2]
        expected_ix = [  0, 0,  1,  1,  1,  1]
        for i in range(x.size):
            with self.subTest(i=i, x=x[i]):
                self.assertAlmostEqual(xd[i], expected_xd[i])
                self.assertEqual(ix[i], expected_ix[i])

                
if __name__ == '__main__':
    unittest.main()
