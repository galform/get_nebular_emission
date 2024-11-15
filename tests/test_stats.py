#python3 -m unittest tests/test_stats.py 

import unittest
import numpy as np

import src.gne_stats as st

class TestPredict(unittest.TestCase):
    def test_locate_interval(self):
        edges = np.array([0,1,2,3,4])
        self.assertEqual(st.locate_interval(0.5,edges),0)
        self.assertEqual(st.locate_interval(0,edges),-1)
        self.assertEqual(st.locate_interval(-0.2,edges),-1)
        self.assertEqual(st.locate_interval(4,edges),3)
        self.assertEqual(st.locate_interval(4.5,edges),4)
        edges = np.array([-3,-2.,-1.])
        self.assertEqual(st.locate_interval(-2.5,edges),0)
        self.assertEqual(st.locate_interval(-3,edges),-1)
        self.assertEqual(st.locate_interval(-3.2,edges),-1)
        self.assertEqual(st.locate_interval(-1,edges),1)
        self.assertEqual(st.locate_interval(0,edges),2)
        
    def test_interpl_weights(self):
        edges = np.array([0,1,2])
        # Test floats
        test_cases = [ # x value, expected (xd, ix)
            (0.5, (0.5, 0)),  # middle of first interval
            (0, (0.0, 0)),      # lower edge
            (2., (1., 1)),       # upper edge
            (-0.2, (0.0, 0)),   # below range
            (2.2, (1., 1)),     # above range
        ]
        for x, (expected_xd, expected_ix) in test_cases:
            with self.subTest(x=x):
                xd, ix = st.interpl_weights(x, edges)
                self.assertAlmostEqual(xd, expected_xd, places=7)
                self.assertEqual(ix, expected_ix)
                
        # Test arrays
        x = np.array([0.5, 0, 2.0, -0.2, 2.2, 1.5])
        xd, ix = st.interpl_weights(x, edges)
        expected_xd = [0.5, 0., 1., 0., 1., 0.5]
        expected_ix = [0, 0, 1, 0, 1, 1]
        for i in range(x.size):
            with self.subTest(i=i, x=x[i]):
                self.assertAlmostEqual(xd[i], expected_xd[i])
                self.assertEqual(ix[i], expected_ix[i])

                
if __name__ == '__main__':
    unittest.main()
