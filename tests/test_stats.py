#python3 -m unittest tests/test_stats.py 

from unittest import TestCase
import numpy as np

import src.gne_stats as st

class TestPredict(TestCase):
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

        
if __name__ == '__main__':
    unittest.main()
