#python3 -m unittest tests/test_stats.py 

import unittest
import numpy as np

import src.gne_stats as st

class TestPredict(unittest.TestCase):
    def test_get_cumulative_2Ddensity(self):
        # Generate a 2D Gaussian
        std1 = 1
        np.random.seed(42) 
        data = np.random.normal(0, std1, size=[10000, 2])
        x = data[:,0]; y = data[:,1]

        # Test for a given grid
        ng = 50
        xout, yout, zout = st.get_cumulative_2Ddensity(x,y, n_grid=ng)
        self.assertEqual(xout.size, ng*ng)
        self.assertEqual(yout.size, ng*ng)
        self.assertEqual(zout.size, ng*ng)

        r = np.sqrt(xout**2 + yout**2)
        mask = np.abs(r - 1.0) < 0.1
        values_at_1sigma = zout[mask]
        expected_value = 1 - np.exp(-0.5)  # ≈ 0.39 (2D Gaussian)
        mean_value = np.mean(values_at_1sigma)
        self.assertAlmostEqual(mean_value, expected_value, delta=0.05)
        
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
            self.assertAlmostEqual(xd[i], expected_xd[i])
            self.assertEqual(ix[i], expected_ix[i])

    def test_bilinear_interpl(self):
        xedges = np.array([14,15,16])
        yedges = np.array([20,21])
        zedges = np.array([[91,162],[210,95],[90,200]])
        zedges_3d = np.array([[[91,100],[162,170]],
                              [[210,220],[95,105]],
                              [[90,95],[200,210]]])
        # Test float and 2D zedges
        val = st.bilinear_interpl(14.5, 20.2, xedges, yedges, zedges)
        self.assertAlmostEqual(val, 146.1, places=7)
        # Test arrays with 2D zedges
        eval_x = np.array([15.5,  15,  14,  13,  16,  17,15.5,15.5,17,14.5,14.5])
        eval_y = np.array([20.8,20.8,20.2,20.2,20.2,20.2,  21,  22,22,  19,  20])
        expval= [148,118,105.2,105.2,112,112,147.5,147.5,200,150.5,150.5]
        val = st.bilinear_interpl(eval_x, eval_y, xedges, yedges, zedges)
        for i in range(eval_x.size):
            self.assertAlmostEqual(val[i], expval[i], places=7)
        # Test float and 3D zedges
        val = st.bilinear_interpl(14.5, 20.2, xedges, yedges, zedges_3d)
        np.testing.assert_array_almost_equal(val,[146.1, 155.5], decimal=7)
        # Test arrays with 3D zedges
        eval_x = np.array([14.5, 15.5])
        eval_y = np.array([20, 20.8])
        expval = [[150.5, 160],[148.0, 157.5]]
        val = st.bilinear_interpl(eval_x, eval_y, xedges, yedges, zedges_3d)
        np.testing.assert_array_almost_equal(val, expval, decimal=7)
        
if __name__ == '__main__':
    unittest.main()
