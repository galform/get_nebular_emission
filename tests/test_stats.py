# python -m unittest tests/test_stats.py 

import unittest
import numpy as np

import gne.gne_stats as st
import gne.gne_const as c

def func(x):
    return x

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

    def test_ensure_2d(self):
        a = np.zeros(3)
        b = st.ensure_2d(a)
        self.assertEqual(np.shape(b),(1,3))        
        b = st.ensure_2d(a,axis=1)
        self.assertEqual(np.shape(b),(3,1))        
        a = np.zeros((2,3))
        b = st.ensure_2d(a)
        self.assertEqual(np.shape(b),(2,3))        

        
    def test_components2tot(self):
        # Test data with one components
        ncomp = 1; xx = np.zeros((3,ncomp))
        xx[0,0]=10;  xx[1,0]=11.2; xx[2,0]=12
        vals = st.components2tot(xx)
        np.testing.assert_allclose(vals,xx, atol=0.001)

        # Generate data with several components
        ncomp = 2; xx = np.zeros((3,ncomp))
        xx[0,0]=10;  xx[1,0]=11.2; xx[2,0]=12
        xx[0,1]=10.5;xx[1,1]=11.8; xx[2,1]=c.notnum

        # Tests with and without log10input
        expected = np.array([10.619,11.897,12.])
        vals = st.components2tot(xx)
        np.testing.assert_allclose(vals,expected, atol=0.001)
        
        expected = np.array([20.5,23.,12.])
        vals = st.components2tot(xx, log10input=False)
        np.testing.assert_allclose(vals,expected, atol=0.001)

    def test_n_gt_x(self):
        # Simple 1D array
        xedges = np.array([0, 1, 2, 3])
        array = np.array([0, 1, 2, 3, 4])
        expected = np.array([4, 3, 2, 1])
        vals = st.n_gt_x(xedges, array)
        np.testing.assert_allclose(vals, expected)
        
        # All elements above all edges
        xedges = np.array([-2, -1])
        array = np.array([0, 1, 2])
        expected = np.array([3, 3])
        vals = st.n_gt_x(xedges, array)
        np.testing.assert_allclose(vals, expected)
        
        # No elements above any edge
        xedges = np.array([10, 20])
        array = np.array([0, 1, 2, 3])
        expected = np.array([0, 0])
        vals = st.n_gt_x(xedges, array)
        np.testing.assert_allclose(vals, expected)
                
        # Array with duplicates
        xedges = np.array([1, 2])
        array = np.array([1, 1, 2, 2, 3])
        expected = np.array([3, 1])
        vals = st.n_gt_x(xedges, array)
        np.testing.assert_allclose(vals, expected)
    
        
    def test_safe_sum(self):
        arrays = [np.array([1,2,c.notnum]),np.array([2,3,4])]
        expected = np.array([3,5,4])
        vals = st.safe_sum_arrays(arrays)
        np.testing.assert_allclose(vals,expected, atol=0.001)

        arrays = [np.array([1,2,c.notnum]),np.array([2,3,c.notnum])]
        expected = np.array([3,5,c.notnum])
        vals = st.safe_sum_arrays(arrays)
        np.testing.assert_allclose(vals,expected, atol=0.001)        
                
    def test_romberg(self):
        expected = 0.5
        val = st.romberg(func,0,1)
        self.assertAlmostEqual(val, expected)        


    def test_vol_sphere(self):
        expected = (4/3)*np.pi
        val = st.vol_sphere(1)
        self.assertEqual(val,expected)

        n=3
        expected = np.zeros(n); expected.fill((4/3)*np.pi)
        r = np.zeros(n); r.fill(1)
        val = st.vol_sphere(r)
        np.testing.assert_allclose(val,expected, atol=0.001)

        
if __name__ == '__main__':
    unittest.main()
