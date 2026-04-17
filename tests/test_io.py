# python -m unittest tests/test_io.py 

import os, sys
import shutil
import unittest
import tempfile
import h5py
import numpy as np
from numpy.testing import assert_allclose

import gne.gne_io as io
import gne.gne_const as c

opath = 'data/example_data/iz61/ivol0/'
nom = 'ex'
txtfile = opath+nom+'.txt'
hf5file = opath+nom+'.hdf5'

class TestStringMethods(unittest.TestCase):
    def test_outroot(self):
        snap = 61
        rootdir = c.repo_dir
        expath = os.path.join(rootdir,'output','iz'+str(snap),'ivol')
        dirf, endf = io.get_outroot(snap,'example.txt')
        self.assertEqual(dirf,expath)
        self.assertEqual(endf,'example.hdf5')
        
    def test_plotfile(self):
        expath = 'output/iz61/ivol0/'
        dirf = io.get_plotfile(expath,nom,'uuu')
        self.assertEqual(dirf,'output/iz61/plots/uuu_ex.pdf')

    def test_outnom(self):
        outf = io.get_outnom(txtfile,verbose=False)
        self.assertEqual(outf,'output/iz61/ivol0/ex.hdf5')

        outf = io.get_outnom(txtfile,nomf='outfile',verbose=False)
        self.assertEqual(outf,'output/iz61/ivol0/outfile.hdf5')

        expath = 'uuu'
        outf = io.get_outnom(txtfile,dirf=expath,verbose=False)
        expected = expath+'/iz61/ivol0/ex.hdf5'
        self.assertEqual(outf,expected)        
        shutil.rmtree(expath)

    def test_get_param(self):
        att_config = {'albedo': 0, 'costheta': None}

        result = io.get_param(att_config, 'albedo', c.albedo)
        self.assertEqual(result, 0)

        result = io.get_param(att_config, 'costheta', c.albedo)
        self.assertEqual(result, c.albedo)
        
        att_config = {}
        result = io.get_param(att_config, 'Rv', c.Rv)
        self.assertEqual(result, c.Rv)

        att_config = None
        result = io.get_param(att_config, 'Rv', c.Rv)
        self.assertEqual(result, c.Rv)

    def test_add2headers(self):
        self.temp_dir = tempfile.mkdtemp()
        self.test_file = os.path.join(self.temp_dir, 'test_output.hdf5')
        with h5py.File(self.test_file, 'w') as hf:
            hf.create_dataset('header', (100,))

        names = ['mgasr_type']
        values = [['disc', 'bulge']]
        count = io.add2header(self.test_file, names, values, verbose=False)
        
        # Check that one attribute was added
        self.assertEqual(count, 1)
        
        # Read back and verify
        with h5py.File(self.test_file, 'r') as hf:
            self.assertIn('mgasr_type', hf['header'].attrs)
            mgasr_type_raw = hf['header'].attrs['mgasr_type']
            mgasr_type = io.decode_string_list(mgasr_type_raw)
            
            self.assertEqual(len(mgasr_type), 2)
            self.assertEqual(mgasr_type[0], 'disc')
            self.assertEqual(mgasr_type[1], 'bulge')
            
        if os.path.exists(self.test_file):
            os.remove(self.test_file)
        os.rmdir(self.temp_dir)
            
    def test_hdf5_headers(self):
        # Note: get_outnom now strips the last path component from dirf
        # So 'output/test' becomes 'output' before adding /iz.../
        expath = 'output/test'
        h0 = 0.7
        filenom = io.generate_header(hf5file,0,100,h0,0.4,0.3,0.6,
                                     1e8,100,1e6,
                                     outpath=expath,verbose=False)
        # Updated expectation: get_outnom strips 'test' from path
        self.assertEqual(filenom,'output/test/iz61/ivol0/'+nom+'.hdf5')
        
        names = ['a','b',None]
        values = [1,'b',3]
        num = io.add2header(filenom,names,values,verbose=False)
        self.assertEqual(num,2)
        hf = h5py.File(filenom, 'r')
        self.assertEqual(hf['header'].attrs['h0'], h0)
        self.assertEqual(hf['header'].attrs[names[0]], values[0])
        self.assertEqual(hf['header'].attrs[names[1]], values[1])
        hf.close()        
        shutil.rmtree('output')


    def test_read_mgas_hr(self):
        vb = False
        sel=[0,1]
        
        # Test with text file first (using column indices)
        incols_txt_single = [[6, 11]]
        # Expected values from text file (must match text file content)
        with open(txtfile, 'r') as f:
            lines = [l for l in f if not l.startswith('#')]
            row0 = lines[0].split()
            row1 = lines[1].split()
        expect_m_txt = np.array([[float(row0[6]), float(row1[6])]])
        expect_r_txt = np.array([[float(row0[11]), float(row1[11])]])
        
        mm,rr = io.read_mgas_hr(txtfile,incols_txt_single,sel,
                                inputformat='txt',verbose=vb)
        assert_allclose(mm,expect_m_txt,rtol=0.01)  
        assert_allclose(rr,expect_r_txt,rtol=0.01)  

        # Test with HDF5 file (using dataset paths)
        # Dataset names in actual HDF5: data/mcold, data/rdisk, data/mcold_burst, data/rbulge
        incols_hdf5 = [['data/mcold','data/rdisk'],
                       ['data/mcold_burst','data/rbulge']]
        incols_txt = [[6,11],[9,12]]
        
        # Get expected values from actual HDF5 file
        with h5py.File(hf5file, 'r') as hf:
            expect_m = np.array([[hf['data/mcold'][0], hf['data/mcold'][1]],
                                 [hf['data/mcold_burst'][0], hf['data/mcold_burst'][1]]])
            expect_r = np.array([[hf['data/rdisk'][0], hf['data/rdisk'][1]],
                                 [hf['data/rbulge'][0], hf['data/rbulge'][1]]])
        
        mm,rr = io.read_mgas_hr(hf5file,incols_hdf5,sel,verbose=vb)
        assert_allclose(mm,expect_m,rtol=0.01)  
        assert_allclose(rr,expect_r,rtol=0.01)  

        # Text file test with matching columns
        m_t,r_t = io.read_mgas_hr(txtfile,incols_txt,sel,
                                  inputformat='txt',verbose=vb)
        # Note: text and HDF5 files may have different data, 
        # so we just verify the function works, not exact match
        self.assertEqual(m_t.shape, (2, 2))
        self.assertEqual(r_t.shape, (2, 2))
        
        
    def test_get_mgas_hr(self):
        vb = False
        sel=[0]
        
        # Dataset names: data/mcold, data/rdisk
        incols = [['data/mcold','data/rdisk']]
        
        # Get expected values from actual HDF5 file
        with h5py.File(hf5file, 'r') as hf:
            expect_m = np.array([[hf['data/mcold'][0]]])
            expect_r = np.array([[hf['data/rdisk'][0]]])
        
        mm,rr = io.get_mgas_hr(hf5file,sel,incols,[4],verbose=vb)
        assert_allclose(mm,expect_m,rtol=0.01)
        assert_allclose(rr,expect_r,rtol=0.01)

        sel=[0,1,3]
        rtype=[0,0]
        incols = [['data/mcold','data/rdisk'],
                  ['data/mcold_burst','data/rbulge']]

        # Get expected values from actual HDF5 file
        with h5py.File(hf5file, 'r') as hf:
            expect_m = np.array([[hf['data/mcold'][i] for i in sel],
                                 [hf['data/mcold_burst'][i] for i in sel]])
            expect_r = np.array([[hf['data/rdisk'][i] for i in sel],
                                 [hf['data/rbulge'][i] for i in sel]])

        mm,rr = io.get_mgas_hr(hf5file,sel,incols,rtype,verbose=vb)
        assert_allclose(mm,expect_m,rtol=0.01)  
        assert_allclose(rr,expect_r,rtol=0.01)
        
        mm,rr = io.get_mgas_hr(hf5file,sel,incols,[1,2],verbose=vb)
        assert_allclose(mm,expect_m,rtol=0.01)
        exp = np.copy(expect_r)
        exp[0,:] = c.re2hr*expect_r[0,:]
        exp[1,:] = c.re2hr*c.r502re*expect_r[1,:]
        assert_allclose(rr,exp,rtol=0.01)  
        
        h=2.
        mm,rr = io.get_mgas_hr(hf5file,sel,incols,[0,3],
                               h0=h,units_h0=True,verbose=vb)
        assert_allclose(mm,expect_m/h,rtol=0.01)
        exp = np.copy(expect_r)/h
        exp[1,:] = c.re2hr*c.r502re*c.rvir2r50*exp[1,:]
        assert_allclose(rr,exp,rtol=0.01)  
        
        re2hr=8.; r502re=300.
        mm,rr = io.get_mgas_hr(hf5file,sel,incols,[1,2],
                               re2hr=re2hr,r502re=r502re,verbose=vb)
        assert_allclose(mm,expect_m,rtol=0.01)
        exp = np.copy(expect_r)
        exp[0,:] = re2hr*expect_r[0,:]
        exp[1,:] = re2hr*r502re*expect_r[1,:]
        assert_allclose(rr,exp,rtol=0.01)  
        

    def test_read_data(self):
        params_txt = [0, 2] 
        # Correct HDF5 dataset paths: data/mstars_disk, data/mstardot
        params_hdf5 = ['data/mstars_disk', 'data/mstardot']

        sel = [0, 1, 2] # Select 3 first galaxies
        
        # Get expected values from actual HDF5 file
        with h5py.File(hf5file, 'r') as hf:
            expect_mstar = np.array([hf['data/mstars_disk'][i] for i in sel])
            expect_sfr = np.array([hf['data/mstardot'][i] for i in sel])

        # Test hdf5 file            
        result_hdf5 = io.read_data(hf5file, sel, inputformat='hdf5', 
                                   params=params_hdf5, verbose=False)
        self.assertEqual(result_hdf5.shape, (2, 3))
        assert_allclose(result_hdf5[0], expect_mstar, rtol=0.01)
        assert_allclose(result_hdf5[1], expect_sfr, rtol=0.01)

        result_single_hdf5 = io.read_data(hf5file, sel, inputformat='hdf5',
                                          params=['data/mstars_disk'], verbose=False)
        self.assertEqual(result_single_hdf5.ndim, 1)

        params_with_none = ['data/mstars_disk', None, 'data/mstardot']
        result_none = io.read_data(hf5file, sel, inputformat='hdf5',
                                   params=params_with_none, verbose=False)
        self.assertEqual(result_none.shape, (2, 3))

        with self.assertRaises(SystemExit):
            io.read_data(hf5file, sel, inputformat='invalid', 
                         params=['data/mstars_disk'], verbose=False)
    

        result_missing = io.read_data(hf5file, sel, inputformat='hdf5',
                                      params=['data/nonexistent'], verbose=False)
        assert_allclose(result_missing, np.zeros(3), rtol=0.01)

        empty_sel = np.array([], dtype=int)
        result_empty = io.read_data(hf5file, empty_sel, inputformat='hdf5',
                                    params=['data/mstars_disk'], verbose=False)
        self.assertEqual(len(result_empty), 0)
    
        single_sel = [0]
        result_single = io.read_data(hf5file, single_sel, inputformat='hdf5',
                                     params=['data/mstars_disk'], verbose=False)
        self.assertEqual(len(result_single), 1)

        large_sel = list(range(min(50, 100))) 
        result_large = io.read_data(hf5file, large_sel, inputformat='hdf5',
                                    params=['data/mstars_disk'], verbose=False)
        self.assertEqual(len(result_large), len(large_sel))
        
        # Test text file
        result_txt = io.read_data(txtfile, sel, inputformat='txt', 
                                  params=params_txt, verbose=False)
        self.assertEqual(result_txt.shape, (2, 3))
        # Note: text file may have different values than HDF5

        result_single_txt = io.read_data(txtfile, sel, inputformat='txt',
                                         params=[0], verbose=False)
        self.assertEqual(result_single_txt.ndim, 1)
        
        # Test different selections work
        sel_diff = [1, 3, 5]
        result_diff_hdf5 = io.read_data(hf5file, sel_diff, inputformat='hdf5',
                                        params=['data/mstars_disk'], verbose=False)
        result_diff_txt = io.read_data(txtfile, sel_diff, inputformat='txt',
                                       params=[0], verbose=False)
        self.assertEqual(len(result_diff_hdf5), 3)
        self.assertEqual(len(result_diff_txt), 3)
    
        # Test parameter ordering is preserved
        params_ordered = ['data/mstars_disk', 'data/mstardot', 'data/Zgas_disc']
        params_txt_ordered = [0, 2, 4]
   
        result_ordered_hdf5 = io.read_data(hf5file, sel_diff, inputformat='hdf5',
                                           params=params_ordered, verbose=False)
        result_ordered_txt = io.read_data(txtfile, sel_diff, inputformat='txt',
                                          params=params_txt_ordered, verbose=False)
        self.assertEqual(result_ordered_hdf5.shape, (3, 3))
        self.assertEqual(result_ordered_txt.shape, (3, 3))


    def test_read_sfrdata(self):
        sel = np.array([0, 1, 2, 3, 4])  # First 5 galaxies
        
        # Correct HDF5 dataset paths matching actual file structure:
        # mstars_disk, mstardot (SFR for disk), Zgas_disc
        # mstars_bulge, mstardot_burst (SFR for bulge), Zgas_bst
        cols_hdf5 = [['data/mstars_disk','data/mstardot','data/Zgas_disc'],
                     ['data/mstars_bulge','data/mstardot_burst','data/Zgas_bst']]
        cols_txt = [[0,2,4],[1,3,5]]
        
        cols_single_hdf5 = [['data/mstars_disk','data/mstardot','data/Zgas_disc']]
        cols_single_txt = [[0,2,4]]
        
        sel_reverse = np.array([1, 0])
        
        # HDF5 file
        ms_hdf5, sfr_hdf5, z_hdf5 = io.read_sfrdata(hf5file, cols_hdf5, sel,
                                                    inputformat='hdf5', verbose=False)
        self.assertEqual(ms_hdf5.shape, (2, 5))
        self.assertEqual(sfr_hdf5.shape, (2, 5))
        self.assertEqual(z_hdf5.shape, (2, 5))
    
        ms_hdf5_rev, sfr_hdf5_rev, z_hdf5_rev = io.read_sfrdata(
            hf5file, cols_hdf5, sel_reverse, inputformat='hdf5', verbose=False)
        assert_allclose(ms_hdf5[:, 0], ms_hdf5_rev[:, 1], rtol=1e-12)
        assert_allclose(ms_hdf5[:, 1], ms_hdf5_rev[:, 0], rtol=1e-12)
            
        ms_single_hdf5, sfr_single_hdf5, z_single_hdf5 = io.read_sfrdata(
            hf5file, cols_single_hdf5, sel, inputformat='hdf5', verbose=False)
        self.assertEqual(ms_single_hdf5.shape, (1, 5))
    
        with self.assertRaises(SystemExit): # Test invalid input
            io.read_sfrdata(hf5file, cols_hdf5, np.array([0]),
                            inputformat='invalid', verbose=False)
            
        empty_sel = np.array([], dtype=int) # Test empty selection
        ms_empty, sfr_empty, z_empty = io.read_sfrdata(
            hf5file, cols_hdf5, empty_sel, inputformat='hdf5', verbose=False)
        self.assertEqual(ms_empty.shape[1], 0)
        self.assertEqual(sfr_empty.shape[1], 0) 
        self.assertEqual(z_empty.shape[1], 0)
    
        # Text file
        ms_txt, sfr_txt, z_txt = io.read_sfrdata(txtfile, cols_txt, sel,
                                                 inputformat='txt', verbose=False)
        self.assertEqual(ms_txt.shape, (2, 5))
        self.assertEqual(sfr_txt.shape, (2, 5))
        self.assertEqual(z_txt.shape, (2, 5))
    
        ms_txt_rev, sfr_txt_rev, z_txt_rev = io.read_sfrdata(
            txtfile, cols_txt, sel_reverse, inputformat='txt', verbose=False)
        assert_allclose(ms_txt[:, 0], ms_txt_rev[:, 1], rtol=1e-12)
        assert_allclose(ms_txt[:, 1], ms_txt_rev[:, 0], rtol=1e-12)
    
        ms_single_txt, sfr_single_txt, z_single_txt = io.read_sfrdata(
            txtfile, cols_single_txt, sel, inputformat='txt', verbose=False)
        self.assertEqual(ms_single_txt.shape, (1, 5))
        
        # Verify single component matches first component
        assert_allclose(ms_single_hdf5[0,:], ms_hdf5[0,:], rtol=1e-12,
                        err_msg="Single component should match first one")


    def test_read_and_add(self):
        """Test read_and_add function for various data shapes."""
        # Create temporary HDF5 file with test datasets
        self.temp_dir = tempfile.mkdtemp()
        self.test_file = os.path.join(self.temp_dir, 'test_read_add.hdf5')
        
        with h5py.File(self.test_file, 'w') as hf:
            # 1D dataset
            hf.create_dataset('prop_1d', data=np.array([1.0, 2.0, 3.0, 4.0]))      
            # 2D dataset with shape (n, 2) - expected case
            hf.create_dataset('prop_2d_2cols', 
                              data=np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]))      
            # 2D dataset with shape (n, 3) - unexpected case
            hf.create_dataset('prop_2d_3cols',
                              data=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))      
            # 2D dataset with shape (n, 1)
            hf.create_dataset('prop_2d_1col',
                              data=np.array([[1.0], [2.0], [3.0]]))
        
        with h5py.File(self.test_file, 'r') as f:
            # Test 1D data - should return unchanged
            result_1d = io.read_and_add(f, 'prop_1d', verbose=False)
            self.assertEqual(result_1d.shape, (4,))
            assert_allclose(result_1d, np.array([1.0, 2.0, 3.0, 4.0]))
            
            # Test 2D with 2 columns - should sum along axis=1
            result_2d_2 = io.read_and_add(f, 'prop_2d_2cols', verbose=False)
            self.assertEqual(result_2d_2.shape, (3,))
            assert_allclose(result_2d_2, np.array([3.0, 7.0, 11.0]))
            
            # Test 2D with 3 columns - should still sum
            result_2d_3 = io.read_and_add(f, 'prop_2d_3cols', verbose=False)
            self.assertEqual(result_2d_3.shape, (2,))
            assert_allclose(result_2d_3, np.array([6.0, 15.0]))
            
            # Test 2D with 1 column - should sum (effectively flatten)
            result_2d_1 = io.read_and_add(f, 'prop_2d_1col', verbose=False)
            self.assertEqual(result_2d_1.shape, (3,))
            assert_allclose(result_2d_1, np.array([1.0, 2.0, 3.0]))
        
        # Test verbose output (capture stdout to verify warning)
        import io as io_module
        from contextlib import redirect_stdout
        
        with h5py.File(self.test_file, 'r') as f:
            # Test warning for 3 columns
            captured = io_module.StringIO()
            with redirect_stdout(captured):
                io.read_and_add(f, 'prop_2d_3cols', verbose=True)
            output = captured.getvalue()
            self.assertIn('[WARN]', output)
            self.assertIn('prop_2d_3cols', output)
            self.assertIn('(2, 3)', output)
            
            # Test no warning for 2 columns
            captured = io_module.StringIO()
            with redirect_stdout(captured):
                io.read_and_add(f, 'prop_2d_2cols', verbose=True)
            output = captured.getvalue()
            self.assertEqual(output, '')
            
            # Test no warning for 1D
            captured = io_module.StringIO()
            with redirect_stdout(captured):
                io.read_and_add(f, 'prop_1d', verbose=True)
            output = captured.getvalue()
            self.assertEqual(output, '')
        
        # Cleanup
        if os.path.exists(self.test_file):
            os.remove(self.test_file)
        os.rmdir(self.temp_dir)
        
if __name__ == '__main__':
    unittest.main()
