# python -m unittest tests/test_slurm.py 

import unittest
import os
import tempfile
import shutil
import sys
from unittest.mock import patch

import gne.gne_const as c
import gne.gne_slurm as su

class TestSlurmUtilsHdf5(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures before running tests in the class"""
        # Create a temporary directory in the system temp location
        cls.test_dir = tempfile.mkdtemp(prefix='gne_test_')
        
        # Create a mock SLURM template file for testing
        cls.template_content = '''#!/bin/sh
#SBATCH --job-name=__GNE_JOB_NAME__
#SBATCH --error=__GNE_LOG_DIR__/%x_ivol%a_%A.err 
#SBATCH --output=__GNE_LOG_DIR__/%x_ivol%a_%A.out
#SBATCH --array=__GNE_VOLS__%30 
#
export SLURM_ID=$SLURM_ARRAY_JOB_ID 
export GNE_SUBVOL_INDEX=$SLURM_ARRAY_TASK_ID
srun python << 'EOF_PYTHON_SCRIPT'
__GNE_PARAM_CONTENT__
EOF_PYTHON_SCRIPT
'''
        # Create a mock parameter file for testing
        cls.param_content = '''import gne.gne_const as const
from gne.gne import gne

verbose = True
testing = False
get_emission_lines = True

outpath = '/home/user/Data/Galform/SU1/'
subvols = 2
root = outpath+'iz87/ivol'
endf = 'gne_input.hdf5'
'''
        # Create a REAL template file in the temp directory
        template_file = os.path.join(cls.test_dir, 'slurm_taurus_template.sh')
        with open(template_file, 'w') as f:
            f.write(cls.template_content)

        # Redirect slurm_temp_dir to temp directory
        cls.dir_patcher = patch.object(c, 'slurm_temp_dir', cls.test_dir)
        cls.dir_patcher.start()

        # Create test parameter file
        cls.param_file_path = os.path.join(cls.test_dir, 'run_gne_TestSim.py')
        with open(cls.param_file_path, 'w') as f:
            f.write(cls.param_content)

        
    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests in the class are finished"""
        # Remove temporary test directory
        if os.path.exists(cls.test_dir):
            try:
                shutil.rmtree(cls.test_dir)
            except (OSError, PermissionError) as e:
                print("Warning: Could not remove test directory {}: {}".format(cls.test_dir, e))
        
        # Stop patcher
        cls.dir_patcher.stop()

    def test_generate_job_name(self):
        job_name = su.generate_job_name(model='sam',snap=87)
        self.assertEqual(job_name, 'gne_sam_iz87')

        job_name = su.generate_job_name(model='sam',snap=87,
                                        job_suffix='custom')
        self.assertEqual(job_name, 'gne_sam_iz87_custom')

    def test_get_slurm_template(self):
        template = su.get_slurm_template('taurus')
        self.assertIsNotNone(template)
        self.assertIn('__GNE_JOB_NAME__', template)
        self.assertIn('__GNE_PARAM_CONTENT__', template)
        self.assertIn('__GNE_VOLS__', template)

        with self.assertRaises(SystemExit):
            su.get_slurm_template('nonexistent_hpc')

    def test_modify_param_file(self):
        """Test modification of parameter file content"""
        modified = su.modify_param_file(self.param_file_path,'TestSim','128')

        self.assertIn('TestSim', modified)
        self.assertIn('iz128', modified)
        self.assertNotIn('iz87', modified)


    def test_create_slurm_script(self):
        verbose = False
        script_path = su.create_slurm_script('taurus',
                                             self.param_file_path,
                                             'TestSim','sam','96','5,46',
                                             logdir=self.test_dir,
                                             verbose=verbose)
        # Check returned values
        self.assertIn('gne_sam_iz96',script_path)
        
        # Check file was created
        self.assertTrue(os.path.exists(script_path))
        
        # Variables should have been replaced
        with open(script_path, 'r') as f:
            content = f.read()
        self.assertNotIn('__GNE_LOG_DIR__', content)
        self.assertNotIn('__GNE_JOB_NAME__', content)
        self.assertNotIn('__GNE_PARAM_CONTENT__', content)
        self.assertNotIn('__GNE_VOLS__', content)


    def test_submit_slurm_job_no_sbatch(self):
        """Test submit_slurm_job when sbatch is not available"""
        # Create a dummy script file
        dummy_script = os.path.join(self.test_dir, 'dummy_script.sh')
        with open(dummy_script, 'w') as f:
            f.write('#!/bin/sh\necho "test"')
    
        with patch('shutil.which', return_value=None):
            with patch('subprocess.Popen', side_effect=FileNotFoundError(2, 'No such file or directory', 'sbatch')):
                try:
                    job_id = su.submit_slurm_job(dummy_script)
                    # Expect None when sbatch unavailable
                    self.assertIsNone(job_id)
                except FileNotFoundError:
                    self.fail("submit_slurm_job raised FileNotFoundError instead of returning None")


    def test_check_jobs_status(self):
        verbose = False

        logdir = os.path.join(self.test_dir, 'logs_check')
        os.makedirs(logdir, exist_ok=True)
        model = 'sam'
        snap = '00'
        
        files = ['gne_sam_iz00_ivol0_12345.err',
                 'gne_sam_iz00_ivol0_12345.out']
        ferr = os.path.join(logdir, files[0])
        fout = os.path.join(logdir, files[1])

         # Clean up any leftover files from previous tests
        for f in (ferr, fout):
            if os.path.exists(f):
                os.remove(f)       
        
        result = su.check_job_status(ferr, verbose=verbose)
        self.assertEqual(result,'not_found')

        with open(ferr, 'w') as f:
            pass
        with open(fout, 'w') as f:
            f.write('dummy content')
        result = su.check_job_status(ferr, verbose=verbose)
        self.assertEqual(result,'incomplete')

        with open(fout, 'w') as f:
            f.write('SUCCESS')
        result = su.check_job_status(ferr, verbose=verbose)
        self.assertEqual(result,'success')

        with open(ferr, 'w') as f:
            f.write('dummy content')
        result = su.check_job_status(ferr, verbose=verbose)
        self.assertEqual(result,'error')

                    
    def test_check_all_jobs(self):
        verbose = True
        logdir = os.path.join(self.test_dir, 'logs_check')
        os.makedirs(logdir, exist_ok=True)
        model = 'sam'
        snap = '00'

        files = ['gne_sam_iz00_ivol0.err',
                 'gne_sam_iz00_ivol0.out']
        ferr = os.path.join(logdir, files[0])
        fout = os.path.join(logdir, files[1])

        with open(ferr, 'w') as f:
            pass
        result = su.check_all_jobs(model,snap,logdir,verbose=verbose)
        self.assertEqual(len(result['not_found']),1)

        with open(fout, 'w') as f:
            f.write('dummy content')
        result = su.check_all_jobs(model,snap,logdir,verbose=verbose)
        self.assertEqual(len(result['incomplete']),1)

        with open(fout, 'w') as f:
            f.write('SUCCESS')
        result = su.check_all_jobs(model, snap, logdir, verbose=verbose)
        self.assertEqual(len(result['success']),1)

        with open(ferr, 'w') as f:
            f.write('dummy content')
        result = su.check_all_jobs(model, snap, logdir, verbose=verbose)
        self.assertEqual(len(result['error']),1)

        files = ['gne_sam_iz00_ivol0_suf.err',
                 'gne_sam_iz00_ivol0_suf.out']
        ferr = os.path.join(logdir, files[0])
        fout = os.path.join(logdir, files[1])
        with open(ferr, 'w') as f:
            pass
        with open(fout, 'w') as f:
            f.write('SUCCESS')
        result = su.check_all_jobs(model, snap, logdir,
                                   job_suffix='suf',verbose=verbose)
        self.assertEqual(len(result['success']),1)
        
                    
    def test_clean_all_jobs(self):
        verbose = False
        # Test only_show=True
        logdir = os.path.join(self.test_dir, 'logs_dry_run')
        os.makedirs(logdir, exist_ok=True)

        model = 'sam'
        snap = '96'
        files_to_create = [
            'gne_sam_iz96_ivol0_12345.err',
            'gne_sam_iz96_ivol0_12345.out',
            'gne_sam_iz96_ivol1_12345.err',
            'gne_sam_iz96_ivol1_12345.out',
            'submit_gne_sam_iz96_01_02_2026_10h_30m_00s000000.sh',
        ]
        for fname in files_to_create:
            with open(os.path.join(logdir, fname), 'w') as f:
                f.write('dummy content')
        deleted_files = su.clean_all_jobs(model, snap, logdir,
                                          only_show=True, verbose=verbose)
        self.assertEqual(len(deleted_files), 5)
        for fname in files_to_create:
            self.assertTrue(os.path.exists(os.path.join(logdir, fname)),
                            f'File {fname} was deleted during dry run')


        # Test only_show=False (actual deletion)
        deleted_files = su.clean_all_jobs(model, snap, logdir,
                                          only_show=False, verbose=verbose)
        self.assertEqual(len(deleted_files), 5)
        for fname in files_to_create:
            self.assertFalse(os.path.exists(os.path.join(logdir, fname)),
                             f'File {fname} was not deleted')

        # Test with job_suffix
        suffixed_files = [
            'gne_sam_iz96_custom_ivol0_12345.err',
            'gne_sam_iz96_custom_ivol0_12345.out',
            'gne_sam_iz96_custom_ivol0_12345.txt',
        ]
        for fname in suffixed_files:
            with open(os.path.join(logdir, fname), 'w') as f:
                f.write('dummy content')

        deleted_files = su.clean_all_jobs(
            model, snap, logdir,
            job_suffix='custom', only_show=True, verbose=verbose
        )
        self.assertEqual(len(deleted_files), 2)
        for f in deleted_files:
            self.assertIn('custom', os.path.basename(f))

        # Test non existing files
        deleted_files = su.clean_all_jobs('nomodel', '999', logdir,
                                          only_show=True, verbose=verbose)

        self.assertEqual(len(deleted_files), 0)
    

        
if __name__ == '__main__':
    unittest.main()
