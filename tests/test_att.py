#python -m unittest tests/test_att.py 

import unittest
import numpy as np
from numpy.testing import assert_allclose

import gne.gne_const as c
import gne.gne_att as att

class TestPredict(unittest.TestCase):
    def test_get_f_saito20(self):
        expect = 1
        self.assertAlmostEqual(att.get_f_saito20(2.9),expect,2)
        for z in [1,2.8]:
            expect=0.44 + 0.2*z
            self.assertAlmostEqual(att.get_f_saito20(z),expect,2)

    def test_get_A_lambda(self):  
        tau = 1.0; costheta = 1.0; albedo = 0.0
        expected = -2.5 * np.log10((1.0 - np.exp(-1.0)) / 1.0)
        result = att.get_A_lambda(tau, costheta=costheta, albedo=albedo)       
        np.testing.assert_almost_equal(result, expected, decimal=10)

        tau = np.array([0.5, 1.0, 1.5, 2.0])
        costheta = 0.6
        albedo = 0.4        
        result = att.get_A_lambda(tau, costheta=costheta, albedo=albedo)
        # Verify output shape matches tau
        assert result.shape == tau.shape
        # Verify all results are finite
        assert np.all(np.isfinite(result))
        # Verify A_lambda increases with tau (more optical depth = more attenuation)
        assert np.all(np.diff(result) > 0)

        tau = np.array([0.5, 1.0, 1.5])
        costheta = np.array([0.5, 0.6, 0.7])
        albedo = np.array([0.2, 0.3, 0.4])
        result = att.get_A_lambda(tau, costheta=costheta, albedo=albedo)
        assert result.shape == tau.shape
        assert np.all(np.isfinite(result))

    def test_find_line_index(self):
        line_names = ['OII3727','Hbeta']
        ind = att.find_line_index('Hbeta',line_names)
        self.assertEqual(1, ind)
        ind = att.find_line_index('OII3726',line_names)
        self.assertEqual(0, ind)
        ind = att.find_line_index('OII3728',line_names)
        self.assertEqual(0, ind)
        ind = att.find_line_index('OII3730',line_names)
        self.assertEqual(None, ind)


    def test_gne_att_ratios(self):
        """Test that input attenuation ratios are conserved in the output
        for the 'ratios' mode, i.e. al = L_att / L_unatt == input ratio."""
        import os, tempfile, h5py
        from unittest.mock import patch

        # Retrieve a valid photoionisation model and its line catalogue
        photmod = list(c.line_names.keys())[0]
        lnames = list(c.line_names[photmod])
        nlines_total = len(lnames)

        # Dimensions
        ncomp = 2
        ngal = 5

        # Choose a subset of lines to provide ratios for
        n_ratio_lines = min(3, nlines_total)
        ratio_line_names = lnames[:n_ratio_lines]

        # Synthetic unattenuated emission-line luminosities (erg/s)
        np.random.seed(42)
        neblines = np.random.rand(nlines_total, ncomp, ngal) * 1e40 + 1e38

        # Attenuation ratios: fraction of light transmitted (0 < al < 1)
        ratios = {
            line: np.random.rand(ngal) * 0.5 + 0.3
            for line in ratio_line_names
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'test_nebular.hdf5')

            # --- Build the input HDF5 file ---
            with h5py.File(filepath, 'w') as f:
                header = f.create_group('header')
                header.attrs['redshift'] = 0.5
                header.attrs['photmod_sfr'] = photmod

                sfr_group = f.create_group('sfr_data')
                for i, line in enumerate(lnames):
                    sfr_group.create_dataset(
                        f'{line}_sfr', data=neblines[i])

                data_group = f.create_group('data')
                for line, r in ratios.items():
                    data_group.create_dataset(f'ratio_{line}', data=r)

            # --- Run attenuation in 'ratios' mode ---
            # Mock get_outnom so the same file is used for output
            with patch.object(att.io, 'get_outnom', return_value=filepath):
                att.gne_att(filepath, attmod='ratios', verbose=False)

            # --- Verify al = L_att / L_unatt equals the input ratios ---
            with h5py.File(filepath, 'r') as f:
                for line in ratio_line_names:
                    # Unattenuated: sum over components (as gne_att does
                    # when ratios are 1-D, i.e. per-galaxy)
                    unatt = f[f'sfr_data/{line}_sfr'][:]
                    if unatt.ndim > 1:
                        unatt_summed = np.sum(unatt, axis=0)
                    else:
                        unatt_summed = unatt

                    # Attenuated luminosity written by gne_att
                    attenuated = f[f'sfr_data/{line}_sfr_att'][:]

                    # Recovered attenuation factor
                    recovered_al = attenuated / unatt_summed

                    np.testing.assert_allclose(
                        recovered_al, ratios[line],
                        rtol=1e-5,
                        err_msg=(
                            f'Attenuation ratio not conserved '
                            f'for line {line}'
                        ),
                    )
        
if __name__ == '__main__':
    unittest.main(verbosity=2)
