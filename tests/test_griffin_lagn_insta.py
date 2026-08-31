"""Tests for Griffin+19 and Bravo+25 L_bol and instantaneous Lagn sampling."""

from unittest.mock import patch

import numpy as np
import pytest
from numpy.testing import assert_allclose

import gne.gne_const as c
import gne.gne_Lagn as agn

MBH_MSUN = np.array([1e8, 1e8])
MDOT_HH = np.array([0.01, 0.1])
MDOT_SB = np.array([0.02, 0.05])
BULGE_INSTA_PARAMS = ['data/rgas_bulge', 'data/vbulge']
R_BULGE_MPC = np.array([0.001])
V_BULGE_KMS = np.array([200.0])


def test_get_Lagn_G19_returns_summed_luminosity():
    """Griffin+19 combines HH and SB as L(mdot_hh) + L(mdot_sb)."""
    l_agn = agn.get_Lagn_G19(MBH_MSUN, MDOT_HH, MDOT_SB, None)
    expected = (
        agn.get_Lbol_from_mdot(MBH_MSUN, MDOT_HH, None)
        + agn.get_Lbol_from_mdot(MBH_MSUN, MDOT_SB, None)
    )
    assert_allclose(l_agn, expected)


def test_get_Lagn_Bravo25_without_weights_returns_window_lbol():
    """Snapshot L_bol uses summed mdot when BOOL weights are disabled."""
    l_agn, l_noinsta = agn.get_Lagn_Bravo25(
        MBH_MSUN, MDOT_HH, MDOT_SB, None, weights=None)
    expected = agn.get_Lbol_from_mdot(MBH_MSUN, MDOT_HH + MDOT_SB, None)
    assert_allclose(l_agn, expected)
    assert l_noinsta is None


def test_get_Lagn_Bravo25_with_weights_is_reproducible():
    """Bernoulli sampling on mdot_sb is deterministic for a fixed RNG seed."""
    weights = np.array([1.0, 0.0])
    args = (MBH_MSUN, MDOT_HH, MDOT_SB, None, weights)
    l_agn_a, l_noinsta_a = agn.get_Lagn_Bravo25(*args)
    l_agn_b, l_noinsta_b = agn.get_Lagn_Bravo25(*args)
    assert_allclose(l_agn_a, l_agn_b)
    assert_allclose(l_noinsta_a, l_noinsta_b)
    l_hh_only = agn.get_Lbol_from_mdot(MBH_MSUN[1:2], MDOT_HH[1:2], None)
    assert_allclose(l_agn_a[1], l_hh_only[0])


@pytest.mark.parametrize('calculate_Lagn_insta', [False, True])
@patch('gne.gne_Lagn._get_weights_insta_Lagn')
@patch('gne.gne_Lagn.read_data')
def test_get_Lagn_griffin19_ignores_instantaneous_flag(
        mock_read_data, mock_weights, calculate_Lagn_insta):
    """Griffin+19 never applies BOOL sampling."""
    mock_read_data.return_value = (MBH_MSUN, MDOT_HH, MDOT_SB)
    l_agn, l_noinsta = agn.get_Lagn(
        infile='dummy.hdf5',
        cut=None,
        Lagn_inputs='Griffin+19',
        calculate_Lagn_insta=calculate_Lagn_insta,
        Lagn_insta_params=BULGE_INSTA_PARAMS,
        redshift=1.0,
        redshift_previous=0.5,
    )
    expected = agn.get_Lagn_G19(MBH_MSUN, MDOT_HH, MDOT_SB, None)
    assert_allclose(l_agn, expected)
    assert l_noinsta is None
    mock_weights.assert_not_called()


@pytest.mark.parametrize('calculate_Lagn_insta', [False, True])
@patch('gne.gne_Lagn._get_weights_insta_Lagn')
@patch('gne.gne_Lagn.read_data')
def test_get_Lagn_bravo25_modes(
        mock_read_data, mock_weights, calculate_Lagn_insta):
    """Bravo+25 returns window L_bol and optional instantaneous L_bol."""
    mock_read_data.return_value = (MBH_MSUN, MDOT_HH, MDOT_SB)
    mock_weights.return_value = np.array([1.0, 1.0])
    l_agn, l_noinsta = agn.get_Lagn(
        infile='dummy.hdf5',
        cut=None,
        Lagn_inputs='Bravo+25',
        calculate_Lagn_insta=calculate_Lagn_insta,
        Lagn_insta_params=BULGE_INSTA_PARAMS,
        redshift=1.0,
        redshift_previous=0.5,
    )
    expected_window = agn.get_Lbol_from_mdot(MBH_MSUN, MDOT_HH + MDOT_SB, None)
    if calculate_Lagn_insta:
        mock_weights.assert_called_once()
        assert l_noinsta is not None
        assert_allclose(l_noinsta, expected_window)
        assert l_agn is not None
        assert np.all(l_agn >= 0.0)
    else:
        mock_weights.assert_not_called()
        assert_allclose(l_agn, expected_window)
        assert l_noinsta is None


@patch('gne.gne_Lagn._get_weights_insta_Lagn')
def test_get_Lagn_insta_zeros_when_duty_cycle_weight_is_zero(mock_weights):
    """Generic BOOL path sets L_bol to zero when the sampling weight is null."""
    mock_weights.return_value = np.array([0.0, 1.0])
    lagn_noinsta = np.array([1e44, 2e44])
    lagn_insta = agn.get_Lagn_insta(
        lagn_noinsta,
        infile='dummy.hdf5',
        cut=np.array([0, 1]),
        redshift=1.0,
        params=BULGE_INSTA_PARAMS,
    )
    assert lagn_insta[0] == 0.0
    assert lagn_insta[1] == lagn_noinsta[1]


@patch('gne.gne_Lagn.age_of_universe')
@patch('gne.gne_Lagn.read_data')
def test_get_weights_insta_lagn_uses_tau_fold_constant(mock_read, mock_age):
    """Duty-cycle weights use c.tau_fold from gne_const."""
    mock_read.return_value = (R_BULGE_MPC, V_BULGE_KMS)
    mock_age.side_effect = lambda redshift: 10.0 if redshift == 1.0 else 8.0
    weights = agn._get_weights_insta_Lagn(
        infile='dummy.hdf5',
        cut=np.array([0]),
        redshift=1.0,
        redshift_previous=0.5,
        params=BULGE_INSTA_PARAMS,
    )
    t_bulge_yr = agn.t_bulge(R_BULGE_MPC[0], V_BULGE_KMS[0]) / c.yr_to_s
    delta_t_window = abs(10.0 - 8.0) * 1e9
    expected = np.minimum(c.tau_fold * t_bulge_yr / delta_t_window, 1.0)
    assert_allclose(weights, [expected])
