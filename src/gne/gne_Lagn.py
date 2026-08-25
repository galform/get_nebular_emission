"""
.. moduleauthor:: Julen Expósito-Márquez <expox7@gmail.com>
.. contributions:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
"""
from typing import Optional
import numpy as np
import gne.gne_const as c
from gne.gne_io import read_data
from gne.gne_cosmology import age_of_universe
from gne.gne_stats import safe_sum_arrays

R_BULGE_COLUMNS = ['data/rgas_bulge', 'data/rstar_bulge', 'data/rbulge'] # Remove rstar_bulge or rgas_bulge
V_BULGE_COLUMNS = ['data/vbulge']
MGAS_BULGE_COLUMNS = ['data/mgas_bulge', "data/mcold_burst"]
MSTARS_BULGE_COLUMNS = ['data/mstars_bulge']


def get_Ledd(Mbh):
    '''
    Calculate the Eddington luminosity (erg/s)
    for a black hole of mass Mbh with Eq.9 in Griffin+2019

    Parameters
    ----------
    Mbh : array of floats
       Mass of the black hole (Msun).
     
    Returns
    -------
    Ledd : array of floats (erg/s)
    '''
    
    Ledd = 1.26e38*Mbh
    return Ledd # erg/s
    

def acc_rate_edd(Mbh):
    '''
    Calculate the Eddington mass accretion rate (Msun/yr)
    of a black hole following Eq.10 in Griffin+2019
    (or 14.9 in the book from Mo, van den Bosch and White)

    Parameters
    ----------
    Mbh : array of floats
         Mass of the black hole (Msun).
     
    Returns
    -------
    acc_rate : array of floats
    '''

    units2Msunyr = 1e-7*c.yr_to_s/c.Msun
    acc_rate = units2Msunyr*get_Ledd(Mbh)/(c.e_r_agn*c.c*c.c)

    return acc_rate # Msun/yr


def t_bulge(r_bulge, v_bulge): # Dynamical timescale of the bulge
    '''
    Given the bulge's half-mass radius and circular velocity at the half-mass radius,
    calculates the dynamical timescale of the bulge.

    Parameters
    ----------
    r_bulge : floats
     Half-mass radius of the bulge (Mpc).
    v_bulge : floats
     Circular velocity at the half-mass radius of the bulge (km/s).
     
    Returns
    -------
    dyn_time : floats
    '''
    
    dyn_time = r_bulge/v_bulge * 1e-5*c.Mpc_to_cm # s
    return dyn_time


def t_bulge_from_mass(r_bulge, m_bulge):
    '''
    Given the bulge's half-mass radius and mass,
    calculates the dynamical timescale of the bulge.

    Parameters
    ----------
    r_bulge : floats
     Half-mass radius of the bulge (Mpc).
    m_bulge : floats
     Mass of the bulge (Msun).
     
    Returns
    -------
    dyn_time : floats
        Dynamical timescale of the bulge (s).
    '''
    dyn_time = np.zeros_like(r_bulge, dtype=float)
    mask = (
        np.isfinite(r_bulge) & (r_bulge > 0) &
        np.isfinite(m_bulge) & (m_bulge > 0)
    )
    if np.any(mask):
        r_bulge_m = r_bulge[mask] * c.Mpc_to_cm /100 # m
        m_bulge_kg = m_bulge[mask] * c.Msun #kg
        dyn_time[mask] = np.sqrt(r_bulge_m**3 / (c.G * m_bulge_kg)) #s
    return dyn_time 


def acc_rate_quasar(M_bulge, r_bulge, v_bulge):
    '''
    Given the mass of the bulge, the half-mass radius and circular velocity 
    at the half-mass radius, calculates the accretion rate of the quasar (starburst) mode.

    Parameters
    ----------
    M_bulge : floats
     Mass of the bulge (Msun)
    r_bulge : floats
     Half-mass radius of the bulge (Mpc).
    v_bulge : floats
     Circular velocity at the half-mass radius of the bulge (km/s).
     
    Returns
    -------
    acc_rate : floats
    '''
    
    acc_rate = M_bulge*c.fbh/(t_bulge(r_bulge, v_bulge)*c.fq) # Msun/s
    return acc_rate


def acc_rate_radio(Mhot, Mbh, kagn=c.kagn, kagn_exp=c.kagn_exp):
    '''
    Given the mass of the hot gas and the mass of the black hole,
    calculates the accretion rate of the radio (hot gas) mode,
    following Eq. 1 in Henriques+2016.
    Default kagn values are derived from a GP20 sample 

    Parameters
    ----------
    Mhot : array of floats
        Mass of the hot gas (Msun).
    Mbh : array of floats
        Mass of the black hole (Msun).
    kagn : float
        Multiplicative factor for units Msun/yr
    kagn_exp : float
        Exponent for the dependence of Mdot with Mbh*Mhot
    
    Returns
    -------
    acc_rate : array of floats
    '''
    
    acc_rate = kagn*(Mhot*Mbh*1e-19)**kagn_exp # Msun/yr

    return acc_rate # Msun/yr


def get_Lagn_H14(Mdot,Mbh):
    '''
    Calculate the AGN bolometric luminosity
    following Sec.4.1 from Hirschmann+2014, and McCarthy+17

    Parameters
    ----------
    Mdot : array of floats
         Accretion rate onto the black hole (Msun/yr).
    Mbh : array of floats
       Mass of the black hole (Msun).
     
    Returns
    -------
    LagnH14 : array of floats
    '''

    Mdot_edd = acc_rate_edd(Mbh)
    fedd = Mdot/Mdot_edd

    units = 1e7*c.Msun/c.yr_to_s
    L_eff = units*Mdot*c.c*c.c*c.e_r_agn/(1-c.e_r_agn)

    L_ineff = 10.*get_Ledd(Mbh)*fedd**2

    LagnH14 = np.where(fedd > 0.1, L_eff, L_ineff)
    if np.ndim(LagnH14) == 0:
        return float(LagnH14)
    return LagnH14 #erg/s


def r_iso(spin):
    Z1 = 1 + (1 - np.abs(spin)**2)**(1/3) * ((1+np.abs(spin))**(1/3) + (1-np.abs(spin))**(1/3))
    Z2 = np.sqrt(3*np.abs(spin)**2 + Z1**2)
    r_iso = 3 + Z2 - np.sign(spin)*np.sqrt((3-Z1)*(3+Z1+2*Z2))
    
    return r_iso


def epsilon_td(spin):
    r = r_iso(spin)
    epsilon_td = 1 - (1 - (2/(3*r)))**(1/2) # Radiative accretion efficiency for a thin accretion disc (general approximation, page 5)
    
    return epsilon_td


def get_Lbol_from_mdot(
    M_bh_msun: np.ndarray,
    mdot_msun_per_yr: np.ndarray, 
    spin: Optional[np.ndarray] = None                  
    ) -> np.ndarray:
    """
    Compute L_bol (erg/s) for given BH mass (Msun), the accretion rate (Msun/yr), 
    and the dimensionless spin parameter (between -1 and 1).
    following eqs. (29)-(31) and (32) of the Griffin+2019.

    Parameters
    ----------
    M_bh_msun: np.ndarray
        BH mass (Msun)
    mdot_msun_per_yr: np.ndarray
        Accretion rate (Msun/yr)
    spin: np.ndarray
        Dimensionless spin parameter (between -1 and 1)
    Returns
    -------
    Lagn_G19: np.ndarray
        Bolometric luminosity (erg/s) following Griffin+2019 eqs. (29)-(31) and (32)
    """
    if np.shape(M_bh_msun) != np.shape(mdot_msun_per_yr):
            raise ValueError(f"M_bh shape {np.shape(M_bh_msun)} != Mdot shape {np.shape(mdot_msun_per_yr)}")
    
    M_bh = np.asarray(M_bh_msun, dtype=float)
    Mdot = np.asarray(mdot_msun_per_yr, dtype=float) * c.Msun * c.kilo / c.yr_to_s # g/s
    
    # Dimensionless mdot ≡ Mdot / Mdot_Edd with Mdot_Edd defined using 0.1 efficiency
    Mdot_edd_value = get_Ledd(M_bh) / (c.eta_acc_eff * c.c_cm * c.c_cm) # g/s
    mdot_dim = Mdot / Mdot_edd_value

    #clip the spin values to stay between -1 and 1 (physical values)
    if spin is None:
        epsTD =0.1
        r_lso_use = 6.0
    else:
        spin = np.clip(np.asarray(spin, dtype=float), -0.999, 0.999)
        r_lso_use = r_iso(spin) 
        epsTD = epsilon_td(spin)
        
    #I need to mask per case, but if no spin r_Iso and eps are a scalar the mask breaks    
    epsTD     = np.broadcast_to(epsTD,     np.shape(mdot_dim))
    r_lso_use = np.broadcast_to(r_lso_use, np.shape(mdot_dim))
    
    # eq. (32): boundary inside ADAF so L is continuous
    #New criteria from Griffin 2020 correction
    mcrit_super = c.eta_edd * (0.1 / epsTD)           
    # Initialize Lbol
    Lbol = np.zeros_like(Mdot, dtype=float)

    # CASE 1: low-Ṁ ADAF (mdot < mcrit_visc) -> eq. (29)
    mask1 = mdot_dim < c.acc_rate_crit_visc
    if np.any(mask1):
        Lbol[mask1] = (
            2e-4 * epsTD[mask1] * Mdot[mask1] * c.c_cm * c.c_cm
            * (c.lambda_adaf/5e-4)
            * ((1 - c.beta)/0.5)
            * (6.0 / r_lso_use[mask1])
        )

    # CASE 2: high-Ṁ ADAF (mcrit_visc ≤ mdot < mcrit_ADAF) -> eq. (30)
    mask2 = (mdot_dim >= c.acc_rate_crit_visc) & (mdot_dim < c.acc_rate_crit_adaf)
    if np.any(mask2):
        Lbol[mask2] = (
            0.2 * epsTD[mask2] * Mdot[mask2] * c.c_cm * c.c_cm
            * (mdot_dim[mask2] / (c.alpha_adaf**2))
            * (c.beta/0.5)
            * (6.0 / r_lso_use[mask2])
        )

    # CASE 3: thin disc (mcrit_ADAF ≤ mdot ≤ mcrit_super) -> eq. (28)
    #mask3 = (mdot_dim >= mcrit_ADAF) & (mdot_dim <= eta_Edd) <--- old limit before 2020 correction
    mask3 = (mdot_dim >= c.acc_rate_crit_adaf) & (mdot_dim <= mcrit_super)
    if np.any(mask3):
        Lbol[mask3] = epsTD[mask3] * Mdot[mask3] * c.c_cm * c.c_cm

    # CASE 4: super-Eddington (mdot > mcrit_super) -> Griffin+2020 eq. (2)
    #mask4 = mdot_dim > eta_Edd<--- old limit before 2020 correction
    mask4 = mdot_dim > mcrit_super
    if np.any(mask4):
        arg = mdot_dim[mask4]/mcrit_super[mask4] #>1 by construction so log should be ok
        Lbol[mask4] = c.eta_edd * (1 + np.log(arg)) * get_Ledd(M_bh[mask4])

    return Lbol # erg/s


# Alternative Griffin+19 implementation (component-wise luminosity sum):
#   L_bol = L(mdot_hh) + L(mdot_sb), with BOOL Bernoulli sampling on L_sb only.
# Because L(mdot) is nonlinear (ADAF / thin disc / super-Eddington), this
# decomposition can be more physically motivated when duty-cycle weights apply
# only to the starburst channel. Production uses get_Lagn_G19 below, which sums
# accretion rates before calling get_Lbol_from_mdot, matching Shark and Galform
# (mdot_tot = mdot_hh + mdot_sb, then L_bol = L(mdot_tot); BOOL turns SB mdot
# on/off for the instantaneous case).
#
# def get_Lagn_G19(
#     Mbh: np.ndarray,
#     mdot_hh: np.ndarray,
#     mdot_sb: np.ndarray,
#     spin: np.ndarray,
#     weights: Optional[np.ndarray] = None
#     ) -> tuple[np.ndarray, Optional[np.ndarray]]:
#     '''
#     Calculate the bolometric luminosity of the AGN following Griffin+2019
#
#     Parameters
#     ----------
#     Mbh: np.ndarray
#         Mass of the black hole (Msun)
#     mdot_hh: np.ndarray
#         Accretion rate onto the black hole (Msun/yr)
#     mdot_sb: np.ndarray
#         Accretion rate onto the black hole (Msun/yr)
#     spin: np.ndarray
#         Dimensionless spin parameter (between -1 and 1)
#     weights: Optional[np.ndarray]
#         Weights for the bolometric luminosity of the BHs at the output snapshot (erg/s)
#
#     Returns
#     -------
#     Lagn: np.ndarray
#         Bolometric luminosity of the AGN (erg/s)
#     Lagn_insta: Optional[np.ndarray]
#         Instantaneous bolometric luminosity of the AGN (erg/s) if weights is not None.
#     '''
#
#     Lbol_hh = get_Lbol_from_mdot(Mbh, mdot_hh, spin)
#     Lbol_sb = get_Lbol_from_mdot(Mbh, mdot_sb, spin)
#     Lbol = Lbol_hh + Lbol_sb
#     if weights is None:
#         return Lbol, None
#
#     rng = np.random.default_rng(42)
#     on_sb = rng.random(weights.shape[0]) < weights
#     Lbol_sb_insta = np.where(on_sb, Lbol_sb, 0.0)   # erg/s
#     Lbol_insta = Lbol_hh + Lbol_sb_insta
#
#     return Lbol, Lbol_insta


def get_Lagn_G19(
    Mbh: np.ndarray, 
    mdot_hh: np.ndarray, 
    mdot_sb: np.ndarray, 
    spin: np.ndarray, 
    weights: Optional[np.ndarray] = None
    ) -> tuple[np.ndarray, Optional[np.ndarray]]:
    '''
    Calculate the bolometric luminosity of the AGN following Griffin+2019.

    Snapshot luminosity: mdot = mdot_hh + mdot_sb, then L_bol = L(mdot).
    Instantaneous luminosity: Bernoulli sampling on mdot_sb, then
    mdot = mdot_hh + where(on_sb, mdot_sb, 0), L_bol_insta = L(mdot).
    This matches Shark/Galform (sum accretion rates before L(mdot)); see the
    commented alternative above (L(mdot_hh) + L(mdot_sb)).

    Parameters
    ----------
    Mbh: np.ndarray
        Mass of the black hole (Msun)
    mdot_hh: np.ndarray
        Accretion rate onto the black hole (Msun/yr)
    mdot_sb: np.ndarray
        Accretion rate onto the black hole (Msun/yr)
    spin: np.ndarray
        Dimensionless spin parameter (between -1 and 1)
    weights: Optional[np.ndarray]
        Weights for the bolometric luminosity of the BHs at the output snapshot (erg/s)

    Returns
    -------
    Lagn_noinsta: np.ndarray
        Bolometric luminosity of the AGN (erg/s)
    Lagn: Optional[np.ndarray]
        Instantaneous bolometric luminosity of the AGN (erg/s) if weights is not None.
    '''
    mdot = mdot_hh + mdot_sb
    Lagn_noinsta = get_Lbol_from_mdot(Mbh, mdot, spin)

    if weights is None:
        return Lagn_noinsta, None

    rng = np.random.default_rng(42)              
    on_sb = rng.random(weights.shape[0]) < weights   
    mdot = mdot_hh + np.where(on_sb, mdot_sb, 0.0)   
    Lagn = get_Lbol_from_mdot(Mbh, mdot, spin) # erg/s

    return Lagn_noinsta, Lagn


def Rsch(Mbh):
    '''
    Schwarzschild radius (Mpc) given the mass of the black hole (Msun)

    Parameters
    ----------
    Mbh : floats
     Mass of the black hole (Msun)
     
    Returns
    -------
    Rs : floats
    '''

    Mbh_kg = Mbh*c.Msun
    Rs = 2*c.G*Mbh_kg/(c.c**2)/(c.mega*c.parsec)
    return Rs

def _get_weights_insta_Lagn(
    infile: str,
    cut: np.ndarray,
    redshift: float,
    redshift_previous: Optional[float] = None,
    h0: float = None,
    units_h0: bool = False,
    inputformat: str = 'hdf5',
    params: str = 'Lagn',
    tau_fold: Optional[float] = None,
    testing: bool = False,
    verbose: bool = True,
    ) -> np.ndarray:
    '''
    Calculate the weights for the bolometric luminosity (erg/s)
    for active AGNs at the output snapshot

    Parameters
    ----------
    infile: str
        Name of the input file
    cut: np.ndarray
        Indexes of selected galaxies for the study
    redshift: float
        Redshift
    redshift_previous: Optional[float]
        Redshift of the previous snapshot
    h0 : float
        Hubble constant h (H0 = 100 h km/s/Mpc). Required when units_h0 is True.
    units_h0: bool
        True if input units use h.
    inputformat: str
        Format of the input file.
    params: str
        Names of the parameters to calculate the AGN emission. 
    tau_fold: Optional[float]
        Ratio of lifetime of AGN episode to bulge dynamical timescale. 
        - The fiducial value used in Shark is 1.0.
        - If is None we use c.fq as the weights.
    testing: bool
        If True only run over few entries for testing purposes
    verbose: bool
        If True print out messages

    Returns
    -------
    weights: np.ndarray
        Weights for the bolometric luminosity of the BHs at the output snapshot (erg/s)
    '''
    r_bulge_col = None
    v_bulge_col = None
    mgas_bulge_col = None
    mstars_bulge_col = None

    for param in params:
        if param in R_BULGE_COLUMNS:
            r_bulge_col = param
        elif param in V_BULGE_COLUMNS:
            v_bulge_col = param
        elif param in MGAS_BULGE_COLUMNS:
            mgas_bulge_col = param
        elif param in MSTARS_BULGE_COLUMNS:
            mstars_bulge_col = param

    params_to_read = [r_bulge_col, v_bulge_col]
    if r_bulge_col is None:
        raise ValueError('r_bulge must be provided')
    if v_bulge_col is None:
        params_to_read = [r_bulge_col, mgas_bulge_col, mstars_bulge_col]
        if mgas_bulge_col is None or mstars_bulge_col is None:
            raise ValueError('v_bulge must be provided if mgas_bulge or mstars_bulge is not provided')

    vals = read_data(infile,cut,inputformat=inputformat,
                     params=params_to_read,testing=testing,verbose=verbose)

    if units_h0 and h0 is None:
        raise ValueError('h0 must be provided when units_h0 is True')

    r_bulge = np.asarray(vals[0], dtype=float)
    if units_h0:
        r_bulge = r_bulge/h0

    if v_bulge_col is not None:
        v_bulge = np.asarray(vals[1], dtype=float)
        t_bulge_val = t_bulge(r_bulge, v_bulge)
    else:
        m_bulge = np.asarray(vals[1], dtype=float) + np.asarray(vals[2], dtype=float)
        if units_h0:
            m_bulge = m_bulge/h0
        t_bulge_val = t_bulge_from_mass(r_bulge, m_bulge)

    t_bulge_val = t_bulge_val / c.yr_to_s # yr
    t_snapshot = age_of_universe(redshift) * 1e9 # Gyr -> yr
    delta_t_window = t_snapshot / 10 # yr

    if redshift_previous is not None:
        delta_t_window = abs(t_snapshot - age_of_universe(redshift_previous) * 1e9) # Gyr -> yr

    if tau_fold is None:
        tau_fold = c.fq

    return np.minimum(tau_fold * t_bulge_val/delta_t_window, 1.0)

def get_Lagn_insta(
    Lagn: np.ndarray,
    infile: str,
    cut: np.ndarray,
    redshift: float,
    redshift_previous: Optional[float] = None,
    h0: float = None,
    units_h0: bool = False,
    inputformat: str = 'hdf5',
    params: str = 'Lagn',
    tau_fold: Optional[float] = None,
    testing: bool = False,
    verbose: bool = True
) -> np.ndarray:
    '''
    Calculate the bolometric luminosity of BHs (erg/s)

    Parameters
    ----------
    Lagn: np.ndarray
        Bolometric luminosity of the BHs (erg/s)
    infile: str
        Name of the input file
    cut: np.ndarray
        Indexes of selected galaxies for the study
    redshift: float
        Redshift
    h0 : float
        Hubble constant h (H0 = 100 h km/s/Mpc). Required when units_h0 is True.
    units_h0: bool
        True if input units use h.
    inputformat: str
        Format of the input file.
    params: str
        Names of the parameters to calculate the AGN emission. 
    tau_fold: Optional[float]
        Ratio of lifetime of AGN episode to bulge dynamical timescale. 
        - The fiducial value used in Shark is 1.0.
        - If is None we use c.fq as the weights.
    testing: bool
        If True only run over few entries for testing purposes
    verbose: bool
        If True print out messages
    fq: float
        Fraction of the mass in the bulge that is converted into the black hole

    Returns
    -------
    Lagn_insta : array of floats
        Bolometric luminosity of the BHs (erg/s)
    '''
    weights = _get_weights_insta_Lagn(
        infile=infile,
        cut=cut,
        redshift=redshift,
        redshift_previous=redshift_previous,
        h0=h0,
        units_h0=units_h0,
        inputformat=inputformat,
        params=params,
        tau_fold=tau_fold,
        testing=testing,
        verbose=verbose)
   
    rng = np.random.default_rng(42)              
    on_sb = rng.random(weights.shape[0]) < weights   
    Lagn_insta = np.where(on_sb, Lagn, 0.0)   # erg/s
    return Lagn_insta


def get_Lagn(infile,cut,inputformat='hdf5',params='Lagn',Lagn_inputs='Lagn',
             h0=None,redshift=None,redshift_previous=None,units_h0=False,units_Gyr=False,units_L=0,
             testing=False,verbose=True, calculate_Lagn_insta=None, Lagn_insta_params=None, tau_fold=None
             ) -> tuple[np.ndarray, Optional[np.ndarray]]:
    '''
    Calculate or get the bolometric luminosity of BHs (erg/s) 

    Parameters
    ----------
    infile : string
        Name of the input file
    cut : array of integers
        Indexes of selected galaxies for the study
    inputformat : string
        Format of the input file.
    params : array of strings
        Names of the parameters to calculate the AGN emission. 
    Lagn_inputs : string
        Type of calculation to obtain Lagn
    h0 : float
        Hubble constant
    redshift : float
        Redshift
    redshift_previous : float
        Redshift of the previous snapshot
    units_h0: boolean
        True if input units with h
    units_Gyr: boolean
        True if input units with */Gyr
    units_L: integer
        0: input units [L]=erg/s  (default);
        1: input units [L]=1e40 h^-2 erg/s
        2: input units [L]=1e40 erg/s
    calculate_Lagn_insta: boolean
        If True calculate the instantaneous Lagn
    Lagn_insta_params: list of strings
        Names of the parameters to calculate the instantaneous Lagn
    tau_fold: Optional[float]
        Ratio of lifetime of AGN episode to bulge dynamical timescale. 
        - The fiducial value used in Shark is 1.0.
        - If is None we use c.fq as the weights.
    testing : boolean
        If True only run over few entries for testing purposes
    verbose : boolean
        If True print out messages
    
    Returns
    -------
    Lagn_noinsta : array of floats
        Bolometric luminosity of the BHs for all AGN that have been active from last snapshot (erg/s)
    Lagn: array of floats
        Instantaneous bolometric luminosity of the BHs (erg/s) if calculate_Lagn_insta is True.
    '''

    vals = read_data(infile,cut,inputformat=inputformat,
                     params=params,
                     testing=testing,verbose=verbose)

    if Lagn_inputs=='Lagn':
        if np.ndim(vals) == 1:
            Lagn = vals
        else:
            Lagn = vals[0]
        if units_L==1:
            Lagn = np.asarray(Lagn, dtype=np.float64)
            cfac = np.float64(1e40)/(h0*h0)
            Lagn = Lagn*cfac
        elif units_L==2:
            cfac = np.float64(1e40)
            Lagn = Lagn*cfac
        elif units_L!=0:
            raise ValueError('units_L must be 0, 1 or 2')
        return Lagn, None # erg/s
    
    elif Lagn_inputs=='Hirschman+14':
        Mdot = vals[0]
        Mbh = vals[1]
        if units_h0:
            Mdot = Mdot/h0
            Mbh = Mbh/h0
        if units_Gyr:
            Mdot = Mdot/1e9

        Lagn = get_Lagn_H14(Mdot,Mbh)
        return Lagn, None # erg/s

    elif Lagn_inputs=='Griffin+19':
        Mbh = vals[0]
        mdot_hh = vals[1]
        mdot_sb = vals[2]
        spin = None

        if units_h0:
            Mbh = Mbh/h0 # Msun
            mdot_hh = mdot_hh/h0 
            mdot_sb = mdot_sb/h0 
        
        if units_Gyr:
            mdot_hh = mdot_hh / 1e9 # Msun/yr
            mdot_sb = mdot_sb / 1e9 # Msun/yr

        if len(vals) > 3:
            spin = vals[3]

        weights = None
        if calculate_Lagn_insta:
            weights = _get_weights_insta_Lagn(
                infile=infile,
                cut=cut,
                redshift=redshift,
                redshift_previous=redshift_previous,
                h0=h0,
                units_h0=units_h0,
                inputformat=inputformat,
                tau_fold=tau_fold,
                params=Lagn_insta_params,
                testing=testing,
                verbose=verbose
            )
        
        Lagn_noinsta, Lagn = get_Lagn_G19(Mbh, mdot_hh, mdot_sb, spin, weights)
        return Lagn_noinsta, Lagn # erg/s

    raise ValueError(f"Invalid Lagn_inputs: {Lagn_inputs}")
