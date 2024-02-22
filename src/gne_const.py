import numpy as np
#from astroML.cosmology import Cosmology

#------------------------------------------------------------------------------------
#   Conversion factors:
#------------------------------------------------------------------------------------
kg_to_Msun=1./(1.989e30)
Msun_to_kg=1.989e30
Mpc_to_cm=3.086e24
yr_to_s=3.154e7
#------------------------------------------------------------------------------------

notnum = -999.

mp = 1.67e-27 # Proton mass, kg
c = 2.998e10 # Light velocity, cm/s
planck = 4.135e-15*1.6e-12 # Planck constant, erg*s
boltzmann = 1.38e-23 * 1e4 * kg_to_Msun/(Mpc_to_cm**2) # Boltzmann constant, Mpc^2 Ms s^-2 K^-1
G = 4.3e-9 # Gravitational constant, km^2 Mpc Ms^-1 s^-2
nH_gal = 100 # cm^-3

zsun = 0.0134 # Asplund 2009
ohsun = 8.69  # Allende Prieto 2001 and Asplund 2009 (12 + log10(O/H))sun

vol_pm = 542.16**3. # Volume of the P-Millenium simulation
vol_sage = 1000**3  # Volume of the UNITSIM1

#------------------------------------------------------------------------------------
#   Orsi et. al. 2014
#------------------------------------------------------------------------------------
Z0_orsi = 0.012
q0_orsi = 2.8e7 # cm/s
#------------------------------------------------------------------------------------

photmods = ['gutkin16', 'feltre16']
mod_lim = {'gutkin16': r"src/nebular_data/gutkin16_tables/limits_gutkin.txt",
           'feltre16': r"src/nebular_data/feltre16_tables/limits_feltre.txt"}

unemods = ['kashino20', 'orsi14', 'panuzzo03']

attmods = ['ratios', 'cardelli89']

inputformats = ['txt','hdf5']

# For Kennicut IMF -> M(Kenn) = corr * M(IMF)
# ---------------------------
# Salpeter 0.47
# Kroupa 0.74
# Chabrier 0.81
# Baldry & Glazebrook 0.85

# SFR(Kenn) = corr * SFR(IMF)
# ---------------------------
# Salpeter 0.94
# Kroupa 1.49
# Chabrier 1.57
# Baldry & Glazebrook 2.26
# Top-heavy IMF (x=1) 3.13

######### Lacey et. al. 2016

IMF_M = {'Kennicut': 1, 'Salpeter': 0.47, 'Kroupa': 0.74, 'Chabrier': 0.81, 
            'Baldry&Glazebrook': 0.85, 'Top-heavy': 1.11}

IMF_SFR = {'Kennicut': 1, 'Salpeter': 0.94, 'Kroupa': 1.49, 'Chabrier': 1.57, 
            'Baldry&Glazebrook': 2.26, 'Top-heavy': 3.13}

phot_to_sfr_kenn = 9.85e52 # phot/s

# FOR CONVERSION FROM LYMANN CONTINUUM PHOTONS TO SFR
# It is assumed that a SFR of 1 Msun/yr produces 9.85 · 10^52 photons/s for Kennicut IMF.
# Reference: Chomiuk & Povich 2011, pag. 2: "According to Kennicutt et al. (1994) and Kennicutt
# (1998a), a SFR of 1 M⊙ yr−1 produces a Lyman continuum photon rate Nc = 9.26 × 1052 photon s−1
# for the Salpeter (1955) IMF (assuming a mass range of 0.1–100 Msun)."
# Reescaled to Kennicut, it gives our number.

#------------------------------------------------------------------------------------

# Mean ratio between global cold gas density at different redshift for GP20

med_to_low = 1.74 # 0.8 to 0
high_to_low = 1.58 # 1.5 to 0

#------------------------------------------------------------------------------------
#    Atenuation:
#------------------------------------------------------------------------------------

# De Barros et. al. 2016 - Halpha, 3.33, z = (0.004, 0.02)
# Holden et. al. 2016 - 1600A, 2.27, z = 0.02
# Calzetti et. al. 2021 - Halpha, 2.54, z = (0.00258, 0.00261)
# Valentino et. al. 2017 - Halpha, 1.75, z = 1.55
# Buat et. al. 2018 - Halpha, 1.85, z = (0.6, 1.6)
# Saito et. al. 2020 - Halpha, 5/(z+2.2), z < 0.46
# Saito et. al. 2020 - OII (3727A, 3729A), 5/(z+2.2), z = (0.48,1.54)

def saito_att(z):
    if z < 2.8:
        return (z+2.2)/5
    else:
        return 1

def line_att_coef_all(z):
    return saito_att(z)

def line_att_coef_func(z=0):
    saito = saito_att(z)
    line_att_coef = {'OII3727' : saito,
                     'Hbeta' : saito,
                     'OIII4959' : saito,
                     'OIII5007' : saito,
                     'OI6300' : saito,
                     'NII6548' : saito,
                     'Halpha' : saito,
                     'NII6584' : saito,
                     'SII6717' : saito,
                     'SII6731' : saito,
                     'NV1240' : saito,
                     'CIV1548' : saito,
                     'CIV1551' : saito,
                     'HeII1640' : saito, 
                     'OIII1661' : saito,
                     'OIII1666' : saito,
                     'SiIII1883' : saito,
                     'SiIII1888' : saito, 
                     'CIII1907' : saito,
                     'CIII1908' : saito,
                     'CIII1910' : saito
        }
    return line_att_coef

#------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
#    AGNs:
#------------------------------------------------------------------------------------

AGNinputs = ['Lagn', 'acc_rate', 'acc_rates', 'radio_mode', 'quasar_mode', 'complete']

# Griffin et. al 2019:
alpha_adaf = 0.1 # Viscosity parameter for ADAFs
alpha_td = 0.1 # Viscosity parameter for TDs
lambda_adaf = 0.2 # Fraction of viscous energy transferred to electrons in ADAF
acc_rate_crit_adaf = 0.01 # Boundary between thin disc and ADAF accretion (in terms of ratio between accretion rate and eddington's)
eta_edd = 4 # Super-Eddington suppression factor
fq = 10 # Ratio of lifetime of AGN episode to bulge dynamical timescale
fbh = 0.005 # Fraction of the mass of stars formed in a starburst accreted onto a black hole

beta = 1 - alpha_adaf/0.55
acc_rate_crit_visc = 0.001*(lambda_adaf/0.0005)*((1-beta)/beta)*alpha_adaf**2 
# Boundary between the two adaf regimes

spin_bh = 0.67 # 0, 0.3, 0.5, 0.67, 0.9, 0.95, 1

# Lagos et. al 2008:
kappa_agn = 5.44e-4#2.46e-4 3.17e-4 5.44e-4
kappa_agn_exp = 0.597 #0.624 #0.460 0.597
# 1e-4  # Efficiency of cold gas accretion onto the BH during gas cooling (1e-4 in Lagos et. al)

# Lacey et. al 2016
epsilon_heat = 0.02 # BH heating efficienty 

nH_AGN = 1000 # cm^-3
radius_NLR = 0.001 # Mpc

#------------------------------------------------------------------------------------

# Relationship between epsilon and M and Z
epsilon_a_sfr = -0.1633
epsilon_b_sfr = 0.3776

#------------------------------------------------------------------------------------

zmet = {
    "gutkin16" : np.array([0.0001,
                         0.0002,
                         0.0005,
                         0.001,
                         0.002,
                         0.004,
                         0.006,
                         0.008,
                         0.010,
                         0.014,
                         0.017,
                         0.020,
                         0.030,
                         0.040]),
    "feltre16" : np.array([0.0001,
                         0.0002,
                         0.0005,
                         0.001,
                         0.002,
                         0.004,
                         0.006,
                         0.008,
                         0.014,
                         0.017,
                         0.020,
                         0.030,
                         0.040,
                         0.050,
                         0.060,
                         0.070])
}

zmet_reduced = {
    "gutkin16" : np.array([0.0001,
                         0.002,
                         0.014,
                         0.030])
}

zmet_str = {
    "gutkin16" : np.array(['0001',
                           '0002',
                           '0005',
                           '001',
                           '002',
                           '004',
                           '006',
                           '008',
                           '010',
                           '014',
                           '017',
                           '020',
                           '030',
                           '040']),
    "feltre16" : np.array(['0001',
                           '0002',
                           '0005',
                           '001',
                           '002',
                           '004',
                           '006',
                           '008',
                           '014',
                           '017',
                           '020',
                           '030',
                           '040',
                           '050',
                           '060',
                           '070'])
}

zmet_reduced = {
    "gutkin16" : np.array([0.0001,
                         0.002,
                         0.014,
                         0.030])
}

lines_model = {
    "gutkin16" : np.array(['OII3727','Hbeta','OIII4959','OIII5007',
                           'NII6548','Halpha','NII6584','SII6717',
                           'SII6731','NV1240','CIV1548','CIV1551',
                           'HeII1640', 'OIII1661','OIII1666','SiIII1883',
                           'SiIII1888', 'CIII1908']),
    "feltre16" : np.array(['OII3727','Hbeta','OIII4959','OIII5007',
                           'OI6300','NII6548','Halpha','NII6584','SII6717',
                           'SII6731','NV1240','CIV1548','CIV1551',
                           'HeII1640', 'OIII1661','OIII1666','SiIII1883',
                           'SiIII1888', 'CIII1907','CIII1910'])
    }

wavelength_model = {
    "gutkin16" : np.array([3727,4861,4959,5007,6548,6563,6584,
                           6717,6731,1240,1548,1551,1640,1661,
                           1666,1883,1888,1908]),
    "feltre16" : np.array([3727,4861,4959,5007,6300,6548,6563,6584,
                           6717,6731,1240,1548,1551,1640,1661,
                           1666,1883,1888,1907,1910])
    }

def coef_att_line_model_func(z=0):
    line_att_coef = line_att_coef_func(z)

    coef_att_line_model = {
        "gutkin16" : np.array([line_att_coef[line] for line in lines_model["gutkin16"]]),
        "feltre16" : np.array([line_att_coef[line] for line in lines_model["feltre16"]])
        }

    return coef_att_line_model

#------------------------------------------------------------------------------------
#   GALFORM:
#------------------------------------------------------------------------------------
# GALFORM tables have [0,1,3,5,6,7,8] of the gutkin16 model.

line_headers = ['L_tot_', 'L_disk_', 'L_bulge_']
att_ext = '_ext'

gal_headers = ['mag_LC_r_disk', 'mag_LC_r_bulge', 'zcold', 'mcold', 
           'zcold_burst', 'mcold_burst', 'mstardot_average', 'L_tot_Halpha',
           'L_tot_NII6583', 'L_tot_Hbeta', 'L_tot_OIII5007', 'mstars_total', 'is_central',
           'mstardot', 'mstardot_burst', 'mstars_bulge', 'L_tot_OII3727', 'L_tot_SII6716',
           'L_tot_SII6731', 'mag_SDSS_r_o_t', 'L_tot_Halpha_ext', 'L_tot_Hbeta_ext', 
           'L_tot_OII3727_ext', 'L_tot_OIII5007_ext', 'L_disk_Halpha', 
           'L_disk_Halpha_ext', 'L_bulge_Halpha', 'L_bulge_Halpha_ext', 
           'L_disk_Hbeta', 'L_disk_Hbeta_ext', 'L_bulge_Hbeta', 'L_bulge_Hbeta_ext',
           'L_disk_OIII5007', 'L_disk_OIII5007_ext', 'L_bulge_OIII5007', 'L_bulge_OIII5007_ext', 
           'L_disk_NII6583', 'L_disk_NII6583_ext', 'L_bulge_NII6583', 'L_bulge_NII6583_ext', 
           'L_disk_OII3727', 'L_disk_OII3727_ext', 'L_bulge_OII3727', 'L_bulge_OII3727_ext', 
           'L_disk_SII6717', 'L_disk_SII6717_ext', 'L_bulge_SII6717', 'L_bulge_SII6717_ext', 
           'L_disk_SII6731', 'L_disk_SII6731_ext', 'L_bulge_SII6731', 'L_bulge_SII6731_ext']
# Comentar en la documentación el tema de añadir líneas aquí.

reff_to_scale_high = 0.0875
halfmass_to_reff = 1/1.67

#------------------------------------------------------------------------------------
#   Plots:
#------------------------------------------------------------------------------------
BPT_lines = {'NII': ['SFR_Composite', 'Composite_AGN', 'LINER_NIIlim', 'LINER_OIIIlim'], 
             'SII': ['SFR_AGN', 'Seyfert_LINER'], 'OI': ['SFR_AGN', 'Seyfert_LINER']}

#------------------------------------------------------------------------------------
#   Cosmology:
#------------------------------------------------------------------------------------
h = 0.704#0.6777 # Hubble constant H = 100 h km s**-1 Mpc**-1
omega0=0.307
omegab = 0.0482
lambda0=0.693
