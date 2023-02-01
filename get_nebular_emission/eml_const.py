import numpy as np
#from astroML.cosmology import Cosmology

notnum = -999.

mp = 1.67e-27 #Proton mass, kg

zsun = 0.0134 # Asplund 2009
# zsum = 0.014 #CMB constant ref ?
ohsun = 8.69  # Allende Prieto 2001 and Asplund 2009 (12 + log10(O/H))sun

vol_pm = 542.16**3. # Volume of the P-Millenium simulation

unemods = ['kashino20']

photmods = ['gutkin16']
mod_lim = {'gutkin16': r"nebular_data/gutkin_tables/limits_gutkin.txt"}

attmods = ['ratios', 'cardelli89']

inputformats = ['textfile','HDF5']

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

IMF_M = {'Kennicut': 1, 'Salpeter': 0.47, 'Kroupa': 0.74, 'Chabrier': 0.81, 
            'Baldry&Glazebrook': 0.85, 'Top-heavy': 1.11}

IMF_SFR = {'Kennicut': 1, 'Salpeter': 0.94, 'Kroupa': 1.49, 'Chabrier': 1.57, 
            'Baldry&Glazebrook': 2.26, 'Top-heavy': 3.13}

# FOR CONVERSION FROM LYMANN CONTINUUM PHOTONS TO SFR
# It is assumed that a SFR of 1 Msun/yr produces 9.85 · 10^52 photons/s for Kennicut IMF.
# Reference: Chomiuk & Povich 2011, pag. 2: "According to Kennicutt et al. (1994) and Kennicutt
# (1998a), a SFR of 1 M⊙ yr−1 produces a Lyman continuum photon rate Nc = 9.26 × 1052 photon s−1
# for the Salpeter (1955) IMF (assuming a mass range of 0.1–100 Msun)."
# Reescaled to Kennicut, it gives our number.

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
                         0.040])
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
                           '040'])
}

lines_model = {
    "gutkin16" : np.array(['OII3727','Hbeta','OIII4959','OIII5007',
                           'NII6548','Halpha','NII6583','SII6717',
                           'SII6731','NV1240','CIV1548','CIV1551',
                           'HeII1640', 'OIII1661','OIII1666','SiIII1883',
                           'SiIII1888', 'CIII1908'])
    }

wavelength_model = {
    "gutkin16" : np.array([3727,4861,4959,5007,6548,6563,6583,
                           6717,6731,1240,1548,1551,1640,1661,
                           1666,1883,1888,1908])
    }

#------------------------------------------------------------------------------------
#   GALFORM:
#------------------------------------------------------------------------------------
# GALFORM tables have [0,1,3,5,6,7,8] of the gutkin16 model.

line_headers = ['L_tot_', 'L_disk_', 'L_bulge_']
att_ext = '_ext'

#------------------------------------------------------------------------------------
#   Cosmology:
#------------------------------------------------------------------------------------
h = 0.6777 # Hubble constant H = 100 h km s**-1 Mpc**-1
omega0=0.307
omegab = 0.0482
lambda0=0.693

#------------------------------------------------------------------------------------
#   Conversion factors:
#------------------------------------------------------------------------------------
kg_to_Msun=1./(1.989e30)
Msun_to_kg=1.989e30
Mpc_to_cm=3.086e24
