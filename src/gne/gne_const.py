import os
import numpy as np

notnum    = -999.
testlimit = 50
eps       = 1e-10

#-------------------Base directories
repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
data_dir = os.path.join(repo_dir, 'data')
obs_data_dir = os.path.join(data_dir, 'observational_data')
slurm_temp_dir = os.path.join(data_dir, 'slurm_templates')

#-------------------Solar constants
Lbolsun = 3.826e33  # erg/s
Msun    = 1.989e30  # kg
zsun = 0.0134       # Asplund 2009
zsunK20 = 0.014     # Kashino 2020
ohsun = 8.69        # Allende Prieto 2001 and Asplund 2009 (12 + log10(O/H))sun
h_nul = 13.6        # Lyman limit h*nu(eV)
parsec  = 3.085677581491367e+16 # m

#--------------------------------------------
#   Conversion factors:
#--------------------------------------------
kg_to_Msun= 1./Msun
Mpc_to_cm = parsec*1e8
yr_to_s   = 365.*24.*60.*60.
kilo      = 1000.0
mega      = 1000000.0
giga      = 1000000000.0
J2erg     = 1e7
eV        = 1.602e-19   #J=kg*m^2*s**-2
#--------------------------------------------
G    = 6.6743e-11          # Gravitational constant, Nm^2/kg^2=m^3/kg/s^2
mp   = 1.67e-27            # Proton mass, kg
c    = 2.998e8             # Light velocity, m/s
h    = 6.62607015e-34      # Planck constant, Js
kB   = 1.380649e-23        # Boltzmann constant, J/K

G_Ms = G*Msun/(kilo*kilo*parsec*mega) # 4.301e-9 km^2*Mpc/Msun/s^-2 
c_cm = c*100.                         # Light velocity, cm/s
h_erg= h*J2erg                        # Planck constant, erg s
kB_Ms= kB/(Msun*(parsec*mega)**2)     # 7.29e-99 Mpc^2*Msun/s^2/K 
#--------------------------------------------
sigma_1Dprobs = [0.682689492137086,    # 1 sigma
                 0.954499736103642,    # 2 sigma
                 0.997300203936740,    # 3 sigma
                 0.999936657516334,    # 4 sigma
                 0.999999426696856,    # 5 sigma
                 0.999999998026825]    # 6 sigma
sigma_2Dprobs = [0.3935,    # 1 sigma
                 0.6321,    # 2 sigma
                 0.7769,    # 3 sigma
                 0.8647,    # 4 sigma
                 0.9179,    # 5 sigma
                 0.9502,    # 6 sigma
                 0.9698,    # 7 sigma
                 0.9817,    # 8 sigma
                 0.9889,    # 9 sigma
                 0.9933]    # 10 sigma
#--------------------------------------------
#   Possible options and models:
#--------------------------------------------
inputformats = ['txt','hdf5']

zeq = ['tremonti2004','tremonti2004b','leblanc']

model_nH_sfr = ['kashino20']
model_U_sfr  = ['kashino20', 'orsi14']

model_spec_agn = 'feltre16'
alpha_NLR_feltre16 = -1.7
xid_NLR_feltre16 = 0.5

model_U_agn    = ['panuzzo03']

photmods = ['gutkin16', 'feltre16']
mod_lim = {'gutkin16': r"nebular_data/gutkin16_tables/limits_gutkin.txt",
           'feltre16': r"nebular_data/feltre16_tables/limits_feltre.txt"}

#--------------------------------------------
#   Orsi et. al. 2014
#--------------------------------------------
Z0_orsi = 0.012
q0_orsi = 2.8e7 # cm/s
gamma_orsi = 1.3
#--------------------------------------------

#--------------------------------------------------------------
#   IMF transformations (Tables B1 and B2 in Lacey et al. 2016) 
#--------------------------------------------------------------
# log10(M1) = log10(IMF_M2/IMF_M1) + log10(M2) 
IMF_M = {'Kennicut': 1, 'Salpeter': 0.47, 'Kroupa': 0.74, 'Chabrier': 0.81, 
            'Baldry&Glazebrook': 0.85, 'Top-heavy': 1.11}

# # log10(SFR1) = log10(IMF_SFR2/IMF_SFR1) + log10(SFR2)
IMF_SFR = {'Kennicut': 1, 'Salpeter': 0.79, 'Kroupa': 1.19,
           'Chabrier': 1.26,'Baldry&Glazebrook': 1.56,
           'Top-heavy': 1.89}
IMF_SFRins = {'Kennicut': 1, 'Salpeter': 0.94, 'Kroupa': 1.49,
              'Chabrier': 1.57, 'Baldry&Glazebrook': 2.26,
              'Top-heavy': 3.13}

phot_to_sfr_kenn = 9.85e52 # phot/s

# FOR CONVERSION FROM LYMANN CONTINUUM PHOTONS TO SFR
# It is assumed that a SFR of 1 Msun/yr produces 9.85 · 10^52 photons/s for Kennicut IMF.
# Reference: Chomiuk & Povich 2011, pag. 2: "According to Kennicutt et al. (1994) and Kennicutt
# (1998a), a SFR of 1 M⊙ yr−1 produces a Lyman continuum photon rate Nc = 9.26 × 1052 photon s−1
# for the Salpeter (1955) IMF (assuming a mass range of 0.1–100 Msun)."
# Reescaled to Kennicut, it gives our number.

#-------------------------------------------
#    Scalelength:
#-------------------------------------------
re2hr    = 1/1.68 # Cole+2000, Leroy+2021
r502re   = 1.     # Lima Neto+1999, Wolf+2010, Huang+2017
rvir2r50 = 0.03   # Huang+2017, Yang+2025

#-------------------------------------------
#    Attenuation:
#-------------------------------------------
attmods = ['ratios', 'favole20']
Rv = 3.1
costheta = 0.6
albedo = 0.8

#-------------------------------------------
#    AGNs:
#-------------------------------------------

Lagn_inputs = ['Lagn', 'acc_rate', 'acc_rates', 'radio_mode', 'quasar_mode', 'complete']

# Griffin et. al 2019:
eta_acc_eff = 0.1 # Nominal efficiency of the accretion. eq 10 in Griffin+2019
alpha_adaf = 0.1 # Viscosity parameter for ADAFs
alpha_td = 0.1 # Viscosity parameter for TDs
lambda_adaf = 0.2 # Fraction of viscous energy transferred to electrons in ADAF
acc_rate_crit_adaf = 0.01 # Boundary between thin disc and ADAF accretion (in terms of ratio between accretion rate and eddington's)
eta_edd = 4 # Super-Eddington suppression factor
fq = 10 # Ratio of lifetime of AGN episode to bulge dynamical timescale
fbh = 0.005 # Fraction of the mass of stars formed in a starburst accreted onto a black hole

beta = 1 - alpha_adaf/0.55 #  Ratio of gas pressure to total pressure. See 3.1 in Griffin+2019.
acc_rate_crit_visc = 0.001*(lambda_adaf/0.0005)*((1-beta)/beta)*alpha_adaf**2 
# Boundary between the two adaf regimes
spin_bh = 0.67 # 0, 0.3, 0.5, 0.67, 0.9, 0.95, 1

# Typical value for accretion efficiency (e.g. Griffin+2019)
e_r_agn = 0.1

# Fit of GP20 data to equation 1 in Henriques et al. 2016: 
kagn = 5.44e-4 
kagn_exp = 0.597

# Typical temperature of ionising regions (Hirschmann+2017):
temp_ionising = 10000  # K

# Typical values from Osterbrock and Ferland, 2006 book:
nH_sfr_cm3 =   100   # cm^-3
nH_NLR_cm3 =  1000   # cm^-3
epsilon_NLR = 0.01
radius_NLR = 0.001   # Mpc
alphaB = {5000: 4.54e-13, 10000: 2.59e-13, 20000: 1.43e-13}  #cm^3/s

#------------------------------------------
nH_bins = {
    "gutkin16" : np.array([10, 100, 1000, 10000]),
    "feltre16" : np.array([100, 1000, 10000])
}

lus_bins = {
    "gutkin16" : np.array([-4., -3.5, -3., -2.5, -2., -1.5, -1.]),
    "feltre16" : np.array([-5., -4.5, -4., -3.5, -3., -2.5, -2.,
                           -1.5, -1.])
}

zmet_str = {
    "gutkin16" : np.array(['0001','0002','0005','001','002','004','006',
                           '008','010','014','017','020','030','040']),
    "feltre16" : np.array(['0001','0002','0005','001','002','004',
                           '006','008','014','017','020','030','040',
                           '050','060','070'])
}

zmet_str_reduced = {
    "gutkin16" : np.array(['0001','002','014','030']),
}

line_names = {
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

line_wavelength = {
    "gutkin16" : np.array([3727,4861,4959,5007,6548,6563,6584,
                           6717,6731,1240,1548,1551,1640,1661,
                           1666,1883,1888,1908]),
    "feltre16" : np.array([3727,4861,4959,5007,6300,6548,6563,6584,
                           6717,6731,1240,1548,1551,1640,1661,
                           1666,1883,1888,1907,1910])
    }

# Limits in h*nu for the piecewise AGN spectral approximation
agn_spec_limits = {
    "feltre16" : np.array([0.1,5,1240])
} 

def coef_att_line_model_func(z=0):
    line_att_coef = line_att_coef_func(z)

    coef_att_line_model = {
        "gutkin16" : np.array([line_att_coef[line] for line in line_names["gutkin16"]]),
        "feltre16" : np.array([line_att_coef[line] for line in line_names["feltre16"]])
        }

    return coef_att_line_model

