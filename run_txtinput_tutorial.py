"""
This is a modular python code to provide model spectral emission 
lines, both from star-forming regions and narrow-line regions of AGNs.
The input of the code are global galactic properties. 

The intrinsic luminosities can be passed through an attenuation model 
to also get the predicted attenuated luminosities.

@authors: viogp, expox7
"""

import gne.gne_const as const
from gne.gne import gne, gne_att, gne_flux
from gne.gne_plots import make_testplots

### RUN the code with the given parameters and/or make plots
testing = False           # If True: use only the first 50 elements
get_emission_lines = True # Obtain nebular emission lines
get_attenuation = True
get_flux = True
plot_tests = True

# Calculate emission from AGNs: AGN = True
AGN = True

###############################################################
### OUTPUT FILES: Default output path is output/
outpath = None

###############################################################
### INPUT FILES: given as a root, ending and number of subvolumes
# Input files are expected to have, AT LEAST:
# Stellar mass (M*) of the galaxy (or disc, SF burst, buldge, etc).
# Star formation rate (SFR) or 12+log(O/H)
# Mean metallicity of the cold gas (Z).
root = 'data/example_data/iz61/ivol'
endf   = '/ex.txt'
subvols = 2

# Redshifts, cosmology and volume of the simulation
redshift = 0.
snapshot = 61
h0     = 0.704
omega0 = 0.307
omegab = 0.0482
lambda0= 0.693
vol    = 12500000.0
mp     = 9.35e8       

### INPUT FORMAT ('txt' for text files; 'hdf5' for HDF5 files)
inputformat = 'txt'

### UNITS: 
# units_h0=False if input units [Mass]=Msun, [Radius]=Mpc (default)
# units_h0=True  if input units [Mass]=Msun/h, [Radius]=Mpc/h
units_h0=True
# units_Gyr=False if input units [SFR,Mdot]=[Mass]/yr (default)
# units_Gyr=True  if input units [SFR,Mdot]=[Mass]/Gyr 
units_Gyr=True 
# Type of input units for luminosity, per component:
# 0: input units [L]=erg/s  (default);
# 1: input units [L]=1e40 h^-2 erg/s
# 2: input units [L]=1e40 erg/s
units_L=1

####################################################
############  Emission from SF regions #############
####################################################

# All available models can be seen in gne_const module.
# NEBULAR model connecting global properties to ionising properties:
# nH: number density of Hydrogen (or electrons); U: ionising parameter
model_nH_sfr='kashino20'
model_U_sfr='kashino20' #'orsi14'   
# PHOTOIONIZATION model for SF regions to get line luminosities
photmod_sfr='gutkin16'

### INPUT PARAMETERS
# m_sfr_z has the location in the input files of the three mandatory parameters:
# M*(units), SFR or 12+log(O/H), and Zgas. 
# m_sfr_z is a list of lists with either the column number
# for each parameters or the name of the HDF5 variable.
# Each list correspond to a different component: 
# m_sfr_z = [[mstar_disk,SFR_disk,Zgas_disk],[mstar_stb,SFR_stb,Zgas_stb]]
# For a single component: m_sfr_z = [[Mstellar,SFR,Zgas]]
#m_sfr_z = [[20,15,6]]
m_sfr_z = [[20,15,6],[18,17,5]]

# mtot2mdisk is True if the stellar mass of discs is calculated 
# from the total and buldge values (False by default)
# mtot2mdisk = True; cols = [[M,SFR,Z],[M_bulge,SFR_bulge,Z_bulge]]
# mtot2mdisk = False; cols = [[M_disk,SFR_disk,Z_disk],[M_bulge,SFR_bulge,Z_bulge]]        
mtot2mdisk = False

# inoh True if the gas metallicity input as log(O/H)+12
#      False if Zgas = MZcold/Mcold (Default option)
inoh = False

### INITIAL MASS FUNCTIONs
# Specify the assumed IMFs for each galaxy component in the input data.
# Example for two components: IMF = ['Kennicut','Kennicut']
IMF = ['Kennicut','Kennicut']

####################################################
#####  Emission from AGN narrow line regions #######
####################################################
# PHOTOIONIZATION model for AGN NLR to get line luminosities
photmod_agn = 'feltre16'

# Columns to read either the central or global metallicity
# If several components are given, they will be added
Zgas_NLR = [6,5]
# Z_correct_grad 
#    False (default) if the central gas metallicity has been provided
#    True to correct a global metallicity with the gradients from Belfiore+2017
Z_correct_grad = True

# Connecting global properties to AGN NLR characteristics:
# Model to calculate the ionising parameter, U
model_U_agn    = 'panuzzo03'

# Panuzzo's model requires the calculation of the filling factor
# epsilon(Mgas, Scalelength, n_NLR, T_NLR, r_NLR)
# n_NLR, T_NLR and r_NLR are taken as constants.
# If mgas_r = None, a fixed volume-filling factor is assumed, otherwise
# mgas_r is a list of lists with either the column number
# for each parameters or the name of the HDF5 variable.
# Each list can correspond to a different component:
# mgas_r = [[mgas_comp1,R_comp1],...]
#mgas_r = [None]
mgas_r = [[11,23],[12,21]]

# If mgas_r given, the type of component needs to be specified
# mgasr_type = ['disc', 'bulge' or None]
mgasr_type = ['disc','bulge']

# Type of radius input, per component:
# 0: scalelength;
# 1: effective radius, Re
# 2: half-mass/light radius, R50 (Re=r502re*R50 with a default r502re=1) 
# 3: radius of the galaxy or host halo
r_type = [1,1]

# spec: model for the spectral distribution of the AGN
model_spec_agn = 'feltre16'

# The AGNs bolometric luminosity, Lagn, is needed.
# This value can be either firectly input or calculated.
# The way of obtaining Lagn is indicated in Lagn_inputs.
# The calcultions require different black hole (BH) parameters.
# Lagn_inputs='Lagn' if Lagn in input
#            in erg/s,h^-2erg/s,1e40erg/s,1e40(h^-2)erg/s
#            Lagn_params=[Lagn, Mbh] 
# Lagn_inputs='Mdot_hh' for a calculation from
#            the mass accretion rate of the BH, Mdot,
#            the BH mass, Mbh,
#            and, as an optional input, the BH spin, Mspin. 
#            Lagn_params=[Mdot,Mbh] or [Mdot,Mbh,Mspin]
# Lagn_inputs='Mdot_stb_hh' for a calculation from
#            the mass accretion rate during the last starburst, Mdot_stb,
#            the hot halo or radio mass accretion, Mdot_hh,
#            the BH mass, Mbh,
#            and, as an optional input, the BH spin, Mspin. 
#            Lagn_params=[Mdot_stb,Mdot_hh,Mbh] or [Mdot_stb,Mdot_hh,Mbh,Mspin]
# Lagn_inputs='radio_mode' for a calculation from
#            the mass of the hot gas, Mhot,
#            the BH mass, Mbh,
#            and, as an optional input, the BH spin, Mspin. 
#            Lagn_params=[Mhot,Mbh] or [Mhot,Mbh,Mspin]
# Lagn_inputs='quasar_mode' for a calculation from
#            the mass of the bulge, Mbulge,
#            the half-mass radius of the bulge, rbulge,
#            the circular velocity of the bulge, vbulge,
#            the BH mass, Mbh,
#            and, as an optional input, the BH spin, Mspin. 
#            Lagn_params=[Mbulge,rbulge,vbulge,Mbh,(Mspin)]
# Lagn_inputs='complete' for a calculation from
#            the mass of the bulge, Mbulge,
#            the half-mass radius of the bulge, rbulge,
#            the circular velocity of the bulge, vbulge,
#            the mass of the hot gas, Mg,
#            the BH mass, Mbh,
#            and, as an optional input, the BH spin, Mspin. 
#            Lagn_params=[Mbulge,rbulge,vbulge,Mhot,Mbh,(Mspin)]
Lagn_inputs = 'Lagn'; Lagn_params=[0,1]
#Lagn_inputs = 'Mdot_hh'; Lagn_params=[16,8,21]
#Lagn_inputs = 'Mdot_stb_hh'; Lagn_params=[15,16,8,21]
#Lagn_inputs = 'radio_mode'; Lagn_params=[9,8]
#Lagn_inputs = 'quasar_mode'; Lagn_params=[25,12,14,21]
#Lagn_inputs = 'complete'; Lagn_params=[25,12,14,9,21]

####################################################
########  Redshift evolution parameters  ###########
####################################################

### HIGH REDSHIFT CORRECTION ###
# Empirical relationships to connect global galaxy properties and nebular
    # properties are often derived from local galaxies. get_emission_lines has
    # a way of evolving the filling factor with redshift. If this correction is to be used,
    # a fixed number of files is needed equal to that at z=0.
    # If local relations are to be used: infiles_z0 = [None]
root_z0 = None

####################################################
##########       Dust attenuation      #############
####################################################

# Continuum and line attenuation calculation. If this option is selected 
    # the output file will have intrinsic AND attenuated values of
    # luminosity for the emission lines.
# att=True to calculate the dust attenuation; False, otherwise
att = False
    
# To use Cardelli's law (following Favole et. al. 2020):
    # attmod = 'cardelli89' (default)
    # att_params = [half-mass radius, cold gas mass, cold gas metallicity]
# To use already available attenuation coefficients: attmod = 'ratios'
    # att_params in this case has the location of the attenuation coefficients
    # for each line for which attenuation is going to be calculated.
    # This mode requieres an extra variable, att_ratio_lines, with the names
    # of the lines corresponding to the coefficients listed in att_params.
# Example:
    # attmod = 'ratios'
    # att_params = [31,32,33,34,35,36,36]
    # att_ratio_lines = ['Halpha','Hbeta','NII6584','OII3727','OIII5007','SII6717','SII6731'] 

attmod='cardelli89'
att_params= []

####################################################
##########      Other calculations     #############
####################################################
# Include other parameters in the output files
extra_params_names = ['type','mh','xgal','ygal','zgal',
                      'vxgal','vygal','vzgal',
                      'magK','magR','M_SMBH']
extra_params_labels = extra_params_names
extra_params = [24,13,29,30,31,26,27,28,10,9,1]

### SELECTION CRITERIA ###
# Cuts can be made on the input file
# In this example, location 7 correspond to the halo mass.
# The dark matter particles of the simulations has a mass of 9.35e8 Msun/h

# Paramter to impose cuts
cutcols = [13]
# List of minimum values. None for no inferior limit.
mincuts = [21*mp]
# List of maximum values. None for no superior limit.
maxcuts = [None]


##################################################################
#############    Run the code and/or make plots   ################
##################################################################

verbose = True
for ivol in range(subvols):
    infile = root+str(ivol)+endf

    infile_z0 = root_z0
    if root_z0 is not None:
        infile_z0 = root_z0+str(ivol)+endf

    if get_emission_lines:  # Obtain nebular emission lines
        gne(infile,redshift,snapshot,h0,omega0,omegab,lambda0,vol,mp,
            inputformat=inputformat,outpath=outpath,
            units_h0=units_h0,units_Gyr=units_Gyr,units_L=units_L,
            model_nH_sfr=model_nH_sfr, model_U_sfr=model_U_sfr,
            photmod_sfr=photmod_sfr,
            m_sfr_z=m_sfr_z,mtot2mdisk=mtot2mdisk,
            inoh=inoh,IMF = IMF,
            AGN=AGN,photmod_agn=photmod_agn,
            Zgas_NLR=Zgas_NLR,Z_correct_grad=Z_correct_grad,
            model_U_agn=model_U_agn,           
            mgas_r_agn=mgas_r,mgasr_type_agn=mgasr_type,r_type_agn=r_type,
            model_spec_agn=model_spec_agn,
            Lagn_inputs=Lagn_inputs, Lagn_params=Lagn_params,
            att=att, line_att=line_att, attmod=attmod,
            att_params=att_params,
            infile_z0=infile_z0, 
            extra_params=extra_params,
            extra_params_names=extra_params_names,
            extra_params_labels=extra_params_labels,
            cutcols=cutcols, mincuts=mincuts, maxcuts=maxcuts,
            testing=testing,verbose=verbose)

if plot_tests:  # Make test plots
    make_testplots(root,endf,snapshot,subvols=subvols,
                   gridplots=False,
                   outpath=outpath,verbose=verbose)
