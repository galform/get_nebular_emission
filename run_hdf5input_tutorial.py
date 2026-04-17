"""
This is a modular python code to provide model spectral emission 
lines, both from star-forming regions and narrow-line regions of AGNs.
The input of the code are global galactic properties. 

The intrinsic luminosities can be passed through an attenuation model 
to also get the predicted attenuated luminosities.

@authors: viogp
"""

import gne.gne_const as const
from gne.gne import gne
from gne.gne_att import gne_att
from gne.gne_flux import gne_flux
from gne.gne_plots import make_testplots
import os, h5py

verbose = True
### RUN the code with the given parameters and/or make plots
testing = False            # If True: use only the first 50 elements
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
subvols = 2
root = 'data/example_data/iz61/ivol'
endf   = 'ex.hdf5'

### INPUT FORMAT ('txt' for text files; 'hdf5' for HDF5 files)
inputformat = 'hdf5'

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
m_sfr_z = [['data/mstars_disk','data/mstardot','data/Zgas_disc'],
           ['data/mstars_bulge','data/mstardot_burst','data/Zgas_bst']]

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
Zgas_NLR = ['data/Zgas_bst','data/Zgas_disc']
# Z_correct_grad 
#    False (default) if the central gas metallicity has been provided
#    True to correct a global metallicity with the gradients from Belfiore+2017
Z_correct_grad = True

# Connecting global properties to AGN NLR characteristics:
# Model to calculate the ionising parameter, U
model_U_agn    = 'panuzzo03'

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
Lagn_inputs = 'Lagn'; Lagn_params=['data/Lbol_AGN','data/mstars_bulge']

###################################################################
########  Filling factor and Cardelli's law parameters  ###########
###################################################################
# Panuzzo's model requires the calculation of the filling factor
# epsilon(Mgas, Scalelength, n_NLR, T_NLR, r_NLR)
# n_NLR, T_NLR and r_NLR are taken as constants.
# If mgas_r = None, a fixed volume-filling factor is assumed, otherwise
# mgas_r is a list of lists with either the column number
# for each parameters or the name of the HDF5 variable.
# Each list can correspond to a different component:
# mgas_r = None
# mgas_r = [[mgas_comp1,R_comp1],...]
mgas_r = [['data/mcold','data/rdisk'],
          ['data/mcold_burst','data/rbulge']]

# If mgas_r given, the type of component needs to be specified
# mgasr_type = 'disc', 'bulge' or None
mgasr_type = ['disc','bulge']

# Type of radius input, per component:
# 0: scalelength;
# 1: effective radius, Re
# 2: half-mass/light radius, R50 (Re=r502re*R50 with a default r502re=1) 
# 3: radius of the galaxy or host halo
r_type = [2,2]

####################################################
##########       Dust attenuation      #############
####################################################
#WARNING: Att w ratios not checked (no adequate ex. data)

# Dust-attenuated luminosities are calculated if get_attenuation=True
# line_att=True to apply the extra line attenuation from Saito+2021
line_att = True

# Available dust attenuation models
# 'favole20' (default)
#    The calculation follows Favole et. al. 2020 and requires
#    the above parameters: mgas_r, mgasr_type and r_type
#    If None is passed, default parameters will be used
attmod = 'favole20'
att_config = {'Rv': None, 'albedo': None, 'costheta': None} 
# 'ratios'
#    The calculation uses already available attenuation coefficients.
#    att_ratios should contain the location of these coefficients, and
#    the names of the lines with available ratios are in att_rlines.
#    Example for
#attmod = 'ratios'
#att_config = ['Halpha', 'Hbeta', 'NII6583', 'OII3727', 'OIII5007', 'SII6716']

####################################################
########  Redshift evolution parameters  ###########
####################################################
# WARNING: Evolution calculation has not been fully tested

### HIGH REDSHIFT CORRECTION ###
# Empirical relationships to connect global galaxy properties and nebular
    # properties are often derived from local galaxies. get_emission_lines has
    # a way of evolving the filling factor with redshift. If this correction is to be used,
    # a fixed number of files is needed equal to that at z=0.
    # If local relations are to be used: infiles_z0 = [None]
root_z0 = None

####################################################
##########      Other calculations     #############
####################################################
# Include other parameters in the output files
# WARNING: magK and magR are the dataset names used
#          for selections in plots (optional)
extra_params_names = ['type','mh','xgal','ygal','zgal',
                      'vxgal','vygal','vzgal',
                      'magK','magR','M_SMBH','index']
extra_params = ['data/type','data/mhhalo',
                'data/xgal','data/ygal','data/zgal',
                'data/vxgal','data/vygal','data/vzgal',
                'data/mag_UKIRT-K_o_tot_ext',
                'data/mag_SDSSz0.1-r_o_tot_ext',
                'data/M_SMBH','data/index']
if attmod == 'ratios':
    for line in att_config:
        extra_params_names.append('ratio_'+line)
        extra_params.append('data/ratio_'+line)
extra_params_labels = extra_params_names

### SELECTION CRITERIA ###
# Cuts can be made on the input file
# In this example, location 7 correspond to the halo mass.
# The dark matter particles of the simulations has a mass of 9.35e8 Msun/h

# Paramter to impose cuts
cutcols = ['data/mhhalo']
# List of minimum values. None for no inferior limit.
mincuts = [21*9.35e8]
# List of maximum values. None for no superior limit.
maxcuts = [None]

##################################################################
#############    Run the code and/or make plots   ################
##################################################################
list_subvols = subvols
if isinstance(subvols, int):
    list_subvols = list(range(subvols))

for ivol in list_subvols:
    infile = os.path.join(root+str(ivol),endf)

    infile_z0 = root_z0
    if root_z0 is not None:
        infile_z0 = os.path.join(root_z0+str(ivol),endf) 

    # Get the redshift, cosmology and volume of the model galaxies
    f = h5py.File(infile) 
    header = f['header'] #; print(list(header.attrs.keys()))
    redshift = header.attrs['redshift']
    snapshot = header.attrs['snapnum']
    boxside = header.attrs['bside_Mpch']
    h0 = header.attrs['h0']
    omega0 = header.attrs['omega0']
    omegab = header.attrs['omegab']
    lambda0 = header.attrs['lambda0']
    mp = header.attrs['mp_Msunh']
    try:
        p = header.attrs['percentage']/100.
    except:
        p = 1
    f.close()
    effvol = p*boxside**3

    if get_emission_lines:  
        # Obtain nebular emission lines
        gne(infile,redshift,snapshot,h0,omega0,omegab,lambda0,
            mp,boxside,effvol,
            inputformat=inputformat,outpath=outpath,
            units_h0=units_h0,units_Gyr=units_Gyr,units_L=units_L,
            model_nH_sfr=model_nH_sfr, model_U_sfr=model_U_sfr,
            photmod_sfr=photmod_sfr,
            m_sfr_z=m_sfr_z,mtot2mdisk=mtot2mdisk,
            inoh=inoh,IMF = IMF,
            AGN=AGN,photmod_agn=photmod_agn,
            Zgas_NLR=Zgas_NLR,Z_correct_grad=Z_correct_grad,
            model_U_agn=model_U_agn,
            mgas_r=mgas_r,mgasr_type=mgasr_type,r_type=r_type,
            model_spec_agn=model_spec_agn,
            Lagn_inputs=Lagn_inputs, Lagn_params=Lagn_params,
            infile_z0=infile_z0, 
            extra_params=extra_params,
            extra_params_names=extra_params_names,
            extra_params_labels=extra_params_labels,
            cutcols=cutcols, mincuts=mincuts, maxcuts=maxcuts,
            testing=testing,verbose=verbose)

    if get_attenuation: # Obtain dust-attenuated luminosities
        gne_att(infile,outpath=outpath,
                attmod=attmod,line_att=line_att,
                att_config=att_config,verbose=verbose)

    if get_flux: # Calculate fluxes from luminosities
        gne_flux(infile,outpath=outpath,
                 line_names=['Halpha','Hbeta','NII6584','OIII5007'],
                 verbose=verbose)

if plot_tests:  # Make test plots
    make_testplots(snapshot,endf,outpath=outpath,
                   subvols=list_subvols,
                   gridplots=False,verbose=verbose)
