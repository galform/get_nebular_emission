"""
This is a modular python code to provide model spectral emission 
lines, both from star-forming regions and narrow-line regions of AGNs.
The input of the code are global galactic properties. 

The intrinsic luminosities can be passed through an attenuation model 
to also get the predicted attenuated luminosities.

@authors: expox7, viogp
"""

import src.gne_const as const
from src.gne import gne
from src.gne_plots import make_testplots

### RUN the code with the given parameters and/or make plots
run_code = True
make_plots = True

# Calculate emission from AGNs: AGN = True
AGN = True

###############################################################
### INPUT FILE(S) and redshift
# Input files are expected to have, AT LEAST:
# Stellar mass (M*) of the galaxy (or disc or buldge).
# Star formation rate (SFR) OR magnitude of Lyman Continuum photons (m_LC).
# Mean metallicity of the cold gas (Z).
infiles = ['src/example_data/GP20_62.5kpc_z0_example.txt']
redshifts = [0.15]

# Cosmology and volume of the simulation
h0     = 0.704
omega0 = 0.307
omegab = 0.0482
lambda0= 0.693
vol    = pow(62.5,3) #Mpc/h

### INPUT FORMAT
# If your input files are text files (.txt, .dat, .csv...): inputformat = 'txt'
# If your input files are HDF5 files: inputformat = 'hdf5'
inputformat = 'txt'

####################################################
############  Emission from SF regions #############
####################################################

### NEBULAR AND PHOTOIONIZATION MODELS for SF regions
# All available models can be seen in gne_const module.
# Model that connects global galactic properties with nebular parameters.
unemod_sfr='kashino20'
# Photoionization model to get line luminosities from nebular parameters.
photmod_sfr='gutkin16'

### INPUT PARAMETERS
# m_sfr_z has the location in the input files of the three mandatory parameters:
# M*(units), SFR or m_LC and Zgas. 
# m_sfr_z is a list of lists with either the column number
# for each parameters or the name of the HDF5 variable.
# Each list correspond to a different component: 
# m_sfr_z = [[M_disk,SFR_disk,Z_disk],[M_bulge,SFR_bulge,Z_bulge]]
# For a single component: m_sfr_z = [[M*,SFR,Zgas]]
# For a HDF5 input file: m_sfr_z = [['Mstellar','SFR','Zgas']]

#m_sfr_z = [[0,2,4]]
m_sfr_z = [[0,2,4],[1,3,5]]

# mtot2mdisk is True if the stellar mass of discs is calculated 
# from the total and buldge values (False by default)
# mtot2mdisk = True; cols = [[M,SFR,Z],[M_bulge,SFR_bulge,Z_bulge]]
# mtot2mdisk = False; cols = [[M_disk,SFR_disk,Z_disk],[M_bulge,SFR_bulge,Z_bulge]]        
mtot2mdisk = False

# LC2sfr is True when Lyman Continuum photons are given  instead of the SFR
# LC2sfr = True; cols = [[M,m_LC,Z]]
# LC2sfr = False; cols = [[M,SFR,Z]] (Default option)      
LC2sfr = False

# inoh True if the gas metallicity input as log(O/H)+12
#      False if Zgas = MZcold/Mcold (Default option)
inoh = False

### INITIAL MASS FUNCTIONs ###
# Specify the assumed IMFs for each galaxy component in the input data.
# Example for two components: IMF = ['Kennicut','Kennicut']
IMF = ['Kennicut','Kennicut']

####################################################
#####  Emission from AGN narrow line regions #######
####################################################

### NEBULAR AND PHOTOIONIZATION MODELS
# All available models can be seen in gne_const module.
# Model that connects global galactic properties with nebular parameters.
unemod_agn='panuzzo03'
# Photoionization model to get line luminosities from nebular parameters.
photmod_agn='feltre16'

# Parameters for calculating emission from AGN NLR:
# Cold gas mass (Mg).
# Baryonic half-mass radius (Rhm).
# For two disk and bulge:epsilon_params = [Mg_disk, Rhm_disk, Mg_bulge, Rhm_bulge]
epsilon_params=[6,11,19,12]
    
# Second, the AGNs bolometric luminosity is needed. This value can be in the
    # input file or it can be estimated from other paremeters. To indicate
    # how are you going to get it you use AGNinputs and then put the location
    # of the corresponding parameters in Lagn_params:
        # Lagn: Put the luminosity value directly.
            # Parameters:
                # Bolometric luminosity of the AGN: Lagn
                # Bulge stellar mass: Ms_bulge
            # Lagn_params = [Lagn, Ms_bulge]
        # acc_rate: Put the mass accretion rate of the black hole.
            # Parameters:
                # Mass accretion rate: M_dot
                # Black hole mass: Mbh
                # Bulge stellar mass: Ms_bulge
            # Lagn_params = [M_dot, Mbh, Ms_bulge]
        # acc_rates: Put the mass accretion rates for the quasar mode and
            # the radio mode.
            # Parameters:
                # Mass accretion rate of the quasar mode: M_dot_stb
                # Mass accretion rate of the radio mode: M_dot_radio
                # Black hole mass: Mbh
                # Bulge stellar mass: Ms_bulge
            # Lagn_params = [M_dot_stb, Mdot_radio, Mbh, Ms_bulge]
        # radio_mode: Put the parameters needed to calculate the mass accretion
            # rate of the radio mode and use it to calculate Lagn.
            # Parameters:
                # Mass of the hot gas: Mhot
                # Black hole mass: Mbh
                # Bulge stellar mass: Ms_bulge
            # Lagn_params = [Mhot, Mbh, Ms_bulge]
AGNinputs = 'Lagn'
Lagn_params=[17,21]

# AGN emission calculation is done assuming that the available metallicity 
    # value is the one corresponding to the NLR, i.e. the metallicity
    # around the center of the galaxy. 
# If Z_central_cor is True, the code estimates the value of metallicity around 
    # the center of the galaxy from the mean value, following the function
    # Zagn from the gne_une module.
# If Z_central_cor is False, it is assumed that the metallicity is 
    # approximatelly uniform in the galaxy.
Z_central_cor=True

####################################################
########  Redshift evolution parameters  ###########
####################################################

### HIGH REDSHIFT CORRECTION ###
# Empirical relationships to connect global galaxy properties and nebular
    # properties are often derived from local galaxies. get_emission_lines has
    # a way of evolving the filling factor with redshift. If this correction is to be used,
    # a fixed number of files is needed equal to that at z=0.
    # If local relations are to be used: infile_z0 = [None]
infile_z0 = [None]


####################################################
##########       Dust attenuation      #############
####################################################

# Continuum and line attenuation calculation. If this option is selected 
    # the output file will have intrinsic AND attenuated values of
    # luminosity for the emission lines.
    
# Calculate attenuation: att = True
# Don't calculate attenuation: att = False
att=True

# attmod defines the attenuation model the code is going to use. 
# Each model requires different parameters. The code also provides the
    # possibility of using attenuation coefficients calculated externally and
    # inserted in the input file.
# att_params is a list with the location of the required parameters. Right now
    # the code only has one attenuation model available, the attenuation curve
    # from Cardelli et. al. 1989 used with the methodology of Favole et. al. 2020.
    
# To use cardelli's attenuation curve: attmod = 'cardelli89'.
    # Parameters:
        # Baryonic half-mass radius (Rhm).
        # Cold gas mass (Mg).
        # Cold gas metallicity (Z).
    # att_params = [Rhm, Mg, Z]

# To use already available attenuation coefficients: attmod = 'ratios'
    # att_params in this case has the location of the attenuation coefficients
    # for each line for which attenuation is going to be calculated.
    # This mode requieres an extra variable, att_ratio_lines, with the names
    # of the lines corresponding to the coefficients listed in att_params.
# Example:
    # Suppose we want to calculate the attenuation por Halpha and OII3727. Our 
    # input file is a text file with the attenuation coefficients of 
    # those two lines in the sixth and eighth columns.
    # So:
        # attmod = 'ratios'
        # att_params = [5,7]
        # att_ratio_lines = ['Halpha','OII3727']

attmod='ratios'
att_params=[31,32,33,34,35,36,36] # It would be [11,6,4] if you want to use Cardelli.
att_ratio_lines=['Halpha','Hbeta','NII6584','OII3727','OIII5007','SII6717','SII6731']

####################################################
##########      Other calculations     #############
####################################################

# Flux calculation.
flux=True

# Include other parameters in the output files
extra_params_names = ['Type','Mbh','Mhalo','Ms_bulge','m_K','m_R','m_R_SDSS','m_I',
                      'Mdot_stb','Mdot_hh','Mhot','Lagn']
extra_params_labels = ['Type of halo (central = 0)',
                       r'Black hole mass ($M_\odot \ h^{-1}$)',
                       r'Halo mass ($M_\odot \ h^{-1}$)',
                       r'Stellar mass of bulge ($M_\odot \ h^{-1}$)',
                       'K band (Apparent magnitude, attenuated)',
                       'R band (Apparent magnitude, attenuated)',
                       'R band SDSS (Apparent magnitude, attenuated)',
                       'I band (Apparent magnitude, attenuated)',
                       'SMBH mass accretion rate (starburst mode)',
                       'SMBH mass accretion rate (hot halo mode)',
                       r'Hot gas mass ($M_\odot \ h^{-1}$)',
                       r'AGN bolometric luminosity (erg $s^{-1}$)']
extra_params = [30,8,7,21,25,27,18,29,15,16,9,17]


### SELECTION CRITERIA ###
# Cuts can be made on the input file
# In this example, location 7 correspond to the stellar mass.
# The dark matter particles of the simulations has a mass of 9.35e8 Msun/h
cutcols = [7]
# List of minimum values. None for no inferior limit.
mincuts = [20*9.35e8]
# List of maximum values. None for no superior limit.
maxcuts = [None]


##################################################################
#############    Run the code and or make plots   ################
##################################################################

for ii, infile in enumerate(infiles):
    zz = redshifts[ii]

    if run_code:
        gne(infile, zz, m_sfr_z,
            h0,omega0,omegab,lambda0,vol,
            infile_z0=infile_z0, inputformat=inputformat,
            cutcols=cutcols, mincuts=mincuts, maxcuts=maxcuts,
            att=att,
            att_params=att_params, att_ratio_lines=att_ratio_lines,
            flux=flux,
            IMF = IMF,
            inoh=inoh, mtot2mdisk=mtot2mdisk, LC2sfr=LC2sfr,
            AGN=AGN,AGNinputs=AGNinputs, Lagn_params=Lagn_params,
            Z_central_cor=Z_central_cor,
            epsilon_params=epsilon_params,
            extra_params=extra_params,extra_params_names=extra_params_names,
            extra_params_labels=extra_params_labels,
            attmod=attmod, unemod_sfr=unemod_sfr, 
            unemod_agn=unemod_agn, photmod_sfr=photmod_sfr,
            photmod_agn=photmod_agn,
            verbose=True)

    if make_plots:
        # Make test plots
        make_testplots(infile,zz,verbose=True)
