"""
This is a modular python code to provide model spectral emission 
lines, both from star-forming regions and narrow-line regions of AGNs.
The input of the code are global galactic properties. 

The intrinsic luminosities can be passed through an attenuation model 
to also get the predicted attenuated luminosities.

@authors: viogp
"""

import src.gne_const as const
from src.gne import gne
from src.gne_plots import make_testplots
import h5py

### RUN the code with the given parameters and/or make plots
testing = True # use only the first 50 elements
run_code = True
make_plots = True

# Calculate emission from AGNs: AGN = True
AGN = True

###############################################################
### INPUT FILES
# Input files are expected to have, AT LEAST:
# Stellar mass (M*) of the galaxy (or disc or buldge).
# Star formation rate (SFR) OR magnitude of Lyman Continuum photons (m_LC).
# Mean metallicity of the cold gas (Z).
infiles = ['data/example_data/GP20_62.5kpc_z0_example.hdf5']

### INPUT FORMAT ('txt' for text files; 'hdf5' for HDF5 files)
inputformat = 'hdf5'

### OUTPUT PATH (Default: output/)
outpath = None  

####################################################
############  Emission from SF regions #############
####################################################

# All available models can be seen in gne_const module.
# NEBULAR model connecting global properties to nebular parameters
unemod_sfr='kashino20'
# PHOTOIONIZATION model for SF regions to get line luminosities
photmod_sfr='gutkin16'


### INPUT PARAMETERS
# m_sfr_z has the location in the input files of the three mandatory parameters:
# M*(units), SFR or m_LC and Zgas. 
# m_sfr_z is a list of lists with either the column number
# for each parameters or the name of the HDF5 variable.
m_sfr_z = [['data/mstar_disk','data/SFR_disk','data/Zgas_disk'],
           ['data/mstar_bulge','data/SFR_bulge','data/Zgas_bulge']]

# mtot2mdisk is True if the stellar mass of discs is calculated 
# from the total and buldge values (False by default)
# mtot2mdisk = True; cols = [[M,SFR,Z],[M_bulge,SFR_bulge,Z_bulge]]
# mtot2mdisk = False; cols = [[M_disk,SFR_disk,Z_disk],[M_bulge,SFR_bulge,Z_bulge]]        
mtot2mdisk = False

# LC2sfr is True when Lyman Continuum photons are given instead of the SFR
# LC2sfr = True; cols = [[M,m_LC,Z]]
# LC2sfr = False; cols = [[M,SFR,Z]] (Default option)      
LC2sfr = False

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

# All available models can be seen in gne_const module.
# NEBULAR model connecting global properties to nebular parameters
unemod_agn = 'panuzzo03'
# PHOTOIONIZATION model for AGN regions to get line luminosities
photmod_agn = 'feltre16'


# mg_r50 has the location of the following parameters:
# Cold gas mass (Mg).
# Baryonic half-mass radius (R50).
# mg_r50 is a list of lists with either the column number
# for each parameters or the name of the HDF5 variable.
# For bulge and disk:
mg_r50 = ['data/mgas_disk','data/rhm_disk',
          'data/mgas_bulge','data/rhm_bulge']
    
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
Lagn_params=['data/lagn','data/mstar_bulge']

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
    # If local relations are to be used: infiles_z0 = [None]
infiles_z0 = [None]


####################################################
##########       Dust attenuation      #############
####################################################

# Continuum and line attenuation calculation. If this option is selected 
    # the output file will have intrinsic AND attenuated values of
    # luminosity for the emission lines.
    
# att=True to calculate the dust attenuation; False, otherwise
att = True
    
# To use Cardelli's law (following Favole et. al. 2020):
    # attmod = 'cardelli89'.
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
att_params=['data/rhm_disk','data/mgas_disk','data/Zgas_disk'] 

####################################################
##########      Other calculations     #############
####################################################

# Flux calculation.
flux=True

# Include other parameters in the output files
extra_params_names = ['mh','magK','magR','type','MBH']
extra_params_labels = extra_params_names
extra_params = ['data/mh','data/magK','data/magR','data/type','data/MBH']

# Make the calculation on a subsample based on selection cuts
# Paramter to impose cuts
cutcols = ['data/mh']
# List of minimum values. None for no inferior limit.
mincuts = [20*9.35e8]
# List of maximum values. None for no superior limit.
maxcuts = [None]

##################################################################
#############    Run the code and or make plots   ################
##################################################################

for ii, infile in enumerate(infiles):
    infile_z0 = infiles_z0[ii]
    
    # Get the redshift, cosmology and volume of the model galaxies
    f = h5py.File(infile)
    header = f['header']
    zz = header.attrs['redshift']
    snapshot = header.attrs['snapnum']
    vol = header.attrs['bside_Mpch']**3
    h0 = header.attrs['h0']
    omega0 = header.attrs['omega0']
    omegab = header.attrs['omegab']
    lambda0 = header.attrs['lambda0']    

    if run_code:
        # Run the code
        gne(infile, zz,h0,omega0,omegab,lambda0,vol,
            inputformat=inputformat, outpath=outpath,
            unemod_sfr=unemod_sfr, photmod_sfr=photmod_sfr,
            m_sfr_z=m_sfr_z,mtot2mdisk=mtot2mdisk, LC2sfr=LC2sfr,
            inoh=inoh,IMF = IMF,
            AGN=AGN,
            unemod_agn=unemod_agn, photmod_agn=photmod_agn,
            mg_r50=mg_r50,
            AGNinputs=AGNinputs, Lagn_params=Lagn_params,
            Z_central_cor=Z_central_cor,
            infile_z0=infile_z0, 
            att=att, attmod=attmod, att_params=att_params,
            flux=flux,
            extra_params=extra_params,extra_params_names=extra_params_names,
            extra_params_labels=extra_params_labels,
            cutcols=cutcols, mincuts=mincuts, maxcuts=maxcuts,
            testing=testing,verbose=True)

    if make_plots:
        # Make test plots
        make_testplots(infile,zz,verbose=True)
