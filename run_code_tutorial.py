#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 08:42:21 2023

@author: expox7
"""

import get_nebular_emission.eml as eml
import get_nebular_emission.eml_const as const
import glob

####################################################
####################################################
####################################################

# This is a modular python code to provide model spectral emission 
# lines, both from star-forming regions and narrow-line regions of AGNs.
# The input of the code are global galactic properties. 

# Emission from star-forming region calculation uses the stellar mass, 
# star formation rate and cold gas metallicity to generate a range of 
# spectral emission lines. 
# Emission from AGNs requires the bolometric 
# luminosity of the AGN, the total mass of cold gas and the effective 
# radius of the galaxy components. 

# These parameters are passed through a nebular model to get the parameters 
# defining the ionizing sources, and then those nebular properties are passed 
# through a photoionization grid to get the intrinsic 
# luminosities of the different emission lines available in the 
# used model by interpolation in such grid. 

# The intrinsic luminosities can be passed through an attenuation model 
# to also get the predicted attenuated luminosities.

####################################################
####################################################
####################################################


####################################################
################       BASICS       ################
####################################################


### INPUT FILES
# It has to be a list with the path to the input files.
# Each file correspond to a sub-volume of the simulation you want to pass
    # through get_nebular_emission.
# If all the galaxy data is in one single input file, it has to be in a list:
    # [input file path]
# The input files are expected to have, AT LEAST, this information 
    # for each galaxy or galaxy component (disk and bulge): 
        # Stellar mass (M*).
        # Star formation rate (SFR) OR magnitude of Lyman Continuum photons (m_LC).
        # Mean metallicity of the cold gas (Z).
infile = glob.glob(r"GP20/GP20_0.0_*.txt")


# OUTPUT FILE
# The path to the output file. It will always create this file and fill it
    # with the results, overwriting it if it already exist.
outfile = r"output_data/emlines_GP20_z0.0_Kashino_test.hdf5"

### REDSHIFT
# Put here the redshift of your galaxy sample.
redshift = 0.15


### INPUT FORMAT
# The code can work with text files (.txt, .dat, .csv...) or with HDF5 files.
# If your input files are text files: inputformat = 'txt'
# If your input files are HDF5 files: inputformat = 'hdf5'
inputformat = 'txt'


### PARAMETERS FOR THE CALCULATION OF THE INTRINSIC EMISSION FROM STAR-FORMING REGIONS
# cols has the location in the input files of the three mandatory parameters:
    # M*, SFR/m_LC and Z. 
    # It is a list of lists with the location of the parameters. Each list
    # correspond to a different component: 
    # cols = [[M_disk,SFR_disk,Z_disk],[M_bulge,SFR_bulge,Z_bulge]]
    # In the case of a single component:
    # cols = [[M,SFR,Z]] 
# Example:
    # For one component galaxies and a text file as input file, if the total
    # stellar mass is in the first column,
    # the SFR is in the second,
    # and the mean metallicity of the cold gas is in the fourth, then:
    # cols = [[0,1,3]]
# Example 2:
    # For two component galaxies and a HDF5 file as input file:
    # cols = [['Ms_disk','SFR_disk','Z_disk'],['Ms_bulge','SFR_bulge','Z_bulge']]
    # Supposing that, for example, 'Ms_disk' is the name of the HDF5 file's 
    # dataset where the values of the stellar mass of the disk are stored.
cols = [[0,2,4],[1,3,5]]
# cols = [[0,2,4]]

# LC2sfr is True when instead of the SFR you input the magnitude of  
    # Lyman Continuum photons.
# So:
    # First option:
        # LC2sfr = True
        # cols = [[M,m_LC,Z]]
    # Second option:
        # LC2sfr = False
        # cols = [[M,SFR,Z]]      
LC2sfr = False

# mtot2mdisk is True when, for galaxy data with two components (disk and bulge),
    # you want to input total values and values for the bulge, instead of 
    # values for the disk and values for the bulge.
# When True, the code internally calculates the values for the disc 
    # from the total and bulb values and then performs the calculations for 
    # the disc and bulb separately.
# So:
    # First option:
        # mtot2mdisk = True
        # cols = [[M,SFR,Z],[M_bulge,SFR_bulge,Z_bulge]]
    # Second option:
        # mtot2mdisk = False
        # cols = [[M_disk,SFR_disk,Z_disk],[M_bulge,SFR_bulge,Z_bulge]]        
mtot2mdisk = False

### NEBULAR AND PHOTOIONIZATION MODELS
# These variables indicate which model are going to be used by the code.
# unemod is the model that connects global galactic properties with 
    # nebular parameters.
# photmod is the photoionization model that calculates emission line
    # luminosities from the nebular parameters.
# sfr correspond to calculation of emission lines from star-forming galaxies
    # and agn to calculation of emission lines from the narrow-line region 
    # of AGNs.
# All available models can be seen in eml_const module.

unemod_sfr='kashino20'
unemod_agn='panuzzo03'

photmod_sfr='gutkin16'
photmod_agn='feltre16'

####################################################
####################################################
####################################################



####################################################
###############       OPTIONAL       ###############
####################################################


####################################################
### HIGH REDSHIFT CORRECTION ###

# Empirical relationships to connect global galaxy properties and nebular
    # properties are often derived from local galaxies. get_emission_lines has
    # a way of partially correcting their results for high-redshift galaxies
    # taking the evolution of the mean gas density of galaxies into 
    # consideration.
    # See section 3.2.3 of ??? for a detailed explanation.

# FILES WITH SIMULATIONS AT REDSHIFT 0.
# Use this if you are calculating emission lines for z > 0 simulations, you
    # are using an empirical relationship for the connection between global
    # galaxy properties and nebular properties taken with local galaxies, 
    # and you want to correct the results following ??? methodology.
# There have to be the same number of files than input files.
# If you don't want to use it, remove it from the eml.eml call or define it as
    # infile_z0 = [None]
    
infile_z0 = [None] #glob.glob(r"GP20/GP20_0.0_*.txt")
####################################################

####################################################
### INITIAL MASS FUNCTION CORRECTION ###

# You can make a correction to the star formation rate and 
    # stellar mass of the semi-analytic models from the IMF assumed in the 
    # simulations to another IMF of your election, following the correlation
    # coefficients from Lagos et. al. 2011.
    # Possible IMFs: 
        # Kennicut, Salpeter, Kroupa, Chabrier, Baldry&Glazebrook, Top-heavy
        
# Note that empirical relationships for nebular models assume a certain IMF!

# IMF_i are the assumed IMFs for each galaxy component in the input data.
# IMF_f are the IMFs to which you want to transform each component.
# Example:
    # For one component:
        # IMF_i = [['Kennicut']]
        # IMF_f = [['Kroupa']]
        # This change the mean stellar mass and SFR of the input data from
        # an assumed Kennicut IMF to a Kroupa IMF.
    # For two components:
        # IMF_i = ['Kennicut','Kennicut']
        # IMF_f = ['Kroupa','Top-heavy'] 
        # This does the same for disk and bulge separately. Note that for the
        # emission line calculations this supposes a Kroupa IMF for the disk
        # and a Top-heavy IMF for the bulge.
        
IMF_i = ['Kennicut','Kennicut']
IMF_f = ['Kroupa','Top-heavy']
####################################################

####################################################
### CALCULATION OF THE INTRINSIC EMISSION FROM THE NARROW-LINE REGION OF AGNS ###

# Along with the calculation of the emission from star-forming regions,
# the code can calculate the emission from AGNs, given the necessary
# parameters.

# Calculate emission from AGNs: AGN = True
# Don't calculate emission from AGNs: AGN = False
AGN=True

# AGN emission can be calculated from different parameters in several levels of
    # assumptions and approximations. 
    
# First, some paramers are necesary to 
    # calculate the volume filling-factor of the NLR:
        # Cold gas mass (Mg).
        # Baryonic half-mass radius (Rhm).
    # For one single component:
        # epsilon_params = [Mg, Rhm]
    # For two components, disk and bulge:
        # epsilon_params = [Mg_disk, Rhm_disk, Mg_bulge, Rhm_bulge]
epsilon_params=[6,11,19,12]
    
# Second, the AGNs bolometric luminosity is needed. This value can be in the
    # input file or can be estimated from other paremeters. To indicate
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
    # Zagn from the eml_une module.
# If Z_central_cor is False, it is assumed that the metallicity is 
    # approximatelly uniform in the galaxy.
Z_central_cor=True
####################################################

####################################################
### ATTENUATION ###

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

# Flux calculation.
flux=True

####################################################

extra_params_names = ['Type','Mbh','Mhalo','Ms_bulge','m_K','m_R','m_R_SDSS','m_I',
                      'Mdot_stb','Mdot_hh','Mhot','Lagn']
extra_params_labels = ['Type of halo (central = 0)',
                       r'Black hole mass ($M_\odot \ h^{-1}$)',
                       r'Bulge mass ($M_\odot \ h^{-1}$)',
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

####################################################
### SELECTION CRITERIA ###

# If you don't want to calculate emission line luminosities for all the 
    # galaxies included in the input files, you5 can make a cut based on
    # values range for any parameter inside the input files.
    
# cutcols is a list with the location of the parameters that will be 
    # considered for the cut.
# mincuts is a list with the minimum value the parameters of cutcols can
    # have. None implies no inferior limit.
# maxcuts is a list with the maximum value the parameters of cutcols can
    # have. None implies no superior limit.

cutcols = [7]
mincuts = [20*9.35e8]
maxcuts = [None]
# In the example, location 7 correspond to the column with the total mass
    # of the galaxy. The dark matter particles of the simulations of 
    # the example has a mass of 9.35 · 10^8 Ms, so we are only taking 
    # galaxies with more than 20 DM particles in the simulations.
####################################################

eml.eml(infile, outfile, m_sfr_z=cols, infile_z0=infile_z0, 
            att=att,
            att_params=att_params, att_ratio_lines=att_ratio_lines,
            flux=flux,
            h0=const.h,
            IMF_i=IMF_i, IMF_f=IMF_f, inputformat=inputformat,
            mtot2mdisk=mtot2mdisk, LC2sfr=LC2sfr,
            redshift=redshift,
            AGN=AGN,flag=0,
            AGNinputs=AGNinputs, Lagn_params=Lagn_params,
            Z_central_cor=Z_central_cor,
            epsilon_params=epsilon_params,
            extra_params=extra_params,extra_params_names=extra_params_names,
            extra_params_labels=extra_params_labels,
            attmod=attmod, unemod_sfr=unemod_sfr, 
            unemod_agn=unemod_agn, photmod_sfr=photmod_sfr,
            photmod_agn=photmod_agn,
            verbose=True)