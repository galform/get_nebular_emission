#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:36:13 2023

@author: expox7
"""

import h5py
import numpy as np
from get_nebular_emission.eml_io import get_nheader, homedir, locate_interval
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_io import check_file
import sys
import warnings
from cosmology import emission_line_flux

#------------------------------------------------------------------------------------
#   Cardelli et al. 1989 extinction laws in FIR and IR/OPT:
#------------------------------------------------------------------------------------
def cardelli(waveA):
    '''
    Given the wavelength, returns Al/Av following Cardelli 1989.

    Parameters
    ----------
    waveA : floats
     Wavelength (A)
     
    Returns
    -------
    Al_Av : floats
    '''
    
    Rv=3.1 
    wave=waveA/10000.
    x=1./wave
    
    if (x < 0.3) or (x > 10):
        print('STOP (eml_photio.cardelli): ',
              'Wavelength out of range.')
        sys.exit()
        return
    elif (x < 1.1): #IR
        ax = 0.574*(x**1.61) 
        bx = -0.527*(x**1.61)
    elif (x < 3.3): #Optical/NIR
        y = x-1.82
        ax = (1.+0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 
        0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7) 
        bx = (1.41338*y+2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 -
        0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7)
    elif (x < 8): #UV
        if (x < 5.9):
            Fa = 0
            Fb = 0
        else: 
            Fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
            Fb = 0.2130*(x-5.9)**2 + 0.1207*(x-5.9)**3
        ax = 1.752-0.316*x - 0.104/((x-4.67)**2 + 0.341) + Fa
        bx = -3.090+1.825*x + 1.206/((x-4.62)**2 + 0.263) + Fb
    else:
        ax = -1.073 - 0.628*(x-8) + 0.137*(x-8)**2 - 0.070*(x-8)**3
        bx = 13.670 + 4.257*(x-8) - 0.420*(x-8)**2 + 0.374*(x-8)**3
    
    Al_Av = ax+bx/Rv
    return Al_Av

def coef_att_cardelli(wavelength, Mcold_disc, rhalf_mass_disc, Z_disc, costheta=0.3, albedo=0.56):
    '''
    Given the wavelength, the cold gas mass, the half-mass radius and the global metallicity of the disk,
    along with the assumed albedo and scattering angle, it gives the attenuation coefficient.

    Parameters
    ----------
    wavelength : floats
     Wavelength (A).
    Mcold_disc : floats
     Cold gas mass (Msun).
    rhalf_mass_disc : floats
     Half-mass radius of the disk (Msun).
    Z_disc : floats
     Disc's global metallicity.
    costheta : float
     Cosine of the assumed cattering angle.
    albedo : float
     Assumed albedo.
     
    Returns
    -------
    coef_att : floats
    '''
    
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    
    a_disc = 1.68
    if wavelength > 2000: #Guiderdoni & Roca-Volmerange 1987 (see De Lucia & Blaizot 2007)
        s = 1.6
    else:
        s = 1.35
    Al_Av = cardelli(wavelength)
    sectheta = 1./costheta

    mean_col_dens_disc_log = (np.log10(Mcold_disc*const.h) + np.log10(const.Msun_to_kg) - 
    np.log10(1.4*const.mp*np.pi)-2.*np.log10(a_disc*rhalf_mass_disc*const.h*const.Mpc_to_cm))
    
    tau_disc = np.log10(Al_Av) + np.log10((Z_disc/const.zsun)**s) + mean_col_dens_disc_log - np.log10(2.1e21)
    tau_disc = 10.**tau_disc
    
    al_disc = (np.sqrt(1.-albedo))*tau_disc
    
    A_lambda = -2.5*np.log10((1.-np.exp(-al_disc*sectheta))/(al_disc*sectheta))
    
    coef_att = 10.**(-0.4*A_lambda)
    
    return coef_att

def coef_att_ratios(infile,cols_notatt,cols_att,cols_photmod,inputformat='HDF5',photmod='gutkin16',verbose=True):
    '''
    It reads luminosities of lines with and without attenuation
    from line emission data and it returns the attenuation coefficients.

    Parameters
    ----------
    infile : string
     Name of the input file. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    cols_notatt : list
     Attenuated flux lines calculated by the semi-analytic model of the input data. 
     Used to calculate attenuation coefficients for the "ratio" attenuation model.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    cols_att : list
     Not attenuated flux lines calculated by the semi-analytic model of the input data. 
     Used to calculate attenuation coefficients for the "ratio" attenuation model.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    cols_photmod : list
     Index in the list of lines of the photoionization model of the lines for which 
     attenuation is going to be calculated in the "ratio" attenuation model.
    inputformat : string
     Format of the input file.
    photmod : string
     Photoionisation model to be used for look up tables.
    verbose : boolean
     If True print out messages.

    Returns
    -------
    coef_att : floats
    '''
    
    check_file(infile, verbose=verbose)
    
    ncomp = len(cols_notatt)
        
    lines = const.lines_model[photmod]
    numlines = len(lines)
    
    cols_att = np.array(cols_att)
    cols_notatt = np.array(cols_notatt)
    cols_photmod = np.array(cols_photmod)
    
    if inputformat=='HDF5':
        with h5py.File(infile, 'r') as f:
            hf = f['data']
            
            coef_att = np.empty((ncomp,numlines,len(hf[cols_notatt[0,0]])))
            coef_att.fill(const.notnum)
            
            for i in range(len(cols_photmod)):
                for comp in range(ncomp):
                    ind = np.where(hf[cols_notatt[comp,i]][:]!=0)
                    ind2 = np.where(hf[cols_notatt[comp,i]][:]==0)
                    coef_att[comp,cols_photmod[i]][ind] = hf[cols_att[comp,i]][ind]/hf[cols_notatt[comp,i]][ind]
                    coef_att[comp,cols_photmod[i]][ind2] = 1
    elif inputformat=='textfile':
        ih = get_nheader(infile)        
        X = np.loadtxt(infile,skiprows=ih).T
        
        coef_att = np.empty((ncomp,numlines,len(X[0])))
        coef_att.fill(const.notnum)
        
        for i in range(len(cols_photmod)):
            for comp in range(ncomp):
                if ncomp!=1:
                    ind = np.where(X[cols_notatt[comp,i]]!=0)
                    ind2 = np.where(X[cols_notatt[comp,i]]==0)
                    coef_att[comp,cols_photmod[i]][ind] = X[cols_att[comp,i]][ind]/X[cols_notatt[comp,i]][ind]
                    coef_att[comp,cols_photmod[i]][ind2] = 1
                else:
                    ind = np.where(X[cols_notatt[comp,i]]!=0)
                    ind2 = np.where(X[cols_notatt[comp,i]]==0)
                    if comp==0:
                        coef_att[comp,cols_photmod[i]][ind] = (X[cols_att[comp,i]][ind]-X[cols_att[1,i]][ind])/(X[cols_notatt[comp,i]][ind]-X[cols_notatt[1,i]][ind])
                        coef_att[comp,cols_photmod[i]][ind2] = 1
                    else:
                        coef_att[comp,cols_photmod[i]][ind] = X[cols_att[comp,i]][ind]/X[cols_notatt[comp,i]][ind]
                        coef_att[comp,cols_photmod[i]][ind2] = 1
                    
                    
        del X        
    
    return coef_att

def attenuation(nebline, att_param=None, att_ratio_lines=None,
                redshift=0, attmod='cardelli89',origin='sfr',
                photmod='gutkin16', cut=None, verbose=True):
    '''
    Get the attenuated emission lines from the raw ones.

    Parameters
    ----------
    nebline : floats
     Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr).
    infile : string
     - Name of the input file. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    att_params : list
     Parameters to look for calculating attenuation. See eml_const to know what each model expects.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    redshift : float
     Redshift of the input data.
    inputformat : string
     Format of the input file.
    attmod : string
     Attenuation model.
    cols_notatt : list
     Attenuated flux lines calculated by the semi-analytic model of the input data. 
     Used to calculate attenuation coefficients for the "ratio" attenuation model.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    cols_att : list
     Not attenuated flux lines calculated by the semi-analytic model of the input data. 
     Used to calculate attenuation coefficients for the "ratio" attenuation model.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    cols_photmod : list
     Index in the list of lines of the photoionization model of the lines for which 
     attenuation is going to be calculated in the "ratio" attenuation model.
    photmod : string
     Photoionisation model to be used for look up tables.
    cut : integers
     Indeces of the not cutted galaxies.
    verbose : boolean
     If True print out messages.

    Returns
    -------
    nebline_att, coef_att : floats
      
    Notes
    -------
    It will ignore the lines for which, for any reason, attenuation can not be calculated.
    '''
    
    ncomp = len(nebline)
    coef_att = np.full(nebline.shape,const.notnum)
    
    if att_param[0][0] != None:
        if attmod not in const.attmods:
            if verbose:
                print('STOP (eml_photio.attenuation): Unrecognised model for attenuation.')
                print('                Possible attmod= {}'.format(const.attmods))
            sys.exit()
        elif attmod=='ratios':
            for i, line in enumerate(const.lines_model[photmod]):
                if line in att_ratio_lines:
                    ind = np.where(np.array(att_ratio_lines)==line)[0]
                else:
                    continue
                
                for comp in range(ncomp):
                    if comp==0:
                        coef_att[comp,i] = att_param[ind] * const.line_att_coef_all(redshift)
                    else:
                        if origin=='sfr':
                            coef_att[comp,i] = coef_att[0,i]
                        else:
                            coef_att[comp,i] = 1.
            coef_att[(coef_att!=coef_att)&(coef_att!=const.notnum)] = 1.
        elif attmod=='cardelli89':
            Rhm = att_param[0]
            Mcold_disc = att_param[1]
            Z_disc = att_param[2]
                
            coef_att = np.full(nebline.shape,const.notnum)
            for comp in range(ncomp):
                for i, line in enumerate(const.lines_model[photmod]):
                    if comp==0:
                        coef_att[comp,i] = coef_att_cardelli(const.wavelength_model[photmod][i], 
                                    Mcold_disc=Mcold_disc, rhalf_mass_disc=Rhm, 
                                    Z_disc=Z_disc, costheta=0.3, albedo=0.56) * const.line_att_coef_all(redshift)
                    else:
                        coef_att[comp,i] = coef_att[0,i]
            coef_att[(coef_att!=coef_att)&(coef_att!=const.notnum)] = 1.
    
    nebline_att = np.full(nebline.shape,const.notnum)
    ind = np.where((coef_att!=const.notnum))
    nebline_att[ind] = nebline[ind]*coef_att[ind]
    
    return nebline_att, coef_att