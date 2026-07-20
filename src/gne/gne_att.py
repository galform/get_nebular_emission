#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:36:13 2023

@Author: expox7, vgp
"""

import h5py
import numpy as np
import gne.gne_io as io
import gne.gne_const as c
import sys
import re
from gne.gne_cosmology import emission_line_flux

def get_f_saito20(z):
    '''
    Calculate the dust attenuation fraction parameter
    from Eq. 5 in Saito+2020, his parameter is defined as
    f = E_stellar(B-V)/E_gas(B-V)
    and is a function of redshift (z)

    Parameters
    ----------
    z : float
      Redshift
     
    Returns
    -------
    f : float
    '''
    if (z <= 2.8):
        f = 0.44 + 0.2*z
    else:
        f = 1
    return f   


def get_AlAv_cardelli89(waveAA,Rv=c.Rv):
    '''
    Given the wavelength, returns Al/Av following
    Cardelli+1989 (doi:10.1086/167900)

    Parameters
    ----------
    waveAA : array of floats
     Wavelength (A)
     
    Returns
    -------
    Al_Av : array floats
    '''
    
    wl=waveAA/10000. #microm
    x=1./wl
    
    if (x < 0.3) or (x > 10):
        print('STOP (gne_att.cardelli): ',
              'Wavelength out of range.')
        sys.exit()
        return
    elif (x <= 1.1): #IR
        ax = 0.574*(x**1.61) 
        bx = -0.527*(x**1.61)
    elif (x <= 3.3): #Optical/NIR
        y = x-1.82
        ax = (1.+0.17699*y - 0.50447*(y**2) - 0.02427*(y**3) +
              0.72085*(y**4) + 0.01979*(y**5) -
              0.77530*(y**6) + 0.32999*(y**7)) 
        bx = (1.41338*y + 2.28305*(y**2) + 1.07233*(y**3) -
              5.38434*(y**4) - 0.62251*(y**5) +
              5.30260*(y**6) - 2.09002*(y**7))
    elif (x <= 8): #UV
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


def get_A_lambda(tau,costheta=c.costheta,albedo=c.albedo):
    """
    Calculate the attenuation coefficient, A_lambda.
    
    Arguments can be floats or arrays of floats.
    """
    # Vectorise
    tau = np.asarray(tau)
    costheta = np.asarray(costheta)
    albedo = np.asarray(albedo)

    # Calculate the attenuation coefficient
    al = tau*np.sqrt(1.0 - albedo)    
    x = al/costheta
    A_lambda = -2.5 * np.log10((1.0 - np.exp(-x)) / x)
    
    return A_lambda


def logNH_delucia07(lmgas,hgas):
    '''
    Calculate the mean hydrogen column density
    following De Lucia and Blaziot 2007 (0606519)
    '''
    logNH = np.zeros(lmgas.shape)

    ind = np.where((hgas > 0) & (lmgas >0))
    if np.shape(ind)[1]>0:
        rcm = c.Mpc_to_cm*hgas[ind]/c.re2hr
        logm_at = lmgas[ind] + np.log10(c.Msun) - np.log10(1.4*c.mp)
        logNH[ind] = logm_at - np.log10(np.pi) - 2*np.log10(1.68*rcm)
    return logNH


def att_favole20(wavelengths,lzgas,lmgas,hgas,Rv=c.Rv,costheta=c.costheta,
                 albedo=c.albedo):
    '''
    Calculate attenuation coefficients following Favole+2020 (1908.05626)

    Parameters
    ----------
    wavelenths : array of floats
        Wavelengths (AA) at which calculate the attenuation coefficients
    lzgas : array of floats
        log10(Z_cold_gas)
    lmgas : array of floats
        log10(M_cold_gas)
    Rv : float
        Slope of the extinction curve in the optical
    costheta : float
        Cosine of the typical dust scattering angle
    albedo : float
        Dust albedo (fraction between 0 and 1)
    Return
    ------
    coeff : array of floats
    '''
    ngal = len(lzgas)
    coeff = np.full((len(wavelengths),ngal),c.notnum)

    # Calculate the mean hydrogen column density
    logNH = logNH_delucia07(lmgas,hgas)

    # For each galaxy get the attenuation for each wavelength
    for ii,wl in enumerate(wavelengths):
        # Use the Cardelli+89 model to obtain A_l/A_V
        alav = get_AlAv_cardelli89(wl,Rv=Rv)

        #Calculate the galaxy optical depth
        if wl > 2000:
            # Following Guiderdoni & Roca-Volmerange 1987 
            s = 1.6
        else:
            s = 1.35

        ind = np.where((lzgas > c.notnum) & (logNH >0))
        if np.shape(ind)[1]>0:
            logz_norm = lzgas[ind] - np.log10(c.zsun)
            logNH_norm = logNH[ind] - np.log10(2.1) - 21.
            tau = alav*10.**(s*logz_norm + logNH_norm)
            
            # For each galaxy calculate the attenuation coefficient
            coeff[ii,ind] =  get_A_lambda(tau,
                                          costheta=costheta,albedo=albedo)
    return coeff                                          


def find_line_index(line, line_names):
    """
    Find the index of a line in line_names, allowing for ±1 in the numeric part.
    
    Parameters:
    -----------
    line : str
        Line name to find (e.g., 'NII6583', 'SII6716')
    line_names : array-like
        List of line names to search in
    
    Returns:
    --------
    int or None
        Index of matching line, or None if not found
    """   
    # Try exact match
    for i, name in enumerate(line_names):
        if line == name:
            return i
    
    # Extract names with letters and numbers
    match = re.match(r'([A-Za-z]+)(\d+)', line)
    if not match: 
        return None    
    prefix = match.group(1)
    number = int(match.group(2))

    # Search for match with ±1 tolerance
    for i, name in enumerate(line_names):
        name_match = re.match(r'([A-Za-z]+)(\d+)', name)
        if name_match:
            name_prefix = name_match.group(1)
            name_number = int(name_match.group(2))
            
            if prefix == name_prefix and abs(number - name_number) <= 1:
                return i   
    return None


def gne_att(infile, outpath=None,out_ending=None,
            attmod='cardelli89',line_att=None,
            att_config=None,verbose=True):
    '''
    Get the neccessary information to calculate
    dust-attenuated luminosities
    
    Parameters
    ----------
    infile : string
        Input file
    outpath : string
        Path to output, default is output/ 
    out_ending : string
        Name root for output file
    attmod : str
        Attenuation model ('cardelli89' or 'ratios')
    att_config : dict
       Model-specific parameters:
       - 'cardelli89': {'s': 1.6, 'costheta': 0.3}
       - 'ratios': {'ratios': array, 'rlines': list}
    verbose : boolean
       If True print out messages.
    '''
    # Get emission lines
    lfile= io.get_outnom(infile,dirf=outpath,nomf=out_ending,
                         verbose=verbose)
    f = h5py.File(lfile, 'r') 
    header = f['header']
    zz = header.attrs['redshift']
    photmod_sfr = header.attrs['photmod_sfr']
    lnames = c.line_names[photmod_sfr]
    nlines = len(lnames)
    group = 'sfr_data'
    ncomp = np.shape(f[group+'/'+lnames[0]+'_sfr'][:])[0]
    ngal = np.shape(f[group+'/'+lnames[0]+'_sfr'][:])[1]
    neblines = np.zeros((nlines, ncomp, ngal))
    for i, line in enumerate(lnames):
        neblines[i, :, :] = f[group+'/'+line+'_sfr'][:]
    if attmod == 'ratios':
        rr = 'ratio_'
        fkeys = f['data'].keys()
        att_rlines = [name.replace(rr, '')
                      for name in fkeys if name.startswith(rr)]
        if not att_rlines:
            if verbose: print('WARNING: No attenuation calculation, as no line ratios, '+rr+'*, found in '+infile)
            return
        # Read the ratios into a dictionary
        att_ratios= {rr+name: f['data/'+rr+name][:]
                     for name in att_rlines}
    else: 
        mgasr_type_attr = header.attrs.get('mgasr_type')
        if mgasr_type_attr is not None:
            mgasr_type = io.decode_string_list(mgasr_type_attr)
        else:
            if verbose: print('WARNING: No attenuation calculation, as no gas information found in '+infile)
            return

        lzgas = f[group+'/lz_sfr'][:]  # log10(Z_cold_gas)
        lm_gas = f['data/lm_gas'][:] # log10(M/Msun)
        h_gas = f['data/h_gas'][:]  #Scalelength(Mpc)

    if 'agn_data' not in f.keys():
        AGN = False
    else:
        AGN = True
        group_agn = 'agn_data'
        photmod_agn = header.attrs['photmod_NLR']
        lnames_agn = c.line_names[photmod_agn]
        nlines = len(lnames_agn)
        ngal = len(f[group_agn+'/'+lnames_agn[0] +'_agn'])
        neblines_agn = np.zeros((nlines, ngal))
        for i, line in enumerate(lnames_agn):
            neblines_agn[i, :] = f[group_agn+'/'+line+'_agn'][:]
        if attmod != 'ratios':
            lzgas_agn = f[group_agn+'/lz_agn'][:] # log10(Z_cold_gas)
    f.close()

    # Add attenuation information to the line file
    nattrs = io.add2header(lfile,['attmod'],[attmod],verbose=verbose)
    
    # Get wavelengths and prep attenuation coefficients
    wavelengths = c.line_wavelength[photmod_sfr]
    coeff = np.full(neblines.shape,c.notnum)
    if AGN:
        wavelengths_agn = c.line_wavelength[photmod_agn]
        coeff_agn = np.full(neblines_agn.shape,c.notnum)
        
    # Get the information needed by the dust attenuation model
    if attmod not in c.attmods:
        if verbose:
            print('STOP (gne_att): Unrecognised dust model.')
            print('                Possible ones= {}'.format(c.attmods))
        sys.exit()
    elif attmod == 'ratios':
        # Add components if attenuation ratios given for the whole
        first_value = next(iter(att_ratios.values()))
        if first_value.ndim == 1:
            neblines = np.sum(neblines, axis=1)
            coeff = np.full(neblines.shape,c.notnum)
        elif first_value.ndim != ncomp:
            if verbose: print('WARNING: No attenuation calculation, different number of components in '+infile)
            return           

        # Assign ratios to coeff matrix, adequately
        valid_lnames = []
        for line in att_rlines:
            ii = find_line_index(line, lnames)
            if ii is not None:
                coeff[ii,:] = att_ratios[rr+line]
                valid_lnames.append(lnames[ii])

        # Keep only valid lines
        valid_indices = np.array([i for i, line in enumerate(lnames)
                                  if line in valid_lnames])
        lnames = lnames[valid_indices]
        neblines = neblines[valid_indices]  # Modifies 1st axis with lines
        coeff = coeff[valid_indices]        # Modifies 1st axis with lines
                
        if AGN:
            valid_lnames = []
            for line in att_rlines:
                ii = find_line_index(line, lnames_agn)
                if ii is not None:
                    coeff_agn[ii,:] = att_ratios[rr+line]
                    valid_lnames.append(lnames_agn[ii])

            # Keep only valid lines
            valid_indices = np.array([i for i, line in enumerate(lnames_agn)
                                      if line in valid_lnames])
            lnames_agn = lnames_agn[valid_indices]
            neblines_agn = neblines_agn[valid_indices] # Modify 1st axis with lines
            coeff_agn = coeff_agn[valid_indices]       # Modify 1st axis with lines
    else:
        Rv = io.get_param(att_config, 'Rv', c.Rv)
        costheta = io.get_param(att_config, 'costheta', c.costheta) 
        albedo = io.get_param(att_config, 'albedo', c.albedo)
        # Add parameters to header
        config_names  = ['Rv', 'costheta', 'albedo']
        config_values = [Rv, costheta, albedo]
        io.add2header(lfile,config_names,config_values,verbose=False)

        if attmod == 'favole20':
            icomp = mgasr_type.index('disc')
            coeffd = att_favole20(wavelengths,lzgas[:,icomp],
                                  lm_gas[:,icomp],h_gas[:,icomp],
                                  Rv=Rv,costheta=costheta,
                                  albedo=albedo)
            coeff[:, :, :] = coeffd[:, np.newaxis, :]
            
        elif attmod == 'favole20_percomponent':
            for icomp in range(ncomp):
                coeff[:,icomp,:] = att_favole20(wavelengths,lzgas[:,icomp],
                                                lm_gas[:,icomp],
                                                h_gas[:,icomp],
                                                Rv=Rv,costheta=costheta,
                                                albedo=albedo)

        if AGN: ####here This to be checked: does it make sense to use the disc?
            icomp = mgasr_type.index('disc')
            coeff_agn = att_favole20(wavelengths_agn,lzgas[:,icomp],
                                  lm_gas[:,icomp],h_gas[:,icomp],
                                  Rv=Rv,costheta=costheta,
                                  albedo=albedo)        
            
    # Lines with valid attenuation coefficients
    ind = (coeff > c.notnum)
    if AGN:
        ind_agn = (coeff_agn > c.notnum)

    # Extra gas attenuation on top of the stellar one
    if line_att and attmod != 'ratios':
        f = get_f_saito20(zz)
        coeff[ind] = coeff[ind]/f
        if AGN:
            coeff_agn[ind_agn] = coeff_agn[ind_agn]/f
            
        #Add information in header
        nattrs = io.add2header(lfile,['att_f_saito20'],[f],verbose=verbose)

    # Calculate dust-attenuated emission lines
    neblines_att = neblines
    if attmod == 'ratios':
        neblines_att[ind] = neblines[ind]*coeff[ind]
    else:
        neblines_att[ind] = neblines[ind]*(10**(-0.4*coeff[ind]))
    lnames_att = [line + '_sfr_att' for line in lnames]
    labels = np.full(np.shape(lnames_att),'erg s^-1')
    io.write_data(lfile,group=group,params=neblines_att,
                  params_names=lnames_att,params_labels=labels)

    if AGN:
        neblines_agn_att = neblines_agn
        if attmod == 'ratios':
            neblines_agn_att[ind_agn] = neblines_agn[ind_agn]*coeff_agn[ind_agn]
        else:
            neblines_agn_att[ind_agn] = neblines_agn[ind_agn]*(
                10**(-0.4*coeff_agn[ind_agn]))
        lnames_att = [line + '_agn_att' for line in lnames_agn]
        labels = np.full(np.shape(lnames_att),'erg s^-1')

        io.write_data(lfile,group=group_agn,params=neblines_agn_att,
                      params_names=lnames_att,params_labels=labels)

    print('SUCCESS (gne_attenuation)')
    return
