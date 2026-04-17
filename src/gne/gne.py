"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
.. contributions:: Olivia Vidal <ovive.pro@gmail.com>
.. contributions:: Julen Expósito-Márquez <expox7@gmail.com>
"""
import time
import numpy as np
import gne.gne_io as io
from gne.gne_model_UnH import get_UnH_sfr, get_UnH_agn
from gne.gne_Z import correct_Z,get_zgasagn
from gne.gne_m_sfr import get_sfrdata
from gne.gne_Lagn import get_Lagn
import gne.gne_const as c
from gne.gne_stats import components2tot
from gne.gne_photio import get_lines, get_limits
from gne.gne_plots import make_testplots

def gne(infile,redshift,snap,h0,omega0,omegab,lambda0,
        mp,boxside,effvol,
        inputformat='hdf5',outpath=None,out_ending=None,
        units_h0=False,units_Gyr=False,units_L=0,
        model_nH_sfr='kashino19',model_U_sfr='kashino19',
        photmod_sfr='gutkin16',nH_sfr=c.nH_sfr_cm3,
        q0=c.q0_orsi, z0=c.Z0_orsi, gamma=c.gamma_orsi,
        T=c.temp_ionising,xid_sfr=0.3,co_sfr=1,
        m_sfr_z=[None],mtot2mdisk=True,inoh=False,
        IMF=['Kroupa','Kroupa'],imf_cut_sfr=100,
        AGN=False,photmod_agn='feltre16',
        Zgas_NLR=None,Z_correct_grad=False,
        model_U_agn='panuzzo03',
        mgas_r=[None],mgasr_type=[None],r_type=[None],
        model_spec_agn='feltre16',
        alpha_NLR=c.alpha_NLR_feltre16,xid_NLR=c.xid_NLR_feltre16,
        nH_NLR=c.nH_NLR_cm3,T_NLR=c.temp_ionising,r_NLR=c.radius_NLR,
        Lagn_inputs='Lagn', Lagn_params=[None],
        zeq=None,infile_z0=None,
        extra_params=[None], extra_params_names=[None],
        extra_params_labels=[None],
        cutcols=[None], mincuts=[None], maxcuts=[None],
        cutlimits=False, 
        testing=False,verbose=True):
    '''
    Calculate emission lines given the properties of model galaxies

    Parameters
    ----------
    infile : string
        Name of the input file.
    redshift : float
        Redshift of the input data.
    m_sfr_z : list of integers (txt input files) or strings (hdf5 files)
        [[component1_stellar_mass,sfr,Z],...]
    inputformat : string
        Format of the input file: 'hdf5' or 'txt'.
    infile_z0 : strings
     List with the name of the input files with the galaxies at redshift 0. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    h0 : float
      If not None: value of h, H0=100h km/s/Mpc.
    snap: integer
        Simulation snapshot number
    cutcols : list
     Parameters to look for cutting the data.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    mincuts : floats
     Minimum value of the parameter of cutcols in the same index. All the galaxies below won't be considered.
    maxcuts : floats
     Maximum value of the parameter of cutcols in the same index. All the galaxies above won't be considered.
    IMF : array of strings
       Assumed IMF for the input data of each component, [[component1_IMF],[component2_IMF],...]
    q0 : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    z0 : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    gamma : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
        If True, corrects Zgas_NLR using gradients from the literature
    extra_params : list
     Parameters from the input files which will be saved in the output file.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    extra_params_names : strings
     Names of the datasets in the output files for the extra parameters.
    extra_params_labels : strings
     Description labels of the datasets in the output files for the extra parameters.
    model_nH_sfr : string
        Model to go from galaxy properties to Hydrogen (or e) number density.
    model_U_sfr : string
        Model to go from galaxy properties to ionising parameter.
    model_spec_agn : string
        Model for the spectral distribution for AGNs.
    model_U_sfr : string
        Model to go from galaxy properties to AGN ionising parameter.
    photmod_sfr : string
        Photoionisation model to be used for look up tables.
    inoh : boolean
       If true, the input is assumed to be 12+log10(O/H), otherwise Zgas
    cutlimits : boolean
     If True the galaxies with U, ne and Z outside the photoionization model's grid limits won't be considered.
    mtot2mdisk : boolean
     If True transform the total mass into the disk mass. disk mass = total mass - bulge mass.
    xid_NLR : float
     Dust-to-metal ratio for the AGN photoionisation model.
    alpha_NLR : array of floats
        Spectral index assumed for the AGN.
    xid_sfr : float
     Dust-to-metal ratio for the SF photoionisation model.
    co_sfr : float
     C/O ratio for the SF photoionisation model.
    imf_cut_sfr : float
     Solar mass high limit for the IMF for the SF photoionisation model.
    AGN : boolean
       If True calculates emission from the narrow-line region of AGNs.
    Lagn_inputs : string
       Type of inputs for AGN's bolometric luminosity calculations.
    Lagn_params : list of integers (text files) or strings (hdf5 files)
       Parameters to obtain the bolometric luminosity.
    Zgas_NLR : list of integer (text file) or strings (hdf5 file)
        Location of the central metallicity in input files
    Z_correct_gradrection : boolean
        True to modify the metallicity by literature gradients
    model_spec_agn : string
        Model for the spectral distribution for AGNs.
    model_U_agn : string
        Model to go from galaxy properties to AGN ionising parameter.
    photmod_agn : string
        Photoionisation model to be used for look up tables.
    nH_NLR : float
        Value assumed for the electron number density in AGN NLR.
    T_NLR : float
        Value assumed for the AGN NLR temperature.
    r_NLR : float
        Value assumed for the radius of the AGN NLR.
    units_h0: boolean
        True if input units with h
    units_Gyr: boolean
        True if input units with */Gyr
    units_L: integer
        0: input units [L]=erg/s  (default);
        1: input units [L]=1e40 h^-2 erg/s
        2: input units [L]=1e40 erg/s
    testing : boolean
        If True only run over few entries for testing purposes
    verbose : boolean
        If True print out messages

    Notes
    -------
    This code returns an .hdf5 file with the mass, specific star formation rate,
    electron density, metallicity, ionization parameter, and the emission lines.

    '''

    # Generate header in the output file from input
    outfile = io.generate_header(infile,redshift,snap,
                                 h0,omega0,omegab,lambda0,
                                 mp,boxside,effvol,
                                 units_h0,outpath=outpath,
                                 out_ending=out_ending,
                                 verbose=verbose)

    #----------------HII region calculation------------------------
    if verbose: print('SF:')        
    # Add relevant constants to header
    names = ['model_nH_sfr','model_U_sfr','photmod_sfr',
             'nH_sfr_cm3','xid_sfr','co_sfr','imf_cut_sfr',
             'q0_orsi','Z0_orsi','gamma_orsi']
    values = [model_nH_sfr,model_U_sfr,photmod_sfr,
              nH_sfr,xid_sfr,co_sfr,imf_cut_sfr,
              q0,z0,gamma]
    nattrs = io.add2header(outfile,names,values,verbose=verbose)
    
    # Number of components
    ncomp = io.get_ncomponents(m_sfr_z)
    
    # Time variables
    start_total_time = time.perf_counter()
    start_time = time.perf_counter()

    # Get indexes for selection
    cut = io.get_selection(infile,outfile,inputformat=inputformat,
                           cutcols=cutcols,mincuts=mincuts,maxcuts=maxcuts,
                           testing=testing,verbose=verbose)

    # Read the input data and correct it to the adequate units, etc.
    lms, lssfr, lzgas = get_sfrdata(infile,m_sfr_z,selection=cut,
                                    h0=h0,units_h0=units_h0,
                                    units_Gyr=units_Gyr,inoh = inoh,
                                    mtot2mdisk=mtot2mdisk,
                                    inputformat=inputformat,
                                    testing=testing,verbose=verbose)
    
    epsilon_param_z0 = [None]
    if infile_z0 is not None:
        epsilon_param_z0 = io.read_data(infile_z0,cut,
                                        inputformat=inputformat,
                                        params=agn_nH_params,
                                        testing=testing,
                                        verbose=verbose)
    
    # Modification of the stellar mass-metallicity relation
    if zeq is not None:
        minZ, maxZ = get_limits(propname='Z', photmod=photmod_sfr)
        ###here this does not work due to Lagn_param, need to change this dependency
        lzgas = correct_Z(zeq,lms,lzgas,minZ,maxZ,Lagn_param)

    # Characterise the HII regions from galaxy global properties
    lu_sfr, lnH_sfr = get_UnH_sfr(lms, lssfr, lzgas, outfile,
                                  q0=q0, z0=z0,gamma=gamma, T=T,
                                  epsilon_param_z0=epsilon_param_z0,
                                  IMF=IMF,model_nH_sfr=model_nH_sfr,
                                  model_U_sfr=model_U_sfr,verbose=verbose)

    if verbose: print(' U and nH calculated.')             
    lu_o_sfr = np.copy(lu_sfr)
    lnH_o_sfr = np.copy(lnH_sfr)
    lzgas_o_sfr = np.copy(lzgas)

    # Obtain spectral emission lines from HII regions
    nebline_sfr = get_lines(lu_sfr.T,lzgas.T,outfile,lnH=lnH_sfr.T,
                            photmod=photmod_sfr,origin='sfr',
                            verbose=verbose)

    # Change units into erg/s 
    if (photmod_sfr == 'gutkin16'):
        # Units: Lbolsun per unit SFR(Msun/yr) for 10^8yr, assuming Chabrier
        sfr = np.zeros(shape=np.shape(lssfr))
        for comp in range(ncomp):
            sfr[:,comp] = 10**(lms[:,comp]+lssfr[:,comp])
            nebline_sfr[comp,:,:] = nebline_sfr[comp,:,:]*c.Lbolsun*sfr[:,comp]

    if verbose:
        print(' Emission lines calculated.')

    # Write output
    extra_param = io.read_data(infile,cut,
                               inputformat=inputformat,
                               params=extra_params,
                               testing=testing,
                               verbose=verbose)        
    io.write_global_data(outfile,lms,lssfr=lssfr,
                      extra_param=extra_param,
                      extra_params_names=extra_params_names,
                      extra_params_labels=extra_params_labels,
                      verbose=verbose)

    io.write_sfr_data(outfile,lu_o_sfr,lnH_o_sfr,lzgas_o_sfr,
                      nebline_sfr,verbose=verbose)
    del lu_sfr, lnH_sfr
    del lu_o_sfr, lnH_o_sfr, lzgas_o_sfr
    del nebline_sfr

    # Get information on Mgas and scalelengths, if provided
    mgas = None; hr = None
    if mgas_r is not None:
        nattrs = io.add2header(outfile,['mgasr_type'],
                               [mgasr_type],verbose=verbose)
        mgas, hr = io.get_mgas_hr(infile,cut,
                                  mgas_r,r_type,
                                  h0=h0,units_h0=units_h0,
                                  inputformat=inputformat,
                                  testing=testing,verbose=verbose)
        lmgas = np.full(mgas.shape,c.notnum)
        ind = (mgas>0.)
        lmgas[ind] = np.log10(mgas[ind])
        io.write_global_data(outfile,lmgas.T,mass_type='gas',
                             extra_param=[hr.T],
                             extra_params_names=['h_gas'],
                             extra_params_labels=['Scalelength(Mpc)'],
                             verbose=verbose)
    
    #----------------NLR AGN calculation------------------------
    if AGN:
        if verbose: print('AGN:')
        # Add relevant constants to header
        names = ['model_spec_NLR','model_U_NLR','photmod_NLR',
                  'nH_NLR_cm3','T_NLR_K','r_NLR_Mpc','alpha_NLR','xid_NLR']
        values = [model_spec_agn,model_U_agn,photmod_agn,
                  nH_NLR,T_NLR,r_NLR,alpha_NLR,xid_NLR]
        nattrs = io.add2header(outfile,names,values,verbose=verbose)

        # Get the central metallicity
        if Z_correct_grad:
            # Get total mass for Z corrections 
            lm_tot = components2tot(lms)
            lzgas_agn = get_zgasagn(infile,Zgas_NLR,selection=cut,inoh=inoh,
                                    Z_correct_grad=True,lm_tot=lm_tot,
                                    inputformat=inputformat,
                                    testing=testing,verbose=verbose)
        else:
            lzgas_agn = get_zgasagn(infile,Zgas_NLR,selection=cut,
                                    inoh=inoh,inputformat=inputformat,
                                    testing=testing,verbose=verbose)

        # Get the AGN bolometric luminosity
        Lagn = get_Lagn(infile,cut,inputformat=inputformat,
                        params=Lagn_params,Lagn_inputs=Lagn_inputs,
                        h0=h0,units_h0=units_h0,
                        units_Gyr=units_Gyr,units_L=units_L,
                        testing=testing,verbose=verbose)

        # Get the ionising parameter, U, (and filling factor)
        lu_agn, epsilon_agn = get_UnH_agn(Lagn, mgas, hr,outfile,
                                          mgasr_type=mgasr_type,
                                          verbose=verbose)
        if verbose: print(' U calculated.')

        # Calculate emission lines in adequate units
        nebline_agn = get_lines(lu_agn,lzgas_agn,outfile,
                                photmod=photmod_agn,origin='NLR',
                                verbose=verbose)

        if (photmod_agn == 'feltre16'):
            # Units: erg/s for a central Lacc=10^45 erg/s
            nebline_agn[0] = nebline_agn[0]*Lagn/1e45
        if verbose: print(' Emission lines calculated.')

        # Write output in a file            
        io.write_agn_data(outfile,Lagn,lu_agn.T,lzgas_agn.T,
                          nebline_agn,
                          epsilon_agn=epsilon_agn,
                          verbose=verbose)             
        del lu_agn, lzgas_agn 
        del nebline_agn
    del lms, lssfr, cut

    if verbose:
        tt = round(time.perf_counter() - start_total_time,2)
        print('* Total time: ', tt, 's.')
    print(f'SUCCESS: {outfile}')
    return
