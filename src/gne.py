import src.gne_io as io
from src.gne_une import get_une, bursttobulge, calculate_ng_hydro_eq, Z_blanc, Z_tremonti, Z_tremonti2, n_ratio
from src.gne_agn import L_agn
import src.gne_const as const
from src.gne_photio import get_lines, get_limits, clean_photarray, calculate_flux
from src.gne_att import attenuation
import time
import numpy as np
from src.gne_plots import make_testplots

def gne(infile, redshift, m_sfr_z,
        h0,omega0,omegab,lambda0,vol,
        inputformat='HDF5',infile_z0=[None], h0units=True, 
        cutcols=[None], mincuts=[None], maxcuts=[None], 
        att=False, att_params=None, att_ratio_lines=None,
        flux=False,
        flag=0,
        IMF_i=['Kroupa', 'Kroupa'], IMF_f=['Kroupa', 'Kroupa'], 
        q0=const.q0_orsi, z0=const.Z0_orsi, gamma=1.3,
        T=10000,
        AGN=False, AGNinputs='Lagn', Lagn_params=None, Z_central_cor=False,
        epsilon_params=None,
        extra_params=None, extra_params_names=None, extra_params_labels=None,
        attmod='cardelli89',
        unemod_sfr='kashino19', unemod_agn='panuzzo03',
        photmod_sfr='gutkin16', photmod_agn='feltre16',
        oh12=False, LC2sfr=False,
        cutlimits=False, mtot2mdisk=True,
        verbose=True, testing=False,
        xid_feltre=0.5,alpha_feltre=-1.7,
        xid_gutkin=0.3,co_gutkin=1,imf_cut_gutkin=100):
    '''
    Calculate emission lines given the properties of model galaxies

    Parameters
    ----------
    infile : strings
     List with the name of the input files. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    outfile : string
     Name of the output file.
    m_sfr_z : list
     - [[component1_stellar_mass,sfr/LC,Z],[component2_stellar_mass,sfr/LC,Z],...]
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    inputformat : string
     Format of the input file.
    infile_z0 : strings
     List with the name of the input files with the galaxies at redshift 0. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    h0 : float
      If not None: value of h, H0=100h km/s/Mpc.
    redshift : float
     Redshift of the input data.
    cutcols : list
     Parameters to look for cutting the data.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    mincuts : floats
     Minimum value of the parameter of cutcols in the same index. All the galaxies below won't be considered.
    maxcuts : floats
     Maximum value of the parameter of cutcols in the same index. All the galaxies above won't be considered.
    att : boolean
     If True calculates attenuated emission.
    att_params : list
     Parameters to look for calculating attenuation. See gne_const to know what each model expects.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    att_ratio_lines : strings
     Names of the lines corresponding to the values in att_params when attmod=ratios.
     They should be written as they are in the selected model (see gne_const).
    flux : boolean
     If True calculates flux of the emission lines based on the given redshift.
    IMF_i : strings
     Assumed IMF in the input data.
     - [[component1_IMF],[component2_IMF],...]
    IMF_f : strings
     Assumed IMF for the luminosity calculation. Please check the assumed IMF of the selected model for calculating U and ne.
     - [[component1_IMF],[component2_IMF],...]
    q0 : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    z0 : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    gamma : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    T : float
     Typical temperature of ionizing regions.
    AGN : boolean
     If True calculates emission from the narrow-line region of AGNs.
    AGNinputs : string
     Type of inputs for AGN's bolometric luminosity calculations.
    Lagn_params : list
     Inputs for AGN's bolometric luminosity calculations.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    Z_central_correction : boolean
     If False, the code supposes the central metallicity of the galaxy to be the mean one.
     If True, the code estimates the central metallicity of the galaxy from the mean one.
    epsilon_params : list
     Inputs for the calculation of the volume-filling factor.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    extra_params : list
     Parameters from the input files which will be saved in the output file.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    extra_params_names : strings
     Names of the datasets in the output files for the extra parameters.
    extra_params_labels : strings
     Description labels of the datasets in the output files for the extra parameters.
    attmod : string
     Attenuation model.
    unemod_sfr : string
     Model to go from galaxy properties to U and ne.
    unemod_agn : string
     Model to go from galaxy properties to U and ne.
    photmod_sfr : string
     Photoionisation model to be used for look up tables.
    photmod_agn : string
     Photoionisation model to be used for look up tables.
    LC2sfr : boolean
     If True magnitude of Lyman Continuum photons expected as input for SFR.
    cutlimits : boolean
     If True the galaxies with U, ne and Z outside the photoionization model's grid limits won't be considered.
    mtot2mdisk : boolean
     If True transform the total mass into the disk mass. disk mass = total mass - bulge mass.
    verbose : boolean
     If True print out messages.
    testing : boolean
     If True only run over few entries for testing purposes.
    xid_feltre : float
     Dust-to-metal ratio for the Feltre et. al. photoionisation model.
    alpha_feltre : float
     Alpha value for the Feltre et. al. photoionisation model.
    xid_gutkin : float
     Dust-to-metal ratio for the Gutkin et. al. photoionisation model.
    co_gutkin : float
     C/O ratio for the Gutkin et. al. photoionisation model.
    imf_cut_gutkin : float
     Solar mass high limit for the IMF for the Gutkin et. al. photoionisation model.
    
    Notes
    -------
    This code returns an .hdf5 file with the mass, specific star formation rate,
    electron density, metallicity, ionization parameter, and the emission lines.

    '''

    # Generate header in the output file from input
    outfile = io.generate_header(infile, redshift,
                                 h0,omega0,omegab,lambda0,vol,
                                 unemod_sfr=unemod_sfr, unemod_agn=unemod_agn,
                                 photmod_sfr=photmod_sfr, photmod_agn=photmod_agn,
                                 attmod=attmod,verbose=verbose)

    # Time and passes variables
    first = True
    start_total_time = time.perf_counter()
    start_time = time.perf_counter()

    # Read the input data and correct it to the adequate units, etc.
    lms, lssfr, loh12, cut = io.get_data(infile, outfile,
                                         m_sfr_z, h0units=h0units,
                                         cutcols=cutcols, mincuts=mincuts,
                                         maxcuts=maxcuts,
                                         inputformat=inputformat,
                                         oh12 = oh12,
                                         LC2sfr=LC2sfr,mtot2mdisk=mtot2mdisk,
                                         IMF_i=IMF_i, IMF_f=IMF_f,
                                         verbose=verbose,testing=testing)

    epsilon_param, epsilon_param_z0, Lagn_param, att_param, extra_param = io.get_secondary_data(
        infile,cut,infile_z0=infile_z0,epsilon_params=epsilon_params,
        extra_params=extra_params,Lagn_params=Lagn_params,
        att_params=att_params,inputformat=inputformat,
        attmod=attmod,verbose=verbose)

    # Modification of the stellar mass-metallicity relation
    # 0 for no correction
    if flag==1:
        loh12 = Z_tremonti(lms,loh12,Lagn_param)[1]
    elif flag==2:
        minZ, maxZ = get_limits(propname='Z', photmod=photmod_sfr)
        loh12 = Z_tremonti2(lms,loh12,minZ,maxZ,Lagn_param)
            
    Q_sfr, lu_sfr, lne_sfr, loh12_sfr, epsilon_sfr, ng_ratio = \
        get_une(lms, lssfr, loh12, outfile,
                q0=q0, z0=z0, T=T, IMF_f=IMF_f, h0units=h0units,
                epsilon_param=epsilon_param,
                epsilon_param_z0=epsilon_param_z0,
                origin='sfr', unemod=unemod_sfr,
                gamma=gamma, verbose=verbose)

    if verbose:
        print('SF:')
        print(' U and ne calculated.')
            
    lu_o_sfr = np.copy(lu_sfr)
    lne_o_sfr = np.copy(lne_sfr)
    loh12_o_sfr = np.copy(loh12_sfr)

    clean_photarray(lms, lssfr, lu_sfr, lne_sfr, loh12_sfr, photmod=photmod_sfr)

    nebline_sfr = get_lines(lu_sfr,lne_sfr,loh12_sfr,photmod=photmod_sfr,
                            verbose=verbose,
                            xid_gutkin=xid_gutkin,co_gutkin=co_gutkin,imf_cut_gutkin=imf_cut_gutkin)

    for comp in range(len(m_sfr_z)):
        nebline_sfr[comp] = nebline_sfr[comp]*3.826e33*10**(lms[:,comp]+lssfr[:,comp])
        
    if verbose:
        print(' Emission calculated.')
            
    if att:
        nebline_sfr_att, coef_sfr_att = attenuation(nebline_sfr, att_param=att_param, 
                                      att_ratio_lines=att_ratio_lines,redshift=redshift,
                                      origin='sfr',
                                      cut=cut, attmod=attmod, photmod=photmod_sfr,verbose=verbose)
        
        if verbose:
            print(' Attenuation calculated.')
    else:
        nebline_sfr_att = np.array(None)

    if flux:
        fluxes_sfr = calculate_flux(nebline_sfr,outfile,
                                    h0units=h0units,origin='sfr')
        fluxes_sfr_att = calculate_flux(nebline_sfr_att,outfile,
                                        h0units=h0units,origin='sfr')
        if verbose:
            print(' Flux calculated.')
    else:
        fluxes_sfr = np.array(None)
        fluxes_sfr_att = np.array(None)

    io.write_sfr_data(outfile,lms,lssfr,lu_o_sfr,lne_o_sfr,loh12_o_sfr,
                      nebline_sfr,nebline_sfr_att,
                      fluxes_sfr,fluxes_sfr_att,
                      extra_param=extra_param,
                      extra_params_names=extra_params_names,
                      extra_params_labels=extra_params_labels,
                      first=first,verbose=verbose)
    
    del lu_sfr, lne_sfr, loh12_sfr
    del lu_o_sfr, lne_o_sfr, loh12_o_sfr
    del nebline_sfr, nebline_sfr_att
        
    if AGN:
        bursttobulge(lms, Lagn_param)
        
        Lagn = L_agn(Lagn_param,AGNinputs=AGNinputs,
                     verbose=verbose)
        
        Q_agn, lu_agn, lne_agn, loh12_agn, epsilon_agn, ng_ratio = \
            get_une(lms,lssfr, loh12, outfile, q0=q0, z0=z0,
                    Z_central_cor=Z_central_cor,Lagn=Lagn, T=T,
                    epsilon_param=epsilon_param,h0units=h0units,
                    IMF_f=IMF_f, origin='agn',
                    unemod=unemod_agn, gamma=gamma, verbose=verbose)
        
        if verbose:
            print('AGN:')
            print(' U and ne calculated.')
        
        lu_o_agn = np.copy(lu_agn)
        lne_o_agn = np.copy(lne_agn)
        loh12_o_agn = np.copy(loh12_agn) 
            
        clean_photarray(lms, lssfr, lu_agn, lne_agn, loh12_agn, photmod=photmod_agn)
            
        nebline_agn = get_lines(lu_agn,lne_agn,loh12_agn,photmod=photmod_agn,verbose=verbose,
                            xid_feltre=xid_feltre,alpha_feltre=alpha_feltre)
        nebline_agn[0] = nebline_agn[0]*Lagn/1e45
        
        if verbose:
            print(' Emission calculated.')

        if att:
            nebline_agn_att, coef_agn_att = attenuation(nebline_agn, att_param=att_param, 
                                          att_ratio_lines=att_ratio_lines,redshift=redshift,
                                          origin='agn',
                                          cut=cut, attmod=attmod, photmod=photmod_agn,verbose=verbose)
            if verbose:
                print(' Attenuation calculated.')     
        else:
            nebline_agn_att = np.array(None)
            
        if flux:
            fluxes_agn = calculate_flux(nebline_agn,outfile,
                                        h0units=h0units,origin='agn')
            fluxes_agn_att = calculate_flux(nebline_agn_att,outfile,
                                            h0units=h0units,origin='agn')
            if verbose:
                print(' Flux calculated.')
        else:
            fluxes_agn = np.array(None)
            fluxes_agn_att = np.array(None)

        io.write_agn_data(outfile,lu_o_agn,lne_o_agn,loh12_o_agn,
                          nebline_agn,nebline_agn_att,
                          fluxes_agn,fluxes_agn_att,
                          epsilon_agn,
                          first=first,verbose=verbose)             
        del lu_agn, lne_agn, loh12_agn 
        del lu_o_agn, lne_o_agn, loh12_o_agn
        del nebline_agn, nebline_agn_att


    del lms, lssfr, cut
    if first:
        first = False
    
    if verbose:
        print('* Total time: ', round(time.perf_counter() - start_total_time,2), 's.')
