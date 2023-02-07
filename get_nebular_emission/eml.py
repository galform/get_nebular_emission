from get_nebular_emission.eml_io import get_data, write_data
from get_nebular_emission.eml_une import get_une
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_photio import get_lines, attenuation, get_limits, clean_photarray
import time
import numpy as np
#import get_nebular_emission.eml_testplots as get_testplot

def eml(infile, outfile, m_sfr_z, h0=None, cutcols=[None], mincuts=[None], 
        maxcuts=[None], att_param=['Rvir','ColdGas'],
        volume = 542.16**3.,inputformat='HDF5',
        IMF_i=['Chabrier', 'Chabrier'], IMF_f=['Kroupa', 'Kroupa'], 
        attmod='cardelli89',unemod='kashino20',photmod='gutkin16',
        gamma=None, LC2sfr=False, cutlimits=False, mtot2mdisk = True,
        verbose=True, Plotting=False, Testing=False):
    '''
    Calculate emission lines given the properties of model galaxies

    Parameters
    ----------
    infile : string
     Name of the input file. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    inputformat : string
     Format of the input file.
    m_sfr_z : list
     - [[component1_stellar_mass,sfr/LC,Z],[component2_stellar_mass,sfr/LC,Z],...]
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    cutcols : list
     Parameters to look for cutting the data.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    mincuts : list
     Minimum value of the parameter of cutcols in the same index. All the galaxies below won't be considered.
    maxcuts : list
     Maximum value of the parameter of cutcols in the same index. All the galaxies above won't be considered.
    att_param : list
     Parameters to look for calculating attenuation. See eml_const to know what each model expects.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    IMF_i : list
     Assumed IMF in the input data.
     - [[component1_IMF],[component2_IMF],...]
    IMF_f : list
     Assumed IMF for the luminosity calculation. Please check the assumed IMF of the selected model for calculating U and ne.
     - [[component1_IMF],[component2_IMF],...]
    cutlimits : boolean
     If True the galaxies with U, ne and Z outside the photoionization model's grid limits won't be considered.
    h0 : float
      If not None: value of h, H0=100h km/s/Mpc.
    volume : float
      Carlton model default value = 500^3 Mpc^3/h^3. If not 500.**3. : value of the simulation volume in Mpc^3/h^3
    unemod : string
      Model to go from galaxy properties to U and ne
    attmod : string
      Attenuation model.
    photmod : string
      Photoionisation model to be used for look up tables.
    LC2sfr : boolean
      If True magnitude of Lyman Continuum photons expected as input for SFR.
    mtot2mdisk : boolean
      If True transform the total mass into the disk mass. disk mass = total mass - bulge mass.
    verbose : boolean
      If True print out messages
    Plotting : boolean
      If True run verification plots with all data.
    Testing : boolean
      If True only run over few entries for testing purposes

    Notes
    -------
    This code returns an .hdf5 file with the mass, specific star formation rate,
    electron density, metallicity, ionization parameter, and the emission lines.

    '''
    
    first = True
    
    start_total_time = time.perf_counter()
    
    for i in range(len(infile)):
        
        start_time = time.perf_counter()
        
        # Read the input data and correct it to the adequate units, etc.
        lms, lssfr, loh12, cut = get_data(infile[i], m_sfr_z, h0=h0,
                                     cutcols=cutcols, mincuts=mincuts, maxcuts=maxcuts,
                                     inputformat=inputformat, LC2sfr=LC2sfr, 
                                     mtot2mdisk=mtot2mdisk,
                                     IMF_i=IMF_i, IMF_f=IMF_f, verbose=verbose, 
                                     Plotting=Plotting, Testing=Testing)
        
        # return lms
        
        if verbose:
            print('Data read.')
        
        # if Testing and Plotting:
        #     lne = np.array([[2,2],[2,2]])
        #     lu = np.array([[-3.2,-3.2],[-1.7,-1.7]])
        #     loh12 = np.array([[np.log10(0.0008),np.log10(0.0008)],[np.log10(0.013),np.log10(0.013)]])
        # else:
        #     lu, lne, loh12 = get_une(lms, lssfr, loh12,
        #                     unemod=unemod, LC2sfr=LC2sfr, verbose=verbose,
        #                     Plotting=Plotting, Testing=Testing)
        
        lu, lne, loh12 = get_une(lms, lssfr, loh12,
                            unemod=unemod, gamma=gamma, LC2sfr=LC2sfr, verbose=verbose,
                            Plotting=Plotting, Testing=Testing)
        
        if verbose:
            print('U and ne calculated.')
        
        lms, lssfr, lu, lne, loh12, limits = clean_photarray(lms, lssfr, lu, lne, loh12, cutlimits=cutlimits, photmod=photmod)
        
        nebline = get_lines(lu,lne,loh12,photmod=photmod,verbose=verbose,
                            Testing=Testing,Plotting=Plotting)
        
        if verbose:
            print('Lines calculated.')
        
        nebline_att, coef_att = attenuation(nebline, m_sfr_z, att_param=att_param, 
                                  infile=infile[i], inputformat=inputformat, 
                                  cut=cut, limits=limits, cutlimits=cutlimits, 
                                  attmod=attmod, verbose=verbose)
        
        if verbose:
            print('Attenuation calculated.')
            
        write_data(lms,lssfr,lu,lne,loh12,nebline,nebline_att,outfile,attmod=attmod,
                   unemod=unemod,photmod=photmod,first=first)
        
        del lms, lssfr, lu, lne, loh12, nebline, nebline_att, cut, limits
        
        time.sleep(1)
        
        if first:
            first = False
            
        if verbose:
            print()
            print('Subvolume', i+1, 'of', len(infile))
            print('Time:', round(time.perf_counter() - start_time,2), 's.')
            print()
            
    if verbose:
        print('Total time: ', round(time.perf_counter() - start_total_time,2), 's.')