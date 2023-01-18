from get_nebular_emission.eml_io import get_data, write_data
from get_nebular_emission.eml_une import get_une
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_photio import get_lines, attenuation
import time
import numpy as np
#import get_nebular_emission.eml_testplots as get_testplot

def eml(infile, outfile, m_sfr_z, h0=None, att_param=['Rvir','ColdGas'],
        volume = 542.16**3.,IMF=['Kennicut','Top-heavy'], inputformat='HDF5',
        attmod='cardelli89',unemod='kashino20',photmod='gutkin16',
        LC2sfr=False, Zloh12=False, mtot2mdisk = True,
        verbose=False, Plotting=False, Testing=False):
    '''
    Calculate emission lines given the properties of model galaxies

    Parameters
    ----------
    infile : string
     - Name of the input file. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    m_sfr_z : list
     - [[component1_stellar_mass,sfr/LC,Z],[component2_stellar_mass,sfr/LC,Z],...]
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    h0 : float
      If not None: value of h, H0=100h km/s/Mpc.
    volume : float
      Carlton model default value = 500^3 Mpc^3/h^3. If not 500.**3. : value of the simulation volume in Mpc^3/h^3
    unemod : string
      Model to go from galaxy properties to U and ne
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
    This code returns an .hdf5 file with the electron density, the metallicity, 
    the ionization parameter, and the emission lines.

    '''
    
    first = True
    
    for i in range(len(infile)):
        
        start_time = time.perf_counter()
        
        # Read the input data and correct it to the adequate units, etc.
        lms, lssfr, loh12, cut = get_data(infile[i], m_sfr_z, h0=h0, 
                                     inputformat=inputformat, LC2sfr=LC2sfr, 
                                     Zloh12=Zloh12, mtot2mdisk=mtot2mdisk,
                                     IMF=IMF, verbose=verbose, 
                                     Plotting=Plotting, Testing=Testing)
        
        print('Data read.')
        
        # lms = lms[:1153192]
        # lssfr = lssfr[:1153192]
        # loh12 = loh12[:1153192]
        
        # if Testing and Plotting:
        #     lne = np.array([[2,2],[2,2]])
        #     lu = np.array([[-3.2,-3.2],[-1.7,-1.7]])
        #     loh12 = np.array([[np.log10(0.0008),np.log10(0.0008)],[np.log10(0.013),np.log10(0.013)]])
        # else:
        #     lu, lne, loh12 = get_une(lms, lssfr, loh12,
        #                     unemod=unemod, LC2sfr=LC2sfr, verbose=verbose,
        #                     Plotting=Plotting, Testing=Testing)
        
        lu, lne, loh12 = get_une(lms, lssfr, loh12,
                            unemod=unemod, LC2sfr=LC2sfr, verbose=verbose,
                            Plotting=Plotting, Testing=Testing)
            
        print('U and ne calculated.')
            
        limits = np.where((lu[:,0]>-4)&(lu[:,0]<-1)&(loh12[:,0]>-4)&(loh12[:,0]<-1.4)&(lu[:,0]!=const.notnum))[0]
        
        lms = lms[limits]
        lssfr = lssfr[limits]
        loh12 = loh12[limits]
        lu = lu[limits]
        lne = lne[limits]
        # From U and ne, obtain the emission lines from HII regions
        nebline = get_lines(lu,lne,loh12,photmod=photmod,verbose=verbose,
                            Testing=Testing,Plotting=Plotting)
        
        print('Lines calculated.')
        
        if Testing and Plotting:
            nebline_att = np.copy(nebline)
        else:
            nebline_att, coef_att, lines = attenuation(nebline, m_sfr_z, att_param=att_param, 
                                infile=infile[i], cut=cut, limits=limits, attmod=attmod, verbose=verbose)
            # nebline_att = np.copy(nebline)
            
        print('Attenuation calculated.')
            
        write_data(lms,lssfr,lu,lne,loh12,nebline[:],nebline_att,outfile,attmod=attmod,
                   unemod=unemod,photmod=photmod,first=first)
        
        if first:
            first = False
            
        if verbose:
            print()
            print('Subvolume', i+1, 'of', len(infile))
            print('Time:', round(time.perf_counter() - start_time,2), 's.')
            print()
            
        return cut, limits