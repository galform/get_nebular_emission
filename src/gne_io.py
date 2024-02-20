"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
"""
import h5py
import sys
import os
import numpy as np
import src.gne_const as const
import math


def stop_if_no_file(infile):
    '''
    It stops the program if a file does not exists

    Parameters
    -------
    infile : string
        Input file
    '''
    
    if (not os.path.isfile(infile)):
        print('STOP: no input file {}'.format(infile)) 
        sys.exit()
    return


def check_file(infile,verbose=False):
    '''
    It checks if a file exists

    Parameters
    -------
    infile : string
        Input file
    verbose : boolean
        If True print out messages

    Returns
    -------
    file_fine : boolean
        True when the file exists.
    '''
    
    file_fine = True  
    if (not os.path.isfile(infile)):
        file_fine = False
        if verbose:
            print('WARNING (gne_io.check_file): file not found {}'.
                  format(infile))

    return file_fine


def create_dir(outdir):
    '''
    Return True if directory already exists or it has been created
    '''
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except:
            print('WARNING (iotools.create_dir): problem creating directory ',outdir)
            return False
    return True


def get_outnom(filenom,ftype='data',verbose=False):
    '''
    Get output from a given filename

    Parameters
    -------
    filenom : string
        Name of file
    ftype : string
        Type of the file, data or plots
    verbose : boolean
        If True print out messages

    Returns
    -------
    outfile : string
        Path to output file
    '''

    nom = os.path.splitext(filenom.split('/')[-1])[0]

    dirf = 'output/' + ftype + '/'
    create_dir(dirf)    

    outfile = dirf + nom + '.hdf5'

    if verbose:
        print(f'* Output {ftype}: {outfile}')
    return outfile


def get_nheader(infile,firstchar=None):
    '''
    Given a text file with a structure: header+data, 
    counts the number of header lines

    Parameters
    -------
    infile : string
        Input file

    Returns
    -------
    ih : integer
        Number of lines with the header text
    '''


    ih = 0
    with open(infile,'r') as ff:
        for line in ff:
            if not line.strip():
                # Count any empty lines in the header
                ih += 1
            else:
                sline = line.strip()
                
                # Check that the first character is not a digit
                char1 = sline[0]
                word1 = sline.split()[0]
                if not firstchar:
                    if (not char1.isdigit()):
                        if (char1 != '-'):
                            ih += 1
                        else:
                            try:
                                float(word1)
                                return ih
                            except:
                                ih += 1
                    else:
                        return ih
                else:
                    if char1 == firstchar:
                        ih+=1
    return ih
        


def get_ncomponents(cols):
    '''
    Get the number of components to estimate the emission lines from

    Parameters
    ----------
    cols : list
      List of columns with M*, SFR, Z, per components

    Returns
    -------
    ncomp : integer
      Number of components (for example 2 for bulge and disk)
    '''
    
    ncomp = 1
    
    try:
        dum = np.shape(cols)[1]
        ncomp = np.shape(cols)[0]
    except:
        ncomp = 1
        print('STOP (gne_io.get_ncomponents): ',
              'Columns should be given as m_sfr_z=[[0,1,2]]')
        sys.exit()
        
    return ncomp


def locate_interval(val, edges):
    '''
    Get the index, i, of the interval with edges.

    Parameters
    ----------
    val : int or float
        Value
    edges : array of int or floats
        Array of the edges of the intervals
        
    Returns
    -------
    ind : integer
        Index of the interval. If outside: index of the limits
    '''

    ind = const.notnum
    low = np.asarray(edges[:-1])
    high = np.asarray(edges[1:])

    if (val >= high[-1] and val >= low[-1]):
        ind = len(high)
    elif (val<=low[0]):
        ind = 0
    else:
        linds = np.where(val >= low)
        hinds = np.where(val < high)

        if (np.shape(linds)[1] > 0 and np.shape(hinds)[1] > 0):
            lind = linds[0]
            hind = hinds[0]
            common = list(set(lind).intersection(hind))

            if (len(common) == 1):
                ind = common[0]

    return ind



def read_data(infile, cols, cutcols=[None], mincuts=[None], maxcuts=[None],
              inputformat='hdf5',testing=False, verbose=True):
    '''
    It reads star masses, star formation rates and metallicities from a file.

    Parameters
    ----------
    infile : string
     - Name of the input file. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    inputformat : string
     Format of the input file.
    cols : list
     - [[component1_stellar_mass,sfr,Z],[component2_stellar_mass,sfr,Z],...]
     - Expected : component1 = total or disk, component2 = bulge
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
    attmod : string
      Model of dust attenuation.
    verbose : boolean
      If True print out messages
    testing : boolean
      If True only run over few entries for testing purposes

    Returns
    -------
    lms, lssfr, loh12 : floats
    cut : integers
    '''
    
    check_file(infile, verbose=verbose)

    ncomp = get_ncomponents(cols)
    
    if testing:
        limit = 50
    else:
        limit = None
        
    if inputformat not in const.inputformats:
        if verbose:
            print('STOP (gne_io): Unrecognised input format.',
                  'Possible input formats = {}'.format(const.inputformats))
        sys.exit()
    elif inputformat=='hdf5':
        with h5py.File(infile, 'r') as f:
            hf = f['data']
            
            cut = np.arange(len(hf[cols[0][0]][:limit]))
            
            for i in range(len(cutcols)):
                if cutcols[i]:
                    param = hf[cutcols[i]][:limit]
                    mincut = mincuts[i]
                    maxcut = maxcuts[i]
                    if mincut and maxcut:
                        cut = np.intersect1d(cut,np.where((mincut<param)&(param<maxcut))[0])
                    elif mincut:
                        cut = np.intersect1d(cut,np.where(mincut<param)[0])
                    elif maxcut:
                        cut = np.intersect1d(cut,np.where(param<maxcut)[0])
                
            for i in range(ncomp):
                if i==0:
                    lms = np.array([hf[cols[i][0]][:limit]])
                    lssfr = np.array([hf[cols[i][1]][:limit]])
                    loh12 = np.array([hf[cols[i][2]][:limit]])
                else:
                    lms = np.append(lms,[hf[cols[i][0]][:limit]],axis=0)
                    lssfr = np.append(lssfr,[hf[cols[i][1]][:limit]],axis=0)
                    loh12 = np.append(loh12,[hf[cols[i][2]][:limit]],axis=0)
    elif inputformat=='txt':
        ih = get_nheader(infile)
        
        cut = np.arange(len(np.loadtxt(infile,usecols=cols[0],skiprows=ih)[:limit]))
        
        if cutcols[0]:
            for i in range(len(cutcols)):
                
                param = np.loadtxt(infile,usecols=cutcols[i],skiprows=ih)[:limit]
                mincut = mincuts[i]
                maxcut = maxcuts[i]
                
                if mincut and maxcut:
                    cut = np.intersect1d(cut,np.where((mincut<param)&(param<maxcut))[0])
                elif mincut:
                    cut = np.intersect1d(cut,np.where(mincut<param)[0])
                elif maxcut:
                    cut = np.intersect1d(cut,np.where(param<maxcut)[0])
                    
        
        for i in range(ncomp):
            X = np.loadtxt(infile,usecols=cols[i],skiprows=ih).T[:,:limit]
            
            if i==0:
                lms = np.array([X[0]])
                lssfr = np.array([X[1]])
                loh12 = np.array([X[2]])
            else:
                lms = np.append(lms,[X[0]],axis=0)
                lssfr = np.append(lssfr,[X[1]],axis=0)
                loh12 = np.append(loh12,[X[2]],axis=0)
    else:
        if verbose:
            print('STOP (gne_io.read_data): ',
                  'Input file has not been found.')
        sys.exit()
        
    lms = lms.T
    lssfr = lssfr.T
    loh12 = loh12.T
            
    return lms[cut], lssfr[cut], loh12[cut], cut



def get_secondary_data(infile, cut, infile_z0=None, epsilon_params=None, 
                       Lagn_params=None, att_params=None, extra_params=None,
                       inputformat='hdf5', attmod='cardelli89', verbose=True):    
    '''
    Get data for epsilon calculation in the adecuate units.
    
    Parameters
    ----------
    infile : string
     - Name of the input file. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    infile_z0 : string
     Name of the files with the galaxies at redshift 0. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    cut : strings
     List of indexes of the selected galaxies from the samples.
    inputformat : string
     Format of the input file.
    epsilon_params : list
     Inputs for epsilon calculation (parameter for Panuzzo 2003 nebular region model).
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    Lagn_params : list
     Inputs for AGN's bolometric luminosity calculations.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    attmod : string
     Attenuation model.
    verbose : boolean
     If True print out messages.
    extra_params : list
     Parameters from the input files which will be saved in the output file.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
     
    Returns
    -------
    epsilon_param, epsilon_param_z0, Lagn_param, att_param, extra_param : floats
    '''
    
    epsilon_param = [[None]]
    epsilon_param_z0 = [[None]]
    Lagn_param = [[None]]
    extra_param = [[None]]
    
    if inputformat not in const.inputformats:
        if verbose:
            print('STOP (gne_io): Unrecognised input format.',
                  'Possible input formats = {}'.format(const.inputformats))
        sys.exit()
    elif inputformat=='hdf5':
        if verbose:
            print('HDF5 not implemented yet for secondary params.')
        sys.exit()
    elif inputformat=='txt':
        ih = get_nheader(infile)
        
        if epsilon_params:
            epsilon_param = np.loadtxt(infile,skiprows=ih,usecols=epsilon_params)[cut].T
            
        if infile_z0[0]:
            epsilon_param_z0 = np.loadtxt(infile_z0,skiprows=ih,usecols=epsilon_params)[cut].T

        if Lagn_params:
            Lagn_param = np.loadtxt(infile,skiprows=ih,usecols=Lagn_params)[cut].T
            
        if extra_params:
            extra_param = np.loadtxt(infile,skiprows=ih,usecols=extra_params)[cut].T
            if len(extra_params)==1:
                extra_param = np.array([extra_param])
        
        if att_params:
                att_param = np.loadtxt(infile,skiprows=ih,usecols=att_params)[cut].T
                
    return epsilon_param, epsilon_param_z0, Lagn_param, att_param, extra_param



def get_data(infile, cols, h0=None, inputformat='hdf5', 
             IMF_i=['Chabrier', 'Chabrier'], IMF_f=['Kroupa', 'Kroupa'], 
             cutcols=None, mincuts=[None], maxcuts=[None],
             attmod='GALFORM', LC2sfr=False, mtot2mdisk=True, 
             verbose=False, testing=False):
    '''
    Get Mstars, sSFR and (12+log(O/H)) in the adecuate units.

    Parameters
    ----------
    infile : strings
     List with the name of the input files. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    inputformat : string
     Format of the input file.
    cols : list
     - [[component1_stellar_mass,sfr,Z],[component2_stellar_mass,sfr,Z],...]
     - Expected : component1 = total or disk, component2 = bulge
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    cutcols : list
     Parameters to look for cutting the data.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    mincuts : strings
     Minimum value of the parameter of cutcols in the same index. All the galaxies below won't be considered.
    maxcuts : strings
     Maximum value of the parameter of cutcols in the same index. All the galaxies above won't be considered.
    attmod : string
     Attenuation model.
    IMF_i : strings
     Assumed IMF in the input data.
     - [[component1_IMF],[component2_IMF],...]
    IMF_f : strings
     Assumed IMF for the luminosity calculation. Please check the assumed IMF of the selected model for calculating U and ne.
     - [[component1_IMF],[component2_IMF],...]
    h0 : float
      If not None: value of h, H0=100h km/s/Mpc.
    LC2sfr : boolean
      If True magnitude of Lyman Continuum photons expected as input for SFR.
    mtot2mdisk : boolean
      If True transform the total mass into the disk mass. disk mass = total mass - bulge mass.
    verbose : boolean
      If True print out messages
    testing : boolean
      If True only run over few entries for testing purposes

    Returns
    -------
    lms, lssfr, loh12 : floats
    cut : integers
    '''
    
    lms,lssfr,loh12,cut = read_data(infile, cols=cols, cutcols=cutcols,
                                    maxcuts=maxcuts, mincuts=mincuts,
                                    inputformat=inputformat, 
                                    testing=testing, verbose=verbose)

    ncomp = get_ncomponents(cols)

    # Set to a default value if negative stellar masses
    ind = np.where(lms<=1.)
    lms[ind] = const.notnum
    lssfr[ind] = const.notnum
    loh12[ind] = const.notnum

    if LC2sfr: # Avoid positives magnitudes of LC photons
        ind = np.where(lssfr>0)
        lssfr[ind] = const.notnum ; loh12[ind] = const.notnum


    else: # Avoid other negative SFR
        ind = np.where(lssfr<=0)
        lssfr[ind] = const.notnum ; loh12[ind] = const.notnum



    ind = np.where(loh12<=0) # Avoid other negative Z
    loh12[ind] = const.notnum

    # Calculate the disk mass if we have only the total and bulge mass

    if mtot2mdisk:
        if ncomp!=2:
            if verbose:
                print('STOP (gne_io.get_data): ',
                      'mtot2mdisk can only be True with two components.')
            sys.exit()
                
        lms_tot = lms[:,0]

        # Calculate the disk mass :
        lmsdisk = lms[:,0] - lms[:,1]
        lms = np.column_stack((lmsdisk,lms[:,1]))     

        # Take the log of the total stellar mass
        ind = np.where(lms_tot > 0.)
        lms_tot[ind] = np.log10(lms_tot[ind])

        # Take the log of the stellar mass:
        ind = np.where(lms > 0.)
        lms[ind] = np.log10(lms[ind])

    else:
        if ncomp!=1:
            lms_tot = np.sum(lms,axis=1)
    
            # Take the log of the total stellar mass:
            ind = np.where(lms_tot > 0.)
            lms_tot[ind] = np.log10(lms_tot[ind])
    
        # Take the log of the stellar mass:
        ind = np.where(lms > 0.)
        lms[ind] = np.log10(lms[ind])

    # Obtain log10(sSFR) in 1/yr and calculate SFR from LC photons if necessary

    if LC2sfr:
        
        for comp in range(ncomp):
            ins = np.zeros(len(lssfr))
            ind = np.where(lssfr[:, comp] != const.notnum)
            ins[ind] = 1.02*(10.**(-0.4*lssfr[ind,comp]-4.))
            ins[ind] = ins[ind]*(const.IMF_SFR[IMF_i[comp]]/const.IMF_M[IMF_i[comp]])*(const.IMF_M[IMF_f[comp]]/const.IMF_SFR[IMF_f[comp]])
            ind = np.where(ins > 0)
            lssfr[ind,comp] = np.log10(ins[ind]) - lms[ind,comp] - 9.
            ind = np.where(ins < 0)
            lssfr[ind,comp] = const.notnum
    
            ind = np.where(lssfr[:, comp] == const.notnum)
            lssfr[ind,comp] = const.notnum
            
        if ncomp!=1:
            lssfr_tot = np.zeros(len(lssfr))
            ssfr = np.zeros(lssfr.shape)
            
            for comp in range(ncomp):
                ind = np.where(lssfr[:,comp] != const.notnum)
                ssfr[ind,comp] = 10. ** lssfr[ind,comp]

            ins = np.sum(ssfr,axis=1)
            ind = np.where(ins > 0)
            lssfr_tot[ind] = np.log10(ins[ind])

    # If instantaneous SFR as input:
    else:
        # Take the log of the ssfr:
        # ind = [np.where(lssfr[:,0] > 0.)[0], np.where(lssfr[:,1] > 0.)[0]]
        for comp in range(ncomp):
            ind = np.where(lssfr[:,comp] > 0.)[0]
            lssfr[ind,comp] = lssfr[ind,comp]*(const.IMF_SFR[IMF_i[comp]]/const.IMF_M[IMF_i[comp]])*(const.IMF_M[IMF_f[comp]]/const.IMF_SFR[IMF_f[comp]])
            lssfr[ind,comp] = np.log10(lssfr[ind,comp]) - lms[ind,comp] - 9.

        if ncomp!=1:
            lssfr_tot = np.zeros(len(lssfr))
            ssfr = np.zeros(lssfr.shape)
            for comp in range(ncomp):
                ind = np.where(lssfr[:,comp]!=const.notnum)
                ssfr[ind,comp] = 10.**(lssfr[ind,comp])

            ins = np.sum(ssfr,axis=1)
            ind = np.where(ins>0)
            lssfr_tot[ind] = np.log10(ins[ind])

    if h0:
        # Correct the units of the stellar mass
        ind = np.where(lms > 0.)
        lms[ind] = lms[ind] - np.log10(h0)
        if ncomp!=1:
            lms_tot = lms_tot - np.log10(h0)

    if ncomp!=1:
        lsfr = lssfr_tot+lms_tot
    else:
        lsfr = lssfr + lms

    # Obtain 12+log10(O/H) from Z=MZcold/Mcold
    ind = np.where(loh12>0)
    loh12[ind] = np.log10(loh12[ind]) #+ const.ohsun - np.log10(const.zsun)

    if ncomp!=1:
        oh12 = np.zeros(loh12.shape)
        loh12_tot = np.zeros(len(loh12))

        for comp in range(ncomp):
            ind = np.where(loh12[:,comp] != const.notnum)
            oh12[ind,comp] = 10. ** (loh12[ind,comp])
    
        ins = np.sum(oh12,axis=1)
        ind = np.where(ins>0)
        loh12_tot[ind] = np.log10(ins[ind])



    if testing: # here : Search more efficient form. Allow more components in the header
        if ncomp==2:
            header1 = 'log(mstars_tot), log(mstars_disk), log(mstars_bulge),' \
                      ' log(SFR_tot), log(sSFR_tot), log(sSFR_disk), log(sSFR_bulge) ' \
                      'log(Z)_tot, log(Z)_disk, log(Z)_bulge'
            datatofile=np.column_stack((lms_tot,lms,lsfr,lssfr_tot,lssfr,loh12_tot,loh12))
        elif ncomp==1:
            header1 = 'log(mstars), ' \
                      'log(SFR), log(sSFR) ' \
                      'log(Z)'
            datatofile=np.column_stack((lms,lsfr,lssfr,loh12))

        if LC2sfr:
            outfil = r"example_data/tmp_LC.dat"
        else:
            outfil = r"example_data/tmp_avSFR.dat"

        with open(outfil, 'w') as outf:
            np.savetxt(outf, datatofile, delimiter=' ', header=header1)
            outf.closed


                    
    return lms,lssfr,loh12,cut

def write_data(lms,lssfr,lu_sfr,lne_sfr,loh12_sfr,
               nebline_sfr,nebline_sfr_att=None,fluxes_sfr=None,fluxes_sfr_att=None,
               extra_param=[[None]],extra_params_names=None,extra_params_labels=None,
               outfile='output.hdf5',attmod='ratios',
               unemod_sfr='kashino20',photmod_sfr='gutkin16',first=True):
    '''
    Create a .hdf5 file from a .dat file.

    Parameters
    ----------
    lms_sfr : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr_sfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lu_sfr : floats
     U of the galaxies per component.
    lne_sfr : floats
     ne of the galaxies per component (cm^-3).
    loh12_sfr : floats
     Metallicity of the galaxies per component (12+log(O/H))
    nebline_sfr : floats
      Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
    nebline_sfr_att : floats
      Array with the luminosity of the attenuated lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)    
    extra_params : list
     Parameters from the input files which will be saved in the output file.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    extra_params_names : strings
     Names of the datasets in the output files for the extra parameters.
    extra_params_labels : strings
     Description labels of the datasets in the output files for the extra parameters.
    outfile : string
      Name of the output file.
    attmod : string
      Attenuation model.
    unemod_sfr : string
      Model to go from galaxy properties to U and ne.
    photmod_sfr : string
      Photoionisation model to be used for look up tables.
    first : boolean
      If True it creates the HDF5 file (first subvolume). If false, it adds elements to the existing one.
    '''
    
    if first: 
        with h5py.File(outfile,'w') as hf:
            head = hf.create_dataset('header',(1,))
            head.attrs[u'HII model'] = unemod_sfr
            head.attrs[u'Lines model for SF'] = photmod_sfr
            head.attrs[u'Attenuation model'] = attmod
    
            # Data
            hfdat = hf.create_group('data')
            
            hfdat.create_dataset('lms', data=lms, maxshape=(None,None))
            hfdat['lms'].dims[0].label = 'log10(M*) (Msun)'
            
            hfdat.create_dataset('lssfr', data=lssfr, maxshape=(None,None))
            hfdat['lssfr'].dims[0].label = 'log10(SFR/M*) (1/yr)'
    
            hfdat.create_dataset('lu_sfr', data=lu_sfr, maxshape=(None,None))
            hfdat['lu_sfr'].dims[0].label = 'log10(U) (dimensionless)'
    
            hfdat.create_dataset('lne_sfr',data=lne_sfr, maxshape=(None,None))
            hfdat['lne_sfr'].dims[0].label = 'log10(nH) (cm**-3)'
    
            hfdat.create_dataset('lz_sfr', data=loh12_sfr, maxshape=(None,None))
            hfdat['lz_sfr'].dims[0].label = 'log10(Z)'
            
            for i in range(len(const.lines_model[photmod_sfr])):           
                hfdat.create_dataset(const.lines_model[photmod_sfr][i] + '_sfr', 
                                     data=nebline_sfr[:,i], maxshape=(None,None))
                hfdat[const.lines_model[photmod_sfr][i] + '_sfr'].dims[0].label = 'Lines units: [Lsun = 3.826E+33egr s^-1 per unit SFR(Mo/yr) for 10^8yr]'
                
                if fluxes_sfr.any():
                    hfdat.create_dataset(const.lines_model[photmod_sfr][i] + '_sfr_flux', 
                                         data=fluxes_sfr[:,i], maxshape=(None,None))
                    hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux'].dims[0].label = 'Lines units: egr s^-1 cm^-2'
                    
                if fluxes_sfr_att.any():
                    hfdat.create_dataset(const.lines_model[photmod_sfr][i] + '_sfr_flux_att', 
                                         data=fluxes_sfr_att[:,i], maxshape=(None,None))
                    hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux_att'].dims[0].label = 'Lines units: egr s^-1 cm^-2'

                
                if nebline_sfr_att.any():
                    if nebline_sfr_att[0,i,0] > 0:
                        hfdat.create_dataset(const.lines_model[photmod_sfr][i] + '_sfr_att', 
                                             data=nebline_sfr_att[:,i], maxshape=(None,None))
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_att'].dims[0].label = 'Lines units: [Lsun = 3.826E+33egr s^-1 per unit SFR(Mo/yr) for 10^8yr]'
    
            if extra_param[0][0] != None:
                for i in range(len(extra_param)):
                    hfdat.create_dataset(extra_params_names[i], data=extra_param[i][:,None], maxshape=(None,None))
                    if extra_params_labels:
                        hfdat[extra_params_names[i]].dims[0].label = extra_params_labels[i]
    
    else:
        with h5py.File(outfile,'a') as hf:
            hfdat = hf['data']
            
            hfdat['lms'].resize((hfdat['lms'].shape[0] + lms.shape[0]),axis=0)
            hfdat['lms'][-lms.shape[0]:] = lms
            
            hfdat['lssfr'].resize((hfdat['lssfr'].shape[0] + lssfr.shape[0]),axis=0)
            hfdat['lssfr'][-lssfr.shape[0]:] = lssfr
            
            hfdat['lu_sfr'].resize((hfdat['lu_sfr'].shape[0] + lu_sfr.shape[0]),axis=0)
            hfdat['lu_sfr'][-lu_sfr.shape[0]:] = lu_sfr
            
            hfdat['lne_sfr'].resize((hfdat['lne_sfr'].shape[0] + lne_sfr.shape[0]),axis=0)
            hfdat['lne_sfr'][-lne_sfr.shape[0]:] = lne_sfr
            
            hfdat['lz_sfr'].resize((hfdat['lz_sfr'].shape[0] + loh12_sfr.shape[0]),axis=0)
            hfdat['lz_sfr'][-loh12_sfr.shape[0]:] = loh12_sfr
            
            
            for i in range(len(const.lines_model[photmod_sfr])): 
                hfdat[const.lines_model[photmod_sfr][i] + '_sfr'].resize((hfdat[const.lines_model[photmod_sfr][i] + '_sfr'].shape[1] + nebline_sfr.shape[2]),axis=1)
                hfdat[const.lines_model[photmod_sfr][i] + '_sfr'][:,-nebline_sfr.shape[2]:] = nebline_sfr[:,i]
                
                if fluxes_sfr.any():
                    hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux'].resize((hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux'].shape[1] + nebline_sfr.shape[2]),axis=1)
                    hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux'][:,-nebline_sfr.shape[2]:] = fluxes_sfr[:,i]
                
                if fluxes_sfr_att.any():
                     hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux_att'].resize((hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux_att'].shape[1] + nebline_sfr.shape[2]),axis=1)
                     hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux_att'][:,-nebline_sfr.shape[2]:] = fluxes_sfr_att[:,i]
                
                if nebline_sfr_att.any():
                    if nebline_sfr_att[0,i,0] > 0:
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_att'].resize((hfdat[const.lines_model[photmod_sfr][i] + '_sfr_att'].shape[1] + nebline_sfr_att.shape[2]),axis=1)
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_att'][:,-nebline_sfr_att.shape[2]:] = nebline_sfr_att[:,i]

            if extra_param[0][0] != None:
                for i in range(len(extra_param)):
                    hfdat[extra_params_names[i]].resize((hfdat[extra_params_names[i]].shape[0] + extra_param[i][:,None].shape[0]),axis=0)
                    hfdat[extra_params_names[i]][-extra_param[i][:,None].shape[0]:] = extra_param[i][:,None]
        

def write_data_AGN(lms,lssfr,lu_sfr,lne_sfr,loh12_sfr,lu_agn,lne_agn,loh12_agn,
               nebline_sfr,nebline_agn,nebline_sfr_att=None,nebline_agn_att=None,
               fluxes_sfr=None,fluxes_agn=None,fluxes_sfr_att=None,fluxes_agn_att=None,
               epsilon_sfr=None,epsilon_agn=None,
               extra_param=[[None]],extra_params_names=None,extra_params_labels=None,
               ew_notatt=None,ew_att=None,outfile='output.hdf5',attmod='ratios',
               unemod_sfr='kashino20',unemod_agn='panuzzo03',photmod_sfr='gutkin16',
               photmod_agn='feltre16',first=True):
    '''
    Create a .hdf5 file from a .dat file.

    Parameters
    ----------
    lms_sfr : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr_sfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lu_sfr : floats
     U of the galaxies per component.
    lne_sfr : floats
     ne of the galaxies per component (cm^-3).
    loh12_sfr : floats
     Metallicity of the galaxies per component (12+log(O/H))
    nebline_sfr : floats
      Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
    nebline_sfr_att : floats
      Array with the luminosity of the attenuated lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
    lms_agn : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr_agn : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lu_agn : floats
     U of the galaxies per component.
    lne_agn : floats
     ne of the galaxies per component (cm^-3).
    loh12_agn : floats
     Metallicity of the galaxies per component (12+log(O/H))
    nebline_agn : floats
      Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
    nebline_agn_att : floats
      Array with the luminosity of the attenuated lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)      
    extra_params : list
     Parameters from the input files which will be saved in the output file.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    extra_params_names : strings
     Names of the datasets in the output files for the extra parameters.
    extra_params_labels : strings
     Description labels of the datasets in the output files for the extra parameters.
    outfile : string
      Name of the output file.
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
    first : boolean
      If True it creates the HDF5 file (first subvolume). If false, it adds elements to the existing one.
    '''
    
    if first: 
        with h5py.File(outfile,'w') as hf:
            head = hf.create_dataset('header',(1,))
            head.attrs[u'HII model'] = unemod_sfr
            head.attrs[u'AGN model'] = unemod_agn
            head.attrs[u'Lines model for SF'] = photmod_sfr
            head.attrs[u'Lines model for AGN'] = photmod_agn
            head.attrs[u'Attenuation model'] = attmod
    
            # Data
            hfdat = hf.create_group('data')
            
            hfdat.create_dataset('lms', data=lms, maxshape=(None,None))
            hfdat['lms'].dims[0].label = 'log10(M*) (Msun)'
            
            hfdat.create_dataset('lssfr', data=lssfr, maxshape=(None,None))
            hfdat['lssfr'].dims[0].label = 'log10(SFR/M*) (1/yr)'
    
            hfdat.create_dataset('lu_sfr', data=lu_sfr, maxshape=(None,None))
            hfdat['lu_sfr'].dims[0].label = 'log10(U) (dimensionless)'
    
            hfdat.create_dataset('lne_sfr',data=lne_sfr, maxshape=(None,None))
            hfdat['lne_sfr'].dims[0].label = 'log10(nH) (cm**-3)'
    
            hfdat.create_dataset('lz_sfr', data=loh12_sfr, maxshape=(None,None))
            hfdat['lz_sfr'].dims[0].label = 'log10(Z)'
            
            hfdat.create_dataset('lu_agn', data=lu_agn, maxshape=(None,None))
            hfdat['lu_agn'].dims[0].label = 'log10(U) (dimensionless)'
    
            hfdat.create_dataset('lne_agn',data=lne_agn, maxshape=(None,None))
            hfdat['lne_agn'].dims[0].label = 'log10(nH) (cm**-3)'
    
            hfdat.create_dataset('lz_agn', data=loh12_agn, maxshape=(None,None))
            hfdat['lz_agn'].dims[0].label = 'log10(Z)'
            
            hfdat.create_dataset('epsilon_agn', data=epsilon_agn[None,:], maxshape=(None,None))
            hfdat['epsilon_agn'].dims[0].label = 'NLRs volume filling factor (dimensionless)'

            for i in range(len(const.lines_model[photmod_sfr])):           
                hfdat.create_dataset(const.lines_model[photmod_sfr][i] + '_sfr', 
                                     data=nebline_sfr[:,i], maxshape=(None,None))
                hfdat[const.lines_model[photmod_sfr][i] + '_sfr'].dims[0].label = 'Lines units: erg s^-1'
                
                if fluxes_sfr.any():
                    hfdat.create_dataset(const.lines_model[photmod_sfr][i] + '_sfr_flux', 
                                         data=fluxes_sfr[:,i], maxshape=(None,None))
                    hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux'].dims[0].label = 'Lines units: egr s^-1 cm^-2'
                    
                if fluxes_sfr_att.any():
                    if fluxes_sfr_att[0,i,0] >= 0:
                        hfdat.create_dataset(const.lines_model[photmod_sfr][i] + '_sfr_flux_att', 
                                             data=fluxes_sfr_att[:,i], maxshape=(None,None))
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux_att'].dims[0].label = 'Lines units: egr s^-1 cm^-2'
                
                if nebline_sfr_att.any():
                    if nebline_sfr_att[0,i,0] >= 0:
                        hfdat.create_dataset(const.lines_model[photmod_sfr][i] + '_sfr_att', 
                                             data=nebline_sfr_att[:,i], maxshape=(None,None))
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_att'].dims[0].label = 'Lines units: erg s^-1'

            for i in range(len(const.lines_model[photmod_agn])):
                hfdat.create_dataset(const.lines_model[photmod_agn][i] + '_agn', 
                                     data=nebline_agn[0,i][None,:], maxshape=(None,None))
                hfdat[const.lines_model[photmod_agn][i] + '_agn'].dims[0].label = 'Lines units: egr s^-1'
                
                if fluxes_agn.any():
                    hfdat.create_dataset(const.lines_model[photmod_agn][i] + '_agn_flux', 
                                         data=fluxes_agn[0,i][None,:], maxshape=(None,None))
                    hfdat[const.lines_model[photmod_agn][i] + '_agn_flux'].dims[0].label = 'Lines units: egr s^-1 cm^-2'
                    
                if fluxes_agn_att.any():
                    if fluxes_agn_att[0,i,0] >= 0:
                        hfdat.create_dataset(const.lines_model[photmod_agn][i] + '_agn_flux_att', 
                                             data=fluxes_agn_att[0,i][None,:], maxshape=(None,None))
                        hfdat[const.lines_model[photmod_agn][i] + '_agn_flux_att'].dims[0].label = 'Lines units: egr s^-1 cm^-2'
                
                if nebline_agn_att.any():
                    if nebline_agn_att[0,i,0] >= 0:
                        hfdat.create_dataset(const.lines_model[photmod_agn][i] + '_agn_att', 
                                             data=nebline_agn_att[0,i][None,:], maxshape=(None,None))
                        hfdat[const.lines_model[photmod_agn][i] + '_agn_att'].dims[0].label = 'Lines units: egr s^-1'

            if extra_param[0][0] != None:
                for i in range(len(extra_param)):
                    hfdat.create_dataset(extra_params_names[i], data=extra_param[i][None,:], maxshape=(None,None))
                    if extra_params_labels:
                        hfdat[extra_params_names[i]].dims[0].label = extra_params_labels[i]
    else:
        with h5py.File(outfile,'a') as hf:
            hfdat = hf['data']
            
            hfdat['lms'].resize((hfdat['lms'].shape[0] + lms.shape[0]),axis=0)
            hfdat['lms'][-lms.shape[0]:] = lms
            
            hfdat['lssfr'].resize((hfdat['lssfr'].shape[0] + lssfr.shape[0]),axis=0)
            hfdat['lssfr'][-lssfr.shape[0]:] = lssfr
            
            hfdat['lu_sfr'].resize((hfdat['lu_sfr'].shape[0] + lu_sfr.shape[0]),axis=0)
            hfdat['lu_sfr'][-lu_sfr.shape[0]:] = lu_sfr
            
            hfdat['lne_sfr'].resize((hfdat['lne_sfr'].shape[0] + lne_sfr.shape[0]),axis=0)
            hfdat['lne_sfr'][-lne_sfr.shape[0]:] = lne_sfr
            
            hfdat['lz_sfr'].resize((hfdat['lz_sfr'].shape[0] + loh12_sfr.shape[0]),axis=0)
            hfdat['lz_sfr'][-loh12_sfr.shape[0]:] = loh12_sfr
            
            hfdat['lu_agn'].resize((hfdat['lu_agn'].shape[0] + lu_agn.shape[0]),axis=0)
            hfdat['lu_agn'][-lu_agn.shape[0]:] = lu_agn
            
            hfdat['lne_agn'].resize((hfdat['lne_agn'].shape[0] + lne_agn.shape[0]),axis=0)
            hfdat['lne_agn'][-lne_agn.shape[0]:] = lne_agn
            
            hfdat['lz_agn'].resize((hfdat['lz_agn'].shape[0] + loh12_agn.shape[0]),axis=0)
            hfdat['lz_agn'][-loh12_agn.shape[0]:] = loh12_agn
            
            hfdat['epsilon_agn'].resize((hfdat['epsilon_agn'].shape[1] + epsilon_agn[None,:].shape[1]),axis=1)
            hfdat['epsilon_agn'][0,-epsilon_agn[None,:].shape[1]:] = epsilon_agn[None,:]
            
            for i in range(len(const.lines_model[photmod_sfr])): 
                hfdat[const.lines_model[photmod_sfr][i] + '_sfr'].resize((hfdat[const.lines_model[photmod_sfr][i] + '_sfr'].shape[1] + nebline_sfr.shape[2]),axis=1)
                hfdat[const.lines_model[photmod_sfr][i] + '_sfr'][:,-nebline_sfr.shape[2]:] = nebline_sfr[:,i]
                
                if fluxes_sfr.any():
                    if fluxes_sfr[0,i,0] >= 0:
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux'].resize((hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux'].shape[1] + nebline_sfr.shape[2]),axis=1)
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux'][:,-nebline_sfr.shape[2]:] = fluxes_sfr[:,i]
                
                if fluxes_sfr_att.any():
                    if fluxes_sfr_att[0,i,0] >= 0:
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux_att'].resize((hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux_att'].shape[1] + nebline_sfr.shape[2]),axis=1)
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_flux_att'][:,-nebline_sfr.shape[2]:] = fluxes_sfr_att[:,i]
                
                if nebline_sfr_att.any():
                    if nebline_sfr_att[0,i,0] >= 0:
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_att'].resize((hfdat[const.lines_model[photmod_sfr][i] + '_sfr_att'].shape[1] + nebline_sfr_att.shape[2]),axis=1)
                        hfdat[const.lines_model[photmod_sfr][i] + '_sfr_att'][:,-nebline_sfr_att.shape[2]:] = nebline_sfr_att[:,i]
                        
            for i in range(len(const.lines_model[photmod_agn])):
                hfdat[const.lines_model[photmod_agn][i] + '_agn'].resize((hfdat[const.lines_model[photmod_agn][i] + '_agn'].shape[1] + nebline_agn.shape[2]),axis=1)
                hfdat[const.lines_model[photmod_agn][i] + '_agn'][:,-nebline_agn.shape[2]:] = nebline_agn[0,i][None,:]
                
                if fluxes_agn.any():
                    hfdat[const.lines_model[photmod_agn][i] + '_agn_flux'].resize((hfdat[const.lines_model[photmod_agn][i] + '_agn_flux'].shape[1] + nebline_agn.shape[2]),axis=1)
                    hfdat[const.lines_model[photmod_agn][i] + '_agn_flux'][:,-nebline_agn.shape[2]:] = fluxes_agn[0,i][None,:]
                
                if fluxes_agn_att.any():
                    if fluxes_agn_att[0,i,0] >= 0:
                        hfdat[const.lines_model[photmod_agn][i] + '_agn_flux_att'].resize((hfdat[const.lines_model[photmod_agn][i] + '_agn_flux_att'].shape[1] + nebline_agn.shape[2]),axis=1)
                        hfdat[const.lines_model[photmod_agn][i] + '_agn_flux_att'][:,-nebline_agn.shape[2]:] = fluxes_agn_att[0,i][None,:]
                
                if nebline_agn_att.any():
                    if nebline_agn_att[0,i,0] >= 0:
                        hfdat[const.lines_model[photmod_agn][i] + '_agn_att'].resize((hfdat[const.lines_model[photmod_agn][i] + '_agn_att'].shape[1] + nebline_agn_att.shape[2]),axis=1)
                        hfdat[const.lines_model[photmod_agn][i] + '_agn_att'][:,-nebline_agn_att.shape[2]:] = nebline_agn_att[0,i][None,:]
            
            if extra_param[0][0] != None:
                for i in range(len(extra_param)):
                    hfdat[extra_params_names[i]].resize((hfdat[extra_params_names[i]].shape[1] + extra_param[i][None,:].shape[1]),axis=1)
                    hfdat[extra_params_names[i]][0,-extra_param[i][None,:].shape[1]:] = extra_param[i][None,:]
        
