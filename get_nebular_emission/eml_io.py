"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
"""
import time
import h5py
import sys
import os
import numpy as np
import get_nebular_emission.eml_const as const
import math
from pathlib import Path

homedir = Path.home()

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
            print('WARNING (eml_io.check_file): file not found {}'.
                  format(infile))

    return file_fine


# def create_dir(outdir):
#     '''
#     Return True if directory has been created
#     '''
    
#     if not os.path.exists(outdir):
#         try:
#             os.makedirs(outdir)
#         except:
#             print('WARNING (eml_io.create_dir): problem creating directory ',
#                   outdir)
#             return False
#     return True


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
        print('STOP (eml_io.get_ncomponents): ',
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
              inputformat='HDF5',Plotting=False, Testing=False, verbose=True):
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
    Plotting : boolean
      If True run verification plots with all data.
    Testing : boolean
      If True only run over few entries for testing purposes

    Returns
    -------
    lms, lssfr, loh12 : floats
    cut : integers
    '''
    
    if verbose:
        print(infile)
    
    check_file(infile, verbose=verbose)

    ncomp = get_ncomponents(cols)
    
    if (Testing and not Plotting):
        limit = 50
    else:
        limit = None
        
    if inputformat not in const.inputformats:
        if verbose:
            print('STOP (eml_io): Unrecognised input format.',
                  'Possible input formats = {}'.format(const.inputformats))
        sys.exit()
    elif inputformat=='HDF5':
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
                    lms = [hf[cols[i][0]][:limit]]
                    lssfr = [hf[cols[i][1]][:limit]]
                    loh12 = [hf[cols[i][2]][:limit]]
                else:
                    lms = np.append(lms,[hf[cols[i][0]][:limit]],axis=0)
                    lssfr = np.append(lssfr,[hf[cols[i][1]][:limit]],axis=0)
                    loh12 = np.append(loh12,[hf[cols[i][2]][:limit]],axis=0)
    elif inputformat=='textfile':
        ih = get_nheader(infile)
        
        cut = np.arange(len(np.loadtxt(infile,usecols=cols[0],skiprows=ih)[:limit]))
        
        if cutcols:
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
                lms = [X[0]]
                lssfr = [X[1]]
                loh12 = [X[2]]
            else:
                lms = np.append(lms,[X[0]],axis=0)
                lssfr = np.append(lssfr,[X[1]],axis=0)
                loh12 = np.append(loh12,[X[2]],axis=0)
    else:
        if verbose:
            print('STOP (eml_io.read_data): ',
                  'Input file has not been found.')
        sys.exit()
        
    lms = lms.T
    lssfr = lssfr.T
    loh12 = loh12.T
            
    return lms[cut], lssfr[cut], loh12[cut], cut


def get_data(infile, cols, h0=None, inputformat='HDF5', 
             IMF_i=['Chabrier', 'Chabrier'], IMF_f=['Kroupa', 'Kroupa'], 
             cutcols=None, mincuts=[None], maxcuts=[None],
             attmod='GALFORM', LC2sfr=False, mtot2mdisk=True, 
             verbose=False, Plotting=False, Testing=False):
    '''
    Get Mstars, sSFR and (12+log(O/H)) in the adecuate units.

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
    h0 : float
      If not None: value of h, H0=100h km/s/Mpc.
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

    Returns
    -------
    lms, lssfr, loh12 : floats
    cut : integers
    '''
    
    lms,lssfr,loh12,cut = read_data(infile, cols=cols, cutcols=cutcols,
                                maxcuts=maxcuts, mincuts=mincuts, inputformat=inputformat, 
                                Plotting=Plotting, Testing=Testing, verbose=verbose)

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
                print('STOP (eml_io.get_data): ',
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
        ind = np.where(lssfr > 0.)
        for comp in range(ncomp):
            lssfr[ind,comp] = lssfr[ind,comp]*(const.IMF_SFR[IMF_i[comp]]/const.IMF_M[IMF_i[comp]])*(const.IMF_M[IMF_f[comp]]/const.IMF_SFR[IMF_f[comp]])
        lssfr[ind] = np.log10(lssfr[ind]) - lms[ind] - 9. 

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



    if Plotting: # here : Search more efficient form. Allow more components in the header
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


def get_reducedfile(infile, outfile, indcol, verbose=False):
    '''
    Get reduced file with only the wanted columns data

    Parameters
    ----------
    infile : string
     - Name of the input file.
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    outfile : string
     - Name of the output file.
     - Text file (*.dat, *txt, *.cat), columns separated by ' '.
    indcol : list
     - [mstars_total, mstars_bulge, mstardot, mstardot_burst, mag_LC_r_disk, mag_LC_r_bulge, zcold, zcold_burst]
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    verbose : boolean
      If True print out messages
    '''

    check_file(infile,verbose=verbose)

    if ('.hdf5' in infile):
        print('Program not set to deal with hdf5 files yet, input a text file')
    else:
        with open(infile,"r") as ff:
            ih = get_nheader(infile)
            with open(outfile, "w") as outf:
                for il,line in enumerate(ff):
                    if il<ih-1:
                        n=line
                        #print(n,il,ih)
                        #exit()
                        outf.write(n)

                    if il==ih-1:
                        headlast=line.split()
                        headlast=headlast[1:]
                        headlast=np.array(headlast)

                        # Here : Better with a matrix form. Future change.
                        # This is necessary to only select the wanted columns of the header too:
                        headlast=np.column_stack(('#',headlast[indcol[0]], headlast[indcol[1]],
                                                  headlast[indcol[2]], headlast[indcol[3]],
                                                  headlast[indcol[4]],headlast[indcol[5]],
                                                  headlast[indcol[6]], headlast[indcol[7]]))
                        #print(headlast)
                        np.savetxt(outf,headlast,fmt='%s')

                    if il>=ih: # Here : Remove the last line. Type string : Number of galaxies in the sample.
                        if ('.csv' in infile):
                            allcols = line.split(',')
                        else:
                            allcols = line.split()


                        mstars_total = np.array(float(allcols[indcol[0]]))
                        mstars_bulge = np.array(float(allcols[indcol[1]]))
                        mstardot = np.array(float(allcols[indcol[2]]))
                        mstardot_burst = np.array(float(allcols[indcol[3]]))
                        mag_LC_r_disk = np.array(float(allcols[indcol[4]]))
                        mag_LC_r_bulge = np.array(float(allcols[indcol[5]]))
                        zcold = np.array(float(allcols[indcol[6]]))
                        zcold_burst = np.array(float(allcols[indcol[7]]))

                        print(mstars_total)


                        datatofile = np.column_stack((mstars_total, mstars_bulge, mstardot,
                                                    mstardot_burst, mag_LC_r_disk, mag_LC_r_bulge,
                                                    zcold, zcold_burst))

                        np.savetxt(outf,datatofile)


            outf.closed
        ff.closed

    #outfile = print(outf, 'Output file :{}'.format(outfile))

    return

def hdf5_from_text(infile, outfile, galform=True, verbose=True):
    '''
    Create a .hdf5 file from a text file.

    Parameters
    ----------
    infile : string
     - Name of the input file.
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    outfile : string
     Name of the output file.
    galform : boolean
     If True, it doesn't read the headers and instead has them hard-coded.
    verbose : boolean
      If True print out messages
    '''

    check_file(infile, verbose=verbose)

    ih = get_nheader(infile)
    
    X = np.loadtxt(infile,skiprows=ih).T
    
    if galform:
        headers = const.gal_headers
    else:
        with open(infile, "r") as ff:
            for il, line in enumerate(ff):
                # Read until the data starts
                if (il<ih-1): continue
            
                # Read data
                if ('.csv' in infile):
                    allcols = line.split(',')
                else:
                    allcols = line.split()
                
                if (il==ih-1):
                    headers = np.copy(allcols[1:])
                    headers = headers.astype('<U30')
                    break
        
    with h5py.File(outfile, 'w') as f:
        f.create_dataset('header',(1,))
        hf = f.create_group('data')
        for i in range(len(headers)):
            hf.create_dataset(headers[i],data=X[i])
            
    return X

def write_data(lms,lssfr,lu,lne,loh12,nebline,nebline_att,outfile,attmod='ratios',
               unemod='kashino20',photmod='gutkin16',first=True):
    '''
    Create a .hdf5 file from a .dat file.

    Parameters
    ----------
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lu : floats
     U of the galaxies per component.
    lne : floats
     ne of the galaxies per component (cm^-3).
    loh12 : floats
     Metallicity of the galaxies per component (12+log(O/H))
    nebline : floats
      Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
    nebline_att : floats
      Array with the luminosity of the attenuated lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
    outfile : string
      Name of the output file.
    unemod : string
      Model to go from galaxy properties to U and ne
    photmod : string
      Photoionisation model to be used for look up tables.
    verbose : boolean
      If True print out messages
    '''
    
    if first: 
        with h5py.File(outfile,'w') as hf:
            head = hf.create_dataset('header',(1,))
            head.attrs[u'HII model'] = unemod
            head.attrs[u'Lines model'] = photmod
            head.attrs[u'Attenuation model'] = attmod
    
            # Data
            hfdat = hf.create_group('data')
            
            hfdat.create_dataset('lms', data=lms, maxshape=(None,None))
            hfdat['lms'].dims[0].label = 'log10(M*) (Msun)'
            
            hfdat.create_dataset('lssfr', data=lssfr, maxshape=(None,None))
            hfdat['lssfr'].dims[0].label = 'log10(SFR/M*) (1/yr)'
    
            hfdat.create_dataset('lu', data=lu, maxshape=(None,None))
            hfdat['lu'].dims[0].label = 'log10(U) (dimensionless)'
    
            hfdat.create_dataset('lne',data=lne, maxshape=(None,None))
            hfdat['lne'].dims[0].label = 'log10(nH) (cm**-3)'
    
            hfdat.create_dataset('lz', data=loh12, maxshape=(None,None))
            hfdat['lz'].dims[0].label = 'log10(Z)'
            
            for i in range(len(const.lines_model[photmod])):           
                hfdat.create_dataset(const.lines_model[photmod][i], 
                                     data=nebline[:,i], maxshape=(None,None))
                hfdat[const.lines_model[photmod][i]].dims[0].label = 'Lines units: [Lsun = 3.826E+33egr s^-1 per unit SFR(Mo/yr) for 10^8yr]'
                
                if nebline_att[0,i,0] >= 0:
                    hfdat.create_dataset(const.lines_model[photmod][i] + '_att', 
                                         data=nebline_att[:,i], maxshape=(None,None))
                    hfdat[const.lines_model[photmod][i] + '_att'].dims[0].label = 'Lines units: [Lsun = 3.826E+33egr s^-1 per unit SFR(Mo/yr) for 10^8yr]'
    else:
        with h5py.File(outfile,'a') as hf:
            hfdat = hf['data']
            
            hfdat['lu'].resize((hfdat['lu'].shape[0] + lu.shape[0]),axis=0)
            hfdat['lu'][-lu.shape[0]:] = lu
            
            hfdat['lne'].resize((hfdat['lne'].shape[0] + lne.shape[0]),axis=0)
            hfdat['lne'][-lne.shape[0]:] = lne
            
            hfdat['lms'].resize((hfdat['lms'].shape[0] + lms.shape[0]),axis=0)
            hfdat['lms'][-lms.shape[0]:] = lms
            
            hfdat['lssfr'].resize((hfdat['lssfr'].shape[0] + lssfr.shape[0]),axis=0)
            hfdat['lssfr'][-lssfr.shape[0]:] = lssfr
            
            hfdat['lz'].resize((hfdat['lz'].shape[0] + loh12.shape[0]),axis=0)
            hfdat['lz'][-loh12.shape[0]:] = loh12
            
            for i in range(len(const.lines_model[photmod])): 
                hfdat[const.lines_model[photmod][i]].resize((hfdat[const.lines_model[photmod][i]].shape[1] + nebline.shape[2]),axis=1)
                hfdat[const.lines_model[photmod][i]][:,-nebline.shape[2]:] = nebline[:,i]
                
                if nebline_att[0,i,0] >= 0:
                    hfdat[const.lines_model[photmod][i] + '_att'].resize((hfdat[const.lines_model[photmod][i] + '_att'].shape[1] + nebline_att.shape[2]),axis=1)
                    hfdat[const.lines_model[photmod][i] + '_att'][:,-nebline_att.shape[2]:] = nebline_att[:,i]
        