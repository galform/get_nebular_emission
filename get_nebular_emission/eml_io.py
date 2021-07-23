"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
"""
import sys
import os
import numpy as np
import get_nebular_emission.eml_const as const

def stop_if_no_file(infile):
    '''
    Stop if the file does not exist
    '''
    if (not os.path.isfile(infile)):
        print('STOP: no input file {}'.format(infile)) 
        sys.exit()
    return


def check_file(infile,verbose=False):
    '''
    Return True if the file exists
    '''
    file_fine = True  
    if (not os.path.isfile(infile)):
        file_fine = False
        if verbose:
            print('WARNING (eml_io.check_file): file not found {}'.
                  format(infile))

    return file_fine


def create_dir(outdir):
    '''
    Return True if directory has been created
    '''
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except:
            print('WARNING (eml_io.create_dir): problem creating directory ',
                  outdir)
            return False
    return True


def get_nheader(infile):
    '''
    Given a file with a structure: header+data, 
    counts the number of header lines

    Args:
    infile: string, input file

    Returns:
    ih: integer, number of lines with the header text
    '''

    ih = 0
    ff = open(infile,'r')
    for line in ff:
        if not line.strip():
            # Count any empty lines in the header
            ih += 1
        else:
            sline = line.strip()
            
            # Check that the first character is not a digit
            char1 = sline[0]
            word1 = sline.split()[0]
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
    return ih
        

def get_zfile(zmet,photmod='gutkin'):
    '''
    Given a metallicity get the name of the corresponding table

    Parameters
    ----------
    zmet : float
       Metallicity value
    photomod : string
       Name of the considered photoionisation model

    Returns
    -------
    zfile : string
       Name of the model file with data for the given metallicity
    '''
    
    dec = str(zmet).split('.')[-1]
    root = 'nebular_data/'+photmod+'_tables/nebular_emission_Z'
    zfile = root+dec+'.txt'

    file_fine = check_file(zfile)
    if (not file_fine):
        zfile = None
    
    return zfile
    

def get_ncomponents(cols):
    '''
    Get the number of components to estimate the emission lines from

    Parameters:
    ----------
    cols : list
      List of columns with M*, SFR, Z, per components

    Return:
    ncomp : int
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

def get_data(infile, cols, h0=None, inoh=False, verbose=False, Testing=False):
    '''
    Get Mstars, sSFR and (12+log(O/H)) in the adecuate units

    Parameters
    ----------
    infile : string
      Name of the input file. 
      In text files (*.dat, *txt, *.cat), columns separated by ' '.
      In csv files (*.csv), columns separated by ','.
    cols : list
      [[component1_stellar_mass,sfr,Z],[component2_stellar_mass,sfr,Z],...]
      For text or csv files: list of integers with column position.
      For hdf5 files: list of data names.
    h0 : float
      If not None: value of h, H0=100h km/s/Mpc.
    inoh : boolean
      If yes, the metallicity has already been provided as 12+log(O/H)
    verbose : boolean
      Yes = print out messages
    Testing : boolean
      Yes = to only run over few entries for testing purposes

    Returns
    -------
    mstars, ssfr, oh : floats
    '''
    
    check_file(infile,verbose=verbose)

    ncomp = get_ncomponents(cols)

    if ('.hdf5' in infile):
        print('Program not set to deal with hdf5 files yet, input a text file')
    else:
        # Text file
        ih = get_nheader(infile)
        # Jump the header and read the provided columns
        with open(infile, "r") as ff:
            for il, line in enumerate(ff):
                # Read until the data starts
                if (il<ih): continue

                # Read data
                if ('.csv' in infile):
                    allcols = line.split(',')
                else:
                    allcols = line.split()

                if (il == ih):
                    lms   = np.array([[float(allcols[cols[ic][0]]) for ic in range(ncomp)]])
                    lssfr = np.array([[float(allcols[cols[ic][1]]) for ic in range(ncomp)]])
                    loh12 = np.array([[float(allcols[cols[ic][2]]) for ic in range(ncomp)]])
                else:
                    comp = np.array([[float(allcols[cols[ic][0]]) for ic in range(ncomp)]])
                    lms = np.append(lms,comp,axis=0)

                    comp = np.array([[float(allcols[cols[ic][1]]) for ic in range(ncomp)]])
                    lssfr= np.append(lssfr,comp,axis=0)

                    comp = np.array([[float(allcols[cols[ic][2]]) for ic in range(ncomp)]])
                    loh12= np.append(loh12,comp,axis=0)

                if (Testing and il>ih+50): break

        # Set to a default value if negative stellar masses
        ind = np.where(lms<=0.)
        lms[ind] = const.notnum
        lssfr[ind] = const.notnum
        loh12[ind] = const.notnum

        ind = np.where(lssfr<=0) # Avoid other negative SFR
        lssfr[ind] = const.notnum ; loh12[ind] = const.notnum

        ind = np.where(loh12<=0) # Avoid other negative Z
        loh12[ind] = const.notnum

        # Take the log of the stellar mass
        ind = np.where(lms>0.)
        lms[ind] = np.log10(lms[ind]) 

        # Obtain log10(sSFR) in 1/yr
        ind = np.where(lssfr>0.)
        lssfr[ind] = np.log10(lssfr[ind]) - 9. - lms[ind]
        #here: allow for input number of ionizing photons
        
        if h0:
            # Correct the units of the stellar mass
            lms = lms - np.log10(h0)

        # Obtain 12+log10(O/H) from Z=MZcold/Mcold
        ind = np.where(loh12>0)
        loh12[ind] = const.ohsun + np.log10(loh12[ind]) - np.log10(const.zsun)
        #here: allow for loh12 direct input
                    
    return lms,lssfr,loh12
