"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
"""
import sys
import os
import numpy as np

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

    return ncomp

def get_data(infile, cols, h0=None, inoh=False, verbose=False):
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

                if (il == ih):#here start arrays
                    
                print(il,allcols[2]); exit()
        #mass, sfr, 
        
    return ih,ih,ih
