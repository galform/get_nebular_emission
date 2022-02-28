"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
"""
import sys
import os
import numpy as np
import get_nebular_emission.eml_const as const
import math

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
        print('STOP (eml_io.get_ncomponents): ',
              'Columns should be given as m_sfr_z=[[0,1,2]]')
        sys.exit()
        
    return ncomp

def get_data(infile, cols, h0=None, inoh=False, LC2sfr=False, mtot2mdisk=True, verbose=False, Plotting=False, Testing=False):
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
      Expected : component1 = total or disk, component2 = bulge
      For text or csv files: list of integers with column position.
      For hdf5 files: list of data names.
    h0 : float
      If not None: value of h, H0=100h km/s/Mpc.
    inoh : boolean
      If yes, the metallicity has already been provided as 12+log(O/H)
    LC2sfr : boolean
      If True magnitude of Lyman Continuum photons expected as input for SFR.
    mtot2mdisk : boolean
      Yes = transform the total mass into the disk mass. disk mass = total mass - bulge mass.
    verbose : boolean
      Yes = print out messages
    Plotting : boolean
      If True run verification plots with all data.
    Testing : boolean
      Yes = to only run over few entries for testing purposes

    Returns
    -------
    mstars, ssfr, oh : floats
    '''

    check_file(infile, verbose=verbose)

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

                if (Testing and not Plotting and il>ih+50): break


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
            lms_tot = lms[:,0] + lms[:,1]

            # Take the log of the total stellar mass:
            ind = np.where(lms_tot > 0.)
            lms_tot[ind] = np.log10(lms_tot[ind])

            # Take the log of the stellar mass:
            ind = np.where(lms > 0.)
            lms[ind] = np.log10(lms[ind])


        # Obtain log10(sSFR) in 1/yr and calculate SFR from LC photons if necessary

        lssfrd = np.zeros(len(lssfr))
        lssfrb = np.zeros(len(lssfr))
        ssfrd = np.zeros(len(lssfr))
        ssfrb = np.zeros(len(lssfr))
        lssfr_tot = np.zeros(len(lssfr))

        if LC2sfr:
            # DISK:
            ins = np.zeros(len(lssfr))
            ind = np.where(lssfr[:, 0] != const.notnum)
            ins[ind] = 1.02*(10.**(-0.4*lssfr[ind, 0]-4.))
            ind = np.where(ins > 0)
            lssfrd[ind] = np.log10(ins[ind]) - 9. - lms[ind, 0]
            ind = np.where(ins < 0)
            lssfrd[ind] = const.notnum

            ind = np.where(lssfr[:, 0] == const.notnum)
            lssfrd[ind] = const.notnum

            # Calculate sSFR disk
            ind = np.where(lssfrd != const.notnum)
            ssfrd[ind] = 10. ** lssfrd[ind]

            # BULGE:
            ind = np.where(lssfr[:, 1] != const.notnum)
            ins[ind] = 0.360*(10.**(-0.4*lssfr[ind, 1]-4.))
            ind = np.where(ins > 0)
            lssfrb[ind] = np.log10(ins[ind]) - 9. - lms[ind, 1]
            ind = np.where(ins < 0)
            lssfrb[ind] = const.notnum

            ind = np.where(lssfr[:,1] == const.notnum)
            lssfrb[ind] = const.notnum

            # Calculate sSFR bulge
            ind = np.where(lssfrb != const.notnum)
            ssfrb[ind] = 10.**lssfrb[ind]

            # Final Vector and Total
            lssfr = np.column_stack((lssfrd,lssfrb))
            ins = ssfrd + ssfrb
            ind = np.where(ins > 0)
            lssfr_tot[ind] = np.log10(ins[ind])

        # If instantaneous SFR as input:
        else:
            # Take the log of the ssfr:
            ind = np.where(lssfr > 0.)
            lssfr[ind] = np.log10(lssfr[ind]) - 9. - lms[ind]

            # Total
            ind = np.where(lssfr[:,0]!=const.notnum)
            ssfrd[ind] = 10.**(lssfr[ind,0])

            ind = np.where(lssfr[:,1]!=const.notnum)
            ssfrb[ind] = 10.**(lssfr[ind,1])

            ins = ssfrd + ssfrb
            ind = np.where(ins>0)
            lssfr_tot[ind] = np.log10(ins[ind])

        if h0:
            # Correct the units of the stellar mass
            lms = lms - np.log10(h0)
            lms_tot = lms_tot - np.log10(h0)

        if Plotting:
            lsfr = lssfr_tot+lms_tot


        # Obtain 12+log10(O/H) from Z=MZcold/Mcold
        ind = np.where(loh12>0)
        loh12[ind] = np.log10(loh12[ind]) + const.ohsun - np.log10(const.zsun)

        oh12d = np.zeros(len(loh12))
        oh12b = np.zeros(len(loh12))
        loh12_tot = np.zeros(len(loh12))

        # Total :
        ind = np.where(loh12[:,0] != const.notnum)
        oh12d[ind] = 10. ** (loh12[ind, 0])

        ind = np.where(loh12[:, 1] != const.notnum)
        oh12b[ind] = 10. ** (loh12[ind, 1])

        ins = oh12d+oh12b
        ind = np.where(ins>0)
        loh12_tot[ind] = np.log10(ins[ind])

        # Here: allow for loh12 direct input



        if Testing and Plotting: # here : Search more efficient form. Allow more components in the header
            header1 = 'log(mstars_tot), log(mstars_disk), log(mstars_bulge),' \
                      ' log(SFR_tot), log(sSFR_tot), log(sSFR_disk), log(sSFR_bulge) ' \
                      '(12 + log(O/H))_tot, (12 + log (O/H))_disk, (12 + log (O/H))_bulge'
            datatofile=np.column_stack((lms_tot,lms,lsfr,lssfr_tot,lssfr,loh12_tot,loh12))

            if LC2sfr:
                outfil = r"example_data/tmp_LC.dat"
            else:
                outfil = r"example_data/tmp_avSFR.dat"

            with open(outfil, 'w') as outf:
                np.savetxt(outf, datatofile, delimiter=' ', header=header1)
                outf.closed


                    
    return lms,lssfr,loh12


def get_reducedfile(infile, outfile, indcol, verbose=False):
    '''
    Get reduced file with only the wanted columns data.

    Parameters
    ----------
    infile : string
      Name of the input file.
      In text files (*.dat, *txt, *.cat), columns separated by ' '.
      In csv files (*.csv), columns separated by ','.
    outfile : string
      Name of the output file.
      Text file (*.dat, *txt, *.cat), columns separated by ' '.
    indcol : list
      [mstars_total, mstars_bulge, mstardot, mstardot_burst, mag_LC_r_disk, mag_LC_r_bulge, zcold, zcold_burst]
      For text or csv files: list of integers with column position.
      For hdf5 files: list of data names.
    verbose : boolean
      Yes = print out messages
    Returns
    -------
    outfile : string
    Path to the reduced output file.
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