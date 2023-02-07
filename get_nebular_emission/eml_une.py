import sys

import numpy as np
import get_nebular_emission.eml_const as const

def get_une_withmodel(lms, lssfr, loh12, unemod, gamma=None, verbose=False):
    '''
    Given log10(Mstar), log10(sSFR) and 12+log(O/H),
    get the ionizing parameter, logU, and the electron density, logne,
    using the specified model.

    Parameters
    ----------
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    loh12 : floats
     Metallicity of the galaxies per component (12+log(O/H)).
    unemod : string
      Model to go from galaxy properties to U and ne
    verbose : boolean
      Yes = print out messages

    Returns
    -------
    lu, lne, loh12 : floats
    '''
    
    lu, lne = [np.full(np.shape(lms), const.notnum) for i in range(2)]

    ind = np.where((lssfr > const.notnum) &
                   (lms > const.notnum) &
                   (loh12 > const.notnum))
    
    if unemod == 'kashino20': #Coefficients in Table 2 from Kashino and Inoue 2019 (https://arxiv.org/pdf/1812.06939.pdf)
        if (np.shape(ind)[1]>1):
            loh12[ind] = loh12[ind] + const.ohsun - np.log10(const.zsun)
            
            lne[ind] = 2.066 + 0.310*(lms[ind]-10) + 0.492*(lssfr[ind] + 9.)
            lu[ind] =  -2.316 - 0.360*(loh12[ind] -8.) -0.292*lne[ind] + 0.428*(lssfr[ind] + 9.)
            # lu[ind] =  -3.073 - 0.137*(lms[ind]-10) + 0.372*(lssfr[ind] + 9.)
            
            loh12[ind] = loh12[ind] - const.ohsun + np.log10(const.zsun) # We leave it in log(Z)
    
    if unemod == 'orsi14':
        if not gamma:
            print('STOP (eml_une.get_une_withmodel): ',
                  'Gamma for Orsi14 model not specified.')
            sys.exit()
        if (np.shape(ind)[1]>1):
            lne[ind] = 2.066 + 0.310*(lms[ind]-10) + 0.492*(lssfr[ind] + 9.) #np.log10(10)
            lu[ind] = np.log10(const.q0*((10**loh12[ind])/const.Z0)**-gamma / const.c)
    

    return lu, lne, loh12
    
def get_une(lms, lssfr, loh12, unemod='kashino20', gamma=None, LC2sfr=False, Testing=False, Plotting=False, verbose=True):
    '''
    Given log10(Mstar), log10(sSFR) and 12+log(O/H),
    get the ionizing parameter, U, and the electron density, ne.

    Parameters
    ----------
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    loh12 : floats
     Metallicity of the galaxies per component (12+log(O/H)).
    unemod : string
      Model to go from galaxy properties to U and ne
    LC2sfr : boolean
      True = Do the change from Lyman Continuum photons to instantaneus SFR
    verbose : boolean
      Yes = print out messages
    Plotting : boolean
      If True run verification plots with all data.
    Testing : boolean
      Yes = to only run over few entries for testing purposes

    Returns
    -------
    lu, lne, loh12 : floats
    '''

    ncomp = len(lms[0])
    
    if unemod not in const.unemods:
        if verbose:
            print('STOP (eml_une): Unrecognised model to get U and ne.')
            print('                Possible unemod= {}'.format(const.unemods))
        sys.exit()
    elif (unemod == 'kashino20'):
        lu, lne, loh12 = get_une_withmodel(lms,lssfr,loh12,unemod,verbose=verbose)
    elif (unemod == 'orsi14'):
        lu, lne, loh12 = get_une_withmodel(lms,lssfr,loh12,unemod,gamma,verbose=verbose)

    if Plotting:
        if ncomp==2:
            header1 = 'log(u_disk),log(u_bulge),log(ne_disk),log(ne_bulge)'
            datatofile = np.append(lu, lne, axis=1)
        elif ncomp==1:
            header1 = 'log(u), log(ne)'
            datatofile = np.append(lu, lne, axis=1)
            
        if LC2sfr:
            outfil = r"example_data/tmp_une_LC_SAGE.dat"
        else:
            outfil = r"example_data/tmp_une_avSFR_SAGE.dat"
        with open(outfil, 'w') as outf:
            np.savetxt(outf, datatofile, delimiter=' ', header=header1)
        
    return lu, lne, loh12
