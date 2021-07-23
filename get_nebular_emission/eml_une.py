import sys
import numpy as np
import get_nebular_emission.eml_const as const

def get_une_kashino20(lms, lssfr, loh12, verbose=False):
    '''
    Given log10(Mstar), log10(sSFR) and 12+log(O/H),
    get the ionizing parameter, logU, and the electron density, logne,
    using equations 10 and 12 
    from Kashino and Inoue 2019 (https://arxiv.org/pdf/1812.06939.pdf)

    Parameters
    ----------
    lms : float
      log10(Mstar/Msun)
    lssfr : float
      log10(sSFR/yr), it should be an instantaneous measurement
    loh12 : float 
      12+log(O/H)
    unemod : string
      Model to go from galaxy properties to U and ne
    verbose : boolean
      Yes = print out messages

    Returns
    -------
    u, ne : floats
    '''
    
    lu, lne = [np.full(np.shape(lms), const.notnum) for i in range(2)]

    ind = np.where((lssfr > const.notnum) &
                   (lms > const.notnum) &
                   (loh12 > const.notnum))
    if (np.shape(ind)[1]>1):
        lne[ind] = (lssfr[ind] + 12.661 + 0.627*(lms[ind]-10))/1.753

        lu[ind] =  -2.316 - 0.360*(loh12[ind] -8.) -\
            0.292*lne[ind] + 0.428*(lssfr[ind] + 9)

    return lu, lne
    
def get_une(lms, lssfr, loh12, unemod='kashino20', verbose=False):
    '''
    Given log10(Mstar), log10(sSFR) and 12+log(O/H),
    get the ionizing parameter, U, and the electron density, ne.

    Parameters
    ----------
    lms : float
      log10(Mstar/Msun)
    lssfr : float
      log10(sSFR/yr), it should be an instantaneous measurement
    loh12 : float 
      12+log(O/H)
    unemod : string
      Model to go from galaxy properties to U and ne
    verbose : boolean
      Yes = print out messages

    Returns
    -------
    lu, lne : floats
    '''

    lu = None ; lne = None
    
    if (unemod == 'kashino20'):
        lu, lne = get_une_kashino20(lms,lssfr,loh12,verbose=verbose)
    else:
        print('STOP (eml_une): Unrecognised model to get U and ne.')
        print('                Possible unemod= {}'.format(const.unemods))
        sys.exit()
        
    return lu, lne
