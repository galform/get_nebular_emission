import sys
import numpy as np

def get_plots(lms, lssfr, loh12, verbose=False):
    '''
       Given log10(Mstar), log10(sSFR) and 12+log(O/H),
       get the 

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

