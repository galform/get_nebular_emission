from get_nebular_emission.eml_io import get_data
from get_nebular_emission.eml_une import get_une
import get_nebular_emission.eml_const as const
import get_nebular_emission.eml_plots as plot

def eml(infile, m_sfr_z=[0,1,2], h0=None,
        unemod='kashino20',photmod='gutkin16',
        LC2sfr=False,
        verbose=False, Testing=False):
    '''
    Calculate emission lines given the properties of model galaxies

    Parameters
    ----------
    infile : string
      Name of the input file. 
      In text files (*.dat, *txt, *.cat), columns separated by ' '.
      In csv files (*.csv), columns separated by ','.
    m_sfr_z : list
      [[component1_stellar_mass,sfr/LC,Z],[component2_stellar_mass,sfr/LC,Z],...]
      For text or csv files: list of integers with column position.
      For hdf5 files: list of data names.
    h0 : float
      If not None: value of h, H0=100h km/s/Mpc.
    unemod : string
      Model to go from galaxy properties to U and ne
    photmod : string
      Photoionisation model to be used for look up tables.
    LC2sfr : boolean
      If True Lyman Continuum photons expected as input for SFR.
    verbose : boolean
      Yes = print out messages

    Returns
    -------

    '''

    # Read the input data and correct it to the adequate units, etc.
    lms, lssfr, loh12 = get_data(infile, m_sfr_z, h0=h0,
                                 verbose=verbose, Testing=Testing)

    print(lssfr.min(), lssfr.max())

    # From the galaxy properties obtain the
    # ionizing parameter, U, and electron density, ne
    u, ne = get_une(lms, lssfr, loh12,
                    unemod=unemod, verbose=verbose)

    # From U and ne, obtain the emission lines from HII regions
    
