'''
# 4 arrays of random numbers: lines. File with 4 columns. Plot BPT to read and plot.
import numpy as np
import os
from matplotlib import pyplot as plt
from get_nebular_emission.eml_style import style1
import get_nebular_emission.eml_const as const
bpt_data = 'C:/Users/Olivia/PRUEBAS/bpt_data.dat' # r"example_data/bpt_data.dat"
header1 = 'line1   line2   line3   line4'

a = np.random.randint(1, 10, size=10)
b = np.random.randint(1, 10, size=10)
c = np.random.randint(1, 10, size=10)
d = np.random.randint(1, 10, size=10)

datatofile=np.column_stack((a, b, c, d))

with open(bpt_data,'w') as svfile:
    np.savetxt(svfile,datatofile,delimiter=' ',header=header1)
    svfile.closed
'''
import numpy as np
from get_nebular_emission.eml_io import get_nheader
import get_nebular_emission.eml_const as const

def get_limits():
    '''


    Returns
    -------
    limits of the parameters of the photoionization model

    '''

    # HERE: Hacer otra función que lea los límites.
    # Como una iteración porque hay tres propiedades n_lim !! VER SLACK propname bucle
    # Si hay una n en la string se considera que son los límites de la n.
    # En esa función va a tomar logaritmos así está todo preparado.
    # Iteración sobre las tres propiedades de forma que solo tenga una cosa.
    propname = ['n', 'Z', 'U']


def get_lines_Gutkin(loh12, lu, lne, verbose=False):
    '''
    Given 12+log(O/H), logU and logne,
    get the interpolations for the emission lines,
    using the tables
    from Gutkin et al. (2016) (https://arxiv.org/pdf/1607.06086.pdf)

    HERE : EXPLAIN HOW TO ADD LIMITS DEPENDING ON PHOTIO MODEL.

    Parameters
    ----------
    loh12 : float
      12+log(O/H)
    lu : float
      log(U)
    lne : float
      log(ne)
    photmods : string
      Model to go from U, Z and ne to emission lines luminosities.
    verbose : boolean
      Yes = print out messages

    Returns
    -------
    emission lines : floats
    '''

    flimits = r"nebular_data/limits_gutkin.txt"
    ih = get_nheader(flimits)

    uplimit = np.loadtxt(flimits, skiprows=ih, usecols=(0),unpack = True)
    lowlimit = np.loadtxt(flimits,skiprows=ih,usecols=(1),unpack = True)


    ind = np.where(lne>uplimit)
    lne[ind] = uplimit
    ind = np.where(lne<lowlimit)
    lne[ind] = lowlimit


    return lines


def get_lines(loh12, lu, lne, photmods='Gutkin16', verbose=False, Testing=False, Plotting=False):
    '''
    Given 12+log(O/H), logU and logne,
    get the interpolations for the emission lines

    Parameters
    ----------
    loh12 : float
      12+log(O/H)
    lu : float
      log(U)
    lne : float
      log(ne)
    photmods : string
      Model to go from U, Z and ne to emission lines luminosities.
    verbose : boolean
      Yes = print out messages

    Returns
    -------
      emission lines : floats
    '''


    if (photmods == 'Gutkin16'):
        lines = get_lines_Gutkin(loh12, lu, lne, verbose=verbose)
    else: #HERE cambiar
        print('STOP (eml_une): Unrecognised model to get emission lines.')
        print('                Possible unemod= {}'.format(const.photmods))
        exit()


    return lines

