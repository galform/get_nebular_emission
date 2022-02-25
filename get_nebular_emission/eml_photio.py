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
from get_nebular_emission.eml_io import check_file

def get_limits(infile, propname):
    '''
    Given a file with a structure: property + lower limit + upper limit,
    gets the limits of the parameters of the photoionization model.

    In the file we must find the properties well specified.
    The header lines have to start with '#'

    Parameters
    -------
    infile: string
        input file

    propname : string
        name of the property that we want

    Returns
    -------

    lower_limit: float
        lower limit of the requested property
    upper_limit: float
        upper limits of the requested property

    Examples
    -------
    infile:

    # Table 3 of Gutkin+2016 (https://arxiv.org/pdf/1607.06086.pdf)
    # Property      Lower_limit     Upper_limit
    Z               0.0001          0.040
    U               -4.0            -1.0
    xid             0.1             0.5
    nH               1               4


    >> get_limits(infile, propname = 'Z'):
        0.0001 0.04

    >> get_limits(infile, propname = 'nH'):
        1  4

    '''

    # HERE: Hacer otra función que lea los límites.
    # Como una iteración porque hay tres propiedades n_lim !! VER SLACK propname bucle
    # Si hay una n en la string se considera que son los límites de la n.
    # En esa función va a tomar logaritmos así está todo preparado.
    # Iteración sobre las tres propiedades de forma que solo tenga una cosa.


    prop = np.loadtxt(infile,dtype=str,comments='#',usecols=(0),unpack=True)
    prop = prop.tolist()
    ind = prop.index(propname)

    # Read the header of '#'
    ih = 0
    ff = open(infile, 'r')
    for line in ff:
        # Count any empty lines in the header
        if not line.strip():
            ih += 1
        else:
            sline = line.strip()
            # Check that the first character is a '#'
            char1 = sline[0]
            if char1 == '#':
                ih += 1

    lower_limit = np.loadtxt(infile, skiprows=ind+ih, max_rows=1, usecols=(1),unpack=True)
    upper_limit = np.loadtxt(infile,skiprows=ind+ih, max_rows=1,usecols=(2),unpack=True)

    return lower_limit,upper_limit

'''
    for line in ff:
        sline = line.strip()
        word1 = sline.split()[0]
        word2 = sline.split()[1]
        word3 = sline.split()[2]

        for ii, name in enumerate(propname):
            if propname[ii]==word1:
                lower_limits.append(float(word2))
                upper_limits.append(float(word3))

    return propname, lower_limits, upper_limits
    '''

def clean_photarray(limfile, infile, col_prop, propname, photmod='Gutkin16', verbose=True):
    '''

    Parameters
    ----------
    limfile
    infile
    cols: list
    propname
    photmod
    verbose

    Returns
    -------

    '''

    # Llama a la función que te dice los limites. Hacer otra función que limpie los límites.
    # Aqui poner un if que si el modelo es tal que haga tal cambio de unidades.
    # Para añadir un nuevo modelo, clean_photarray necesitará añadir cosas en el if respecto al cambio de unidades.
    # Aquí da array limpios, array 5 => array 5, con los elementos en los límites.

    lowerl,upperl = get_limits(infile=limfile, propname=propname)


    # Read the data file:
    check_file(infile, verbose=verbose)

    if ('.hdf5' in infile):
        print('Program not set to deal with hdf5 files yet, input a text file')
    else:
        # Text file
        ih = get_nheader(infile)
        if ('.cvs' in infile):
            deli = ','
        else:
            deli = None

        
        prop = np.loadtxt(infile, delimiter=deli, skiprows=ih, usecols = col_prop, unpack=True)
        #Z = np.loadtxt(infile,delimiter=deli, skiprows=ih,usecols=cols[1],unpack=True)
        #nH = np.loadtxt(infile,delimiter=deli, skiprows=ih,usecols=cols[2],unpack=True)

        #U = U.tolist()
        print(prop)

    ind = np.where(prop>upperl)
    prop[ind] = upperl
    ind = np.where(prop<lowerl)
    prop[ind] = lowerl

    # Change of units
    #if photmod == 'Gutkin16':
     #   prop = np.log10(prop) # Here: it works with nH but no with U, U is in log.



    '''
    ind = np.where(Z>upperl)
    Z[ind] = upperl
    ind = np.where(Z<lowerl)
    Z[ind] = lowerl

    ind = np.where(nH>upperl)
    nH[ind] = upperl
    ind = np.where(nH<lowerl)
    nH[ind] = lowerl
    '''

    return prop



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

    #low_limit,up_limit = get_photlimits(photmod=,prop=,)

    ind = np.where(lne>uplimit)
    lne[ind] = uplimit
    ind = np.where(lne<lowlimit)
    lne[ind] = lowlimit


    #return lines


def get_lines(in_loh12, in_lu, in_lne, photmods='Gutkin16', verbose=False, Testing=False, Plotting=False):
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

# Hacer un loop sobre las tres arrays y que te las limpie el código.
# loh12 = clean_photarrray(in_loh12=,photmod=,prop=loh12) Hacer para los tres
'''
    if (photmods == 'Gutkin16'):
        lines = get_lines_Gutkin(loh12, lu, lne, verbose=verbose)
    else: #HERE cambiar
        print('STOP (eml_une): Unrecognised model to get emission lines.')
        print('                Possible unemod= {}'.format(const.photmods))
        exit()


    return lines
'''
