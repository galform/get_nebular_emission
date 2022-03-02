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


# DICCIONARIO DE MODELOS DE PHOTOIO: Nombre del modelo y nombre del fichero con los límites
# Si lía mucho matriz de strings

def get_zfile(zmet, photmod='gutkin'):

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
    root = 'nebular_data/' + photmod + '_tables/nebular_emission_Z'
    zfile = root + dec + '.txt'

    file_fine = check_file(zfile)
    if (not file_fine):
        zfile = None

    return zfile


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
    if propname not in prop:
        print('STOP: property {} '.format(propname)+'not found in the limits file {}'.format(infile))
        exit()
    else:
        ind = prop.index(propname)

        ih = get_nheader(infile,firstchar='#')
        lower_limit = np.loadtxt(infile, skiprows=ind+ih, max_rows=1, usecols=(1),unpack=True)
        upper_limit = np.loadtxt(infile,skiprows=ind+ih, max_rows=1,usecols=(2),unpack=True)
        return lower_limit,upper_limit

def clean_photarray(limfile, infile, col_prop, propname, photmod='Gutkin16', verbose=True):
    # Quitar limfile porque con el diccionario con Gutkin lo va a entender
    '''
    Given a file with a structure: property + lower limit + upper limit,
    a file with the data, the name of the property that we want and its column in the file well specify,
    gets the array of the data property with the necessary changes due to the limits.


    Parameters
    ----------
    limfile: string
        file with the limits
    infile: string
        file with the data
    col_prop: float
        column of the property in the infile
    propname: string
        name of the property
    photmod:  string
        Model to get the luminosity lines for the interpolations
    verbose : boolean
        Yes = print out messages


    Returns
    -------
    prop : array
        array of the property choose with all the data in the limits.
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

    # Change of units
    #if photmod == 'Gutkin16':
     #   prop = np.log10(prop)  # Here: it works with nH but no with U, U is in log.
    # It must be consistent the limits.txt and the data.

    ind = np.where(prop>upperl)
    prop[ind] = upperl
    ind = np.where(prop<lowerl)
    prop[ind] = lowerl

    return prop


def get_lines_Gutkin(infile, verbose=False):
    '''
    Given a file with the limits of the Gutkin model and a file with data,
    get 12+log(O/H), logU and logne to
    get the interpolations for the emission lines,
    using the tables
    from Gutkin et al. (2016) (https://arxiv.org/pdf/1607.06086.pdf)


    Parameters
    ----------
    infile: string
        file with the data
    photmods : string
      Model to go from U, Z and ne to emission lines luminosities.
    verbose : boolean
      Yes = print out messages

    Returns
    -------
    emission lines : floats
    '''

    limfile = r"/nebular_data/gutkin_tables/limits_gutkin.txt" # Estara en el diccionario, dado arriba
    # loh12 = clean_photarray(limfile, infile, col_prop= it is in the name of the file,popname= 'Z', photmod='Gutkin16', verbose=True)
    lu = clean_photarray(limfile, infile, col_prop=(0), photmod='Gutkin16', verbose=True)
    lne = clean_photarray(limfile, infile, col_prop=(2), photmod='Gutkin16', verbose=True)

    lines = 'interpolations done'

    return lines


def get_lines(infile, photomod='Gutkin16',verbose=False, Testing=False, Plotting=False):
    #(in_loh12, in_lu, in_lne, photmods='Gutkin16', verbose=False, Testing=False, Plotting=False):
    '''
    Given 12+log(O/H), logU and logne,
    get the interpolations for the emission lines

    Parameters
    ----------
    infile: string
        file with the data
    photmods : string
      Model to go from U, Z and ne to emission lines luminosities.
    verbose : boolean
      Yes = print out messages
    Plotting : boolean
      If True run verification plots with all data.
    Testing : boolean
      Yes = to only run over few entries for testing purposes

    Returns
    -------
      emission lines : floats
    '''

    # Hacer un loop sobre las tres arrays y que te las limpie el código.
    # loh12 = clean_photarrray(in_loh12=,photmod=,prop=loh12) Hacer para los tres

    if (photomod == 'Gutkin16'):
        lines = get_lines_Gutkin(infile, verbose=verbose)
    else:
        print('STOP (eml_photio): Unrecognised model to get emission lines.')
        print('                Possible unemod= {}'.format(const.photmods))
        exit()


    return lines

