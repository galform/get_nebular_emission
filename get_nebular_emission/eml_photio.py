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
import h5py
import numpy as np
from get_nebular_emission.eml_io import get_nheader
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_io import check_file


# DICCIONARIO DE MODELOS DE PHOTOIO: Nombre del modelo y nombre del fichero con los límites
# Si lía mucho matriz de strings

mod_lim = {'gutkin16': r"nebular_data/gutkin_tables/limits_gutkin.txt"}
# print(mod_lim.keys()) # To show the keys

def get_zfile(zmet, photmod='gutkin16'):

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


def get_limits(propname, photmod = 'gutkin16',verbose=True):
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

    try:
        infile = mod_lim[photmod]
    except KeyError:
        print('STOP (eml_photio): the {}'.format(photmod) + ' model is an unrecognised model in the dictionary mod_lim')
        print('                  Possible photmod= {}'.format(mod_lim.keys()))
        exit()

    # Check if the data file exists:
    check_file(infile, verbose=verbose)

    prop = np.loadtxt(infile,dtype=str,comments='#',usecols=(0),unpack=True)
    prop = prop.tolist()
    if propname not in prop:
        print('STOP (eml_photio): property {} '.format(propname)+'not found in the limits file {}'.format(infile))
        print('                   In the limits file we must find the properties written as: U, Z and nH')
        exit()
    else:
        ind = prop.index(propname)

        ih = get_nheader(infile,firstchar='#')
        lower_limit = np.loadtxt(infile, skiprows=ind+ih, max_rows=1, usecols=(1),unpack=True)
        upper_limit = np.loadtxt(infile,skiprows=ind+ih, max_rows=1,usecols=(2),unpack=True)
        return lower_limit,upper_limit

def clean_photarray(photmod='gutkin16', verbose=True):
    #clean_photarray(infile, in_Z, in_U, in_ne, photmod = 'gutkin16',verbose=True)
    #infile: siempre va a ser el que salga de eml_une, se podría poner el enlace y ya
    # lo que pasa es que solo tendríamos un archivo. Se podría introducir en el diccionario y así tener varios guardados
    # dependiendo del modelo. La idea además es que sea hdf5, siempre.
    #
    # in_Z, in_U, in_ne: Los nombres de las propiedades para buscarlas directamente en el archivo de límites. Siempre
    # van a ser las mismas columnas en el infile. Quizá ponerlo en un array: [Z, U, ne] = ['Z', 'U', 'nH']

    # No estoy teniendo en cuenta disco y bulbo: ¿Cómo hacer?
    # ¿Guardo en arrays loh12 = [disk, bulge], lu = [disk, bulge], lne = [disk, bulge]?

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



    # Read the data file:
    infile = r"output_data/U_ne_loh12.hdf5"

    check_file(infile, verbose=verbose)

    f = h5py.File(infile,'r')
    header = f['header']
    data = f['data']

    lu = data['lu'][:]
    lud=[]; lub=[]

    for ii, u in enumerate(lu):
        u1 = lu[ii]
        u1 = u1.tolist()
        lud.append(u1[0])
        lub.append((u1[1]))
    lud = np.array(lud)
    lub = np.array(lub)


    lne = data['lne'][:]
    lned=[];lneb=[]
    for ii, ne in enumerate(lne):
        ne1 = lne[ii]
        ne1 = ne1.tolist()
        lned.append(ne1[0])
        lneb.append((ne1[1]))
    lned = np.array(lned)
    lneb = np.array(lneb)


    loh12 = data['loh12'][:]
    loh12d=[];loh12b=[]
    for ii, oh12 in enumerate(loh12):
        oh121 = loh12[ii]
        oh121 = oh121.tolist()
        loh12d.append(oh121[0])
        loh12b.append((oh121[1]))
    loh12d = np.array(loh12d)
    loh12b = np.array(loh12b)


    lowerl, upperl = get_limits(propname='U', photmod=photmod)

    '''
    for ii, u in enumerate(lu):
        for iv in range(len(lu[ii])):
            ind = np.where(lu[ii,iv]>upperl)
            lu[ii, ind] = upperl
            ind = np.where(lu[ii,iv]<lowerl)
            lu[ii,ind] = lowerl
    print(lu)
    '''
    ind = np.where(lud > upperl)
    lud[ind] = upperl
    ind = np.where(lud < lowerl)
    lud[ind] = lowerl

    ind = np.where(lub > upperl)
    lub[ind] = upperl
    ind = np.where(lub < lowerl)
    lub[ind] = lowerl

    lu = np.stack((lud, lub))

    lowerl, upperl = get_limits(propname='nH', photmod=photmod)
    ind = np.where(lned > upperl)
    lned[ind] = upperl
    ind = np.where(lned < lowerl)
    lned[ind] = lowerl

    ind = np.where(lneb > upperl)
    lneb[ind] = upperl
    ind = np.where(lneb < lowerl)
    lneb[ind] = lowerl

    lne = np.stack((lned, lneb))

    lowerl, upperl = get_limits(propname='Z', photmod=photmod)
    ind = np.where(loh12d > upperl)
    loh12d[ind] = upperl
    ind = np.where(loh12d < lowerl)
    loh12d[ind] = lowerl

    ind = np.where(loh12b > upperl)
    loh12b[ind] = upperl
    ind = np.where(loh12b < lowerl)
    loh12b[ind] = lowerl

    loh12 = np.stack((loh12d, loh12b))

    return lu, lne, loh12
'''
    else:
        # Text file
        ih = get_nheader(infile)
        if ('.cvs' in infile):
            deli = ','
        else:
            deli = None

        
        prop = np.loadtxt(infile, delimiter=deli, skiprows=ih,usecols = col_prop, unpack=True)

    # Change of units
    #if photmod == 'Gutkin16':
     #   prop = np.log10(prop)  # Here: it works with nH but no with U, U is in log.
    # It must be consistent the limits.txt and the data.
'''



def get_lines_Gutkin(infile,in_loh12, in_U, in_ne, verbose=True):
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

    loh12 = clean_photarray(infile, propname= in_loh12, col_prop=(4),photmod='gutkin16', verbose=True)
    lu = clean_photarray(infile, propname= in_U,col_prop=(0), photmod='gutkin16', verbose=True)
    lne = clean_photarray(infile, propname= in_ne,col_prop=(2), photmod='gutkin16', verbose=True)
    #print(loh12,lu,lne)

    lines = 'interpolations done'

    return lines


def get_lines(infile, photmod='Gutkin16',verbose=False, Testing=False, Plotting=False):
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

    if (photmod == 'Gutkin16'):
        lines = get_lines_Gutkin(infile, verbose=verbose)
    else:
        print('STOP (eml_photio): Unrecognised model to get emission lines.')
        print('                Possible photmod= {}'.format(const.photmods))
        exit()


    return lines

