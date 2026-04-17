"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
.. contributions:: Olivia Vidal <ovive.pro@gmail.com>
.. contributions:: Julen Expósito-Márquez <expox7@gmail.com>
"""
import h5py
import sys
import os
import glob
import numpy as np
import gne.gne_const as c

def stop_if_no_file(infile):
    '''
    It stops the program if a file does not exists

    Parameters
    -------
    infile : string
        Input file
    '''
    
    if (not os.path.isfile(infile)):
        print('STOP: no input file {}'.format(infile)) 
        sys.exit()
    return


def check_file(infile,verbose=False):
    '''
    It checks if a file exists

    Parameters
    -------
    infile : string
        Input file
    verbose : boolean
        If True print out messages

    Returns
    -------
    file_fine : boolean
        True when the file exists.
    '''
    
    file_fine = True  
    if (not os.path.isfile(infile)):
        file_fine = False
        if verbose:
            print('WARNING (gne_io.check_file): file not found {}'.
                  format(infile))

    return file_fine


def create_dir(outdir):
    '''
    Return True if directory already exists or it has been created
    '''
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except:
            print('WARNING (iotools.create_dir): problem creating directory ',outdir)
            return False
    return True


def get_outroot(snap,ending,outpath=None,verbose=False):
    '''
    Get path to output line data and the file root name

    Parameters
    -------
    root : string
        Root for input files
    ending : string
        Ending for input files
    snap: integer
        Simulation snapshot number
    outpath : string
        Path to output
    verbose : boolean
        If True print out messages

    Returns
    -------
    root, endf : string
        Root and ending of output line data files
    '''
    endf = ending
    if ending.endswith('.txt'):
        endf = ending[:-4] + '.hdf5'
    elif not ending.endswith('.hdf5'):
        endf = ending + '.hdf5'

    if outpath is None:
        opath = os.path.join(c.repo_dir, 'output')
    else:
        opath = outpath
    root = os.path.join(opath, 'iz' + str(snap), 'ivol')
        
    if not glob.glob(root + '*'):
        print(f'STOP: no adequate output directories {root}*') 
        sys.exit()

    if verbose:
        print(f'* Root to output: {root}')
    return root, endf


def get_metadata(filenom, verbose=True):
    """
    Extract cosmology and metadata from line file.
    
    Parameters
    ----------
    filenom : string
        Path to line file
    verbose : boolean
        If True print out messages
        
    Returns
    -------
    dict
        Dictionary containing:
        - redshift, omega0, omegab, lambda0, h0
        - vol
        - photmod_sfr, photmod_agn (if AGN exists)
        - AGN : boolean flag
    """
    f = h5py.File(filenom, 'r')
    header = f['header']
    boxside = header.attrs['boxside_Mpc']
    veff = header.attrs['eff_vol_Mpc3']

    metadata = {
        'redshift': header.attrs['redshift'],
        'omega0': header.attrs['omega0'],
        'omegab': header.attrs['omegab'],
        'lambda0': header.attrs['lambda0'],
        'h0': header.attrs['h0'],
        'vol_eff': veff,
        'photmod_sfr': header.attrs['photmod_sfr'],
        'AGN': 'agn_data' in f.keys(),
        'att': 'attmod' in header.attrs,
    }
    # AGN info
    if metadata['AGN']:
        metadata['photmod_agn'] = header.attrs['photmod_NLR']
    # Attenuation info
    if metadata['att']:
        metadata['attmod'] = header.attrs['attmod']
    # Flux info
    fsfr = f['sfr_data']
    required_flux_datasets =['Halpha_sfr_flux','Hbeta_sfr_flux',
                             'NII6584_sfr_flux','OIII5007_sfr_flux']
    metadata['flux'] = any(dataset in fsfr for dataset
                           in required_flux_datasets)
    f.close()

    # Message if reduced samples at input 
    eff_boxside = pow(veff, 1./3.)
    effdiff = abs(eff_boxside - boxside)
    if effdiff>1e-5 and verbose:
        effp = 100*veff/(boxside**3)
        print(f'    Side of the effective box = {eff_boxside:.1f} Mpc;') 
        print(f'     out of the original {boxside:.1f} Mpc ({effp:.1f}%)\n')
    
    return metadata


def get_plotfile(root,ending,plot_type):
    '''
    Get the path and name to a plot

    Parameters
    -------
    root : string
        Root to data to use for plotting
    plot_type : string
        Type of plot: bpt, etc

    Returns
    -------
    plotfile : string
        Path to plots
    '''
    plotpath = root.split('ivol')[0]+'plots/'
    create_dir(plotpath)

    nom = ending.split('.hdf5')[0]
    if nom.startswith('lines_'):
        nom = nom.split('lines_')[1]
    
    plotfile = plotpath + plot_type + '_' + nom + '.pdf'
    return plotfile


def get_outnom(filenom,dirf=None,nomf=None,verbose=False):
    '''
    Get output from a given filename

    Parameters
    -------
    filenom : string
        Name of file
    dirf : string
        Path to output
    nomf : string
        Name root for output file
    verbose : boolean
        If True print out messages

    Returns
    -------
    outfile : string
        Path to output file
    '''
    path, nomf1 = filenom.rsplit('/', 1)
    afteriz = path.split('iz')[-1]

    if nomf is None:
        nom = nomf1.split('.')[0]
    else:
        nom = nomf

    if dirf is None:
        dirf = 'output'
    #dirf = dirf.rsplit('/', 1)[0] + '/iz' + afteriz + '/'
    root = os.path.join(dirf,'iz'+afteriz)
    create_dir(root)

    outfile = os.path.join(root, nom+'.hdf5')
    if outfile == filenom:
        outfile = os.path.join(root, nom+'_output.hdf5')
        print(f'WARNING: outfile is the same as filenom {filenom}. Renaming outfile to {outfile}')

    if verbose:
        print(f'* Output: {outfile}')
    return outfile


def get_param(config, key, default):
    param = default
    if config:
        try:
            value = config.get(key)
            if value is not None: 
                param = value 
        except:
            param = default
    return param


def print_h5attr(infile,inhead='header'):
    """
    Print out the group attributes of a hdf5 file

    Parameters
    ----------
    infile : string
      Name of input file (this should be a hdf5 file)
    inhead : string
      Name of the group to read the attributes from

    Example
    -------
    >>> import h2s_io as io
    >>> infile = '/hpcdata0/simulations/BAHAMAS/AGN_TUNED_nu0_L100N256_WMAP9/Data/Snapshots/snapshot_026/snap_026.27.hdf5'
    >>> io.print_h5attr(infile,inhead='Units')
    """

    filefine = check_file(infile) #print(filefine)
    if (not filefine):
        print('WARNING (h2s_io.printh5attr): Check that the file provided is correct')
        return ' '
    
    f = h5py.File(infile, 'r')
    header = f[inhead]
    for hitem in list(header.attrs.items()): 
        print(hitem)
    f.close()

    return ' '


def get_nheader(infile,firstchar=None):
    '''
    Given a text file with a structure: header+data, 
    counts the number of header lines

    Parameters
    -------
    infile : string
        Input file

    Returns
    -------
    ih : integer
        Number of lines with the header text
    '''


    ih = 0
    with open(infile,'r') as ff:
        for line in ff:
            if not line.strip():
                # Count any empty lines in the header
                ih += 1
            else:
                sline = line.strip()
                
                # Check that the first character is not a digit
                char1 = sline[0]
                word1 = sline.split()[0]
                if not firstchar:
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
                else:
                    if char1 == firstchar:
                        ih+=1
    return ih
        


def generate_header(infile,redshift,snap,
                    h0,omega0,omegab,lambda0,
                    mp,boxside,effvol,
                    units_h0=False,outpath=None,
                    out_ending=None,verbose=True):
    """
    Generate the header of the file with the line data

    Parameters
    -----------
    infile : string
        Path to input
    zz: float
        Redshift of the simulation snapshot
    snap: integer
        Simulation snapshot number
    h0 : float
        Hubble constant divided by 100
    omega0 : float
        Matter density at z=0
    omegab : float
        Baryonic density at z=0
    lambda0 : float
        Cosmological constant z=0
    mp : float
        Simulation resolution, particle mass
    boxside : float
        Side of the simulation box
    effvol : float
        Simulation effective volume (if downsampled in input)
    units_h0: boolean
        True if input units with h
    outpath : string
        Path to output
    out_ending : string
        Name for output file, if different from root
    verbose : bool
        True for messages
 
    Returns
    -----
    filenom : string
       Full path to the output file
    """

    # Get the file name
    filenom = get_outnom(infile,dirf=outpath,nomf=out_ending,
                         verbose=verbose)

    # Change units if required
    if units_h0:
        mp = mp/h0
        boxside = boxside/h0
        effvol = effvol/(h0*h0*h0)

    # Generate the output file (the file is rewrtitten)
    hf = h5py.File(filenom, 'w')

    # Generate a header
    headnom = 'header'
    head = hf.create_dataset(headnom,(100,))
    head.attrs[u'redshift'] = redshift
    head.attrs[u'snapnum'] = snap
    head.attrs[u'h0'] = h0
    head.attrs[u'omega0'] = omega0
    head.attrs[u'omegab'] = omegab
    head.attrs[u'lambda0'] = lambda0
    head.attrs[u'mp_Msun'] = mp
    head.attrs[u'boxside_Mpc'] = boxside
    head.attrs[u'eff_vol_Mpc3'] = effvol
    hf.close()

    return filenom


def decode_string_list(raw_list):
    """
    Helper function to decode a list of strings from HDF5 attribute.
    Handles both byte strings (older h5py) and
    regular strings (newer h5py).
    """
    result = []
    for item in raw_list:
        if isinstance(item, bytes):
            result.append(item.decode('utf-8'))
        else:
            result.append(str(item))
    return result


def add2header(filenom,names,values,verbose=True):
    """
    Add attributes to header

    Parameters
    -----------
    filenom : string
        Path to file 
    names : list of strings
        Atribute names
    values: list
        Values of attributes
    verbose : bool
        True for messages
    """
    
    # Open the file header
    hf = h5py.File(filenom, 'a')
    head = hf['header']
    
    # Append attributes
    count = 0
    for ii, nom in enumerate(names):
        if nom is not None:
            head.attrs[nom] = values[ii]
            count += 1
    hf.close()

    if verbose:
        print(f'* gne_io.add2header: Appended {count} attributes out of {len(names)}')
    
    return count


def get_selection(infile, outfile, inputformat='hdf5',
                  cutcols=None, mincuts=[None], maxcuts=[None],
                  testing=False,verbose=False):
    '''
    Get indexes of selected galaxies

    Parameters
    ----------
    infile : strings
     List with the name of the input files. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    inputformat : string
     Format of the input file.
    cutcols : list
     Parameters to look for cutting the data.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    mincuts : strings
     Minimum value of the parameter of cutcols in the same index. All the galaxies below won't be considered.
    maxcuts : strings
     Maximum value of the parameter of cutcols in the same index. All the galaxies above won't be considered.
    verbose : boolean
      If True print out messages
    testing : boolean
      If True only run over few entries for testing purposes

    Returns
    -------
    selection : array of integers
    '''

    selection = None

    if cutcols is None:
        cutcols = [None]
    
    check_file(infile, verbose=verbose)

    if testing:
        limit = c.testlimit
    else:
        limit = None    

    if inputformat not in c.inputformats:
        if verbose:
            print('STOP (gne_io): Unrecognised input format.',
                  'Possible input formats = {}'.format(c.inputformats))
        sys.exit()
    elif inputformat=='hdf5':
        with h5py.File(infile, 'r') as hf:
            if len(cutcols)==0 or cutcols==[None]:
                key = [k for k in hf["data"].keys() if k != 'header'][0]
                ind = np.arange(len(hf["data"][key]))
            else:
                ind = np.arange(len(hf[cutcols[0]][:]))

            for i in range(len(cutcols)):
                if cutcols[i]:
                    param = hf[cutcols[i]][:]
                    mincut = mincuts[i]
                    maxcut = maxcuts[i]

                    if mincut and maxcut:
                        ind = np.intersect1d(ind,np.where((mincut<param)&(param<maxcut))[0])
                    elif mincut:
                        ind = np.intersect1d(ind,np.where(mincut<param)[0])
                    elif maxcut:
                        ind = np.intersect1d(ind,np.where(param<maxcut)[0])
            selection = ind[:limit]
    elif inputformat=='txt':
        ih = get_nheader(infile)

        if len(cutcols) == 0 or cutcols == [None]:
            ind = np.arange(len(np.loadtxt(infile, usecols=0, skiprows=ih)))
        else:
            ind = np.arange(len(np.loadtxt(infile,usecols=cutcols[0],skiprows=ih)))

        for i in range(len(cutcols)):
            if cutcols[i]:
                param = np.loadtxt(infile,usecols=cutcols[i],skiprows=ih)
                mincut = mincuts[i]
                maxcut = maxcuts[i]

                if mincut and maxcut:
                    ind = np.intersect1d(ind,np.where((mincut<param)&(param<maxcut))[0])
                elif mincut:
                    ind = np.intersect1d(ind,np.where(mincut<param)[0])
                elif maxcut:
                    ind = np.intersect1d(ind,np.where(param<maxcut)[0])
        selection = ind[:limit]
    else:
        if verbose:
            print('STOP (gne_io.get_selection): ',
                  'Input file has not been found.')
        sys.exit()

    return selection


def read_data(infile, cut, inputformat='hdf5', params=[None],
              testing=False, verbose=True):    
    '''
    Read input data per column/dataset
    
    Parameters
    ----------
    infile : string
       Name of the input file. 
    cut : array of integers
       List of indexes of the selected galaxies from the samples.
    inputformat : string
       Format of the input file.
    params : list of either integers or strings
       Inputs columns for text files or dataset name for hdf5 files.
    testing : boolean
       If True only run over few entries for testing purposes
    verbose : boolean
       If True print out messages.
     
    Returns
    -------
    outparams : array of floats
    '''

    check_file(infile, verbose=verbose)

    if inputformat not in c.inputformats:
        if verbose:
            print('STOP (gne_io): Unrecognised input format.',
                  'Possible input formats = {}'.format(c.inputformats))
        sys.exit()
    elif inputformat=='hdf5':
        with h5py.File(infile, 'r') as hf:
            ii = 0
            for nomparam in params:
                if (nomparam is not None):
                    prop = np.zeros(np.shape(cut))
                    try:
                        ###here what if Pos/vel as matrix?
                        prop = hf[nomparam][cut]
                    except:
                        if verbose:
                            print('\n WARNING (gne_io): no {} found in {}'.format(
                                nomparam,infile))

                    if (ii == 0):
                        outparams = prop
                    else:
                        outparams = np.vstack((outparams,prop))
                    ii += 1

    elif inputformat=='txt':
        ih = get_nheader(infile)
        outparams = np.loadtxt(infile,skiprows=ih,usecols=params)[cut].T

    return outparams


def read_and_add(f, prop_name, verbose=True):
    """
    Read property from HDF5, add if shape is (n, 2).
    
    Parameters
    ----------
    f : h5py.File
        Open HDF5 file object.
    prop_name : str
        Name of the dataset to read.
    verbose : bool
        If True, print warnings about unexpected shapes.
        
    Returns
    -------
    numpy.ndarray
        1D array of summed values or original data if already 1D.
    """
    data = f[prop_name][:]
    result = data
    
    # Check if it's 2D and has exactly 2 columns
    if data.ndim == 2:
        result = np.sum(data, axis=1)
        if data.shape[1] != 2:
            # Unexpected case: 2D but not 2 columns
            if verbose:
                print(f"  [WARN] '{prop_name}': Detected shape {data.shape}. Expected (n, 2) for summation.")        
        
    return result


def get_ncomponents(cols):
    '''
    Get the number of components to estimate the emission lines from

    Parameters
    ----------
    cols : list
      List of columns with properties per components

    Returns
    -------
    ncomp : integer
      Number of components (for example 2 for bulge and disk)
    '''
    
    ncomp = 1

    try:
        dum = np.shape(cols)[1]
        ncomp = np.shape(cols)[0]
    except:
        ncomp = 1
        print('STOP (gne_io.get_ncomponents): ',
              'Columns should be given as [[0,1,..][...]]')
        sys.exit()
        
    return ncomp


def read_sfrdata(infile, cols, cut, inputformat='hdf5',
                 testing=False, verbose=True):    
    '''
    Read input M, SFR and Z data for each component
    
    Parameters
    ----------
    infile : string
       Name of the input file.
    cols : list of either integers or strings
       Inputs columns for text files or dataset name for hdf5 files.
    cut : array of integers
       List of indexes of the selected galaxies from the samples.
    inputformat : string
       Format of the input file.
    testing : boolean
       If True only run over few entries for testing purposes
    verbose : boolean
       If True print out messages.
     
    Returns
    -------
    outms, outssfr, outzgas : array of floats
    '''

    check_file(infile, verbose=verbose)
    
    ncomp = get_ncomponents(cols)
    # Initialise output
    outms, outssfr, outzgas = [np.zeros((ncomp,cut.size)) for i in range(3)]

    # Read input data
    if inputformat not in c.inputformats:
        if verbose:
            print('STOP (gne_io): Unrecognised input format.',
                  'Possible input formats = {}'.format(c.inputformats))
        sys.exit()
    elif inputformat=='hdf5':
        with h5py.File(infile, 'r') as hf:
            for i in range(ncomp):
                if i==0:
                    ms = np.array([hf[cols[i][0]][:]])
                    ssfr = np.array([hf[cols[i][1]][:]])
                    zgas = np.array([hf[cols[i][2]][:]])
                else:
                    ms = np.append(ms,[hf[cols[i][0]][:]],axis=0)
                    ssfr = np.append(ssfr,[hf[cols[i][1]][:]],axis=0)
                    zgas = np.append(zgas,[hf[cols[i][2]][:]],axis=0)
    elif inputformat=='txt':
        ih = get_nheader(infile)            
        for i in range(ncomp):
            X = np.loadtxt(infile,usecols=cols[i],skiprows=ih).T
            
            if i==0:
                ms = np.array([X[0]])
                ssfr = np.array([X[1]])
                zgas = np.array([X[2]])
            else:
                ms = np.append(ms,[X[0]],axis=0)
                ssfr = np.append(ssfr,[X[1]],axis=0)
                zgas = np.append(zgas,[X[2]],axis=0)

    for i in range(ncomp):
        outms[i,:] = ms[i,cut]
        outssfr[i,:] = ssfr[i,cut]
        outzgas[i,:] = zgas[i,cut]

    return outms, outssfr, outzgas        



def read_mgas_hr(infile, cols, selection, inputformat='hdf5',
                 testing=False, verbose=True):    
    '''
    Read input Mgas and scalelenght for each component
    
    Parameters
    ----------
    infile : string
       Name of the input file.
    cols : list of either integers or strings
       Inputs columns for text files or dataset name for hdf5 files.
    selection : array of integers
       List of indexes of the selected galaxies from the samples.
    inputformat : string
       Format of the input file.
    testing : boolean
       If True only run over few entries for testing purposes
    verbose : boolean
       If True print out messages.
     
    Returns
    -------
    mgas, hr : array of floats
    '''

    check_file(infile, verbose=verbose)

    # Initialise output
    ncomp = get_ncomponents(cols)
    outm, outr = [np.zeros((ncomp,len(selection))) for i in range(2)]

    # Read input data
    if inputformat not in c.inputformats:
        if verbose:
            print('STOP (gne_io): Unrecognised input format.',
                  'Possible input formats = {}'.format(c.inputformats))
        sys.exit()
    elif inputformat=='hdf5':
        with h5py.File(infile, 'r') as hf:
            for i in range(ncomp):
                if i==0:
                    mgas = np.array([hf[cols[i][0]][:]])
                    hr = np.array([hf[cols[i][1]][:]])
                else:
                    mgas = np.append(mgas,[hf[cols[i][0]][:]],axis=0)
                    hr = np.append(hr,[hf[cols[i][1]][:]],axis=0)
    elif inputformat=='txt':
        ih = get_nheader(infile)            
        for i in range(ncomp):
            X = np.loadtxt(infile,usecols=cols[i],skiprows=ih).T
            
            if i==0:
                mgas = np.array([X[0]])
                hr = np.array([X[1]])
            else:
                mgas = np.append(mgas,[X[0]],axis=0)
                hr = np.append(hr,[X[1]],axis=0)

    for i in range(ncomp):
        outm[i,:] = mgas[i,selection]
        outr[i,:] = hr[i,selection]

    return outm, outr        


          
def get_mgas_hr(infile,selection,cols,r_type,
                h0=None,units_h0=False,
                re2hr=c.re2hr,r502re=c.r502re,rvir2r50=c.rvir2r50,
                inputformat='hdf5',
                testing=False,verbose=False):
    '''
    Get Mgas and scalelength in the adecuate units.

    Parameters
    ----------
    infile : string
       Name of the input file.
    cols : list of either integers or strings
       Inputs columns for text files or dataset name for hdf5 files.
    r_type : list of integers per component
       0 for scalelength; 1 for R50 or Reff; and 2 for a full radius
    selection : array of integers
       List of indexes of the selected galaxies from the samples.
    h0 : float
       Hubble constant divided by 100 (only needed if units_h0 is True)
    units_h0: boolean
       True if input units with h
    re2hr: float
       Constant, a, to get the scalength, hr, from
       an effective radius (3D), Re, as hr=a*Re
    r502re: float
       Constant, a, to get the effective radius (3D), Re, from
       a half-mass(light) radius (2D), R50, as Re=a*R50
    rvir2r50: float
       Constant, a, to get the half-mass(light) radius (2D), R50, from
       the halo radius, Rvir, as R50=a*Rvir
    inputformat : string
       Format of the input file.
    testing : boolean
       If True only run over few entries for testing purposes
    verbose : boolean
       If True print out messages.

    Returns
    -------
    mgas, hr : array of floats (Msun, Mpc)
    '''

    # Read Mgas and hr
    mgas, hr = read_mgas_hr(infile,cols,selection,
                            inputformat=inputformat,
                            testing=testing,verbose=verbose)
    if units_h0:
        mgas = mgas/h0
        hr = hr/h0

    # Initialise output with change of units
    outm = np.copy(mgas)
    outr = np.copy(hr)

    # Check that rtype is adequate
    ncomp = np.shape(outr)[0]
    if ((min(r_type)<0 or max(r_type)>3 or len(r_type)!=ncomp) and verbose):
        print('WARNING! Input r_type should be 0, 1, 2 or 3, per component.')

    # Correct scalelenght for each component
    for i in range(ncomp):
        if r_type[i] == 1:
            # Get the scalelenght from an effective (2D, projected) radius
            outr[i,:] = re2hr*hr[i,:]
        elif r_type[i] == 2:
            # Get the scalelenght from the half-mass(light) 3D radius
            outr[i,:] = re2hr*r502re*hr[i,:]
        elif r_type[i] == 3:
            # Get the scalelenght from the radius of the halo
            outr[i,:] = re2hr*r502re*rvir2r50*hr[i,:]

    return outm, outr


def write_data(filenom,group=None,params=None,params_names=None,
               params_labels=None,verbose=True):
    '''
    Write global data to output file

    Parameters
    ----------
    filenom : string
       Full path to the output file
    group : string
       Group in hdf5 file to write into
    params : list of strings (hdf5 files)
       Datasets to be written
    params_names : array of strings
       Names of the datasets to be written
    params_labels : array of strings
       Units and maybe descriptions for the parameters
    '''     
    # Open the file 
    hf = h5py.File(filenom, 'a')
    if group is not None:
        # Create the output group if needed
        if group not in hf.keys():
            gdat = hf.create_group(group)
        else:
            gdat = hf[group]

    if params is not None:
        if params_labels is None:
            params_labels = params_names
        
        for i in range(len(params)):
            nom = params_names[i]
            data = params[i]
            maxshape = tuple(None for _ in data.shape)

            if nom in gdat.keys():
                # Remove dataset if already present
                del gdat[nom]
            gdat.create_dataset(nom,data=data,maxshape=maxshape)
            gdat[nom].dims[0].label = params_labels[i]

    hf.close()
    return 


def write_global_data(filenom,lmass,mass_type='s',
                      lssfr=None,scalelength=None,
                      extra_param=None,extra_params_names=None,
                      extra_params_labels=None,verbose=True):
    '''
    Write global data to output file

    Parameters
    ----------
    filenom : string
       Full path to the output file
    mass : array of floats
       Masses (log10(M) (Msun)).
    mass_type : string
       's' for stellar, 'g' for gas, etc. indicating the type of mass
    ssfr : array of floats
       sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    scalelength : array of floats
        Scalelenth of the galaxies per component (log10(SFR/M*) (1/yr)).
    extra_params : list of integers (txt files) or strings (hdf5 files)
        Column number or names of the extra parameters to be saved.
    extra_params_names : array of strings
        Names of the datasets in the output files for the extra parameters.
    extra_params_labels : array of strings
        Descriptions (expected name and units) for the extra parameters.
    '''     
    # Open the file 
    hf = h5py.File(filenom, 'a')
    # Check if there is a group for global data and create it otherwise
    if 'data' not in hf.keys():
        gdat = hf.create_group('data')
    else:
        gdat = hf['data']

    maxshape = tuple(None for _ in lmass.shape)
    nom = 'lm_'+mass_type
    gdat.create_dataset(nom, data=lmass, maxshape=maxshape)
    gdat[nom].dims[0].label = 'log10(M/Msun)'

    if lssfr is not None:
        maxshape = tuple(None for _ in lssfr.shape)
        gdat.create_dataset('lssfr', data=lssfr, maxshape=maxshape)
        gdat['lssfr'].dims[0].label = 'log10(SFR/M*/yr)'

    if extra_param is not None:
        if extra_params_labels is None:
            extra_param_labels = extra_param_names
        
        for i in range(len(extra_param)):
            nom = extra_params_names[i]
            data = extra_param[i]
            maxshape = tuple(None for _ in data.shape)
            
            gdat.create_dataset(nom,data=data,maxshape=maxshape)
            gdat[nom].dims[0].label = extra_params_labels[i]

    hf.close()
    return 


def write_sfr_data(filenom,lu_sfr,lnH_sfr,lzgas_sfr,nebline_sfr,
                   verbose=True):
    '''
    Write line data from star forming regions

    Parameters
    ----------
    filenom : string
       Full path to the output file
    lu_sfr : floats
     U of the galaxies per component.
    lnH_sfr : floats
     nH of the galaxies per component (cm^-3).
    lzgas_sfr : floats
     Metallicity of the galaxies per component (12+log(O/H))
    nebline_sfr : floats
      Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
    '''

    # Read information on models
    f = h5py.File(filenom, 'r')   
    header = f['header']
    photmod_sfr = header.attrs['photmod_sfr']
    f.close()

    # Output data
    with h5py.File(filenom,'a') as hf:
        # SF data
        hfdat = hf.create_group('sfr_data')
        hfdat.create_dataset('lu_sfr', data=lu_sfr, maxshape=(None,None))
        hfdat['lu_sfr'].dims[0].label = 'log10(U) (dimensionless)'
    
        hfdat.create_dataset('lnH_sfr',data=lnH_sfr, maxshape=(None,None))
        hfdat['lnH_sfr'].dims[0].label = 'log10(nH) (cm**-3)'
    
        hfdat.create_dataset('lz_sfr', data=lzgas_sfr, maxshape=(None,None))
        hfdat['lz_sfr'].dims[0].label = 'log10(Z_cold_gas) (dimensionless)'

        for i in range(len(c.line_names[photmod_sfr])):           
            hfdat.create_dataset(c.line_names[photmod_sfr][i] + '_sfr', 
                                 data=nebline_sfr[:,i], maxshape=(None,None))
            hfdat[c.line_names[photmod_sfr][i] + '_sfr'].dims[0].label = \
                'erg s^-1 (Photoionisation model: Lsun = 3.826e+33 erg/s per unit SFR(Mo/yr) for 10^8yr)'
                
    return 


def write_agn_data(filenom,Lagn,lu_agn,lzgas_agn,
                   nebline_agn,epsilon_agn=None,
                   ew_notatt=None,ew_att=None,
                   verbose=True):
    '''
    Write line data from AGNs in output file

    Parameters
    ----------
    filenom : string
       Name of the output file.
    Lagn : numpy array
       Bolometric luminosity (erg/s)
    lu_agn : floats
     U of the galaxies per component.
    lzgas_agn : floats
     Metallicity of the galaxies per component (12+log(O/H))
    nebline_agn : array of floats
       Luminosities (erg/s)
    '''

    # Read information on models
    f = h5py.File(filenom, 'r')   
    header = f['header']
    photmod_agn = header.attrs['photmod_NLR']
    f.close()

    with h5py.File(filenom,'a') as hf:
        # AGN data
        hfdat = hf.create_group('agn_data')

        hfdat.create_dataset('Lagn', data=Lagn, maxshape=(None))
        hfdat['Lagn'].dims[0].label = 'L_bol (erg/s)'

        lu_agn_val = lu_agn.reshape(-1)
        hfdat.create_dataset('lu_agn', data=lu_agn_val, maxshape=(None))
        hfdat['lu_agn'].dims[0].label = 'log10(U) (dimensionless)'

        lzgas_agn_val = lzgas_agn.reshape(-1)
        hfdat.create_dataset('lz_agn', data=lzgas_agn_val, maxshape=(None))
        hfdat['lz_agn'].dims[0].label = 'log10(Z)'

        if epsilon_agn is not None:
            epsilon_agn_val = epsilon_agn.reshape(-1)
            hfdat.create_dataset('epsilon_NLR', data=epsilon_agn_val,
                                 maxshape=(None))
            hfdat['epsilon_NLR'].dims[0].label = \
                'AGN NLRs volume filling factor (dimensionless)'

        for i in range(len(c.line_names[photmod_agn])):
            ndata = nebline_agn[0,i,:]
            hfdat.create_dataset(c.line_names[photmod_agn][i] + '_agn',
                                 data=ndata, maxshape=(None))
            hfdat[c.line_names[photmod_agn][i] + '_agn'].dims[0].label = \
                'erg s^-1'
    return

