"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
.. contributions:: Olivia Vidal <ovive.pro@gmail.com>
.. contributions:: Julen Expósito-Márquez <expox7@gmail.com>
"""
import h5py
import numpy as np
import os,sys
import gne.gne_io as io
import gne.gne_const as c
import gne.gne_stats as st
from gne.gne_cosmology import emission_line_flux, logL2flux, set_cosmology

def get_zfile(zmet_str, photmod='gutkin16'):
    '''
    Given a metallicity string get the name of the corresponding table

    Parameters
    ----------
    zmet_str : string
        Metallicity name in files.
    photomod : string
        Name of the considered photoionisation model.

    Returns
    -------
    zfile : string
        Name of the model file with data for the given metallicity.
    '''
    root = os.path.join(c.data_dir, 'nebular_data',
                        photmod + '_tables', 'nebular_emission_Z')

    if len(zmet_str)<3:
        zmet_str = zmet_str+'0'
    zfile = root + zmet_str + '.txt'
    # Si son 2 numeros que le añada un cero
    file_fine = io.check_file(zfile)
    if (not file_fine):
        zfile = None

    return zfile
                

def get_limits(propname, photmod='gutkin16',verbose=True):
    '''
    Given a file with a structure: property + lower limit + upper limit,
    gets the limits of the parameters of the photoionization model.

    In the file we must find the properties well specified i.e U, Z and nH.
    The header lines have to start with '#'

    Parameters
    -------
    propname : string
        name of the property that we want
    photomod : string
        Name of the considered photoionisation model

    Returns
    -------
    lower_limit: float
        lower limit of the requested property
    upper_limit: float
        upper limits of the requested property

    Examples
    -------
    >>> get_limits(propname = 'nH', photmod = 'gutkin16')
        10  10000
    '''

    try:
        infile = os.path.join(c.data_dir,c.mod_lim[photmod])
    except KeyError:
        print('STOP (gne_photio.get_limits): the {}'.format(photmod) +
              ' is an unrecognised model in the dictionary mod_lim')
        print('                  Possible photmod= {}'.format(c.mod_lim.keys()))
        sys.exit()

    # Check if the limits file exists:
    io.check_file(infile, verbose=verbose)

    # Properties
    prop = np.loadtxt(infile,dtype=str,comments='#',usecols=(0),unpack=True)
    prop = prop.tolist()
    if propname not in prop:
        print('STOP (gne_photio.get_limits): property {} '.format(propname)+
              'not found in the limits file {}'.format(infile))
        sys.exit()
    else:
        ind = prop.index(propname)

        ih = io.get_nheader(infile,firstchar='#')
        lower_limit = np.loadtxt(infile, skiprows=ind+ih, max_rows=1, usecols=(1),unpack=True)
        upper_limit = np.loadtxt(infile,skiprows=ind+ih, max_rows=1,usecols=(2),unpack=True)
        return float(lower_limit),float(upper_limit)


def get_Zgrid(zgrid_str):
    '''
    Get the metallicity values from the file names

    Params
    -------
    zgrid_str : array or list of strings
        String with decimal part of the metallicity grid
      
    Returns
    -------
    nz : integer
        Number of elements of zgrid_str
    zgrid, lzgrid : array of floats
        Numerical value of the metallicity grid, Z, and log10(Z)
    '''
    nz = len(zgrid_str)

    zgrid = np.full(nz,c.notnum)
    zgrid = np.array([float('0.' + zz) for zz in zgrid_str])

    lzgrid = np.full(nz, c.notnum)
    ind = np.where(zgrid > 0.)
    if (np.shape(ind)[1] > 0):
        lzgrid[ind] = np.log10(zgrid[ind])

    return nz,zgrid,lzgrid


def read_gutkin16_grids(xid_phot, co_phot, imf_cut_phot):
    """
    Read the photoionisation model tables from Gutkin+16
    into matrices per each value of nH.
    
    Parameters
    ----------
    xid_phot : float
       Dust-to-metal ratio
    co_phot : float
       C/O ratio
    imf_cut_phot : float
       Solar mass high limit for the IMF
        
    Returns
    -------
    emline_grid1 : matrix, shape (nZ_reduced,nU,nELines)
        Grid for nH = 10
    emline_grid2 : matrix, shape (nZ,nU,nELines)
        Grid for nH = 100
    emline_grid3 : matrix, shape (nZ,nU,nELines)
        Grid for nH = 1000
    emline_grid4 : matrix, shape (nZ_reduced,nU,nELines)
        Grid for nH = 10000
    """
    
    photmod = 'gutkin16'

    # Read line names
    line_names = c.line_names[photmod]
    nemline = len(line_names)

    # Read grid of Zs
    zmet_str = c.zmet_str[photmod]
    nzmet, zmets, zedges = get_Zgrid(zmet_str)

    # Read reduced grid of Zs and map its indeces to the full grid
    zmet_str_reduced = c.zmet_str_reduced[photmod]
    nzmet_reduced, zmets_reduced, zredges = get_Zgrid(zmet_str_reduced)
    k_to_kred = {0: 0, 4: 1, 9: 2, 12: 3}

    # Read grid of Us
    uedges = c.lus_bins[photmod]
    nu = len(uedges)
    
    # Read grid of nHs
    nHbins = c.nH_bins[photmod]
    nHedges = np.array([np.log10(val) for val in nHbins])
    nnH = len(nHbins)

    # Store the photoionisation model tables into matrices
    emline_grid1 = np.zeros((nzmet_reduced,nu,nemline))
    emline_grid2 = np.zeros((nzmet,nu,nemline))
    emline_grid3 = np.zeros((nzmet,nu,nemline))
    emline_grid4 = np.zeros((nzmet_reduced,nu,nemline))    
    for k, zname in enumerate(zmet_str):
        infile = get_zfile(zname, photmod=photmod)
        io.check_file(infile, verbose=True)
        ih = io.get_nheader(infile)
        
        # Read all lines after header and extract columns
        data = np.loadtxt(infile, skiprows=ih)
        u = data[:, 0]          # log(Us)
        xid = data[:, 1]        # xid
        nH = data[:, 2]         # nh
        co = data[:, 3]         # (C/O)/(C/O)sol
        imf_cut = data[:, 4]    # mup
        emission_lines = data[:, 5:]  # All emission line data
        
        # Create mask for matching conditions
        mask = (xid == xid_phot) & (co == co_phot) & (imf_cut == imf_cut_phot)
        filtered_indices = np.where(mask)[0]
        
        # Process filtered data
        for idx in filtered_indices:
            if nH[idx] not in nHbins:
                continue
                
            # Find index for the read u value
            l = np.where(uedges == u[idx])[0]
            if len(l) == 0:
                continue
            l = l[0]
            
            # Get emission line values
            em_values = emission_lines[idx]
            
            # Fill appropriate grid based on nH value
            if nH[idx] == 10:
                kred = k_to_kred.get(k)
                emline_grid1[kred, l, :] = em_values
            elif nH[idx] == 100:
                emline_grid2[k, l, :] = em_values
            elif nH[idx] == 1000:
                emline_grid3[k, l, :] = em_values
            elif nH[idx] == 10000:
                kred = k_to_kred.get(k)
                emline_grid4[kred, l, :] = em_values
                
    return emline_grid1, emline_grid2, emline_grid3, emline_grid4


def get_lines_gutkin16(lu, lzgas, filenom, lnH=None,
                       origin='sfr', verbose=True):
    '''
    Get emission line luminosities, using tables from
    Gutkin et al. (2016) (https://arxiv.org/pdf/1607.06086.pdf)
    Units: Lbolsun per unit SFR(Msun/yr) for 10^8yr, assuming Chabrier
    
    Parameters
    ----------
    lu : array of floats
       Ionising parameter of the ionising region(s), log10(Us)
    lzgas : array of floats
       Metallicity of the ionising region(s), log10(Zgas)
    filenome : string
        Name of the file with relevant information    
    lnH : array of floats (or None)
       Hydrogen or electron density, log10(nH/cm^-3)
    origin : string
        Type of ionising region.    
    verbose : boolean
       If True print out messages
      
    Returns
    -------
    nebline : array of floats      
    '''

    # Read relevant constanst from file
    f = h5py.File(filenom, 'r')
    header = f['header']
    photmod  = header.attrs['photmod_'+origin]
    xid_phot = header.attrs['xid_'+origin]
    co_phot  = header.attrs['co_'+origin]
    imf_cut  = header.attrs['imf_cut_'+origin]
    if lnH is None: # Get constant value
        nH = header.attrs['nH_'+origin+'_cm3']
        lnH = np.zeros(lu.shape); lnH.fill(np.log10(nH))
    f.close()

    if (photmod != 'gutkin16'):
        if verbose:
            print('STOP (gne_lines_gutkin16): Photoionisation model mismatch.')
        return None

    emline_grid1, emline_grid2, \
        emline_grid3, emline_grid4 = read_gutkin16_grids(
            xid_phot, co_phot, imf_cut)

    # Initialize the matrix to store the emission lines
    ndat = lu.shape[1]
    ncomp = lu.shape[0]
    nemline = emline_grid1.shape[2]
    nebline = np.zeros((ncomp,nemline,ndat))

    # Edges of Z, U and nH grids
    zmet_str = c.zmet_str[photmod]
    nzmet, zmets, zedges = get_Zgrid(zmet_str)
    zmet_str_reduced = c.zmet_str_reduced[photmod]
    nzmet_reduced, zmets_reduced, zredges = get_Zgrid(zmet_str_reduced)
    uedges = c.lus_bins[photmod]
    nHbins = c.nH_bins[photmod]
    nHedges = np.array([np.log10(val) for val in nHbins])

    # Interpolate in all three grids: logUs, logZ, nH
    for comp in range(ncomp):
        ucomp = lu[comp,:]; zcomp=lzgas[comp,:]; nHcomp = lnH[comp,:]
        ind = np.where((ucomp > c.notnum) &
                       (zcomp > c.notnum) &
                       (nHcomp > c.notnum))[0]
        if (ind.size < 1):
            print('WARNING (get_lines_gutkin16):',
                  'no adequate log(Us)+ found for component',comp)
            return nebline
        
        # Initialize matrices with interpolated values on Z and U
        uu = ucomp[ind]
        zz = zcomp[ind]
        ngal = uu.size
        int1_zu, int2_zu, int3_zu, int4_zu = [np.zeros((ngal,nemline)) for i in range(4)]

        # Interpolate over Zgas and U
        int1_zu = st.bilinear_interpl(zz,uu,zredges,uedges,emline_grid1)
        int4_zu = st.bilinear_interpl(zz,uu,zredges,uedges,emline_grid4)

        int2_zu = st.bilinear_interpl(zz,uu,zedges,uedges,emline_grid2)
        int3_zu = st.bilinear_interpl(zz,uu,zedges,uedges,emline_grid3)

        # Interpolate over nH
        nn = nHcomp[ind]
        nHd, inH = st.interpl_weights(nn,nHedges)
        if verbose:
            print('Comp ',comp,': <logUs>=',np.average(uu),
                  '; <logZ>=',np.average(zz),
                  '; <lognH>=',np.average(nn))        

        c0,c1 = [np.zeros((ngal,nemline)) for i in range(2)]

        mask0 = (inH == 0)
        mask1 = (inH == 1)
        mask2 = (inH == 2)
        c0[mask0,:] = int1_zu[mask0,:]; c1[mask0,:] = int2_zu[mask0,:]
        c0[mask1,:] = int2_zu[mask1,:]; c1[mask1,:] = int3_zu[mask1,:]
        c0[mask2,:] = int3_zu[mask2,:]; c1[mask2,:] = int4_zu[mask2,:] 

        nebline[comp,:,ind] = c0*(1-nHd[:, np.newaxis]) + c1*nHd[:, np.newaxis]
        
    return nebline


def read_feltre16_grids(xid_phot, alpha_phot):
    """
    Read the photoionisation model tables from Feltre+16
    into matrices per each value of nH.
    Units: erg/s for a central Lacc=10^45 erg/s
    
    Parameters
    ----------
    xid_phot : float
     Dust-to-metal ratio for the Feltre et. al. photoionisation model.
    alpha_phot : float
     Alpha value for the Feltre et. al. photoionisation model.
        
    Returns
    -------
    emline_grid1 : matrix, shape (nZ,nU,nELines)
        Grid for nH = 100
    emline_grid2 : matrix, shape (nZ,nU,nELines)
        Grid for nH = 1000
    emline_grid3 : matrix, shape (nZ,nU,nELines)
        Grid for nH = 10000
    """
    photmod = 'feltre16'
    
    # Read line names
    line_names = c.line_names[photmod]
    nemline = len(line_names)

    # Read grid of Zs
    zmet_str = c.zmet_str[photmod]
    nzmet, zmets, zedges = get_Zgrid(zmet_str)

    # Read grid of Us
    uedges = c.lus_bins[photmod]
    nu = len(uedges)

    # Read grid of nHs
    nHbins = c.nH_bins[photmod]
    nHedges = np.array([np.log10(val) for val in nHbins])
    nnH = len(nHbins)
    
    # Store the photoionisation model tables into matrices
    emline_grid1 = np.zeros((nzmet,nu,nemline))
    emline_grid2 = np.zeros((nzmet,nu,nemline))
    emline_grid3 = np.zeros((nzmet,nu,nemline))

    for k, zname in enumerate(zmet_str):
        infile = get_zfile(zname,photmod=photmod)
        io.check_file(infile,verbose=True)
        ih = io.get_nheader(infile)

        # Read all lines after header and extract columns
        data = np.loadtxt(infile, skiprows=ih)
        u = data[:, 0]          # log(Us)
        xid = data[:, 1]        # xid
        nH = data[:, 2]         # nh
        alpha = data[:, 3]      # alpha
        emission_lines = data[:, 4:]  # All emission line data
        
        # Create mask for matching conditions
        mask = (xid == xid_phot) & (alpha == alpha_phot)
        filtered_indices = np.where(mask)[0]

        # Process filtered data
        for idx in filtered_indices:
            if nH[idx] not in nHbins: continue

            # Find index for the read u value
            l = np.where(uedges == u[idx])[0]
            if len(l) == 0: continue
            l = l[0] 
            
            # Get emission line values
            em_values = emission_lines[idx]

            # Determine grid and index based on nH value
            if nH[idx] == 100:
                emline_grid1[k, l, :] = em_values
            elif nH[idx] == 1000:
                emline_grid2[k, l, :] = em_values
            elif nH[idx] == 10000:
                emline_grid3[k, l, :] = em_values
    return emline_grid1, emline_grid2, emline_grid3


def get_lines_feltre16(lu, lzgas, filenom, lnH=None,
                       origin='NLR',verbose=True):
    '''
    Get emission line luminosities, using tables from
    Feltre+2016 (https://arxiv.org/pdf/1511.08217)
    Units: Lsun for L_AGN = 10^45 erg/s
    
    Parameters
    ----------
    lu : array of floats
       Ionising parameter of the ionising region(s), log10(Us)
    lzgas : array of floats
       Metallicity of the ionising region(s), log10(Zgas)
    filenom : string
        Name of the file with relevant information    
    lnH : array of floats (or None)
       Hydrogen or electron density, log10(nH/cm^-3)
    origin : string
        Type of ionising region.    
    verbose : boolean
       If True print out messages
      
    Returns
    -------
    nebline : array of floats      
    '''
    # Read relevant constanst from file
    f = h5py.File(filenom, 'r')
    header = f['header']
    photmod    = header.attrs['photmod_'+origin]
    alpha_phot = header.attrs['alpha_'+origin]    
    xid_phot   = header.attrs['xid_'+origin]
    if lnH is None: # Get constant value
        nH = header.attrs['nH_'+origin+'_cm3']
        lnH = np.zeros(lu.shape); lnH.fill(np.log10(nH))
    f.close()

    if (photmod != 'feltre16'):
        if verbose:
            print('STOP (gne_lines_feltre16): Photoionisation model mismatch.')
        return None
    
    emline_grid1, emline_grid2, emline_grid3 = read_feltre16_grids(
            xid_phot, alpha_phot)

    # Initialize the matrix to store the emission lines
    ndat = lu.shape[1]
    ncomp = lu.shape[0]
    nemline = emline_grid1.shape[2]
    nebline = np.zeros((ncomp,nemline,ndat))

    # Edges of Z, U and nH grids
    zmet_str = c.zmet_str[photmod]
    nzmet, zmets, zedges = get_Zgrid(zmet_str)
    uedges = c.lus_bins[photmod]
    nHbins = c.nH_bins[photmod]
    nHedges = np.array([np.log10(val) for val in nHbins])

    # Interpolate in all three grids: logUs, logZ, nH
    for comp in range(ncomp):
        ucomp = lu[comp,:]; zcomp=lzgas[comp,:]; nHcomp = lnH[comp,:]
        ind = np.where((ucomp > c.notnum) &
                       (zcomp > c.notnum) &
                       (nHcomp > c.notnum))[0]
        if (ind.size < 1):
            print('WARNING (get_lines_feltre16):',
                  'no adequate log(Us)+ found for component',comp)
            return nebline

        # Initialize matrices with interpolated values on Z and U
        uu = ucomp[ind]
        zz = zcomp[ind]
        ngal = uu.size
        int1_zu, int2_zu, int3_zu = [np.zeros((ngal,nemline)) for i in range(3)]

        # Interplate over Zgas and U
        int1_zu = st.bilinear_interpl(zz,uu,zedges,uedges,emline_grid1)
        int2_zu = st.bilinear_interpl(zz,uu,zedges,uedges,emline_grid2)
        int3_zu = st.bilinear_interpl(zz,uu,zedges,uedges,emline_grid3)    
         
        # Interpolate over nH
        nn = nHcomp[ind]
        nHd, inH = st.interpl_weights(nn,nHedges)
        if verbose:
            print(comp,'<logUs>=',np.average(uu),
                  '<logZ>=',np.average(zz),
                  '<lognH>=',np.average(nn))        

        c0,c1 = [np.zeros((ngal,nemline)) for i in range(2)]
        
        mask0 = (inH == 0)
        mask1 = (inH == 1)
        c0[mask0,:] = int1_zu[mask0,:]; c1[mask0,:] = int2_zu[mask0,:]
        c0[mask1,:] = int2_zu[mask1,:]; c1[mask1,:] = int3_zu[mask1,:]
        
        nebline[comp,:,ind] = c0*(1-nHd[:, np.newaxis]) + c1*nHd[:, np.newaxis]

    return nebline


def get_lines(lu, lzgas, filenom, lnH=None, origin='sfr',
              photmod='gutkin16', verbose=True):
    '''
    Get the emission line luminosity per ionising region.
    Units depend on the photoionisation model.

    Parameters
    ----------
    lu : array of floats
       Ionising parameter of the ionising region(s), log10(Us)
    lzgas : array of floats
       Metallicity of the ionising region(s), log10(Zgas)
    filenome : string
        Name of the file with relevant information    
    lnH : array of floats (or None)
       Hydrogen or electron density, log10(nH/cm^-3)
    origin : string
        Type of ionising region.    
    photomod : string
       Name of the considered photoionisation model
    verbose : boolean
       If True print out messages

    Returns
    -------
    nebline : array of floats
    '''

    if photmod not in c.photmods:
        if verbose:
            print('STOP (gne_photio.get_lines): Unrecognised model to get emission lines.')
            print('                Possible photmod= {}'.format(c.photmods))
        sys.exit()
    elif (photmod == 'gutkin16'):
        nebline = get_lines_gutkin16(lu,lzgas, filenom, lnH=lnH,
                                     origin=origin,verbose=verbose)
    elif (photmod == 'feltre16'):
        nebline = get_lines_feltre16(lu,lzgas, filenom, lnH=lnH,
                                     origin=origin,verbose=verbose)
    return nebline

