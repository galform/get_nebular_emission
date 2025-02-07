"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
.. contributions:: Olivia Vidal <ovive.pro@gmail.com>
.. contributions:: Julen Expósito-Márquez <expox7@gmail.com>
"""
import h5py
import numpy as np
import sys
#import warnings
import src.gne_io as io
import src.gne_const as c
import src.gne_stats as st
from src.gne_cosmology import emission_line_flux, logL2flux, set_cosmology


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

    root = 'data/nebular_data/' + photmod + '_tables/nebular_emission_Z'
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
        infile = c.mod_lim[photmod]
    except KeyError:
        print('STOP (gne_photio.get_limits): the {}'.format(photmod) +
              ' is an unrecognised model in the dictionary mod_lim')
        print('                  Possible photmod= {}'.format(c.mod_lim.keys()))
        sys.exit()

    # Check if the limits file exists:
    io.check_file(infile, verbose=verbose)
    #print(infile)

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


    
def calculate_flux(nebline,filenom,origin='sfr'):
    '''
    Get the fluxes for the emission lines given the luminosity and redshift.

    Params
    -------
    nebline : array of floats
        Luminosities of the lines per component.
        Lsun for L_AGN = 10^45 erg/s
    filenom : string
        Name of file with output
    origin : string
        Emission source (star-forming region or AGN).
      
    Returns
    -------
    fluxes : floats
        Array with the fluxes of the lines per component.
    '''
    
    if nebline.any():
        # Read redshift and cosmo
        f = h5py.File(filenom, 'r')
        header = f['header']
        redshift = header.attrs['redshift']
        h0 = header.attrs['h0']
        omega0 = header.attrs['omega0']
        omegab = header.attrs['omegab']
        lambda0 = header.attrs['lambda0']
        f.close()
        
        set_cosmology(omega0=omega0, omegab=omegab,lambda0=lambda0,h0=h0)
        
        luminosities = np.zeros(nebline.shape)
        luminosities[nebline>0] = np.log10(nebline[nebline>0]*h0**2)
        if (origin=='agn') and (luminosities.shape[0]==2):
            luminosities[1] = 0
            
        fluxes = np.zeros(luminosities.shape)
        for comp in range(luminosities.shape[0]):
            for i in range(luminosities.shape[1]):
                for j in range(luminosities.shape[2]):
                    if luminosities[comp,i,j] == 0:  
                        fluxes[comp,i,j] = 0
                    else:
                        fluxes[comp,i,j] = logL2flux(luminosities[comp,i,j],redshift)
    else:
        fluxes = np.copy(nebline)
            
    return fluxes


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


def get_lines_gutkin16(lu, lnH, lzgas, xid_phot=0.3,
                       co_phot=1,imf_cut_phot=100,verbose=True):
    '''
    Get the interpolations for the emission lines,
    using the tables
    from Gutkin et al. (2016) (https://arxiv.org/pdf/1607.06086.pdf)

    Parameters
    ----------
    lu : floats
     U of the galaxies per component.
    lnH : floats
     ne of the galaxies per component (cm^-3).
    lzgas : floats
     Metallicity of the galaxies per component (log10(Z))
    xid_phot : float
       Dust-to-metal ratio
    co_phot : float
       C/O ratio
    imf_cut_phot : float
       Solar mass high limit for the IMF
    verbose : boolean
       If True print out messages
      
    Returns
    -------
    nebline : array of floats
       Line luminosity per component
       Units: Lbolsun per unit SFR(Msun/yr) for 10^8yr, assuming Chabrier
    '''
    
    photmod = 'gutkin16'

    # Read line names
    line_names = c.line_names[photmod]
    nemline = len(line_names)

    # Initialize the matrix to store the emission lines
    ndat = lu.shape[1]
    ncomp = lu.shape[0]
    nebline = np.zeros((ncomp,nemline,ndat)); nebline.fill(c.notnum)

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
        infile = get_zfile(zname,photmod=photmod)
        io.check_file(infile,verbose=True)
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
            if nH[idx] not in nHbins: continue

            # Find index for the read u value
            l = np.where(uedges == u[idx])[0]
            if len(l) == 0: continue
            l = l[0] 
            
            # Get emission line values
            em_values = emission_lines[idx]

            # Determine grid and index based on nH value
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

    # Interpolate in all three grids: logUs, logZ, nH
    for comp in range(ncomp):
        ucomp = lu[comp,:]; zcomp=lzgas[comp,:]; nHcomp = lnH[comp,:]
        ind = np.where((ucomp > c.notnum) &
                       (zcomp > c.notnum) &
                       (nHcomp > c.notnum))[0]
        if (ind.size < 1):
            print('WARNING (get_lines_gutkin16): no adequate log(Us) found')
            return nebline
        
        # Initialize matrices with interpolated values on Z and U
        uu = ucomp[ind]
        zz = zcomp[ind]
        ngal = uu.size
        int1_zu, int2_zu, int3_zu, int4_zu = [np.zeros((ngal,nemline)) for i in range(4)]

        # Interplate over Zgas and U
        int1_zu = st.bilinear_interpl(zz,uu,zredges,uedges,emline_grid1)
        int4_zu = st.bilinear_interpl(zz,uu,zredges,uedges,emline_grid4)

        int2_zu = st.bilinear_interpl(zz,uu,zedges,uedges,emline_grid2)
        int3_zu = st.bilinear_interpl(zz,uu,zedges,uedges,emline_grid3)

        # Interpolate over nH
        nn = nHcomp[ind]
        nHd, inH = st.interpl_weights(nn,nHedges)

        c0,c1 = [np.zeros((ngal,nemline)) for i in range(2)]

        mask0 = (inH == 0)
        mask1 = (inH == 1)
        mask2 = (inH == 2)
        c0[mask0,:] = int1_zu[mask0,:]; c1[mask0,:] = int2_zu[mask0,:]
        c0[mask1,:] = int2_zu[mask1,:]; c1[mask1,:] = int3_zu[mask1,:]
        c0[mask2,:] = int3_zu[mask2,:]; c1[mask2,:] = int4_zu[mask2,:] 

        nebline[comp,:,ind] = c0*(1-nHd[:, np.newaxis]) + c1*nHd[:, np.newaxis]
        
    return nebline


def get_lines_feltre16(lu, lnH, lzgas, xid_phot=0.5,
                     alpha_phot=-1.7,verbose=True):
    '''
    Get the interpolations for the emission lines,
    using the tables from Feltre+2016 (https://arxiv.org/pdf/1511.08217)
    
    lnH : floats
     ne of the galaxies per component (cm^-3).
    lzgas : floats
     Metallicity of the galaxies per component (log10(Z))
    xid_phot : float
     Dust-to-metal ratio for the Feltre et. al. photoionisation model.
    alpha_phot : float
     Alpha value for the Feltre et. al. photoionisation model.
    verbose : boolean
      If True print out messages
      
    Returns
    -------
    nebline : array of floats
       Line luminosities per galaxy component.
       Units: Lsun for L_AGN = 10^45 erg/s
    '''

    photmod = 'feltre16'
    
    # Read line names
    line_names = c.line_names[photmod]
    nemline = len(line_names)

    # Initialize the matrix to store the emission lines
    ndat = lu.shape[1]
    ncomp = lu.shape[0]
    nebline = np.zeros((ncomp,nemline,ndat))
    
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

    # Interpolate in all three grids: logUs, logZ, nH
    for comp in range(ncomp):
        ucomp = lu[comp,:]; zcomp=lzgas[comp,:]; nHcomp = lnH[comp,:]
        ind = np.where((ucomp > c.notnum) &
                       (zcomp > c.notnum) &
                       (nHcomp > c.notnum))[0]
        if (ind.size < 1):
            print('WARNING (get_lines_feltre16): no adequate log(Us) found')
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
        
        c0,c1 = [np.zeros((ngal,nemline)) for i in range(2)]
        
        mask0 = (inH == 0)
        mask1 = (inH == 1)
        c0[mask0,:] = int1_zu[mask0,:]; c1[mask0,:] = int2_zu[mask0,:]
        c0[mask1,:] = int2_zu[mask1,:]; c1[mask1,:] = int3_zu[mask1,:]
        
        nebline[comp,:,ind] = c0*(1-nHd[:, np.newaxis]) + c1*nHd[:, np.newaxis]
                        
    return nebline


def get_lines(lu, lnH, lzgas, photmod='gutkin16',xid_phot=0.3,
              co_phot=1,imf_cut_phot=100,alpha_phot=-1.7, verbose=True):
    '''
    Get the emission lines

    Parameters
    ----------
    lu : floats
       U of the galaxies per component.
    lnH : floats
       ne of the galaxies per component (cm^-3).
    lzgas : floats
       Metallicity of the galaxies per component (log10(Z))
    photomod : string
       Name of the considered photoionisation model.
    xid_phot : float
       Dust-to-metal ratio for the photoionisation model
    co_phot : float
       C/O ratio  for the photoionisation model
    imf_cut_phot : float
       Solar mass high limit for the IMF  for the photoionisation model
    alpha_phot : float
       Alpha value for the AGN photoionisation model.
    verbose : boolean
       If True print out messages

    Returns
    -------
    nebline : array of floats
        Line luminosity per galaxy component, if relevant.
        Units depend on the photoionisation model.
    '''

    if photmod not in c.photmods:
        if verbose:
            print('STOP (gne_photio.get_lines): Unrecognised model to get emission lines.')
            print('                Possible photmod= {}'.format(c.photmods))
        sys.exit()
    elif (photmod == 'gutkin16'):
        nebline = get_lines_gutkin16(lu,lnH,lzgas,xid_phot=xid_phot,
                                     co_phot=co_phot,imf_cut_phot=imf_cut_phot,
                                     verbose=verbose)
    elif (photmod == 'feltre16'):
        nebline = get_lines_feltre16(lu,lnH,lzgas,xid_phot=xid_phot,
                                   alpha_phot=alpha_phot,verbose=verbose)

    return nebline

