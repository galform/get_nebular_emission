import os
import numpy as np
import gne.gne_cosmology as cosmo
import gne.gne_stats as st
import gne.gne_const as c
import gne.gne_io as io

def lines_BPT(x, BPT, line):
    '''
    
    Boundary lines for the distinction of ELG types in BPT diagrams.
    It assummes OIII/Hb on the y axis.
 
    Parameters
    ----------
    
    x : floats
       Array of points on the x axis to define the lines. 
       It should correspond to the wanted BPT.
    BPT : string
       Key corresponding to the wanted x axis for the BPT.
    line : string
       Key corresponding to the wanted boundary line for the BPT.
    
    Returns
    -------
    boundary : floats
    Values of the boundary line in the desired range.
    '''

    boundary = np.zeros(len(x)); boundary.fill(-999.)
    
    if BPT=='NII':
        if line=='Kauffmann2003':
            x0 = 0.05
            boundary[x<x0] = 0.61/(x[x<x0] - x0) + 1.3
        elif line=='Kewley2001':
            x0 = 0.47
            boundary[x<x0] = 0.61/(x[x<x0] - x0) + 1.19
        elif line=='LINER_NIIlim':
            boundary = np.log10(0.6) # Kauffmann 2003
        elif line=='LINER_OIIIlim':
            boundary = np.log10(3) # Kauffmann 2003
    elif BPT=='SII':
        if line=='Kewley2001':
            x0 = 0.32
            boundary[x<x0] = 0.72/(x[x<x0] - x0) + 1.3
        elif line=='Kewley2006':
            boundary = 1.89*x + 0.76
    else:
        print('STOP (gne_plots.lines_BPT): ',
              'BPT plot not recognized.')
        return None
            
    return boundary


def get_obs_bpt(redshift,bpt):
    '''
    Get observational data for BPT diagrams at a given redshift
    
    Parameters
    ----------
    redshift : float
       Redshift of interest
    bpt: string
        Type of BPT diagram: 'NII' (OIII/Hbeta vs N2/Ha)
        or 'SII' (OIII/Hbeta vs S2/Ha)

    Returns
    -------
    xobs, yobs : array of floats
       Ratios for each observed spectral emission line
    obsdata : boolean
       True if there is any observational data at the given redshift
    '''

    xobs = -999.; yobs = -999.; obsdata = False

    # Use different data sets for different redshifts
    if redshift <= 0.2:
        obsdata = True
        obsfile = os.path.join(c.obs_data_dir,'favole2024.txt')
        l1,l2 = np.loadtxt(obsfile,skiprows=1,usecols=(15,9),unpack=True)
        xx, yy = [np.zeros(len(l1)) for i in range(2)]
        ind = np.where((l1>0.) & (l2>0.))
        if (np.shape(ind)[1]>0): #O3/Hb
            yy[ind] = np.log10(l1[ind]/l2[ind])
            
        if bpt=='NII': #N2/Ha
            l1,l2 = np.loadtxt(obsfile,skiprows=1,usecols=(18,6),unpack=True)
            ind = np.where((l1>0.) & (l2>0.))
            if (np.shape(ind)[1]>0):
                xx[ind] = np.log10(l1[ind]/l2[ind]) 
        elif bpt=='SII': #S2/Ha
            l1,l2 = np.loadtxt(obsfile,skiprows=1,usecols=(21,6),unpack=True)
            ind = np.where((l1>0.) & (l2>0.))
            if (np.shape(ind)[1]>0):
                xx[ind] = np.log10(l1[ind]/l2[ind]) 

    elif 1.45 <= redshift <= 1.75:
        obsdata = True
        if bpt=='NII':
            obsfile = os.path.join(c.obs_data_dir,'NII_Kashino.txt')
            yy = np.loadtxt(obsfile,skiprows=18,usecols=(6)) #O3/Hb
            xx = np.loadtxt(obsfile,skiprows=18,usecols=(3)) #N2/Ha
                
        elif bpt=='SII':
            obsfile = os.path.join(c.obs_data_dir,'SII_Kashino.txt')
            yy = np.loadtxt(obsfile,skiprows=18,usecols=(6)) #O3/Hb
            xx = np.loadtxt(obsfile,skiprows=18,usecols=(3)) #N2/Ha

    if obsdata:
        ind = np.where((xx>c.notnum) & (yy>c.notnum))
        if (np.shape(ind)[1]>0):
            xobs = xx[ind]
            yobs = yy[ind]
        else:
            obsdata = False

    return xobs,yobs,obsdata


def get_pozzetti(metadata=None,outpath=None,verbose=False):
    xobs = -999.; yobs = -999.; pozzetti_mod = False
    nomtab = 'pozzetti_tab3.txt'
    pozzetti_table = os.path.join(c.obs_data_dir,nomtab)
    
    # Read table 3 form Pozzetti+2018
    data = np.loadtxt(pozzetti_table)
    zmin = data[:,0]
    zmax = data[:,1]
    nzz  = len(zmin)

    mods = data[:,2:]
    nmod = np.shape(mods)[1]

    # Find if redshift within the intervals
    redshift = metadata['redshift']
    edges = np.insert(zmax,0,zmin[0])
    izz = st.locate_interval(redshift,edges,side='left')
    if (izz<0) or (izz>nzz-1):
        return xobs,yobs,pozzetti_mod
    pozzetti_mod = True

    # Name of Pozzetti's table in Mpc for a cosmology
    omega0 = metadata['omega0']
    omegab = metadata['omegab']
    lambda0 = metadata['lambda0']
    h0 = metadata['h0']

    om = str(omega0)
    ob = str(omegab)
    l0 = str(lambda0)
    hh = str(h0)
    cosmo_str = ('_m'+om+'_b'+ob+'_l'+l0+'_h'+hh).replace('.','p')
    nomtabMpc = 'pozzettiMpc'+cosmo_str+'.txt'

    # Check if Pozzeti's table in Mpc exists for this cosmology
    if outpath is None:
        opath = os.path.join(c.repo_dir, 'output')
    else:
        opath = outpath
    pozzettiMpc = os.path.join(opath,nomtabMpc)
    file_exists = io.check_file(pozzettiMpc)

    # Generate the table in Mpc if needed
    if (not file_exists):
        cosmo.set_cosmology(omega0=omega0, omegab=omegab,
                            lambda0=lambda0,h0=h0)
        for j in range(nmod):
            for i in range(nzz):
                if (mods[i,j]>0):
                    n_Mpch = cosmo.ndeg2nV(mods[i,j],zmin[i],zmax[i])
                    data[i,j+2] = n_Mpch*h0*h0*h0 #Mpc-3

        # Read header lines from original table
        with open(pozzetti_table, 'r') as fin:
            header_lines = [line for line in fin if line.startswith('#')]
            header_lines = [line.replace('deg-2', 'Mpc-3') for line in header_lines]

        with open(pozzettiMpc,'w') as f:
            # Write cosmology info
            header = f'# omega0={omega0}, omegab={omegab},lambda0={lambda0}, h0={h0}'
            f.write(header + '\n') 
            # Write table info
            for line in header_lines:
                f.write(line)
            # Write out the converted values
            np.savetxt(f, data, fmt=['%.1f']*2+['%.5e']*nmod)  
        if verbose:
            print('Table 3 in Pozzetti+2016 converted to (Mpc/h)^3:\n',pozzettiMpc)

    # Flux limits, units of 1e-16 erg cm-2 s-1
    xx = np.array([0.5, 1, 2, 3, 5])
    xx *= 1e-16

    # Read the data for model 3    
    data = np.loadtxt(pozzettiMpc)
    yy = data[izz,12:] # Model 3 (for 1)2:7, 2)7:12)

    ind = np.where(yy>0)
    if np.shape(ind)[1]>1:
        xobs = np.log10(xx[ind])
        yobs = np.log10(yy[ind])
    else:
        pozzetti_mod = False
    
    return xobs,yobs,pozzetti_mod
