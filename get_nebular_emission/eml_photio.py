import h5py
import numpy as np
from get_nebular_emission.eml_io import get_nheader, homedir, locate_interval
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_io import check_file
import sys
import warnings


mod_lim = {'gutkin16': r"nebular_data/gutkin_tables/limits_gutkin.txt"}

# print(mod_lim.keys()) # To show the keys

#------------------------------------------------------------------------------------
#   Cardelli et al. 1989 extinction laws in FIR and IR/OPT:
#------------------------------------------------------------------------------------
def cardelli(waveA):
    Rv=3.1 
    wave=waveA/10000.
    x=1./wave
    
    if (x < 0.3) or (x > 10):
        print('STOP (eml_photio.cardelli): ',
              'Wavelength out of range.')
        sys.exit()
        return
    elif (x < 1.1): #IR
        ax = 0.574*(x**1.61) 
        bx = -0.527*(x**1.61)
    elif (x < 3.3): #Optical/NIR
        y = x-1.82
        ax = (1.+0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 
        0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7) 
        bx = (1.41338*y+2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 -
        0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7)
    elif (x < 8): #UV
        if (x < 5.9):
            Fa = 0
            Fb = 0
        else: 
            Fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
            Fb = 0.2130*(x-5.9)**2 + 0.1207*(x-5.9)**3
        ax = 1.752-0.316*x - 0.104/((x-4.67)**2 + 0.341) + Fa
        bx = -3.090+1.825*x + 1.206/((x-4.62)**2 + 0.263) + Fb
    else:
        ax = -1.073 - 0.628*(x-8) + 0.137*(x-8)**2 - 0.070*(x-8)**3
        bx = 13.670 + 4.257*(x-8) - 0.420*(x-8)**2 + 0.374*(x-8)**3
        
    return ax+bx/Rv

def coef_att_cardelli(wavelength, Mcold_disc, rhalf_mass_disc, Z_disc, costheta=0.3, albedo=0.56):
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    
    a_disc = 1.68
    if wavelength > 2000: #Guiderdoni & Roca-Volmerange 1987 (see De Lucia & Blaizot 2007)
        s = 1.6
    else:
        s = 1.35
    Al_Av = cardelli(wavelength)
    sectheta = 1./costheta

    mean_col_dens_disc_log = (np.log10(Mcold_disc) + np.log10(const.Msun_to_kg) - 
    np.log10(1.4*const.mp*np.pi)-2.*np.log10(a_disc*rhalf_mass_disc*const.Mpc_to_cm))
    
    tau_disc = np.log10(Al_Av) + np.log10((Z_disc/const.zsun)**s) + mean_col_dens_disc_log - np.log10(2.1e21)
    tau_disc = 10.**tau_disc
    
    al_disc = (np.sqrt(1.-albedo))*tau_disc
    
    A_lambda = -2.5*np.log10((1.-np.exp(-al_disc*sectheta))/(al_disc*sectheta))
    
    coef_att = 10.**(-0.4*A_lambda)
    
    return coef_att

def coef_att_galform(infile,cols,lines=[0,1,3,5,6,7,8],photmod='gutkin16',verbose=True):
    '''
    It reads luminosities of lines with and without attenuation
    from GALFORM data and it returns the attenuation coefficients.
    It supposes two components.

    Parameters
    ----------
    infile : string
      Name of the input file. Expects an hdf5 file.
    verbose : boolean
      If True print out messages

    Returns
    -------
    coef_att : floats
    '''
    
    check_file(infile, verbose=verbose)
    
    ncomp = len(cols)
    
    if ('.hdf5' not in infile):
        infile = infile[:-4] + '.hdf5' #This expects a .txt or a .dat
    
    lines = const.lines_model[photmod][lines]
    
    numlines = len(lines)
    
    with h5py.File(infile, 'r') as f:
        hf = f['data']
        
        coef_att = np.empty((ncomp,numlines,len(hf['L_disk_' + lines[0]])))
        coef_att.fill(1)
    
        if ncomp==2:
            for i in range(numlines):
                indd = np.where(hf['L_disk_' + lines[i]][:] > 0.)[0]
                indb = np.where(hf['L_bulge_' + lines[i]][:] > 0.)[0]
                
                coef_att[0,i,indd] = hf['L_disk_' + lines[i] + '_ext'][indd]/hf['L_disk_' + lines[i]][indd]
                coef_att[1,i,indb] = hf['L_bulge_' + lines[i] + '_ext'][indb]/hf['L_bulge_' + lines[i]][indb]
                
                coef_att[0,i][(coef_att[0,i]>1.)&(coef_att[0,i]-1 < 0.01)] = 1.
                coef_att[1,i][(coef_att[1,i]>1.)&(coef_att[1,i]-1 < 0.01)] = 1.
        elif ncomp==1:
            for i in range(numlines):
                ind = np.where(hf['L_tot_' + lines[i]][:] > 0.)[0]
                
                coef_att[0,i,ind] = hf['L_tot_' + lines[i] + '_ext'][ind]/hf['L_tot_' + lines[i]][ind]
                
                coef_att[0,i][(coef_att[0,i]>1.)&(coef_att[0,i]-1 < 0.01)] = 1.
        else:
            if verbose:
                print('STOP (eml_photio.coef_att_galform): ',
                      'GALFORM does not have more than two components.')
                sys.exit()
    
    return coef_att


def get_zfile(zmet_str, photmod='gutkin16'):

    '''
    Given a metallicity string get the name of the corresponding table

    Parameters
    ----------
    zmet_str : string
        Metallicity name in files
    photomod : string
        Name of the considered photoionisation model

    Returns
    -------
    zfile : string
        Name of the model file with data for the given metallicity
    '''

    root = 'nebular_data/' + photmod + '_tables/nebular_emission_Z'
    if len(zmet_str)<3:
        zmet_str = zmet_str+'0'
    zfile = root + zmet_str + '.txt'
    # Si son 2 numeros que le añada un cero
    file_fine = check_file(zfile)
    if (not file_fine):
        zfile = None

    return zfile

def clean_photarray(lu, lne, loh12, photmod='gutkin16', verbose=True):

    '''
    Given the model, take the values outside the limits and give them the apropriate
    value inside the limits depending on the model.

    Parameters
    ----------
    photomod : string
     Name of the considered photoionisation model

    verbose : boolean
     If True print out messages

    Returns
    -------
    lu,lne,loh12 : arrays
     - Array of the properties with all the data in the limits.
     - Two components [disk, bulge]
    '''


    # Read the data file:
    # infile = r"output_data/U_ne_loh12.hdf5"

    # check_file(infile, verbose=verbose)

    # f = h5py.File(infile,'r')
    # header = f['header']
    # data = f['data']

    # lu = data['lu'][:]
    lud=[]; lub=[]

    for ii, u in enumerate(lu):
        u1 = lu[ii]
        u1 = u1.tolist()
        lud.append(u1[0])
        lub.append((u1[1]))
    lud = np.array(lud)
    lub = np.array(lub)


    # lne = data['lne'][:]
    lned=[];lneb=[]
    for ii, ne in enumerate(lne):
        ne1 = lne[ii]
        ne1 = ne1.tolist()
        lned.append(ne1[0])
        lneb.append((ne1[1]))
    lned = np.array(lned)
    lneb = np.array(lneb)


    # loh12 = data['loh12'][:]
    loh12d=[];loh12b=[]
    for ii, oh12 in enumerate(loh12):
        oh121 = loh12[ii]
        oh121 = oh121.tolist()
        loh12d.append(oh121[0])
        loh12b.append((oh121[1]))
    loh12d = np.array(loh12d)
    loh12b = np.array(loh12b)


    lowerl, upperl = get_limits(propname='U', photmod=photmod)

    ind = np.where((lud > upperl)&(lud != const.notnum))
    lud[ind] = upperl
    ind = np.where((lud < lowerl)&(lud != const.notnum))
    lud[ind] = lowerl

    ind = np.where((lub > upperl)&(lub != const.notnum))
    lub[ind] = upperl
    ind = np.where((lub < lowerl)&(lub != const.notnum))
    lub[ind] = lowerl

    lu = np.stack((lud, lub))

    lowerl, upperl = get_limits(propname='nH', photmod=photmod)
    
    ind = np.where((lned > upperl)&(lned != const.notnum))
    lned[ind] = upperl
    ind = np.where((lned < lowerl)&(lned != const.notnum))
    lned[ind] = lowerl

    ind = np.where((lneb > upperl)&(lneb != const.notnum))
    lneb[ind] = upperl
    ind = np.where((lneb < lowerl)&(lneb != const.notnum))
    lneb[ind] = lowerl

    lne = np.stack((lned, lneb))

    lowerl, upperl = get_limits(propname='Z', photmod=photmod)
    
    ind = np.where((loh12d > upperl)&(loh12d != const.notnum))
    loh12d[ind] = upperl
    ind = np.where((loh12d < lowerl)&(loh12d != const.notnum))
    loh12d[ind] = lowerl

    ind = np.where((loh12b > upperl)&(loh12b != const.notnum))
    loh12b[ind] = upperl
    ind = np.where((loh12b < lowerl)&(loh12b != const.notnum))
    loh12b[ind] = lowerl

    loh12 = np.stack((loh12d, loh12b))

    return lu.T, lne.T, loh12.T

def get_limits(propname, photmod = 'gutkin16',verbose=True):
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
    >>> # Table 3 of Gutkin+2016 (https://arxiv.org/pdf/1607.06086.pdf)
    >>> # Property      Lower_limit     Upper_limit
    >>> Z               0.0001          0.040
    >>> U               -4.0            -1.0
    >>> xid             0.1             0.5
    >>> nH               1               4 
    >>> get_limits(propname = 'Z', photmod = 'gutkin16')
        0.0001 0.04
    >>> get_limits(propname = 'nH', photmod = 'gutkin16')
        1  4

    '''

    try:
        infile = mod_lim[photmod]
    except KeyError:
        print('STOP (eml_photio): the {}'.format(photmod) + ' model is an unrecognised model in the dictionary mod_lim')
        print('                  Possible photmod= {}'.format(mod_lim.keys()))
        exit()

    # Check if the limits file exists:
    check_file(infile, verbose=verbose)
    # print(infile)

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

def get_lines_Gutkin(lu, lne, loh12, Testing=False, Plotting=False, verbose=True):
    '''
    Get the interpolations for the emission lines,
    using the tables
    from Gutkin et al. (2016) (https://arxiv.org/pdf/1607.06086.pdf)

    Parameters
    ----------
    lu : floats
      Ionization parameter, log10(U) (dimensionless) [disk, bulge]
    lne : floats
      Electron density, log10(nH) (cm**-3) [disk, bulge]
    loh12 : floats
      Metallicity, 12+log(O/H) [disk, bulge]
    Plotting : boolean
      If True run verification plots with all data.
    Testing : boolean
      If True to only run over few entries for testing purposes
    verbose : boolean
      If True print out messages
      
    Returns
    -------
    nebline_disk, nebline_bulge : floats
      Fluxes of the lines for the disk and bulge components.
    '''
    
    zmet_str = const.zmet_str['gutkin16']
    zmets = np.full(len(zmet_str),const.notnum)
    zmets = np.array([float('0.' + zmet) for zmet in zmet_str])

    logubins = [-4., -3.5, -3., -2.5, -2., -1.5, -1.]
    
    nemline = 18
    ndat = lu.shape[0]
    ncomp = lu.shape[1]

    nzmet = 14
    nu = 7
    nzmet_reduced = 4
    zmets_reduced = const.zmet_reduced['gutkin16']

    emline_grid1 = np.zeros((nzmet_reduced,nu,nemline)) # From slower to faster
    emline_grid2 = np.zeros((nzmet,nu,nemline))
    emline_grid3 = np.zeros((nzmet,nu,nemline))
    emline_grid4 = np.zeros((nzmet_reduced,nu,nemline))

    l = 0
    kred = 0
    nn = 0

    for k, zname in enumerate(zmets):
        infile = get_zfile(zmet_str[k],photmod='gutkin')
        check_file(infile,verbose=True)
        #print(k,infile)
        ih = get_nheader(infile)

        with open(infile,'r') as ff:
            iline = -1.
            for line in ff:
                iline += 1

                if iline<ih:continue

                data = np.array((line.split()))
                u = float(data[0])
                xid = float(data[1])
                nH = float(data[2])
                co = float(data[3])
                imf_cut = float(data[4])

                if xid==0.3 and co==1.and imf_cut==100:
                    if u == -4.:
                        l = 0
                    if u == -3.5:
                        l = 1
                    if u == -3.:
                        l = 2
                    if u == -2.5:
                        l = 3
                    if u == -2.:
                        l = 4
                    if u == -1.5:
                        l = 5
                    if u == -1.:
                        l = 6


                    if nH==10 or nH==100 or nH==1000 or nH==10000:

                        if nH==10 or nH==10000:
                            if k==0:
                                kred = 0
                            if k==4:
                                kred = 1
                            if k==9:
                                kred = 2
                            if k==12:
                                kred = 3
                        for j in range(nemline):
                            if nH == 10:
                                emline_grid1[kred,l,j] = float(data[j+5])
                            if nH == 100:
                                emline_grid2[k,l,j] = float(data[j+5])
                            if nH == 1000:
                                emline_grid3[k,l,j] = float(data[j+5])
                            if nH == 10000:
                                emline_grid4[kred,l,j] = float(data[j+5])
        ff.close()
    # From Carlton's code, I do not understand from line 201 to line 212:
    '''
    for kred in range(nzmet_reduced):
        for l in range(nu):
            emline_grid1(kred,l,1)
            emline_grid1(kred,l,3)
            emline_grid1(kred,l,5)
            emline_grid1(kred,l,6)
    
    for k in range(nzmet):
        for l in range(nu):
            emline_grid2(k,l,1)
            emline_grid2(k,l,3)
            emline_grid2(k,l,5)
            emline_grid2(k,l,6)
        
            emline_grid3(k,l,1)
            emline_grid3(k,l,3)
            emline_grid3(k,l,5)
            emline_grid3(k,l,6)
    '''

    # log metallicity bins ready for interpolation:

    lzmets_reduced = np.full(len(zmets_reduced), const.notnum)
    ind = np.where(zmets_reduced > 0.)
    if (np.shape(ind)[1]) > 0:
        lzmets_reduced[ind] = np.log10(zmets_reduced[ind])


    lzmets = np.full(len(zmets), const.notnum)
    ind = np.where(zmets > 0.)
    if (np.shape(ind)[1] > 0):
        lzmets[ind] = np.log10(zmets[ind])

    nebline = np.zeros((ncomp,nemline,ndat))

    # Interpolate in all three ne grids to start with u-grid first, since the same for all grids
    
    for comp in range(ncomp):

        emline_int1 = np.zeros((nemline,ndat))
        emline_int2 = np.zeros((nemline, ndat))
        emline_int3 = np.zeros((nemline, ndat))
        emline_int4 = np.zeros((nemline, ndat))
    
        # Interpolate over ionisation parameter
        du = []
        j = []
        for logu in lu[:,comp]:
            j1 = locate_interval(logu,logubins)
            if j1 == 0:
                du.append(0.0)
                j.append(0)
                #du = 0.0
                j1 = 0
            elif j1 == nu - 1:
                du.append(1.0)
                j.append(nu-2)
                #du = 1.0
                j1 = nu - 2
            else:
                d = (logu - logubins[j1]) / (logubins[j1 + 1] - logubins[j1])
                du.append(d)
                j.append(j1)
    
        # Interpolate over disk gas metallicity loh12[comp]
        dz = []
        i = []
        for logz in loh12[:,comp]:
            i1 = locate_interval(logz,lzmets_reduced)
    
            if i1==0:
                dz.append(0.0)
                #dz = 0.0
                i1 = 0
                i.append(0)
            elif i1 == nzmet_reduced-1:
                dz.append(1.0)
                #dz = 1.0
                i1 = nzmet_reduced-2
                i.append(nzmet_reduced-2)
            else:
                d = (logz - lzmets_reduced[i1])/(lzmets_reduced[i1+1]-lzmets_reduced[i1])
                dz.append(d)
                i.append(i1)
    
    
        for k in range(nemline):
            for ii in range(ndat):
                #emline_grid1 = np.zeros((nzmet_reduced, nu, nemline))
                #print(emline_grid1[i[ii]][j[ii]][k])
                emline_int1[k][ii] = (1.-dz[ii])*(1.-du[ii])*emline_grid1[i[ii]][j[ii]][k]+\
                                     dz[ii]*(1-du[ii])*emline_grid1[i[ii]+1][j[ii]][k]+\
                                     (1.-dz[ii])*du[ii]*emline_grid1[i[ii]][j[ii]+1][k]+\
                                     dz[ii]*du[ii]*emline_grid1[i[ii]+1][j[ii]+1][k]
    
                emline_int4[k][ii] = (1.-dz[ii])*(1.-du[ii])*emline_grid4[i[ii]][j[ii]][k]+\
                                     dz[ii]*(1-du[ii])*emline_grid4[i[ii]+1][j[ii]][k]+\
                                     (1.-dz[ii])*du[ii]*emline_grid4[i[ii]][j[ii]+1][k]+\
                                     dz[ii]*du[ii]*emline_grid4[i[ii]+1][j[ii]+1][k]
    
    
        # full metallicity grid for emlines_grid2 ne=100 and emlines_grid3 ne=1000
    
        dz = []
        i = []
    
        for logz in loh12[:,comp]:
            i1 = locate_interval(logz, lzmets)
            if i1 == 0:
                dz.append(0.0)
                # dz = 0.0
                i1 = 0
                i.append(0)
            elif i1 == nzmet - 1:
                dz.append(1.0)
                # dz = 1.0
                i1 = nzmet - 2
                i.append(nzmet - 2)
            else:
                d = (logz - lzmets[i1]) / (lzmets[i1 + 1] - lzmets[i1])
                dz.append(d)
                i.append(i1)
    
    
        for k in range(nemline):
            for ii in range(ndat):
                emline_int2[k][ii] = (1.-dz[ii])*(1.-du[ii])*emline_grid2[i[ii]][j[ii]][k]+\
                                     dz[ii]*(1-du[ii])*emline_grid2[i[ii]+1][j[ii]][k]+\
                                     (1.-dz[ii])*du[ii]*emline_grid2[i[ii]][j[ii]+1][k]+\
                                     dz[ii]*du[ii]*emline_grid2[i[ii]+1][j[ii]+1][k]
    
                emline_int3[k][ii] = (1.-dz[ii])*(1.-du[ii])*emline_grid3[i[ii]][j[ii]][k]+\
                                     dz[ii]*(1-du[ii])*emline_grid3[i[ii]+1][j[ii]][k]+\
                                     (1.-dz[ii])*du[ii]*emline_grid3[i[ii]][j[ii]+1][k]+\
                                     dz[ii]*du[ii]*emline_grid3[i[ii]+1][j[ii]+1][k]
    
        # Interpolate over ne
        # use gas density in disk logned
        for n in range(ndat):
            if (lne[:,comp][n] > 2. and lne[:,comp][n] <= 3.):
                dn = (lne[:,comp][n] -2.)/(3. - 2.)
                for k in range(nemline):
                    nebline[comp][k][n] = (1.-dn)*emline_int2[k][n] + (dn)*emline_int3[k][n]
    
            elif (lne[:,comp][n] > 1. and lne[:,comp][n] <= 2.):
                dn = (lne[:,comp][n] -1.)/(2. - 1.)
                for k in range(nemline):
                    nebline[comp][k][n] = (1.-dn)*emline_int1[k][n] + (dn)*emline_int2[k][n]
    
            elif (lne[:,comp][n] > 3. and lne[:,comp][n]<=4.):
                dn = (lne[:,comp][n] - 3.)/(4. - 3.)
                for k in range(nemline):
                    nebline[comp][k][n] = (1. - dn) * emline_int3[k][n] + (dn) * emline_int4[k][n]
                print('hay mayor que 3')
    
            elif (lne[:,comp][n] <= 1.):
                for k in range(nemline):
                    nebline[comp][k][n] = emline_int1[k][n]
            elif (lne[:,comp][n] > 4.):
                for k in range(nemline):
                    nebline[comp][k][n] = emline_int4[k][n]
            else:
                print('log(ne)disk out of limits','log(ne)disk = {}'.format(lne[:,comp][n]))

    # check if there is an ongoing bulge with the masses

    return nebline


def get_lines(lu, lne, loh12, photmod='gutkin16', verbose=True, Testing=False, Plotting=False):
    '''
    Get the emission lines

    Parameters
    ----------
    lu : floats
      Ionization parameter, log10(U) (dimensionless) [disk, bulge]
    lne : floats
      Electron density, log10(nH) (cm**-3) [disk, bulge]
    loh12 : floats
      Metallicity, 12+log(O/H) [disk, bulge]
    photomod : string
      Name of the considered photoionisation model
    verbose : boolean
      If True print out messages
    Plotting : boolean
      If True run verification plots with all data.
    Testing : boolean
      If True to only run over few entries for testing purposes

    Returns
    -------
    nebline_disk, nebline_bulge : floats
      Fluxes of the lines for the disk and bulge components.
    '''
    
    #lu,lne,loh12 = clean_photarray(lu,lne,loh12,photmod=photmod,verbose=verbose)

    if photmod not in const.photmods:
        if verbose:
            print('STOP (eml_photio.get_lines): Unrecognised model to get emission lines.')
            print('                Possible photmod= {}'.format(const.photmods))
        sys.exit()
    elif (photmod == 'gutkin16'):
        nebline = get_lines_Gutkin(lu,lne,loh12, Testing=Testing, 
                Plotting=Plotting, verbose=verbose)

    return nebline

def attenuation(nebline, m_sfr_z, att_param=['Rvir','ColdGas'], attmod='GALFORM', 
                photmod='gutkin16', infile=None, cut=None, limits=None, verbose=True):
    '''
    Get the attenuated emission lines from the raw ones.

    Parameters
    ----------
    nebline_disk, nebline_bulge : floats
      Fluxes of the lines for the disk and bulge components.

    Returns
    -------
    nebline_att : floats
      Fluxes of the attenuated lines for the disk and bulge components.
      
    lines : floats
      Array with indexes of the lines whose fluxes will be stored in the output file.
      
    Notes
    -------
    In general, lines will cover all the lines in the model. It will only eliminate the lines
    for which, for any reason, attenuation can not be calculated.
    '''
    
    ncomp = len(nebline)
    
    if ('.hdf5' not in infile):
        infile = infile[:-4] + '.hdf5' #This expects a .txt or a .dat 
        
    if attmod not in const.attmods:
        if verbose:
            print('STOP (eml_photio.attenuation): Unrecognised model for attenuation.')
            print('                Possible attmod= {}'.format(const.attmods))
        sys.exit()
    elif attmod=='ratios':
        lines = []
        with h5py.File(infile,'r') as f:
            hf = f['data']
            for i, line in enumerate(const.lines_model[photmod]):
                if any(line in key for key in hf.keys()):
                    lines.append(i)
        coef_att = coef_att_galform(infile,m_sfr_z, lines=lines,photmod=photmod,verbose=verbose)
        
        coef_att = coef_att[:,:,cut]
        coef_att = coef_att[:,:,limits]
    elif attmod=='cardelli89':
        lines = np.arange(len(nebline[0,:]))
        with h5py.File(infile,'r') as f:
            Rvir = f['data'][att_param[0]][:]
            Mcold_disc = f['data'][att_param[1]][:]
            Z_disc = f['data'][m_sfr_z[0][2]][:]
            
            Rvir = Rvir[cut]
            Mcold_disc = Mcold_disc[cut]
            Z_disc = Z_disc[cut]
            
            Rvir = Rvir[:1153192]
            Mcold_disc = Mcold_disc[:1153192]
            Z_disc = Z_disc[:1153192]
            
            Rvir = Rvir[limits]
            Mcold_disc = Mcold_disc[limits]
            Z_disc = Z_disc[limits]
            
            coef_att = np.empty(nebline.shape)
        for i, line in enumerate(const.lines_model[photmod]):
            for comp in range(ncomp):
                if comp==0:
                    coef_att[comp,i] = coef_att_cardelli(const.wavelength_model[photmod][i], 
                                Mcold_disc=Mcold_disc, rhalf_mass_disc=Rvir*const.h/2, 
                                Z_disc=Z_disc, costheta=0.3, albedo=0.56)
                else:
                    coef_att[comp,i] = coef_att[0,i]
        coef_att[coef_att!=coef_att] = 1.
    else:
        if verbose:
            print('STOP (eml_photio): Unrecognised model to attenuate lines.')
            print('                Possible attmod= {}'.format(const.attmods))
        sys.exit()
        
    nebline_att = nebline[:,lines]*coef_att
    
    
    return nebline_att, coef_att, lines

