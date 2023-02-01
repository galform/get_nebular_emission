import h5py
import numpy as np
from get_nebular_emission.eml_io import get_nheader, homedir, locate_interval
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_io import check_file
import sys
import warnings
from cosmology import emission_line_flux

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

def coef_att_galform(infile,cols,inputformat='HDF5',photmod='gutkin16',verbose=True):
    '''
    It reads luminosities of lines with and without attenuation
    from GALFORM data and it returns the attenuation coefficients.
    It supposes two components.

    Parameters
    ----------
    infile : string
     Name of the input file. Expects an hdf5 file.
    inputformat : string
     Format of the input file.
    cols : list
     - [[component1_stellar_mass,sfr,Z],[component2_stellar_mass,sfr,Z],...]
     - Expected : component1 = total or disk, component2 = bulge
     - For text or csv files: list of integers with column position.
    photmod : string
     Photoionisation model to be used for look up tables.
    verbose : boolean
     If True print out messages

    Returns
    -------
    coef_att : floats
    '''
    
    check_file(infile, verbose=verbose)
    
    ncomp = len(cols)
    
    if ncomp>2:
        if verbose:
            print('STOP (eml_photio.coef_att_galform): ',
                  'GALFORM does not have more than two components.')
        sys.exit()
        
    lines = const.lines_model[photmod]
    numlines = len(lines)
    
    if inputformat=='HDF5':
        with h5py.File(infile, 'r') as f:
            hf = f['data']
            
            coef_att = np.empty((ncomp,numlines,len(hf[const.line_headers[0] + lines[0]])))
            coef_att.fill(1)
        
            for i in range(numlines):
                
                if ncomp==2:
                    line_disk = const.line_headers[1] + lines[i]
                    line_bulge = const.line_headers[2] + lines[i]
                    att_disk = line_disk + const.att_ext
                    att_bulge = line_bulge + const.att_ext
                    try:
                        indd = np.where(hf[line_disk][:] > 0.)[0] 
                        indb = np.where(hf[line_bulge][:] > 0.)[0]
                        
                        coef_att[0,i,indd] = hf[att_disk][indd]/hf[line_disk][indd]
                        coef_att[1,i,indb] = hf[att_bulge][indb]/hf[line_bulge][indb]
                        
                        coef_att[0,i][(coef_att[0,i]>1.)&(coef_att[0,i]-1 < 0.01)] = 1.
                        coef_att[1,i][(coef_att[1,i]>1.)&(coef_att[1,i]-1 < 0.01)] = 1.
                    except:
                        coef_att[0,i] = const.notnum
                        coef_att[1,i] = const.notnum
                else:
                    line = const.line_headers[0] + lines[i]
                    att = line + const.att_ext
                    try:
                        ind = np.where(hf[line][:] > 0.)[0]
                        
                        coef_att[0,i,ind] = hf[att][ind]/hf[line][ind]
                        
                        coef_att[0,i][(coef_att[0,i]>1.)&(coef_att[0,i]-1 < 0.01)] = 1.
                    except:
                        coef_att[0,i] = const.notnum
    elif inputformat=='textfile':
        headers = ['mag_LC_r_disk', 'mag_LC_r_bulge', 'zcold', 'mcold', 
                   'zcold_burst', 'mcold_burst', 'mstardot_average', 'L_tot_Halpha',
                   'L_tot_NII6583', 'L_tot_Hbeta', 'L_tot_OIII5007', 'mstars_total', 'is_central',
                   'mstardot', 'mstardot_burst', 'mstars_bulge', 'L_tot_OII3727', 'L_tot_SII6716',
                   'L_tot_SII6731', 'mag_SDSS_r_o_t', 'L_tot_Halpha_ext', 'L_tot_Hbeta_ext', 
                   'L_tot_OII3727_ext', 'L_tot_OIII5007_ext', 'L_disk_Halpha', 
                   'L_disk_Halpha_ext', 'L_bulge_Halpha', 'L_bulge_Halpha_ext', 
                   'L_disk_Hbeta', 'L_disk_Hbeta_ext', 'L_bulge_Hbeta', 'L_bulge_Hbeta_ext',
                   'L_disk_OIII5007', 'L_disk_OIII5007_ext', 'L_bulge_OIII5007', 'L_bulge_OIII5007_ext', 
                   'L_disk_NII6583', 'L_disk_NII6583_ext', 'L_bulge_NII6583', 'L_bulge_NII6583_ext', 
                   'L_disk_OII3727', 'L_disk_OII3727_ext', 'L_bulge_OII3727', 'L_bulge_OII3727_ext', 
                   'L_disk_SII6717', 'L_disk_SII6717_ext', 'L_bulge_SII6717', 'L_bulge_SII6717_ext', 
                   'L_disk_SII6731', 'L_disk_SII6731_ext', 'L_bulge_SII6731', 'L_bulge_SII6731_ext']
        
        ih = get_nheader(infile)        
        X = np.loadtxt(infile,skiprows=ih).T
        
        coef_att = np.empty((ncomp,numlines,len(X[0])))
        coef_att.fill(1)
    
        for i in range(numlines):
            
            if ncomp==2:
                line_disk = const.line_headers[1] + lines[i]
                line_bulge = const.line_headers[2] + lines[i]
                att_disk = line_disk + const.att_ext
                att_bulge = line_bulge + const.att_ext
                
                for j, elem in enumerate(headers):
                    if elem == line_disk:
                        line_disk = j
                    elif elem == line_bulge:
                        line_bulge = j
                    elif elem == att_disk:
                        att_disk = j
                    elif elem == att_bulge:
                        att_bulge = j
                try:
                    indd = np.where(X[line_disk] > 0.)[0]
                    indb = np.where(X[line_bulge] > 0.)[0]
                    
                    coef_att[0,i,indd] = X[att_disk][indd]/X[line_disk][indd]
                    coef_att[1,i,indb] = X[att_bulge][indb]/X[line_bulge][indb]
                    
                    coef_att[0,i][(coef_att[0,i]>1.)&(coef_att[0,i]-1 < 0.01)] = 1.
                    coef_att[1,i][(coef_att[1,i]>1.)&(coef_att[1,i]-1 < 0.01)] = 1.
                except:
                    coef_att[0,i] = const.notnum
                    coef_att[1,i] = const.notnum
            else:
                line = const.line_headers[0] + lines[i]
                att = line + const.att_ext
                
                for j, elem in enumerate(headers):
                    if elem == line:
                        line = j
                    elif elem == att:
                        att = j
                try:
                    ind = np.where(X[line] > 0.)[0]
                    
                    coef_att[0,i,ind] = X[att][ind]/X[line][ind]
                    
                    coef_att[0,i][(coef_att[0,i]>1.)&(coef_att[0,i]-1 < 0.01)] = 1.
                except:
                    coef_att[0,i] = const.notnum
                    
        del X
        
    
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

def clean_photarray(lms, lssfr, lu, lne, loh12, cutlimits=False, photmod='gutkin16', verbose=True):

    '''
    Given the model, take the values outside the limits and give them the apropriate
    value inside the limits depending on the model.

    Parameters
    ----------
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lu : floats
     U of the galaxies per component.
    lne : floats
     ne of the galaxies per component (cm^-3).
    loh12 : floats
     Metallicity of the galaxies per component (12+log(O/H))
    cutlimits : boolean
     If True the galaxies with U, ne and Z outside the photoionization model's grid limits won't be considered.
    photomod : string
     Name of the considered photoionisation model
    verbose : boolean
     If True print out messages

    Returns
    -------
    lms,lssfr,lu,lne,loh12 : floats
    limits : integers
    '''

    minU, maxU = get_limits(propname='U', photmod=photmod)
    minnH, maxnH = get_limits(propname='nH', photmod=photmod)
    minZ, maxZ = get_limits(propname='Z', photmod=photmod)
    
    limits = np.where((lu[:,0]>minU)&(lu[:,0]<maxU)&(loh12[:,0]>np.log10(minZ))&(loh12[:,0]<np.log10(maxZ))&
                      (lne[:,0]>np.log10(minnH))&(lne[:,0]<np.log10(maxnH))&(lu[:,0]!=const.notnum))[0]
    
    if cutlimits:
        lms = lms[limits]
        lssfr = lssfr[limits]
        loh12 = loh12[limits]
        lu = lu[limits]
        lne = lne[limits]
    else:
        for i in range(lu.shape[1]):        
                lu[:,i][(lu[:,i] > maxU)&(lu[:,i] != const.notnum)] = maxU
                lu[:,i][(lu[:,i] < minU)&(lu[:,i] != const.notnum)] = minU
                
                lne[:,i][(lne[:,i] > np.log10(maxnH))&(lne[:,i] != const.notnum)] = np.log10(maxnH)
                lne[:,i][(lne[:,i] < np.log10(minnH))&(lne[:,i] != const.notnum)] = np.log10(minnH)
                
                loh12[:,i][(loh12[:,i] > np.log10(maxZ))&(loh12[:,i] != const.notnum)] = np.log10(maxZ)
                loh12[:,i][(loh12[:,i] < np.log10(minZ))&(loh12[:,i] != const.notnum)] = np.log10(minZ)
                
    return lms, lssfr, lu, lne, loh12, limits

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
        infile = const.mod_lim[photmod]
    except KeyError:
        print('STOP (eml_photio): the {}'.format(photmod) + ' model is an unrecognised model in the dictionary mod_lim')
        print('                  Possible photmod= {}'.format(const.mod_lim.keys()))
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
     U of the galaxies per component.
    lne : floats
     ne of the galaxies per component (cm^-3).
    loh12 : floats
     Metallicity of the galaxies per component (12+log(O/H))
    Plotting : boolean
      If True run verification plots with all data.
    Testing : boolean
      If True to only run over few entries for testing purposes
    verbose : boolean
      If True print out messages
      
    Returns
    -------
    nebline : floats
     Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
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
        
        ind = np.where(lu[:,comp] != const.notnum)[0]

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
            for ii in ind:
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
            for ii in ind:
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
        for n in ind:
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
     U of the galaxies per component.
    lne : floats
     ne of the galaxies per component (cm^-3).
    loh12 : floats
     Metallicity of the galaxies per component (12+log(O/H))
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
    nebline : floats
     Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
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

def attenuation(nebline, m_sfr_z, att_param=['Rvir','ColdGas'], inputformat='HDF5', attmod='GALFORM', 
                photmod='gutkin16', infile=None, cut=None, limits=None, cutlimits=False, verbose=True):
    '''
    Get the attenuated emission lines from the raw ones.

    Parameters
    ----------
    nebline : floats
     Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr).
    infile : string
     - Name of the input file. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    m_sfr_z : list
     - [[component1_stellar_mass,sfr/LC,Z],[component2_stellar_mass,sfr/LC,Z],...]
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    att_param : list
     Parameters to look for calculating attenuation. See eml_const to know what each model expects.
     - For text or csv files: list of integers with column position.
     - For hdf5 files: list of data names.
    inputformat : string
     Format of the input file.
    attmod : string
     Attenuation model.
    photmod : string
     Photoionisation model to be used for look up tables.
    cut : integers
     Indeces of the cutted galaxies.
    limits : integers
     Indeces of the galaxies with U, ne and Z outside the limits of the photoionization model
     after the initial cut.
    cutlimits : boolean
     If True the galaxies with U, ne and Z outside the photoionization model's grid limits won't be considered.
    verbose : boolean
     If True print out messages.

    Returns
    -------
    nebline_att : floats
      Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
    lines : floats
      Array with indexes of the lines whose fluxes will be stored in the output file.
      
    Notes
    -------
    It will ignore the lines for which, for any reason, attenuation can not be calculated.
    '''
    
    ncomp = len(nebline)
    
    if attmod not in const.attmods:
        if verbose:
            print('STOP (eml_photio.attenuation): Unrecognised model for attenuation.')
            print('                Possible attmod= {}'.format(const.attmods))
        sys.exit()
    elif attmod=='ratios':
        coef_att = coef_att_galform(infile,m_sfr_z,inputformat=inputformat,photmod=photmod,verbose=verbose)
        
        coef_att = coef_att[:,:,cut]        
        if cutlimits:
            coef_att = coef_att[:,:,limits]
    elif attmod=='cardelli89':
        with h5py.File(infile,'r') as f:
            Rvir = f['data'][att_param[0]][:]
            Mcold_disc = f['data'][att_param[1]][:]
            Z_disc = f['data'][m_sfr_z[0][2]][:]
            
            Rvir = Rvir[cut]
            Mcold_disc = Mcold_disc[cut]
            Z_disc = Z_disc[cut]
            
            if cutlimits: 
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
        
    nebline_att = nebline*coef_att
    
    
    return nebline_att, coef_att

