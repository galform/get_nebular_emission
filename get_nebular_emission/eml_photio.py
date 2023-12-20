import h5py
import numpy as np
from get_nebular_emission.eml_io import get_nheader, homedir, locate_interval
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_io import check_file
import sys
import warnings
from cosmology import emission_line_flux

from cosmology import logL2flux, set_cosmology

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

    root = 'nebular_data/' + photmod + '_tables/nebular_emission_Z'
    if len(zmet_str)<3:
        zmet_str = zmet_str+'0'
    zfile = root + zmet_str + '.txt'
    # Si son 2 numeros que le añada un cero
    file_fine = check_file(zfile)
    if (not file_fine):
        zfile = None

    return zfile

def clean_photarray(lms, lssfr, lu, lne, loh12, photmod='gutkin16', verbose=True):

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
     Metallicity of the galaxies per component (log10(Z)).
    photomod : string
     Name of the considered photoionisation model.
    verbose : boolean
     If True print out messages.

    Returns
    -------
    lms,lssfr,lu,lne,loh12 : floats
    '''

    minU, maxU = get_limits(propname='U', photmod=photmod)
    minnH, maxnH = get_limits(propname='nH', photmod=photmod)
    minZ, maxZ = get_limits(propname='Z', photmod=photmod)
    
    limits = np.where((lu[:,0]>minU)&(lu[:,0]<maxU)&(loh12[:,0]>np.log10(minZ))&(loh12[:,0]<np.log10(maxZ))&
                      (lne[:,0]>np.log10(minnH))&(lne[:,0]<np.log10(maxnH))&(lu[:,0]!=const.notnum))[0]
    
    for i in range(lu.shape[1]):        
        lu[:,i][(lu[:,i] > maxU)&(lu[:,i] != const.notnum)] = maxU
        lu[:,i][(lu[:,i] < minU)&(lu[:,i] != const.notnum)] = minU
        
        lne[:,i][(lne[:,i] > np.log10(maxnH))&(lne[:,i] != const.notnum)] = np.log10(maxnH)
        lne[:,i][(lne[:,i] < np.log10(minnH))&(lne[:,i] != const.notnum)] = np.log10(minnH)
        
        loh12[:,i][(loh12[:,i] > np.log10(maxZ))&(loh12[:,i] != const.notnum)] = np.log10(maxZ)
        loh12[:,i][(loh12[:,i] < np.log10(minZ))&(loh12[:,i] != const.notnum)] = np.log10(minZ)
                
    return lms, lssfr, lu, lne, loh12

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
    
def calculate_flux(nebline,redshift,h0=const.h,origin='sfr'):
    
    if nebline.any():
        set_cosmology(omega0=const.omega0, omegab=const.omegab,lambda0=const.lambda0,h0=h0)
        
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
            
def get_lines_Feltre(lu, lne, loh12, verbose=True, 
                     xid_feltre=0.5,alpha_feltre=-1.7):
    '''
    Get the interpolations for the emission lines,
    using the tables
    from Feltre et al. (2016) (https://arxiv.org/pdf/1511.08217utkingalaxies per component.
    lne : floats
     ne of the galaxies per component (cm^-3).
    loh12 : floats
     Metallicity of the galaxies per component (log10(Z))
    verbose : boolean
      If True print out messages
      
    Returns
    -------
    nebline : floats
     Array with the luminosity of the lines per component. (Lsun for L_AGN = 10^45 erg/s)
    '''
    
    minU, maxU = get_limits(propname='U', photmod='feltre16')
    minnH, maxnH = get_limits(propname='nH', photmod='feltre16')
    minZ, maxZ = get_limits(propname='Z', photmod='feltre16')
    minZ, maxZ = np.log10(minZ), np.log10(maxZ)
    
    zmet_str = const.zmet_str['feltre16']
    zmets = np.full(len(zmet_str),const.notnum)
    zmets = np.array([float('0.' + zmet) for zmet in zmet_str])

    logubins = [-5., -4.5, -4., -3.5, -3., -2.5, -2., -1.5, -1.]
    
    nemline = 20
    ndat = lu.shape[0]
    ncomp = lu.shape[1]

    nzmet = 16
    nu = 9

    emline_grid1 = np.zeros((nzmet,nu,nemline))
    emline_grid2 = np.zeros((nzmet,nu,nemline))
    emline_grid3 = np.zeros((nzmet,nu,nemline))

    l = 0
    for k, zname in enumerate(zmets):
        infile = get_zfile(zmet_str[k],photmod='feltre16')
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
                alpha = float(data[3])

                if xid==xid_feltre and alpha==alpha_feltre:
                    if u == -5.:
                        l = 0
                    if u == -4.5:
                        l = 1
                    if u == -4.:
                        l = 2
                    if u == -3.5:
                        l = 3
                    if u == -3.:
                        l = 4
                    if u == -2.5:
                        l = 5
                    if u == -2.:
                        l = 6
                    if u == -1.5:
                        l = 7
                    if u == -1.:
                        l = 8


                    if nH==100 or nH==1000 or nH==10000:
                        for j in range(nemline):
                            if nH == 100:
                                emline_grid1[k,l,j] = float(data[j+4])
                            if nH == 1000:
                                emline_grid2[k,l,j] = float(data[j+4])
                            if nH == 10000:
                                emline_grid3[k,l,j] = float(data[j+4])
        ff.close()

    # log metallicity bins ready for interpolation:

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
    
        # Interpolate over ionisation parameter
        du = []
        j = []
        for logu in lu[:,comp]:
            j1 = locate_interval(logu,logubins)
            if logu<minU:
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

        dz = []
        i = []
    
        for logz in loh12[:,comp]:
            i1 = locate_interval(logz, lzmets)
            if logz<minZ:
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
                emline_int1[k][ii] = (1.-dz[ii])*(1.-du[ii])*emline_grid1[i[ii]][j[ii]][k]+\
                                     dz[ii]*(1-du[ii])*emline_grid1[i[ii]+1][j[ii]][k]+\
                                     (1.-dz[ii])*du[ii]*emline_grid1[i[ii]][j[ii]+1][k]+\
                                     dz[ii]*du[ii]*emline_grid1[i[ii]+1][j[ii]+1][k]
                                     
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
                    nebline[comp][k][n] = (1.-dn)*emline_int1[k][n] + (dn)*emline_int2[k][n]
    
            elif (lne[:,comp][n] > 3. and lne[:,comp][n] <= 4.):
                dn = (lne[:,comp][n] - 3.)/(4. - 3.)
                for k in range(nemline):
                    nebline[comp][k][n] = (1. - dn) * emline_int2[k][n] + (dn) * emline_int3[k][n]
                # print('hay mayor que 3')
    
            elif (lne[:,comp][n] <= 2.):
                for k in range(nemline):
                    nebline[comp][k][n] = emline_int1[k][n]
            elif (lne[:,comp][n] > 4.):
                for k in range(nemline):
                    nebline[comp][k][n] = emline_int3[k][n]
            else:
                print('log(ne)disk out of limits','log(ne)disk = {}'.format(lne[:,comp][n]))
                
    return nebline

def get_lines_Gutkin(lu, lne, loh12, verbose=True,
                     xid_gutkin=0.3,co_gutkin=1,imf_cut_gutkin=100):
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
     Metallicity of the galaxies per component (log10(Z))
    verbose : boolean
      If True print out messages
      
    Returns
    -------
    nebline : floats
     Array with the luminosity of the lines per component. (Lsun per unit SFR(Mo/yr) for 10^8yr)
    '''
    
    minU, maxU = get_limits(propname='U', photmod='feltre16')
    minnH, maxnH = get_limits(propname='nH', photmod='feltre16')
    minZ, maxZ = get_limits(propname='Z', photmod='feltre16')
    minZ, maxZ = np.log10(minZ), np.log10(maxZ)
    
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
        infile = get_zfile(zmet_str[k],photmod='gutkin16')
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
                
                if xid==xid_gutkin and co==co_gutkin and imf_cut==imf_cut_gutkin:
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
            if logu<minU:
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
    
            if logz<minZ:
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
            if logz<minZ:
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
                # print('hay mayor que 3')
    
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


def get_lines(lu, lne, loh12, photmod='gutkin16', verbose=True,
              xid_gutkin=0.3,co_gutkin=1,imf_cut_gutkin=100,
              xid_feltre=0.5,alpha_feltre=-1.7):
    '''
    Get the emission lines

    Parameters
    ----------
    lu : floats
     U of the galaxies per component.
    lne : floats
     ne of the galaxies per component (cm^-3).
    loh12 : floats
     Metallicity of the galaxies per component (log10(Z))
    photomod : string
      Name of the considered photoionisation model.
    verbose : boolean
      If True print out messages

    Returns
    -------
    nebline : floats
     Array with the luminosity of the lines per component. Units depends on photmod.
    '''
    
    if photmod not in const.photmods:
        if verbose:
            print('STOP (eml_photio.get_lines): Unrecognised model to get emission lines.')
            print('                Possible photmod= {}'.format(const.photmods))
        sys.exit()
    elif (photmod == 'gutkin16'):
        nebline = get_lines_Gutkin(lu,lne,loh12,
                verbose=verbose,
                xid_gutkin=xid_gutkin,co_gutkin=co_gutkin,imf_cut_gutkin=imf_cut_gutkin)
    elif (photmod == 'feltre16'):
        nebline = get_lines_Feltre(lu,lne,loh12,
                verbose=verbose,
                xid_feltre=xid_feltre,alpha_feltre=alpha_feltre)

    return nebline

