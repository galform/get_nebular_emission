"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
.. contributions:: Olivia Vidal <ovive.pro@gmail.com>
.. contributions:: Julen Expósito-Márquez <expox7@gmail.com>
"""
import sys
import numpy as np
import h5py
import src.gne_const as const
from src.gne_stats import perc_2arrays
    
def Z_blanc(logM_or):
    logZ = np.zeros(logM_or.shape)
    logM = logM_or - 9.35
    for comp in range(logM_or.shape[1]):
        for i in range(logM_or.shape[0]):
            if logM[i,comp]<(8.7-9.35):
                logZ[i,comp] = 8.37
            elif logM[i,comp]<(9.5-9.35):
                logZ[i,comp] = 8.37 + 0.14*logM[i,comp] - 0.14*(8.7-9.35)
            elif logM[i,comp]<(10.5-9.35):
                logZ[i,comp] = 8.37 + 0.14*(9.5-9.35) - 0.14*(8.7-9.35) + 0.37*logM[i,comp] - 0.37*(9.5-9.35)
            else:
                logZ[i,comp] = 8.37 + 0.14*(9.5-9.35) - 0.14*(8.7-9.35) + 0.37*(10.5-9.35) - 0.37*(9.5-9.35) + 0.03*logM[i,comp] - 0.03*(10.5-9.35)
            
    logZ = logZ - const.ohsun + np.log10(const.zsun) # We leave it in log(Z)
    
    return logZ


def get_Ztremonti(logM,logZ,Lagn_param=[[None],[None]]):
    # Ms and Z scale relation from Tremonti et. al. 2004
    
    try:
        if logZ.shape[1] > 1:
            if Lagn_param[-1][0] != None:
                logMt = np.log10(10**logM[:,0] + Lagn_param[-1])
            else:
                logMt = np.log10(10**logM[:,0] + 10**logM[:,1])
            logZ[:,0] = -1.492 + 1.847*logMt - 0.08026*logMt**2
            # logZ[:,1] = -1.492 + 1.847*logMt - 0.08026*logMt**2
    except:
        logZ = -1.492 + 1.847*logM - 0.08026*logM**2
    
    logZ = logZ - const.ohsun + np.log10(const.zsun) # We leave it in log(Z)
    
    return logMt, logZ

def get_Ztremonti2(logM,logZ,minZ,maxZ,Lagn_param=[[None],[None]]):
    # Correction in bins to Z values using the Ms and Z scale relation from Tremonti et. al. 2004
    
    logZt = np.copy(logZ)    

    logMtot, logZt = get_Ztremonti(logM,logZt,Lagn_param)
    # logZt = Z_blanc(logM)
    
    logZtot = logZ[:,0]
    logZt = logZt[:,0]
    
    # ind_lims = np.where((logZtot > np.log10(minZ)) & (logZtot < np.log10(maxZ)))[0]
    ind_lims = np.where((logZtot > -900) & (logZtot < 900))[0]
    
    smin = 7
    smax = 12
    ds = 0.05
    sbins = np.arange(smin, (smax + ds), ds)
    sbinsH = np.arange(smin, smax, ds)
    shist = sbinsH + ds * 0.5

    median = perc_2arrays(sbins, logMtot[ind_lims], logZtot[ind_lims], 0.5)
    median_t = perc_2arrays(sbins, logMtot[ind_lims], logZt[ind_lims], 0.5)
    ind_med = np.where(median != -999.)[0]

    shist = shist[ind_med]
    median = median[ind_med]
    median_t = median_t[ind_med]
    
    final_bin = sbins[ind_med[-1]+1]
    sbins = sbins[ind_med]
    sbins = np.append(sbins,final_bin)
    
    dif = median_t - median

    for i in range(len(sbins)-1):
        ind = np.where((logMtot>sbins[i])&(logMtot<sbins[i+1]))
        logZ[:,0][ind] = logZ[:,0][ind] + dif[i]
        logZ[:,1][ind] = logZ[:,1][ind] + dif[i]
        
    # smin = 7
    # smax = 12
    # ds = 0.05
    # sbins = np.arange(smin, (smax + ds), ds)
    # sbinsH = np.arange(smin, smax, ds)
    # shist = sbinsH + ds * 0.5

    # median2 = perc_2arrays(sbins, logMtot[ind_lims], logZ[:,0][ind_lims], 0.5)
    # ind_med = np.where(median2 != -999.)[0]
    # median2 = median2[ind_med]
        
    # print(median2)
    # print()
    # print(median_t)
    
    return logZ
    

def alpha_B(T):
    '''
    Given the temperature of the ionizing region, it interpolates from Osterbrock & Ferland 2006
    tables the corresponding value of the hydrogen case B recombination coefficient.

    Parameters
    ----------
    T : floats
     Typical temperature considered for nebular regions.
     
    Returns
    -------
    alphaB : floats
    '''
    
    temps = [5000, 10000, 20000]
    values = [4.54e-13, 2.59e-13, 1.43e-13] # Osterbrock & Ferland 2006, Second Edition, page 22
    
    alphaB = np.interp(T,temps,values)
    
    return alphaB


def surface_density(x,M,reff,profile='exponential',verbose=True):
    '''
    Given the mass of a disk, M, its effective radius, reff,
    the surface density is calculated at a given distance from the center, x.

    Parameters
    ----------
    x : array of floats
     Distance to the center in which surface density is going to be calculated (Mpc).
    M : array of floats
     Mass of the desired component of the galaxy (Msun).
    reff : floats
     Effective radius of the galaxy (Mpc)
    profile : string
     Assumed density profile form for the surface density.
    verbose : boolean
     If True print out messages.
     
    Returns
    -------
    surf_den : floats
    '''
    
    profiles = ['exponential']
    
    if profile not in profiles:
        if verbose:
            print('STOP (gne_une): Unrecognised profile for the surface density.')
            print('                Possible profiles= {}'.format(profiles))
        sys.exit()
    elif profile=='exponential':
        central = M/(2*np.pi*(reff**2))
        
        surf_den = central*np.exp(-x/reff)
        return surf_den

    
def enclosed_mass_disk(x,M,reff,profile='exponential',verbose=True):
    '''
    Given the mass of the desired component of the galaxy, the disk effective radius
    and a distance to the center, it calculates the surface density at that distance.

    Parameters
    ----------
    x : floats
     Distance to the center in which surface density is going to be calculated (Mpc).
    M : floats
     Mass of the desired component of the galaxy (Msun).
    reff : floats
     Effective radius of the galaxy (Mpc)
    profile : string
     Assumed density profile form for the surface density.
    verbose : boolean
     If True print out messages.
     
    Returns
    -------
    surf_den : floats
    '''
    
    profiles = [profile]
    
    if profile not in profiles:
        if verbose:
            print('STOP (gne_une): Unrecognised profile for the surface density.')
            print('                Possible profiles= {}'.format(profiles))
        sys.exit()
    elif profile=='exponential':
        ind = np.where((M>1e-5)&(reff>1e-5))[0]  ###here limitis to 0
        mass_enclosed = np.zeros(M.shape)
        if len(x) > 1:
            mass_enclosed[ind] = (M[ind]/reff[ind])*(reff[ind] - np.exp(-x[ind]/reff[ind])*reff[ind] - x[ind]*np.exp(-x[ind]/reff[ind]))
        else:
            mass_enclosed[ind] = (M[ind]/reff[ind])*(reff[ind] - np.exp(-x[0]/reff[ind])*reff[ind] - x[0]*np.exp(-x[0]/reff[ind]))
            
        return mass_enclosed
    

    
def enclosed_mass_sphere(x,M,reff,profile='exponential',verbose=True):
    '''
    Given the mass of the desired component of the galaxy, the disk effective radius
    and a distance to the center, it calculates the surface density at that distance.

    Parameters
    ----------
    x : floats
     Distance to the center in which surface density is going to be calculated (Mpc).
    M : floats
     Mass of the desired component of the galaxy (Msun).
    reff : floats
     Effective radius of the galaxy (Mpc)
    profile : string
     Assumed density profile form for the surface density.
    verbose : boolean
     If True print out messages.
     
    Returns
    -------
    surf_den : floats
    '''
    
    profiles = ['exponential']
    
    if profile not in profiles:
        if verbose:
            print('STOP (gne_une): Unrecognised profile for the surface density.')
            print('                Possible profiles= {}'.format(profiles))
        sys.exit()
    elif profile=='exponential':
        ind = np.where((M>1e-5)&(reff>1e-5))[0] ###here limits, shouldn't be 0?
        mass_enclosed = np.zeros(M.shape)
        if len(x) > 1:
            mass_enclosed[ind] = (M[ind]/(2*reff[ind]**3))*(2*reff[ind]**3 - np.exp(-x[ind]/reff[ind])*reff[ind]*(2**reff[ind]**2 + 2*reff[ind]*x[ind] + x[ind]**2))
        else:
            mass_enclosed[ind] = (M[ind]/(2*reff[ind]**3))*(2*reff[ind]**3 - np.exp(-x[0]/reff[ind])*reff[ind]*(2**reff[ind]**2 + 2*reff[ind]*x[0] + x[0]**2)) ###here different eqs from mine
        
        return mass_enclosed

    
def vol_sphere(r):
    '''
    Given the radius of a sphere, returns its value.

    Parameters
    ----------
    r : float
     Radius of the sphere.
     
    Returns
    -------
    V : float
    '''
    
    V = (4./3.)*np.pi*r**3
    return V


def mean_density(x,M,r_hm,profile='exponential',bulge=False,verbose=True):
    '''
    Given the mass of the desired component of the galaxy, the disk effective radius
    and a distance to the center, it calculates the particle density at that distance.

    Parameters
    ----------
    x : floats
     Distance to the center in which surface density is going to be calculated (Mpc).
    Ms : floats
     Stellar mass of the galaxy (Msun).
    Mg : floats
     Cold gas mass of the galaxy (Msun).
    r_hm : floats
     Half-mass radius of the galaxy (Mpc)
    profile : string
     Assumed density profile form for the surface density.
    bulge : boolean
     True if the calculation is being applied to a bulge.
     False if the calculation is being applied to a disk.
    verbose : boolean
     If True print out messages.
     
    Returns
    -------
    n : floats
    '''
    
    profiles = ['exponential']
    
    if profile not in profiles:
        if verbose:
            print('STOP (gne_une): Unrecognised profile for the surface density.')
            print('                Possible profiles= {}'.format(profiles))
        sys.exit()
    elif profile=='exponential':
        reff = const.halfmass_to_reff*r_hm # GALFORM ###here not to be hardwired
        
        if bulge:
            M_enclosed = enclosed_mass_sphere(x,M,reff,profile=profile,verbose=verbose)
        else:
            M_enclosed = enclosed_mass_disk(x,M,reff,profile=profile,verbose=verbose)

        n = np.zeros(M_enclosed.shape)
        
        if len(x) > 1:
            ind = np.where((M_enclosed>0)&(x>0))[0]
            n[ind] = M_enclosed[ind]/vol_sphere(x[ind]) / (const.mp*const.kg_to_Msun*const.Mpc_to_cm**3)
        else:
            ind = np.where((M_enclosed>0))[0]
            n[ind] = M_enclosed[ind]/vol_sphere(x[0]) / (const.mp*const.kg_to_Msun*const.Mpc_to_cm**3)
        
        return n

    
def particle_density(x,M,r_hm,T=10000,profile='exponential',verbose=True):
    '''
    Given the mass of the desired component of the galaxy, the disk effective radius
    and a distance to the center, it calculates the particle density at that distance.

    Parameters
    ----------
    x : floats
      Distance to the center in which surface density is going to be calculated (Mpc).
    Ms : floats
      Stellar mass of the galaxy (Msun).
    Mg : floats
      Cold gas mass of the galaxy (Msun).
    r_hm : floats
      Half-mass radius of the galaxy (Mpc)
    T : float
     Typical temperature of ionizing regions.
    profile : string
      Assumed density profile form for the surface density.
    verbose : boolean
      If True print out messages.
     
    Returns
    -------
    n : floats
    '''
    
    reff = const.halfmass_to_reff*r_hm 
    
    den_gas = surface_density(x,M,reff,profile=profile,verbose=verbose)
    
    # h_star = const.reff_to_scale_high*reff
    # den_star = surface_density(x,Ms,reff,profile=profile,verbose=verbose)
    # gamma_gas = gamma_gas_func()
    # gamma_star = gamma_star_func(h_star,den_star)
    # Pext = 0.5*np.pi*const.G*den_gas*(den_gas + (gamma_gas/gamma_star)*den_star) * 1e10/(const.Mpc_to_cm**2)
    
    Pext = 0.5*np.pi*const.G*den_gas**2 * 1e10/(const.Mpc_to_cm**2)
    
    # P = nkT
    n = Pext/(T*const.boltzmann) / const.Mpc_to_cm**3 # cm^-3
    
    return n


def gamma_gas_func():
    '''
    Calculates the velocity dispersion of the gas component (see Lagos et. al. 2011).
     
    Returns
    -------
    gamma_gas : float
    '''
    gamma_gas = 10 #km s^-1, Lagos et al. 2011
    
    return gamma_gas


def gamma_star_func(h_star,den_star):
    '''
    Calculates the velocity disparsion of the star component (see Lagos et. al. 2011).
    
    Parameters
    ----------
    h_star : float
     Stellar scaleheight.
    den_star : float
     Stellar surface density.
     
    Returns
    -------
    gamma_gas : float
    '''
    
    gamma_star = np.sqrt(np.pi*const.G*h_star*den_star) # GALFORM
    
    return gamma_star
    

def mean_density_hydro_eq(max_r,M,r_hm,profile='exponential',verbose=True):
    '''
    Given the mass of the desired component of the galaxy, the disk effective radius
    and a distance to the center, it calculates the mean particle density within that distance.

    Parameters
    ----------
    max_r : floats
      Distance to the center within the surface density is going to be calculated (Mpc).
    Ms : floats
      Stellar mass of the galaxy (Msun).
    Mg : floats
      Cold gas mass of the galaxy (Msun).
    r_hm : floats
      Half-mass radius of the galaxy (Mpc)
    profile : string
      Assumed density profile form for the surface density.
    verbose : boolean
      If True print out messages.
     
    Returns
    -------
    n : floats
    '''
    n = 0
    
    if max_r>0:
        x = np.arange(0,max_r,max_r/100)
        for xbin in x:
            n += particle_density(xbin,M,r_hm,profile=profile,verbose=verbose)
        n = n/len(x)
        return n # cm^-3
    else:
        return n # cm^-3

    
def calculate_ng_hydro_eq(max_r,M,r_hm,profile='exponential',verbose=True):
    '''
    Given the mass of the desired component of the galaxy, the disk effective radius
    and a distance to the center, it calculates the mean particle density within that distance.

    Parameters
    ----------
    max_r : floats
      Distance to the center within the surface density is going to be calculated (Mpc).
    Ms : floats
      Stellar mass of the galaxy (Msun).
    Mg : floats
      Cold gas mass of the galaxy (Msun).
    r_hm : floats
      Half-mass radius of the galaxy (Mpc)
    profile : string
      Assumed density profile form for the surface density.
    verbose : boolean
      If True print out messages.
     
    Returns
    -------
    n : floats
    '''
    
    ng = np.zeros(M.shape)
    if len(max_r) > 1:
        for i in range(len(M)):
            ng[i] = mean_density_hydro_eq(max_r[i],M[i],r_hm[i],profile=profile,verbose=verbose)
    else:
        for i in range(len(M)):
            ng[i] = mean_density_hydro_eq(max_r[0],M[i],r_hm[i],profile=profile,verbose=verbose)
            
    return ng # cm^-3
        

def epsilon_simplemodel(max_r,Mg,r_hm,nH=1000,profile='exponential',bulge=False,verbose=True):
    '''
    Given the mass of the desired component of the galaxy, the disk effective radius
    and a distance to the center, it calculates the volume filling-factor within that distance.

    Parameters
    ----------
    max_r : floats
     Distance to the center within the surface density is going to be calculated (Mpc).
    Ms : floats
     Stellar mass of the galaxy (Msun).
    Mg : floats
     Cold gas mass of the galaxy (Msun).
    r_hm : floats
     Half-mass radius of the galaxy (Mpc).
    nH : float
     Assumed hydrogen density in the ionizing regions.
    profile : string
     Assumed density profile form for the surface density.
    bulge : boolean
     True if the calculation is being applied to a bulge.
     False if the calculation is being applied to a disk.
    verbose : boolean
     If True print out messages.
     
    Returns
    -------
    epsilon : floats
    '''
    
    n = mean_density(max_r,Mg,r_hm,profile=profile,bulge=bulge,verbose=verbose)
    epsilon = n/nH
    
    return n, epsilon


def calculate_epsilon(epsilon_param,max_r,filenom,nH=const.nH_AGN,
                      profile='exponential',verbose=True):
    '''
    It reads the relevant parameters in the input file and calculates 
    the volume filling-factor within that distance.

    Parameters
    ----------
    epsilon_param : array of floats
       Parameters for epsilon calculation.
    max_r : array of floats
       Distance to the center within the surface density is going to be calculated (Mpc).
    filenom : string
       File with output
    nH : float
     Assumed hydrogen density in the ionizing regions.
    profile : string
     Assumed density profile form for the surface density.
    verbose : boolean
     If True print out messages.
     
    Returns
    -------
    epsilon : array of floats
    '''

    if epsilon_param.shape[0] == 2: #2
        Mg, r = epsilon_param
        # Mg = Mg + Mg_bulge
        ind_epsilon = np.where((Mg>5e-5)&(r>5e-5)) ###here why this arbitrary values?
        epsilon = np.zeros(Mg.shape)
        ng = np.zeros(Mg.shape)
        if len(max_r) > 1:
            max_r = max_r[ind_epsilon]
        ng[ind_epsilon], epsilon[ind_epsilon]=epsilon_simplemodel(max_r,
                                                                  Mg[ind_epsilon],r[ind_epsilon],nH=nH,verbose=verbose)
    else:
        Mg, r, Mg_bulge, r_bulge = epsilon_param
        ind_epsilon = np.where((Mg>5e-5)&(r>5e-5))
        epsilon = np.zeros(Mg.shape)
        ng = np.zeros(Mg.shape)
        if len(max_r) > 1:
            max_r = max_r[ind_epsilon]
        ng_disk, ep_disk = epsilon_simplemodel(max_r,
                                               Mg[ind_epsilon],r[ind_epsilon],nH=nH,verbose=verbose)
        ng_bulge, ep_bulge = epsilon_simplemodel(max_r,
                                                 Mg_bulge[ind_epsilon],r_bulge[ind_epsilon],nH=nH,
                                                 bulge=True,verbose=verbose)
        epsilon[ind_epsilon]= ep_disk + ep_bulge
        ng[ind_epsilon]= ng_disk + ng_bulge
    
    epsilon[epsilon>1] = 1
    return epsilon


def n_ratio(n,n_z0):
    '''
    Estimates the metallicity of the AGN from the global metallicity.

    Parameters
    ----------
    n : floats
     Particle number density.
    n_z0 : floats
     Particle number density of the galaxies from the sample at redshift 0.
     
    Returns
    -------
    ratio : floats
    '''
    
    ratio = np.full(n.shape,1.)
    ind = np.where((n>0)&(n_z0>0))[0]
    
    mean = np.mean(n[ind])
    mean_0 = np.mean(n_z0[ind])
    
    ratio = mean/mean_0
    
    return ratio
    

def get_Zagn(logMin,logz):
    '''
    Estimates the metallicity of the AGN from the global metallicity.

    Parameters
    ----------
    logMin : array of floats
       Stellar mass of the galaxy or its components (log10(M*/Msun ???? )).
    logz : array of floats
       Cold gas metallicity (log10(Z)).
     
    Returns
    -------
    lzout : array of floats
        Metallicity of the AGN
    '''

    rows = logMin.shape[0]
    logM = np.zeros(shape=rows)
    lzagn = np.zeros(shape=rows)
    
    if logz.shape[1] >= 2:
        ind = np.where(logz[:,1]>const.notnum)
        lzagn[ind] = np.copy(logz[ind,1])
    
        Ms = 10**logMin
        Ms = np.sum(Ms,axis=1)
        ind = np.where(Ms>0)
        logM[ind] = np.log10(Ms[ind])

    else:
        logM = logMin

    ###here where is this justified?
    lzagn[(logM>9.5) & (logM<=10)] = lzagn[(logM>9.5) & (logM<=10)] + 0.1
    lzagn[(logM>10)  & (logM<=10.5)] = lzagn[(logM>10)  & (logM<=10.5)] + 0.3
    lzagn[(logM>10.5)& (logM<=11)] = lzagn[(logM>10.5)& (logM<=11)] + 0.1

    ###here input should directly be zbuldge and M* total
    lzout = np.zeros(shape=(np.shape(logz))); lzout.fill(const.notnum)
    lzout[:,0] = lzagn

    return lzout

######################

def phot_rate(lssfr=None, lms=None, IMF=None, Lagn=None, origin='sfr'):
    '''
    Given log10(Mstar), log10(sSFR), log10(Z), Lagn and the assumed IMF,
    get the rate of ionizing photons in photon per second.

    Parameters
    ----------
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lzgas : floats
     Metallicity of the galaxies per component (log10(Z)).
    IMF : array of strings
     Assumed IMF for the input data of each component.
    Lagn : floats
     Bolometric luminosity of the AGN (Lsun).
    origin : string
     Source of the emission (star-forming region or AGN)
     
    Returns
    -------
    Q : floats
    '''
    
    if origin=='sfr':
        Q = np.zeros(np.shape(lssfr))
        for comp in range(Q.shape[1]):
            ###here ref. missing
            Q[:,comp] = 10**(lssfr[:,comp] + lms[:,comp]) * const.IMF_SFR[IMF[comp]] * const.phot_to_sfr_kenn
            # lssfr[:,comp] = np.log10(Q[:,comp]/(const.IMF_SFR[IMF[comp]] * const.phot_to_sfr_kenn)) - lms[:,comp]
            
    if origin=='agn':
        Q = np.zeros(np.shape(lssfr))
        ind = np.where(Lagn>0)[0]
        # Q[ind,0] = Lagn[ind]*2.3e10 # Panda 2022
        Q[ind,0] = Lagn[ind]*((3.28e15)**-1.7)/(1.7*8.66e-11*const.planck) # Feltre 2016
        # This implies that Lion = Lbol/5 aprox.
            
    return Q


def get_une_kashino20(Q, lms, lssfr, lzgas, T, ng_ratio, IMF):
    '''
    Characterise the SF ionising region from global galactic properties,
    using the model from
    Kashino & Inoue 2019 (https://arxiv.org/abs/1812.06939).

    Parameters
    ----------
    Q : floats
     Rate of ionizing photons (phot/s).
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lzgas : floats
     Metallicity of the galaxies per component (log10(Z)).
    T : float
     Typical temperature of ionizing regions.
    ng_ratio : floats
     Ratio between the mean particle number density of the cold gas of the 
     input sample and the sample at redshift 0.
    IMF : array of strings
     Assumed IMF for the input data of each component.

    Returns
    -------
    lu, lne, lzgas : floats
    '''

    ###here missing transformation to the IMF assumed by Kashino
    lu, lne = [np.full(np.shape(lms), const.notnum) for i in range(2)]
    
    lssfr_new = np.full(np.shape(lssfr),const.notnum)
    for comp in range(lssfr.shape[1]):
        for i in range(lssfr.shape[0]):
            if Q[i,comp] == 0:
                continue
            else:
                ###here why do we need this?
                lssfr_new[i,comp] = np.log10(Q[i,comp]/(const.IMF_SFR[IMF[comp]] * const.phot_to_sfr_kenn)) - lms[i,comp]

    ind = np.where((lssfr_new > const.notnum) &
                   (lms > 0) &
                   (lzgas > const.notnum))
    
    
    if (np.shape(ind)[1]>1):
        # Transform the metallicity into log10(O/H)+12
        lzgas[ind] = lzgas[ind] + const.ohsun - np.log10(const.zsun)

        # Apply equation
        lne[ind] = 2.066 + 0.310*(lms[ind]-10) + 0.492*(lssfr_new[ind] + 9.)

        # Eq. 12 from Kashino & Inoue 2019 ####here esta eq está mal
        lu[ind] =  -2.316 - 0.360*(lzgas[ind] -8.) -0.292*lne[ind] + 0.428*(lssfr_new[ind] + 9.)
        # lu[ind] =  -3.073 - 0.137*(lms[ind]-10) + 0.372*(lssfr[ind] + 9.)
        lzgas[ind] = lzgas[ind] - const.ohsun + np.log10(const.zsun) # We leave it in log(Z)
        
    ind = np.where((lssfr > const.notnum) &
                   (lms > 0) &
                   (lzgas > const.notnum))
    
    ind_comp = []   
    for comp in range(len(Q[0])):
        ind_comp.append(np.where((lssfr[:,comp] > const.notnum) &
                       (lms[:,comp] > 0) &
                       (lzgas[:,comp] > const.notnum) &
                       (Q[:,comp] > 0))[0])
        
    epsilon = np.full(np.shape(lssfr),const.notnum)
    cte = np.zeros(np.shape(lssfr))
    
    for comp in range(len(Q[0])):
        epsilon[:,comp][ind_comp[comp]] = ((1/alpha_B(T)) * ((4*const.c_cm*(10**lu[:,comp][ind_comp[comp]]))/3)**(3/2) * 
                              ((4*np.pi)/(3*Q[:,comp][ind_comp[comp]]*(10**lne[:,comp][ind_comp[comp]])))**(1/2))
        
        if ng_ratio != None:
            epsilon[:,comp][ind_comp[comp]] = epsilon[:,comp][ind_comp[comp]] * ng_ratio
        
        cte[:,comp][ind_comp[comp]] = 3*(alpha_B(T)**(2/3)) * (3*epsilon[:,comp][ind_comp[comp]]**2*(10**lne[:,comp][ind_comp[comp]])/(4*np.pi))**(1/3) / (4*const.c_cm)    
    
    lu[ind] = np.log10(cte[ind] * Q[ind]**(1/3))

    return lu, lne, lzgas


def get_une_orsi14(Q, lms, lssfr, lzgas, T, q0, z0, gamma, ng_ratio):
    '''
    Given log10(Mstar), log10(sSFR), log10(Z) and the values for the free parameters,
    get the ionizing parameter, logU, and the electron density, logne,
    using the model from Orsi 2014.

    Parameters
    ----------
    Q : floats
     Rate of ionizing photons (phot/s).
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lzgas : floats
     Metallicity of the galaxies per component (log10(Z)).
    T : float
     Typical temperature of ionizing regions.
    q0 : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    z0 : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    gamma : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    ng_ratio : floats
     Ratio between the mean particle number density of the cold gas of the 
     input sample and the sample at redshift 0.

    Returns
    -------
    lu, lne, lzgas : floats
    '''
    
    lu, lne = [np.full(np.shape(lms), const.notnum) for i in range(2)]

    ind = np.where((lssfr > const.notnum) &
                   (lms > 0) &
                   (lzgas > const.notnum))
    
    if (np.shape(ind)[1]>1):
        lne[ind] = 2.066 + 0.310*(lms[ind]-10) + 0.492*(lssfr[ind] + 9.)
        lu[ind] = np.log10(q0*((10**lzgas[ind])/z0)**-gamma / const.c_cm)
        
    ind = np.where((lssfr > const.notnum) &
                   (lms > 0) &
                   (lzgas > const.notnum))
    
    ind_comp = []   
    for comp in range(len(Q[0])):
        ind_comp.append(np.where((lssfr[:,comp] > const.notnum) &
                       (lms[:,comp] > 0) &
                       (lzgas[:,comp] > const.notnum) &
                       (Q[:,comp] > 0))[0])
        
    epsilon = np.full(np.shape(lssfr),const.notnum)
    cte = np.zeros(np.shape(lssfr))
    
    for comp in range(len(Q[0])):
        epsilon[:,comp][ind_comp[comp]] = ((1/alpha_B(T)) * ((4*const.c_cm*(10**lu[:,comp][ind_comp[comp]]))/3)**(3/2) * 
                              ((4*np.pi)/(3*Q[:,comp][ind_comp[comp]]*(10**lne[:,comp][ind_comp[comp]])))**(1/2))
        
        if ng_ratio != None:
            epsilon[:,comp][ind_comp[comp]] = epsilon[:,comp][ind_comp[comp]] * ng_ratio
        
        cte[:,comp][ind_comp[comp]] = 3*(alpha_B(T)**(2/3)) * (3*epsilon[:,comp][ind_comp[comp]]**2*(10**lne[:,comp][ind_comp[comp]])/(4*np.pi))**(1/3) / (4*const.c_cm)    
    
    lu[ind] = np.log10(cte[ind] * Q[ind]**(1/3))
    

    return lu, lne, lzgas

# def get_une_carton17(lms, lssfr, lzgas):
#     '''
#     Given log10(Mstar), log10(sSFR), log10(Z),
#     get the ionizing parameter, logU, and the electron density, logne,
#     using the model from Carton 2017.

#     Parameters
#     ----------
#     lms : floats
#      Masses of the galaxies per component (log10(M*) (Msun)).
#     lssfr : floats
#      sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
#     lzgas : floats
#      Metallicity of the galaxies per component (log10(Z)).

#     Returns
#     -------
#     lu, lne, lzgas : floats
#     '''
    
#     lu, lne = [np.full(np.shape(lms), const.notnum) for i in range(2)]

#     ind = np.where((lssfr > const.notnum) &
#                    (lms > 0) &
#                    (lzgas > const.notnum))
    
#     if (np.shape(ind)[1]>1):
#         lne[ind] = 2.066 + 0.310*(lms[ind]-10) + 0.492*(lssfr[ind] + 9.)
#         lu[ind] = -0.8*np.log10(10**lzgas[ind]/const.zsun) - 3.58   

#     return lu, lne, lzgas


def get_une_panuzzo03(Q, lms, lssfr, lzgas, T, epsilon0, ng_ratio, origin, IMF):
    '''
    Given the rate of ionizing photons, log10(Mstar), log10(sSFR), log10(Z),
    the assumed temperature for the ionizing regions, the volume filling-factor
    and the assumed IMF,
    get the ionizing parameter, logU, and the electron density, logne,
    using the model from Panuzzo 2003.

    Parameters
    ----------
    Q : floats
     Rate of ionizing photons (phot/s)
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lzgas : floats
     Metallicity of the galaxies per component (log10(Z)).
    T : float
     Typical temperature of ionizing regions.
    epsilon0 : floats
     Volume filling-factor of the galaxies.
    ng_ratio : floats
     Ratio between the mean particle number density of the cold gas of the 
     input sample and the sample at redshift 0.
    origin : string
     Source of the ionizing photons.
    IMF : array of strings
     Assumed IMF for the input data of each component.

    Returns
    -------
    lu, lne, lzgas : floats
    '''
    
    lzgas_all = np.copy(lzgas)
    
    lu, lne, lzgas = [np.full(np.shape(lms), const.notnum) for i in range(3)]

    ind = np.where((lssfr > const.notnum) &
                   (lms > 0) &
                   (lzgas_all > const.notnum) &
                   (Q > 0))
    
    # ind1 = np.where((lssfr[:,0] > const.notnum) &
    #                (lms[:,0] > 0) &
    #                (lzgas[:,0] > const.notnum) &
    #                (Q[:,0] > 0))[0]
    # ind2 = np.where((lssfr[:,1] > const.notnum) &
    #                (lms[:,1] > 0) &
    #                (lzgas[:,1] > const.notnum) &
    #                (Q[:,1] > 0))[0]
    
    # ind_comp = [ind1,ind2]
    
    ind_comp = []   
    for comp in range(len(Q[0])):
        ind_comp.append(np.where((lssfr[:,comp] > const.notnum) &
                       (lms[:,comp] > 0) &
                       (lzgas_all[:,comp] > const.notnum) &
                       (Q[:,comp] > 0))[0])
    
    if (np.shape(ind)[1]>1):
        
        epsilon = np.full(np.shape(lssfr),const.notnum)
        cte = np.zeros(np.shape(lssfr))
        
        if origin=='sfr':
            # lu, lne, lzgas = get_une_orsi14(Q, lms, lssfr, lzgas, T, q0=const.q0_orsi, z0=const.Z0_orsi, gamma=1.3)
            lu, lne, lzgas = get_une_kashino20(Q,lms,lssfr,lzgas_all,T,ng_ratio,IMF)
            
            for comp in range(len(Q[0])):
                epsilon[:,comp][ind_comp[comp]] = ((1/alpha_B(T)) * ((4*const.c_cm*(10**lu[:,comp][ind_comp[comp]]))/3)**(3/2) * 
                                      ((4*np.pi)/(3*Q[:,comp][ind_comp[comp]]*(10**lne[:,comp][ind_comp[comp]])))**(1/2))
                
                if ng_ratio != None:
                    epsilon[:,comp][ind_comp[comp]] = epsilon[:,comp][ind_comp[comp]] * ng_ratio
                
                cte[:,comp][ind_comp[comp]] = 3*(alpha_B(T)**(2/3)) * (3*epsilon[:,comp][ind_comp[comp]]**2*(10**lne[:,comp][ind_comp[comp]])/(4*np.pi))**(1/3) / (4*const.c_cm)    
            
            lu[ind] = np.log10(cte[ind] * Q[ind]**(1/3))
        
        if origin=='agn':
            lne[ind] = 3
            lzgas[ind] = lzgas_all[ind]
            
            for comp in range(len(Q[0])):
                
                epsilon[:,comp][ind_comp[comp]] = epsilon0[ind_comp[comp]]
                
                cte[:,comp][ind_comp[comp]] = ( (3*(alpha_B(T)**(2/3)) / (4*const.c_cm)) 
                 * (3*epsilon[:,comp][ind_comp[comp]]**2*(10**lne[:,comp][ind_comp[comp]])/(4*np.pi))**(1/3) )
                
            cte[cte==0] = 1e-50
            lu[ind] = np.log10(cte[ind] * Q[ind]**(1/3) / 3)
            lu[cte==1e-50] = const.notnum
    

    return lu, lne, lzgas


def get_une(lms_o, lssfr_o, lzgas_o,filenom,
            q0=const.q0_orsi, z0=const.Z0_orsi, Lagn=None, ng_ratio=None,
            Z_central_cor=False,
            gamma=1.3, T=10000, epsilon_param=[None], epsilon_param_z0=[None],
            epsilon=0.01,IMF=['Kennicut','Kennicut'],
            unemod='kashino20', origin='sfr', verbose=True):
    '''
    Given the global properties of a galaxy or a region
    (log10(Mstar), log10(sSFR) and 12+log(O/H)),
    get the characteristics of the ionising region
    (ionizing parameter, U, and the electron density, ne).

    Parameters
    ----------
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lzgas : floats
     Metallicity of the galaxies per component (log10(Z)).
    q0 : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    z0 : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    gamma : float
     Ionization parameter constant to calibrate Orsi 2014 model for nebular regions. q0(z/z0)^-gamma
    Lagn : floats
     Bolometric luminosity of the AGNs (erg/s).
    T : float
     Typical temperature of ionizing regions.
    epsilon_param : floats
     Parameters for epsilon calculation.
    epsilon_param_z0 : floats
     Parameters for epsilon calculation in the sample of galaxies at redshift 0.
    epsilon : floats
     Volume filling-factor of the galaxy.
    IMF : array of strings
     Assumed IMF for the input data of each component.
    unemod : string
     Model to go from galaxy properties to U and ne.
    origin : string
     Source of the ionizing photons.
    verbose : boolean
     Yes = print out messages.

    Returns
    -------
    Q, lu, lne, lzgas : floats
    '''

    # Read redshift
    f = h5py.File(filenom, 'r')
    header = f['header']
    redshift = header.attrs['redshift']
    f.close()
    
    # ncomp = len(lms[0])
    Q = phot_rate(lssfr=lssfr_o,lms=lms_o,IMF=IMF,Lagn=Lagn,origin=origin)
    
    epsilon = None
    if origin=='agn' and epsilon_param is not None:
        epsilon = calculate_epsilon(epsilon_param,[const.radius_NLR],
                                    filenom,nH=const.nH_AGN,
                                    profile='exponential',verbose=verbose)
    if origin=='sfr' and epsilon_param_z0 is not None:
        # ng = calculate_ng_hydro_eq(2*epsilon_param[1],epsilon_param[0],epsilon_param[1],profile='exponential',verbose=True)
        # epsilon = ng/const.nH_gal
        # epsilon[epsilon>1] = 1
            
        
        # ng_z0 = calculate_ng_hydro_eq(2*epsilon_param_z0[1],epsilon_param_z0[0],epsilon_param_z0[1],profile='exponential',verbose=True)
        # ng_ratio = n_ratio(ng,ng_z0)
        if redshift==0.8:
            ng_ratio = const.med_to_low
        elif redshift==1.5:
            ng_ratio = const.high_to_low
        else:
            ng_ratio = 1.
                    
    if Z_central_cor and origin=='agn':
        lzgas = get_Zagn(lms_o,lzgas_o)
    else:
        lzgas = np.copy(lzgas_o)
    
    if unemod not in const.unemods:
        if verbose:
            print('STOP (gne_une): Unrecognised model to get U and ne.')
            print('                Possible unemod= {}'.format(const.unemods))
        sys.exit()
    elif (unemod == 'kashino20'):
        lu, lne, lzgas = get_une_kashino20(Q,lms_o,lssfr_o,lzgas,T,ng_ratio,IMF)
    elif (unemod == 'orsi14'):
        lu, lne, lzgas = get_une_orsi14(Q,lms_o,lssfr_o,lzgas,T,q0,z0,gamma,ng_ratio)
    elif (unemod == 'panuzzo03'):
        lu, lne, lzgas = get_une_panuzzo03(Q,lms_o,lssfr_o,lzgas,T,epsilon,ng_ratio,origin,IMF)
        
    return Q, lu, lne, lzgas, epsilon, ng_ratio
