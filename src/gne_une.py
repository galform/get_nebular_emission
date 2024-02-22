import sys

import numpy as np
import h5py
import src.gne_const as const
from scipy.stats import gaussian_kde
from src.gne_stats import perc_2arrays
from src.gne_io import get_nheader


def bursttobulge(lms,Lagn_param):
    '''
    Changes the bulge component of the stellar mass from the mass of the starburst
    to the total mass of the bulge.
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    Lagn_params : floats
     Parameters to calculate the AGN emission. 
     The last one is always the stellar mass of the bulge.
    '''
    ind = np.where(Lagn_param[-1]>0)
    lms[:,1] = const.notnum
    lms[:,1][ind] = np.log10(Lagn_param[-1][ind])

    
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


def Z_tremonti(logM,logZ,Lagn_param=[[None],[None]]): # Ms and Z scale relation from Tremonti et. al. 2004
    
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

def Z_tremonti2(logM,logZ,minZ,maxZ,Lagn_param=[[None],[None]]): # Correction in bins to Z values using the Ms and Z scale relation from Tremonti et. al. 2004
    
    logZt = np.copy(logZ)    

    logMtot, logZt = Z_tremonti(logM,logZt,Lagn_param)
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
    

def Ledd(Mbh): # Eddington luminosity
    '''
    Given the mass of the black hole, calculates the Eddington luminosity.

    Parameters
    ----------
    Mbh : floats
     Mass of the black hole (Msun).
     
    Returns
    -------
    Ledd : floats
    '''
    
    Ledd = 1.26e38*Mbh # erg s^-1
    return Ledd
    

def acc_rate_edd(Mbh): # Eddington mass accretion rate
    '''
    Given the mass of the black hole, calculates the eddington accretion rate.

    Parameters
    ----------
    Mbh : floats
     Mass of the black hole (Msun).
     
    Returns
    -------
    acc_rate : floats
    '''
    
    acc_rate = Ledd(Mbh)/(0.1*const.c**2) * (const.kg_to_Msun/1000) # Msun/s
    return acc_rate


def t_bulge(r_bulge, v_bulge): # Dynamical timescale of the bulge
    '''
    Given the bulge's half-mass radius and circular velocity at the half-mass radius,
    calculates the dynamical timescale of the bulge.

    Parameters
    ----------
    r_bulge : floats
     Half-mass radius of the bulge (Mpc).
    v_bulge : floats
     Circular velocity at the half-mass radius of the bulge (km/s).
     
    Returns
    -------
    dyn_time : floats
    '''
    
    dyn_time = r_bulge/v_bulge * 1e-5*const.Mpc_to_cm # s
    return dyn_time


def acc_rate_quasar(M_bulge, r_bulge, v_bulge):
    '''
    Given the mass of the bulge, the half-mass radius and circular velocity 
    at the half-mass radius, calculates the accretion rate of the quasar (starburst) mode.

    Parameters
    ----------
    M_bulge : floats
     Mass of the bulge (Msun)
    r_bulge : floats
     Half-mass radius of the bulge (Mpc).
    v_bulge : floats
     Circular velocity at the half-mass radius of the bulge (km/s).
     
    Returns
    -------
    acc_rate : floats
    '''
    
    acc_rate = M_bulge*const.fbh/(t_bulge(r_bulge, v_bulge)*const.fq) # Msun/s
    return acc_rate

def acc_rate_radio(Mhot, Mbh):
    '''
    Given the mass of the hot gas and the mass of the black hole, calculates
    the accretion rate of the radio (hot gas) mode.

    Parameters
    ----------
    Mhot : floats
     Mass of the hot gas (Msun).
    Mbh : floats
     Mass of the black hole (Msun).
     
    Returns
    -------
    acc_rate : floats
    '''
    
    acc_rate = const.kappa_agn*(Mhot*Mbh*1e-19)**const.kappa_agn_exp * (1/3.154e7) # Msun/s
    return acc_rate


def r_iso(spin):
    Z1 = 1 + (1 - np.abs(spin)**2)**(1/3) * ((1+np.abs(spin))**(1/3) + (1-np.abs(spin))**(1/3))
    Z2 = np.sqrt(3*np.abs(spin)**2 + Z1**2)
    r_iso = 3 + Z2 - np.sign(spin)*np.sqrt((3-Z1)*(3+Z1+2*Z2))
    
    return r_iso
    
def epsilon_td(spin):
    r = r_iso(spin)
    epsilon_td = 1 - (1 - (2/(3*r)))**(1/2) # Radiative accretion efficiency for a thin accretion disc (general approximation, page 5)
    
    return epsilon_td

# AGNinputs = ['Lagn', 'acc_rate', 'acc_rates', 'radio_mode', 'quasar_mode', 'complete']


########################

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
    
    profiles = ['exponential']
    
    if profile not in profiles:
        if verbose:
            print('STOP (gne_une): Unrecognised profile for the surface density.')
            print('                Possible profiles= {}'.format(profiles))
        sys.exit()
    elif profile=='exponential':
        ind = np.where((M>1e-5)&(reff>1e-5))[0]
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
        ind = np.where((M>1e-5)&(reff>1e-5))[0]
        mass_enclosed = np.zeros(M.shape)
        if len(x) > 1:
            mass_enclosed[ind] = (M[ind]/(2*reff[ind]**3))*(2*reff[ind]**3 - np.exp(-x[ind]/reff[ind])*reff[ind]*(2**reff[ind]**2 + 2*reff[ind]*x[ind] + x[ind]**2))
        else:
            mass_enclosed[ind] = (M[ind]/(2*reff[ind]**3))*(2*reff[ind]**3 - np.exp(-x[0]/reff[ind])*reff[ind]*(2**reff[ind]**2 + 2*reff[ind]*x[0] + x[0]**2))
        
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
        reff = const.halfmass_to_reff*r_hm # GALFORM
        
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


def calculate_epsilon(epsilon_param,max_r,h0=const.h,nH=1000,profile='exponential',verbose=True):
    '''
    It reads the relevant parameters in the input file and calculates 
    the volume filling-factor within that distance.

    Parameters
    ----------
    infile : string
     - Name of the input file. 
     - In text files (*.dat, *txt, *.cat), columns separated by ' '.
     - In csv files (*.csv), columns separated by ','.
    max_r : floats
     Distance to the center within the surface density is going to be calculated (Mpc).
    epsilon_param : floats
     Parameters for epsilon calculation.
    cut : integers
     Indeces of the galaxies which have survived the specified cut in the main program. 
    h0 : float
     If not None: value of h, H0=100h km/s/Mpc.
    nH : float
     Assumed hydrogen density in the ionizing regions.
    profile : string
     Assumed density profile form for the surface density.
    verbose : boolean
     If True print out messages.
     
    Returns
    -------
    epsilon : floats
    '''
    
    if epsilon_param.shape[0] == 2: #2
        Mg, r = epsilon_param
        # Mg = Mg + Mg_bulge
        ind_epsilon = np.where((Mg>5e-5)&(r>5e-5))
        epsilon = np.zeros(Mg.shape)
        ng = np.zeros(Mg.shape)
        if len(max_r) > 1:
            max_r = max_r[ind_epsilon]
        if h0:
            ng[ind_epsilon], epsilon[ind_epsilon]=epsilon_simplemodel(max_r,
                    Mg[ind_epsilon]/const.h,r[ind_epsilon]/const.h,nH=nH,verbose=verbose)
        else:
            ng[ind_epsilon], epsilon[ind_epsilon]=epsilon_simplemodel(max_r,
                    Mg[ind_epsilon],r[ind_epsilon],nH=nH,verbose=verbose)
    else:
        Mg, r, Mg_bulge, r_bulge = epsilon_param
        ind_epsilon = np.where((Mg>5e-5)&(r>5e-5))
        epsilon = np.zeros(Mg.shape)
        ng = np.zeros(Mg.shape)
        if len(max_r) > 1:
            max_r = max_r[ind_epsilon]
        if h0:
            ng_disk, ep_disk = epsilon_simplemodel(max_r,
                    Mg[ind_epsilon]/const.h,r[ind_epsilon]/const.h,nH=nH,verbose=verbose)
            ng_bulge, ep_bulge = epsilon_simplemodel(max_r,
                        Mg_bulge[ind_epsilon]/const.h,r_bulge[ind_epsilon]/const.h,nH=nH,
                        bulge=True,verbose=verbose)
            epsilon[ind_epsilon]= ep_disk + ep_bulge
            ng[ind_epsilon]= ng_disk + ng_bulge
        else:
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
    

def Zagn(logM,logz):
    '''
    Estimates the metallicity of the AGN from the global metallicity.

    Parameters
    ----------
    logM : string
     Stellar mass of the galaxy (log10(Msun)).
    logz : floats
     Cold gas metallicity (log10(Z)).
     
    Returns
    -------
    logz : floats
    '''
    
    if logz.shape[1] >= 2:
        ind = np.where(logz[:,1]>const.notnum)
        logz[ind,0] = np.copy(logz[ind,1])
        logz[:,1] = const.notnum
    
    Ms = 10**logM
    Ms = np.sum(Ms,axis=1)
    
    ind = np.where(Ms>0)
    Ms[ind] = np.log10(Ms[ind])
    
    # logz = logz + 0.1
    
    for i in range(len(Ms)):
        if Ms[i]<9.5:
            continue
        elif Ms[i]<10:
            logz[i] = logz[i] + 0.1
        elif Ms[i]<10.5:
            logz[i] = logz[i] + 0.3
        elif Ms[i]<11:
            logz[i] = logz[i] + 0.1
        else:
            continue
    
    return logz

######################

def phot_rate(lssfr=None, lms=None, IMF_f=None, Lagn=None, origin='sfr'):
    '''
    Given log10(Mstar), log10(sSFR), log10(Z), Lagn and the assumed IMF,
    get the rate of ionizing photons in photon per second.

    Parameters
    ----------
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    loh12 : floats
     Metallicity of the galaxies per component (log10(Z)).
    IMF_f : string
     Assumed IMF.
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
            Q[:,comp] = 10**(lssfr[:,comp] + lms[:,comp]) * const.IMF_SFR[IMF_f[comp]] * const.phot_to_sfr_kenn
            # lssfr[:,comp] = np.log10(Q[:,comp]/(const.IMF_SFR[IMF_f[comp]] * const.phot_to_sfr_kenn)) - lms[:,comp]
    if origin=='agn':
        Q = np.zeros(np.shape(lssfr))
        ind = np.where(Lagn>0)[0]
        # Q[ind,0] = Lagn[ind]*2.3e10 # Panda 2022
        Q[ind,0] = Lagn[ind]*((3.28e15)**-1.7)/(1.7*8.66e-11*const.planck) # Feltre 2016
        # This implies that Lion = Lbol/5 aprox.
            
    return Q


def get_une_kashino20(Q, lms, lssfr, loh12, T, ng_ratio, IMF_f):
    '''
    Given log10(Mstar), log10(sSFR), log10(Z) and the assumed IMF,
    get the ionizing parameter, logU, and the electron density, logne,
    using the model from Kashino 2020.

    Parameters
    ----------
    Q : floats
     Rate of ionizing photons (phot/s).
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    loh12 : floats
     Metallicity of the galaxies per component (log10(Z)).
    T : float
     Typical temperature of ionizing regions.
    ng_ratio : floats
     Ratio between the mean particle number density of the cold gas of the 
     input sample and the sample at redshift 0.
    IMF_f : string
     Assumed IMF.

    Returns
    -------
    lu, lne, loh12 : floats
    '''
    
    lu, lne = [np.full(np.shape(lms), const.notnum) for i in range(2)]
    
    lssfr_new = np.full(np.shape(lssfr),const.notnum)
    for comp in range(lssfr.shape[1]):
        for i in range(lssfr.shape[0]):
            if Q[i,comp] == 0:
                continue
            else:
                lssfr_new[i,comp] = np.log10(Q[i,comp]/(const.IMF_SFR[IMF_f[comp]] * const.phot_to_sfr_kenn)) - lms[i,comp]

    ind = np.where((lssfr_new > const.notnum) &
                   (lms > 0) &
                   (loh12 > const.notnum))
    
    
    if (np.shape(ind)[1]>1):
        loh12[ind] = loh12[ind] + const.ohsun - np.log10(const.zsun)
        
        lne[ind] = 2.066 + 0.310*(lms[ind]-10) + 0.492*(lssfr_new[ind] + 9.)
        # lne[ind] = 2
        lu[ind] =  -2.316 - 0.360*(loh12[ind] -8.) -0.292*lne[ind] + 0.428*(lssfr_new[ind] + 9.)
        # lu[ind] =  -3.073 - 0.137*(lms[ind]-10) + 0.372*(lssfr[ind] + 9.)
        loh12[ind] = loh12[ind] - const.ohsun + np.log10(const.zsun) # We leave it in log(Z)
        
    ind = np.where((lssfr > const.notnum) &
                   (lms > 0) &
                   (loh12 > const.notnum))
    
    ind_comp = []   
    for comp in range(len(Q[0])):
        ind_comp.append(np.where((lssfr[:,comp] > const.notnum) &
                       (lms[:,comp] > 0) &
                       (loh12[:,comp] > const.notnum) &
                       (Q[:,comp] > 0))[0])
        
    epsilon = np.full(np.shape(lssfr),const.notnum)
    cte = np.zeros(np.shape(lssfr))
    
    for comp in range(len(Q[0])):
        epsilon[:,comp][ind_comp[comp]] = ((1/alpha_B(T)) * ((4*const.c*(10**lu[:,comp][ind_comp[comp]]))/3)**(3/2) * 
                              ((4*np.pi)/(3*Q[:,comp][ind_comp[comp]]*(10**lne[:,comp][ind_comp[comp]])))**(1/2))
        
        if ng_ratio != None:
            epsilon[:,comp][ind_comp[comp]] = epsilon[:,comp][ind_comp[comp]] * ng_ratio
        
        cte[:,comp][ind_comp[comp]] = 3*(alpha_B(T)**(2/3)) * (3*epsilon[:,comp][ind_comp[comp]]**2*(10**lne[:,comp][ind_comp[comp]])/(4*np.pi))**(1/3) / (4*const.c)    
    
    lu[ind] = np.log10(cte[ind] * Q[ind]**(1/3))

    return lu, lne, loh12


def get_une_orsi14(Q, lms, lssfr, loh12, T, q0, z0, gamma, ng_ratio):
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
    loh12 : floats
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
    lu, lne, loh12 : floats
    '''
    
    lu, lne = [np.full(np.shape(lms), const.notnum) for i in range(2)]

    ind = np.where((lssfr > const.notnum) &
                   (lms > 0) &
                   (loh12 > const.notnum))
    
    if (np.shape(ind)[1]>1):
        lne[ind] = 2.066 + 0.310*(lms[ind]-10) + 0.492*(lssfr[ind] + 9.)
        lu[ind] = np.log10(q0*((10**loh12[ind])/z0)**-gamma / const.c)
        
    ind = np.where((lssfr > const.notnum) &
                   (lms > 0) &
                   (loh12 > const.notnum))
    
    ind_comp = []   
    for comp in range(len(Q[0])):
        ind_comp.append(np.where((lssfr[:,comp] > const.notnum) &
                       (lms[:,comp] > 0) &
                       (loh12[:,comp] > const.notnum) &
                       (Q[:,comp] > 0))[0])
        
    epsilon = np.full(np.shape(lssfr),const.notnum)
    cte = np.zeros(np.shape(lssfr))
    
    for comp in range(len(Q[0])):
        epsilon[:,comp][ind_comp[comp]] = ((1/alpha_B(T)) * ((4*const.c*(10**lu[:,comp][ind_comp[comp]]))/3)**(3/2) * 
                              ((4*np.pi)/(3*Q[:,comp][ind_comp[comp]]*(10**lne[:,comp][ind_comp[comp]])))**(1/2))
        
        if ng_ratio != None:
            epsilon[:,comp][ind_comp[comp]] = epsilon[:,comp][ind_comp[comp]] * ng_ratio
        
        cte[:,comp][ind_comp[comp]] = 3*(alpha_B(T)**(2/3)) * (3*epsilon[:,comp][ind_comp[comp]]**2*(10**lne[:,comp][ind_comp[comp]])/(4*np.pi))**(1/3) / (4*const.c)    
    
    lu[ind] = np.log10(cte[ind] * Q[ind]**(1/3))
    

    return lu, lne, loh12

# def get_une_carton17(lms, lssfr, loh12):
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
#     loh12 : floats
#      Metallicity of the galaxies per component (log10(Z)).

#     Returns
#     -------
#     lu, lne, loh12 : floats
#     '''
    
#     lu, lne = [np.full(np.shape(lms), const.notnum) for i in range(2)]

#     ind = np.where((lssfr > const.notnum) &
#                    (lms > 0) &
#                    (loh12 > const.notnum))
    
#     if (np.shape(ind)[1]>1):
#         lne[ind] = 2.066 + 0.310*(lms[ind]-10) + 0.492*(lssfr[ind] + 9.)
#         lu[ind] = -0.8*np.log10(10**loh12[ind]/const.zsun) - 3.58   

#     return lu, lne, loh12


def get_une_panuzzo03(Q, lms, lssfr, loh12, T, epsilon0, ng_ratio, origin, IMF_f):
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
    loh12 : floats
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
    IMF_f : string
     Assumed IMF.

    Returns
    -------
    lu, lne, loh12 : floats
    '''
    
    loh12_all = np.copy(loh12)
    
    lu, lne, loh12 = [np.full(np.shape(lms), const.notnum) for i in range(3)]

    ind = np.where((lssfr > const.notnum) &
                   (lms > 0) &
                   (loh12_all > const.notnum) &
                   (Q > 0))
    
    # ind1 = np.where((lssfr[:,0] > const.notnum) &
    #                (lms[:,0] > 0) &
    #                (loh12[:,0] > const.notnum) &
    #                (Q[:,0] > 0))[0]
    # ind2 = np.where((lssfr[:,1] > const.notnum) &
    #                (lms[:,1] > 0) &
    #                (loh12[:,1] > const.notnum) &
    #                (Q[:,1] > 0))[0]
    
    # ind_comp = [ind1,ind2]
    
    ind_comp = []   
    for comp in range(len(Q[0])):
        ind_comp.append(np.where((lssfr[:,comp] > const.notnum) &
                       (lms[:,comp] > 0) &
                       (loh12_all[:,comp] > const.notnum) &
                       (Q[:,comp] > 0))[0])
    
    if (np.shape(ind)[1]>1):
        
        epsilon = np.full(np.shape(lssfr),const.notnum)
        cte = np.zeros(np.shape(lssfr))
        
        if origin=='sfr':
            # lu, lne, loh12 = get_une_orsi14(Q, lms, lssfr, loh12, T, q0=const.q0_orsi, z0=const.Z0_orsi, gamma=1.3)
            lu, lne, loh12 = get_une_kashino20(Q,lms,lssfr,loh12_all,T,ng_ratio,IMF_f)
            
            for comp in range(len(Q[0])):
                epsilon[:,comp][ind_comp[comp]] = ((1/alpha_B(T)) * ((4*const.c*(10**lu[:,comp][ind_comp[comp]]))/3)**(3/2) * 
                                      ((4*np.pi)/(3*Q[:,comp][ind_comp[comp]]*(10**lne[:,comp][ind_comp[comp]])))**(1/2))
                
                if ng_ratio != None:
                    epsilon[:,comp][ind_comp[comp]] = epsilon[:,comp][ind_comp[comp]] * ng_ratio
                
                cte[:,comp][ind_comp[comp]] = 3*(alpha_B(T)**(2/3)) * (3*epsilon[:,comp][ind_comp[comp]]**2*(10**lne[:,comp][ind_comp[comp]])/(4*np.pi))**(1/3) / (4*const.c)    
            
            lu[ind] = np.log10(cte[ind] * Q[ind]**(1/3))
        
        if origin=='agn':
            lne[ind] = 3
            loh12[ind] = loh12_all[ind]
            
            for comp in range(len(Q[0])):
                
                epsilon[:,comp][ind_comp[comp]] = epsilon0[ind_comp[comp]]
                
                cte[:,comp][ind_comp[comp]] = ( (3*(alpha_B(T)**(2/3)) / (4*const.c)) 
                 * (3*epsilon[:,comp][ind_comp[comp]]**2*(10**lne[:,comp][ind_comp[comp]])/(4*np.pi))**(1/3) )
                
            cte[cte==0] = 1e-50
            lu[ind] = np.log10(cte[ind] * Q[ind]**(1/3) / 3)
            lu[cte==1e-50] = const.notnum
    

    return lu, lne, loh12


def get_une(lms_o, lssfr_o, loh12_o,
            q0=const.q0_orsi, z0=const.Z0_orsi, Lagn=None, ng_ratio=None,
            Z_central_cor=False,
            gamma=1.3, T=10000, epsilon_param=[[None]], epsilon_param_z0=[[None]],
            epsilon=0.01, h0=None, IMF_f=['Kroupa','Kroupa'],
            redshift=0.1,
            unemod='kashino20', origin='sfr', verbose=True):
    '''
    Given log10(Mstar), log10(sSFR) and 12+log(O/H),
    get the ionizing parameter, U, and the electron density, ne.

    Parameters
    ----------
    lms : floats
     Masses of the galaxies per component (log10(M*) (Msun)).
    lssfr : floats
     sSFR of the galaxies per component (log10(SFR/M*) (1/yr)).
    loh12 : floats
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
    h0 : float
     If not None: value of h, H0=100h km/s/Mpc.
    IMF_f : string
     Assumed IMF.
    redshift : float
     Redshift of the input data. 
    unemod : string
     Model to go from galaxy properties to U and ne.
    origin : string
     Source of the ionizing photons.
    verbose : boolean
     Yes = print out messages.

    Returns
    -------
    Q, lu, lne, loh12 : floats
    '''

    # ncomp = len(lms[0])
    Q = phot_rate(lssfr=lssfr_o,lms=lms_o,IMF_f=IMF_f,Lagn=Lagn,origin=origin)
    
    epsilon = None
    if epsilon_param[0][0] != None:
        if origin=='agn':
            epsilon = calculate_epsilon(epsilon_param,[const.radius_NLR],h0=h0,
                              nH=const.nH_AGN,profile='exponential',verbose=verbose)
        if origin=='sfr':
            # ng = calculate_ng_hydro_eq(2*epsilon_param[1],epsilon_param[0],epsilon_param[1],profile='exponential',verbose=True)
            # epsilon = ng/const.nH_gal
            # epsilon[epsilon>1] = 1
            
            if epsilon_param_z0[0][0] != None:
                # ng_z0 = calculate_ng_hydro_eq(2*epsilon_param_z0[1],epsilon_param_z0[0],epsilon_param_z0[1],profile='exponential',verbose=True)
                # ng_ratio = n_ratio(ng,ng_z0)
                if redshift==0.8:
                    ng_ratio = const.med_to_low
                elif redshift==1.5:
                    ng_ratio = const.high_to_low
                else:
                    ng_ratio = 1.
                    
    if Z_central_cor and origin=='agn':
        loh12 = Zagn(lms_o,loh12_o)
    else:
        loh12 = np.copy(loh12_o)
    
    if unemod not in const.unemods:
        if verbose:
            print('STOP (gne_une): Unrecognised model to get U and ne.')
            print('                Possible unemod= {}'.format(const.unemods))
        sys.exit()
    elif (unemod == 'kashino20'):
        lu, lne, loh12 = get_une_kashino20(Q,lms_o,lssfr_o,loh12,T,ng_ratio,IMF_f)
    elif (unemod == 'orsi14'):
        lu, lne, loh12 = get_une_orsi14(Q,lms_o,lssfr_o,loh12,T,q0,z0,gamma,ng_ratio)
    elif (unemod == 'panuzzo03'):
        lu, lne, loh12 = get_une_panuzzo03(Q,lms_o,lssfr_o,loh12,T,epsilon,ng_ratio,origin,IMF_f)
        
    return Q, lu, lne, loh12, epsilon, ng_ratio
