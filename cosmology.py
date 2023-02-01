#! /usr/bin/env python

"""
Cosmology + some useful functions
This module contains various functions to compute distances and 
times in a Universe with a given cosmology.
List of functions:
  set_cosmology(): lets user specify a cosmology.
  cosmology_set(): determines wheter an input cosmology
                   has been specfied.
  report_cosmology(): report back parameters for specified 
                      cosmology.
  comoving_distance(): calculates the comoving distance at 
                       redshift, z (Mpc/h).
  redshift_at_distance(): calculates the redshift at comoving 
                          disance, r.
  age_of_universe(): calculates the age of the Universe at 
                     redshift, z (Gyr).
  lookback_time(): calculates lookback time to given redshift, 
                   z (Gyr).
  angular_diamater_distance(): calculates the angular diameter
                               distance a redshift, z (Mpc/h).
  angular_scale(): calculates the angular scale at redshift, z.
  luminosity_distance(): calculates the luminosity distance at 
                         redshift, z (Mpc/h).
  comoving_volume(z): calculates the comoving volume contained
                    within a sphere extending out to redshift,
                    z ((Mpc/h)^3).
  cv_survey(z1,z2,area): calculates the comoving volume of a survey ((Mpc/h)^3).
  dVdz() :  calculates dV/dz at redshift, z  (Mpc/h)^3.
  distance_modulus(z): DM = 5log10(Dl/10) - 5logh; Dl expected in Mpc/h
  band_corrected_distance_modulus(z): DM for Galform
  H(): return Hubble constant as measured at redshift, z.
  tHubble(): returns Hubble time at redshift z (Gyr).
  E(): returns Peebles' E(z) function at redshift, z, for
       specified cosmology.
  Hubble(): return Hubble parameter for specified cosmology.
  Omega_M(): return matter density for specified cosmology.
  Omega_b(): return baryon density for specified cosmology.
  Omega_V(): return vacuum density for specified cosmology.
  Omega_r(): return radiation density for specified cosmology.
  Omega_k(): return curvature for specfied cosmology.
  omegam(z): return matter density at z.
  omegab(z): return baryon density at z.
  omegav(z): return vacuum density at z.
  omegar(z): return radiation density at z.
  ndeg2nV(ndeg,z1,z2,area,verbose=False): Transforms number of objects per deg2 to number density (N/V). 
  emission_line_flux(luminosity_data,z): returns flux from luminosity
  emission_line_luminosity(flux_data,z): returns luminosity from flux
  polar2cartesians(ra,dec,zz): returns cartesian coordinates (cMpc/h)
NOTE: this module requires the numpy and scipy libraries to be
      available for import!
Based upon the 'Cosmology Calculator' (Wright, 2006, PASP, 
118, 1711) and Fortran 90 code written by John Helly.
"""

import sys
import numpy as np
import scipy as sp
from scipy.constants import c,constants
from scipy.integrate import romberg
from astropy.constants import M_sun

WM = None
WV = None
WB = None
WR = None
WK = None
h = None

kmpersec_to_mpchpergyr = None

zmax = 20.0
dz = 0.0005
nzmax = int(zmax/dz)
r_comoving = np.zeros(nzmax)
redshift = np.arange(0.0,zmax,dz)
inv_dz = 1.0/dz

zlow_lim = 0.001

Mpc = constants.mega*constants.parsec #10**6*3.08567758131e+16
H100 = 100.0*constants.kilo/Mpc # in h/s units
Gyr = constants.giga*constants.year
invH0 = (Mpc/(100.0*constants.kilo))/Gyr
DH = c/1000./100. # Hubble Distance in Mpc/h (c is in m/s)
Mpc2cm = constants.mega*constants.parsec*100.
zlow = 0.00001 ; dlz = np.log(zmax)/float(nzmax) 
lredshift = np.arange(np.log(zlow),np.log(zmax),dz)

asky =  4.0*np.pi*(180/np.pi)**2

def set_cosmology(omega0=None,omegab=None,lambda0=None,h0=None, \
                      universe="Flat",include_radiation=False):
    """
    set_cosmology(): Sets the cosmological parameters and evaluates
                     the comoving distance relation as a function
                     of redshift.
                     
    USAGE: set_cosmology([Omega_M],[Omega_b],[Omega_V],[h],
                         [universe=Flat],[include_radion=True])
                         
           Omega_M: total matter density at z=0
           Omega_b: baryon matter density at z=0
           Omega_V: vacuum energy density at z=0
                    (default value depends on type of universe;
                    1-(Omega_M+Omega_R) for flat universe, 0 for 
                    open universe)
           h: Hubble parameter at z=0
           universe: specifies desired geomerty of Universe -- only
                     considered if Omega_V not specified
                     (can be "Flat" or "Open"; default value is "Flat")
           include_radiation: include Omega_R in calculations?
                     -- can be "True" (i.e. set Omega_R = 4.165e-5/(h*h))
                     of "False" (i.e. set Omega_R = 0.0)
                     (default value is False)
          Default values: Planck18
    """

    global WM, WV, WB, WR, WK, h
    if(h0 is None):
        h = 0.674
    else:
        h = h0    
    if(include_radiation):
        WR = 8.985075e-5
    else:
        WR = 0.0
    if(omegab is None):
        WB = 0.0224/(h*h)
    else:
        WB = omegab
    if(omega0 is None):
        WM = 0.315
    else:
        WM = omega0
    if(lambda0 is None):
        if(universe in ("Flat","F","flat","f")):
            WV = 1.0 - (WM + WR)
        if(universe in ("Open","O","open","o")):
            WV = 0
    else:
        WV = lambda0
    WK = 1.0 - (WM + WV + WR)

    global r_comoving, redshift
    for i in range(1,len(redshift)):
        z1 = redshift[i-1]
        z2 = redshift[i]
        r_comoving[i] = r_comoving[i-1] + romberg(f,z1,z2)

    global kmpersec_to_mpchpergyr
    kmpersec_to_mpchpergyr = constants.kilo * (Gyr/Mpc) * h

    return


def set_Millennium():
    set_cosmology(0.25,0.045,0.75,0.73)
    return


def set_MR7():
    set_cosmology(0.272,0.0455,0.728,0.704)
    return


def set_bahamasW9():
    set_cosmology(0.2793, 0.0463, 0.7207, 0.700)
    return


def set_bahamasP():
    set_cosmology(0.3175, 0.0490, 0.6825, 0.6711)
    return


def set_Planck13():
    set_cosmology(0.3086,0.0483,0.6914,0.6777)
    return


def set_Planck15(): #UNIT
    set_cosmology(0.3089,0.0486,0.6911,0.6774)
    return


def cosmology_set():
    """
    cosmology_set(): determines whether an input cosmology
                     has been specified (yes ==> TRUE,
                     no ==> FALSE).
    USAGE: cosmology_set()
    """
    if(WM is None):
        return False
    else:
        return True
    

def report_cosmology():
    """
    report_cosmology(): reports parameters for inputted cosmology
    USAGE: report_comology()
    """
    print("***********************")
    print("COSMOLOGY:")
    print("   Omega_M = {0:5.3f}".format(WM))
    print("   Omega_b = {0:5.3f}".format(WB))
    print("   Omega_V = {0:5.3f}".format(WV))
    print("   h       = {0:5.3f}".format(h))
    print("   Omega_R = {0:5.3e}".format(WR))
    print("   Omega_k = {0:5.3f}".format(WK))
    print("***********************")
    return


def f(z):
    """
    f(z): Function relating comoving distance to redshift.
          Integrating f(z)dz from 0 to z' gives comoving
          distance r(z'). Result is in Mpc/h.
          
          Note: uses global cosmology variables.          
    """
    a = 1.0/(1.0+z)
    result = WK*np.power(a,-2) + WV + \
        WM*np.power(a,-3) + WR*np.power(a,-4)
    result = DH/np.sqrt(result)
    return result

def E(z):
    """
    E(z): Peebles' E(z) function.
              
          Note: uses global cosmology variables.  
    """
    a = 1.0/(1.0+z)
    result = WK*np.power(a,-2) + WV + \
             WM*np.power(a,-3) + WR*np.power(a,-4)
    return np.sqrt(result)


def rez(lz):
    """
    E(ln_z): Function relating comoving distance to redshift.
          Integrating rez(z)d(ln_z) from zlow to z' gives comoving
          distance r(z'). Result is in Mpc/h.
          
          Note: uses global cosmology variables.          
    """
    
    z = exp(lz)
    a = 1.0/(1.0+z)
    result = WK*np.power(a,-2) + WV + \
        WM*np.power(a,-3) + WR*np.power(a,-4)

    result = DH/np.sqrt(result)
    return result



def H(z):
    """
    H(z): Function to return the Hubble parameter as measured
           by an observer at redshift, z.
    """
    cosmology_set()
    result = 100.0*E(z)
    return result

def tHubble(z):
    """
    tHubble(z): Function to return the Hubble time at z in Gyr.
    """
    cosmology_set()
    result = 1./(H(z)*kmpersec_to_mpchpergyr)
    return result

def comoving_distance(z):
    """
    comoving_distance(): returns the comoving distance (in Mpc/h)
                         corresponding to redshift, z.
                       
    USAGE: r = comoving_distance(z)
    NOTE: requires that a cosmology must first have been
          set using set_cosmology()
    """
    cosmology_set()
    r = np.interp(z,redshift,r_comoving)
    return r

def redshift_at_distance(r):
    """
    redshift_at_distance(): returns the redshift corresponding
                            to comoving distance, r (in Mpc/h).
    USAGE: z = redshift_at_distance(z)
    NOTE: requires that a cosmology must first have been
          set using set_cosmology()
    """
    cosmology_set()
    z = np.interp(r,r_comoving,redshift)
    return z


def age_of_universe(z):
    """
    age_of_universe(): returns the age of the Universe (in Gyr) at
                       a redshift, z, for the given cosmology.
    USAGE: age = age_of_universe(z)
    NOTE: requires that a cosmology must first have been
          set using set_cosmology()
    """
    cosmology_set()
    a = 1.0/(1.0+z)
    if(WM >= 0.99999): # Einstein de Sitter Universe
        result = invH0*2.0*np.sqrt(a)/(3.0*h)
    else:
        if(WV <= 0.0): # Open Universe
            zplus1 = 1.0/a
            result1 = (WM/(2.0*h*np.power(1-WM,1.5)))
            result2 = 2.0*np.sqrt(1.0-WM)*np.sqrt(WM*(zplus1-1.0)+1.0)
            result3 = np.arccosh((WM*(zplus1-1.0)-WM+2.0)/(WM*zplus1))
            result = invH0*result1*(result2/result3)
        else: # Flat Universe with non-zero Cosmological Constant
            result1 = (2.0/(3.0*h*np.sqrt(1.0-WM)))
            result2 = np.arcsinh(np.sqrt((1.0/WM-1.0)*a)*a)
            result = invH0*result1*result2
    return result


def lookback_time(z):
    """
    lookback_time(): returns the lookback time (in Gyr) to 
                     redshift, z.
    USAGE: t = lookback_time(z)
    
    NOTE: requires that a cosmology must first have been
          set using set_cosmology()    
    """
    cosmology_set()
    t = age_of_universe(0.0) - age_of_universe(z)
    return t


def angular_diameter_distance(z):
    """
    angular_diameter_distance(): returns the angular diameter
                                 distance (in Mpc/h) corresponding
                                 to redshift, z.
                                 Da = size/rad
    USAGE: dA = angular_diameter_distance(z)    
    NOTE: requires that a cosmology must first have been
          set using set_cosmology()    
    """
    cosmology_set()
    dr = comoving_distance(z)*Mpc/(c/H100) #Unitless
    x = np.sqrt(np.abs(WK))*dr
    if np.ndim(x) > 0:
        ratio = np.ones_like(x)*-1.00
        mask = (x > 0.1)
        y = x[np.where(mask)]
        if(WK > 0.0):
            np.place(ratio,mask,0.5*(np.exp(y)-np.exp(-y))/y)
        else:
            np.place(ratio,mask,np.sin(y)/y)
        mask = (x <= 0.1)
        y = np.power(x[np.where(mask)],2)
        if(WK < 0.0): 
            y = -y
        np.place(ratio,mask,1.0 + y/6.0 + np.power(y,2)/120.0)
    else:        
        ratio = -1.0
        if(x > 0.1):
            if(WK > 0.0):
                ratio = 0.5*(np.exp(x)-np.exp(-x))/x
            else:
                ratio = np.sin(x)/x
        else:
            y = np.power(x,2)
            if(WK < 0.0): 
                y = -y
            ratio = 1.0 + y/6.0 + np.power(y,2)/120.0
    dt = ratio*dr/(1.0+z)
    dA = (c/H100)*dt/Mpc
    return dA


def angular_scale(z):
    """
    angular_scale(): returns the angular scale (in kpc/h/arcsec)
                     corresponding to redshift, z.
                   
    USAGE: a = angular_scale(z)
    NOTE: requires that a cosmology must first have been
          set using set_cosmology()    
    """
    cosmology_set()
    da = angular_diameter_distance(z) #Mpc/h/rad
    a = da/206.26480 # 1 rad = 206265 arcsec
    return a
    

def luminosity_distance(z):
    """
    luminosity_distance(): returns the luminosity distance
                           (in Mpc/h) corresponding to a
                           redshift, z.
    USAGE: dL = luminosity_distance(z)
    NOTE: requires that a cosmology must first have been
          set using set_cosmology()    
    """
    dL = angular_diameter_distance(z)*(1.0+z)**2
    return dL
    

def comoving_volume(z, verbose=False):
    """
    comoving_volume(): returns the comoving volume (in (Mpc/h)^3)
                       contained within a sphere extending out
                       to redshift, z.
    Example:
    import Cosmology as cosmo
    cosmo.set_Planck15()
    cosmo.comoving_volume(0.9,verbose=True)
    > cV (z=0.9) = 4.03e+10 (Mpc/h)^-3
    """
    cosmology_set()

    if (z<zlow_lim):
        return 0.0
    
    dr = comoving_distance(z)*Mpc/(c/H100) #Unitless: DC/DH
    x = np.sqrt(np.abs(WK))*dr
    if np.ndim(z) > 0:
        ratio = np.ones_like(z)*-1.0
        mask = (x > 0.1)
        y = x[np.where(mask)]
        if(WK > 0.0):
            rat = (0.125*(np.exp(2.0*y)-np.exp(-2.0*y))-y/2.0)
        else:
            rat = (y/2.0 - np.sin(2.0*y)/4.0)
        np.place(ratio,mask,rat/(np.power(y,3)/3.0))
        mask = (x <= 0.1)
        y = np.power(x[np.where(mask)],2)
        if(WK < 0.0): 
            y = -y
        np.place(ratio,mask,1.0 + y/5.0 + np.power(y,2)*(2.0/105.0))
    else:  
        ratio = -1.0
        if(x > 0.1):
            if(WK > 0.0):
                ratio = (0.125*(np.exp(2.0*x)-np.exp(-2.0*x))-x/2.0)
            else:
                ratio = (x/2.0 - np.sin(2.0*x)/4.0)
            ratio = ratio/(np.power(x,3)/3.0)
        else:
            y = np.power(x,2)
            if(WK < 0.0): 
                y = -y
            ratio = 1.0 + y/5.0 + np.power(y,2)*(2.0/105.0)

    vol = 4.0*np.pi*ratio*np.power((c/H100)*dr/Mpc,3)/3.0

    if verbose:
        print('cV (z={:.1f}) = {:.4e} (Mpc/h)^-3'.format(z,vol))
    return vol


def cv_survey(z1,z2,area,verbose=False):
    '''
    Get the volume of a survey.
    A cosmology needs to have been set.
    Parameters:
    z1 : float, lower redshift limit of the survey
    z2 : float, higher redshift limit of the survey
    area : float, total area of the survey (deg2)
    Returns:
    vsurvey : survey's comoving volume [(Mpc/h)^-3]
    Example:
    import Cosmology as cosmo
    cosmo.set_bahamasW9()
    cosmo.cv_survey(0.9,1.0,14000,verbose=True)
    > V survey (dz=0.1) = 3.9e+09 (Mpc/h)^-3
    '''

    if (z1<zlow_lim):
        dV = comoving_volume(z2)
    else:
        dV = comoving_volume(z2) - comoving_volume(z1)
        
    vsurvey = dV*area/asky
    
    if verbose:
        print('V survey (dz={:.1f}) = {:.4e} (Mpc/h)^-3'.format(z2-z1,vsurvey))    
    return vsurvey


def dVdz(z):
    """
    dVdz() : returns the comoving volume element dV/dz
             at redshift, z, for all sky (Mpc/h)^3.
             dV = (c/H100)*(1+z)**2*D_A**2/E(z) dz dOmega
             f(z) = (c/H100)/E(z)
             ==> dV/dz(z,all sky) = 4*PI*f(z)*(1+z)**2*D_A**2
             
    USAGE: dVdz = dVdz(z)
    NOTE: requires that a cosmology must first have been
          set using set_cosmology()         
    """
    cosmology_set()
    dA = angular_diameter_distance(z)
    return f(z)*np.power(dA,2)*np.power(1.0+z,2)*4.0*np.pi
    

def distance_modulus(z):
    '''
    Dinstance modulus 5log10(Dl/10) - 5logh
    Dl is expected in Mpc/h
    '''
    cosmology_set()
    if (z < zlow_lim):
        dm = 0.0
    else:
        dL = luminosity_distance(z)
        dm = 5.0*np.log10(dL) + 25
    return dm


def band_corrected_distance_modulus(z):
    """
    band_corrected_distance_modulus(): returns the Band Corrected
            Distance Modulus (BCDM) at redshift, z.
    USAGE: bcdm = band_corrected_distance_modulus(z)
            
    NOTE: requires that a cosmology must first have been
          set using set_cosmology()
    FURTHER INFORMATION:
    There is no h dependence as we work always in length units of Mpc/h 
    such that our absolute magnitudes are really Mabs-5logh and no 
    additional h dependence is needed here to get apparent magnitudes 
    that are h independent.
    
    In Galform versions 2.5.1 onwards the additional -2.5 * log10(1.0+z)
    is needed to convert from absolute to apparent magnitude as the 
    definition of absolute magnitude in the Galform code has been changed
    by a factor of (1+z). With the new definition a galaxy with a SED in 
    which f_nu is a constant will, quite sensibly, have the same AB 
    absolute magnitude independent wave band range (including whether it 
    is rest or observer frame) and independent of redshift. 
    
    One way of thinking about this is that while the standard luminosity
    distance and corresponding distance modulus applies to bolometric 
    luminosities, for a filter of finite width the flux depends on the
    band width of the filter in the galaxy's rest frame and it is this 
    that we are taking into account when defining this "band corrected"
    distance modulus. 
    """
    cosmology_set()
    if (z < zlow_lim):
        bcdm = 0.0
    else:
        dref = 10.0/constants.mega # 10pc in Mpc
        dL = luminosity_distance(z)
        bcdm = 5.0*np.log10(dL/dref) - 2.5*np.log10(1.0+z)
    return bcdm


def Hubble():
    """
    Hubble(): returns h for the specified cosmology
              (will exit if no cosmology has been set)
    
    USAGE: h = Hubble()
    """
    cosmology_set()
    return h

def Omega_M():
    """
    Omega_M(): returns Omega_M for the specified cosmology
              (will exit if no cosmology has been set)
    
    USAGE: wm = Omega_M()
    """
    cosmology_set()
    return WM


def Omega_b():
    """
    Omega_b(): returns Omega_b for the specified cosmology
               (will exit if no cosmology has been set)
                   
    USAGE: wb = Omega_b()
    """
    cosmology_set()
    return WB

def Omega_V():
    """
    Omega_V(): returns Omega_V for the specified cosmology
               (will exit if no cosmology has been set)
                   
    USAGE: wv = Omega_V()
    """
    cosmology_set()
    return WV

def Omega_r():
    """
    Omega_r(): returns Omega_r for the specified cosmology
               (will exit if no cosmology has been set)
    USAGE: wr = Omega_r()
    """
    cosmology_set()
    return WR

    
def Omega_k():
    """
    Omega_k(): returns Omega_k for the specified cosmology
               (will exit if no cosmology has been set)
                   
    USAGE: wk = Omega_k()
    """
    cosmology_set()
    return WK


def omegam(z):
    """
    Matter density at z
              
          Note: uses global cosmology variables.  
    """
    a = 1.0/(1.0+z)
    omegam=WM*np.power(a,-3)/(E(z)**2)
    return omegam

def omegab(z):
    """
    Baryonic density at z
              
          Note: uses global cosmology variables.  
    """
    a = 1.0/(1.0+z)
    omegab=WB*np.power(a,-3)/(E(z)**2)
    return omegab


def omegav(z):
    """
    Vacuum/Cosmological constant density at z
              
          Note: uses global cosmology variables.  
    """
    a = 1.0/(1.0+z)
    omegav=WV/(E(z)**2)
    return omegav


def omegar(z):
    """
    Radiation density at z
              
          Note: uses global cosmology variables.  
    """
    a = 1.0/(1.0+z)
    omegar=WR*np.power(a,-4)/(E(z)**2)
    return omegar


def kaiser_factor(z,bias,gamma=None):
    """
    Calculate the Kaiser Factor from either a linear bias value
    or an array bias 
    """
    if (gamma is None):
        gamma = 0.55
    
    omb = np.power(omegam(z),gamma)/bias 
    kaiser_factor = 1. + (2./3.)*omb + (1./5.)*omb**2.

    return kaiser_factor


def ndeg2nV(ndeg,z1,z2,verbose=False):
    '''
    Transforms number of objects per deg2 per dz to 
    number density (N/V). 
    A cosmology needs to have been set.
    Parameters:
    ndeg : float, number of objects per deg2 per dz
    z1 : float, lower redshift limit of the survey
    z2 : float, higher redshift limit of the survey
    Returns:
    nV : Number density (NV) [(Mpc/h)^-3]
    Example:
    import Cosmology as cosmo
    cosmo.set_bahamasW9()
    cosmo.ndeg2nV(2400,0.6,1.6)
    > 0.0002671477226063551
    '''

    dV = comoving_volume(z2) - comoving_volume(z1)
    dz = z2-z1
    nV = ndeg*dz*asky/dV
    
    if verbose:
        print('n (dz={:.1f}) = {:.4e} (Mpc/h)^-3'.format(dz,nV))
    
    return nV

    
def logL2flux(log10luminosity,z):
    """
    Returns flux in units of erg/s/cm^2 from input of
    log10(Luminosity in units of h-2erg/s)
    and corresponding redshifts.
    """
    if log10luminosity>-9.:
        # Luminosity distance in cm/h
        d_L = max(luminosity_distance(z),10.**-5)*Mpc2cm

        # Luminosities are in h-2 erg/s units
        den = 4.0*np.pi*(d_L**2)
        emission_line_flux = log10luminosity - np.log10(den)
        # Flux in erg/s/cm^2
        emission_line_flux = 10**(emission_line_flux)
    else:
        emission_line_flux = 0.

    return emission_line_flux


def emission_line_flux(luminosity_data,z):
    """
    Returns flux in units of erg/s/cm^2 from input of 
    luminosity_data in units of E+40*h-2erg/s and corresponding redshifts.
    """

    if luminosity_data>0.:
        # Luminosity distance in cm/h
        d_L = max(luminosity_distance(z),10.**-5)*Mpc2cm

        # Luminosities are in 10^40 h-2 erg/s units
        den = 4.0*np.pi*(d_L**2)
        emission_line_flux = np.log10(luminosity_data/den) + 40.
        # Flux in erg/s/cm^2
        emission_line_flux = 10**(emission_line_flux)
    else:
        emission_line_flux = 0.
        
    return emission_line_flux

def emission_line_luminosity(flux_data,z):
    """
    Returns luminosity in units of E+40*h-2erg/s from input of 
    flux_data in units of erg/s/cm^2 and corresponding redshifts.
    """

    if flux_data>0.:
        # Luminosity distance in cm/h
        d_L = max(luminosity_distance(z),10.**-5)*Mpc2cm
        
        #emission_line_luminosity = np.log10(4.0*np.pi*(d_L**2)*flux_data) - 40.
        print('WARNING: unsure (1+z) F to L')
        emission_line_luminosity = np.log10((1+z)*4.0*np.pi*(d_L**2)*flux_data) - 40. 
        emission_line_luminosity = 10**(emission_line_luminosity)
    else:
        emission_line_luminosity = 0.
        
    return emission_line_luminosity


def polar2cartesians(ra,dec,zz):
    """ Returns cartesian coordinates given polar ones in degrees"""
    dz = comoving_distance(zz)
    cx = dz*np.cos(dec*(np.pi/180.))*np.cos(ra*(np.pi/180.))
    cy = dz*np.cos(dec*(np.pi/180.))*np.sin(ra*(np.pi/180.))
    cz = dz*np.sin(dec*(np.pi/180.))
    return cx,cy,cz

if __name__== "__main__":
    #set_Planck13()
    set_cosmology(0.2,0.0483,0.8,0.6777)
    zz = 1.
    print(distance_modulus(zz))