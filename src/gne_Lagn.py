"""
.. moduleauthor:: Julen Expósito-Márquez <expox7@gmail.com>
.. contributions:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
"""
import numpy as np
import src.gne_const as const
from src.gne_io import read_data

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



def Rsch(Mbh):
    '''
    Given the mass of the black hole, it calculates the Schwarzschild radius.

    Parameters
    ----------
    Mbh : floats
     Mass of the black hole (Msun)
     
    Returns
    -------
    Rs : floats
    '''
    
    Rs = 2*const.G*Mbh/(const.c**2) * 1e10 * 3.086e19 #km
    return Rs


def get_Lagn(infile,cut,inputformat='hdf5',params='Lagn',AGNinputs='Lagn',
             h0=None,units_h0=False,units_Gyr=False,units_L40h2=False,
             testing=False,verbose=True):
    '''
    Calculate or get the bolometric luminosity of BHs (erg/s) 

    Parameters
    ----------
    infile : string
        Name of the input file
    cut : array of integers
        Indexes of selected galaxies for the study
    inputformat : string
        Format of the input file.
    params : array of strings
        Names of the parameters to calculate the AGN emission. 
    AGNinputs : string
        Type of calculation to obtain Lagn
    units_h0: boolean
        True if input units with h
    units_Gyr: boolean
        True if input units with */Gyr
    units_L40h2: boolean
        True if input units with 1e40erg/s
    testing : boolean
        If True only run over few entries for testing purposes
    verbose : boolean
        If True print out messages
    
    Returns
    -------
    Lagn : array of floats
        Bolometric luminosity of the BHs (erg/s)
    '''

    vals = read_data(infile,cut,inputformat=inputformat,
                                 params=params,
                                 testing=testing,verbose=verbose)
    
    if AGNinputs=='Lagn':
        Lagn = vals[0]
        if units_L40h2:
            Lagn = Lagn*1e40/h0/h0

        return Lagn # erg s^-1
    
    elif AGNinputs=='acc_rate':
        Mdot = vals[0]
        if units_h0:
            Mdot = Mdot/h0
        if units_Gyr:
            Mdot = Mdot/1e9
        Mdot = Mdot/const.yr_to_s
        
        Mbh = vals[1]
        if units_h0:
            Mbh = Mbh/h0

        if len(vals) > 2:
            spin = vals[2]
        else:
            spin = np.full(Mbh.shape,const.spin_bh)

    elif AGNinputs=='acc_stb':
        Mdot = vals[0] + vals[1]
        if units_h0:
            Mdot = Mdot/h0
        if units_Gyr:
            Mdot = Mdot/1e9
        Mdot = Mdot/const.yr_to_s
        
        Mbh = vals[2]
        if units_h0:
            Mbh = Mbh/h0

        if len(vals) > 3:
            spin = vals[3]
        else:
            spin = np.full(Mbh.shape,const.spin_bh)

    elif AGNinputs=='radio_mode':
        Mhot = vals[0]
        Mbh = vals[1]
        if units_h0:
            Mhot = Mhot/h0
            Mbh = Mbh/h0            
        
        Mdot = np.zeros(Mhot.shape)
        ind_radio = np.where((Mhot!=0)&(Mbh!=0))[0]
        Mdot[ind_radio] = acc_rate_radio(Mhot[ind_radio], 
                                         Mbh[1][ind_radio])
        
        if len(vals) > 2:
            spin = vals[2]
        else:
            spin = np.full(Mbh.shape,const.spin_bh)

    elif AGNinputs=='quasar_mode':
        M_b = vals[0]
        r_b = vals[1]
        v_b = vals[2]
        Mbh = vals[3]
        if units_h0:
            M_b = M_b/h0
            r_b = r_b/h0
            Mbh = Mbh/h0

        Mdot = np.zeros(M_b.shape)
        ind_quasar = np.where(M_b!=0)[0]
        Mdot[ind_quasar] = acc_rate_quasar(M_b[ind_quasar], 
                                           r_b[ind_quasar], v_b[ind_quasar])

        if len(vals) > 4:
            spin = vals[4]
        else:
            spin = np.full(Mbh.shape,const.spin_bh)

    elif AGNinputs=='complete': #Mbulg, rbulg, vbulg, Mhot, Mbh
        M_b = vals[0]
        r_b = vals[1]
        v_b = vals[2]
        Mhot = vals[3]
        Mbh = vals[4]
        if units_h0:
            M_b = M_b/h0
            r_b = r_b/h0
            Mhot = Mhot/h0
            Mbh = Mbh/h0

        Mdot_quasar = np.zeros(Mbh.shape)
        ind_quasar = np.where(M_b!=0)[0]
        Mdot_quasar[ind_quasar] = acc_rate_quasar(M_b[ind_quasar], 
                                                  r_b[ind_quasar],
                                                  v_b[ind_quasar])
        
        Mdot_radio = np.zeros(Mbh.shape)        
        ind_radio = np.where((Mhot!=0)&(Mbh!=0))[0]
        Mdot_radio[ind_radio] = acc_rate_radio(Mhot[ind_radio], 
                                               Mbh[ind_radio])
        
        Mdot = Mdot_radio + Mdot_quasar
        
        if len(vals) > 5:
            spin = vals[5]
        else:
            spin = np.full(Mbh.shape,const.spin_bh)
    
    Mdot_edd = acc_rate_edd(Mbh)
    
    mdot = np.zeros(Mdot.shape)
    bh = np.where(Mdot_edd!=0)[0]
    mdot[bh] = Mdot[bh]/Mdot_edd[bh]
    
    Lagn = np.zeros(np.shape(mdot))
    # n1, n2, n3, n4, n0 = [],[],[],[],[]
    for i in range(len(Lagn)):
        if mdot[i] == 0:
            # n0.append(i)
            continue
        
        # Lagn[i] = epsilon_td(const.spin_bh)*Mdot[i]*const.c**2 * (1000/const.kg_to_Msun)
        
        if mdot[i] < const.acc_rate_crit_visc:
            # n1.append(i)
            Lagn[i] = 0.2*epsilon_td(spin[i])*Mdot[i]*const.c**2*(mdot[i]/const.alpha_adaf**2)*((1-const.beta)/0.5)*(6/r_iso(spin[i])) * (1000/const.kg_to_Msun)
        elif mdot[i] < const.acc_rate_crit_adaf:
            # n2.append(i)
            Lagn[i] = 0.2*epsilon_td(spin[i])*Mdot[i]*const.c**2*(mdot[i]/const.alpha_adaf**2)*(const.beta/0.5)*(6/r_iso(spin[i])) * (1000/const.kg_to_Msun)
        elif mdot[i] < const.eta_edd:
            # n3.append(i)
            Lagn[i] = epsilon_td(spin[i])*Mdot[i]*const.c**2 * (1000/const.kg_to_Msun)
        else:
            # n4.append(i)
            Lagn[i] = const.eta_edd*Ledd(Mbh[i])*(1 + np.log(mdot[i]/const.eta_edd))
            
    # logLagn = np.log10(Lagn)
    # print(len(n1),len(n2),len(n3),len(n4))
    # print(np.mean(logLagn[n1]),np.mean(logLagn[n2]),np.mean(logLagn[n3]))
        
    return Lagn # erg s^-1
