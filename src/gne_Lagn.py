"""
.. moduleauthor:: Julen Expósito-Márquez <expox7@gmail.com>
.. contributions:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
"""
import numpy as np
import src.gne_const as c
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

    lms[:,1] = c.notnum
    lms[:,1][ind] = np.log10(Lagn_param[-1][ind])


def get_Ledd(Mbh): # Eddington luminosity
    '''
    Calculate the Eddington luminosity for a black hole of mass Mbh,
    following Eq.3 in Griffin+2019

    Parameters
    ----------
    Mbh : array of floats
       Mass of the black hole (Msun).
     
    Returns
    -------
    Ledd : array of floats
    '''
    
    Ledd = 1.26e38*Mbh # erg/s
    return Ledd # erg/s
    

def acc_rate_edd(Mbh): # Eddington mass accretion rate
    '''
    Calculate the Eddington mass accretion rate of a black hole,
    following Eq.4 in Griffin+2019
    (or 14.9 in the book from Mo, van den Bosch and White)

    Parameters
    ----------
    Mbh : array of floats
         Mass of the black hole (Msun).
     
    Returns
    -------
    acc_rate : array of floats
    '''

    acc_rate = (1e-7*c.yr_to_s/c.Msun) * get_Ledd(Mbh)/(c.e_r_agn*c.c*c.c)

    return acc_rate # Msun/yr


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
    
    dyn_time = r_bulge/v_bulge * 1e-5*c.Mpc_to_cm # s
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
    
    acc_rate = M_bulge*c.fbh/(t_bulge(r_bulge, v_bulge)*c.fq) # Msun/s
    return acc_rate


def acc_rate_radio(Mhot, Mbh, kagn=c.kagn, kagn_exp=c.kagn_exp):
    '''
    Given the mass of the hot gas and the mass of the black hole,
    calculates the accretion rate of the radio (hot gas) mode,
    following Eq. 1 in Henriques+2016.
    Default kagn values are derived from a GP20 sample 

    Parameters
    ----------
    Mhot : array of floats
        Mass of the hot gas (Msun).
    Mbh : array of floats
        Mass of the black hole (Msun).
    kagn : float
        Multiplicative factor for units Msun/yr
    kagn_exp : float
        Exponent for the dependence of Mdot with Mbh*Mhot
    
    Returns
    -------
    acc_rate : array of floats
    '''
    
    acc_rate = kagn*(Mhot*Mbh*1e-19)**kagn_exp # Msun/yr

    return acc_rate # Msun/yr


def get_Lagn_M16(Mdot):
    '''
    Calculate the AGN bolometric luminosity
    following Sec.7 and Table 1 in McCarthy+2016

    Parameters
    ----------
    Mdot : array of floats
         Accretion rate onto the black hole (Msun/yr).
     
    Returns
    -------
    LagnM16 : array of floats
    '''

    LagnM16 = c.e_r_agn*(1-c.e_f_agn)*Mdot*c.c*c.c * (1e7*c.Msun/c.yr_to_s)

    return LagnM16 #erg/s


def get_Lagn_H14(Mdot,Mbh):
    '''
    Calculate the AGN bolometric luminosity
    following Sec.4.1 from Hirschmann+2014, and McCarthy+16

    Parameters
    ----------
    Mdot : array of floats
         Accretion rate onto the black hole (Msun/yr).
    Mbh : array of floats
       Mass of the black hole (Msun).
     
    Returns
    -------
    LagnH14 : array of floats
    '''

    Mdot_edd = acc_rate_edd(Mbh)
    fedd = Mdot/Mdot_edd

    if isinstance(fedd, (float, int)): # Floats
        if fedd>0.1:
            LagnH14 = get_Lagn_M16(Mdot)
        else:
            LagnH14 = 10.*get_Ledd(Mbh)*(Mdot/Mdot_edd)**2
    
    else: # Arrays
        LagnH14 = 10.*get_Ledd(Mbh)*(Mdot/Mdot_edd)**2
        LagnH14[fedd>0.1] = get_Lagn_M16(Mdot[fedd>0.1])

    return LagnH14 #erg/s



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
    
    Rs = 2*c.G*Mbh/(c.c_cm**2) * 1e10 * 3.086e19 #km
    return Rs


def get_Lagn(infile,cut,inputformat='hdf5',params='Lagn',AGNinputs='Lagn',
             h0=None,units_h0=False,units_Gyr=False,units_L40h2=False,
             kagn=c.kagn,kagn_exp=c.kagn_exp,testing=False,verbose=True):
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
    kagn : float
        Multiplicative factor for units Msun/yr, radio mode
    kagn_exp : float
        Exponent for the dependence of Mdot with Mbh*Mhot, radio mode
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
        return Lagn # erg/s
    
    elif AGNinputs=='Mdot_hh':
        Mdot = vals[0]
        Mbh = vals[1]
        if units_h0:
            Mdot = Mdot/h0
            Mbh = Mbh/h0
        if units_Gyr:
            Mdot = Mdot/1e9

        if len(vals) > 2:
            spin = vals[2]
        else:
            #spin = np.full(Mbh.shape,c.spin_bh)
            Lagn = get_Lagn_H14(Mdot,Mbh)
            return Lagn # erg/s

    elif AGNinputs=='Mdot_stb_hh':
        Mdot = vals[0] + vals[1]
        Mbh = vals[2]
        if units_h0:
            Mdot = Mdot/h0
            Mbh = Mbh/h0
        if units_Gyr:
            Mdot = Mdot/1e9
        
        if len(vals) > 3:
            spin = vals[3]
        else:
            #spin = np.full(Mbh.shape,c.spin_bh)
            Lagn = get_Lagn_H14(Mdot,Mbh)
            return Lagn # erg/s

    elif AGNinputs=='radio_mode':
        Mhot = vals[0]
        Mbh = vals[1]
        if units_h0:
            Mhot = Mhot/h0
            Mbh = Mbh/h0            

        Mdot = np.zeros(Mhot.shape)
        ind = np.where((Mhot!=0)&(Mbh!=0))
        if (np.shape(ind)[1]>0):
            Mdot[ind] = acc_rate_radio(Mhot[ind],Mbh[ind],
                                       kagn=kagn, kagn_exp=kagn_exp)

        if len(vals) > 2:
            spin = vals[2]
        else:
            Lagn = get_Lagn_H14(Mdot,Mbh)
            return Lagn # erg/s
            
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
        Mdot[M_b>0] = acc_rate_quasar(M_b[M_b>0],r_b[M_b>0], v_b[M_b>0])

        if len(vals) > 4:
            spin = vals[4]
        else:
            #spin = np.full(Mbh.shape,c.spin_bh)
            Lagn = get_Lagn_H14(Mdot,Mbh)
            return Lagn # erg/s

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
        Mdot_quasar[M_b>0] = acc_rate_quasar(M_b[M_b>0], r_b[M_b>0],
                                             v_b[M_b>0])
        
        Mdot_radio = np.zeros(Mbh.shape)
        ind = np.where((Mhot>0)&(Mbh>0))
        if (np.shape(ind)[1]>0):
            Mdot_radio[ind] = acc_rate_radio(Mhot[ind],Mbh[ind])
        
        Mdot = Mdot_radio + Mdot_quasar
        
        if len(vals) > 5:
            spin = vals[5]
        else:
            #spin = np.full(Mbh.shape,c.spin_bh)
            Lagn = get_Lagn_H14(Mdot,Mbh)
            return Lagn # erg/s

    Mdot = Mdot/c.yr_to_s ####here to check units w spin
    Mdot_edd = acc_rate_edd(Mbh) ###here units in what follow might be bad as now this is in Msun/yr
    
    mdot = np.zeros(Mdot.shape)
    bh = np.where(Mdot_edd!=0)[0]
    mdot[bh] = Mdot[bh]/Mdot_edd[bh]
    
    Lagn = np.zeros(np.shape(mdot))
    # n1, n2, n3, n4, n0 = [],[],[],[],[]
    for i in range(len(Lagn)):
        if mdot[i] == 0:
            # n0.append(i)
            continue
        
        # Lagn[i] = epsilon_td(c.spin_bh)*Mdot[i]*c.c_cm**2 * (1000/c.kg_to_Msun)
        
        if mdot[i] < c.acc_rate_crit_visc:
            # n1.append(i)
            Lagn[i] = 0.2*epsilon_td(spin[i])*Mdot[i]*c.c_cm**2*(mdot[i]/c.alpha_adaf**2)*((1-c.beta)/0.5)*(6/r_iso(spin[i])) * (1000/c.kg_to_Msun)
        elif mdot[i] < c.acc_rate_crit_adaf:
            # n2.append(i)
            Lagn[i] = 0.2*epsilon_td(spin[i])*Mdot[i]*c.c_cm**2*(mdot[i]/c.alpha_adaf**2)*(c.beta/0.5)*(6/r_iso(spin[i])) * (1000/c.kg_to_Msun)
        elif mdot[i] < c.eta_edd:
            # n3.append(i)
            Lagn[i] = epsilon_td(spin[i])*Mdot[i]*c.c_cm**2 * (1000/c.kg_to_Msun)
        else:
            # n4.append(i)
            Lagn[i] = c.eta_edd*get_Ledd(Mbh[i])*(1 + np.log(mdot[i]/c.eta_edd))
            
    # logLagn = np.log10(Lagn)
    # print(len(n1),len(n2),len(n3),len(n4))
    # print(np.mean(logLagn[n1]),np.mean(logLagn[n2]),np.mean(logLagn[n3]))
        
    return Lagn # erg s^-1
