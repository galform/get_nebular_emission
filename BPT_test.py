#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 10:02:26 2023

@author: expox7
"""

import numpy as np
import h5py
from matplotlib import pyplot as plt
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_photio import get_limits
import sys
from cosmology import set_cosmology
import mpl_style

# File selection
infile = ["output_data/emlines_GP20_z0.0_Kashino_test.hdf5"]
file = 0

# Flux, Luminosity, Mass, Metallicity, Ratio, sSFR
color_labels = [r'$\log H_\alpha \ [\rm erg/s]$', r'$\log L[H_\alpha] \ [\rm erg/s]$',
                r'$\log M_* \ [M_\odot]$', r'$\log Z$',
                r'$\frac{L_{\rm H_\alpha, AGN}}{L_{\rm H_\alpha}}$', r'$\log sSFR \ [yr^{-1}]$']
color = 4

# NII, SII
plot = 0

# Redshift of the galaxies from the sample. Selection criteria stablished for
# z = 0.1, 0.8, 1.5, following Exposito-Marquez et. al. 2024.
redshift = 0.1

# Minimum ratio between luminosity in Halpha from AGNs and total luminosity.
ratio_min = -0.1

# Consider attenuation?
att = False

##########################################

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
    
    if BPT not in const.BPT_lines:
        print('STOP (eml_plots.lines_BPT): ',
              'BPT plot not recognized.')
        sys.exit()
        return
    elif BPT=='NII':
        if line=='SFR_Composite':
            boundary = 0.61/(x-0.05) + 1.3 # Kauffmann 2003
        elif line=='Composite_AGN':
            boundary = 0.61/(x-0.47) + 1.19 # Kewley 2001
        elif line=='LINER_NIIlim':
            boundary = np.log10(0.6) # Kauffmann 2003
        elif line=='LINER_OIIIlim':
            boundary = np.log10(3) # Kauffmann 2003
    elif BPT=='SII':
        if line=='SFR_AGN':
            boundary = 0.72/(x-0.32) + 1.3 # Kewley 2001
        elif line=='Seyfert_LINER':
            boundary = 1.89*x + 0.76 # Kewley 2006
    elif BPT=='SII':
        if line=='SFR_AGN':
            boundary = 0.73/(x-0.59) + 1.33 # Kewley 2001
        elif line=='Seyfert_LINER':
            boundary = 1.18*x + 1.3 # Kewley 2006
            
    return boundary

plt.style.use(mpl_style.style1)

set_cosmology(omega0=const.omega0, omegab=const.omegab,lambda0=const.lambda0,h0=const.h)

titles = []
for i in range(100):
    titles.append(str(i))
    
X = np.loadtxt('observational_data/Main-ELG.txt',skiprows=1).T

Ha_favole = X[6]
Hb_favole = X[9]
NII6548_favole = X[18]
OII3727_favole = X[12]
OIII5007_favole = X[15]
SII6731_favole = X[21]

#############################################################

plt.figure(figsize=(15,15))

f = h5py.File(infile[file], 'r')
data = f['data']

lu_disk = data['lu_sfr'][:,0]
lne_disk = data['lne_sfr'][:,0]
loh12_disk = data['lz_sfr'][:,0]

lu_agn = data['lu_agn'][:,0]

minU, maxU = get_limits(propname='U', photmod='gutkin16')
minnH, maxnH = get_limits(propname='nH', photmod='gutkin16')
minZ, maxZ = get_limits(propname='Z', photmod='gutkin16')

minU, maxU = float(minU), float(maxU)
minnH, maxnH = float(minnH), float(maxnH)
minZ, maxZ = float(minZ), float(maxZ)

minU_agn, maxU_agn = get_limits(propname='U', photmod='feltre16')

minU_agn, maxU_agn = float(minU_agn), float(maxU_agn)

if att:
    Hbeta = np.sum(data['Hbeta_sfr_att'],axis=0)
    OIII5007 = np.sum(data['OIII5007_sfr_att'],axis=0)
    NII6548 = np.sum(data['NII6584_sfr_att'],axis=0)
    Halpha = np.sum(data['Halpha_sfr_att'],axis=0)
    SII6731 = np.sum(data['SII6717_sfr_att'],axis=0) + np.sum(data['SII6731_sfr_att'],axis=0)
    OII3727 = np.sum(data['OII3727_sfr_att'],axis=0)
else:
    Hbeta = np.sum(data['Hbeta_sfr'],axis=0)
    OIII5007 = np.sum(data['OIII5007_sfr'],axis=0)
    NII6548 = np.sum(data['NII6584_sfr'],axis=0)
    Halpha = np.sum(data['Halpha_sfr'],axis=0)
    SII6731 = np.sum(data['SII6717_sfr'],axis=0) + np.sum(data['SII6731_sfr'],axis=0)
    OII3727 = np.sum(data['OII3727_sfr'],axis=0)

lz = data['lz_sfr'][:,0]

lu_sfr = data['lu_sfr'][:,0]
lne = data['lne_sfr'][:,0]

lssfr = data['lssfr'][:,0]

lms = np.log10(10**data['lms'][:,0] + 10**data['lms'][:,1])

lu_agn = data['lu_agn'][:,0]

Lagn = data['Lagn'][0]

Halpha_sfr = np.copy(Halpha)

if att:
    Halpha_agn = data['Halpha_agn_att'][0,:]
else:
    Halpha_agn = data['Halpha_agn'][0,:]

Halpha2 = np.copy(Halpha)

if att:
    Hbeta += data['Hbeta_agn_att'][0,:]
    OIII5007 += data['OIII5007_agn_att'][0,:]
    NII6548 += data['NII6584_agn_att'][0,:]
    Halpha += data['Halpha_agn_att'][0,:]
    SII6731 += data['SII6717_agn_att'][0,:] + data['SII6731_agn_att'][0,:]
    OII3727 += data['OII3727_agn_att'][0,:]
else:
    Hbeta += data['Hbeta_agn'][0,:]
    OIII5007 += data['OIII5007_agn'][0,:]
    NII6548 += data['NII6584_agn'][0,:]
    Halpha += data['Halpha_agn'][0,:]
    SII6731 += data['SII6717_agn'][0,:] + data['SII6731_agn'][0,:]
    OII3727 += data['OII3727_agn'][0,:] 

r = data['m_R'][0]
k = data['m_K'][0]

if att:
    Ha_flux_sfr = np.sum(data['Halpha_sfr_flux_att'],axis=0)
    Ha_flux_agn = np.sum(data['Halpha_agn_flux_att'],axis=0)
    Ha_flux = Ha_flux_sfr + Ha_flux_agn
    
    Hb_flux_sfr = np.sum(data['Hbeta_sfr_flux_att'],axis=0)
    Hb_flux_agn = np.sum(data['Hbeta_agn_flux_att'],axis=0)
    Hb_flux = Hb_flux_sfr + Hb_flux_agn
    
    NII6548_flux_sfr = np.sum(data['NII6584_sfr_flux_att'],axis=0)
    NII6548_flux_agn = np.sum(data['NII6584_agn_flux_att'],axis=0)
    NII6548_flux = NII6548_flux_sfr + NII6548_flux_agn
    
    OII3727_flux_sfr = np.sum(data['OII3727_sfr_flux_att'],axis=0)
    OII3727_flux_agn = np.sum(data['OII3727_agn_flux_att'],axis=0)
    OII3727_flux = OII3727_flux_sfr + OII3727_flux_agn
    
    OIII5007_flux_sfr = np.sum(data['OIII5007_sfr_flux_att'],axis=0)
    OIII5007_flux_agn = np.sum(data['OIII5007_agn_flux_att'],axis=0)
    OIII5007_flux = OIII5007_flux_sfr + OIII5007_flux_agn
    
    SII6731_flux_sfr = np.sum(data['SII6731_sfr_flux_att'],axis=0)
    SII6731_flux_agn = np.sum(data['SII6731_agn_flux_att'],axis=0)
    SII6731_flux = SII6731_flux_sfr + SII6731_flux_agn
    
    SII6717_flux_sfr = np.sum(data['SII6717_sfr_flux_att'],axis=0)
    SII6717_flux_agn = np.sum(data['SII6717_agn_flux_att'],axis=0)
    SII6731_flux = SII6731_flux + SII6717_flux_sfr + SII6717_flux_agn
else:
    Ha_flux_sfr = np.sum(data['Halpha_sfr_flux'],axis=0)
    Ha_flux_agn = np.sum(data['Halpha_agn_flux'],axis=0)
    Ha_flux = Ha_flux_sfr + Ha_flux_agn
    
    Hb_flux_sfr = np.sum(data['Hbeta_sfr_flux'],axis=0)
    Hb_flux_agn = np.sum(data['Hbeta_agn_flux'],axis=0)
    Hb_flux = Hb_flux_sfr + Hb_flux_agn
    
    NII6548_flux_sfr = np.sum(data['NII6584_sfr_flux'],axis=0)
    NII6548_flux_agn = np.sum(data['NII6584_agn_flux'],axis=0)
    NII6548_flux = NII6548_flux_sfr + NII6548_flux_agn
    
    OII3727_flux_sfr = np.sum(data['OII3727_sfr_flux'],axis=0)
    OII3727_flux_agn = np.sum(data['OII3727_agn_flux'],axis=0)
    OII3727_flux = OII3727_flux_sfr + OII3727_flux_agn
    
    OIII5007_flux_sfr = np.sum(data['OIII5007_sfr_flux'],axis=0)
    OIII5007_flux_agn = np.sum(data['OIII5007_agn_flux'],axis=0)
    OIII5007_flux = OIII5007_flux_sfr + OIII5007_flux_agn
    
    SII6731_flux_sfr = np.sum(data['SII6731_sfr_flux'],axis=0)
    SII6731_flux_agn = np.sum(data['SII6731_agn_flux'],axis=0)
    SII6731_flux = SII6731_flux_sfr + SII6731_flux_agn
    
    SII6717_flux_sfr = np.sum(data['SII6717_sfr_flux'],axis=0)
    SII6717_flux_agn = np.sum(data['SII6717_agn_flux'],axis=0)
    SII6731_flux = SII6731_flux + SII6717_flux_sfr + SII6717_flux_agn

ind2 = np.where((Hbeta>0)&(OIII5007>0)&(NII6548>0)&(Halpha>0)&
    (SII6731>0)&(OII3727>0)&(lu_sfr>minU)&(lu_sfr<maxU)&
    (lz>np.log10(minZ))&(lz<np.log10(maxZ)))[0]

f.close()

Halpha_sfr = Halpha_sfr[ind2]
Halpha_agn = Halpha_agn[ind2]

Hbeta = Hbeta[ind2]
OIII5007 = OIII5007[ind2]
NII6548 = NII6548[ind2]
Halpha = Halpha[ind2]
SII6731 = SII6731[ind2]
OII3727 = OII3727[ind2]

Halpha_ratio = Halpha_agn/(Halpha2[ind2])

lz = lz[ind2]
lssfr = lssfr[ind2]
lms = lms[ind2]

lu_sfr = lu_sfr[ind2]
lne = lne[ind2]

lu_agn = lu_agn[ind2]
Lagn = Lagn[ind2]

Ha_flux = Ha_flux[ind2]
Hb_flux = Hb_flux[ind2]
NII6548_flux = NII6548_flux[ind2]
OII3727_flux = OII3727_flux[ind2]
OIII5007_flux = OIII5007_flux[ind2]
SII6731_flux = SII6731_flux[ind2]

r = r[ind2]
k = k[ind2]

bpt_x = ['log$_{10}$([NII]$\\lambda$6584/H$\\alpha$)',
         'log$_{10}$([SII]$\\lambda$6724/H$\\alpha$)']
my_x = np.array([np.log10(NII6548 / Halpha),np.log10(SII6731 / Halpha)])

bpt_y = ['log$_{10}$([OIII]$\\lambda$5007/H$\\beta$)',
         'log$_{10}$([OIII]$\\lambda$5007/H$\\beta$)']
my_y = np.array([np.log10(OIII5007 / Hbeta),np.log10(OIII5007 / Hbeta)])


my_x_obs = [np.log10(NII6548_favole / Ha_favole),np.log10(SII6731_favole / Ha_favole)]
my_y_obs = [np.log10(OIII5007_favole / Hb_favole),np.log10(OIII5007_favole / Hb_favole)]

lum_Halpha = np.log10(Halpha*const.h**2) 

z = [np.log10(Ha_flux), lum_Halpha, lms, lz, Halpha_ratio, lssfr]
vmin = [-16,                39,     8,  -3,     0,         -12]
vmax = [-13,                42,    12,  -1,     1,          -7]


if plot==0:
    xmin=-1.9
    xmax=0.8
    ymin=-1.5
    ymax=1.5
    
    x = np.arange(xmin, xmax+0.1, 0.03)
    
    SFR_Composite = lines_BPT(x,'NII','SFR_Composite')
    Composite_AGN = lines_BPT(x,'NII','Composite_AGN')
    LINER_NIIlim = lines_BPT(x,'NII','LINER_NIIlim')
    LINER_OIIIlim = lines_BPT(x,'NII','LINER_OIIIlim')
    
    plt.plot(x[x<0.05],SFR_Composite[x<0.05],'k--',markersize=3)
    plt.plot(x[x<0.47],Composite_AGN[x<0.47],'k.',markersize=3)
    plt.vlines(LINER_NIIlim,ymin,LINER_OIIIlim,'k',linestyles='dashdot')
    plt.hlines(LINER_OIIIlim,LINER_NIIlim,xmax,'k',linestyles='dashdot')
elif plot==1:
    xmin=-1.9
    xmax=0.9
    ymin=-2.1
    ymax=1.6
    
    x = np.arange(xmin, xmax+0.1, 0.03)
    
    SFR_AGN = lines_BPT(x,'SII','SFR_AGN')
    Seyfert_LINER = lines_BPT(x,'SII','Seyfert_LINER')
    
    plt.plot(x[x<0.32], SFR_AGN[x<0.32], 'k.', markersize=3)
    
    plt.plot(x[(Seyfert_LINER>SFR_AGN)|(x>=0.32)], Seyfert_LINER[(Seyfert_LINER>SFR_AGN)|(x>=0.32)], 'k.', markersize=3)   
    
if redshift==0.1:
    flux = 2e-16
    ind3 = (np.where((Ha_flux>flux)&(Hb_flux>flux)&(OIII5007_flux>flux)&
                      (NII6548_flux>flux)&(SII6731_flux>flux)&(r<17.77)&(Halpha_ratio>ratio_min))[0])
elif redshift==0.8:
    flux = 1e-16
    ind3 = (np.where((Ha_flux>flux)&(r<24.1)&(Halpha_ratio>ratio_min))[0])
elif redshift==1.5:
    flux = 5e-17
    ind3 = (np.where((Ha_flux>5e-17)&(k<23.5)&(Halpha_ratio>ratio_min))[0])
else:
    ind3 = np.where(r == r)[0]

plt.scatter(my_x_obs[plot], my_y_obs[plot], s=20, c='darkgrey', alpha=0.8)

plt.scatter(my_x[plot][ind3], my_y[plot][ind3], c=z[color][ind3], s=50, marker='o',cmap='jet',vmin=vmin[color],vmax=vmax[color],zorder=1)
cbar = plt.colorbar()
cbar.set_label(color_labels[color], rotation=270, labelpad =60, size=30)
cbar.ax.tick_params(labelsize=30)

xs_att = my_x[plot][ind3]
ys_att = my_y[plot][ind3]


plt.xlabel(bpt_x[plot],size=30)
plt.ylabel(bpt_y[plot],size=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

plt.xlim((xmin,xmax))
plt.ylim((ymin,ymax))
plt.grid()

plt.savefig('plots/BPT.png')
