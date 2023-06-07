#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:53:20 2023

@author: expox7
"""

import get_nebular_emission.eml as eml
import get_nebular_emission.eml_const as const
import time
import glob

# infile = glob.glob(r"red_GP14_z0.115_125Mpc.txt")
# infile = glob.glob(r"red_SAGE_z0.115_250Mpc.txt")
# infile = glob.glob(r"red_SAG_z0.115_250Mpc.txt")
# infile = glob.glob(r"red_L16_z0.115_125Mpc.txt")
# infile = glob.glob(r"red_GALACTICUS_z0.117_100Mpc.txt")
# infile = glob.glob(r"red_Henriques2015I_z0.114_250Mpc.txt")
# infile = glob.glob(r"GP20/iz61/ivol0/GP20_z0.0_0.txt")
infile = glob.glob(r"*/*/*/GP20_z0.0_*.txt")
# infile = glob.glob(r"example_data/emlines_lc16_PMILL_iz245_ivol*.dat")

infile_z0 = [None] #glob.glob(r"*/*/*/GP20_z0.0_*.txt")

redshift = 0

# outfile = r"output_data/emlines_L16_radio.hdf5"
# outfile = r"output_data/emlines_GAL_radio.hdf5"
# outfile = r"output_data/emlines_GP14_radio.hdf5"
# outfile = r"output_data/emlines_H15_radio.hdf5"
# outfile = r"output_data/emlines_SAGE.hdf5"
# outfile = r"output_data/emlines_Carlton.hdf5"
# outfile = r"output_data/emlines_SAG_2.hdf5"

# outfile = r"output_data/emlines_GP20_manual_expfit_nw.hdf5"
# outfile = r"output_data/emlines_GP20_manual_expfit.hdf5"

# outfile = r"output_data/emlines_GP20_z0.0_Lagn.hdf5"
# outfile = r"output_data/emlines_GP20_Lagn_Ztremonti2.hdf5"
# outfile = r"output_data/emlines_GP20_Lagn_Ztremonti2.hdf5"
# outfile = r"output_data/emlines_GP20_manual_no_starburst.hdf5"

# outfile = r"output_data/emlines_GP20_Lagn.hdf5"
# outfile = r"output_data/emlines_GP20_Lagn_starburst.hdf5"
# outfile = r"output_data/emlines_GP20_Lagn_no_starburst.hdf5"
# outfile = r"output_data/emlines_GP20_manual_starburst.hdf5"
# outfile = r"output_data/emlines_GP20_manual_no_starburst_expfit.hdf5"

outfile = r"output_data/emlines_GP20_z0.0_Lagn.hdf5"
# outfile = r"output_data/emlines_GP20_z0.8_Lagn.hdf5"
# outfile = r"output_data/emlines_GP20_z1.5_Lagn.hdf5"

# outfile = r"output_data/emlines_GP20_z0.0_Lagn.hdf5"
# outfile = r"output_data/emlines_GP20_z0.8_Lagn_kashino.hdf5"
# outfile = r"output_data/emlines_GP20_z1.5_Lagn_kashino.hdf5"

# Calculate the nebular emission lines
# The columns where the stellar mass, SFR and Zcold=Mzcold/Mcold are
# to be found in the input file, should be specified within the argument m_sfr_z
# Several components are allowed (e.g. disk and bulge).

# Mass: Msun/h
# SFR: Msun/(Gyr*h)
# Zcold: MZcold/Mcold

model = 'GP20'
# GP14, H15, GP20, GALACTICUS, SAG, SAGE

if model=='CARLTON':
    inputformat = 'textfile'
    LC2sfr = True
    
    IMF_i = ['Kennicut','Kennicut']
    
    cols = [[11,0,2],[15,1,4]]
    mtot2mdisk=False
    
    att_params=None
    
    epsilon_params=None
    
    AGNinputs = 'Lagn' #complete: Mburst, rbulg, vbulg, Mhot, Mbh
    Lagn_params=None #complete: [1,12,14,9,8]
    
    cutcols = [19]
    mincuts = [None]
    maxcuts = [19.8 - 38.034]

if model=='GP20':
    inputformat = 'textfile'
    LC2sfr = False
    
    IMF_i = ['Kennicut','Kennicut'] 
    
    cols = [[0,2,4],[1,3,5]]
    mtot2mdisk=False
    
    att_params=[11,6,4]
    
    epsilon_params=[6,11,19,12]
    
    AGNinputs = 'Lagn' #complete: Mburst, rbulg, vbulg, Mhot, Mbh, Spin
    Lagn_params=[17,8,20]#,15,16,9] #complete: [1,12,14,9,8,20]
    
    cutcols = [7,18,15]
    mincuts = [50*9.35e8,None,None]
    maxcuts = [None,19.8,None]

if model=='SAG':
    inputformat = 'textfile'
    LC2sfr = False
    
    IMF_i = ['Kennicut','Kennicut']
    
    cols = [[0,2,5],[1,3,7]]
    mtot2mdisk=False
    
    att_params=[12,4,5]
    
    epsilon_params=[4,12]
    
    AGNinputs = 'radio_mode' #complete: Mburst, rbulg, vbulg, Mhot, Mbh
    Lagn_params=[10,9] #complete: [1,12,14,9,8]
    
    cutcols = [8]
    mincuts = [50*1.51e9]
    maxcuts = [None]
    
if model=='SAGE':
    inputformat = 'textfile'
    LC2sfr = False
    
    IMF_i = ['Kennicut','Kennicut']
    
    cols = [[0,2,5],[1,3,5]]
    mtot2mdisk=False
    
    att_params=[10,4]
    
    epsilon_params=[4,10]
    
    AGNinputs = 'radio_mode' #complete: Mburst, rbulg, vbulg, Mhot, Mbh
    Lagn_params=[8,7] #complete: [1,12,14,9,8]
    
    cutcols = [6]
    mincuts = [50*1.51e9]
    maxcuts = [None]

if model=='GALACTICUS':
    inputformat = 'textfile'
    LC2sfr = False
    
    IMF_i = ['Kennicut','Kennicut']
    
    cols = [[0,2,5],[1,3,7]]
    mtot2mdisk=False
    
    att_params=[12,4]
    
    epsilon_params=[4,12,6,13]
    
    AGNinputs = 'radio_mode' #complete: Mburst, rbulg, vbulg, Mhot, Mbh
    Lagn_params=[10,9] #complete: [1,12,14,9,8]
    
    cutcols = [8]
    mincuts = [50*1.51e9]
    maxcuts = [None]

if model=='GP14':
    inputformat = 'textfile'
    LC2sfr = False
    
    IMF_i = ['Kennicut','Kennicut']
    
    cols = [[0,2,4],[1,3,4]]
    mtot2mdisk=True
    
    cols_att = [[10,14,18,22,26,30,34,38],[8,12,16,20,24,28,32,36]] # Total, bulge
    cols_notatt = [[11,15,19,23,27,31,35,39],[9,13,17,21,25,29,33,37]]
    
    cols_photmod = [0,1,2,3,5,6,7,8]
    # 'OII3727','Hbeta','OIII4959','OIII5007',
    # 'Halpha','NII6583','SII6717','SII6731'
    
    att_params=[cols_att,cols_notatt,cols_photmod]
    
    epsilon_params=[6,60]
    
    AGNinputs = 'radio_mode'
    Lagn_params=[58,57] #complete: [1,61,63,58,57]
    
    cutcols = [56,65]
    mincuts = [50*9.35e8,None]
    maxcuts = [None,19.8]
    
if model=='H15':
    inputformat = 'textfile'
    LC2sfr = False
    
    IMF_i = ['Kennicut','Kennicut']
    
    cols = [[0,2,5],[1,3,5]]
    mtot2mdisk=True
    
    att_params=[10,4]
    
    epsilon_params=[4,10]
    
    AGNinputs = 'radio_mode' #complete: Mbulg, rbulg, vbulg, Mhot, Mbh
    Lagn_params=[8,7] #[1,11,13,8,7] #[14,15,7]
    
    cutcols = [6]
    mincuts = [50*9.35e8]
    maxcuts = [None]

IMF_f = ['Kennicut','Top-heavy'] # IMF to which transform each component
# Kennicut, Salpeter, Kroupa, Chabrier, Baldry&Glazebrook, Top-heavy
# Kashino asumes Kroupa.

gamma = 1.3 # Orsi 2014

attmod='cardelli89'

unemod_sfr='panuzzo03'
unemod_agn='panuzzo03'

photmod_sfr='gutkin16'
photmod_agn='feltre16'

# lms_H15, lssfr_H15, lz_H15, Lagn_H15, n_H15, n2_H15 = 
eml.eml(infile, outfile, m_sfr_z=cols, infile_z0=infile_z0, 
                  cutcols=cutcols, mincuts=mincuts, 
                  maxcuts=maxcuts, att_params=att_params, h0=const.h, 
                  IMF_i=IMF_i, IMF_f=IMF_f, inputformat=inputformat,
                  mtot2mdisk=mtot2mdisk, LC2sfr=LC2sfr,
                  gamma=gamma,
                  AGNinputs=AGNinputs, Lagn_params=Lagn_params,
                  epsilon_params=epsilon_params,
                  attmod=attmod, unemod_sfr=unemod_sfr, 
                  unemod_agn=unemod_agn, photmod_sfr=photmod_sfr,
                  photmod_agn=photmod_agn,
                  verbose=True, Plotting=False, Testing=False)

# n, lms, lssfr, lz = np.copy(n_GAL), np.copy(lms_GAL[:,0]), np.copy(lssfr_GAL[:,0]), np.copy(lz_GAL[:,0])
# plt.figure(figsize=(10,10))
# plt.hist(np.log10(n[(lssfr<-11)]),range=(-3,3),bins=500,color='r',alpha=0.5,label=r'$sSFR < 10^{-11} M_\odot$')
# plt.hist(np.log10(n[(lssfr>-11)&(lssfr<-10)]),range=(-3,3),bins=500,color='g',alpha=0.5,label=r'$sSFR > 10^{-11} M_\odot$')
# plt.hist(np.log10(n[(lssfr>-10)&(lssfr<-9)]),range=(-3,3),bins=500,color='k',alpha=0.5,label=r'$sSFR > 10^{-10} M_\odot$')
# plt.hist(np.log10(n[(lssfr>-9)]),range=(-3,3),bins=500,color='b',alpha=0.5,label=r'$sSFR > 10^{-9} M_\odot$')
# plt.xlabel(r'$n_{\rm gas}$ [part/cm$^3$]',size=20)
# plt.ylabel('Number of galaxies',size=20)
# plt.xticks(size=20)
# plt.yticks(size=20)
# plt.legend(fontsize=20)
# plt.title('GALACTICUS',size=20)
# plt.grid()