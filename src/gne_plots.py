"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
.. contributions:: Olivia Vidal <ovive.pro@gmail.com>
.. contributions:: Julen Expósito-Márquez <expox7@gmail.com>
"""
import src.gne_io as io
import numpy as np
#import os.path
import h5py
## matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
#import matplotlib.gridspec as gridspec
#import src.gne_style as style
#from src.gne_stats import perc_2arrays
import src.gne_const as const
#from src.gne_io import get_nheader, check_file
from src.gne_photio import get_limits #,get_lines_Gutkin
#from numpy import random
#from scipy.stats import gaussian_kde
#import sys
#from cosmology import logL2flux, set_cosmology
from src.gne_cosmology import set_cosmology,emission_line_luminosity
import src.gne_style
plt.style.use(src.gne_style.style1)

cmap = 'jet'

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

    boundary = np.zeros(len(x)); boundary.fill(-999.)
    
    if BPT=='NII':
        if line=='Kauffmann2003':
            x0 = 0.05
            boundary[x<x0] = 0.61/(x[x<x0] - x0) + 1.3
        elif line=='Kewley2001':
            x0 = 0.47
            boundary[x<x0] = 0.61/(x[x<x0] - x0) + 1.19
        elif line=='LINER_NIIlim':
            boundary = np.log10(0.6) # Kauffmann 2003
        elif line=='LINER_OIIIlim':
            boundary = np.log10(3) # Kauffmann 2003
    elif BPT=='SII':
        if line=='Kewley2001':
            x0 = 0.32
            boundary[x<x0] = 0.72/(x[x<x0] - x0) + 1.3
        elif line=='Kewley2006':
            boundary = 1.89*x + 0.76
    else:
        print('STOP (gne_plots.lines_BPT): ',
              'BPT plot not recognized.')
        return None
            
    return boundary



#def test_sfrf(inputdata, outplot, obsSFR=None, obsGSM=None, colsSFR=[0,1,2,3],
#              colsGSM=[0,1,2,3], labelObs=None, specific=False, h0=const.h, volume=const.vol_pm, verbose=False):
#
#    '''
#    
#    Given log10(Mstar) and log10(sSFR) get the plots to compare log10(SFR) vs log10(Mstar).
#    Get the GSMF and the SFRF plots. 
#    Given the observations, compare the plots with the observations too.
# 
#    Parameters
#    ----------
# 
#    obsSFR : string
#      - Name of the input file for the SFR data observed.
#      - In text files (*.dat, *txt, *.cat), columns separated by ' '.
#      - In csv files (*.csv), columns separated by ','.
#      - Expected histogram mode:
#       - A column with the low value of the bin,
#       - A column with the high value of the bin,
#       - A column with the frequency in the bin,
#       - A column with the error. 
# 
#    obsGSM : string
#      - Name of the input file for the GSM data observed.
#      - In text files (*.dat, *txt, *.cat), columns separated by ' '.
#      - In csv files (*.csv), columns separated by ','.
#      - Expected histogram mode:
#       - A column with the low value of the bin,
#       - A column with the high value of the bin,
#       - A column with the frequency in the bin,
#       - A column with the error.
# 
#    colsSFR : list
#      - Columns with the data required to do the observational histogram of the SFR.
#      - Expected: [ind_column1, ind_column2, ind_column3, ind_column4]
#       - column1 is the column with the low values of the bins, in Msun/yr,
#       - column2 with the high values of the bins, in Msun/yr,
#       - column3 with the frequency, in Mpc^-3 dex^-1
#       - column4 with the error, in Mpc^-3 dex^-1
#       
#    colsGSM : list
#      - Columns with the data required to do the observational histogram of the GSM.
#      - Expected: [ind_column1, ind_column2, ind_column3, ind_column4]
#       - column1 is the column with the low values of the bins, in h^-2Msun,
#       - column2 with the high values of the bins, in h^-2Msun,
#       - column3 with the frequency, in h^-3 Mpc^-3,
#       - column4 with the error, in h^-3 Mpc^-3.
# 
#    labelObs : list of strings
#      - For the legend, add the name to cite the observational data source.
#      - ['GSM observed', 'SFR observed']
# 
#    outplot : string
#      - Name of the output file.
#      - Image-type files (*.pdf, *.jpg, ...)
#      
#    specific : boolean
#      If True it makes the plots with the sSFR. Otherwise, it makes the plots with the SFR.
# 
#    h0 : float
#      If not None: value of h, H0=100h km/s/Mpc.
#      
#    volume : float
#      - Carlton model default value = 542.16^3 Mpc^3/h^3.
#      - table 1: https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4922B/abstract
#      - If not 542.16**3. : valume of the simulation volume in Mpc^3/h^3
#    verbose : boolean
#      If True print out messages
# 
#    Notes
#    -------
#    It makes plot(log10(SFR),log10(Mstar)), plot GSMF and plot SFRF,
#    all three in one grid and saves it in the outplot path.
#    '''
#
#
#
#    # Define a class that forces representation of float to look a certain way
#    # This remove trailing zero so '1.0' becomes '1'
#    class nf(float):
#        def __repr__(self):
#            str = '%.1f' % (self.__float__(),)
#            if str[-1] == '0':
#                return '%.0f' % self.__float__()
#            else:
#                return '%.1f' % self.__float__()
#    # -----------------------------------------------------
#
#
#    # Correct the units of the simulation volume to Mpc^3:
#    if h0:
#        volume=volume/(h0**3)
#
#    #Prepare the plot
#    lsty = ['-',(0,(2,3))] # Line form
#
#    nds = np.array([-2., -3., -4., -5.]) # Contours values
#    al = np.sort(nds)
#
#    cm = plt.get_cmap('tab10')  # Colour map to draw colours from
#    color = []
#    for ii in range(0, 10):
#        col = cm(ii)
#        color.append(col)  # col change for each iteration
#
#
#    # Initialize GSMF (Galaxy Cosmological Mass Function)
#    mmin = 8 #10.3 # mass resolution 2.12 * 10**9 h0 M_sun (Baugh 2019)
#    mmax = 15 
#    dm = 0.2
#    mbins = np.arange(mmin, mmax, dm)
#    mhist = mbins + dm * 0.5
#    gsmf = np.zeros((len(mhist)))
#
#    # Initialize SSFRF
#    smin = -4.4
#    smax = 3
#    ds = 0.2
#    sbins = np.arange(smin, smax, ds)
#    shist = sbins + ds * 0.5
#    ssfrf = np.zeros((len(shist)))
#
#    # Initialize SFR vs M function
#    lenm = len(mhist)
#    lens = len(shist)
#    smf = np.zeros((lens,lenm))
#
#    # Plots limits and style
#    fig = plt.figure(figsize=(8.5, 9.))
#    gs = gridspec.GridSpec(3, 3)
#    gs.update(wspace=0., hspace=0.)
#    ax = plt.subplot(gs[1:, :-1])
#
#    # Fig. sSFR vs M
#    xtit = "log$_{10}(\\rm M_{*}$ [M$_\odot$])"
#    if specific:
#        ytit = "log$_{10}(\\rm sSFR/Gyr^{-1})$"
#    else:
#        ytit = "log$_{10}(\\rm SFR$ [M$_\odot$ yr$^{-1}$])"
#    xmin = 8.5; xmax = 12.25; ymin = smin;  ymax = smax
#    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
#    ax.set_xlabel(xtit); ax.set_ylabel(ytit)
#
#    # GSMF
#    axm = plt.subplot(gs[0, :-1],sharex=ax)
#    ytit="log$_{10}(\Phi(M_*))$" ; axm.set_ylabel(ytit)
#    axm.set_autoscale_on(False) ;  axm.minorticks_on()
#    axm.set_ylim(-5.5,-1)
#    plt.setp(axm.get_xticklabels(), visible=False)
#
#    # SSFRF
#    axs = plt.subplot(gs[1:, 2], sharey=ax)
#    if specific:
#        xtit = "log$_{10}(\Phi(sSFR))$"; axs.set_xlabel(xtit)
#    else:
#        xtit = "log$_{10}(\Phi(SFR))$"; axs.set_xlabel(xtit)
#    axs.set_autoscale_on(False); axs.minorticks_on()
#    axs.set_xlim(-5.5, 0.0)
#    start, end = axs.get_xlim()
#    axs.xaxis.set_ticks(np.arange(-4., end, 1.))
#    plt.setp(axs.get_yticklabels(), visible=False)
#
#    # Data Observations
#
#    # SFR observed
#
#    if obsSFR:
#        ih = get_nheader(obsSFR)
#
#        dataSFR = [0]*len(colsSFR)
#
#        for ii, col in enumerate(colsSFR):
#            #print(ii,col,colsSFR[ii])
#            data = np.loadtxt(obsSFR,skiprows=ih, usecols=col, unpack=True)
#            dataSFR[ii] = np.array(data)
#
#        dex = dataSFR[1]-dataSFR[0]
#        histSFR = dataSFR[1]-0.5*dex
#        errorSFR = dataSFR[3]
#
#    # GSM observed
#    if obsGSM:
#        ih = get_nheader(obsGSM)
#
#        dataGSM = [0]*len(colsGSM)
#
#        for ii, col in enumerate(colsGSM):
#            data = np.loadtxt(obsGSM,skiprows=ih, usecols=col, unpack=True)
#            dataGSM[ii] = np.array(data)
#
#        dex = dataGSM[1] - dataGSM[0]
#
#        # Change the units from h^-2 Msun to Msun.
#        histGSM = dataGSM[1] - 2*np.log10(h0) - 0.5*dex
#
#        # Change the units from h^3 Mpc^-3 to Mpc^-3
#        freqGSM = np.log10((dataGSM[2])) + 3 * np.log10(h0)
#        
#        lowGSM = np.log10(dataGSM[2]-dataGSM[3]) + 3 * np.log10(h0)
#        
#        lowGSM = abs(lowGSM - freqGSM)
#
#    for ii in range(len(inputdata)):
#
#        with h5py.File(inputdata[ii],'r') as file:
#            data = file['data']          
#            lms = np.log10((10**data['lms'][:,0])/const.IMF_M['Chabrier']+10**data['lms'][:,1]*const.IMF_M['Top-heavy']/const.IMF_M['Chabrier']) #+ np.log10(h0)
#            if specific:
#                lsfr = np.log10(10**data['lssfr'][:,0]+10**data['lssfr'][:,1]) + 9
#            else: 
#                lsfr = np.log10(10**data['lssfr'][:,0]+10**data['lssfr'][:,1]) + lms
#                lsfr = lsfr/const.IMF_SFR['Chabrier']
#            # lms = lms + np.log10(h0)     
#            del data
#
#
#        # Make the histograms
#
#        H, bins_edges = np.histogram(lms, bins=np.append(mbins, mmax))
#        gsmf = H / volume / dm  # In Mpc^3/h^3
#
#        H, bins_edges = np.histogram(lsfr, bins=np.append(sbins, smax))
#        sfrf = H / volume / ds # / const.h**-3
#
#        H, xedges, yedges = np.histogram2d(lsfr, lms,
#                                           bins=([np.append(sbins, smax),
#                                                  np.append(mbins, mmax)]))
#        smf = H / volume / dm / ds
#
#
#        # Plot SMF vs SFR
#
#        matplotlib.rcParams['contour.negative_linestyle'] = lsty[ii]
#        zz = np.zeros(shape=(len(shist), len(mhist))); zz.fill(const.notnum)
#        ind = np.where(smf > 0.)
#        zz[ind] = np.log10(smf[ind])
#        
#        # print(zz[ind])
#
#        ind = np.where(zz > const.notnum)
#
#        if (np.shape(ind)[1] > 1):
#
#            # Contours
#            xx, yy = np.meshgrid(mbins, sbins)
#            # Here: How to find the levels of the data?
#            cs = ax.contour(xx, yy, zz, levels=al, colors=color[ii])
#            ax.clabel(cs, inline=1, fontsize=10)
#
#        # Plot GSMF
#        py = gsmf; ind = np.where(py > 0.)
#        x = mhist[ind]; y = np.log10(py[ind])
#        ind = np.where(y < 0.)
#        axm.plot(x[ind], y[ind], color=color[ii])
#
#        # Plot observations GSMF
#        if obsGSM and ii==0:
#            axm.errorbar(histGSM, freqGSM, yerr=lowGSM, marker='o', color=color[ii + 2],
#                             label=''+ labelObs[0] +'')
#                
#            leg2 = axm.legend(bbox_to_anchor=(0.025, -0.87, 1.5, 1.5), fontsize='small',
#                              handlelength=1.2, handletextpad=0.4)
#            leg2.get_texts()
#            leg2.draw_frame(False)
#        
#        # Plot SFRF
#        px = sfrf; ind = np.where(px > 0.)
#        y = shist[ind]; x = np.log10(px[ind])
#        ind = np.where(x < 0.)
#        axs.plot(x[ind], y[ind], color=color[ii], label='Model')
#            
#        # Plot observations SFRF
#        if obsSFR and ii==0:
#            axs.errorbar(dataSFR[2], histSFR, xerr=errorSFR, marker='o', color=color[ii + 3],
#                          label=''+ labelObs[1] +'')
#
#        leg = axs.legend(bbox_to_anchor=(-0.47, 0.1, 1.5, 1.38), fontsize='small',
#                          handlelength=1.2, handletextpad=0.4)
#        leg.get_texts()
#        leg.draw_frame(False)
#
#    plotf = outplot
#
#    # Save figures
#    print('Plot: {}'.format(plotf))
#    fig.savefig(plotf)
#
#def test_medians(infile, outplot, lines_cut=0, r_cut=999 ,i_cut=999 ,k_cut=999, verbose=False):
#    '''
#    Given U and ne calculated from the Mstar and the SFR in gne_une.
#    get the plot of the medians of these quantities in masses bins.
#    Compare the plots between the data calculated from the average SFR and
#    from the SFR calculated from the LC photons.
#
#    Parameters
#    ----------
#    outplot : string
#     Path to the folder where save the plot.
#    verbose : boolean
#     If True print out messages.
#    lines_cut : float
#     Minimum flux for emision lines (Ha, Hb, OIII, NII and SII; Units: erg s^-1 cm^-2).
#    r_cut : float
#     Maximum r magnitude.
#    i_cut : float
#     Maximum i magnitude.
#    k_cut : float
#     Maximum k magnitude.
#
#    Notes
#    -------
#    Medians U and ne for each bin of mass in two differents plots.
#
#    '''
#    
#    set_cosmology(omega0=const.omega0, omegab=const.omegab,lambda0=const.lambda0,h0=const.h)
#    U_ne = ['u', 'ne']
#    cm = plt.get_cmap('tab10')  # Colour map to draw colours from
#    color = []
#
#    # Prepare the bins
#    mmin = 9.25
#    mmax = 12.25
#    dm = 0.2
#    mbins = np.arange(mmin, (mmax + dm), dm)
#    mbinsH = np.arange(mmin, mmax, dm)
#    mhist = mbinsH + dm * 0.5
#
#    for iu, une in enumerate(U_ne):
#        
#        if iu==0:
#            mmin = 9.5
#            mmax = 11.5
#            dm = 0.2
#            mbins = np.arange(mmin, (mmax + dm), dm)
#            mbinsH = np.arange(mmin, mmax, dm)
#            mhist = mbinsH + dm * 0.5
#            
#            fig, ax = plt.subplots(2, 1, figsize=(10,20), sharex=True, sharey=True)
#
#            for axis in ax.flat:
#                axis.tick_params(labelsize=35)
#                axis.set_xlim((mmin,mmax))
#                axis.set_ylim((-4.75,-1.25))
#                # axis.locator_params(axis='both', nbins=6)
#                axis.locator_params(axis='x', nbins=6)
#                axis.locator_params(axis='y', nbins=12)
#                axis.grid()
#                
#            ax[1].set_xlabel(r'$\log M_*$ [M$_\odot$]',size=35)
#            ylabels = [r'$\log U_{\rm SF}$',r'$\log U_{\rm SF}$',r'$\log U_{\rm AGN}$']
#        else:
#            mmin = 9.25
#            mmax = 12.25
#            dm = 0.2
#            mbins = np.arange(mmin, (mmax + dm), dm)
#            mbinsH = np.arange(mmin, mmax, dm)
#            mhist = mbinsH + dm * 0.5
#            
#            plt.figure()
#            plt.xlabel(r'$\log M_*$ [M$_\odot$]',size=15)
#            plt.ylim((1.25,2.75))
#            plt.ylabel(r'$\log n_e$ [cm$^{-3}$]',size=15)
#
#        plotf = outplot + '/test_medians_'+ une+'.pdf'
#        col = cm(iu)
#        color.append(col)# col change for each iteration
#        
#        if iu==0:
#            plt.ylim((-4.75,-1.25))
#            ylabels = [r'$\log U_{\rm SF}$',r'$\log U_{\rm SF}$',r'$\log U_{\rm AGN}$']
#        else:
#            plt.ylim((1.25,2.75))
#            plt.ylabel(r'$\log n_e$ [cm$^{-3}$]',size=15)
#
#        with h5py.File(infile,'r') as file:
#            f = file['data']
#            
#            lu = f['lu_sfr'][:,0]
#            lne = f['lne_sfr'][:,0]
#            
#            lus = [f['lu_sfr'][:,0],f['lu_sfr'][:,1],f['lu_agn'][:,0]]
#            lnes = [f['lne_sfr'][:,0],f['lne_sfr'][:,1]]
#            lms = np.log10(10**f['lms'][:,0] + 10**f['lms'][:,1])
#            
#            Ha_flux_sfr = np.sum(f['Halpha_sfr_flux'],axis=0)
#            Ha_flux_agn = np.sum(f['Halpha_agn_flux'],axis=0)
#            Ha_flux = Ha_flux_sfr + Ha_flux_agn
#            
#            Hb_flux_sfr = np.sum(f['Hbeta_sfr_flux'],axis=0)
#            Hb_flux_agn = np.sum(f['Hbeta_agn_flux'],axis=0)
#            Hb_flux = Hb_flux_sfr + Hb_flux_agn
#            
#            NII6548_flux_sfr = np.sum(f['NII6584_sfr_flux'],axis=0)
#            NII6548_flux_agn = np.sum(f['NII6584_agn_flux'],axis=0)
#            NII6548_flux = NII6548_flux_sfr + NII6548_flux_agn
#            
#            OII3727_flux_sfr = np.sum(f['OII3727_sfr_flux'],axis=0)
#            OII3727_flux_agn = np.sum(f['OII3727_agn_flux'],axis=0)
#            OII3727_flux = OII3727_flux_sfr + OII3727_flux_agn
#            
#            OIII5007_flux_sfr = np.sum(f['OIII5007_sfr_flux'],axis=0)
#            OIII5007_flux_agn = np.sum(f['OIII5007_agn_flux'],axis=0)
#            OIII5007_flux = OIII5007_flux_sfr + OIII5007_flux_agn
#            
#            SII6731_flux_sfr = np.sum(f['SII6731_sfr_flux'],axis=0)
#            SII6731_flux_agn = np.sum(f['SII6731_agn_flux'],axis=0)
#            SII6731_flux = SII6731_flux_sfr + SII6731_flux_agn
#            
#            SII6717_flux_sfr = np.sum(f['SII6717_sfr_flux'],axis=0)
#            SII6717_flux_agn = np.sum(f['SII6717_agn_flux'],axis=0)
#            SII6731_flux = SII6731_flux + SII6717_flux_sfr + SII6717_flux_agn
#            
#            r = f['m_R'][0]
#            k = f['m_K'][0]
#            I = f['m_I'][0]
#            
#        for num in range(3):
#            
#            if iu==1 and num==2:
#                break
#            
#            if iu==0:
#                lu = lus[num]
#            else:
#                lu = lnes[num]
#            
#            cut = (np.where((lu!=const.notnum)&(lne!=const.notnum)&
#                             (Ha_flux>lines_cut)&(Hb_flux>lines_cut)&(OIII5007_flux>lines_cut)&
#                              (NII6548_flux>lines_cut)&(SII6731_flux>lines_cut)&
#                              (r<r_cut)&(I<i_cut)&(k<k_cut))[0])
#            
#            # MEDIANS:
#            median = perc_2arrays(mbins, lms[cut], lu[cut], 0.5)
#        
#            ind = np.where(median>const.notnum)
#            
#            median = median[ind]
#        
#            # QUARTILES:
#            up_qu = perc_2arrays(mbins, lms[cut], lu[cut], 0.75)[ind]
#            low_qu = perc_2arrays(mbins, lms[cut], lu[cut], 0.25)[ind]    
#            qu = np.append([median-low_qu],[up_qu-median],axis=0)
#            
#            if iu==0:
#                if num==0:
#                    eb1=ax[0].errorbar(mhist[ind],median,marker='o',yerr=qu,elinewidth=0.5, capsize=5, 
#                                 color='black')
#                    eb1[-1][0].set_linestyle('-')
#                elif num==1:
#                    eb2=ax[0].errorbar(mhist[ind],median,marker='o',ls='--',yerr=qu,elinewidth=0.5, capsize=5, 
#                                 color='black')
#                    eb2[-1][0].set_linestyle('--')
#                else:
#                    ax[1].errorbar(mhist[ind],median,marker='o',yerr=qu,elinewidth=0.5, capsize=5, 
#                             color='black')
#
#                if num==0:
#                    ax[0].set_ylabel(ylabels[num],size=35)
#                elif num==2:
#                    ax[1].set_ylabel(ylabels[num],size=35)
#            else:
#                if num==0:
#                    eb1=plt.errorbar(mhist[ind],median,marker='o',ls='-',yerr=qu,elinewidth=0.5, capsize=5, 
#                             color='black')
#                    eb1[-1][0].set_linestyle('-')
#                else:
#                    eb2=plt.errorbar(mhist[ind],median,marker='o',ls='--',yerr=qu,elinewidth=0.5, capsize=5, 
#                             color='black')
#                    eb2[-1][0].set_linestyle('--')
#        
#        if iu==1:
#            plt.xlim((mmin,mmax))
#            plt.grid()
#        else:
#            fig.subplots_adjust(wspace=0,hspace=0)
#        plt.savefig(plotf)
#        # plt.close()
#
#

#def test_interpolation(infile, zz, verbose=True):
#    '''
#    Run a test of the interpolations done in gne_photio.
#    Two plots, one to verify the U interpolation and the other one to verify the Z interpolation
#    
#    Parameters
#    ----------
#    infile : string
#     Name of the input file. 
#    outplot : string
#     Path to the folder plot.
#    photmod : string
#      Photoionisation model to be used for look up tables.
#    plot_phot : boolean
#     If True it plots points from the photoionization tables.
#    create_file : boolean
#     If True it creates textfiles to read the photoionization tables.
#    file_folder : string
#     Folder where the textfiles to read the tables will be/are stored.
#    verbose : boolean
#     If True print out messages.
#
#    Notes
#    -------
#    Plot of several BPT diagrams.
#    '''
#    
#    set_cosmology(omega0=const.omega0, omegab=const.omegab,lambda0=const.lambda0,h0=const.h)
#    
#    for num in range(len(infile)):
#    
#        check_file(infile[num], verbose=True)
#        f = h5py.File(infile[num], 'r')
#        data = f['data']
#    
#        lu_disk = data['lu'][:,0]
#        lne_disk = data['lne'][:,0]
#        lzgas_disk = data['lz'][:,0]
#        
#        minU, maxU = get_limits(propname='U', photmod=photmod)
#        minnH, maxnH = get_limits(propname='nH', photmod=photmod)
#        minZ, maxZ = get_limits(propname='Z', photmod=photmod)
#        
#        ignore = True
#        if ignore:
#            ind = np.where((lu_disk!=minU)&(lu_disk!=maxU)&(lzgas_disk!=np.log10(minZ))&(lzgas_disk!=np.log10(maxZ))&
#                       (lne_disk!=np.log10(minnH))&(lne_disk!=np.log10(maxnH)))[0]
#        else:
#            ind = np.arange(len(lu_disk))
#
#        Hbeta = np.sum(data['Hbeta'],axis=0)[ind]
#        OIII5007 = np.sum(data['OIII5007'],axis=0)[ind]
#        NII6548 = np.sum(data['NII6583'],axis=0)[ind]
#        Halpha = np.sum(data['Halpha'],axis=0)[ind]
#        SII6717_6731 = np.sum(data['SII6731'],axis=0)[ind]
#        OII3727 = np.sum(data['OII3727'],axis=0)[ind]
#        
#        lz = data['lz'][:,0]
#        lz = lz[ind]
#        
#        lssfr = data['lssfr'][:,0]
#        lssfr = lssfr[ind]
#        
#        lms = np.log10(10**data['lms'][:,0] + 10**data['lms'][:,1])
#        lms = lms[ind]
#        
#        ind2 = np.where((Hbeta>0)&(OIII5007>0)&(NII6548>0)&(Halpha>0)&(SII6717_6731>0)&(OII3727>0))[0]
#        
#        print(len(ind),len(ind2))
#        
#        Hbeta = Hbeta[ind2]
#        OIII5007 = OIII5007[ind2]
#        NII6548 = NII6548[ind2]
#        Halpha = Halpha[ind2]
#        SII6717_6731 = SII6717_6731[ind2]
#        OII3727 = OII3727[ind2]
#        
#        lz = lz[ind2]
#        lssfr = lssfr[ind2]
#        lms = lms[ind2]
#        
#        bpt_x = ['log$_{10}$([NII]$\\lambda$6584/H$\\alpha$)',
#                 'log$_{10}$([SII]$\\lambda$6731/H$\\alpha$)',
#                 'log$_{10}$([NII]$\\lambda$6584/[OII]$\\lambda$3727)',
#                 'log$_{10}$([NII]$\\lambda$6584/H$\\alpha$)']
#        my_x = [np.log10(NII6548 / Halpha),np.log10(SII6717_6731 / Halpha),np.log10(NII6548 / OII3727)]#,np.log10(NII6548 / Halpha)]
#        
#        bpt_y = ['log$_{10}$([OIII]$\\lambda$5007/H$\\beta$)',
#                 'log$_{10}$([OIII]$\\lambda$5007/H$\\beta$)',
#                 'log$_{10}$([OIII]$\\lambda$5007/[OII]$\\lambda$3727)',
#                 'log$_{10}$(EW(H$\\alpha$)/$\dot{A}$)']
#        my_y = [np.log10(OIII5007 / Hbeta),np.log10(OIII5007 / Hbeta),np.log10(OIII5007 / OII3727)]#,np.log10(EW_Halpha)]
#        
#        if not plot_phot:
#            for i in range(4):
#                plt.figure(figsize=(15,15))
#                
#                # X1, Y1 = np.mgrid[xmin:xmax:68j, ymin:ymax:68j]
#                # positions = np.vstack([X1.ravel(), Y1.ravel()])
#                # values = np.vstack([my_x[i], my_y[i]])
#                # kernel = stats.gaussian_kde(values,0.75)
#                # BPT = np.reshape(kernel(positions).T, X1.shape)
#                # plt.imshow(BPT, cmap=plt.cm.gist_earth_r,extent=[xmin, xmax, ymin, ymax],aspect=(xmax-xmin)/(ymax-ymin))#,vmin=0,vmax=1)
#                
#                if i==0:
#                    xmin=-2.2
#                    xmax=1
#                    ymin=-2
#                    ymax=2
#                    
#                    x = np.arange(xmin, xmax+0.1, 0.03)
#                    
#                    SFR_Composite = lines_BPT(x,'NII','SFR_Composite')
#                    Composite_AGN = lines_BPT(x,'NII','Composite_AGN')
#                    LINER_NIIlim = lines_BPT(x,'NII','LINER_NIIlim')
#                    LINER_OIIIlim = lines_BPT(x,'NII','LINER_OIIIlim')
#                    
#                    plt.plot(x[x<0.05],SFR_Composite[x<0.05],'k--',markersize=3)
#                    plt.plot(x[x<0.47],Composite_AGN[x<0.47],'k.',markersize=3)
#                    plt.vlines(LINER_NIIlim,ymin,LINER_OIIIlim,'k',linestyles='dashdot')
#                    plt.hlines(LINER_OIIIlim,LINER_NIIlim,xmax,'k',linestyles='dashdot')
#                elif i==1:
#                    xmin=-2.6 #-1.6
#                    xmax=0.2
#                    ymin=-1.9
#                    ymax=1.5
#                    
#                    x = np.arange(xmin, xmax+0.1, 0.03)
#                    
#                    SFR_AGN = lines_BPT(x,'SII','SFR_AGN')
#                    Seyfert_LINER = lines_BPT(x,'SII','Seyfert_LINER')
#                    
#                    plt.plot(x[x<0.32], SFR_AGN[x<0.32], 'k.', markersize=3)
#                    
#                    plt.plot(x[(Seyfert_LINER>SFR_AGN)|(x>=0.32)], Seyfert_LINER[(Seyfert_LINER>SFR_AGN)|(x>=0.32)], 'k.', markersize=3)
#                elif i==2:
#                    xmin=-1.9
#                    xmax=0.9
#                    ymin=-2.1
#                    ymax=1.6
#                elif i==3:
#                    xmin=-2.2
#                    xmax=1.2
#                    ymin=-1
#                    ymax=3
#                    
#                    # x = np.arange(xmin, xmax+0.1, 0.03)
#            
#                # xy = np.vstack([my_x[i], my_y[i]])
#                # z = gaussian_kde(xy)(xy)
#                # z = z/np.amax(z)
#                # np.save('density_galform_o_g1.3_ratios_' + str(i),z)
#                
#                # z = np.load('density_galform_kashino_ratios_' + str(i) + '.npy')
#                # z = np.log10(z)
#                
#                # Ha_flux = np.zeros(Halpha.shape)
#                # for j in range(len(Halpha)):
#                #     Ha_flux[j] = logL2flux(Halpha[j],0.131)
#                    
#                # ind = np.where((Ha_flux>2e-15))
#                
#                z = Halpha
#                
#                vmin = 40.5
#                vmax = 43
#
#                plt.scatter(my_x[i][ind], my_y[i][ind], c=z[ind], s=1, marker='o',cmap='jet',vmin=vmin, vmax=vmax)
#                cbar = plt.colorbar()
#                cbar.set_label(r'$\log H_\alpha \ [\rm erg/s]$', rotation=270, labelpad =40, size=30)
#                cbar.ax.tick_params(labelsize=30)
#                
#                #'$\log \bar{n}_p$'
#                #'$\log M_* \ [M_\odot]$'
#                #'$\log Z$'
#                #'$\log SFR \ [M_\odot/yr]$'
#                #'$\log H_\alpha \ [\rm erg/s]$'
#                
#                plt.xlabel(bpt_x[i],size=30)
#                plt.ylabel(bpt_y[i],size=30)
#                plt.xticks(fontsize=30)
#                plt.yticks(fontsize=30)
#                
#                plt.xlim((xmin,xmax))
#                plt.ylim((ymin,ymax))
#                plt.grid()
#                
#                plotnom = outplot + '/BPTplot_' + str(i) + '_' + str(num) + '_k.png'
#                
#                # np.save(outplot + '/BPTplot_' + str(i) + '_' + str(num), np.array([my_x,my_y]))
#            
#                plt.savefig(plotnom)
#                # plt.close()
#                
#                print(str(i+1) + ' de 4.')     
#        
#        if plot_phot:       
#            if photmod not in const.photmods:
#                if verbose:
#                    print('STOP (gne_photio.test_bpt): Unrecognised model to get emission lines.')
#                    print('                Possible photmod= {}'.format(const.photmods))
#                sys.exit()
#            elif (photmod == 'gutkin16'):
#                
#                Z = ['0001', '0002', '0005', '001', '002', '004', '006', '008', '010', '014', '017', '020', '030', '040']
#            
#                zz = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.014, 0.017, 0.02, 0.03, 0.04]
#            
#                uu = [-1., -1.5, -2., -2.5, -3., -3.5, -4.]
#            
#                ne = ['100']  # ['10', '100', '1000','10000']
#            
#                cm = plt.get_cmap('tab20') # Colour map to draw colours from
#            
#                if create_file:
#                    for iz, zname in enumerate(Z):
#                        infile = r"nebular_data/gutkin16_tables/nebular_emission_Z" + zname + ".txt"
#                
#                        ih = get_nheader(infile)
#                
#                        datane = np.loadtxt(infile, skiprows=ih, usecols=(2), unpack=True)
#                        datalu = np.loadtxt(infile, skiprows=ih, usecols=(0), unpack=True)
#                
#                        OIII5007_model = np.loadtxt(infile, skiprows=ih, usecols=(8), unpack=True)
#                        Hb_model = np.loadtxt(infile, skiprows=ih, usecols=(6), unpack=True)
#                        NII6548_model = np.loadtxt(infile, skiprows=ih, usecols=(9), unpack=True)
#                        Ha_model = np.loadtxt(infile, skiprows=ih, usecols=(10), unpack=True)
#                        SII6717_6731_model = np.loadtxt(infile, skiprows=ih, usecols=(12), unpack=True) + np.loadtxt(infile, skiprows=ih, usecols=(12), unpack=True)
#                
#                        for ii, nh in enumerate(ne):
#                            outfile = r"output_data/Gutkinfile_n_" + nh + ".txt"
#                            if iz==0 and os.path.exists(outfile):
#                                os.remove(outfile)
#                
#                            header1 = 'Z, U, NII6584/Ha  OIII5007/Hb, SII(6717+6731)/Ha'
#                
#                            ind = np.where(datane == float(nh))
#                            x = np.log10(NII6548_model[ind] / Ha_model[ind])
#                            y = np.log10(OIII5007_model[ind] / Hb_model[ind])
#                            p = np.log10(SII6717_6731_model[ind] / Ha_model[ind])
#                            u = datalu[ind]
#                            z = np.full(np.shape(u), zz[iz])
#                
#                            tofile = np.column_stack((z, u, x, y, p))
#                
#                            with open(outfile, 'a') as outf:
#                                if iz == 0:
#                                    np.savetxt(outf, tofile, delimiter=' ', header=header1)
#                                else:
#                                    np.savetxt(outf, tofile, delimiter=' ')
#                                outf.closed
#                else:
#                    for ii, nh in enumerate(ne):
#                        outfile = r"output_data/Gutkinfile_n_" + nh + ".txt"
#                        if not os.path.exists(outfile):
#                            print('STOP (gne_photio.test_bpt): Textfiles for table reading dont exist.')
#                            print('Create them with create_file = True.')
#            
#                cols = []
#                for iz, lz in enumerate(zz):
#                    col = cm(iz)
#                    cols.append(col)
#            
#                for ii, nh in enumerate(ne):
#                    infile = r"output_data/Gutkinfile_n_" + nh + ".txt"
#            
#                    ih = get_nheader(infile)
#            
#                    z = np.loadtxt(infile, skiprows=ih, usecols=(0), unpack=True)
#                    u = np.loadtxt(infile, skiprows=ih, usecols=(1), unpack=True)
#                    x = np.loadtxt(infile, skiprows=ih, usecols=(2), unpack=True)
#                    y = np.loadtxt(infile, skiprows=ih, usecols=(3), unpack=True)
#                    p = np.loadtxt(infile, skiprows=ih, usecols=(4), unpack=True)
#                    
#                    comp_x = [x,p]
#                    comp_y = [y,y]
#            
#                    # DIFERENTS COLORS FOR U:
#                        
#                    for i in range(2):
#                        plt.figure(figsize=(15,15))
#            
#                        for iu, lu in enumerate(uu):
#                            ind2 = np.where(u == uu[iu])
#                            plt.plot(comp_x[i][ind2], comp_y[i][ind2], marker='.', linewidth=0, color=cols[iu], label='U = ' + str(lu) + '')
#                
#                        labelsU = []
#                        for elem in lu_disk: labelsU.append('U = {}'.format(np.round(elem,2)))
#                        
#                        plt.plot(my_x[i], my_y[i], marker='o', markersize=2, linewidth=0, color='black')
#                        
#                        plt.xlabel(bpt_x[i],size=30)
#                        plt.ylabel(bpt_y[i],size=30)
#                        plt.xticks(fontsize=30)
#                        plt.yticks(fontsize=30)
#                        plt.grid()
#                        plt.legend()
#                        
#                        plotnom = outplot + '/BPTplot_U_' + str(i) + '_' + str(num) + '.png'
#                        
#                        print('U', str(i))
#                    
#                        plt.savefig(plotnom)
#                        plt.close()
#            
#                    # DIFFERENTS COLORS FOR Z:
#                        
#                    for i in range(2):
#                        plt.figure(figsize=(15,15))
#            
#                        for iz, lz in enumerate(zz):
#                            ind2 = np.where(z == zz[iz])
#                            plt.plot(comp_x[i][ind2], comp_y[i][ind2], marker='.', linewidth=0, color=cols[iz], label='Z = ' + str(lz) + '')
#                
#                        labelsZ = []
#                        for elem in lzgas_disk: labelsZ.append('Z = {:.4f}'.format(10 ** (elem)))
#                        
#                        plt.plot(my_x[i], my_y[i], marker='o', markersize=2, linewidth=0, color='black')
#                        
#                        plt.xlabel(bpt_x[i],size=30)
#                        plt.ylabel(bpt_y[i],size=30)
#                        plt.xticks(fontsize=30)
#                        plt.yticks(fontsize=30)
#                        plt.grid()
#                        plt.legend()
#                        
#                        plotnom = outplot + '/BPTplot_Z_' + str(i) + '_' + str(num) + '.png'
#                        
#                        print('Z', str(i))
#                    
#                        plt.savefig(plotnom)
#                        plt.close()
#


def test_bpts(infile, zz, snap, verbose=True):
    '''
    Make the 2 BPT diagrams without attenuation
    
    Parameters
    ----------
    infile : string
       Name of the input file. 
    zz : float
       Redshift
    snap: integer
        Simulation snapshot number    
    verbose : boolean
       If True print out messages.
    '''

    bptnoms = ['','']
    xmins = [-1.9,-1.9]
    xmaxs = [0.8,0.9]
    ymins = [-1.5,-2.1]
    ymaxs = [1.5,1.6]

    # Read header and SF information from file
    filenom = io.get_outnom(infile,snap,ftype='line_data')

    f = h5py.File(filenom, 'r') 
    header = f['header']
    redshift = header.attrs['redshift']
    omega0 = header.attrs['omega0']
    omegab = header.attrs['omegab']
    lambda0 = header.attrs['lambda0']
    h0 = header.attrs['h0']
    photmod_sfr = header.attrs['photmod_sfr']
    photmod_agn = header.attrs['photmod_agn']

    lu_sfr = f['sfr_data/lu_sfr'][:,0]
    lz_sfr = f['sfr_data/lz_sfr'][:,0]
    Ha_sfr = np.sum(f['sfr_data/Halpha_sfr'],axis=0)
    Hb_sfr = np.sum(f['sfr_data/Hbeta_sfr'],axis=0)
    NII6548_sfr = np.sum(f['sfr_data/NII6584_sfr'],axis=0)
    OII3727_sfr = np.sum(f['sfr_data/OII3727_sfr'],axis=0)
    OIII5007_sfr = np.sum(f['sfr_data/OIII5007_sfr'],axis=0)
    SII6731_sfr = np.sum(f['sfr_data/SII6731_sfr'],axis=0)
    SII6717_sfr = np.sum(f['sfr_data/SII6717_sfr'],axis=0)

    # Read AGN information if it exists
    AGN = True
    if 'agn_data' not in f.keys(): AGN = False
    
    if AGN:
        lu_agn = f['agn_data/lu_agn'][:,0]
        lz_agn = f['agn_data/lz_agn'][:,0]
        Ha_agn = np.sum(f['agn_data/Halpha_agn'],axis=0)
        Hb_agn = np.sum(f['agn_data/Hbeta_agn'],axis=0)
        NII6548_agn = np.sum(f['agn_data/NII6584_agn'],axis=0)
        OII3727_agn = np.sum(f['agn_data/OII3727_agn'],axis=0)
        OIII5007_agn = np.sum(f['agn_data/OIII5007_agn'],axis=0)
        SII6731_agn = np.sum(f['agn_data/SII6731_agn'],axis=0)
        SII6717_agn = np.sum(f['agn_data/SII6717_agn'],axis=0)

    # Magnitudes for cuts
    ismagr = True
    try:
        magr = f['data/magR'][:,0]
    except:
        ismagr = False

    ismagk = True
    try:
        magk = f['data/magK'][:,0]
    except:
        ismagk = False

    f.close()

    # Set the cosmology from the simulation
    set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0,h0=h0)
    
    # Line information
    if AGN:
        Ha = Ha_sfr + Ha_agn
        Hb = Hb_sfr + Hb_agn
        NII = NII6548_sfr + NII6548_agn
        OII = OII3727_sfr + OII3727_agn
        OIII = OIII5007_sfr + OIII5007_agn
        SII = SII6731_sfr + SII6731_agn +\
            SII6717_sfr + SII6717_agn
    else:
        Ha = Ha_sfr
        Hb = Hb_sfr
        NII = NII6548_sfr
        OII = OII3727_sfr
        OIII = OIII5007_sfr
        SII = SII6731_sfr + SII6717_sfr

    minU, maxU = get_limits(propname='U', photmod=photmod_sfr)
    minZ, maxZ = get_limits(propname='Z', photmod=photmod_sfr)

    minU, maxU = float(minU), float(maxU)
    minZ, maxZ = float(minZ), float(maxZ)

    ind = np.where((Ha>0)   & (Hb>0)  & 
                    (NII>0)  & (OII>0) &
                    (OIII>0) & (SII>0) &
                    (lu_sfr>minU)&(lu_sfr<maxU)&
                    (lz_sfr>np.log10(minZ))&(lz_sfr<np.log10(maxZ)))
    if (np.shape(ind)[1] < 1):
        print('STOP BPT plots: not enough adequate data')
        return None, None

    # For colourbar
    if AGN:
        Halpha_ratio = Ha_agn[ind]/Ha[ind]
    else:
        Halpha_ratio = Ha[ind]
    
    Ha = Ha[ind]
    Hb = Hb[ind]
    NII = NII[ind]
    OII = OII[ind]
    OIII = OIII[ind]
    SII = SII[ind]

    O3Hb =np.log10(OIII) - np.log10(Hb)
    N2Ha =np.log10(NII) - np.log10(Ha)
    S2Ha =np.log10(SII) - np.log10(Ha)

    if ismagr:
        mag_r = magr[ind]
    if ismagk:
        mag_k = magk[ind]

    sel = np.copy(ind)
    for ii, bpt in enumerate(['NII','SII']):
        # Set figure
        plt.figure(figsize=(15,15))
        ax = plt.subplot()
        ytit = 'log$_{10}$([OIII]$\\lambda$5007/H$\\beta$)'
        ax.set_xlim(xmins[ii], xmaxs[ii])
        ax.set_ylim(ymins[ii], ymaxs[ii])
        
        # Add obs. data and further cuts if adequate
        if redshift <= 0.2:
            obsplot = True
            obsfile = 'data/observational_data/favole2024.txt'
            data = np.loadtxt(obsfile,skiprows=1,usecols=(15,9))
            yobs = np.log10(data[:,0]/data[:,1]) #O3/Hb
            if bpt=='NII':
                data = np.loadtxt(obsfile,skiprows=1,usecols=(18,6))
                xobs = np.log10(data[:,0]/data[:,1]) #N2/Ha
            elif bpt=='SII':
                data = np.loadtxt(obsfile,skiprows=1,usecols=(21,6))
                xobs = np.log10(data[:,0]/data[:,1]) #S2/Ha

            flux = 2e-16 # erg/s/cm^2 Favole+2024
            Lmin = emission_line_luminosity(flux,redshift)*1e40 #erg/s

            if ismagr:
                sel = np.where((Ha> Lmin) & (Hb> Lmin) &
                               (OIII> Lmin) & (NII> Lmin) &
                               (SII> Lmin)&(mag_r<17.77))
            else:
                sel = np.where((Ha> Lmin) & (Hb> Lmin) &
                               (OIII> Lmin) & (NII> Lmin) &
                               (SII> Lmin))

        elif 0.7 <= redshift <= 0.9:            
            flux = 1e-16  # erg/s/cm^2 Kashino+2019
            Lmin = emission_line_luminosity(flux,redshift)*1e40 #erg/s
            
            if ismagr:
                sel = np.where((Ha> Lmin) & (mag_r<124.1))
            else:
                sel = np.where(Ha> Lmin)
                
        elif 1.45 <= redshift <= 1.75:
            obsplot = True
            if bpt=='NII':
                obsfile = 'data/observational_data/NII_Kashino.txt'
                yobs = np.loadtxt(obsfile,skiprows=18,usecols=(6)) #O3/Hb
                xobs = np.loadtxt(obsfile,skiprows=18,usecols=(3)) #N2/Ha
            elif bpt=='SII':
                obsfile = 'data/observational_data/SII_Kashino.txt'
                yobs = np.loadtxt(obsfile,skiprows=18,usecols=(6)) #O3/Hb
                xobs = np.loadtxt(obsfile,skiprows=18,usecols=(3)) #N2/Ha

            flux = 5e-17  # erg/s/cm^2 Kashino+2019
            Lmin = emission_line_luminosity(flux,redshift)*1e40 #erg/s
            
            if ismagk:
                sel = np.where((Ha> Lmin) & (mag_k<23.5))
            else:
                sel = np.where(Ha > Lmin)
                
        if obsplot:
            ax.scatter(xobs,yobs, s=20, c='darkgrey', alpha=0.8)
            if (np.shape(sel)[1]<1):
                sel1 = np.arange(0,len(Ha),1)
                sel = np.expand_dims(sel1, axis=0)
                print('WARNING: due to low numbers, using minimal cuts')

        # Model lines
        yobs = O3Hb[sel] #O3/Hb
        cha = Halpha_ratio[sel]
        
        if bpt=='NII':
            xobs = N2Ha[sel] #N2/Ha
        elif bpt=='SII':
            xobs = S2Ha[sel] #S2/Ha
        ax.scatter(xobs,yobs, c=cha,s=50, marker='o', cmap=cmap)
        
        # Add colorbar
        sm = ScalarMappable(cmap=cmap) # Create ScalarMappable
        sm.set_array(cha)
        cbar = plt.colorbar(sm, ax=ax, cmap=cmap)
        if AGN:
            collabel = r'$F_{\rm H_{\alpha}, AGN}/F_{\rm H_{\alpha}, tot}$'
        else:
            collabel = r'$F_{\rm H_{\alpha}}$'
        cbar.set_label(collabel,rotation=270,labelpad=60)        
            
        # Lines
        xline = np.arange(xmins[ii],xmaxs[ii]+0.1, 0.03)
        if bpt=='NII':
            xtit = 'log$_{10}$([NII]$\\lambda$6584/H$\\alpha$)'

            yline = lines_BPT(xline,bpt,'Kauffmann2003')
            ax.plot(xline,yline,'k.')

            yline = lines_BPT(xline,bpt,'Kewley2001')
            ax.plot(xline,yline,'k-')
        elif bpt=='SII':
            xtit = 'log$_{10}$([SII]$\\lambda$6584/H$\\alpha$)'

            yline = lines_BPT(xline,bpt,'Kewley2001')
            ax.plot(xline,yline,'k-')

            ylinel = lines_BPT(xline,bpt,'Kewley2006')
            ax.plot(xline[ylinel>yline],ylinel[ylinel>yline],'k--')

        ax.set_xlabel(xtit); ax.set_ylabel(ytit)
        
        # Output files
        bptnoms[ii] = io.get_outnom(infile,snap,ftype='plots',
                                    ptype=bpt+'bpt',verbose=verbose)
        plt.savefig(bptnoms[ii])

        if verbose:
            if ismagr and ismagk:
                magmsg = '(R and K mag. used for selection)'
            elif ismagr:
                magmsg = '(R mag. used for selection)'
            elif ismagk:
                magmsg = '(K mag. used for selection)'
            else:
                magmsg = ''
            print('    {} gal. for {}-BPT plot at z={} {}\n'.format(
                np.shape(sel)[1],bpt,redshift,magmsg))
        
    return bptnoms


def make_testplots(fnom,zz,snap,verbose=True):    
    # Get output file for BPT plot
    nbpt, sbpt = test_bpts(fnom,zz,snap,verbose=verbose)
    
    return ' '
