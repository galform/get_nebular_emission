"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
.. contributions:: Olivia Vidal <ovive.pro@gmail.com>
.. contributions:: Julen Expósito-Márquez <expox7@gmail.com>
"""

import os.path
import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
import src.gne_io as io
from src.gne_stats import perc_2arrays,calculate_sigma_levels
import src.gne_const as c
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
#              colsGSM=[0,1,2,3], labelObs=None, specific=False, h0=c.h, volume=c.vol_pm, verbose=False):
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
#            lms = np.log10((10**data['lms'][:,0])/c.IMF_M['Chabrier']+10**data['lms'][:,1]*c.IMF_M['Top-heavy']/c.IMF_M['Chabrier']) #+ np.log10(h0)
#            if specific:
#                lsfr = np.log10(10**data['lssfr'][:,0]+10**data['lssfr'][:,1]) + 9
#            else: 
#                lsfr = np.log10(10**data['lssfr'][:,0]+10**data['lssfr'][:,1]) + lms
#                lsfr = lsfr/c.IMF_SFR['Chabrier']
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
#        sfrf = H / volume / ds # / c.h**-3
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
#        zz = np.zeros(shape=(len(shist), len(mhist))); zz.fill(c.notnum)
#        ind = np.where(smf > 0.)
#        zz[ind] = np.log10(smf[ind])
#        
#        # print(zz[ind])
#
#        ind = np.where(zz > c.notnum)
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
#    set_cosmology(omega0=c.omega0, omegab=c.omegab,lambda0=c.lambda0,h0=c.h)
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
#        minU, maxU = get_limits(propname='logUs', photmod=photmod)
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
#            if photmod not in c.photmods:
#                if verbose:
#                    print('STOP (gne_photio.test_bpt): Unrecognised model to get emission lines.')
#                    print('                Possible photmod= {}'.format(c.photmods))
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


def plot_umz(root, subvols=1, outpath=None, verbose=True):
    '''
    Make plots of the ionizing parameter versus stellar mass and Zgas,
    for SF regions and AGNs.
    
    Parameters
    ----------
    root : string
       Path to input files. 
    subvols: integer
        Number of subvolumes to be considered
    outpath : string
        Path to output, default is output/ 
    verbose : boolean
       If True print out messages.
    '''

    # Get redshift and model information from data
    filenom = root+'0.hdf5'
    f = h5py.File(filenom, 'r') 
    header = f['header']
    redshift = header.attrs['redshift']
    photmod_sfr = header.attrs['photmod_sfr']
    photmod_agn = header.attrs['photmod_agn']
    mp = header.attrs['mp_Msun']
    lms = f['data/lms'][:,:]
    f.close()

    # Get number of components
    ncomp = np.shape(lms)[1]
    
    # Prep plots
    fig = plt.figure(figsize=(30, 30))
    gs = fig.add_gridspec(2, 2, hspace=0, wspace=0)
    (axsz, axsm), (axaz, axam) = gs.subplots(sharex='col', sharey='row')
    
    # Read limits for properties and photoionisation models
    minUs, maxUs = get_limits(propname='logUs', photmod=photmod_sfr)
    umins = minUs - 0.5; umaxs = maxUs +0.5
    minZs1, maxZs1 = get_limits(propname='Z', photmod=photmod_sfr)
    minZs = np.log10(minZs1); maxZs = np.log10(maxZs1);
    
    minUa, maxUa = get_limits(propname='logUs', photmod=photmod_agn)
    umina = minUa - 0.5; umaxa = maxUa +0.5
    minZa1, maxZa1 = get_limits(propname='Z', photmod=photmod_agn)
    minZa = np.log10(minZa1); maxZa = np.log10(maxZa1);

    zmin = min(minZs,minZa)-0.5
    zmax = max(maxZs,maxZa)+0.5
    
    xtit = 'log$_{10}(Z)$'
    axaz.set_xlabel(xtit)
    ytit = 'log$_{10}U_{\\rm AGN}$'
    axaz.set_ylabel(ytit)
    axaz.set_xlim(zmin, zmax)
    axaz.set_ylim(umina, umaxa)

    ytit = 'log$_{10}U_{\\rm SF}$'
    axsz.set_ylabel(ytit)
    axsz.set_ylim(umins, umaxs)

    xtit = 'log$_{10}(M_*/M_{\\odot})$'
    axam.set_xlabel(xtit)

    # Squares with limits
    axsz.add_patch(plt.Rectangle((minZs, minUs), maxZs-minZs,
                                 maxUs-minUs,ec="r",fc="none"))
    axaz.add_patch(plt.Rectangle((minZa, minUa), maxZa-minZa,
                                 maxUa-minUa,ec="r",fc="none"))

    # Initialise counters per component
    tots, tota, ina, ins = [np.zeros(ncomp) for i in range(4)]
    
    # Read data in each subvolume
    for ivol in range(subvols):
        filenom = root+str(ivol)+'.hdf5' #; print(filenom); exit()
        f = h5py.File(filenom, 'r')
        lms1 = f['data/lms'][:]

        # Read SF information from file
        lusfr1 = f['sfr_data/lu_sfr'][:]
        lzsfr1 = f['sfr_data/lz_sfr'][:]

        # Read AGN information if it exists
        AGN = True
        if 'agn_data' not in f.keys(): AGN = False
        if AGN: 
            luagn1 = f['agn_data/lu_agn'][:]
            lzagn1 = f['agn_data/lz_agn'][:]
        f.close()

        if ivol == 0:
            lms = lms1; lusfr = lusfr1; lzsfr = lzsfr1
            if AGN: luagn = luagn1; lzagn = lzagn1
        else:
            lms = np.append(lms,lms1,axis=0)
            lusfr = np.append(lusfr,lusfr1,axis=0)
            lzsfr = np.append(lzsfr,lzsfr1,axis=0)            
            if AGN:
                luagn = np.append(luagn,luagn1,axis=0)
                lzagn = np.append(lzagn,lzagn1,axis=0)            
            
        # Check number of galaxies within model limits
        if (len(lusfr) != len(lzsfr)):
            print('WARNING plots.umz, SFR: different length arrays U and Z')
        if AGN:
            if (len(luagn) != len(lzagn)):
                print('WARNING plots.umz, AGN: different length arrays U and Z')

        # Loop over components 
        for i in range(ncomp):
            u = lusfr[:,i]
            z = lzsfr[:,i]
            tots[i] = tots[i] + len(u)
            ind = np.where((u>=minUs) & (u<=maxUs) &
                           (z>=minZs) & (z<=maxZs))
            ins[i] = ins[i] + np.shape(ind)[1]
            
            if AGN:
                u = luagn[:,i]
                z = lzagn[:,i]
                tota[i] = tota[i] + len(u)
                ind = np.where((u>=minUa) & (u<=maxUa) &
                               (z>=minZa) & (z<=maxZa))
                ina[i] = ina[i] + np.shape(ind)[1]

    # Set limits for the mass
    
    # Bins for percentiles
    dx = 0.2

    zbins = np.arange(zmin, (zmax + dx), dx)
    zhist = zbins + dx * 0.5

    mmin = max(np.log10(mp),np.min(lms[lms>0]))
    mmax = np.max(lms[lms>0])
    mbins = np.arange(mmin, (mmax + dx), dx)
    mhist = mbins + dx * 0.5
    axam.set_xlim(mmin, mmax)
    
    # Find percentages within photoionisation model limits
    for i in range(ncomp):
        col = np.array([plt.cm.tab20(float(i)/ncomp)])
        
        leg = "{} components {}, {:.1f}% in".format(int(tots[i]),i,
                                                 ins[i]*100./tots[i])
        m = lms[:,i]
        u = lusfr[:,i]
        z = lzsfr[:,i]

        ind = np.where((u>c.notnum) & (z>c.notnum))
        if (np.shape(ind)[1]>len(zbins)):
            for isub in range(2):
                arr = z[ind]
                if isub == 1: arr = m[ind]

                med = perc_2arrays(zbins, arr, u[ind], 0.5)        
                upq = perc_2arrays(zbins, arr, u[ind], 0.84)
                low = perc_2arrays(zbins, arr, u[ind], 0.16)
                jnd = np.where(med>c.notnum)
                if (np.shape(jnd)[1]>1):
                    el = med[jnd] - low[jnd]
                    eh = upq[jnd] - med[jnd]
                    if isub == 0:
                        axsz.errorbar(zhist[jnd],med[jnd],
                                      yerr=[el,eh],c=col,label=leg)
                    elif isub == 1:
                        axsm.errorbar(mhist[jnd],med[jnd],yerr=[el,eh],c=col)
        else:
            axsz.scatter(z[ind],u[ind],c=col,label=leg)
            axsm.scatter(m[ind],u[ind],c=col)                
            
        if AGN:
            leg = "{} components {}, {:.1f}% in".format(int(tota[i]),i,
                                                     ina[i]*100./tota[i])
            u = luagn[:,i]
            z = lzagn[:,i]

            ind = np.where((u>c.notnum) & (z>c.notnum))
            if (np.shape(ind)[1]>len(zbins)):
                for isub in range(2):
                    arr = z[ind]
                    if isub == 1: arr = m[ind]

                    med = perc_2arrays(zbins, arr, u[ind], 0.5)        
                    upq = perc_2arrays(zbins, arr, u[ind], 0.84)
                    low = perc_2arrays(zbins, arr, u[ind], 0.16)
                    jnd = np.where(med>c.notnum)
                    if (np.shape(jnd)[1]>1):
                        el = med[jnd] - low[jnd]
                        eh = upq[jnd] - med[jnd]
                        if isub == 0:
                            axaz.errorbar(zhist[jnd],med[jnd],
                                          yerr=[el,eh],c=col,label=leg)
                        elif isub == 1:
                            axam.errorbar(mhist[jnd],med[jnd],
                                          yerr=[el,eh],c=col)
            else:
                axaz.scatter(z[ind],u[ind],c=col,label=leg)
                axam.scatter(m[ind],u[ind],c=col)                
            

    # Legend
    leg = axsz.legend(loc=1); leg.draw_frame(False)
    leg = axaz.legend(loc=1); leg.draw_frame(False)
   
    # Output
    pltpath = io.get_plotpath(root)
    plotnom = pltpath+'umz.pdf'
    plt.savefig(plotnom)
    if verbose:
         print(f'* U plots: {plotnom}')
    
    return plotnom


def get_obs_bpt(redshift,bpt):
    '''
    Get observational data for BPT diagrams at a given redshift
    
    Parameters
    ----------
    redshift : float
       Redshift of interest
    bpt: string
        Type of BPT diagram: 'NII' (OIII/Hbeta vs N2/Ha)
        or 'SII' (OIII/Hbeta vs S2/Ha)

    Returns
    -------
    xobs, yobs : array of floats
       Ratios for each observed spectral emission line
    obsdata : boolean
       True if there is any observational data at the given redshift
    '''

    xobs = -999.; yobs = -999.; obsdata = False
    
    # Use different data sets for different redshifts
    if redshift <= 0.2:
        obsdata = True
        obsfile = 'data/observational_data/favole2024.txt'
        l1,l2 = np.loadtxt(obsfile,skiprows=1,usecols=(15,9),unpack=True)
        xx, yy = [np.zeros(len(l1)) for i in range(2)]
        ind = np.where((l1>0.) & (l2>0.))
        if (np.shape(ind)[1]>0): #O3/Hb
            yy[ind] = np.log10(l1[ind]/l2[ind])
            
        if bpt=='NII': #N2/Ha
            l1,l2 = np.loadtxt(obsfile,skiprows=1,usecols=(18,6),unpack=True)
            ind = np.where((l1>0.) & (l2>0.))
            if (np.shape(ind)[1]>0):
                xx[ind] = np.log10(l1[ind]/l2[ind]) 
        elif bpt=='SII': #S2/Ha
            l1,l2 = np.loadtxt(obsfile,skiprows=1,usecols=(21,6),unpack=True)
            ind = np.where((l1>0.) & (l2>0.))
            if (np.shape(ind)[1]>0):
                xx[ind] = np.log10(l1[ind]/l2[ind]) 

    elif 1.45 <= redshift <= 1.75:
        obsdata = True
        if bpt=='NII':
            obsfile = 'data/observational_data/NII_Kashino.txt'
            yy = np.loadtxt(obsfile,skiprows=18,usecols=(6)) #O3/Hb
            xx = np.loadtxt(obsfile,skiprows=18,usecols=(3)) #N2/Ha
                
        elif bpt=='SII':
            obsfile = 'data/observational_data/SII_Kashino.txt'
            yy = np.loadtxt(obsfile,skiprows=18,usecols=(6)) #O3/Hb
            xx = np.loadtxt(obsfile,skiprows=18,usecols=(3)) #N2/Ha

    if obsdata:
        ind = np.where((xx>c.notnum) & (yy>c.notnum))
        if (np.shape(ind)[1]>0):
            xobs = xx[ind]
            yobs = yy[ind]
        else:
            obsdata = False

    return xobs,yobs,obsdata


def plot_bpts(root, subvols=1, outpath=None, verbose=True):
    '''
    Make the 2 BPT diagrams without attenuation
    
    Parameters
    ----------
    root : string
       Path to input files. 
    subvols: integer
        Number of subvolumes to be considered
    outpath : string
        Path to output, default is output/ 
    verbose : boolean
       If True print out messages.
    '''

    # Get redshift and cosmology from data
    filenom = root+'0.hdf5'
    f = h5py.File(filenom, 'r') 
    header = f['header']
    redshift = header.attrs['redshift']
    omega0 = header.attrs['omega0']
    omegab = header.attrs['omegab']
    lambda0 = header.attrs['lambda0']
    h0 = header.attrs['h0']
    photmod_sfr = header.attrs['photmod_sfr']
    photmod_agn = header.attrs['photmod_agn']

    f.close()
    
    # Set the cosmology from the simulation
    set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0,h0=h0)

    # Read limits for properties and photoionisation models
    minU, maxU = get_limits(propname='logUs', photmod=photmod_sfr)
    minZ, maxZ = get_limits(propname='Z', photmod=photmod_sfr)

    # Prep plots
    fig, (axn, axs) = plt.subplots(1, 2, figsize=(30, 15),
                                   layout='constrained')
    ytit = 'log$_{10}$([OIII]$\\lambda$5007/H$\\beta$)'
    xmins = [-1.9,-1.9]
    xmaxs = [0.8,0.9]
    ymins = [-1.5,-2.1]
    ymaxs = [1.5,1.6]
    for ii, bpt in enumerate(['NII','SII']):
        if bpt=='NII':
            xtit = 'log$_{10}$([NII]$\\lambda$6584/H$\\alpha$)'
            axn.set_xlim(xmins[ii], xmaxs[ii])
            axn.set_ylim(ymins[ii], ymaxs[ii])
            axn.set_xlabel(xtit); axn.set_ylabel(ytit)
        elif bpt=='SII':
            xtit = 'log$_{10}$([SII]$\\lambda$6584/H$\\alpha$)'
            axs.set_xlim(xmins[ii], xmaxs[ii])
            axs.set_ylim(ymins[ii], ymaxs[ii])
            axs.set_xlabel(xtit); axs.set_ylabel(ytit)

        xobs, yobs, obsdata = get_obs_bpt(redshift,bpt)
        if obsdata and bpt=='NII':
            #xi,yi,z,levels = calculate_sigma_levels(xobs,yobs,n_grid=100,n_levels=3)
            #contour = axn.contourf(xi, yi, z, levels=levels,
            #           colors=['darkgrey']*3)
            #           #alpha=[0.8, 0.5, 0.2])
            axn.scatter(xobs,yobs, s=20, c='darkgrey', alpha=0.8)                       
        elif obsdata and bpt=='SII':
            axs.scatter(xobs,yobs, s=20, c='darkgrey', alpha=0.8)

    # Read data in each subvolume and add data to plots
    seltot = 0
    for ivol in range(subvols):
        filenom = root+str(ivol)+'.hdf5'
        f = h5py.File(filenom, 'r')
        
        # Read SF information from file
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
            # Read AGN information from file
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
        
        # Combine luminosities
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

        ind = np.where((Ha>0)   & (Hb>0)  & 
                       (NII>0)  & (OII>0) &
                       (OIII>0) & (SII>0) &
                       (lu_sfr>minU)&(lu_sfr<maxU)&
                       (lz_sfr>np.log10(minZ))&(lz_sfr<np.log10(maxZ)))
        if (np.shape(ind)[1] < 1):
            print('STOP BPT plots: not enough adequate data')
            return None

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
        # Add further cuts if adequate
        if redshift <= 0.2:
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
            flux = 5e-17  # erg/s/cm^2 Kashino+2019
            Lmin = emission_line_luminosity(flux,redshift)*1e40 #erg/s
            
            if ismagk:
                sel = np.where((Ha> Lmin) & (mag_k<23.5))
            else:
                sel = np.where(Ha > Lmin)
                
        if (np.shape(sel)[1]<1):
            continue
        seltot = seltot + np.shape(sel)[1]
                
        # Model spectral line ratios
        yy = O3Hb[sel] #O3/Hb
        cha = Halpha_ratio[sel]
        
        xx = N2Ha[sel] #N2/Ha
        axn.scatter(xx,yy, c=cha,s=50, marker='o', cmap=cmap)

        xx = S2Ha[sel] #S2/Ha
        axs.scatter(xx,yy, c=cha,s=50, marker='o', cmap=cmap)

        # Join all data for the colourbar
        if ivol == 0:
            chatot = cha
        else:
            chatot = np.append(chatot,cha)

    # Add colorbar
    sm = cm.ScalarMappable(cmap=cmap) # Create ScalarMappable
    sm.set_array(chatot)    
    cbar = plt.colorbar(sm, ax=axs, cmap=cmap, location='right')
    if AGN:
        collabel = r'$L_{\rm H_{\alpha}, AGN}/L_{\rm H_{\alpha}, tot}$'
    else:
        collabel = r'$L_{\rm H_{\alpha}}$'
    cbar.set_label(collabel,rotation=270,labelpad=60)        

    if verbose:
        if ismagr and ismagk:
            magmsg = '(R and K mag. used for selection)'
        elif ismagr:
            magmsg = '(R mag. used for selection)'
        elif ismagk:
            magmsg = '(K mag. used for selection)'
        else:
            magmsg = ''
        print('    {} gal. for BPT plots at z={} {}\n'.format(
            seltot,redshift,magmsg))

    for ii, bpt in enumerate(['NII','SII']):
        # Lines
        xline = np.arange(xmins[ii],xmaxs[ii]+0.1, 0.03)
        if bpt=='NII':
            yline = lines_BPT(xline,bpt,'Kauffmann2003')
            axn.plot(xline,yline,'k.')

            yline = lines_BPT(xline,bpt,'Kewley2001')
            axn.plot(xline,yline,'k-')
            
        elif bpt=='SII':
            yline = lines_BPT(xline,bpt,'Kewley2001')
            axs.plot(xline,yline,'k-')

            ylinel = lines_BPT(xline,bpt,'Kewley2006')
            axs.plot(xline[ylinel>yline],ylinel[ylinel>yline],'k--')

    # Output
    pltpath = io.get_plotpath(root)
    bptnom = pltpath+'bpts.pdf'
    plt.savefig(bptnom)
    if verbose:
         print(f'* BPT plots: {bptnom}')
    
    return bptnom


def make_testplots(rootf,snap,subvols=1,outpath=None,verbose=True):
    '''
    Make test plots
    
    Parameters
    ----------
    rootf : string
       Path and root to input files
    snap: integer
        Simulation snapshot number
    subvols: integer
        Number of subvolumes to be considered
    outpath : string
        Path to output, default is output/ 
    verbose : boolean
       If True print out messages.
    '''

    root = io.get_outroot(rootf,snap,outpath=outpath,verbose=True)

    # Get output file for BPT plot
    umz = plot_umz(root,subvols=subvols,verbose=verbose)
    
    # Get output file for BPT plot
    bpt = plot_bpts(root,subvols=subvols,verbose=verbose)
    
    return
