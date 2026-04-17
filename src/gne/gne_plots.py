"""
.. moduleauthor:: Violeta Gonzalez-Perez <violetagp@protonmail.com>
.. contributions:: Olivia Vidal <ovive.pro@gmail.com>
.. contributions:: Julen Expósito-Márquez <expox7@gmail.com>
"""

import os.path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.colors as mcol
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.offsetbox as moffbox

import warnings
warnings.filterwarnings('ignore', message='Input line .* contained no data')

import gne.gne_const as c
import gne.gne_io as io
import gne.gne_stats as st
from gne.gne_stats import n_gt_x
from gne.gne_photio import get_limits,read_gutkin16_grids,read_feltre16_grids
from gne.gne_cosmology import set_cosmology
from gne.gne_flux import flux2L
import gne.gne_style
plt.style.use(gne.gne_style.style1)

cmap = 'jet'
n4contour = 1000
min_Lbol = 42 # Based on Griffin+2020, fig 14
max_Lbol = 50
min_Ms = 8    # To be obtained from sim. res. ###here
max_Ms = 12   # To be obtained from sim. res. ###here

markers = ['o','^', 's', '*','D', 'p', 'h', 'H', '+', 'x', 'v', '<', '>', '|', '_']

def contour2Dsigma(n_levels=None,color='darkgrey'):
    '''
    Get levels following the standard deviation numbers expected for
    a 2D-Gaussian distribution. Generate colours varying in intensity.
    '''
    if n_levels is not None:
        levels=c.sigma_2Dprobs[0:n_levels]
    else:
        levels=c.sigma_2Dprobs.copy()
        
    nl = len(levels); levels.insert(0,0)
    alphas = np.linspace(0.2, 1, nl)[::-1].tolist()
    colors = [(*mcol.to_rgba(color, alpha=a),)
              for a in alphas]

    return levels,colors


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


def plot_comp_contour(ax, xx, yy, tots, ins, cm=plt.cm.tab20):
    """
    Plot components as a contour or scatter plot,
    on given axis and return legend elements.
    
    Parameters:
    -----------
    ax : matplotlib axis
        Axis to plot on
    xx : ndarray
        X vayyes for each component
    yy : ndarray
        Y vayyes for each component
    tots : ndarray
        Total number of elements for each component
    ins : ndarray
        Number within limits for each component
    cm : matplotlib colormap, optional
        Colormap to use
        
    Returns:
    --------
    proxies : list
        List of proxy artists for legend
    labels : list
        List of labels for legend
    """
    proxies = []; labels = []

    n_comp = np.shape(xx)[1]
    for i in range(n_comp):
        ind = np.where((xx[:, i] > c.notnum) & (yy[:, i] > c.notnum))[0]
        if len(ind) > 0:
            x = xx[ind, i]
            y = yy[ind, i]
            col = np.array([cm(float(i) / n_comp)])
            
            if len(ind) > n4contour:
                xc, yc, zc = st.get_cumulative_2Ddensity(x, y, n_grid=100)
                levels, colors = contour2Dsigma(color=col)
                contour = ax.contourf(xc, yc, zc, levels=levels, colors=colors)
                proxies.append(plt.Rectangle((0, 0), 1, 1, fc=col[0]))
            else:
                scatter = ax.scatter(x, y, c=col)
                proxies.append(scatter)
            
            leg = "{} component {} ({:.1f}% in)".format(
                int(tots[i]), i, ins[i]*100./tots[i])
            labels.append(leg)
    
    return proxies, labels


def plot_comp_quartiles(ax, xx, yy, xmin, xmax, tots, ins, cm=plt.cm.tab20):
    """
    Plot quartiles for components and return legend elements.
    
    Parameters:
    -----------
    ax : matplotlib axis
        Axis to plot on
    xin : ndarray
        X values for each component
    yin : ndarray
        Y values for each component
    n_comp : int
        Number of components
    xmin, xmax: float
        Limits for the x-axis
    tots : ndarray
        Total number of elements for each component
    ins : ndarray
        Number within limits for each component
    cm : matplotlib colormap, optional
        Colormap to use
        
    Returns:
    --------
    proxies : list
        List of proxy artists for legend
    labels : list
        List of labels for legend
    """
    proxies = []; labels = []
    
    dx = 0.2    
    xbins = np.arange(xmin,xmax + dx, dx)
    xhist = xbins + dx * 0.5
    ax.set_xlim(xmin, xmax)

    n_comp = np.shape(xx)[1]
    for i in range(n_comp):
        ind = np.where((xx[:, i] > c.notnum) & (yy[:, i] > c.notnum))[0]
        if len(ind) > 0:
            x = xx[ind, i]
            y = yy[ind, i]
            col = np.array([cm(float(i) / n_comp)])

            med = st.perc_2arrays(xbins, x, y, 0.5)        
            upq = st.perc_2arrays(xbins, x, y, 0.84)
            low = st.perc_2arrays(xbins, x, y, 0.16)
            jnd = np.where(med>c.notnum)
            if (np.shape(jnd)[1]>1):
                el = med[jnd] - low[jnd]
                eh = upq[jnd] - med[jnd]
                quart = ax.errorbar(xhist[jnd],med[jnd],yerr=[el,eh],c=col)
                proxies.append(quart)
                
                leg = "{} component {} ({:.1f}% in)".format(
                    int(tots[i]), i, ins[i]*100./tots[i])
                labels.append(leg)
                
    return proxies, labels


def plot_uzn(root, endf, subvols=1, outpath=None, verbose=True):
    '''
    Make plots of the ionizing parameter versus Zgas,
    and electron density as a function of stellar mass,
    for SF regions and NLR AGNs, if calculated.
    
    Parameters
    ----------
    root : string
       Path to files with calculated data (lines, etc)
    endf : string
       Ending of input files. 
    subvols: integer or list of integers
        Number of subvolumes to be considered
    outpath : string
        Path to output, default is output/ 
    verbose : boolean
       If True print out messages.
    '''
    # Get redshift and model information from data
    filenom = os.path.join(root+'0',endf)
    f = h5py.File(filenom, 'r') 
    header = f['header']
    redshift = header.attrs['redshift']
    photmod_sfr = header.attrs['photmod_sfr']
    mp = header.attrs['mp_Msun']
    lu = f['sfr_data/lu_sfr'][:]
    # Read AGN information if it exists
    if 'agn_data' not in f.keys():
        AGN = False
    else:
        AGN = True
        photmod_agn = header.attrs['photmod_NLR']
        lua = f['agn_data/lu_agn'][:]

    # Get number of components
    ncomp = np.shape(lu)[1]

    # Prep plots
    side = 15
    if AGN:
        fig, ((axu,axn),(axua,axna)) = plt.subplots(2, 2,
                                                    figsize=(2*side, 2*side),
                                                    layout='constrained')
    else:
        fig, (axu,axn) = plt.subplots(1, 2, figsize=(2*side, side),
                                 layout='constrained')

    axu.set_ylabel('log$_{10}U_{\\rm SF}$')
    axu.set_xlabel('log$_{10}(Z_{\\rm gas})$')
    axn.set_ylabel('log$_{10}(n_{H, \\rm SFR})$')
    axn.set_xlabel('log$_{10}(M_{*})$')    
    if AGN:
        axua.set_ylabel('log$_{10}U_{\\rm AGN}$')
        axua.set_xlabel('log$_{10}(Z_{\\rm gas})$')
        axna.set_ylabel('log$_{10}(n_{H, \\rm AGN})$')
        axna.set_xlabel('log$_{10}(L_{\\rm AGN}/{\\rm erg/s})$')    

    # Read limits for photoionisation models
    pad = 0.5
    umin, umax = get_limits(propname='logUs', photmod=photmod_sfr)
    axu.set_ylim(umin-pad, umax+pad)
    
    zmin, zmax = np.log10(get_limits(propname='Z', photmod=photmod_sfr))
    axu.set_xlim(zmin-pad, zmax+pad)

    axu.add_patch(plt.Rectangle((zmin, umin), zmax-zmin,umax-umin,
                                ec="gray",ls='--',lw=10,fc="none"))

    nmin, nmax = np.log10(get_limits(propname='nH', photmod=photmod_sfr))
    axn.set_ylim(nmin-pad, nmax+pad)

    mmin = min_Ms-pad; mmax = max_Ms+pad 
    axn.plot([mmin, mmax], [nmin, nmin], 'gray', ls='--', lw=10)
    axn.plot([mmin, mmax], [nmax, nmax], 'gray', ls='--', lw=10)

    if AGN:
        uamin, uamax = get_limits(propname='logUs', photmod=photmod_agn)
        axua.set_ylim(uamin-pad, uamax+pad)
    
        zamin, zamax = np.log10(get_limits(propname='Z', photmod=photmod_agn))
        axua.set_xlim(zamin-pad, zamax+pad)

        axua.add_patch(plt.Rectangle((zamin, uamin), zamax-zamin,uamax-uamin,
                                     ec="gray",ls='--',lw=10,fc="none"))

        namin, namax = np.log10(get_limits(propname='nH', photmod=photmod_agn))
        axna.set_ylim(namin-pad, namax+pad)

        mmin = min_Lbol-pad; mmax = max_Lbol+pad 
        axna.plot([mmin, mmax], [namin, namin], 'gray', ls='--', lw=10)
        axna.plot([mmin, mmax], [namax, namax], 'gray', ls='--', lw=10)

    # Initialise counters per component
    tots, ins, inns = [np.zeros(ncomp) for i in range(3)]
    if AGN:
        tota, ina, inna = [np.zeros(1) for i in range(3)]

    # Read data in each subvolume
    list_subvols = subvols
    if isinstance(subvols, int):
        list_subvols = list(range(subvols))

    first_vol = True

    for ivol in list_subvols:
        filenom = os.path.join(root+'0',endf) #; print(filenom); exit()
        f = h5py.File(filenom, 'r'); header = f['header']

        # Read information from file
        lms1 = f['data/lms'][:]
        lzsfr1 = f['sfr_data/lz_sfr'][:]
        lusfr1 = f['sfr_data/lu_sfr'][:]
        lnsfr1 = f['sfr_data/lnH_sfr'][:]
        if AGN:
            Lagn1  = f['agn_data/Lagn'][:]
            lzagn1 = f['agn_data/lz_agn'][:]
            luagn1 = f['agn_data/lu_agn'][:]
            if 'epsilon_NLR' in header.attrs:
                epsilon_is_constant = True
                epsilon1 = header.attrs['epsilon_NLR']
            else:
                epsilon_is_constant = False
                epsilon1 = f['agn_data/epsilon_NLR'][:]
        f.close()

        if first_vol:
            lusfr = lusfr1; lzsfr = lzsfr1
            lnsfr = lnsfr1; lms = lms1
            if AGN:
                luagn = luagn1; lzagn = lzagn1
                Lagn = Lagn1
                if not epsilon_is_constant:
                    epsilon = epsilon1    
            first_vol = False
        else:
            lusfr = np.append(lusfr,lusfr1,axis=0)
            lzsfr = np.append(lzsfr,lzsfr1,axis=0)
            lnsfr = np.append(lnsfr,lnsfr1,axis=0)
            lms = np.append(lms,lms1,axis=0)
            if AGN:
                luagn = np.append(luagn,luagn1,axis=0)
                lzagn = np.append(lzagn,lzagn1,axis=0)
                Lagn = np.append(Lagn,Lagn1,axis=0)
                if not epsilon_is_constant:
                    epsilon = np.append(epsilon,epsilon1,axis=0)
                else:
                    epsilon = np.zeros(Lagn.shape); epsilon.fill(epsilon1)

        # Check number of galaxies within model limits
        if (len(lusfr) != len(lzsfr)):
            print('WARNING plots.uzn, SFR: different length arrays U and Z')
        if AGN:
            if (len(luagn) != len(lzagn)):
                print('WARNING plots.uzn, AGN: different length arrays U and Z')
    
        # Count parameters within the limits of photoionising models
        for i in range(ncomp):
            mask = (lusfr[:,i] > c.notnum) & (lzsfr[:,i] > c.notnum)
            u = lusfr[mask,i]
            z = lzsfr[mask,i]
            tots[i] = tots[i] + len(u)
            ind = np.where((u>=umin) & (u<=umax) &
                           (z>=zmin) & (z<=zmax))
            ins[i] = ins[i] + np.shape(ind)[1]        

            mask = (lnsfr[:,i] > c.notnum)
            nn = lnsfr[mask,i]            
            ind = np.where((nn>=nmin) & (nn<=nmax))
            inns[i] = inns[i] + np.shape(ind)[1]        

        if AGN:
            mask = (luagn[:] > c.notnum) & (lzagn[:] > c.notnum)
            u = luagn[mask]
            z = lzagn[mask]
            tota = tota + len(u)
            ind = np.where((u>=uamin) & (u<=uamax) &
                           (z>=zamin) & (z<=zamax))
            ina = ina + np.shape(ind)[1]

    ###here to check as not working
    ## Plot per component U versus Z
    #proxies, labels = plot_comp_contour(axu, lzsfr, lusfr, tots, ins)
    #if AGN:
    #    aproxies, alabels = plot_comp_contour(axua, lzagn, luagn, tota, ina)
    #
    ## Legend for U vs Z
    #leg = axu.legend(proxies, labels, loc=0); leg.draw_frame(False)
    #if AGN:
    #    leg = axua.legend(aproxies, alabels, loc=0); leg.draw_frame(False)
    #print(leg); exit()
    ###here end
    ## Plot per component nH versus M* (or Lagn)
    #proxies, labels = plot_comp_quartiles(axn, lms, lnsfr,
    #                                      min_Ms, max_Ms, tots, inns)
    #if AGN:
    #    col = np.zeros(Lagn[:,0].shape); col.fill(c.notnum)
    #    mask = Lagn[:,0] > 0
    #    col[mask] = np.log10(Lagn[mask,0])
    #    lLagn = np.repeat(col[:, np.newaxis], nacomp, axis=1)
    #
    #    aproxies, alabels = plot_comp_quartiles(axna, lLagn, lnagn,
    #                                            min_Lbol, max_Lbol, tota, inna)

    ## Legend for nH plots
    #leg = axn.legend(proxies,labels, loc=0); leg.draw_frame(False)
    ##if AGN:
    ##    leg = axna.legend(aproxies,alabels, loc=0); leg.draw_frame(False)
        
    # Output
    plotnom = io.get_plotfile(root,endf,'uzn')
    plt.savefig(plotnom)
    if verbose:
         print(f'* U plots: {plotnom}')

#    pltpath = io.get_plotpath(root)
#    plotnom = pltpath+'uzn.pdf'
#    plt.savefig(plotnom)
#    if verbose:
#         print(f'* U plots: {plotnom}')
    
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
        obsfile = os.path.join(c.obs_data_dir,'favole2024.txt')
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
            obsfile = os.path.join(c.obs_data_dir,'NII_Kashino.txt')
            yy = np.loadtxt(obsfile,skiprows=18,usecols=(6)) #O3/Hb
            xx = np.loadtxt(obsfile,skiprows=18,usecols=(3)) #N2/Ha
                
        elif bpt=='SII':
            obsfile = os.path.join(c.obs_data_dir,'SII_Kashino.txt')
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


def plot_model_bpt_grids(photmod='gutkin16',xid=0.3,co=1,imf_cut=100,
                         alpha=-1.7,verbose=True):
    '''
    Plot photoionisation grids on 2 BPT diagrams.
    
    Parameters
    ----------
    photmod : string
       Name of the photoionisation model to be plotted
    verbose : boolean
       If True print out messages.

    Return
    ------
    outpath : string
       Name of output plot within output/photoio_grids
    '''
    # Prep plots
    fig, (axn, axs) = plt.subplots(1, 2, figsize=(32, 17),
                                   layout='constrained')
    #plt.subplots_adjust(right=0.85, top=0.9) 
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
            xtit = 'log$_{10}$([SII]$\\lambda\\lambda$6717,6731/H$\\alpha$)'
            axs.set_xlim(xmins[ii], xmaxs[ii])
            axs.set_ylim(ymins[ii], ymaxs[ii])
            axs.set_xlabel(xtit); axs.set_ylabel(ytit)

        xobs, yobs, obsdata = get_obs_bpt(0.,bpt)
        if obsdata and bpt=='NII':
            x,y,z = st.get_cumulative_2Ddensity(xobs,yobs,n_grid=100)
            levels,colors= contour2Dsigma()
            contour = axn.contourf(x, y, z, levels=levels,colors=colors)
        elif obsdata and bpt=='SII':
            x,y,z = st.get_cumulative_2Ddensity(xobs,yobs,n_grid=100)
            levels,colors= contour2Dsigma()
            contour = axs.contourf(x, y, z, levels=levels,colors=colors)

    for ii, bpt in enumerate(['NII','SII']):
        # Lines
        xline = np.arange(xmins[ii],xmaxs[ii]+0.1, 0.03)
        if bpt=='NII':
            yline = lines_BPT(xline,bpt,'Kauffmann2003')
            axn.plot(xline,yline,'k--')

            yline = lines_BPT(xline,bpt,'Kewley2001')
            axn.plot(xline,yline,'k-')
            
        elif bpt=='SII':
            yline = lines_BPT(xline,bpt,'Kewley2001')
            axs.plot(xline,yline,'k-')

            ylinel = lines_BPT(xline,bpt,'Kewley2006')
            axs.plot(xline[ylinel>yline],ylinel[ylinel>yline],'k-.')


    # Read grids of photoionisation models
    if photmod == 'gutkin16':
        grid1,grid2,grid3,grid4 = read_gutkin16_grids(xid, co, imf_cut)
        grids = [grid1, grid2, grid3, grid4]
        nl = 5
        col_ha = 10-nl
        col_hb = 6-nl 
        col_o3 = 8-nl     #[OIII]5007
        col_n2 = 11-nl    #[NII]6584
        col_s2_a = 12-nl  #[SII]6717
        col_s2_b = 13-nl  #[SII]6731        
    elif photmod == 'feltre16':
        grid1,grid2,grid3 = read_feltre16_grids(xid, alpha)
        grids = [grid1, grid2, grid3]
        nl = 4
        col_ha = 10-nl
        col_hb = 5-nl 
        col_o3 = 7-nl     #[OIII]5007
        col_n2 = 11-nl    #[NII]6584  
        col_s2_a = 12-nl  #[SII]6717
        col_s2_b = 13-nl  #[SII]6731

    for i, grid in enumerate(grids):
        nz = grid.shape[0]
        nu = grid.shape[1]

        size = 100 + i * 30
        for iz in range(nz):
            color = cm.tab20(iz % 20)
            
            for iu in range(nu):
                marker = markers[iu % len(markers)]
                
                el = grid[iz,iu,:]
                y = np.log10(el[col_o3]/el[col_hb])
                for ii, bpt in enumerate(['NII','SII']):
                    if bpt=='NII':
                        x = np.log10(el[col_n2]/el[col_ha])
                        axn.scatter(x,y,s=size,c=[color],marker=marker)
                    elif bpt=='SII':
                        s2 = el[col_s2_a] + el[col_s2_b]
                        x = np.log10(s2/el[col_ha])
                        axs.scatter(x,y,s=size,c=[color],marker=marker)

    # Add legend on Z values
    zvals = c.zmet_str[photmod]
    z_handles = []
    for iz, z_val in enumerate(zvals):
        color = cm.tab20(iz % 20)
        z_handle = mpatches.Patch(color=color, label=f'0.{z_val}')
        z_handles.append(z_handle)
    legend1 = fig.legend(handles=z_handles, loc='center left',
                         bbox_to_anchor=(1.02, 0.5), 
                         title="Z$_{gas}$", frameon=False)
    fig.add_artist(legend1)

    # Add legend on U values
    uvals = c.lus_bins[photmod]
    u_handles = []
    for iu, u_val in enumerate(uvals):
        marker = markers[iu % len(markers)]
        u_handle = mlines.Line2D([], [], color='black',
                                 marker=marker, ls='None',
                                 markersize=10, label=f'{u_val}')
        u_handles.append(u_handle)
    legend2 = fig.legend(handles=u_handles, 
                         loc='upper center',
                         bbox_to_anchor=(0.5, 1.1),
                         title='log$_{10}$ U',
                         ncol=min(5, len(u_handles)),
                         frameon=False)
    legend2._legend_box.align = "left"
    legend2._legend_box.sep = 7 
    plt.setp(legend2.get_title(), ha='right')
    legend2._legend_box = moffbox.VPacker(
        pad=0, sep=0, align="left",
        children=[
            moffbox.HPacker(
                pad=0, sep=5, align="center",
                children=[legend2._legend_box.get_children()[0],
                          legend2._legend_box.get_children()[1]])])
    fig.add_artist(legend2)
    
    # Add legend on model information
    if photmod == 'gutkin16':
        legend_model = (f'Gutkin+16\n'
                        f'$\\xi_d$ = {xid}\n'
                        f'C/O = {co} (C/O)$_\\odot$\n'
                        f'M(IMF)$<{imf_cut}$ M$_\\odot$')
    elif photmod == 'feltre16':
        legend_model = (f'Feltre+16\n'
                        f'$\\xi_d$ = {xid}\n'
                        f'$\\alpha$ = {alpha}\n')
    axn.text(0.05, 0.97, legend_model, transform=axn.transAxes,
             verticalalignment='top')

    # Output
    pltpath = 'output/plots/photoio_grids/'
    io.create_dir(pltpath)
    if (photmod == 'gutkin16'):
        bptnom = pltpath+photmod+'_xi'+str(xid)+'_bpts.pdf'
    elif (photmod == 'feltre16'):
        bptnom = pltpath+photmod+'_alpha'+str(abs(alpha))+'_bpts.pdf'

    plt.savefig(bptnom)
    if verbose:
         print(f'* Photoionisation model grids on BPT plots: {bptnom}')
    
    return bptnom


def plot_bpts(root, endf, subvols=[0], outpath=None,
              metadata=None,verbose=True):
    '''
    Make the 2 BPT diagrams without attenuation
    
    Parameters
    ----------
    root : string
       Path to input files. 
    endf : string
       Ending of input files. 
    subvols: List of integers
        List of subvolumes to be considered
    outpath : string
        Path to output, default is output/
    metadata : dictionary
        Cosmology and other metadata information
    verbose : boolean
       If True print out messages.
    '''
    
    # Get metadata
    photmod_sfr = metadata['photmod_sfr']
    AGN = metadata['AGN']
    if AGN:
        photmod_agn = metadata['photmod_agn']
    redshift = metadata['redshift']

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
            xtit = 'log$_{10}$([SII]$\\lambda\\lambda$6717,6731/H$\\alpha$)'
            axs.set_xlim(xmins[ii], xmaxs[ii])
            axs.set_ylim(ymins[ii], ymaxs[ii])
            axs.set_xlabel(xtit); axs.set_ylabel(ytit)

        xobs, yobs, obsdata = get_obs_bpt(redshift,bpt)
        if obsdata and bpt=='NII':
            x,y,z = st.get_cumulative_2Ddensity(xobs,yobs,n_grid=100)
            levels,colors= contour2Dsigma()
            contour = axn.contourf(x, y, z, levels=levels,colors=colors)
        elif obsdata and bpt=='SII':
            x,y,z = st.get_cumulative_2Ddensity(xobs,yobs,n_grid=100)
            levels,colors= contour2Dsigma()
            contour = axs.contourf(x, y, z, levels=levels,colors=colors)

    # Read data in each subvolume and add data to plots
    seltot = 0
    chatot = None
    for ivol in subvols:
        filenom = os.path.join(root+str(ivol),endf)
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
        if AGN:
            # Read AGN information from file
            lu_agn = f['agn_data/lu_agn'][:]
            lz_agn = f['agn_data/lz_agn'][:]
            Ha_agn = f['agn_data/Halpha_agn'][:]
            Hb_agn = f['agn_data/Hbeta_agn'][:]
            NII6548_agn = f['agn_data/NII6584_agn'][:]
            OII3727_agn = f['agn_data/OII3727_agn'][:]
            OIII5007_agn= f['agn_data/OIII5007_agn'][:]
            SII6731_agn = f['agn_data/SII6731_agn'][:]
            SII6717_agn = f['agn_data/SII6717_agn'][:]
        
        # Magnitudes for cuts
        ismagr = True
        try:
            magr = f['data/magR'][:]
        except:
            ismagr = False
        
        ismagk = True
        try:
            magk = f['data/magK'][:]
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
    
        O3Hb = np.log10(OIII) - np.log10(Hb)
        N2Ha = np.log10(NII) - np.log10(Ha)
        S2Ha = np.log10(SII) - np.log10(Ha)

        if ismagr:
            mag_r = magr[ind]
        if ismagk:
            mag_k = magk[ind]

        sel = (np.arange(len(O3Hb)),)
        # Add further cuts if adequate
        if redshift <= 0.2:
            flux = 2e-16 # erg/s/cm^2 Favole+2024
            Lmin = flux2L(flux,redshift) #erg/s

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
            Lmin = flux2L(flux,redshift) #erg/s
            
            if ismagr:
                sel = np.where((Ha> Lmin) & (mag_r<124.1))
            else:
                sel = np.where(Ha> Lmin)
        elif 1.45 <= redshift <= 1.75:
            flux = 5e-17  # erg/s/cm^2 Kashino+2019
            Lmin = flux2L(flux,redshift) #erg/s
            
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
        if chatot is None:
            chatot = cha
        else:
            chatot = np.append(chatot,cha)

    # Add colorbar
    sm = cm.ScalarMappable(cmap=cmap) # Create ScalarMappable
    sm.set_array(chatot)    
    cbar = plt.colorbar(sm, ax=axs, cmap=cmap, location='right')
    if AGN:
        collabel = (r'$L_{\rm H_{\alpha}, AGN}/L_{\rm H_{\alpha}, tot}$'+
                    f' (z={redshift:.1f})')
    else:
        collabel = r'$L_{\rm H_{\alpha}}$'+f' (z={redshift:.1f})'
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
        print(f'    {seltot} gal. for BPT plots at z={redshift:.1f} {magmsg}\n')

    for ii, bpt in enumerate(['NII','SII']):
        # Lines
        xline = np.arange(xmins[ii],xmaxs[ii]+0.1, 0.03)
        if bpt=='NII':
            yline = lines_BPT(xline,bpt,'Kauffmann2003')
            axn.plot(xline,yline,'k--')

            yline = lines_BPT(xline,bpt,'Kewley2001')
            axn.plot(xline,yline,'k-')
            
        elif bpt=='SII':
            yline = lines_BPT(xline,bpt,'Kewley2001')
            axs.plot(xline,yline,'k-')

            ylinel = lines_BPT(xline,bpt,'Kewley2006')
            axs.plot(xline[ylinel>yline],ylinel[ylinel>yline],'k-.')

    # Output
    bptnom = io.get_plotfile(root,endf,'bpt')
    plt.savefig(bptnom)
    if verbose:
         print(f'* BPT plots: {bptnom}')
    
    return bptnom


def plot_lf(root, endf, subvols=[0], outpath=None,
            outnom = 'masses',
            props=['data/mh','data/lm_s','data/lm_gas'],
            prop_labels=[r'M$_{\rm h}(M_{\odot})$',
                         r'M$_{\rm *}(M_{\odot})$',
                         r'M$_{\rm gas}(M_{\odot})$'],
            xmin=9,xmax=14.,dx=0.1,
            metadata=None,verbose=True):
    '''
    Make a (luminosity) function plot
    
    Parameters
    ----------
    root : string
       Path to input files. 
    endf : string
       Ending of input files. 
    subvols: integer or list of integers
       Number of subvolumes to be considered
    outpath : string
       Path to output, default is output/
    outnom : string
       Root name for output plot
    props : array of strings
       Dataset names to be plotted
    prop_labels : array of strings
       Dataset names for plot labels
    xmin : float
       Minimum value for the histogram
    xmax : float
       Maximum value for the histogram
    dx : float
       Histogram bin size
    metadata : dictionary
        Cosmology and other metadata information    
    verbose : boolean
       If True print out messages.
    '''
    # Get metadata
    vol_eff = metadata['vol_eff']
    redshift = metadata['redshift']    
    photmod_sfr = metadata['photmod_sfr']
    AGN = metadata['AGN']
    if AGN:
        photmod_agn = metadata['photmod_agn']
    att = metadata['att']
    if att:
        attmod = metadata['attmod']

    # Initialise histogram bins for luminosity functions
    xbins = np.arange(xmin, xmax, dx)
    xhist = xbins + dx*0.5

    # Initialise LF arrays
    nprops = len(props)
    yhist = np.zeros((nprops, len(xhist)))

    # Read data from each subvolume
    for ivol in subvols:
        # Read information from file and get histogram
        filenom = os.path.join(root+str(ivol),endf)
        f = h5py.File(filenom, 'r')
        for iprop, prop_name in enumerate(props):
            yall = io.read_and_add(f,prop_name,verbose=verbose)
            ind = np.where(yall > 0) 
            if np.shape(ind)[1] > 0:
                yy = np.log10(yall[ind])
                H, dum = np.histogram(yy,bins=np.append(xbins,xmax))
                yhist[iprop, :] += H
        f.close()

    # Normalize by bin size and volume
    yhist = yhist/dx/vol_eff

    # Plot settings
    plt.figure(figsize=(18, 18.))
    xtit = r'$\log_{10}$'+'(Property)' 
    plt.xlabel(xtit)
    ytit = r'$\log_{10}(\Phi/\mathrm{Mpc}^{-3}\,\mathrm{dex}^{-1})$'
    plt.ylabel(ytit)

    # Plot each property
    for iprop, prop_name in enumerate(props):
        leg = prop_labels[iprop]
        ilf = yhist[iprop, :]
        ind = np.where(ilf > 0)
        if len(ind[0]) > 0:
            x = xhist[ind]
            y = ilf[ind]
            indy = np.where(y > 0)
            if len(indy[0]) > 0:
                logy = np.log10(y[indy])
                plt.plot(x[indy], logy,label=leg)
    plt.legend(loc='best',frameon=False)
    
    # Output
    outf = outnom+'f'
    nom = io.get_plotfile(root,endf,outf)
    plt.savefig(nom)
    if verbose:
         print(f'* Function plot: {nom}')    
    return nom



def plot_line_lfs(root, endf, subvols=[0], outpath=None,
             metadata=None,verbose=True):
    '''
    Make line luminosity function plots
    
    Parameters
    ----------
    root : string
       Path to input files. 
    endf : string
       Ending of input files. 
    subvols: integer or list of integers
        Number of subvolumes to be considered
    outpath : string
        Path to output, default is output/
    metadata : dictionary
        Cosmology and other metadata information
    verbose : boolean
       If True print out messages.
    '''

    # Get metadata
    vol_eff = metadata['vol_eff']
    redshift = metadata['redshift']    
    photmod_sfr = metadata['photmod_sfr']
    AGN = metadata['AGN']
    if AGN:
        photmod_agn = metadata['photmod_agn']
    att = metadata['att']
    if att:
        attmod = metadata['attmod']

    # Read limits for properties and photoionisation models
    minU, maxU = get_limits(propname='logUs', photmod=photmod_sfr)
    minZ, maxZ = get_limits(propname='Z', photmod=photmod_sfr)

    # Input lines
    line_names = ['Halpha', 'Hbeta', 'NII6584', 'OII3727',
                  'OIII5007', 'SII6731', 'SII6717']

    # Define emission lines to plot and initialise LF arrays
    line_labels = [r'H$_{\alpha}$', r'H$_{\beta}$',
                   r'[OII]$\lambda\lambda 3727$', 
                   r'[OIII]$\lambda 5007$', r'[NII]$\lambda 6584$', 
                   r'[SII]$\lambda\lambda 6724$']

    # Initialise histogram bins for luminosity functions
    lmin = 38.0
    lmax = 46.0
    dl = 0.1
    lbins = np.arange(lmin, lmax, dl)
    lhist = lbins + dl * 0.5

    # Initialise LF arrays
    nlines = len(line_labels)
    lf = np.zeros((nlines, len(lhist)))
    lf_att = np.zeros((nlines, len(lhist)))

    # Read data from each subvolume
    for ivol in subvols:
        filenom = os.path.join(root+str(ivol),endf)
        f = h5py.File(filenom, 'r')

        # Read SF information from file
        lu_sfr = f['sfr_data/lu_sfr'][:,0]
        lz_sfr = f['sfr_data/lz_sfr'][:,0]

        ldims = f['sfr_data/Halpha_sfr'][:].ndim
        if ldims > 1:
            sfr_data = {line: np.sum(f[f'sfr_data/{line}_sfr'], axis=0)
                        for line in line_names}
        else:
            sfr_data = {line: f[f'sfr_data/{line}_sfr'][:]
                        for line in line_names}

        if att:
            # Initialize 
            ngal = sfr_data[line_names[0]].shape[0]
            sfr_data_att = {line: np.full(ngal, c.notnum) for line in line_names}

            for line in line_names: # Fill in available data
                key = f'sfr_data/{line}_sfr_att'
                if key in f:
                    ldims = f[key].ndim
                    if ldims > 1:
                        sfr_data_att[line] = np.sum(f[key], axis=0)
                    else:
                        sfr_data_att[line] = f[key][:]
        if AGN:
            # Read AGN information if it exists
            agn_data = {line: f[f'agn_data/{line}_agn'][:]
                        for line in line_names}
            if att:
                ngal = agn_data[line_names[0]].shape[0]
                agn_data_att = {line: np.full(ngal, c.notnum) for line in line_names}
                for line in line_names: # Fill in available data
                    key = f'agn_data/{line}_agn_att'
                    if key in f:
                        agn_data_att[line] = f[key][:]
        f.close()

        # Combine luminosities if adequate
        combined = {}
        combined_att = {}
        line_mapping = {
            'Ha': ['Halpha'],
            'Hb': ['Hbeta'],
            'NII': ['NII6584'],
            'OII': ['OII3727'],
            'OIII': ['OIII5007'],
            'SII': ['SII6731', 'SII6717']
        }

        for out_name, in_lines in line_mapping.items():
            arrays = [sfr_data[line] for line in in_lines]
            if AGN:
                arrays.extend([agn_data[line] for line in in_lines])
            combined[out_name] = st.safe_sum_arrays(arrays)
    
            if att:
                arrays_att = [sfr_data_att[line] for line in in_lines]
                if AGN:
                    arrays_att.extend([agn_data_att[line] for line in in_lines])
                combined_att[out_name] = st.safe_sum_arrays(arrays_att)

        # Calculate histograms for each line
        for iline, line in enumerate(line_mapping.keys()):
            # Intrinsic luminosity function
            lums = combined[line]
            ind = np.where(lums > 0) ###here more cuts like in bpt?
            if np.shape(ind)[1] > 0:
                ll = np.log10(lums[ind])
                H, dum = np.histogram(ll,bins=np.append(lbins,lmax))
                lf[iline, :] += H

            if att:
                # Dust attenuated luminosity function
                lums = combined_att[line]
                ind = np.where(lums > 0) ###here more cuts like in bpt?
                if np.shape(ind)[1] > 0:
                    ll = np.log10(lums[ind])
                    H, dum = np.histogram(ll,bins=np.append(lbins,lmax))
                    lf_att[iline, :] += H

    # Normalize by bin size and volume
    lf = lf/dl/vol_eff
    if att:
        lf_att = lf_att/dl/vol_eff

    # Plot settings
    fig, axes = plt.subplots(2, 3, figsize=(30,21))
    axes = axes.flatten()
    ytit = r'$\log_{10}(\Phi/\mathrm{Mpc}^{-3}\,\mathrm{dex}^{-1})$'
    xmin = 39.0
    xmax = 44.0
    ymin = -5.5
    ymax = -1.0

    # Plot each line
    for iline in range(nlines):
        ax = axes[iline]
        xtit = r'$\log_{10}$(L' + line_labels[iline] +\
            r'$/\mathrm{erg\,s^{-1}})$' 
        # Plot intrinsic LF (dotted line)
        ilf = lf[iline, :]
        ind = np.where(ilf > 0)
        if len(ind[0]) > 0:
            x = lhist[ind]
            y = ilf[ind]
            indy = np.where(y > 0)
            if len(indy[0]) > 0:
                logy = np.log10(y[indy])
                ax.plot(x[indy], logy, 'r:',
                        label=f'Intrinsic (z={redshift:.1f})')

        if att:
            # Plot dust-attenuated LF (solid line)
            ilf = lf_att[iline, :]
            ind = np.where(ilf > 0)
            if len(ind[0]) > 0:
                x = lhist[ind]
                y = ilf[ind]
                indy = np.where(y > 0)
                if len(indy[0]) > 0:
                    logy = np.log10(y[indy])
                    ax.plot(x[indy], logy, 'b-',
                            label='Dust-attenuated')
            
        # Set axis properties
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        ax.minorticks_on()
        ax.set_xlabel(xtit); ax.set_ylabel(ytit)
        if (iline==0) and len(ind[0]) > 0:
            ax.legend(loc='best',frameon=False)
    
    plt.tight_layout()
    
    # Output
    nom = io.get_plotfile(root,endf,'lf')
    plt.savefig(nom)
    if verbose:
         print(f'* LFs plots: {nom}')
    
    return nom


def plot_ncumu_flux(root, endf, subvols=[0], outpath=None,
                    metadata=None,verbose=True):
    '''
    Make plots with the cumulative numbers as a function of flux
    
    Parameters
    ----------
    root : string
       Path to input files. 
    endf : string
       Ending of input files. 
    subvols: integer or list of integers
        Number of subvolumes to be considered
    outpath : string
        Path to output, default is output/
    metadata : dictionary
        Cosmology and other metadata information
    verbose : boolean
       If True print out messages.
    '''

    # Get metadata
    vol_eff = metadata['vol_eff']
    redshift = metadata['redshift']    
    photmod_sfr = metadata['photmod_sfr']
    AGN = metadata['AGN']
    if AGN:
        photmod_agn = metadata['photmod_agn']
    att = metadata['att']
    if att:
        attmod = metadata['attmod']

    # Read limits for properties and photoionisation models
    minU, maxU = get_limits(propname='logUs', photmod=photmod_sfr)
    minZ, maxZ = get_limits(propname='Z', photmod=photmod_sfr)

    # Define emission lines to plot and initialise arrays
    line_labels = [r'H$_{\alpha}$',r'H$_{\alpha}$+N[II]','O[III]',r'O[III]+H$_{\beta}$']

    # Initialise histogram bins for luminosity functions
    fmin = -18
    fmax = -14
    df = 0.2
    fbins = np.arange(fmin, fmax, df)

    # Initialise arrays
    nlines = len(line_labels)
    ncum = np.zeros((nlines, len(fbins)))
    ncum_att = np.zeros((nlines, len(fbins)))

    # Read data from each subvolume
    for ivol in subvols:
        filenom = os.path.join(root+str(ivol),endf)
        f = h5py.File(filenom, 'r')

        # Read SF information from file
        lu_sfr = f['sfr_data/lu_sfr'][:,0]
        lz_sfr = f['sfr_data/lz_sfr'][:,0]
        Ha_sfr = np.sum(f['sfr_data/Halpha_sfr_flux'],axis=0)
        Hb_sfr = np.sum(f['sfr_data/Hbeta_sfr_flux'],axis=0)
        NII_sfr = np.sum(f['sfr_data/NII6584_sfr_flux'],axis=0)
        OIII_sfr = np.sum(f['sfr_data/OIII5007_sfr_flux'],axis=0)

        if att:
            Ha_sfr_att = np.sum(f['sfr_data/Halpha_sfr_att_flux'],axis=0)
            Hb_sfr_att = np.sum(f['sfr_data/Hbeta_sfr_att_flux'],axis=0)
            NII_sfr_att = np.sum(f['sfr_data/NII6584_sfr_att_flux'],axis=0)
            OIII_sfr_att = np.sum(f['sfr_data/OIII5007_sfr_att_flux'],axis=0)

        if AGN:
            # Read AGN information if it exists
            Ha_agn = f['agn_data/Halpha_agn_flux'][:]
            Hb_agn = f['agn_data/Hbeta_agn_flux'][:]
            OIII_agn = f['agn_data/OIII5007_agn_flux'][:]
            NII_agn = f['agn_data/NII6584_agn_flux'][:]

            if att:
                Ha_agn_att = f['agn_data/Halpha_agn_att_flux'][:]
                Hb_agn_att = f['agn_data/Hbeta_agn_att_flux'][:]
                OIII_agn_att = f['agn_data/OIII5007_agn_att_flux'][:]
                NII_agn_att = f['agn_data/NII6584_agn_att_flux'][:]
        f.close()

        if AGN:
            # Combine luminosities (SFR + AGN if available)
            Ha = Ha_sfr + Ha_agn
            HaN2 = Ha + NII_sfr + NII_agn
            O3 = OIII_sfr + OIII_agn
            O3Hb = O3 + Hb_sfr + Hb_agn 
            if att:
                Ha_att = Ha_sfr_att + Ha_agn_att
                HaN2_att = Ha_att + NII_sfr_att + NII_agn_att
                O3_att = OIII_sfr_att + OIII_agn_att
                O3Hb_att = O3_att + Hb_sfr_att + Hb_agn_att 
        else:
            Ha = Ha_sfr
            HaN2 = Ha + NII_sfr
            O3 = OIII_sfr
            O3Hb = O3 + Hb_sfr
            if att:
                Ha_att = Ha_sfr_att
                HaN2_att = Ha_att + NII_sfr_att
                O3_att = OIII_sfr_att
                O3Hb_att = O3_att + Hb_sfr_att

        flux_int = [Ha, HaN2, O3, O3Hb]
        if att:
            flux_att = [Ha_att, HaN2_att, O3_att, O3Hb_att]
        
        # Calculate the cumulative numbers for each line
        for iline in range(nlines):
            # Intrinsic flux
            flux = flux_int[iline]
            ind = np.where(flux > 0) ###here more cuts like in bpt?
            if np.shape(ind)[1] > 0:
                ff = np.log10(flux[ind])
                H = n_gt_x(fbins,ff)
                ncum[iline,:] = ncum[iline,:] + H

            if att: # Attenuated flux
                flux = flux_att[iline]
                ind = np.where(flux > 0) ###here more cuts like in bpt?
                if np.shape(ind)[1] > 0:
                    ff = np.log10(flux[ind])
                    H = n_gt_x(fbins,ff)
                    ncum_att[iline,:] = ncum_att[iline,:] + H
                    
    # Get number per volume
    ncum = ncum/vol_eff
    if att:
        ncum_att = ncum_att/vol_eff

    # Plot settings
    nfigs = 2
    fig, axes = plt.subplots(1, 2, figsize=(30,21))
    axes = axes.flatten()
    ytit = r'$\log_{10}(n_{\rm gal}(>F_{\rm lim})/\mathrm{Mpc}^{-3}\,\mathrm{dex}^{-1})$'
    xmin = fmin
    xmax = fmax
    ymin = -2
    ymax = -1.0

    line = -2
    for ifig in range(nfigs):
        ax = axes[ifig]
        xtit = r'$\log_{10}(F_{\rm lim}/\mathrm{erg\,s^{-1}\,cm^{-2}})$'
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        ax.minorticks_on()
        ax.set_xlabel(xtit); ax.set_ylabel(ytit)
        line += 2
        for iline in [line,line+1]:
            color = plt.cm.tab10(iline % 10)
            # Plot intrinsic n
            yy = ncum[iline, :]
            ind = np.where(yy > 0)
            if len(ind[0]) > 0:
                x = fbins[ind]
                y = np.log10(yy[ind])
                ll = line_labels[iline]+'(int.)'
                ax.plot(x, y, '-',color=color,label=ll)

            if att: # Dust-attenuated 
                yy = ncum_att[iline, :]
                ind = np.where(yy > 0)
                if len(ind[0]) > 0:
                    x = fbins[ind]
                    y = np.log10(yy[ind])
                    ll = line_labels[iline]+'(att.)'
                    ax.plot(x, y,'--',color=color,label=ll)
            
        # Legend
        if len(ind[0]) > 0:
            ax.legend(loc='best',frameon=False)
    fig.suptitle(f'z = {redshift}')
    plt.tight_layout()
    
    # Output
    nom = io.get_plotfile(root,endf,'flux_ncumu')
    plt.savefig(nom)
    if verbose:
         print(f'* Cumulative numbers vs flux plots: {nom}')
    
    return nom


def make_gridplots(xid_sfr=0.3,co_sfr=1,imf_cut_sfr=100,
                   xid_NLR=0.5,alpha_NLR=-1.7,verbose=True):
    '''
    Make plots for photoionisation tables
    
    Parameters
    ----------
    verbose : boolean
       If True print out messages.
    '''

    # Plot photoionisation grids on BPT diagrams
    grids_sfr = plot_model_bpt_grids(photmod='gutkin16',xid=xid_sfr,
                                     co=co_sfr,imf_cut=imf_cut_sfr,
                                     verbose=verbose)
    grids_agn = plot_model_bpt_grids(photmod='feltre16',
                                     xid=xid_NLR,alpha=alpha_NLR,
                                     verbose=verbose)
    
    return


def make_testplots(snap,ending,outpath=None,
                   subvols=[0],gridplots=False,verbose=True):
    '''
    Make test plots
    
    Parameters
    ----------
    snap: integer
        Simulation snapshot number
    out_ending : string
       End name of input files
    outpath : string
       Path to input files
    subvols: list of integers
        List of subvolumes to be considered
    outpath : string
        Path to output, default is output/ 
    verbose : boolean
       If True print out messages.
    '''
    root, endf = io.get_outroot(snap,ending,outpath=outpath,
                                verbose=verbose)
    plots_dir = os.path.join(os.path.dirname(root), 'plots')

    # Get metadata from line files
    ivol0 = str(subvols[0])
    filenom = os.path.join(root+ivol0,endf)
    metadata = io.get_metadata(filenom, verbose=verbose)

    # Set cosmology only once
    set_cosmology(omega0 = metadata['omega0'],
                  omegab = metadata['omegab'],
                  lambda0 = metadata['lambda0'],
                  h0 = metadata['h0'])  
    

    ### Photoionisation plots
    #if gridplots:
    #    make_gridplots() ###here work in progress

    #### Characterisation of global properties
    #if (metadata['AGN']):
    #    lbol_lf = plot_lf(root,endf,subvols=subvols,outpath=outpath,
    #                      metadata=metadata,verbose=verbose)
        
    ### Characterisation of properties of ionising regions
    # U vs Z
    #uzn = plot_uzn(root,endf,subvols=subvols,verbose=verbose) 
        
    # Line plots
    # Make NII and SII bpt plots
    bpt = plot_bpts(root,endf,subvols=subvols,outpath=outpath,
                    metadata=metadata,verbose=verbose)
    
    # Make line LFs
    lfs = plot_line_lfs(root,endf,subvols=subvols,outpath=outpath,
                   metadata=metadata,verbose=verbose)
    
    # Cumulative numbers with flux limits (if possible) 
    if (metadata['flux'] and metadata['redshift']>0):
        ncumu_flux = plot_ncumu_flux(root,endf,subvols=subvols,
                                     outpath=outpath,metadata=metadata,
                                     verbose=verbose)
    else:
        if verbose:
            if (metadata['flux']):
                print(f'WARNING: Skipping cumulative flux plot at z=0.')
            else:
                print(f'WARNING: No flux data found in {filenom}.')

    print(f'SUCCESS: plots in {plots_dir}')
    return
