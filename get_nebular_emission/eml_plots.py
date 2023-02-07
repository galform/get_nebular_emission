import numpy as np
import matplotlib
import os.path
import h5py
# matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import get_nebular_emission.eml_style as style
from get_nebular_emission.stats import perc_2arrays
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_io import get_nheader, check_file
from get_nebular_emission.eml_photio import get_lines_Gutkin, get_limits
from numpy import random
from scipy import stats
import sys

plt.style.use(style.style1)

def test_sfrf(inputdata, outplot, obsSFR=None, obsGSM=None, colsSFR=[0,1,2,3],
              colsGSM=[0,1,2,3], labelObs=None, specific=False, h0=const.h, volume=const.vol_pm, verbose=False):

    '''
    
    Given log10(Mstar) and log10(sSFR) get the plots to compare log10(SFR) vs log10(Mstar).
    Get the GSMF and the SFRF plots. 
    Given the observations, compare the plots with the observations too.
 
    Parameters
    ----------
 
    obsSFR : string
      - Name of the input file for the SFR data observed.
      - In text files (*.dat, *txt, *.cat), columns separated by ' '.
      - In csv files (*.csv), columns separated by ','.
      - Expected histogram mode:
       - A column with the low value of the bin,
       - A column with the high value of the bin,
       - A column with the frequency in the bin,
       - A column with the error. 
 
    obsGSM : string
      - Name of the input file for the GSM data observed.
      - In text files (*.dat, *txt, *.cat), columns separated by ' '.
      - In csv files (*.csv), columns separated by ','.
      - Expected histogram mode:
       - A column with the low value of the bin,
       - A column with the high value of the bin,
       - A column with the frequency in the bin,
       - A column with the error.
 
    colsSFR : list
      - Columns with the data required to do the observational histogram of the SFR.
      - Expected: [ind_column1, ind_column2, ind_column3, ind_column4]
       - column1 is the column with the low values of the bins, in Msun/yr,
       - column2 with the high values of the bins, in Msun/yr,
       - column3 with the frequency, in Mpc^-3 dex^-1
       - column4 with the error, in Mpc^-3 dex^-1
       
    colsGSM : list
      - Columns with the data required to do the observational histogram of the GSM.
      - Expected: [ind_column1, ind_column2, ind_column3, ind_column4]
       - column1 is the column with the low values of the bins, in h^-2Msun,
       - column2 with the high values of the bins, in h^-2Msun,
       - column3 with the frequency, in h^-3 Mpc^-3,
       - column4 with the error, in h^-3 Mpc^-3.
 
    labelObs : list of strings
      - For the legend, add the name to cite the observational data source.
      - ['GSM observed', 'SFR observed']
 
    outplot : string
      - Name of the output file.
      - Image-type files (*.pdf, *.jpg, ...)
      
    specific : boolean
      If True it makes the plots with the sSFR. Otherwise, it makes the plots with the SFR.
 
    h0 : float
      If not None: value of h, H0=100h km/s/Mpc.
      
    volume : float
      - Carlton model default value = 542.16^3 Mpc^3/h^3.
      - table 1: https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4922B/abstract
      - If not 542.16**3. : valume of the simulation volume in Mpc^3/h^3
    verbose : boolean
      If True print out messages
 
    Notes
    -------
    It makes plot(log10(SFR),log10(Mstar)), plot GSMF and plot SFRF,
    all three in one grid and saves it in the outplot path.
    '''



    # Define a class that forces representation of float to look a certain way
    # This remove trailing zero so '1.0' becomes '1'
    class nf(float):
        def __repr__(self):
            str = '%.1f' % (self.__float__(),)
            if str[-1] == '0':
                return '%.0f' % self.__float__()
            else:
                return '%.1f' % self.__float__()
    # -----------------------------------------------------


    # Correct the units of the simulation volume to Mpc^3:
    if h0:
        volume=volume/(h0**3)

    #Prepare the plot
    lsty = ['-',(0,(2,3))] # Line form

    nds = np.array([-2., -3., -4.]) # Contours values
    al = np.sort(nds)

    SFR = ['avSFR']
    labels = ['SFR']

    cm = plt.get_cmap('tab10')  # Colour map to draw colours from
    color = []
    for ii in range(0, 10):
        col = cm(ii)
        color.append(col)  # col change for each iteration


    # Initialize GSMF (Galaxy Cosmological Mass Function)
    mmin = 8 #10.3 # mass resolution 2.12 * 10**9 h0 M_sun (Baugh 2019)
    mmax = 15 
    dm = 0.1
    mbins = np.arange(mmin, mmax, dm)
    mhist = mbins + dm * 0.5
    gsmf = np.zeros((len(mhist)))

    # Initialize SSFRF
    smin = -5
    smax = 2.5
    ds = 0.5
    sbins = np.arange(smin, smax, ds)
    shist = sbins + ds * 0.5
    ssfrf = np.zeros((len(shist)))

    # Initialize SFR vs M function
    lenm = len(mhist)
    lens = len(shist)
    smf = np.zeros((lens,lenm))

    # Plots limits and style
    fig = plt.figure(figsize=(8.5, 9.))
    gs = gridspec.GridSpec(3, 3)
    gs.update(wspace=0., hspace=0.)
    ax = plt.subplot(gs[1:, :-1])

    # Fig. sSFR vs M
    xtit = "log$_{10}(\\rm M_{*}/M_{\odot})$"
    if specific:
        ytit = "log$_{10}(\\rm sSFR Gyr^{-1})$"
    else:
        ytit = "log$_{10}(\\rm SFR/M_{\odot} yr^{-1})$"
    xmin = mmin; xmax = 12; ymin = smin;  ymax = smax
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xtit); ax.set_ylabel(ytit)

    # GSMF
    axm = plt.subplot(gs[0, :-1],sharex=ax)
    ytit="log$_{10}(\Phi(M_*))$" ; axm.set_ylabel(ytit)
    axm.set_autoscale_on(False) ;  axm.minorticks_on()
    axm.set_ylim(-8,-1)
    plt.setp(axm.get_xticklabels(), visible=False)

    # SSFRF
    axs = plt.subplot(gs[1:, 2], sharey=ax)
    xtit = "log$_{10}(\Phi(SFR))$"; axs.set_xlabel(xtit)
    axs.set_autoscale_on(False); axs.minorticks_on()
    axs.set_xlim(-5.5, 0.0)
    start, end = axs.get_xlim()
    axs.xaxis.set_ticks(np.arange(-4., end, 1.))
    plt.setp(axs.get_yticklabels(), visible=False)

    # Data Observations

    # SFR observed

    if obsSFR:
        ih = get_nheader(obsSFR)
    
        dataSFR = [0]*len(colsSFR)
    
        for ii, col in enumerate(colsSFR):
            #print(ii,col,colsSFR[ii])
            data = np.loadtxt(obsSFR,skiprows=ih, usecols=col, unpack=True)
            dataSFR[ii] = np.array(data)
    
        dex = dataSFR[1]-dataSFR[0]
        histSFR = dataSFR[1]-0.5*dex
        errorSFR = dataSFR[3]

    # GSM observed
    if obsGSM:
        ih = get_nheader(obsGSM)
    
        dataGSM = [0]*len(colsGSM)
    
        for ii, col in enumerate(colsGSM):
            data = np.loadtxt(obsGSM,skiprows=ih, usecols=col, unpack=True)
            dataGSM[ii] = np.array(data)
    
        dex = dataGSM[1] - dataGSM[0]
    
        # Change the units from h^-2 Msun to Msun.
        histGSM = dataGSM[1] - 2*np.log10(h0) - 0.5*dex
    
        # Change the units from h^3 Mpc^-3 to Mpc^-3
        freqGSM = np.log10((dataGSM[2])) + 3 * np.log10(h0)
        
        lowGSM = np.log10(dataGSM[2]-dataGSM[3]) + 3 * np.log10(h0)
        
        lowGSM = abs(lowGSM - freqGSM)

    for ii in range(len(inputdata)):
        with h5py.File(inputdata[ii],'r') as file:
            data = file['data']          
            lms = np.log10(10**data['lms'][:,0]+10**data['lms'][:,1])
            if specific:
                lsfr = np.log10(10**data['lssfr'][:,0]+10**data['lssfr'][:,1]) + 9
            else: 
                lsfr = np.log10(10**data['lssfr'][:,0]+10**data['lssfr'][:,1]) + lms
            # lms = lms + np.log10(h0)     
            del data


        # Make the histograms

        H, bins_edges = np.histogram(lms, bins=np.append(mbins, mmax))
        gsmf = H / volume / dm  # In Mpc^3/h^3

        H, bins_edges = np.histogram(lsfr, bins=np.append(sbins, smax))
        sfrf = H / volume / ds

        H, xedges, yedges = np.histogram2d(lsfr, lms,
                                           bins=([np.append(sbins, smax),
                                                  np.append(mbins, mmax)]))
        smf = H / volume / dm / ds


        # Plot SMF vs SFR

        matplotlib.rcParams['contour.negative_linestyle'] = lsty[ii]
        zz = np.zeros(shape=(len(shist), len(mhist))); zz.fill(const.notnum)
        ind = np.where(smf > 0.)
        zz[ind] = np.log10(smf[ind])
        
        # print(zz[ind])

        ind = np.where(zz > const.notnum)

        if (np.shape(ind)[1] > 1):

            # Contours
            xx, yy = np.meshgrid(mbins, sbins)
            # Here: How to find the levels of the data?
            cs = ax.contour(xx, yy, zz, levels=al, colors=color[ii])
            ax.clabel(cs, inline=1, fontsize=10)
            for i in range(len(labels)):
                cs.collections[i].set_label(labels[i])

        # Plot GSMF
        py = gsmf; ind = np.where(py > 0.)
        x = mhist[ind]; y = np.log10(py[ind])
        ind = np.where(y < 0.)

        axm.plot(x[ind], y[ind], color=color[ii],
                 linestyle=lsty[ii], label=labels[ii])
    
            # Plot observations GSMF
        if obsGSM:
            axm.errorbar(histGSM, freqGSM, yerr=lowGSM, marker='o', color=color[ii + 1],
                             label=''+ labelObs[0] +'')
                
            leg2 = axm.legend(bbox_to_anchor=(0.135, -0.34, 1.5, 1.4), fontsize='small',
                              handlelength=1.2, handletextpad=0.4)
            leg2.get_texts()
            leg2.draw_frame(False)


        # Plot SFRF
        px = sfrf; ind = np.where(px > 0.)
        y = shist[ind]; x = np.log10(px[ind])
        ind = np.where(x < 0.)
        axs.plot(x[ind], y[ind], color=color[ii],
                 linestyle=lsty[ii])
            
            # Plot observations SFRF
        if obsSFR:
            axs.errorbar(dataSFR[2], histSFR, xerr=errorSFR, marker='o', color=color[ii + 2],
                          label=''+ labelObs[1] +'')

            leg = axs.legend(bbox_to_anchor=(1.5, 1.4), fontsize='small',
                              handlelength=1.2, handletextpad=0.4)
            leg.get_texts()
            leg.draw_frame(False)

    plotf = outplot

    # Save figures
    print('Plot: {}'.format(plotf))
    fig.savefig(plotf)

def test_medians(infile, outplot, verbose=False):
    '''
    Given U and ne calculated from the Mstar and the SFR in eml_une.
    get the plot of the medians of these quantities in masses bins.
    Compare the plots between the data calculated from the average SFR and
    from the SFR calculated from the LC photons.

    Parameters
    ----------
    outplot : string
     Path to the folder where save the plot.
    verbose : boolean
     If True print out messages


    Notes
    -------
    Medians U and ne for each bin of mass in two differents plots.

    '''

    # Prepare the plots
    
    SFR = ['avSFR']
    U_ne = ['u', 'ne']
    cm = plt.get_cmap('tab10')  # Colour map to draw colours from
    color = []

    # Prepare the bins
    mmin = 8.5
    mmax = 11.5
    dm = 0.2
    mbins = np.arange(mmin, (mmax + dm), dm)
    mbinsH = np.arange(mmin, mmax, dm)
    mhist = mbinsH + dm * 0.5

    for iu, une in enumerate(U_ne):

        # Prepare the figures
        plt.figure()
        plt.xlabel(r'$\log M_*/h^{-1} (M_\odot)$',size=15)
        plotf = outplot + '/test_medians_'+ une+'.pdf'
        col = cm(iu)
        color.append(col)  # col change for each iteration
        
        if iu==0:
            plt.ylim((-5,-2))
            plt.ylabel(r'$\log U$',size=15)
        else:
            plt.ylim((-0.5,2.5))
            plt.ylabel(r'$\log n_e (cm^{-3})$',size=15)

        for ii,sfr in enumerate(SFR):
            
            with h5py.File(infile,'r') as file:
                print(infile)
                f = file['data']
                lu = f['lu'][:,0]
                lne = f['lne'][:,0]   
                lms = f['lms'][:,0]  + np.log10(const.h)
            
            data = np.append([lu], [lne], axis=0)
            print(data.shape)
            
            cut = np.where((lms>7)&(lu!=const.notnum)&(lne!=const.notnum))
            
            # MEDIANS:
            median = perc_2arrays(mbins, lms[cut], data[iu][cut], 0.5) #data[iu]

            ind = np.where(median>const.notnum)
            
            median = median[ind]

            # QUARTILES:
            up_qu = perc_2arrays(mbins, lms[cut], data[iu][cut], 0.75)[ind]
            #qup[ii] = up_qu
            low_qu = perc_2arrays(mbins, lms[cut], data[iu][cut], 0.25)[ind]
            #qlow[ii] = low_qu
            
            qu = np.append([median-low_qu],[up_qu-median],axis=0)

            #plt.plot(mhist[ind],median[ind],'o', color=col, label='Calculated from the ' + SFR[ii] + '')
            #plt.plot(mhist[ind],up_qu[ind],'o', color='r')
            #plt.plot(mhist[ind],low_qu[ind],'o', color='r')
            plt.errorbar(mhist[ind],median,marker='o',yerr=qu,elinewidth=0.5, color=col)#, label='Calculated from the ' + SFR[ii] + '')
        # plt.legend()
        plt.savefig(plotf)
        plt.close()



def test_bpt(infile, outplot, photmod='gutkin16', plot_phot=False, create_file=False, file_folder='output_data', verbose=False):
    '''
    Run a test of the interpolations done in eml_photio.
    Two plots, one to verify the U interpolation and the other one to verify the Z interpolation
    
    Parameters
    ----------
    infile : string
     Name of the input file. 
    outplot : string
     Path to the folder plot.
    photmod : string
      Photoionisation model to be used for look up tables.
    plot_phot : boolean
     If True it plots points from the photoionization tables.
    create_file : boolean
     If True it creates textfiles to read the photoionization tables.
    file_folder : string
     Folder where the textfiles to read the tables will be/are stored.
    verbose : boolean
     If True print out messages.

    Notes
    -------
    Plot of several BPT diagrams.
    '''
    
    check_file(infile, verbose=True)
    f = h5py.File(infile, 'r')
    data = f['data']

    lu_disk = data['lu'][:,0]
    lne_disk = data['lne'][:,0]
    loh12_disk = data['lz'][:,0]
    
    minU, maxU = get_limits(propname='U', photmod=photmod)
    minnH, maxnH = get_limits(propname='nH', photmod=photmod)
    minZ, maxZ = get_limits(propname='Z', photmod=photmod)
    
    ind = np.where((lu_disk!=minU)&(lu_disk!=maxU)&(loh12_disk!=np.log10(minZ))&(loh12_disk!=np.log10(maxZ))&
                      (lne_disk!=np.log10(minnH))&(lne_disk!=np.log10(maxnH))) #We ignore the galaxies in the limits.
    
    Hbeta = np.sum(data['Hbeta_att'],axis=0)[ind]
    OIII5007 = np.sum(data['OIII5007_att'],axis=0)[ind]
    NII6548 = np.sum(data['NII6583_att'],axis=0)[ind]
    Halpha = np.sum(data['Halpha_att'],axis=0)[ind]
    SII6717_6731 = np.sum(data['SII6717_att'],axis=0)[ind] + np.sum(data['SII6731_att'],axis=0)[ind]
    
    # Hbeta = Hbeta[::10000]
    # OIII5007 = OIII5007[::10000]
    # NII6548 = NII6548[::10000]
    # Halpha = Halpha[::10000]
    # SII6717_6731 = SII6717_6731[::10000]
    
    ind2 = np.where((Hbeta!=0)&(OIII5007!=0)&(NII6548!=0)&(Halpha!=0)&(SII6717_6731!=0))
    
    print(len(ind[0]),len(ind2[0]))
    
    Hbeta = Hbeta[ind2]
    OIII5007 = OIII5007[ind2]
    NII6548 = NII6548[ind2]
    Halpha = Halpha[ind2]
    SII6717_6731 = SII6717_6731[ind2]
    
    bpt_x = ['log$_{10}$([NII]$\\lambda$6584/H$\\alpha$)','log$_{10}$([SII]$\\lambda$(6717+6731)/H$\\alpha$)']
    my_x = [np.log10(NII6548 / Halpha),np.log10(SII6717_6731 / Halpha)]
    
    bpt_y = ['log$_{10}$([OIII]$\\lambda$5007/H$\\beta$)','log$_{10}$([OIII]$\\lambda$5007/H$\\beta$)']
    my_y = [np.log10(OIII5007 / Hbeta),np.log10(OIII5007 / Hbeta)]
    
    if not plot_phot:
        for i in range(2):
            plt.figure(figsize=(15,15))
            
            xmin=-2.5
            xmax=1
            ymin=-2
            ymax=2
            
            # X1, Y1 = np.mgrid[xmin:xmax:68j, ymin:ymax:68j]
            # positions = np.vstack([X1.ravel(), Y1.ravel()])
            # values = np.vstack([my_x[i], my_y[i]])
            # kernel = stats.gaussian_kde(values,0.75)
            # BPT = np.reshape(kernel(positions).T, X1.shape)
            # plt.imshow(BPT, cmap=plt.cm.gist_earth_r,extent=[xmin, xmax, ymin, ymax],aspect=(xmax-xmin)/(ymax-ymin))#,vmin=0,vmax=1)
            
            plt.plot(my_x[i], my_y[i], marker='o', markersize=1, linewidth=0, color='black')
            
            plt.xlabel(bpt_x[i],size=30)
            plt.ylabel(bpt_y[i],size=30)
            plt.xticks(fontsize=30)
            plt.yticks(fontsize=30)
            
            plt.xlim((xmin,xmax))
            plt.ylim((ymin,ymax))
            plt.grid()
            
            plotnom = outplot + '/BPTplot2_' + str(i) + '.png'
        
            plt.savefig(plotnom)
            plt.close()
    
    if plot_phot:       
        if photmod not in const.photmods:
            if verbose:
                print('STOP (eml_photio.test_bpt): Unrecognised model to get emission lines.')
                print('                Possible photmod= {}'.format(const.photmods))
            sys.exit()
        elif (photmod == 'gutkin16'):
            
            Z = ['0001', '0002', '0005', '001', '002', '004', '006', '008', '010', '014', '017', '020', '030', '040']
        
            zz = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.014, 0.017, 0.02, 0.03, 0.04]
        
            uu = [-1., -1.5, -2., -2.5, -3., -3.5, -4.]
        
            ne = ['100']  # ['10', '100', '1000','10000']
        
            cm = plt.get_cmap('tab20') # Colour map to draw colours from
        
            if create_file:
                for iz, zname in enumerate(Z):
                    infile = r"nebular_data/gutkin_tables/nebular_emission_Z" + zname + ".txt"
            
                    ih = get_nheader(infile)
            
                    datane = np.loadtxt(infile, skiprows=ih, usecols=(2), unpack=True)
                    datalu = np.loadtxt(infile, skiprows=ih, usecols=(0), unpack=True)
            
                    OIII5007_model = np.loadtxt(infile, skiprows=ih, usecols=(8), unpack=True)
                    Hb_model = np.loadtxt(infile, skiprows=ih, usecols=(6), unpack=True)
                    NII6548_model = np.loadtxt(infile, skiprows=ih, usecols=(9), unpack=True)
                    Ha_model = np.loadtxt(infile, skiprows=ih, usecols=(10), unpack=True)
                    SII6717_6731_model = np.loadtxt(infile, skiprows=ih, usecols=(12), unpack=True) + np.loadtxt(infile, skiprows=ih, usecols=(12), unpack=True)
            
                    for ii, nh in enumerate(ne):
                        outfile = r"output_data/Gutkinfile_n_" + nh + ".txt"
                        if iz==0 and os.path.exists(outfile):
                            os.remove(outfile)
            
                        header1 = 'Z, U, NII6584/Ha  OIII5007/Hb, SII(6717+6731)/Ha'
            
                        ind = np.where(datane == float(nh))
                        x = np.log10(NII6548_model[ind] / Ha_model[ind])
                        y = np.log10(OIII5007_model[ind] / Hb_model[ind])
                        p = np.log10(SII6717_6731_model[ind] / Ha_model[ind])
                        u = datalu[ind]
                        z = np.full(np.shape(u), zz[iz])
            
                        tofile = np.column_stack((z, u, x, y, p))
            
                        with open(outfile, 'a') as outf:
                            if iz == 0:
                                np.savetxt(outf, tofile, delimiter=' ', header=header1)
                            else:
                                np.savetxt(outf, tofile, delimiter=' ')
                            outf.closed
            else:
                for ii, nh in enumerate(ne):
                    outfile = r"output_data/Gutkinfile_n_" + nh + ".txt"
                    if not os.path.exists(outfile):
                        print('STOP (eml_photio.test_bpt): Textfiles for table reading dont exist.')
                        print('Create them with create_file = True.')
        
            cols = []
            for iz, lz in enumerate(zz):
                col = cm(iz)
                cols.append(col)
        
            for ii, nh in enumerate(ne):
                infile = r"output_data/Gutkinfile_n_" + nh + ".txt"
        
                ih = get_nheader(infile)
        
                z = np.loadtxt(infile, skiprows=ih, usecols=(0), unpack=True)
                u = np.loadtxt(infile, skiprows=ih, usecols=(1), unpack=True)
                x = np.loadtxt(infile, skiprows=ih, usecols=(2), unpack=True)
                y = np.loadtxt(infile, skiprows=ih, usecols=(3), unpack=True)
                p = np.loadtxt(infile, skiprows=ih, usecols=(4), unpack=True)
                
                comp_x = [x,p]
                comp_y = [y,y]
        
                # DIFERENTS COLORS FOR U:
                    
                for i in range(2):
                    plt.figure(figsize=(15,15))
        
                    for iu, lu in enumerate(uu):
                        ind2 = np.where(u == uu[iu])
                        plt.plot(comp_x[i][ind2], comp_y[i][ind2], marker='.', linewidth=0, color=cols[iu], label='U = ' + str(lu) + '')
            
                    labelsU = []
                    for elem in lu_disk: labelsU.append('U = {}'.format(np.round(elem,2)))
                    
                    plt.plot(my_x[i], my_y[i], marker='o', markersize=2, linewidth=0, color='black')
                    
                    plt.xlabel(bpt_x[i],size=30)
                    plt.ylabel(bpt_y[i],size=30)
                    plt.xticks(fontsize=30)
                    plt.yticks(fontsize=30)
                    plt.grid()
                    plt.legend()
                    
                    plotnom = outplot + '/BPTplot_U_' + str(i) + '.png'
                    
                    print('U', str(i))
                
                    plt.savefig(plotnom)
                    plt.close()
        
                # DIFFERENTS COLORS FOR Z:
                    
                for i in range(2):
                    plt.figure(figsize=(15,15))
        
                    for iz, lz in enumerate(zz):
                        ind2 = np.where(z == zz[iz])
                        plt.plot(comp_x[i][ind2], comp_y[i][ind2], marker='.', linewidth=0, color=cols[iz], label='Z = ' + str(lz) + '')
            
                    labelsZ = []
                    for elem in loh12_disk: labelsZ.append('Z = {:.4f}'.format(10 ** (elem)))
                    
                    plt.plot(my_x[i], my_y[i], marker='o', markersize=2, linewidth=0, color='black')
                    
                    plt.xlabel(bpt_x[i],size=30)
                    plt.ylabel(bpt_y[i],size=30)
                    plt.xticks(fontsize=30)
                    plt.yticks(fontsize=30)
                    plt.grid()
                    plt.legend()
                    
                    plotnom = outplot + '/BPTplot_Z_' + str(i) + '.png'
                    
                    print('Z', str(i))
                
                    plt.savefig(plotnom)
                    plt.close()