import numpy as np
import matplotlib
import os.path

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import get_nebular_emission.eml_style as style
from get_nebular_emission.stats import perc_2arrays
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_io import get_nheader

plt.style.use(style.style1)

def test_sfrf(obsSFR, obsGSM, colsSFR,colsGSM,labelObs, outplot, h0=None, volume=const.vol_pm, verbose=False):

    '''
       Given log10(Mstar) and log10(sSFR)
       get the plots to compare log10(SFR) vs log10(Mstar).
       Get the GSMF and the SFRF plots.

       Given the observations, compare the plots with the observations too.


       Parameters
       ----------

       obsSFR : string
         Name of the input file for the SFR data observed.
         Expected histogram mode:
         with a column with the low value of the bin,
         a column with the high value of the bin,
         a column with the frequency in the bin,
         and a column with the error.
         These columns must be specify in the colsSFR variable.

         In text files (*.dat, *txt, *.cat), columns separated by ' '.
         In csv files (*.csv), columns separated by ','.

       obsGSM : string
         Name of the input file for the GSM data observed.

         Expected histogram mode:
         with a column with the low value of the bin,
         a column with the high value of the bin,
         a column with the frequency in the bin,
         and a column with the error.
         These columns must be specify in the colsGSM variable.

         In text files (*.dat, *txt, *.cat), columns separated by ' '.
         In csv files (*.csv), columns separated by ','.

       colsSFR : list
         Columns with the data required to do the observational histogram of the SFR.
         Expected: [ind_column1, ind_column2, ind_column3, ind_column4]
         where: column1 is the column with the low values of the bins, in Msun/yr,
                column2 with the high values of the bins, in Msun/yr,
                column3 with the frequency, in Mpc^-3 dex^-1
                column4 with the error, in Mpc^-3 dex^-1
       colsGSM :
         Columns with the data required to do the observational histogram of the GSM.
         Expected: [ind_column1, ind_column2, ind_column3, ind_column4]
         where: column1 is the column with the low values of the bins, in h^-2Msun,
                column2 with the high values of the bins, in h^-2Msun,
                column3 with the frequency, in h^-3 Mpc^-3,
                column4 with the error, in h^-3 Mpc^-3.

       labelObs : list of strings
         For the legend, add the name to cite the observational data source.
         ['GSM observed', 'SFR observed']

       outplot : string
         Name of the output file.
         Image-type files (*.pdf, *.jpg, ...)

       h0 : float
         If not None: value of h, H0=100h km/s/Mpc.
       volume : float
         Carlton model default value = 542.16^3 Mpc^3/h^3.
         table 1: https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4922B/abstract
         If not 542.16**3. : value of the simulation volume in Mpc^3/h^3
       verbose : boolean
         Yes = print out messages

       Returns
       -------
       plot(log10(SFR),log10(Mstar)), plot GSMF and plot SFRF,
       all three in one grid.
       Save it in the outplot path.
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
        volume=volume*(h0**3)

    volume = volume/200 # There is only 1 sub-volume

    #Prepare the plot
    lsty = ['-',(0,(2,3))] # Line form

    nds = np.array([-2., -3., -4.]) # Contours values
    al = np.sort(nds)

    SFR = ['LC', 'avSFR']
    labels = ['average SFR', 'SFR from LC photons']

    cm = plt.get_cmap('tab10')  # Colour map to draw colours from
    color = []
    for ii in range(0, 10):
        col = cm(ii)
        color.append(col)  # col change for each iteration


    # Initialize GSMF (Galaxy Cosmological Mass Function)
    mmin = 10.3 # mass resolution 2.12 * 10**9 h0 M_sun (Baugh 2019)
    mmax = 15. #15.
    dm = 0.1
    mbins = np.arange(mmin, mmax, dm)
    mhist = mbins + dm * 0.5
    gsmf = np.zeros((len(mhist)))

    # Initialize SFRF
    smin = -6.
    smax = 3.5
    ds = 0.1
    sbins = np.arange(smin, smax, ds)
    shist = sbins + ds * 0.5
    sfrf = np.zeros((len(shist)))

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
    ytit = "log$_{10}(\\rm SFR/M_{\odot}yr^{-1})$"
    xmin = mmin; xmax = 11.6; ymin = smin;  ymax = smax
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xtit); ax.set_ylabel(ytit)

    # GSMF
    axm = plt.subplot(gs[0, :-1],sharex=ax)
    ytit="log$_{10}(\Phi(M_*))$" ; axm.set_ylabel(ytit)
    axm.set_autoscale_on(False) ;  axm.minorticks_on()
    axm.set_ylim(-4.5,-2.)
    plt.setp(axm.get_xticklabels(), visible=False)

    # SFRF
    axs = plt.subplot(gs[1:, 2], sharey=ax)
    xtit = "log$_{10}(\Phi(SFR))$"; axs.set_xlabel(xtit)
    axs.set_autoscale_on(False); axs.minorticks_on()
    axs.set_xlim(-6.4, 0.0)
    start, end = axs.get_xlim()
    axs.xaxis.set_ticks(np.arange(-6., end, 1.))
    plt.setp(axs.get_yticklabels(), visible=False)

    # Data Observations

    # SFR observed

    ih = get_nheader(obsSFR)

    dataSFR = [0]*len(colsSFR)

    for ii, col in enumerate(colsSFR):
        #print(ii,col,colsSFR[ii])
        data = np.loadtxt(obsSFR,skiprows=ih, usecols=col, unpack=True)
        dataSFR[ii] = np.array(data)

    dex = dataSFR[1]-dataSFR[0]
    histSFR = dataSFR[1]-0.5*dex

    # GSM observed

    ih = get_nheader(obsGSM)

    dataGSM = [0]*len(colsGSM)

    for ii, col in enumerate(colsGSM):
        data = np.loadtxt(obsGSM,skiprows=ih, usecols=col, unpack=True)
        dataGSM[ii] = np.array(data)

    dex = dataGSM[1] - dataGSM[0]

    # Change the units from h^-2 Msun to Msun.
    histGSM = dataGSM[1] + 2*np.log10(h0) - 0.5*dex

    # Change the units from h^3 Mpc^-3 to Mpc^-3
    freqGSM = np.log10((dataGSM[2]))- 3 * np.log10(h0)


    for ii, sfr in enumerate(SFR):
        tempfile = r"example_data/tmp_"+sfr+".dat"
        if not os.path.isfile(tempfile): continue

        ih = get_nheader(tempfile) # Number of lines in header

        # Jump the header and read the provided columns
        lms = np.loadtxt(tempfile, skiprows=ih, usecols=(0), unpack=True)
        lsfr = np.loadtxt(tempfile, skiprows=ih, usecols=(3), unpack=True)
        # loh12 = np.loadtxt(tempfile, skiprows=ih, usecols=(6), unpack=True). Not necessary in this plot


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
        axm.plot(histGSM, freqGSM, 'o', color=color[ii + 1])


        # Plot SFRF
        px = sfrf; ind = np.where(px > 0.)
        y = shist[ind]; x = np.log10(px[ind])
        ind = np.where(x < 0.)
        axs.plot(x[ind], y[ind], color=color[ii],
                 linestyle=lsty[ii], label=labels[ii])

        # Plot observations SFRF
        axs.plot(dataSFR[2], histSFR, 'o', color=color[ii + 2],
                 label=''+ labelObs[ii] +'')

    leg = axs.legend(bbox_to_anchor=(1.5, 1.4), fontsize='small',
                     handlelength=1.2, handletextpad=0.4)
    # for item in leg.legendHandles:
    # item.set_visible(True)
    leg.get_texts()
    leg.draw_frame(False)

    # for col,text in zip(color,leg.get_texts()):
    #   text.set_color(color)
    #  leg.draw_frame(False)

    plotf = outplot

    # Save figures
    print('Plot: {}'.format(plotf))
    fig.savefig(plotf)

    # os.remove(r"example_data/tmp_LC.dat")
    # os.remove(r"example_data/tmp_avSFR.dat")

def test_zm(obsZF, obsGSM, colsZ ,colsGSM,labelObs,outplot,h0=None, volume=const.vol_pm, verbose=False):


        '''
       Given log10(Mstar) and (12 + log(O/H))
       get the plots to compare (12 + log(O/H)) vs log10(Mstar).
       Get the GSMF and the ZF plots. (If that exists)

       Given the observations, compare the plots with the observations too. (If that exists)


       Parameters
       ----------
       (HERE : IF THAT EXISTS)
       obsZF : string
         Name of the input file for the Z data observed.
         Expected histogram mode:
         with a column with the low value of the bin,
         a column with the high value of the bin,
         a column with the frequency in the bin,
         and a column with the error.
         These columns must be specify in the colsZ variable.

         In text files (*.dat, *txt, *.cat), columns separated by ' '.
         In csv files (*.csv), columns separated by ','.

       obsGSM : string
         Name of the input file for the GSM data observed.

         Expected histogram mode:
         with a column with the low value of the bin,
         a column with the high value of the bin,
         a column with the frequency in the bin,
         and a column with the error.
         These columns must be specify in the colsGSM variable.

         In text files (*.dat, *txt, *.cat), columns separated by ' '.
         In csv files (*.csv), columns separated by ','.

       colsZ : list (HERE: IF THAT EXISTS, change everything)
         Columns with the data required to do the observational histogram of the SFR.
         Expected: [ind_column1, ind_column2, ind_column3, ind_column4]
         where: column1 is the column with the low values of the bins, in Msun/yr,
                column2 with the high values of the bins, in Msun/yr,
                column3 with the frequency, in Mpc^-3 dex^-1
                column4 with the error, in Mpc^-3 dex^-1
       colsGSM :
         Columns with the data required to do the observational histogram of the GSM.
         Expected: [ind_column1, ind_column2, ind_column3, ind_column4]
         where: column1 is the column with the low values of the bins, in h^-2Msun,
                column2 with the high values of the bins, in h^-2Msun,
                column3 with the frequency, in h^-3 Mpc^-3,
                column4 with the error, in h^-3 Mpc^-3.

       labelObs : list of strings
         For the legend, add the name to cite the observational data source.
         ['GSM observed', 'Z observed']

       outplot : string
         Name of the output file.
         Image-type files (*.pdf, *.jpg, ...)

       h0 : float
         If not None: value of h, H0=100h km/s/Mpc.
       volume : float
         Carlton model default value = 542.16^3 Mpc^3/h^3.
         table 1: https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4922B/abstract
         If not 542.16**3. : value of the simulation volume in Mpc^3/h^3
       verbose : boolean
         Yes = print out messages

       Returns
       -------
       plot((12 + log(O/H)),log10(Mstar)), plot GSMF and plot ZF,
       all three in one grid.
       Save it in the outplot path.
    '''

        # HERE: All like above if Z observational exists. I am not sure, I don't think so.

def test_medians(outplot, volume=542.16 ** 3, verbose=False):
    '''
        HERE: CHANGE ALL THE DESCRIPTION

        Given log10(Mstar), log10(sSFR) and 12+log(O/H),
        get the plot 12+log(O/H) vs log10(Mstar) when Plotting

        plot

        Parameters
        ----------
        cols : list
            [[component1_stellar_mass,sfr,Z],[component2_stellar_mass,sfr,Z],...]
            For text or csv files: list of integers with column position.
            For hdf5 files: list of data names.
        h0 : float
            If not None: value of h, H0=100h km/s/Mpc.
        volume : float
            Carlton model default value = 542.16^3 Mpc^3/h^3.
            table 1: https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4922B/abstract
            If not 542.16**3. : value of the simulation volume in Mpc^3/h^3
        verbose : boolean
            Yes = print out messages
        Plotting : boolean
            If True run verification plots with all data.

        Returns
        -------
        plot(log10(sSFR),log10(Mstar)), plot(12+log(O/H),log10(Mstar)) : plot #Change these names later
           '''

    # Prepare the plots

    SFR = ['avSFR', 'LC']
    U_ne = ['u', 'ne']
    col_une = [0,2] # [U, ne]
    cm = plt.get_cmap('tab10')  # Colour map to draw colours from
    color = []





    # Read the data
    tempfile_io = r"example_data/tmp_LC.dat"
    '''
    plotf1 = 'C:/Users/Olivia/PRUEBAS/pruebaplot_u.pdf'
    plotf2 = 'C:/Users/Olivia/PRUEBAS/pruebaplot_ne.pdf'
    '''
    ih = get_nheader(tempfile_io)
    lms = np.loadtxt(tempfile_io, skiprows=ih, usecols=(0), unpack=True)

    # Prepare the bins
    mmin = 8.5
    mmax = 11.5
    dm = 0.1
    mbins = np.arange(mmin, (mmax + dm), dm)
    mbinsH = np.arange(mmin, mmax, dm)
    mhist = mbinsH + dm * 0.5

    for iu, une in enumerate(U_ne):

        # Prepare the figures
        plotf = outplot + '/test_medians_'+ une+'.pdf'
        col = cm(iu)
        color.append(col)  # col change for each iteration

        med=[0]*2 # HERE: Search an efficient form
        qlow = med
        qup = med


        for ii,sfr in enumerate(SFR):
            tempfile_une = r"example_data/tmp_une_" + sfr + ".dat"
            # ih = get_nheader(tempfile_une) : Not necessary, it is the same that tempfile_io always.
            data = np.loadtxt(tempfile_une,skiprows=ih,usecols=col_une,unpack=True) # data[0] = lu, data[1] = lne

        # data = np.loadtxt(tempfile_une,skiprows=ih,usecols=(0,2),unpack=True)

        #lu = np.loadtxt(tempfile_une, skiprows=ih, usecols=(0), unpack=True)

        # lne = np.loadtxt(tempfile_une,skiprows=ih,usecols=(2),unpack=True)

            # MEDIANS:
            median = perc_2arrays(mbins, lms, data[iu], 0.5)
            # median2 = perc_2arrays(mbins,lms,lne,0.5)
            #print(data[ii],median,ii,une)
            med[ii]=median

            # QUARTILES:
            up_qu = perc_2arrays(mbins, lms, data[iu], 0.75)
            qup[ii] = up_qu
            low_qu = perc_2arrays(mbins, lms, data[iu], 0.25)
            qlow[ii] = low_qu

            ind = np.where(median>const.notnum)
            #print(ind,median[ind])
            #print(med[ii].item(0))
            plt.plot(mhist[ind],median[ind],'o', color=col, label='Calculated from the ' + SFR[ii] + '')
            plt.legend()
        plt.savefig('C:/Users/Olivia/PRUEBAS/pruebaplot'+une+'.png')
        #print('hecho')

        '''
        # PLOT
            plt.xlabel('log$_{10}$(M$_*$/M$_{\odot}$)')
            plt.ylabel('log$_{10}$ (U)')
            plt.xlim(8.5, 11.2)
            plt.ylim(-5, -2)
            ind = np.where((median > const.notnum) & (low_qu > const.notnum) & (up_qu > const.notnum))
            plt.plot(mhist[ind], median[ind], 'o', color=col, label='Calculated from the ' + sfr + '')
            plt.plot(mhist[ind], low_qu[ind], '_', color=col)
            plt.plot(mhist[ind], up_qu[ind], '_', color=col)
            plt.vlines(mhist[ind], low_qu[ind], up_qu[ind], color=col, linestyles='dashed')

            plt.legend()


        plt.savefig(plotf)
        print('Plot: {}'.format(plotf))
        
        '''

        '''
        # QUARTILES:

        up_qu = perc_2arrays(mbins, lms, lu, 0.75)
        low_qu = perc_2arrays(mbins, lms, lu, 0.25)
        # up_qne  = perc_2arrays(mbins,lms,lne,0.75)
        # low_qne = perc_2arrays(mbins,lms,lne,0.25)

        # COLOR:
        col = cm(ii)  # cm(1.*ii/len(SFR));
        color.append(col)  # col va cambiando en cada iteración

        plt.style.use(style.style1)
        plt.xlabel('log$_{10}$(M$_*$/M$_{\odot}$)')
        plt.ylabel('log$_{10}$ (U)')
        plt.xlim(8.5, 11.2)
        plt.ylim(-5, -2)
        ind = np.where((up_qu > const.notnum) & (low_qu > const.notnum))
        plt.plot(mhist[ind], median1[ind], 'o', color=col, label='Calculated from the ' + sfr + '')
        plt.plot(mhist[ind], low_qu[ind], '_', color=col)
        plt.plot(mhist[ind], up_qu[ind], '_', color=col)
        plt.vlines(mhist[ind], low_qu[ind], up_qu[ind], color=col, linestyles='dashed')

    plt.legend()
    plt.savefig(plotf1)
    plt.show()

    # fig.savefig(plotf)

    # Save figures
    print('Plot: {}'.format(plotf1))

    for ii, sfr in enumerate(SFR):
        tempfile_une = 'C:/Users/Olivia/get_nebular_emission/example_data/tmp_une_' + sfr + '.dat'

        # data = np.loadtxt(tempfile_une,skiprows=ih,usecols=(0,2),unpack=True)
        # lu = np.loadtxt(tempfile_une,skiprows=ih,usecols=(0),unpack=True)
        lne = np.loadtxt(tempfile_une, skiprows=ih, usecols=(2), unpack=True)

        # MEDIANS:
        # median1 = perc_2arrays(mbins,lms,lu,0.5)
        median2 = perc_2arrays(mbins, lms, lne, 0.5)

        # QUARTILES:

        # up_qu   = perc_2arrays(mbins,lms,lu,0.75)
        # low_qu  = perc_2arrays(mbins,lms,lu,0.25)
        up_qne = perc_2arrays(mbins, lms, lne, 0.75)
        low_qne = perc_2arrays(mbins, lms, lne, 0.25)

        # COLOR:
        col = cm(ii)  # cm(1.*ii/len(SFR));
        color.append(col)  # col change in each iteration

        plt.style.use(style.style1)
        plt.xlabel('log$_{10}$(M$_*$/M$_{\odot}$)')
        plt.ylabel('log$_{10}$ (n$_e$/cm$^{-3}$)')
        plt.xlim(8.5, 11.2)
        plt.ylim(-1., 2.)
        ind = np.where((up_qne > const.notnum) & (low_qne > const.notnum))
        plt.plot(mhist[ind], median2[ind], 'o', color=col, label='Calculated from the ' + sfr + '')
        plt.plot(mhist[ind], low_qne[ind], '_', color=col)
        plt.plot(mhist[ind], up_qne[ind], '_', color=col)
        plt.vlines(mhist[ind], low_qne[ind], up_qne[ind], color=col, linestyles='dashed')

    plt.legend()
    plt.savefig(plotf2)
    plt.show()

    # fig.savefig(plotf)

    # Save figures
    print('Plot: {}'.format(plotf2))
'''

def plot_bpt(cols, h0=None, volume=542.16 ** 3, verbose=False):
    '''
        Given log10(Mstar), log10(sSFR) and 12+log(O/H),
        get the plot 12+log(O/H) vs log10(Mstar) when Plotting

        plot

        Parameters
        ----------
        cols : list
            [[component1_stellar_mass,sfr,Z],[component2_stellar_mass,sfr,Z],...]
            For text or csv files: list of integers with column position.
            For hdf5 files: list of data names.
        h0 : float
            If not None: value of h, H0=100h km/s/Mpc.
        volume : float
            Carlton model default value = 542.16^3 Mpc^3/h^3.
            table 1: https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4922B/abstract
            If not 542.16**3. : value of the simulation volume in Mpc^3/h^3
        verbose : boolean
            Yes = print out messages
        Plotting : boolean
            If True run verification plots with all data.

        Returns
        -------
        plot(log10(sSFR),log10(Mstar)), plot(12+log(O/H),log10(Mstar)) : plot #Change these names later
           '''