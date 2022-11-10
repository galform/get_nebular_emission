import numpy as np
import matplotlib
import os.path
import h5py
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import get_nebular_emission.eml_style as style
from get_nebular_emission.stats import perc_2arrays
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_io import get_nheader, check_file
from get_nebular_emission.eml_photio import get_lines_Gutkin

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
       Get the GSMF and the ZF plots.

       Given the observations, compare the plots with the observations too.


       Parameters
       ----------

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

       colsZ : list
         Columns with the data required to do the observational histogram of the Z.
         Expected: [ind_column1, ind_column2, ind_column3, ind_column4]
         where: column1 is the column with the low values of the bins,
                column2 with the high values of the bins,
                column3 with the frequency,
                column4 with the error,
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

def test_medians(outplot, verbose=False):
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
            Yes = print out messages


        Returns
        -------
        Medians U and ne for each bin of mass in two differents plots.

    '''

    # Prepare the plots
    
    SFR = ['avSFR', 'LC']
    U_ne = ['u', 'ne']
    col_une = [0,2] # [U, ne]
    cm = plt.get_cmap('tab10')  # Colour map to draw colours from
    #color = []

    # Prepare the bins
    mmin = 8.5
    mmax = 11.5
    dm = 0.1
    mbins = np.arange(mmin, (mmax + dm), dm)
    mbinsH = np.arange(mmin, mmax, dm)
    mhist = mbinsH + dm * 0.5

    for iu, une in enumerate(U_ne):

        # Prepare the figures
        plt.figure()
        plotf = outplot + '/test_medians_'+ une+'.pdf'

        for ii,sfr in enumerate(SFR):
            
            # Read the data
            tempfile_io = r"example_data/tmp_" + sfr + ".dat"

            ih = get_nheader(tempfile_io)
            lms = np.loadtxt(tempfile_io, skiprows=ih, usecols=(0), unpack=True)
            
            col = cm(ii)
            #color.append(col)  # col change for each iteration
            tempfile_une = r"example_data/tmp_une_" + sfr + ".dat"
            # ih = get_nheader(tempfile_une) : Not necessary, it is the same that tempfile_io always.
            data = np.loadtxt(tempfile_une,skiprows=ih,usecols=col_une,unpack=True) # data[0] = lu, data[1] = lne

            # MEDIANS:
            median = perc_2arrays(mbins, lms, data[iu], 0.5)


            # QUARTILES:
            up_qu = perc_2arrays(mbins, lms, data[iu], 0.55)
            #qup[ii] = up_qu
            low_qu = perc_2arrays(mbins, lms, data[iu], 0.45)
            #qlow[ii] = low_qu
            
            qu = np.append([low_qu],[up_qu],axis=0)

            ind = np.where(median>const.notnum)

            plt.plot(mhist[ind],median[ind],'o', color=col, label='Calculated from the ' + SFR[ii] + '')
            #plt.errorbar(mhist[ind],median[ind],yerr=qu,marker='o',elinewidth=0.5, color=col, label='Calculated from the ' + SFR[ii] + '')
        plt.legend()
        plt.savefig(plotf)
        plt.close()



def test_bpt(outplot, photmod='gutkin16',verbose=False):
    '''
        Run a test of the interpolations done in eml_photio.
        Two plots, one to verify the U interpolation and the other one to verify the Z interpolation
        ----------
        outplot : string
            Path to the folder plot.

        verbose : boolean
            Yes = print out messages

        Returns
        -------
        plot of BPT diagram with the photoionisation model data and two points of our interpolation
           '''

    if photmod == 'gutkin16':
        print(get_lines_Gutkin(Testing=True,Plotting=True,verbose=True))
    else:
        print('STOP (eml_plots): Unrecognised model to get the interpolations ({})'.format(photmod))
        print('                  Possible photmod = {}'.format(const.photmods))
        exit()

    Z = ['0001', '0002', '0005', '001'] #, '002', '004', '006', '008', '010', '014', '017', '020', '030', '040']

    zz = [0.0001, 0.0002, 0.0005, 0.001] #, 0.002, 0.004, 0.006, 0.008, 0.01, 0.014, 0.017, 0.02, 0.03, 0.04]

    uu = [-1., -1.5, -2., -2.5, -3., -3.5, -4.]

    ne = ['100']  # ['10', '100', '1000','10000']

    cm = plt.get_cmap('tab20') # Colour map to draw colours from

    for iz, zname in enumerate(Z):
        infile = r"nebular_data/gutkin_tables/nebular_emission_Z" + zname + ".txt"

        ih = get_nheader(infile)

        datane = np.loadtxt(infile, skiprows=ih, usecols=(2), unpack=True)
        datalu = np.loadtxt(infile, skiprows=ih, usecols=(0), unpack=True)

        OIII5007 = np.loadtxt(infile, skiprows=ih, usecols=(8), unpack=True)
        Hb = np.loadtxt(infile, skiprows=ih, usecols=(6), unpack=True)
        NII6548 = np.loadtxt(infile, skiprows=ih, usecols=(9), unpack=True)
        Ha = np.loadtxt(infile, skiprows=ih, usecols=(10), unpack=True)

        for ii, nh in enumerate(ne):
            outfile = r"output_data/Gutkinfile_n_" + nh + ".txt"

            header1 = 'Z, U, NII6584[ind]/Ha[ind]  OIII5007[ind]/Hb[ind]'

            ind = np.where(datane == float(nh))
            x = np.log10(NII6548[ind] / Ha[ind])
            y = np.log10(OIII5007[ind] / Hb[ind])
            u = datalu[ind]
            z = np.full(np.shape(u), zz[iz])

            tofile = np.column_stack((z, u, x, y))

            with open(outfile, 'a') as outf:
                if iz == 0:
                    np.savetxt(outf, tofile, delimiter=' ', header=header1)
                else:
                    np.savetxt(outf, tofile, delimiter=' ')
                outf.closed

    cols = []
    for iz, lz in enumerate(uu):
        col = cm(iz)
        cols.append(col)

    for ii, nh in enumerate(ne):
        infile = r"output_data/Gutkinfile_n_" + nh + ".txt"

        ih = get_nheader(infile)

        z = np.loadtxt(infile, skiprows=ih, usecols=(0), unpack=True)
        u = np.loadtxt(infile, skiprows=ih, usecols=(1), unpack=True)
        x = np.loadtxt(infile, skiprows=ih, usecols=(2), unpack=True)
        y = np.loadtxt(infile, skiprows=ih, usecols=(3), unpack=True)

        #os.remove(infile)

        file = r"output_data/output_kashino20_test.hdf5"
        check_file(file, verbose=True)
        f = h5py.File(file, 'r')
        header = f['header']
        data = f['data']

        # print(list(data.keys()))

        lu_disk = data['lud'][:]
        lu_bulge = data['lub'][:]
        lne_disk = data['lned'][:]
        lne_bulge = data['lneb'][:]
        loh12_disk = (data['loh12d'][:])
        loh12_bulge = (data['loh12b'][:])
        neblines_disk = data['nebline_disk'][:]
        neblines_bulge = data['nebline_bulge'][:]

        # print(np.shape(neblines_disk)) ###### SHAPE neblines_disk >>> (18,n) 18 lines n data for each line

        # Take the galaxies with ne >= 2:
        ind = np.where(lne_disk >= 2)

        Hbeta = neblines_disk[1][ind]

        OIII5007 = neblines_disk[3][ind]
        NII6548 = neblines_disk[4][ind]
        Halpha = neblines_disk[5][ind]

        lne_disk = lne_disk[ind]
        lne_bulge = lne_bulge[ind]
        loh12_disk = loh12_disk[ind]
        loh12_bulge = loh12_bulge[ind]

        lu_disk = lu_disk[ind]

        lu_bulge = lu_bulge[ind]

        my_x = np.log10(NII6548 / Halpha)
        my_y = np.log10(OIII5007 / Hbeta)

        # DIFERENTS COLORS FOR U:

        plt.figure()
        for iu, lu in enumerate(uu):
            ind = np.where(u == uu[iu])
            plt.plot(x[ind], y[ind], marker='.', linewidth=0, color=cols[iu], label='U = ' + str(lu) + '')


        labelsU = ['U = {}'.format(lu_disk[0]), 'U = {}'.format(lu_disk[1])]

        plt.plot(my_x, my_y, marker='o', linewidth=0, color='black')  # ,label = 'U = {}'.format(lu_disk))
        plt.xlabel('log$_{10}$([NII]$\\lambda$6584/H$\\alpha$)')
        plt.ylabel('log$_{10}$([OIII]$\\lambda$5007/H$\\beta$)')

        for i, label in enumerate(labelsU):
            plt.annotate(label, (my_x[i], my_y[i]))
        plotnom = outplot + '/BPTplot_n_' + nh + '_U_test.png'
        plt.legend()

        plt.savefig(plotnom)
        plt.close()
        #plt.show()

        # DIFFERENTS COLORS FOR Z:

        plt.figure()
        for iz, lz in enumerate(zz):
            ind = np.where(z == zz[iz])
            plt.plot(x[ind], y[ind], marker='.', linewidth=0, color=cols[iz], label='Z = ' + str(lz) + '')

        labelsZ = ['Z = {:.1E}'.format(10 ** (loh12_disk[0])), 'Z = {:.1E}'.format(10 ** (loh12_disk[1]))]
        plt.plot(my_x, my_y, marker='o', linewidth=0, color='black')

        plt.xlabel('log$_{10}$([NII]$\\lambda$6584/H$\\alpha$)')
        plt.ylabel('log$_{10}$([OIII]$\\lambda$5007/H$\\beta$)')
        for i, label in enumerate(labelsZ):
            plt.annotate(label, (my_x[i], my_y[i]))
        plotnom = outplot + '/BPTplot_n_' + nh + '_Z_test.png'
        plt.legend()

        plt.savefig(plotnom)
        plt.close()
        
        f.close()
        #plt.show()
        #os.remove(file)