import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import get_nebular_emission.eml_style as style
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_io import get_ncomponents
from get_nebular_emission.eml_io import get_nheader
import time
inicio = time.time()

plt.style.use(style.style1)

def test_ssfr(cols, h0=None, volume=const.vol_pm, verbose=False): # Poner volumen como constante con las constantes

    '''
       Given log10(Mstar), log10(sSFR) and 12+log(O/H),
       get the plots log10(sSFR) vs log10(Mstar)
       and 12+log(O/H) vs log10(Mstar) when Plotting OTRA FUNCIÓN PARA ESTO:  test_zm
       test_medians
       plot_bpt

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

    #Here: We start with sSFR vs M. Add later Z vs M.

    # We call the file with the data
    tempfile = r"example_data/tmp2.dat"

    # Correct the units of the simulation volume to Mpc^3:
    if h0:
        volume=volume*(h0**3)

    #Prepare the plot

    #color = ['navy', 'royalblue', 'lightsteelblue']
    lsty = ['-', '--', ':']

    # Initialize GSMF (Galaxy Cosmological Mass Function)
    mmin = 8.5
    mmax = 15.
    dm = 0.1
    mbins = np.arange(mmin, mmax, dm)
    mhist = mbins + dm * 0.5

    # Initialize SSFR
    smin = 3.
    smax = 13.
    ds = 0.1
    sbins = np.arange(smin, smax, ds)
    shist = sbins + ds * 0.5

    # Plots limits and style
    fig = plt.figure(figsize=(8.5, 9.))
    gs = gridspec.GridSpec(3, 3)
    gs.update(wspace=0., hspace=0.)
    ax = plt.subplot(gs[1:, :-1])

    # Fig. sSFR vs M

    xtit = "$log_{10}(\\rm M_{*}/M_{\odot})$"#h^{-1})$"
    ytit = "$log_{10}(\\rm SFR/M_{\odot}yr^{-1})$"#h^{-1}
    xmin = mmin; xmax = mmax; ymin = smin;  ymax = smax
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xtit); ax.set_ylabel(ytit)

    # GSMF
    axm = plt.subplot(gs[0, :-1],sharex=ax)
    ytit="$log_{10}(\Phi)$" ; axm.set_ylabel(ytit)
    axm.set_autoscale_on(False) ;  axm.minorticks_on()
    axm.set_ylim(-5.5,-1.)
    plt.setp(axm.get_xticklabels(), visible=False)

    # SFRF
    axs = plt.subplot(gs[1:, 2], sharey=ax)
    xtit = "$log_{10}(\Phi)$"; axs.set_xlabel(xtit)
    axs.set_autoscale_on(False); axs.minorticks_on()
    axs.set_xlim(-4.4, 0.0)
    start, end = axs.get_xlim()
    axs.xaxis.set_ticks(np.arange(-4., end, 1.))
    plt.setp(axs.get_yticklabels(), visible=False)

    #print(range(len(lssfr[0])-1))
    #exit()
    ncomp = get_ncomponents(cols)

    #lms = np.array([clms[:,0]])
    #print(lssfr)
    #exit()

    cm = plt.get_cmap('tab10')  # Colour map to draw colours from
    color = []


    #for ii in range(ncomp):

    ih = get_nheader(tempfile)

    # Jump the header and read the provided columns
    with open(tempfile, "r") as ff:
        for il, line in enumerate(ff):
            # Read until the data starts:
            if (il<ih): continue

            # Read data
            allcols = line.split()

            lms = np.array([[float(allcols[0])]])
            comp = np.array([[float(allcols[1])]])
            lms = np.append(lms, comp, axis=0)

            lssfr = np.array([[float(allcols[2])]])
            comp = np.array([[float(allcols[3])]])
            lssfr = np.append(lssfr, comp, axis=0)

            loh12 = np.array([[float(allcols[4])]])
            comp = np.array([[float(allcols[5])]])
            loh12 = np.append(loh12, comp, axis=0)
            for ii in range(ncomp):
                # GSMF
                H, bins_edges = np.histogram(lms[ii], bins=np.append(mbins, mmax))
                gsmf = H / volume / dm  # In Mpc^3/h^3
                #if mmin > H.min() and mmax < H.max():
                #   print('WARNING (eml_testplots): mass out of range for the plots')

                # sSFR
                H, bins_edges = np.histogram(lssfr[ii],bins=np.append(sbins,smax))
                sfrf = H/volume/ds
                #print(sfrf.max(), sfrf.min())

                # sSFR-GSMF
                H, xedges, yedges = np.histogram2d(lssfr[ii], lms[ii],
                                        bins=[np.append(sbins,smax),
                                        np.append(mbins,mmax)])
                smf = H/volume/dm/ds

                col = cm(ii); color.append(col) #col va cambiando en cada iteración

                #Plot SMF vs SFR
                matplotlib.rcParams['contour.negative_linestyle'] = lsty[ii]
                zz = np.zeros(shape=(len(shist), len(mhist))) ; zz.fill(const.notnum)
                ind = np.where(smf>0.)
                zz[ind] = np.log10(smf[ind])

                ind = np.where(zz>const.notnum)
                if (np.shape(ind)[1] > 1):
                    # Contours
                    xx, yy = np.meshgrid(mbins, sbins)
                    cs = ax.contour(xx, yy, zz, colors=col)

                # GSMF
                py = gsmf; ind = np.where(py>0.)
                x = mhist[ind]; y = np.log10(py[ind])

                ind = np.where(y < 0.)
                axm.plot(x[ind], y[ind], color=col,
                        linestyle=lsty[ii])


                # SFRF
                px = sfrf; ind = np.where(px>0.)
                y = shist[ind]; x = np.log10(px[ind])
                ind = np.where(x < 0.)
                axs.plot(x[ind], y[ind], color=col,
                                    linestyle=lsty[ii])#,label='$n_{\\rm gal}=10^{'+str(iic)+'}{\\rm Mpc}^{-3}h^{-3}$')

        plotf = 'C:/Users/Olivia/PRUEBAS/' + 'pruebaplot.pdf'
        # Save figures
        print('Plot: {}'.format(plotf))
        fig.savefig(plotf)

def test_zm(cols, h0=None, volume=542.16 ** 3., verbose=False):
    '''
          Given log10(Mstar), log10(sSFR) and 12+log(O/H),
          get the plot 12+log(O/H) vs log10(Mstar) when Plotting
          test_medians
          plot_bpt

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

def test_medians(cols, h0=None, volume=542.16 ** 3, verbose=False):
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