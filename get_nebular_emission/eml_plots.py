import numpy as np
import matplotlib
import os.path

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import get_nebular_emission.eml_style as style
import get_nebular_emission.eml_const as const
from get_nebular_emission.eml_io import get_ncomponents
from get_nebular_emission.eml_io import get_nheader

plt.style.use(style.style1)

def test_sfrf(cols, h0=None, volume=const.vol_pm, verbose=False): # Poner volumen como constante con las constantes

    '''
       Given log10(Mstar), log10(sSFR) and 12 + log(O/H),
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


    # Correct the units of the simulation volume to Mpc^3:
    if h0:
        volume=volume*(h0**3)

    #Prepare the plot
    lsty = ['-',(0,(2,3))] # Line form

    nds = np.array([-2., -3., -4.]) # Contours values
    al = np.sort(nds)

    SFR = ['LC', 'avSFR']
    labels = ['average SFR', 'SFR from LC photons']
    obs_labels = ['mstars_total', 'SFR_total'] # Here: Add Bibliography
    cm = plt.get_cmap('tab10')  # Colour map to draw colours from
    color = []
    for ii in range(0, 10):
        col = cm(ii)  # cm(1.*ii/len(SFR));
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
    xtit = "$log_{10}(\\rm M_{*}/M_{\odot})$"
    ytit = "$log_{10}(\\rm SFR/M_{\odot}yr^{-1})$"
    xmin = mmin; xmax = 11.6; ymin = smin;  ymax = smax
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xtit); ax.set_ylabel(ytit)

    # GSMF
    axm = plt.subplot(gs[0, :-1],sharex=ax)
    ytit="$log_{10}(\Phi)$" ; axm.set_ylabel(ytit)
    axm.set_autoscale_on(False) ;  axm.minorticks_on()
    axm.set_ylim(-4.5,-2.)
    plt.setp(axm.get_xticklabels(), visible=False)

    # SFRF
    axs = plt.subplot(gs[1:, 2], sharey=ax)
    xtit = "$log_{10}(\Phi)$"; axs.set_xlabel(xtit)
    axs.set_autoscale_on(False); axs.minorticks_on()
    axs.set_xlim(-6.4, 0.0)
    start, end = axs.get_xlim()
    axs.xaxis.set_ticks(np.arange(-6., end, 1.))
    plt.setp(axs.get_yticklabels(), visible=False)

    # Data Observations

    # Here: Ask how to change these paths. Put these data in gitHub ?
    fobs_sfrf = 'C:/Users/Olivia/Desktop/Prácticas externas/Proyecto/Código Carlton/get_nebular_emission' \
                '/cmb/gruppioni_2015_z2.0-2.5_cha.txt'
    #if not os.path.isfile(fobs_sfrf): continue

    fobs_gsmf = 'C:/Users/Olivia/Desktop/Prácticas externas/Proyecto/Código Carlton/get_nebular_emission' \
                '/cmb/henriques_2014_z2_cha.txt'
    # if not os.path.isfile(fobs_gsmf):continue

    # Here : Search an efficient form to do that.
    SFRobs_Low = np.loadtxt(fobs_sfrf, skiprows=4, usecols=(0), unpack=True)
    SFRobs_High = np.loadtxt(fobs_sfrf, skiprows=4, usecols=(1), unpack=True)
    freqObs_sfr = np.loadtxt(fobs_sfrf, skiprows=4, usecols=(2), unpack=True)
    errorObs_sfr = np.loadtxt(fobs_sfrf, skiprows=4, usecols=(3), unpack=True)

    dexObsSFR = SFRobs_High - SFRobs_Low
    ghistObsSFR = SFRobs_High - 0.5 * dexObsSFR

    lmsobs_Low = np.loadtxt(fobs_gsmf, skiprows=6, usecols=(0), unpack=True) + 2 * np.log10(h0)
    lmsobs_High = np.loadtxt(fobs_gsmf, skiprows=6, usecols=(1), unpack=True) + 2 * np.log10(h0)
    freqObs_lms = np.log10(np.loadtxt(fobs_gsmf, skiprows=6, usecols=(2), unpack=True)) - 3 * np.log10(h0)
    errorObs_lms = np.log10(np.loadtxt(fobs_gsmf, skiprows=6, usecols=(3), unpack=True)) - 3 * np.log10(h0)

    dexObslms = lmsobs_High - lmsobs_Low
    ghistObslms = lmsobs_High - 0.5 * dexObslms

    for ii, sfr in enumerate(SFR):
        tempfile = r"example_data/tmp_"+sfr+".dat"
        if not os.path.isfile(tempfile): continue

        ih = get_nheader(tempfile) # Number of lines of header

        # Jump the header and read the provided columns
        lms = np.loadtxt(tempfile, skiprows=ih, usecols=(0), unpack=True)
        lsfr = np.loadtxt(tempfile, skiprows=ih, usecols=(3), unpack=True)
        loh12 = np.loadtxt(tempfile, skiprows=ih, usecols=(6), unpack=True)

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
        zz = np.zeros(shape=(len(shist), len(mhist)));
        zz.fill(const.notnum)
        ind = np.where(smf > 0.)
        zz[ind] = np.log10(smf[ind])

        ind = np.where(zz > const.notnum)
        if (np.shape(ind)[1] > 1):

            # Contours
            xx, yy = np.meshgrid(mbins, sbins)
            cs = ax.contour(xx, yy, zz, levels=al, colors=color[ii])

            # labels = ['average SFR', 'SFR from LC photons']
            ax.clabel(cs, inline=1, fontsize=10)
            # for i in range(len(labels)):
            #    cs.collections[i].set_label(labels[i])

        # Plot GSMF
        py = gsmf; ind = np.where(py > 0.)
        x = mhist[ind]; y = np.log10(py[ind])
        ind = np.where(y < 0.)

        axm.plot(x[ind], y[ind], color=color[ii],
                 linestyle=lsty[ii], label=labels[ii])
        axm.plot(ghistObslms, freqObs_lms, 'o', color=color[ii + 1],
                 label='log$_{10}$(mstars_total) observed' + SFR[ii] + '')


        # Plot SFRF
        px = sfrf; ind = np.where(px > 0.)
        y = shist[ind]; x = np.log10(px[ind])
        ind = np.where(x < 0.)
        axs.plot(x[ind], y[ind], color=color[ii],
                 linestyle=lsty[ii], label=labels[ii])

        axs.plot(freqObs_sfr, ghistObsSFR, 'o', color=color[ii + 2],
                 label='log$_{10}$(' + obs_labels[ii] + ') observed')
                # Here: Change to the real bibliography

    leg = axs.legend(bbox_to_anchor=(1.5, 1.4), fontsize='small',
                     handlelength=1.2, handletextpad=0.4)
    # for item in leg.legendHandles:
    # item.set_visible(True)
    leg.get_texts()
    leg.draw_frame(False)

    # for col,text in zip(color,leg.get_texts()):
    #   text.set_color(color)
    #  leg.draw_frame(False)

    plotf = 'C:/Users/Olivia/PRUEBAS/' + 'pruebaplot.pdf' # Here: Change path

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