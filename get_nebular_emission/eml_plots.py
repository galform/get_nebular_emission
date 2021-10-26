import numpy as np
import os.path, sys
import subprocess
from Cosmology import *
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import ndimage
from stats import *
import get_nebular_emission.eml_style as style
plt.style.use(style.style1)
import get_nebular_emission.eml_const as const

def get_plots(lms, lssfr, loh12, verbose=False):
    '''
       Given log10(Mstar), log10(sSFR) and 12+log(O/H),
       get the plots log10(sSFR) vs log10(Mstar)
       and 12+log(O/H) vs log10(Mstar) when Testing

       Parameters
       ----------
       lms : float
         log10(Mstar/Msun)
       lssfr : float
         log10(sSFR/yr), it should be an instantaneous measurement
       loh12 : float
         12+log(O/H)
       verbose : boolean
         Yes = print out messages

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

    #Prepare the plot
    volume = (500. ** 3.) #Here: Add later to constant file.
                          #Here: Allow for differents volumes
    cuts = ['lo2', 'sfr', 'm']
    lsty = ['-', '--', ':']

    # Initialize GSMF
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
    xtit = "$log_{10}(\\rm M_{*}/M_{\odot})$"
    ytit = "$log_{10}(\\rm sSFR/yr)$"
    xmin = mmin; xmax = 11.9; ymin = 6.;  ymax = smax
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

    for iic,cut in enumerate(cuts):
        # GSMF
        H, bins_edges = np.histogram(lms, bins=np.append(mbins, mmax))
        gsmf = H / volume / dm  # In Mpc^3/h^3

        # sSFR
        H, bins_edges = np.histogram(lssfr,bins=np.append(sbins,smax))
        sfrf = H/volume/ds

        # sSFR-GSMF
        H, xedges, yedges = np.histogram2d(lssfr,lms,
                                       bins=[np.append(sbins,smax),
                                             np.append(mbins,mmax)])
        smf = H/volume/dm/ds

        #Plot SMF vs SFR
        matplotlib.rcParams['contour.negative_linestyle'] = lsty[iic]
        zz = np.zeros(shape=(len(shist),len(mhist))) ; zz.fill(const.notnum)
        ind = np.where(smf>0.)
        zz[ind] = np.log10(smf[ind])

        ind = np.where(zz > -999.)
        if (np.shape(ind)[1] > 1):
            # Contours
            xx, yy = np.meshgrid(mbins, sbins)
            # al = nds[iin] ; print(al)
            cs = ax.contour(xx, yy, zz, levels=al,colors=cols[iin])
            # cs.levels = [nf(val) for val in cs.levels]
            # ax.clabel(cs, cs.levels, inline=1,inline_spacing=0,\
            #  fontsize=10,fmt='%r')#fmt='%r %%')

        # GSMF
        py = gsmf; ind = np.where(py > 0.)
        x = mhist[ind]; y = np.log10(py[ind])
        ind = np.where(y < 0.)
        axm.plot(x[ind], y[ind], color=cols[iin], linestyle=lsty[iic])

        # SFRF
        px = sfrf ; ind = np.where(px>0.)
        y = shist[ind] ; x = np.log10(px[ind])
        ind = np.where(x < 0.)
        if (iic == 1):
          inleg = '$n_{\\rm gal}=10^{'+str(nds[iin])+'}{\\rm Mpc}^{-3}h^{-3}$'
          axs.plot(x[ind],y[ind],color=cols[iin], linestyle=lsty[iic],label=inleg)
        else:
            if (iin == 2 and iic == 0):
                axs.plot([],[],' ',
                         label=survey+', z='+zz_list[iiz])

            axs.plot(x[ind],y[ind],color=cols[iin],
                     linestyle=lsty[iic])