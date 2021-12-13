import numpy as np
import matplotlib
import os.path
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import get_nebular_emission.eml_const as const
from forplots import eml_style as style


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

# Here: We start with sSFR vs M. Add later Z vs M.

# We call the file with the data
'''
LC2sfr = False
if LC2sfr:
    tempfile = 'C:/Users/Olivia/get_nebular_emission/example_data/tmp_LC.dat'
else:
    tempfile = 'C:/Users/Olivia/get_nebular_emission/example_data/tmp_avSFR.dat'
'''

# Correct the units of the simulation volume to Mpc^3
volume = (542.16**3)/200.
h0 = 0.6777

# Prepare the plot

# color = ['navy', 'royalblue', 'lightsteelblue']
lsty = ['-',(0,(2,3))]


# Initialize GSMF (Galaxy Cosmological Mass Function)
mmin = 10.3 # Buscar resolucion en masa
mmax = 12. #15.
dm = 0.1
mbins = np.arange(mmin, mmax, dm)
mhist = mbins + dm * 0.5
gsmf = np.zeros((len(mhist)))


# Initialize SFRF
smin = -6. #3.
smax = 3.5 #13.
ds = 0.1
sbins = np.arange(smin, smax, ds)
shist = sbins + ds * 0.5
sfrf = np.zeros((len(shist)))

# Initialize SFR vs M function
lens= len(shist)
lenm= len(mhist)
smf = np.zeros((lens,lenm))


# Plots limits and style
fig = plt.figure(figsize=(8.5, 9.))
gs = gridspec.GridSpec(3, 3)
gs.update(wspace=0., hspace=0.)
ax = plt.subplot(gs[1:, :-1])

# Fig. SFR vs M
xtit = "$log_{10}(\\rm M_{*}/M_{\odot})$"  # h^{-1})$"
ytit = "$log_{10}(\\rm SFR/M_{\odot}yr^{-1})$"  # h^{-1}
xmin = mmin;
xmax = mmax;#11.9;
ymin = smin;#6;
ymax = smax;
ax.set_xlim(xmin, xmax);
ax.set_ylim(ymin, ymax)
ax.set_xlabel(xtit);
ax.set_ylabel(ytit)

# GSMF
axm = plt.subplot(gs[0, :-1], sharex=ax)
ytit = "$log_{10}(\Phi(M_*))$";
axm.set_ylabel(ytit)
axm.set_autoscale_on(False);
axm.minorticks_on()
axm.set_ylim(-5.5, -2.)
plt.setp(axm.get_xticklabels(), visible=False)

# SFRF
axs = plt.subplot(gs[1:, 2], sharey=ax)
xtit = "$log_{10}(\Phi(SFR))$";
axs.set_xlabel(xtit)
axs.set_autoscale_on(False);
axs.minorticks_on()
axs.set_xlim(-6.4, 0.)
start, end = axs.get_xlim()
axs.xaxis.set_ticks(np.arange(-6., end, 1.))
plt.setp(axs.get_yticklabels(), visible=False)


fobs_sfrf = 'C:/Users/Olivia/Desktop/Prácticas externas/Proyecto/Código Carlton/get_nebular_emission' \
            '/cmb/gruppioni_2015_z2.0-2.5_cha.txt'
#if not os.path.isfile(fobs_sfrf):
 #   continue

fobs_gsmf = 'C:/Users/Olivia/Desktop/Prácticas externas/Proyecto/Código Carlton/get_nebular_emission' \
            '/cmb/henriques_2014_z2_cha.txt'
#if not os.path.isfile(fobs_gsmf):
 #   continue

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

ih = 1

nds = np.array([-2., -3., -4.])
al = np.sort(nds)


SFR = ['LC','avSFR']
labels = ['average SFR', 'SFR from LC photons']
obs_labels = ['mstars_total','SFR_total']
cm = plt.get_cmap('tab10')  # Colour map to draw colours from
color = []
for ii in range(0,10):
    col = cm(ii)  # cm(1.*ii/len(SFR));
    color.append(col)  # col va cambiando en cada iteración

for ii, sfr in enumerate(SFR):
    tempfile = 'C:/Users/Olivia/get_nebular_emission/example_data/tmp_'+sfr+'.dat'
    if not os.path.isfile(tempfile):
        continue

    lms   = np.loadtxt(tempfile, skiprows=ih, usecols=(0), unpack=True)
    lsfr = np.loadtxt(tempfile, skiprows=ih, usecols=(3), unpack=True)
    loh12 = np.loadtxt(tempfile, skiprows=ih, usecols=(6), unpack=True)

    H, bins_edges = np.histogram(lms, bins=np.append(mbins, mmax))
    gsmf = H / volume / dm  # In Mpc^3/h^3

    H, bins_edges = np.histogram(lsfr, bins=np.append(sbins, smax))
    sfrf = H / volume / ds

    H, xedges, yedges = np.histogram2d(lsfr, lms,
                                    bins=([np.append(sbins, smax),
                                     np.append(mbins, mmax)]))
    smf = H / volume / dm / ds

    plt.style.use(style.style1)

    # Plot SMF vs SFR
    matplotlib.rcParams['contour.negative_linestyle'] = lsty[ii]
    zz = np.zeros(shape=(len(shist), len(mhist))); zz.fill(const.notnum)
    ind = np.where(smf > 0.)
    zz[ind] = np.log10(smf[ind])

    ind = np.where(zz > const.notnum)
    if (np.shape(ind)[1] > 1):
        # Contours
        xx, yy = np.meshgrid(mbins, sbins)
        cs = ax.contour(xx, yy, zz, levels=al, colors=color[ii])

        #labels = ['average SFR', 'SFR from LC photons']
        ax.clabel(cs,inline=1,fontsize=10)
        #for i in range(len(labels)):
        #    cs.collections[i].set_label(labels[i])

    # GSMF
    py = gsmf;
    ind = np.where(py > 0.)
    x = mhist[ind];
    y = np.log10(py[ind])
    ind = np.where(y < 0.)

    axm.plot(x[ind], y[ind], color=color[ii],
             linestyle=lsty[ii], label=labels[ii])
    axm.plot(ghistObslms,freqObs_lms,'o',color=color[ii+1],label='log$_{10}$(mstars_total) observed'+SFR[ii]+'')
    #plt.legend()

    # SFRF
    px = sfrf;
    ind = np.where(px > 0.)
    y = shist[ind];
    x = np.log10(px[ind])
    ind = np.where(x < 0.)
    axs.plot(x[ind], y[ind], color=color[ii],
            linestyle=lsty[ii],label=labels[ii])  # ,label='$n_{\\rm gal}=10^{'+str(iic)+'}{\\rm Mpc}^{-3}h^{-3}$')

    axs.plot(freqObs_sfr, ghistObsSFR, 'o', color=color[ii + 2],label='log$_{10}$('+obs_labels[ii]+') observed')


leg = axs.legend(bbox_to_anchor=(1.5, 1.4),fontsize='small', \
                handlelength=1.2, handletextpad=0.4)
#for item in leg.legendHandles:
    #item.set_visible(True)
leg.get_texts()
leg.draw_frame(False)

#for col,text in zip(color,leg.get_texts()):
 #   text.set_color(color)
  #  leg.draw_frame(False)

plotf = 'C:/Users/Olivia/PRUEBAS/' + 'pruebaplot.pdf'
# Save figures
print('Plot: {}'.format(plotf))
fig.savefig(plotf)


'''
# Jump the header and read the provided columns
with open(tempfile, "r") as ff:
    for il, line in enumerate(ff):
        # Read until the data starts:
        if (il < ih): continue

        # Read data
        allcols = line.split()

        lms1 = np.array([[float(allcols[0])]])
        lms2 = np.array([[float(allcols[1])]])
        lms = np.append(lms1,lms2,axis=0)
        #print(lms[1])
    
        lssfr1 = np.array([[float(allcols[2])]])
        lssfr2 = np.array([[float(allcols[3])]])
        lssfr = np.append(lssfr1,lssfr2,axis=0)

        loh121 = np.array([[float(allcols[4])]])
        loh122 = np.array([[float(allcols[5])]])
        loh12 = np.append(loh121,loh122,axis=0)

        for ii in range(ncomp):
            #if mmin > lms[ii] or mmax < lms[ii]:
             #   print('WARNING (eml_testplots): mass out of range for the plots', lms[ii])
            #if smin > lssfr[ii] or smax < lssfr[ii]:
             #   print('WARNING (eml_testplots): SFR out of range for the plots', lssfr[ii])
            #if mmin > loh12[ii] or mmax < loh12[ii]:
            #   print('WARNING (eml_testplots): mass out of range for the plots')


            # GSMF
            H, bins_edges = np.histogram(lms[ii], bins=np.append(mbins, mmax))
            H = H / volume / dm  # In Mpc^3/h^3
            gsmf[:,ii] += H
        print(lms)
    #print(gsmf)
           
            # sSFR
            H, bins_edges = np.histogram(lssfr[ii], bins=np.append(sbins, smax))
            H = H / volume / ds
            #H = np.log10(H,where=0<H)
            sfrf[:,ii] += H

            # sSFR-GSMF
            H, xedges, yedges = np.histogram2d(lssfr[ii], lms[ii],
                                            bins=[np.append(sbins, smax),
                                                    np.append(mbins, mmax)])
            H = H / volume / dm / ds
            H = np.log10(H,where=0<H)
            smf[ii] +=H

    #print(H)
    #smf = H / volume / dm / ds
    
        
        for ii in range(ncomp):
        col = cm(ii);
            color.append(col)  # col va cambiando en cada iteración

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
                cs = ax.contour(xx, yy, zz, colors=col)

            # GSMF
            py = gsmf;
            ind = np.where(py > 0.)
            x = mhist[ind];
            y = np.log10(py[ind])

            ind = np.where(y < 0.)
            axm.plot(x[ind], y[ind], color=col,
                     linestyle=lsty[ii])

            # SFRF
            px = sfrf;
            ind = np.where(px > 0.)
            y = shist[ind];
            x = np.log10(px[ind])
            ind = np.where(x < 0.)
            axs.plot(x[ind], y[ind], color=col,
                     linestyle=lsty[ii])  # ,label='$n_{\\rm gal}=10^{'+str(iic)+'}{\\rm Mpc}^{-3}h^{-3}$')

    plotf = 'C:/Users/Olivia/PRUEBAS/' + 'pruebaplot.pdf'
    # Save figures
    print('Plot: {}'.format(plotf))
    fig.savefig(plotf)
'''