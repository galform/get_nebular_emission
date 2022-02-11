import matplotlib.pyplot as plt
import numpy as np
from get_nebular_emission import eml_const as const
from get_nebular_emission.stats import perc_2arrays
import get_nebular_emission.eml_style as style

SFR = ['avSFR','LC']
variable = ['u','ne']
cm = plt.get_cmap('tab10')  # Colour map to draw colours from
color = []
plt.style.use(style.style1)
fig1=plt.figure(num=1)
fig2=plt.figure(num=2)
tempfile_io = 'C:/Users/Olivia/get_nebular_emission/example_data/tmp_LC.dat'
plotf1 = 'C:/Users/Olivia/PRUEBAS/pruebaplot_u.pdf'
plotf2 = 'C:/Users/Olivia/PRUEBAS/pruebaplot_ne.pdf'

verbose = False
ih = 1

lms   = np.loadtxt(tempfile_io, skiprows=ih, usecols=(0), unpack=True)
# Prepare the bins
mmin = 8.5
mmax = 12
dm = 0.1
mbins = np.arange(mmin, (mmax + dm), dm)
mbinsH = np.arange(mmin, mmax, dm)
mhist = mbinsH + dm * 0.5

for ii,sfr in enumerate(SFR):
    tempfile_une = 'C:/Users/Olivia/get_nebular_emission/example_data/tmp_une_'+sfr+'.dat'

    #data = np.loadtxt(tempfile_une,skiprows=ih,usecols=(0,2),unpack=True)
    lu = np.loadtxt(tempfile_une,skiprows=ih,usecols=(0),unpack=True)
    #lne = np.loadtxt(tempfile_une,skiprows=ih,usecols=(2),unpack=True)


    # MEDIANS:
    median1 = perc_2arrays(mbins,lms,lu,0.5)
    #median2 = perc_2arrays(mbins,lms,lne,0.5)

    # QUARTILES:

    up_qu   = perc_2arrays(mbins,lms,lu,0.75)
    low_qu  = perc_2arrays(mbins,lms,lu,0.25)
    #up_qne  = perc_2arrays(mbins,lms,lne,0.75)
    #low_qne = perc_2arrays(mbins,lms,lne,0.25)


    # COLOR:
    col = cm(ii)  # cm(1.*ii/len(SFR));
    color.append(col)  # col va cambiando en cada iteración

    plt.style.use(style.style1)
    plt.xlabel('log$_{10}$(M$_*$/M$_{\odot}$)')
    plt.ylabel('log$_{10}$ (U)')
    plt.xlim(8.5,11.2)
    plt.ylim(-5,-2)
    ind = np.where((up_qu > const.notnum) & (low_qu > const.notnum))
    plt.plot(mhist[ind],median1[ind],'o',color=col,label='Calculated from the '+sfr+'')
    plt.plot(mhist[ind],low_qu[ind],'_',color=col)
    plt.plot(mhist[ind],up_qu[ind],'_',color=col)
    plt.vlines(mhist[ind],low_qu[ind],up_qu[ind],color=col,linestyles='dashed')


plt.legend()
plt.savefig(plotf1)
plt.show()

    #fig.savefig(plotf)



# Save figures
print('Plot: {}'.format(plotf1))

for ii,sfr in enumerate(SFR):
    tempfile_une = 'C:/Users/Olivia/get_nebular_emission/example_data/tmp_une_'+sfr+'.dat'

    #data = np.loadtxt(tempfile_une,skiprows=ih,usecols=(0,2),unpack=True)
    #lu = np.loadtxt(tempfile_une,skiprows=ih,usecols=(0),unpack=True)
    lne = np.loadtxt(tempfile_une,skiprows=ih,usecols=(2),unpack=True)


    # MEDIANS:
    #median1 = perc_2arrays(mbins,lms,lu,0.5)
    median2 = perc_2arrays(mbins,lms,lne,0.5)

    # QUARTILES:

    #up_qu   = perc_2arrays(mbins,lms,lu,0.75)
    #low_qu  = perc_2arrays(mbins,lms,lu,0.25)
    up_qne  = perc_2arrays(mbins,lms,lne,0.75)
    low_qne = perc_2arrays(mbins,lms,lne,0.25)


    # COLOR:
    col = cm(ii)  # cm(1.*ii/len(SFR));
    color.append(col)  # col change in each iteration

    plt.style.use(style.style1)
    plt.xlabel('log$_{10}$(M$_*$/M$_{\odot}$)')
    plt.ylabel('log$_{10}$ (n$_e$/cm$^{-3}$)')
    plt.xlim(8.5,11.2)
    plt.ylim(-1.,2.)
    ind = np.where((up_qne > const.notnum) & (low_qne > const.notnum))
    plt.plot(mhist[ind],median2[ind],'o',color=col,label='Calculated from the '+sfr+'')
    plt.plot(mhist[ind],low_qne[ind],'_',color=col)
    plt.plot(mhist[ind],up_qne[ind],'_',color=col)
    plt.vlines(mhist[ind],low_qne[ind],up_qne[ind],color=col,linestyles='dashed')


plt.legend()
plt.savefig(plotf2)
plt.show()

    #fig.savefig(plotf)



# Save figures
print('Plot: {}'.format(plotf2))


'''
meds = np.zeros(len(mbins)) ; meds.fill(const.notnum)
errp = np.zeros(len(mbins)) ; #meds.fill(const.notnum)
errn = np.zeros(len(mbins)) ; #meds.fill(const.notnum)


for ii, bin in enumerate(mbins):
    dd = []
    for im, m in enumerate(lms):
        if ii==0:
            if mbins[ii]-dm<=m<=mbins[ii]:
                d = lu[im]
                dd.append(d)
        else:
            if mbins[ii-1] < m <= bin:# Here : Problem with the first bin, I do not know how to count in the first bin
                d = lu[im]
                dd.append(d)
    #print(len(dd),im, ii,mbins[ii-1],bin,mbins[0],mbins[0]-dm)


    #ind = np.where(dd[ii] > const.notnum)
    median = percentiles(0.5, dd)
    errorpos = percentiles(0.5+0.68/2,dd)
    errorneg = percentiles(0.5-0.68/2,dd)

    #ind = np.where(median > const.notnum)
    #if median > const.notnum:

    meds[ii] = median
    errp[ii] = errorpos-median
    errn[ii] = errorneg-median
    #else: meds[ii] = 0


plt.style.use(style.style1)
ind=np.where(meds>const.notnum)
plt.plot(mhist[ind],meds[ind],'ob')
ind = np.where((errp>const.notnum) & (meds>const.notnum))
#ax = plt.subplots()
ind = np.where((errn>const.notnum) & (errp>const.notnum) & (meds>const.notnum))
plt.errorbar(mhist[ind],meds[ind],yerr=[errp[ind],errn[ind]],xerr=None,fmt='.b',elinewidth=0.8,capsize=3)
#
#plt.errorbar(mhist[ind],meds[ind],yerr=errn[ind])

#plt.errorbar(mhist[ind],meds[ind],yerr=err[ind],xerr=None,fmt='.b')
plt.show()




#print(medians)




for iv in range(val0):

    dd= data[iv, :]
    ind = np.where(dd > const.notnum)
    median = percentiles(0.5, dd[ind])
    meds[iv] = median
    mtext = mtext + ' ' + str(median)




                for ii in range(val0):
                    dd = data[ii,:]
                    ind = np.where(dd>-999.)
                    median = percentiles(0.5,dd[ind])
                    meds[ii] = median
                    mtext = mtext+' '+str(median)

                # Write output file
                with open(ofile, 'a') as of:
                    of.write(mtext+' \n')

print('File with medians: {}'.format(ofile))
'''
