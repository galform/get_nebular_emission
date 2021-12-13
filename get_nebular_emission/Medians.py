import sys, os.path
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from get_nebular_emission import eml_const as const
from forplots.stats import percentiles
from forplots.stats import perc_2arrays
import forplots.eml_style as style

SFR = ['avSFR','LC']

for ii,sfr in enumerate(SFR):
    tempfile_une = 'C:/Users/Olivia/get_nebular_emission/example_data/tmp_une_'+sfr+'.dat'
    tempfile_io = 'C:/Users/Olivia/get_nebular_emission/example_data/tmp_LC.dat'

    verbose = False

    ih = 1

    data = np.loadtxt(tempfile_une,skiprows=ih,usecols=(0,2),unpack=True)
    lu = np.loadtxt(tempfile_une,skiprows=ih,usecols=(0),unpack=True)
    lne = np.loadtxt(tempfile_une,skiprows=ih,usecols=(2),unpack=True)
    lms   = np.loadtxt(tempfile_io, skiprows=ih, usecols=(0), unpack=True)

    val0 = np.shape(data)[0] ; val1 = np.shape(data)[1]

    # Prepare the bins
    mmin = 8.5
    mmax = 12
    dm = 0.1
    mbins = np.arange(mmin,(mmax+dm),dm)
    mbinsH = np.arange(mmin, mmax, dm)
    mhist = mbinsH + dm * 0.5

    medians = np.zeros((len(mhist)))



    median1 = perc_2arrays(mbins,lms,lu,0.5)

    #error = np.zeros((2, len(median1)))
    up_q = perc_2arrays(mbins,lms,lu,0.75)
    low_q =perc_2arrays(mbins,lms,lu,0.25)
    ind = np.where((up_q != const.notnum) & (low_q != const.notnum))
    error = np.row_stack((low_q[ind],up_q[ind]))

    plt.xlim(8.5,11.2)
    plt.ylim(-5,-2.5)
    plt.plot(mhist[ind],median1[ind],'ob')
    plt.plot(mhist[ind],low_q[ind],'_')
    plt.plot(mhist[ind],up_q[ind],'_')
    plt.vlines(mhist[ind],low_q[ind],up_q[ind])

    #plt.errorbar(mhist[ind],median1[ind],yerr=error,fmt='.b',elinewidth=0.8,capsize=3.)
    plt.show()

    #errorpos1 = perc_2arrays(mbins,lms,lu,(0.5+0.68/2))
    #errorneg1 = perc_2arrays(mbins,lms,lu,(0.5-0.68/2))


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
