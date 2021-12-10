# 4 arrays de numeros aleatorios: lineas. Fichero con 4 columnas. Plot bpt que lea y que pinte.
import numpy as np
import os
from matplotlib import pyplot as plt
from get_nebular_emission.forplots.eml_style import style1
import get_nebular_emission.eml_const as const
bpt_data = 'C:/Users/Olivia/PRUEBAS/bpt_data.dat' # r"example_data/bpt_data.dat"
header1 = 'line1   line2   line3   line4'

a = np.random.randint(1, 10, size=10)
b = np.random.randint(1, 10, size=10)
c = np.random.randint(1, 10, size=10)
d = np.random.randint(1, 10, size=10)

datatofile=np.column_stack((a, b, c, d))

with open(bpt_data,'w') as svfile:
    np.savetxt(svfile,datatofile,delimiter=' ',header=header1)
    svfile.closed

line1 = np.loadtxt(bpt_data,skiprows=1,usecols=(0),unpack=True)
line2 = np.loadtxt(bpt_data,skiprows=1,usecols=(1),unpack=True)
line3 = np.loadtxt(bpt_data,skiprows=1,usecols=(2),unpack=True)
line4 = np.loadtxt(bpt_data,skiprows=1,usecols=(3),unpack=True)

ax = line1-line2 # As they are logarithms
ay = line3-line4


plt.style.use(style1)
ind=np.where((ax > const.notnum) & (ay > const.notnum))
plt.plot(ax[ind],ay[ind],'ob',label = 'GALFORM')
plt.xlabel('log$_{10}$(line$_1$/line$_2$)')
plt.ylabel('log$_{10}$(line$_3$/line$_4$)')
plt.legend()


plt.show()

os.remove(bpt_data)
