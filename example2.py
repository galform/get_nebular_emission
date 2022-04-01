import numpy as np
import get_nebular_emission.eml_const as const


'''REDUCED FILE'''
from get_nebular_emission.eml_io import get_reducedfile as reducedfile

infile = r"PRUEBAS/emlines_lc16_PMILL_iz245_ivol8.cat"
outfile = r"example_data/emlines_lc16_PMILL_iz245_ivol8_reduced.dat"


# reducedfile(infile,outfile,indcol=[11,15,13,14,0,1,2,4],verbose=True)

'''TEST PLOT SFRF'''

from get_nebular_emission.eml_plots import test_sfrf as test1

plot2file = r"PRUEBAS/" + "pruebaplotfeb.pdf"
SFRfile = r"get_nebular_emission" \
                "/cmb/gruppioni_2015_z2.0-2.5_cha.txt"
GMSfile = r"get_nebular_emission" \
                "/cmb/henriques_2014_z2_cha.txt"
obs_labels = ['Henriques+2014, z = 2.0', 'Gruppioni+2015, z = 2.0-2.5'] # (gsmf observed, sfrf observed)

#test1(obsGSM=GMSfile,obsSFR=SFRfile, colsSFR=[0,1,2,3],colsGSM=[0,1,2,3],
     # labelObs=obs_labels, outplot=plot2file, h0=0.6777,volume = const.vol_pm,verbose=True)


'''TEST PLOT MEDIANS'''
from get_nebular_emission.eml_plots import test_medians as testmed
#plot2file = r"PRUEBAS"
#testmed(outplot=plot2file, verbose=True)

'''TEST BPT DIAGRAM'''
from get_nebular_emission.eml_plots import test_bpt
plot = r"PRUEBAS"
test_bpt(outplot=plot,photmod='gutkin16',verbose=True)








