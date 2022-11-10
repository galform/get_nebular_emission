import numpy as np
import get_nebular_emission.eml_const as const


'''REDUCED FILE'''
from get_nebular_emission.eml_io import get_reducedfile as reducedfile

# infile = r"example_data/emlines_lc16_PMILL_iz245_ivol0.dat"
# outfile = r"example_data/emlines_lc16_PMILL_iz245_ivol0_reduced.dat"

# # infile = r"example_data/example_data.dat"
# # outfile = r"example_data/example_data_reduced.dat"


# reducedfile(infile,outfile,indcol=[11,15,13,14,0,1,2,4],verbose=True)

'''TEST PLOT SFRF'''
def plot_sfrf(plot2file = r'pruebaplotfeb.pdf',
              SFRfile = r'sfrf/gruppioni_2015_z2.0-2.5_cha.txt',
              GMSfile = r'gsmf/henriques_2014_z2_cha.txt',
              obs_labels = ['Henriques+2014, z = 2.0', 'Gruppioni+2015, z = 2.0-2.5']): # (gsmf observed, sfrf observed)
    
    from get_nebular_emission.eml_plots import test_sfrf
    
    test_sfrf(obsGSM=GMSfile,obsSFR=SFRfile, colsSFR=[0,1,2,3],colsGSM=[0,1,2,3],
          labelObs=obs_labels, outplot=plot2file, h0=0.6777,volume = const.vol_pm,verbose=True)
    
    


'''TEST PLOT MEDIANS'''
def plot_medians(plot2file=r'.'):
    from get_nebular_emission.eml_plots import test_medians as testmed
    testmed(outplot=plot2file, verbose=True)

'''TEST BPT DIAGRAM'''
def plot_bpt(plot2file=r'.'):
    from get_nebular_emission.eml_plots import test_bpt
    test_bpt(outplot=plot2file,photmod='gutkin16',verbose=True)







