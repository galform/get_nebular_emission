import numpy as np
import get_nebular_emission.eml_const as const

'''TEST PLOT SFRF'''
def plot_sfrf(inputdata=[r'output_data/emlines_GP20_z0.0_Lagn.hdf5'],
              plot2file = r'plots/sfrf.pdf',
               SFRfile = r'sfrf/gruppioni_2015_z0.0-0.3_cha.txt',
              GMSfile = r'gsmf/henriques_2014_z0_cha.txt',
              obs_labels = ['Henriques+2014, z = 0', 'Gruppioni+2015, z = 0.0-0.3'],
              volume = 5 * 125**3, #1*(542.16**3)/200,
              specific=False):
    
    from get_nebular_emission.eml_plots import test_sfrf
    
    test_sfrf(inputdata=inputdata, obsGSM=GMSfile,obsSFR=SFRfile,colsSFR=[0,1,2,3],colsGSM=[0,1,2,3],
          labelObs=obs_labels, outplot=plot2file, specific=specific, h0=const.h,volume = volume,verbose=True)


'''TEST PLOT MEDIANS'''
def plot_medians(infile=r"output_data/emlines_GALFORM_p_e0.001.hdf5",plot2folder=r'./plots'):
    from get_nebular_emission.eml_plots import test_medians as testmed
    testmed(infile=infile,outplot=plot2folder, verbose=True)

'''TEST BPT DIAGRAM'''
def plot_bpt(infile=["output_data/emlines_GALFORM_k.hdf5"],plot2folder=r'./plots'):
    from get_nebular_emission.eml_plots import test_bpt
    test_bpt(infile=infile, outplot=plot2folder,photmod='gutkin16',plot_phot=False,create_file=False,verbose=True)