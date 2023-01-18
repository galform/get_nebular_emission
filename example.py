import get_nebular_emission.eml as eml
import time
import glob
#from preparenovonix.compare import plot_vct

# The r" " are needed to handle Windows paths
# otherwise ' ' can be enough to include the path.
infile = glob.glob(r"example_data/emlines_lc16_PMILL_iz245_ivol8.dat")
infile = glob.glob(r"example_data/SAGE_z0.142_*.hdf5")
outfile = r"output_data/emlines_SAGE_test.hdf5"

# Calculate the nebular emission lines
# The columns where the stellar mass, SFR and Zcold=Mzcold/Mcold are
# to be found in the input file, should be specified within the argument m_sfr_z
# Several components are allowed (e.g. disk and bulge).

inputformat = 'HDF5'
LC2sfr = False
if LC2sfr:
    if inputformat=='HDF5':
        cols=[['mstars_total','mag_LC_r_disk','zcold'],['mstars_bulge','mag_LC_r_bulge','zcold_burst']]
    else:
        cols = [[11,0,2],[15,1,4]]
else:
    if inputformat=='HDF5':
        cols=[['mstars_total','mstardot','zcold'],['mstars_bulge','mstardot_burst','zcold_burst']]
    else:
        cols = [[11,13,2],[15,14,4]]
    
IMF = ['Kennicut','Top-heavy'] 
# Kennicut, Salpeter, Kroupa, Chabrier, Baldry&Glazebrook, Top-heavy

cols=[['StellarMass','SfrDisk','zcold'],['BulgeMass','SfrBulge','zcold']]
att_param=['Rvir','ColdGas']

start_time = time.perf_counter()

eml.eml(infile, outfile, m_sfr_z=cols, att_param=att_param, h0=0.6777, 
                  IMF=IMF, inputformat=inputformat,LC2sfr=LC2sfr, 
                  Zloh12=False, mtot2mdisk=False, 
                  attmod='cardelli89',unemod='kashino20',photmod='gutkin16',
                  verbose=True, Plotting=True, Testing=False)

print('Total time: ', round(time.perf_counter() - start_time,2), 's.')

