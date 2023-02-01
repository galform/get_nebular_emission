import get_nebular_emission.eml as eml
import get_nebular_emission.eml_const as const
import time
import glob

# The r" " are needed to handle Windows paths
# otherwise ' ' can be enough to include the path.
infile = glob.glob(r"example_data/emlines_lc16_PMILL_iz245_ivol*.hdf5")
# infile = glob.glob(r"example_data/SAGE_z0.142_*.hdf5")

outfile = r"output_data/emlines_GALFORM.hdf5"

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
# Mass: Msun/h
# SFR: Msun/(yr*h)
# Zcold: MZcold/Mcold

IMF_i = ['Kennicut','Kennicut'] # Input IMF
IMF_f = ['Kennicut','Top-heavy'] # IMF to which transform each component
# Kennicut, Salpeter, Kroupa, Chabrier, Baldry&Glazebrook, Top-heavy
# Kashino asumes Kroupa.

if 'SAGE' in infile[0]:
    cols=[['StellarMass','SfrDisk','zcold'],['BulgeMass','SfrBulge','zcold']]
att_param=['Rvir','ColdGas']

if 'emlines' in infile[0]:
    if inputformat=='HDF5':
        cutcols = ['mag_SDSS_r_o_t']
    else:
        cutcols = [19]
    mincuts = [None]
    maxcuts = [19.8 - 38.034]
else:
    cutcols = ['Mvir']
    mincuts = [20*1.25*1e9]
    maxcuts = [None]
    
eml.eml(infile, outfile, m_sfr_z=cols, cutcols=cutcols, mincuts=mincuts, 
                  maxcuts=maxcuts, att_param=att_param, h0=const.h, 
                  IMF_i=IMF_i, IMF_f=IMF_f, inputformat=inputformat,
                  cutlimits=False, mtot2mdisk=True, LC2sfr=LC2sfr,
                  attmod='ratios', unemod='kashino20', photmod='gutkin16',
                  verbose=True, Plotting=True, Testing=False)

