import get_nebular_emission.eml as eml
#from preparenovonix.compare import plot_vct

# The r" " are needed to handle Windows paths
# otherwise ' ' can be enough to include the path.
infile = r"example_data/example_data.dat"

# Calculate the nebular emission lines
eml.eml(infile, m_sfr_z=[[13,12,4]],h0=0.6777, verbose=True, Testing=True)
#eml.eml(infile, m_sfr_z=[[13,12,4],[13,12,4]],h0=0.6777, verbose=True, Testing=True)

## Compare the prepared file with the original one
#plot_vct(infile, first_loop=0, plot_type='png', plot_show=True)
