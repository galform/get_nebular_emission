import preparenovonix.novonix_prep as prep
from preparenovonix.compare import plot_vct

# The r" " are needed to handle Windows paths
# otherwise ' ' can be enough to include the path.
infile = r"example_data/example_data.csv"

# Prepare the Novonix data file
prep.prepare_novonix(infile, addstate=True, lprotocol=True,
                     overwrite=False, verbose=True)

# Compare the prepared file with the original one
plot_vct(infile, first_loop=0, plot_type='png', plot_show=True)
