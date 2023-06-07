|docs|

.. inclusion-marker-do-not-remove

This code is in beta. Please report issues and be patient with the incomplete documentation.
   
Get nebular emission
======================

**get_nebular_emission** is a Python package that given the metallicity of cold gas, the stellar mass and either the specific star formation or the number of ionizing photons, calculates the intensity of nebular emission lines from star forming HII regions.

The current version of the code is based on the Fortran90 one used for the paper Baugh, Lacey, Gonzalez-Perez and Manzoni 2022 (https://arxiv.org/abs/2112.00129).

The code reads galaxy properties from an input text file with values separated by spaces. This can have a header, as long as the header lines start with characters or signs well specify (different from the minus sign). The file can contain several columns with properties for different galactic components. For each component, the column number (starting from 0) should be provided for the following properties: the stellar mass (Msun), the (instantaneous) star formation rate (SFR Msun/Gyr), the metallicity of the cold gas (defined as the ratio MZcold/Mcold).

If a h0 value is specified, then the units will be assumed to be: stellar mass (Msun/h), SFR (Msun/h/Gyr).

Tests
-----

The test can be run using unittest:
:code:`python3 -m unittest discover -s tests`.

Requirements and Installation
-----------------------------

This code has been developed in Python 3.7.1 and it is compatible with Python above 3.5 versions. The code has been tested to run in Linux and Windows operating systems. 

This code uses numpy as specified in docs/requirements.txt.

The code can be run directly from a cloned GitHub `repository`_ and then makeing a local installation from the main directory (where you can see `setup.py`:
:code:`python3 setup.py install`.


The functions in the package can be used after importing novonix_add, for example as follows:



.. _pyversion: https://uk.mathworks.com/help/matlab/getting-started-with-python.html

.. _package: https://pypi.org/project/get_nebular_emission/

.. _repository: https://github.com/galform/get_nebular_emission
	     
.. |docs| image:: https://readthedocs.org/projects/get_nebular_emission/badge/?version=latest
   :target: https://get_nebular_emission.readthedocs.io/en/latest/
   :alt: Documentation Status

ISSUES
------

**eml_photio**

· The interpolations can be made more efficient.

· If the limits file does not exist, allow the program to continue the interpolations.

· Program set to do linear interpolations. See what happends with other types of interpolations.

· Test properly Cardelli's model for dust attenuation.

**eml_plots**

· In test_sfrf set the program to allow several observation data and automate for every redshift.

· In test_sfrf there is a problem with the contours, the plot is not well in my opinion. Set the program to allow automate the levels of the contours. One way to do it: In SFRF search the bin with at least 100 galaxies (before divide the data by the volume) to be the minimum level, get the maximum value of phi to be the maximum level and take another one between these two. 

· Add test_zm. Similar to test_sfrf but with the metallicity.

· In test_sfrf verify the limits for the SFR and the mass. 

· Put some warnings if we have data out of the limits to do the plots. 

· Put more warnings to notify the errors better.

· In the future, for verification, we could compare the mstardot_average versus the mstardot given by GALFORM and the mstardot_average versus the (mstardot + mstardot_burst). We have assumed that mstardot is the SFR averaged of the disc, so this is a way to verify this assumption.

**Others**

· Allow the escape fraction vary: (100-escf)*LyC



