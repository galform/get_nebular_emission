|build| |coverage| |docs| |pypi| |zenodo| 

.. inclusion-marker-do-not-remove

Get nebular emission
======================

**get_nebular_emission** is a Python package that given the metallicity of cold gas, the stellar mass and either the specific star formation or the number of ionizing photons, calculates the intensity of nebular emission lines from star forming HII regions.


The code will read galaxy properties from an input text file with values separated by spaces. This can have a header, as long as the header lines start with characters (different from the minus sign). The file can contain several columns with properties for different galactic components. For each component, the column number (starting from 0) should be provided for the following properties: the stellar mass (Msun), the (instantaneous) star formation rate (SFR Msun/Gyr), the metallicity of the cold gas (defined as the ratio MZcold/Mcold).

If a h0 value is specified, then the units will be assumed to be: stellar mass (Msun/h), SFR (Msun/h/Gyr).


TO DO
-----
[] Calculate Bulge and Disk separately and then add the fluxes from the 2 contributions.
[] Adding AGN lines. From Cedric: For Galform, we should try to do this in a way that the spectrum and luminosity of the ionizing radiation assumed in the emission lines calculation is consistent with the AGN SED assumed in Galform. Currently we use an AGN template spectrum from Marconi, for which the shape depends on the AGN bolometric luminosity. Andrew updated the Galform code so that it calculated observer-frame magnitudes for any bands of interest, along with observer-frame luminosities in hard and soft X-ray bands, for outputs at snapshots, based on the Marconi template (which depends on bolometric luminosity). We should check that this is properly included in the main version of Galform on the archive, and that it works correctly. Andrew was still working on calculating AGN luminosities correctly on lightcones (the issue being that AGN luminosities do not vary smoothly between snapshots) when he left, so I don't think there is code for that yet.


Example
-------

The **example.py** runs over the given example data, producing a new file and a plot that compares the original and the prepared data. To run this
example, simply type: :code:`python example.py`.

Tests
-----

The test can be run using unittest:
:code:`python3 -m unittest discover -s tests`.

Requirements and Installation
-----------------------------

This code has been developed in Python 3.7.1 and it is compatible with Python above 3.5 versions. The code has been tested to run in Linux operating systems. 

This code uses numpy as specified in docs/requirements.txt. The ploting routine from the *example.py* also requires the use of matplotlib.

The code can be run directly from a cloned GitHub `repository`_ and then makeing a local installation from the main directory (where you can see `setup.py`:
:code:`python3 setup.py install`.


The functions in the package can be used after importing novonix_add, for example as follows:



.. _pyversion: https://uk.mathworks.com/help/matlab/getting-started-with-python.html

.. _package: https://pypi.org/project/get_nebular_emission/

.. _repository: https://github.com/galform/get_nebular_emission

.. |build| image:: https://travis-ci.org/galform/get_nebular_emission.svg?branch=master
    :target: https://travis-ci.org/galform/get_nebular_emission

.. |coverage| image:: https://codecov.io/gh/galform/get_nebular_emission/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/galform/get_nebular_emission
	     
.. |docs| image:: https://readthedocs.org/projects/get_nebular_emission/badge/?version=latest
   :target: https://get_nebular_emission.readthedocs.io/en/latest/
   :alt: Documentation Status

.. |pypi| image:: https://img.shields.io/pypi/v/get_nebular_emission.svg
    :target: https://pypi.org/project/get_nebular_emissioin/
	 
.. |zenodo| image:: https://zenodo.org/badge/186994865.svg
   :target: https://zenodo.org/badge/latestdoi/186994865
