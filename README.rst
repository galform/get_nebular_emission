|docs|

.. inclusion-marker-do-not-remove

This code is in beta. Please report issues and be patient with the incomplete documentation.
   
Get nebular emission
======================

**get_nebular_emission** is a Python package that given the metallicity of cold gas, the stellar mass and either the specific star formation or the number of ionizing photons, calculates the intensity of nebular emission lines from star forming HII regions. It also calculates the emission lines from narrow-line regions of AGNs from the AGN luminosity and the size and mass of galactic components. 

The current version of the code is based on the Fortran90 one used for the paper Baugh, Lacey, Gonzalez-Perez and Manzoni 2022 (https://arxiv.org/abs/2112.00129).

The code reads galaxy properties from an input text or HDF5 file with values separated by spaces. This can have a header, as long as the header lines start with characters or signs well specify (different from the minus sign). The file can contain several columns with properties for different galactic components. For each component, the column number (starting from 0) should be provided for the following properties: the stellar mass (Msun), the (instantaneous) star formation rate (SFR Msun/Gyr), the metallicity of the cold gas (defined as the ratio MZcold/Mcold).

If a h0 value is specified, then the units will be assumed to be: stellar mass (Msun/h), SFR (Msun/h/Gyr).

Requirements and Installation
-----------------------------

This code has been developed in Python 3.7.1 and it is compatible with Python above 3.5 versions. The code has been tested to run in Linux and Windows operating systems. 

This code uses numpy as specified in docs/requirements.txt.

The code can be run directly from a cloned GitHub `repository`_ and then makeing a local installation from the main directory (where you can see `setup.py`:
:code:`python3 setup.py install`.

Tutorial
-----------------------------

|NII|

.. _pyversion: https://uk.mathworks.com/help/matlab/getting-started-with-python.html

.. _package: https://pypi.org/project/get_nebular_emission/

.. _repository: https://github.com/galform/get_nebular_emission
	     
.. |docs| image:: https://readthedocs.org/projects/get-nebular-emission/badge/?version=latest
   :target: https://get-nebular-emission.readthedocs.io/en/latest/
   :alt: Documentation Status
   
.. |NII| image:: https://i.ibb.co/xSxQr58/NII-test.png



