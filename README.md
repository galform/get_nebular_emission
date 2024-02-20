|docs|

.. inclusion-marker-do-not-remove

This code is in beta. Please report issues and be patient with the incomplete documentation.
   
Get nebular emission
======================

**get_nebular_emission** is a Python package that given the metallicity of cold gas, the stellar mass and either the specific star formation or the number of ionizing photons, calculates the intensity of nebular emission lines from star forming HII regions. It also calculates the emission lines from narrow-line regions of AGNs from the bolometric luminosity (or properties to derive this, such as the mass accretion rate), the total mass of cold gas and the effective radius of the galactic components. 

The code expects text or HDF5 files with data of global galactic properties from one or more galaxies. These properties can correspond to any number of galaxy component (for example, disk and bulge), for which the code calculates emission lines independently. Finally, it outputs a HDF5 file with all the results. As long as the needed input parameters are available, the code can be coupled with model galaxies produced with different codes: hydrodynamical simulations, semi-analytical models, empirical ones, etc. Although the current code uses global galactic properties (or bulge and disc ones) as input, it can be easily adapted to use galactic regions or particles as input.  

The code considers two possible sources of emission: star-forming or HII regions [Baugh+2022](https://arxiv.org/abs/2112.00129) and the narrow-line regions (NLR) of AGNs. Although the code structure is the same for both ionising origins, different inputs are required. First, global galactic properties are connected to the ionising source, either HII regions or the NLR. Then, a photoionisation model is used to obtain the emission line luminosities given the characteristics of the ionising sources. This is done interpolating the values from predefined photoionisation grids. Each part of the procedure is handled by a different module, specified in the flowchart:

|flowchart|

Requirements and Installation
-----------------------------

This code has been developed in Python 3.7.1 and it is compatible with Python above 3.5 versions. The code has been tested to run in Linux and Windows operating systems. 

This code uses numpy as specified in docs/requirements.txt.

The code can be run directly from a cloned GitHub repository and then run following the **run_code_tutorial.py** as a normal python3 program.  


Tutorial and testing
-----------------------------

A brief tutorial going through all the main options of the code is available at `run_code_tutorial.py`. A more in-depth review of all the options in the code can be found at Expósito-Márquez et. al. 2024 in prep. The tutorial code is prepared to run over model galaxies produced at z = 0 by the semi-analytical model described in [Gonzalez-Perez et al. 2020](https://academic.oup.com/mnras/article/498/2/1852/5894931), stored in the *example_data* directory. The repository also has a program to plot line-ratio diagrams from the output files of the code, with several predefined selection criteria. 

A BPT with the results for the example file without considering attenuation and with the selection criteria for local galaxies from [Favole et al. 2023](https://arxiv.org/abs/2303.11031) can be seen below for testing purposes.


|NII|

.. _pyversion: https://uk.mathworks.com/help/matlab/getting-started-with-python.html

.. _repository: https://github.com/galform/get_nebular_emission
	     
.. |docs| image:: https://readthedocs.org/projects/get-nebular-emission/badge/?version=latest
   :target: https://get-nebular-emission.readthedocs.io/en/latest/
   :alt: Documentation Status
   
.. |NII| image:: https://i.ibb.co/xSxQr58/NII-test.png

.. |flowchart| image:: https://i.ibb.co/CsdZjgm/flow-chart.png



