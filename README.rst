|docs|

.. inclusion-marker-do-not-remove

This code is in beta. Please report issues and be patient with the incomplete documentation.
   
Get nebular emission
======================

**get_nebular_emission** is a Python package that given the metallicity of cold gas, the stellar mass and either the specific star formation or the number of ionizing photons, calculates the intensity of nebular emission lines from star forming HII regions. It also calculates the emission lines from narrow-line regions of AGNs from the bolometric luminosity (or properties to derive this, such as the mass accretion rate), the total mass of cold gas and the effective radius of the galactic components. 

The code expects text or HDF5 files with data of global galactic properties from one or more galaxies. These properties can correspond to any number of galaxy component (for example, disk and bulge), for which the code calculates emission lines independently. Finally, it outputs a HDF5 file with all the results. As long as the needed input parameters are available, the code can be coupled with model galaxies produced with different codes: hydrodynamical simulations, semi-analytical models, empirical ones, etc. Although the current code uses global galactic properties (or bulge and disc ones) as input, it can be easily adapted to use galactic regions or particles as input.  

The code considers two possible sources of emission: star-forming or HII regions [Baugh et al. 2022](https://arxiv.org/abs/2112.00129) and the narrow-line regions (NLR) of AGNs. Although the code structure is the same for both ionising origins, different inputs are required. First, global galactic properties are connected to the ionising source, either HII regions or the NLR. Then, a photoionisation model is used to obtain the emission line luminosities given the characteristics of the ionising sources. This is done interpolating the values from predefined photoionisation grids. Each part of the procedure is handled by a different module, specified in the flowchart:

|flowchart|

Requirements and Installation
-----------------------------

This code has been developed in Python 3.7.1 and it is compatible with Python above 3.5 versions. The code has been tested to run in Linux and Windows operating systems. 

This code uses numpy as specified in docs/requirements.txt.

To install this code, clone the GitHub repository:

```
git clone git@github.com:galform/get_nebular_emission.git
```

Running the code
-----------------------------
The code can be run following the **run_txtinput_tutorial.py** or **run_hdf5input_tutorial.py** as a normal python3 programs.  
```
python3 run_txtinput_tutorial.py
```


Tutorials and testing
-----------------------------

Brief tutorial going through all the main options of the code are available for a input text file, **run_txtinput_tutorial.py**. A more in-depth review of all the options in the code can be found at Expósito-Márquez et. al. 2024 in prep. The tutorials run over model galaxies stored in the *example_data* directory. These are produced at z = 0 by the semi-analytical model described in [Gonzalez-Perez et al. 2020](https://academic.oup.com/mnras/article/498/2/1852/5894931).

The example code make selections within the data, calculate emission line luminosities and make test plots. A line ratio diagram with the results from the example text file, without considering dust attenuation, can be seen below together with the local data from [Favole et al. 2023](https://arxiv.org/abs/2303.11031).


|NII|
|SII|


Citing
-----------------------
As you use get_emission_lines, cite at least the first of the following papers:

* [Expósito-Márquez et al 2024a](in prep)
* [Expósito-Márquez et al 2024b](in prep)
* [Baugh et al. 2022](https://arxiv.org/abs/2112.00129)
  ```
  @ARTICLE{2022MNRAS.510.1880B,
       author = {{Baugh}, C.~M. and {Lacey}, Cedric G. and {Gonzalez-Perez}, Violeta and {Manzoni}, Giorgio},
        title = "{Modelling emission lines in star-forming galaxies}",
      journal = {\mnras},
     keywords = {methods: numerical, H II regions, galaxies: formation, Astrophysics - Astrophysics of Galaxies},
         year = 2022,
        month = feb,
       volume = {510},
       number = {2},
        pages = {1880-1893},
          doi = {10.1093/mnras/stab3506},
archivePrefix = {arXiv},
       eprint = {2112.00129},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.1880B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
  ```	

.. _pyversion: https://uk.mathworks.com/help/matlab/getting-started-with-python.html

.. _repository: https://github.com/galform/get_nebular_emission
	     
.. |docs| image:: https://readthedocs.org/projects/get-nebular-emission/badge/?version=latest
   :target: https://get-nebular-emission.readthedocs.io/en/latest/
   :alt: Documentation Status
   
.. |NII| image:: src/example_data/NIIbpt_GP20_62.5kpc_z0_example.pdf

.. |SII| image:: src/example_data/SIIbpt_GP20_62.5kpc_z0_example.pdf
		 
.. |flowchart| image:: https://i.ibb.co/CsdZjgm/flow-chart.png



