|build| |coverage| |docs| |pypi| |zenodo| 

.. inclusion-marker-do-not-remove

Get nebular emission
======================

**get_nebular_emission** is a Python package that given the metallicity of cold gas, the stellar mass and either the specific star formation or the number of ionizing photons, calculates the intensity of nebular emission lines from star forming HII regions.


Example
-------

The **example.py** runs over the given example data, producing a new file and a plot that compares the original and the prepared data. To run this
example, simply type: :code:`python example.py`.

Requirements and Installation
-----------------------------

This code has been developed in Python 3.7.1 and it is compatible with Python above 3.5 versions. The code has been tested to run in Linux operating systems. 

This code uses numpy as specified in docs/requirements.txt. The ploting routine from the *example.py* also requires the use of matplotlib.

The code can be run directly from a cloned GitHub `repository`_ or it can also be installed as a python `package`_ through pip:

.. code::

   pip install get_nebular_emission

The functions in the package can be used after importing novonix_add, for example as follows:



.. _pyversion: https://uk.mathworks.com/help/matlab/getting-started-with-python.html

.. _package: https://pypi.org/project/preparenovonix/

.. _repository: https://github.com/BatLabLancaster/preparenovonix

.. |build| image:: https://travis-ci.org/BatLabLancaster/preparenovonix.svg?branch=master
    :target: https://travis-ci.org/BatLabLancaster/preparenovonix

.. |coverage| image:: https://codecov.io/gh/BatLabLancaster/preparenovonix/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/BatLabLancaster/preparenovonix
	     
.. |docs| image:: https://readthedocs.org/projects/prepare-novonix-data/badge/?version=latest
   :target: https://prepare-novonix-data.readthedocs.io/en/latest/
   :alt: Documentation Status

.. |pypi| image:: https://img.shields.io/pypi/v/preparenovonix.svg
    :target: https://pypi.org/project/preparenovonix/
	 
.. |zenodo| image:: https://zenodo.org/badge/186994865.svg
   :target: https://zenodo.org/badge/latestdoi/186994865
