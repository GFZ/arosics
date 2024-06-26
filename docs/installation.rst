============
Installation
============

Using Anaconda or Miniconda (recommended)
-----------------------------------------

Using conda_ (latest version recommended), AROSICS is installed as follows:


1. Create virtual environment for arosics (optional but recommended):

   .. code-block:: bash

    $ conda create -c conda-forge --name arosics python=3
    $ conda activate arosics


2. Then install AROSICS itself:

   .. code-block:: bash

    $ conda install -c conda-forge 'arosics>=1.3.0'


This is the preferred method to install AROSICS, as it always installs the most recent stable release and
automatically resolves all the dependencies.


Using pip (not recommended)
---------------------------

There is also a `pip`_ installer for AROSICS. However, please note that AROSICS depends on some
open source packages that may cause problems when installed with pip. Therefore, we strongly recommend
to resolve the following dependencies before the pip installer is run:

    * cartopy
    * gdal
    * geopandas
    * joblib >=1.3.0
    * matplotlib
    * numpy
    * pandas
    * pykrige
    * pyproj >2.2.0
    * scikit-image >=0.21.0
    * shapely

.. note::

    The gdal library must be installed before numpy, otherwise do a
    `re-installation <https://gis.stackexchange.com/a/465888/140483>`_
    in order to respect the build isolation.


Then, the pip installer can be run by:

   .. code-block:: bash

    $ pip install arosics

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.



.. note::

    AROSICS has been tested with Python 3.8+. It should be fully compatible to all Python versions
    from 3.8 onwards. Python 2.7 support was dropped in AROSICS 1.3 due to its end of life status.


.. _pip: https://pip.pypa.io
.. _Python installation guide: https://docs.python-guide.org/en/latest/starting/installation/
.. _conda: https://docs.conda.io
