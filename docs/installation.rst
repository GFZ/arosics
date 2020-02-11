============
Installation
============

AROSICS depends on some open source packages which are usually installed without problems by the automatic install
routine. However, for some projects, we strongly recommend resolving the dependency before the automatic installer
is run. This approach avoids problems with conflicting versions of the same software.
Using conda_, the recommended approach is:


1. Create virtual environment for arosics (optional but recommended):

   .. code-block:: bash

    $ conda create --name arosics python=3
    $ source activate arosics


2. Install some libraries needed for AROSICS:

   .. code-block:: bash

    $ conda install -c conda-forge numpy gdal scikit-image matplotlib pyproj "shapely<=1.6.4" geopandas pandas cmocean


3. Install optional libraries for AROSICS (only needed for some specific functions):

   .. code-block:: bash

    $ conda install -c conda-forge basemap pykrige
    $ conda install -c conda-forge pyfftw  # Linux and MacOS
    $ conda install -c jesserobertson pyfftw  # Windows


4. Then install AROSICS using the pip installer:

   .. code-block:: bash

    $ pip install arosics


This is the preferred method to install arosics, as it will always install the most recent stable release.

.. note::

    AROSICS has been tested with Python 3.4+ and Python 2.7. It should be fully compatible to all Python versions
    from 2.7 onwards. However, we will continously drop the support for Python 2.7 in future.


If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _conda: https://conda.io/docs
