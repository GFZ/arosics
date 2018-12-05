.. highlight:: shell

============
Installation
============

AROSICS depends on some open source packages which are usually installed without problems by the automatic install
routine. However, for some projects, we strongly recommend resolving the dependency before the automatic installer
is run. This approach avoids problems with conflicting versions of the same software.
Using conda_, the recommended approach is:

 .. code-block:: console

    # create virtual environment for arosics, this is optional
    conda create --name arosics python=3
    source activate arosics
    conda install -c conda-forge numpy gdal scikit-image matplotlib pyproj rasterio shapely geopandas cmocean

    # optional libraries:
    conda install -c conda-forge basemap pykrige
    conda install -c conda-forge pyfftw  # Linux and MacOS
    conda install -c jesserobertson pyfftw  # Windows


Stable release
--------------

To install AROSICS, use the pip installer:

 .. code-block:: console

    pip install arosics


Or clone the repository via GIT and update the PATH environment variable:

 .. code-block:: console

    cd /your/installation/folder
    git clone https://gitext.gfz-potsdam.de/danschef/arosics.git
    git clone https://gitext.gfz-potsdam.de/danschef/geoarray.git
    git clone https://gitext.gfz-potsdam.de/danschef/py_tools_ds.git
    PATH=$PATH:/path/to/your/installation/folder/arosics:/path/to/your/installation/folder/geoarray:/path/to/your/installation/folder/py_tools_ds


This is the preferred method to install arosics, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for arosics can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/danschef/arosics

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://gitext.gfz-potsdam.de/danschef/arosics/repository/archive.tar.gz?ref=master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://gitext.gfz-potsdam.de/danschef/arosics
.. _tarball: https://gitext.gfz-potsdam.de/danschef/arosics/repository/archive.tar.gz?ref=master
.. _conda: https://conda.io/docs/
