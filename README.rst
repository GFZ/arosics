.. figure:: http://danschef.gitext.gfz-potsdam.de/arosics/images/arosics_logo.png
        :target: https://gitext.gfz-potsdam.de/danschef/arosics

An Automated and Robust Open-Source Image Co-Registration Software for Multi-Sensor Satellite Data


* Free software: GNU General Public License v3
* Documentation: http://danschef.gitext.gfz-potsdam.de/arosics/doc/


Status
------

.. .. image:: https://img.shields.io/travis/danschef/arosics.svg
        :target: https://travis-ci.org/danschef/arosics

.. .. image:: https://readthedocs.org/projects/arosics/badge/?version=latest
        :target: https://arosics.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. .. image:: https://pyup.io/repos/github/danschef/arosics/shield.svg
     :target: https://pyup.io/repos/github/danschef/arosics/
     :alt: Updates


.. image:: https://gitext.gfz-potsdam.de/danschef/arosics/badges/master/build.svg
        :target: https://gitext.gfz-potsdam.de/danschef/arosics/commits/master
.. image:: https://gitext.gfz-potsdam.de/danschef/arosics/badges/master/coverage.svg
        :target: http://danschef.gitext.gfz-potsdam.de/arosics/coverage/
.. image:: https://img.shields.io/pypi/v/arosics.svg
        :target: https://pypi.python.org/pypi/arosics

See also the latest coverage_ report and the nosetests_ HTML report.


Features
--------

* Detection and correction of local or global geometric displacements between two input images.


Installation
------------

AROSICS depends on some open source packages which are usually installed without problems by the automatic install
routine. However, for some projects, we strongly recommend resolving the dependency before the automatic installer
is run. This approach avoids problems with conflicting versions of the same software.
Using conda_, the recommended approach is:

 .. code-block:: console

    # create virtual environment for arosics, this is optional
    conda create -y -q --name arosics python=3
    source activate arosics
    conda install -y -q -c conda-forge numpy gdal scikit-image matplotlib pyproj rasterio shapely
    conda install -y -q -c conda-forge pyfftw basemap pykrige  # these libraries are optional


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

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.
The test data represent modified Copernicus Sentinel data (2016).

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _coverage: http://danschef.gitext.gfz-potsdam.de/arosics/coverage/
.. _nosetests: http://danschef.gitext.gfz-potsdam.de/arosics/nosetests_reports/nosetests.html
.. _conda: https://conda.io/docs/

