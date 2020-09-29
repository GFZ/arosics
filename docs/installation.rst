============
Installation
============

Using conda_, AROSICS is installed as follows:


1. Create virtual environment for arosics (optional but recommended):

   .. code-block:: bash

    $ conda create -c conda-forge --name arosics python=3
    $ source activate arosics


2. Then install AROSICS itself:

   .. code-block:: bash

    $ conda install -c conda-forge arosics


This is the preferred method to install AROSICS, as it will always installs the most recent stable release and
automatically resolves all the dependencies.


If you don't use conda_, you may use the pip installer of AROSICS. However, please note that AROSICS depends on some
open source packages that may cause problems when installed with pip. Therefore, we strongly recommend resolving the
`dependencies of AROSICS`_ before the automatic installer is run.


To run the pip installer of AROSICS use:

   .. code-block:: bash

    $ pip install arosics


.. note::

    AROSICS has been tested with Python 3.4+ and Python 2.7. It should be fully compatible to all Python versions
    from 2.7 onwards. However, we will continously drop the support for Python 2.7 in future.


If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _conda: https://conda.io/docs
.. _`dependencies of AROSICS`: https://gitext.gfz-potsdam.de/danschef/arosics/-/blob/master/requirements.txt
