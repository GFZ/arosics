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
.. image:: https://img.shields.io/pypi/l/arosics.svg
        :target: https://gitext.gfz-potsdam.de/danschef/arosics/blob/master/LICENSE
.. image:: https://img.shields.io/pypi/pyversions/arosics.svg
        :target: https://img.shields.io/pypi/pyversions/arosics.svg

See also the latest coverage_ report and the nosetests_ HTML report.


Feature overview
----------------

Detection and correction of local or global geometric displacements between two input images:

* Local co-registration:

  A dense grid of tie points is automatically computed, whereas tie points are subsequently validated using a
  multistage workflow. Only those tie points not marked as false-positives are used to compute the parameters of an
  affine transformation. Warping of the target image is done using an appropriate resampling technique
  (cubic by default).

|

* Global co-registration:

  Only a global X/Y translation is computed within a small subset of the input images (window position is adjustable).
  This allows very fast co-registration but only corrects for translational (global) X/Y shifts.
  The calculated subpixel-shifts are (by default) applied to the geocoding information of the output image.
  No spatial resampling is done automatically as long as both input images have the same projection.
  If you need the output image to be aligned to the reference image coordinate grid
  (by using an appropriate resampling algorithm), use the '-align_grids' option.


Credits
-------

AROSICS was developed within the context of the GeoMultiSens project funded by the German Federal Ministry
of Education and Research (project grant code: 01 IS 14 010 A-C).

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.
The test data represent modified Copernicus Sentinel data (2016).

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _coverage: http://danschef.gitext.gfz-potsdam.de/arosics/coverage/
.. _nosetests: http://danschef.gitext.gfz-potsdam.de/arosics/nosetests_reports/nosetests.html
.. _conda: https://conda.io/docs/

