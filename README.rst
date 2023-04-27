.. figure:: https://danschef.git-pages.gfz-potsdam.de/arosics/images/arosics_logo.png
    :target: https://git.gfz-potsdam.de/danschef/arosics
    :align: center

==================================================================================================
An Automated and Robust Open-Source Image Co-Registration Software for Multi-Sensor Satellite Data
==================================================================================================

* Free software: Apache 2.0
* **Documentation:** https://danschef.git-pages.gfz-potsdam.de/arosics/doc/
* The (open-access) **paper** corresponding to this software repository can be found here:
  `Scheffler et al. 2017 <https://www.mdpi.com/2072-4292/9/7/676>`__
  (cite as: Scheffler D, Hollstein A, Diedrich H, Segl K, Hostert P. AROSICS: An Automated and Robust Open-Source
  Image Co-Registration Software for Multi-Sensor Satellite Data. Remote Sensing. 2017; 9(7):676).
* Information on how to **cite the AROSICS Python package** can be found in the
  `CITATION <https://git.gfz-potsdam.de/danschef/arosics/-/blob/main/CITATION>`__ file.
* Submit feedback by filing an issue `here <https://git.gfz-potsdam.de/danschef/arosics/issues>`__
  or join our chat here: |Gitter|

.. |Gitter| image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/arosics/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link
    :alt: https://gitter.im/arosics/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link

Status
------

.. image:: https://git.gfz-potsdam.de/danschef/arosics/badges/main/pipeline.svg
        :target: https://git.gfz-potsdam.de/danschef/arosics/commits/main
.. image:: https://git.gfz-potsdam.de/danschef/arosics/badges/main/coverage.svg
        :target: https://danschef.git-pages.gfz-potsdam.de/arosics/coverage/
.. image:: https://img.shields.io/pypi/v/arosics.svg
        :target: https://pypi.python.org/pypi/arosics
.. image:: https://img.shields.io/conda/vn/conda-forge/arosics.svg
        :target: https://anaconda.org/conda-forge/arosics
.. image:: https://img.shields.io/pypi/l/arosics.svg
        :target: https://git.gfz-potsdam.de/danschef/arosics/blob/main/LICENSE
.. image:: https://img.shields.io/pypi/pyversions/arosics.svg
        :target: https://img.shields.io/pypi/pyversions/arosics.svg
.. image:: https://img.shields.io/pypi/dm/arosics.svg
        :target: https://pypi.python.org/pypi/arosics
.. image:: https://zenodo.org/badge/253474603.svg
        :target: https://zenodo.org/badge/latestdoi/253474603

See also the latest coverage_ report and the pytest_ HTML report.

Feature overview
----------------

AROSICS is a python package to perform **automatic subpixel co-registration** of two satellite image datasets
based on an image matching approach working in the frequency domain, combined with a multistage workflow for
effective detection of false-positives.

It detects and corrects **local as well as global misregistrations** between two input images in the subpixel scale,
that are often present in satellite imagery. The algorithm is robust against the typical difficulties of
multi-sensoral/multi-temporal images. Clouds are automatically handled by the implemented outlier detection algorithms.
The user may provide user-defined masks to exclude certain image areas from tie point creation. The image overlap area
is automatically detected. AROSICS supports a wide range of input data formats and can be used from the command
line (without any Python experience) or as a normal Python package.


Global co-registration - fast but only for static X/Y-shifts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Only a global X/Y translation is computed within a small subset of the input images (window position is adjustable).
This allows very fast co-registration but only corrects for translational (global) X/Y shifts.
The calculated subpixel-shifts are (by default) applied to the geocoding information of the output image.
No spatial resampling is done automatically as long as both input images have the same projection. However, AROSICS
also allows to align the output image to the reference image coordinate grid if needed.

Here is an example of a Landsat-8 / Sentinel-2 image pair before and after co-registration using AROSICS:

.. image:: https://git.gfz-potsdam.de/danschef/arosics/raw/main/docs/images/animation_testcase1_zoom_L8_S2_global_coreg_before_after_900x456.gif


Local co-registration - for spatially variable shifts but a bit slower
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A dense grid of tie points is automatically computed, whereas tie points are subsequently validated using a
multistage workflow. Only those tie points not marked as false-positives are used to compute the parameters of an
affine transformation. Warping of the target image is done using an appropriate resampling technique
(cubic by default).

Here is an example of the computed shift vectors after filtering false-positives
(mainly due to clouds in the target image):

.. image:: https://git.gfz-potsdam.de/danschef/arosics/raw/main/docs/images/shift_vectors_testcase1__900x824.gif


For further details check out the `documentation <https://danschef.git-pages.gfz-potsdam.de/arosics/doc/>`__!


History / Changelog
-------------------

You can find the protocol of recent changes in the AROSICS package
`here <https://git.gfz-potsdam.de/danschef/arosics/-/blob/main/HISTORY.rst>`__.


Credits
-------

AROSICS was developed by Daniel Scheffler (German Research Centre of Geosciences) within the context of the
`GeoMultiSens <http://www.geomultisens.gfz-potsdam.de/>`__ project funded by the German Federal Ministry of Education and Research
(project grant code: 01 IS 14 010 A-C).

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.
The test data represent modified Copernicus Sentinel-2 data (ESA 2016). The input data for the figures in the
documentation have been provided by NASA (Landsat-8) and ESA (Sentinel-2).

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _coverage: https://danschef.git-pages.gfz-potsdam.de/arosics/coverage/
.. _pytest: https://danschef.git-pages.gfz-potsdam.de/arosics/test_reports/report.html
.. _conda: https://docs.conda.io/

