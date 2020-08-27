Requirements to your input data
*******************************

Compatible image formats
~~~~~~~~~~~~~~~~~~~~~~~~

The input images can have any GDAL compatible image format. You can find a list here:
http://www.gdal.org/formats_list.html


Geocoding
~~~~~~~~~

Your target image must be approximately geocoded to your reference image.
In case of ENVI files, this means they must have a 'map info' and a 'coordinate system string' as attributes of their
header file.

.. note::

    AROSICS also allows to compute the misregistration between two input images without any geocoding. In this case,
    it is assumed that both images have the same spatial resolution and their upper-left coordinates approximately
    represents the same map coordinate. The computed misregistration is then returned in image coordinate units.


Supported projections
~~~~~~~~~~~~~~~~~~~~~

AROSICS provides full support for UTM projections and geographic coordinates. Providing support for other projections
is currently work in progress (see `here <https://gitext.gfz-potsdam.de/danschef/arosics/-/issues/37>`__ for the
status) and may lead to some incompatibities.


Geographic overlap
~~~~~~~~~~~~~~~~~~

The input images must have a geographic overlap but clipping them to same geographical extent is NOT neccessary.
The image overlap area is automatically calculated. Thereby, no-data regions within the images are automatically
respected.


Spatial resolution
~~~~~~~~~~~~~~~~~~

The input images may have different spatial resolutions. Any needed resampling of the data is done automatically.

.. attention::

    Please try avoid any spatial resampling of the input images before running AROSICS. It might affect
    the accuracy of the computed misregistration.


Orthorectified datasets
~~~~~~~~~~~~~~~~~~~~~~~

Please use ortho-rectified input data in order to minimize local shifts in the input images.


No-data values
~~~~~~~~~~~~~~

The no-data value of each image is automatically derived from the image corners. However, this may fail if the actual
no-data value is not present within a 3x3 matrix at the image corners. User provided no-data values will speed up the
computation and avoid wrongly derived values.


Actual image corner coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Providing the map coordinates of the actual image corners lets you save some computation time,
because in this case the implemented automatic algorithm can be skipped.


Image masks / areas to be excluded from tie-point creation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The user may provide user-defined masks to exclude certain image areas from tie point creation. This is useful for
example in case of cloud areas or moving objects. However, the outlier detection algorithms of AROSICS will filter out
tie points with large differences to the surrounding area.


Unequal sensors and image acquisition conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AROSICS is designed to robustly handle the typical difficulties of multi-sensoral/multi-temporal images.
