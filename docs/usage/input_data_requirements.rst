Requirements to your input data
*******************************

Compatible image formats
~~~~~~~~~~~~~~~~~~~~~~~~

The input images can have any GDAL compatible image format. You can find a list here:
https://gdal.org/drivers/raster/index.html


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

AROSICS was initially written with support for UTM and geographic coordinates only. Full support for any other
projection was added in version 1.4.0. However, make sure your input images have the same projection. Different
projections for the reference and target image are currently not supported.

AROSICS can also be applied to images without any projection and geocoding information. In this case, however,
the input images need to have the same pixel size and must cover more or less the same spatial area
(with a shift a few pixels at most).


Geographic overlap
~~~~~~~~~~~~~~~~~~

The input images must have a geographic overlap but clipping them to same geographical extent is NOT neccessary
The image overlap area is automatically calculated, given that your input images have valid geocoding and projection
information). Thereby, no-data regions within the images are automatically respected.


Spatial resolution
~~~~~~~~~~~~~~~~~~

The input images may have different spatial resolutions. Any needed resampling of the data is done automatically.

.. attention::

    Please try to avoid any spatial resampling of the input images before running AROSICS. It might affect
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
