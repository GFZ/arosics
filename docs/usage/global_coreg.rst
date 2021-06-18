Global image co-registration
****************************


Use the class :class:`arosics.COREG` to detect and correct global spatial shifts between a reference and a target image.
It computes a global X-/Y-shift based on a single matching window with customizable position.


Using the Python API
--------------------

calculate spatial shifts - with input data on disk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from arosics import COREG

    >>> im_reference = '/path/to/your/ref_image.bsq'
    >>> im_target    = '/path/to/your/tgt_image.bsq'

    >>> CR = COREG(im_reference, im_target, wp=(354223, 5805559), ws=(256,256))
    >>> CR.calculate_spatial_shifts()

    Calculating actual data corner coordinates for reference image...
    Corner coordinates of reference image:
        [[319090.0, 5790510.0], [351800.0, 5899940.0], [409790.0, 5900040.0], [409790.0, 5790250.0], [319090.0, 5790250.0]]
    Calculating actual data corner coordinates for image to be shifted...
    Corner coordinates of image to be shifted:
        [[319460.0, 5790510.0], [352270.0, 5900040.0], [409790.0, 5900040.0], [409790.0, 5790250.0], [319460.0, 5790250.0]]
    Matching window position (X,Y): 354223/5805559
    Detected integer shifts (X/Y):       0/-2
    Detected subpixel shifts (X/Y):      0.357885632465/0.433837319984
    Calculated total shifts in fft pixel units (X/Y):         0.357885632465/-1.56616268002
    Calculated total shifts in reference pixel units (X/Y):   0.357885632465/-1.56616268002
    Calculated total shifts in target pixel units (X/Y):      0.357885632465/-1.56616268002
    Calculated map shifts (X,Y):				  3.578856324660592 15.661626799963415
    Original map info: ['UTM', 1, 1, 300000.0, 5900040.0, 10.0, 10.0, 33, 'North', 'WGS-84']
    Updated map info:  ['UTM', 1, 1, '300003.57885632466', '5900055.6616268', 10.0, 10.0, 33, 'North', 'WGS-84']


calculate spatial shifts - without any disk access
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, create some example input images for AROSICS in-memory
(we use instances of :class:`GeoArray<geoarray.GeoArray>` for that):

.. code-block:: python

    >>> from geoarray import GeoArray
    >>> from arosics import COREG

    >>> im_reference = '/path/to/your/ref_image.bsq'
    >>> im_target    = '/path/to/your/tgt_image.bsq'

    # get a sample numpy array with corresponding geoinformation as reference image
    >>> geoArr  = GeoArray(im_reference)

    >>> ref_ndarray = geoArr[:]            # numpy.ndarray with shape (10980, 10980)
    >>> ref_gt      = geoArr.geotransform  # GDAL geotransform: (300000.0, 10.0, 0.0, 5900040.0, 0.0, -10.0)
    >>> ref_prj     = geoArr.projection    # projection as WKT string ('PROJCS["WGS 84 / UTM zone 33N....')

    # get a sample numpy array with corresponding geoinformation as target image
    >>> geoArr  = GeoArray(im_target)

    >>> tgt_ndarray = geoArr[:]            # numpy.ndarray with shape (10980, 10980)
    >>> tgt_gt      = geoArr.geotransform  # GDAL geotransform: (300000.0, 10.0, 0.0, 5900040.0, 0.0, -10.0)
    >>> tgt_prj     = geoArr.projection    # projection as WKT string ('PROJCS["WGS 84 / UTM zone 33N....')

    # create in-memory instances of GeoArray from the numpy array data, the GDAL geotransform tuple and the WKT
    # projection string
    >>> geoArr_reference = GeoArray(ref_ndarray, ref_gt, ref_prj)
    >>> geoArr_target    = GeoArray(tgt_ndarray, tgt_gt, tgt_prj)


Now pass these in-memory :class:`GeoArray<geoarray.GeoArray>` instances to :class:`arosics.COREG`
and calculate spatial shifts:

.. code-block:: python

    >>> CR = COREG(geoArr_reference, geoArr_target, wp=(354223, 5805559), ws=(256,256))
    >>> CR.calculate_spatial_shifts()

    Calculating actual data corner coordinates for reference image...
    Corner coordinates of reference image:
        [[300000.0, 5848140.0], [409790.0, 5848140.0], [409790.0, 5790250.0], [300000.0, 5790250.0]]
    Calculating actual data corner coordinates for image to be shifted...
    Corner coordinates of image to be shifted:
        [[300000.0, 5847770.0], [409790.0, 5847770.0], [409790.0, 5790250.0], [300000.0, 5790250.0]]
    Matching window position (X,Y): 354223/5805559
    Detected integer shifts (X/Y):                            0/-2
    Detected subpixel shifts (X/Y):                           0.357885632465/0.433837319984
    Calculated total shifts in fft pixel units (X/Y):         0.357885632465/-1.56616268002
    Calculated total shifts in reference pixel units (X/Y):   0.357885632465/-1.56616268002
    Calculated total shifts in target pixel units (X/Y):      0.357885632465/-1.56616268002
    Calculated map shifts (X,Y):				  3.578856324660592/15.661626799963415
    Calculated absolute shift vector length in map units:     16.065328089207995
    Calculated angle of shift vector in degrees from North:   192.8717191970359
    Original map info: ['UTM', 1, 1, 300000.0, 5900040.0, 10.0, 10.0, 33, 'North', 'WGS-84']
    Updated map info:  ['UTM', 1, 1, '300003.57885632466', '5900055.6616268', 10.0, 10.0, 33, 'North', 'WGS-84']

    'success'


correct shifts
~~~~~~~~~~~~~~

:meth:`CR.correct_shifts() <arosics.COREG.correct_shifts>` returns an
:class:`OrderedDict<collections.OrderedDict>` containing the co-registered
numpy array and its corresponding geoinformation.

.. code-block:: python

    >>> CR.correct_shifts()

    OrderedDict([('band', None),
                 ('is shifted', True),
                 ('is resampled', False),
                 ('updated map info',
                  ['UTM',
                   1,
                   1,
                   300003.57885632466,
                   5900025.6616268,
                   10.0,
                   10.0,
                   33,
                   'North',
                   'WGS-84']),
                 ('updated geotransform',
                  [300000.0, 10.0, 0.0, 5900040.0, 0.0, -10.0]),
                 ('updated projection',
                  'PROJCS["WGS 84 / UTM zone 33N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32633"]]'),
                 ('arr_shifted', array([[   0,    0,    0, ...,  953,  972, 1044],
                         [   0,    0,    0, ..., 1001,  973, 1019],
                         [   0,    0,    0, ...,  953,  985, 1020],
                         ...,
                         [   0,    0,    0, ...,  755,  763,  773],
                         [   0,    0,    0, ...,  760,  763,  749],
                         [9999, 9999, 9999, ..., 9999, 9999, 9999]], dtype=uint16)),
                 ('GeoArray_shifted',
                  <geoarray.GeoArray at 0x7f6c5a1cabe0>)])


To write the coregistered image to disk, the :class:`arosics.COREG` class needs to be instanced with a filepath given to
keyword 'path_out'. The output raster format can be any format supported by GDAL.
Find a list of supported formats here: https://gdal.org/drivers/raster/index.html


apply detected shifts to multiple images
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes it can be useful to apply the same shifts to multiple images, e.g., to different mask images derived from
the same satellite dataset. For this purpose you can calculate spatial shifts using the :class:`arosics.COREG` class
(see above) and then apply the calculated shifts to mulitple images using the :class:`arosics.DESHIFTER` class.
Take a look at the keyword arguments of the :class:`arosics.DESHIFTER` class when you need further adjustments
(e.g. output paths for the corrected images; aligned output grid, ...).

.. code-block:: python

    >>> from arosics import DESHIFTER

    >>> DESHIFTER(im_target1, CR.coreg_info).correct_shifts()
    >>> DESHIFTER(im_target2, CR.coreg_info).correct_shifts()

    OrderedDict([('band', None),
                 ('is shifted', True),
                 ('is resampled', False),
                 ('updated map info',
                  ['UTM',
                   1,
                   1,
                   300003.57885632466,
                   5900025.6616268,
                   10.0,
                   10.0,
                   33,
                   'North',
                   'WGS-84']),
                 ('updated geotransform',
                  [300000.0, 10.0, 0.0, 5900040.0, 0.0, -10.0]),
                 ('updated projection',
                  'PROJCS["WGS 84 / UTM zone 33N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32633"]]'),
                 ('arr_shifted', array([[   0,    0,    0, ...,  953,  972, 1044],
                         [   0,    0,    0, ..., 1001,  973, 1019],
                         [   0,    0,    0, ...,  953,  985, 1020],
                         ...,
                         [   0,    0,    0, ...,  755,  763,  773],
                         [   0,    0,    0, ...,  760,  763,  749],
                         [9999, 9999, 9999, ..., 9999, 9999, 9999]], dtype=uint16)),
                 ('GeoArray_shifted',
                  <geoarray.GeoArray at 0x7f6c5a1caa58>)])


----


Using the Shell console
-----------------------

The help instructions of the console interface can be accessed like this:

.. code-block:: bash

    $ arosics -h

Follow these instructions to run AROSICS from a shell console. For example, the most simple call for a global
co-registration would look like this:

.. code-block:: bash

    $ arosics global /path/to/your/ref_image.bsq /path/to/your/tgt_image.bsq
