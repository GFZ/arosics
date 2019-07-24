Local image co-registration
***************************

This local co-registration module of AROSICS has been designed to detect and correct geometric shifts present locally
in your input image. The class :class:`~arosics.COREG_LOCAL` calculates a grid of spatial shifts with points spread
over the whole overlap area of the input images. Based on this grid a correction of local shifts can be performed.


Using the Python API
--------------------

detect and correct local shifts - with input data on disk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from arosics import COREG_LOCAL

    >>> im_reference = '/path/to/your/ref_image.bsq'
    >>> im_target    = '/path/to/your/tgt_image.bsq'
    >>> kwargs = {
    >>>     'grid_res'     : 200,
    >>>     'window_size'  : (64,64),
    >>>     'path_out'     : 'auto',
    >>>     'projectDir'   : 'my_project',
    >>>     'q'            : False,
    >>> }

    >>> CRL = COREG_LOCAL(im_reference,im_target,**kwargs)
    >>> CRL.correct_shifts()

    Calculating actual data corner coordinates for reference image...
    Corner coordinates of reference image:
        [[319090.0, 5790510.0], [351800.0, 5899940.0], [409790.0, 5900040.0], [409790.0, 5790250.0], [319090.0, 5790250.0]]
    Calculating actual data corner coordinates for image to be shifted...
    Corner coordinates of image to be shifted:
        [[319460.0, 5790510.0], [352270.0, 5900040.0], [409790.0, 5900040.0], [409790.0, 5790250.0], [319460.0, 5790250.0]]
    Matching window position (X,Y): 372220.10753674706/5841066.947109019
    Calculating tie point grid (1977 points) in mode 'multiprocessing'...
        progress: |==================================================| 100.0% [1977/1977] Complete 9.75 sek
    Found 1144 valid GCPs.
    Correcting geometric shifts...
    Translating progress |==================================================| 100.0% Complete
    Warping progress     |==================================================| 100.0% Complete
    Writing GeoArray of size (10979, 10979) to /home/gfz-fe/scheffler/jupyter/arosics_jupyter/my_project/S2A_OPER_MSI_L1C_TL_SGS__20160608T153121_A005024_T33UUU_B03__shifted_to__S2A_OPER_MSI_L1C_TL_SGS__20160529T153631_A004881_T33UUU_B03.bsq.


    OrderedDict([('band', None),
                 ('is shifted', True),
                 ('is resampled', True),
                 ('updated map info',
                  ['UTM',
                   1,
                   1,
                   300000.0,
                   5900030.0,
                   10.0,
                   10.0,
                   33,
                   'North',
                   'WGS-84']),
                 ('updated geotransform',
                  [300000.0, 10.0, 0.0, 5900030.0, 0.0, -10.0]),
                 ('updated projection',
                  'PROJCS["WGS 84 / UTM zone 33N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32633"]]'),
                 ('arr_shifted', array([[   0,    0,    0, ..., 1034,  996, 1001],
                         [   0,    0,    0, ..., 1046, 1114, 1124],
                         [   0,    0,    0, ..., 1021, 1126, 1148],
                         ...,
                         [   0,    0,    0, ...,  760,  769,  805],
                         [   0,    0,    0, ...,  762,  755,  765],
                         [   0,    0,    0, ...,    0,    0,    0]], dtype=uint16)),
                 ('GeoArray_shifted',
                  <geoarray.GeoArray at 0x7f451ac14a90>)])


detect and correct local shifts - without any disk access
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All you have to do is to instanciate :class:`arosics.COREG_LOCAL` with two instances of the :class:`geoarray.GeoArray`
class as described above.


.. code-block:: python

    >>> from geoarray import GeoArray

    >>> CRL = COREG_LOCAL(GeoArray(ref_ndarray, ref_gt, ref_prj),
    >>>                   GeoArray(tgt_ndarray, tgt_gt, tgt_prj),
    >>>                   **kwargs)
    >>> CRL.correct_shifts()


visualize tie point grid with INITIAL shifts present in your input target image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the method :meth:`CRL.view_CoRegPoints()<arosics.COREG_LOCAL.view_CoRegPoints>` to visualize the tie point grid with
the calculated absolute lenghts of the shift vectors (the unit corresponds to the input projection - UTM in the shown
example, thus the unit is 'meters'.).

.. note::

    A calculation of reliable shifts above cloud covered areas is not possible.
    In the current version of AROSICS these areas are not masked. A proper masking is planned.


.. code-block:: python

    >>> CRL.view_CoRegPoints(figsize=(15,15), backgroundIm='ref')

    Note: array has been downsampled to 1000 x 1000 for faster visualization.

.. image:: ../images/output_40_1.png


The output figure shows the calculated absolute lenghts of the shift vectors - in this case with shifts up to ~25 meters.


visualize tie point grid with shifts present AFTER shift correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The remaining shifts after local correction can be calculated and visualized by instanciating the
:class:`arosics.COREG_LOCAL` with the output path of the above instance of :class:`COREG_LOCAL<arosics.COREG_LOCAL>`.

.. code-block:: python

    >>> CRL_after_corr = COREG_LOCAL(im_reference, CRL.path_out, **kwargs)
    >>> CRL_after_corr.view_CoRegPoints(figsize=(15,15),backgroundIm='ref')

    Calculating actual data corner coordinates for reference image...
    Corner coordinates of reference image:
        [[319090.0, 5790510.0], [351800.0, 5899940.0], [409790.0, 5900040.0], [409790.0, 5790250.0], [319090.0, 5790250.0]]
    Calculating actual data corner coordinates for image to be shifted...
    Corner coordinates of image to be shifted:
        [[319460.0, 5790540.0], [352270.0, 5900030.0], [409780.0, 5900030.0], [409780.0, 5790260.0], [322970.0, 5790250.0], [319460.0, 5790280.0]]
    Matching window position (X,Y): 372216.38593955856/5841068.390957352
    Note: array has been downsampled to 1000 x 1000 for faster visualization.
    Calculating tie point grid (1977 points) in mode 'multiprocessing'...
        progress: |==================================================| 100.0% [1977/1977] Complete 10.78 sek

.. image:: ../images/output_44_1.png


The output figure shows a significant reduction of geometric shifts.


show the points table of the calculated tie point grid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    Point records where no valid match has been found are filled with -9999.

.. code-block:: python

    >>> CRL.CoRegPoints_table


.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>geometry</th>
          <th>POINT_ID</th>
          <th>X_IM</th>
          <th>Y_IM</th>
          <th>X_UTM</th>
          <th>Y_UTM</th>
          <th>X_WIN_SIZE</th>
          <th>Y_WIN_SIZE</th>
          <th>X_SHIFT_PX</th>
          <th>Y_SHIFT_PX</th>
          <th>X_SHIFT_M</th>
          <th>Y_SHIFT_M</th>
          <th>ABS_SHIFT</th>
          <th>ANGLE</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>POINT (352000 5898040)</td>
          <td>81</td>
          <td>5200</td>
          <td>200</td>
          <td>352000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>POINT (354000 5898040)</td>
          <td>82</td>
          <td>5400</td>
          <td>200</td>
          <td>354000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.372470</td>
          <td>-0.285500</td>
          <td>3.724704</td>
          <td>2.855005</td>
          <td>4.693024</td>
          <td>232.529646</td>
        </tr>
        <tr>
          <th>2</th>
          <td>POINT (356000 5898040)</td>
          <td>83</td>
          <td>5600</td>
          <td>200</td>
          <td>356000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.260948</td>
          <td>-0.293539</td>
          <td>2.609479</td>
          <td>2.935389</td>
          <td>3.927580</td>
          <td>221.636201</td>
        </tr>
        <tr>
          <th>3</th>
          <td>POINT (358000 5898040)</td>
          <td>84</td>
          <td>5800</td>
          <td>200</td>
          <td>358000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>4</th>
          <td>POINT (360000 5898040)</td>
          <td>85</td>
          <td>6000</td>
          <td>200</td>
          <td>360000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>5</th>
          <td>POINT (362000 5898040)</td>
          <td>86</td>
          <td>6200</td>
          <td>200</td>
          <td>362000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>6</th>
          <td>POINT (364000 5898040)</td>
          <td>87</td>
          <td>6400</td>
          <td>200</td>
          <td>364000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.141693</td>
          <td>0.187036</td>
          <td>1.416935</td>
          <td>-1.870360</td>
          <td>2.346476</td>
          <td>322.853405</td>
        </tr>
        <tr>
          <th>7</th>
          <td>POINT (366000 5898040)</td>
          <td>88</td>
          <td>6600</td>
          <td>200</td>
          <td>366000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-0.230941</td>
          <td>0.121139</td>
          <td>-2.309409</td>
          <td>-1.211389</td>
          <td>2.607841</td>
          <td>62.320969</td>
        </tr>
        <tr>
          <th>8</th>
          <td>POINT (368000 5898040)</td>
          <td>89</td>
          <td>6800</td>
          <td>200</td>
          <td>368000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>9</th>
          <td>POINT (370000 5898040)</td>
          <td>90</td>
          <td>7000</td>
          <td>200</td>
          <td>370000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-0.035693</td>
          <td>0.084596</td>
          <td>-0.356928</td>
          <td>-0.845957</td>
          <td>0.918172</td>
          <td>22.875994</td>
        </tr>
        <tr>
          <th>10</th>
          <td>POINT (372000 5898040)</td>
          <td>91</td>
          <td>7200</td>
          <td>200</td>
          <td>372000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>11</th>
          <td>POINT (374000 5898040)</td>
          <td>92</td>
          <td>7400</td>
          <td>200</td>
          <td>374000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>12</th>
          <td>POINT (376000 5898040)</td>
          <td>93</td>
          <td>7600</td>
          <td>200</td>
          <td>376000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>13</th>
          <td>POINT (378000 5898040)</td>
          <td>94</td>
          <td>7800</td>
          <td>200</td>
          <td>378000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>14</th>
          <td>POINT (380000 5898040)</td>
          <td>95</td>
          <td>8000</td>
          <td>200</td>
          <td>380000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>15</th>
          <td>POINT (382000 5898040)</td>
          <td>96</td>
          <td>8200</td>
          <td>200</td>
          <td>382000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>16</th>
          <td>POINT (384000 5898040)</td>
          <td>97</td>
          <td>8400</td>
          <td>200</td>
          <td>384000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>17</th>
          <td>POINT (386000 5898040)</td>
          <td>98</td>
          <td>8600</td>
          <td>200</td>
          <td>386000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>18</th>
          <td>POINT (388000 5898040)</td>
          <td>99</td>
          <td>8800</td>
          <td>200</td>
          <td>388000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.656098</td>
          <td>2.533985</td>
          <td>6.560977</td>
          <td>-25.339852</td>
          <td>26.175457</td>
          <td>345.483797</td>
        </tr>
        <tr>
          <th>19</th>
          <td>POINT (390000 5898040)</td>
          <td>100</td>
          <td>9000</td>
          <td>200</td>
          <td>390000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>20</th>
          <td>POINT (392000 5898040)</td>
          <td>101</td>
          <td>9200</td>
          <td>200</td>
          <td>392000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>21</th>
          <td>POINT (394000 5898040)</td>
          <td>102</td>
          <td>9400</td>
          <td>200</td>
          <td>394000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>22</th>
          <td>POINT (396000 5898040)</td>
          <td>103</td>
          <td>9600</td>
          <td>200</td>
          <td>396000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>23</th>
          <td>POINT (398000 5898040)</td>
          <td>104</td>
          <td>9800</td>
          <td>200</td>
          <td>398000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>24</th>
          <td>POINT (400000 5898040)</td>
          <td>105</td>
          <td>10000</td>
          <td>200</td>
          <td>400000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-0.147210</td>
          <td>-0.223871</td>
          <td>-1.472098</td>
          <td>2.238708</td>
          <td>2.679344</td>
          <td>146.672433</td>
        </tr>
        <tr>
          <th>25</th>
          <td>POINT (402000 5898040)</td>
          <td>106</td>
          <td>10200</td>
          <td>200</td>
          <td>402000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>26</th>
          <td>POINT (404000 5898040)</td>
          <td>107</td>
          <td>10400</td>
          <td>200</td>
          <td>404000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>27</th>
          <td>POINT (406000 5898040)</td>
          <td>108</td>
          <td>10600</td>
          <td>200</td>
          <td>406000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.249318</td>
          <td>0.214416</td>
          <td>2.493182</td>
          <td>-2.144158</td>
          <td>3.288369</td>
          <td>310.695805</td>
        </tr>
        <tr>
          <th>28</th>
          <td>POINT (408000 5898040)</td>
          <td>109</td>
          <td>10800</td>
          <td>200</td>
          <td>408000.0</td>
          <td>5898040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.372511</td>
          <td>-1.410450</td>
          <td>3.725107</td>
          <td>14.104504</td>
          <td>14.588127</td>
          <td>194.794441</td>
        </tr>
        <tr>
          <th>29</th>
          <td>POINT (352000 5896040)</td>
          <td>136</td>
          <td>5200</td>
          <td>400</td>
          <td>352000.0</td>
          <td>5896040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>1947</th>
          <td>POINT (350000 5792040)</td>
          <td>2995</td>
          <td>5000</td>
          <td>10800</td>
          <td>350000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.209144</td>
          <td>-1.750348</td>
          <td>2.091443</td>
          <td>17.503485</td>
          <td>17.627992</td>
          <td>186.813809</td>
        </tr>
        <tr>
          <th>1948</th>
          <td>POINT (352000 5792040)</td>
          <td>2996</td>
          <td>5200</td>
          <td>10800</td>
          <td>352000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.367216</td>
          <td>-1.643834</td>
          <td>3.672159</td>
          <td>16.438337</td>
          <td>16.843505</td>
          <td>192.592548</td>
        </tr>
        <tr>
          <th>1949</th>
          <td>POINT (354000 5792040)</td>
          <td>2997</td>
          <td>5400</td>
          <td>10800</td>
          <td>354000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.288332</td>
          <td>-1.711756</td>
          <td>2.883320</td>
          <td>17.117562</td>
          <td>17.358700</td>
          <td>189.561275</td>
        </tr>
        <tr>
          <th>1950</th>
          <td>POINT (356000 5792040)</td>
          <td>2998</td>
          <td>5600</td>
          <td>10800</td>
          <td>356000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.349523</td>
          <td>-1.629551</td>
          <td>3.495229</td>
          <td>16.295510</td>
          <td>16.666142</td>
          <td>192.105965</td>
        </tr>
        <tr>
          <th>1951</th>
          <td>POINT (358000 5792040)</td>
          <td>2999</td>
          <td>5800</td>
          <td>10800</td>
          <td>358000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>1952</th>
          <td>POINT (360000 5792040)</td>
          <td>3000</td>
          <td>6000</td>
          <td>10800</td>
          <td>360000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.356829</td>
          <td>-1.353932</td>
          <td>3.568290</td>
          <td>13.539322</td>
          <td>14.001641</td>
          <td>194.764576</td>
        </tr>
        <tr>
          <th>1953</th>
          <td>POINT (362000 5792040)</td>
          <td>3001</td>
          <td>6200</td>
          <td>10800</td>
          <td>362000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.332107</td>
          <td>-1.475567</td>
          <td>3.321073</td>
          <td>14.755674</td>
          <td>15.124795</td>
          <td>192.684252</td>
        </tr>
        <tr>
          <th>1954</th>
          <td>POINT (364000 5792040)</td>
          <td>3002</td>
          <td>6400</td>
          <td>10800</td>
          <td>364000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.260931</td>
          <td>-1.235276</td>
          <td>2.609308</td>
          <td>12.352761</td>
          <td>12.625339</td>
          <td>191.927413</td>
        </tr>
        <tr>
          <th>1955</th>
          <td>POINT (366000 5792040)</td>
          <td>3003</td>
          <td>6600</td>
          <td>10800</td>
          <td>366000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>1956</th>
          <td>POINT (368000 5792040)</td>
          <td>3004</td>
          <td>6800</td>
          <td>10800</td>
          <td>368000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.230095</td>
          <td>-1.258021</td>
          <td>2.300948</td>
          <td>12.580208</td>
          <td>12.788901</td>
          <td>190.364959</td>
        </tr>
        <tr>
          <th>1957</th>
          <td>POINT (370000 5792040)</td>
          <td>3005</td>
          <td>7000</td>
          <td>10800</td>
          <td>370000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-0.096170</td>
          <td>-0.463691</td>
          <td>-0.961701</td>
          <td>4.636910</td>
          <td>4.735589</td>
          <td>168.282899</td>
        </tr>
        <tr>
          <th>1958</th>
          <td>POINT (372000 5792040)</td>
          <td>3006</td>
          <td>7200</td>
          <td>10800</td>
          <td>372000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.194545</td>
          <td>0.126613</td>
          <td>1.945447</td>
          <td>-1.266134</td>
          <td>2.321176</td>
          <td>303.056848</td>
        </tr>
        <tr>
          <th>1959</th>
          <td>POINT (374000 5792040)</td>
          <td>3007</td>
          <td>7400</td>
          <td>10800</td>
          <td>374000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
          <td>-9999.000000</td>
        </tr>
        <tr>
          <th>1960</th>
          <td>POINT (376000 5792040)</td>
          <td>3008</td>
          <td>7600</td>
          <td>10800</td>
          <td>376000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-0.192273</td>
          <td>-0.410461</td>
          <td>-1.922730</td>
          <td>4.104609</td>
          <td>4.532627</td>
          <td>154.900105</td>
        </tr>
        <tr>
          <th>1961</th>
          <td>POINT (378000 5792040)</td>
          <td>3009</td>
          <td>7800</td>
          <td>10800</td>
          <td>378000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.411476</td>
          <td>-1.231980</td>
          <td>4.114758</td>
          <td>12.319801</td>
          <td>12.988792</td>
          <td>198.469086</td>
        </tr>
        <tr>
          <th>1962</th>
          <td>POINT (380000 5792040)</td>
          <td>3010</td>
          <td>8000</td>
          <td>10800</td>
          <td>380000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.262658</td>
          <td>-0.490337</td>
          <td>2.626580</td>
          <td>4.903369</td>
          <td>5.562549</td>
          <td>208.176553</td>
        </tr>
        <tr>
          <th>1963</th>
          <td>POINT (382000 5792040)</td>
          <td>3011</td>
          <td>8200</td>
          <td>10800</td>
          <td>382000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.186922</td>
          <td>-1.105403</td>
          <td>1.869221</td>
          <td>11.054032</td>
          <td>11.210959</td>
          <td>189.597841</td>
        </tr>
        <tr>
          <th>1964</th>
          <td>POINT (384000 5792040)</td>
          <td>3012</td>
          <td>8400</td>
          <td>10800</td>
          <td>384000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-0.267606</td>
          <td>0.342886</td>
          <td>-2.676062</td>
          <td>-3.428858</td>
          <td>4.349526</td>
          <td>37.970358</td>
        </tr>
        <tr>
          <th>1965</th>
          <td>POINT (386000 5792040)</td>
          <td>3013</td>
          <td>8600</td>
          <td>10800</td>
          <td>386000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.368027</td>
          <td>-1.232417</td>
          <td>3.680269</td>
          <td>12.324169</td>
          <td>12.861941</td>
          <td>196.626786</td>
        </tr>
        <tr>
          <th>1966</th>
          <td>POINT (388000 5792040)</td>
          <td>3014</td>
          <td>8800</td>
          <td>10800</td>
          <td>388000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.405260</td>
          <td>-0.790863</td>
          <td>4.052597</td>
          <td>7.908634</td>
          <td>8.886509</td>
          <td>207.131823</td>
        </tr>
        <tr>
          <th>1967</th>
          <td>POINT (390000 5792040)</td>
          <td>3015</td>
          <td>9000</td>
          <td>10800</td>
          <td>390000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.372675</td>
          <td>-1.224506</td>
          <td>3.726746</td>
          <td>12.245065</td>
          <td>12.799619</td>
          <td>196.927457</td>
        </tr>
        <tr>
          <th>1968</th>
          <td>POINT (392000 5792040)</td>
          <td>3016</td>
          <td>9200</td>
          <td>10800</td>
          <td>392000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.386730</td>
          <td>-1.438051</td>
          <td>3.867297</td>
          <td>14.380515</td>
          <td>14.891447</td>
          <td>195.052215</td>
        </tr>
        <tr>
          <th>1969</th>
          <td>POINT (394000 5792040)</td>
          <td>3017</td>
          <td>9400</td>
          <td>10800</td>
          <td>394000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.433132</td>
          <td>-1.209992</td>
          <td>4.331321</td>
          <td>12.099919</td>
          <td>12.851785</td>
          <td>199.695480</td>
        </tr>
        <tr>
          <th>1970</th>
          <td>POINT (396000 5792040)</td>
          <td>3018</td>
          <td>9600</td>
          <td>10800</td>
          <td>396000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.410025</td>
          <td>-0.784237</td>
          <td>4.100254</td>
          <td>7.842365</td>
          <td>8.849563</td>
          <td>207.602090</td>
        </tr>
        <tr>
          <th>1971</th>
          <td>POINT (398000 5792040)</td>
          <td>3019</td>
          <td>9800</td>
          <td>10800</td>
          <td>398000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.376237</td>
          <td>-1.138838</td>
          <td>3.762373</td>
          <td>11.388382</td>
          <td>11.993777</td>
          <td>198.281973</td>
        </tr>
        <tr>
          <th>1972</th>
          <td>POINT (400000 5792040)</td>
          <td>3020</td>
          <td>10000</td>
          <td>10800</td>
          <td>400000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.071339</td>
          <td>-0.964923</td>
          <td>0.713385</td>
          <td>9.649233</td>
          <td>9.675568</td>
          <td>184.228288</td>
        </tr>
        <tr>
          <th>1973</th>
          <td>POINT (402000 5792040)</td>
          <td>3021</td>
          <td>10200</td>
          <td>10800</td>
          <td>402000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.246210</td>
          <td>-1.129963</td>
          <td>2.462099</td>
          <td>11.299628</td>
          <td>11.564754</td>
          <td>192.292166</td>
        </tr>
        <tr>
          <th>1974</th>
          <td>POINT (404000 5792040)</td>
          <td>3022</td>
          <td>10400</td>
          <td>10800</td>
          <td>404000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>-0.263890</td>
          <td>-0.903314</td>
          <td>-2.638901</td>
          <td>9.033142</td>
          <td>9.410709</td>
          <td>163.715048</td>
        </tr>
        <tr>
          <th>1975</th>
          <td>POINT (406000 5792040)</td>
          <td>3023</td>
          <td>10600</td>
          <td>10800</td>
          <td>406000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.239090</td>
          <td>-1.235482</td>
          <td>2.390904</td>
          <td>12.354817</td>
          <td>12.584034</td>
          <td>190.952493</td>
        </tr>
        <tr>
          <th>1976</th>
          <td>POINT (408000 5792040)</td>
          <td>3024</td>
          <td>10800</td>
          <td>10800</td>
          <td>408000.0</td>
          <td>5792040.0</td>
          <td>64</td>
          <td>64</td>
          <td>0.272772</td>
          <td>-0.964375</td>
          <td>2.727717</td>
          <td>9.643754</td>
          <td>10.022098</td>
          <td>195.793451</td>
        </tr>
      </tbody>
    </table>
    <p>1977 rows × 14 columns</p>
    </div>


export tie point grid to an ESRI point shapefile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> CRL.tiepoint_grid.to_PointShapefile(path_out='/path/to/your/output_shapefile.shp')


----


Using the Shell console
-----------------------

Follow these instructions to run AROSICS from a shell console. For example, the most simple call for a local
co-registration would look like this:

.. code-block:: bash

    $ python arosics_cli.py local /path/to/your/ref_image.bsq /path/to/your/tgt_image.bsq 50
