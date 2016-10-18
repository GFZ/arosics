
# Description

Perform subpixel coregistration of two satellite image datasets using Fourier Shift Theorem proposed by Foroosh et al. 2002: Foroosh, H., Zerubia, J. B., & Berthod, M. (2002). Extension of phase correlation to subpixel registration. IEEE Transactions on Image Processing, 11(3), 188-199. doi:10.1109/83.988953); Python implementation by Daniel Scheffler (daniel.scheffler [at] gfz-potsdam [dot] de).

The program detects and corrects global X/Y-shifts between two input images in the subpixel scale, that are often present in satellite imagery. It does not correct scaling or rotation issues and will not apply any higher grade transformation. Therefore it will also not correct for shifts that are locally present in the input images.

Prerequisites and hints:
The input images can have any GDAL compatible image format (http://www.gdal.org/formats_list.html). Both of them must be georeferenced. In case of ENVI files, this means they must have a 'map info' and a 'coordinate system string' as attributes of their header file. Different projection systems are currently not supported. The input images must have a geographic overlap but clipping them to same geographical extent is NOT neccessary. Please do not perform any spatial resampling of the input images before applying this algorithm. Any needed resampling of the data is done automatically. Thus the input images can have different spatial resolutions. The current algorithm will not perform any ortho-rectification. So please use ortho-rectified input data in order to prevent local shifts in the output image. By default the calculated subpixel-shifts are applied to the header file of the output image. No spatial resampling is done automatically as long as the both input images have the same projection. If you need the output image to be aligned to the reference image coordinate grid (by using an appropriate resampling algorithm), use the '-align_grids' option. The image overlap area is automatically calculated. Thereby no-data regions within the images are standardly respected. Providing the map coordinates of the actual data corners lets you save some calculation time, because in this case the automatic algorithm can be skipped. The no-data value of each image is automatically derived from the image corners. The verbose program mode gives some more output about the interim results, shows some figures and writes the used footprint and overlap polygons to disk. The figures must be manually closed in in order to continue the processing.


# Install

For now, there is no automatic install script. Just clone the repository, install the dependencies and add the root directory of CoReg_Sat to your PATH environment variable.

CoReg_Sat has been tested with Python 3.5. It is not completely compatible with Python 2.7 (at least at the moment).

The following non-standard Python libraries are required:
    - gdal, osr, ogr
    - geopandas
    - pykrige
    - argparse
    - shapely
    - pyfftw is optional but will speed up calculation
    
In addition clone the repository "py_tools_ds" and add its root directory to your PATH environment variable:

    https://gitext.gfz-potsdam.de/danschef/py_tools_ds
Since py_tools_ds is not a public repository right now, contact Daniel Scheffler if you can not access it.

# Modules

## CoReg

This module calculates spatial shifts and performs a global correction (based on a single matching window).

### Python Interface

#### calculate spatial shifts - with input data on disk


```python
from CoReg_Sat import COREG

im_reference = '/path/to/your/ref_image.bsq'
im_target    = '/path/to/your/tgt_image.bsq'

CR = COREG(im_reference, im_target, wp=(354223, 5805559), ws=(256,256))
CR.calculate_spatial_shifts()
```

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


#### calculate spatial shifts - without any disk access


```python
from py_tools_ds.ptds import GeoArray
from CoReg_Sat import COREG

im_reference = '/path/to/your/ref_image.bsq'
im_target    = '/path/to/your/tgt_image.bsq'

# get a sample numpy array with corresponding geoinformation as reference image
geoArr  = GeoArray(im_reference)

ref_ndarray = geoArr[:]            # numpy.ndarray with shape (10980, 10980)
ref_gt      = geoArr.geotransform  # GDAL geotransform: (300000.0, 10.0, 0.0, 5900040.0, 0.0, -10.0)
ref_prj     = geoArr.projection    # projection as WKT string ('PROJCS["WGS 84 / UTM zone 33N....')

# get a sample numpy array with corresponding geoinformation as target image
geoArr  = GeoArray(im_target)

tgt_ndarray = geoArr[:]            # numpy.ndarray with shape (10980, 10980)
tgt_gt      = geoArr.geotransform  # GDAL geotransform: (300000.0, 10.0, 0.0, 5900040.0, 0.0, -10.0)
tgt_prj     = geoArr.projection    # projection as WKT string ('PROJCS["WGS 84 / UTM zone 33N....')

# pass an instance of GeoArray to COREG and calculate spatial shifts
geoArr_reference = GeoArray(ref_ndarray, ref_gt, ref_prj)
geoArr_target    = GeoArray(tgt_ndarray, tgt_gt, tgt_prj)

CR = COREG(geoArr_reference, geoArr_target, wp=(354223, 5805559), ws=(256,256))
CR.calculate_spatial_shifts()
```

    Calculating actual data corner coordinates for reference image...
    Corner coordinates of reference image:
    	[[300000.0, 5848140.0], [409790.0, 5848140.0], [409790.0, 5790250.0], [300000.0, 5790250.0]]
    Calculating actual data corner coordinates for image to be shifted...
    Corner coordinates of image to be shifted:
    	[[300000.0, 5847770.0], [409790.0, 5847770.0], [409790.0, 5790250.0], [300000.0, 5790250.0]]
    Matching window position (X,Y): 354223/5805559
    Detected integer shifts (X/Y):       0/-2
    Detected subpixel shifts (X/Y):      0.357885632465/0.433837319984
    Calculated total shifts in fft pixel units (X/Y):         0.357885632465/-1.56616268002
    Calculated total shifts in reference pixel units (X/Y):   0.357885632465/-1.56616268002
    Calculated total shifts in target pixel units (X/Y):      0.357885632465/-1.56616268002
    Calculated map shifts (X,Y):				  3.578856324660592 15.661626799963415
    Original map info: ['UTM', 1, 1, 300000.0, 5900040.0, 10.0, 10.0, 33, 'North', 'WGS-84']
    Updated map info:  ['UTM', 1, 1, '300003.57885632466', '5900055.6616268', 10.0, 10.0, 33, 'North', 'WGS-84']


#### correct shifts

CR.correct_shifts() returns an an OrderedDict containing the coregistered numpy array and its corresponding geoinformation.


```python
CR.correct_shifts()
```

    Writing GeoArray of size (10978, 10978) to /home/gfz-fe/scheffler/temp.





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
                 ('updated projection',
                  'PROJCS["WGS 84 / UTM zone 33N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32633"]]'),
                 ('arr_shifted', array([[   0,    0,    0, ...,  953,  972, 1044],
                         [   0,    0,    0, ..., 1001,  973, 1019],
                         [   0,    0,    0, ...,  953,  985, 1020],
                         ..., 
                         [   0,    0,    0, ...,  755,  763,  773],
                         [   0,    0,    0, ...,  760,  763,  749],
                         [9999, 9999, 9999, ..., 9999, 9999, 9999]], dtype=uint16))])



To write the coregistered image to disk, the COREG class needs to be instanced with a filepath given to keyword 'path_out'.

### Shell console interface

The help instructions of the console interface can be accessed like this:


```python
cd /path/to/CoReg_Sat/
python ./dsc__CoReg_Sat_FourierShiftTheorem.py -h
```

Follow these instructions to run CoReg_Sat from a shell console. For example, the most simple call would be like this:


```python
python ./dsc__CoReg_Sat_FourierShiftTheorem.py /path/to/your/ref_image.bsq /path/to/your/tgt_image.bsq
```

## Geometric quality grid

This module calculates a grid of spatial shifts with points spread over the whole overlap area of the input images. At the moment, a shift correction using all of these points is planned but not yet implemented.

### Python interface

#### calculate geometric quality grid - with input data on disk


```python
from CoReg_Sat import Geom_Quality_Grid

im_reference = '/path/to/your/ref_image.bsq'
im_target    = '/path/to/your/tgt_image.bsq'
kwargs = {
    'grid_res'     : 100,
    'window_size'  : (256,256),
    'calc_corners' : True,
    'dir_out'      : '/path/to/your/output/',
    'projectName'  : 'my_project',
    'q'            : True,
}

GQG = Geom_Quality_Grid(im_reference,im_target,**kwargs)
GQG.get_quality_grid()
```




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
      <th>WIN_SIZE</th>
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
      <th>162</th>
      <td>POINT (352000 5899040)</td>
      <td>162</td>
      <td>5200</td>
      <td>100</td>
      <td>352000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>163</th>
      <td>POINT (353000 5899040)</td>
      <td>163</td>
      <td>5300</td>
      <td>100</td>
      <td>353000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>164</th>
      <td>POINT (354000 5899040)</td>
      <td>164</td>
      <td>5400</td>
      <td>100</td>
      <td>354000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>165</th>
      <td>POINT (355000 5899040)</td>
      <td>165</td>
      <td>5500</td>
      <td>100</td>
      <td>355000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-0.247889</td>
      <td>-0.251688</td>
      <td>-2.478889</td>
      <td>-2.516880</td>
      <td>3.532644</td>
      <td>135.435699</td>
    </tr>
    <tr>
      <th>166</th>
      <td>POINT (356000 5899040)</td>
      <td>166</td>
      <td>5600</td>
      <td>100</td>
      <td>356000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>167</th>
      <td>POINT (357000 5899040)</td>
      <td>167</td>
      <td>5700</td>
      <td>100</td>
      <td>357000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>168</th>
      <td>POINT (358000 5899040)</td>
      <td>168</td>
      <td>5800</td>
      <td>100</td>
      <td>358000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>169</th>
      <td>POINT (359000 5899040)</td>
      <td>169</td>
      <td>5900</td>
      <td>100</td>
      <td>359000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>170</th>
      <td>POINT (360000 5899040)</td>
      <td>170</td>
      <td>6000</td>
      <td>100</td>
      <td>360000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>171</th>
      <td>POINT (361000 5899040)</td>
      <td>171</td>
      <td>6100</td>
      <td>100</td>
      <td>361000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>172</th>
      <td>POINT (362000 5899040)</td>
      <td>172</td>
      <td>6200</td>
      <td>100</td>
      <td>362000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>173</th>
      <td>POINT (363000 5899040)</td>
      <td>173</td>
      <td>6300</td>
      <td>100</td>
      <td>363000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>0.185042</td>
      <td>-26.788517</td>
      <td>1.850425</td>
      <td>-267.885169</td>
      <td>267.891560</td>
      <td>180.395766</td>
    </tr>
    <tr>
      <th>174</th>
      <td>POINT (364000 5899040)</td>
      <td>174</td>
      <td>6400</td>
      <td>100</td>
      <td>364000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>175</th>
      <td>POINT (365000 5899040)</td>
      <td>175</td>
      <td>6500</td>
      <td>100</td>
      <td>365000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>176</th>
      <td>POINT (366000 5899040)</td>
      <td>176</td>
      <td>6600</td>
      <td>100</td>
      <td>366000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>177</th>
      <td>POINT (367000 5899040)</td>
      <td>177</td>
      <td>6700</td>
      <td>100</td>
      <td>367000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-0.255553</td>
      <td>0.188893</td>
      <td>-2.555527</td>
      <td>1.888925</td>
      <td>3.177854</td>
      <td>53.529930</td>
    </tr>
    <tr>
      <th>178</th>
      <td>POINT (368000 5899040)</td>
      <td>178</td>
      <td>6800</td>
      <td>100</td>
      <td>368000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>179</th>
      <td>POINT (369000 5899040)</td>
      <td>179</td>
      <td>6900</td>
      <td>100</td>
      <td>369000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>180</th>
      <td>POINT (370000 5899040)</td>
      <td>180</td>
      <td>7000</td>
      <td>100</td>
      <td>370000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>181</th>
      <td>POINT (371000 5899040)</td>
      <td>181</td>
      <td>7100</td>
      <td>100</td>
      <td>371000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>182</th>
      <td>POINT (372000 5899040)</td>
      <td>182</td>
      <td>7200</td>
      <td>100</td>
      <td>372000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>183</th>
      <td>POINT (373000 5899040)</td>
      <td>183</td>
      <td>7300</td>
      <td>100</td>
      <td>373000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>184</th>
      <td>POINT (374000 5899040)</td>
      <td>184</td>
      <td>7400</td>
      <td>100</td>
      <td>374000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>1.687416</td>
      <td>0.206791</td>
      <td>16.874158</td>
      <td>2.067906</td>
      <td>17.000396</td>
      <td>276.986685</td>
    </tr>
    <tr>
      <th>185</th>
      <td>POINT (375000 5899040)</td>
      <td>185</td>
      <td>7500</td>
      <td>100</td>
      <td>375000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>0.715882</td>
      <td>-1.562482</td>
      <td>7.158821</td>
      <td>-15.624824</td>
      <td>17.186734</td>
      <td>204.615818</td>
    </tr>
    <tr>
      <th>186</th>
      <td>POINT (376000 5899040)</td>
      <td>186</td>
      <td>7600</td>
      <td>100</td>
      <td>376000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>187</th>
      <td>POINT (377000 5899040)</td>
      <td>187</td>
      <td>7700</td>
      <td>100</td>
      <td>377000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>188</th>
      <td>POINT (378000 5899040)</td>
      <td>188</td>
      <td>7800</td>
      <td>100</td>
      <td>378000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>189</th>
      <td>POINT (379000 5899040)</td>
      <td>189</td>
      <td>7900</td>
      <td>100</td>
      <td>379000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>0.195506</td>
      <td>42.665374</td>
      <td>1.955060</td>
      <td>426.653741</td>
      <td>426.658220</td>
      <td>359.737455</td>
    </tr>
    <tr>
      <th>190</th>
      <td>POINT (380000 5899040)</td>
      <td>190</td>
      <td>8000</td>
      <td>100</td>
      <td>380000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
      <td>-9999.000000</td>
    </tr>
    <tr>
      <th>191</th>
      <td>POINT (381000 5899040)</td>
      <td>191</td>
      <td>8100</td>
      <td>100</td>
      <td>381000.0</td>
      <td>5899040.0</td>
      <td>228</td>
      <td>-0.210588</td>
      <td>0.042504</td>
      <td>-2.105879</td>
      <td>0.425039</td>
      <td>2.148345</td>
      <td>78.589040</td>
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
    </tr>
    <tr>
      <th>12070</th>
      <td>POINT (380000 5791040)</td>
      <td>12070</td>
      <td>8000</td>
      <td>10900</td>
      <td>380000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.423103</td>
      <td>-1.284985</td>
      <td>4.231027</td>
      <td>-12.849851</td>
      <td>13.528498</td>
      <td>198.224989</td>
    </tr>
    <tr>
      <th>12071</th>
      <td>POINT (381000 5791040)</td>
      <td>12071</td>
      <td>8100</td>
      <td>10900</td>
      <td>381000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.295017</td>
      <td>-1.258081</td>
      <td>2.950171</td>
      <td>-12.580810</td>
      <td>12.922086</td>
      <td>193.197273</td>
    </tr>
    <tr>
      <th>12072</th>
      <td>POINT (382000 5791040)</td>
      <td>12072</td>
      <td>8200</td>
      <td>10900</td>
      <td>382000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.351648</td>
      <td>-1.200180</td>
      <td>3.516481</td>
      <td>-12.001803</td>
      <td>12.506354</td>
      <td>196.330375</td>
    </tr>
    <tr>
      <th>12073</th>
      <td>POINT (383000 5791040)</td>
      <td>12073</td>
      <td>8300</td>
      <td>10900</td>
      <td>383000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.404020</td>
      <td>-1.261354</td>
      <td>4.040197</td>
      <td>-12.613543</td>
      <td>13.244798</td>
      <td>197.760587</td>
    </tr>
    <tr>
      <th>12074</th>
      <td>POINT (384000 5791040)</td>
      <td>12074</td>
      <td>8400</td>
      <td>10900</td>
      <td>384000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.440928</td>
      <td>-1.271966</td>
      <td>4.409277</td>
      <td>-12.719659</td>
      <td>13.462223</td>
      <td>199.118903</td>
    </tr>
    <tr>
      <th>12075</th>
      <td>POINT (385000 5791040)</td>
      <td>12075</td>
      <td>8500</td>
      <td>10900</td>
      <td>385000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.448053</td>
      <td>-1.264458</td>
      <td>4.480532</td>
      <td>-12.644577</td>
      <td>13.414936</td>
      <td>199.511483</td>
    </tr>
    <tr>
      <th>12076</th>
      <td>POINT (386000 5791040)</td>
      <td>12076</td>
      <td>8600</td>
      <td>10900</td>
      <td>386000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.368041</td>
      <td>-1.225450</td>
      <td>3.680410</td>
      <td>-12.254503</td>
      <td>12.795244</td>
      <td>196.716656</td>
    </tr>
    <tr>
      <th>12077</th>
      <td>POINT (387000 5791040)</td>
      <td>12077</td>
      <td>8700</td>
      <td>10900</td>
      <td>387000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.421414</td>
      <td>-1.248958</td>
      <td>4.214138</td>
      <td>-12.489578</td>
      <td>13.181370</td>
      <td>198.645029</td>
    </tr>
    <tr>
      <th>12078</th>
      <td>POINT (388000 5791040)</td>
      <td>12078</td>
      <td>8800</td>
      <td>10900</td>
      <td>388000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.390062</td>
      <td>-1.179842</td>
      <td>3.900622</td>
      <td>-11.798421</td>
      <td>12.426487</td>
      <td>198.294167</td>
    </tr>
    <tr>
      <th>12079</th>
      <td>POINT (389000 5791040)</td>
      <td>12079</td>
      <td>8900</td>
      <td>10900</td>
      <td>389000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.415643</td>
      <td>-1.148986</td>
      <td>4.156429</td>
      <td>-11.489856</td>
      <td>12.218539</td>
      <td>199.887474</td>
    </tr>
    <tr>
      <th>12080</th>
      <td>POINT (390000 5791040)</td>
      <td>12080</td>
      <td>9000</td>
      <td>10900</td>
      <td>390000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.404408</td>
      <td>-1.165013</td>
      <td>4.044082</td>
      <td>-11.650125</td>
      <td>12.332073</td>
      <td>199.143308</td>
    </tr>
    <tr>
      <th>12081</th>
      <td>POINT (391000 5791040)</td>
      <td>12081</td>
      <td>9100</td>
      <td>10900</td>
      <td>391000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.425620</td>
      <td>-1.187113</td>
      <td>4.256201</td>
      <td>-11.871128</td>
      <td>12.611064</td>
      <td>199.724473</td>
    </tr>
    <tr>
      <th>12082</th>
      <td>POINT (392000 5791040)</td>
      <td>12082</td>
      <td>9200</td>
      <td>10900</td>
      <td>392000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.412727</td>
      <td>-1.140732</td>
      <td>4.127267</td>
      <td>-11.407319</td>
      <td>12.131004</td>
      <td>199.890562</td>
    </tr>
    <tr>
      <th>12083</th>
      <td>POINT (393000 5791040)</td>
      <td>12083</td>
      <td>9300</td>
      <td>10900</td>
      <td>393000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.420770</td>
      <td>-1.208876</td>
      <td>4.207698</td>
      <td>-12.088759</td>
      <td>12.800110</td>
      <td>199.191321</td>
    </tr>
    <tr>
      <th>12084</th>
      <td>POINT (394000 5791040)</td>
      <td>12084</td>
      <td>9400</td>
      <td>10900</td>
      <td>394000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.411952</td>
      <td>-1.213787</td>
      <td>4.119517</td>
      <td>-12.137869</td>
      <td>12.817889</td>
      <td>198.746894</td>
    </tr>
    <tr>
      <th>12085</th>
      <td>POINT (395000 5791040)</td>
      <td>12085</td>
      <td>9500</td>
      <td>10900</td>
      <td>395000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.442422</td>
      <td>-0.854839</td>
      <td>4.424218</td>
      <td>-8.548387</td>
      <td>9.625416</td>
      <td>207.363827</td>
    </tr>
    <tr>
      <th>12086</th>
      <td>POINT (396000 5791040)</td>
      <td>12086</td>
      <td>9600</td>
      <td>10900</td>
      <td>396000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.410208</td>
      <td>-1.201714</td>
      <td>4.102077</td>
      <td>-12.017142</td>
      <td>12.697982</td>
      <td>198.847450</td>
    </tr>
    <tr>
      <th>12087</th>
      <td>POINT (397000 5791040)</td>
      <td>12087</td>
      <td>9700</td>
      <td>10900</td>
      <td>397000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.394172</td>
      <td>-1.189054</td>
      <td>3.941716</td>
      <td>-11.890545</td>
      <td>12.526858</td>
      <td>198.340361</td>
    </tr>
    <tr>
      <th>12088</th>
      <td>POINT (398000 5791040)</td>
      <td>12088</td>
      <td>9800</td>
      <td>10900</td>
      <td>398000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.372462</td>
      <td>-1.105048</td>
      <td>3.724623</td>
      <td>-11.050478</td>
      <td>11.661299</td>
      <td>198.626667</td>
    </tr>
    <tr>
      <th>12089</th>
      <td>POINT (399000 5791040)</td>
      <td>12089</td>
      <td>9900</td>
      <td>10900</td>
      <td>399000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.222491</td>
      <td>-1.150124</td>
      <td>2.224915</td>
      <td>-11.501243</td>
      <td>11.714471</td>
      <td>190.948626</td>
    </tr>
    <tr>
      <th>12090</th>
      <td>POINT (400000 5791040)</td>
      <td>12090</td>
      <td>10000</td>
      <td>10900</td>
      <td>400000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.228518</td>
      <td>-1.205668</td>
      <td>2.285181</td>
      <td>-12.056680</td>
      <td>12.271332</td>
      <td>190.732335</td>
    </tr>
    <tr>
      <th>12091</th>
      <td>POINT (401000 5791040)</td>
      <td>12091</td>
      <td>10100</td>
      <td>10900</td>
      <td>401000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.179267</td>
      <td>-1.124293</td>
      <td>1.792668</td>
      <td>-11.242928</td>
      <td>11.384950</td>
      <td>189.059463</td>
    </tr>
    <tr>
      <th>12092</th>
      <td>POINT (402000 5791040)</td>
      <td>12092</td>
      <td>10200</td>
      <td>10900</td>
      <td>402000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.249695</td>
      <td>-1.089078</td>
      <td>2.496946</td>
      <td>-10.890779</td>
      <td>11.173353</td>
      <td>192.913120</td>
    </tr>
    <tr>
      <th>12093</th>
      <td>POINT (403000 5791040)</td>
      <td>12093</td>
      <td>10300</td>
      <td>10900</td>
      <td>403000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.244407</td>
      <td>-0.743836</td>
      <td>2.444071</td>
      <td>-7.438364</td>
      <td>7.829607</td>
      <td>198.189303</td>
    </tr>
    <tr>
      <th>12094</th>
      <td>POINT (404000 5791040)</td>
      <td>12094</td>
      <td>10400</td>
      <td>10900</td>
      <td>404000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.271931</td>
      <td>-0.778958</td>
      <td>2.719314</td>
      <td>-7.789583</td>
      <td>8.250592</td>
      <td>199.243899</td>
    </tr>
    <tr>
      <th>12095</th>
      <td>POINT (405000 5791040)</td>
      <td>12095</td>
      <td>10500</td>
      <td>10900</td>
      <td>405000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.227180</td>
      <td>-0.774121</td>
      <td>2.271795</td>
      <td>-7.741213</td>
      <td>8.067679</td>
      <td>196.355253</td>
    </tr>
    <tr>
      <th>12096</th>
      <td>POINT (406000 5791040)</td>
      <td>12096</td>
      <td>10600</td>
      <td>10900</td>
      <td>406000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.225044</td>
      <td>-0.771500</td>
      <td>2.250443</td>
      <td>-7.715001</td>
      <td>8.036525</td>
      <td>196.261809</td>
    </tr>
    <tr>
      <th>12097</th>
      <td>POINT (407000 5791040)</td>
      <td>12097</td>
      <td>10700</td>
      <td>10900</td>
      <td>407000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.307129</td>
      <td>-0.782714</td>
      <td>3.071287</td>
      <td>-7.827140</td>
      <td>8.408146</td>
      <td>201.424520</td>
    </tr>
    <tr>
      <th>12098</th>
      <td>POINT (408000 5791040)</td>
      <td>12098</td>
      <td>10800</td>
      <td>10900</td>
      <td>408000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.284867</td>
      <td>-0.810308</td>
      <td>2.848671</td>
      <td>-8.103080</td>
      <td>8.589228</td>
      <td>199.369333</td>
    </tr>
    <tr>
      <th>12099</th>
      <td>POINT (409000 5791040)</td>
      <td>12099</td>
      <td>10900</td>
      <td>10900</td>
      <td>409000.0</td>
      <td>5791040.0</td>
      <td>207</td>
      <td>0.304630</td>
      <td>-0.859407</td>
      <td>3.046297</td>
      <td>-8.594072</td>
      <td>9.118004</td>
      <td>199.517633</td>
    </tr>
  </tbody>
</table>
<p>8035 rows × 13 columns</p>
</div>



#### calculate geometric quality grid - without any disk access

All you have to do is to instanciate Geom_Quality_Grid with two instances of the GeoArray class as described above.


```python
GQG = Geom_Quality_Grid(GeoArray(ref_ndarray, ref_gt, ref_prj),GeoArray(tgt_ndarray, tgt_gt, tgt_prj),**kwargs)
GQG.get_quality_grid()
```

#### export geometric quality grid to an ESRI point shapefile


```python
GQG.quality_grid_to_PointShapefile()
```

### Shell console interface

By far, there is no shell console interface for this module.
