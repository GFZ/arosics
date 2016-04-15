import math
import os
import re
import warnings
import numpy  as np
import pandas as pd

try:
    import gdal
    import osr
    import ogr
except ImportError:
    from osgeo import gdal
    from osgeo import osr
    from osgeo import ogr
import rasterio
from rasterio.warp import RESAMPLING
from rasterio.warp import calculate_default_transform as rio_calc_transform
from rasterio.warp import reproject as rio_reproject
from shapely.geometry import Polygon, shape

from .utilities import get_dtypeStr
from .io        import gdal_ReadAsArray_mp



def get_corner_coordinates(fPath=None, gt=None, rows=None,cols=None):
    """Return (UL, LL, LR, UR)"""
    if fPath is None: assert None not in [gt,rows,cols],"'gt', 'rows' and 'cols' must be given if gdal_ds is missing."
    if fPath:
        gdal_ds = gdal.Open(fPath)
        gdal_ds_GT = gdal_ds.GetGeoTransform()
        rows, cols = gdal_ds.RasterYSize, gdal_ds.RasterXSize
        gdal_ds = None
    else:
        gdal_ds_GT = gt
        rows, cols = rows,cols
    ext=[]
    xarr=[0,cols]
    yarr=[0,rows]
    for px in xarr:
        for py in yarr:
            x=gdal_ds_GT[0]+(px*gdal_ds_GT[1])+(py*gdal_ds_GT[2])
            y=gdal_ds_GT[3]+(px*gdal_ds_GT[4])+(py*gdal_ds_GT[5])
            ext.append([x,y])
        yarr.reverse()
    return ext


def corner_coord_to_minmax(corner_coords):
    """Return (UL, LL, LR, UR)"""
    x_vals = [int(i[0]) for i in corner_coords] # FIXME int?
    y_vals = [int(i[1]) for i in corner_coords] # FIXME int?
    xmin,xmax,ymin,ymax = min(x_vals),max(x_vals),min(y_vals),max(y_vals)
    return xmin,xmax,ymin,ymax


def get_footprint_polygon(CornerLonLat):
    outpoly = Polygon(CornerLonLat)
    assert outpoly. is_valid, 'The given coordinates result in an invalid polygon. Check coordinate order.'
    return outpoly


def get_overlap_polygon(poly1, poly2,v=0):
    overlap_poly = poly1.intersection(poly2)
    if not overlap_poly.is_empty:
        overlap_percentage = 100 * shape(overlap_poly).area / shape(poly2).area
        if v: print('%.2f percent of the image to be shifted is covered by the reference image.' %overlap_percentage)
        return {'overlap poly':overlap_poly, 'overlap percentage':overlap_percentage, 'overlap area':overlap_poly.area}
    else:
        return {'overlap poly':None, 'overlap percentage':0, 'overlap area':0}


def get_true_corner_lonlat(fPath,bandNr=1,noDataVal=None,mp=1,v=0):
    gdal_ds   = gdal.Open(fPath)
    fPath     = gdal_ds.GetDescription()
    rows,cols = gdal_ds.RasterYSize, gdal_ds.RasterXSize
    gt,prj    = gdal_ds.GetGeoTransform(),gdal_ds.GetProjection()
    gdal_ds = None

    mask_1bit = np.zeros((rows,cols),dtype='uint8') # zeros -> image area later overwritten by ones

    if noDataVal is None:
        mask_1bit[:,:] = 1
    elif noDataVal=='unclear':
        warnings.warn("No data value could not be automatically detected. Thus the matching window used for shift "
              "calulation had to be centered in the middle of the overlap center without respecting no data values. "
              "To avoid this provide the correct no data values for reference and shift image via '-nodata'")
        mask_1bit[:,:] = 1
    else:
        if mp:
            band_data = gdal_ReadAsArray_mp(fPath,bandNr)
        else:
            gdal.Open(fPath)
            band_data = gdal_ds.GetRasterBand(bandNr).ReadAsArray()
            gdal_ds = None
        mask_1bit[band_data!=noDataVal] = 1

    if v: print('detected no data value',noDataVal)

    cols_containing_data = np.any(mask_1bit, axis=0)
    rows_containing_data = np.any(mask_1bit, axis=1)

    first_dataCol = list(cols_containing_data).index(True)
    last_dataCol  = cols - (list(reversed(cols_containing_data)).index(True) + 1)
    first_dataRow = list(rows_containing_data).index(True)
    last_dataRow  = rows - (list(reversed(rows_containing_data)).index(True) + 1)

    if first_dataCol == 0 and first_dataRow == 0 and last_dataCol == cols-1 and last_dataRow == rows-1 \
            and mask_1bit[0, 0] == True and mask_1bit[0, cols - 1] == True \
            and mask_1bit[rows - 1, 0] == True and mask_1bit[rows - 1, cols - 1] == True:
        UL, UR, LL, LR = (0, 0), (0, cols - 1), (rows - 1, 0), (rows - 1, cols - 1)
    else:
        StartStopRows_in_first_dataCol = [list(mask_1bit[:, first_dataCol]).index(True),
                                         (rows - list(reversed(mask_1bit[:, first_dataCol])).index(True) + 1)]
        StartStopRows_in_last_dataCol  = [list(mask_1bit[:, last_dataCol]).index(True),
                                         (rows - list(reversed(mask_1bit[:, last_dataCol])).index(True) + 1)]
        StartStopCols_in_first_dataRow = [list(mask_1bit[first_dataRow, :]).index(True),
                                         (cols - list(reversed(mask_1bit[first_dataRow, :])).index(True) + 1)]
        StartStopCols_in_last_dataRow  = [list(mask_1bit[last_dataRow, :]).index(True),
                                         (cols - list(reversed(mask_1bit[last_dataRow, :])).index(True) + 1)]

        if True in [abs(np.diff(i)[0]) > 10 for i in [StartStopRows_in_first_dataCol, StartStopRows_in_last_dataCol,
                                                      StartStopCols_in_first_dataRow, StartStopCols_in_last_dataRow]]:
            ''' In case of cut image corners (e.g. ALOS AVNIR-2 Level-1B2):
             Calculation of trueDataCornerPos outside of the image.'''
            def find_line_intersection_point(line1, line2):
                xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
                ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

                det = lambda a,b: a[0] * b[1] - a[1] * b[0]
                div = det(xdiff, ydiff)
                if div == 0: # lines do not intersect
                    return None, None
                d = (det(*line1), det(*line2))
                x = det(d, xdiff) / div
                y = det(d, ydiff) / div
                return x, y

            line_N = ((StartStopCols_in_first_dataRow[1], first_dataRow),
                      (last_dataCol, StartStopRows_in_last_dataCol[0]))
            line_S = ((first_dataCol, StartStopRows_in_first_dataCol[1]),
                      (StartStopCols_in_last_dataRow[0], last_dataRow))
            line_W = ((StartStopCols_in_first_dataRow[0], first_dataRow),
                      (first_dataCol, StartStopRows_in_first_dataCol[0]))
            line_O = ((last_dataCol, StartStopRows_in_last_dataCol[1]),
                      (StartStopCols_in_last_dataRow[1], last_dataRow))
            corners = [list(reversed(find_line_intersection_point(line_N, line_W))),
                       list(reversed(find_line_intersection_point(line_N, line_O))),
                       list(reversed(find_line_intersection_point(line_S, line_W))),
                       list(reversed(find_line_intersection_point(line_S, line_O)))]
        else:
            dataRow_in_first_dataCol = np.mean(StartStopRows_in_first_dataCol)
            dataRow_in_last_dataCol  = np.mean(StartStopRows_in_last_dataCol)
            dataCol_in_first_dataRow = np.mean(StartStopCols_in_first_dataRow)
            dataCol_in_last_dataRow  = np.mean(StartStopCols_in_last_dataRow)
            corners = [(first_dataRow, dataCol_in_first_dataRow), (dataRow_in_last_dataCol, last_dataCol),
                       (last_dataRow, dataCol_in_last_dataRow), (dataRow_in_first_dataCol, first_dataCol)]

        if not [None, None] in corners:
            distUL = [np.sqrt((0 - c[0]) ** 2 + (0 - c[1]) ** 2) for c in corners]       # distance of each corner to UL
            UL     = corners[distUL.index(min(distUL))]
            distUR = [np.sqrt((0 - c[0]) ** 2 + (cols - c[1]) ** 2) for c in corners]    # distance of each corner to UR
            UR     = corners[distUR.index(min(distUR))]
            distLL = [np.sqrt((rows - c[0]) ** 2 + (0 - c[1]) ** 2) for c in corners]    # distance of each corner to LL
            LL     = corners[distLL.index(min(distLL))]
            distLR = [np.sqrt((rows - c[0]) ** 2 + (cols - c[1]) ** 2) for c in corners] # distance of each corner to LR
            LR     = corners[distLR.index(min(distLR))]
        else:
            UL, UR, LL, LR = (0, 0), (0, cols - 1), (rows - 1, 0), (rows - 1, cols - 1)
            print('')
            warnings.warn('Automatic detection of the actual image corners of %s failed, because not all of the four '
                          'corners could be found within the image. The current detection algorithm assumes that the '
                          'actual image corners are within the image. Using the outer image corners instead.\n You can '
                          'use -cor0 and -cor1 to manually set the image corner coordinates or -v to check if the '
                            'the matching window is positioned correcly.\n' %os.path.basename(fPath))

    # check if all points are unique
    all_coords_are_unique = len([UL, UR, LL, LR]) == len(pd.DataFrame([UL, UR, LL, LR]).drop_duplicates().values)
    UL, UR, LL, LR = (UL, UR, LL, LR) if all_coords_are_unique else ((0, 0), (0, cols-1), (rows-1, 0), (rows-1, cols-1))

    get_mapYX = lambda YX: pixelToMapYX(list(reversed(YX)),geotransform=gt,projection=prj)[0]
    corner_pos_XY = [list(reversed(i)) for i in [get_mapYX(YX) for YX in [UL, UR, LR, LL]]]
    return corner_pos_XY


mapYX2imYX = lambda mapYX,gt: ((mapYX[0]-gt[3])/gt[5],  (mapYX[1]-gt[0])/gt[1])
imYX2mapYX = lambda imYX ,gt: (gt[3]-(imYX[0]*abs(gt[5])), gt[0]+(imYX[1]*gt[1]))


def get_smallest_boxImYX_that_contains_boxMapYX(box_mapYX, gt_im):
    xmin,ymin,xmax,ymax = Polygon([(i[1],i[0]) for i in box_mapYX]).bounds # map-bounds box_mapYX
    (ymin,xmin),(ymax,xmax) = mapYX2imYX([ymin,xmin],gt_im), mapYX2imYX([ymax,xmax],gt_im) # image coord bounds
    xmin,ymin,xmax,ymax = int(xmin),math.ceil(ymin), math.ceil(xmax), int(ymax) # round min coords off and max coords on
    return (ymax,xmin),(ymax,xmax),(ymin,xmax),(ymin,xmin) # UL_YX,UR_YX,LR_YX,LL_YX


def get_subset_GeoTransform(gt_fullArr,subset_box_imYX):
    gt_subset = list(gt_fullArr[:]) # copy
    gt_subset[3],gt_subset[0] = imYX2mapYX(subset_box_imYX[0],gt_fullArr)
    return gt_subset


def get_gdalReadInputs_from_boxImYX(boxImYX):
    """Returns row_start,col_start,rows_count,cols_count and assumes boxImYX as [UL_YX,UR_YX,LR_YX,LL_YX)"""
    rS, cS = boxImYX[0]
    clip_sz_x = abs(boxImYX[1][1]-boxImYX[0][1]) # URx-ULx
    clip_sz_y = abs(boxImYX[0][0]-boxImYX[2][0]) # ULy-LLy
    return cS, rS, clip_sz_x,clip_sz_y


def write_shp(path_out, shapely_geom, prj=None,attrDict=None):
    shapely_geom = [shapely_geom] if not isinstance(shapely_geom,list) else shapely_geom
    attrDict     = [attrDict] if not isinstance(attrDict,list) else attrDict
    assert len(shapely_geom) == len(attrDict), "'shapely_geom' and 'attrDict' must have the same length."

    print('Writing %s ...' %path_out)
    if os.path.exists(path_out): os.remove(path_out)
    ds = ogr.GetDriverByName("Esri Shapefile").CreateDataSource(path_out)

    if prj is not None:
        srs = osr.SpatialReference()
        srs.ImportFromWkt(prj)
    else:
        srs = None

    geom_type = list(set([gm.type for gm in shapely_geom]))
    assert len(geom_type)==1, 'All shapely geometries must belong to the same type. Got %s.' %geom_type

    layer = ds.CreateLayer('', srs, ogr.wkbPoint)      if geom_type[0] == 'Point'      else \
            ds.CreateLayer('', srs, ogr.wkbLineString) if geom_type[0] == 'LineString' else \
            ds.CreateLayer('', srs, ogr.wkbPolygon)    if geom_type[0] == 'Polygon'    else \
            None # FIXME

    if isinstance(attrDict[0],dict):
        for attr in attrDict[0].keys():
            assert len(attr)<=10, "ogr does not support fieldnames longer than 10 digits. '%s' is too long" %attr
            DTypeStr  = get_dtypeStr(attrDict[0][attr])
            FieldType = ogr.OFTInteger if DTypeStr.startswith('int') else ogr.OFTReal     if DTypeStr.startswith('float') else \
                        ogr.OFTString  if DTypeStr.startswith('str') else ogr.OFTDateTime if DTypeStr.startswith('date') else None
            FieldDefn = ogr.FieldDefn(attr, FieldType)
            if DTypeStr.startswith('float'):
                FieldDefn.SetPrecision(6)
            layer.CreateField(FieldDefn)  # Add one attribute

    for i in range(len(shapely_geom)):
        # Create a new feature (attribute and geometry)
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometry(ogr.CreateGeometryFromWkb(shapely_geom[i].wkb)) # Make a geometry, from Shapely object

        list_attr2set = attrDict[0].keys() if isinstance(attrDict[0],dict) else []

        for attr in list_attr2set:
            val       = attrDict[i][attr]
            DTypeStr  = get_dtypeStr(val)
            val = int(val) if DTypeStr.startswith('int') else float(val) if DTypeStr.startswith('float') else \
                  str(val) if DTypeStr.startswith('str') else val
            feat.SetField(attr, val)

        layer.CreateFeature(feat)
        feat.Destroy()

    # Save and close everything
    ds = layer = None


def get_proj4info(ds=None,proj=None):
    assert ds is not None or proj is not None, "Specify at least one of the arguments 'ds' or 'proj'"
    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjection() if ds is not None else proj)
    return srs.ExportToProj4()


def wkt2epsg(wkt, epsg=os.environ['GDAL_DATA'].replace('/gdal','/proj/epsg')):
    """ Transform a WKT string to an EPSG code
    :param wkt:  WKT definition
    :param epsg: the proj.4 epsg file (defaults to '/usr/local/share/proj/epsg')
    :returns:    EPSG code
    http://gis.stackexchange.com/questions/20298/is-it-possible-to-get-the-epsg-value-from-an-osr-spatialreference-class-using-th
    """
    code = None #default
    p_in = osr.SpatialReference()
    s = p_in.ImportFromWkt(wkt)
    if s == 5:  # invalid WKT
        raise Exception('Invalid WKT.')
    if p_in.IsLocal():
        raise Exception('The given WKT is a local coordinate system.')
    cstype = 'GEOGCS' if p_in.IsGeographic() else 'PROJCS'
    an = p_in.GetAuthorityName(cstype)
    assert an in [None,'EPSG'], "No EPSG code found. Found %s instead." %an
    ac = p_in.GetAuthorityCode(cstype)
    if ac is None:  # try brute force approach by grokking proj epsg definition file
        p_out = p_in.ExportToProj4()
        if p_out:
            with open(epsg) as f:
                for line in f:
                    if line.find(p_out) != -1:
                        m = re.search('<(\\d+)>', line)
                        if m:
                            code = m.group(1)
                            break
                code = code if code else None # match or no match
    else:
        code = ac
    return code


def pixelToLatLon(pixelPairs, path_im=None, geotransform=None, projection=None):
    """The following method translates given pixel locations into latitude/longitude locations on a given GEOTIF
    This method does not take into account pixel size and assumes a high enough
    image resolution for pixel size to be insignificant
    :param pixelPairs:      The pixel pairings to be translated in the form [[x1,y1],[x2,y2]]
    :param path_im:         The file location of the input image
    :param geotransform:    GDAL GeoTransform
    :param projection:      GDAL Projection
    :returns:               The lat/lon translation of the pixel pairings in the form [[lat1,lon1],[lat2,lon2]]
    """
    assert path_im is not None or geotransform is not None and projection is not None, \
        "GEOP.pixelToLatLon: Missing argument! Please provide either 'path_im' or 'geotransform' AND 'projection'."
    if geotransform is not None and projection is not None:
        gt,proj = geotransform, projection
    else:
        ds       = gdal.Open(path_im)
        gt, proj = ds.GetGeoTransform(), ds.GetProjection()

    # Create a spatial reference object for the dataset
    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj)
    # Set up the coordinate transformation object
    srsLatLong = srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(srs, srsLatLong)
    # Go through all the point pairs and translate them to pixel pairings
    latLonPairs = []
    for point in pixelPairs:
        # Translate the pixel pairs into untranslated points
        ulon = point[0] * gt[1] + gt[0]
        ulat = point[1] * gt[5] + gt[3]
        # Transform the points to the space
        (lon, lat, holder) = ct.TransformPoint(ulon, ulat)
        # Add the point to our return array
        latLonPairs.append([lat, lon])

    return latLonPairs


def pixelToMapYX(pixelCoords, path_im=None, geotransform=None, projection=None):
    assert path_im is not None or geotransform is not None and projection is not None, \
        "pixelToMapYX: Missing argument! Please provide either 'path_im' or 'geotransform' AND 'projection'."
    if geotransform is not None and projection is not None:
        gt,proj = geotransform, projection
    else:
        ds       = gdal.Open(path_im)
        gt, proj = ds.GetGeoTransform(), ds.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj)
    # Set up the coordinate transformation object
    ct  = osr.CoordinateTransformation(srs, srs)

    mapYmapXPairs = []
    pixelCoords   = [pixelCoords] if not type(pixelCoords[0]) in [list,tuple] else pixelCoords

    for point in pixelCoords:
        # Translate the pixel pairs into untranslated points
        u_mapX = point[0] * gt[1] + gt[0]
        u_mapY = point[1] * gt[5] + gt[3]
        # Transform the points to the space
        (mapX, mapY, holder) = ct.TransformPoint(u_mapX, u_mapY)
        # Add the point to our return array
        mapYmapXPairs.append([mapY, mapX])

    return mapYmapXPairs


def isProjectedOrGeographic(prj):
    srs = osr.SpatialReference()
    if prj.startswith('EPSG:'):
        srs.ImportFromEPSG(int(prj.split(':')[1]))
    elif prj.startswith('+proj='):
        srs.ImportFromProj4(prj)
    elif prj.startswith('GEOGCS') or prj.startswith('PROJCS'):
        srs.ImportFromWkt(prj)
    else:
        raise Exception('Unknown input projection.')
    return 'projected' if srs.IsProjected() else 'geographic' if srs.IsGeographic() else None


def warp_ndarray(ndarray,in_gt, in_prj, out_prj, out_gt=None, outRowsCols=None, outUL=None,
                 out_res=None,out_extent=None,out_dtype=None,rsp_alg=0,in_nodata=None,out_nodata=None):
    """Reproject / warp a numpy array with given geo information to target coordinate system.
    :param ndarray:     numpy.ndarray [rows,cols,bands]
    :param in_gt:       input gdal GeoTransform
    :param in_prj:      input projection as WKT string
    :param out_prj:     output projection as WKT string
    :param out_gt:      output gdal GeoTransform as float tuple in the source coordinate system (optional)
    :param outUL:       [X,Y] output upper left coordinates as floats in the source coordinate system
                        (requires outRowsCols)
    :param outRowsCols: [rows, cols] (optional)
    :param out_res:     output resolution as tuple of floats (x,y) in the TARGET coordinate system
    :param out_extent:  [left, bottom, right, top] as floats in the source coordinate system
    :param out_dtype:   output data type as numpy data type
    :param rsp_alg:     Resampling method to use. One of the following (int, default is 0):
                        0 = nearest neighbour, 1 = bilinear, 2 = cubic, 3 = cubic spline, 4 = lanczos,
                        5 = average, 6 = mode
    :param in_nodata    no data value of the input image
    :param out_nodata   no data value of the output image
    :return out_arr:    warped numpy array
    :return out_gt:     warped gdal GeoTransform
    :return out_prj:    warped projection as WKT string
    """

    if not ndarray.flags['OWNDATA']:
        temp    = np.empty_like(ndarray)
        temp[:] = ndarray
        ndarray = temp # deep copy: converts view to its own array in oder to avoid wrong output

    with rasterio.drivers():
        if out_gt is None:
            assert out_extent is not None, 'Provide eighter out_gt or out_extent!'
            if outUL is not None:
                assert outRowsCols is not None, 'outRowsCols must be given if outUL is given.'
            outUL = [in_gt[0],in_gt[3]] if outUL is None else outUL

        inEPSG, outEPSG = [wkt2epsg(prj) for prj in [in_prj, out_prj]]
        assert inEPSG  is not None, 'Could not derive input EPSG code.'
        assert outEPSG is not None, 'Could not derive output EPSG code.'

        src_crs = {'init': 'EPSG:%s' %inEPSG}
        dst_crs = {'init': 'EPSG:%s' %outEPSG}

        if len(ndarray.shape)==3:
            # convert input array axis order to rasterio axis order
            ndarray = np.swapaxes(np.swapaxes(ndarray,0,2),1,2)
            bands,rows,cols = ndarray.shape
            rows,cols       = outRowsCols if outRowsCols else rows,cols
        else:
            bands = 0
            rows,cols = ndarray.shape if outRowsCols is None else outRowsCols

        out_dtype = ndarray.dtype if out_dtype is None else out_dtype
        gt2bounds = lambda gt,r,c: [gt[0], gt[3] + r*gt[5], gt[0] + c*gt[1], gt[3]] # left, bottom, right, top

        #get dst_transform
        if out_gt is None and out_extent is None:
            if outRowsCols:
                outUL = [in_gt[0],in_gt[3]] if outUL is None else outUL
                ulRC2bounds = lambda ul,r,c: [ul[0], ul[1] + r*in_gt[5], ul[0] + c*in_gt[1], ul[1]] # left, bottom, right, top
                left, bottom, right, top = ulRC2bounds(outUL,rows,cols)
            else: # outRowsCols is None and outUL is None: use in_gt
                left, bottom, right, top = gt2bounds(in_gt,rows,cols)
                    #,im_xmax,im_ymin,im_ymax = corner_coord_to_minmax(get_corner_coordinates(self.ds_im2shift))
        elif out_extent:
            left, bottom, right, top = out_extent
            #input array is only a window of the actual input array
            assert left >= in_gt[0] and right <= (in_gt[0]+(cols+1)*in_gt[1]) and \
                   bottom >= in_gt[3]+(rows+1)*in_gt[5] and top <= in_gt[3], \
                   "out_extent has to be completely within the input image bounds."
        else: # out_gt is given
            left, bottom, right, top = gt2bounds(in_gt,rows,cols)

        if out_res is None:
            # get pixel resolution in target coord system
            prj_in_out = (isProjectedOrGeographic(in_prj),isProjectedOrGeographic(out_prj))
            assert None not in prj_in_out, 'prj_in_out contains None.'
            if prj_in_out[0] == prj_in_out[1]:
                out_res = (in_gt[1],abs(in_gt[5]))
            elif prj_in_out == ('geographic','projected'):
                raise NotImplementedError
            else:             #('projected','geographic')
                px_size_LatLon = np.array(pixelToLatLon([1,1],geotransform=in_gt,projection=in_prj)) - \
                                 np.array(pixelToLatLon([0,0],geotransform=in_gt,projection=in_prj))
                out_res = tuple(reversed(abs(px_size_LatLon)))
                print('OUT_RES NOCHMAL CHECKEN: ', out_res)

        dst_transform,out_cols,out_rows = rio_calc_transform(
            src_crs, dst_crs, cols, rows, left, bottom, right, top, resolution=out_res, densify_pts=21)

        aff = list(dst_transform)
        if out_gt is None:
            out_gt = (aff[2],aff[0],aff[1],aff[5],aff[3],aff[4])

        #src_transform = rasterio.transform.from_origin(in_gt[0], in_gt[3], in_gt[1], in_gt[5])

        dict_rspInt_rspAlg = {0:RESAMPLING.nearest, 1:RESAMPLING.bilinear, 2:RESAMPLING.cubic,
                              3:RESAMPLING.cubic_spline, 4:RESAMPLING.lanczos, 5:RESAMPLING.average, 6:RESAMPLING.mode}

        out_arr = np.zeros((bands,out_rows,out_cols), out_dtype)\
            if len(ndarray.shape)==3 else np.zeros((out_rows,out_cols), out_dtype)

        # FIXME direct passing of src_transform and dst_transform results in a wrong output image. Maybe a rasterio-bug?
        #rio_reproject(ndarray, out_arr, src_transform=src_transform, src_crs=src_crs, dst_transform=dst_transform,
        #    dst_crs=dst_crs, resampling=dict_rspInt_rspAlg[rsp_alg])
        # FIXME indirect passing causes Future warning
        with warnings.catch_warnings():
            warnings.simplefilter('ignore') # FIXME supresses: FutureWarning: GDAL-style transforms are deprecated and will not be supported in Rasterio 1.0.
            rio_reproject(ndarray, out_arr,
                      src_transform=in_gt, src_crs=src_crs, dst_transform=out_gt, dst_crs=dst_crs,
                      resampling=dict_rspInt_rspAlg[rsp_alg], src_nodata=in_nodata, dst_nodata=out_nodata)

        if len(ndarray.shape)==3:
            # convert output array axis order to GMS axis order [rows,cols,bands]
            out_arr = np.swapaxes(np.swapaxes(out_arr,0,1),1,2)

    return out_arr, out_gt, out_prj


def find_noDataVal(path_im,bandNr=1,sz=3,is_vrt=0):
    """tries to derive no data value from homogenious corner pixels within 3x3 windows (by default)
    :param gdal_ds:
    :param bandNr:
    :param sz: window size in which corner pixels are analysed
    """
    gdal_ds = gdal.Open(path_im) if not is_vrt else gdal.OpenShared(path_im)
    get_mean_std = lambda corner_subset: {'mean':np.mean(corner_subset), 'std':np.std(corner_subset)}
    rows,cols = gdal_ds.RasterYSize,gdal_ds.RasterXSize
    UL = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(0,0,sz,sz))
    UR = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(cols-sz,0,sz,sz))
    LR = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(cols-sz,rows-sz,sz,sz))
    LL = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(0,rows-sz,sz,sz))
    gdal_ds = None
    possVals  = [i['mean'] for i in [UL,UR,LR,LL] if i['std']==0]
    # possVals==[]: all corners are filled with data; np.std(possVals)==0: noDataVal clearly identified
    return None if possVals==[] else possVals[0] if np.std(possVals)==0 else 'unclear'