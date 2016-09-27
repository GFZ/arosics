# -*- coding: utf-8 -*-
import math
import os
import re
import warnings

import numpy  as np
import pyproj

try:
    import gdal
    import osr
    import ogr
except ImportError:
    from osgeo import gdal
    from osgeo import osr
    from osgeo import ogr

from geopandas import GeoDataFrame


# custom
from py_tools_ds.ptds.geo.vector.geometry  import boxObj
from py_tools_ds.ptds.geo.vector.geometry  import round_shapelyPoly_coords
from py_tools_ds.ptds.geo.vector.topology  import get_footprint_polygon, get_overlap_polygon, \
                                                  find_line_intersection_point, get_largest_onGridPoly_within_poly, \
                                                  get_smallest_boxImYX_that_contains_boxMapYX, \
                                                  get_smallest_shapelyImPolyOnGrid_that_contains_shapelyImPoly
from py_tools_ds.ptds.geo.coord_calc       import get_corner_coordinates, corner_coord_to_minmax
from py_tools_ds.ptds.geo.coord_grid       import move_shapelyPoly_to_image_grid, find_nearest_grid_coord
from py_tools_ds.ptds.geo.coord_trafo      import transform_utm_to_wgs84, pixelToLatLon, pixelToMapYX, imYX2mapYX
from py_tools_ds.ptds.geo.projection       import get_proj4info, get_UTMzone, WKT2EPSG, isProjectedOrGeographic
from py_tools_ds.ptds.geo.map_info         import geotransform2mapinfo
from py_tools_ds.ptds.geo.raster.reproject import warp_ndarray


from . import io as IO


__all__=['boxObj',
         'round_shapelyPoly_coords',
         'get_footprint_polygon',
         'get_overlap_polygon',
         'find_line_intersection_point',
         'get_largest_onGridPoly_within_poly',
         'get_smallest_boxImYX_that_contains_boxMapYX',
         'get_smallest_shapelyImPolyOnGrid_that_contains_shapelyImPoly',
         'get_corner_coordinates',
         'corner_coord_to_minmax,',
         'move_shapelyPoly_to_image_grid',
         'find_nearest_grid_coord',
         'transform_utm_to_wgs84',
         'pixelToLatLon',
         'pixelToMapYX',
         'imYX2mapYX',
         'get_proj4info',
         'get_UTMzone',
         'WKT2EPSG',
         'isProjectedOrGeographic',
         'geotransform2mapinfo',
         'warp_ndarray']





def angle_to_north(XY):
    """Calculates the angle between the lines [origin:[0,0],north:[0,1]] and
     [origin:[0,0],pointXY:[X,Y]] in clockwise direction. Returns values between 0 and 360 degrees."""
    XY    = np.array(XY)
    XYarr = XY if len(XY.shape)==2 else XY.reshape((1,2))
    return np.abs(np.degrees(np.arctan2(XYarr[:,1],XYarr[:,0])-np.pi/2)%360)


def get_true_corner_lonlat(fPath,bandNr=1,noDataVal=None,mp=1,v=0): # TODO implement shapely algorithm
    gdal_ds   = gdal.Open(fPath)
    fPath     = gdal_ds.GetDescription()
    rows,cols = gdal_ds.RasterYSize, gdal_ds.RasterXSize
    gt,prj    = gdal_ds.GetGeoTransform(),gdal_ds.GetProjection()
    gdal_ds   = None

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
            band_data = IO.gdal_ReadAsArray_mp(fPath, bandNr)
        else:
            gdal_ds   = gdal.Open(fPath)
            band_data = gdal_ds.GetRasterBand(bandNr).ReadAsArray()
            gdal_ds   = None
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
    all_coords_are_unique = len([UL, UR, LL, LR]) == len(GeoDataFrame([UL, UR, LL, LR]).drop_duplicates().values)
    UL, UR, LL, LR = (UL, UR, LL, LR) if all_coords_are_unique else ((0, 0), (0, cols-1), (rows-1, 0), (rows-1, cols-1))

    get_mapYX = lambda YX: pixelToMapYX(list(reversed(YX)),geotransform=gt,projection=prj)[0]
    corner_pos_XY = [list(reversed(i)) for i in [get_mapYX(YX) for YX in [UL, UR, LR, LL]]]
    return corner_pos_XY


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


def find_noDataVal(path_im,bandNr=1,sz=3,is_vrt=0):
    """tries to derive no data value from homogenious corner pixels within 3x3 windows (by default)
    :param path_im:
    :param bandNr:
    :param sz: window size in which corner pixels are analysed
    :param is_vrt:
    """
    gdal_ds      = gdal.Open(path_im) if not is_vrt else gdal.OpenShared(path_im)
    get_mean_std = lambda corner_subset: {'mean':np.mean(corner_subset), 'std':np.std(corner_subset)}
    rows,cols    = gdal_ds.RasterYSize,gdal_ds.RasterXSize
    UL = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(0,0,sz,sz))
    UR = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(cols-sz,0,sz,sz))
    LR = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(cols-sz,rows-sz,sz,sz))
    LL = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(0,rows-sz,sz,sz))
    gdal_ds  = None
    possVals = [i['mean'] for i in [UL,UR,LR,LL] if i['std']==0]
    # possVals==[]: all corners are filled with data; np.std(possVals)==0: noDataVal clearly identified
    return None if possVals==[] else possVals[0] if np.std(possVals)==0 else 'unclear'