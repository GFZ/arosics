# -*- coding: utf-8 -*-

# AROSICS - Automated and Robust Open-Source Image Co-Registration Software
#
# Copyright (C) 2017-2020  Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
#
# This software was developed within the context of the GeoMultiSens project funded
# by the German Federal Ministry of Education and Research
# (project grant code: 01 IS 14 010 A-C).
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

import warnings
import sys

# custom
import numpy as np
from geopandas import GeoDataFrame

# internal modules
from py_tools_ds.geo.coord_calc import calc_FullDataset_corner_positions
from py_tools_ds.geo.coord_trafo import pixelToMapYX, imYX2mapYX
from geoarray import GeoArray

__author__ = 'Daniel Scheffler'


def angle_to_north(XY):
    """Calculate the angle in degrees of a given line to north in clockwise direction.

    Angle definition: between [origin:[0,0],north:[0,1]] and [origin:[0,0],pointXY:[X,Y]].
    """
    XY = np.array(XY)
    XYarr = XY if len(XY.shape) == 2 else XY.reshape((1, 2))
    return np.abs(np.degrees(np.array(np.arctan2(XYarr[:, 1], XYarr[:, 0]) - np.pi / 2)) % 360)


def get_true_corner_mapXY(fPath_or_geoarray, band=0, noDataVal=None, mp=1, v=0, q=0):  # pragma: no cover
    """Return the actual map corner coordinates of a given image file or GeoArray instance.

    :param fPath_or_geoarray:
    :param band:                <int> index of the band to be used (starting with 0)
    :param noDataVal:
    :param mp:
    :param v:
    :param q:
    :return:
    """
    # FIXME this function is not used anymore
    warnings.warn('This function is not in use anymore. Use it on your own risk!', DeprecationWarning)
    geoArr = GeoArray(fPath_or_geoarray) if not isinstance(fPath_or_geoarray, GeoArray) else fPath_or_geoarray

    rows, cols = geoArr.shape[:2]
    gt, prj = geoArr.geotransform, geoArr.projection

    assert gt and prj, 'GeoTransform an projection must be given for calculation of LonLat corner coordinates.'

    mask_1bit = np.zeros((rows, cols), dtype='uint8')  # zeros -> image area later overwritten by ones

    if noDataVal is None:
        mask_1bit[:, :] = 1
    elif noDataVal == 'ambiguous':
        warnings.warn("No data value could not be automatically detected. Thus the matching window used for shift "
                      "calculation had to be centered in the middle of the overlap area without respecting no data "
                      "values. To avoid this provide the correct no data values for reference and shift image via "
                      "'-nodata'")
        mask_1bit[:, :] = 1
    else:
        band_data = geoArr[band]  # TODO implement gdal_ReadAsArray_mp (reading in multiprocessing)
        mask_1bit[band_data != noDataVal] = 1

    if v:
        print('detected no data value', noDataVal)

    try:
        corner_coords_YX = calc_FullDataset_corner_positions(mask_1bit, assert_four_corners=False, algorithm='shapely')
    except Exception:
        if v:
            warnings.warn("\nCalculation of corner coordinates failed within algorithm 'shapely' (Exception: %s)."
                          " Using algorithm 'numpy' instead." % sys.exc_info()[1])
        # FIXME numpy algorithm returns wrong values for S2A_OPER_MSI_L1C_TL_SGS__20160608T153121_A005024_T33UUU_B03.jp2
        # FIXME (Hannes)
        corner_coords_YX = \
            calc_FullDataset_corner_positions(mask_1bit, assert_four_corners=False, algorithm='numpy')

    if len(corner_coords_YX) == 4:  # this avoids shapely self intersection
        corner_coords_YX = list(np.array(corner_coords_YX)[[0, 1, 3, 2]])  # UL, UR, LL, LR => UL, UR, LR, LL

    # check if enough unique coordinates have been found
    if not len(GeoDataFrame(corner_coords_YX).drop_duplicates().values) >= 3:
        if not q:
            warnings.warn('\nThe algorithm for automatically detecting the actual image coordinates did not find '
                          'enough unique corners. Using outer image corner coordinates instead.')
        corner_coords_YX = ((0, 0), (0, cols - 1), (rows - 1, 0), (rows - 1, cols - 1))

    # check if all points are unique
    # all_coords_are_unique = len([UL, UR, LL, LR]) == len(GeoDataFrame([UL, UR, LL, LR]).drop_duplicates().values)
    # UL, UR, LL, LR = \
    #   (UL, UR, LL, LR) if all_coords_are_unique else ((0, 0), (0, cols-1), (rows-1, 0), (rows-1, cols-1))

    def get_mapYX(YX): return pixelToMapYX(list(reversed(YX)), geotransform=gt, projection=prj)[0]
    corner_pos_XY = [list(reversed(i)) for i in [get_mapYX(YX) for YX in corner_coords_YX]]
    return corner_pos_XY


def get_subset_GeoTransform(gt_fullArr, subset_box_imYX):
    gt_subset = list(gt_fullArr[:])  # copy
    gt_subset[3], gt_subset[0] = imYX2mapYX(subset_box_imYX[0], gt_fullArr)
    return gt_subset


def get_gdalReadInputs_from_boxImYX(boxImYX):
    """Return row_start,col_start,rows_count,cols_count and assumes boxImYX as [UL_YX,UR_YX,LR_YX,LL_YX)."""
    rS, cS = boxImYX[0]
    clip_sz_x = abs(boxImYX[1][1] - boxImYX[0][1])  # URx-ULx
    clip_sz_y = abs(boxImYX[0][0] - boxImYX[3][0])  # ULy-LLy
    return cS, rS, clip_sz_x, clip_sz_y


def get_GeoArrayPosition_from_boxImYX(boxImYX):
    """Return row_start,row_end,col_start,col_end and assumes boxImYX as [UL_YX,UR_YX,LR_YX,LL_YX)."""
    rS, cS = boxImYX[0]  # UL
    rE, cE = boxImYX[2]  # LR
    return rS, rE - 1, cS, cE - 1  # -1 because boxImYX represents outer box and includes the LR corner of LR pixel
