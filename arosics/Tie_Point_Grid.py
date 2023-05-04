# -*- coding: utf-8 -*-

# AROSICS - Automated and Robust Open-Source Image Co-Registration Software
#
# Copyright (C) 2017-2023
# - Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
# - Helmholtz Centre Potsdam - GFZ German Research Centre for Geosciences Potsdam,
#   Germany (https://www.gfz-potsdam.de/)
#
# This software was developed within the context of the GeoMultiSens project funded
# by the German Federal Ministry of Education and Research
# (project grant code: 01 IS 14 010 A-C).
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import collections
import multiprocessing
import os
import warnings
from time import time, sleep
from typing import Optional, Union

# custom
from osgeo import gdal
import numpy as np
from geopandas import GeoDataFrame
from pandas import DataFrame, Series
from shapely.geometry import Point
from matplotlib import pyplot as plt
from scipy.interpolate import RBFInterpolator, RegularGridInterpolator

# internal modules
from .CoReg import COREG
from py_tools_ds.geo.projection import isLocal
from py_tools_ds.io.pathgen import get_generic_outpath
from py_tools_ds.processing.progress_mon import ProgressBar
from py_tools_ds.geo.vector.conversion import points_to_raster
from py_tools_ds.io.vector.writer import write_shp
from geoarray import GeoArray

from .CoReg import GeoArray_CoReg  # noqa F401  # flake8 issue

__author__ = 'Daniel Scheffler'

global_shared_imref: Optional[Union[GeoArray, str]] = None
global_shared_im2shift: Optional[Union[GeoArray, str]] = None


def mp_initializer(imref, imtgt):
    """Declare global variables needed for self._get_spatial_shifts().

    :param imref:   reference image
    :param imtgt:   target image
    """
    global global_shared_imref, global_shared_im2shift
    global_shared_imref = imref
    global_shared_im2shift = imtgt


class Tie_Point_Grid(object):
    """
    The 'Tie_Point_Grid' class applies the algorithm to detect spatial shifts to the overlap area of the input images.

    Spatial shifts are calculated for each point in grid of which the parameters can be adjusted using keyword
    arguments. Shift correction performs a polynomial transformation using te calculated shifts of each point in the
    grid as GCPs. Thus, 'Tie_Point_Grid' can be used to correct for locally varying geometric distortions of the target
    image.

    See help(Tie_Point_Grid) for documentation!
    """

    def __init__(self,
                 COREG_obj: COREG,
                 grid_res: float,
                 max_points: int = None,
                 outFillVal: int = -9999,
                 resamp_alg_calc: str = 'cubic',
                 tieP_filter_level: int = 3,
                 outlDetect_settings: dict = None,
                 dir_out: str = None,
                 CPUs: int = None,
                 progress: bool = True,
                 v: bool = False,
                 q: bool = False):
        """Get an instance of the 'Tie_Point_Grid' class.

        :param COREG_obj:
            an instance of COREG class

        :param grid_res:
            grid resolution in pixels of the target image (x-direction)

        :param max_points:
            maximum number of points used to find coregistration tie points

            NOTE: Points are selected randomly from the given point grid (specified by 'grid_res'). If the point does
            not provide enough points, all available points are chosen.

        :param outFillVal:
            if given the generated tie points grid is filled with this value in case no match could be found during
            co-registration (default: -9999)

        :param resamp_alg_calc:
            the resampling algorithm to be used for all warping processes during calculation of spatial shifts
            (valid algorithms: nearest, bilinear, cubic, cubic_spline, lanczos, average, mode, max, min, med, q1, q3)
            default: cubic (highly recommended)

        :param tieP_filter_level:
            filter tie points used for shift correction in different levels (default: 3).
            NOTE: lower levels are also included if a higher level is chosen

            - Level 0: no tie point filtering
            - Level 1: Reliablity filtering
                       - filter all tie points out that have a low reliability according to internal tests
            - Level 2: SSIM filtering
                       - filters all tie points out where shift correction does not increase image similarity within
                         matching window (measured by mean structural similarity index)
            - Level 3: RANSAC outlier detection

        :param outlDetect_settings:
            a dictionary with the settings to be passed to arosics.TiePointGrid.Tie_Point_Refiner.
            Available keys: min_reliability, rs_max_outlier, rs_tolerance, rs_max_iter, rs_exclude_previous_outliers,
            rs_timeout, rs_random_state, q. See documentation there.

        :param dir_out:
            output directory to be used for all outputs if nothing else is given to the individual methods

        :param CPUs:
            number of CPUs to use during calculation of tie points grid
            (default: None, which means 'all CPUs available')

        :param progress:
            show progress bars (default: True)

        :param v:
            verbose mode (default: False)

        :param q:
            quiet mode (default: False)
        """
        if not isinstance(COREG_obj, COREG):
            raise ValueError("'COREG_obj' must be an instance of COREG class.")

        self.COREG_obj = COREG_obj  # type: COREG
        self.grid_res = grid_res
        self.max_points = max_points
        self.outFillVal = outFillVal
        self.rspAlg_calc = resamp_alg_calc
        self.tieP_filter_level = tieP_filter_level
        self.outlDetect_settings = outlDetect_settings or dict()
        self.dir_out = dir_out
        self.CPUs = CPUs
        self.v = v
        self.q = q if not v else False  # overridden by v
        self.progress = progress if not q else False  # overridden by q

        if 'q' not in self.outlDetect_settings:
            self.outlDetect_settings['q'] = self.q

        self.ref = self.COREG_obj.ref  # type: GeoArray_CoReg
        self.shift = self.COREG_obj.shift  # type: GeoArray_CoReg

        self.XY_points, self.XY_mapPoints = self._get_imXY__mapXY_points(self.grid_res)
        self._CoRegPoints_table = None  # set by self.CoRegPoints_table
        self._GCPList = None  # set by self.to_GCPList()

    @property
    def mean_x_shift_px(self):
        return self.CoRegPoints_table['X_SHIFT_PX'][self.CoRegPoints_table['X_SHIFT_PX'] != self.outFillVal].mean()

    @property
    def mean_y_shift_px(self):
        return self.CoRegPoints_table['Y_SHIFT_PX'][self.CoRegPoints_table['Y_SHIFT_PX'] != self.outFillVal].mean()

    @property
    def mean_x_shift_map(self):
        return self.CoRegPoints_table['X_SHIFT_M'][self.CoRegPoints_table['X_SHIFT_M'] != self.outFillVal].mean()

    @property
    def mean_y_shift_map(self):
        return self.CoRegPoints_table['Y_SHIFT_M'][self.CoRegPoints_table['Y_SHIFT_M'] != self.outFillVal].mean()

    @property
    def CoRegPoints_table(self):
        """Return a GeoDataFrame containing all the results from coregistration for all points in the tie point grid.

        Columns of the GeoDataFrame: 'geometry','POINT_ID','X_IM','Y_IM','X_MAP','Y_MAP','X_WIN_SIZE', 'Y_WIN_SIZE',
                                     'X_SHIFT_PX','Y_SHIFT_PX', 'X_SHIFT_M', 'Y_SHIFT_M', 'ABS_SHIFT' and 'ANGLE'
        """
        if self._CoRegPoints_table is not None:
            return self._CoRegPoints_table
        else:
            self._CoRegPoints_table = self.get_CoRegPoints_table()
            return self._CoRegPoints_table

    @CoRegPoints_table.setter
    def CoRegPoints_table(self, CoRegPoints_table):
        self._CoRegPoints_table = CoRegPoints_table

    @property
    def GCPList(self):
        """Return a list of GDAL compatible GCP objects."""
        if self._GCPList:
            return self._GCPList
        else:
            self._GCPList = self.to_GCPList()
            return self._GCPList

    @GCPList.setter
    def GCPList(self, GCPList):
        self._GCPList = GCPList

    def _get_imXY__mapXY_points(self, grid_res):
        """Return a numpy array containing possible positions for coregistration tie points.

        NOTE: The returned positions are dependent from the given grid resolution.

        :param grid_res:
        :return:
        """
        if not self.q:
            print('Initializing tie points grid...')

        Xarr, Yarr = np.meshgrid(np.arange(0, self.shift.shape[1] + grid_res, grid_res),
                                 np.arange(0, self.shift.shape[0] + grid_res, grid_res))

        mapXarr = np.full_like(Xarr, self.shift.gt[0], dtype=np.float64) + Xarr * self.shift.gt[1]
        mapYarr = np.full_like(Yarr, self.shift.gt[3], dtype=np.float64) - Yarr * abs(self.shift.gt[5])

        XY_points = np.empty((Xarr.size, 2), Xarr.dtype)
        XY_points[:, 0] = Xarr.flat
        XY_points[:, 1] = Yarr.flat

        XY_mapPoints = np.empty((mapXarr.size, 2), mapXarr.dtype)
        XY_mapPoints[:, 0] = mapXarr.flat
        XY_mapPoints[:, 1] = mapYarr.flat

        assert XY_points.shape == XY_mapPoints.shape

        return XY_points, XY_mapPoints

    def _exclude_bad_XYpos(self, GDF):
        """Exclude all points outside of the image overlap area and where the bad data mask is True (if given).

        :param GDF:     <geopandas.GeoDataFrame> must include the columns 'X_MAP' and 'Y_MAP'
        :return:
        """
        from skimage.measure import points_in_poly  # import here to avoid static TLS ImportError

        # exclude all points outside of overlap area
        inliers = points_in_poly(self.XY_mapPoints,
                                 np.swapaxes(np.array(self.COREG_obj.overlap_poly.exterior.coords.xy), 0, 1))
        GDF = GDF[inliers].copy()
        # GDF = GDF[GDF['geometry'].within(self.COREG_obj.overlap_poly.simplify(tolerance=15))] # works but much slower

        assert not GDF.empty, 'No coregistration point could be placed within the overlap area. Check your input data!'

        # exclude all points where bad data mask is True (e.g. points on clouds etc.)
        orig_len_GDF = len(GDF)  # length of GDF after dropping all points outside the overlap polygon
        mapXY = np.array(GDF.loc[:, ['X_MAP', 'Y_MAP']])
        GDF['REF_BADDATA'] = self.COREG_obj.ref.mask_baddata.read_pointData(mapXY) \
            if self.COREG_obj.ref.mask_baddata is not None else False
        GDF['TGT_BADDATA'] = self.COREG_obj.shift.mask_baddata.read_pointData(mapXY) \
            if self.COREG_obj.shift.mask_baddata is not None else False
        GDF = GDF[(~GDF['REF_BADDATA']) & (~GDF['TGT_BADDATA'])]
        if self.COREG_obj.ref.mask_baddata is not None or self.COREG_obj.shift.mask_baddata is not None:
            if not self.q:
                if not GDF.empty:
                    print('With respect to the provided bad data mask(s) %s points of initially %s have been excluded.'
                          % (orig_len_GDF - len(GDF), orig_len_GDF))
                else:
                    warnings.warn('With respect to the provided bad data mask(s) no coregistration point could be '
                                  'placed within an image area usable for coregistration.')

        return GDF

    @staticmethod
    def _get_spatial_shifts(coreg_kwargs):
        # unpack
        pointID = coreg_kwargs['pointID']
        fftw_works = coreg_kwargs['fftw_works']
        del coreg_kwargs['pointID'], coreg_kwargs['fftw_works']

        # assertions
        assert global_shared_imref is not None
        assert global_shared_im2shift is not None

        # run CoReg
        CR = COREG(global_shared_imref, global_shared_im2shift, CPUs=1, **coreg_kwargs)
        CR.fftw_works = fftw_works
        CR.calculate_spatial_shifts()

        # fetch results
        last_err = CR.tracked_errors[-1] if CR.tracked_errors else None
        win_sz_y, win_sz_x = CR.matchBox.imDimsYX if CR.matchBox else (None, None)
        CR_res = [win_sz_x, win_sz_y, CR.x_shift_px, CR.y_shift_px, CR.x_shift_map, CR.y_shift_map,
                  CR.vec_length_map, CR.vec_angle_deg, CR.ssim_orig, CR.ssim_deshifted, CR.ssim_improved,
                  CR.shift_reliability, last_err]

        return [pointID] + CR_res

    def _get_coreg_kwargs(self, pID, wp):
        return dict(
            pointID=pID,
            fftw_works=self.COREG_obj.fftw_works,
            wp=wp,
            ws=self.COREG_obj.win_size_XY,
            resamp_alg_calc=self.rspAlg_calc,
            footprint_poly_ref=self.COREG_obj.ref.poly,
            footprint_poly_tgt=self.COREG_obj.shift.poly,
            r_b4match=self.ref.band4match + 1,  # band4match is internally saved as index, starting from 0
            s_b4match=self.shift.band4match + 1,  # band4match is internally saved as index, starting from 0
            max_iter=self.COREG_obj.max_iter,
            max_shift=self.COREG_obj.max_shift,
            nodata=(self.COREG_obj.ref.nodata, self.COREG_obj.shift.nodata),
            force_quadratic_win=self.COREG_obj.force_quadratic_win,
            binary_ws=self.COREG_obj.bin_ws,
            v=False,  # otherwise this would lead to massive console output
            q=True,  # otherwise this would lead to massive console output
            ignore_errors=True
        )

    def get_CoRegPoints_table(self):
        assert self.XY_points is not None and self.XY_mapPoints is not None

        # create a dataframe containing 'geometry','POINT_ID','X_IM','Y_IM','X_MAP','Y_MAP'
        # (convert imCoords to mapCoords
        XYarr2PointGeom = np.vectorize(lambda X, Y: Point(X, Y), otypes=[Point])
        geomPoints = np.array(XYarr2PointGeom(self.XY_mapPoints[:, 0], self.XY_mapPoints[:, 1]))

        crs = self.COREG_obj.shift.prj if not isLocal(self.COREG_obj.shift.prj) else None

        GDF = GeoDataFrame(index=range(len(geomPoints)),
                           crs=crs,
                           columns=['geometry', 'POINT_ID', 'X_IM', 'Y_IM', 'X_MAP', 'Y_MAP'])
        GDF['geometry'] = geomPoints
        GDF['POINT_ID'] = range(len(geomPoints))
        GDF[['X_IM', 'Y_IM']] = self.XY_points
        GDF[['X_MAP', 'Y_MAP']] = self.XY_mapPoints

        # exclude offsite points and points on bad data mask
        GDF = self._exclude_bad_XYpos(GDF)
        if GDF.empty:
            self.CoRegPoints_table = GDF
            return self.CoRegPoints_table

        # choose a random subset of points if a maximum number has been given
        if self.max_points and len(GDF) > self.max_points:
            GDF = GDF.sample(self.max_points).copy()

        # equalize pixel grids in order to save warping time
        if len(GDF) > 100:
            # NOTE: actually grid res should be also changed here because self.shift.xgsd changes and grid res is
            # connected to that
            self.COREG_obj.equalize_pixGrids()
            self.ref = self.COREG_obj.ref
            self.shift = self.COREG_obj.shift

        # validate reference and target image inputs
        assert self.ref.footprint_poly  # this also checks for mask_nodata and nodata value
        assert self.shift.footprint_poly

        # ensure the input arrays for CoReg are in memory -> otherwise the code will get stuck in multiprocessing if
        # neighboured matching windows overlap during reading from disk!!
        self.ref.cache_array_subset(
            [self.COREG_obj.ref.band4match])  # only sets geoArr._arr_cache; does not change number of bands
        self.shift.cache_array_subset([self.COREG_obj.shift.band4match])

        # get all variations of kwargs for coregistration
        list_coreg_kwargs = (self._get_coreg_kwargs(i, self.XY_mapPoints[i]) for i in GDF.index)  # generator

        # run co-registration for whole grid
        if self.CPUs is None or self.CPUs > 1:
            if not self.q:
                cpus = self.CPUs if self.CPUs is not None else multiprocessing.cpu_count()
                print("Calculating tie point grid (%s points) using %s CPU cores..." % (len(GDF), cpus))

            with multiprocessing.Pool(self.CPUs, initializer=mp_initializer, initargs=(self.ref, self.shift)) as pool:
                if self.q or not self.progress:
                    results = pool.map(self._get_spatial_shifts, list_coreg_kwargs)
                else:
                    results = pool.map_async(self._get_spatial_shifts, list_coreg_kwargs, chunksize=1)
                    bar = ProgressBar(prefix='\tprogress:')
                    while True:
                        sleep(.1)
                        # this does not really represent the remaining tasks but the remaining chunks
                        # -> thus chunksize=1
                        # noinspection PyProtectedMember
                        numberDone = len(GDF) - results._number_left
                        if self.progress:
                            bar.print_progress(percent=numberDone / len(GDF) * 100)
                        if results.ready():
                            # <= this is the line where multiprocessing can freeze if an exception appears within
                            # COREG and is not raised
                            results = results.get()
                            break
                pool.close()  # needed to make coverage work in multiprocessing
                pool.join()

        else:
            # declare global variables needed for self._get_spatial_shifts()
            global global_shared_imref, global_shared_im2shift
            global_shared_imref = self.ref
            global_shared_im2shift = self.shift

            if not self.q:
                print("Calculating tie point grid (%s points) on 1 CPU core..." % len(GDF))
            results = np.empty((len(geomPoints), 14), object)
            bar = ProgressBar(prefix='\tprogress:')
            for i, coreg_kwargs in enumerate(list_coreg_kwargs):
                if self.progress:
                    bar.print_progress((i + 1) / len(GDF) * 100)
                results[i, :] = self._get_spatial_shifts(coreg_kwargs)

        # merge results with GDF
        # NOTE: We use a pandas.DataFrame here because the geometry column is missing.
        #       GDF.astype(...) fails with geopandas>0.6.0 if the geometry columns is missing.
        records = DataFrame(results,
                            columns=['POINT_ID', 'X_WIN_SIZE', 'Y_WIN_SIZE', 'X_SHIFT_PX', 'Y_SHIFT_PX', 'X_SHIFT_M',
                                     'Y_SHIFT_M', 'ABS_SHIFT', 'ANGLE', 'SSIM_BEFORE', 'SSIM_AFTER',
                                     'SSIM_IMPROVED', 'RELIABILITY', 'LAST_ERR'])

        # merge DataFrames (dtype must be equal to records.dtypes; We need np.object due to None values)
        GDF = GDF.astype(object).merge(records.astype(object), on='POINT_ID', how="inner")
        GDF = GDF.replace([np.nan, None], int(self.outFillVal))  # fillna fails with geopandas==0.6.0
        GDF.crs = crs  # gets lost when using GDF.astype(np.object), so we have to reassign that

        if not self.q:
            print("Found %s matches." % len(GDF[GDF.LAST_ERR == int(self.outFillVal)]))

        # filter tie points according to given filter level
        if self.tieP_filter_level > 0:
            if not self.q:
                print('Performing validity checks...')
            TPR = Tie_Point_Refiner(GDF[GDF.ABS_SHIFT != self.outFillVal], **self.outlDetect_settings)
            GDF_filt, new_columns = TPR.run_filtering(level=self.tieP_filter_level)
            GDF = GDF.merge(GDF_filt[['POINT_ID'] + new_columns], on='POINT_ID', how="outer")

        GDF = GDF.replace([np.nan, None], int(self.outFillVal))  # fillna fails with geopandas==0.6.0

        self.CoRegPoints_table = GDF

        if not self.q:
            if GDF.empty:
                warnings.warn('No valid GCPs could by identified.')
            else:
                if self.tieP_filter_level > 0:
                    print("%d valid tie points remain after filtering." % len(GDF[GDF.OUTLIER.__eq__(False)]))

        return self.CoRegPoints_table

    def calc_rmse(self, include_outliers: bool = False) -> float:
        """Calculate root mean square error of absolute shifts from the tie point grid.

        :param include_outliers:    whether to include tie points that have been marked as false-positives (if present)
        """
        if self.CoRegPoints_table.empty:
            raise RuntimeError('Cannot compute the RMSE because no tie points were found at all.')

        tbl = self.CoRegPoints_table
        tbl = tbl if include_outliers else tbl[tbl['OUTLIER'] == 0].copy() if 'OUTLIER' in tbl.columns else tbl

        if not include_outliers and tbl.empty:
            raise RuntimeError('Cannot compute the RMSE because all tie points are flagged as false-positives.')

        shifts = np.array(tbl['ABS_SHIFT'])
        shifts_sq = [i * i for i in shifts if i != self.outFillVal]

        return np.sqrt(sum(shifts_sq) / len(shifts_sq))

    def calc_overall_ssim(self,
                          include_outliers: bool = False,
                          after_correction: bool = True
                          ) -> float:
        """Calculate the median value of all SSIM values contained in tie point grid.

        :param include_outliers:    whether to include tie points that have been marked as false-positives
        :param after_correction:    whether to compute median SSIM before correction or after
        """
        if self.CoRegPoints_table.empty:
            raise RuntimeError('Cannot compute the overall SSIM because no tie points were found at all.')

        tbl = self.CoRegPoints_table
        tbl = tbl if include_outliers else tbl[tbl['OUTLIER'] == 0].copy()

        if not include_outliers and tbl.empty:
            raise RuntimeError('Cannot compute the overall SSIM because all tie points are flagged as false-positives.')

        ssim_col = np.array(tbl['SSIM_AFTER' if after_correction else 'SSIM_BEFORE'])
        ssim_col = [i * i for i in ssim_col if i != self.outFillVal]

        return float(np.median(ssim_col))

    def calc_overall_stats(self, include_outliers: bool = False) -> dict:
        """Calculate statistics like RMSE, MSE, MAE, ... from the tie point grid.

        Full list of returned statistics:

        - N_TP:                 number of tie points
        - N_VALID_TP:           number of valid tie points
        - N_INVALID_TP:         number of invalid tie points (false-positives)
        - PERC_VALID_TP:        percentage of valid tie points
        - RMSE_M:               root mean squared error of absolute shift vector length in map units
        - RMSE_X_M:             root mean squared error of shift vector length in x-direction in map units
        - RMSE_Y_M:             root mean squared error of shift vector length in y-direction in map units
        - RMSE_X_PX:            root mean squared error of shift vector length in x-direction in pixel units
        - RMSE_Y_PX:            root mean squared error of shift vector length in y-direction in pixel units
        - MSE_M:                mean squared error of absolute shift vector length in map units
        - MSE_X_M:              mean squared error of shift vector length in x-direction in map units
        - MSE_Y_M:              mean squared error of shift vector length in y-direction in map units
        - MSE_X_PX:             mean squared error of shift vector length in x-direction in pixel units
        - MSE_Y_PX:             mean squared error of shift vector length in y-direction in pixel units
        - MAE_M:                mean absolute error of absolute shift vector length in map units
        - MAE_X_M:              mean absolute error of shift vector length in x-direction in map units
        - MAE_Y_M:              mean absolute error of shift vector length in y-direction in map units
        - MAE_X_PX:             mean absolute error of shift vector length in x-direction in pixel units
        - MAE_Y_PX:             mean absolute error of shift vector length in y-direction in pixel units
        - MEAN_ABS_SHIFT:       mean absolute shift vector length in map units
        - MEAN_X_SHIFT_M:       mean shift vector length in x-direction in map units
        - MEAN_Y_SHIFT_M:       mean shift vector length in y-direction in map units
        - MEAN_X_SHIFT_PX:      mean shift vector length in x-direction in pixel units
        - MEAN_Y_SHIFT_PX:      mean shift vector length in y-direction in pixel units
        - MEAN_ANGLE:           mean direction of the shift vectors in degrees from north
        - MEAN_SSIM_BEFORE:     mean structural similatity index within each matching window before co-registration
        - MEAN_SSIM_AFTER:      mean structural similatity index within each matching window after co-registration
        - MEAN_RELIABILITY:     mean tie point reliability in percent
        - MEDIAN_ABS_SHIFT:     median absolute shift vector length in map units
        - MEDIAN_X_SHIFT_M:     median shift vector length in x-direction in map units
        - MEDIAN_Y_SHIFT_M:     median shift vector length in y-direction in map units
        - MEDIAN_X_SHIFT_PX:    median shift vector length in x-direction in pixel units
        - MEDIAN_Y_SHIFT_PX:    median shift vector length in y-direction in pixel units
        - MEDIAN_ANGLE:         median direction of the shift vectors in degrees from north
        - MEDIAN_SSIM_BEFORE:   median structural similatity index within each matching window before co-registration
        - MEDIAN_SSIM_AFTER:    median structural similatity index within each matching window after co-registration
        - MEDIAN_RELIABILITY:   median tie point reliability in percent
        - STD_ABS_SHIFT:        standard deviation of absolute shift vector length in map units
        - STD_X_SHIFT_M:        standard deviation of shift vector length in x-direction in map units
        - STD_Y_SHIFT_M:        standard deviation of shift vector length in y-direction in map units
        - STD_X_SHIFT_PX:       standard deviation of shift vector length in x-direction in pixel units
        - STD_Y_SHIFT_PX:       standard deviation of shift vector length in y-direction in pixel units
        - STD_ANGLE:            standard deviation of direction of the shift vectors in degrees from north
        - STD_SSIM_BEFORE:      standard deviation of structural similatity index within each matching window before
                                co-registration
        - STD_SSIM_AFTER:       standard deviation of structural similatity index within each matching window after
                                co-registration
        - STD_RELIABILITY:      standard deviation of tie point reliability in percent
        - MIN_ABS_SHIFT:        minimal absolute shift vector length in map units
        - MIN_X_SHIFT_M:        minimal shift vector length in x-direction in map units
        - MIN_Y_SHIFT_M:        minimal shift vector length in y-direction in map units
        - MIN_X_SHIFT_PX:       minimal shift vector length in x-direction in pixel units
        - MIN_Y_SHIFT_PX:       minimal shift vector length in y-direction in pixel units
        - MIN_ANGLE:            minimal direction of the shift vectors in degrees from north
        - MIN_SSIM_BEFORE:      minimal structural similatity index within each matching window before co-registration
        - MIN_SSIM_AFTER:       minimal structural similatity index within each matching window after co-registration
        - MIN_RELIABILITY:      minimal tie point reliability in percent
        - MIN_ABS_SHIFT:        maximal absolute shift vector length in map units
        - MAX_X_SHIFT_M:        maximal shift vector length in x-direction in map units
        - MAX_Y_SHIFT_M:        maximal shift vector length in y-direction in map units
        - MAX_X_SHIFT_PX:       maximal shift vector length in x-direction in pixel units
        - MAX_Y_SHIFT_PX:       maximal shift vector length in y-direction in pixel units
        - MAX_ANGLE:            maximal direction of the shift vectors in degrees from north
        - MAX_SSIM_BEFORE:      maximal structural similatity index within each matching window before co-registration
        - MAX_SSIM_AFTER:       maximal structural similatity index within each matching window after co-registration
        - MAX_RELIABILITY:      maximal tie point reliability in percent

        :param include_outliers:    whether to include tie points that have been marked as false-positives (if present)
        """
        if self.CoRegPoints_table.empty:
            raise RuntimeError('Cannot compute overall statistics because no tie points were found at all.')

        tbl = self.CoRegPoints_table

        n_tiepoints = sum(tbl['ABS_SHIFT'] != self.outFillVal)
        n_outliers = sum(tbl['OUTLIER'] == 1)

        tbl = tbl if include_outliers else tbl[tbl['OUTLIER'] == 0].copy() if 'OUTLIER' in tbl.columns else tbl
        tbl = tbl.copy().replace(self.outFillVal, np.nan)

        if not include_outliers and tbl.empty:
            raise RuntimeError('Cannot compute overall statistics '
                               'because all tie points are flagged as false-positives.')

        def RMSE(shifts):
            shifts_sq = shifts ** 2
            return np.sqrt(sum(shifts_sq) / len(shifts_sq))

        def MSE(shifts):
            shifts_sq = shifts ** 2
            return sum(shifts_sq) / len(shifts_sq)

        def MAE(shifts):
            shifts_abs = np.abs(shifts)
            return sum(shifts_abs) / len(shifts_abs)

        abs_shift, x_shift_m, y_shift_m, x_shift_px, y_shift_px, angle, ssim_before, ssim_after, reliability = \
            [tbl[k].dropna().values for k in ['ABS_SHIFT', 'X_SHIFT_M', 'Y_SHIFT_M', 'X_SHIFT_PX', 'Y_SHIFT_PX',
                                              'ANGLE', 'SSIM_BEFORE', 'SSIM_AFTER', 'RELIABILITY']]

        stats = dict(
            N_TP=n_tiepoints,
            N_VALID_TP=len(abs_shift),
            N_INVALID_TP=n_outliers,
            PERC_VALID_TP=(n_tiepoints - n_outliers) / n_tiepoints * 100,

            RMSE_M=RMSE(abs_shift),
            RMSE_X_M=RMSE(x_shift_m),
            RMSE_Y_M=RMSE(y_shift_m),
            RMSE_X_PX=RMSE(x_shift_px),
            RMSE_Y_PX=RMSE(y_shift_px),

            MSE_M=MSE(abs_shift),
            MSE_X_M=MSE(x_shift_m),
            MSE_Y_M=MSE(y_shift_m),
            MSE_X_PX=MSE(x_shift_px),
            MSE_Y_PX=MSE(y_shift_px),

            MAE_M=MAE(abs_shift),
            MAE_X_M=MAE(x_shift_m),
            MAE_Y_M=MAE(y_shift_m),
            MAE_X_PX=MAE(x_shift_px),
            MAE_Y_PX=MAE(y_shift_px),
        )

        for stat, func in zip(['mean', 'median', 'std', 'min', 'max'],
                              [np.mean, np.median, np.std, np.min, np.max]):
            for n in ['abs_shift', 'x_shift_m', 'y_shift_m', 'x_shift_px', 'y_shift_px',
                      'angle', 'ssim_before', 'ssim_after', 'reliability']:

                vals = locals()[n]
                stats[f'{stat}_{n}'.upper()] = func(vals)

        return stats

    def plot_shift_distribution(self,
                                include_outliers: bool = True,
                                unit: str = 'm',
                                interactive: bool = False,
                                figsize: tuple = None,
                                xlim: list = None,
                                ylim: list = None,
                                fontsize: int = 12,
                                title: str = 'shift distribution',
                                savefigPath: str = '',
                                savefigDPI: int = 96,
                                showFig: bool = True,
                                return_fig: bool = False
                                ) -> tuple:
        """Create a 2D scatterplot containing the distribution of calculated X/Y-shifts.

        :param include_outliers:    whether to include tie points that have been marked as false-positives
        :param unit:                'm' for meters or 'px' for pixels (default: 'm')
        :param interactive:         whether to use interactive mode (uses plotly for visualization)
        :param figsize:             (xdim, ydim)
        :param xlim:                [xmin, xmax]
        :param ylim:                [ymin, ymax]
        :param fontsize:            size of all used fonts
        :param title:               the title to be plotted above the figure
        :param savefigPath:         path where to save the figure
        :param savefigDPI:          DPI resolution of the output figure when saved to disk
        :param showFig:             whether to show or to hide the figure
        :param return_fig:          whether to return the figure and axis objects
        """
        from matplotlib import pyplot as plt

        if unit not in ['m', 'px']:
            raise ValueError("Parameter 'unit' must have the value 'm' (meters) or 'px' (pixels)! Got %s." % unit)

        if self.CoRegPoints_table.empty:
            raise RuntimeError('Shift distribution cannot be plotted because no tie points were found at all.')

        tbl = self.CoRegPoints_table
        tbl = tbl[tbl['ABS_SHIFT'] != self.outFillVal]
        tbl_il = tbl[tbl['OUTLIER'] == 0].copy() if 'OUTLIER' in tbl.columns else tbl
        tbl_ol = tbl[tbl['OUTLIER']].copy() if 'OUTLIER' in tbl.columns else None
        x_attr = 'X_SHIFT_M' if unit == 'm' else 'X_SHIFT_PX'
        y_attr = 'Y_SHIFT_M' if unit == 'm' else 'Y_SHIFT_PX'
        rmse = self.calc_rmse(include_outliers=False)  # always exclude outliers when calculating RMSE
        figsize = figsize if figsize else (10, 10)

        if interactive:
            from plotly.offline import iplot, init_notebook_mode
            import plotly.graph_objs as go
            # FIXME outliers are not plotted

            init_notebook_mode(connected=True)

            # Create a trace
            trace = go.Scatter(
                x=tbl_il[x_attr],
                y=tbl_il[y_attr],
                mode='markers'
            )

            data = [trace]

            # Plot and embed in ipython notebook!
            iplot(data, filename='basic-scatter')

            return None, None

        else:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)

            if include_outliers and 'OUTLIER' in tbl.columns:
                ax.scatter(tbl_ol[x_attr], tbl_ol[y_attr], marker='+', c='r', label='false-positives')
            ax.scatter(tbl_il[x_attr], tbl_il[y_attr], marker='+', c='g', label='valid tie points')

            # set axis limits
            if not xlim:
                xmax = np.abs(tbl_il[x_attr]).max()
                xlim = [-xmax, xmax]
            if not ylim:
                ymax = np.abs(tbl_il[y_attr]).max()
                ylim = [-ymax, ymax]
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            # add text box containing RMSE of plotted shifts
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
            ax.text(xlim[1] - (xlim[1] / 20), -ylim[1] + (ylim[1] / 20),
                    'RMSE:  %s m / %s px' % (np.round(rmse, 2), np.round(rmse / self.shift.xgsd, 2)),
                    ha='right', va='bottom', fontsize=fontsize, bbox=dict(facecolor='w', pad=None, alpha=0.8))

            # add grid and increase linewidth of middle line
            ax.grid(visible=True)
            ax.spines["right"].set_visible(True)
            ax.spines["top"].set_visible(True)
            xgl = ax.get_xgridlines()
            middle_xgl = xgl[int(np.median(np.array(range(len(xgl)))))]
            middle_xgl.set_linewidth(2)
            middle_xgl.set_linestyle('-')
            ygl = ax.get_ygridlines()
            middle_ygl = ygl[int(np.median(np.array(range(len(ygl)))))]
            middle_ygl.set_linewidth(2)
            middle_ygl.set_linestyle('-')

            # set title and adjust tick labels
            ax.set_title(title, fontsize=fontsize)
            [tick.label1.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
            [tick.label1.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
            ax.set_xlabel('x-shift [%s]' % 'meters' if unit == 'm' else 'pixels', fontsize=fontsize)
            ax.set_ylabel('y-shift [%s]' % 'meters' if unit == 'm' else 'pixels', fontsize=fontsize)

            # add legend with labels in the right order
            handles, labels = ax.get_legend_handles_labels()
            leg = plt.legend(reversed(handles), reversed(labels), fontsize=fontsize, loc='upper right', scatterpoints=3)
            leg.get_frame().set_edgecolor('black')

            # remove white space around the figure
            fig.subplots_adjust(top=.94, bottom=.06, right=.96, left=.09)

            if savefigPath:
                fig.savefig(savefigPath, dpi=savefigDPI, pad_inches=0.3, bbox_inches='tight')

            if return_fig:
                return fig, ax

            if showFig and not self.q:
                plt.show(block=True)
            else:
                plt.close(fig)

    def dump_CoRegPoints_table(self, path_out=None):
        if self.CoRegPoints_table.empty:
            raise RuntimeError('Cannot dump tie points table because it is empty.')

        path_out = path_out if path_out else \
            get_generic_outpath(dir_out=self.dir_out,
                                fName_out="CoRegPoints_table_grid%s_ws(%s_%s)__T_%s__R_%s.pkl"
                                          % (self.grid_res, self.COREG_obj.win_size_XY[0],
                                             self.COREG_obj.win_size_XY[1], self.shift.basename,
                                             self.ref.basename))
        if not self.q:
            print('Writing %s ...' % path_out)
        self.CoRegPoints_table.to_pickle(path_out)

    def to_GCPList(self):
        # get copy of tie points grid without no data
        try:
            GDF = self.CoRegPoints_table.loc[self.CoRegPoints_table.ABS_SHIFT != self.outFillVal, :].copy()
        except AttributeError:
            # self.CoRegPoints_table has no attribute 'ABS_SHIFT' because all points have been excluded
            return []

        if getattr(GDF, 'empty'):  # GDF.empty returns AttributeError
            return []
        else:
            # exclude all points flagged as outliers
            if 'OUTLIER' in GDF.columns:
                GDF = GDF[GDF.OUTLIER.__eq__(False)].copy()
            avail_TP = len(GDF)

            if not avail_TP:
                # no point passed all validity checks
                return []

            if avail_TP > 7000:
                GDF = GDF.sample(7000)
                warnings.warn('By far not more than 7000 tie points can be used for warping within a limited '
                              'computation time (due to a GDAL bottleneck). Thus these 7000 points are randomly chosen '
                              'out of the %s available tie points.' % avail_TP)

            # calculate GCPs
            GDF['X_MAP_new'] = GDF.X_MAP + GDF.X_SHIFT_M
            GDF['Y_MAP_new'] = GDF.Y_MAP + GDF.Y_SHIFT_M
            GDF['GCP'] = GDF.apply(lambda GDF_row: gdal.GCP(GDF_row.X_MAP_new,
                                                            GDF_row.Y_MAP_new,
                                                            0,
                                                            GDF_row.X_IM,
                                                            GDF_row.Y_IM),
                                   axis=1)
            self.GCPList = GDF.GCP.tolist()

            return self.GCPList

    def test_if_singleprocessing_equals_multiprocessing_result(self):
        # RANSAC filtering always produces different results because it includes random sampling
        self.tieP_filter_level = 1

        self.CPUs = None
        dataframe = self.get_CoRegPoints_table()
        mp_out = np.empty_like(dataframe.values)
        mp_out[:] = dataframe.values
        self.CPUs = 1
        dataframe = self.get_CoRegPoints_table()
        sp_out = np.empty_like(dataframe.values)
        sp_out[:] = dataframe.values

        return np.array_equal(sp_out, mp_out)

    def _get_line_by_PID(self, PID):
        return self.CoRegPoints_table.loc[PID, :]

    def _get_lines_by_PIDs(self, PIDs):
        assert isinstance(PIDs, list)
        lines = np.zeros((len(PIDs), self.CoRegPoints_table.shape[1]))
        for i, PID in enumerate(PIDs):
            lines[i, :] = self.CoRegPoints_table[self.CoRegPoints_table['POINT_ID'] == PID]
        return lines

    def to_PointShapefile(self,
                          path_out: str = None,
                          skip_nodata: bool = True,
                          skip_nodata_col: str = 'ABS_SHIFT'
                          ) -> None:
        """Write the calculated tie point grid to a point shapefile (e.g., for visualization by a GIS software).

        NOTE: The shapefile uses Tie_Point_Grid.CoRegPoints_table as attribute table.

        :param path_out:        <str> the output path. If not given, it is automatically defined.
        :param skip_nodata:     <bool> whether to skip all points where no valid match could be found
        :param skip_nodata_col: <str> determines which column of Tie_Point_Grid.CoRegPoints_table is used to
                                identify points where no valid match could be found
        """
        if self.CoRegPoints_table.empty:
            raise RuntimeError('Cannot save a point shapefile because no tie points were found at all.')

        GDF = self.CoRegPoints_table

        if skip_nodata:
            GDF2pass = GDF[GDF[skip_nodata_col] != self.outFillVal].copy()
        else:
            GDF2pass = GDF
            GDF2pass.LAST_ERR = GDF2pass.apply(lambda GDF_row: repr(GDF_row.LAST_ERR), axis=1)

        # replace boolean values (cannot be written)
        GDF2pass = GDF2pass.replace(False, 0).copy()  # replace booleans where column dtype is not bool but np.object
        GDF2pass = GDF2pass.replace(True, 1).copy()
        for col in GDF2pass.columns:
            if GDF2pass[col].dtype == bool:
                GDF2pass[col] = GDF2pass[col].astype(int)

        path_out = path_out if path_out else \
            get_generic_outpath(dir_out=os.path.join(self.dir_out, 'CoRegPoints'),
                                fName_out="CoRegPoints_grid%s_ws(%s_%s)__T_%s__R_%s.shp"
                                          % (self.grid_res, self.COREG_obj.win_size_XY[0],
                                             self.COREG_obj.win_size_XY[1], self.shift.basename, self.ref.basename))
        if not self.q:
            print('Writing %s ...' % path_out)
        GDF2pass.to_file(path_out)

    def _to_PointShapefile(self, skip_nodata=True, skip_nodata_col='ABS_SHIFT'):  # pragma: no cover
        warnings.warn(DeprecationWarning(
            "'_tiepoints_grid_to_PointShapefile' is deprecated."  # TODO delete if other method validated
            " 'tiepoints_grid_to_PointShapefile' is much faster."))
        if self.CoRegPoints_table.empty:
            raise RuntimeError('Cannot save a point shapefile because no tie points were found at all.')

        GDF = self.CoRegPoints_table
        GDF2pass = \
            GDF if not skip_nodata else \
            GDF[GDF[skip_nodata_col] != self.outFillVal]
        shapely_points = GDF2pass['geometry'].values.tolist()
        attr_dicts = [collections.OrderedDict(zip(GDF2pass.columns,
                                                  GDF2pass.loc[i].values))
                      for i in GDF2pass.index]

        fName_out = "CoRegPoints_grid%s_ws%s.shp" \
                    % (self.grid_res, self.COREG_obj.win_size_XY)
        path_out = os.path.join(self.dir_out, fName_out)
        write_shp(path_out,
                  shapely_points,
                  prj=self.COREG_obj.shift.prj,
                  attrDict=attr_dicts)

    def to_vectorfield(self, path_out: str = None, fmt: str = None, mode: str = 'md') -> GeoArray:
        """Save the calculated X-/Y-shifts to a 2-band raster file that can be used to visualize a vectorfield.

        NOTE: For example ArcGIS is able to visualize such 2-band raster files as a vectorfield.

        :param path_out:    the output path. If not given, it is automatically defined.
        :param fmt:         output raster format string
        :param mode:        The mode how the output is written ('uv' or 'md'; default: 'md')
                            - 'uv': outputs X-/Y shifts
                            - 'md': outputs magnitude and direction
        """
        assert mode in ['uv', 'md'], "'mode' must be either 'uv' (outputs X-/Y shifts) or 'md' " \
                                     "(outputs magnitude and direction)'. Got %s." % mode
        attr_b1 = 'X_SHIFT_M' if mode == 'uv' else 'ABS_SHIFT'
        attr_b2 = 'Y_SHIFT_M' if mode == 'uv' else 'ANGLE'

        if self.CoRegPoints_table.empty:
            raise RuntimeError('Cannot save the vector field because no tie points were found at all.')

        xshift_arr, gt, prj = points_to_raster(points=self.CoRegPoints_table['geometry'],
                                               values=self.CoRegPoints_table[attr_b1],
                                               tgt_res=self.shift.xgsd * self.grid_res,
                                               prj=self.CoRegPoints_table.crs.to_wkt(),
                                               fillVal=self.outFillVal)

        yshift_arr, gt, prj = points_to_raster(points=self.CoRegPoints_table['geometry'],
                                               values=self.CoRegPoints_table[attr_b2],
                                               tgt_res=self.shift.xgsd * self.grid_res,
                                               prj=self.CoRegPoints_table.crs.to_wkt(),
                                               fillVal=self.outFillVal)

        out_GA = GeoArray(np.dstack([xshift_arr, yshift_arr]), gt, prj, nodata=self.outFillVal)

        path_out = path_out if path_out else \
            get_generic_outpath(dir_out=os.path.join(self.dir_out, 'CoRegPoints'),
                                fName_out="CoRegVectorfield%s_ws(%s_%s)__T_%s__R_%s.tif"
                                          % (self.grid_res, self.COREG_obj.win_size_XY[0],
                                             self.COREG_obj.win_size_XY[1], self.shift.basename, self.ref.basename))

        out_GA.save(path_out, fmt=fmt if fmt else 'Gtiff')

        return out_GA

    def to_interpolated_raster(self,
                               metric: str = 'ABS_SHIFT',
                               method: str = 'RBF',
                               plot_result: bool = False,
                               lowres_spacing: int = 5,
                               v: bool = False
                               ) -> np.ndarray:
        """Interpolate the point data of the given metric into space.

        :param metric:          metric name to interpolate, i.e., one of the column names of
                                Tie_Point_Grid.CoRegPoints_table, e.g., 'ABS_SHIFT'.
        :param method:          interpolation algorithm
                                - 'RBF' (Radial Basis Function)
                                - 'GPR' (Gaussian Process Regression; equivalent to Simple Kriging)
                                - 'Kriging' (Ordinary Kriging based on pykrige)
        :param plot_result:     plot the result to assess the interpolation quality
        :param lowres_spacing:  by default, RBF, GPR, and Kriging run a lower resolution which is then linearly
                                interpolated to the full output image resolution. lowres_spacing defines the number of
                                pixels between the low resolution grid points
                                (higher values are faster but less accurate, default: 5)
        :param v:               enable verbose mode
        :return:    interpolation result as numpy array in the X/Y dimension of the target image of the co-registration
        """
        TPGI = Tie_Point_Grid_Interpolator(self, v=v)

        return TPGI.interpolate(metric=metric, method=method, plot_result=plot_result, lowres_spacing=lowres_spacing)

    @staticmethod
    def to_Raster_using_Kriging(*args, **kwargs):
        raise NotImplementedError('This method was removed in arosics version 1.8.2. To interpolate tie points into '
                                  'space, you may now use Tie_Point_Grid.to_interpolated_raster() instead.')


class Tie_Point_Refiner(object):
    """A class for performing outlier detection."""

    def __init__(self, GDF,
                 min_reliability=60,
                 rs_max_outlier: float = 10,
                 rs_tolerance: float = 2.5,
                 rs_max_iter: int = 15,
                 rs_exclude_previous_outliers: bool = True,
                 rs_timeout: float = 20,
                 rs_random_state: Optional[int] = 0,
                 q: bool = False):
        """Get an instance of Tie_Point_Refiner.

        :param GDF:                             GeoDataFrame like TiePointGrid.CoRegPoints_table containing all tie
                                                points to be filtered and the corresponding metadata
        :param min_reliability:                 minimum threshold for previously computed tie X/Y shift
                                                reliability (default: 60%)
        :param rs_max_outlier:                  RANSAC: maximum percentage of outliers to be detected
                                                (default: 10%)
        :param rs_tolerance:                    RANSAC: percentage tolerance for max_outlier_percentage
                                                (default: 2.5%)
        :param rs_max_iter:                     RANSAC: maximum iterations for finding the best RANSAC threshold
                                                (default: 15)
        :param rs_exclude_previous_outliers:    RANSAC: whether to exclude points that have been flagged as
                                                outlier by earlier filtering (default:True)
        :param rs_timeout:                      RANSAC: timeout for iteration loop in seconds (default: 20)
        :param rs_random_state:                 RANSAC random state (an integer corresponds to a fixed/pseudo-random
                                                state, None randomizes the result)

        :param q:
        """
        self.GDF = GDF.copy()
        self.min_reliability = min_reliability
        self.rs_max_outlier_percentage = rs_max_outlier
        self.rs_tolerance = rs_tolerance
        self.rs_max_iter = rs_max_iter
        self.rs_exclude_previous_outliers = rs_exclude_previous_outliers
        self.rs_timeout = rs_timeout
        self.rs_random_state = rs_random_state
        self.q = q
        self.new_cols = []
        self.ransac_model_robust = None

    def run_filtering(self, level=3):
        """Filter tie points used for shift correction.

        :param level:   tie point filter level (default: 3).
                        NOTE: lower levels are also included if a higher level is chosen

                        - Level 0: no tie point filtering
                        - Level 1: Reliablity filtering
                                   - filter all tie points out that have a low reliability according to internal tests
                        - Level 2: SSIM filtering
                                   - filters all tie points out where shift correction does not increase image
                                     similarity within matching window (measured by mean structural similarity index)
                        - Level 3: RANSAC outlier detection

        :return:
        """
        # TODO catch empty GDF

        # RELIABILITY filtering
        if level > 0:
            marked_recs = self._reliability_thresholding()  # type: Series
            self.GDF['L1_OUTLIER'] = marked_recs
            self.new_cols.append('L1_OUTLIER')

            n_flagged = len(marked_recs[marked_recs])
            perc40 = np.percentile(self.GDF.RELIABILITY, 40)

            if n_flagged / len(self.GDF) > .7:
                warnings.warn(r"More than 70%% of the found tie points have a reliability lower than %s%% and are "
                              r"therefore marked as false-positives. Consider relaxing the minimum reliability "
                              r"(parameter 'min_reliability') to avoid that. For example min_reliability=%d would only "
                              r"flag 40%% of the tie points in case of your input data."
                              % (self.min_reliability, perc40))

            if not self.q:
                print('%s tie points flagged by level 1 filtering (reliability).'
                      % n_flagged)

        # SSIM filtering
        if level > 1:
            marked_recs = self._SSIM_filtering()
            self.GDF['L2_OUTLIER'] = marked_recs  # type: Series
            self.new_cols.append('L2_OUTLIER')

            if not self.q:
                print('%s tie points flagged by level 2 filtering (SSIM).' % (len(marked_recs[marked_recs])))

        # RANSAC filtering
        if level > 2:
            # exclude previous outliers
            ransacInGDF = self.GDF[~self.GDF[self.new_cols].any(axis=1)].copy() \
                    if self.rs_exclude_previous_outliers else self.GDF

            if len(ransacInGDF) > 4:
                # running RANSAC with less than four tie points makes no sense

                marked_recs = self._RANSAC_outlier_detection(ransacInGDF)  # type: Series
                # we need to join a list here because otherwise it's merged by the 'index' column
                self.GDF['L3_OUTLIER'] = marked_recs.tolist()

                if not self.q:
                    print('%s tie points flagged by level 3 filtering (RANSAC)'
                          % (len(marked_recs[marked_recs])))
            else:
                print('RANSAC skipped because too less valid tie points have been found.')
                self.GDF['L3_OUTLIER'] = False

            self.new_cols.append('L3_OUTLIER')

        self.GDF['OUTLIER'] = self.GDF[self.new_cols].any(axis=1)
        self.new_cols.append('OUTLIER')

        return self.GDF, self.new_cols

    def _reliability_thresholding(self):
        """Exclude all records where estimated reliability of the calculated shifts is below the given threshold."""
        return self.GDF.RELIABILITY < self.min_reliability

    def _SSIM_filtering(self):
        """Exclude all records where SSIM decreased."""
        # ssim_diff  = np.median(self.GDF['SSIM_AFTER']) - np.median(self.GDF['SSIM_BEFORE'])

        # self.GDF.SSIM_IMPROVED = \
        #     self.GDF.apply(lambda GDF_row: GDF_row['SSIM_AFTER']>GDF_row['SSIM_BEFORE'] + ssim_diff, axis=1)

        return ~self.GDF.SSIM_IMPROVED

    def _RANSAC_outlier_detection(self, inGDF):
        """Detect geometric outliers between point cloud of source and estimated coordinates using RANSAC algorithm."""
        # from skimage.transform import PolynomialTransform  # import here to avoid static TLS ImportError

        src_coords = np.array(inGDF[['X_MAP', 'Y_MAP']])
        xyShift = np.array(inGDF[['X_SHIFT_M', 'Y_SHIFT_M']])
        est_coords = src_coords + xyShift

        for co, n in zip([src_coords, est_coords], ['src_coords', 'est_coords']):
            assert co.ndim == 2 and co.shape[1] == 2, "'%s' must have shape [Nx2]. Got shape %s." % (n, co.shape)

        if not 0 < self.rs_max_outlier_percentage < 100:
            raise ValueError
        min_inlier_percentage = 100 - self.rs_max_outlier_percentage

        # class PolyTF_1(PolynomialTransform):  # pragma: no cover
        #     def estimate(*data):
        #         return PolynomialTransform.estimate(*data, order=1)

        # robustly estimate affine transform model with RANSAC
        # eliminates not more than the given maximum outlier percentage of the tie points

        model_robust, inliers = None, None
        count_inliers = None
        th = 5  # start RANSAC threshold
        th_checked = {}  # dict of thresholds that already have been tried + calculated inlier percentage
        th_substract = 2
        count_iter = 0
        time_start = time()
        ideal_count = min_inlier_percentage * src_coords.shape[0] / 100

        # optimize RANSAC threshold so that it marks not much more or less than the given outlier percentage
        while True:
            if th_checked:
                th_too_strict = count_inliers < ideal_count  # True if too less inliers remaining

                # calculate new theshold using old increment
                # (but ensure th_new>0 by adjusting increment if needed)
                th_new = 0
                while th_new <= 0:
                    th_new = th + th_substract if th_too_strict else th - th_substract
                    if th_new <= 0:
                        th_substract /= 2

                # check if calculated new threshold has been used before
                th_already_checked = th_new in th_checked.keys()

                # if yes, decrease increment and recalculate new threshold
                th_substract = \
                    th_substract if not th_already_checked else \
                    th_substract / 2
                th = th_new if not th_already_checked else \
                    (th + th_substract if th_too_strict else
                     th - th_substract)

            ###############
            # RANSAC call #
            ###############

            # model_robust, inliers = ransac((src, dst),
            #                                PolynomialTransform,
            #                                min_samples=3)
            if src_coords.size and \
               est_coords.size and \
               src_coords.shape[0] > 6:
                # import here to avoid static TLS ImportError
                from skimage.measure import ransac
                from skimage.transform import AffineTransform

                model_robust, inliers = \
                    ransac((src_coords, est_coords),
                           AffineTransform,
                           min_samples=6,
                           residual_threshold=th,
                           max_trials=2000,
                           stop_sample_num=int(
                               (min_inlier_percentage - self.rs_tolerance) /
                               100 * src_coords.shape[0]
                           ),
                           stop_residuals_sum=int(
                               (self.rs_max_outlier_percentage - self.rs_tolerance) /
                               100 * src_coords.shape[0]),
                           random_state=self.rs_random_state
                           )
            else:
                warnings.warn('RANSAC filtering could not be applied '
                              'because there were too few tie points to fit a model.')
                inliers = np.array([])
                break

            count_inliers = np.count_nonzero(inliers)

            th_checked[th] = count_inliers / src_coords.shape[0] * 100
            # print(th,'\t', th_checked[th], )

            if min_inlier_percentage - self.rs_tolerance <\
               th_checked[th] <\
               min_inlier_percentage + self.rs_tolerance:
                # print('in tolerance')
                break

            if count_iter > self.rs_max_iter or \
               time() - time_start > self.rs_timeout:
                break  # keep last values and break while loop

            count_iter += 1

        outliers = inliers.__eq__(False) if inliers is not None and inliers.size else np.array([])

        if inGDF.empty or outliers is None or \
           (isinstance(outliers, list) and not outliers) or \
           (isinstance(outliers, np.ndarray) and not outliers.size):
            outseries = Series([False] * len(self.GDF))

        elif len(inGDF) < len(self.GDF):
            inGDF['outliers'] = outliers
            fullGDF = GeoDataFrame(self.GDF['POINT_ID'])
            fullGDF = fullGDF.merge(inGDF[['POINT_ID', 'outliers']],
                                    on='POINT_ID',
                                    how="outer")
            # fullGDF.outliers.copy()[~fullGDF.POINT_ID.isin(GDF.POINT_ID)] = False
            fullGDF = fullGDF.fillna(False)  # NaNs are due to exclude_previous_outliers
            outseries = fullGDF['outliers']

        else:
            outseries = Series(outliers)

        assert len(outseries) == len(self.GDF), \
            'RANSAC output validation failed.'

        self.ransac_model_robust = model_robust

        return outseries


class Tie_Point_Grid_Interpolator(object):
    """Class to interpolate tie point data into space."""

    def __init__(self, tiepointgrid: Tie_Point_Grid, v: bool = False) -> None:
        """Get an instance of Tie_Point_Grid_Interpolator.

        :param tiepointgrid:    instance of Tie_Point_Grid after computing spatial shifts
        :param v:               enable verbose mode
        """
        self.tpg = tiepointgrid
        self.v = v

    def interpolate(self,
                    metric: str,
                    method: str = 'RBF',
                    plot_result: bool = False,
                    lowres_spacing: int = 5
                    ) -> np.array:
        """Interpolate the point data of the given metric into space.

        :param metric:          metric name to interpolate, i.e., one of the column names of
                                Tie_Point_Grid.CoRegPoints_table, e.g., 'ABS_SHIFT'.
        :param method:          interpolation algorithm
                                - 'RBF' (Radial Basis Function)
                                - 'GPR' (Gaussian Process Regression; equivalent to Simple Kriging)
                                - 'Kriging' (Ordinary Kriging based on pykrige)
        :param plot_result:     plot the result to assess the interpolation quality
        :param lowres_spacing:  by default, RBF, GPR, and Kriging run a lower resolution which is then linearly
                                interpolated to the full output image resolution. lowres_spacing defines the number of
                                pixels between the low resolution grid points
                                (higher values are faster but less accurate, default: 5)
        :return:    interpolation result as numpy array in the X/Y dimension of the target image of the co-registration
        """
        t0 = time()

        rows, cols, data = self._get_pointdata(metric)
        nrows_out, ncols_out = self.tpg.shift.shape[:2]

        rows_lowres = np.arange(0, nrows_out + lowres_spacing, lowres_spacing)
        cols_lowres = np.arange(0, ncols_out + lowres_spacing, lowres_spacing)
        args = rows, cols, data, rows_lowres, cols_lowres

        if method == 'RBF':
            data_lowres = self._interpolate_via_rbf(*args)
        elif method == 'GPR':
            data_lowres = self._interpolate_via_gpr(*args)
        elif method == 'Kriging':
            data_lowres = self._interpolate_via_kriging(*args)
        else:
            raise ValueError(method)

        if lowres_spacing > 1:
            rows_full = np.arange(nrows_out)
            cols_full = np.arange(ncols_out)
            data_full = self._interpolate_regulargrid(rows_lowres, cols_lowres, data_lowres, rows_full, cols_full)
        else:
            data_full = data_lowres

        if self.v:
            print('interpolation runtime: %.2fs' % (time() - t0))
        if plot_result:
            self._plot_interpolation_result(data_full, rows, cols, data, metric)

        return data_full

    def _get_pointdata(self, metric: str):
        """Get the point data for the given metric from Tie_Point_Grid.CoRegPoints_table while ignoring outliers."""
        tiepoints = self.tpg.CoRegPoints_table[self.tpg.CoRegPoints_table.OUTLIER.__eq__(False)].copy()

        rows = np.array(tiepoints.Y_IM)
        cols = np.array(tiepoints.X_IM)
        data = np.array(tiepoints[metric])

        return rows, cols, data

    @staticmethod
    def _plot_interpolation_result(data_full: np.ndarray,
                                   rows: np.ndarray,
                                   cols: np.ndarray,
                                   data: np.ndarray,
                                   metric: str
                                   ):
        """Plot the interpolation result together with the input point data."""
        plt.figure(figsize=(7, 7))
        im = plt.imshow(data_full)
        plt.colorbar(im)
        plt.scatter(cols, rows, c=data, edgecolors='black')
        plt.title(metric)
        plt.show()

    @staticmethod
    def _interpolate_regulargrid(rows: np.ndarray,
                                 cols: np.ndarray,
                                 data: np.ndarray,
                                 rows_full: np.ndarray,
                                 cols_full: np.ndarray
                                 ):
        """Run linear regular grid interpolation."""
        RGI = RegularGridInterpolator(points=[cols, rows],
                                      values=data.T,  # must be in shape [x, y]
                                      method='linear',
                                      bounds_error=False)
        data_full = RGI(np.dstack(np.meshgrid(cols_full, rows_full)))
        return data_full

    @staticmethod
    def _interpolate_via_rbf(rows: np.ndarray,
                             cols: np.ndarray,
                             data: np.ndarray,
                             rows_full: np.ndarray,
                             cols_full: np.ndarray
                             ):
        """Run Radial Basis Function (RBF) interpolation.

        -> https://github.com/agile-geoscience/xlines/blob/master/notebooks/11_Gridding_map_data.ipynb
        -> documents the legacy scipy.interpolate.Rbf
        """
        rbf = RBFInterpolator(
            np.column_stack([cols, rows]), data,
            kernel="linear",
            # kernel="thin_plate_spline",
        )
        cols_grid, rows_grid = np.meshgrid(cols_full, rows_full)
        data_full = \
            rbf(np.column_stack([cols_grid.flat, rows_grid.flat]))\
            .reshape(rows_grid.shape)

        return data_full

    @staticmethod
    def _interpolate_via_gpr(rows: np.ndarray,
                             cols: np.ndarray,
                             data: np.ndarray,
                             rows_full: np.ndarray,
                             cols_full: np.ndarray
                             ):
        """Run Gaussian Process Regression (GPR) interpolation.

        -> https://stackoverflow.com/questions/24978052/interpolation-over-regular-grid-in-python
        """
        try:
            import sklearn  # noqa F401
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                "GPR interpolation requires the optional package 'scikit-learn' to be installed. You may install it "
                "with Conda (conda install -c conda-forge scikit-learn) or Pip (pip install scikit-learn)."
            )

        from sklearn.gaussian_process.kernels import RBF
        from sklearn.gaussian_process import GaussianProcessRegressor

        gp = GaussianProcessRegressor(
            normalize_y=False,
            alpha=0.001,  # Larger values imply more input noise and result in smoother grids; default: 1e-10
            kernel=RBF(length_scale=100))
        gp.fit(np.column_stack([cols, rows]), data.T)

        cols_grid, rows_grid = np.meshgrid(cols_full, rows_full)
        data_full = \
            gp.predict(np.column_stack([cols_grid.flat, rows_grid.flat]))\
            .reshape(rows_grid.shape)

        return data_full

    @staticmethod
    def _interpolate_via_kriging(rows: np.ndarray,
                                 cols: np.ndarray,
                                 data: np.ndarray,
                                 rows_full: np.ndarray,
                                 cols_full: np.ndarray
                                 ):
        """Run Ordinary Kriging interpolation based on pykrige.

        Reference: P.K. Kitanidis, Introduction to Geostatistics: Applications in Hydrogeology,
                   (Cambridge University Press, 1997) 272 p.
        """
        try:
            import pykrige  # noqa F401
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                "Ordinary Kriging requires the optional package 'pykrige' to be installed. You may install it with "
                "Conda (conda install -c conda-forge pykrige) or Pip (pip install pykrige)."
            )

        from pykrige.ok import OrdinaryKriging

        OK = OrdinaryKriging(cols.astype(float), rows.astype(float), data.astype(float),
                             variogram_model='spherical',
                             verbose=False)

        data_full, sigmasq = \
            OK.execute('grid',
                       cols_full.astype(float),
                       rows_full.astype(float),
                       backend='C',
                       # n_closest_points=12
                       )

        return data_full
