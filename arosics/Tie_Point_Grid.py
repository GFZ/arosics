# -*- coding: utf-8 -*-

# AROSICS - Automated and Robust Open-Source Image Co-Registration Software
#
# Copyright (C) 2017-2021  Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
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

import collections
import multiprocessing
import os
import warnings
import time
from typing import Optional, Union

# custom
from osgeo import gdal
import numpy as np
from geopandas import GeoDataFrame
from pandas import DataFrame, Series
from shapely.geometry import Point

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
    grid as GCPs. Thus 'Tie_Point_Grid' can be used to correct for locally varying geometric distortions of the target
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
            rs_timeout, q. See documentation there.

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
        self.kriged = None  # set by Raster_using_Kriging()

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
        GDF.loc[:, ['X_IM', 'Y_IM']] = self.XY_points
        GDF.loc[:, ['X_MAP', 'Y_MAP']] = self.XY_mapPoints

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
                        time.sleep(.1)
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

        ssim_col = np.array(tbl['SSIM_AFTER' if after_correction else 'SSIM_BEFORE'])
        ssim_col = [i * i for i in ssim_col if i != self.outFillVal]

        return float(np.median(ssim_col))

    def plot_shift_distribution(self,
                                include_outliers: bool = True,
                                unit: str = 'm',
                                interactive: bool = False,
                                figsize: tuple = None,
                                xlim: list = None,
                                ylim: list = None,
                                fontsize: int = 12,
                                title: str = 'shift distribution'
                                ) -> tuple:
        """Create a 2D scatterplot containing the distribution of calculated X/Y-shifts.

        :param include_outliers:    whether to include tie points that have been marked as false-positives
        :param unit:                'm' for meters or 'px' for pixels (default: 'm')
        :param interactive:         interactive mode uses plotly for visualization
        :param figsize:             (xdim, ydim)
        :param xlim:                [xmin, xmax]
        :param ylim:                [ymin, ymax]
        :param fontsize:            size of all used fonts
        :param title:               the title to be plotted above the figure
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
            plt.text(xlim[1] - (xlim[1] / 20), -ylim[1] + (ylim[1] / 20),
                     'RMSE:  %s m / %s px' % (np.round(rmse, 2), np.round(rmse / self.shift.xgsd, 2)),
                     ha='right', va='bottom', fontsize=fontsize, bbox=dict(facecolor='w', pad=None, alpha=0.8))

            # add grid and increase linewidth of middle line
            plt.grid()
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
            plt.xlabel('x-shift [%s]' % 'meters' if unit == 'm' else 'pixels', fontsize=fontsize)
            plt.ylabel('y-shift [%s]' % 'meters' if unit == 'm' else 'pixels', fontsize=fontsize)

            # add legend with labels in the right order
            handles, labels = ax.get_legend_handles_labels()
            leg = plt.legend(reversed(handles), reversed(labels), fontsize=fontsize, loc='upper right', scatterpoints=3)
            leg.get_frame().set_edgecolor('black')

            plt.show()

            return fig, ax

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
        """Write the calculated tie points grid to a point shapefile (e.g., for visualization by a GIS software).

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

    def to_Raster_using_Kriging(self, attrName, skip_nodata=1, skip_nodata_col='ABS_SHIFT', outGridRes=None,
                                fName_out=None, tilepos=None, tilesize=500, mp=None):

        mp = False if self.CPUs == 1 else True
        self._Kriging_sp(attrName, skip_nodata=skip_nodata, skip_nodata_col=skip_nodata_col,
                         outGridRes=outGridRes, fName_out=fName_out, tilepos=tilepos)

        # if mp:
        #     tilepositions = UTL.get_image_tileborders([tilesize,tilesize],self.tgt_shape)
        #     args_kwargs_dicts=[]
        #     for tp in tilepositions:
        #         kwargs_dict = {'skip_nodata':skip_nodata,'skip_nodata_col':skip_nodata_col,'outGridRes':outGridRes,
        #                        'fName_out':fName_out,'tilepos':tp}
        #         args_kwargs_dicts.append({'args':[attrName],'kwargs':kwargs_dict})
        #     # self.kriged=[]
        #     # for i in args_kwargs_dicts:
        #     #     res = self.Kriging_mp(i)
        #     #     self.kriged.append(res)
        #     #     print(res)
        #
        #     with multiprocessing.Pool() as pool:
        #        self.kriged = pool.map(self.Kriging_mp,args_kwargs_dicts)
        # else:
        #     self.Kriging_sp(attrName,skip_nodata=skip_nodata,skip_nodata_col=skip_nodata_col,
        #                     outGridRes=outGridRes,fName_out=fName_out,tilepos=tilepos)
        res = self.kriged if mp else None
        return res

    def _Kriging_sp(self, attrName, skip_nodata=1, skip_nodata_col='ABS_SHIFT', outGridRes=None,
                    fName_out=None, tilepos=None):
        GDF = self.CoRegPoints_table
        GDF2pass = GDF if not skip_nodata else GDF[GDF[skip_nodata_col] != self.outFillVal]

        X_coords, Y_coords, ABS_SHIFT = GDF2pass['X_MAP'], GDF2pass['Y_MAP'], GDF2pass[attrName]

        xmin, ymin, xmax, ymax = GDF2pass.total_bounds

        grid_res = outGridRes if outGridRes else int(min(xmax - xmin, ymax - ymin) / 250)
        grid_x, grid_y = np.arange(xmin, xmax + grid_res, grid_res), np.arange(ymax, ymin - grid_res, -grid_res)

        # Reference: P.K. Kitanidis, Introduction to Geostatistcs: Applications in Hydrogeology,
        #            (Cambridge University Press, 1997) 272 p.
        from pykrige.ok import OrdinaryKriging
        OK = OrdinaryKriging(X_coords, Y_coords, ABS_SHIFT, variogram_model='spherical', verbose=False)
        zvalues, sigmasq = OK.execute('grid', grid_x, grid_y, backend='C', n_closest_points=12)

        if self.CPUs is None or self.CPUs > 1:
            fName_out = fName_out if fName_out else \
                "Kriging__%s__grid%s_ws%s_%s.tif" % (attrName, self.grid_res, self.COREG_obj.win_size_XY, tilepos)
        else:
            fName_out = fName_out if fName_out else \
                "Kriging__%s__grid%s_ws%s.tif" % (attrName, self.grid_res, self.COREG_obj.win_size_XY)
        path_out = get_generic_outpath(dir_out=self.dir_out, fName_out=fName_out)
        # add a half pixel grid points are centered on the output pixels
        xmin = xmin - grid_res / 2
        # ymin = ymin - grid_res / 2
        # xmax = xmax + grid_res / 2
        ymax = ymax + grid_res / 2

        GeoArray(zvalues,
                 geotransform=(xmin, grid_res, 0, ymax, 0, -grid_res),
                 projection=self.COREG_obj.shift.prj).save(path_out)

        return zvalues

    def _Kriging_mp(self, args_kwargs_dict):
        args = args_kwargs_dict.get('args', [])
        kwargs = args_kwargs_dict.get('kwargs', [])

        return self._Kriging_sp(*args, **kwargs)


class Tie_Point_Refiner(object):
    """A class for performing outlier detection."""

    def __init__(self, GDF,
                 min_reliability=60,
                 rs_max_outlier: float = 10,
                 rs_tolerance: float = 2.5,
                 rs_max_iter: int = 15,
                 rs_exclude_previous_outliers: bool = True,
                 rs_timeout: float = 20,
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

        :param q:
        """
        self.GDF = GDF.copy()
        self.min_reliability = min_reliability
        self.rs_max_outlier_percentage = rs_max_outlier
        self.rs_tolerance = rs_tolerance
        self.rs_max_iter = rs_max_iter
        self.rs_exclude_previous_outliers = rs_exclude_previous_outliers
        self.rs_timeout = rs_timeout
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
        time_start = time.time()
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
                               100 * src_coords.shape[0])
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
               time.time() - time_start > self.rs_timeout:
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
