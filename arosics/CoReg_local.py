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

import warnings
import os
from copy import copy
from typing import Tuple, Union, Optional
from collections import OrderedDict

# custom
from osgeo import gdal
try:
    import pyfftw
except ImportError:
    pyfftw = None
import numpy as np

from matplotlib import pyplot as plt  # noqa F401
from geopandas import GeoDataFrame  # noqa F401

from .Tie_Point_Grid import Tie_Point_Grid
from .CoReg import COREG
from .DeShifter import DESHIFTER
from .geometry import has_metaRotation, remove_metaRotation
from py_tools_ds.geo.coord_trafo import transform_any_prj, reproject_shapelyGeometry
from py_tools_ds.geo.map_info import geotransform2mapinfo
from geoarray import GeoArray

__author__ = 'Daniel Scheffler'


class COREG_LOCAL(object):
    """
    COREG_LOCAL applies the algorithm to detect spatial shifts to the whole overlap area of the input images.

    Spatial shifts are calculated for each point in grid of which the parameters can be adjusted using keyword
    arguments. Shift correction performs a polynomial transformation using the calculated shifts of each point in the
    grid as GCPs. Thus this class can be used to correct for locally varying geometric distortions of the target image.

    See help(COREG_LOCAL) for documentation.
    """

    def __init__(self,
                 im_ref: Union[GeoArray, str],
                 im_tgt: Union[GeoArray, str],
                 grid_res: float,
                 max_points: int = None,
                 window_size: Tuple[int, int] = (256, 256),
                 path_out: str = None,
                 fmt_out: str = 'ENVI',
                 out_crea_options: list = None,
                 projectDir: str = None,
                 r_b4match: int = 1,
                 s_b4match: int = 1,
                 max_iter: int = 5,
                 max_shift: int = 5,
                 tieP_filter_level: int = 3,
                 min_reliability: float = 60,
                 rs_max_outlier: float = 10,
                 rs_tolerance: float = 2.5,
                 align_grids: bool = True,
                 match_gsd: bool = False,
                 out_gsd: float = None,
                 target_xyGrid=None,
                 resamp_alg_deshift: str = 'cubic',
                 resamp_alg_calc: str = 'cubic',
                 footprint_poly_ref: str = None,
                 footprint_poly_tgt: str = None,
                 data_corners_ref: list = None,
                 data_corners_tgt: list = None,
                 outFillVal: int = -9999,
                 nodata: Tuple[int, int] = (None, None),
                 calc_corners: bool = True,
                 binary_ws: bool = True,
                 force_quadratic_win: bool = True,
                 mask_baddata_ref: Union[GeoArray, str] = None,
                 mask_baddata_tgt: Union[GeoArray, str] = None,
                 CPUs: int = None,
                 progress: bool = True,
                 v: bool = False,
                 q: bool = False,
                 ignore_errors: bool = True
                 ) -> None:
        """
        Get an instance of COREG_LOCAL.

        :param im_ref:
            source path of reference image (any GDAL compatible image format is supported)

        :param im_tgt:
            source path of image to be shifted (any GDAL compatible image format is supported)

        :param grid_res:
            tie point grid resolution in pixels of the target image (x-direction)

        :param max_points:
            maximum number of points used to find coregistration tie points

            NOTE: Points are selected randomly from the given point grid (specified by 'grid_res').
                  If the point grid does not provide enough points, all available points are chosen.

        :param window_size:
            custom matching window size [pixels] (default: (256,256))

        :param path_out:
            target path of the coregistered image
            - if None (default), no output is written to disk
            - if 'auto': /dir/of/im1/<im1>__shifted_to__<im0>.bsq

        :param fmt_out:
            raster file format for output file. ignored if path_out is None. Can be any GDAL compatible raster file
            format (e.g. 'ENVI', 'GTIFF'; default: ENVI). Refer to https://gdal.org/drivers/raster/index.html to get a
            full list of supported formats.

        :param out_crea_options:
            GDAL creation options for the output image, e.g. ["QUALITY=80", "REVERSIBLE=YES", "WRITE_METADATA=YES"]

        :param projectDir:
            name of a project directory where to store all the output results. If given, name is inserted into all
            automatically generated output paths.

        :param r_b4match:
            band of reference image to be used for matching (starts with 1; default: 1)

        :param s_b4match:
            band of shift image to be used for matching (starts with 1; default: 1)

        :param max_iter:
            maximum number of iterations for matching (default: 5)

        :param max_shift:
            maximum shift distance in reference image pixel units (default: 5 px)

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

        :param min_reliability:
            Tie point filtering: minimum reliability threshold, below which tie points are marked as false-positives
            (default: 60%)
            - accepts values between 0% (no reliability) and 100 % (perfect reliability)
            HINT: decrease this value in case of poor signal-to-noise ratio of your input data

        :param rs_max_outlier:
            RANSAC tie point filtering: proportion of expected outliers (default: 10%)

        :param rs_tolerance:
            RANSAC tie point filtering: percentage tolerance for max_outlier_percentage (default: 2.5%)

        :param align_grids:
            True: align the input coordinate grid to the reference (does not affect the output pixel size as long as
            input and output pixel sizes are compatible (5:30 or 10:30 but not 4:30), default = True

        :param match_gsd:
            True: match the input pixel size to the reference pixel size,
            default = False

        :param out_gsd:
            output pixel size in units of the reference coordinate system (default = pixel size of the input array),
            given values are overridden by match_gsd=True

        :param target_xyGrid:
            a list with a target x-grid and a target y-grid like [[15,45], [15,45]]
            This overrides 'out_gsd', 'align_grids' and 'match_gsd'.

        :param resamp_alg_deshift:
            the resampling algorithm to be used for shift correction (if neccessary)
            valid algorithms: nearest, bilinear, cubic, cubic_spline, lanczos, average, mode, max, min, med, q1, q3
            (default: cubic)

        :param resamp_alg_calc:
            the resampling algorithm to be used for all warping processes during calculation of spatial shifts
            valid algorithms: nearest, bilinear, cubic, cubic_spline, lanczos, average, mode, max, min, med, q1, q3
            (default: cubic (highly recommended))

        :param footprint_poly_ref:
            footprint polygon of the reference image (WKT string or shapely.geometry.Polygon),
            e.g. 'POLYGON ((299999 6000000, 299999 5890200, 409799 5890200, 409799 6000000, 299999 6000000))'

        :param footprint_poly_tgt:
            footprint polygon of the image to be shifted (WKT string or shapely.geometry.Polygon)
            e.g. 'POLYGON ((299999 6000000, 299999 5890200, 409799 5890200, 409799 6000000, 299999 6000000))'

        :param data_corners_ref:
            map coordinates of data corners within reference image. ignored if footprint_poly_ref is given.

        :param data_corners_tgt:
            map coordinates of data corners within image to be shifted. ignored if footprint_poly_tgt is given.

        :param outFillVal:
            if given the generated tie point grid is filled with this value in case no match could be found during
            co-registration (default: -9999)

        :param nodata:
            no data values for reference image and image to be shifted

        :param calc_corners:
            calculate true positions of the dataset corners in order to get a useful matching window position within
            the actual image overlap
            (default: True; deactivated if 'data_corners_im0' and 'data_corners_im1' are given)

        :param binary_ws:
            use binary X/Y dimensions for the matching window (default: True)

        :param force_quadratic_win:
            force a quadratic matching window (default: 1)

        :param mask_baddata_ref:
            path to a 2D boolean mask file (or an instance of BadDataMask) for the reference image where all bad data
            pixels (e.g. clouds) are marked with True and the remaining pixels with False. Must have the same
            geographic extent and projection like 'im_ref'. The mask is used to check if the chosen matching window
            position is valid in the sense of useful data. Otherwise this window position is rejected.

        :param mask_baddata_tgt:
            path to a 2D boolean mask file (or an instance of BadDataMask) for the image to be shifted where all bad
            data pixels (e.g. clouds) are marked with True and the remaining pixels with False. Must have the same
            geographic extent and projection like 'im_ref'. The mask is used to check if the chosen matching window
            position is valid in the sense of useful data. Otherwise this window position is rejected.

        :param CPUs:
            number of CPUs to use during calculation of tie point grid (default: None, which means 'all CPUs available')

        :param progress:
            show progress bars (default: True)

        :param v:
            verbose mode (default: False)

        :param q:
            quiet mode (default: False)

        :param ignore_errors:
            Useful for batch processing. (default: False)
        """
        # assertions / input validation
        assert gdal.GetDriverByName(fmt_out), "'%s' is not a supported GDAL driver." % fmt_out
        if match_gsd and out_gsd:
            warnings.warn("'-out_gsd' is ignored because '-match_gsd' is set.\n")
        if out_gsd:
            assert isinstance(out_gsd, list) and len(out_gsd) == 2, 'out_gsd must be a list with two values.'

        self.params = dict([x for x in locals().items() if x[0] != "self" and not x[0].startswith('__')])

        # NOTE: self.imref and self.im2shift are handled completely independent from self.COREG_obj.ref and
        #       self.COREG_obj.shift. self.COREG_obj.ref and self.COREG_obj.shift are used for shift calculation and
        #       correction is applied to self.im2shift.
        self.imref = GeoArray(im_ref, nodata=nodata[0], progress=progress, q=q)
        self.im2shift = GeoArray(im_tgt, nodata=nodata[1], progress=progress, q=q)
        self.path_out = path_out  # updated by self.set_outpathes
        self.fmt_out = fmt_out
        self.out_creaOpt = out_crea_options
        self._projectDir = projectDir
        self.grid_res = grid_res
        self.max_points = max_points
        self.window_size = window_size
        self.max_shift = max_shift
        self.max_iter = max_iter
        self.tieP_filter_level = tieP_filter_level
        self.min_reliability = min_reliability
        self.rs_max_outlier = rs_max_outlier
        self.rs_tolerance = rs_tolerance
        self.align_grids = align_grids
        self.match_gsd = match_gsd
        self.out_gsd = out_gsd
        self.target_xyGrid = target_xyGrid
        self.rspAlg_DS = resamp_alg_deshift  # TODO convert integers to strings
        self.rspAlg_calc = resamp_alg_calc
        self.calc_corners = calc_corners
        self.nodata = nodata
        self.outFillVal = outFillVal
        self.bin_ws = binary_ws
        self.force_quadratic_win = force_quadratic_win
        self.CPUs = CPUs
        self.path_verbose_out = ''  # TODO
        self.v = v
        self.q = q if not v else False  # overridden by v
        self.progress = progress if not q else False  # overridden by v
        self.ignErr = ignore_errors  # FIXME this is not yet implemented for COREG_LOCAL

        assert self.tieP_filter_level in range(4), 'Invalid tie point filter level.'
        assert isinstance(self.imref, GeoArray) and isinstance(self.im2shift, GeoArray), \
            'Something went wrong with the creation of GeoArray instances for reference or target image. The created ' \
            'instances do not seem to belong to the GeoArray class. If you are working in Jupyter Notebook, reset ' \
            'the kernel and try again.'

        COREG.__dict__['_set_outpathes'](self, self.imref, self.im2shift)
        # make sure that the output directory of coregistered image is the project directory if a project directory is
        # given
        if path_out and projectDir and os.path.basename(self.path_out):
            self.path_out = os.path.join(self.projectDir, os.path.basename(self.path_out))

        gdal.AllRegister()

        # resample input data in case there is a metadata rotation (not handled by AROSICS)
        self._check_and_handle_metaRotation()

        try:
            # ignore_errors must be False because in case COREG init fails, coregistration for the whole scene fails
            self.COREG_obj = COREG(self.imref, self.im2shift,
                                   ws=window_size,
                                   footprint_poly_ref=footprint_poly_ref,
                                   footprint_poly_tgt=footprint_poly_tgt,
                                   data_corners_ref=data_corners_ref,
                                   data_corners_tgt=data_corners_tgt,
                                   resamp_alg_calc=self.rspAlg_calc,
                                   calc_corners=calc_corners,
                                   r_b4match=r_b4match,
                                   s_b4match=s_b4match,
                                   max_iter=max_iter,
                                   max_shift=max_shift,
                                   nodata=nodata,
                                   mask_baddata_ref=None,  # see below
                                   mask_baddata_tgt=None,
                                   CPUs=self.CPUs,
                                   force_quadratic_win=self.force_quadratic_win,
                                   binary_ws=self.bin_ws,
                                   progress=self.progress,
                                   v=v,
                                   q=q,
                                   ignore_errors=False)
        except Exception:
            warnings.warn('\nFirst attempt to check if functionality of co-registration failed. Check your '
                          'input data and parameters. The following error occurred:', stacklevel=3)
            raise

        if pyfftw:
            self.check_if_fftw_works()

        # add bad data mask
        # (mask is not added during initialization of COREG object in order to avoid bad data area errors there)
        if mask_baddata_ref is not None:
            self.COREG_obj.ref.mask_baddata = mask_baddata_ref
        if mask_baddata_tgt is not None:
            self.COREG_obj.shift.mask_baddata = mask_baddata_tgt

        self._tiepoint_grid = None  # set by self.tiepoint_grid
        self._CoRegPoints_table = None  # set by self.CoRegPoints_table
        self._coreg_info = None  # set by self.coreg_info
        self.deshift_results = None  # set by self.correct_shifts()
        self._success = None  # set by self.success property

    def _check_and_handle_metaRotation(self):
        """Check if the provided input data have a metadata rotation and if yes, correct it AND equalize grids.

        In case there is a rotation, the GDAL GeoTransform is not 0 at positions 2 or 4. So far, AROSICS does not
        handle such rotations, so the resampling is needed to make things work. The pixel grid equalization is also
        done here to avoid a double-resampling (grid would be equalized by COREG.equalize_pixGrids() otherwise).
        """
        grid2use = 'ref' if self.im2shift.xgsd <= self.imref.xgsd else 'shift'

        if has_metaRotation(self.imref) or has_metaRotation(self.im2shift):

            msg = 'The %s image needs to be resampled because it has a row/column rotation in '\
                  'its map info which is not handled by AROSICS.'

            if grid2use == 'ref':
                if has_metaRotation(self.imref):
                    warnings.warn(msg % 'reference')
                    self.imref = remove_metaRotation(self.imref)

                # resample target to reference image
                if not self.q:
                    print('Adapting the target image pixel grid to the one of the reference image for shift detection.')
                self.im2shift.reproject_to_new_grid(prototype=self.imref, CPUs=self.CPUs)

            else:
                # remove any metadata rotation (a rotation that only exists in the map info)
                if has_metaRotation(self.im2shift):
                    warnings.warn(msg % 'target')
                    self.im2shift = remove_metaRotation(self.im2shift)

                # resample reference to target image
                print('Adapting the reference image pixel grid to the one of the target image for shift detection.')
                self.imref.reproject_to_new_grid(prototype=self.im2shift, CPUs=self.CPUs)

    def check_if_fftw_works(self) -> None:
        """Assign the attribute 'fftw_works' to self.COREG_obj by executing shift calculation once with muted output."""
        # calculate global shift once in order to check is fftw works
        try:
            self.COREG_obj.q = True
            self.COREG_obj.v = False
            self.COREG_obj.calculate_spatial_shifts()
        except RuntimeError:
            if self.COREG_obj.fftw_works is not None:
                pass
            else:
                warnings.warn('\nFirst attempt to check if functionality of co-registration failed. Check your '
                              'input data and parameters. The following error occurred:', stacklevel=3)
                raise

        self.COREG_obj.q = self.q
        self.COREG_obj.v = self.v

    @property
    def projectDir(self) -> str:
        if self._projectDir:
            if len(os.path.split(self._projectDir)) == 1:
                return os.path.abspath(os.path.join(os.path.curdir, self._projectDir))
            else:
                return os.path.abspath(self._projectDir)
        else:
            # return a project name that not already has a corresponding folder on disk
            root_dir = os.path.dirname(self.im2shift.filePath) if self.im2shift.filePath else os.path.curdir
            fold_name = 'UntitledProject_1'

            while os.path.isdir(os.path.join(root_dir, fold_name)):
                fold_name = '%s_%s' % (fold_name.split('_')[0], int(fold_name.split('_')[-1]) + 1)

            self._projectDir = os.path.join(root_dir, fold_name)
            return self._projectDir

    @property
    def tiepoint_grid(self) -> Tie_Point_Grid:
        if self._tiepoint_grid:
            return self._tiepoint_grid
        else:
            self.calculate_spatial_shifts()
            return self._tiepoint_grid

    @property
    def CoRegPoints_table(self) -> GeoDataFrame:
        """Return a GeoDataFrame containing all the results from coregistration for all points in the tie point grid.

        Columns of the GeoDataFrame: 'geometry','POINT_ID','X_IM','Y_IM','X_MAP','Y_MAP','X_WIN_SIZE', 'Y_WIN_SIZE',
                                     'X_SHIFT_PX','Y_SHIFT_PX', 'X_SHIFT_M', 'Y_SHIFT_M', 'ABS_SHIFT' and 'ANGLE'
        """
        return self.tiepoint_grid.CoRegPoints_table

    @property
    def success(self) -> bool:
        self._success = self.tiepoint_grid.GCPList != []
        return self._success

    def calculate_spatial_shifts(self) -> None:
        self._tiepoint_grid = \
            Tie_Point_Grid(self.COREG_obj, self.grid_res,
                           max_points=self.max_points,
                           outFillVal=self.outFillVal,
                           resamp_alg_calc=self.rspAlg_calc,
                           tieP_filter_level=self.tieP_filter_level,
                           outlDetect_settings=dict(
                               min_reliability=self.min_reliability,
                               rs_max_outlier=self.rs_max_outlier,
                               rs_tolerance=self.rs_tolerance),
                           dir_out=self.projectDir,
                           CPUs=self.CPUs,
                           progress=self.progress,
                           v=self.v,
                           q=self.q)
        self._tiepoint_grid.get_CoRegPoints_table()

        if self.v:
            print('Visualizing CoReg points grid...')
            self.view_CoRegPoints(figsize=(10, 10))

    def show_image_footprints(self):
        """Show a web map containing the calculated footprints and overlap area of the input images.

        NOTE: This method is intended to be called from Jupyter Notebook.
        """
        return self.COREG_obj.show_image_footprints()

    def view_CoRegPoints(self,
                         shapes2plot: str = 'points',
                         attribute2plot: str = 'ABS_SHIFT',
                         cmap: plt.cm = None,
                         exclude_fillVals: bool = True,
                         backgroundIm: str = 'tgt',
                         hide_filtered: bool = True,
                         figsize: tuple = None,
                         title: str = '',
                         vector_scale: float = 1.,
                         savefigPath: str = '',
                         savefigDPI: int = 96,
                         showFig: bool = True,
                         vmin: float = None,
                         vmax: float = None,
                         return_map: bool = False
                         ) -> Optional[Tuple]:
        """
        Show a map of the calculated tie point grid with the target image as background.

        :param shapes2plot:         'points': plot points representing values of 'attribute2plot' onto the map
                                    'vectors': plot shift vectors onto the map
        :param attribute2plot:      the attribute of the tie point grid to be shown (default: 'ABS_SHIFT')
        :param cmap:                a custom color map to be applied to the plotted grid points (default: 'RdYlGn_r')
        :param exclude_fillVals:    whether to exclude those points of the grid where spatial shift detection failed
        :param backgroundIm:        whether to use the target or the reference image as map background. Possible
                                    options are 'ref' and 'tgt' (default: 'tgt')
        :param hide_filtered:       hide all points that have been filtered out according to tie point filter level
        :param figsize:             size of the figure to be viewed, e.g. (10, 10)
        :param title:               plot title
        :param vector_scale:        scale factor for shift vector length (default: 1 -> no scaling)
        :param savefigPath:
        :param savefigDPI:
        :param showFig:             whether to show or to hide the figure
        :param vmin:
        :param vmax:
        :param return_map:
        :return:
        """
        from matplotlib import pyplot as plt  # noqa
        from matplotlib.offsetbox import AnchoredText
        from cartopy.crs import PlateCarree
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        # get a map showing the reference or target image
        if backgroundIm not in ['tgt', 'ref']:
            raise ValueError('backgroundIm')
        backgroundIm = self.im2shift if backgroundIm == 'tgt' else self.imref
        fig, ax = backgroundIm.show_map(figsize=figsize,
                                        nodataVal=self.nodata[1],
                                        return_map=True,
                                        band=self.COREG_obj.shift.band4match)

        # make sure the output figure has a reasonable size, also if figsize is not given
        if not figsize:
            w, h = fig.get_size_inches()
            fig.set_size_inches(w * 1.6, h * 1.6)

        # set figure title
        dict_attr_title = dict(
            X_WIN_SIZE='size of the matching window in x-direction [pixels]',
            Y_WIN_SIZE='size of the matching window in y-direction [pixels]',
            X_SHIFT_PX='absolute shifts in x-direction [pixels]',
            Y_SHIFT_PX='absolute shifts in y-direction [pixels]',
            X_SHIFT_M='absolute shifts in x-direction [map units]',
            Y_SHIFT_M='absolute shifts in y-direction [map units]',
            ABS_SHIFT='absolute shift vector length [map units]',
            ANGLE='shift vector direction [angle in degrees]',
            SSIM_BEFORE='structural similarity index before co-registration',
            SSIM_AFTER='structural similarity index after co-registration',
            SSIM_IMPROVED='structural similarity index improvement through co-registration [yes/no]',
            RELIABILITY='reliability of the computed shift vector'
        )
        if title:
            ax.set_title(title)
        elif attribute2plot in dict_attr_title:
            ax.set_title(dict_attr_title[attribute2plot], pad=20)
        elif attribute2plot in self.CoRegPoints_table.columns:
            ax.set_title(attribute2plot)
        else:
            raise ValueError(attribute2plot, "Invalid value for 'attribute2plot'. Valid values are: %s."
                             % ", ".join(self.CoRegPoints_table.columns))

        if not self.CoRegPoints_table.empty:
            # get GeoDataFrame containing everything needed for plotting
            outlierCols = [c for c in self.CoRegPoints_table.columns if 'OUTLIER' in c]
            attr2include = ['geometry', attribute2plot] + outlierCols + ['X_SHIFT_M', 'Y_SHIFT_M']
            GDF = self.CoRegPoints_table.loc[self.CoRegPoints_table.X_SHIFT_M != self.outFillVal, attr2include].copy()\
                if exclude_fillVals else self.CoRegPoints_table.loc[:, attr2include]

            # get LonLat coordinates for all points
            XY = np.array([geom.coords.xy for geom in GDF.geometry]).reshape(-1, 2)
            lon, lat = transform_any_prj(self.im2shift.projection, 4326, XY[:, 0], XY[:, 1])
            GDF['Lon'], GDF['Lat'] = lon, lat

            # get colors for all points
            palette = cmap if cmap is not None else plt.cm.get_cmap('RdYlGn_r')
            if cmap is None and attribute2plot == 'ANGLE':
                import cmocean
                palette = getattr(cmocean.cm, 'delta')

            if hide_filtered:
                if self.tieP_filter_level > 0:
                    GDF = GDF[GDF.L1_OUTLIER.__eq__(False)].copy()
                if self.tieP_filter_level > 1:
                    GDF = GDF[GDF.L2_OUTLIER.__eq__(False)].copy()
                if self.tieP_filter_level > 2:
                    GDF = GDF[GDF.L3_OUTLIER.__eq__(False)].copy()
            else:
                marker = 'o' if len(GDF) < 10000 else '.'
                common_kw = dict(marker=marker, alpha=1.0, transform=PlateCarree())
                if self.tieP_filter_level > 0:
                    # flag level 1 outliers
                    GDF_filt = GDF[GDF.L1_OUTLIER.__eq__(True)].copy()
                    ax.scatter(GDF_filt['Lon'], GDF_filt['Lat'], c='b', s=250, label='reliability',
                               **common_kw)
                if self.tieP_filter_level > 1:
                    # flag level 2 outliers
                    GDF_filt = GDF[GDF.L2_OUTLIER.__eq__(True)].copy()
                    ax.scatter(GDF_filt['Lon'], GDF_filt['Lat'], c='r', s=150, label='SSIM',
                               **common_kw)
                if self.tieP_filter_level > 2:
                    # flag level 3 outliers
                    GDF_filt = GDF[GDF.L3_OUTLIER.__eq__(True)].copy()
                    ax.scatter(GDF_filt['Lon'], GDF_filt['Lat'], c='y', s=250, label='RANSAC',
                               **common_kw)

                if self.tieP_filter_level > 0:
                    ax.legend(loc=0, scatterpoints=1)

            # plot all points or vectors on top
            if not GDF.empty:
                vmin_auto, vmax_auto = \
                    (np.percentile(GDF[attribute2plot], 0),
                     np.percentile(GDF[attribute2plot], 98)) \
                    if attribute2plot != 'ANGLE' else (0, 360)
                vmin = vmin if vmin is not None else vmin_auto
                vmax = vmax if vmax is not None else vmax_auto

                if shapes2plot == 'vectors':
                    # plot shift vectors
                    # doc: https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.quiver.html
                    mappable = ax.quiver(
                        GDF['Lon'].values, GDF['Lat'].values,
                        -GDF['X_SHIFT_M'].values,
                        -GDF['Y_SHIFT_M'].values,  # invert absolute shifts to make arrows point to tgt
                        GDF[attribute2plot].clip(vmin, vmax),  # sets the colors
                        scale=1200 / vector_scale,  # larger values decrease the arrow length
                        width=.0015,  # arrow width (in relation to plot width)
                        # linewidth=1, # maybe use this to mark outliers instead of scatter points
                        cmap=palette,
                        pivot='middle',  # position the middle point of the arrows onto the tie point location
                        transform=PlateCarree()
                        )

                elif shapes2plot == 'points':
                    # plot tie points
                    mappable = ax.scatter(
                        GDF['Lon'], GDF['Lat'],
                        c=GDF[attribute2plot],
                        lw=0,
                        cmap=palette,
                        marker='o' if len(GDF) < 10000 else '.',
                        s=50,
                        alpha=1.0,
                        vmin=vmin,
                        vmax=vmax,
                        transform=PlateCarree())
                    pass
                else:
                    raise ValueError("The parameter 'shapes2plot' must be set to 'vectors' or 'points'. "
                                     "Received %s." % shapes2plot)

                # add colorbar
                divider = make_axes_locatable(ax)
                cax = divider.new_vertical(size="2%", pad=0.4, pack_start=True,
                                           axes_class=plt.Axes  # needed because ax is a GeoAxis instance
                                           )
                fig.add_axes(cax)
                fig.colorbar(mappable, cax=cax, orientation="horizontal")

                # hack to enlarge the figure on the top to avoid cutting off the title (everthing else has no effect)
                divider.new_vertical(size="2%", pad=0.4, pack_start=False, axes_class=plt.Axes)

            else:
                msg = "The map does not contain any tie points \n" \
                      "because all the found tie points were flagged as false-positives."
                ax.add_artist(AnchoredText(msg, loc='lower center', prop=dict(c='r')))

                if not self.q:
                    warnings.warn(msg)

        else:
            msg = "The map does not contain any tie points because no tie points were found at all."
            ax.add_artist(AnchoredText(msg, loc='lower center', prop=dict(c='r')))

            if not self.q:
                warnings.warn(msg)

        if savefigPath:
            fig.savefig(savefigPath, dpi=savefigDPI)

        if return_map:
            return fig, ax

        if showFig and not self.q:
            plt.show(block=True)
        else:
            plt.close(fig)

    def view_CoRegPoints_folium(self, attribute2plot: str = 'ABS_SHIFT'):
        warnings.warn(UserWarning('This function is still under construction and may not work as expected!'))
        assert self.CoRegPoints_table is not None, 'Calculate tie point grid first!'

        import folium
        import geojson
        from folium.raster_layers import ImageOverlay

        lon_min, lat_min, lon_max, lat_max = \
            reproject_shapelyGeometry(self.im2shift.box.mapPoly, self.im2shift.projection, 4326).bounds
        center_lon, center_lat = (lon_min + lon_max) / 2, (lat_min + lat_max) / 2

        # get image to plot
        image2plot = self.im2shift[:, :, 0]  # FIXME hardcoded band

        from py_tools_ds.geo.raster.reproject import warp_ndarray
        image2plot, gt, prj = \
            warp_ndarray(image2plot, self.im2shift.geotransform, self.im2shift.projection,
                         in_nodata=self.nodata[1], out_nodata=self.nodata[1], out_XYdims=(1000, 1000), q=True,
                         out_prj='epsg:3857')  # image must be transformed into web mercator projection

        # create map
        map_osm = folium.Map(location=[center_lat, center_lon])  # ,zoom_start=3)
        # import matplotlib
        ImageOverlay(
            colormap=lambda x: (1, 0, 0, x),  # TODO a colormap must be given
            # colormap=matplotlib.cm.gray, # does not work
            image=image2plot, bounds=[[lat_min, lon_min], [lat_max, lon_max]],
        ).add_to(map_osm)

        points_values = self.CoRegPoints_table[['geometry', attribute2plot]]
        points_values.geometry.crs = points_values.crs
        folium.GeoJson(points_values).add_to(map_osm)

        # add overlap polygon
        overlapPoly = reproject_shapelyGeometry(self.COREG_obj.overlap_poly, self.im2shift.prj, 4326)
        gjs = geojson.Feature(geometry=overlapPoly, properties={})
        folium.GeoJson(gjs).add_to(map_osm)

        return map_osm

    def _get_updated_map_info_meanShifts(self) -> list:
        """Return the updated map info of the target image, shifted on the basis of the mean X/Y shifts."""
        original_map_info = geotransform2mapinfo(self.im2shift.gt, self.im2shift.prj)
        updated_map_info = copy(original_map_info)
        updated_map_info[3] = str(float(original_map_info[3]) + self.tiepoint_grid.mean_x_shift_map)
        updated_map_info[4] = str(float(original_map_info[4]) + self.tiepoint_grid.mean_y_shift_map)
        return updated_map_info

    @property
    def coreg_info(self) -> dict:
        """Return a dictionary containing everthing to correct the detected local displacements of the target image."""
        if self._coreg_info:
            return self._coreg_info
        else:
            if not self._tiepoint_grid:
                self.calculate_spatial_shifts()

            TPG = self._tiepoint_grid

            self._coreg_info = {
                'GCPList': TPG.GCPList,
                'mean_shifts_px': {'x': TPG.mean_x_shift_px if TPG.GCPList else None,
                                   'y': TPG.mean_y_shift_px if TPG.GCPList else None},
                'mean_shifts_map': {'x': TPG.mean_x_shift_map if TPG.GCPList else None,
                                    'y': TPG.mean_y_shift_map if TPG.GCPList else None},
                'updated map info means': self._get_updated_map_info_meanShifts() if TPG.GCPList else None,
                'original map info': geotransform2mapinfo(self.imref.gt, self.imref.prj),
                'reference projection': self.imref.prj,
                'reference geotransform': self.imref.gt,
                'reference grid': [[self.imref.gt[0], self.imref.gt[0] + self.imref.gt[1]],
                                   [self.imref.gt[3], self.imref.gt[3] + self.imref.gt[5]]],
                'reference extent': {'cols': self.imref.xgsd, 'rows': self.imref.ygsd},  # FIXME not needed anymore
                'success': self.success
            }
            return self.coreg_info

    def correct_shifts(self,
                       max_GCP_count: int = None,
                       cliptoextent: bool = False,
                       min_points_local_corr: int = 5
                       ) -> OrderedDict:
        """Perform a local shift correction using all points from the previously calculated tie point grid.

        NOTE: Only valid matches are used as GCP points.

        :param max_GCP_count: maximum number of GCPs to use
        :param cliptoextent:  whether to clip the output image to its real extent
        :param min_points_local_corr:   number of valid tie points, below which a global shift correction is
                                        performed instead of a local correction (global X/Y shift is then computed as
                                        the mean shift of the remaining points)(default: 5 tie points)
        :return:
        """
        if not self._tiepoint_grid:
            self.calculate_spatial_shifts()

        if self.tiepoint_grid.GCPList:
            if max_GCP_count:
                self.coreg_info['GCPList'] = self.coreg_info['GCPList'][:max_GCP_count]

            # make sure the correction is applied to the original target image
            im2shift = GeoArray(self.params['im_tgt'], nodata=self.nodata[1], progress=self.progress, q=self.q)

            if has_metaRotation(im2shift):
                # resample the target image because (so far) the computed shifts cannot be applied to a dataset with
                # a metadata rotation (GDAL GeoTransform not 0 at positons 2 and 4)
                im2shift = remove_metaRotation(im2shift)

            # apply the correction
            DS = DESHIFTER(im2shift, self.coreg_info,
                           path_out=self.path_out,
                           fmt_out=self.fmt_out,
                           out_crea_options=self.out_creaOpt,
                           align_grids=self.align_grids,
                           match_gsd=self.match_gsd,
                           out_gsd=self.out_gsd,
                           target_xyGrid=self.target_xyGrid,
                           min_points_local_corr=min_points_local_corr,
                           resamp_alg=self.rspAlg_DS,
                           cliptoextent=cliptoextent,
                           # clipextent=self.im2shift.box.boxMapYX,
                           progress=self.progress,
                           v=self.v,
                           q=self.q)

            self.deshift_results = DS.correct_shifts()
            return self.deshift_results
        else:
            if not self.q:
                warnings.warn('Correction of geometric shifts failed because the input GCP list is empty!')
