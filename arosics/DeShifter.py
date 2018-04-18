# -*- coding: utf-8 -*-

import collections
import time
import warnings

# internal modules
from geoarray import GeoArray
from py_tools_ds.geo.map_info import mapinfo2geotransform, geotransform2mapinfo
from py_tools_ds.geo.coord_grid import is_coord_grid_equal
from py_tools_ds.geo.projection import prj_equal
from py_tools_ds.geo.raster.reproject import warp_ndarray
from py_tools_ds.numeric.vector import find_nearest

__author__ = 'Daniel Scheffler'

_dict_rspAlg_rsp_Int = {'nearest': 0, 'bilinear': 1, 'cubic': 2, 'cubic_spline': 3, 'lanczos': 4, 'average': 5,
                        'mode': 6, 'max': 7, 'min': 8, 'med': 9, 'q1': 10, 'q2': 11,
                        0: 'nearest', 1: 'bilinear', 2: 'cubic', 3: 'cubic_spline', 4: 'lanczos', 5: 'average',
                        6: 'mode', 7: 'max', 8: 'min', 9: 'med', 10: 'q1', 11: 'q2'}


class DESHIFTER(object):
    """See help(DESHIFTER) for documentation!"""

    def __init__(self, im2shift, coreg_results, **kwargs):
        """
        Deshift an image array or one of its products by applying the coregistration info calculated by COREG class.

        :param im2shift:            <path,GeoArray> path of an image to be de-shifted or alternatively a GeoArray object
        :param coreg_results:       <dict> the results of the co-registration as given by COREG.coreg_info or
                                    COREG_LOCAL.coreg_info respectively

        :Keyword Arguments:
            - path_out(str):        /output/directory/filename for coregistered results
            - fmt_out (str):        raster file format for output file. ignored if path_out is None. can be any GDAL
                                    compatible raster file format (e.g. 'ENVI', 'GeoTIFF'; default: ENVI)
            - out_crea_options(list): GDAL creation options for the output image,
                                    e.g. ["QUALITY=20", "REVERSIBLE=YES", "WRITE_METADATA=YES"]
            - band2process (int):   The index of the band to be processed within the given array (starts with 1),
                                    default = None (all bands are processed)
            - nodata(int, float):   no data value of the image to be de-shifted
            - out_gsd (float):      output pixel size in units of the reference coordinate system (default = pixel size
                                    of the input array), given values are overridden by match_gsd=True
            - align_grids (bool):   True: align the input coordinate grid to the reference (does not affect the
                                    output pixel size as long as input and output pixel sizes are compatible
                                    (5:30 or 10:30 but not 4:30), default = False
            - match_gsd (bool):     True: match the input pixel size to the reference pixel size,
                                    default = False
            - target_xyGrid(list):  a list with an x-grid and a y-grid like [[15,45], [15,45]].
                                    This overrides 'out_gsd', 'align_grids' and 'match_gsd'.
            - min_points_local_corr number of valid tie points, below which a global shift correction is performed
                                    instead of a local correction (global X/Y shift is then computed as the mean shift
                                    of the remaining points)(default: 5 tie points)
            - resamp_alg(str)       the resampling algorithm to be used if neccessary
                                    (valid algorithms: nearest, bilinear, cubic, cubic_spline, lanczos, average, mode,
                                                       max, min, med, q1, q3)
            - cliptoextent (bool):  True: clip the input image to its actual bounds while deleting possible no data
                                    areas outside of the actual bounds, default = False
            - clipextent (list):    xmin, ymin, xmax, ymax - if given the calculation of the actual bounds is skipped.
                                    The given coordinates are automatically snapped to the output grid.
            - CPUs(int):            number of CPUs to use (default: None, which means 'all CPUs available')
            - progress(bool):       show progress bars (default: True)
            - v(bool):              verbose mode (default: False)
            - q(bool):              quiet mode (default: False)

        """

        # private attributes
        self._grids_alignable = None

        # store args / kwargs
        self.init_args = dict([x for x in locals().items() if x[0] != "self" and not x[0].startswith('__')])
        self.init_kwargs = self.init_args['kwargs']

        # unpack args
        self.im2shift = im2shift if isinstance(im2shift, GeoArray) else GeoArray(im2shift)
        self.GCPList = coreg_results['GCPList'] if 'GCPList' in coreg_results else None
        self.ref_gt = coreg_results['reference geotransform']
        self.ref_grid = coreg_results['reference grid']
        self.ref_prj = coreg_results['reference projection']

        # unpack kwargs
        self.path_out = kwargs.get('path_out', None)
        self.fmt_out = kwargs.get('fmt_out', 'ENVI')
        self.out_creaOpt = kwargs.get('out_crea_options', [])
        self.band2process = kwargs.get('band2process', None)  # starts with 1 # FIXME why?
        self.band2process = \
            self.band2process - 1 if self.band2process is not None else None  # internally handled as band index
        self.nodata = kwargs.get('nodata', self.im2shift.nodata)
        self.align_grids = kwargs.get('align_grids', False)
        self.min_points_local_corr = kwargs.get('min_points_local_corr', 5)
        self.rspAlg = kwargs.get('resamp_alg', 'cubic')  # TODO accept also integers
        self.cliptoextent = kwargs.get('cliptoextent', False)
        self.clipextent = kwargs.get('clipextent', None)
        self.CPUs = kwargs.get('CPUs', None)
        self.v = kwargs.get('v', False)
        self.q = kwargs.get('q', False) if not self.v else False  # overridden by v
        self.progress = kwargs.get('progress', True) if not self.q else False  # overridden by q

        self.im2shift.nodata = kwargs.get('nodata', self.im2shift.nodata)
        self.im2shift.q = self.q
        self.shift_prj = self.im2shift.projection
        self.shift_gt = list(self.im2shift.geotransform)

        # in case of local shift correction and local coreg results contain less points than min_points_local_corr:
        # force global correction based on mean X/Y shifts
        if 'GCPList' in coreg_results and len(coreg_results['GCPList']) < self.min_points_local_corr:
            warnings.warn('Only %s valid tie point(s) could be identified. A local shift correction is therefore not '
                          'reasonable and could cause artifacts in the output image. The target image is '
                          'corrected globally with the mean X/Y shift of %.3f/%.3f pixels.'
                          % (len(self.GCPList), coreg_results['mean_shifts_px']['x'],
                             coreg_results['mean_shifts_px']['y']))
            self.GCPList = None
            coreg_results['updated map info'] = coreg_results['updated map info means']

        # in case of global shift correction -> the updated map info from coreg_results already has the final map info
        # BUT: this will updated in correct_shifts() if clipextent is given or warping is needed
        if not self.GCPList:
            mapI = coreg_results['updated map info']
            self.updated_map_info = mapI if mapI else geotransform2mapinfo(self.shift_gt, self.shift_prj)
            self.updated_gt = mapinfo2geotransform(self.updated_map_info) if mapI else self.shift_gt
            self.original_map_info = coreg_results['original map info']
        self.updated_projection = self.ref_prj

        self.out_grid = self._get_out_grid()  # needs self.ref_grid, self.im2shift
        self.out_gsd = [abs(self.out_grid[0][1] - self.out_grid[0][0]),
                        abs(self.out_grid[1][1] - self.out_grid[1][0])]  # xgsd, ygsd

        # assertions
        assert self.rspAlg in _dict_rspAlg_rsp_Int.keys(), \
            "'%s' is not a supported resampling algorithm." % self.rspAlg
        if self.band2process is not None:
            assert self.im2shift.bands - 1 >= self.band2process >= 0, \
                "The %s '%s' has %s %s. So 'band2process' must be %s%s. Got %s." \
                % (self.im2shift.__class__.__name__, self.im2shift.basename, self.im2shift.bands,
                   'bands' if self.im2shift.bands > 1 else 'band', 'between 1 and ' if self.im2shift.bands > 1 else '',
                   self.im2shift.bands, self.band2process + 1)

        # set defaults for general class attributes
        self.is_shifted = False  # this is not included in COREG.coreg_info
        self.is_resampled = False  # this is not included in COREG.coreg_info
        self.tracked_errors = []
        self.arr_shifted = None  # set by self.correct_shifts
        self.GeoArray_shifted = None  # set by self.correct_shifts

    def _get_out_grid(self):
        # parse given params
        out_gsd = self.init_kwargs.get('out_gsd', None)
        match_gsd = self.init_kwargs.get('match_gsd', False)
        out_grid = self.init_kwargs.get('target_xyGrid', None)

        # assertions
        assert out_grid is None or (isinstance(out_grid, (list, tuple)) and len(out_grid) == 2)
        assert out_gsd is None or (isinstance(out_gsd, (int, tuple, list)) and len(out_gsd) == 2)

        ref_xgsd, ref_ygsd = (self.ref_grid[0][1] - self.ref_grid[0][0], abs(self.ref_grid[1][1] - self.ref_grid[1][0]))

        def get_grid(gt, xgsd, ygsd): return [[gt[0], gt[0] + xgsd], [gt[3], gt[3] - ygsd]]

        # get out_grid
        if out_grid:
            # output grid is given
            pass

        elif out_gsd:
            out_xgsd, out_ygsd = [out_gsd, out_gsd] if isinstance(out_gsd, int) else out_gsd

            if match_gsd and (out_xgsd, out_ygsd) != (ref_xgsd, ref_ygsd):
                warnings.warn("\nThe parameter 'match_gsd is ignored because another output ground sampling distance "
                              "was explicitly given.")
            if self.align_grids and \
               self._are_grids_alignable(self.im2shift.xgsd, self.im2shift.ygsd, out_xgsd, out_ygsd):
                # use grid of reference image with the given output gsd
                out_grid = get_grid(self.ref_gt, out_xgsd, out_ygsd)
            else:  # no grid alignment
                # use grid of input image with the given output gsd
                out_grid = get_grid(self.im2shift.geotransform, out_xgsd, out_ygsd)

        elif match_gsd:
            if self.align_grids:
                # use reference grid
                out_grid = self.ref_grid
            else:
                # use grid of input image and reference gsd
                out_grid = get_grid(self.im2shift.geotransform, ref_xgsd, ref_ygsd)

        else:
            if self.align_grids and \
               self._are_grids_alignable(self.im2shift.xgsd, self.im2shift.ygsd, ref_xgsd, ref_ygsd):
                # use origin of reference image and gsd of input image
                out_grid = get_grid(self.ref_gt, self.im2shift.xgsd, self.im2shift.ygsd)
            else:
                # use input image grid
                out_grid = get_grid(self.im2shift.geotransform, self.im2shift.xgsd, self.im2shift.ygsd)

        return out_grid

    @property
    def warping_needed(self):
        """Returns True if image warping is needed in consideration of the input parameters of DESHIFTER."""

        assert self.out_grid, 'Output grid must be calculated before.'
        equal_prj = prj_equal(self.ref_prj, self.shift_prj)
        return \
            False if (equal_prj and not self.GCPList and is_coord_grid_equal(self.updated_gt, *self.out_grid)) else True

    def _are_grids_alignable(self, in_xgsd, in_ygsd, out_xgsd, out_ygsd):
        """Checks if the input image pixel grid is alignable to the output grid.

        :param in_xgsd:
        :param in_ygsd:
        :param out_xgsd:
        :param out_ygsd:
        :return:
        """
        if self._grids_alignable is None:
            def is_alignable(gsd1, gsd2):
                """Check if pixel sizes are divisible."""
                return max(gsd1, gsd2) % min(gsd1, gsd2) == 0

            self._grids_alignable = \
                False if (not is_alignable(in_xgsd, out_xgsd) or not is_alignable(in_ygsd, out_ygsd)) else True

            if self._grids_alignable is False and not self.q:
                warnings.warn("\nThe coordinate grid of %s cannot be aligned to the desired grid because their pixel "
                              "sizes are not exact multiples of each other (input [X/Y]: %s/%s; desired [X/Y]: %s/%s). "
                              "Therefore the original grid is chosen for the resampled output image. If you don´t like "
                              "that you can use the 'out_gsd' or 'match_gsd' parameters to set an appropriate output "
                              "pixel size or to allow changing the pixel size.\n"
                              % (self.im2shift.basename, in_xgsd, in_ygsd, out_xgsd, out_ygsd))

        return self._grids_alignable

    def _get_out_extent(self):
        if self.clipextent is None:
            # no clip extent has been given
            if self.cliptoextent:
                # use actual image corners as clip extent
                self.clipextent = self.im2shift.footprint_poly.envelope.bounds
            else:
                # use outer bounds of the image as clip extent
                xmin, xmax, ymin, ymax = self.im2shift.box.boundsMap
                self.clipextent = xmin, ymin, xmax, ymax

        # snap clipextent to output grid
        # (in case of odd input coords the output coords are moved INSIDE the input array)
        xmin, ymin, xmax, ymax = self.clipextent
        xmin = find_nearest(self.out_grid[0], xmin, roundAlg='on', extrapolate=True)
        ymin = find_nearest(self.out_grid[1], ymin, roundAlg='on', extrapolate=True)
        xmax = find_nearest(self.out_grid[0], xmax, roundAlg='off', extrapolate=True)
        ymax = find_nearest(self.out_grid[1], ymax, roundAlg='off', extrapolate=True)
        return xmin, ymin, xmax, ymax

    def correct_shifts(self):
        # type: () -> collections.OrderedDict

        if not self.q:
            print('Correcting geometric shifts...')

        t_start = time.time()

        if not self.warping_needed:
            """NO RESAMPLING NEEDED"""

            self.is_shifted = True
            self.is_resampled = False
            xmin, ymin, xmax, ymax = self._get_out_extent()

            if self.cliptoextent:
                # TODO validate results
                # TODO -> output extent does not seem to be the requested one! (only relevant if align_grids=False)
                # get shifted array
                shifted_geoArr = GeoArray(self.im2shift[:], tuple(self.updated_gt), self.shift_prj)

                # clip with target extent
                #  NOTE: get_mapPos() does not perform any resampling as long as source and target projection are equal
                self.arr_shifted, self.updated_gt, self.updated_projection = \
                    shifted_geoArr.get_mapPos((xmin, ymin, xmax, ymax),
                                              self.shift_prj,
                                              fillVal=self.nodata,
                                              band2get=self.band2process)

                self.updated_map_info = geotransform2mapinfo(self.updated_gt, self.updated_projection)

            else:
                # array keeps the same; updated gt and prj are taken from coreg_info
                self.arr_shifted = self.im2shift[:, :, self.band2process] \
                    if self.band2process is not None else self.im2shift[:]

            out_geoArr = GeoArray(self.arr_shifted, self.updated_gt, self.updated_projection, q=self.q)
            out_geoArr.nodata = self.nodata  # equals self.im2shift.nodata after __init__()
            out_geoArr.metadata = self.im2shift.metadata[[self.band2process]] \
                if self.band2process is not None else self.im2shift.metadata

            self.GeoArray_shifted = out_geoArr

            if self.path_out:
                out_geoArr.save(self.path_out, fmt=self.fmt_out)

        else:  # FIXME equal_prj==False ist noch NICHT implementiert
            """RESAMPLING NEEDED"""
            # FIXME avoid reading the whole band if clip_extent is passed

            in_arr = self.im2shift[:, :, self.band2process] \
                if self.band2process is not None and self.im2shift.ndim == 3 else self.im2shift[:]

            if not self.GCPList:
                # apply XY-shifts to input image gt 'shift_gt' in order to correct the shifts before warping
                self.shift_gt[0], self.shift_gt[3] = self.updated_gt[0], self.updated_gt[3]

            # get resampled array
            out_arr, out_gt, out_prj = \
                warp_ndarray(in_arr, self.shift_gt, self.shift_prj, self.ref_prj,
                             rspAlg=_dict_rspAlg_rsp_Int[self.rspAlg],
                             in_nodata=self.nodata,
                             out_nodata=self.nodata,
                             out_gsd=self.out_gsd,
                             out_bounds=self._get_out_extent(),  # always returns an extent snapped to the target grid
                             gcpList=self.GCPList,
                             # polynomialOrder=str(3),
                             # options='-refine_gcps 500 1.9',
                             # warpOptions=['-refine_gcps 500 1.9'],
                             # options='-wm 10000',# -order 3',
                             # options=['-order 3'],
                             # options=['GDAL_CACHEMAX 800 '],
                             # warpMemoryLimit=125829120, # 120MB
                             CPUs=self.CPUs,
                             progress=self.progress,
                             q=self.q)

            out_geoArr = GeoArray(out_arr, out_gt, out_prj, q=self.q)
            out_geoArr.nodata = self.nodata  # equals self.im2shift.nodata after __init__()
            out_geoArr.metadata = self.im2shift.metadata[[self.band2process]] \
                if self.band2process is not None else self.im2shift.metadata

            self.arr_shifted = out_arr
            self.updated_gt = out_gt
            self.updated_projection = out_prj
            self.updated_map_info = geotransform2mapinfo(out_gt, out_prj)
            self.GeoArray_shifted = out_geoArr
            self.is_shifted = True
            self.is_resampled = True

            if self.path_out:
                out_geoArr.save(self.path_out, fmt=self.fmt_out, creationOptions=self.out_creaOpt)

        # validation
        if not is_coord_grid_equal(self.updated_gt, *self.out_grid, tolerance=1.e8):
            raise RuntimeError('DESHIFTER output dataset has not the desired target pixel grid. Target grid '
                               'was %s. Output geotransform is %s.' % (str(self.out_grid), str(self.updated_gt)))
        # TODO to be continued (extent, map info, ...)

        if self.v:
            print('Time for shift correction: %.2fs' % (time.time() - t_start))
        return self.deshift_results

    @property
    def deshift_results(self):
        deshift_results = collections.OrderedDict()
        deshift_results.update({
            'band': self.band2process,
            'is shifted': self.is_shifted,
            'is resampled': self.is_resampled,
            'updated map info': self.updated_map_info,
            'updated geotransform': self.updated_gt,
            'updated projection': self.updated_projection,
            'arr_shifted': self.arr_shifted,
            'GeoArray_shifted': self.GeoArray_shifted
        })
        return deshift_results


def deshift_image_using_coreg_info(im2shift, coreg_results, path_out=None, fmt_out='ENVI', q=False):
    """Corrects a geometrically distorted image using previously calculated coregistration info. This function can be
    used for example to corrects spatial shifts of mask files using the same transformation parameters that have been
    used to correct their source images.

    :param im2shift:      <path,GeoArray> path of an image to be de-shifted or alternatively a GeoArray object
    :param coreg_results: <dict> the results of the co-registration as given by COREG.coreg_info or
                          COREG_LOCAL.coreg_info respectively
    :param path_out:      /output/directory/filename for coregistered results. If None, no output is written - only
                          the shift corrected results are returned.
    :param fmt_out:       raster file format for output file. ignored if path_out is None. can be any GDAL
                                        compatible raster file format (e.g. 'ENVI', 'GeoTIFF'; default: ENVI)
    :param q:             quiet mode (default: False)
    :return:
    """
    deshift_results = DESHIFTER(im2shift, coreg_results).correct_shifts()

    if path_out:
        deshift_results['GeoArray_shifted'].save(path_out, fmt_out=fmt_out, q=q)

    return deshift_results
