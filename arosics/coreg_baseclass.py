# -*- coding: utf-8 -*-
__author__='Daniel Scheffler'


import warnings
import os

# custom
try:
    import gdal
except ImportError:
    from osgeo import gdal
try:
    import pyfftw
except ImportError:
    pyfftw = None
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .Tie_Point_Grid import Tie_Point_Grid
from .CoReg import COREG
from .DeShifter import DESHIFTER, _dict_rspAlg_rsp_Int
from py_tools_ds.geo.coord_trafo import transform_any_prj, reproject_shapelyGeometry
from geoarray import GeoArray



class CoReg(object):
    def __init__(self, im_ref, im_tgt, path_out=None, fmt_out='ENVI', out_crea_options=None, r_b4match=1, s_b4match=1,
                 win_size=(256, 256), max_iter=5, max_shift=5, align_grids=True, match_gsd=False, out_gsd=None,
                 target_xyGrid=None,
                 resamp_alg_deshift='cubic', resamp_alg_calc='cubic', footprint_poly_ref=None, footprint_poly_tgt=None,
                 data_corners_ref=None, data_corners_tgt=None, nodata=(None, None), calc_corners=True,
                 binary_ws=True, force_quadratic_win=True, mask_baddata_ref=None, mask_baddata_tgt=None, CPUs=None,
                 progress=True, v=False, path_verbose_out=None, q=False, ignore_errors=True):

        """CoReg Baseclass
        
        
        :param im_ref(str, GeoArray):   source path of reference image (any GDAL compatible image format is supported)
        :param im_tgt(str, GeoArray):   source path of image to be shifted (any GDAL compatible image format is supported)
        :param path_out(str):           target path of the coregistered image
                                            - if None (default), the method correct_shifts() does not write to disk
                                            - if 'auto': /dir/of/im1/<im1>__shifted_to__<im0>.bsq
        :param fmt_out(str):            raster file format for output file. ignored if path_out is None. can be any GDAL
                                        compatible raster file format (e.g. 'ENVI', 'GeoTIFF'; default: ENVI). Refer to
                                        http://www.gdal.org/formats_list.html to get a full list of supported formats.
        :param out_crea_options(list):  GDAL creation options for the output image,
                                        e.g. ["QUALITY=80", "REVERSIBLE=YES", "WRITE_METADATA=YES"]
        :param r_b4match(int):          band of reference image to be used for matching (starts with 1; default: 1)
        :param s_b4match(int):          band of shift image to be used for matching (starts with 1; default: 1)
        :param win_size(tuple):         custom matching window size [pixels] (default: (256,256))
        :param max_iter(int):           maximum number of iterations for matching (default: 5)
        :param max_shift(int):          maximum shift distance in reference image pixel units (default: 5 px)
        :param align_grids(bool):       align the coordinate grids of the image to be and the reference image (default: 0)
        :param match_gsd(bool):         match the output pixel size to pixel size of the reference image (default: 0)
        :param out_gsd(tuple):          xgsd ygsd: set the output pixel size in map units
                                        (default: original pixel size of the image to be shifted)
        :param target_xyGrid(list):     a list with a target x-grid and a target y-grid like [[15,45], [15,45]]
                                        This overrides 'out_gsd', 'align_grids' and 'match_gsd'.
        :param resamp_alg_deshift(str)  the resampling algorithm to be used for shift correction (if neccessary)
                                        valid algorithms: nearest, bilinear, cubic, cubic_spline, lanczos, average, mode,
                                                          max, min, med, q1, q3
                                        default: cubic
        :param resamp_alg_calc(str)     the resampling algorithm to be used for all warping processes during calculation
                                        of spatial shifts
                                        (valid algorithms: nearest, bilinear, cubic, cubic_spline, lanczos, average, mode,
                                                       max, min, med, q1, q3)
                                        default: cubic (highly recommended)
        :param footprint_poly_ref(str): footprint polygon of the reference image (WKT string or shapely.geometry.Polygon),
                                        e.g. 'POLYGON ((299999 6000000, 299999 5890200, 409799 5890200, 409799 6000000,
                                                        299999 6000000))'
        :param footprint_poly_tgt(str): footprint polygon of the image to be shifted (WKT string or shapely.geometry.Polygon)
                                        e.g. 'POLYGON ((299999 6000000, 299999 5890200, 409799 5890200, 409799 6000000,
                                                        299999 6000000))'
        :param data_corners_ref(list):  map coordinates of data corners within reference image.
                                        ignored if footprint_poly_ref is given.
        :param data_corners_tgt(list):  map coordinates of data corners within image to be shifted.
                                        ignored if footprint_poly_tgt is given.
        :param nodata(tuple):           no data values for reference image and image to be shifted
        :param calc_corners(bool):      calculate true positions of the dataset corners in order to get a useful
                                        matching window position within the actual image overlap
                                        (default: 1; deactivated if '-cor0' and '-cor1' are given
        :param binary_ws(bool):         use binary X/Y dimensions for the matching window (default: 1)
        :param mask_baddata_ref(str, GeoArray): path to a 2D boolean mask file (or an instance of GeoArray) for the
                                                reference image where all bad data pixels (e.g. clouds) are marked with
                                                True and the remaining pixels with False. Must have the same geographic
                                                extent and projection like 'im_ref'. The mask is used to check if the
                                                chosen matching window position is valid in the sense of useful data.
                                                Otherwise this window position is rejected.
        :param mask_baddata_tgt(str, GeoArray): path to a 2D boolean mask file (or an instance of GeoArray) for the
                                                image to be shifted where all bad data pixels (e.g. clouds) are marked
                                                with True and the remaining pixels with False. Must have the same
                                                geographic extent and projection like 'im_ref'. The mask is used to
                                                check if the chosen matching window position is valid in the sense of
                                                useful data. Otherwise this window position is rejected.
        :param CPUs(int):               number of CPUs to use during pixel grid equalization
                                        (default: None, which means 'all CPUs available')
        :param force_quadratic_win(bool):   force a quadratic matching window (default: 1)
        :param progress(bool):          show progress bars (default: True)
        :param v(bool):                 verbose mode (default: False)
        :param path_verbose_out(str):   an optional output directory for intermediate results
                                        (if not given, no intermediate results are written to disk)
        :param q(bool):                 quiet mode (default: False)
        :param ignore_errors(bool):     Useful for batch processing. (default: False)
                                        In case of error COREG.success == False and COREG.x_shift_px/COREG.y_shift_px
                                        is None
        """

        self.params = dict([x for x in locals().items() if x[0] != "self"])

        # assertions
        assert gdal.GetDriverByName(fmt_out), "'%s' is not a supported GDAL driver." % fmt_out
        if match_gsd and out_gsd: warnings.warn("'-out_gsd' is ignored because '-match_gsd' is set.\n")
        if out_gsd:  assert isinstance(out_gsd, list) and len(out_gsd) == 2, 'out_gsd must be a list with two values.'
        if data_corners_ref and not isinstance(data_corners_ref[0], list):  # group if not [[x,y],[x,y]..] but [x,y,x,y,]
            data_corners_ref = [data_corners_ref[i:i + 2] for i in range(0, len(data_corners_ref), 2)]
        if data_corners_tgt and not isinstance(data_corners_tgt[0], list):  # group if not [[x,y],[x,y]..]
            data_corners_tgt = [data_corners_tgt[i:i + 2] for i in range(0, len(data_corners_tgt), 2)]
        if nodata: assert isinstance(nodata, tuple) and len(nodata) == 2, "'nodata' must be a tuple with two values." \
                                                                          "Got %s with length %s." % (
                                                                          type(nodata), len(nodata))
        for rspAlg in [resamp_alg_deshift, resamp_alg_calc]:
            assert rspAlg in _dict_rspAlg_rsp_Int.keys(), "'%s' is not a supported resampling algorithm." % rspAlg
        if resamp_alg_calc == 'average' and (v or not q):
            warnings.warn("The resampling algorithm 'average' causes sinus-shaped patterns in fft images that will "
                          "affect the precision of the calculated spatial shifts! It is highly recommended to "
                          "choose another resampling algorithm.")


        self.path_out = path_out  # updated by self.set_outpathes
        self.fmt_out = fmt_out
        self.out_creaOpt = out_crea_options
        self.win_size_XY         = win_size                  # updated by self.get_opt_winpos_winsize()
        self.max_iter            = max_iter
        self.max_shift           = max_shift
        self.align_grids         = align_grids
        self.match_gsd           = match_gsd
        self.out_gsd             = out_gsd
        self.target_xyGrid       = target_xyGrid
        self.rspAlg_DS           = resamp_alg_deshift
        self.rspAlg_calc         = resamp_alg_calc
        self.calc_corners        = calc_corners
        self.CPUs                = CPUs
        self.bin_ws              = binary_ws
        self.force_quadratic_win = force_quadratic_win
        self.v                   = v
        self.path_verbose_out    = path_verbose_out
        self.q                   = q if not v else False # overridden by v
        self.progress            = progress if not q else False  # overridden by q

        self.ignErr              = ignore_errors
        self.max_win_sz_changes  = 3                   # TODO: änderung der window size, falls nach max_iter kein valider match gefunden