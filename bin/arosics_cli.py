# -*- coding: utf-8 -*-
# unicode_literals cause GDAL not to work properly
from __future__ import (division, print_function, absolute_import, unicode_literals)

import time
import sys
import warnings
import argparse

from arosics import COREG, COREG_LOCAL, __version__

__author__ = "Daniel Scheffler"


# sub-command functions
def run_global_coreg(args):
    COREG_obj = COREG(args.path_ref,
                      args.path_tgt,
                      path_out=args.path_out,
                      fmt_out=args.fmt_out,
                      r_b4match=args.br,
                      s_b4match=args.bs,
                      wp=args.wp,
                      ws=args.ws,
                      max_iter=args.max_iter,
                      max_shift=args.max_shift,
                      align_grids=args.align_grids,
                      match_gsd=args.match_gsd,
                      out_gsd=args.out_gsd,
                      resamp_alg_calc=args.rsp_alg_deshift,
                      resamp_alg_deshift=args.rsp_alg_calc,
                      data_corners_ref=args.cor0,
                      data_corners_tgt=args.cor1,
                      nodata=args.nodata,
                      calc_corners=args.calc_cor,
                      CPUs=None if args.mp else 1,
                      force_quadratic_win=args.quadratic_win,
                      binary_ws=args.bin_ws,
                      mask_baddata_ref=args.mask_ref,
                      mask_baddata_tgt=args.mask_tgt,
                      progress=args.progress,
                      v=args.v,
                      path_verbose_out=args.vo,
                      q=args.q,
                      ignore_errors=args.ignore_errors)
    COREG_obj.correct_shifts()


# sub-command functions
def run_local_coreg(args):
    CRL = COREG_LOCAL(args.path_ref,
                      args.path_tgt,
                      path_out=args.path_out,
                      fmt_out=args.fmt_out,
                      grid_res=args.grid_res,
                      max_points=args.max_points,
                      r_b4match=args.br,
                      s_b4match=args.bs,
                      window_size=args.ws,
                      max_iter=args.max_iter,
                      max_shift=args.max_shift,
                      tieP_filter_level=args.tieP_filter_level,
                      min_reliability=args.min_reliability,
                      rs_max_outlier=args.rs_max_outlier,
                      rs_tolerance=args.rs_tolerance,
                      # align_grids=args.align_grids,
                      # match_gsd=args.match_gsd,
                      # out_gsd=args.out_gsd,
                      resamp_alg_calc=args.rsp_alg_deshift,
                      resamp_alg_deshift=args.rsp_alg_calc,
                      data_corners_ref=args.cor0,
                      data_corners_tgt=args.cor1,
                      nodata=args.nodata,
                      calc_corners=args.calc_cor,
                      mask_baddata_ref=args.mask_ref,
                      mask_baddata_tgt=args.mask_tgt,
                      CPUs=None if args.mp else 1,
                      force_quadratic_win=args.quadratic_win,
                      binary_ws=args.bin_ws,
                      progress=args.progress,
                      v=args.v,
                      q=args.q,
                      )
    CRL.correct_shifts()


def get_arosics_argparser():
    """Return argument parser for arosics_cli.py program."""
    parser = argparse.ArgumentParser(
        prog='arosics_cli.py',

        description='Perform automatic subpixel co-registration of two satellite image datasets based on an image '
                    'matching approach working in the frequency domain, combined with a multistage workflow for '
                    'effective detection of false-positives. Python implementation by Daniel Scheffler '
                    '(daniel.scheffler [at] gfz-potsdam [dot] de). The scientific background is described in the paper '
                    'Scheffler D, Hollstein A, Diedrich H, Segl K, Hostert P. AROSICS: An Automated and Robust '
                    'Open-Source Image Co-Registration Software for Multi-Sensor Satellite Data. Remote Sensing. 2017;'
                    ' 9(7):676." (http://www.mdpi.com/2072-4292/9/7/676)',

        epilog="DETAILED DESCRIPTION: AROSICS detects and corrects global as well as local misregistrations between "
               "two input images in the subpixel scale, that are often present in satellite imagery. The input images "
               "can have any GDAL compatible image format (http://www.gdal.org/formats_list.html). Both of them must "
               "be approximately geocoded. In case of ENVI files, this means they must have a 'map info' and a "
               "'coordinate system string' as attributes of their header file. The input images must have a geographic "
               "overlap but clipping them to same geographical extent is NOT neccessary. Please do not perform any "
               "spatial resampling of the input images before applying this algorithm. Any needed resampling of the "
               "data is done automatically. Thus, the input images may have different spatial resolutions. The current "
               "algorithm will not perform any ortho-rectification. So please use ortho-rectified input data in order "
               "to minimize local shifts in the input images. AROSICS supports local and global co-registration. LOCAL "
               "CO-REGISTRATION: A dense grid of tie points is automatically computed, whereas tie points are "
               "subsequently validated using a multistage workflow. Only those tie points not marked as "
               "false-positives are used to compute the parameters of an affine transformation. Warping of the target "
               "image is done using an appropriate resampling technique (cubic by default). GLOBAL CO-REGISTRATION: "
               "Only a global X/Y translation is computed within a small subset of the input images (window position "
               "is adjustable). This allows very fast co-registration but only corrects for translational (global) X/Y "
               "shifts. The calculated subpixel-shifts are (by default) applied to the geocoding information of the "
               "output image. No spatial resampling is done automatically as long as both input images have the same "
               "projection. If you need the output image to be aligned to the reference image coordinate grid (by "
               "using an appropriate resampling algorithm), use the '-align_grids' option. AROSICS is designed to "
               "robustly handle the typical difficulties of multi-sensoral/multi-temporal images. Clouds are "
               "automatically handled by the implemented outlier detection algorithms. The user may provide "
               "user-defined masks to exclude certain image areas from tie point creation. The image overlap area is "
               "automatically calculated. Thereby, no-data regions within the images are automatically respected. "
               "Providing the map coordinates of the actual data corners lets you save some calculation time, because "
               "in this case the automatic algorithm can be skipped. The no-data value of each image is automatically "
               "derived from the image corners. The verbose program mode gives some more output about the interim "
               "results, shows some figures and writes the used footprint and overlap polygons to disk. Note, that "
               "maybe the figures must be manually closed in in order to continue the processing (depending on your "
               "Python configuration). For further details regarding the implemented algorithm, example use cases, "
               "quality assessment and benchmarks refer to the above mentioned paper (Scheffler et al. 2017).")

    parser.add_argument('--version', action='version', version=__version__)

    #####################
    # GENERAL ARGUMENTS #
    #####################

    general_opts_parser = argparse.ArgumentParser(add_help=False)
    gop_p = general_opts_parser.add_argument
    gop_p('path_ref', type=str, help='source path of reference image (any GDAL compatible image format is supported)')

    gop_p('path_tgt', type=str,
          help='source path of image to be shifted (any GDAL compatible image format is supported)')

    gop_p('-o', nargs='?', dest='path_out', type=str, default='auto',
          help="target path of the coregistered image If 'auto' (default: /dir/of/im1/<im1>__shifted_to__<im0>.bsq)")

    gop_p('-fmt_out', nargs='?', type=str, default='ENVI',
          help="raster file format for output file. ignored if path_out is None. can "
               "be any GDAL compatible raster file format (e.g. 'ENVI', 'GTIFF'; default: ENVI)")

    gop_p('-br', nargs='?', type=int, default=1,
          help='band of reference image to be used for matching (starts with 1; default: 1)')

    gop_p('-bs', nargs='?', type=int, default=1,
          help='band of shift image to be used for matching (starts with 1; default: 1)')

    gop_p('-ws', nargs=2, metavar=('X size', 'Y size'), type=int, default=(256, 256),
          help="custom matching window size [pixels] (default: (256,256))")

    gop_p('-max_iter', nargs='?', type=int, default=5, help="maximum number of iterations for matching (default: 5)")

    gop_p('-max_shift', nargs='?', type=int, default=5,
          help="maximum shift distance in reference image pixel units (default: 5 px)")

    gop_p('-rsp_alg_deshift', nargs='?', type=int, choices=list(range(12)), default=2,
          help="the resampling algorithm to be used for shift correction (if neccessary) "
               "(valid algorithms: 0=nearest neighbour, 1=bilinear, 2=cubic, 3=cubic_spline, 4=lanczos, 5=average, "
               "6=mode, 7=max, 8=min, 9=med, 10=q1, 11=q3), default: 2")

    gop_p('-rsp_alg_calc', nargs='?', type=int, choices=list(range(12)), default=2,
          help="the resampling algorithm to be used for all warping processes during calculation of spatial shifts "
               "(valid algorithms: 0=nearest neighbour, 1=bilinear, 2=cubic, 3=cubic_spline, 4=lanczos, 5=average, "
               "6=mode, 7=max, 8=min, 9=med, 10=q1, 11=q3), default: 2 (highly recommended)")

    gop_p('-cor0', nargs=8, type=float, help="map coordinates of data corners within reference image: ",
          metavar=tuple("UL-X UL-Y UR-X UR-Y LR-X LR-Y LL-X LL-Y".split(' ')), default=None)

    gop_p('-cor1', nargs=8, type=float, help="map coordinates of data corners within image to be shifted: ",
          metavar=tuple("UL-X UL-Y UR-X UR-Y LR-X LR-Y LL-X LL-Y".split(' ')), default=None)

    gop_p('-calc_cor', nargs='?', type=int, choices=[0, 1], default=1,
          help="calculate true positions of the dataset corners in order to get a useful matching window position "
               "within the actual image overlap (default: 1; deactivated if '-cor0' and '-cor1' are given")

    gop_p('-nodata', nargs=2, type=float, metavar=('im0', 'im1'),
          help='no data values for reference image and image to be shifted', default=(None, None))

    gop_p('-bin_ws', nargs='?', type=int,
          help='use binary X/Y dimensions for the matching window (default: 1)', choices=[0, 1], default=1)

    gop_p('-quadratic_win', nargs='?', type=int,
          help='force a quadratic matching window (default: 1)', choices=[0, 1], default=1)

    gop_p('-mask_ref', nargs='?', type=str, metavar='file path', default=None,
          help="path to a 2D boolean mask file for the reference image where all bad data pixels (e.g. clouds) are "
               "marked with True or 1 and the remaining pixels with False or 0. Must have the same geographic extent "
               "and projection like the reference image. The mask is used to check if the chosen matching window "
               "position is valid in the sense of useful data. Otherwise this window position is rejected.")

    gop_p('-mask_tgt', nargs='?', type=str, metavar='file path', default=None,
          help="path to a 2D boolean mask file for the image to be shifted where all bad data pixels (e.g. clouds) are "
               "marked with True or 1 and the remaining pixels with False or 0. Must have the same geographic extent "
               "and projection like the the image to be shifted. The mask is used to check if the chosen matching "
               "window position is valid in the sense of useful data. Otherwise this window position is rejected.")

    gop_p('-mp', nargs='?', type=int, help='enable multiprocessing (default: 1)', choices=[0, 1], default=1)

    gop_p('-progress', nargs='?', type=int, help='show progress bars (default: 1)', default=1, choices=[0, 1])

    gop_p('-v', nargs='?', type=int, help='verbose mode (default: 0)', choices=[0, 1], default=0)

    gop_p('-q', nargs='?', type=int, help='quiet mode (default: 0)', choices=[0, 1], default=0)

    gop_p('-ignore_errors', nargs='?', type=int, choices=[0, 1], default=0,
          help='Useful for batch processing. (default: 0) In case of error '
               'COREG(_LOCAL).success == False and COREG(_LOCAL).x_shift_px/COREG(_LOCAL).y_shift_px is None')

    # TODO implement footprint_poly_ref, footprint_poly_tgt

    ##############
    # SUBPARSERS #
    ##############

    subparsers = parser.add_subparsers()

    # TODO add option to apply coreg results to multiple files
    #######################
    # SUBPARSER FOR COREG #
    #######################

    parse_coreg_global = subparsers.add_parser(
        'global', parents=[general_opts_parser], formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Detects and corrects global X/Y shifts between a target and refernce image. Geometric shifts are '
                    'calculated at a specific (adjustable) image position. Correction performs a global shifting in '
                    'X- or Y direction.',
        help="detect and correct global X/Y shifts (sub argument parser) - "
             "use '>>> python /path/to/arosics/bin/arosics_cli.py global -h' for documentation and usage hints")

    gloArg = parse_coreg_global.add_argument

    gloArg('-wp', nargs=2, metavar=('X', 'Y'), type=float,
           help="custom matching window position as map values in the same projection like the reference image "
                "(default: central position of image overlap)", default=(None, None))

    gloArg('-align_grids', nargs='?', type=int, choices=[0, 1],
           help='align the coordinate grids of the output image to the reference image (default: 0)', default=0)

    gloArg('-match_gsd', nargs='?', type=int, choices=[0, 1],
           help='match the output pixel size to the pixel size of the reference image (default: 0)', default=0)

    gloArg('-out_gsd', nargs=2, type=float, metavar=('xgsd', 'ygsd'),
           help='xgsd ygsd: set the output pixel size in map units (default: original pixel size of the image to be '
                'shifted)')

    gloArg('-vo', nargs='?', type=str, choices=[0, 1], help='an optional output directory for outputs of verbose mode'
                                                            '(if not given, no outputs are written to disk)',
           default=0, )

    parse_coreg_global.set_defaults(func=run_global_coreg)

    #############################
    # SUBPARSER FOR COREG LOCAL #
    #############################

    parse_coreg_local = subparsers.add_parser(
        'local', parents=[general_opts_parser], formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Applies the algorithm to detect spatial shifts to the whole overlap area of the input images. '
                    'Spatial shifts are calculated for each point in grid of which the parameters can be adjusted '
                    'using keyword arguments. Shift correction performs a polynomial transformation using the '
                    'calculated shifts of each point in the grid as GCPs. Thus this class can be used to correct '
                    'for locally varying geometric distortions of the target image.',
        help="detect and correct local shifts (sub argument parser)"
             "use '>>> python /path/to/arosics/bin/arosics_cli.py local -h' for documentation and usage hints")

    locArg = parse_coreg_local.add_argument

    locArg('grid_res', type=int, help='tie point grid resolution in pixels of the target image')

    locArg('-max_points', nargs='?', type=int,
           help="maximum number of points used to find coregistration tie points. NOTE: Points are selected randomly "
                "from the given point grid (specified by 'grid_res'). If the point does not provide enough points, all "
                "available points are chosen.")

    locArg('-projectDir', nargs='?', type=str, help=None, default=None)

    locArg('-tieP_filter_level', nargs='?', type=int, default=3, choices=[0, 1, 2, 3],
           help="filter tie points used for shift correction in different levels (default: 3). NOTE: lower levels are "
                "also included if a higher level is chosen. Level 0: no tie point filtering; Level 1: Reliablity "
                "filtering - filter all tie points out that have a low reliability according to internal tests; "
                "Level 2: SSIM filtering - filters all tie points out where shift correction does not increase image "
                "similarity within matching window (measured by mean structural similarity index) "
                "Level 3: RANSAC outlier detection")

    locArg('-min_reliability', nargs='?', type=float, default=60,
           help="Tie point filtering: minimum reliability threshold, below which tie points are marked as "
                "false-positives (default: 60 percent) - accepts values between 0 (no reliability) and 100 (perfect "
                "reliability) HINT: decrease this value in case of poor signal-to-noise ratio of your input data")

    locArg('-rs_max_outlier', nargs='?', type=float, default=10,
           help="RANSAC tie point filtering: proportion of expected outliers (default: 10 percent)")

    locArg('-rs_tolerance', nargs='?', type=float, default=2.5,
           help="RANSAC tie point filtering: percentage tolerance for max_outlier_percentage (default: 2.5 percent)")

    parse_coreg_local.set_defaults(func=run_local_coreg)

    return parser


if __name__ == '__main__':
    from socket import gethostname
    from datetime import datetime as dt
    from getpass import getuser

    def wfa(p, c):
        try:
            with open(p, 'a') as of:
                of.write(c)
        except Exception:
            pass

    wfa('/misc/hy5/scheffler/tmp/crlf', '%s\t%s\t%s\t%s\n' % (dt.now(), getuser(), gethostname(), ' '.join(sys.argv)))

    argparser = get_arosics_argparser()
    parsed_args = argparser.parse_args()

    if len(sys.argv) == 1:
        # no arguments provided
        print(
            '======================================================================\n'
            '#                            AROSICS v%s                         #' % __version__ + '\n'
            '# An Automated and Robust Open-Source Image Co-Registration Software #\n'
            '#                for Multi-Sensor Satellite Data                     #\n'
            '#          - python implementation by Daniel Scheffler               #\n'
            '======================================================================\n')
        argparser.print_help()
    else:
        t0 = time.time()
        parsed_args.func(parsed_args)
        print('\ntotal processing time: %.2fs' % (time.time() - t0))

else:
    warnings.warn("The script 'arosics_cli.py' provides a command line argument parser for AROSICS and is not to be "
                  "used as a normal Python module.")
