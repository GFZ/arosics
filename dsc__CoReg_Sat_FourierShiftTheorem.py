# -*- coding: utf-8 -*-
from __future__ import (division, print_function, absolute_import, unicode_literals) #unicode_literals cause GDAL not to work properly

__author__ = "Daniel Scheffler"

import time
import sys

# custom
try:
    import pyfftw
except ImportError:
    print('PYFFTW library is missing. However, coregistration works. But in some cases it can be much slower.')
try:
    import gdal
    import osr
    import ogr
except ImportError:
    from osgeo import gdal
    from osgeo import osr
    from osgeo import ogr

# internal modules
from CoReg_Sat import COREG, __version__



if __name__ == '__main__':
    import argparse
    from socket import gethostname
    from datetime import datetime as dt
    from getpass import getuser
    from components.io import wfa

    wfa('/misc/hy5/scheffler/tmp/crlf', '%s\t%s\t%s\t%s\n' % (dt.now(), getuser(), gethostname(), ' '.join(sys.argv)))
    parser = argparse.ArgumentParser(
        prog='dsc__CoReg_Sat_FourierShiftTheorem.py',

        description='Perform subpixel coregistration of two satellite image datasets ' \
                    'using Fourier Shift Theorem proposed by Foroosh et al. 2002: ' \
                    'Foroosh, H., Zerubia, J. B., & Berthod, M. (2002). Extension of phase correlation to subpixel ' \
                    'registration. IEEE Transactions on Image Processing, 11(3), 188-199. doi:10.1109/83.988953); '\
                    'Python implementation by Daniel Scheffler (daniel.scheffler@gfz-potsdam.de)',

        epilog="DETAILED DESCRIPTION: The program detects and corrects global X/Y-shifts between two input images in "\
        "the subpixel scale, that are often present in satellite imagery. It does not correct scaling or rotation "\
        "issues and will not apply any higher grade transformation. Therefore it will also not correct for shifts "\
        "that are locally present in the input images. "
        "Prerequisites and hints: The input images can have any GDAL compatible image format "\
        "(http://www.gdal.org/formats_list.html). Both of them must be georeferenced. In case of ENVI files, this "\
        "means they must have a 'map info' and a 'coordinate system string' as attributes of their header file. "\
        "Different projection systems are currently not supported. The input images must have a geographic overlap "\
        "but clipping them to same geographical extent is NOT neccessary. Please do not perform any spatial "\
        "resampling of the input images before applying this algorithm. Any needed resampling of the data is done "\
        "automatically. Thus the input images can have different spatial resolutions. The current algorithm will not "\
        "perform any ortho-rectification. So please use ortho-rectified input data in order to prevent local shifts "\
        "in the output image. By default the calculated subpixel-shifts are applied to the header file of the output "\
        "image. No spatial resampling is done automatically as long as the both input images have the same "\
        "projection. If you need the output image to be aligned to the reference image coordinate grid (by using an "\
        "appropriate resampling algorithm), use the '-align_grids' option. The image overlap area is automatically "\
        "calculated. Thereby no-data regions within the images are standardly respected. Providing the map "\
        "coordinates of the actual data corners lets you save some calculation time, because in this case the "\
        "automatic algorithm can be skipped. The no-data value of each image is automatically derived from the image "\
        "corners. The verbose program mode gives some more output about the interim results, shows some figures and "\
        "writes the used footprint and overlap polygons to disk. The figures must be manually closed in in order to "\
        "continue the processing."
        "The following non-standard Python libraries are required: gdal, osr, ogr, geopandas, rasterio, pykrige, "\
        "argparse and shapely. pyfftw is optional but will speed up calculation.")

    parser.add_argument('path_im0',
                        type=str,
                        help='source path of reference image (any GDAL compatible image format is supported)')

    parser.add_argument('path_im1',
                        type=str,
                        help='source path of image to be shifted (any GDAL compatible image format is supported)')

    parser.add_argument('-o', nargs='?', dest='path_out', type=str,
                        help="target path of the coregistered image (default: /dir/of/im1/<im1>__shifted_to__<im0>.bsq)",
                        default='.')

    parser.add_argument('-br', nargs='?',type=int, help='band of reference image to be used for matching '\
                       '(starts with 1; default: 1)', default=1)

    parser.add_argument('-bs', nargs='?',type=int, help='band of shift image to be used for matching '\
                       '(starts with 1; default: 1)', default=1)

    parser.add_argument('-wp', nargs=2, metavar=('X', 'Y'),type=float,help="custom matching window position as map "\
                       "values in the same projection like the reference image "\
                       "(default: central position of image overlap)", default=(None,None))

    parser.add_argument('-ws', nargs=2, metavar=('X size', 'Y size'),type=float,
                       help="custom matching window size [pixels] (default: (512,512))", default=(512,512))

    parser.add_argument('-max_iter', nargs='?', type=int,help="maximum number of iterations for matching (default: 5)",
                        default=5)

    parser.add_argument('-max_shift', nargs='?', type=int,help="maximum shift distance in reference image pixel units "\
                       "(default: 5 px)", default=5)

    parser.add_argument('-align_grids', nargs='?',type=int, help='align the coordinate grids of the image to be '\
           'and the reference image (default: 0)', default=0, choices=[0,1])

    parser.add_argument('-match_gsd', nargs='?',type=int, help='match the output pixel size to pixel size of the '\
        'reference image (default: 0)', default=0, choices=[0,1])

    parser.add_argument('-out_gsd', nargs=2,type=float, help='xgsd ygsd: set the output pixel size in map units'\
        '(default: original pixel size of the image to be shifted)',metavar=('xgsd','ygsd'))

    parser.add_argument('-cor0', nargs=8,type=float, help="map coordinates of data corners within reference image: ",
        metavar=tuple("UL-X UL-Y UR-X UR-Y LR-X LR-Y LL-X LL-Y".split(' ')),default=None)

    parser.add_argument('-cor1', nargs=8,type=float,help="map coordinates of data corners within image to be shifted: ",
        metavar=tuple("UL-X UL-Y UR-X UR-Y LR-X LR-Y LL-X LL-Y".split(' ')),default=None)

    parser.add_argument('-nodata', nargs=2,type=float, metavar=('im0','im1'),
                        help='no data values for reference image and image to be shifted',default=(None,None))

    parser.add_argument('-calc_cor', nargs='?',type=int, choices=[0,1],default=1, help="calculate true positions of "\
                        "the dataset corners in order to get a useful matching window position within the actual "\
                        "image overlap (default: 1; deactivated if '-cor0' and '-cor1' are given")

    parser.add_argument('-mp', nargs='?', type=int, help='enable multiprocessing (default: 1)', default=1, choices=[0,1])

    parser.add_argument('-bin_ws', nargs='?', type=int, help='use binary X/Y dimensions for the matching window '
                                                             '(default: 1)', default=1, choices=[0, 1])

    parser.add_argument('-quadratic_win', nargs='?', type=int, help='force a quadratic matching window (default: 1)',
                        default=1, choices=[0, 1])

    parser.add_argument('-v', nargs='?',type=int, help='verbose mode (default: 0)', default=0, choices=[0,1])

    parser.add_argument('-vo', nargs='?',type=str, help='an optional output directory for intermediate results '
                        '(if not given, no intermediate results are written to disk)', default=0, choices=[0,1])

    parser.add_argument('-q', nargs='?',type=int, help='quiet mode (default: 0)', default=0, choices=[0,1])

    parser.add_argument('-ignore_errors', nargs='?',type=int, help='Useful for batch processing. (default: 0) '
                        'In case of error COREG.success == False and COREG.x_shift_px/COREG.y_shift_px is None',
                        default=0, choices=[0,1])

    parser.add_argument('--version', action='version', version=__version__)

    args = parser.parse_args()

    print('==================================================================\n'
          '#          SUBPIXEL COREGISTRATION FOR SATELLITE IMAGERY         #\n'
          '#          - algorithm proposed by Foroosh et al. 2002           #\n'
          '#          - python implementation by Daniel Scheffler           #\n'
          '==================================================================\n')

    t0 = time.time()
    COREG_obj = COREG(args.path_im0,
                      args.path_im1,
                      path_out         = args.path_out,
                      r_b4match        = args.br,
                      s_b4match        = args.bs,
                      wp               = args.wp,
                      ws               = args.ws,
                      max_iter         = args.max_iter,
                      max_shift        = args.max_shift,
                      align_grids      = args.align_grids,
                      match_gsd        = args.match_gsd,
                      out_gsd          = args.out_gsd,
                      data_corners_im0 = args.cor0,
                      data_corners_im1 = args.cor1,
                      nodata           = args.nodata,
                      calc_corners     = args.calc_cor,
                      multiproc        = args.mp,
                      binary_ws        = args.bin_ws,
                      v                = args.v,
                      path_verbose_out = args.vo,
                      q                = args.q,
                      ignore_errors    = args.ignore_errors)
    COREG_obj.calculate_spatial_shifts()
    COREG_obj.correct_shifts()
    print('\ntotal processing time: %.2fs' %(time.time()-t0))

else:
    pass # TODO