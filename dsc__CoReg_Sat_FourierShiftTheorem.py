# -*- coding: utf-8 -*-
from __future__ import (division, print_function,absolute_import) #unicode_literals cause GDAL not to work properly

__author__ = "Daniel Scheffler"
__version__= "2016-09-26_01"

import collections
import multiprocessing
import os
import re
import shutil
import subprocess
import warnings
import time
import sys
import numpy as np
import tempfile

from shapely.geometry import Point, box
from geopandas  import GeoDataFrame
from pykrige.ok import OrdinaryKriging
import rasterio

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


from py_tools_ds.ptds import GeoArray


if __name__=='__main__':
    from components import geometry  as GEO
    from components import plotting  as PLT
    from components import utilities as UTL
    from components import io        as IO

else:
    from .components import geometry  as GEO
    from .components import plotting  as PLT
    from .components import utilities as UTL
    from .components import io        as IO



class imParamObj(object):
    def __init__(self, CoReg_params, imID):
        assert imID in ['ref', 'shift']
        self.imName = 'reference image' if imID == 'ref' else 'image to be shifted'

        # set GeoArray
        get_geoArr = lambda p: GeoArray(p) if not isinstance(p,GeoArray) else p
        self.GeoArray = get_geoArr(CoReg_params['im_ref']) if imID == 'ref' else get_geoArr(CoReg_params['im_tgt'])

        # set title to be used in plots
        self.title = os.path.basename(self.GeoArray.filePath) if self.GeoArray.filePath else self.imName

        # set params
        self.prj   = self.GeoArray.projection
        self.gt    = self.GeoArray.geotransform
        self.xgsd  = abs(self.gt[1])
        self.ygsd  = abs(self.gt[5])
        shape      = self.GeoArray.shape
        self.rows  = shape[0]
        self.cols  = shape[1]
        self.bands = shape[2] if len(shape)==3 else 1

        # validate params
        assert self.prj, 'The %s has no projection.' % self.imName
        assert not re.search('LOCAL_CS', self.prj), 'The %s is not georeferenced.' % self.imName
        assert self.gt, 'The %s has no map information.' % self.imName

        # set band4match
        self.band4match = CoReg_params['r_b4match'] if imID == 'ref' else CoReg_params['s_b4match']
        assert self.bands >= self.band4match >= 1, "The %s has %s %s. So its band number to match must be %s%s." \
            % (self.imName, self.bands, 'bands' if self.bands > 1 else 'band', 'between 1 and '
                if self.bands > 1 else '', self.bands)

        # set nodata
        if CoReg_params['nodata'][0] and CoReg_params['nodata'][1]:
            self.nodata = CoReg_params['nodata'][0 if imID == 'ref' else 1]
        elif self.GeoArray.filePath:
            self.nodata = GEO.find_noDataVal(self.GeoArray.filePath)
        else: #FIXME
            warnings.warn('Automatic nodata value detection for numpy array is currently not implemented. '
                          'Please provide a no data value. Otherwise it is set to None.')

        # set corner coords
        given_corner_coord = CoReg_params['data_corners_%s' % ('im0' if imID == 'ref' else 'im1')]
        if given_corner_coord is None:
            if CoReg_params['calc_corners']:
                if not CoReg_params['q']:
                    print('Calculating actual data corner coordinates for %s...' % self.imName)
                self.corner_coord = GEO.get_true_corner_mapXY(
                    self.GeoArray, self.band4match, self.nodata, CoReg_params['multiproc'])
            else:
                self.corner_coord = GEO.get_corner_coordinates(gt=self.GeoArray.geotransform,
                                                               cols=self.cols,rows=self.rows)
        else:
            self.corner_coord = given_corner_coord

        # set footprint polygon
        self.poly = GEO.get_footprint_polygon(self.corner_coord)
        for XY in self.corner_coord:
            assert self.poly.contains(Point(XY)) or self.poly.touches(Point(XY)), \
                "The corner position '%s' is outside of the %s." % (XY, self.imName)

        if not CoReg_params['q']: print('Corner coordinates of %s:\n\t%s' % (self.imName, self.corner_coord))



class COREG(object):
    def __init__(self, im_ref, im_tgt, path_out='.', r_b4match=1, s_b4match=1, wp=(None,None), ws=(512, 512),
                 max_iter=5, max_shift=5, align_grids=False, match_gsd=False, out_gsd=None, data_corners_im0=None,
                 data_corners_im1=None, nodata=(None,None), calc_corners=True, multiproc=True, binary_ws=True,
                 force_quadratic_win=True, v=False, q=False, ignore_errors=False): # FIXME path_im0/1 can also be a GeoArray
        """
        :param im_ref(str, GeoArray):   source path of reference image (any GDAL compatible image format is supported)
        :param im_ref(str, GeoArray):   source path of image to be shifted (any GDAL compatible image format is supported)
        :param path_out(str):           target path of the coregistered image
                                        (default: /dir/of/im1/<im1>__shifted_to__<im0>.bsq)
        :param r_b4match(int):          band of reference image to be used for matching (starts with 1; default: 1)
        :param s_b4match(int):          band of shift image to be used for matching (starts with 1; default: 1)
        :param wp(tuple):               custom matching window position as map values in the same projection like the
                                        reference image (default: central position of image overlap)
        :param ws(tuple):               custom matching window size [pixels] (default: (512,512))
        :param max_iter(int):           maximum number of iterations for matching (default: 5)
        :param max_shift(int):          maximum shift distance in reference image pixel units (default: 5 px)
        :param align_grids(bool):       align the coordinate grids of the image to be and the reference image (default: 0)
        :param match_gsd(bool):         match the output pixel size to pixel size of the reference image (default: 0)
        :param out_gsd(tuple):          xgsd ygsd: set the output pixel size in map units
                                        (default: original pixel size of the image to be shifted)
        :param resamp_alg(str)          the resampling algorithm to be used if neccessary
        :param data_corners_im0(list):  map coordinates of data corners within reference image
        :param data_corners_im1(list):  map coordinates of data corners within image to be shifted
        :param nodata(tuple):           no data values for reference image and image to be shifted
        :param calc_corners(bool):      calculate true positions of the dataset corners in order to get a useful
                                        matching window position within the actual image overlap
                                        (default: 1; deactivated if '-cor0' and '-cor1' are given
        :param multiproc(bool):         enable multiprocessing (default: 1)
        :param binary_ws(bool):         use binary X/Y dimensions for the matching window (default: 1)
        :param force_quadratic_win(bool):   force a quadratic matching window (default: 1)
        :param v(bool):                 verbose mode (default: 0)
        :param q(bool):                 quiet mode (default: 0)
        :param ignore_errors(bool):     Useful for batch processing. (default: 0)
                                        In case of error COREG.success == False and COREG.x_shift_px/COREG.y_shift_px
                                        is None
        """

        # FIXME add resamp_alg  # FIXME resamp_alg: "average" causes sinus-shaped patterns in fft images
        self.params              = dict([x for x in locals().items() if x[0] != "self"])

        if match_gsd and out_gsd: warnings.warn("'-out_gsd' is ignored because '-match_gsd' is set.\n")
        if out_gsd:  assert isinstance(out_gsd, list) and len(out_gsd) == 2, 'out_gsd must be a list with two values.'
        if data_corners_im0 and not isinstance(data_corners_im0[0],list): # group if not [[x,y],[x,y]..] but [x,y,x,y,]
            data_corners_im0 = [data_corners_im0[i:i+2] for i in range(0, len(data_corners_im0), 2)]
        if data_corners_im1 and not isinstance(data_corners_im1[0],list): # group if not [[x,y],[x,y]..]
            data_corners_im1 = [data_corners_im1[i:i+2] for i in range(0, len(data_corners_im1), 2)]
        if nodata: assert isinstance(nodata, tuple) and len(nodata) == 2, "'nodata' must be a tuple with two values." \
                                                                          "Got %s with length %s." %(type(nodata),len(nodata))

        self.path_out            = path_out            # updated by self.set_outpathes
        self.win_pos_XY          = wp                  # updated by self.get_opt_winpos_winsize()
        self.win_size_XY         = ws                  # updated by self.get_opt_winpos_winsize()
        self.max_iter            = max_iter
        self.max_shift           = max_shift
        self.align_grids         = align_grids
        self.match_gsd           = match_gsd
        self.out_gsd             = out_gsd
        self.calc_corners        = calc_corners
        self.mp                  = multiproc
        self.bin_ws              = binary_ws
        self.force_quadratic_win = force_quadratic_win
        self.v                   = v
        self.q                   = q if not v else False
        self.ignErr              = ignore_errors
        self.verbose_out         = None                # set by self.set_outpathes
        self.max_win_sz_changes  = 3                   # TODO: änderung der window size, falls nach max_iter kein valider match gefunden
        self.ref                 = None                # set by self.get_image_params
        self.shift               = None                # set by self.get_image_params
        self.matchWin            = None                # set by self.get_clip_window_properties()
        self.otherWin            = None                # set by self.get_clip_window_properties()
        self.imfft_gsd           = None                # set by self.get_clip_window_properties()
        self.fftw_win_size_YX    = None                # set by calc_shifted_cross_power_spectrum

        self.x_shift_px          = None                # always in shift image units (image coords) # set by calculate_spatial_shifts()
        self.y_shift_px          = None                # always in shift image units (image coords) # set by calculate_spatial_shifts()
        self.x_shift_map         = None                # set by self.get_updated_map_info()
        self.y_shift_map         = None                # set by self.get_updated_map_info()
        self.updated_map_info    = None                # set by self.get_updated_map_info()

        self.tracked_errors      = []                  # expanded each time an error occurs
        self.success             = False               # default

        gdal.AllRegister()
        self.get_image_params()
        self.set_outpathes(im_ref, im_tgt)
        self.grid2use                 = 'ref' if self.shift.xgsd <= self.ref.xgsd else 'shift'
        if self.v: print('resolutions: ', self.ref.xgsd, self.shift.xgsd)

        overlap_tmp                   = GEO.get_overlap_polygon(self.ref.poly, self.shift.poly, self.v)
        self.overlap_poly             = overlap_tmp['overlap poly'] # has to be the reference projection
        assert self.overlap_poly, 'The input images have no spatial overlap.'
        self.overlap_percentage       = overlap_tmp['overlap percentage']
        self.overlap_area             = overlap_tmp['overlap area']

        if self.v: IO.write_shp(os.path.join(self.verbose_out, 'poly_imref.shp'),    self.ref.poly,     self.ref.prj)
        if self.v: IO.write_shp(os.path.join(self.verbose_out, 'poly_im2shift.shp'), self.shift.poly,   self.shift.prj)
        if self.v: IO.write_shp(os.path.join(self.verbose_out, 'overlap_poly.shp'),  self.overlap_poly, self.ref.prj)

        ### FIXME: transform_mapPt1_to_mapPt2(im2shift_center_map, ds_imref.GetProjection(), ds_im2shift.GetProjection()) # später basteln für den fall, dass projektionen nicht gleich sind

        # get_clip_window_properties
        self.get_opt_winpos_winsize()
        if not self.q: print('Matching window position (X,Y): %s/%s' % (self.win_pos_XY[0], self.win_pos_XY[1]))
        self.get_clip_window_properties()

        if self.v and self.matchWin.mapPoly:
            IO.write_shp(os.path.join(self.verbose_out, 'poly_matchWin.shp'), self.matchWin.mapPoly, self.matchWin.prj)

        self.success = None if self.matchWin.boxMapYX else False


    def set_outpathes(self, im_ref, im_tgt):
        path_im_ref         = (im_ref if isinstance(im_ref, GeoArray) else GeoArray(im_ref)).filePath
        path_im_tgt         = (im_tgt if isinstance(im_tgt, GeoArray) else GeoArray(im_tgt)).filePath
        dir_out, fName_out  = os.path.split(self.path_out)
        if not dir_out:
            if not path_im_ref:
                raise #FIXME
            else:
                dir_out = os.path.dirname(path_im_ref)
        get_baseN           = lambda path: os.path.splitext(os.path.basename(path))[0]
        fName_out           = fName_out if not fName_out in ['.',''] else '%s__shifted_to__%s.bsq' \
                                %(get_baseN(path_im_tgt), get_baseN(path_im_ref))
        self.path_out       = os.path.abspath(os.path.join(dir_out,fName_out))
        assert ' ' not in self.path_out, \
            "The path of the output image contains whitespaces. This is not supported by GDAL."
        self.verbose_out    = os.path.abspath(os.path.join(dir_out, 'CoReg_verboseOut__%s__shifted_to__%s' \
                                %(get_baseN(path_im_tgt), get_baseN(path_im_ref))))
        if self.v and not os.path.isdir(self.verbose_out): os.makedirs(self.verbose_out)


    def get_image_params(self):
        self.ref   = imParamObj(self.params,'ref')
        self.shift = imParamObj(self.params,'shift')
        assert GEO.get_proj4info(proj=self.ref.prj) == GEO.get_proj4info(proj=self.shift.prj), \
            'Input projections are not equal. Different projections are currently not supported.'


    def get_opt_winpos_winsize(self):
        # type: (tuple,tuple) -> tuple,tuple
        """Calculates optimal window position and size in reference image units according to DGM, cloud_mask and
        trueCornerLonLat."""
        # dummy algorithm: get center position of overlap instead of searching ideal window position in whole overlap
        # TODO automatischer Algorithmus zur Bestimmung der optimalen Window Position

        wp = tuple(self.win_pos_XY)
        assert type(self.win_pos_XY) in [tuple,list,np.ndarray],\
            'The window position must be a tuple of two elements. Got %s with %s elements.' %(type(wp),len(wp))
        wp = tuple(wp)

        if None in wp:
            overlap_center_pos_x, overlap_center_pos_y = self.overlap_poly.centroid.coords.xy
            wp = (wp[0] if wp[0] else overlap_center_pos_x[0]), (wp[1] if wp[1] else overlap_center_pos_y[0])

        assert self.overlap_poly.contains(Point(wp)), 'The provided window position is outside of the overlap area of '\
                                                      'the two input images. Check the coordinates.'
        self.win_pos_XY  = wp
        self.win_size_XY = (int(self.win_size_XY[0]), int(self.win_size_XY[1])) if self.win_size_XY else (512,512)


    def get_clip_window_properties(self):
        """Calculate all properties of the matching window and the other window. These windows are used to read the
        corresponding image positions in the reference and the target image.
        hint: Even if X- and Y-dimension of the target window is equal, the output window can be NOT quadratic!
        """

        wpX,wpY             = self.win_pos_XY
        wsX,wsY             = self.win_size_XY
        ref_wsX, ref_wsY    = (wsX*self.ref.xgsd  , wsY*self.ref.ygsd)   # image units -> map units
        shift_wsX,shift_wsY = (wsX*self.shift.xgsd, wsY*self.shift.ygsd) # image units -> map units
        ref_box_kwargs      = {'wp':(wpX,wpY),'ws':(ref_wsX,ref_wsY)    ,'gt':self.ref.gt  }
        shift_box_kwargs    = {'wp':(wpX,wpY),'ws':(shift_wsX,shift_wsY),'gt':self.shift.gt}
        matchWin            = GEO.boxObj(**ref_box_kwargs)   if self.grid2use=='ref' else GEO.boxObj(**shift_box_kwargs)
        otherWin            = GEO.boxObj(**shift_box_kwargs) if self.grid2use=='ref' else GEO.boxObj(**ref_box_kwargs)
        overlapWin          = GEO.boxObj(mapPoly=self.overlap_poly,gt=self.ref.gt)

        # clip matching window to overlap area
        matchWin.mapPoly = matchWin.mapPoly.intersection(overlapWin.mapPoly)

        # move matching window to imref grid or im2shift grid
        mW_rows, mW_cols = (self.ref.rows, self.ref.cols) if self.grid2use == 'ref' else (self.shift.rows, self.shift.cols)
        matchWin.mapPoly = GEO.move_shapelyPoly_to_image_grid(matchWin.mapPoly, matchWin.gt, mW_rows, mW_cols, 'NW')

        # check, ob durch Verschiebung auf Grid das matchWin außerhalb von overlap_poly geschoben wurde
        if not matchWin.mapPoly.within(overlapWin.mapPoly):
            # matchPoly weiter verkleinern # 1 px buffer reicht, weil window nur auf das Grid verschoben wurde
            xLarger,yLarger = matchWin.is_larger_DimXY(overlapWin.boundsIm)
            matchWin.buffer_imXY(-1 if xLarger else 0, -1 if yLarger else 0)

        # matching_win direkt auf grid2use (Rundungsfehler bei Koordinatentrafo beseitigen)
        matchWin.imPoly = GEO.round_shapelyPoly_coords(matchWin.imPoly, precision=0, out_dtype=int)

        # Check, ob match Fenster größer als anderes Fenster
        if not (matchWin.mapPoly.within(otherWin.mapPoly) or matchWin.mapPoly==otherWin.mapPoly):
            # dann für anderes Fenster kleinstes Fenster finden, das match-Fenster umgibt
            otherWin.boxImYX = GEO.get_smallest_boxImYX_that_contains_boxMapYX(matchWin.boxMapYX,otherWin.gt)

        # evtl. kann es sein, dass bei Shift-Fenster-Vergrößerung das shift-Fenster zu groß für den overlap wird
        while not otherWin.mapPoly.within(overlapWin.mapPoly):
            # -> match Fenster verkleinern und neues anderes Fenster berechnen
            xLarger, yLarger = otherWin.is_larger_DimXY(overlapWin.boundsIm)
            matchWin.buffer_imXY(-1 if xLarger else 0, -1 if yLarger else 0)
            otherWin.boxImYX = GEO.get_smallest_boxImYX_that_contains_boxMapYX(matchWin.boxMapYX,otherWin.gt)

        # check results
        assert matchWin.mapPoly.within(otherWin.mapPoly)
        assert otherWin.mapPoly.within(overlapWin.mapPoly)

        self.imfft_gsd              = self.ref.xgsd       if self.grid2use =='ref' else self.shift.xgsd
        self.ref.win,self.shift.win = (matchWin,otherWin) if self.grid2use =='ref' else (otherWin,matchWin)
        self.matchWin,self.otherWin = matchWin, otherWin
        self.ref.  win.size_YX      = tuple([int(i) for i in self.ref.  win.imDimsYX])
        self.shift.win.size_YX      = tuple([int(i) for i in self.shift.win.imDimsYX])
        match_win_size_XY           = tuple(reversed([int(i) for i in matchWin.imDimsYX]))
        if match_win_size_XY != self.win_size_XY:
            print('Target window size %s not possible due to too small overlap area or window position too close '
                  'to an image edge. New matching window size: %s.' %(self.win_size_XY,match_win_size_XY))
        #IO.write_shp('/misc/hy5/scheffler/Temp/matchMapPoly.shp', matchWin.mapPoly,matchWin.prj)
        #IO.write_shp('/misc/hy5/scheffler/Temp/otherMapPoly.shp', otherWin.mapPoly,otherWin.prj)


    def get_image_windows_to_match(self):
        """Reads the matching window and the other window using subset read, and resamples the other window to the
        resolution and the pixel grid of the matching window. The result consists of two images with the same
        dimensions and exactly the same corner coordinates."""

        is_avail_rsp_average = int(gdal.VersionInfo()[0]) >= 2 #FIXME move to output image generation
        if not is_avail_rsp_average:
            warnings.warn("The GDAL version on this server does not yet support the resampling algorithm "
                          "'average'. This can affect the correct detection of subpixel shifts. To avoid this "
                          "please update GDAL to a version above 2.0.0!")

        self.matchWin.imParams = self.ref   if self.grid2use=='ref' else self.shift
        self.otherWin.imParams = self.shift if self.grid2use=='ref' else self.ref

        # matchWin per subset-read einlesen -> self.matchWin.data
        rS, rE, cS, cE = GEO.get_GeoArrayPosition_from_boxImYX(self.matchWin.boxImYX)
        assert np.array_equal(np.abs(np.array([rS,rE,cS,cE])), np.array([rS,rE,cS,cE])), \
            'Got negative values in gdalReadInputs for %s.' %self.matchWin.imParams.imName
        self.matchWin.data = self.matchWin.imParams.GeoArray[rS:rE,cS:cE, self.matchWin.imParams.band4match-1]

        # otherWin per subset-read einlesen
        rS, rE, cS, cE = GEO.get_GeoArrayPosition_from_boxImYX(self.otherWin.boxImYX)
        assert np.array_equal(np.abs(np.array([rS,rE,cS,cE])), np.array([rS,rE,cS,cE])), \
            'Got negative values in gdalReadInputs for %s.' %self.otherWin.imParams.imName
        self.otherWin.data = self.otherWin.imParams.GeoArray[rS:rE, cS:cE, self.otherWin.imParams.band4match - 1]

        if self.v:
            print('Original matching windows:')
            ref_data, shift_data =  (self.matchWin.data, self.otherWin.data) if self.grid2use=='ref' else \
                                    (self.otherWin.data, self.matchWin.data)
            PLT.subplot_imshow([ref_data, shift_data],[self.ref.title,self.shift.title], grid=True)

        otherWin_subgt = GEO.get_subset_GeoTransform(self.otherWin.gt, self.otherWin.boxImYX)

        # resample otherWin.data to the resolution of matchWin AND make sure the pixel edges are identical
        # (in order to make each image show the same window with the same coordinates)
        # rsp_algor = 5 if is_avail_rsp_average else 2 # average if possible else cubic # OLD
        # TODO replace cubic resampling by PSF resampling - average resampling leads to sinus like distortions in the fft image that make a precise coregistration impossible. Thats why there is currently no way around cubic resampling.
        tgt_xmin,tgt_xmax,tgt_ymin,tgt_ymax = self.matchWin.boundsMap
        self.otherWin.data = GEO.warp_ndarray(self.otherWin.data, otherWin_subgt, self.otherWin.imParams.prj,
                                              self.matchWin.imParams.prj, out_gsd=(self.imfft_gsd, self.imfft_gsd),
                                              out_bounds=([tgt_xmin, tgt_ymin, tgt_xmax, tgt_ymax]),
                                              rspAlg='cubic', in_nodata=self.otherWin.imParams.nodata)[0]

        if self.matchWin.data.shape != self.otherWin.data.shape:
            self.tracked_errors.append(
                RuntimeError('Bad output of get_image_windows_to_match. Reference image shape is %s whereas shift '
                             'image shape is %s.' % (self.matchWin.data.shape, self.otherWin.data.shape)))
            raise self.tracked_errors[-1]
        rows, cols = [i if i % 2 == 0 else i - 1 for i in self.matchWin.data.shape]
        self.matchWin.data, self.otherWin.data = self.matchWin.data[:rows, :cols], self.otherWin.data[:rows, :cols]

        assert self.matchWin.data is not None and self.otherWin.data is not None, 'Creation of matching windows failed.'


    @staticmethod
    def shrink_winsize_to_binarySize(win_shape_YX, target_size=None):
        # type: (tuple, tuple, int , int) -> tuple
        """Shrinks a given window size to the closest binary window size (a power of 2) -
        separately for X- and Y-dimension.

        :param win_shape_YX:    <tuple> source window shape as pixel units (rows,colums)
        :param target_size:     <tuple> source window shape as pixel units (rows,colums)
        """

        binarySizes   = [2**i for i in range(3,14)] # [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
        possibSizes_X = [i for i in binarySizes if i <= win_shape_YX[1]]
        possibSizes_Y = [i for i in binarySizes if i <= win_shape_YX[0]]
        if possibSizes_X and possibSizes_Y:
            tgt_size_X,tgt_size_Y = target_size if target_size else (max(possibSizes_X),max(possibSizes_Y))
            closest_to_target_X = int(min(possibSizes_X, key=lambda x:abs(x-tgt_size_X)))
            closest_to_target_Y = int(min(possibSizes_Y, key=lambda y:abs(y-tgt_size_Y)))
            return closest_to_target_Y,closest_to_target_X
        else:
            return None


    def calc_shifted_cross_power_spectrum(self,im0=None,im1=None,precision=np.complex64):
        """Calculates shifted cross power spectrum for quantifying x/y-shifts.

            :param precision:   to be quantified as a datatype
            :param im0:         reference image
            :param im1:         subject image to shift
            :return:            2D-numpy-array of the shifted cross power spectrum
        """

        im0 = im0 if im0 is not None else self.ref.win.data
        im1 = im1 if im1 is not None else self.shift.win.data
        assert im0.shape == im1.shape, 'The reference and the target image must have the same dimensions.'
        if im0.shape[0]%2!=0: warnings.warn('Odd row count in one of the match images!')
        if im1.shape[1]%2!=0: warnings.warn('Odd column count in one of the match images!')

        wsYX = self.shrink_winsize_to_binarySize(im0.shape) if self.bin_ws              else im0.shape
        wsYX = ((min(wsYX),) * 2                            if self.force_quadratic_win else wsYX) if wsYX else None

        if wsYX:
            time0 = time.time()
            if self.v: print('final window size: %s/%s (X/Y)' % (wsYX[1], wsYX[0]))
            center_YX = np.array(im0.shape)/2
            xmin,xmax,ymin,ymax = int(center_YX[1]-wsYX[1]/2), int(center_YX[1]+wsYX[1]/2),\
                                  int(center_YX[0]-wsYX[0]/2), int(center_YX[0]+wsYX[0]/2)
            in_arr0  = im0[ymin:ymax,xmin:xmax].astype(precision)
            in_arr1  = im1[ymin:ymax,xmin:xmax].astype(precision)

            if self.v:
                PLT.subplot_imshow([in_arr0.astype(np.float32), in_arr1.astype(np.float32)],
                               ['FFTin '+self.ref.title,'FFTin '+self.shift.title], grid=True)

            if 'pyfft' in globals():
                fft_arr0 = pyfftw.FFTW(in_arr0,np.empty_like(in_arr0), axes=(0,1))()
                fft_arr1 = pyfftw.FFTW(in_arr1,np.empty_like(in_arr1), axes=(0,1))()
            else:
                fft_arr0 = np.fft.fft2(in_arr0)
                fft_arr1 = np.fft.fft2(in_arr1)
            if self.v: print('forward FFTW: %.2fs' %(time.time() -time0))

            eps = np.abs(fft_arr1).max() * 1e-15
            # cps == cross-power spectrum of im0 and im2

            temp = (fft_arr0 * fft_arr1.conjugate()) / (np.abs(fft_arr0) * np.abs(fft_arr1) + eps)

            time0 = time.time()
            if 'pyfft' in globals():
                ifft_arr = pyfftw.FFTW(temp,np.empty_like(temp), axes=(0,1),direction='FFTW_BACKWARD')()
            else:
                ifft_arr = np.fft.ifft2(temp)
            if self.v: print('backward FFTW: %.2fs' %(time.time() -time0))

            cps = np.abs(ifft_arr)
            # scps = shifted cps
            scps = np.fft.fftshift(cps)
            if self.v:
                PLT.subplot_imshow([in_arr0.astype(np.uint16), in_arr1.astype(np.uint16), fft_arr0.astype(np.uint8),
                                fft_arr1.astype(np.uint8), scps], titles=['matching window im0', 'matching window im1',
                                "fft result im0", "fft result im1", "cross power spectrum"], grid=True)
                PLT.subplot_3dsurface(scps.astype(np.float32))
        else:
            self.tracked_errors.append(
                RuntimeError('The matching window became too small for calculating a reliable match. Matching failed.'))
            if self.ignErr:
                scps = None
            else:
                raise self.tracked_errors[-1]

        self.fftw_win_size_YX = wsYX
        return scps


    @staticmethod
    def get_peakpos(scps):
        max_flat_idx = np.argmax(scps)
        return np.unravel_index(max_flat_idx, scps.shape)


    @staticmethod
    def get_shifts_from_peakpos(peakpos,arr_shape):
        y_shift = peakpos[0]-arr_shape[0]//2
        x_shift = peakpos[1]-arr_shape[1]//2
        return x_shift,y_shift


    @staticmethod
    def clip_image(im,center_YX,winSzYX):
        get_bounds = lambda YX,wsY,wsX: (int(YX[1]-(wsX/2)),int(YX[1]+(wsX/2)),int(YX[0]-(wsY/2)),int(YX[0]+(wsY/2)))
        wsY,wsX    = winSzYX
        xmin,xmax,ymin,ymax = get_bounds(center_YX,wsY,wsX)
        return im[ymin:ymax,xmin:xmax]


    def get_grossly_deshifted_images(self,im0,im1,x_intshift,y_intshift):
        # get_grossly_deshifted_im0
        old_center_YX = np.array(im0.shape)/2
        new_center_YX = [old_center_YX[0]+y_intshift, old_center_YX[1]+x_intshift]

        x_left  = new_center_YX[1]
        x_right = im0.shape[1]-new_center_YX[1]
        y_above = new_center_YX[0]
        y_below = im0.shape[0]-new_center_YX[0]
        maxposs_winsz = 2*min(x_left,x_right,y_above,y_below)

        gdsh_im0 = self.clip_image(im0,new_center_YX,[maxposs_winsz,maxposs_winsz])

        # get_corresponding_im1_clip
        crsp_im1  = self.clip_image(im1,np.array(im1.shape)/2,gdsh_im0.shape)

        if self.v:
            PLT.subplot_imshow([self.clip_image(im0, old_center_YX, gdsh_im0.shape), crsp_im1],
                           titles=['reference original', 'target'], grid=True)
            PLT.subplot_imshow([gdsh_im0, crsp_im1], titles=['reference virtually shifted', 'target'], grid=True)
        return gdsh_im0,crsp_im1


    @staticmethod
    def find_side_maximum(scps, v=0):
        centerpos = [scps.shape[0]//2, scps.shape[1]//2]
        profile_left  = scps[ centerpos [0]  ,:centerpos[1]+1]
        profile_right = scps[ centerpos [0]  , centerpos[1]:]
        profile_above = scps[:centerpos [0]+1, centerpos[1]]
        profile_below = scps[ centerpos [0]: , centerpos[1]]

        if v:
            max_count_vals = 10
            PLT.subplot_2dline([[range(len(profile_left)) [-max_count_vals:], profile_left[-max_count_vals:]],
                                [range(len(profile_right))[:max_count_vals] , profile_right[:max_count_vals]],
                                [range(len(profile_above))[-max_count_vals:], profile_above[-max_count_vals:]],
                                [range(len(profile_below))[:max_count_vals:], profile_below[:max_count_vals]]],
                                titles =['Profile left', 'Profile right', 'Profile above', 'Profile below'],
                                shapetuple=(2,2))

        get_sidemaxVal_from_profile = lambda pf: np.array(pf)[::-1][1] if pf[0]<pf[-1] else np.array(pf)[1]
        sm_dicts_lr  = [{'side':si, 'value': get_sidemaxVal_from_profile(pf)} \
                        for pf,si in zip([profile_left,profile_right],['left','right'])]
        sm_dicts_ab  = [{'side':si, 'value': get_sidemaxVal_from_profile(pf)} \
                        for pf,si in zip([profile_above,profile_below],['above','below'])]
        sm_maxVal_lr = max([i['value'] for i in sm_dicts_lr])
        sm_maxVal_ab = max([i['value'] for i in sm_dicts_ab])
        sidemax_lr   = [sm for sm in sm_dicts_lr if sm['value'] is sm_maxVal_lr][0]
        sidemax_ab   = [sm for sm in sm_dicts_ab if sm['value'] is sm_maxVal_ab][0]
        sidemax_lr['direction_factor'] = {'left':-1, 'right':1} [sidemax_lr['side']]
        sidemax_ab['direction_factor'] = {'above':-1,'below':1} [sidemax_ab['side']]

        if v:
            print('Horizontal side maximum found %s. value: %s' %(sidemax_lr['side'],sidemax_lr['value']))
            print('Vertical side maximum found %s. value: %s' %(sidemax_ab['side'],sidemax_ab['value']))

        return sidemax_lr, sidemax_ab


    def calc_subpixel_shifts(self,scps):
        sidemax_lr, sidemax_ab = self.find_side_maximum(scps,self.v)
        x_subshift = (sidemax_lr['direction_factor']*sidemax_lr['value'])/(np.max(scps)+sidemax_lr['value'])
        y_subshift = (sidemax_ab['direction_factor']*sidemax_ab['value'])/(np.max(scps)+sidemax_ab['value'])
        return x_subshift, y_subshift


    @staticmethod
    def get_total_shifts(x_intshift,y_intshift,x_subshift,y_subshift):
        return x_intshift+x_subshift, y_intshift+y_subshift


    def calculate_spatial_shifts(self):
        if self.success is False: return None,None

        if self.q:  warnings.simplefilter('ignore')
        self.get_image_windows_to_match()
        #
        im0,im1 = self.ref.win.data, self.shift.win.data
        #write_envi(im0, '/misc/hy5/scheffler/Projekte/2016_03_Geometrie_Paper/v2_imref__%s_%s__ws%s.hdr' %(self.win_pos_XY,self.ref_win_size_YX[0]))
        #write_envi(im1, '/misc/hy5/scheffler/Projekte/2016_03_Geometrie_Paper/v2_imtgt__%s_%s__ws%s.hdr' %(self.win_pos_XY,self.shift_win_size_YX[0]))
        if self.v:
            print('Matching windows with equalized spatial resolution:')
            PLT.subplot_imshow([im0, im1], [self.ref.title, self.shift.title], grid=True)

        gsd_factor = self.imfft_gsd/self.shift.xgsd

        if self.v: print('gsd_factor',         gsd_factor)
        if self.v: print('imfft_gsd_mapvalues',self.imfft_gsd)

        scps = self.calc_shifted_cross_power_spectrum()
        x_shift_px,y_shift_px = None,None # defaults
        if scps is None:
            self.success = False
        else:
            peakpos = self.get_peakpos(scps)
            x_intshift, y_intshift  = [], []
            x_tempshift,y_tempshift = self.get_shifts_from_peakpos(peakpos,scps.shape)
            x_intshift.append(x_tempshift)
            y_intshift.append(y_tempshift)
            if not self.q: print('Detected integer shifts (X/Y):       %s/%s' %(x_tempshift,y_tempshift))

            count_iter = 1
            while [x_tempshift, y_tempshift] != [0,0]:
                count_iter +=1
                if not self.q: print('No clear match found yet. Jumping to iteration %s...' %count_iter)
                if len(x_intshift)>1:
                    x_tempshift = x_tempshift + x_intshift[-2]
                    y_tempshift = y_tempshift + y_intshift[-2]
                if not self.q: print('input shifts: ', x_tempshift, y_tempshift)

                if scps is None: self.success = False; break # FIXME sinnlos
                gdsh_im0,crsp_im1 = self.get_grossly_deshifted_images(im0, im1,x_tempshift,y_tempshift)
                scps              = self.calc_shifted_cross_power_spectrum(gdsh_im0,crsp_im1)

                if scps is None: self.success = False; break
                peakpos = self.get_peakpos(scps)
                x_tempshift,y_tempshift = self.get_shifts_from_peakpos(peakpos,scps.shape)
                x_intshift.append(x_tempshift)
                y_intshift.append(y_tempshift)
                if not self.q: print('Detected integer shifts (X/Y):       %s/%s' %(x_tempshift,y_tempshift))

                if count_iter > self.max_iter:
                    self.success = False
                    self.tracked_errors.append(RuntimeError('No match found in the given window.'))
                    if not self.ignErr:
                        raise self.tracked_errors[-1]
                    else:
                        warnings.warn('No match found in the given window.'); break

            if not self.success==False:
                print('COREG 642', x_intshift, y_intshift)
                x_intshift, y_intshift      = sum(x_intshift), sum(y_intshift) ## FIXME SUMME??
                x_subshift, y_subshift      = self.calc_subpixel_shifts(scps)
                print(x_intshift, y_intshift)
                x_totalshift, y_totalshift  = self.get_total_shifts(x_intshift,y_intshift,x_subshift,y_subshift)
                x_shift_px, y_shift_px      = x_totalshift*gsd_factor, y_totalshift*gsd_factor
                if not self.q:
                    print('Detected subpixel shifts (X/Y):      %s/%s' %(x_subshift,y_subshift))
                    print('Calculated total shifts in fft image units (X/Y):         %s/%s' %(x_totalshift,y_totalshift))
                    print('Calculated total shifts in reference image units (X/Y):   %s/%s' %(x_totalshift,y_totalshift))
                    print('Calculated total shifts in target image units (X/Y):      %s/%s' %(x_shift_px,y_shift_px))

                if max([abs(x_totalshift),abs(y_totalshift)]) > self.max_shift:
                    self.success = False
                    self.tracked_errors.append(
                        RuntimeError("The calculated shift (X: %s px / Y: %s px) is recognized as too large to "
                                     "be valid. If you know that it is valid, just set the '-max_shift' "
                                     "parameter to an appropriate value. Otherwise try to use a different window "
                                     "size for matching via the '-ws' parameter or define the spectral bands "
                                     "to be used for matching manually ('-br' and '-bs')."
                                     % (self.x_shift_px, self.y_shift_px)))
                    if not self.ignErr:
                        raise self.tracked_errors[-1]
                else:
                    self.success = True

        self.x_shift_px, self.y_shift_px = (x_shift_px,y_shift_px) if self.success else (None,None)
        if self.x_shift_px or self.y_shift_px:
            self.get_updated_map_info()

        warnings.simplefilter('default')


    def get_updated_map_info(self):
        original_map_info = GEO.geotransform2mapinfo(self.shift.gt, self.shift.prj)
        if not self.q: print('Original map info:', original_map_info)
        new_originY, new_originX = GEO.pixelToMapYX([self.x_shift_px,self.y_shift_px],
                                                     geotransform=self.shift.gt, projection=self.shift.prj)[0]
        self.x_shift_map = new_originX - self.shift.gt[0]
        self.y_shift_map = new_originY - self.shift.gt[3]
        if not self.q: print('Calculated map shifts (X,Y): ', self.x_shift_map,self.y_shift_map)

        self.updated_map_info    = original_map_info.copy()
        self.updated_map_info[3] = str(float(original_map_info[3]) + self.x_shift_map)
        self.updated_map_info[4] = str(float(original_map_info[4]) + self.y_shift_map)
        if not self.q: print('Updated map info:',self.updated_map_info)


    @property
    def coreg_info(self):
        return {
            'corrected_shifts_px'   : {'x':self.x_shift_px,  'y':self.y_shift_px },
            'corrected_shifts_map'  : {'x':self.x_shift_map, 'y':self.y_shift_map},
            'original map info'     : GEO.geotransform2mapinfo(self.shift.gt, self.shift.prj),
            'updated map info'      : self.updated_map_info,
            'reference projection'  : self.ref.prj,
            'reference geotransform': self.ref.gt,
            'reference grid'        : [ [self.ref.gt[0], self.ref.gt[0]+self.ref.gt[1]],
                                        [self.ref.gt[3], self.ref.gt[3]+self.ref.gt[5]] ],
            'reference extent'      : {'cols':self.ref.xgsd, 'rows':self.ref.ygsd}, # FIXME not needed anymore
            'success'               : self.success}


    def correct_shifts(self):
        #DS = DESHIFTER(self.shift.GeoArray, self.coreg_info, band2process=2)
        #get_deshift_results = lambda band:\
        #    DESHIFTER(self.shift.GeoArray, self.coreg_info, band2process=band).correct_shifts()
        #with multiprocessing.Pool() as pool:
        #    deshift_results = pool.map(get_deshift_results, range(self.shift.GeoArray.bands))
        DS = DESHIFTER(self.shift.GeoArray, self.coreg_info, out_gsd=self.out_gsd, align_grids=self.align_grids,
                       match_gsd=self.match_gsd, v=self.v, q=self.q)
        return DS.correct_shifts()


    def correct_shifts_OLD(self):
        if self.success:
            if not os.path.exists(os.path.dirname(self.path_out)): os.makedirs(os.path.dirname(self.path_out))
            equal_prj = GEO.get_proj4info(proj=self.ref.prj) == GEO.get_proj4info(proj=self.shift.prj)

            if equal_prj and not self.align_grids and not self.match_gsd and \
                self.out_gsd in [None,[self.shift.xgsd,self.shift.ygsd]]:
                self.shift_image_by_updating_map_info()
            elif equal_prj and self.align_grids: # match_gsd and out_gsd are respected
                self.align_coordinate_grids()
            else: # match_gsd and out_gsd are respected ### TODO: out_proj implementieren
                self.resample_without_grid_aligning()
        else:
            warnings.warn('No result written because detection of image displacements failed.')


    def shift_image_by_updating_map_info(self):
        if not self.q: print('\nWriting output...')
        ds_im2shift = gdal.Open(self.shift.path)
        if not ds_im2shift.GetDriver().ShortName == 'ENVI': # FIXME laaangsam
            if self.mp:
                IO.convert_gdal_to_bsq__mp(self.shift.path, self.path_out)
            else:
                os.system('gdal_translate -of ENVI %s %s' %(self.shift.path, self.path_out))
            file2getHdr = self.path_out
        else:
            shutil.copy(self.shift.path,self.path_out)
            file2getHdr = self.shift.path
        ds_im2shift = None

        path_hdr = '%s.hdr' %os.path.splitext(file2getHdr)[0]
        path_hdr = '%s.hdr' %file2getHdr if not os.path.exists(path_hdr) else path_hdr
        path_hdr = None if not os.path.exists(path_hdr) else path_hdr

        assert path_hdr, 'No header file found for %s. Applying shifts failed.' %file2getHdr

        new_originY, new_originX = GEO.pixelToMapYX([self.x_shift_px, self.y_shift_px],
                                                    geotransform=self.shift.gt, projection=self.shift.prj)[0]
        map_xshift,  map_yshift  = new_originX - self.shift.gt[0], new_originY - self.shift.gt[3]
        if not self.q: print('Calculated map shifts (X,Y): %s, %s' %(map_xshift,map_yshift))

        with open(path_hdr,'r') as inF:
            content = inF.read()
            map_info_line   = [i for i in content.split('\n') if re.search('map info [\w+]*=',i,re.I)][0]
            map_info_string = re.search('{([\w+\s\S]*)}',map_info_line,re.I).group(1)
            map_info        = [i.strip() for i in map_info_string.split(',')]
        if not self.q: print('Original map info: %s' %map_info)
        new_map_info = map_info
        new_map_info[3] = str(float(map_info[3]) + map_xshift)
        new_map_info[4] = str(float(map_info[4]) + map_yshift)
        if not self.q: print('Updated map info:  %s' %new_map_info)
        new_map_info_line = 'map info = { %s }' %' , '.join(new_map_info)
        new_content = content.replace(map_info_line,new_map_info_line)
        outpath_hdr = '%s.hdr' %os.path.splitext(self.path_out)[0]
        with open(outpath_hdr,'w') as outF:
            outF.write(new_content)

        if not self.q: print('\nCoregistered image written to %s.' %self.path_out)


    def align_coordinate_grids(self):
        xgsd, ygsd = (self.ref.xgsd, self.ref.ygsd) if self.match_gsd else self.out_gsd if self.out_gsd \
                        else (self.shift.xgsd,self.shift.ygsd) # self.match_gsd overrides self.out_gsd in __init__

        if not self.q:
            print('Resampling shift image in order to align coordinate grid of shift image to reference image. ')

        is_alignable = lambda gsd1,gsd2: max(gsd1,gsd2)%min(gsd1,gsd2)==0 # checks if pixel sizes are divisible
        if not is_alignable(self.ref.xgsd,xgsd) or not is_alignable(self.ref.ygsd,ygsd) and not self.q:
            print("\nWARNING: The coordinate grids of the reference image and the image to be shifted cannot be "
                  "aligned because their pixel sizes are not exact multiples of each other (ref [X/Y]: "
                  "%s %s; shift [X/Y]: %s %s). Therefore the pixel size of the reference image is chosen for the "
                  "resampled output image. If you don´t like that you can use the '-out_gsd' parameter to set an "
                  "appropriate output pixel size.\n" %(self.ref.xgsd,self.ref.ygsd,xgsd,ygsd))
            xgsd, ygsd = self.ref.xgsd, self.ref.ygsd

        r_xmin,r_xmax,r_ymin,r_ymax =  GEO.corner_coord_to_minmax(GEO.get_corner_coordinates(self.ref.path))
        imref_xgrid = np.arange(r_xmin,r_xmax,self.ref.xgsd)
        imref_ygrid = np.arange(r_ymin,r_ymax,self.ref.ygsd)
        s_xmin,s_xmax,s_ymin,s_ymax =  GEO.corner_coord_to_minmax(GEO.get_corner_coordinates(self.shift.path))
        xmin,xmax = UTL.find_nearest(imref_xgrid, s_xmin, 'off',extrapolate=1), UTL.find_nearest(imref_xgrid, s_xmax, 'on',extrapolate=1)
        ymin,ymax = UTL.find_nearest(imref_ygrid, s_ymin, 'off',extrapolate=1), UTL.find_nearest(imref_ygrid, s_ymax, 'on',extrapolate=1)

        new_originY, new_originX = GEO.pixelToMapYX([self.x_shift_px, self.y_shift_px],
                                                    geotransform=self.shift.gt, projection=self.shift.prj)[0]
        map_xshift,  map_yshift  = new_originX - self.shift.gt[0], new_originY - self.shift.gt[3]
        if not self.q: print('Calculated map shifts (X,Y): %s, %s' %(map_xshift,map_yshift))
        gt_shifted                   = list(self.shift.gt[:])
        gt_shifted[0], gt_shifted[3] = new_originX, new_originY

        pathVRT = os.path.splitext(os.path.basename(self.shift.path))[0] + '.vrt'
        ds_im2shift = gdal.Open(self.shift.path)
        dst_ds  = gdal.GetDriverByName('VRT').CreateCopy(pathVRT, ds_im2shift, 0)
        dst_ds.SetGeoTransform(gt_shifted)
        dst_ds = ds_im2shift  = None

        cmd = "gdalwarp -r cubic -tr %s %s -t_srs '%s' -of ENVI %s %s -overwrite -te %s %s %s %s" \
              %(xgsd,ygsd,self.ref.prj,pathVRT,self.path_out, xmin,ymin,xmax,ymax)
        output = subprocess.check_output(cmd, shell=True)
        if output!=1 and os.path.exists(self.path_out):
            if not self.q: print('Coregistered and resampled image written to %s.' %self.path_out)
        else:
            print(output)
            self.tracked_errors.append(RuntimeError('Resampling failed.'))
            raise self.tracked_errors[-1]


    def resample_without_grid_aligning(self):
        xgsd, ygsd = (self.ref.xgsd, self.ref.ygsd) if self.match_gsd else self.out_gsd if self.out_gsd \
                        else (self.shift.xgsd,self.shift.ygsd) # self.match_gsd overrides self.out_gsd in __init__

        new_originY, new_originX = GEO.pixelToMapYX([self.x_shift_px, self.y_shift_px],
                                                    geotransform=self.shift.gt, projection=self.shift.prj)[0]
        map_xshift,  map_yshift  = new_originX - self.shift.gt[0], new_originY - self.shift.gt[3]
        if not self.q: print('Calculated map shifts (X,Y): %s, %s' %(map_xshift,map_yshift))
        gt_shifted                   = list(self.shift.gt[:])
        gt_shifted[0], gt_shifted[3] = new_originX, new_originY

        pathVRT     = os.path.splitext(os.path.basename(self.shift.path))[0] + '.vrt'
        ds_im2shift = gdal.Open(self.shift.path)
        dst_ds  = (lambda drv: drv.CreateCopy(pathVRT, ds_im2shift, 0))(gdal.GetDriverByName('VRT'))
        dst_ds.SetGeoTransform(gt_shifted)
        dst_ds = ds_im2shift = None

        if not self.q: print('Resampling shift image without aligning of coordinate grids. ')
        cmd = "gdalwarp -r cubic -tr %s %s -t_srs '%s' -of ENVI %s %s -overwrite" \
              %(xgsd,ygsd,self.ref.prj,pathVRT,self.path_out)
        output = subprocess.check_output(cmd, shell=True)
        if output!=1 and os.path.exists(self.path_out):
            if not self.q: print('Coregistered and resampled image written to %s.' %self.path_out)
        else:
            print(output)
            self.tracked_errors.append(RuntimeError('Resampling failed.'))
            raise self.tracked_errors[-1]


class DESHIFTER(object):
    dict_rspAlg_rsp_Int = {'nearest': 0, 'bilinear': 1, 'cubic': 2,
                           'cubic_spline': 3, 'lanczos': 4, 'average': 5, 'mode': 6}

    def __init__(self, im2shift, coreg_results, **kwargs):
        """
        Deshift an image array or one of its products by applying the coregistration info calculated by COREG class.

        :param dict_GMS_obj:        the copied dictionary of a GMS object, containing the attribute 'coreg_info'
        :param attrname2deshift:    attribute name of the GMS object containing the array to be shifted (or a its path)

        :Keyword Arguments:
            - band2process (int):   The index of the band to be processed within the given array (starts with 1),
                                    default = None (all bands are processed)
            - out_gsd (float):      output pixel size in units of the reference coordinate system (default = pixel size
                                    of the input array), given values are overridden by match_gsd=True
            - align_grids (bool):   True: align the input coordinate grid to the reference (does not affect the
                                    output pixel size as long as input and output pixel sizes are compatible
                                    (5:30 or 10:30 but not 4:30), default = False
            - match_gsd (bool):     True: match the input pixel size to the reference pixel size,
                                    default = False
            - target_xyGrid(list):  a list with an x-grid and a y-grid like [[15,45], [15,45]]
            - resamp_alg(str)       the resampling algorithm to be used if neccessary
                                    (valid algorithms: nearest, bilinear, cubic, cubic_spline, lanczos, average, mode)
            - warp_alg(str):        'GDAL' or 'rasterio' (default = 'rasterio')
            - cliptoextent (bool):  True: clip the input image to its actual bounds while deleting possible no data
                                    areas outside of the actual bounds, default = True
            - clipextent (list):    xmin, ymin, xmax, ymax - if given the calculation of the actual bounds is skipped
            - tempDir(str):         directory to be used for tempfiles (default: /dev/shm/)

        """
        # FIXME add v and q, mp?
        # unpack args
        self.im2shift           = im2shift if isinstance(im2shift, GeoArray) else GeoArray(im2shift)
        self.shift_prj          = im2shift.projection
        self.shift_gt           = list(im2shift.geotransform)
        self.nodata             = GEO.get_outFillZeroSaturated(self.im2shift.dtype)[0]
        self.updated_map_info   = coreg_results['updated map info']
        self.original_map_info  = coreg_results['original map info']
        self.updated_gt         = GEO.mapinfo2geotransform(self.updated_map_info)
        self.ref_gt             = coreg_results['reference geotransform']
        self.ref_grid           = coreg_results['reference grid']
        self.ref_prj            = coreg_results['reference projection']
        self.updated_projection = self.ref_prj

        # unpack kwargs
        self.band2process = kwargs.get('band2process', None) # starts with 1 # FIXME warum?
        tempAsENVI        = kwargs.get('tempAsENVI'  , False)
        self.outFmt       = 'VRT' if not tempAsENVI else 'ENVI'
        self.rspAlg       = kwargs.get('resamp_alg'  , 'cubic')
        self.warpAlg      = kwargs.get('warp_alg'    , 'GDAL_lib')
        self.cliptoextent = kwargs.get('cliptoextent', True)
        self.clipextent   = kwargs.get('clipextent'  , None)
        self.tempDir      = kwargs.get('tempDir'     ,'/dev/shm/')
        self.out_grid     = self.get_out_grid(kwargs) # needs self.ref_grid, self.im2shift
        self.out_gsd      = [abs(self.out_grid[0][1]-self.out_grid[0][0]), abs(self.out_grid[1][1]-self.out_grid[1][0])]  # xgsd, ygsd

        # assertions
        assert self.rspAlg  in self.dict_rspAlg_rsp_Int.keys()
        assert self.warpAlg in ['GDAL_cmd', 'GDAL_lib']

        # set defaults for general class attributes
        self.is_shifted       = False # this is not included in COREG.coreg_info
        self.is_resampled     = False # this is not included in COREG.coreg_info
        self.tracked_errors   = []
        self.arr_shifted      = None  # set by self.correct_shifts


    def get_out_grid(self, init_kwargs):
        # parse given params
        out_gsd     = init_kwargs.get('out_gsd'      , None)
        align_grids = init_kwargs.get('align_grids'  , False)
        match_gsd   = init_kwargs.get('match_gsd'    , False)
        out_grid    = init_kwargs.get('target_xyGrid', None)

        # assertions
        assert out_grid is None or (isinstance(out_grid,(list, tuple)) and len(out_grid)==2)
        assert out_gsd  is None or (isinstance(out_gsd, (int, list))   and len(out_gsd) ==2)

        ref_xgsd, ref_ygsd = (self.ref_grid[0][1]-self.ref_grid[0][0],self.ref_grid[1][1]-self.ref_grid[1][0])
        get_grid           = lambda gt, xgsd, ygsd: [[gt[0], gt[0] + xgsd], [gt[3], gt[3] - ygsd]]

        # get out_grid
        if out_grid:
            # output grid is given
            return out_grid

        elif out_gsd:
            out_xgsd, out_ygsd = [out_gsd, out_gsd] if isinstance(out_gsd, int) else out_gsd

            if match_gsd and (out_xgsd, out_ygsd)!=(ref_xgsd, ref_ygsd):
                warnings.warn("\nThe parameter 'match_gsd is ignored because another output ground sampling distance "
                              "was explicitly given.")
            if align_grids and self.grids_alignable(self.im2shift.xgsd, self.im2shift.ygsd, out_xgsd, out_ygsd):
                # use grid of reference image with the given output gsd
                return get_grid(self.ref_gt, out_xgsd, out_ygsd)
            else: # no grid alignment
                # use grid of input image with the given output gsd
                return get_grid(self.im2shift.geotransform, out_xgsd, out_ygsd)

        elif match_gsd:
            if align_grids:
                # use reference grid
                return self.ref_grid
            else:
                # use grid of input image and reference gsd
                return get_grid(self.im2shift.geotransform, ref_xgsd, ref_ygsd)

        else:
            if align_grids and self.grids_alignable(self.im2shift.xgsd, self.im2shift.ygsd, ref_xgsd, ref_ygsd):
                # use origin of reference image and gsd of input image
                return get_grid(self.ref_gt, self.im2shift.xgsd, self.im2shift.ygsd)
            else:
                # use input image grid
                return get_grid(self.im2shift.geotransform, self.im2shift.xgsd, self.im2shift.ygsd)


    @staticmethod
    def grids_alignable(in_xgsd, in_ygsd, out_xgsd, out_ygsd):
        is_alignable = lambda gsd1, gsd2: max(gsd1, gsd2) % min(gsd1, gsd2) == 0  # checks if pixel sizes are divisible
        if not is_alignable(in_xgsd, out_xgsd) or not is_alignable(in_ygsd, out_ygsd):
            warnings.warn("\nThe targeted output coordinate grid is not alignable with the image to be shifted because "
                          "their pixel sizes are not exact multiples of each other (input [X/Y]: "
                          "%s %s; output [X/Y]: %s %s). Therefore the targeted output grid is "
                          "chosen for the resampled output image. If you don´t like that you can use the '-out_gsd' "
                          "parameter to set an appropriate output pixel size.\n"
                          % (in_xgsd, in_ygsd, out_xgsd, out_ygsd))
            return False
        else:
            return True


    def get_out_extent(self):
        if self.cliptoextent and self.clipextent is None:
            # calculate actual corner coords (input image projection)
            trueCorner_mapXY       = GEO.get_true_corner_mapXY(self.im2shift, bandNr=1, noDataVal=self.nodata, mp=1,v=0)
            xmin, xmax, ymin, ymax = GEO.corner_coord_to_minmax(trueCorner_mapXY)
            self.clipextent        = box(xmin, ymin, xmax, ymax).bounds
        else:
            self.clipextent = GEO.corner_coord_to_minmax(GEO.get_corner_coordinates(
                                    gt=self.shift_gt, cols=self.im2shift.cols, rows=self.im2shift.rows))

        # snap clipextent to output grid (in case of odd input coords the output coords are moved INSIDE the input array)
        xmin, ymin, xmax, ymax = self.clipextent
        xmin = UTL.find_nearest(self.out_grid[0], xmin, roundAlg='on' , extrapolate=True)
        ymin = UTL.find_nearest(self.out_grid[1], ymin, roundAlg='on' , extrapolate=True)
        xmax = UTL.find_nearest(self.out_grid[0], xmax, roundAlg='off', extrapolate=True)
        ymax = UTL.find_nearest(self.out_grid[0], ymax, roundAlg='off', extrapolate=True)
        return xmin, ymin, xmax, ymax


    def correct_shifts(self):
        # type: (L2A_P.DESHIFTER) -> collections.OrderedDict

        t_start = time.time()
        equal_prj = GEO.prj_equal(self.ref_prj,self.shift_prj)

        if equal_prj and GEO.is_coord_grid_equal(self.shift_gt, *self.out_grid):
            """NO RESAMPLING NEEDED"""
            self.is_shifted     = True
            self.is_resampled   = False
            xmin,ymin,xmax,ymax = self.get_out_extent()

            if self.cliptoextent:
                # clip with target extent
                self.arr_shifted, self.updated_gt, self.updated_projection = \
                    self.im2shift.get_mapPos((xmin,ymin,xmax,ymax), self.shift_prj, fillVal=self.nodata)
                self.updated_map_info = GEO.geotransform2mapinfo(self.updated_gt, self.updated_projection)
            else:
                self.arr_shifted = self.im2shift[:]

            # apply shifts by updating gt
            self.updated_projection = self.updated_projection if self.updated_projection else self.shift_prj  # self.updated_projection is taken from ref_prj which is None if target image don't has to be coregistered
            self.updated_gt[3], self.updated_gt[0] = \
                GEO.pixelToMapYX([xmin, ymin], geotransform=self.updated_gt, projection=self.updated_projection)[0]
            self.updated_map_info = GEO.geotransform2mapinfo(self.updated_gt, self.updated_projection)

        else: # FIXME equal_prj==False ist noch NICHT implementiert
            """RESAMPLING NEEDED"""
            if self.warpAlg=='GDAL_cmd':
                warnings.warn('This method has not been tested in its current state!')
                # FIXME nicht multiprocessing-fähig, weil immer kompletter array gewarpt wird und sich ergebnisse gegenseitig überschreiben
                # create tempfile
                fd, path_tmp = tempfile.mkstemp(prefix='CoReg_Sat', suffix=self.outFmt, dir=self.tempDir)
                os.close(fd)

                t_extent   = " -te %s %s %s %s" %self.get_out_extent()
                xgsd, ygsd = self.out_gsd
                cmd = "gdalwarp -r %s -tr %s %s -t_srs '%s' -of %s %s %s -srcnodata %s -dstnodata %s -overwrite%s"\
                      %(self.rspAlg, xgsd,ygsd,self.ref_prj,self.outFmt,self.im2shift.filePath,
                        path_tmp, self.nodata, self.nodata, t_extent)
                out, exitcode, err = UTL.subcall_with_output(cmd)

                if exitcode!=1 and os.path.exists(path_tmp):
                    """update map info, arr_shifted, geotransform and projection"""
                    ds_shifted = gdal.OpenShared(path_tmp) if self.outFmt == 'VRT' else gdal.Open(path_tmp)
                    self.shift_gt, self.shift_prj = ds_shifted.GetGeoTransform(), ds_shifted.GetProjection()
                    self.updated_map_info         = GEO.geotransform2mapinfo(self.shift_gt,self.shift_prj)

                    print('reading from', ds_shifted.GetDescription())
                    if self.band2process is None:
                        dim2RowsColsBands = lambda A: np.swapaxes(np.swapaxes(A,0,2),0,1) # rasterio.open(): [bands,rows,cols]
                        self.arr_shifted  = dim2RowsColsBands(rasterio.open(path_tmp).read())
                    else:
                        self.arr_shifted  = rasterio.open(path_tmp).read(self.band2process)

                    self.is_shifted       = True
                    self.is_resampled     = True

                    ds_shifted            = None
                    [gdal.Unlink(p) for p in [path_tmp] if os.path.exists(p)] # delete tempfiles
                else:
                    print("\n%s\nCommand was:  '%s'" %(err.decode('utf8'),cmd))
                    [gdal.Unlink(p) for p in [path_tmp] if os.path.exists(p)] # delete tempfiles
                    self.tracked_errors.append(RuntimeError('Resampling failed.'))
                    raise self.tracked_errors[-1]

            elif self.warpAlg=='GDAL_lib':
                # apply XY-shifts to shift_gt
                in_arr = self.im2shift[self.band2process] if self.band2process else self.im2shift[:]
                self.shift_gt[0], self.shift_gt[3] = self.updated_gt[0], self.updated_gt[3]

                # get resampled array
                out_arr, out_gt, out_prj = \
                    GEO.warp_ndarray(in_arr,self.shift_gt,self.shift_prj,self.ref_prj,
                                     rspAlg            = self.dict_rspAlg_rsp_Int[self.rspAlg],
                                     in_nodata          = self.nodata,
                                     out_nodata         = self.nodata,
                                     out_gsd            = self.out_gsd,
                                     out_bounds         = self.get_out_extent(),
)

                self.updated_projection = self.ref_prj
                self.arr_shifted        = out_arr
                self.updated_map_info   = GEO.geotransform2mapinfo(out_gt,out_prj)
                self.shift_gt           = GEO.mapinfo2geotransform(self.updated_map_info)
                self.is_shifted         = True
                self.is_resampled       = True

        print('Time for shift correction: %.2fs' %(time.time()-t_start))
        return self.deshift_results


    @property
    def deshift_results(self):
        deshift_results = collections.OrderedDict()
        deshift_results.update({'band'              :self.band2process})
        deshift_results.update({'is shifted'        :self.is_shifted})
        deshift_results.update({'is resampled'      :self.is_resampled})
        deshift_results.update({'updated map info'  :self.updated_map_info})
        deshift_results.update({'updated projection':self.updated_projection})
        deshift_results.update({'arr_shifted'       :self.arr_shifted})
        return deshift_results



class Geom_Quality_Grid(object):
    def __init__(self, path_im0, path_im1, grid_res, window_size=256, dir_out=None, projectName=None, multiproc=True,
                 r_b4match=1, s_b4match=1, max_iter=5, max_shift=5, data_corners_im0=None,
                 data_corners_im1=None, outFillVal=-9999, nodata=None, calc_corners=True, binary_ws=True,
                 v=False, q=False, ignore_errors=False):
        self.path_imref    = path_im0
        self.path_im2shift = path_im1
        self.dir_out       = dir_out
        self.grid_res      = grid_res
        self.window_size   = window_size
        self.mp            = multiproc
        self.max_shift     = max_shift
        self.max_iter      = max_iter
        self.r_b4match     = r_b4match
        self.s_b4match     = s_b4match
        self.calc_corners  = calc_corners
        self.nodata        = nodata
        self.outFillVal    = outFillVal
        self.bin_ws        = binary_ws
        self.v             = v
        self.q             = q
        self.ignore_errors = ignore_errors

        self.projectName   = projectName if projectName else 'UntitledProject_1'
        while projectName is None and os.path.isdir(os.path.join(self.dir_out,self.projectName)):
            self.projectName = '%s_%s' %(self.projectName.split('_')[0],int(self.projectName.split('_')[1])+1)
        self.dir_out = os.path.join(self.dir_out,self.projectName)
        if not os.path.exists(self.dir_out): os.makedirs(self.dir_out)

        gdal.AllRegister()

        self.COREG_obj = COREG(path_im0, path_im1, data_corners_im0=data_corners_im0,
                               data_corners_im1=data_corners_im1, calc_corners=calc_corners, r_b4match=r_b4match,
                               s_b4match=s_b4match, max_iter=max_iter, max_shift=max_shift, nodata=nodata,
                               multiproc=multiproc, binary_ws=self.bin_ws, v=v, q=q, ignore_errors=ignore_errors)
        self.ref_shape                = [self.COREG_obj.ref.rows,self.COREG_obj.ref.cols]
        self.tgt_shape                = [self.COREG_obj.shift.rows,self.COREG_obj.shift.cols]
        self.corner_coord_imref       = self.COREG_obj.ref.corner_coord
        self.corner_coord_im2shift    = self.COREG_obj.shift.corner_coord
        self.overlap_poly             = self.COREG_obj.overlap_poly

        self.XY_points, self.XY_mapPoints = self._get_imXY__mapXY_points(self.grid_res)
        self.quality_grid                 = None # set by self.get_quality_grid()

    def _get_imXY__mapXY_points(self,grid_res):
        Xarr,Yarr       = np.meshgrid(np.arange(0,self.tgt_shape[1],grid_res),
                                      np.arange(0,self.tgt_shape[0],grid_res))
        ULmapYX,LRmapYX = GEO.pixelToMapYX([[0,0],[self.tgt_shape[1],self.tgt_shape[0]]], path_im=self.path_im2shift)

        mapXarr,mapYarr = np.meshgrid(np.arange(ULmapYX[1],LRmapYX[1],     self.grid_res*self.COREG_obj.shift.xgsd),
                                      np.arange(ULmapYX[0],LRmapYX[0],-abs(self.grid_res*self.COREG_obj.shift.ygsd)))

        XY_points      = np.empty((Xarr.size,2),Xarr.dtype)
        XY_points[:,0] = Xarr.flat
        XY_points[:,1] = Yarr.flat

        XY_mapPoints      = np.empty((mapXarr.size,2),mapXarr.dtype)
        XY_mapPoints[:,0] = mapXarr.flat
        XY_mapPoints[:,1] = mapYarr.flat

        return XY_points,XY_mapPoints

    def _get_spatial_shifts(self, args):
        pointID,winpos_XY,cornerCoords_l8,cornerCoords_tgt = args
        COREG_obj = COREG(self.path_imref, self.path_im2shift, wp=winpos_XY, ws=self.window_size,
                          data_corners_im0=cornerCoords_l8, data_corners_im1=cornerCoords_tgt, r_b4match=self.r_b4match,
                          s_b4match=self.s_b4match, max_iter=self.max_iter, max_shift=self.max_shift,
                          nodata=self.nodata, binary_ws=self.bin_ws, v=self.v, q=self.q, ignore_errors=self.ignore_errors)
        COREG_obj.calculate_spatial_shifts()
        res = pointID, COREG_obj.ref.win.size_YX[0], COREG_obj.x_shift_px, COREG_obj.y_shift_px
        return res

    def get_quality_grid(self,exclude_outliers=1,dump_values=1):
        assert self.XY_points is not None and self.XY_mapPoints is not None

        #ref_ds,tgt_ds = gdal.Open(self.path_imref),gdal.Open(self.path_im2shift)
        #ref_pathTmp, tgt_pathTmp = None,None
        #if ref_ds.GetDriver().ShortName!='ENVI':
        #    ref_pathTmp = IO.get_tempfile(ext='.bsq')
        #    IO.convert_gdal_to_bsq__mp(self.path_imref,ref_pathTmp)
        #    self.path_imref = ref_pathTmp
        #if tgt_ds.GetDriver().ShortName!='ENVI':
        #    tgt_pathTmp = IO.get_tempfile(ext='.bsq')
        #    IO.convert_gdal_to_bsq__mp(self.path_im2shift,tgt_pathTmp)
        #    self.path_im2shift = tgt_pathTmp
        #ref_ds=tgt_ds=None

        XYarr2PointGeom = np.vectorize(lambda X,Y: Point(X,Y), otypes=[Point])
        geomPoints      = np.array(XYarr2PointGeom(self.XY_mapPoints[:,0],self.XY_mapPoints[:,1]))

        if GEO.isProjectedOrGeographic(self.COREG_obj.shift.prj)=='geographic':
            crs = dict(ellps='WGS84', datum='WGS84', proj='longlat')
        elif GEO.isProjectedOrGeographic(self.COREG_obj.shift.prj)=='projected':
            UTMzone = abs(GEO.get_UTMzone(prj=self.COREG_obj.shift.prj))
            south   = GEO.get_UTMzone(prj=self.COREG_obj.shift.prj)<0
            crs     = dict(ellps='WGS84', datum='WGS84', proj='utm', zone=UTMzone,south=south,units='m', no_defs=True)
            if not south: del crs['south']
        else:
            crs = None

        GDF        = GeoDataFrame(index=range(len(geomPoints)),crs=crs,
                                  columns=['geometry','POINT_ID','X_IM','Y_IM','X_UTM','Y_UTM','WIN_SIZE',
                                           'X_SHIFT_PX','Y_SHIFT_PX','X_SHIFT_M','Y_SHIFT_M','ABS_SHIFT','ANGLE'])
        GDF       ['geometry']       = geomPoints
        GDF       ['POINT_ID']       = range(len(geomPoints))
        GDF.loc[:,['X_IM','Y_IM']]   = self.XY_points
        GDF.loc[:,['X_UTM','Y_UTM']] = self.XY_mapPoints
        GDF = GDF if not exclude_outliers else GDF[GDF['geometry'].within(self.overlap_poly)]
        GDF.loc[:,['WIN_SIZE','X_SHIFT_PX','Y_SHIFT_PX','X_SHIFT_M','Y_SHIFT_M','ABS_SHIFT','ANGLE']] = self.outFillVal # Fehlwert

        args = [[i,self.XY_mapPoints[i],self.corner_coord_imref,self.corner_coord_im2shift] for i in GDF.index]#[:1000]

        if self.mp:
            print('multiprocessing')
            with multiprocessing.Pool() as pool:
                results = pool.map(self._get_spatial_shifts, args)
        else:
            print('singleprocessing')
            results = []
            for i,argset in enumerate(args):
                #print(argset[1])
                #if not 0<i<10: continue
                #if i>300 or i<100: continue
                #if i!=127: continue
                if i%100==0: print('Point #%s, ID %s' %(i,argset[0]))
                res = self._get_spatial_shifts(argset)
                results.append(res)

        #if ref_pathTmp and os.path.exists(ref_pathTmp): os.remove(ref_pathTmp)
        #if tgt_pathTmp and os.path.exists(tgt_pathTmp): os.remove(tgt_pathTmp)

        for res in results:
            pointID                      = res[0]
            GDF.loc[pointID,'WIN_SIZE']   = res[1] if res[1] is not None else self.outFillVal
            GDF.loc[pointID,'X_SHIFT_PX'] = res[2] if res[2] is not None else self.outFillVal
            GDF.loc[pointID,'Y_SHIFT_PX'] = res[3] if res[3] is not None else self.outFillVal

        s2_gsd = self.COREG_obj.shift.xgsd
        GDF.loc[:,'X_SHIFT_M']  = np.where(GDF['X_SHIFT_PX']==self.outFillVal,self.outFillVal,GDF['X_SHIFT_PX']*s2_gsd)
        GDF.loc[:,'Y_SHIFT_M']  = np.where(GDF['Y_SHIFT_PX']==self.outFillVal,self.outFillVal,GDF['Y_SHIFT_PX']*s2_gsd)
        GDF.loc[:,'ABS_SHIFT']  = np.where(GDF['X_SHIFT_M'] ==self.outFillVal,self.outFillVal,
                                          np.sqrt((GDF['X_SHIFT_M']**2 + GDF['Y_SHIFT_M']**2).astype(np.float64)))
        GDF.loc[:,'ANGLE']      = np.where(GDF['X_SHIFT_PX']==self.outFillVal,self.outFillVal,
                                          GEO.angle_to_north(GDF.loc[:,['X_SHIFT_PX','Y_SHIFT_PX']]))

        if dump_values:
            fName_out = "CoRegMatrix_grid%s_ws%s__T_%s__R_%s.pkl" %(self.grid_res, self.window_size, os.path.splitext(
                os.path.basename(self.path_im2shift))[0], os.path.splitext(os.path.basename(self.path_imref))[0])
            path_out  = os.path.join(self.dir_out, 'CoRegMatrix', fName_out)
            if not os.path.exists(os.path.dirname(path_out)): os.makedirs(os.path.dirname(path_out))
            GDF.to_pickle(path_out)

        self.quality_grid = GDF
        return GDF

    def test_if_singleprocessing_equals_multiprocessing_result(self):
        self.mp = 1
        dataframe = self.get_quality_grid()
        mp_out    = np.empty_like(dataframe.values)
        mp_out[:] = dataframe.values
        self.mp = 0
        dataframe = self.get_quality_grid()
        sp_out    = np.empty_like(dataframe.values)
        sp_out[:] = dataframe.values

        return np.array_equal(sp_out,mp_out)

    def get_line_by_PID(self,PID):
        assert self.quality_grid, 'Calculate quality grid first!'
        return self.quality_grid.loc[PID,:]

    def get_lines_by_PIDs(self,PIDs):
        assert self.quality_grid, 'Calculate quality grid first!'
        assert isinstance(PIDs,list)
        lines = np.zeros((len(PIDs),self.quality_grid.shape[1]))
        for i,PID in enumerate(PIDs):
            lines[i,:] = self.quality_grid[self.quality_grid['POINT_ID']==PID]
        return lines

    def quality_grid_to_PointShapefile(self,skip_nodata=1,skip_nodata_col = 'ABS_SHIFT'):
        GDF            = self.quality_grid
        GDF2pass       = GDF if not skip_nodata else GDF[GDF[skip_nodata_col]!=self.outFillVal]

        fName_out = "CoRegPoints_grid%s_ws%s__T_%s__R_%s.shp" %(self.grid_res, self.window_size, os.path.splitext(
            os.path.basename(self.path_im2shift))[0], os.path.splitext(os.path.basename(self.path_imref))[0])
        path_out  = os.path.join(self.dir_out, 'CoRegPoints', fName_out)
        if not os.path.exists(os.path.dirname(path_out)): os.makedirs(os.path.dirname(path_out))
        print('Writing %s ...' %path_out)
        GDF2pass.to_file(path_out)

    def _quality_grid_to_PointShapefile(self,skip_nodata=1,skip_nodata_col = 'ABS_SHIFT'):
        warnings.warn(DeprecationWarning("'_quality_grid_to_PointShapefile' deprecated."
                                         " 'quality_grid_to_PointShapefile' is much faster."))
        GDF            = self.quality_grid
        GDF2pass       = GDF if not skip_nodata else GDF[GDF[skip_nodata_col]!=self.outFillVal]
        shapely_points = GDF2pass['geometry'].values.tolist()
        attr_dicts     = [collections.OrderedDict(zip(GDF2pass.columns,GDF2pass.loc[i].values)) for i in GDF2pass.index]


        fName_out = "CoRegPoints_grid%s_ws%s.shp" %(self.grid_res, self.window_size)
        path_out = os.path.join(self.dir_out, fName_out)
        IO.write_shp(path_out, shapely_points, prj=self.COREG_obj.shift.prj, attrDict=attr_dicts)

    def quality_grid_to_Raster_using_KrigingOLD(self,attrName,skip_nodata=1,skip_nodata_col='ABS_SHIFT',outGridRes=None,
                                             fName_out=None,tilepos=None):
        GDF             = self.quality_grid
        GDF2pass        = GDF if not skip_nodata else GDF[GDF[skip_nodata_col]!=self.outFillVal]

        # subset if tilepos is given
        rows,cols = tilepos if tilepos else self.tgt_shape
        GDF2pass        = GDF2pass.loc[(GDF2pass['X_IM']>=cols[0])&(GDF2pass['X_IM']<=cols[1])&
                                       (GDF2pass['Y_IM']>=rows[0])&(GDF2pass['Y_IM']<=rows[1])]


        X_coords,Y_coords,ABS_SHIFT = GDF2pass['X_UTM'], GDF2pass['Y_UTM'],GDF2pass[attrName]

        xmin,ymin,xmax,ymax = GDF2pass.total_bounds

        grid_res            = outGridRes if outGridRes else int(min(xmax-xmin,ymax-ymin)/250)
        grid_x,grid_y       = np.arange(xmin, xmax+grid_res, grid_res), np.arange(ymax, ymin-grid_res, -grid_res)

        # Reference: P.K. Kitanidis, Introduction to Geostatistcs: Applications in Hydrogeology,
        #            (Cambridge University Press, 1997) 272 p.
        OK = OrdinaryKriging(X_coords, Y_coords, ABS_SHIFT, variogram_model='spherical',verbose=False)
        zvalues, sigmasq = OK.execute('grid', grid_x, grid_y)#,backend='C',)

        fName_out = fName_out if fName_out else \
            "Kriging__%s__grid%s_ws%s.tif" %(attrName,self.grid_res, self.window_size)
        path_out  = os.path.join(self.dir_out, fName_out)
        print('Writing %s ...' %path_out)
        # add a half pixel grid points are centered on the output pixels
        xmin,ymin,xmax,ymax = xmin-grid_res/2,ymin-grid_res/2,xmax+grid_res/2,ymax+grid_res/2
        IO.write_numpy_to_image(zvalues, path_out, gt=(xmin, grid_res, 0, ymax, 0, -grid_res), prj=self.COREG_obj.shift.prj)

        return zvalues


    def quality_grid_to_Raster_using_Kriging(self,attrName,skip_nodata=1,skip_nodata_col='ABS_SHIFT',outGridRes=None,
                                             fName_out=None,tilepos=None,tilesize=500,mp=None):

        mp = mp if mp else self.mp
        self.Kriging_sp(attrName,skip_nodata=skip_nodata,skip_nodata_col=skip_nodata_col,
                outGridRes=outGridRes,fName_out=fName_out,tilepos=tilepos)

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

    def Kriging_sp(self,attrName,skip_nodata=1,skip_nodata_col='ABS_SHIFT',outGridRes=None,
                                             fName_out=None,tilepos=None):
        GDF             = self.quality_grid
        GDF2pass        = GDF if not skip_nodata else GDF[GDF[skip_nodata_col]!=self.outFillVal]

#         # subset if tilepos is given
# #        overlap_factor =
#         rows,cols = tilepos if tilepos else self.tgt_shape
#         xvals, yvals = np.sort(GDF2pass['X_IM'].values.flat),np.sort(GDF2pass['Y_IM'].values.flat)
#         cS,cE = UTL.find_nearest(xvals,cols[0],'off',1), UTL.find_nearest(xvals,cols[1],'on',1)
#         rS,rE = UTL.find_nearest(yvals,rows[0],'off',1), UTL.find_nearest(yvals,rows[1],'on',1)
#         # GDF2pass        = GDF2pass.loc[(GDF2pass['X_IM']>=cols[0])&(GDF2pass['X_IM']<=cols[1])&
#         #                                (GDF2pass['Y_IM']>=rows[0])&(GDF2pass['Y_IM']<=rows[1])]
#         GDF2pass        = GDF2pass.loc[(GDF2pass['X_IM']>=cS)&(GDF2pass['X_IM']<=cE)&
#                                        (GDF2pass['Y_IM']>=rS)&(GDF2pass['Y_IM']<=rE)]

        X_coords,Y_coords,ABS_SHIFT = GDF2pass['X_UTM'], GDF2pass['Y_UTM'],GDF2pass[attrName]

        xmin,ymin,xmax,ymax = GDF2pass.total_bounds

        grid_res            = outGridRes if outGridRes else int(min(xmax-xmin,ymax-ymin)/250)
        grid_x,grid_y       = np.arange(xmin, xmax+grid_res, grid_res), np.arange(ymax, ymin-grid_res, -grid_res)

        # Reference: P.K. Kitanidis, Introduction to Geostatistcs: Applications in Hydrogeology,
        #            (Cambridge University Press, 1997) 272 p.
        OK = OrdinaryKriging(X_coords, Y_coords, ABS_SHIFT, variogram_model='spherical',verbose=False)
        zvalues, sigmasq = OK.execute('grid', grid_x, grid_y,backend='C',n_closest_points=12)

        if self.mp:
            fName_out = fName_out if fName_out else \
                "Kriging__%s__grid%s_ws%s_%s.tif" %(attrName,self.grid_res, self.window_size,tilepos)
        else:
            fName_out = fName_out if fName_out else \
                "Kriging__%s__grid%s_ws%s.tif" %(attrName,self.grid_res, self.window_size)
        path_out  = os.path.join(self.dir_out, fName_out)
        print('Writing %s ...' %path_out)
        # add a half pixel grid points are centered on the output pixels
        xmin,ymin,xmax,ymax = xmin-grid_res/2,ymin-grid_res/2,xmax+grid_res/2,ymax+grid_res/2
        IO.write_numpy_to_image(zvalues, path_out, gt=(xmin, grid_res, 0, ymax, 0, -grid_res), prj=self.COREG_obj.shift.prj)

        return zvalues

    def Kriging_mp(self,args_kwargs_dict):
        args   = args_kwargs_dict.get('args',[])
        kwargs = args_kwargs_dict.get('kwargs',[])

        return self.Kriging_sp(*args,**kwargs)




if __name__ == '__main__':
    import argparse
    from socket import gethostname
    from datetime import datetime as dt
    from getpass import getuser
    from components.io import wfa
    wfa('/misc/hy5/scheffler/tmp/crlf', '%s\t%s\t%s\t%s\n' % (dt.now(), getuser(), gethostname(), ' '.join(sys.argv)))
    parser = argparse.ArgumentParser(prog='dsc__CoReg_Sat_FourierShiftTheorem.py',
        description='Perform subpixel coregistration of two satellite image datasets ' \
        'using Fourier Shift Theorem proposed by Foroosh et al. 2002: ' \
        'Foroosh, H., Zerubia, J. B., & Berthod, M. (2002). Extension of phase correlation to subpixel registration. '\
        'IEEE Transactions on Image Processing, 11(3), 188–199. doi:10.1109/83.988953); Python implementation by '\
        'Daniel Scheffler (daniel.scheffler@gfz-potsdam.de)',
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
    parser.add_argument('path_im0', type=str, help='source path of reference image (any GDAL compatible image format is supported)')
    parser.add_argument('path_im1', type=str, help='source path of image to be shifted (any GDAL compatible image format is supported)')
    parser.add_argument('-o', nargs='?', dest='path_out', type=str,
                        help="target path of the coregistered image (default: /dir/of/im1/<im1>__shifted_to__<im0>.bsq)", default='.')
    parser.add_argument('-br', nargs='?',type=int, help='band of reference image to be used for matching '\
                        '(starts with 1; default: 1)', default=1)
    parser.add_argument('-bs', nargs='?',type=int, help='band of shift image to be used for matching '\
                        '(starts with 1; default: 1)', default=1)
    parser.add_argument('-wp', nargs=2, metavar=('X', 'Y'),type=float,help="custom matching window position as map "\
                        "values in the same projection like the reference image "\
                        "(default: central position of image overlap)", default=(None,None))
    parser.add_argument('-ws', nargs=2, metavar=('X size', 'Y size'),type=float,
                        help="custom matching window size [pixels] (default: (512,512))", default=(512,512))
    parser.add_argument('-max_iter', nargs='?', type=int,help="maximum number of iterations for matching (default: 5)", default=5)
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
    parser.add_argument('-calc_cor', nargs=1,type=int, choices=[0,1],default=1, help="calculate true positions of "\
                        "the dataset corners in order to get a useful matching window position within the actual "\
                        "image overlap (default: 1; deactivated if '-cor0' and '-cor1' are given")
    parser.add_argument('-mp', nargs='?', type=int, help='enable multiprocessing (default: 1)', default=1, choices=[0,1])
    parser.add_argument('-bin_ws', nargs='?', type=int, help='use binary X/Y dimensions for the matching window '
                                                             '(default: 1)', default=1, choices=[0, 1])
    parser.add_argument('-quadratic_win', nargs='?', type=int, help='force a quadratic matching window (default: 1)',
                        default=1, choices=[0, 1])
    parser.add_argument('-v', nargs='?',type=int, help='verbose mode (default: 0)', default=0, choices=[0,1])
    parser.add_argument('-q', nargs='?',type=int, help='quiet mode (default: 0)', default=0, choices=[0,1])
    parser.add_argument('-ignore_errors', nargs='?',type=int, help='Useful for batch processing. (default: 0) '
                        'In case of error COREG.success == False and COREG.x_shift_px/COREG.y_shift_px is None',
                        default=0, choices=[0,1])
    parser.add_argument('--version', action='version', version='%(prog)s 2016-09-26_01')
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
                      q                = args.q,
                      ignore_errors    = args.ignore_errors)
    COREG_obj.calculate_spatial_shifts()
    COREG_obj.correct_shifts()
    print('\ntotal processing time: %.2fs' %(time.time()-t0))