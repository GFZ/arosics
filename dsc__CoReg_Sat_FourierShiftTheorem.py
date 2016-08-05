# -*- coding: utf-8 -*-
from __future__ import (division, print_function,absolute_import) #unicode_literals cause GDAL not to work properly

__author__ = "Daniel Scheffler"
__version__= "2016-08-02_03"

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

from shapely.geometry import Point, Polygon, box
from geopandas  import GeoDataFrame
from pykrige.ok import OrdinaryKriging

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


class CoReg(object):
    def __init__(self, path_im0,path_im1,path_out='.',r_b4match=1,s_b4match=1,wp=None,ws=(512,512),max_iter=5,max_shift=5,
                 align_grids=0,match_gsd=0,out_gsd=None,data_corners_im0=None,data_corners_im1=None,nodata=None,
                 calc_corners=1,multiproc=1,binary_ws=1,v=0,q=0,ignore_errors=0):
        self.mp                       = multiproc
        self.v                        = v
        self.q                        = q if not v else 0
        self.ignErr                   = ignore_errors

        assert os.path.exists(path_im0), "Reference image not found at '%s'. Please check the path."     %path_im0
        assert os.path.exists(path_im1), "Image to be shifted not found at '%s'. Please check the path." %path_im1
        for p,n in zip([path_im0,path_im1,path_out],['reference image','image to be shifted','output image']):
            assert ' ' not in p, "The path of the %s contains whitespaces. This is not supported by GDAL." %n
        self.path_imref               = path_im0
        self.path_im2shift            = path_im1
        if not self.q: print('reference image:', path_im0)
        if not self.q: print('shift image:    ', path_im1, '\n')
        dir_out, fName_out            = os.path.split(path_out)
        dir_out                       = dir_out if not dir_out  ==''  else os.path.dirname(path_im1)
        get_baseN                     = lambda path: os.path.splitext(os.path.basename(path))[0]
        fName_out                     = fName_out if not fName_out in ['.',''] else '%s__shifted_to__%s.bsq' \
                                         %(get_baseN(path_im1), get_baseN(path_im0))
        self.path_out                 = os.path.abspath(os.path.join(dir_out,fName_out))
        self.verbose_out              = os.path.abspath(os.path.join(dir_out,'CoReg_verboseOut__%s__shifted_to__%s' \
                                         %(get_baseN(path_im1), get_baseN(path_im0)) ))
        self.max_iter                 = max_iter
        self.max_shift                = max_shift
        self.max_win_sz_changes       = 3  ### TODO: änderung der window size, falls nach max_iter kein valider match gefunden
        self.align_grids              = align_grids
        self.match_gsd                = match_gsd
        self.out_gsd                  = out_gsd
        if match_gsd and out_gsd: warnings.warn("'-out_gsd' is ignored because '-match_gsd' is set.\n")
        if out_gsd:  assert isinstance(out_gsd,list) and len(out_gsd)==2, 'out_gsd must be a list with two values.'
        self.corner_coord_imref       = data_corners_im0
        self.corner_coord_im2shift    = data_corners_im1
        if data_corners_im0 and not isinstance(data_corners_im0[0],list): # group if not [[x,y],[x,y]..] but [x,y,x,y,]
            self.corner_coord_imref   = [data_corners_im0[i:i+2] for i in range(0, len(data_corners_im0), 2)]
        if data_corners_im1 and not isinstance(data_corners_im1[0],list): # group if not [[x,y],[x,y]..]
            self.corner_coord_im2shift= [data_corners_im1[i:i+2] for i in range(0, len(data_corners_im1), 2)]
        if self.v and not os.path.isdir(self.verbose_out): os.makedirs(self.verbose_out)

        if nodata: assert isinstance(nodata,list) and len(nodata)==2, 'nodata must be a list with two values.'
        self.ref_nodata, self.shift_nodata = \
            nodata if nodata else [GEO.find_noDataVal(self.path_imref), GEO.find_noDataVal(self.path_im2shift)]
        self.bin_ws = binary_ws

        gdal.AllRegister()

        self.ds_imref                  = gdal.Open(self.path_imref,0)    # 0 = GA_ReadOnly
        self.ds_im2shift               = gdal.Open(self.path_im2shift,0) # 0 = GA_ReadOnly
        assert self.ds_imref,    'Reference image can not be read by GDAL. You can try another image format.'
        assert self.ds_im2shift, 'The image to be shifted can not be read by GDAL. You can try another image format.'
        self.ref_prj,  self.shift_prj   = self.ds_imref.GetProjection(),   self.ds_im2shift.GetProjection()
        self.ref_gt ,  self.shift_gt    = self.ds_imref.GetGeoTransform(), self.ds_im2shift.GetGeoTransform()
        self.ref_xgsd, self.shift_xgsd  = abs(self.ref_gt[1]), abs(self.shift_gt[1])
        self.ref_ygsd, self.shift_ygsd  = abs(self.ref_gt[5]), abs(self.shift_gt[5])
        self.ref_rows, self.shift_rows  = self.ds_imref.RasterYSize, self.ds_im2shift.RasterYSize
        self.ref_cols, self.shift_cols  = self.ds_imref.RasterXSize, self.ds_im2shift.RasterXSize
        self.ref_bands,self.shift_bands = self.ds_imref.RasterCount, self.ds_im2shift.RasterCount
        assert self.ref_prj                             , 'Reference image has no projection.'
        assert not re.search('LOCAL_CS', self.ref_prj)  , 'Reference image is not georeferenced.'
        assert self.ref_gt                              , 'Reference image has no map information.'
        assert self.shift_prj                           , 'The image to be shifted has no projection.'
        assert not re.search('LOCAL_CS', self.shift_prj), 'The image to be shifted is not georeferenced.'
        assert self.shift_gt                            , 'The image to be shifted has no map information.'
        assert GEO.get_proj4info(self.ds_imref) == GEO.get_proj4info(self.ds_im2shift),\
            'Input projections are not equal. Different projections are currently not supported.'
        self.ds_imref = self.ds_im2shift = None # close GDAL datasets


        self.imref_band4match        = r_b4match
        self.im2shift_band4match     = s_b4match
        assert self.ref_bands >= self.imref_band4match >= 1, "The reference image has %s %s. So '-rb' must be %s%s." \
                % (self.ref_bands,'bands' if self.ref_bands>1 else 'band',
                   'between 1 and ' if self.ref_bands>1 else '', self.ref_bands)
        assert self.shift_bands >= self.im2shift_band4match >= 1, 'The image to be shifted '\
                "has %s %s. So '-sb' must be %s%s." % (self.shift_bands,'bands' if self.shift_bands>1 else 'band',
                'between 1 and ' if self.shift_bands>1 else '', self.shift_bands)

        if calc_corners:
            if self.corner_coord_imref is None:
                if not self.q: print('Calculating actual data corner coordinates for reference image...')
                self.corner_coord_imref    = GEO.get_true_corner_lonlat(self.path_imref, r_b4match,
                                                                        self.ref_nodata,self.mp)
            if self.corner_coord_im2shift is None:
                if not self.q: print('Calculating actual data corner coordinates for image to be shifted...')
                self.corner_coord_im2shift = GEO.get_true_corner_lonlat(self.path_im2shift, s_b4match,
                                                                        self.shift_nodata,self.mp)
        else:
            self.corner_coord_imref    = GEO.get_corner_coordinates(self.path_imref)
            self.corner_coord_im2shift = GEO.get_corner_coordinates(self.path_im2shift)
        if not self.q: print('Corner coordinates of reference image: %s'     %self.corner_coord_imref)
        if not self.q: print('Corner coordinates of image to be shifted: %s' %self.corner_coord_im2shift)

        self.poly_imref              = GEO.get_footprint_polygon(self.corner_coord_imref)
        self.poly_im2shift           = GEO.get_footprint_polygon(self.corner_coord_im2shift)
        overlap_tmp                  = GEO.get_overlap_polygon(self.poly_imref, self.poly_im2shift, self.v)
        self.overlap_poly            = overlap_tmp['overlap poly'] # has to be the reference projection

        for r_XY,s_XY in zip(self.corner_coord_imref,self.corner_coord_im2shift):
            assert self.poly_imref.contains(Point(r_XY))    or self.poly_imref.touches(Point(r_XY)),\
                "The corner position '%s' is outside of the reference image." %r_XY
            assert self.poly_im2shift.contains(Point(s_XY)) or self.poly_im2shift.touches(Point(s_XY)), \
                "The corner position '%s' is outside of the image to be shifted." %s_XY
        assert self.overlap_poly, 'The input images have no spatial overlap.'
        self.overlap_percentage      = overlap_tmp['overlap percentage']
        self.overlap_area            = overlap_tmp['overlap area']
        if self.v: IO.write_shp(os.path.join(self.verbose_out, 'poly_imref.shp')   , self.poly_imref,    self.ref_prj)
        if self.v: IO.write_shp(os.path.join(self.verbose_out, 'poly_im2shift.shp'), self.poly_im2shift, self.shift_prj)
        if self.v: IO.write_shp(os.path.join(self.verbose_out, 'overlap_poly.shp') , self.overlap_poly,  self.ref_prj)

        self.win_pos, opt_win_size_XY = self.get_opt_winpos_winsize(wp,ws)
        if not self.q: print('Matching window position (X,Y): %s/%s' %(self.win_pos[1],self.win_pos[0]))

        ### FIXME: transform_mapPt1_to_mapPt2(im2shift_center_map, ds_imref.GetProjection(), ds_im2shift.GetProjection()) # später basteln für den fall, dass projektionen nicht gleich sind

        # get_clip_window_properties
        self.ref_box_YX, self.shift_box_YX, self.ref_win_size_YX, self.shift_win_size_YX, \
        self.box_mapYX, self.poly_matchWin, self.poly_matchWin_prj, self.imfft_gsd \
            = self.get_clip_window_properties(opt_win_size_XY)


        if self.v and self.poly_matchWin:
            IO.write_shp(os.path.join(self.verbose_out, 'poly_matchWin.shp'), self.poly_matchWin, self.poly_matchWin_prj)

        self.imref_matchWin    = None # set by get_image_windows_to_match
        self.im2shift_matchWin = None # set by get_image_windows_to_match
        self.fftw_win_size     = None # set by get_opt_fftw_winsize

        self.x_shift_px        = None # always in shift image units (image coords) # set by calculate_spatial_shifts()
        self.y_shift_px        = None # always in shift image units (image coords) # set by calculate_spatial_shifts()
        self.success           = None if self.box_mapYX else False

    def get_opt_winpos_winsize(self,wp,ws):
        # type: (tuple,tuple) -> tuple,tuple
        """Calculates optimal window position and size in reference image units according to DGM, cloud_mask and
        trueCornerLonLat."""
        # dummy algorithm: get center position of overlap instead of searching ideal window position in whole overlap
        # TODO automatischer Algorithmus zur Bestimmung der optimalen Window Position
        if wp is None:
            overlap_center_pos_x, overlap_center_pos_y = self.overlap_poly.centroid.coords.xy
            win_pos = (overlap_center_pos_y[0], overlap_center_pos_x[0])
        else:
            wp = wp if not isinstance(wp,np.ndarray) else wp.tolist()
            assert isinstance(wp,list), 'The window position must be a list of two elements. Got %s.' %wp
            assert len(wp)==2,          'The window position must be a list of two elements. Got %s.' %wp
            assert self.overlap_poly.contains(Point(wp)), 'The provided window position is ' \
                'outside of the overlap area of the two input images. Check the coordinates.'
            win_pos = tuple(reversed(wp))

        # TODO automatischer Algorithmus zur Bestimmung der optimalen Window Position
        gained_win_size_XY = (int(ws[0]), int(ws[1])) if ws else (512,512)
        return win_pos, gained_win_size_XY

    @staticmethod
    def get_winPoly(wp_imYX,ws,gt,match_grid=0):
        # type: (tuple, tuple, list, bool) -> shapely.Polygon, tuple, tuple
        """Creates a shapely polygon from a given set of image cordinates, window size and geotransform.
        :param wp_imYX:
        :param ws:          <tuple> X/Y window size
        :param gt:
        :param match_grid:  <bool> """

        ws = (ws,ws) if isinstance(ws,int) else ws
        xmin,xmax,ymin,ymax = (wp_imYX[1]-ws[0]/2, wp_imYX[1]+ws[0]/2, wp_imYX[0]+ws[1]/2, wp_imYX[0]-ws[1]/2)
        if match_grid:
            xmin,xmax,ymin,ymax = [int(i) for i in [xmin,xmax,ymin,ymax]]
        box_YX      = (ymax,xmin),(ymax,xmax),(ymin,xmax),(ymin,xmin) # UL,UR,LR,LL
        UL,UR,LR,LL = [GEO.imYX2mapYX(imYX, gt) for imYX in box_YX]
        box_mapYX   = UL,UR,LR,LL
        return Polygon([(UL[1],UL[0]),(UR[1],UR[0]),(LR[1],LR[0]),(LL[1],LR[0])]), box_mapYX, box_YX

    def get_clip_window_properties_OLD(self, dst_ws_XY): # TODO alles auf shapely polygone umstellen?
        """TODO
        hint: Even if X- and Y-dimension of the target window is equal, the output window can be NOT quadratic!"""

        # Fenster-Positionen in Bildkoordinaten in imref und im2shift ausrechnen
        refYX, shiftYX = GEO.mapYX2imYX(self.win_pos, self.ref_gt), GEO.mapYX2imYX(self.win_pos, self.shift_gt)
        # maximale Fenster-Größen in imref und im2shift ausrechnen
        max_ref_clipWS   = 2*int(min(list(refYX)   + [self.ref_cols  -refYX[1],  self.ref_rows  -refYX[0]])  ) # TODO auf XY-dim anpassen
        max_shift_clipWS = 2*int(min(list(shiftYX) + [self.shift_cols-shiftYX[1],self.shift_rows-shiftYX[0]])) # TODO auf XY-dim anpassen
        print(max_ref_clipWS)
        print(max_shift_clipWS)
        max_ref_clipWS = 20 # FIXME

        for maxWs,imName in zip([max_ref_clipWS, max_shift_clipWS], ['reference image', 'image to be shifted']):
            if maxWs<8:
                if not self.ignErr:
                    raise RuntimeError('The matching window is too close to the borders of the %s. '
                                       'Matching failed.' %imName)
                else: return [None]*8

        gsdFact = max([self.ref_xgsd,self.shift_xgsd])/min([self.ref_xgsd,self.shift_xgsd]) # immer positiv

        # Fenster-Größen in imref und im2shift anhand Ziel-Größe und maximal möglicher Größe setzen

        if self.shift_xgsd <= self.ref_xgsd:
            """Matching-Fall 1: Shift-Auflösung >= Ref-Auflösung"""
            imfft_gsd_mapvalues = self.ref_xgsd
            ref_clipWS_tmp   = dst_ws_XY if dst_ws_XY <= max_ref_clipWS   else max_ref_clipWS
            shift_clipWS_tmp = dst_ws_XY * gsdFact if dst_ws_XY * gsdFact <= max_shift_clipWS else max_shift_clipWS

            ref_winPoly,  ref_box_mapYX,  ref_box_YX   = self.get_winPoly(refYX,  ref_clipWS_tmp,  self.ref_gt, match_grid=1)
            shift_winPoly,shift_box_mapYX,shift_box_YX = self.get_winPoly(shiftYX,shift_clipWS_tmp,self.shift_gt)
#            IO.write_shp('/misc/hy5/scheffler/tmp/ref_win.shp',ref_winPoly,prj=self.ref_prj)
#            IO.write_shp('/misc/hy5/scheffler/tmp/shift_win.shp',shift_winPoly,prj=self.shift_prj)

            if ref_winPoly.within(shift_winPoly):
                # wenn Ref Fenster innerhalb des Shift-Fensters:
                # dann direkt übernehmen und shift_box_XY an ref_box_mapYX anpassen
                ref_box_YX, box_mapYX = ref_box_YX, ref_box_mapYX
                #shift_box_YX = [mapYX2imYX(YX,self.shift_gt) for YX in box_mapYX]
                shift_box_YX = GEO.get_smallest_boxImYX_that_contains_boxMapYX(box_mapYX, self.shift_gt)
            else: # wenn Ref Fenster zu groß für Shift-Bild: dann Fenster vom Shift-Bild als Ref-Fenster verwenden
                ref_box_YX_flt = [GEO.mapYX2imYX(YX, self.ref_gt) for YX in shift_box_mapYX]
                ref_box_YX     = [[int(YX[0]),int(YX[1])] for YX in ref_box_YX_flt]
                box_mapYX      = [GEO.imYX2mapYX(YX, self.ref_gt) for YX in ref_box_YX]
                shift_box_YX   = GEO.get_smallest_boxImYX_that_contains_boxMapYX(box_mapYX, self.shift_gt)

        else:
            """Matching-Fall 2: Ref-Auflösung > Shift-Auflösung"""
            imfft_gsd_mapvalues = self.shift_xgsd # fft-Bild bekommt Auflösung des warp-Bilds
            ref_clipWS_tmp   = dst_ws_XY * gsdFact if dst_ws_XY * gsdFact <= max_ref_clipWS   else max_ref_clipWS
            shift_clipWS_tmp = dst_ws_XY         if dst_ws_XY <= max_shift_clipWS else max_shift_clipWS
            ref_winPoly,  ref_box_mapYX,  ref_box_YX   = self.get_winPoly(refYX,  ref_clipWS_tmp,  self.ref_gt )
            shift_winPoly,shift_box_mapYX,shift_box_YX = self.get_winPoly(shiftYX,shift_clipWS_tmp,self.shift_gt,match_grid=1)

            if shift_winPoly.within(ref_winPoly):
                # wenn Shift Fenster innerhalb des Ref-Fensters:
                # dann direkt übernehmen und ref_box_XY an shift_box_mapYX anpassen
                shift_box_YX, box_mapYX = shift_box_YX, shift_box_mapYX
                #ref_box_YX              = [mapYX2imYX(YX,self.ref_gt) for YX in box_mapYX]
                ref_box_YX              = GEO.get_smallest_boxImYX_that_contains_boxMapYX(box_mapYX, self.ref_gt)
            else: # wenn Shift Fenster zu groß für Ref-Bild: dann Fenster vom Ref-Bild als Shift-Fenster verwenden
                shift_box_YX_flt = [GEO.mapYX2imYX(YX, self.shift_gt) for YX in ref_box_mapYX] # bildkoordinaten von mapYX
                shift_box_YX     = [[int(YX[0]),int(YX[1])] for YX in shift_box_YX_flt]   # auf grid schieben
                box_mapYX        = [GEO.imYX2mapYX(YX, self.shift_gt) for YX in shift_box_YX]  # mapYX aus shift-Bildkoord.
                ref_box_YX       = ref_box_YX
                ref_box_YX       = GEO.get_smallest_boxImYX_that_contains_boxMapYX(box_mapYX, self.ref_gt)

        ## Test, ob zusätzlicher Buffer der Fenstergröße zum Umprojizieren gebraucht wird
        #buffFact = 1.0 if get_proj4info(self.ref_prj)==get_proj4info(self.shift_prj) else 1.5

        poly_matchWin     = GEO.get_footprint_polygon([(YX[1], YX[0]) for YX in box_mapYX])
        poly_matchWin_prj = self.ref_prj if self.shift_xgsd <= self.ref_xgsd else self.shift_prj
        ref_box_YX        = [[int(i[0]),int(i[1])] for i in ref_box_YX]
        shift_box_YX      = [[int(i[0]),int(i[1])] for i in shift_box_YX]
        ref_win_size_YX   = abs(ref_box_YX  [2][0]-ref_box_YX  [0][0]), abs(ref_box_YX  [1][1]-ref_box_YX  [0][1]) # LLy-ULy, URx-ULx
        shift_win_size_YX = abs(shift_box_YX[2][0]-shift_box_YX[0][0]), abs(shift_box_YX[1][1]-shift_box_YX[0][1]) # LLy-ULy, URx-ULx

        return ref_box_YX, shift_box_YX, ref_win_size_YX, shift_win_size_YX, \
               box_mapYX, poly_matchWin, poly_matchWin_prj, imfft_gsd_mapvalues

    def move_shapelyPoly_to_image_grid(self,shapelyPoly,gt,rows,cols,moving_dir='NW'):
        polyULxy       = (min(shapelyPoly.exterior.coords.xy[0]), max(shapelyPoly.exterior.coords.xy[1]))
        polyLRxy       = (max(shapelyPoly.exterior.coords.xy[0]), min(shapelyPoly.exterior.coords.xy[1]))
        UL, LL, LR, UR = GEO.get_corner_coordinates(gt=gt,rows=rows,cols=cols) # (x,y) tuples
        round_x        = {'NW':'off','NO':'on','SW':'off','SE':'on'}[moving_dir]
        round_y        = {'NW':'on','NO':'on','SW':'off','SE':'off'}[moving_dir]
        tgt_xgrid      = np.arange(UL[0],UR[0]+gt[1]     , gt[1])
        tgt_ygrid      = np.arange(LL[1],UL[1]+abs(gt[5]), abs(gt[5]))
        tgt_xmin       = UTL.find_nearest(tgt_xgrid, polyULxy[0], round=round_x)
        tgt_xmax       = UTL.find_nearest(tgt_xgrid, polyLRxy[0], round=round_x)
        tgt_ymin       = UTL.find_nearest(tgt_ygrid, polyLRxy[1], round=round_y)
        tgt_ymax       = UTL.find_nearest(tgt_ygrid, polyULxy[1], round=round_y)
        from shapely.geometry import box
        return box(tgt_xmin,tgt_ymin,tgt_xmax,tgt_ymax)

    def get_clip_window_properties(self, dst_ws_XY):
            """TODO
            hint: Even if X- and Y-dimension of the target window is equal, the output window can be NOT quadratic!"""

            wp_Y,wp_X = self.win_pos
            ws_X,ws_Y = dst_ws_XY[0]*self.ref_xgsd, dst_ws_XY[1]*self.ref_ygsd # image units -> map units
            tgt_matchPoly = box(wp_X-ws_X/2, wp_Y-ws_Y/2, wp_X+ws_X/2, wp_Y+ws_Y/2) # minx,miny,maxx,maxy
            matchPoly     = tgt_matchPoly.intersection(self.overlap_poly) # clip matching window to overlap area

            # move matching window to imref grid or im2shift grid
            matchPoly = self.move_shapelyPoly_to_image_grid(matchPoly,self.ref_gt,self.ref_rows,self.ref_cols,'NW') \
                if self.shift_xgsd <= self.ref_xgsd else \
                self.move_shapelyPoly_to_image_grid(matchPoly, self.shift_gt, self.shift_rows, self.shift_cols, 'NW')

            # check, ob durch Vershiebung auf Grid das matchWin außerhalb von overlap_poly geschoben wurde
            if not matchPoly.within(self.overlap_poly):
                # matchPoly weiter verkleinern
                matchPoly = matchPoly.buffer(-self.ref_xgsd)

            ref_winImPoly   = Polygon(GEO.get_imageCoords_from_shapelyPoly(matchPoly, self.ref_gt)) # asserts equal projections
            shift_winImPoly = Polygon(GEO.get_imageCoords_from_shapelyPoly(matchPoly, self.shift_gt)) # asserts equal projections

            ref_winMapPoly   = GEO.shapelyImPoly_to_shapelyMapPoly(ref_winImPoly, self.ref_gt, self.ref_prj)
            shift_winMapPoly = GEO.shapelyImPoly_to_shapelyMapPoly(shift_winImPoly, self.shift_gt, self.shift_prj)

            if self.shift_xgsd <= self.ref_xgsd:
                """Matching-Fall 1: Shift-Auflösung >= Ref-Auflösung => Grid von imref übernehmen, im2shift anpassen"""
                imfft_gsd_mapvalues = self.ref_xgsd
                # matching_win und ref_winImPoly direkt auf imref grid
                ref_winImPoly = GEO.round_shapelyPoly_coords(ref_winImPoly, precision=0, out_dtype=int) # Rundungsfehler bei Koordinatentrafo beseitigen

                if ref_winMapPoly.within(shift_winMapPoly) or ref_winMapPoly==shift_winMapPoly:
                    # wenn Ref Fenster innerhalb des Shift-Fensters: dann direkt übernehmen
                    pass
                else:
                    # wenn Ref Fenster größer als Shift-Fenster: dann kleinstes Shift-Fenster finden, das Ref-Fenster umgibt
                    ref_winMapPoly = GEO.shapelyImPoly_to_shapelyMapPoly(ref_winImPoly, self.ref_gt, self.ref_prj)
                    ref_boxMapYX = GEO.shapelyBox2BoxYX(ref_winMapPoly,coord_type='map')
                    shift_box_YX = GEO.get_smallest_boxImYX_that_contains_boxMapYX(ref_boxMapYX,self.shift_gt)

                    # update shift_winImPoly
                    shift_box_XY = [(i[1],i[0]) for i in shift_box_YX]
                    xmin, xmax, ymin, ymax = GEO.corner_coord_to_minmax(shift_box_XY)
                    shift_winImPoly = box(xmin, ymin, xmax, ymax)

                    overlapImPoly  = Polygon(GEO.get_imageCoords_from_shapelyPoly(self.overlap_poly, self.shift_gt)) # asserts equal projections

                    while not shift_winImPoly.within(overlapImPoly): # evtl. kann es sein, dass bei Shift-Fenster-Vergrößerung das shift-Fenster zu groß für den overlap wird -> ref-Fenster verkleinern und neues shift-fenster berechnen
                        ref_winImPoly  = ref_winImPoly.buffer(-1)
                        ref_winMapPoly = GEO.shapelyImPoly_to_shapelyMapPoly(ref_winImPoly, self.ref_gt, self.ref_prj)
                        ref_boxMapYX   = GEO.shapelyBox2BoxYX(ref_winMapPoly, coord_type='map')
                        shift_box_YX   = GEO.get_smallest_boxImYX_that_contains_boxMapYX(ref_boxMapYX, self.shift_gt)
                        shift_box_XY   = [(i[1], i[0]) for i in shift_box_YX]
                        xmin, xmax, ymin, ymax = GEO.corner_coord_to_minmax(shift_box_XY)
                        shift_winImPoly = box(xmin, ymin, xmax, ymax)
                    shift_winMapPoly = GEO.shapelyImPoly_to_shapelyMapPoly(shift_winImPoly, self.shift_gt, self.shift_prj)

                assert ref_winMapPoly.within(shift_winMapPoly)
                matchPoly = GEO.shapelyImPoly2shapelyMapPoly(ref_winImPoly, self.ref_gt)

            else:
                """Matching-Fall 2: Ref-Auflösung > Shift-Auflösung => Grid von im2shift übernehmen, imref anpassen"""
                imfft_gsd_mapvalues = self.shift_xgsd  # fft-Bild bekommt Auflösung des warp-Bilds
                #  matching_win und shift_winImPoly direkt auf im2shift grid
                shift_winImPoly = GEO.round_shapelyPoly_coords(shift_winImPoly, precision=0, out_dtype=int)  # Rundungsfehler bei Koordinatentrafo beseitigen

                if shift_winMapPoly.within(ref_winMapPoly) or shift_winMapPoly == ref_winMapPoly:
                    # wenn Ref Fenster innerhalb des Shift-Fensters: dann direkt übernehmen
                    pass
                else:
                    # wenn Shift Fenster größer als Ref-Fenster: dann kleinstes Ref-Fenster finden, das shift-Fenster umgibt
                    shift_winMapPoly = GEO.shapelyImPoly_to_shapelyMapPoly(shift_winImPoly, self.shift_gt, self.shift_prj)
                    shift_boxMapYX = GEO.shapelyBox2BoxYX(shift_winMapPoly, coord_type='map')
                    ref_box_YX = GEO.get_smallest_boxImYX_that_contains_boxMapYX(shift_boxMapYX, self.ref_gt)

                    # update ref_winImPoly
                    ref_box_XY = [(i[1], i[0]) for i in ref_box_YX]
                    xmin, xmax, ymin, ymax = GEO.corner_coord_to_minmax(ref_box_XY)
                    ref_winImPoly = box(xmin, ymin, xmax, ymax)

                    overlapImPoly = Polygon(GEO.get_imageCoords_from_shapelyPoly(self.overlap_poly, self.ref_gt))  # asserts equal projections

                    while not ref_winImPoly.within(overlapImPoly):  # evtl. kann es sein, dass bei Ref-Fenster-Vergrößerung das Ref-Fenster zu groß für den overlap wird -> shift-Fenster verkleinern und neues Ref-fenster berechnen
                        shift_winImPoly = shift_winImPoly.buffer(-1)
                        shift_winMapPoly = GEO.shapelyImPoly_to_shapelyMapPoly(shift_winImPoly, self.shift_gt, self.shift_prj)
                        shift_boxMapYX = GEO.shapelyBox2BoxYX(shift_winMapPoly, coord_type='map')
                        ref_box_YX = GEO.get_smallest_boxImYX_that_contains_boxMapYX(shift_boxMapYX, self.ref_gt)
                        ref_box_XY = [(i[1], i[0]) for i in ref_box_YX]
                        xmin, xmax, ymin, ymax = GEO.corner_coord_to_minmax(ref_box_XY)
                        ref_winImPoly = box(xmin, ymin, xmax, ymax)
                    ref_winMapPoly = GEO.shapelyImPoly_to_shapelyMapPoly(ref_winImPoly, self.ref_gt,
                                                                           self.ref_prj)

                assert shift_winMapPoly.within(ref_winMapPoly)
                matchPoly = GEO.shapelyImPoly2shapelyMapPoly(shift_winImPoly, self.shift_gt)

            ref_box_YX        = [[int(i[0]),int(i[1])] for i in GEO.shapelyBox2BoxYX(ref_winImPoly  ,coord_type='image')]
            shift_box_YX      = [[int(i[0]),int(i[1])] for i in GEO.shapelyBox2BoxYX(shift_winImPoly,coord_type='image')]
            box_mapYX         = GEO.shapelyBox2BoxYX(matchPoly,coord_type='map')
            poly_matchWin_prj = self.ref_prj if self.shift_xgsd <= self.ref_xgsd else self.shift_prj
            ref_win_size_YX   = abs(ref_box_YX[2][0]  -ref_box_YX[0][0]),   abs(ref_box_YX[1][1]  -ref_box_YX[0][1])  # LLy-ULy, URx-ULx
            shift_win_size_YX = abs(shift_box_YX[2][0]-shift_box_YX[0][0]), abs(shift_box_YX[1][1]-shift_box_YX[0][1])# LLy-ULy, URx-ULx
            match_win_size_XY = tuple(reversed(ref_win_size_YX if self.shift_xgsd <= self.ref_xgsd else shift_win_size_YX))
            if ref_win_size_YX != match_win_size_XY:
                print('Target window size %s not possible due to too small overlap area or window position too close '
                      'to an image edge. New matching window size: %s.' %(dst_ws_XY,match_win_size_XY))

            #print(ref_win_size_YX)
            #print(shift_win_size_YX)
            #print('RWP', ref_winImPoly)
            #print('SWP', shift_winImPoly)
            #print(ref_winMapPoly.bounds)
            #print(shift_winMapPoly.bounds)
            #print(ref_box_YX)
            #print(shift_box_YX)
            # IO.write_shp('/misc/hy5/scheffler/Temp/ref_winMapPoly.shp', ref_winMapPoly,self.ref_prj)
            # IO.write_shp('/misc/hy5/scheffler/Temp/shift_winMapPoly.shp', shift_winMapPoly, self.shift_prj)

            return ref_box_YX, shift_box_YX, ref_win_size_YX, shift_win_size_YX, \
                   box_mapYX, matchPoly, poly_matchWin_prj, imfft_gsd_mapvalues

    def get_image_windows_to_match(self):
        if self.v: print('resolutions: ', self.ref_xgsd, self.shift_xgsd)
        box_mapXvals = [i[1] for i in self.box_mapYX] # UTM-RW
        box_mapYvals = [i[0] for i in self.box_mapYX] # UTM-HW

        is_avail_rsp_average = int(gdal.VersionInfo()[0]) >= 2
        if not is_avail_rsp_average:
            warnings.warn("The GDAL version on this server does not yet support the resampling algorithm "
                          "'average'. This can affect the correct detection of subpixel shifts. To avoid this "
                          "please update GDAL to a version above 2.0.0!")

        if self.shift_xgsd <= self.ref_xgsd: # Matching-Fall 1: Shift-Auflösung >= Ref-Auflösung
            """falls referenz die niedrigere auflösung hat oder auflösungen gleich sind:
            im2shift auf imref downsamplen"""
            # 1. Eckkoordinaten Warp-Bild -> per Subset-Read referenz-Bild mit vorgegebener Fenstergröße (wie in späterer fft, um Rechenzeit zu sparen) einlesen -> zwischenspeichern in dev/shm/
            # 2. warp-Bild (bzw. subset je nach fft-fenstergröße) downsamplen auf Auflösung des referenz-bildes per pixel-aggregate -> dev/shm/    ## http://gis.stackexchange.com/questions/110769/gdal-python-aggregate-raster-into-lower-resolution
            # 3. referenz und warp-bild exakt in der fft-fenstergröße einlesen

            # imref per subset-read einlesen -> im0
            gdalReadInputs  = GEO.get_gdalReadInputs_from_boxImYX(self.ref_box_YX)
            assert np.array_equal(np.abs(np.array(gdalReadInputs)),np.array(gdalReadInputs)),\
                'Got negative values in gdalReadInputs for refernce image.'
            ds_imref        = gdal.Open(self.path_imref)
            imref_clip_data = ds_imref.GetRasterBand(self.imref_band4match).ReadAsArray(*gdalReadInputs)
            ds_imref        = None

            # im2shift per subset-read einlesen
            gdalReadInputs  = GEO.get_gdalReadInputs_from_boxImYX(self.shift_box_YX)
            assert np.array_equal(np.abs(np.array(gdalReadInputs)),np.array(gdalReadInputs)),\
                'Got negative values in gdalReadInputs for image to be shifted.'
            ds_im2shift             = gdal.Open(self.path_im2shift)
            im2shift_clip_data_orig = ds_im2shift.GetRasterBand(self.im2shift_band4match).ReadAsArray(*gdalReadInputs)
            ds_im2shift             = None

            if self.v:
                print('Original matching windows:')
                PLT.subplot_imshow([imref_clip_data,im2shift_clip_data_orig],
                                   [os.path.basename(self.path_imref), os.path.basename(self.path_im2shift)], grid=True)

            im2shift_clip_gt = GEO.get_subset_GeoTransform(self.shift_gt, self.shift_box_YX)

            # resample im2shift to the resolution of imref AND make sure the pixel edges are identical
            # (in order to make each image show the same window with the same coordinates)
            #rsp_algor = 5 if is_avail_rsp_average else 2 # average if possible else cubic # OLD
            rsp_algor = 2 # TODO replace cubic resampling by PSF resampling - average resampling leads to sinus like distortions in the fft image that make a precise coregistration impossible. Thats why there is currently no way around cubic resampling.
            im2shift_clip_data = GEO.warp_ndarray(im2shift_clip_data_orig, im2shift_clip_gt, self.shift_prj,
                                                  self.ref_prj, out_res=(self.imfft_gsd, self.imfft_gsd),
                                                  out_extent=([min(box_mapXvals), min(box_mapYvals),
                                                               max(box_mapXvals), max(box_mapYvals)]),
                                                  rsp_alg=rsp_algor, in_nodata=self.shift_nodata) [0]
        else:
            """falls referenz die höhere auflösung hat"""
            """imref auf im2shift downsamplen"""  ## funktioniert perfekt (getestet an UTM_Z43:5m - UTM_Z43:30m)

            # 1. shift-window ausclippen
            # 2. ref-window per gdalwarp anhand clip-koordinaten von shift ausschneiden, ggf. umprojizieren und downsamplen

            # im2shift per subset-read einlesen -> im0
            gdalReadInputs     = GEO.get_gdalReadInputs_from_boxImYX(self.shift_box_YX)
            assert np.array_equal(np.abs(np.array(gdalReadInputs)),np.array(gdalReadInputs)),\
                'Got negative values in gdalReadInputs for refernce image.'
            ds_im2shift        = gdal.Open(self.path_im2shift)
            im2shift_clip_data = ds_im2shift.GetRasterBand(self.im2shift_band4match).ReadAsArray(*gdalReadInputs)
            ds_im2shift        = None

            # imref per subset-read einlesen
            gdalReadInputs       = GEO.get_gdalReadInputs_from_boxImYX(self.ref_box_YX)
            assert np.array_equal(np.abs(np.array(gdalReadInputs)),np.array(gdalReadInputs)),\
                'Got negative values in gdalReadInputs for image to be shifted.'
            ds_imref             = gdal.Open(self.path_imref)
            imref_clip_data_orig = ds_imref.GetRasterBand(self.imref_band4match).ReadAsArray(*gdalReadInputs)
            ds_imref             = None

            #subplot_imshow([imref_clip_data_orig, im2shift_clip_data],grid=True)
            imref_clip_gt = GEO.get_subset_GeoTransform(self.ref_gt, self.ref_box_YX)
            #print(get_corner_coordinates(gt=imref_clip_gt,rows=clip_sz,cols=clip_sz))

#            rsp_algor = 5 if is_avail_rsp_average else 2 # average if possible else cubic
            rsp_algor = 2  # TODO replace cubic resampling by PSF resampling - average resampling leads to sinus like distortions in the fft image that make a precise coregistration impossible. Thats why there is currently no way around cubic resampling.
            imref_clip_data = GEO.warp_ndarray(imref_clip_data_orig, imref_clip_gt, self.ref_prj, self.shift_prj,
                                               out_res=(self.imfft_gsd, self.imfft_gsd),
                                               out_extent=([min(box_mapXvals), min(box_mapYvals),
                                                            max(box_mapXvals),max(box_mapYvals)]),
                                               rsp_alg=rsp_algor, in_nodata=self.ref_nodata) [0]
            #subplot_imshow([imref_clip_data,im2shift_clip_data],grid=True)

            # #########################################################################################################
            # # 1. Eckkoordinaten des Clip-Windows für imref berechnen (anhand optimaler fft-Fenstergröße und Zentrumspunkt der Bildüberlappung)
            # # 2. imref ausschneiden und per pixel-aggregate auf auflösung von im2shift downsamplen -> dev/shm/ -> einlesen als array ## http://gis.stackexchange.com/questions/110769/gdal-python-aggregate-raster-into-lower-resolution
            # # 3. warp-bild anhand des Clip-Windows des Referenz-bildes ausschneiden und auflösung unverändert lassen
            #
            # # imref per subset-read einlesen -> im0_orig
            # rS, cS  = self.ref_box_YX[0] # UL
            # clip_sz = abs(self.ref_box_YX[1][0]-self.ref_box_YX[0][0])
            #
            # ds_imref             = gdal.Open(self.path_imref)
            # imref_clip_data_orig = ds_imref.GetRasterBand(self.imref_band4match).ReadAsArray(cS, rS, clip_sz, clip_sz)
            # ds_imref             = None
            #
            # # imref mit imref_box ausclippen und per pixel aggregate auf auflösung von im2shift downsamplen
            # imref_clip_gt = list(self.ref_gt[:])
            # imref_clip_gt[3],imref_clip_gt[0] = \
            #     pixelToMapYX([int(cS),int(rS)],geotransform=self.ref_gt,projection=self.ref_prj)[0]
            # rsp_algor = 5 # average
            # if gdal.VersionInfo().startswith('1'):
            #     warnings.warn("The GDAL version on this server does not yet support the resampling algorithm "
            #                   "'average'. This can affect the correct detection of subpixel shifts. To avoid this "
            #                   "please update GDAL to a version above 2.0.0!")
            #     rsp_algor = 2 # cubic
            # imref_clip_data = warp_ndarray(imref_clip_data_orig,imref_clip_gt,self.ref_prj,self.ref_prj,
            #     out_res=(self.imfft_gsd, self.imfft_gsd),out_extent=([min(box_mapXvals), min(box_mapYvals),
            #     max(box_mapXvals),max(box_mapYvals)]), rsp_alg=rsp_algor, in_nodata=self.ref_nodata) [0]
            #     # target extent srs is not neccessary because -t_srs = imref_proj and output extent is derived from imref
            #
            # # im2shift per subset-read einlesen
            # rS, cS  = self.shift_box_YX[0] # UL
            # rS, cS  = rS-1, cS-1 # Input-Bild muss etwas größer sein als out_extent
            # clip_sz = abs(self.shift_box_YX[1][0]-self.shift_box_YX[0][0])+2
            # ds_im2shift             = gdal.Open(self.path_im2shift)
            # im2shift_clip_data_orig = \
            #     ds_im2shift.GetRasterBand(self.im2shift_band4match).ReadAsArray(cS, rS, clip_sz, clip_sz)
            # ds_im2shift             = None
            #
            # # im2shift mit imref_box_map ausclippen (cubic, da pixelgrenzen verschoben werden müssen aber auflösung
            # # gleich bleibt)
            # im2shift_clip_gt = list(self.shift_gt[:])
            # im2shift_clip_gt[3],im2shift_clip_gt[0] = \
            #     pixelToMapYX([int(cS),int(rS)],geotransform=self.shift_gt,projection=self.shift_prj)[0]
            # im2shift_clip_data = warp_ndarray(im2shift_clip_data_orig,im2shift_clip_gt,self.shift_prj,self.ref_prj,
            #     out_res=(self.imfft_gsd, self.imfft_gsd),out_extent=([min(box_mapXvals), min(box_mapYvals),
            #     max(box_mapXvals),max(box_mapYvals)]), rsp_alg=2, in_nodata=self.shift_nodata) [0]

        if imref_clip_data.shape != im2shift_clip_data.shape:
            raise RuntimeError('Bad output of get_image_windows_to_match. Reference image shape is %s whereas shift'
                               ' image shape is %s.' %(imref_clip_data.shape, im2shift_clip_data.shape))
        rows,cols = [i if i%2==0 else i-1 for i in imref_clip_data.shape]
        imref_clip_data, im2shift_clip_data = imref_clip_data[:rows,:cols], im2shift_clip_data[:rows,:cols]

        assert imref_clip_data is not None and im2shift_clip_data is not None, 'Creation of matching windows failed.'

        self.imref_matchWin    = imref_clip_data
        self.im2shift_matchWin = im2shift_clip_data

        return imref_clip_data, im2shift_clip_data

    def get_opt_fftw_winsize(self,im_shape, target_size=None,v=0,ignErr=0,bin_ws=1):
        binarySizes = [2**i for i in range(5,14)] # [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
        possibSizes = [i for i in binarySizes if i <= min(im_shape)] if bin_ws else list(range(8192))
        if not possibSizes:
            if not ignErr:  raise RuntimeError('The matching window became too small for calculating a reliable match. '
                                               'Matching failed.')
            else:           return None
        target_size = max(possibSizes) if target_size is None else target_size
        closest_to_target = int(min(possibSizes, key=lambda x:abs(x-target_size)))
        if v: print('final window size: %s' %closest_to_target)
        self.fftw_win_size = closest_to_target
        return closest_to_target

    def calc_shifted_cross_power_spectrum(self,im0,im1,window_size=1024,precision=np.complex64, v=0, ignErr=0, bin_ws=1):
        """Calculates shifted cross power spectrum for quantifying x/y-shifts.
            :param im0:         reference image
            :param im1:         subject image to shift
            :param window_size: size of image area to be processed
            :param precision:   to be quantified as a datatype
            :param v:           verbose
            :return:            2D-numpy-array of the shifted cross power spectrum
        """
        window_size = self.get_opt_fftw_winsize(im0.shape, target_size=window_size,v=v,ignErr=ignErr,bin_ws=bin_ws)
        assert im0.shape == im1.shape, 'The reference and the target image must have the same dimensions.'
        if im0.shape[0]%2!=0: warnings.warn('Odd row count in one of the match images!')
        if im1.shape[1]%2!=0: warnings.warn('Odd column count in one of the match images!')

        #print('fft_size',window_size)
        if window_size:
            t0 = time.time()

            center_YX = np.array(im0.shape)/2
            xmin,xmax,ymin,ymax = int(center_YX[1]-window_size/2), int(center_YX[1]+window_size/2),\
                                  int(center_YX[0]-window_size/2), int(center_YX[0]+window_size/2)
            in_arr0  = im0[ymin:ymax,xmin:xmax].astype(precision)
            in_arr1  = im1[ymin:ymax,xmin:xmax].astype(precision)

            if self.v:
                PLT.subplot_imshow([in_arr0.astype(np.float32), in_arr1.astype(np.float32)],
                               ['FFTin_'+os.path.basename(self.path_imref),
                                'FFTin_'+os.path.basename(self.path_im2shift)], grid=True)

            if 'pyfft' in globals():
                fft_arr0 = pyfftw.FFTW(in_arr0,np.empty_like(in_arr0), axes=(0,1))()
                fft_arr1 = pyfftw.FFTW(in_arr1,np.empty_like(in_arr1), axes=(0,1))()
            else:
                fft_arr0 = np.fft.fft2(in_arr0)
                fft_arr1 = np.fft.fft2(in_arr1)
            if v: print('forward FFTW: %.2fs' %(time.time() -t0))

            eps = np.abs(fft_arr1).max() * 1e-15
            # cps == cross-power spectrum of im0 and im2

            temp = (fft_arr0 * fft_arr1.conjugate()) / (np.abs(fft_arr0) * np.abs(fft_arr1) + eps)

            t0 = time.time()
            if 'pyfft' in globals():
                ifft_arr = pyfftw.FFTW(temp,np.empty_like(temp), axes=(0,1),direction='FFTW_BACKWARD')()
            else:
                ifft_arr = np.fft.ifft2(temp)
            if v: print('backward FFTW: %.2fs' %(time.time() -t0))

            cps = np.abs(ifft_arr)
            # scps = shifted cps
            scps = np.fft.fftshift(cps)
            if v:
                PLT.subplot_imshow([in_arr0.astype(np.uint16), in_arr1.astype(np.uint16), fft_arr0.astype(np.uint8),
                                fft_arr1.astype(np.uint8), scps], titles=['matching window im0', 'matching window im1',
                                "fft result im0", "fft result im1", "cross power spectrum"], grid=True)
                PLT.subplot_3dsurface(scps.astype(np.float32))
        else:
            scps = None
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
        wsY,wsX   = winSzYX
        xmin,xmax,ymin,ymax = get_bounds(center_YX,wsY,wsX)
        return im[ymin:ymax,xmin:xmax]

    def get_grossly_deshifted_images(self,im0,im1,x_intshift,y_intshift,v=0):
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

        if v:
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

    def calc_subpixel_shifts(self,scps,v=0):
        sidemax_lr, sidemax_ab = self.find_side_maximum(scps,v)
        x_subshift = (sidemax_lr['direction_factor']*sidemax_lr['value'])/(np.max(scps)+sidemax_lr['value'])
        y_subshift = (sidemax_ab['direction_factor']*sidemax_ab['value'])/(np.max(scps)+sidemax_ab['value'])
        return x_subshift, y_subshift

    @staticmethod
    def get_total_shifts(x_intshift,y_intshift,x_subshift,y_subshift):
        return x_intshift+x_subshift, y_intshift+y_subshift

    def calculate_spatial_shifts(self):
        if self.success is False: return None,None

        if self.q:  warnings.simplefilter('ignore')
        im0, im1 = self.get_image_windows_to_match()
        #
        #write_envi(im0, '/misc/hy5/scheffler/Projekte/2016_03_Geometrie_Paper/v2_imref__%s_%s__ws%s.hdr' %(self.win_pos[1],self.win_pos[0],self.ref_win_size_YX[0]))
        #write_envi(im1, '/misc/hy5/scheffler/Projekte/2016_03_Geometrie_Paper/v2_imtgt__%s_%s__ws%s.hdr' %(self.win_pos[1],self.win_pos[0],self.shift_win_size_YX[0]))
        if self.v:
            print('Matching windows with equalized spatial resolution:')
            PLT.subplot_imshow([im0, im1], [os.path.basename(self.path_imref),
                                            os.path.basename(self.path_im2shift)], grid=True)

        gsd_factor = self.imfft_gsd/self.shift_xgsd


        if self.v: print('gsd_factor',         gsd_factor)
        if self.v: print('imfft_gsd_mapvalues',self.imfft_gsd)

        scps = self.calc_shifted_cross_power_spectrum(im0,im1,v=self.v,ignErr=self.ignErr,bin_ws=self.bin_ws)
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
                gdsh_im0,crsp_im1 =self.get_grossly_deshifted_images(im0,im1,x_tempshift,y_tempshift,v=self.v)
                scps = self.calc_shifted_cross_power_spectrum(gdsh_im0,crsp_im1,v=self.v,ignErr=self.ignErr,
                                                              bin_ws=self.bin_ws)

                if scps is None: self.success = False; break
                peakpos = self.get_peakpos(scps)
                x_tempshift,y_tempshift = self.get_shifts_from_peakpos(peakpos,scps.shape)
                x_intshift.append(x_tempshift)
                y_intshift.append(y_tempshift)
                if not self.q: print('Detected integer shifts (X/Y):       %s/%s' %(x_tempshift,y_tempshift))

                if count_iter > self.max_iter:
                    self.success = False
                    if not self.ignErr:
                        raise RuntimeError('No match found in the given window.')
                    else:
                        warnings.warn('No match found in the given window.'); break

            if not self.success==False:
                x_intshift, y_intshift      = sum(x_intshift), sum(y_intshift)
                x_subshift, y_subshift      = self.calc_subpixel_shifts(scps,v=self.v)
                x_totalshift, y_totalshift  = self.get_total_shifts(x_intshift,y_intshift,x_subshift,y_subshift)
                x_shift_px, y_shift_px      = x_totalshift*gsd_factor, y_totalshift*gsd_factor
                if not self.q:
                    print('Detected subpixel shifts (X/Y):      %s/%s' %(x_subshift,y_subshift))
                    print('Calculated total shifts in fft image units (X/Y):         %s/%s' %(x_totalshift,y_totalshift))
                    print('Calculated total shifts in reference image units (X/Y):   %s/%s' %(x_totalshift,y_totalshift))
                    print('Calculated total shifts in target image units (X/Y):      %s/%s' %(x_shift_px,y_shift_px))

                if max([abs(x_totalshift),abs(y_totalshift)]) > self.max_shift:
                    self.success = False
                    if not self.ignErr:
                        raise RuntimeError("The calculated shift is recognized as too large to be valid. "
                        "If you know that it is valid, just set the '-max_shift' parameter to an appropriate value. "
                        "Otherwise try to use a different window size for matching via the '-ws' parameter or define "
                        "the spectral bands to be used for matching manually ('-br' and '-bs'.)")
                else:
                    self.success = True

        self.x_shift_px, self.y_shift_px = (x_shift_px,y_shift_px) if self.success else (None,None)

        warnings.simplefilter('default')

    def correct_shifts(self):
        if self.success:
            if not os.path.exists(os.path.dirname(self.path_out)): os.makedirs(os.path.dirname(self.path_out))
            equal_prj = GEO.get_proj4info(proj=self.ref_prj) == GEO.get_proj4info(proj=self.shift_prj)

            if equal_prj and not self.align_grids and not self.match_gsd and \
                self.out_gsd in [None,[self.shift_xgsd,self.shift_ygsd]]:
                self.shift_image_by_updating_map_info()
            elif equal_prj and self.align_grids: # match_gsd and out_gsd are respected
                self.align_coordinate_grids()
            else: # match_gsd and out_gsd are respected ### TODO: out_proj implementieren
                self.resample_without_grid_aligning()
        else:
            warnings.warn('No result written because detection of image displacements failed.')

    def shift_image_by_updating_map_info(self):
        if not self.q: print('\nWriting output...')
        ds_im2shift = gdal.Open(self.path_im2shift)
        if not ds_im2shift.GetDriver().ShortName == 'ENVI': # FIXME laaangsam
            if self.mp:
                IO.convert_gdal_to_bsq__mp(self.path_im2shift,self.path_out)
            else:
                os.system('gdal_translate -of ENVI %s %s' %(self.path_im2shift, self.path_out))
            file2getHdr = self.path_out
        else:
            shutil.copy(self.path_im2shift,self.path_out)
            file2getHdr = self.path_im2shift
        ds_im2shift = None

        path_hdr = '%s.hdr' %os.path.splitext(file2getHdr)[0]
        path_hdr = '%s.hdr' %file2getHdr if not os.path.exists(path_hdr) else path_hdr
        path_hdr = None if not os.path.exists(path_hdr) else path_hdr

        assert path_hdr, 'No header file found for %s. Applying shifts failed.' %file2getHdr

        new_originY, new_originX = GEO.pixelToMapYX([self.x_shift_px, self.y_shift_px],
                                                    geotransform=self.shift_gt, projection=self.shift_prj)[0]
        map_xshift,  map_yshift  = new_originX - self.shift_gt[0], new_originY - self.shift_gt[3]
        if not self.q: print('Calculated map shifts (X,Y): %s, %s' %(map_xshift,map_yshift))

        with open(path_hdr,'r') as inF:
            content = inF.read()
            map_info_line   = [i for i in content.split('\n') if re.search('map info [\w+]*=',i,re.I) is not None][0]
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
        xgsd, ygsd = (self.ref_xgsd, self.ref_ygsd) if self.match_gsd else self.out_gsd if self.out_gsd is not None \
                        else (self.shift_xgsd,self.shift_ygsd) # self.match_gsd overrides self.out_gsd in __init__

        if not self.q:
            print('Resampling shift image in order to align coordinate grid of shift image to reference image. ')

        is_alignable = lambda gsd1,gsd2: max(gsd1,gsd2)%min(gsd1,gsd2)==0 # checks if pixel sizes are divisible
        if not is_alignable(self.ref_xgsd,xgsd) or not is_alignable(self.ref_ygsd,ygsd) and not self.q:
            print("\nWARNING: The coordinate grids of the reference image and the image to be shifted cannot be "
                  "aligned because their pixel sizes are not exact multiples of each other (ref [X/Y]: "
                  "%s %s; shift [X/Y]: %s %s). Therefore the pixel size of the reference image is chosen for the "
                  "resampled output image. If you don´t like that you can use the '-out_gsd' parameter to set an "
                  "appropriate output pixel size.\n" %(self.ref_xgsd,self.ref_ygsd,xgsd,ygsd))
            xgsd, ygsd = self.ref_xgsd, self.ref_ygsd

        r_xmin,r_xmax,r_ymin,r_ymax =  GEO.corner_coord_to_minmax(GEO.get_corner_coordinates(self.path_imref))
        imref_xgrid = np.arange(r_xmin,r_xmax,self.ref_xgsd)
        imref_ygrid = np.arange(r_ymin,r_ymax,self.ref_ygsd)
        s_xmin,s_xmax,s_ymin,s_ymax =  GEO.corner_coord_to_minmax(GEO.get_corner_coordinates(self.path_im2shift))
        xmin,xmax = UTL.find_nearest(imref_xgrid, s_xmin, 'off',extrapolate=1), UTL.find_nearest(imref_xgrid, s_xmax, 'on',extrapolate=1)
        ymin,ymax = UTL.find_nearest(imref_ygrid, s_ymin, 'off',extrapolate=1), UTL.find_nearest(imref_ygrid, s_ymax, 'on',extrapolate=1)

        new_originY, new_originX = GEO.pixelToMapYX([self.x_shift_px, self.y_shift_px],
                                                    geotransform=self.shift_gt, projection=self.shift_prj)[0]
        map_xshift,  map_yshift  = new_originX - self.shift_gt[0], new_originY - self.shift_gt[3]
        if not self.q: print('Calculated map shifts (X,Y): %s, %s' %(map_xshift,map_yshift))
        gt_shifted                   = list(self.shift_gt[:])
        gt_shifted[0], gt_shifted[3] = new_originX, new_originY

        pathVRT = os.path.splitext(os.path.basename(self.path_im2shift))[0] + '.vrt'
        ds_im2shift = gdal.Open(self.path_im2shift)
        dst_ds  = gdal.GetDriverByName('VRT').CreateCopy(pathVRT, ds_im2shift, 0)
        dst_ds.SetGeoTransform(gt_shifted)
        dst_ds = ds_im2shift  = None

        cmd = "gdalwarp -r cubic -tr %s %s -t_srs '%s' -of ENVI %s %s -overwrite -te %s %s %s %s" \
              %(xgsd,ygsd,self.ref_prj,pathVRT,self.path_out, xmin,ymin,xmax,ymax)
        output = subprocess.check_output(cmd, shell=True)
        if output!=1 and os.path.exists(self.path_out):
            if not self.q: print('Coregistered and resampled image written to %s.' %self.path_out)
        else:
            print(output)
            raise RuntimeError('Resampling failed.')

    def resample_without_grid_aligning(self):
        xgsd, ygsd = (self.ref_xgsd, self.ref_ygsd) if self.match_gsd else self.out_gsd if self.out_gsd is not None \
                        else (self.shift_xgsd,self.shift_ygsd) # self.match_gsd overrides self.out_gsd in __init__

        new_originY, new_originX = GEO.pixelToMapYX([self.x_shift_px, self.y_shift_px],
                                                    geotransform=self.shift_gt, projection=self.shift_prj)[0]
        map_xshift,  map_yshift  = new_originX - self.shift_gt[0], new_originY - self.shift_gt[3]
        if not self.q: print('Calculated map shifts (X,Y): %s, %s' %(map_xshift,map_yshift))
        gt_shifted                   = list(self.shift_gt[:])
        gt_shifted[0], gt_shifted[3] = new_originX, new_originY

        pathVRT = os.path.splitext(os.path.basename(self.path_im2shift))[0] + '.vrt'
        dst_ds  = (lambda drv: drv.CreateCopy(pathVRT, self.ds_im2shift, 0))(gdal.GetDriverByName('VRT'))
        dst_ds.SetGeoTransform(gt_shifted)
        dst_ds  = None

        if not self.q: print('Resampling shift image without aligning of coordinate grids. ')
        cmd = "gdalwarp -r cubic -tr %s %s -t_srs '%s' -of ENVI %s %s -overwrite" \
              %(xgsd,ygsd,self.ref_prj,pathVRT,self.path_out)
        output = subprocess.check_output(cmd, shell=True)
        if output!=1 and os.path.exists(self.path_out):
            if not self.q: print('Coregistered and resampled image written to %s.' %self.path_out)
        else:
            print(output)
            raise RuntimeError('Resampling failed.')



class Geom_Quality_Grid(object):
    def __init__(self, path_im0, path_im1, grid_res, window_size=256, dir_out=None, projectName=None, multiproc=True,
                 r_b4match=1, s_b4match=1, max_iter=5, max_shift=5, data_corners_im0=None,
                 data_corners_im1=None, outFillVal=-9999, nodata=None, calc_corners=1, binary_ws=1,
                 v=0, q=0, ignore_errors=0):
        self.path_imref    = path_im0
        self.path_im2shift = path_im1
        self.dir_out      = dir_out
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

        self.COREG_obj = CoReg(path_im0, path_im1, data_corners_im0=data_corners_im0,
                               data_corners_im1=data_corners_im1, calc_corners=calc_corners, r_b4match=r_b4match,
                               s_b4match=s_b4match, max_iter=max_iter, max_shift=max_shift, nodata=nodata,
                               multiproc=multiproc, binary_ws=self.bin_ws, v=v, q=q, ignore_errors=ignore_errors)
        self.ref_shape                = [self.COREG_obj.ref_rows,self.COREG_obj.ref_cols]
        self.tgt_shape                = [self.COREG_obj.shift_rows,self.COREG_obj.shift_cols]
        self.corner_coord_imref       = self.COREG_obj.corner_coord_imref
        self.corner_coord_im2shift    = self.COREG_obj.corner_coord_im2shift
        self.overlap_poly             = self.COREG_obj.overlap_poly

        self.XY_points, self.XY_mapPoints = self._get_imXY__mapXY_points(self.grid_res)
        self.quality_grid                 = None # set by self.get_quality_grid()

    def _get_imXY__mapXY_points(self,grid_res):
        Xarr,Yarr       = np.meshgrid(np.arange(0,self.tgt_shape[1],grid_res),
                                      np.arange(0,self.tgt_shape[0],grid_res))
        ULmapYX,LRmapYX = GEO.pixelToMapYX([[0,0],[self.tgt_shape[1],self.tgt_shape[0]]], path_im=self.path_im2shift)

        mapXarr,mapYarr = np.meshgrid(np.arange(ULmapYX[1],LRmapYX[1],     self.grid_res*self.COREG_obj.shift_xgsd),
                                      np.arange(ULmapYX[0],LRmapYX[0],-abs(self.grid_res*self.COREG_obj.shift_ygsd)))

        XY_points      = np.empty((Xarr.size,2),Xarr.dtype)
        XY_points[:,0] = Xarr.flat
        XY_points[:,1] = Yarr.flat

        XY_mapPoints      = np.empty((mapXarr.size,2),mapXarr.dtype)
        XY_mapPoints[:,0] = mapXarr.flat
        XY_mapPoints[:,1] = mapYarr.flat

        return XY_points,XY_mapPoints

    def _get_spatial_shifts(self, args):
        pointID,winpos_XY,cornerCoords_l8,cornerCoords_tgt = args
        COREG_obj = CoReg(self.path_imref, self.path_im2shift, wp=winpos_XY, ws=self.window_size,
                          data_corners_im0=cornerCoords_l8, data_corners_im1=cornerCoords_tgt, r_b4match=self.r_b4match,
                          s_b4match=self.s_b4match, max_iter=self.max_iter, max_shift=self.max_shift,
                          nodata=self.nodata, binary_ws=self.bin_ws, v=self.v, q=self.q, ignore_errors=self.ignore_errors)
        COREG_obj.calculate_spatial_shifts()
        res = pointID, COREG_obj.ref_win_size_YX[0], COREG_obj.x_shift_px, COREG_obj.y_shift_px
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
        geomPoints = np.array(XYarr2PointGeom(self.XY_mapPoints[:,0],self.XY_mapPoints[:,1]))

        if GEO.isProjectedOrGeographic(self.COREG_obj.shift_prj)=='geographic':
            crs = dict(ellps='WGS84', datum='WGS84', proj='longlat')
        elif GEO.isProjectedOrGeographic(self.COREG_obj.shift_prj)=='projected':
            UTMzone = abs(GEO.get_UTMzone(prj=self.COREG_obj.shift_prj))
            south   = GEO.get_UTMzone(prj=self.COREG_obj.shift_prj)<0
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

        s2_gsd = self.COREG_obj.shift_xgsd
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
        IO.write_shp(path_out, shapely_points, prj=self.COREG_obj.shift_prj, attrDict=attr_dicts)

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
        IO.write_numpy_to_image(zvalues,path_out,gt=(xmin,grid_res,0,ymax,0,-grid_res),prj=self.COREG_obj.shift_prj)

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
        IO.write_numpy_to_image(zvalues,path_out,gt=(xmin,grid_res,0,ymax,0,-grid_res),prj=self.COREG_obj.shift_prj)

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
    parser.add_argument('-wp', nargs=2, metavar=('X', 'Y'),type=float,help="custom matching window position as  map "\
                        "values in the same projection like the reference image "\
                        "(default: central position of image overlap)", default=None)
    parser.add_argument('-ws', nargs=2, metavar=('X size', 'Y size'),type=float,help="custom matching window size [pixels] (default: 512)", default=512)
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
                        help='no data values for reference image and image to be shifted',default=[None,None])
    parser.add_argument('-calc_cor', nargs=1,type=int, choices=[0,1],default=1, help="calculate true positions of "\
                        "the dataset corners in order to get a useful matching window position within the actual "\
                        "image overlap (default: 1; deactivated if '-cor0' and '-cor1' are given")
    parser.add_argument('-mp', nargs='?', type=int, help='enable multiprocessing (default: 1)', default=1, choices=[0,1])
    parser.add_argument('-bin_ws', nargs='?', type=int, help='use binary X/Y dimensions for the matching window '
                                                             '(default: 1)', default=1, choices=[0, 1])
    parser.add_argument('-v', nargs='?',type=int, help='verbose mode (default: 0)', default=0, choices=[0,1])
    parser.add_argument('-q', nargs='?',type=int, help='quiet mode (default: 0)', default=0, choices=[0,1])
    parser.add_argument('-ignore_errors', nargs='?',type=int, help='Useful for batch processing. (default: 0) '
                        'In case of error COREG.success == False and COREG.x_shift_px/COREG.y_shift_px is None',
                        default=0, choices=[0,1])
    parser.add_argument('--version', action='version', version='%(prog)s 2016-08-02_03')
    args = parser.parse_args()

    print('==================================================================\n'
          '#          SUBPIXEL COREGISTRATION FOR SATELLITE IMAGERY         #\n'
          '#          - algorithm proposed by Foroosh et al. 2002           #\n'
          '#          - python implementation by Daniel Scheffler           #\n'
          '==================================================================\n')

    t0 = time.time()
    COREG_obj = CoReg(args.path_im0,
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
    #COREG_obj.correct_shifts()
    print('\ntotal processing time: %.2fs' %(time.time()-t0))