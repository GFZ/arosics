# -*- coding: utf-8 -*-
from __future__ import (division, print_function,absolute_import) #unicode_literals cause GDAL not to work properly
__author__ = "Daniel Scheffler"
__version__= "2016-03-17_01"


import numpy as np
import os
import shutil
import time
import re
import subprocess
import warnings
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



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

from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import shape

class COREG(object):
    def __init__(self, path_im0,path_im1,path_out='.',r_b4match=1,s_b4match=1,wp=None,ws=512,max_iter=5,max_shift=5,
                 align_grids=0,match_gsd=0,out_gsd=None,data_corners_im0=None,data_corners_im1=None,nodata=[None,None],
                 calc_corners=1,v=0):
        assert os.path.exists(path_im0), "Reference image not found at '%s'. Please check the path."     %path_im0
        assert os.path.exists(path_im1), "Image to be shifted not found at '%s'. Please check the path." %path_im1
        for p,n in zip([path_im0,path_im1,path_out],['reference image','image to be shifted','output image']):
            assert ' ' not in p, "The path of the %s contains whitespaces. This is not supported by GDAL." %n
        self.path_imref               = path_im0
        self.path_im2shift            = path_im1
        print('reference image:', path_im0)
        print('shift image:    ', path_im1, '\n')
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
        if match_gsd and out_gsd is not None: print("WARNING:'-out_gsd' is ignored because '-match_gsd' is set.\n")
        if out_gsd is not None:
            assert isinstance(out_gsd,list) and len(out_gsd)==2, 'out_gsd must be a list with two values.'
        self.corner_coord_imref       = data_corners_im0
        self.corner_coord_im2shift    = data_corners_im1
        if data_corners_im0 is not None and not isinstance(data_corners_im0[0],list): # group if not [[x,y],[x,y]..]
            self.corner_coord_imref = [data_corners_im0[i:i+n] for i in range(0, len(data_corners_im0), 2)]
        if data_corners_im1 is not None and not isinstance(data_corners_im1[0],list): # group if not [[x,y],[x,y]..]
            self.corner_coord_im2shift = [data_corners_im1[i:i+n] for i in range(0, len(data_corners_im1), 2)]
        self.nodata                   = nodata
        self.v                        = v
        if self.v and not os.path.isdir(self.verbose_out): os.makedirs(self.verbose_out)

        self.ds_imref                 = gdal.Open(self.path_imref)
        assert self.ds_imref is not None, 'Reference image can not be read by GDAL. You can try another image format.'
        self.ds_im2shift              = gdal.Open(self.path_im2shift)
        assert self.ds_im2shift is not None, \
            'The image to be shifted can not be read by GDAL. You can try another image format.'
        self.ref_prj, self.shift_prj  = self.ds_imref.GetProjection(),   self.ds_im2shift.GetProjection()
        self.ref_gt , self.shift_gt   = self.ds_imref.GetGeoTransform(), self.ds_im2shift.GetGeoTransform()
        self.ref_xgsd,self.shift_xgsd = abs(self.ref_gt[1]), abs(self.shift_gt[1])
        self.ref_ygsd,self.shift_ygsd = abs(self.ref_gt[5]), abs(self.shift_gt[5])
        assert self.ref_prj                         is not '',   'Reference image has no projection.'
        assert re.search('LOCAL_CS',self.ref_prj)   is None,     'Reference image is not georeferenced.'
        assert self.ref_gt                          is not None, 'Reference image has no map information.'
        assert self.shift_prj                       is not '',   'The image to be shifted has no projection.'
        assert re.search('LOCAL_CS',self.shift_prj) is None,     'The image to be shifted is not georeferenced.'
        assert self.shift_gt                        is not None, 'The image to be shifted has no map information.'
        assert get_proj4info(self.ds_imref) == get_proj4info(self.ds_im2shift),\
            'Input projections are not equal. Different projections are currently not supported.'


        self.imref_band4match        = r_b4match
        self.im2shift_band4match     = s_b4match
        ref_bands, shift_bands = self.ds_imref.RasterCount, self.ds_im2shift.RasterCount
        assert ref_bands >= self.imref_band4match >= 1, "The reference image has %s %s. So '-rb' must be %s%s." \
                % (ref_bands,'bands' if ref_bands>1 else 'band', 'between 1 and ' if ref_bands>1 else '', ref_bands)
        assert shift_bands >= self.im2shift_band4match >= 1, 'The image to be shifted '\
                "has %s %s. So '-sb' must be %s%s." % (shift_bands,'bands' if shift_bands>1 else 'band',
                'between 1 and ' if shift_bands>1 else '', shift_bands)
        if calc_corners:
            if self.corner_coord_imref is None:
                print('Calculating actual data corner coordinates for reference image...')
                self.corner_coord_imref    = get_true_corner_lonlat(self.ds_imref,r_b4match,self.nodata[0])
            if self.corner_coord_im2shift is None:
                print('Calculating actual data corner coordinates for image to be shifted...')
                self.corner_coord_im2shift = get_true_corner_lonlat(self.ds_im2shift,s_b4match,self.nodata[1])
        else:
            self.corner_coord_imref    = get_corner_coordinates(self.ds_imref)
            self.corner_coord_im2shift = get_corner_coordinates(self.ds_im2shift)
        if self.v: print('Corner coordinates of reference image: %s' %self.corner_coord_imref)
        if self.v: print('Corner coordinates of image to be shifted: %s' %self.corner_coord_im2shift)
        self.poly_imref              = get_footprint_polygon(self.corner_coord_imref)
        self.poly_im2shift           = get_footprint_polygon(self.corner_coord_im2shift)

        overlap_tmp                  = get_overlap_polygon(self.poly_imref , self.poly_im2shift, self.v)
        self.overlap_poly            = overlap_tmp['overlap poly'] # has to be the reference projection
        assert self.overlap_poly is not None, 'The input images have no spatial overlap.'
        self.overlap_percentage      = overlap_tmp['overlap percentage']
        self.overlap_area            = overlap_tmp['overlap area']

        if self.v: write_shp(self.poly_imref,   os.path.join(self.verbose_out,'poly_imref.shp'), self.ref_prj)
        if self.v: write_shp(self.poly_im2shift,os.path.join(self.verbose_out,'poly_im2shift.shp'), self.shift_prj)
        if self.v: write_shp(self.overlap_poly, os.path.join(self.verbose_out,'overlap_poly.shp'), self.ref_prj)

        self.win_pos, self.win_size = self.get_opt_winpos_winsize()         if wp is None and ws == 512 else \
                                     [self.get_opt_winpos_winsize()[0],ws]  if wp is None and ws != 512 else \
                                     [list(reversed(wp)),ws] # wp is not None
        assert self.overlap_poly.contains(Point(list(reversed(self.win_pos)))), 'The provided window position is '\
            ' outside of the overlap area of the two input images. Check the coordinates.'
        ### FIXME: transform_mapPt1_to_mapPt2(im2shift_center_map, ds_imref.GetProjection(), ds_im2shift.GetProjection()) # später basteln für den fall, dass projektionen nicht gleich sind

        self.x_shift_px              = None
        self.y_shift_px              = None

    def get_opt_winpos_winsize(self):
        """according to DGM, cloud_mask, trueCornerLonLat"""
        # dummy algorithm: get center position of overlap instead of searching ideal window position in whole overlap
        overlap_center_pos_x, overlap_center_pos_y = self.overlap_poly.centroid.coords.xy
        win_pos = [overlap_center_pos_y[0], overlap_center_pos_x[0]]
        win_size = 512 ### TODO automatischer Algorithmus zur Bestimmung der optimalen Window Position und Groesse
        return win_pos, win_size

    def calculate_spatial_shifts(self):
        im0, im1, imfft_gsd_mapvalues, gsd_factor = get_image_windows_to_match(self.ds_imref,self.ds_im2shift,
            self.win_pos, self.win_size, self.imref_band4match, self.im2shift_band4match, v=self.v)
        if self.v: print('gsd_factor',         gsd_factor)
        if self.v: print('imfft_gsd_mapvalues',imfft_gsd_mapvalues)

        def write_envi(arr,outpath):
            from spectral.io import envi as envi
            shape = (arr.shape[0],arr.shape[1],1) if len(arr.shape)==3 else arr.shape
            out    = envi.create_image(outpath,shape=shape,dtype=arr.dtype,interleave='bsq',ext='.bsq', force=True) # 8bit for muliple masks in one file
            out_mm = out.open_memmap(writable=True)
            out_mm[:,:,0] = arr

        write_envi(im0, '/misc/hy5/scheffler/Projekte/2016_03_Geometrie_Paper/v1_imref__%s_%s__ws%s.hdr' %(self.win_pos[1],self.win_pos[0],'--'))
        write_envi(im1, '/misc/hy5/scheffler/Projekte/2016_03_Geometrie_Paper/v1_imtgt__%s_%s__ws%s.hdr' %(self.win_pos[1],self.win_pos[0],'--'))


        scps = calc_shifted_cross_power_spectrum(im0,im1,v=self.v)
        peakpos = get_peakpos(scps)
        x_intshift, y_intshift  = [], []
        x_tempshift,y_tempshift = get_shifts_from_peakpos(peakpos,scps.shape)
        x_intshift.append(x_tempshift)
        y_intshift.append(y_tempshift)
        print('Detected integer shifts (X/Y):       %s/%s' %(x_tempshift,y_tempshift))

        count_iter = 1
        while [x_tempshift, y_tempshift] != [0,0]:
            count_iter +=1
            print('No clear match found yet. Jumping to iteration %s...' %count_iter)
            if len(x_intshift)>1:
                x_tempshift = x_tempshift + x_intshift[-2]
                y_tempshift = y_tempshift + y_intshift[-2]
            print('input shifts: ', x_tempshift, y_tempshift)
            gdsh_im0,crsp_im1 = get_grossly_deshifted_images(im0,im1,x_tempshift,y_tempshift,scps.shape, v=self.v)
            scps = calc_shifted_cross_power_spectrum(gdsh_im0,crsp_im1,v=self.v)
            peakpos = get_peakpos(scps)
            x_tempshift,y_tempshift = get_shifts_from_peakpos(peakpos,scps.shape)
            x_intshift.append(x_tempshift)
            y_intshift.append(y_tempshift)
            print('Detected integer shifts (X/Y):       %s/%s' %(x_tempshift,y_tempshift))
            if count_iter > self.max_iter:
                raise RuntimeError('No match found in the given window.')
        x_intshift, y_intshift = sum(x_intshift), sum(y_intshift)

        x_subshift, y_subshift = calc_subpixel_shifts(scps,v=self.v)
        print('Detected subpixel shifts (X/Y):      %s/%s' %(x_subshift,y_subshift))

        x_totalshift, y_totalshift = get_total_shifts(x_intshift,y_intshift,x_subshift,y_subshift)
        print('Calculated total shifts in fft image units (X/Y):         %s/%s' %(x_totalshift,y_totalshift))

        self.x_shift_px, self.y_shift_px = x_totalshift*gsd_factor, y_totalshift*gsd_factor
        print('Calculated total shifts in reference image units (X/Y):   %s/%s' %(x_totalshift,y_totalshift))
        if max([x_totalshift,y_totalshift]) > self.max_shift:
            raise RuntimeError("The calculated shift is recognized as too large to be valid. "
            "If you know that it is valid, just set the '-max_shift' parameter to an appropriate value. Otherwise "
            "try to use a different window size for matching via the '-ws' parameter or define the spectral bands "
            "to be used for matching manually ('-br' and '-bs'.)")

    def correct_shifts(self):
        equal_prj = get_proj4info(proj=self.ref_prj)==get_proj4info(proj=self.shift_prj)

        if equal_prj and not self.align_grids and not self.match_gsd and \
            self.out_gsd in [None,[self.shift_xgsd,self.shift_ygsd]]:
            self.shift_image_by_updating_map_info()
        elif equal_prj and self.align_grids: # match_gsd and out_gsd are respected
            self.align_coordinate_grids()
        else: # match_gsd and out_gsd are respected ### TODO: out_proj implementieren
            self.resample_without_grid_aligning()

    def shift_image_by_updating_map_info(self):
        print('\nWriting output...')
        if not self.ds_im2shift.GetDriver().ShortName == 'ENVI':
            os.system('gdal_translate -of ENVI %s %s' %(self.path_im2shift, self.path_out))
            file2getHdr = self.path_out

        else:
            shutil.copy(self.path_im2shift,self.path_out)
            file2getHdr = self.path_im2shift

        path_hdr = '%s.hdr' %os.path.splitext(file2getHdr)[0]
        path_hdr = '%s.hdr' %file2getHdr if not os.path.exists(path_hdr) else path_hdr
        path_hdr = None if not os.path.exists(path_hdr) else path_hdr

        assert path_hdr is not None, 'No header file found for %s. Applying shifts failed.' %file2getHdr

        #path_hdr = '%s.hdr' %os.path.splitext(self.path_im2shift)[0]
        #hdr_dict = (lambda SpFH: SpFH.metadata)(envi.open(path_hdr))
        #with open(path_hdr,'r') as inF:
        #    re_res       = [re.search('([\w+\s\S]*)=',i,re.I) for i in inF.read().split('\n')]
        #    attr_ordered = [i.group(1).strip() for i in re_res if i is not None]
        #print('Original map info:', hdr_dict['map info'])

        new_originY, new_originX = pixelToMapYX([self.x_shift_px,self.y_shift_px],
                                    geotransform=self.shift_gt, projection=self.shift_prj)[0]
        map_xshift,  map_yshift  = new_originX - self.shift_gt[0], new_originY - self.shift_gt[3]
        print('Calculated map shifts (X,Y): %s, %s' %(map_xshift,map_yshift))

        #hdr_dict['map info'][3] = str(float(hdr_dict['map info'][3]) + map_xshift)
        #hdr_dict['map info'][4] = str(float(hdr_dict['map info'][4]) + map_yshift)
        #print('Updated map info:',hdr_dict['map info'])

        #outpath_hdr = '%s.hdr' %os.path.splitext(self.path_out)[0]
        #if not os.path.isdir(os.path.dirname(outpath_hdr)): os.makedirs(os.path.dirname(outpath_hdr))
        #get_valStr = lambda k,v: \
        #    '{\n%s}' % '\n'.join(['  ' + line for line in v.split('\n')])  if k.lower() == 'description' else \
        #    '{ %s }' %(' , '.join([str(i).replace(',', '-') for i in v]),) if not is_str(v) and hasattr(k, '__len__') \
        #    else str(v)
        #with open(outpath_hdr,'w') as outF:
        #    outF.write('ENVI\n')
        #    [outF.write('%s = %s\n' %(k,get_valStr(k,hdr_dict[k]))) for k in attr_ordered]

        #shutil.copy(self.path_im2shift,self.path_out)
        #print('\nCoregistered image written to %s.' %self.path_out)

        with open(path_hdr,'r') as inF:
            content = inF.read()
            map_info_line   = [i for i in content.split('\n') if re.search('map info [\w+]*=',i,re.I) is not None][0]
            map_info_string = re.search('{([\w+\s\S]*)}',map_info_line,re.I).group(1)
            map_info        = [i.strip() for i in map_info_string.split(',')]
        print('Original map info: %s' %map_info)
        new_map_info = map_info
        new_map_info[3] = str(float(map_info[3]) + map_xshift)
        new_map_info[4] = str(float(map_info[4]) + map_yshift)
        print('Updated map info:  %s' %new_map_info)
        new_map_info_line = 'map info = { %s }' %' , '.join(new_map_info)
        new_content = content.replace(map_info_line,new_map_info_line)
        outpath_hdr = '%s.hdr' %os.path.splitext(self.path_out)[0]
        with open(outpath_hdr,'w') as outF:
            outF.write(new_content)

        print('\nCoregistered image written to %s.' %self.path_out)

    def align_coordinate_grids(self):
        xgsd, ygsd = (self.ref_xgsd, self.ref_ygsd) if self.match_gsd else self.out_gsd if self.out_gsd is not None \
                        else (self.shift_xgsd,self.shift_ygsd) # self.match_gsd overrides self.out_gsd in __init__

        print('Resampling shift image in order to align coordinate grid of shift image to reference image. ')

        is_alignable = lambda gsd1,gsd2: max(gsd1,gsd2)%min(gsd1,gsd2)==0 # checks if pixel sizes are divisible
        if not is_alignable(self.ref_xgsd,xgsd) or not is_alignable(self.ref_ygsd,ygsd):
            print("\nWARNING: The coordinate grids of the reference image and the image to be shifted cannot be "\
                    "aligned because their pixel sizes are not exact multiples of each other (ref [X/Y]: "\
                    "%s %s; shift [X/Y]: %s %s). Therefore the pixel size of the reference image is chosen for the "\
                    "resampled output image. If you don´t like that you can use the '-out_gsd' parameter to set an "\
                    "appropriate output pixel size.\n" %(self.ref_xgsd,self.ref_ygsd,xgsd,ygsd))
            xgsd, ygsd = self.ref_xgsd, self.ref_ygsd

        r_xmin,r_xmax,r_ymin,r_ymax =  corner_coord_to_minmax(get_corner_coordinates(self.ds_imref))
        imref_xgrid = np.arange(r_xmin,r_xmax,self.ref_xgsd)
        imref_ygrid = np.arange(r_ymin,r_ymax,self.ref_ygsd)
        s_xmin,s_xmax,s_ymin,s_ymax =  corner_coord_to_minmax(get_corner_coordinates(self.ds_im2shift))
        xmin,xmax = find_nearest(imref_xgrid,s_xmin,'off'),find_nearest(imref_xgrid,s_xmax,'on')
        ymin,ymax = find_nearest(imref_ygrid,s_ymin,'off'),find_nearest(imref_ygrid,s_ymax,'on')

        new_originY, new_originX = pixelToMapYX([self.x_shift_px,self.y_shift_px],
                                    geotransform=self.shift_gt, projection=self.shift_prj)[0]
        map_xshift,  map_yshift  = new_originX - self.shift_gt[0], new_originY - self.shift_gt[3]
        print('Calculated map shifts (X,Y): %s, %s' %(map_xshift,map_yshift))
        gt_shifted                   = list(self.shift_gt[:])
        gt_shifted[0], gt_shifted[3] = new_originX, new_originY

        pathVRT = os.path.splitext(os.path.basename(self.path_im2shift))[0] + '.vrt'
        dst_ds  = (lambda drv: drv.CreateCopy(pathVRT, self.ds_im2shift, 0))(gdal.GetDriverByName('VRT'))
        dst_ds.SetGeoTransform(gt_shifted)
        dst_ds  = None

        cmd = "gdalwarp -r cubic -tr %s %s -t_srs '%s' -of ENVI %s %s -overwrite -te %s %s %s %s" \
              %(xgsd,ygsd,self.ref_prj,pathVRT,self.path_out, xmin,ymin,xmax,ymax)
        output = subprocess.check_output(cmd, shell=True)
        if output!=1 and os.path.exists(self.path_out):
            print('Coregistered and resampled image written to %s.' %self.path_out)
        else:
            print(output)
            raise RuntimeError('Resampling failed.')

    def resample_without_grid_aligning(self):
        xgsd, ygsd = (self.ref_xgsd, self.ref_ygsd) if self.match_gsd else self.out_gsd if self.out_gsd is not None \
                        else (self.shift_xgsd,self.shift_ygsd) # self.match_gsd overrides self.out_gsd in __init__

        new_originY, new_originX = pixelToMapYX([self.x_shift_px,self.y_shift_px],
                                    geotransform=self.shift_gt, projection=self.shift_prj)[0]
        map_xshift,  map_yshift  = new_originX - self.shift_gt[0], new_originY - self.shift_gt[3]
        print('Calculated map shifts (X,Y): %s, %s' %(map_xshift,map_yshift))
        gt_shifted                   = list(self.shift_gt[:])
        gt_shifted[0], gt_shifted[3] = new_originX, new_originY

        pathVRT = os.path.splitext(os.path.basename(self.path_im2shift))[0] + '.vrt'
        dst_ds  = (lambda drv: drv.CreateCopy(pathVRT, self.ds_im2shift, 0))(gdal.GetDriverByName('VRT'))
        dst_ds.SetGeoTransform(gt_shifted)
        dst_ds  = None

        print('Resampling shift image without aligning of coordinate grids. ')
        cmd = "gdalwarp -r cubic -tr %s %s -t_srs '%s' -of ENVI %s %s -overwrite" \
              %(xgsd,ygsd,self.ref_prj,pathVRT,self.path_out)
        output = subprocess.check_output(cmd, shell=True)
        if output!=1 and os.path.exists(self.path_out):
            print('Coregistered and resampled image written to %s.' %self.path_out)
        else:
            print(output)
            raise RuntimeError('Resampling failed.')

def get_corner_coordinates(gdal_ds):
    gdal_ds_GT = gdal_ds.GetGeoTransform()
    ext=[]
    xarr=[0,gdal_ds.RasterXSize]
    yarr=[0,gdal_ds.RasterYSize]
    for px in xarr:
        for py in yarr:
            x=gdal_ds_GT[0]+(px*gdal_ds_GT[1])+(py*gdal_ds_GT[2])
            y=gdal_ds_GT[3]+(px*gdal_ds_GT[4])+(py*gdal_ds_GT[5])
            ext.append([x,y])
        yarr.reverse()
    return ext

def corner_coord_to_minmax(corner_coords):
    x_vals = [int(i[0]) for i in corner_coords]
    y_vals = [int(i[1]) for i in corner_coords]
    xmin,xmax,ymin,ymax = min(x_vals),max(x_vals),min(y_vals),max(y_vals)
    return xmin,xmax,ymin,ymax

def get_footprint_polygon(CornerLonLat):
    outpoly = Polygon(CornerLonLat)
    assert outpoly. is_valid, 'The given coordinates result in an invalid polygon. Check coordinate order.'
    return outpoly

def get_true_corner_lonlat(gdal_ds,bandNr=1,noDataVal=None,v=0):
    rows,cols = gdal_ds.RasterYSize, gdal_ds.RasterXSize
    mask_1bit = np.zeros((rows,cols),dtype='uint8') # zeros -> image area later overwritten by ones

    noDataVal = find_noDataVal(gdal_ds,bandNr) if noDataVal is None else noDataVal
    if noDataVal is None:
        mask_1bit[:,:] = 1
    elif noDataVal=='unclear':
        print("WARNING: No data value could not be automatically detected. Thus the matching window used for shift "
              "calulation had to be centered in the middle of the overlap center without respecting no data values. "
              "To avoid this provide the correct no data values for reference and shift image via '-nodata'")
        mask_1bit[:,:] = 1
    else:
        band_data = gdal_ds.GetRasterBand(bandNr).ReadAsArray()
        mask_1bit[band_data!=noDataVal] = 1

    if v: print('detected no data value',noDataVal)

    cols_containing_data = np.any(mask_1bit, axis=0)
    rows_containing_data = np.any(mask_1bit, axis=1)

    first_dataCol = list(cols_containing_data).index(True)
    last_dataCol  = cols - (list(reversed(cols_containing_data)).index(True) + 1)
    first_dataRow = list(rows_containing_data).index(True)
    last_dataRow  = rows - (list(reversed(rows_containing_data)).index(True) + 1)

    if first_dataCol == 0 and first_dataRow == 0 and last_dataCol == cols-1 and last_dataRow == rows-1 \
            and mask_1bit[0, 0] == True and mask_1bit[0, cols - 1] == True \
            and mask_1bit[rows - 1, 0] == True and mask_1bit[rows - 1, cols - 1] == True:
        UL, UR, LL, LR = (0, 0), (0, cols - 1), (rows - 1, 0), (rows - 1, cols - 1)
    else:
        StartStopRows_in_first_dataCol = [list(mask_1bit[:, first_dataCol]).index(True),
                                         (rows - list(reversed(mask_1bit[:, first_dataCol])).index(True) + 1)]
        StartStopRows_in_last_dataCol  = [list(mask_1bit[:, last_dataCol]).index(True),
                                         (rows - list(reversed(mask_1bit[:, last_dataCol])).index(True) + 1)]
        StartStopCols_in_first_dataRow = [list(mask_1bit[first_dataRow, :]).index(True),
                                         (cols - list(reversed(mask_1bit[first_dataRow, :])).index(True) + 1)]
        StartStopCols_in_last_dataRow  = [list(mask_1bit[last_dataRow, :]).index(True),
                                         (cols - list(reversed(mask_1bit[last_dataRow, :])).index(True) + 1)]

        if True in [abs(np.diff(i)[0]) > 10 for i in [StartStopRows_in_first_dataCol, StartStopRows_in_last_dataCol,
                                                      StartStopCols_in_first_dataRow, StartStopCols_in_last_dataRow]]:
            ''' In case of cut image corners (e.g. ALOS AVNIR-2 Level-1B2): Calculation of trueDataCornerPos outside of the image.'''
            def find_line_intersection_point(line1, line2):
                xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
                ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

                det = lambda a,b: a[0] * b[1] - a[1] * b[0]
                div = det(xdiff, ydiff)
                if div == 0: # lines do not intersect
                    return None, None
                d = (det(*line1), det(*line2))
                x = det(d, xdiff) / div
                y = det(d, ydiff) / div
                return x, y

            line_N = ((StartStopCols_in_first_dataRow[1], first_dataRow),
                      (last_dataCol, StartStopRows_in_last_dataCol[0]))
            line_S = ((first_dataCol, StartStopRows_in_first_dataCol[1]),
                      (StartStopCols_in_last_dataRow[0], last_dataRow))
            line_W = ((StartStopCols_in_first_dataRow[0], first_dataRow),
                      (first_dataCol, StartStopRows_in_first_dataCol[0]))
            line_O = ((last_dataCol, StartStopRows_in_last_dataCol[1]),
                      (StartStopCols_in_last_dataRow[1], last_dataRow))
            corners = [list(reversed(find_line_intersection_point(line_N, line_W))),
                       list(reversed(find_line_intersection_point(line_N, line_O))),
                       list(reversed(find_line_intersection_point(line_S, line_W))),
                       list(reversed(find_line_intersection_point(line_S, line_O)))]
        else:
            dataRow_in_first_dataCol = np.mean(StartStopRows_in_first_dataCol)
            dataRow_in_last_dataCol  = np.mean(StartStopRows_in_last_dataCol)
            dataCol_in_first_dataRow = np.mean(StartStopCols_in_first_dataRow)
            dataCol_in_last_dataRow  = np.mean(StartStopCols_in_last_dataRow)
            corners = [(first_dataRow, dataCol_in_first_dataRow), (dataRow_in_last_dataCol, last_dataCol),
                       (last_dataRow, dataCol_in_last_dataRow), (dataRow_in_first_dataCol, first_dataCol)]

        if not [None, None] in corners:
            distUL = [np.sqrt((0 - c[0]) ** 2 + (0 - c[1]) ** 2) for c in corners]        # distance of each corner to UL
            UL     = corners[distUL.index(min(distUL))]
            distUR = [np.sqrt((0 - c[0]) ** 2 + (cols - c[1]) ** 2) for c in corners]     # distance of each corner to UR
            UR     = corners[distUR.index(min(distUR))]
            distLL = [np.sqrt((rows - c[0]) ** 2 + (0 - c[1]) ** 2) for c in corners]     # distance of each corner to LL
            LL     = corners[distLL.index(min(distLL))]
            distLR = [np.sqrt((rows - c[0]) ** 2 + (cols - c[1]) ** 2) for c in corners]  # distance of each corner to LR
            LR     = corners[distLR.index(min(distLR))]
        else:
            UL, UR, LL, LR = (0, 0), (0, cols - 1), (rows - 1, 0), (rows - 1, cols - 1)
            print('')
            warnings.warn('Automatic detection of the actual image corners of %s failed, because not all of the four '
                          'corners could be found within the image. The current detection algorithm assumes that the '
                          'actual image corners are within the image. Using the outer image corners instead.\n You can '
                          'use -cor0 and -cor1 to manually set the image corner coordinates or -v to check if the '
                            'the matching window is positioned correcly.\n' %os.path.basename(gdal_ds.GetDescription()))

    get_mapYX = lambda YX: pixelToMapYX(list(reversed(YX)),geotransform=gdal_ds.GetGeoTransform(),projection=gdal_ds.GetProjection())[0]
    corner_pos_XY = [list(reversed(i)) for i in [get_mapYX(YX) for YX in [UL, UR, LR, LL]]]
    print(corner_pos_XY)
    return corner_pos_XY

def find_noDataVal(gdal_ds,bandNr=1,sz=3):
    """tries to derive no data value from homogenious corner pixels within 3x3 windows (by default)
    :param gdal_ds:
    :param bandNr:
    :param sz:
    """
    get_mean_std = lambda corner_subset: {'mean':np.mean(corner_subset), 'std':np.std(corner_subset)}
    rows,cols = gdal_ds.RasterYSize,gdal_ds.RasterXSize
    UL = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(0,0,sz,sz))
    UR = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(cols-sz,0,sz,sz))
    LR = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(cols-sz,rows-sz,sz,sz))
    LL = get_mean_std(gdal_ds.GetRasterBand(bandNr).ReadAsArray(0,rows-sz,sz,sz))
    possVals  = [i['mean'] for i in [UL,UR,LR,LL] if i['std']==0]
    # possVals==[]: all corners are filled with data; np.std(possVals)==0: noDataVal clearly identified
    return None if possVals==[] else possVals[0] if np.std(possVals)==0 else 'unclear'

def find_nearest(array,value,round='off'):
    """finds the value of an array nearest to a another single value
    :param array:
    :param value:
    :param round:
    """
    assert round in ['on','off']
    idx = (np.abs(np.array(array)-value)).argmin()
    if round=='off' and array[idx]>value and idx!=0: idx -= 1
    if round=='on'  and array[idx]<value and idx!=len(array)-1: idx += 1
    return array[idx]

def get_overlap_polygon(poly1, poly2,v=0):
    overlap_poly = poly1.intersection(poly2)
    if not overlap_poly.is_empty:
        overlap_percentage = 100 * shape(overlap_poly).area / shape(poly2).area
        if v: print('%.2f percent of the image to be shifted is covered by the reference image.' %overlap_percentage)
        return {'overlap poly':overlap_poly, 'overlap percentage':overlap_percentage, 'overlap area':overlap_poly.area}
    else:
        return {'overlap poly':None, 'overlap percentage':0, 'overlap area':0}

def pixelToMapYX(pixelCoords, path_im=None, geotransform=None, projection=None):
    assert path_im is not None or geotransform is not None and projection is not None, \
        "pixelToMapYX: Missing argument! Please provide either 'path_im' or 'geotransform' AND 'projection'."
    if geotransform is not None and projection is not None:
        gt,proj = geotransform, projection
    else:
        ds       = gdal.Open(path_im)
        gt, proj = ds.GetGeoTransform(), ds.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj)
    # Set up the coordinate transformation object
    ct  = osr.CoordinateTransformation(srs, srs)

    mapYmapXPairs = []
    pixelCoords   = [pixelCoords] if not type(pixelCoords[0]) in [list,tuple] else pixelCoords

    for point in pixelCoords:
        # Translate the pixel pairs into untranslated points
        u_mapX = point[0] * gt[1] + gt[0]
        u_mapY = point[1] * gt[5] + gt[3]
        # Transform the points to the space
        (mapX, mapY, holder) = ct.TransformPoint(u_mapX, u_mapY)
        # Add the point to our return array
        mapYmapXPairs.append([mapY, mapX])

    return mapYmapXPairs

def write_shp(shapely_poly,path_out,prj=None):
    print('Writing %s ...' %path_out)
    if os.path.exists(path_out): os.remove(path_out)
    ds = (lambda drv: drv.CreateDataSource(path_out))(ogr.GetDriverByName("Esri Shapefile"))
    if prj is not None:
        srs = osr.SpatialReference()
        srs.ImportFromWkt(prj)
        layer = ds.CreateLayer('', srs, ogr.wkbPolygon)
    else:
        layer = ds.CreateLayer('', None, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))  # Add one attribute
    defn = layer.GetLayerDefn()

    # Create a new feature (attribute and geometry)
    feat = ogr.Feature(defn)
    feat.SetField('id', 123)

    # Make a geometry, from Shapely object
    geom = ogr.CreateGeometryFromWkb(shapely_poly.wkb)
    feat.SetGeometry(geom)
    layer.CreateFeature(feat)

    # Save and close everything
    ds = layer = feat = geom = None

def get_proj4info(ds=None,proj=None):
    assert ds is not None or proj is not None, "Specify at least one of the arguments 'ds' or 'proj'"
    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjection() if ds is not None else proj)
    return srs.ExportToProj4()

def get_image_windows_to_match(ds_imref,ds_im2shift, win_pos, win_sz, rb, sb, tempAsENVI=0,v=0):
    path_imref         = ds_imref.GetDescription()
    path_im2shift      = ds_im2shift.GetDescription()
    if not tempAsENVI:
        path_imref_clip    = os.path.splitext(os.path.basename(path_imref))[0]+'_clip.vrt'
        path_im2shift_clip = os.path.splitext(os.path.basename(path_im2shift))[0]+'_clip.vrt'
    else:
        tempdir = '/dev/shm/' if os.path.exists('/dev/shm/') else os.path.dirname(path_im2shift)
        path_imref_clip    = os.path.join(tempdir, os.path.splitext(os.path.basename(path_imref))[0] + '_clip.bsq')
        path_im2shift_clip = os.path.join(tempdir, os.path.splitext(os.path.basename(path_im2shift))[0]+'_clip.bsq')
    outFmt = 'VRT' if not tempAsENVI else 'ENVI'

    # auflösungsfaktor
    imref_gt   = ds_imref.GetGeoTransform()
    imref_proj = ds_imref.GetProjection()
    imref_gsd, im2shift_gsd = imref_gt[1], ds_im2shift.GetGeoTransform()[1]
    if v: print('resolutions: ', imref_gsd, im2shift_gsd)

    # win_pos in imref-bildkoordinaten umrechnen
    mapY,mapX,refULy,refULx,refGSDy,refGSDx = win_pos[0],win_pos[1],imref_gt[3],imref_gt[0],imref_gt[5],imref_gt[1]
    overlap_center_in_imref = {'row':(mapY-refULy)/refGSDy,'col':(mapX-refULx)/refGSDx} #mapX = refULx+refULx*cols
    overlap_center_in_imref['row'] = int(overlap_center_in_imref['row']) # int(), um clip-grenze auf pixelgrenzen der referenz zu schieben
    overlap_center_in_imref['col'] = int(overlap_center_in_imref['col'])

    if imref_gsd >= im2shift_gsd:
        """falls referenz die niedrigere auflösung hat:
        im2shift auf imref downsamplen"""
        # 1. Eckkoordinaten Warp-Bild -> per Subset-Read referenz-Bild mit vorgegebener Fenstergröße (wie in späterer fft, um Rechenzeit zu sparen) einlesen -> zwischenspeichern in dev/shm/
        # 2. warp-Bild (bzw. subset je nach fft-fenstergröße) downsamplen auf Auflösung des referenz-bildes per pixel-aggregate -> dev/shm/    ## http://gis.stackexchange.com/questions/110769/gdal-python-aggregate-raster-into-lower-resolution
        # 3. referenz und warp-bild exakt in der fft-fenstergröße einlesen
        imfft_gsd_mapvalues = imref_gsd # fft-Bild bekommt Auflösung des ref-Bilds
        win_sz_in_imref, win_sz_in_im2shift = win_sz, win_sz*imref_gsd/im2shift_gsd
        get_maxposs_sz = lambda ds: min([ds.RasterXSize,ds.RasterYSize])
        if win_sz_in_imref > get_maxposs_sz(ds_imref):
            print('WARNING: The given target window size exceeds the extent of the reference image. '
                  'It has been reduced to the maximum possible size.')
            win_sz = get_maxposs_sz(ds_imref)
        if win_sz_in_im2shift > get_maxposs_sz(ds_im2shift):
            print('WARNING: The given target window size equals %s px in the image to be shifted which exceeds '\
                  'its extent. Thus it has been reduced to the maximum possible size.' %win_sz_in_im2shift)
            win_sz = int(get_maxposs_sz(ds_im2shift)/(imref_gsd/im2shift_gsd))
        clip_sz = win_sz

        # um overlap_center_in_imref eine box erstellen und imref per subset-read einlesen -> im0
        imref_box = [[overlap_center_in_imref['row']-clip_sz//2, overlap_center_in_imref['col']-clip_sz//2],
                     [overlap_center_in_imref['row']+clip_sz//2, overlap_center_in_imref['col']+clip_sz//2]] # [UL,LR]
        imref_clip_rowS, imref_clip_colS = imref_box[0] # UL
        imref_clip_data = ds_imref.GetRasterBand(rb).ReadAsArray(imref_clip_colS, imref_clip_rowS, clip_sz, clip_sz)

        # map-koordinaten der clip-box von imref ausrechnen
        imref_box_xy  = [[int(i[1]),int(i[0])] for i in imref_box] # swap X and Y values
        imref_box_map = pixelToMapYX(imref_box_xy, projection=imref_proj, geotransform=imref_gt) # [[4455560.0, 299060.0], [4450440.0, 293940.0]] [[HW,RW],[HW,RW]]
        map_xvals = [i[1] for i in imref_box_map] # UTM-RW
        map_yvals = [i[0] for i in imref_box_map] # UTM-HW

        # gdalwarp für im2shift mit pixel aggregate
        rsp_algor = 'average'
        if gdal.VersionInfo().startswith('1'):
            print("WARNING: The GDAL version on this server does not yet support the resampling algorithm "\
                  "'average'. This can affect the correct detection of subpixel shifts. To avoid this "\
                  "please update GDAL to a version above 2.0.0!")
            rsp_algor = 'cubic'
        cmd = "gdalwarp -r %s -tr %s %s -t_srs '%s' -of %s -te %s %s %s %s %s %s -overwrite%s" \
               %(rsp_algor,imfft_gsd_mapvalues, imfft_gsd_mapvalues, imref_proj, outFmt, min(map_xvals),\
               min(map_yvals),max(map_xvals),max(map_yvals),path_im2shift,path_im2shift_clip,' -q' if not v else '')
                # te_srs is not neccessary because -t_srs = imref_proj and output extent is derived from imref
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        if output==1: print(output) # in case of error

        # im1 daten in array einlesen
        ds_im2shift_clip  = gdal.OpenShared(path_im2shift_clip) if not tempAsENVI else gdal.Open(path_im2shift_clip)
        im2shift_clip_data= ds_im2shift_clip.GetRasterBand(sb).ReadAsArray()
        [os.remove(p) for p in [path_im2shift_clip, os.path.splitext(path_im2shift_clip)[0]+'.hdr', \
            path_im2shift_clip +'.aux.xml'] if tempAsENVI and os.path.exists(p)]

    else:
        """falls referenz die höhere auflösung hat"""
        """imref auf im2shift downsamplen"""  ## funktioniert perfekt (getestet an UTM_Z43:5m - UTM_Z43:30m)
        # 1. Eckkoordinaten des Clip-Windows für imref berechnen (anhand optimaler fft-Fenstergröße und Zentrumspunkt der Bildüberlappung)
        # 2. imref ausschneiden und per pixel-aggregate auf auflösung von im2shift downsamplen -> dev/shm/ -> einlesen als array ## http://gis.stackexchange.com/questions/110769/gdal-python-aggregate-raster-into-lower-resolution
        # 3. warp-bild anhand des Clip-Windows des refernz-bildes ausschneiden und auflösung unverändert lassen

        imfft_gsd_mapvalues  = im2shift_gsd # fft-Bild bekommt Auflösung des warp-Bilds
        im2shift_max_poss_sz = min([ds_im2shift.RasterXSize, ds_im2shift.RasterYSize])
        win_sz               = win_sz if not win_sz > im2shift_max_poss_sz else im2shift_max_poss_sz
        # clip muss größer sein als Ziel-fft_bild, weil clip erst ausgeschnitten und dann gedreht wird
        clip_sz              = int(1.5 * win_sz * im2shift_gsd/imref_gsd) # clip-size in referenz (1.5 wegen evtl. verdrehung bei unterschiedlichen projektionen)
        imref_max_poss_sz = min([ds_imref.RasterXSize, ds_imref.RasterYSize])
        clip_sz = clip_sz if not clip_sz > imref_max_poss_sz else win_sz * im2shift_gsd/imref_gsd
        if clip_sz > imref_max_poss_sz:
            print('WARNING: The given target window size equals %s px in the image to be shifted which exceeds '\
                  'its extent. Thus it has been reduced to the maximum possible size.' %clip_sz)
            clip_sz = imref_max_poss_sz

        # um overlap_center_in_imref eine box erstellen und map-koordinaten der clip-box von imref ausrechnen
        # (um downsampling per gdalwarp rechnen zu können -> map koordinaten nötig)
        imref_box = [[overlap_center_in_imref['row']-clip_sz//2, overlap_center_in_imref['col']-clip_sz//2],
                     [overlap_center_in_imref['row']+clip_sz//2, overlap_center_in_imref['col']+clip_sz//2]] # [UL,LR]
        imref_box_xy  = [[int(i[1]),int(i[0])] for i in imref_box] # swap X and Y values
        imref_box_map = pixelToMapYX(imref_box_xy, projection=imref_proj, geotransform=imref_gt) # [[4455560.0, 299060.0], [4450440.0, 293940.0]] [[HW,RW],[HW,RW]]
        map_xvals = [i[1] for i in imref_box_map] # UTM-RW
        map_yvals = [i[0] for i in imref_box_map] # UTM-HW

        # imref mit imref_box ausclippen und per pixel aggregate auf auflösung von im2shift downsamplen
        rsp_algor = 'average'
        if gdal.VersionInfo().startswith('1'):
            print("WARNING: The GDAL version on this server does not yet support the resampling algorithm "\
                      "'average'. This can affect the correct detection of subpixel shifts. To avoid this "\
                      "please update GDAL to a version above 2.0.0!")
            rsp_algor = 'cubic'
        cmd = "gdalwarp -r %s -tr %s %s -t_srs '%s' -of %s -te %s %s %s %s %s %s -overwrite%s" \
                %(rsp_algor,imfft_gsd_mapvalues, imfft_gsd_mapvalues, imref_proj, outFmt, min(map_xvals),\
                  min(map_yvals),max(map_xvals),max(map_yvals),path_imref,path_imref_clip,' -q' if not v else '')
                # te_srs is not neccessary because -t_srs = imref_proj and output extent is derived from imref
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        if output==1: print(output) # in case of error

        # imref_clip einlesen
        ds_im2ref_clip = gdal.OpenShared(path_imref_clip) if not tempAsENVI else gdal.Open(path_imref_clip)
        imref_clip_data= ds_im2ref_clip.GetRasterBand(sb).ReadAsArray()
        ds_im2ref_clip.FlushCache()
        ds_im2ref_clip = None
        [os.remove(p) for p in [path_imref_clip, os.path.splitext(path_imref_clip)[0]+'.hdr', \
            path_imref_clip +'.aux.xml'] if tempAsENVI and os.path.exists(p)]

        # im2shift mit imref_box_map ausclippen (cubic, da pixelgrenzen verschoben werden müssen aber auflösung
        # gleich bleibt)
        cmd = "gdalwarp -r cubic -tr %s %s -t_srs '%s' -of %s -te %s %s %s %s %s %s -overwrite%s" \
                %(imfft_gsd_mapvalues, imfft_gsd_mapvalues, imref_proj, outFmt, min(map_xvals),min(map_yvals), \
                  max(map_xvals),max(map_yvals),path_im2shift,path_im2shift_clip,' -q' if not v else '')
                # te_srs is not neccessary because -t_srs = imref_proj and output extent is derived from imref
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        if output==1: print(output) # in case of error

        # im2shift_clip einlesen
        ds_im2shift_clip  = gdal.OpenShared(path_im2shift_clip) if not tempAsENVI else gdal.Open(path_im2shift_clip)
        im2shift_clip_data= ds_im2shift_clip.GetRasterBand(sb).ReadAsArray()
        [os.remove(p) for p in [path_im2shift_clip, os.path.splitext(path_im2shift_clip)[0]+'.hdr', \
            path_im2shift_clip +'.aux.xml'] if tempAsENVI and os.path.exists(p)]


    ds_im2shift_clip.FlushCache()
    ds_im2shift_clip  = None
    gsd_factor        = imfft_gsd_mapvalues/im2shift_gsd
    return imref_clip_data, im2shift_clip_data, imfft_gsd_mapvalues, gsd_factor

def get_opt_fftw_winsize(im_shape, target_size=None,v=0):
    binarySizes = [2**i for i in range(5,14)] # [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
    possibSizes = [i for i in binarySizes if i <= min(im_shape)]
    assert possibSizes != [], 'The matching window became too small for calculating a reliable match. Matching failed.'
    target_size = max(possibSizes) if target_size is None else target_size
    closest_to_target = min(possibSizes, key=lambda x:abs(x-target_size))
    if v: print('final window size: %s' %closest_to_target)
    return closest_to_target

def calc_shifted_cross_power_spectrum(im0,im1,window_size=1024,precision=np.complex64, v=0):
    """Calculates shifted cross power spectrum for quantifying x/y-shifts.
        :param im0:         reference image
        :param im1:         subject image to shift
        :param window_size: size of image area to be processed
        :param precision:   to be quantified as a datatype
        :param v:           verbose
        :return:            2D-numpy-array of the shifted cross power spectrum
    """
    window_size = get_opt_fftw_winsize(im0.shape, target_size=window_size,v=v)

    t0 = time.time()
    in_arr0  = im0[:window_size,:window_size].astype(precision)
    in_arr1  = im1[:window_size,:window_size].astype(precision)
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
        subplot_imshow([in_arr0.astype(np.uint16),in_arr1.astype(np.uint16),fft_arr0.astype(np.uint8),
                        fft_arr1.astype(np.uint8), scps], titles=['matching window im0', 'matching window im1',
                        "fft result im0", "fft result im1", "cross power spectrum"])
        subplot_3dsurface(scps.astype(np.float32))
    return scps

def get_peakpos(scps):
    max_flat_idx = np.argmax(scps)
    return np.unravel_index(max_flat_idx, scps.shape)

def wfa(p,c):
    try:
        with open(p,'a') as of: of.write(c)
    except: pass

def get_shifts_from_peakpos(peakpos,arr_shape):
    y_shift = peakpos[0]-arr_shape[0]//2
    x_shift = peakpos[1]-arr_shape[1]//2
    return x_shift,y_shift

def get_grossly_deshifted_im0(im0,x_intshift,y_intshift,shape_power_spectrum):
    old_center_coord = [i//2 for i in im0.shape]
    new_center_coord = [old_center_coord[0]+y_intshift, old_center_coord[1]+x_intshift]
    x_left  = new_center_coord[1]
    x_right = im0.shape[1]-new_center_coord[1]
    y_above = new_center_coord[0]
    y_below = im0.shape[0]-new_center_coord[0]
    maxposs_winsz = 2*min(x_left,x_right,y_above,y_below)
    wsz = maxposs_winsz if maxposs_winsz <=min(shape_power_spectrum) and shape_power_spectrum is not None \
        else min(shape_power_spectrum)
    #wsz = get_opt_fftw_winsize(maxposs_winsz,(maxposs_winsz,maxposs_winsz))
    return im0[new_center_coord[0]-wsz//2:new_center_coord[0]+wsz//2+1, \
           new_center_coord[1]-wsz//2:new_center_coord[1]+wsz//2+1]

def get_corresponding_im1_clip(im1, shape_im0):
    center_coord = [i//2 for i in im1.shape]
    return im1[center_coord[0]-shape_im0[0]//2:center_coord[0]+shape_im0[0]//2+1, \
           center_coord[1]-shape_im0[1]//2:center_coord[1]+shape_im0[1]//2+1]

def get_grossly_deshifted_images(im0,im1,x_intshift,y_intshift,shape_power_spectrum=None,v=0):
    gdsh_im0 = get_grossly_deshifted_im0(im0,x_intshift,y_intshift,shape_power_spectrum)
    crsp_im1 = get_corresponding_im1_clip(im1,gdsh_im0.shape)
    if v:
        subplot_imshow([get_corresponding_im1_clip(im0,gdsh_im0.shape),get_corresponding_im1_clip(im1,gdsh_im0.shape)],
                       titles=['reference original', 'target'], grid=True)
        subplot_imshow([gdsh_im0, crsp_im1],titles=['reference virtually shifted', 'target'], grid=True)
    return gdsh_im0,crsp_im1

def find_side_maximum(scps, v=0):
    centerpos = [scps.shape[0]//2, scps.shape[1]//2]
    profile_left  = scps[ centerpos [0]  ,:centerpos[1]+1]
    profile_right = scps[ centerpos [0]  , centerpos[1]:]
    profile_above = scps[:centerpos [0]+1, centerpos[1]]
    profile_below = scps[ centerpos [0]: , centerpos[1]]

    if v:
        max_count_vals = 10
        subplot_2dline([[range(len(profile_left)) [-max_count_vals:], profile_left[-max_count_vals:]],
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

def calc_subpixel_shifts(scps,v=0):
    sidemax_lr, sidemax_ab = find_side_maximum(scps,v)
    x_subshift = (sidemax_lr['direction_factor']*sidemax_lr['value'])/(np.max(scps)+sidemax_lr['value'])
    y_subshift = (sidemax_ab['direction_factor']*sidemax_ab['value'])/(np.max(scps)+sidemax_ab['value'])
    return x_subshift, y_subshift

def get_total_shifts(x_intshift,y_intshift,x_subshift,y_subshift):
    return x_intshift+x_subshift, y_intshift+y_subshift

def subplot_2dline(XY_tuples, titles=None, shapetuple=None, grid=False):
    shapetuple = (1,len(XY_tuples)) if shapetuple is None else shapetuple
    assert titles is None or len(titles)==len(XY_tuples), \
        'List in titles keyword must have the same length as the passed XY_tuples.'
    norm = lambda array, normto: [float(i)*(normto/max(array)) for i in array]
    fig = plt.figure(figsize=norm(plt.figaspect(shapetuple[0]/shapetuple[1]*1.), 10))
    for i,XY in enumerate(XY_tuples):
        ax    = fig.add_subplot(shapetuple[0], shapetuple[1], i+1)
        X,Y = XY
        ax.plot(X,Y,linestyle='-')
        if titles is not None: ax.set_title(titles[i])
        if grid: ax.grid(which='major', axis='both', linestyle='-')
    plt.tight_layout()
    plt.show(); time.sleep(5); plt.close()

def subplot_imshow(ims, titles=None, shapetuple=None, grid=False):
    ims        = [ims] if not isinstance(ims,list) else ims
    assert titles is None or len(titles)==len(ims), 'Error: Got more or less titles than images.'
    shapetuple = (1,len(ims)) if shapetuple is None else shapetuple
    norm       = lambda array, normto: [float(i)*(normto/max(array)) for i in array]
    fig,axes   = plt.subplots(shapetuple[0],shapetuple[1], figsize=norm(plt.figaspect(shapetuple[0]/shapetuple[1]*1.), 20))
    [axes[i].imshow(im,cmap='binary',interpolation='none', vmin=np.percentile(im,2), vmax = np.percentile(im,98)) \
        for i,im in enumerate(ims)]
    if titles is not None: [axes[i].set_title(titles[i]) for i in range(len(ims))]
    if grid: [axes[i].grid(which='major', axis='both', linestyle='-') for i in range(len(ims))]
    plt.tight_layout()
    plt.show(); time.sleep(5); plt.close()

def subplot_3dsurface(ims,shapetuple=None):
    ims        = [ims] if not isinstance(ims,list) else ims
    shapetuple = (1,len(ims)) if shapetuple is None else shapetuple
    norm       = lambda array, normto: [float(i)*(normto/max(array)) for i in array]
    fig        = plt.figure(figsize=norm(plt.figaspect((shapetuple[0]/2.)/shapetuple[1]*1.), 20))
    for i,im in enumerate(ims):
        ax    = fig.add_subplot(shapetuple[0], shapetuple[1], i+1, projection='3d')
        x     = np.arange(0, im.shape[0], 1)
        y     = np.arange(0, im.shape[1], 1)
        X, Y  = np.meshgrid(x, y)
        Z     = im.reshape(X.shape)
        ax.plot_surface(X, Y, Z, cmap = plt.cm.hot)
        ax.contour(X, Y, Z, zdir='x', cmap=plt.cm.coolwarm, offset=0)
        ax.contour(X, Y, Z, zdir='y', cmap=plt.cm.coolwarm, offset=im.shape[1])
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    plt.tight_layout()
    plt.show();time.sleep(5);plt.close()

if __name__ == '__main__':
    import argparse
    from socket import gethostname; from datetime import datetime as dt; from getpass import getuser; import sys
    wfa('/misc/hy5/scheffler/tmp/crlf','%s\t%s\t%s\t%s\n' %(dt.now(), getuser(),gethostname(),' '.join(sys.argv)))
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
        "The following non-standard Python libraries are required: gdal, osr, ogr "\
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
    parser.add_argument('-ws', nargs='?', type=int,help="custom matching window size [pixels] (default: 512)", default=512)
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
    parser.add_argument('-v', nargs='?',type=int, help='verbose mode (default: 0)', default=0, choices=[0,1])
    parser.add_argument('--version', action='version', version='%(prog)s 2016-03-17_01')
    args = parser.parse_args()

    print('==================================================================\n'
          '#          SUBPIXEL COREGISTRATION FOR SATELLITE IMAGERY         #\n'
          '#          - algorithm proposed by Foroosh et al. 2002           #\n'
          '#          - python implementation by Daniel Scheffler           #\n'
          '==================================================================\n')

    t0 = time.time()
    COREG_obj = COREG(args.path_im0,args.path_im1,args.path_out,args.br,args.bs, args.wp,args.ws,args.max_iter,
                      args.max_shift,args.align_grids,args.match_gsd,args.out_gsd,args.cor0,args.cor1,args.nodata,
                      args.calc_cor,args.v)
    COREG_obj.calculate_spatial_shifts()
    COREG_obj.correct_shifts()
    print('\ntotal processing time: %.2fs' %(time.time()-t0))


