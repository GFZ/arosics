# -*- coding: utf-8 -*-
__author__='Daniel Scheffler'

import warnings
import os

# custom
try:
    import gdal
except ImportError:
    from osgeo import gdal
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .Geom_Quality_Grid import Geom_Quality_Grid
from .CoReg             import COREG
from .DeShifter         import DESHIFTER
from py_tools_ds.ptds.geo.coord_trafo import transform_any_prj, reproject_shapelyGeometry
from py_tools_ds.ptds                 import GeoArray



class COREG_LOCAL(object):
    def __init__(self, im_ref, im_tgt, grid_res, window_size=(256,256), path_out=None, fmt_out='ENVI', projectDir=None,
                 r_b4match=1, s_b4match=1, max_iter=5, max_shift=5, data_corners_im0=None,
                 data_corners_im1=None, outFillVal=-9999, nodata=(None,None), calc_corners=True, binary_ws=True,
                 CPUs=None, progress=True, v=False, q=False):

        """Applies the algorithm to detect spatial shifts to the whole overlap area of the input images. Spatial shifts
        are calculated for each point in grid of which the parameters can be adjusted using keyword arguments. Shift
        correction performs a polynomial transformation using the calculated shifts of each point in the grid as GCPs.
        Thus this class can be used to correct for locally varying geometric distortions of the target image.

        :param im_ref(str, GeoArray):   source path of reference image (any GDAL compatible image format is supported)
        :param im_tgt(str, GeoArray):   source path of image to be shifted (any GDAL compatible image format is supported)
        :param grid_res:                quality grid resolution in pixels of the target image
        :param window_size(tuple):      custom matching window size [pixels] (default: (256,256))
        :param path_out(str):           target path of the coregistered image
                                            - if None (default), no output is written to disk
                                            - if 'auto': /dir/of/im1/<im1>__shifted_to__<im0>.bsq
        :param fmt_out(str):            raster file format for output file. ignored if path_out is None. can be any GDAL
                                        compatible raster file format (e.g. 'ENVI', 'GeoTIFF'; default: ENVI)
        :param projectDir:
        :param r_b4match(int):          band of reference image to be used for matching (starts with 1; default: 1)
        :param s_b4match(int):          band of shift image to be used for matching (starts with 1; default: 1)
        :param max_iter(int):           maximum number of iterations for matching (default: 5)
        :param max_shift(int):          maximum shift distance in reference image pixel units (default: 5 px)
        :param data_corners_im0(list):  map coordinates of data corners within reference image
        :param data_corners_im1(list):  map coordinates of data corners within image to be shifted
        :param outFillVal(int):         if given the generated geometric quality grid is filled with this value in case
                                        no match could be found during co-registration (default: -9999)
        :param nodata(tuple):           no data values for reference image and image to be shifted
        :param calc_corners(bool):      calculate true positions of the dataset corners in order to get a useful
                                        matching window position within the actual image overlap
                                        (default: True; deactivated if 'data_corners_im0' and 'data_corners_im1' are given
        :param binary_ws(bool):         use binary X/Y dimensions for the matching window (default: True)
        :param CPUs(int):               number of CPUs to use during calculation of geometric quality grid
                                        (default: None, which means 'all CPUs available')
        :param progress(bool):          show progress bars (default: True)
        :param v(bool):                 verbose mode (default: False)
        :param q(bool):                 quiet mode (default: False)
        """
        # TODO outgsd, matchgsd
        # assertions
        assert fmt_out

        self.params = dict([x for x in locals().items() if x[0] != "self" and not x[0].startswith('__')])

        self.imref        = im_ref if isinstance(im_ref, GeoArray) else GeoArray(im_ref)
        self.im2shift     = im_tgt if isinstance(im_tgt, GeoArray) else GeoArray(im_tgt)

        self.path_out     = path_out  # updated by self.set_outpathes
        self.fmt_out      = fmt_out
        self._projectDir  = projectDir
        self.grid_res     = grid_res
        self.window_size  = window_size
        self.max_shift    = max_shift
        self.max_iter     = max_iter
        self.calc_corners = calc_corners
        self.nodata       = nodata
        self.outFillVal   = outFillVal
        self.bin_ws       = binary_ws
        self.CPUs         = CPUs
        self.path_verbose_out = '' # TODO
        self.v            = v
        self.q            = q if not v else False        # overridden by v
        self.progress     = progress if not q else False # overridden by v

        COREG.__dict__['_set_outpathes'](self, self.imref, self.im2shift)
        # make sure that the output directory of coregistered image is the project directory if a project directory is given
        if path_out and projectDir and os.path.basename(self.path_out):
            self.path_out = os.path.join(self.projectDir, os.path.basename(self.path_out))

        gdal.AllRegister()

        self.COREG_obj = COREG(self.imref, self.im2shift,
                               ws               = window_size,
                               data_corners_im0 = data_corners_im0,
                               data_corners_im1 = data_corners_im1,
                               calc_corners     = calc_corners,
                               r_b4match        = r_b4match,
                               s_b4match        = s_b4match,
                               max_iter         = max_iter,
                               max_shift        = max_shift,
                               nodata           = nodata,
                               multiproc        = self.CPUs is None or self.CPUs > 1,
                               binary_ws        = self.bin_ws,
                               v                = v,
                               q                = q,
                               ignore_errors    = True)

        self._quality_grid      = None # set by self.quality_grid
        self._CoRegPoints_table = None # set by self.CoRegPoints_table
        self._coreg_info        = None # set by self.coreg_info
        self.deshift_results    = None # set by self.correct_shifts()
        self._success           = None # set by self.success property


    @property
    def projectDir(self):
        if self._projectDir:
            if len(os.path.split(self._projectDir))==1:
                return os.path.abspath(os.path.join(os.path.curdir, self._projectDir))
            else:
                return os.path.abspath(self._projectDir)
        else:
            # return a project name that not already has a corresponding folder on disk
            projectDir = os.path.join(os.path.dirname(self.im2shift.filePath), 'UntitledProject_1')
            while os.path.isdir(projectDir):
                projectDir = '%s_%s' % (projectDir.split('_')[0], int(projectDir.split('_')[1]) + 1)
            self._projectDir = projectDir
            return self._projectDir


    @property
    def quality_grid(self):
        if self._quality_grid:
            return self._quality_grid
        else:
            self._quality_grid = Geom_Quality_Grid(self.COREG_obj, self.grid_res, outFillVal=self.outFillVal,
                                                   dir_out=self.projectDir, CPUs=self.CPUs, progress=self.progress,
                                                   v=self.v, q=self.q)
            if self.v:
                self.view_CoRegPoints(figsize=(10,10))
            return self._quality_grid


    @property
    def CoRegPoints_table(self):
        return self.quality_grid.CoRegPoints_table


    @property
    def success(self):
        self._success = self.quality_grid.GCPList != []
        if not self._success and not self.q:
            warnings.warn('No valid GCPs could by identified.')
        return self._success


    def view_CoRegPoints(self, attribute2plot='ABS_SHIFT', cmap=None, exclude_fillVals=True, backgroundIm='tgt',
                         figsize=None, savefigPath='', savefigDPI=96):
        """Shows a map of the calculated quality grid with the target image as background.

        :param attribute2plot:      <str> the attribute of the quality grid to be shown (default: 'ABS_SHIFT')
        :param cmap:                <plt.cm.<colormap>> a custom color map to be applied to the plotted grid points
                                                        (default: 'RdYlGn_r')
        :param exclude_fillVals:    <bool> whether to exclude those points of the grid where spatial shift detection failed
        :param backgroundIm:        <str> whether to use the target or the reference image as map background. Possible
                                          options are 'ref' and 'tgt' (default: 'tgt')
        :param figsize:             <tuple> size of the figure to be viewed, e.g. (10,10)
        :param savefigPath:
        :param savefigDPI:
        :return:
        """

        # get a map showing target image
        if backgroundIm not in ['tgt','ref']: raise ValueError('backgroundIm')
        backgroundIm      = self.im2shift if backgroundIm=='tgt' else self.imref
        fig, ax, map2show = backgroundIm.show_map(figsize=figsize, nodataVal=self.nodata[1], return_map=True)
        # fig, ax, map2show = backgroundIm.show_map_utm(figsize=(20,20), nodataVal=self.nodata[1], return_map=True)
        plt.title(attribute2plot)

        # transform all points of quality grid to LonLat
        GDF = self.CoRegPoints_table.loc[self.CoRegPoints_table.X_SHIFT_M != self.outFillVal, ['geometry', attribute2plot]].copy() \
                if exclude_fillVals else self.CoRegPoints_table.loc[:, ['geometry', attribute2plot]]

        # get LonLat coordinates for all points
        get_LonLat    = lambda X, Y: transform_any_prj(self.im2shift.projection, 4326, X, Y)
        GDF['LonLat'] = list(GDF['geometry'].map(lambda geom: get_LonLat(*tuple(np.array(geom.coords.xy)[:,0]))))

        # get colors for all points
        #vmin = min(GDF[GDF[attribute2plot] != self.outFillVal][attribute2plot])
        #vmax = max(GDF[GDF[attribute2plot] != self.outFillVal][attribute2plot])
        #norm = mpl_normalize(vmin=vmin, vmax=vmax)
        palette = cmap if cmap else plt.cm.RdYlGn_r
        #GDF['color'] = [*GDF[attribute2plot].map(lambda val: palette(norm(val)))]

        # add quality grid to map
        #plot_point = lambda row: ax.plot(*map2show(*row['LonLat']), marker='o', markersize=7.0, alpha=1.0, color=row['color'])
        #GDF.apply(plot_point, axis=1)
        GDF['plt_XY'] = list(GDF['LonLat'].map(lambda ll: map2show(*ll)))
        GDF['plt_X']  = list(GDF['plt_XY'].map(lambda XY: XY[0]))
        GDF['plt_Y']  = list(GDF['plt_XY'].map(lambda XY: XY[1]))
        vmin, vmax = np.percentile(GDF[attribute2plot], 0), np.percentile(GDF[attribute2plot], 95)
        points = plt.scatter(GDF['plt_X'],GDF['plt_Y'], c=GDF[attribute2plot],
                             cmap=palette, marker='o' if len(GDF)<10000 else '.', s=50, alpha=1.0,
                             vmin=vmin, vmax=vmax)

        # add colorbar
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size="2%", pad=0.1) # create axis on the right; size =2% of ax; padding = 0.1 inch
        plt.colorbar(points, cax=cax)

        plt.show(block=True)

        if savefigPath:
            fig.savefig(savefigPath, dpi=savefigDPI)


    def view_CoRegPoints_folium(self, attribute2plot='ABS_SHIFT', cmap=None, exclude_fillVals=True):
        warnings.warn(UserWarning('This function is still under construction and may not work as expected!'))
        assert self.CoRegPoints_table is not None, 'Calculate quality grid first!'

        try:
            import folium, geojson
            from folium import plugins
        except ImportError:
            folium, geojson, plugins = [None]*3
        if not folium or not geojson:
            raise ImportError("This method requires the libraries 'folium' and 'geojson'. They can be installed with "
                              "the shell command 'pip install folium geojson'.")

        lon_min, lat_min, lon_max, lat_max = \
            reproject_shapelyGeometry(self.im2shift.box.mapPoly, self.im2shift.projection, 4326).bounds
        center_lon, center_lat = (lon_min+lon_max)/2, (lat_min+lat_max)/2

        # get image to plot
        image2plot = self.im2shift[0] # FIXME hardcoded band

        from py_tools_ds.ptds.geo.raster.reproject import warp_ndarray
        image2plot, gt, prj = warp_ndarray(image2plot, self.im2shift.geotransform, self.im2shift.projection,
                                           in_nodata=self.nodata[1], out_nodata=self.nodata[1],
                                           out_XYdims=(1000, 1000), q=True, out_prj='epsg:3857') # image must be transformed into web mercator projection

        # create map
        map_osm = folium.Map(location=[center_lat, center_lon])#,zoom_start=3)
        import matplotlib
        plugins.ImageOverlay(
            colormap=lambda x: (1, 0, 0, x), # TODO a colormap must be given
            # colormap=matplotlib.cm.gray, # does not work
            image=image2plot, bounds=[[lat_min, lon_min], [lat_max, lon_max]],
            ).add_to(map_osm)

        folium.GeoJson(self.CoRegPoints_table.loc[:, ['geometry', attribute2plot]]).add_to(map_osm)

        # add overlap polygon
        overlapPoly = reproject_shapelyGeometry(self.COREG_obj.overlap_poly, self.im2shift.epsg, 4326)
        gjs         = geojson.Feature(geometry=overlapPoly, properties={})
        folium.GeoJson(gjs).add_to(map_osm)


        return map_osm


    @property
    def coreg_info(self):
        if self._coreg_info:
            return self._coreg_info
        else:
            self._coreg_info = {
                'GCPList'               : self.quality_grid.GCPList,
                'reference projection'  : self.imref.prj,
                'reference geotransform': self.im2shift.gt,
                'reference grid'        : [ [self.imref.gt[0], self.imref.gt[0]+self.imref.gt[1]],
                                            [self.imref.gt[3], self.imref.gt[3]+self.imref.gt[5]] ],
                'reference extent'      : {'cols':self.imref.xgsd, 'rows':self.imref.ygsd}, # FIXME not needed anymore
                'success'               : self.success
            }
            return self.coreg_info


    def correct_shifts(self, max_GCP_count=None, cliptoextent=False):
        """Performs a local shift correction using all points from the previously calculated geometric quality grid
        that contain valid matches as GCP points.

        :param max_GCP_count: <int> maximum number of GCPs to use
        :param cliptoextent:  <bool> whether to clip the output image to its real extent
        :return:
        """

        coreg_info = self.coreg_info

        if self.quality_grid.GCPList:
            if max_GCP_count:
                coreg_info['GCPList'] = coreg_info['GCPList'][:max_GCP_count] # TODO should be a random sample

            DS = DESHIFTER(self.im2shift, coreg_info,
                           path_out     = self.path_out,
                           fmt_out      = self.fmt_out,
                           out_gsd      = (self.im2shift.xgsd,self.im2shift.ygsd),
                           align_grids  = True,
                           cliptoextent = cliptoextent,
                           #clipextent   = self.im2shift.box.boxMapYX,
                           progress     = self.progress,
                           v            = self.v,
                           q            = self.q)

            self.deshift_results = DS.correct_shifts()
            return self.deshift_results
        else:
            if not self.q:
                warnings.warn('Correction of geometric shifts failed because the input GCP list is empty!')
