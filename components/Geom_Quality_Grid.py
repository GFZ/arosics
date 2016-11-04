# -*- coding: utf-8 -*-
__author__='Daniel Scheffler'

import collections
import multiprocessing
import os
import warnings
import time

# custom
try:
    import gdal
except ImportError:
    from osgeo import gdal
import numpy as np
from geopandas         import GeoDataFrame
from pykrige.ok        import OrdinaryKriging
from shapely.geometry  import Point

# internal modules
from .CoReg  import COREG
from .       import io    as IO
from py_tools_ds.ptds.geo.projection          import isProjectedOrGeographic, get_UTMzone
from py_tools_ds.ptds.io.pathgen              import get_generic_outpath
from py_tools_ds.ptds.processing.progress_mon import printProgress



global_shared_imref    = None
global_shared_im2shift = None


class Geom_Quality_Grid(object):
    """See help(Geom_Quality_Grid) for documentation!"""

    def __init__(self, COREG_obj, grid_res, outFillVal=-9999, dir_out=None, CPUs=None, progress=True, v=False, q=False):

        """Applies the algorithm to detect spatial shifts to the whole overlap area of the input images. Spatial shifts
        are calculated for each point in grid of which the parameters can be adjusted using keyword arguments. Shift
        correction performs a polynomial transformation using te calculated shifts of each point in the grid as GCPs.
        Thus 'Geom_Quality_Grid' can be used to correct for locally varying geometric distortions of the target image.

        :param COREG_obj(object):       an instance of COREG class
        :param grid_res:                grid resolution in pixels of the target image
        :param outFillVal(int):         if given the generated geometric quality grid is filled with this value in case
                                        no match could be found during co-registration (default: -9999)
        :param dir_out(str):            output directory to be used for all outputs if nothing else is given
                                        to the individual methods
        :param CPUs(int):               number of CPUs to use during calculation of geometric quality grid
                                        (default: None, which means 'all CPUs available')
        :param progress(bool):          show progress bars (default: True)
        :param v(bool):                 verbose mode (default: 0)
        :param q(bool):                 quiet mode (default: 0)
        """

        if not isinstance(COREG_obj, COREG): raise ValueError("'COREG_obj' must be an instance of COREG class.")

        self.COREG_obj  = COREG_obj
        self.grid_res   = grid_res
        self.dir_out    = dir_out
        self.outFillVal = outFillVal
        self.CPUs       = CPUs
        self.v          = v
        self.q          = q        if not v else False # overridden by v
        self.progress   = progress if not q else False # overridden by q

        self.ref        = self.COREG_obj.ref  .GeoArray
        self.shift      = self.COREG_obj.shift.GeoArray

        self.XY_points, self.XY_mapPoints = self._get_imXY__mapXY_points(self.grid_res)
        self._CoRegPoints_table           = None # set by self.CoRegPoints_table
        self._GCPList                     = None # set by self.to_GCPList()
        self.kriged                       = None # set by Raster_using_Kriging()


    @property
    def CoRegPoints_table(self):
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
        if self._GCPList:
            return self._GCPList
        else:
            self._GCPList = self.to_GCPList()
            return self._GCPList


    @GCPList.setter
    def GCPList(self, GCPList):
        self._GCPList = GCPList


    def _get_imXY__mapXY_points(self, grid_res):
        Xarr,Yarr       = np.meshgrid(np.arange(0,self.shift.shape[1],grid_res),
                                      np.arange(0,self.shift.shape[0],grid_res))

        ULmapYX, URmapYX, LRmapYX, LLmapYX = self.shift.box.boxMapYX

        mapXarr,mapYarr = np.meshgrid(np.arange(ULmapYX[1],LRmapYX[1],     self.grid_res*self.COREG_obj.shift.xgsd),
                                      np.arange(ULmapYX[0],LRmapYX[0],-abs(self.grid_res*self.COREG_obj.shift.ygsd)))

        XY_points      = np.empty((Xarr.size,2),Xarr.dtype)
        XY_points[:,0] = Xarr.flat
        XY_points[:,1] = Yarr.flat

        XY_mapPoints      = np.empty((mapXarr.size,2),mapXarr.dtype)
        XY_mapPoints[:,0] = mapXarr.flat
        XY_mapPoints[:,1] = mapYarr.flat

        return XY_points,XY_mapPoints


    @staticmethod
    def _get_spatial_shifts(coreg_kwargs):
        pointID = coreg_kwargs['pointID']
        del coreg_kwargs['pointID']

        #for im in [global_shared_imref, global_shared_im2shift]:
        #    imX, imY = mapXY2imXY(coreg_kwargs['wp'], im.gt)
        #    if im.GeoArray[int(imY), int(imX), im.band4match]==im.nodata,\
        #        return
        assert global_shared_imref    is not None
        assert global_shared_im2shift is not None
        CR = COREG(global_shared_imref, global_shared_im2shift, multiproc=False, **coreg_kwargs)
        CR.calculate_spatial_shifts()

        CR_res = [int(CR.matchWin.imDimsYX[1]), int(CR.matchWin.imDimsYX[0]),
                  CR.x_shift_px, CR.y_shift_px, CR.x_shift_map, CR.y_shift_map,
                  CR.vec_length_map, CR.vec_angle_deg]

        return [pointID]+CR_res


    def get_CoRegPoints_table(self, exclude_outliers=1):
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

        if isProjectedOrGeographic(self.COREG_obj.shift.prj)=='geographic':
            crs = dict(ellps='WGS84', datum='WGS84', proj='longlat')
        elif isProjectedOrGeographic(self.COREG_obj.shift.prj)=='projected':
            UTMzone = abs(get_UTMzone(prj=self.COREG_obj.shift.prj))
            south   = get_UTMzone(prj=self.COREG_obj.shift.prj)<0
            crs     = dict(ellps='WGS84', datum='WGS84', proj='utm', zone=UTMzone,south=south,units='m', no_defs=True)
            if not south: del crs['south']
        else:
            crs = None

        GDF        = GeoDataFrame(index=range(len(geomPoints)),crs=crs,
                                  columns=['geometry','POINT_ID','X_IM','Y_IM','X_UTM','Y_UTM'])
        GDF       ['geometry']       = geomPoints
        GDF       ['POINT_ID']       = range(len(geomPoints))
        GDF.loc[:,['X_IM','Y_IM']]   = self.XY_points
        GDF.loc[:,['X_UTM','Y_UTM']] = self.XY_mapPoints

        if exclude_outliers:
            GDF = GDF[GDF['geometry'].within(self.COREG_obj.overlap_poly)]

            #not_within_nodata = \
            #    lambda r: np.array(self.ref[r.Y_IM,r.X_IM,self.COREG_obj.ref.band4match]!=self.COREG_obj.ref.nodata and \
            #              self.shift[r.Y_IM,r.X_IM, self.COREG_obj.shift.band4match] != self.COREG_obj.shift.nodata)[0,0]

            #GDF['not_within_nodata'] = GDF.apply(lambda GDF_row: not_within_nodata(GDF_row),axis=1)
            #GDF = GDF[GDF.not_within_nodata]

        # declare global variables needed for self._get_spatial_shifts()
        global global_shared_imref,global_shared_im2shift
        assert self.ref  .footprint_poly # this also checks for mask_nodata and nodata value
        assert self.shift.footprint_poly
        if not self.ref  .is_inmem: self.ref.cache_array_subset(self.ref  [self.COREG_obj.ref  .band4match])
        if not self.shift.is_inmem: self.ref.cache_array_subset(self.shift[self.COREG_obj.shift.band4match])
        global_shared_imref    = self.ref
        global_shared_im2shift = self.shift

        # get all variations of kwargs for coregistration
        get_coreg_kwargs = lambda pID, wp: {
            'pointID'         : pID,
            'wp'              : wp,
            'ws'              : self.COREG_obj.win_size_XY,
            'data_corners_im0': self.COREG_obj.ref.corner_coord,
            'data_corners_im1': self.COREG_obj.shift.corner_coord,
            'r_b4match'       : self.COREG_obj.ref.band4match+1,   # band4match is internally saved as index, starting from 0
            's_b4match'       : self.COREG_obj.shift.band4match+1, # band4match is internally saved as index, starting from 0
            'max_iter'        : self.COREG_obj.max_iter,
            'max_shift'       : self.COREG_obj.max_shift,
            'nodata'          : (self.COREG_obj.ref.nodata, self.COREG_obj.shift.nodata),
            'binary_ws'       : self.COREG_obj.bin_ws,
            'v'               : False, # otherwise this would lead to massive console output
            'q'               : True,  # otherwise this would lead to massive console output
            'ignore_errors'   : True
        }
        list_coreg_kwargs = (get_coreg_kwargs(i, self.XY_mapPoints[i]) for i in GDF.index) # generator

        # run co-registration for whole grid
        if self.CPUs is None or self.CPUs>1:
            if not self.q:
                print("Calculating geometric quality grid (%s points) in mode 'multiprocessing'..." %len(GDF))

            t0 = time.time()
            with multiprocessing.Pool(self.CPUs) as pool:
                if self.q or not self.progress:
                    results = pool.map(self._get_spatial_shifts, list_coreg_kwargs)
                else:
                    results = pool.map_async(self._get_spatial_shifts, list_coreg_kwargs, chunksize=1)

                    while True:
                        time.sleep(.1)
                        numberDone = len(GDF)-results._number_left # this does not really represent the remaining tasks but the remaining chunks -> thus chunksize=1
                        printProgress(percent=numberDone/len(GDF)*100, barLength=50, prefix='\tprogress:',
                                      suffix='[%s/%s] Complete %.2f sek'%(numberDone,len(GDF), time.time()-t0))
                        if results.ready():
                            results = results.get() # FIXME in some cases the code hangs here ==> x
                            break
        else:
            if not self.q:
                print("Calculating geometric quality grid (%s points) in mode 'singleprocessing'..." %len(GDF))
            results = np.empty((len(geomPoints),9))
            for i,coreg_kwargs in enumerate(list_coreg_kwargs):
                if not self.q:
                    printProgress(percent=(i+1)/len(GDF)*100, prefix='\tprogress:',
                                  suffix='[%s/%s] Complete'%((i+1),len(GDF)), decimals=1, barLength=50)
                results[i,:] = self._get_spatial_shifts(coreg_kwargs)
            # FIXME in some cases the code hangs here ==> x

         # merge results with GDF
        records = GeoDataFrame(np.array(results, np.object), columns=['POINT_ID', 'X_WIN_SIZE', 'Y_WIN_SIZE',
                                                                      'X_SHIFT_PX','Y_SHIFT_PX', 'X_SHIFT_M',
                                                                      'Y_SHIFT_M', 'ABS_SHIFT', 'ANGLE'])
        GDF = GDF.merge(records, on='POINT_ID', how="inner")
        GDF = GDF.fillna(int(self.outFillVal))

        self.CoRegPoints_table = GDF
        return GDF


    def dump_CoRegPoints_table(self, path_out=None):
        path_out = path_out if path_out else get_generic_outpath(dir_out=self.dir_out,
            fName_out="CoRegPoints_table_grid%s_ws(%s_%s)__T_%s__R_%s.pkl" % (self.grid_res, self.COREG_obj.win_size_XY[0],
                        self.COREG_obj.win_size_XY[1], self.shift.basename, self.ref.basename))
        if not self.q:
            print('Writing %s ...' % path_out)
        self.CoRegPoints_table.to_pickle(path_out)


    def to_GCPList(self):
        # get copy of quality grid without no data
        GDF = self.CoRegPoints_table.loc[self.CoRegPoints_table.X_SHIFT_M != self.outFillVal, :].copy()

        if getattr(GDF,'empty'): # GDF.empty returns AttributeError
            return []
        else:
            # calculate GCPs
            GDF['X_UTM_new'] = GDF.X_UTM + GDF.X_SHIFT_M
            GDF['Y_UTM_new'] = GDF.Y_UTM + GDF.Y_SHIFT_M
            GDF['GCP']       = GDF.apply(lambda GDF_row:    gdal.GCP(GDF_row.X_UTM_new, GDF_row.Y_UTM_new, 0,
                                                                     GDF_row.X_IM, GDF_row.Y_IM), axis=1)
            self.GCPList = GDF.GCP.tolist()

            if not self.q:
                print('Found %s valid GCPs.' %len(self.GCPList))

            return self.GCPList


    def test_if_singleprocessing_equals_multiprocessing_result(self):
        self.CPUs = None
        dataframe = self.get_CoRegPoints_table()
        mp_out    = np.empty_like(dataframe.values)
        mp_out[:] = dataframe.values
        self.CPUs = 1
        dataframe = self.get_CoRegPoints_table()
        sp_out    = np.empty_like(dataframe.values)
        sp_out[:] = dataframe.values

        return np.array_equal(sp_out,mp_out)


    def _get_line_by_PID(self, PID):
        return self.CoRegPoints_table.loc[PID, :]


    def _get_lines_by_PIDs(self, PIDs):
        assert isinstance(PIDs,list)
        lines = np.zeros((len(PIDs),self.CoRegPoints_table.shape[1]))
        for i,PID in enumerate(PIDs):
            lines[i,:] = self.CoRegPoints_table[self.CoRegPoints_table['POINT_ID'] == PID]
        return lines


    def to_PointShapefile(self, skip_nodata=1, skip_nodata_col ='ABS_SHIFT', path_out=None):
        GDF            = self.CoRegPoints_table
        GDF2pass       = GDF if not skip_nodata else GDF[GDF[skip_nodata_col]!=self.outFillVal]

        path_out = path_out if path_out else get_generic_outpath(dir_out=os.path.join(self.dir_out, 'CoRegPoints'),
            fName_out="CoRegPoints_grid%s_ws(%s_%s)__T_%s__R_%s.shp" % (self.grid_res, self.COREG_obj.win_size_XY[0],
                        self.COREG_obj.win_size_XY[1], self.shift.basename, self.ref.basename))
        if not self.q:
            print('Writing %s ...' %path_out)
        GDF2pass.to_file(path_out)


    def _to_PointShapefile(self, skip_nodata=1, skip_nodata_col ='ABS_SHIFT'):
        warnings.warn(DeprecationWarning("'_quality_grid_to_PointShapefile' is deprecated." # TODO delete if other method validated
                                         " 'quality_grid_to_PointShapefile' is much faster."))
        GDF            = self.CoRegPoints_table
        GDF2pass       = GDF if not skip_nodata else GDF[GDF[skip_nodata_col]!=self.outFillVal]
        shapely_points = GDF2pass['geometry'].values.tolist()
        attr_dicts     = [collections.OrderedDict(zip(GDF2pass.columns,GDF2pass.loc[i].values)) for i in GDF2pass.index]


        fName_out = "CoRegPoints_grid%s_ws%s.shp" %(self.grid_res, self.COREG_obj.win_size_XY)
        path_out = os.path.join(self.dir_out, fName_out)
        IO.write_shp(path_out, shapely_points, prj=self.COREG_obj.shift.prj, attrDict=attr_dicts)


    def to_Raster_using_KrigingOLD(self, attrName, skip_nodata=1, skip_nodata_col='ABS_SHIFT', outGridRes=None,
                                   path_out=None, tilepos=None):
        GDF             = self.CoRegPoints_table
        GDF2pass        = GDF if not skip_nodata else GDF[GDF[skip_nodata_col]!=self.outFillVal]

        # subset if tilepos is given
        rows,cols = tilepos if tilepos else self.shift.shape
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

        path_out = path_out if path_out else get_generic_outpath(dir_out=os.path.join(self.dir_out, 'CoRegPoints'),
            fName_out="Kriging__%s__grid%s_ws(%s_%s).tif"  % (attrName, self.grid_res, self.COREG_obj.win_size_XY[0],
                        self.COREG_obj.win_size_XY[1]))
        print('Writing %s ...' %path_out)
        # add a half pixel grid points are centered on the output pixels
        xmin,ymin,xmax,ymax = xmin-grid_res/2,ymin-grid_res/2,xmax+grid_res/2,ymax+grid_res/2
        IO.write_numpy_to_image(zvalues, path_out, gt=(xmin, grid_res, 0, ymax, 0, -grid_res), prj=self.COREG_obj.shift.prj)

        return zvalues


    def Raster_using_Kriging(self, attrName, skip_nodata=1, skip_nodata_col='ABS_SHIFT', outGridRes=None,
                             fName_out=None, tilepos=None, tilesize=500, mp=None):

        mp = False if self.CPUs==1 else True
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
        GDF             = self.CoRegPoints_table
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

        if self.CPUs is None or self.CPUs>1:
            fName_out = fName_out if fName_out else \
                "Kriging__%s__grid%s_ws%s_%s.tif" %(attrName,self.grid_res, self.COREG_obj.win_size_XY,tilepos)
        else:
            fName_out = fName_out if fName_out else \
                "Kriging__%s__grid%s_ws%s.tif" %(attrName,self.grid_res, self.COREG_obj.win_size_XY)
        path_out  = get_generic_outpath(dir_out=self.dir_out, fName_out=fName_out)
        print('Writing %s ...' %path_out)
        # add a half pixel grid points are centered on the output pixels
        xmin,ymin,xmax,ymax = xmin-grid_res/2,ymin-grid_res/2,xmax+grid_res/2,ymax+grid_res/2
        IO.write_numpy_to_image(zvalues, path_out, gt=(xmin, grid_res, 0, ymax, 0, -grid_res), prj=self.COREG_obj.shift.prj)

        return zvalues


    def _Kriging_mp(self, args_kwargs_dict):
        args   = args_kwargs_dict.get('args',[])
        kwargs = args_kwargs_dict.get('kwargs',[])

        return self._Kriging_sp(*args, **kwargs)
