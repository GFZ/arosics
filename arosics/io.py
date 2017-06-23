import ctypes
import multiprocessing
import os
import time
import numpy as np
import ogr
import osr

try:
    import gdal
except ImportError:
    from osgeo import gdal
from spectral.io import envi

# internal modules
from .utilities import get_image_tileborders, convertGdalNumpyDataType
from py_tools_ds.geo.map_info import geotransform2mapinfo
from py_tools_ds.geo.projection import EPSG2WKT
from py_tools_ds.dtypes.conversion import get_dtypeStr



def wait_if_used(path_file,lockfile, timeout=100, try_kill=0):
    globs = globals()
    same_gdalRefs = [k for k,v in globs.items() if isinstance(globs[k],gdal.Dataset) and globs[k].GetDescription()==path_file]
    t0 = time.time()
    update_same_gdalRefs = lambda sRs: [sR for sR in sRs if sR in globals() and globals()[sR] is not None]
    while same_gdalRefs != [] or os.path.exists(lockfile):
        if os.path.exists(lockfile): continue
        if time.time()-t0 > timeout:
            if try_kill:
                for sR in same_gdalRefs:
                    globals()[sR] = None
                    print('had to kill %s' %sR)
            else:
                if os.path.exists(lockfile): os.remove(lockfile)
                raise TimeoutError('The file %s is permanently used by another variable.' %path_file)
        same_gdalRefs = update_same_gdalRefs(same_gdalRefs)


def write_envi(arr,outpath,gt=None,prj=None):
    if gt or prj: assert gt and prj, 'gt and prj must be provided together or left out.'
    meta = {'map info':geotransform2mapinfo(gt,prj),'coordinate system string':prj} if gt else None
    shape = (arr.shape[0],arr.shape[1],1) if len(arr.shape)==3 else arr.shape
    out    = envi.create_image(outpath,metadata=meta,shape=shape,dtype=arr.dtype,interleave='bsq',ext='.bsq',
                               force=True) # 8bit for multiple masks in one file
    out_mm = out.open_memmap(writable=True)
    out_mm[:,:,0] = arr


def wfa(p,c):
    try:
        with open(p,'a') as of: of.write(c)
    except:
        pass


shared_array = None
def init_SharedArray_in_globals(dims):
    rows, cols = dims
    global shared_array
    shared_array_base = multiprocessing.Array(ctypes.c_double, rows*cols)
    shared_array      = np.ctypeslib.as_array(shared_array_base.get_obj())
    shared_array      = shared_array.reshape(rows, cols)


def gdal_read_subset(fPath,pos,bandNr):
    (rS,rE),(cS,cE) = pos
    ds = gdal.Open(fPath)
    data = ds.GetRasterBand(bandNr).ReadAsArray(cS,rS,cE-cS+1,rE-rS+1)
    ds=None
    return data


def fill_arr(argDict, def_param=shared_array):
    pos    = argDict.get('pos')
    func   = argDict.get('func2call')
    args   = argDict.get('func_args',[])
    kwargs = argDict.get('func_kwargs',{})

    (rS,rE),(cS,cE) = pos
    shared_array[rS:rE+1,cS:cE+1] = func(*args,**kwargs)


def gdal_ReadAsArray_mp(fPath,bandNr,tilesize=1500):
    ds = gdal.Open(fPath)
    rows,cols = ds.RasterYSize, ds.RasterXSize
    ds = None

    init_SharedArray_in_globals((rows,cols))

    tilepos = get_image_tileborders([tilesize,tilesize],(rows,cols))
    fill_arr_argDicts = [{'pos':pos,'func2call':gdal_read_subset, 'func_args':(fPath,pos,bandNr)} for pos in tilepos]

    with multiprocessing.Pool() as pool: pool.map(fill_arr, fill_arr_argDicts)

    return shared_array


def write_shp(path_out, shapely_geom, prj=None,attrDict=None):
    shapely_geom = [shapely_geom] if not isinstance(shapely_geom,list) else shapely_geom
    attrDict     = [attrDict] if not isinstance(attrDict,list) else attrDict
    # print(len(shapely_geom))
    # print(len(attrDict))
    assert len(shapely_geom) == len(attrDict), "'shapely_geom' and 'attrDict' must have the same length."
    assert os.path.exists(os.path.dirname(path_out)), 'Directory %s does not exist.' % os.path.dirname(path_out)

    print('Writing %s ...' %path_out)
    if os.path.exists(path_out): os.remove(path_out)
    ds = ogr.GetDriverByName("Esri Shapefile").CreateDataSource(path_out)

    if prj is not None:
        prj = prj if not isinstance(prj, int) else EPSG2WKT(prj)
        srs = osr.SpatialReference()
        srs.ImportFromWkt(prj)
    else:
        srs = None

    geom_type = list(set([gm.type for gm in shapely_geom]))
    assert len(geom_type)==1, 'All shapely geometries must belong to the same type. Got %s.' %geom_type

    layer = ds.CreateLayer('', srs, ogr.wkbPoint)      if geom_type[0] == 'Point'      else \
            ds.CreateLayer('', srs, ogr.wkbLineString) if geom_type[0] == 'LineString' else \
            ds.CreateLayer('', srs, ogr.wkbPolygon)    if geom_type[0] == 'Polygon'    else \
            None # FIXME

    if isinstance(attrDict[0],dict):
        for attr in attrDict[0].keys():
            assert len(attr)<=10, "ogr does not support fieldnames longer than 10 digits. '%s' is too long" %attr
            DTypeStr  = get_dtypeStr(attrDict[0][attr])
            FieldType = ogr.OFTInteger if DTypeStr.startswith('int') else ogr.OFTReal     if DTypeStr.startswith('float') else \
                        ogr.OFTString  if DTypeStr.startswith('str') else ogr.OFTDateTime if DTypeStr.startswith('date') else None
            FieldDefn = ogr.FieldDefn(attr, FieldType)
            if DTypeStr.startswith('float'):
                FieldDefn.SetPrecision(6)
            layer.CreateField(FieldDefn)  # Add one attribute

    for i in range(len(shapely_geom)):
        # Create a new feature (attribute and geometry)
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometry(ogr.CreateGeometryFromWkb(shapely_geom[i].wkb)) # Make a geometry, from Shapely object

        list_attr2set = attrDict[0].keys() if isinstance(attrDict[0],dict) else []

        for attr in list_attr2set:
            val       = attrDict[i][attr]
            DTypeStr  = get_dtypeStr(val)
            val = int(val) if DTypeStr.startswith('int') else float(val) if DTypeStr.startswith('float') else \
                  str(val) if DTypeStr.startswith('str') else val
            feat.SetField(attr, val)

        layer.CreateFeature(feat)
        feat.Destroy()

    # Save and close everything
    ds = layer = None


def write_numpy_to_image(array,path_out,outFmt='GTIFF',gt=None,prj=None):
    rows,cols,bands = list(array.shape)+[1] if len(array.shape)==2 else array.shape
    gdal_dtype      = gdal.GetDataTypeByName(convertGdalNumpyDataType(array.dtype))
    outDs           = gdal.GetDriverByName(outFmt).Create(path_out,cols, rows, bands,gdal_dtype)
    for b in range(bands):
        band  = outDs.GetRasterBand(b+1)
        arr2write = array if len(array.shape)==2 else array[:,:,b]
        band.WriteArray(arr2write)
        band=None
    if gt:  outDs.SetGeoTransform(gt)
    if prj: outDs.SetProjection(prj)
    outDs=None


#def get_tempfile(ext=None,prefix=None,tgt_dir=None):
#    """Returns the path to a tempfile.mkstemp() file that can be passed to any function that expects a physical path.
#    The tempfile has to be deleted manually.
#    :param ext:     file extension (None if None)
#    :param prefix:  optional file prefix
#    :param tgt_dir: target directory (automatically set if None)
#     """
#    prefix   = 'danschef__CoReg__' if prefix is None else prefix
#    fd, path = tempfile.mkstemp(prefix=prefix,suffix=ext,dir=tgt_dir)
#    os.close(fd)
#    return path


shared_array_on_disk__memmap = None
def init_SharedArray_on_disk(out_path,dims,gt=None,prj=None):
    global shared_array_on_disk__memmap
    global shared_array_on_disk__path
    path = out_path if not os.path.splitext(out_path)[1]=='.bsq' else \
                            os.path.splitext(out_path)[0]+'.hdr'
    Meta={}
    if gt and prj:
        Meta['map info']                 = geotransform2mapinfo(gt,prj)
        Meta['coordinate system string'] = prj
    shared_array_on_disk__obj    = envi.create_image(path,metadata=Meta,shape=dims,dtype='uint16',
                                                     interleave='bsq',ext='.bsq', force=True)
    shared_array_on_disk__memmap = shared_array_on_disk__obj.open_memmap(writable=True)


def fill_arr_on_disk(argDict):
    pos       = argDict.get('pos')
    in_path   = argDict.get('in_path')
    band      = argDict.get('band')

    (rS,rE),(cS,cE) = pos
    ds   = gdal.Open(in_path)
    band = ds.GetRasterBand(band)
    data = band.ReadAsArray(cS,rS,cE-cS+1,rE-rS+1)
    shared_array_on_disk__memmap[rS:rE+1,cS:cE+1,0] = data
    ds = band = None


def convert_gdal_to_bsq__mp(in_path,out_path,band=1):
    """
    Usage:
        ref_ds,tgt_ds = gdal.Open(self.path_imref),gdal.Open(self.path_im2shift)
        ref_pathTmp, tgt_pathTmp = None,None
        if ref_ds.GetDriver().ShortName!='ENVI':
            ref_pathTmp = IO.get_tempfile(ext='.bsq')
            IO.convert_gdal_to_bsq__mp(self.path_imref,ref_pathTmp)
            self.path_imref = ref_pathTmp
        if tgt_ds.GetDriver().ShortName!='ENVI':
            tgt_pathTmp = IO.get_tempfile(ext='.bsq')
            IO.convert_gdal_to_bsq__mp(self.path_im2shift,tgt_pathTmp)
            self.path_im2shift = tgt_pathTmp
        ref_ds=tgt_ds=None

    :param in_path:
    :param out_path:
    :param band:
    :return:
    """


    ds   = gdal.Open(in_path)
    dims = (ds.RasterYSize,ds.RasterXSize)
    gt,prj = ds.GetGeoTransform(), ds.GetProjection()
    ds   = None
    init_SharedArray_on_disk(out_path,dims,gt,prj)
    positions = get_image_tileborders([512,512],dims)
    argDicts = [{'pos':pos,'in_path':in_path,'band':band} for pos in positions]
    with multiprocessing.Pool() as pool:  pool.map(fill_arr_on_disk,argDicts)
