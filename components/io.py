import ctypes
import multiprocessing
import os
import time
import numpy as np

try:
    import gdal
except ImportError:
    from osgeo import gdal

from .utilities import get_image_tileborders




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


def write_envi(arr,outpath):
    from spectral.io import envi as envi
    shape = (arr.shape[0],arr.shape[1],1) if len(arr.shape)==3 else arr.shape
    out    = envi.create_image(outpath,shape=shape,dtype=arr.dtype,interleave='bsq',ext='.bsq', force=True) # 8bit for muliple masks in one file
    out_mm = out.open_memmap(writable=True)
    out_mm[:,:,0] = arr


def wfa(p,c):
    try:
        with open(p,'a') as of: of.write(c)
    except: pass


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