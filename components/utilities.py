import numpy as np
import datetime
import collections
import bisect


def find_nearest(array, value, roundAlg='auto', extrapolate=False, exclude_val=False):
    """finds the value of an array nearest to a another single value
    NOTE: In case of extrapolation an EQUALLY INCREMENTED array (like a coordinate grid) is assumed!

    :param array:
    :param value:
    :param roundAlg:
    :param extrapolate: extrapolate the given array if the given value is outside the array
    :param exclude_val:
    """
    assert roundAlg in ['auto', 'on', 'off']
    assert isinstance(array, list) or (isinstance(array, np.ndarray) and len(array.shape) == 1)
    array = sorted(list(array))

    if exclude_val and value in array:
        array.remove(value)

    if extrapolate:
        increment = array[1] - array[0]
        if value > max(array):  # expand array until value
            array = np.arange(min(array), value + increment, increment)
        if value < min(array):  # negatively expand array until value
            array = (np.arange(-max(array), value + increment, increment) * -1)[::-1]
    elif value < min(array) or value > max(array):
        raise ValueError('Value %s is outside of the given array.' % value)

    if roundAlg == 'auto':
        diffs = np.abs(np.array(array) - value)
        minDiff = diffs.min()
        minIdx = diffs.argmin()
        isMiddleVal = collections.Counter(diffs)[minDiff] > 1
        out = array[minIdx] if not isMiddleVal else array[bisect.bisect_left(array, value)]
    elif roundAlg == 'off':
        out = array[bisect.bisect_left(array, value) - 1]
    else: # roundAlg == 'on'
        out = array[bisect.bisect_left(array, value)]

    return out

def get_dtypeStr(val):
    is_numpy = 'numpy' in str(type(val))
    DType = str(np.dtype(val)) if is_numpy else 'int' if isinstance(val,int) else 'float' if isinstance(val,float) else \
                'str' if isinstance(val,str) else 'complex' if isinstance(val,complex) else \
                'date' if isinstance(val,datetime.datetime) else None
    assert DType, 'data type not understood'
    return DType

def get_image_tileborders(target_tileSize,shape_fullArr):
    rows,cols     = shape_fullArr[:2]
    row_bounds=[0]
    while row_bounds[-1]+target_tileSize[0] < rows:
        row_bounds.append(row_bounds[-1] + target_tileSize[0]-1)
        row_bounds.append(row_bounds[-2] + target_tileSize[0])
    else:
        row_bounds.append(rows-1)

    col_bounds=[0]
    while col_bounds[-1]+target_tileSize[1] < cols:
        col_bounds.append(col_bounds[-1] + target_tileSize[1]-1)
        col_bounds.append(col_bounds[-2] + target_tileSize[1])
    else:
        col_bounds.append(cols-1)
    return [[tuple([row_bounds[r], row_bounds[r+1]]), tuple([col_bounds[c], col_bounds[c+1]])] \
            for r in range(0, len(row_bounds), 2) for c in range(0, len(col_bounds), 2)]

def cornerPoints_to_listOfXYPairs(corYX,out_shape,out_resXY,shrinkSize=None):
    """
    :param corYX: list of XY pairs in the order UL,UR,LR,LL
    """
    Xarr         = np.zeros(out_shape,np.float64)
    Xarr[None,:] = np.arange(corYX[0][1],corYX[1][1],out_resXY[0])
    Xarr         = Xarr
    Yarr         = np.zeros(list(reversed(out_shape)),np.float64)
    out_resY     = out_resXY[1] if corYX[0][0]<corYX[2][0] else -out_resXY[1]
    Yarr[None,:] = np.arange(corYX[0][0],corYX[2][0],out_resY)
    Yarr         = Yarr.T

    Xarr = Xarr[shrinkSize:-shrinkSize,shrinkSize:-shrinkSize] if shrinkSize else Xarr
    Yarr = Yarr[shrinkSize:-shrinkSize,shrinkSize:-shrinkSize] if shrinkSize else Yarr

    XYarr = np.empty((Xarr.size,2),np.float64)
    XYarr[:,0] = Xarr.flat
    XYarr[:,1] = Yarr.flat
    return XYarr

def get_coord_grid(ULxy,LRxy,out_resXY):
    X_vec = np.arange(ULxy[0],LRxy[0],out_resXY[0])
    Y_vec = np.arange(ULxy[1],LRxy[1],out_resXY[1])
    return np.meshgrid(X_vec,Y_vec)

def convertGdalNumpyDataType(dType):
    """convertGdalNumpyDataType
    :param dType: GDALdataType string or numpy dataType
    :return: corresponding dataType
    """
    # dictionary to translate GDAL data types (strings) in corresponding numpy data types
    dTypeDic = {"Byte": np.uint8, "UInt16": np.uint16, "Int16": np.int16, "UInt32": np.uint32, "Int32": np.int32,
                "Float32": np.float32, "Float64": np.float64, "GDT_UInt32": np.uint32}
    outdType = None

    if dType in dTypeDic:
        outdType = dTypeDic[dType]
    elif dType in dTypeDic.values():
        for i in dTypeDic.items():
            if dType == i[1]:
                outdType = i[0]
    elif dType in [np.int8, np.int64, np.int]:
        outdType = "Int32"
        print(">>>  Warning: %s is converted to GDAL_Type 'Int_32'\n" % dType)
    elif dType in [np.bool, np.bool_]:
        outdType = "Byte"
        print(">>>  Warning: %s is converted to GDAL_Type 'Byte'\n" % dType)
    elif dType in [np.float]:
        outdType = "Float32"
        print(">>>  Warning: %s is converted to GDAL_Type 'Float32'\n" % dType)
    elif dType in [np.float16]:
        outdType = "Float32"
        print(">>>  Warning: %s is converted to GDAL_Type 'Float32'\n" % dType)
    else:
        raise Exception('GEOP.convertGdalNumpyDataType: Unexpected input data type %s.' %dType)
    return outdType
