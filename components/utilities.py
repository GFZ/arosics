import numpy as np
import datetime


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

def get_dtypeStr(val):
    is_numpy = 'numpy' in str(type(val))
    DType = str(np.dtype(val)) if is_numpy else 'int' if isinstance(val,int) else 'float' if isinstance(val,float) else \
                'str' if isinstance(val,str) else 'complex' if isinstance(val,complex) else \
                'date' if isinstance(val,datetime.datetime) else None
    assert DType is not None, 'data type not understood'
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