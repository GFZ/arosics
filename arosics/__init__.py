# -*- coding: utf-8 -*-

"""Top-level package for arosics."""

import warnings
import os

if 'GDAL_DATA' not in os.environ:
    raise EnvironmentError("Please ensure that the GDAL_DATA environment variable is set and try again!")

from arosics.CoReg import COREG
from arosics.CoReg_local import COREG_LOCAL
from arosics.DeShifter import DESHIFTER
from arosics.Tie_Point_Grid import Tie_Point_Grid

__author__ = """Daniel Scheffler"""
__email__ = 'daniel.scheffler@gfz-potsdam.de'
__version__ = '0.4.26'
__versionalias__ = '2017-09-07_01'


# check optional dependencies
try:
    import pyfftw
except ImportError:
    pyfftw = None
    warnings.warn('PYFFTW library is missing. However, coregistration works. But in some cases it can be much slower.')

del warnings, pyfftw
