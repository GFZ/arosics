# -*- coding: utf-8 -*-

"""Top-level package for arosics."""

import warnings

from arosics.CoReg import COREG
from arosics.CoReg_local import COREG_LOCAL
from arosics.DeShifter import DESHIFTER
from arosics.Tie_Point_Grid import Tie_Point_Grid

__author__ = """Daniel Scheffler"""
__email__ = 'daniel.scheffler@gfz-potsdam.de'
__version__ = '0.4.23'
__versionalias__ = '2017-09-05_03'


# check optional dependencies
try:
    import pyfftw
except ImportError:
    pyfftw = None
    warnings.warn('PYFFTW library is missing. However, coregistration works. But in some cases it can be much slower.')

del warnings, pyfftw
