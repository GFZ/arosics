# -*- coding: utf-8 -*-

"""Top-level package for arosics."""

import warnings as _warnings
from pkgutil import find_loader as _find_loader

from arosics.CoReg import COREG
from arosics.CoReg_local import COREG_LOCAL
from arosics.DeShifter import DESHIFTER
from arosics.Tie_Point_Grid import Tie_Point_Grid

from .version import __version__, __versionalias__   # noqa (E402 + F401)


__author__ = """Daniel Scheffler"""
__email__ = 'daniel.scheffler@gfz-potsdam.de'
__all__ = ['COREG',
           'COREG_LOCAL',
           'DESHIFTER',
           'Tie_Point_Grid']


# check optional dependencies
if not _find_loader('pyfftw'):
    _warnings.warn('PYFFTW library is missing. However, coregistration works. But in some cases it can be much slower.')
