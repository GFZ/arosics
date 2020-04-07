# -*- coding: utf-8 -*-

# AROSICS - Automated and Robust Open-Source Image Co-Registration Software
#
# Copyright (C) 2017-2020  Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
#
# This software was developed within the context of the GeoMultiSens project funded
# by the German Federal Ministry of Education and Research
# (project grant code: 01 IS 14 010 A-C).
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

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
