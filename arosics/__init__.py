# -*- coding: utf-8 -*-

# AROSICS - Automated and Robust Open-Source Image Co-Registration Software
#
# Copyright (C) 2017-2023
# - Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
# - Helmholtz Centre Potsdam - GFZ German Research Centre for Geosciences Potsdam,
#   Germany (https://www.gfz-potsdam.de/)
#
# This software was developed within the context of the GeoMultiSens project funded
# by the German Federal Ministry of Education and Research
# (project grant code: 01 IS 14 010 A-C).
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Top-level package for arosics."""

import warnings as _warnings
from pkgutil import find_loader as _find_loader
import os as __os

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


# $PROJ_LIB was renamed to $PROJ_DATA in proj=9.1.1, which leads to issues with fiona>=1.8.20,<1.9
# https://github.com/conda-forge/pyproj-feedstock/issues/130
# -> fix it by setting PROJ_DATA
if 'GDAL_DATA' in __os.environ and 'PROJ_DATA' not in __os.environ and 'PROJ_LIB' not in __os.environ:
    __os.environ['PROJ_DATA'] = __os.path.join(__os.path.dirname(__os.environ['GDAL_DATA']), 'proj')
