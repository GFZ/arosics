#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AROSICS - Automated and Robust Open-Source Image Co-Registration Software
#
# Copyright (C) 2019  Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
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

import os

# custom
import arosics

tests_path = os.path.abspath(os.path.join(arosics.__file__, "..", "..", "tests"))

# define test data pathes
test_cases = dict(
    INTER1=dict(
        ref_path=os.path.join(tests_path, 'data', 'testcase_inter1_S2A_S2A', 'ref_S2A_20160608T153121_T33UUU_sub.tif'),
        tgt_path=os.path.join(tests_path, 'data', 'testcase_inter1_S2A_S2A', 'tgt_S2A_20160529T153631_T33UUU_sub.tif'),
        kwargs_global=dict(
            path_out=os.path.join(tests_path, 'output', 'testcase_inter1_S2A_S2A/'
                                              'tgt_S2A_20160529T153631_T33UUU_sub_CR_global.bsq'),
            footprint_poly_ref='POLYGON ((340870 5862000, 354000 5862000, 354000 5830000, 331320 5830000, '
                               '340870 5862000))',
            footprint_poly_tgt='POLYGON ((341890 5866490, 356180 5866490, 356180 5834970, 335440 5834970, '
                               '335490 5845270, 341890 5866490))',
            progress=False,
            v=False),
        wp_inside=(344720, 5848485),  # inside of overlap
        wp_covering_nodata=(339611, 5856426),  # close to the image edge of the input images -> win>64px covers nodata
        wp_close_to_edge=(353810, 5840516),  # close to the image edge of the input images -> win>64px covers nodata
        wp_cloudy=(353308, 5859404),  # at a cloudy position of the target image
        wp_outside=(349533, 5818862),  # outside of overlap
        kwargs_local=dict(
            grid_res=100,
            path_out=os.path.join(tests_path, 'output', 'testcase_inter1_S2A_S2A',
                                              'tgt_S2A_20160529T153631_T33UUU_sub_CR_local.bsq'),
            progress=False)
    )
)
