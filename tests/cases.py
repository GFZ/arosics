#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AROSICS - Automated and Robust Open-Source Image Co-Registration Software
#
# Copyright (C) 2017-2021
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
