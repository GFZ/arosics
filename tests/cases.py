#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

# custom
import arosics

tests_path = os.path.abspath(os.path.join(arosics.__file__,"../../tests"))

# define test data pathes
test_cases = dict(
    INTER1=dict(
        ref_path = os.path.join(tests_path, 'data/testcase_inter1_S2A_S2A/ref_S2A_20160608T153121_T33UUU_sub.jp2'),
        tgt_path = os.path.join(tests_path, 'data/testcase_inter1_S2A_S2A/tgt_S2A_20160529T153631_T33UUU_sub.jp2'),
        kwargs_global = dict(
            path_out = os.path.join(tests_path, 'output/testcase_inter1_S2A_S2A/'
                                                'tgt_S2A_20160529T153631_T33UUU_sub_CR_global.bsq'),
            progress = False,
            v = 0),
        wp_inside  = (344720, 5848485),# (344932, 5842974), # inside of overlap
        wp_outside = (349533, 5818862), # outside of overlap
        kwargs_local = dict(
            grid_res = 100,
            path_out=os.path.join(tests_path, 'output/testcase_inter1_S2A_S2A/'
                                              'tgt_S2A_20160529T153631_T33UUU_sub_CR_local.bsq'),
            progress = False)
    )
)
