#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AROSICS - Automated and Robust Open-Source Image Co-Registration Software
#
# Copyright (C) 2017-2024
# - Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz.de)
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
#   https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for the module arosics.Tie_Point_Grid."""

import unittest
import tempfile
import os
from importlib.util import find_spec
import shutil
import warnings
import struct

# custom
import pytest
import numpy as np

from .cases import test_cases
from arosics import COREG_LOCAL, Tie_Point_Grid


class Test_Tie_Point_Grid(unittest.TestCase):

    @classmethod
    def setUp(cls):
        CRL = COREG_LOCAL(test_cases['INTER1']['ref_path'], test_cases['INTER1']['tgt_path'],
                          **test_cases['INTER1']['kwargs_local'])
        cls.TPG = Tie_Point_Grid(CRL.COREG_obj, CRL.grid_res,
                                 max_points=100,  # limit to 100 to reduce computational load
                                 outFillVal=CRL.outFillVal,
                                 resamp_alg_calc=CRL.rspAlg_calc,
                                 tieP_filter_level=CRL.tieP_filter_level,
                                 tieP_random_state=CRL.tieP_random_state,
                                 outlDetect_settings=dict(
                                     min_reliability=CRL.min_reliability,
                                     rs_max_outlier=CRL.rs_max_outlier,
                                     rs_tolerance=CRL.rs_tolerance),
                                 dir_out=CRL.projectDir,
                                 CPUs=CRL.CPUs,
                                 progress=CRL.progress,
                                 v=CRL.v,
                                 q=CRL.q)

    def tearDown(self):
        if os.path.isdir(self.TPG.dir_out):
            shutil.rmtree(self.TPG.dir_out)

    def test_mean_shifts(self):
        assert isinstance(self.TPG.mean_x_shift_px, float)
        assert isinstance(self.TPG.mean_y_shift_px, float)
        assert isinstance(self.TPG.mean_x_shift_map, float)
        assert isinstance(self.TPG.mean_y_shift_map, float)

    def test_get_CoRegPoints_table(self):
        self.TPG.get_CoRegPoints_table()

    def test_calc_rmse(self):
        self.TPG.calc_rmse(include_outliers=False)
        self.TPG.calc_rmse(include_outliers=True)

    def test_calc_overall_ssim(self):
        self.TPG.calc_overall_ssim(include_outliers=False, after_correction=True)
        self.TPG.calc_overall_ssim(include_outliers=True, after_correction=False)

    def test_calc_overall_stats(self):
        stats_noOL = self.TPG.calc_overall_stats(include_outliers=False)
        stats_OL = self.TPG.calc_overall_stats(include_outliers=True)

        assert stats_noOL
        assert stats_OL
        assert isinstance(stats_noOL, dict)
        assert isinstance(stats_OL, dict)
        assert stats_noOL != stats_OL

    def test_plot_shift_distribution(self):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', category=UserWarning, message='Matplotlib is currently using agg, '
                                                        'which is a non-GUI backend, so cannot show the figure.')
            self.TPG.plot_shift_distribution()

    def test_dump_CoRegPoints_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, 'CoRegPoints_table.pkl')
            self.TPG.dump_CoRegPoints_table(outpath)
            assert os.path.isfile(outpath)

    def test_to_GCPList(self):
        self.TPG.to_GCPList()

    def test_to_PointShapefile(self):
        with (pytest.warns(UserWarning, match='.*Column names longer than 10 characters will be truncated.*'),
              warnings.catch_warnings()
              ):
            warnings.filterwarnings("ignore", message=".*recognized as too large to be valid.*")

            tbl = self.TPG.CoRegPoints_table

            n_all_points = len(tbl)
            n_nodata = sum(tbl['ABS_SHIFT'] == self.TPG.outFillVal)
            n_outliers = sum(tbl['OUTLIER'].__eq__(True))

            def _get_n_records(filepath: str):
                with open(filepath, 'rb') as inF:
                    header = inF.read(32)
                    return struct.unpack('<I', header[4:8])[0]

            with tempfile.TemporaryDirectory() as tmpdir:
                outpath = os.path.join(tmpdir, 'test_out_shapefile_excl_nodata.shp')
                self.TPG.to_PointShapefile(outpath, skip_nodata=True)
                assert os.path.isfile(outpath)
                assert _get_n_records(f'{outpath[:-4]}.dbf') == n_all_points - n_nodata

            with tempfile.TemporaryDirectory() as tmpdir:
                outpath = os.path.join(tmpdir, 'test_out_shapefile_incl_nodata.shp')
                self.TPG.to_PointShapefile(outpath, skip_nodata=False)
                assert os.path.isfile(outpath)
                assert _get_n_records(f'{outpath[:-4]}.dbf') == n_all_points

            with tempfile.TemporaryDirectory() as tmpdir:
                outpath = os.path.join(tmpdir, 'test_out_shapefile_excl_outliers.shp')
                self.TPG.to_PointShapefile(outpath, skip_nodata=False, skip_outliers=True)
                assert os.path.isfile(outpath)
                assert _get_n_records(f'{outpath[:-4]}.dbf') == n_all_points - n_outliers

            with tempfile.TemporaryDirectory() as tmpdir:
                outpath = os.path.join(tmpdir, 'test_out_shapefile_excl_nodata_outliers.shp')
                self.TPG.to_PointShapefile(outpath, skip_nodata=True, skip_outliers=True)
                assert os.path.isfile(outpath)
                assert _get_n_records(f'{outpath[:-4]}.dbf') == n_all_points - n_nodata - n_outliers

    def test_to_vectorfield(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, 'test_vectorfield.bsq')
            self.TPG.to_vectorfield(outpath, fmt='ENVI', mode='md')
            assert os.path.isfile(outpath)
            self.TPG.to_vectorfield(outpath, fmt='ENVI', mode='uv')
            assert os.path.isfile(outpath)

    def test_interpolate_to_raster_rbf(self):
        arr_interp = self.TPG.to_interpolated_raster('ABS_SHIFT', 'RBF', plot_result=True)

        assert isinstance(arr_interp, np.ndarray)

    def test_interpolate_to_raster_gpr(self):
        if find_spec('sklearn'):
            arr_interp = self.TPG.to_interpolated_raster('ABS_SHIFT', 'GPR', plot_result=True)

            assert isinstance(arr_interp, np.ndarray)

    def test_interpolate_to_raster_kriging(self):
        if find_spec('pykrige.ok'):
            arr_interp = self.TPG.to_interpolated_raster('ABS_SHIFT', 'Kriging', plot_result=True)

            assert isinstance(arr_interp, np.ndarray)

    def test_random_state(self):
        self.TPG.tieP_random_state = None
        point_ids = [self.TPG.get_CoRegPoints_table()['POINT_ID'] for _ in range(2)]
        assert not np.array_equal(point_ids[0], point_ids[1]), \
            "Samples should not be identical when random state is None"

        self.TPG.tieP_random_state = 0
        point_ids = [self.TPG.get_CoRegPoints_table()['POINT_ID'] for _ in range(2)]
        assert np.array_equal(point_ids[0], point_ids[1]), \
            "Samples should be identical when random state is fixed"


if __name__ == '__main__':
    pytest.main()
