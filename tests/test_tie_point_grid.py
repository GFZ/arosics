#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the module arosics.Tie_Point_Grid."""

import unittest
import tempfile
import os
from pkgutil import find_loader
import shutil

# custom
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
        self.assertIsInstance(self.TPG.mean_x_shift_px, float)
        self.assertIsInstance(self.TPG.mean_y_shift_px, float)
        self.assertIsInstance(self.TPG.mean_x_shift_map, float)
        self.assertIsInstance(self.TPG.mean_y_shift_map, float)

    def test_get_CoRegPoints_table(self):
        self.TPG.get_CoRegPoints_table()

    def test_calc_rmse(self):
        self.TPG.calc_rmse(include_outliers=False)
        self.TPG.calc_rmse(include_outliers=True)

    def test_calc_overall_ssim(self):
        self.TPG.calc_overall_ssim(include_outliers=False, after_correction=True)
        self.TPG.calc_overall_ssim(include_outliers=True, after_correction=False)

    def test_plot_shift_distribution(self):
        self.TPG.plot_shift_distribution()

    def test_dump_CoRegPoints_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, 'CoRegPoints_table.pkl')
            self.TPG.dump_CoRegPoints_table(outpath)
            self.assertTrue(os.path.isfile(outpath))

    def test_to_GCPList(self):
        self.TPG.to_GCPList()

    def test_to_PointShapefile(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, 'test_out_shapefile.shp')
            self.TPG.to_PointShapefile(outpath)
            self.assertTrue(os.path.isfile(outpath))

        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, 'test_out_shapefile_incl_nodata.shp')
            self.TPG.to_PointShapefile(outpath, skip_nodata=False)
            self.assertTrue(os.path.isfile(outpath))

    def test_to_vectorfield(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, 'test_vectorfield.bsq')
            self.TPG.to_vectorfield(outpath, fmt='ENVI', mode='md')
            self.assertTrue(os.path.isfile(outpath))
            self.TPG.to_vectorfield(outpath, fmt='ENVI', mode='uv')
            self.assertTrue(os.path.isfile(outpath))

    def test_to_Raster_using_Kriging(self):
        if find_loader('pykrige.ok'):
            with tempfile.TemporaryDirectory() as tmpdir:
                outpath = os.path.join(tmpdir, 'X_SHIFT_M__interpolated.bsq')
                self.TPG.to_Raster_using_Kriging(attrName='X_SHIFT_M', fName_out=outpath)
                self.assertTrue(os.path.isfile(outpath))
