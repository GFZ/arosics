#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the local co-registration module of AROSICS."""


import unittest
import os

# custom
import arosics
from arosics import COREG_LOCAL
from geoarray import GeoArray

tests_path = os.path.abspath(os.path.join(arosics.__file__,"../../tests"))

# define test data pathes
test_cases = dict(
    INTER1=dict(
        ref_path = os.path.join(tests_path, 'data/testcase_inter1_S2A_S2A/ref_S2A_20160608T153121_T33UUU_sub.jp2'),
        tgt_path = os.path.join(tests_path, 'data/testcase_inter1_S2A_S2A/tgt_S2A_20160529T153631_T33UUU_sub.jp2'),
        kwargs = dict(
            grid_res = 100,
            progress = False)
    )
)



class COREG_LOCAL_init(unittest.TestCase):
    """Test case on object initialization of COREG_LOCAL."""

    def setUp(self):
        self.ref_path = test_cases['INTER1']['ref_path']
        self.tgt_path = test_cases['INTER1']['tgt_path']
        self.coreg_kwargs = test_cases['INTER1']['kwargs']

    def test_coreg_init_from_disk(self):
        self.CRL = COREG_LOCAL(self.ref_path, self.tgt_path, **self.coreg_kwargs)

    def test_coreg_init_from_inMem_GeoArray(self):
        # get GeoArray instances
        self.ref_gA = GeoArray(self.ref_path)
        self.tgt_gA = GeoArray(self.tgt_path)

        # assure the raster data are in-memory
        self.ref_gA.to_mem()
        self.tgt_gA.to_mem()

        # get instance of COREG_LOCAL object
        self.CRL = COREG_LOCAL(self.ref_gA, self.tgt_gA, **self.coreg_kwargs)


class CompleteWorkflow_INTER1_S2A_S2A(unittest.TestCase):
    """Test case for the complete workflow of local co-registration based on two Sentinel-2 datasets, one with
    ~25% cloud cover, the other one without any clouds. The subsets cover the S2A tiles only partly (nodata areas
    are present).
    """

    def setUp(self):
        self.ref_path = test_cases['INTER1']['ref_path']
        self.tgt_path = test_cases['INTER1']['tgt_path']
        self.coreg_kwargs = test_cases['INTER1']['kwargs']

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_calculation_of_tie_point_grid(self):
        # get instance of COREG_LOCAL object
        CRL = COREG_LOCAL(self.ref_path, self.tgt_path, **self.coreg_kwargs)

        # use the getter of the CoRegPoints_table to calculate tie point grid
        TPG = CRL.CoRegPoints_table

        # test tie point grid visualization
        #CRL.view_CoRegPoints() # only works if basemap is installed

        # test shift correction and output writer
        CRL.correct_shifts()
