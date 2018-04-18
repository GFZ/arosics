#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the local co-registration module of AROSICS."""

import unittest
import shutil
import os
from importlib import util

# custom
from .cases import test_cases
from arosics import COREG_LOCAL
from geoarray import GeoArray


class COREG_LOCAL_init(unittest.TestCase):
    """Test case on object initialization of COREG_LOCAL."""

    def setUp(self):
        self.ref_path = test_cases['INTER1']['ref_path']
        self.tgt_path = test_cases['INTER1']['tgt_path']
        self.coreg_kwargs = test_cases['INTER1']['kwargs_local']

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
        self.coreg_kwargs = test_cases['INTER1']['kwargs_local']

    def tearDown(self):
        """Delete output."""
        dir_out = os.path.dirname(self.coreg_kwargs['path_out'])
        if os.path.isdir(dir_out):
            shutil.rmtree(dir_out)

    def test_calculation_of_tie_point_grid(self):
        # get instance of COREG_LOCAL object
        CRL = COREG_LOCAL(self.ref_path, self.tgt_path, **self.coreg_kwargs)

        # use the getter of the CoRegPoints_table to calculate tie point grid
        # noinspection PyStatementEffect
        CRL.CoRegPoints_table

        # test tie point grid visualization
        if util.find_spec('mpl_toolkits.basemap'):  # only works if basemap is installed
            CRL.view_CoRegPoints(hide_filtered=True)
            CRL.view_CoRegPoints(hide_filtered=False)
            CRL.view_CoRegPoints(shapes2plot='vectors')

        if util.find_spec('folium') and util.find_spec('geojson'):
            CRL.view_CoRegPoints_folium()

        # test shift correction and output writer
        CRL.correct_shifts()

        self.assertTrue(os.path.exists(self.coreg_kwargs['path_out']),
                        'Output of local co-registration has not been written.')

# if __name__ == '__main__':
#     unittest.main(argv=['first-arg-is-ignored'],exit=False, verbosity=2)
#
#      suite = unittest.TestLoader().loadTestsFromTestCase(eval("CompleteWorkflow_INTER1_S2A_S2A"))
#     alltests = unittest.TestSuite(suite)
#
#      # Part 2: Saving the results of each testsuite and the query for the job.status in individual variables.
#      testResult = unittest.TextTestRunner(verbosity=2).run(alltests)
