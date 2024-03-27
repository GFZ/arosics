#!/usr/bin/env python
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

"""Tests for the local co-registration module of AROSICS."""

import unittest
import shutil
import os
import warnings

# custom
from .cases import test_cases
from arosics import COREG_LOCAL
from geoarray import GeoArray

class COREG_LOCAL_init(unittest.TestCase):
    """Test case on object initialization of COREG_LOCAL."""

    def setUp(self):
        self.ref_path = test_cases["INTER1"]["ref_path"]
        self.tgt_path = test_cases["INTER1"]["tgt_path"]
        self.coreg_kwargs = test_cases["INTER1"]["kwargs_local"]

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
        self.ref_path = test_cases["INTER1"]["ref_path"]
        self.tgt_path = test_cases["INTER1"]["tgt_path"]
        self.coreg_kwargs = test_cases["INTER1"]["kwargs_local"]

    def tearDown(self):
        """Delete output."""
        dir_out = os.path.dirname(self.coreg_kwargs["path_out"])
        if os.path.isdir(dir_out):
            shutil.rmtree(dir_out)

    def test_calculation_of_tie_point_grid(self):
        # get instance of COREG_LOCAL object
        CRL = COREG_LOCAL(self.ref_path, self.tgt_path, **self.coreg_kwargs)

        # calculate tie point grid
        CRL.calculate_spatial_shifts()

        # test tie point grid visualization
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=UserWarning,
                message="Matplotlib is currently using agg, " "which is a non-GUI backend, so cannot show the figure.",
            )
            CRL.view_CoRegPoints(hide_filtered=True)
            CRL.view_CoRegPoints(hide_filtered=False)
            CRL.view_CoRegPoints(shapes2plot="vectors")
            CRL.view_CoRegPoints_folium()

        # test shift correction and output writer
        CRL.correct_shifts()

        self.assertTrue(
            os.path.exists(self.coreg_kwargs["path_out"]), "Output of local co-registration has not been written."
        )

    def test_calculation_of_tie_point_grid_float_coords(self):
        # NOTE: This does not test against unequaly sized output of get_image_windows_to_match().

        # overwrite gt and prj
        ref = GeoArray(self.ref_path)
        ref.to_mem()
        ref.filePath = None
        tgt = GeoArray(self.tgt_path)
        tgt.to_mem()
        tgt.filePath = None

        ref.gt = [330000.19999996503, 10.00000001, 0.0, 5862000.7999997628, 0.0, -10.00000001]
        # ref.gt = [330000.1, 10.1, 0.0, 5862000.1, 0.0, -10.1]
        tgt.gt = [335440.19999996503, 10.00000001, 0.0, 5866490.7999997628, 0.0, -10.00000001]
        # tgt.gt = [330000.1, 10.1, 0.0, 5862000.1, 0.0, -10.1]

        # get instance of COREG_LOCAL object
        CRL = COREG_LOCAL(ref, tgt, **self.coreg_kwargs)
        CRL.calculate_spatial_shifts()
        # CRL.view_CoRegPoints()

    def test_calculation_of_tie_point_grid_noepsg(self):
        """Test local coregistration with a proj. other than LonLat and UTM and a WKT which has no EPSG code (FORCE)."""
        wkt_noepsg = """
            PROJCRS["BU MEaSUREs Lambert Azimuthal Equal Area - SA - V01",
                BASEGEOGCRS["WGS 84",
                    DATUM["World Geodetic System 1984",
                        ELLIPSOID["WGS 84",6378137,298.257223563,
                            LENGTHUNIT["metre",1]]],
                    PRIMEM["Greenwich",0,
                        ANGLEUNIT["degree",0.0174532925199433]],
                    ID["EPSG",4326]],
                CONVERSION["unnamed",
                    METHOD["Lambert Azimuthal Equal Area",
                        ID["EPSG",9820]],
                    PARAMETER["Latitude of natural origin",-15,
                        ANGLEUNIT["degree",0.0174532925199433],
                        ID["EPSG",8801]],
                    PARAMETER["Longitude of natural origin",-60,
                        ANGLEUNIT["degree",0.0174532925199433],
                        ID["EPSG",8802]],
                    PARAMETER["False easting",0,
                        LENGTHUNIT["metre",1],
                        ID["EPSG",8806]],
                    PARAMETER["False northing",0,
                        LENGTHUNIT["metre",1],
                        ID["EPSG",8807]]],
                CS[Cartesian,2],
                    AXIS["easting",east,
                        ORDER[1],
                        LENGTHUNIT["metre",1]],
                    AXIS["northing",north,
                        ORDER[2],
                        LENGTHUNIT["metre",1]]]
            """
        wkt_noepsg = " ".join(wkt_noepsg.split())

        # overwrite prj
        ref = GeoArray(self.ref_path)
        ref.to_mem()
        ref.filePath = None
        tgt = GeoArray(self.tgt_path)
        tgt.to_mem()
        tgt.filePath = None

        ref.prj = wkt_noepsg
        tgt.prj = wkt_noepsg

        # get instance of COREG_LOCAL object
        CRL = COREG_LOCAL(ref, tgt, **self.coreg_kwargs)
        CRL.calculate_spatial_shifts()
        # CRL.view_CoRegPoints()

    def test_calculation_of_tie_point_grid_with_metaRotation(self):
        """Test with default parameters - should compute X/Y shifts properly and write the de-shifted target image."""

        # overwrite gt and prj
        ref = GeoArray(self.ref_path)
        ref.to_mem()
        ref.filePath = None
        ref.gt = [330000, 10, 0.00001, 5862000, 0.00001, -10]
        tgt = GeoArray(self.tgt_path)
        tgt.to_mem()
        tgt.filePath = None
        # tgt.gt = [335440, 5.8932, 0.0, 5866490, 0.0, -10.1]
        tgt.gt = [335440, 10, 0.00001, 5866490, 0.00001, -10]

        # get instance of COREG_LOCAL object
        CRL = COREG_LOCAL(ref, tgt, **self.coreg_kwargs)
        CRL.calculate_spatial_shifts()
        # CRL.view_CoRegPoints()

        self.assertTrue(CRL.success)


if __name__ == "__main__":
    import pytest

    pytest.main()
