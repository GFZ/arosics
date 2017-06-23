#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `arosics` package."""


import unittest

from arosics import COREG
from arosics import COREG_LOCAL


class TestArosics(unittest.TestCase):
    """Tests for `arosics` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        try:
            from spectral.io import envi
            print('Import worked during test!')
        except ImportError:
            print('Could not import spectral!!')

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_000_something(self):
        """Test something."""
