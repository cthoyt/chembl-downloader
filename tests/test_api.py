"""Tests for the API."""

import unittest

import chembl_downloader


class TestApi(unittest.TestCase):
    """Tests for the API."""

    def test_latest_version(self):
        """Test getting the latest version."""
        latest_version = chembl_downloader.latest()
        self.assertIsInstance(latest_version, str)
