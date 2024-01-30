"""Tests for querying the target."""

import unittest

import chembl_downloader


class TestQueries(unittest.TestCase):
    """Tests for the target."""

    def test_get_target_sql(self):
        """Test getting the target sql."""
        target = chembl_downloader.queries.get_target_sql()
        self.assertIsInstance(target, str)

# if __name__ == '__main__':
#     unittest.main()