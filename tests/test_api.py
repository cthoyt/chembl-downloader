"""Tests for the API."""

import unittest
from textwrap import dedent

import chembl_downloader


class TestApi(unittest.TestCase):
    """Tests for the API."""

    def test_latest_version(self):
        """Test getting the latest version."""
        latest_version = chembl_downloader.latest()
        self.assertIsInstance(latest_version, str)

    def test_get_target_sql(self):
        """Test getting the target sql."""
        target = chembl_downloader.queries.get_target_sql("CHEMBL3467")
        expected = dedent(
            """\
        SELECT
            ASSAYS.chembl_id              AS assay_chembl_id,
            TARGET_DICTIONARY.target_type,
            TARGET_DICTIONARY.tax_id,
            COMPOUND_STRUCTURES.canonical_smiles,
            MOLECULE_DICTIONARY.chembl_id AS molecule_chembl_id,
            ACTIVITIES.standard_type,
            ACTIVITIES.pchembl_value
        FROM TARGET_DICTIONARY
             JOIN ASSAYS ON TARGET_DICTIONARY.tid == ASSAYS.tid
             JOIN ACTIVITIES ON ASSAYS.assay_id == ACTIVITIES.assay_id
             JOIN MOLECULE_DICTIONARY ON MOLECULE_DICTIONARY.molregno == ACTIVITIES.molregno
             JOIN COMPOUND_STRUCTURES ON MOLECULE_DICTIONARY.molregno == COMPOUND_STRUCTURES.molregno
        WHERE TARGET_DICTIONARY.chembl_id = 'CHEMBL3467'
            AND ACTIVITIES.pchembl_value IS NOT NULL
        """
        ).strip()
        self.assertEqual(expected, target)

    def test_get_target_with_max_phase(self):
        """Test getting target sql with a max phase."""
        target = chembl_downloader.queries.get_target_sql("CHEMBL3467", max_phase=True)
        expected = dedent(
            """\
                SELECT
                    ASSAYS.chembl_id              AS assay_chembl_id,
                    TARGET_DICTIONARY.target_type,
                    TARGET_DICTIONARY.tax_id,
                    COMPOUND_STRUCTURES.canonical_smiles,
                    MOLECULE_DICTIONARY.chembl_id AS molecule_chembl_id,
                    MOLECULE_DICTIONARY.max_phase,
                    ACTIVITIES.standard_type,
                    ACTIVITIES.pchembl_value
                FROM TARGET_DICTIONARY
                     JOIN ASSAYS ON TARGET_DICTIONARY.tid == ASSAYS.tid
                     JOIN ACTIVITIES ON ASSAYS.assay_id == ACTIVITIES.assay_id
                     JOIN MOLECULE_DICTIONARY ON MOLECULE_DICTIONARY.molregno == ACTIVITIES.molregno
                     JOIN COMPOUND_STRUCTURES ON MOLECULE_DICTIONARY.molregno == COMPOUND_STRUCTURES.molregno
                WHERE TARGET_DICTIONARY.chembl_id = 'CHEMBL3467'
                    AND ACTIVITIES.pchembl_value IS NOT NULL
                """
        ).strip()
        self.assertEqual(expected, target)
