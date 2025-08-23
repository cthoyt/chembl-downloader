"""Tests for the API."""

import shutil
import sqlite3
import tempfile
import unittest
from pathlib import Path
from textwrap import dedent

import pandas as pd
from pystow.constants import PYSTOW_HOME_ENVVAR
from pystow.utils import mock_envvar, open_tarfile

import chembl_downloader

TV = "26"


class TestApi(unittest.TestCase):
    """Tests for the API."""

    def setUp(self) -> None:
        """Set up the test case."""
        rows = [(1,), (2,), (3,)]
        columns = ["activity_id"]
        self.directory_obj = tempfile.TemporaryDirectory()
        self.directory = Path(self.directory_obj.name)

        self.version_directory = self.directory.joinpath(TV)
        self.version_directory.mkdir(exist_ok=True)

        data_path = self.directory.joinpath("test.db")
        self.tarfile_path = self.version_directory.joinpath(f"chembl_{TV}_sqlite.tar.gz")

        df = pd.DataFrame(rows, columns=columns)
        with sqlite3.connect(data_path) as conn:
            df.to_sql("activities", conn)

        inner_path = f"data/chembl_{TV}/chembl_{TV}_sqlite/chembl_{TV}.db"
        with (
            open_tarfile(self.tarfile_path, inner_path, operation="write") as inner_tarfile,
            data_path.open("rb") as source_file,
        ):
            shutil.copyfileobj(source_file, inner_tarfile)

    def tearDown(self) -> None:
        """Tear down the test case."""
        self.directory_obj.cleanup()

    def test_query_scalar(self) -> None:
        """Test querying."""
        with mock_envvar(PYSTOW_HOME_ENVVAR, self.directory.as_posix()):
            self.assertEqual(
                3,
                chembl_downloader.query_scalar(
                    "SELECT COUNT(activity_id) FROM activities",
                    prefix=[],
                    version=TV,
                ),
            )

    def test_cursor(self) -> None:
        """Test querying via a cursor."""
        with mock_envvar(PYSTOW_HOME_ENVVAR, self.directory.as_posix()):
            with chembl_downloader.cursor(prefix=[], version=TV) as cursor:
                cursor.execute("SELECT COUNT(activity_id) FROM activities")
                self.assertEqual(3, cursor.fetchone()[0])

    def test_latest_version(self) -> None:
        """Test getting the latest version."""
        latest_version = chembl_downloader.latest()
        self.assertIsInstance(latest_version, str)

    def test_get_target_sql(self) -> None:
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
             JOIN MOLECULE_DICTIONARY
                ON MOLECULE_DICTIONARY.molregno == ACTIVITIES.molregno
             JOIN COMPOUND_STRUCTURES
                ON MOLECULE_DICTIONARY.molregno == COMPOUND_STRUCTURES.molregno
        WHERE TARGET_DICTIONARY.chembl_id = 'CHEMBL3467'
            AND ACTIVITIES.pchembl_value IS NOT NULL
        """
        ).strip()
        self.assertEqual(expected, target)

    def test_get_target_with_max_phase(self) -> None:
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
                     JOIN MOLECULE_DICTIONARY
                        ON MOLECULE_DICTIONARY.molregno == ACTIVITIES.molregno
                     JOIN COMPOUND_STRUCTURES
                        ON MOLECULE_DICTIONARY.molregno == COMPOUND_STRUCTURES.molregno
                WHERE TARGET_DICTIONARY.chembl_id = 'CHEMBL3467'
                    AND ACTIVITIES.pchembl_value IS NOT NULL
                """
        ).strip()
        self.assertEqual(expected, target)
