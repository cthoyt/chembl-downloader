# -*- coding: utf-8 -*-

"""API for :mod:`chembl_downloader`."""

import gzip
import logging
import os
import sqlite3
import tarfile
from contextlib import closing, contextmanager
from pathlib import Path
from typing import Optional, Sequence, TYPE_CHECKING, Tuple

import pystow

if TYPE_CHECKING:
    import pandas

__all__ = [
    "latest",
    "download_sqlite",
    "download_sdf",
    "download_extract_sqlite",
    "connect",
    "cursor",
    "query",
    "supplier",
]

logger = logging.getLogger(__name__)

#: The default path inside the :mod:`pystow` directory
PYSTOW_PARTS = ["chembl"]


def latest() -> str:
    """Get the latest version of ChEBML as a string.

    :returns: The latest version string of ChEBML
    :raises ImportError: If :mod:`bioversions` is not installed.
    """
    try:
        import bioversions
    except ImportError:
        raise ImportError("Could not import `bioversions`. Install with `pip install bioversions`.")
    else:
        return bioversions.get_version("chembl")


def _download_helper(
    suffix: str, version: Optional[str] = None, prefix: Optional[Sequence[str]] = None
) -> Tuple[str, Path]:
    """Ensure the latest ChEMBL file with the given suffix is downloaded.

    :param suffix: The suffix of the file
    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: A pair of the version and the path to the downloaded tar.gz file
    """
    if version is None:
        version = latest()
    url = f"ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{version}/chembl_{version}{suffix}"
    return version, pystow.ensure(*(prefix or PYSTOW_PARTS), version, url=url)


def download_sqlite(
    version: Optional[str] = None, prefix: Optional[Sequence[str]] = None
) -> Tuple[str, Path]:
    """Ensure the latest ChEMBL SQLite dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: A pair of the version and the local file path to the downloaded *.tar.gz file
    """
    return _download_helper(suffix="_sqlite.tar.gz", version=version, prefix=prefix)


def download_extract_sqlite(
    version: Optional[str] = None, prefix: Optional[Sequence[str]] = None
) -> Path:
    """Ensure the latest ChEMBL SQLite dump is downloaded and extracted.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: The path to the extract ChEMBL SQLite database file
    :raises FileNotFoundError: If no database file could be found in the
        extracted directories
    """
    version, path = download_sqlite(version=version, prefix=prefix)

    # Extraction will be done in the same directory as the download.
    # All ChEMBL SQLite dumps have the same internal folder structure,
    # so assume there's going to be a directory here
    directory = path.parent.joinpath("data")
    if not directory.is_dir():
        logger.info("unarchiving %s to %s", path, directory)
        with tarfile.open(path, mode="r", encoding="utf-8") as tar_file:
            tar_file.extractall(directory)
    else:
        logger.debug("did not re-unarchive %s to %s", path, directory)

    # Since the structure of the zip changes from version to version,
    # it's better to just walk through the unarchived folders recursively
    # and find the DB file
    for root, _dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".db"):
                return Path(root).joinpath(file)

    raise FileNotFoundError("could not find a .db file in the ChEMBL archive")


@contextmanager
def connect(version: Optional[str] = None, prefix: Optional[Sequence[str]] = None):
    """Ensure and connect to the database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :yields: The SQLite connection object.

    Example:
    .. code-block:: python

        import chembl_downloader

        with chembl_downloader.connect() as conn:
            with closing(conn.cursor()) as cursor:
                cursor.execute(...)
    """
    path = download_extract_sqlite(version=version, prefix=prefix)
    with closing(sqlite3.connect(path.as_posix())) as conn:
        yield conn


@contextmanager
def cursor(version: Optional[str] = None, prefix: Optional[Sequence[str]] = None):
    """Ensure, connect, and get a cursor from the database to the database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :yields: The SQLite cursor object.

    Example:
    .. code-block:: python

        import chembl_downloader

        with chembl_downloader.cursor() as cursor:
            cursor.execute(...)
    """
    with connect(version=version, prefix=prefix) as conn:
        with closing(conn.cursor()) as yv:
            yield yv


def query(
    sql: str, version: Optional[str] = None, prefix: Optional[Sequence[str]] = None, **kwargs
) -> "pandas.DataFrame":
    """Ensure the data is available, run the query, then put the results in a dataframe.

    :param sql: A SQL query string or table name
    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param kwargs: keyword arguments to pass through to :func:`pandas.read_sql`, such as
        ``index_col``.
    :return: A dataframe

    Example:
    .. code-block:: python

        import chembl_downloader
        from chembl_downloader.queries import ID_NAME_QUERY_EXAMPLE

        df = chembl_downloader.query(ID_NAME_QUERY_EXAMPLE)
    """
    import pandas as pd

    with connect(version=version, prefix=prefix) as con:
        return pd.read_sql(sql, con=con, **kwargs)


def download_sdf(
    version: Optional[str] = None, prefix: Optional[Sequence[str]] = None
) -> Tuple[str, Path]:
    """Ensure the latest ChEMBL SDF dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: A pair of the version and the local file path to the downloaded *.sdf.gz file
    """
    return _download_helper(suffix=".sdf.gz", version=version, prefix=prefix)


@contextmanager
def supplier(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    **kwargs,
):
    """Get a :class:`rdkit.Chem.ForwardSDMolSupplier` for the given version of ChEMBL.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param kwargs: keyword arguments to pass through to :class:`rdkit.Chem.ForwardSDMolSupplier`, such as
        ``sanitize`` and ``removeHs``.

    Example:
    .. code-block:: python

        from rdkit import Chem

        import chembl_downloader

        data = []
        with chembl_downloader.supplier() as suppl:
            for i, mol in enumerate(suppl):
                if mol is None or mol.GetNumAtoms() > 50:
                    continue
                fp = Chem.PatternFingerprint(mol, fpSize=1024, tautomerFingerprints=True)
                smi = Chem.MolToSmiles(mol)
                data.append((smi, fp))
    """
    from rdkit import Chem

    _, path = download_sdf(version=version, prefix=prefix)
    with gzip.open(path) as file:
        yield Chem.ForwardSDMolSupplier(file, **kwargs)
