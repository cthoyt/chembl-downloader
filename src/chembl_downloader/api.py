# -*- coding: utf-8 -*-

"""API for :mod:`chembl_downloader`."""

import logging
import os
import sqlite3
import tarfile
from contextlib import closing, contextmanager
from pathlib import Path
from typing import Optional, Sequence, Tuple

import pystow

__all__ = [
    "download",
    "connect",
    "cursor",
]

logger = logging.getLogger(__name__)

#: The default path inside the :mod:`pystow` directory
PYSTOW_PARTS = ["chembl"]


def _download_helper(
    version: Optional[str] = None, prefix: Optional[Sequence[str]] = None
) -> Tuple[str, Path]:
    """Ensure the latest ChEMBL SQLite dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: A pair of the version and the path to the downloaded tar.gz file
    """
    if version is None:
        import bioversions

        version = bioversions.get_version("chembl")
    url = f"ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{version}/chembl_{version}_sqlite.tar.gz"
    return version, pystow.ensure(*(prefix or PYSTOW_PARTS), version, url=url)


def download(version: Optional[str] = None, prefix: Optional[Sequence[str]] = None) -> Path:
    """Get a connection as a context to the ChEMBL database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: The path to the extract ChEMBL SQLite database file
    :raises FileNotFoundError: If no database file could be found in the
        extracted directories
    """
    version, path = _download_helper(version=version, prefix=prefix)

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
    path = download(version=version, prefix=prefix)
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
