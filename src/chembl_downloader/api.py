# -*- coding: utf-8 -*-

"""API for :mod:`chembl_downloader`."""

import sqlite3
import tarfile
from contextlib import closing, contextmanager
from pathlib import Path
from typing import Optional, Sequence

import pystow

__all__ = [
    "download",
    "connect",
    "cursor",
]

PYSTOW_PARTS = "pyobo", "raw", "chembl.compound"


def _download_helper(version: Optional[str] = None, prefix: Optional[Sequence[str]] = None) -> Path:
    """Ensure the latest ChEMBL SQLite dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: The path to the downloaded tar.gz file
    """
    if version is None:
        import bioversions

        version = bioversions.get_version("chembl")
    url = f"ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{version}/chembl_{version}_sqlite.tar.gz"
    return pystow.ensure(*(prefix or PYSTOW_PARTS), version, url=url)


def download(version: Optional[str] = None, prefix: Optional[Sequence[str]] = None):
    """Get a connection as a context to the ChEMBL database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: The path to the extract ChEMBL SQLite database file
    """
    path = _download_helper(version=version, prefix=prefix)
    directory = path.parent.joinpath(f"chembl_{version}")
    rv = directory.joinpath(f"chembl_{version}_sqlite", f"chembl_{version}.db")
    if path.parent.joinpath(f"chembl_{version}").is_dir():
        return rv
    with tarfile.open(path, mode="r", encoding="utf-8") as tar_file:
        tar_file.extractall(path.parent)
    if not rv.is_file():
        raise FileNotFoundError(rv.as_posix())
    return rv


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
