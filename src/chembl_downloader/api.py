# -*- coding: utf-8 -*-

"""API for :mod:`chembl_downloader`."""

import sqlite3
import tarfile
from contextlib import closing, contextmanager
from pathlib import Path
from typing import Optional

import pystow
import requests_ftp

__all__ = [
    'ensure',
    'ensure_extract',
    'get_connection',
    'get_cursor',
]

requests_ftp.monkeypatch_session()
PYSTOW_PARTS = 'pyobo', 'raw', 'chembl.compound'


def ensure(version: Optional[str] = None) -> Path:
    """Ensure the latest ChEMBL SQLite dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :return: The path to the downloaded tar.gz file
    """
    if version is None:
        import bioversions
        version = bioversions.get_version('chembl')
    url = f'ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{version}/chembl_{version}_sqlite.tar.gz'
    return pystow.ensure(*PYSTOW_PARTS, version, url=url)


def ensure_extract(version: Optional[str] = None):
    """Get a connection as a context to the ChEMBL database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :return: The path to the extract ChEMBL SQLite database file
    """
    path = ensure(version=version)
    directory = path.parent.joinpath(f'chembl_{version}')
    rv = directory.joinpath(f"chembl_{version}_sqlite", f"chembl_{version}.db")
    if path.parent.joinpath(f'chembl_{version}').is_dir():
        return rv
    with tarfile.open(path, mode="r", encoding="utf-8") as tar_file:
        tar_file.extractall(path.parent)
    if not rv.is_file():
        raise FileNotFoundError(rv.as_posix())
    return rv


@contextmanager
def get_connection(version: Optional[str] = None):
    """Ensure and connect to the database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :yields: The SQLite connection object.

    Example:

    .. code-block:: python

        from chembl_downloader import get_connection

        with get_connection() as conn:
            with closing(conn.cursor()) as cursor:
                cursor.execute(...)
    """
    path = ensure_extract(version=version)
    with closing(sqlite3.connect(path.as_posix())) as conn:
        yield conn


@contextmanager
def get_cursor(version: Optional[str] = None):
    """Ensure, connect, and get a cursor from the database to the database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :yields: The SQLite cursor object.

    Example:

    .. code-block:: python

        from chembl_downloader import cursor

        with get_cursor() as cursor:
            cursor.execute(...)
    """
    with get_connection(version=version) as conn:
        with closing(conn.cursor()) as cursor:
            yield cursor
