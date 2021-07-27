# -*- coding: utf-8 -*-

"""Download, open, and query ChEMBL through SQLite."""

import sqlite3
import tarfile
from pathlib import Path
from typing import Optional

import click
import pystow
from more_click import verbose_option


def ensure(version: Optional[str] = None) -> Path:
    """Ensure the latest ChEMBL SQLite dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :return: The path to the downloaded tar.gz file
    """
    if version is None:
        import bioversions
        version = bioversions.get_version('chembl')
    url = f'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{version}/chembl_{version}_sqlite.tar.gz'
    return pystow.ensure('chembl', version, url=url)


def connect(version: Optional[str] = None):
    """Get a connection as a context to the ChEMBL database."""
    path = ensure(version=version)
    with tarfile.open(path) as tar_file:
        with tar_file.extractfile('inner_path') as file:
            return sqlite3.connect(file)


@click.command()
@verbose_option
def main():
    """Test the connection."""
    con = connect()
    click.echo(str(con))


if __name__ == '__main__':
    main()
