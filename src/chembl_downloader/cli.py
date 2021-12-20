# -*- coding: utf-8 -*-

"""CLI for :mod:`chembl_downloader`."""

from typing import Optional

import click
from more_click import verbose_option

from .api import download_extract_sqlite, get_substructure_library, query
from .queries import ACTIVITIES_QUERY, ID_NAME_QUERY

__all__ = [
    "main",
]

version_option = click.option("--version", help="The ChEMBL version to use.")


@click.group()
def main():
    """Test the connection."""


@main.command()
@version_option
@verbose_option
def download(version: Optional[str]):
    """Download the data."""
    click.echo(download_extract_sqlite(version=version))


@main.command()
@version_option
@verbose_option
def test(version: Optional[str]):
    """Run test queries."""
    click.secho("ID to Name Query\n", fg="green")
    df = query(ID_NAME_QUERY + "\nLIMIT 5", version=version)
    click.echo(df.to_markdown(index=False))

    click.secho("\n\nActivity Query\n", fg="green")
    df = query(ACTIVITIES_QUERY + "\nLIMIT 5", version=version)
    click.echo(df.to_markdown(index=False))


@main.command()
@version_option
@verbose_option
def substructure(version: Optional[str]):
    """Build a substructure library."""
    get_substructure_library(version=version)


if __name__ == "__main__":
    main()
