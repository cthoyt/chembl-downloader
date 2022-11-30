# -*- coding: utf-8 -*-

"""CLI for :mod:`chembl_downloader`."""

import sys
from typing import Optional

import click
from more_click import verbose_option
from tqdm import tqdm

from .api import (
    download_extract_sqlite,
    get_date,
    get_substructure_library,
    query,
    versions,
)
from .queries import ACTIVITIES_QUERY, ID_NAME_QUERY

__all__ = [
    "main",
]

version_option = click.option("--version", help="The ChEMBL version to use. Defaults to latest.")


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


@main.command()
def history():
    """Generate a history command."""
    try:
        from tabulate import tabulate
    except ImportError:
        click.secho("Could not import `tabulate`. Please run `python -m pip install tabulate`")
        sys.exit(1)
    rows = [
        (version, get_date(version=version), _count_compounds(version=version))
        for version in tqdm(versions())
    ]
    click.echo(
        tabulate(
            rows,
            tablefmt="github",
            headers=["ChEMBL Version", "Release Date", "Total Named Compounds *from SQLite*"],
        )
    )


def _count_compounds(version: str) -> str:
    """Test downloader for specific ChEMBL version."""
    from .queries import COUNT_QUERY_SQL

    try:
        total_compounds = query(COUNT_QUERY_SQL, version=version)["count"][0]
    except Exception:
        return "-"
    return f"{total_compounds:,}"


if __name__ == "__main__":
    main()
