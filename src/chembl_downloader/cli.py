# -*- coding: utf-8 -*-

"""CLI for :mod:`chembl_downloader`."""

import sys
from typing import Optional

import click
from more_click import verbose_option

from .api import download_extract_sqlite, query
from .queries import ACTIVITIES_QUERY, ID_NAME_QUERY

__all__ = [
    "main",
]


@click.command()
@verbose_option
@click.option("--version")
@click.option("--test", is_flag=True, help="Run a test query")
def main(version: Optional[str], test: bool):
    """Test the connection."""
    if not test:
        click.echo(download_extract_sqlite(version=version))
        sys.exit(0)

    click.secho("ID to Name Query\n", fg="green")
    df = query(ID_NAME_QUERY + "\nLIMIT 5", version=version)
    click.echo(df.to_markdown(index=False))

    click.secho("\n\nActivity Query\n", fg="green")
    df = query(ACTIVITIES_QUERY + "\nLIMIT 5", version=version)
    click.echo(df.to_markdown(index=False))


if __name__ == "__main__":
    main()
