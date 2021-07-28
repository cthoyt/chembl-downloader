# -*- coding: utf-8 -*-

"""CLI for :mod:`chembl_downloader`."""

from typing import Optional

import click
from more_click import verbose_option

from .api import query
from .queries import ID_NAME_QUERY

__all__ = [
    "main",
]


@click.command()
@verbose_option
@click.option("--version")
def main(version: Optional[str]):
    """Test the connection."""
    click.secho("ID to Name Query\n", fg='green')
    df = query(ID_NAME_QUERY + "\nLIMIT 5", version=version)
    click.echo(df.to_markdown(index=False))


if __name__ == "__main__":
    main()
