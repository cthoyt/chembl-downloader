# -*- coding: utf-8 -*-

"""CLI for :mod:`chembl_downloader`."""

from typing import Optional

import click
from more_click import verbose_option

from .api import query
from .queries import ID_NAME_QUERY_EXAMPLE

__all__ = [
    "main",
]


@click.command()
@verbose_option
@click.option("--version")
def main(version: Optional[str]):
    """Test the connection."""
    df = query(ID_NAME_QUERY_EXAMPLE, columns=["chembl_id", "name"], version=version)
    click.echo(df.to_markdown())


if __name__ == "__main__":
    main()
