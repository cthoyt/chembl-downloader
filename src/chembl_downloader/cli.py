# -*- coding: utf-8 -*-

"""CLI for :mod:`chembl_downloader`."""

from typing import Optional

import click
from more_click import verbose_option
from tabulate import tabulate

from .api import cursor
from .queries import ID_NAME_QUERY_EXAMPLE

__all__ = [
    "main",
]


@click.command()
@verbose_option
@click.option("--version")
def main(version: Optional[str]):
    """Test the connection."""
    with cursor(version=version) as c:
        c.execute(ID_NAME_QUERY_EXAMPLE)
        click.echo(
            tabulate(
                c.fetchall(),
                headers=["chembl_id", "name"],
            )
        )


if __name__ == "__main__":
    main()
