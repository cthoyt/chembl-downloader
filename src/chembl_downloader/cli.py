# -*- coding: utf-8 -*-

"""CLI for :mod:`chembl_downloader`."""

import click
from more_click import verbose_option

from .api import connect


@click.command()
@verbose_option
def main():
    """Test the connection."""
    con = connect()
    click.echo(str(con))


if __name__ == '__main__':
    main()
