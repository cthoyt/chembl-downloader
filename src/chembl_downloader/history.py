"""Generate a version history table to include in the README."""

import click
from tabulate import tabulate
from tqdm import tqdm

from chembl_downloader.api import download_readme, get_date, versions


@click.command()
def _main():
    rows = [(version, get_date(version=version)) for version in tqdm(versions())]
    click.echo(tabulate(rows, tablefmt="github", headers=["ChEMBL Version", "Release Date"]))


if __name__ == "__main__":
    _main()
