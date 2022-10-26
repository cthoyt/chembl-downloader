"""Generate a version history table to include in the README."""

import click
from tabulate import tabulate
from tqdm import tqdm

from chembl_downloader.api import download_readme, versions


@click.command()
def _main():
    rows = []
    for version in tqdm(versions()):
        path = download_readme(version=version)
        try:
            date_p = (
                next(line for line in path.read_text().splitlines() if line.startswith("* Date:"))
                .removeprefix("* Date:")
                .lstrip()
            )
            day, month, year = date_p.split("/")
            date = f"{year}-{month}-{day}"
        except StopIteration:
            date = ""  # happens on 22.1 and 24.1
        rows.append((version, date))
    click.echo(tabulate(rows, tablefmt="github", headers=["ChEMBL Version", "Release Date"]))


if __name__ == "__main__":
    _main()
