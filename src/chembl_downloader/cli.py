"""CLI for :mod:`chembl_downloader`."""

import sys

import click
from more_click import verbose_option
from tqdm import tqdm

from .api import (
    SummaryTuple,
    VersionHint,
    download_extract_sqlite,
    get_substructure_library,
    query,
    query_scalar,
    summarize,
    versions,
)
from .queries import (
    ACTIVITIES_QUERY,
    COUNT_ACTIVITIES_SQL,
    COUNT_ASSAYS_SQL,
    COUNT_COMPOUNDS_SQL,
    ID_NAME_QUERY,
)

__all__ = [
    "main",
]

version_option = click.option("--version", help="The ChEMBL version to use. Defaults to latest.")


@click.group()
def main() -> None:
    """Test the connection."""


@main.command()
@version_option
@verbose_option  # type:ignore
def download(version: str | None) -> None:
    """Download the data."""
    click.echo(download_extract_sqlite(version=version))


@main.command()
@version_option
@verbose_option  # type:ignore
def test(version: str | None) -> None:
    """Run test queries."""
    click.secho("Number of Activities", fg="green")
    count = query_scalar(COUNT_ACTIVITIES_SQL, version=version)
    click.echo(f"{count:,}")

    click.secho("\nNumber of Compounds", fg="green")
    count = query_scalar(COUNT_COMPOUNDS_SQL, version=version)
    click.echo(f"{count:,}")

    click.secho("\nNumber of Assays", fg="green")
    count = query_scalar(COUNT_ASSAYS_SQL, version=version)
    click.echo(f"{count:,}\n")

    click.secho("ID to Name Query\n", fg="green")
    df = query(ID_NAME_QUERY + "\nLIMIT 5", version=version)
    click.echo(df.to_markdown(index=False))

    click.secho("\nActivity Query\n", fg="green")
    df = query(ACTIVITIES_QUERY + "\nLIMIT 5", version=version)
    click.echo(df.to_markdown(index=False))


@main.command()
@version_option
@verbose_option  # type:ignore
def substructure(version: str | None) -> None:
    """Build a substructure library."""
    get_substructure_library(version=version)


@main.command()
def history() -> None:
    """Generate a history command."""
    import csv

    import pystow

    try:
        from tabulate import tabulate
    except ImportError:
        click.secho("Could not import `tabulate`. Please run `python -m pip install tabulate`")
        sys.exit(1)

    versions_: list[VersionHint] = list(versions())
    rows = [summarize(version) for version in tqdm(versions_)]
    columns = SummaryTuple._fields

    with pystow.join("chembl", name="summary.tsv").open("w") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)

    click.echo(
        tabulate(
            rows,
            tablefmt="github",
            headers=columns,
            intfmt=",",
        )
    )


if __name__ == "__main__":
    main()
