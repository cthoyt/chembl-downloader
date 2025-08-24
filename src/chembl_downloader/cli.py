"""CLI for :mod:`chembl_downloader`."""

import sys

import click
from more_click import verbose_option
from tqdm import tqdm

from .api import (
    SummaryTuple,
    _get_version_info,
    download_extract_sqlite,
    download_readme,
    get_substructure_library,
    latest,
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
@click.option("--delete-old", is_flag=True)
def history(delete_old: bool) -> None:
    """Generate a history command."""
    import csv

    import pystow

    try:
        from tabulate import tabulate
    except ImportError:
        click.secho("Could not import `tabulate`. Please run `python -m pip install tabulate`")
        sys.exit(1)

    latest_version_info = latest(full=True)
    versions_: list[str] = versions(full=False)

    summary_path = pystow.join("chembl", name="summary.tsv")
    columns = SummaryTuple._fields
    rows = []
    for version in tqdm(versions_, desc="Summarizing ChEMBL over time", unit="versions"):
        version_info = _get_version_info(version)

        rows.append(summarize(version_info))
        if delete_old and version_info.version != latest_version_info.version:
            tqdm.write(f"[v{version_info.version}] cleaning up")

            download_readme(version=version_info, return_version=False).unlink()
            db_path = download_extract_sqlite(version=version_info, return_version=False)
            db_path.unlink()

            # if the parent directory is empty, remove it too
            version_directory = db_path.parent
            if not any(version_directory.iterdir()):
                version_directory.rmdir()

        # write on every iteration to make monitoring possible
        with summary_path.open("w") as file:
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


@main.command()
def history_draw() -> None:
    """Draw charts after running the ``history`` CLI command."""
    import matplotlib.pyplot as plt
    import pandas as pd
    import pystow
    import seaborn as sns
    from humanize import intcomma

    count_columns = ["compounds", "assays", "activities", "named_compounds"]

    summary_path = pystow.join("chembl", name="summary.tsv")
    df = pd.read_csv(summary_path, sep="\t")[::-1]

    # do this before parsing dates because it looks nicer
    chart_markdown_path = pystow.join("chembl", name="summary.md")
    df_copy = df.copy()
    for column in count_columns:
        df_copy[column] = df_copy[column].map(intcomma)
    df_copy.to_markdown(chart_markdown_path, tablefmt="github", index=False)

    df["date"] = pd.to_datetime(df["date"])

    fig, axes = plt.subplots(2, 2, figsize=(10, 5))
    for column, ax in zip(count_columns, axes.ravel(), strict=False):
        sns.lineplot(df, x="date", y=column, ax=ax)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title(column.replace("_", " ").title())
        ax.set_yscale("log")

    fig.suptitle("ChEMBL Statistics over Time")
    fig.tight_layout()

    chart_png_path = pystow.join("chembl", name="summary.png")
    chart_svg_path = pystow.join("chembl", name="summary.svg")
    fig.savefig(chart_png_path, dpi=450)
    fig.savefig(chart_svg_path)

    click.echo(f"output chart to {chart_png_path}")


if __name__ == "__main__":
    main()
