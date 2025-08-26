"""CLI for :mod:`chembl_downloader`."""

from functools import partial
from pathlib import Path

import click
from more_click import verbose_option
from tqdm import tqdm
from tqdm.contrib.concurrent import thread_map

from .api import (
    SummaryTuple,
    VersionHint,
    VersionInfo,
    _get_version_info,
    download_extract_sqlite,
    download_readme,
    get_substructure_library,
    latest,
    query,
    summarize,
    versions,
)
from .queries import ACTIVITIES_QUERY, ID_NAME_QUERY

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


@main.command(name="summarize")
@version_option
@verbose_option  # type:ignore
def print_summary(version: str | None) -> None:
    """Run test queries and print results to the console."""
    summary = summarize(version=version)

    click.secho("Number of Activities", fg="green")
    click.echo(f"{summary.activities:,}")

    click.secho("\nNumber of Compounds", fg="green")
    click.echo(f"{summary.compounds:,}")

    click.secho("\nNumber of Assays", fg="green")
    click.echo(f"{summary.assays:,}\n")

    click.echo(str(summary))

    click.secho("ID to Name Query", fg="green")
    try:
        df = query(ID_NAME_QUERY + "\nLIMIT 5", version=version)
    except OSError:
        click.secho("failed", fg="red")
    else:
        click.echo("\n" + df.to_markdown(index=False))

    click.secho("\nActivity Query", fg="green")
    try:
        df = query(ACTIVITIES_QUERY + "\nLIMIT 5", version=version)
    except OSError:
        click.secho("failed", fg="red")
    else:
        click.echo("\n" + df.to_markdown(index=False))


@main.command()
@version_option
@verbose_option  # type:ignore
def substructure(version: str | None) -> None:
    """Build a substructure library."""
    get_substructure_library(version=version)


TEST_VERSIONS = ["1", "19", "25", "35"]


@main.command()
@click.option(
    "--delete-old",
    is_flag=True,
    help="Delete old versions as you go, useful if you have limited hard disk "
    "space but also more time intensive on re-run.",
)
@click.option("--test", is_flag=True, help=f"Run on test versions: {', '.join(TEST_VERSIONS)}")
@click.option(
    "--max-workers",
    type=int,
    default=1,
    help="Use threading to download/process multiple at the same time. In practice, "
    "this doesn't work because the ChEMBL server is limited in how fast it can serve the data.",
)
@click.pass_context
def history(ctx: click.Context, delete_old: bool, test: bool, max_workers: int) -> None:
    """Generate a historical analysis of ChEMBL."""
    import csv

    import pystow

    latest_version_info = latest(full=True)
    if test:
        versions_: list[str] = TEST_VERSIONS
    else:
        versions_ = versions(full=False)

    summary_path = pystow.join("chembl", name="summary.tsv")
    columns = SummaryTuple._fields

    f = partial(_help_summarize, latest_version_info=latest_version_info, delete_old=delete_old)
    rows = []
    tqdm_kwargs = {
        "desc": "Summarizing ChEMBL over time",
        "unit": "versions",
    }
    if max_workers > 1:
        row_it = thread_map(f, versions_, max_workers=2, **tqdm_kwargs)
    else:
        row_it = map(f, tqdm(versions_, **tqdm_kwargs))

    for row in row_it:
        rows.append(row)

        # write on every iteration to make monitoring possible
        with summary_path.open("w") as file:
            writer = csv.writer(file, delimiter="\t")
            writer.writerow(columns)
            writer.writerows(rows)

    ctx.invoke(history_draw)


def _help_summarize(
    version: VersionHint, latest_version_info: VersionInfo, delete_old: bool
) -> SummaryTuple:
    version_info = _get_version_info(version)

    rv = summarize(version_info)
    if (
        delete_old
        and version_info.version != latest_version_info.version
        and version_info.version not in TEST_VERSIONS
    ):
        tqdm.write(f"[v{version_info.version}] cleaning up")

        download_readme(version=version_info, return_version=False).unlink()
        db_path = download_extract_sqlite(version=version_info, return_version=False)
        db_path.unlink()

        # if the parent directory is empty, remove it too
        version_directory = db_path.parent
        if not any(version_directory.iterdir()):
            version_directory.rmdir()
    return rv


@main.command()
def history_draw() -> None:
    """Draw charts after running the ``history`` CLI command."""
    import matplotlib.pyplot as plt
    import pandas as pd
    import pystow
    import seaborn as sns
    from humanize import intcomma

    here = Path(__file__).parent.resolve()
    data_dir = here.parent.parent.joinpath("docs", "_data")

    count_columns = SummaryTuple._fields[2:]

    summary_path = pystow.join("chembl", name="summary.tsv")
    df = pd.read_csv(summary_path, sep="\t")
    if data_dir:
        df.to_csv(data_dir.joinpath("summary.tsv"), sep="\t", index=False)

    # do this before parsing dates because it looks nicer
    chart_markdown_path = pystow.join("chembl", name="summary.md")
    df_copy = df.copy()
    for column in count_columns:
        df_copy[column] = df_copy[column].map(intcomma)
    df_copy.to_markdown(chart_markdown_path, tablefmt="github", index=False)

    df["date"] = pd.to_datetime(df["date"])

    n_rows = len(count_columns) // 2
    figsize = (8.5, 1.8 * n_rows + 0.5)

    fig, axes = plt.subplots(n_rows, 2, figsize=figsize, sharex=True)
    for column, ax in zip(count_columns, axes.ravel(), strict=False):
        sns.lineplot(df[df[column] > 0], x="date", y=column, ax=ax)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title(column.replace("_", " ").title())

    fig.suptitle("ChEMBL Statistics over Time")
    fig.tight_layout()

    chart_stub = pystow.join("chembl", name="summary")
    fig.savefig(chart_stub.with_suffix(".png"), dpi=450)
    fig.savefig(chart_stub.with_suffix(".svg"))
    if data_dir:
        fig.savefig(data_dir.joinpath("summary.svg"))

    click.echo(f"output chart to {chart_stub}.svg")

    df_diff = df.set_index("version")[::-1].diff().reset_index()
    df_diff = df_diff[df_diff["date"].notna()]
    for col in count_columns:
        df_diff[col] = df_diff[col].astype(int)

    fig, axes = plt.subplots(n_rows, 2, figsize=figsize)
    for column, ax in zip(count_columns, axes.ravel(), strict=False):
        sns.lineplot(
            df_diff,
            x=df_diff["version"].map(lambda s: float(s.replace("_", "."))),
            y=column,
            ax=ax,
        )
        ax.set_xlabel("Version")
        ax.set_ylabel("")
        ax.set_yscale("symlog")
        ax.set_title(column.replace("_", " ").title())
        ax.axhline(0.0)

    chart_diff_stub = pystow.join("chembl", name="summary-diff")
    fig.suptitle("ChEMBL Statistics (Discrete Derivative)")
    fig.tight_layout()
    fig.savefig(chart_diff_stub.with_suffix(".png"), dpi=450)
    fig.savefig(chart_diff_stub.with_suffix(".svg"))
    if data_dir:
        df_diff[::-1].to_csv(data_dir.joinpath("summary-diff.tsv"), sep="\t", index=False)
        fig.savefig(data_dir.joinpath("summary-diff.svg"))


if __name__ == "__main__":
    main()
