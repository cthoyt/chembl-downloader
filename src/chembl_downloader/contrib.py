# -*- coding: utf-8 -*-

"""Extended functionality not in main scope of chembl-downloader."""

from pathlib import Path
from typing import TYPE_CHECKING, Optional

import click
from scipy import stats

from chembl_downloader.api import query
from chembl_downloader.queries import get_target_sql

if TYPE_CHECKING:
    import pandas

# TODO see https://github.com/PatWalters/datafiles/blob/main/CHEMBL4550.smi


def get_target_smi_df(
    target_id: str,
    *,
    version: Optional[str] = None,
    aggregate: Optional[str] = "gmean",
    **kwargs,
) -> "pandas.DataFrame":
    """Get a SMI."""
    sql = get_target_sql(target_id)
    df = query(sql=sql, version=version, **kwargs)
    if aggregate is not None:
        group_object = df.groupby(["canonical_smiles", "molecule_chembl_id"])
        if aggregate == "gmean":
            df = group_object.agg(stats.gmean)["pchembl_value"]
        elif aggregate == "mean":
            df = group_object.mean()["pchembl_value"]
        else:
            raise ValueError
    return df


def write_target_smi_file(
    target_id: str, path: Path, *, version: Optional[str] = None, sep: str = ",", **kwargs
) -> None:
    """Write SMI file for the given target."""
    df = get_target_smi_df(target_id=target_id, version=version, **kwargs)
    df.to_csv(path, sep=sep, index=False, header=False)


@click.command()
def _main():
    path = Path(__file__).parent.resolve().joinpath("CHEMBL4550.smi")
    write_target_smi_file("CHEMBL4550", path=path, version="31", aggregate=None)
    click.echo(f"output to {path}")


if __name__ == "__main__":
    _main()
