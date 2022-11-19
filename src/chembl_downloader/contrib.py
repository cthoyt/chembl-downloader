# -*- coding: utf-8 -*-

"""Extended functionality not in main scope of chembl-downloader."""

from pathlib import Path
from typing import TYPE_CHECKING, Optional

import click
import pystow

from chembl_downloader.api import latest, query
from chembl_downloader.queries import (
    get_assay_sql,
    get_document_molecule_sql,
    get_target_sql,
)

if TYPE_CHECKING:
    import pandas

__all__ = [
    "get_target_smi_df",
    "write_target_smi_file",
    "get_assay_smi_df",
    "write_assay_smi_file",
    "get_document_smi_df",
    "write_document_smi_file",
]


def get_target_smi_df(
    target_id: str,
    *,
    version: Optional[str] = None,
    aggregate: Optional[str] = "mean",
    refresh: bool = False,
    **kwargs,
) -> "pandas.DataFrame":
    """Get a SMI."""
    import pandas as pd

    if version is None:
        version = latest()

    path = pystow.join("chembl", version, "targets", name=f"{target_id}.smi")
    if path.is_file() and not refresh:
        df = pd.read_csv(path)
    else:
        sql = get_target_sql(target_id)
        df = query(sql=sql, version=version, **kwargs)
        df.to_csv(path, index=False)

    if aggregate is not None:
        group_object = df.groupby(["canonical_smiles", "molecule_chembl_id"])
        if aggregate == "gmean":
            from scipy import stats

            df = group_object.agg(stats.gmean)["pchembl_value"]
        elif aggregate == "mean":
            df = group_object.mean()["pchembl_value"]
        else:
            raise ValueError(f"unknown aggregate: {aggregate}")
    return df


def write_target_smi_file(
    target_id: str, path: Path, *, version: Optional[str] = None, sep: str = ",", **kwargs
) -> None:
    """Write SMI file for the given target."""
    df = get_target_smi_df(target_id=target_id, version=version, **kwargs)
    df.to_csv(path, sep=sep, index=False, header=False)


def get_assay_smi_df(
    assay_chembl_id: str,
    *,
    version: Optional[str] = None,
    refresh: bool = False,
    **kwargs,
) -> "pandas.DataFrame":
    """Get a SMI."""
    import pandas as pd

    if version is None:
        version = latest()

    path = pystow.join("chembl", version, "assays", name=f"{assay_chembl_id}.smi")
    if path.is_file() and not refresh:
        df = pd.read_csv(path)
    else:
        sql = get_assay_sql(assay_chembl_id)
        df = query(sql=sql, version=version, **kwargs)
        df.to_csv(path, index=False)

    return df


def write_assay_smi_file(
    assay_chembl_id: str, path: Path, *, version: Optional[str] = None, sep: str = ",", **kwargs
) -> None:
    """Write SMI file for the given assay."""
    df = get_assay_smi_df(assay_chembl_id=assay_chembl_id, version=version, **kwargs)
    df.to_csv(path, sep=sep, index=False, header=False)


def get_document_smi_df(
    document_chembl_id: str,
    *,
    version: Optional[str] = None,
    refresh: bool = False,
    **kwargs,
) -> "pandas.DataFrame":
    """Get a SMI."""
    import pandas as pd

    if version is None:
        version = latest()

    path = pystow.join("chembl", version, "documents", name=f"{document_chembl_id}.smi")
    if path.is_file() and not refresh:
        df = pd.read_csv(path)
    else:
        sql = get_document_molecule_sql(document_chembl_id)
        df = query(sql=sql, version=version, **kwargs)
        df.to_csv(path, index=False)

    return df


def write_document_smi_file(
    document_chembl_id: str, path: Path, *, version: Optional[str] = None, sep: str = ",", **kwargs
) -> None:
    """Write SMI file for the given document."""
    df = get_document_smi_df(document_chembl_id=document_chembl_id, version=version, **kwargs)
    df.to_csv(path, sep=sep, index=False, header=False)


@click.command()
def _main():
    path = Path(__file__).parent.resolve().joinpath("CHEMBL3098111.smi")
    write_document_smi_file("CHEMBL3098111", path=path, version="31", refresh=True)
    click.echo(f"output to {path}")


if __name__ == "__main__":
    _main()
