# -*- coding: utf-8 -*-

"""Extended functionality not in main scope of chembl-downloader."""

from pathlib import Path
from typing import TYPE_CHECKING, Optional

import click
import pystow
from typing_extensions import Literal

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
    refresh: bool = False,
    standard_relation: Optional[str] = None,
    standard_type: Optional[str] = None,
    aggregate: Optional[Literal["mean", "gmean"]] = "mean",
    **kwargs,
) -> "pandas.DataFrame":
    """Geta dataframe for activities of compounds against the given target.

    :param target_id: ChEMBL identifier for the target.
        For example, use CHEMBL1867 for the human A2A receptor.
    :param version: The version of ChEMBL to use. If not given, uses the latest version.
    :param aggregate:
        The aggregation to use (either "mean" or "gmean" for geometric mean).
        If none, do not do aggregation.
    :param refresh:
        If true, rebuild the cached file.
    :param standard_relation:
        Relation type filter, applied before aggregation. For example, can be "="
    :param standard_type:
        Assay type filter, applied before aggregation. For example, can be "IC50"
    :param kwargs:
        Remaining keyword arguments to pass through to :func:`get_target_sql`
    :return:
        A dataframe
    :raises ValueError:
        If an unknown ``aggregate`` value is given

    Note, this caches the unfiltered, unaggregated data as a SMI file for later reuse.
    """
    import pandas as pd

    if version is None:
        version = latest()

    path = pystow.join("chembl", version, "targets", name=f"{target_id}.smi")
    if path.is_file() and not refresh:
        df = pd.read_csv(path)
    else:
        sql = get_target_sql(
            target_id,
            standard_relation=standard_relation,
            standard_type=standard_type,
            **kwargs,
        )
        df = query(sql=sql, version=version)
        df.to_csv(path, index=False)

    if aggregate is not None:
        group_object = df[["canonical_smiles", "molecule_chembl_id", "pchembl_value"]].groupby(
            ["canonical_smiles", "molecule_chembl_id"]
        )
        if aggregate == "gmean":
            from scipy import stats

            df = group_object.agg(stats.gmean)["pchembl_value"].reset_index()
        elif aggregate == "mean":
            df = group_object.mean(numeric_only=True)["pchembl_value"].reset_index()
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
