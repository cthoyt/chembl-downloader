"""Extended functionality not in main scope of chembl-downloader."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import click

from chembl_downloader.api import VersionHint, _get_version_info, query
from chembl_downloader.queries import (
    get_assay_sql,
    get_document_molecule_sql,
    get_target_sql,
)

if TYPE_CHECKING:
    import pandas

__all__ = [
    "get_assay_smi_df",
    "get_document_smi_df",
    "get_target_smi_df",
    "write_assay_smi_file",
    "write_document_smi_file",
    "write_target_smi_file",
]


def get_target_smi_df(
    target_id: str,
    *,
    version: VersionHint | None = None,
    prefix: Sequence[str] | None = None,
    refresh: bool = False,
    standard_relation: str | None = None,
    standard_type: str | None = None,
    aggregate: Literal["mean", "gmean"] | None = "mean",
    **kwargs: Any,
) -> pandas.DataFrame:
    """Get bioactivities for compounds against the given target (from :func:`get_target_sql`).

    :param target_id: ChEMBL identifier for the target. For example, use CHEMBL1867 for
        the human A2A receptor.
    :param version: The version of ChEMBL to use. If not given, uses the latest version.
    :param aggregate: The aggregation to use (either "mean" or "gmean" for geometric
        mean). If none, do not do aggregation.
    :param refresh: If true, rebuild the cached file.
    :param standard_relation: Relation type filter, applied before aggregation. For
        example, can be "="
    :param standard_type: Assay type filter, applied before aggregation. For example,
        can be "IC50"
    :param kwargs: Remaining keyword arguments to pass through to :func:`get_target_sql`

    :returns: A dataframe

    :raises ValueError: If an unknown ``aggregate`` value is given

    Note, this caches the unfiltered, unaggregated data as a SMI file for later reuse.
    """
    import pandas as pd

    version_info = _get_version_info(version=version, prefix=prefix)

    path = version_info.module.join("targets", name=f"{target_id}.smi")
    if path.is_file() and not refresh:
        df = pd.read_csv(path)
    else:
        sql = get_target_sql(
            target_id,
            standard_relation=standard_relation,
            standard_type=standard_type,
            **kwargs,
        )
        df = query(sql=sql, version=version_info)
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
    target_id: str,
    path: Path,
    *,
    version: VersionHint | None = None,
    prefix: Sequence[str] | None = None,
    sep: str = ",",
    **kwargs: Any,
) -> None:
    """Write SMI file for the given target."""
    df = get_target_smi_df(target_id=target_id, version=version, prefix=prefix, **kwargs)
    df.to_csv(path, sep=sep, index=False, header=False)


def get_assay_smi_df(
    assay_chembl_id: str,
    *,
    version: VersionHint | None = None,
    prefix: Sequence[str] | None = None,
    refresh: bool = False,
    **kwargs: Any,
) -> pandas.DataFrame:
    """Get a dataframe for bioactivties in a given assay (from :func:`get_assay_sql`)."""
    import pandas as pd

    version_info = _get_version_info(version=version, prefix=prefix)

    path = version_info.module.join("assays", name=f"{assay_chembl_id}.smi")
    if path.is_file() and not refresh:
        df = pd.read_csv(path)
    else:
        sql = get_assay_sql(assay_chembl_id)
        df = query(sql=sql, version=version_info, **kwargs)
        df.to_csv(path, index=False)

    return df


def write_assay_smi_file(
    assay_chembl_id: str,
    path: Path,
    *,
    version: VersionHint | None = None,
    prefix: Sequence[str] | None = None,
    sep: str = ",",
    **kwargs: Any,
) -> None:
    """Write SMI file for the given assay."""
    df = get_assay_smi_df(assay_chembl_id=assay_chembl_id, version=version, prefix=prefix, **kwargs)
    df.to_csv(path, sep=sep, index=False, header=False)


def get_document_smi_df(
    document_chembl_id: str,
    *,
    version: VersionHint | None = None,
    prefix: Sequence[str] | None = None,
    refresh: bool = False,
    **kwargs: Any,
) -> pandas.DataFrame:
    """Get bioactivties in a given document (from :func:`get_document_molecule_sql`)."""
    import pandas as pd

    version_info = _get_version_info(version=version, prefix=prefix)

    path = version_info.module.join("documents", name=f"{document_chembl_id}.smi")
    if path.is_file() and not refresh:
        df = pd.read_csv(path)
    else:
        sql = get_document_molecule_sql(document_chembl_id)
        df = query(sql=sql, version=version_info, **kwargs)
        df.to_csv(path, index=False)

    return df


def write_document_smi_file(
    document_chembl_id: str,
    path: Path,
    *,
    version: VersionHint | None = None,
    prefix: Sequence[str] | None = None,
    sep: str = ",",
    **kwargs: Any,
) -> None:
    """Write SMI file for bioactivities in the given document."""
    df = get_document_smi_df(
        document_chembl_id=document_chembl_id, version=version, prefix=prefix, **kwargs
    )
    df.to_csv(path, sep=sep, index=False, header=False)


@click.command()
def _main() -> None:
    path = Path(__file__).parent.resolve().joinpath("CHEMBL3098111.smi")
    write_document_smi_file("CHEMBL3098111", path=path, version="31", refresh=True)
    click.echo(f"output to {path}")


if __name__ == "__main__":
    _main()
