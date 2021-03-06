# -*- coding: utf-8 -*-

"""API for :mod:`chembl_downloader`."""

import gzip
import logging
import os
import pickle
import sqlite3
import tarfile
from contextlib import closing, contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Sequence, Tuple, Union, cast
from xml.etree import ElementTree

import pystow
from tqdm import tqdm

if TYPE_CHECKING:
    import pandas

__all__ = [
    "latest",
    # Database
    "download_sqlite",
    "download_extract_sqlite",
    "connect",
    "cursor",
    "query",
    # SDF
    "download_sdf",
    "supplier",
    "get_substructure_library",
    # Chemreps
    "download_chemreps",
    "get_chemreps_df",
    # Fingerprints
    "download_fps",
    "chemfp_load_fps",
    # Monomers
    "download_monomer_library",
    "get_monomer_library_root",
]

logger = logging.getLogger(__name__)

#: The default path inside the :mod:`pystow` directory
PYSTOW_PARTS = ["chembl"]


def latest() -> str:
    """Get the latest version of ChEBML as a string.

    :returns: The latest version string of ChEBML
    :raises ImportError: If :mod:`bioversions` is not installed.
    """
    try:
        import bioversions
    except ImportError:
        raise ImportError("Could not import `bioversions`. Install with `pip install bioversions`.")
    else:
        return bioversions.get_version("chembl")


def _download_helper(
    suffix: str,
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    *,
    return_version: bool,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL file with the given suffix is downloaded.

    :param suffix: The suffix of the file
    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded file. Otherwise, just return the path.
    """
    if version is None:
        version = latest()
    url = f"ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{version}/chembl_{version}{suffix}"
    path = pystow.ensure(*(prefix or PYSTOW_PARTS), version, url=url)
    if return_version:
        return version, path
    else:
        return path


def download_sqlite(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL SQLite dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded ``*.tar.gz`` file. Otherwise, just
        return the path.
    """
    return _download_helper(
        suffix="_sqlite.tar.gz", version=version, prefix=prefix, return_version=return_version
    )


def download_extract_sqlite(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL SQLite dump is downloaded and extracted.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded ChEMBLSQLite database file. Otherwise,
        just return the path.
    :raises FileNotFoundError: If no database file could be found in the
        extracted directories
    """
    version, path = cast(
        Tuple[str, Path], download_sqlite(version=version, prefix=prefix, return_version=True)
    )

    # Extraction will be done in the same directory as the download.
    # All ChEMBL SQLite dumps have the same internal folder structure,
    # so assume there's going to be a directory here
    directory = path.parent.joinpath("data")
    if not directory.is_dir():
        logger.info("unarchiving %s to %s", path, directory)
        with tarfile.open(path, mode="r", encoding="utf-8") as tar_file:
            tar_file.extractall(directory)
    else:
        logger.debug("did not re-unarchive %s to %s", path, directory)

    # Since the structure of the zip changes from version to version,
    # it's better to just walk through the unarchived folders recursively
    # and find the DB file
    for root, _dirs, files in os.walk(directory):
        for file in files:
            if not file.endswith(".db"):
                continue
            rv = Path(root).joinpath(file)
            if return_version:
                return version, rv
            else:
                return rv

    raise FileNotFoundError("could not find a .db file in the ChEMBL archive")


@contextmanager
def connect(version: Optional[str] = None, prefix: Optional[Sequence[str]] = None):
    """Ensure and connect to the database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :yields: The SQLite connection object.

    Example:
    .. code-block:: python

        import chembl_downloader

        with chembl_downloader.connect() as conn:
            with closing(conn.cursor()) as cursor:
                cursor.execute(...)
    """
    path = cast(Path, download_extract_sqlite(version=version, prefix=prefix, return_version=False))
    with closing(sqlite3.connect(path.as_posix())) as conn:
        yield conn


@contextmanager
def cursor(version: Optional[str] = None, prefix: Optional[Sequence[str]] = None):
    """Ensure, connect, and get a cursor from the database to the database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :yields: The SQLite cursor object.

    Example:
    .. code-block:: python

        import chembl_downloader

        with chembl_downloader.cursor() as cursor:
            cursor.execute(...)
    """
    with connect(version=version, prefix=prefix) as conn:
        with closing(conn.cursor()) as yv:
            yield yv


def query(
    sql: str, version: Optional[str] = None, prefix: Optional[Sequence[str]] = None, **kwargs
) -> "pandas.DataFrame":
    """Ensure the data is available, run the query, then put the results in a dataframe.

    :param sql: A SQL query string or table name
    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param kwargs: keyword arguments to pass through to :func:`pandas.read_sql`, such as
        ``index_col``.
    :return: A dataframe

    Example:
    .. code-block:: python

        import chembl_downloader
        from chembl_downloader.queries import ID_NAME_QUERY_EXAMPLE

        df = chembl_downloader.query(ID_NAME_QUERY_EXAMPLE)
    """
    import pandas as pd

    with connect(version=version, prefix=prefix) as con:
        return pd.read_sql(sql, con=con, **kwargs)


def download_fps(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL fingerprints file is downloaded.

    This file contains 2048 bit radius 2 morgan fingerprints.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded ``*.fps.gz`` file. Otherwise,
        just return the path.
    """
    return _download_helper(
        suffix=".fps.gz", version=version, prefix=prefix, return_version=return_version
    )


def chemfp_load_fps(
    version: Optional[str] = None, prefix: Optional[Sequence[str]] = None, **kwargs
):
    """Ensure the ChEMBL fingerprints file is downloaded and open with :func:`chemfp.load_fingerprints`.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param kwargs: Remaining keyword arguments are passed into :func:`chemfp.load_fingerprints`.
    :rtype: chemfp.arena.FingerprintArena
    """
    import chemfp

    path = download_fps(version=version, prefix=prefix, return_version=False)
    return chemfp.load_fingerprints(path, **kwargs)


def download_chemreps(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL chemical representations file is downloaded.

    This file is tab-separated and has four columns:

    1. ``chembl_id``
    2. ``canonical_smiles``
    3. ``standard_inchi``
    4. ``standard_inchi_key``

    If you want to directly parse it with :mod:`pandas`, use :func:`get_chemreps_df`.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded ``*_chemreps.txt.gz`` file. Otherwise,
        just return the path.
    """
    return _download_helper(
        suffix="_chemreps.txt.gz ", version=version, prefix=prefix, return_version=return_version
    )


def get_chemreps_df(
    version: Optional[str] = None, prefix: Optional[Sequence[str]] = None
) -> "pandas.DataFrame":
    """Download and parse the latest ChEMBL chemical representations file.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: A dataframe with four columns:
        1. ``chembl_id``
        2. ``canonical_smiles``
        3. ``standard_inchi``
        4. ``standard_inchi_key``
    """
    import pandas

    path = cast(Path, download_chemreps(version=version, prefix=prefix, return_version=False))
    df = pandas.read_csv(path, sep="\t", compression="gzip")
    return df


def download_sdf(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL SDF dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded ``*.sdf.gz`` file. Otherwise,
        just return the path.
    """
    return _download_helper(
        suffix=".sdf.gz", version=version, prefix=prefix, return_version=return_version
    )


def download_monomer_library(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL monomer library is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded ``*_monomer_library.xml`` file. Otherwise,
        just return the path.
    """
    return _download_helper(
        suffix="_monomer_library.xml", version=version, prefix=prefix, return_version=return_version
    )


def get_monomer_library_root(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
) -> ElementTree.Element:
    """Ensure the latest ChEMBL monomer library is downloaded and parse its root with :mod:`xml`.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: Return the root of the monomers XML tree, parsed
    """
    monomers_path = cast(Path, download_monomer_library(version=version, prefix=prefix, return_version=False))
    tree = ElementTree.parse(monomers_path)
    return tree.getroot()


@contextmanager
def supplier(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    **kwargs,
):
    """Get a :class:`rdkit.Chem.ForwardSDMolSupplier` for the given version of ChEMBL.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param kwargs: keyword arguments to pass through to :class:`rdkit.Chem.ForwardSDMolSupplier`, such as
        ``sanitize`` and ``removeHs``.

    Example:
    .. code-block:: python

        from rdkit import Chem

        import chembl_downloader

        data = []
        with chembl_downloader.supplier() as suppl:
            for i, mol in enumerate(suppl):
                if mol is None or mol.GetNumAtoms() > 50:
                    continue
                fp = Chem.PatternFingerprint(mol, fpSize=1024, tautomerFingerprints=True)
                smi = Chem.MolToSmiles(mol)
                data.append((smi, fp))
    """
    from rdkit import Chem

    path = cast(Path, download_sdf(version=version, prefix=prefix, return_version=False))
    with gzip.open(path) as file:
        yield Chem.ForwardSDMolSupplier(file, **kwargs)


def get_substructure_library(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    max_heavy: int = 75,
    **kwargs,
):
    """Get the ChEMBL substructure library.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`bioversions.get_version` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param max_heavy: The largest number of heavy atoms that are considered before skipping the molecule.
    :param kwargs: keyword arguments to pass through to :class:`rdkit.Chem.ForwardSDMolSupplier`, such as
        ``sanitize`` and ``removeHs`` via :func:`supplier`.
    :returns: A substructure library object
    :rtype: rdkit.Chem.rdSubstructLibrary.SubstructLibrary

    .. seealso::

        https://greglandrum.github.io/rdkit-blog/tutorial/substructure/2021/12/20/substructlibrary-search-order.html
    """
    # Requires minimum version of v2021.09
    from rdkit.Chem.rdSubstructLibrary import (
        CachedTrustedSmilesMolHolder,
        KeyFromPropHolder,
        SubstructLibrary,
        TautomerPatternHolder,
    )

    if version is None:
        version = latest()

    path = pystow.join(*(prefix or PYSTOW_PARTS), version, name="ssslib.pkl")
    if path.is_file():
        logger.info("loading substructure library from pickle: %s", path)
        with path.open("rb") as file:
            return pickle.load(file)

    molecule_holder = CachedTrustedSmilesMolHolder()
    tautomer_pattern_holder = TautomerPatternHolder()
    key_from_prop_holder = KeyFromPropHolder()
    library = SubstructLibrary(molecule_holder, tautomer_pattern_holder, key_from_prop_holder)
    with supplier(version=version, prefix=prefix, **kwargs) as suppl:
        for mol in tqdm(
            suppl, unit="molecule", unit_scale=True, desc="Building substructure library"
        ):
            if mol is None:
                continue
            if mol.GetNumHeavyAtoms() > max_heavy:  # skip huge molecules
                continue
            library.AddMol(mol)
    with path.open("wb") as file:
        pickle.dump(library, file, protocol=pickle.HIGHEST_PROTOCOL)
    return library
