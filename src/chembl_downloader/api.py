# -*- coding: utf-8 -*-

"""API for :mod:`chembl_downloader`."""

import ftplib
import gzip
import io
import logging
import os
import pickle
import sqlite3
import tarfile
from contextlib import closing, contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, List, Optional, Sequence, Tuple, Union, cast
from xml.etree import ElementTree

import pystow
from tqdm import tqdm

if TYPE_CHECKING:
    import pandas

__all__ = [
    "latest",
    "versions",
    "download_readme",
    "get_date",
    # Database
    "download_sqlite",
    "download_extract_sqlite",
    "connect",
    "cursor",
    "query",
    # SDF
    "download_sdf",
    "supplier",
    "iterate_smiles",
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
    # UniProt mappings
    "download_uniprot_mapping",
    "get_uniprot_mapping_df",
]

logger = logging.getLogger(__name__)

#: The default path inside the :mod:`pystow` directory
PYSTOW_PARTS = ["chembl"]
RELEASE_PREFIX = "* Release:"
DATE_PREFIX = "* Date:"


def _removeprefix(s: str, prefix: str) -> str:
    if s.startswith(prefix):
        return s[len(prefix) :]
    return s


def latest() -> str:
    """Get the latest version of ChEMBL as a string.

    :returns: The latest version string of ChEMBL
    :raises ValueError: If the latest README can not be parsed
    """
    bio = io.BytesIO()
    with ftplib.FTP("ftp.ebi.ac.uk") as ftp:
        ftp.login()
        ftp.retrbinary("RETR pub/databases/chembl/ChEMBLdb/latest/README", bio.write)
    bio.seek(0)
    for line in bio.read().decode("utf-8").split("\n"):
        if line.startswith(RELEASE_PREFIX):
            return _removeprefix(_removeprefix(line, RELEASE_PREFIX).strip(), "chembl_")
    raise ValueError("could not find latest ChEMBL version")


def versions() -> List[str]:
    """Get all versions of ChEMBL."""
    version_list = [str(i).zfill(2) for i in range(1, int(latest()) + 1)]
    # Side version in ChEMBL
    version_list.extend(["22_1", "24_1"])
    return sorted(version_list, reverse=True)


def _download_helper(
    suffix: str,
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    *,
    return_version: bool,
    filename_repeats_version: bool = True,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL file with the given suffix is downloaded.

    :param suffix: The suffix of the file
    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :param filename_repeats_version:
        True if filename contains ``chembl_<version>`` in the beginning. Set to false
        to allow downloading arbitrarily named files.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded file. Otherwise, just return the path.
    :raises ValueError:
        If file could not be downloaded
    """
    if version is None:
        version = latest()

    # for versions 22.1 and 24.1, it's important to canonicalize the version number
    # for versions < 10 it's important to left pad with a zero
    fmt_version = version.replace(".", "_").zfill(2)

    base = f"ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{fmt_version}"
    if filename_repeats_version:
        filename = f"chembl_{fmt_version}{suffix}"
    else:
        filename = suffix
    for url in [
        f"{base}/{filename}",
        f"{base}/archived/{filename}",
    ]:
        try:
            path = pystow.ensure(*(prefix or PYSTOW_PARTS), fmt_version, url=url)
        except IOError:
            continue
        if return_version:
            return version, path
        else:
            return path
    raise ValueError(f"could not find {filename} in data for ChEMBL {fmt_version} in {base}")


def download_sqlite(
    version: Optional[str] = None,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL SQLite dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
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
    *,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL SQLite dump is downloaded and extracted.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded ChEMBLSQLite database file. Otherwise,
        just return the path.
    :raises FileNotFoundError: If no database file could be found in the
        extracted directories
    """
    if version is not None:
        _directory = pystow.join(*(prefix or PYSTOW_PARTS), version)
        if _directory.is_dir():
            rv = _find_sqlite_file(_directory)
            if rv:
                if return_version:
                    return version, rv
                return rv

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
            tar_file.extractall(directory)  # noqa:S202
    else:
        logger.debug("did not re-unarchive %s to %s", path, directory)

    rv = _find_sqlite_file(directory)
    if rv is None:
        raise FileNotFoundError("could not find a .db file in the ChEMBL archive")
    elif return_version:
        return version, rv
    else:
        return rv


def _find_sqlite_file(directory: Union[str, Path]) -> Optional[Path]:
    # Since the structure of the zip changes from version to version,
    # it's better to just walk through the unarchived folders recursively
    # and find the DB file
    for root, _dirs, files in os.walk(directory):
        for file in files:
            if not file.endswith(".db"):
                continue
            rv = Path(root).joinpath(file)
            return rv
    return None


@contextmanager
def connect(version: Optional[str] = None, *, prefix: Optional[Sequence[str]] = None):
    """Ensure and connect to the database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
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
def cursor(version: Optional[str] = None, *, prefix: Optional[Sequence[str]] = None):
    """Ensure, connect, and get a cursor from the database to the database.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
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
    sql: str, version: Optional[str] = None, *, prefix: Optional[Sequence[str]] = None, **kwargs
) -> "pandas.DataFrame":
    """Ensure the data is available, run the query, then put the results in a dataframe.

    :param sql: A SQL query string or table name
    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
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
    *,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL fingerprints file is downloaded.

    This file contains 2048 bit radius 2 morgan fingerprints.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
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
    version: Optional[str] = None, *, prefix: Optional[Sequence[str]] = None, **kwargs
):
    """Ensure the ChEMBL fingerprints file is downloaded and open with :func:`chemfp.load_fingerprints`.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param kwargs: Remaining keyword arguments are passed into :func:`chemfp.load_fingerprints`.
    :return: A fingerpring arena object
    :rtype: chemfp.arena.FingerprintArena
    """
    import chemfp

    path = download_fps(version=version, prefix=prefix, return_version=False)
    return chemfp.load_fingerprints(path, **kwargs)


def download_chemreps(
    version: Optional[str] = None,
    *,
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
        :func:`latest` to look up the latest.
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
    version: Optional[str] = None, *, prefix: Optional[Sequence[str]] = None
) -> "pandas.DataFrame":
    """Download and parse the latest ChEMBL chemical representations file.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
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
    *,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL SDF dump is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
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
    *,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL monomer library is downloaded.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
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
    *,
    prefix: Optional[Sequence[str]] = None,
) -> ElementTree.Element:
    """Ensure the latest ChEMBL monomer library is downloaded and parse its root with :mod:`xml`.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: Return the root of the monomers XML tree, parsed
    """
    monomers_path = cast(
        Path, download_monomer_library(version=version, prefix=prefix, return_version=False)
    )
    tree = ElementTree.parse(monomers_path)
    return tree.getroot()


@contextmanager
def supplier(
    version: Optional[str] = None,
    *,
    prefix: Optional[Sequence[str]] = None,
    **kwargs,
):
    """Get a :class:`rdkit.Chem.ForwardSDMolSupplier` for the given version of ChEMBL.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param kwargs: keyword arguments to pass through to :class:`rdkit.Chem.ForwardSDMolSupplier`, such as
        ``sanitize`` and ``removeHs``.
    :yields: A supplier to be used in a context manager

    In the following example, a supplier is used to get fingerprints and SMILES.

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


def iterate_smiles(
    version: Optional[str] = None,
    *,
    prefix: Optional[Sequence[str]] = None,
    **kwargs,
) -> Iterable[str]:
    """Iterate over SMILES via RDKit."""
    from rdkit import Chem

    with supplier(version=version, prefix=prefix, **kwargs) as suppl:
        for mol in suppl:
            if mol is None:
                continue
            smiles = Chem.MolToSmiles(mol)
            if smiles:
                yield smiles


def get_substructure_library(
    version: Optional[str] = None,
    *,
    max_heavy: int = 75,
    prefix: Optional[Sequence[str]] = None,
    **kwargs,
):
    """Get the ChEMBL substructure library.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param max_heavy: The largest number of heavy atoms that are considered before skipping the molecule.
    :param prefix: The directory inside :mod:`pystow` to use
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


def download_readme(
    version: Optional[str] = None,
    *,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
) -> Union[Path, Tuple[str, Path]]:
    """Ensure the latest ChEMBL README.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded ``*.sdf.gz`` file. Otherwise,
        just return the path.
    """
    return _download_helper(
        "README",
        version=version,
        prefix=prefix,
        return_version=return_version,
        filename_repeats_version=False,
    )


def get_date(version: str, **kwargs) -> str:
    """Get the date of a given version."""
    path = cast(Path, download_readme(version=version, **kwargs))
    try:
        date_p = _removeprefix(
            next(line for line in path.read_text().splitlines() if line.startswith("* Date:")),
            DATE_PREFIX,
        ).lstrip()
    except StopIteration:
        return ""  # happens on 22.1 and 24.1
    else:
        day, month, year = date_p.split("/")
        return f"{year}-{month}-{day}"


def download_uniprot_mapping(
    version: Optional[str] = None,
    *,
    prefix: Optional[Sequence[str]] = None,
    return_version: bool = False,
):
    """Ensure the latest ChEMBL-UniProt target mapping TSV file.

    :param version: The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :param return_version: Should the version get returned? Turn this to true
        if you're looking up the latest version and want to reduce redundant code.
    :return: If ``return_version`` is true, return a pair of the version and the
        local file path to the downloaded ``*.txt`` file. Otherwise,
        just return the path.
    """
    return _download_helper(
        "chembl_uniprot_mapping.txt",
        version=version,
        prefix=prefix,
        return_version=return_version,
        filename_repeats_version=False,
    )


def get_uniprot_mapping_df(
    version: Optional[str] = None,
    *,
    prefix: Optional[Sequence[str]] = None,
) -> "pandas.DataFrame":
    """Download and parse the latest ChEMBL-UniProt target mapping TSV file.

    :param version:
        The version number of ChEMBL to get. If none specified, uses
        :func:`latest` to look up the latest.
    :param prefix: The directory inside :mod:`pystow` to use
    :return: A dataframe with four columns:

        1. ``uniprot_id``
        2. ``chembl_target_id``
        3. ``name``, the name from ChEMBL
        4. ``type``, which can have one of the following values:

           - ``CHIMERIC PROTEIN``
           - ``NUCLEIC-ACID``
           - ``PROTEIN COMPLEX``
           - ``PROTEIN COMPLEX GROUP``
           - ``PROTEIN FAMILY``
           - ``PROTEIN NUCLEIC-ACID COMPLEX``
           - ``PROTEIN-PROTEIN INTERACTION``
           - ``SELECTIVITY GROUP``
           - ``SINGLE PROTEIN``

    """
    import pandas as pd

    path = cast(
        Path, download_uniprot_mapping(version=version, prefix=prefix, return_version=False)
    )
    df = pd.read_csv(
        path,
        sep="\t",
        skiprows=1,
        header=None,
        names=["uniprot_id", "chembl_target_id", "name", "type"],
    )
    return df
