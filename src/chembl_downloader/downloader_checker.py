# -*- coding: utf-8 -*-

"""Codebase for logging and checking the working of chembl_downloader."""

from typing import Tuple

from tabulate import tabulate
from tqdm import tqdm

import chembl_downloader
from chembl_downloader import versions
from chembl_downloader.queries import COUNT_QUERY_SQL


def check_downloader(version_name: str) -> Tuple[str, str]:
    """Test downloader for specific ChEMBL version."""
    try:
        chembl_downloader.download_extract_sqlite(version=version_name)
    except Exception:
        return version_name, "-"

    total_compounds = chembl_downloader.query(COUNT_QUERY_SQL, version=version_name)["count"][0]
    return version_name, total_compounds


def main():
    """Run ChEMBL downloader on all versions fo ChEMBL."""
    chembl_versions = versions()

    headers = ["ChEMBL Version", "Total compounds"]
    table = [check_downloader(version_name=version) for version in tqdm(chembl_versions)]

    with open("../../README.md", "a") as file:
        file.writelines("\n \n")
        file.writelines("## ChEMBL Downloader Status \n")
        file.writelines(
            "> **_NOTE:_** `chembl_downloader` runs only on SQLite distributions of ChEMBL. \
            If the run was unsuccessful, it is because earlier versions of ChEMBL were provided \
            on a MongoDB distribution.",
        )
        file.writelines("\n \n")
        file.writelines(tabulate(table, headers=headers, tablefmt="github"))
        file.writelines("\n")


if __name__ == "__main__":
    main()
