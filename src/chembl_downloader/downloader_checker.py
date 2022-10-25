# -*- coding: utf-8 -*-

"""Codebase for logging and checking the working of chembl_downloader."""

from typing import List, Tuple
from tqdm import tqdm
import chembl_downloader
import bioversions
import ftplib

from tabulate import tabulate


def get_version() -> List[str]:
    """Get all version of ChEMBL."""

    try:
        latest_chembl_version = bioversions.get_version("chembl")
    except ftplib.error_temp:  # error due to too many connected users
        latest_chembl_version = 31

    version_list = [
        str(i).zfill(2)
        for i in (range(1, int(latest_chembl_version) + 1))
    ]

    # Side version in ChEMBL
    version_list.extend(['22.1', '24.1'])
    version_list.sort()
    return version_list


def check_downloader(version_name: str) -> Tuple[str, str, str]:
    try:
        path = chembl_downloader.download_extract_sqlite(version=version_name)
    except Exception as e:
        return version_name, 'No', str(e)

    _count_sql = """
    SELECT COUNT( MOLECULE_DICTIONARY.chembl_id)
    FROM MOLECULE_DICTIONARY
    JOIN COMPOUND_STRUCTURES ON MOLECULE_DICTIONARY.molregno == COMPOUND_STRUCTURES.molregno
    WHERE molecule_dictionary.pref_name IS NOT NULL
"""
    total_compound = chembl_downloader.query(_count_sql)
    return version_name, 'Yes', ''


def main():
    chembl_versions = get_version()

    table = [
        ['ChEMBL Version', 'Downloader working', 'Error']
    ]

    for version in tqdm(chembl_versions):
        data = check_downloader(version_name=version)
        table.append(list(data))

    with open('../README.md', 'a') as file:
        print('\n \n', file=file)
        print('## ChEMBL Downloader Status', file=file)
        print('\n \n', file=file)
        print(tabulate(table, headers="firstrow", tablefmt='pipe'), file=file)
        print('\n', file=file)


if __name__ == '__main__':
    main()
