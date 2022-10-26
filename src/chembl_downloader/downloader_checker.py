# -*- coding: utf-8 -*-

"""Codebase for logging and checking the working of chembl_downloader."""

from typing import List, Tuple
from tqdm import tqdm
import chembl_downloader
import bioversions
import ftplib

from tabulate import tabulate
from chembl_downloader.queries import COUNT_QUERY


def get_version() -> List[str]:
    """Get all version of ChEMBL."""

    try:
        latest_chembl_version = bioversions.get_version("chembl")
    except ftplib.error_temp:  # error due to too many connected users
        latest_chembl_version = 31

    version_list = [
        str(i)
        for i in range(1, int(latest_chembl_version) + 1)
    ]

    # Side version in ChEMBL
    version_list.extend(['22.1', '24.1'])
    return sorted(version_list, reverse=True)


def check_downloader(version_name: str) -> Tuple[str, str, str]:
    try:
        path = chembl_downloader.download_extract_sqlite(version=version_name)
    except Exception as e:
        return version_name, 'No', str(e)

    total_compound = chembl_downloader.query(COUNT_QUERY)
    return version_name, 'Yes', ''


def main():
    chembl_versions = get_version()

    headers = ['ChEMBL Version', 'Downloader working', 'Error']
    table = [
        check_downloader(version_name=version)
        for version in tqdm(chembl_versions)
    ]
    with open('../README.md', 'a') as file:
        print('\n \n', file=file)
        print('## ChEMBL Downloader Status', file=file)
        print('\n \n', file=file)
        print(tabulate(table, headers=headers, tablefmt='github'), file=file)
        print('\n', file=file)


if __name__ == '__main__':
    main()
