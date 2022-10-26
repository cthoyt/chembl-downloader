# -*- coding: utf-8 -*-

"""Codebase for logging and checking the working of chembl_downloader."""

from typing import Tuple
from tqdm import tqdm
import chembl_downloader

from tabulate import tabulate
from chembl_downloader.queries import COUNT_QUERY
from chembl_downloader import versions


def check_downloader(version_name: str) -> Tuple[str, str, str]:
    try:
        path = chembl_downloader.download_extract_sqlite(version=version_name)
    except Exception as e:
        return version_name, 'No', str(e)

    total_compound = chembl_downloader.query(COUNT_QUERY)
    return version_name, 'Yes', ''


def main():
    chembl_versions = versions()

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
