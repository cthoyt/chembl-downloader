<h1 align="center">
    chembl_downloader
</h1>

<p align="center">
    <a href="https://pypi.org/project/chembl_downloader">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/chembl_downloader" />
    </a>
    <a href="https://pypi.org/project/chembl_downloader">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/chembl_downloader" />
    </a>
    <a href="https://github.com/cthoyt/chembl_downloader/blob/main/LICENSE">
        <img alt="PyPI - License" src="https://img.shields.io/pypi/l/chembl_downloader" />
    </a>
    <a href="https://zenodo.org/badge/latestdoi/390113187">
        <img src="https://zenodo.org/badge/390113187.svg" alt="DOI" />
    </a>
    <a href="https://github.com/psf/black">
        <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Code style: black" />
    </a>
</p>

Don't worry about downloading/extracting ChEMBL or versioning - just use ``chembl_downloader`` to write code that knows
how to download it and use it automatically.

Install with:

```bash
$ pip install chembl-downloader
```

## Database Usage

### Download A Specific Version

```python
import chembl_downloader

path = chembl_downloader.download_extract_sqlite(version='28')
```

After it's been downloaded and extracted once, it's smart and does not need to download again. It gets stored
using [`pystow`](https://github.com/cthoyt/pystow) automatically in the `~/.data/chembl`
directory.

We'd like to implement something such that it could load directly into SQLite from the archive, but it appears this is
a [paid feature](https://sqlite.org/purchase/zipvfs).

### Download the Latest Version

First, you'll have to install [`bioversions`](https://github.com/cthoyt/bioversions)
with `pip install bioversions`, whose job it is to look up the latest version of many databases. Then, you can modify
the previous code slightly by omitting the `version` keyword argument:

```python
import chembl_downloader

path = chembl_downloader.download_extract_sqlite()
```

The `version` keyword argument is available for all functions in this package (e.g., including
`connect()`, `cursor()`, and `query()`), but will be omitted below for brevity.

### Automate Connection

Inside the archive is a single SQLite database file. Normally, people manually untar this folder then do something with
the resulting file. Don't do this, it's not reproducible!
Instead, the file can be downloaded and a connection can be opened automatically with:

```python
import chembl_downloader

with chembl_downloader.connect() as conn:
    with conn.cursor() as cursor:
        cursor.execute(...)  # run your query string
        rows = cursor.fetchall()  # get your results
```

The `cursor()` function provides a convenient wrapper around this operation:

```python
import chembl_downloader

with chembl_downloader.cursor() as cursor:
    cursor.execute(...)  # run your query string
    rows = cursor.fetchall()  # get your results
```

### Run a query and get a pandas DataFrame

The most powerful function is `query()` which builds on the previous `connect()` function in combination
with [`pandas.read_sql`](https://pandas.pydata.org/docs/reference/api/pandas.read_sql.html)
to make a query and load the results into a pandas DataFrame for any downstream use.

```python
import chembl_downloader

sql = """
SELECT
    MOLECULE_DICTIONARY.chembl_id,
    MOLECULE_DICTIONARY.pref_name
FROM MOLECULE_DICTIONARY
JOIN COMPOUND_STRUCTURES ON MOLECULE_DICTIONARY.molregno == COMPOUND_STRUCTURES.molregno
WHERE molecule_dictionary.pref_name IS NOT NULL
LIMIT 5
"""

df = chembl_downloader.query(sql)
df.to_csv(..., sep='\t', index=False)
```

Suggestion 1: use `pystow` to make a reproducible file path that's portable to other people's machines
(e.g., it doesn't have your username in the path).

Suggestion 2: RDKit is now pip-installable with `pip install rdkit-pypi`, which means most users don't have to muck
around with complicated conda environments and configurations. One of the powerful but understated tools in RDKit is
the [rdkit.Chem.PandasTools](https://rdkit.org/docs/source/rdkit.Chem.PandasTools.html)
module.

### Access an RDKit supplier over entries in the SDF dump

This example is a bit more fit-for-purpose than the last two. The `supplier()` function makes sure that the latest SDF
dump is downloaded and loads it from the gzip file into a `rdkit.Chem.ForwardSDMolSupplier`
using a context manager to make sure the file doesn't get closed until after parsing is done. Like the previous
examples, it can also explicitly take a `version`.

```python
from rdkit import Chem

import chembl_downloader

with chembl_downloader.supplier() as suppl:
    data = []
    for i, mol in enumerate(suppl):
        if mol is None or mol.GetNumAtoms() > 50:
            continue
        fp = Chem.PatternFingerprint(mol, fpSize=1024, tautomerFingerprints=True)
        smi = Chem.MolToSmiles(mol)
        data.append((smi, fp))
```

This example was adapted from Greg Landrum's RDKit blog post
on [generalized substructure search](https://greglandrum.github.io/rdkit-blog/tutorial/substructure/2021/08/03/generalized-substructure-search.html).

## SDF Usage

### Get an RDKit substructure library

Building on the `supplier()` function, the `get_substructure_library()`
makes the preparation of a [substructure library](https://www.rdkit.org/docs/cppapi/classRDKit_1_1SubstructLibrary.html)
automated and reproducible. Additionally, it caches the results of the build,
which takes on the order of tens of minutes, only has to be done once and future
loading from a pickle object takes on the order of seconds.

The implementation was inspired by Greg Landrum's RDKit blog post,
[Some new features in the SubstructLibrary](https://greglandrum.github.io/rdkit-blog/tutorial/substructure/2021/12/20/substructlibrary-search-order.html).
The following example shows how it can be used to accomplish some of the first
tasks presented in the post:

```python
from rdkit import Chem

import chembl_downloader

library = chembl_downloader.get_substructure_library()
query = Chem.MolFromSmarts('[O,N]=C-c:1:c:c:n:c:c:1')
matches = library.GetMatches(query)
```

## Morgan Fingerprints Usage

### Get the Morgan Fingerprint file

ChEMBL makes a file containing pre-computed 2048 bit radius 2 morgan
fingerprints for each molecule available. It can be downloaded using:

```python
import chembl_downloader

path = chembl_downloader.download_fps()
```

The `version` and other keyword arguments are also valid for this function.

### Load fingerprints with [`chemfp`](https://chemfp.com/)

The following wraps the `download_fps` function with `chemfp`'s fingerprint
loader:

```python
import chembl_downloader

arena = chembl_downloader.chemfp_load_fps()
```

The `version` and other keyword arguments are also valid for this function.
More information on working with the `arena` object can be found
[here](https://chemfp.readthedocs.io/en/latest/using-api.html#working-with-a-fingerprintarena).

## Extras

### Store in a Different Place

If you want to store the data elsewhere using `pystow` (e.g., in [`pyobo`](https://github.com/pyobo/pyobo)
I also keep a copy of this file), you can use the `prefix` argument.

```python
import chembl_downloader

# It gets downloaded/extracted to 
# ~/.data/pyobo/raw/chembl/29/chembl_29/chembl_29_sqlite/chembl_29.db
path = chembl_downloader.download_extract_sqlite(prefix=['pyobo', 'raw', 'chembl'])
```

See the `pystow` [documentation](https://github.com/cthoyt/pystow#%EF%B8%8F-configuration) on configuring the storage
location further.

The `prefix` keyword argument is available for all functions in this package (e.g., including
`connect()`, `cursor()`, and `query()`).

### Download via CLI

After installing, run the following CLI command to ensure it and send the path to stdout

```bash
$ chembl_downloader
```

Use `--test` to show two example queries

```bash
$ chembl_downloader --test
```

## Contributing

Please read the contribution guidelines in [CONTRIBUTING.md](.github/CONTRIBUTING.md).

If you'd like to contribute, there's a submodule called `chembl_downloader.queries`
where you can add a useful SQL queries along with a description of what it does for easy
importing and reuse.

## Statistics

|   ChEMBL Version | Release Date   |
|------------------|----------------|
|               31 | 2022-07-12     |
|               30 | 2022-02-22     |
|               29 | 2021-07-01     |
|               28 | 2021-01-15     |
|               27 | 2020-05-18     |
|               26 | 2020-02-14     |
|               25 | 2019-02-01     |
|             24_1 | 2018-05-01     |
|               24 |                |
|               23 | 2017-05-18     |
|             22_1 | 2016-11-17     |
|               22 |                |
|               21 | 2015-02-12     |
|               20 | 2015-02-03     |
|               19 | 2014-07-2333   |
|               18 | 2014-04-02     |
|               17 | 2013-09-16     |
|               16 | 2013-055555-15 |
|               15 | 2013-01-30     |
|               14 | 2012 -07-18    |
|               13 | 2012-02-29     |
|               12 | 2011-11-30     |
|               11 | 2011-06-07     |
|               10 | 2011-06-07     |
|               09 | 2011-01-04     |
|               08 | 2010-11-05     |
|               07 | 2010-09-03     |
|               06 | 2010-09-03     |
|               05 | 2010-06-07     |
|               04 | 2010-05-26     |
|               03 | 2010-04-30     |
|               02 | 2009-12-07     |
|               01 | 2009-10-28     |

 

## ChEMBL Downloader Status

 

|   ChEMBL Version | Downloader working   | Error                                                                                                                                        | Total compounds   |
|------------------|----------------------|----------------------------------------------------------------------------------------------------------------------------------------------|-------------------|
|               31 | Yes                  |                                                                                                                                              | 41585             |
|               30 | Yes                  |                                                                                                                                              | 41585             |
|               29 | No                   | could not find chembl_29_sqlite.tar.gz in data for ChEMBL 29 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_29       | -                 |
|               28 | No                   | could not find chembl_28_sqlite.tar.gz in data for ChEMBL 28 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_28       | -                 |
|               27 | No                   | could not find chembl_27_sqlite.tar.gz in data for ChEMBL 27 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_27       | -                 |
|               26 | No                   | could not find chembl_26_sqlite.tar.gz in data for ChEMBL 26 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_26       | -                 |
|               25 | No                   | could not find chembl_25_sqlite.tar.gz in data for ChEMBL 25 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_25       | -                 |
|             24_1 | No                   | could not find chembl_24_1_sqlite.tar.gz in data for ChEMBL 24_1 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_24_1 | -                 |
|               24 | No                   | could not find chembl_24_sqlite.tar.gz in data for ChEMBL 24 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_24       | -                 |
|               23 | No                   | could not find chembl_23_sqlite.tar.gz in data for ChEMBL 23 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_23       | -                 |
|             22_1 | No                   | could not find chembl_22_1_sqlite.tar.gz in data for ChEMBL 22_1 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_22_1 | -                 |
|               22 | No                   | could not find chembl_22_sqlite.tar.gz in data for ChEMBL 22 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_22       | -                 |
|               21 | No                   | could not find chembl_21_sqlite.tar.gz in data for ChEMBL 21 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_21       | -                 |
|               20 | No                   | could not find chembl_20_sqlite.tar.gz in data for ChEMBL 20 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_20       | -                 |
|               19 | No                   | could not find chembl_19_sqlite.tar.gz in data for ChEMBL 19 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_19       | -                 |
|               18 | No                   | could not find chembl_18_sqlite.tar.gz in data for ChEMBL 18 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_18       | -                 |
|               17 | No                   | could not find chembl_17_sqlite.tar.gz in data for ChEMBL 17 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_17       | -                 |
|               16 | No                   | could not find chembl_16_sqlite.tar.gz in data for ChEMBL 16 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_16       | -                 |
|               15 | No                   | could not find chembl_15_sqlite.tar.gz in data for ChEMBL 15 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_15       | -                 |
|               14 | No                   | could not find chembl_14_sqlite.tar.gz in data for ChEMBL 14 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_14       | -                 |
|               13 | No                   | could not find chembl_13_sqlite.tar.gz in data for ChEMBL 13 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_13       | -                 |
|               12 | No                   | could not find chembl_12_sqlite.tar.gz in data for ChEMBL 12 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_12       | -                 |
|               11 | No                   | could not find chembl_11_sqlite.tar.gz in data for ChEMBL 11 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_11       | -                 |
|               10 | No                   | could not find chembl_10_sqlite.tar.gz in data for ChEMBL 10 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_10       | -                 |
|               09 | No                   | could not find chembl_09_sqlite.tar.gz in data for ChEMBL 09 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_09       | -                 |
|               08 | No                   | could not find chembl_08_sqlite.tar.gz in data for ChEMBL 08 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_08       | -                 |
|               07 | No                   | could not find chembl_07_sqlite.tar.gz in data for ChEMBL 07 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_07       | -                 |
|               06 | No                   | could not find chembl_06_sqlite.tar.gz in data for ChEMBL 06 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_06       | -                 |
|               05 | No                   | could not find chembl_05_sqlite.tar.gz in data for ChEMBL 05 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_05       | -                 |
|               04 | No                   | could not find chembl_04_sqlite.tar.gz in data for ChEMBL 04 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_04       | -                 |
|               03 | No                   | could not find chembl_03_sqlite.tar.gz in data for ChEMBL 03 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_03       | -                 |
|               02 | No                   | could not find chembl_02_sqlite.tar.gz in data for ChEMBL 02 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_02       | -                 |
|               01 | No                   | could not find chembl_01_sqlite.tar.gz in data for ChEMBL 01 in https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_01       | -                 |


