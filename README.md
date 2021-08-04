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

## Installation

```bash
$ pip install chembl-downloader
```

## Usage

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

If you'd like to contribute, there's a submodule called `chembl_downloader.queries`
where you can add an SQL query along with a description of what it does for easy importing.
