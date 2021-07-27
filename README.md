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
    <!--<a href="https://zenodo.org/badge/latestdoi/XXXX">
        <img src="https://zenodo.org/badge/XXXX.svg" alt="DOI" />
    </a>-->
</p>

Don't worry about downloading/extracting ChEMBL or versioning - just use ``chembl_downloader`` to write code that knows
how to download it and use it automatically.

## Installation

```bash
$ pip install chembl-downloader
```

## Download A Specific Version

```python
import chembl_downloader

path = chembl_downloader.download(version='28')
```

After it's been downloaded and extracted once, it's smart and doesn't need to download again. It gets stored
using [`pystow`](https://github.com/cthoyt/pystow) automatically in the `~/.data/chembl`
directory.

We'd like to implement something such that it could load directly into SQLite from the archive, but it appears this is
a [paid feature](https://sqlite.org/purchase/zipvfs).

## Download the Latest Version

First, you'll have to install [`bioversions`](https://github.com/cthoyt/bioversions)
with `pip install bioversions`, whose job it is to look up the latest version of many databases. Then, you can modify
the previous code slightly by omitting the `version` keyword argument:

```python
import chembl_downloader

path = chembl_downloader.download()
```

The `version` keyword argument is available for all functions in this package, but like the username and password will
be omitted for brevity.

## Automate Connection

Inside the archive is a single SQLite database file. Normally, people manually untar this folder then do something with
the resulting file. Don't do this, it's not reproducible!
Instead, the file can be downloaded and a connection can be opened automatically with:

```python
import chembl_downloader

with chembl_downloader.cursor(version='28') as cursor:
    cursor.execute(...)  # run your query string
    rows = cursor.fetchall()  # get your results
```

You now know everything I can teach you. Please use these tools to do re-usable, reproducible science!

## Store in a Different Place

If you want to store the data elsewhere using `pystow` (e.g., in [`pyobo`](https://github.com/pyobo/pyobo)
I also keep a copy of this file), you can use the `prefix` argument.

```python
import chembl_downloader

# It gets downloaded/extracted to 
# ~/.data/pyobo/raw/chembl/29/chembl_29/chembl_29_sqlite/chembl_29.db
path = chembl_downloader.download(prefix=['pyobo', 'raw', 'chembl'])
```

See the `pystow` [documentation](https://github.com/cthoyt/pystow#%EF%B8%8F-configuration) on configuring the storage
location further.

## Download via CLI

After installing, run the following CLI command to ensure it and send the path to stdout

```bash
$ chembl_downloader
```
