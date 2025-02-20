<!--
<p align="center">
  <img src="https://github.com/cthoyt/chembl-downloader/raw/main/docs/source/logo.png" height="150">
</p>
-->

<h1 align="center">
  ChEMBL Downloader
</h1>

<p align="center">
    <a href="https://github.com/cthoyt/chembl-downloader/actions/workflows/tests.yml">
        <img alt="Tests" src="https://github.com/cthoyt/chembl-downloader/actions/workflows/tests.yml/badge.svg" /></a>
    <a href="https://pypi.org/project/chembl_downloader">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/chembl_downloader" /></a>
    <a href="https://pypi.org/project/chembl_downloader">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/chembl_downloader" /></a>
    <a href="https://github.com/cthoyt/chembl-downloader/blob/main/LICENSE">
        <img alt="PyPI - License" src="https://img.shields.io/pypi/l/chembl_downloader" /></a>
    <a href='https://chembl_downloader.readthedocs.io/en/latest/?badge=latest'>
        <img src='https://readthedocs.org/projects/chembl_downloader/badge/?version=latest' alt='Documentation Status' /></a>
    <a href="https://codecov.io/gh/cthoyt/chembl-downloader/branch/main">
        <img src="https://codecov.io/gh/cthoyt/chembl-downloader/branch/main/graph/badge.svg" alt="Codecov status" /></a>  
    <a href="https://github.com/cthoyt/cookiecutter-python-package">
        <img alt="Cookiecutter template from @cthoyt" src="https://img.shields.io/badge/Cookiecutter-snekpack-blue" /></a>
    <a href="https://github.com/astral-sh/ruff">
        <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json" alt="Ruff" style="max-width:100%;"></a>
    <a href="https://github.com/cthoyt/chembl-downloader/blob/main/.github/CODE_OF_CONDUCT.md">
        <img src="https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg" alt="Contributor Covenant"/></a>
    <a href="https://zenodo.org/badge/latestdoi/390113187">
        <img src="https://zenodo.org/badge/390113187.svg" alt="DOI"></a>
</p>

Reproducibly download, open, parse, and query ChEMBL.

Don't worry about downloading/extracting ChEMBL or versioning - just use
`chembl_downloader` to write code that knows how to download it and use it
automatically.

## 💪 Getting Started

Download an extract the SQLite dump using the following:

```python
import chembl_downloader

path = chembl_downloader.download_extract_sqlite(version='28')
```

After it's been downloaded and extracted once, it's smart and does not need to
download again. It gets stored using
[`pystow`](https://github.com/cthoyt/pystow) automatically in the
`~/.data/chembl` directory.

Full technical documentation can be found on
[ReadTheDocs](https://chembl-downloader.readthedocs.io). Tutorials can be found
in Jupyter notebooks in the [notebooks/](notebooks/) directory of the
repository.

### Download the Latest Version

You can modify the previous code slightly by omitting the `version` keyword
argument to automatically find the latest version of ChEMBL:

```python
import chembl_downloader

path = chembl_downloader.download_extract_sqlite()
```

The `version` keyword argument is available for all functions in this package
(e.g., including `connect()`, `cursor()`, and `query()`), but will be omitted
below for brevity.

### Automatically Connect to SQLite

Inside the archive is a single SQLite database file. Normally, people manually
untar this folder then do something with the resulting file. Don't do this, it's
not reproducible! Instead, the file can be downloaded and a connection can be
opened automatically with:

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

### Run a Query and Get a Pandas DataFrame

The most powerful function is `query()` which builds on the previous `connect()`
function in combination with
[`pandas.read_sql`](https://pandas.pydata.org/docs/reference/api/pandas.read_sql.html)
to make a query and load the results into a pandas DataFrame for any downstream
use.

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

Suggestion 1: use `pystow` to make a reproducible file path that's portable to
other people's machines (e.g., it doesn't have your username in the path).

Suggestion 2: RDKit is now pip-installable with `pip install rdkit-pypi`, which
means most users don't have to muck around with complicated conda environments
and configurations. One of the powerful but understated tools in RDKit is the
[rdkit.Chem.PandasTools](https://rdkit.org/docs/source/rdkit.Chem.PandasTools.html)
module.

### SDF Usage

#### Access an RDKit supplier over entries in the SDF dump

This example is a bit more fit-for-purpose than the last two. The `supplier()`
function makes sure that the latest SDF dump is downloaded and loads it from the
gzip file into a `rdkit.Chem.ForwardSDMolSupplier` using a context manager to
make sure the file doesn't get closed until after parsing is done. Like the
previous examples, it can also explicitly take a `version`.

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

This example was adapted from Greg Landrum's RDKit blog post on
[generalized substructure search](https://greglandrum.github.io/rdkit-blog/tutorial/substructure/2021/08/03/generalized-substructure-search.html).

#### Iterate over SMILES

This example uses the `supplier()` method and RDKit to get SMILES strings from
molecules in ChEMBL's SDF file. If you want direct access to the RDKit molecule
objects, use `supplier()`.

```python
import chembl_downloader

for smiles in chembl_downloader.iterate_smiles():
    print(smiles)
```

### Get an RDKit substructure library

Building on the `supplier()` function, the `get_substructure_library()` makes
the preparation of a
[substructure library](https://www.rdkit.org/docs/cppapi/classRDKit_1_1SubstructLibrary.html)
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

### Morgan Fingerprints Usage

#### Get the Morgan Fingerprint file

ChEMBL makes a file containing pre-computed 2048 bit radius 2 morgan
fingerprints for each molecule available. It can be downloaded using:

```python
import chembl_downloader

path = chembl_downloader.download_fps()
```

The `version` and other keyword arguments are also valid for this function.

#### Load fingerprints with [`chemfp`](https://chemfp.com/)

The following wraps the `download_fps` function with `chemfp`'s fingerprint
loader:

```python
import chembl_downloader

arena = chembl_downloader.chemfp_load_fps()
```

The `version` and other keyword arguments are also valid for this function. More
information on working with the `arena` object can be found
[here](https://chemfp.readthedocs.io/en/latest/using-api.html#working-with-a-fingerprintarena).

### Command Line Interface

After installing, run the following CLI command to ensure it and send the path
to stdout

```console
$ chembl_downloader
```

Use `--test` to show two example queries

```console
$ chembl_downloader --test
```

## Configuration

If you want to store the data elsewhere using `pystow` (e.g., in
[`pyobo`](https://github.com/pyobo/pyobo) I also keep a copy of this file), you
can use the `prefix` argument.

```python
import chembl_downloader

# It gets downloaded/extracted to
# ~/.data/pyobo/raw/chembl/29/chembl_29/chembl_29_sqlite/chembl_29.db
path = chembl_downloader.download_extract_sqlite(prefix=['pyobo', 'raw', 'chembl'])
```

See the `pystow`
[documentation](https://github.com/cthoyt/pystow#%EF%B8%8F-configuration) on
configuring the storage location further.

The `prefix` keyword argument is available for all functions in this package
(e.g., including `connect()`, `cursor()`, and `query()`).

## 🚀 Installation

The most recent release can be installed from
[PyPI](https://pypi.org/project/chembl_downloader/) with uv:

```console
$ uv pip install chembl_downloader
```

or with pip:

```console
$ python3 -m pip install chembl_downloader
```

The most recent code and data can be installed directly from GitHub with uv:

```console
$ uv --preview pip install git+https://github.com/cthoyt/chembl-downloader.git
```

or with pip:

```console
$ UV_PREVIEW=1 python3 -m pip install git+https://github.com/cthoyt/chembl-downloader.git
```

Note that this requires setting `UV_PREVIEW` mode enabled until the uv build
backend becomes a stable feature.

## Users

See
[who's using `chembl-downloader`](https://github.com/search?q=chembl_downloader+-user%3Acthoyt+-is%3Afork&type=code).

## Statistics and Compatibility

`chembl-downloader` is compatible with all versions of ChEMBL. However, some
files are not available for all versions. For example, the SQLite version of the
database was first added in release 21 (2015-02-12).

| ChEMBL Version | Release Date   | Total Named Compounds _from SQLite_ |
| -------------- | -------------- | ----------------------------------: |
| 31             | 2022-07-12     |                              41,585 |
| 30             | 2022-02-22     |                              41,549 |
| 29             | 2021-07-01     |                              41,383 |
| 28             | 2021-01-15     |                              41,049 |
| 27             | 2020-05-18     |                              40,834 |
| 26             | 2020-02-14     |                              40,822 |
| 25             | 2019-02-01     |                              39,885 |
| 24_1           | 2018-05-01     |                              39,877 |
| 24             |                |                                     |
| 23             | 2017-05-18     |                              39,584 |
| 22_1           | 2016-11-17     |                                     |
| 22             |                |                              39,422 |
| 21             | 2015-02-12     |                              39,347 |
| 20             | 2015-02-03     |                                   - |
| 19             | 2014-07-2333   |                                   - |
| 18             | 2014-04-02     |                                   - |
| 17             | 2013-09-16     |                                   - |
| 16             | 2013-055555-15 |                                   - |
| 15             | 2013-01-30     |                                   - |
| 14             | 2012 -07-18    |                                   - |
| 13             | 2012-02-29     |                                   - |
| 12             | 2011-11-30     |                                   - |
| 11             | 2011-06-07     |                                   - |
| 10             | 2011-06-07     |                                   - |
| 09             | 2011-01-04     |                                   - |
| 08             | 2010-11-05     |                                   - |
| 07             | 2010-09-03     |                                   - |
| 06             | 2010-09-03     |                                   - |
| 05             | 2010-06-07     |                                   - |
| 04             | 2010-05-26     |                                   - |
| 03             | 2010-04-30     |                                   - |
| 02             | 2009-12-07     |                                   - |
| 01             | 2009-10-28     |                                   - |

## 👐 Contributing

Contributions, whether filing an issue, making a pull request, or forking, are
appreciated. See
[CONTRIBUTING.md](https://github.com/cthoyt/chembl-downloader/blob/master/.github/CONTRIBUTING.md)
for more information on getting involved.

## 👋 Attribution

### ⚖️ License

The code in this package is licensed under the MIT License.

<!--
### 📖 Citation

Citation goes here!
-->

<!--
### 🎁 Support

This project has been supported by the following organizations (in alphabetical order):

- [Biopragmatics Lab](https://biopragmatics.github.io)

-->

<!--
### 💰 Funding

This project has been supported by the following grants:

| Funding Body  | Program                                                      | Grant Number |
|---------------|--------------------------------------------------------------|--------------|
| Funder        | [Grant Name (GRANT-ACRONYM)](https://example.com/grant-link) | ABCXYZ       |
-->

### 🍪 Cookiecutter

This package was created with
[@audreyfeldroy](https://github.com/audreyfeldroy)'s
[cookiecutter](https://github.com/cookiecutter/cookiecutter) package using
[@cthoyt](https://github.com/cthoyt)'s
[cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack)
template.

## 🛠️ For Developers

<details>
  <summary>See developer instructions</summary>

The final section of the README is for if you want to get involved by making a
code contribution.

### Development Installation

To install in development mode, use the following:

```console
$ git clone git+https://github.com/cthoyt/chembl-downloader.git
$ cd chembl-downloader
$ uv --preview pip install -e .
```

Alternatively, install using pip:

```console
$ UV_PREVIEW=1 python3 -m pip install -e .
```

Note that this requires setting `UV_PREVIEW` mode enabled until the uv build
backend becomes a stable feature.

### Updating Package Boilerplate

This project uses `cruft` to keep boilerplate (i.e., configuration, contribution
guidelines, documentation configuration) up-to-date with the upstream
cookiecutter package. Install cruft with either `uv tool install cruft` or
`python3 -m pip install cruft` then run:

```console
$ cruft update
```

More info on Cruft's update command is available
[here](https://github.com/cruft/cruft?tab=readme-ov-file#updating-a-project).

### 🥼 Testing

After cloning the repository and installing `tox` with
`uv tool install tox --with tox-uv` or `python3 -m pip install tox tox-uv`, the
unit tests in the `tests/` folder can be run reproducibly with:

```console
$ tox -e py
```

Additionally, these tests are automatically re-run with each commit in a
[GitHub Action](https://github.com/cthoyt/chembl-downloader/actions?query=workflow%3ATests).

### 📖 Building the Documentation

The documentation can be built locally using the following:

```console
$ git clone git+https://github.com/cthoyt/chembl-downloader.git
$ cd chembl-downloader
$ tox -e docs
$ open docs/build/html/index.html
```

The documentation automatically installs the package as well as the `docs` extra
specified in the [`pyproject.toml`](pyproject.toml). `sphinx` plugins like
`texext` can be added there. Additionally, they need to be added to the
`extensions` list in [`docs/source/conf.py`](docs/source/conf.py).

The documentation can be deployed to [ReadTheDocs](https://readthedocs.io) using
[this guide](https://docs.readthedocs.io/en/stable/intro/import-guide.html). The
[`.readthedocs.yml`](.readthedocs.yml) YAML file contains all the configuration
you'll need. You can also set up continuous integration on GitHub to check not
only that Sphinx can build the documentation in an isolated environment (i.e.,
with `tox -e docs-test`) but also that
[ReadTheDocs can build it too](https://docs.readthedocs.io/en/stable/pull-requests.html).

#### Configuring ReadTheDocs

1. Log in to ReadTheDocs with your GitHub account to install the integration at
   https://readthedocs.org/accounts/login/?next=/dashboard/
2. Import your project by navigating to https://readthedocs.org/dashboard/import
   then clicking the plus icon next to your repository
3. You can rename the repository on the next screen using a more stylized name
   (i.e., with spaces and capital letters)
4. Click next, and you're good to go!

### 📦 Making a Release

#### Configuring Zenodo

[Zenodo](https://zenodo.org) is a long-term archival system that assigns a DOI
to each release of your package.

1. Log in to Zenodo via GitHub with this link:
   https://zenodo.org/oauth/login/github/?next=%2F. This brings you to a page
   that lists all of your organizations and asks you to approve installing the
   Zenodo app on GitHub. Click "grant" next to any organizations you want to
   enable the integration for, then click the big green "approve" button. This
   step only needs to be done once.
2. Navigate to https://zenodo.org/account/settings/github/, which lists all of
   your GitHub repositories (both in your username and any organizations you
   enabled). Click the on/off toggle for any relevant repositories. When you
   make a new repository, you'll have to come back to this

After these steps, you're ready to go! After you make "release" on GitHub (steps
for this are below), you can navigate to
https://zenodo.org/account/settings/github/repository/cthoyt/chembl-downloader
to see the DOI for the release and link to the Zenodo record for it.

#### Registering with the Python Package Index (PyPI)

You only have to do the following steps once.

1. Register for an account on the
   [Python Package Index (PyPI)](https://pypi.org/account/register)
2. Navigate to https://pypi.org/manage/account and make sure you have verified
   your email address. A verification email might not have been sent by default,
   so you might have to click the "options" dropdown next to your address to get
   to the "re-send verification email" button
3. 2-Factor authentication is required for PyPI since the end of 2023 (see this
   [blog post from PyPI](https://blog.pypi.org/posts/2023-05-25-securing-pypi-with-2fa/)).
   This means you have to first issue account recovery codes, then set up
   2-factor authentication
4. Issue an API token from https://pypi.org/manage/account/token

#### Configuring your machine's connection to PyPI

You have to do the following steps once per machine.

```console
$ uv tool install keyring
$ keyring set https://upload.pypi.org/legacy/ __token__
$ keyring set https://test.pypi.org/legacy/ __token__
```

Note that this deprecates previous workflows using `.pypirc`.

#### Uploading to PyPI

After installing the package in development mode and installing `tox` with
`uv tool install tox --with tox-uv` or `python3 -m pip install tox tox-uv`, run
the following from the console:

```console
$ tox -e finish
```

This script does the following:

1. Uses [bump-my-version](https://github.com/callowayproject/bump-my-version) to
   switch the version number in the `pyproject.toml`, `CITATION.cff`,
   `src/chembl_downloader/version.py`, and
   [`docs/source/conf.py`](docs/source/conf.py) to not have the `-dev` suffix
2. Packages the code in both a tar archive and a wheel using
   [`uv build`](https://docs.astral.sh/uv/guides/publish/#building-your-package)
3. Uploads to PyPI using
   [`uv publish`](https://docs.astral.sh/uv/guides/publish/#publishing-your-package).
4. Push to GitHub. You'll need to make a release going with the commit where the
   version was bumped.
5. Bump the version to the next patch. If you made big changes and want to bump
   the version by minor, you can use `tox -e bumpversion -- minor` after.

#### Releasing on GitHub

1. Navigate to https://github.com/cthoyt/chembl-downloader/releases/new to draft
   a new release
2. Click the "Choose a Tag" dropdown and select the tag corresponding to the
   release you just made
3. Click the "Generate Release Notes" button to get a quick outline of recent
   changes. Modify the title and description as you see fit
4. Click the big green "Publish Release" button

This will trigger Zenodo to assign a DOI to your release as well.

</details>
