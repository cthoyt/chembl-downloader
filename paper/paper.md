---
title:
  Improving reproducibility of cheminformatics workflows with
  `chembl-downloader`

tags:
  - Python
  - chemistry
  - computational chemistry
  - cheminformatics
  - pharmacology
  - toxicology
  - systems pharmacology
  - systems toxicology
  - systems medicine
  - systems biology

authors:
  - name: Charles Tapley Hoyt
    orcid: 0000-0003-4423-4370
    affiliation: 1
    email: cthoyt@gmail.com
    corresponding: true

affiliations:
  - name: RWTH Aachen University, Institute of Inorganic Chemistry
    index: 1
    ror: 04xfq0f34

date: 7 July 2025
bibliography: paper.bib
repository: cthoyt/chembl-downloader
---

# Statement of need

Many modern cheminformatics workflows derive datasets from ChEMBL [@Gaulton2017;
@Zdrazil2023], but few of these datasets are published with accompanying code
for their generation. Consequently, their methodologies (e.g., selection,
filtering, aggregation) are opaque, reproduction is difficult, and
interpretation of results therefore lacks important context. Further, such
static datasets quickly become out-of-date. For example, the current version of
ChEMBL is v35 (as of December 2024), but ExCAPE-DB [@Sun2017] uses v20, Deep
Confidence [@Cortes-Ciriano2019] uses v23, the consensus dataset from
@Isigkeit2022 uses v28, and Papyrus [@Béquignon2023] uses v30. Therefore, there
is a need for tools that provide reproducible bulk access to the latest (or a
given) version of ChEMBL in order to enable researchers to make their derived
datasets more transparent, updatable, and trustworthy.

# State of the field

ChEMBL is typically accessed through its
[application programming interface (API)](https://www.ebi.ac.uk/chembl/api/data/docs),
through its Python client [@Davies2015], through its RDF platform [@Jupp2014],
or in bulk through its
[file transfer protocol (FTP) server](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases).
However, APIs and their respective wrapper libraries are generally not efficient
for querying and processing data in bulk due to dependency on network
connection, remote server uptime and load, rate limits, and the need to paginate
over results. Alternatively, bulk access is currently cumbersome for most
potential users due to the need to download, set up, and connect to databases
locally. Finally, third-party software (e.g., Pipeline Pilot, KNIME) that
provide access to ChEMBL are often inflexible or inextensible due to being
proprietary, closed source, or lacking approachable documentation.

# Summary

This article introduces `chembl-downloader`, a Python package for the
reproducible acquisition, access, and manipulation of ChEMBL data through its
FTP server.

At a low-level, it uses a combination of the
[`pystow`](https://github.com/cthoyt/pystow) Python software package and custom
logic for the reproducible acquisition and pre-processing (e.g., uncompressing)
of either the latest or a given version of most resources in the ChEMBL FTP
server. These include relational database dumps (e.g., in SQLite), molecule
lists (e.g., in SDF, TSV), pre-computed molecular fingerprints (e.g., in
binary), a monomer library (e.g., in XML), and UniProt target mappings (e.g., in
TSV).

At a mid-level, it provides utilities for accessing these files through useful
data structures and functions such as querying the SQLite database with
combination of Python's
[`sqlite3`](https://docs.python.org/3/library/sqlite3.html) library and `pandas`
[@Mckinney2010], parsing SDF files with RDKit [@rdkit], parsing TSVs with
`pandas`, loading fingerprints with `chemfp` [@Dalke2019], and parsing the
monomer library with Python's `xml` library. The low- and mid-level utilities
are kept small and simple such that they can be arbitrarily extended by users.

At a high-level, it maintains a small number of task-specific utilities. It
contains several pre-formatted SQL queries to retrieve the bioactivities
associated with a given assay or target, to retrieve the compounds mentioned in
a publication or patent, to retrieve the names of all compounds, etc.

# Case studies

A first case study demonstrates how the high-level utilities in
`chembl-downloader` can be used to reproduce the dataset generation from
@Cortes-Ciriano2019, highlight some of the methodological controversies (e.g.,
using arithmetic mean instead of geometric mean for pIC~50~ values), show the
immense variability introduced by new datapoints in later versions of ChEMBL,
and highlight the number of additional compounds added to each of its 24
target-specific datasets. @Landrum2024 further investigated the impact of such
aggregation. See the corresponding
[Jupyter notebook](https://github.com/cthoyt/chembl-downloader/blob/v0.5.2/notebooks/cortes-ciriano-refresh.ipynb).

A second case study demonstrates the value of having a reproducible script for
identifying missing identifier mappings between molecules in ChEMBL and ChEBI
[@Hastings2016] _via_ lexical mappings produced by Gilda [@Gyori2022], which
identified 4,266 potential mappings for curation e.g., in a workflow such as
Biomappings [@Hoyt2022]. See the corresponding
[Jupyter notebook](https://github.com/cthoyt/chembl-downloader/blob/v0.5.2/notebooks/chebi-mappings.ipynb).

A final case study demonstrates the utility of `chembl-downloader` by making
pull requests to the code repositories corresponding to three popular
cheminformatics blogs
([the RDKit Blog](https://greglandrum.github.io/rdkit-blog/),
[Practical Cheminformatics](https://practicalcheminformatics.blogspot.com), and
[Is Life Worth Living?](https://iwatobipen.wordpress.com/)) to make the code
more reproducible (where the source data was not available) and
ChEMBL-version-agnostic:

- [greglandrum/rdkit_blog (#5)](https://github.com/greglandrum/rdkit_blog/pull/5)
  for generating a substructure library
- [PatWalters/sfi (#11)](https://github.com/PatWalters/sfi/pull/11) for
  calculating compounds' solubility forecast index, originally proposed by
  @Hill2010
- [PatWalters/jcamd_model_comparison (#1)](https://github.com/PatWalters/jcamd_model_comparison/pull/1)
  for comparing classification models
- [iwatobipen/playground (#4)](https://github.com/iwatobipen/playground/pull/4)
  for parsing and using ChEMBL's monomer library
- [iwatobipen/playground (#5)](https://github.com/iwatobipen/playground/pull/5)
  for analyzing chemical space of molecules from a set of patents

Additional external use cases can be found by
[searching GitHub](https://github.com/search?q=%22import%20chembl_downloader%22%20OR%20%22from%20chembl_downloader%20import%22%20language%3APython%20NOT%20is%3Afork%20-owner%3Acthoyt&type=code).
Finally, several scholarly articles, including from the ChEMBL group itself,
have used `chembl-downloader` in their associated code [@Domingo-Fernández2023;
@Gadiya2023; @Nisonoff2023; @Zdrazil2023; @Gorostiola2024; @Zhang2024;
@Schoenmaker2025].

# Availability and usage

`chembl-downloader` is available as a package on
[PyPI](https://pypi.org/project/chembl-downloader) with the source code
available at
[https://github.com/cthoyt/chembl-downloader](https://github.com/cthoyt/chembl-downloader)
and documentation available at
[https://chembl-downloader.readthedocs.io](https://chembl-downloader.readthedocs.io).
The repository also contains an interactive Jupyter notebook tutorial and
notebooks for the case studies described above.

# Acknowledgements

The author would like to thank Yojana Gadiya and Jennifer HY Lin for helpful
discussions and the NFDI4Chem Consortium (https://www.nfdi4chem.de) for support.

# References
