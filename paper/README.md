# ChEMBL Downloader Paper

This folder contains source material for a
[Journal of Open Source Software](https://joss.theoj.org/) submission:

1. [paper.md](paper.md): the main manuscript
2. [paper.bib](paper.bib): BibTex citations for the manuscript

## Abstract

_adapted from the statement of need section_

Many modern cheminformatics workflows derive datasets from ChEMBL, but few of
these datasets are published with accompanying code for their generation.
Consequently, their methodologies (e.g., selection, filtering, aggregation) are
opaque, reproduction is difficult, and interpretation of results therefore lacks
important context. Further, such static datasets quickly become out-of-date. For
example, the current version of ChEMBL is v35 (as of December 2024), but
ExCAPE-DB uses v20, Deep Confidence uses v23, the consensus dataset from
Isigkeit _et al._ (2022) uses v28, and Papyrus uses v30. Therefore, there is a
need for tools that provide reproducible bulk access to the latest (or a given)
version of ChEMBL in order to enable researchers to make their derived datasets
more transparent, updatable, and trustworthy. This article introduces
`chembl-downloader`, a Python package for the reproducible acquisition, access,
and manipulation of ChEMBL data through its FTP server. It can be downloaded
under the MIT license from https://github.com/cthoyt/chembl-downloader and
installed from PyPI with `pip install chembl-downloader.`

## Build

Follow the instructions at
https://joss.readthedocs.io/en/latest/submitting.html#docker:

```console
$ sh build.sh
```

## Linting

```console
$ npx prettier --prose-wrap always --check "**/*.md" --write
```

## Submission

Follow the instructions at https://joss.theoj.org/papers/new.

## License

The source files and outputs for the manuscript are licensed under CC-BY-4.0.
