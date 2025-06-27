# ChEMBL Downloader Paper

This folder contains source material for a
[Journal of Open Source Software](https://joss.theoj.org/) submission:

1. [paper.md](paper.md): the main manuscript
2. [paper.bib](paper.bib): BibTex citations for the manuscript

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
