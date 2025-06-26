# Extras

the Bioactivity-Explorer [@Liang2019] uses ChEMBL 24 (2018-05-01)

None of these datasets include the code used to produce the dataset. Some
include other code that consumes the dataset (e.g. [@Liang2019] links code for
web service, [@Bequignon2021] links to code for analysis). It's difficult to say
exactly what's going on in here, but people are not happy to share what they
did. Therefore, there's a need for simple, reproducible access to ChEMBL that
works in bulk.

- ExCAPE-DB [@Sun2017] uses ChEMBL 20 and max aggregation
- Deep Confidence [@Cortes-Ciriano2019] uses ChEMBL 23 and arithmetic mean
  aggregation
- Bioactivity-explorer [@Liang2019] uses ChEMBL 24 and uses
  min/max/median/arithmetic aggregation
- The consensus dataset from [@Isigkeit2022] uses ChEMBL 28 but does not give an
  explanation of how they aggregated
- Papyrus [@Bequignon2021] uses ChEMBL 29 and arithmetic mean
- Data inside various analyses contributed to RDKit (e.g.,
  [Free Wilson analysis examples](https://github.com/rdkit/rdkit/tree/master/Contrib/FreeWilson/data))

Potential CDK integration [@Willighagen2017] or Bacting

## Third case study

A third case study demonstrates using the query interface during the
construction of the INDRA CoGEx knowledge graph. Interestingly, no published
knowledge graphs have yet to include more generic bioactivity databases like
ChEMBL while many have focused on drug-specific databases such as DrugBank and
Drug Central [@Bonner2022].

> This case study won't be used since there is no reference publication for
> INDRA CoGEx
