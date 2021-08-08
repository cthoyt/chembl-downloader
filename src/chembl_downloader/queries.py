# -*- coding: utf-8 -*-

"""A collection of query strings for ChEBML."""

from textwrap import dedent

#: This query yields tuples of ChEMBL identifiers and their preferred names, omitting
#: pairs where there is no name or if there's no structure
ID_NAME_QUERY = """
SELECT
    MOLECULE_DICTIONARY.chembl_id,
    MOLECULE_DICTIONARY.pref_name
FROM MOLECULE_DICTIONARY
JOIN COMPOUND_STRUCTURES ON MOLECULE_DICTIONARY.molregno == COMPOUND_STRUCTURES.molregno
WHERE molecule_dictionary.pref_name IS NOT NULL
"""

#: This query returns five examples from the molecule dictionary table
ID_NAME_QUERY_EXAMPLE = ID_NAME_QUERY + "\nLIMIT 5"

ACTIVITIES_QUERY = """
SELECT
    COMPOUND_STRUCTURES.canonical_smiles,
    MOLECULE_DICTIONARY.chembl_id,
    ACTIVITIES.BAO_ENDPOINT,
    ACTIVITIES.STANDARD_RELATION,
    ACTIVITIES.STANDARD_VALUE,
    ACTIVITIES.STANDARD_UNITS,
    ASSAYS.ASSAY_TYPE,
    TARGET_DICTIONARY.organism,
    ASSAYS.TID,
    TARGET_DICTIONARY.chembl_id as target_chembl_id,
    TARGET_DICTIONARY.pref_name
FROM MOLECULE_DICTIONARY
JOIN COMPOUND_STRUCTURES ON MOLECULE_DICTIONARY.molregno == COMPOUND_STRUCTURES.molregno
JOIN ACTIVITIES ON MOLECULE_DICTIONARY.molregno == ACTIVITIES.molregno
JOIN ASSAYS ON ACTIVITIES.ASSAY_ID == ASSAYS.ASSAY_ID
JOIN TARGET_DICTIONARY ON ASSAYS.TID == TARGET_DICTIONARY.TID
WHERE
    assay_type in ('B','F')
    and ACTIVITIES.standard_value is not null
    and ACTIVITIES.standard_units = 'nM'
    and ACTIVITIES.standard_relation is not null
    and ACTIVITIES.standard_type = 'IC50'
    and ACTIVITIES.standard_relation = '='
"""


def get_assay_sql(assay_chembl_id: str) -> str:
    """Get the SQL for the given assay."""
    return dedent(
        f"""\
        SELECT
            COMPOUND_STRUCTURES.canonical_smiles,
            MOLECULE_DICTIONARY.chembl_id,
            ACTIVITIES.STANDARD_TYPE,
            ACTIVITIES.STANDARD_RELATION,
            ACTIVITIES.STANDARD_VALUE,
            ACTIVITIES.STANDARD_UNITS
        FROM MOLECULE_DICTIONARY
        JOIN COMPOUND_STRUCTURES ON MOLECULE_DICTIONARY.molregno == COMPOUND_STRUCTURES.molregno
        JOIN ACTIVITIES ON MOLECULE_DICTIONARY.molregno == ACTIVITIES.molregno
        JOIN ASSAYS ON ACTIVITIES.ASSAY_ID == ASSAYS.ASSAY_ID
        WHERE
            ASSAYS.chembl_id = '{assay_chembl_id}'
            and ACTIVITIES.standard_value is not null
            and ACTIVITIES.standard_relation is not null
            and ACTIVITIES.standard_relation = '='
    """  # noqa: S608
    )
