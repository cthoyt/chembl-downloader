# -*- coding: utf-8 -*-

"""A collection of query strings for ChEBML."""

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
