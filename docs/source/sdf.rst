SDF Usage
=========

This example is a bit more fit-for-purpose than the last two. The
:func:`chembl_downloader.supplier` function makes sure that the latest SDF dump is
downloaded and loads it from the gzip file into a `rdkit.Chem.ForwardSDMolSupplier`
using a context manager to make sure the file doesn't get closed until after parsing is
done. Like the previous examples, it can also explicitly take a ``version``.

.. code-block:: python

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

This example was adapted from Greg Landrum's RDKit blog post on `generalized
substructure search
<https://greglandrum.github.io/rdkit-blog/tutorial/substructure/2021/08/03/generalized-substructure-search.html>`_.

Iterate over SMILES
-------------------

This example uses the :func:`chembl_downloader.supplier` method and RDKit to get
SMILES strings from molecules in ChEMBL's SDF file. If you want direct access to the
RDKit molecule objects, use :func:`chembl_downloader.supplier`.

.. code-block:: python

    import chembl_downloader

    for smiles in chembl_downloader.iterate_smiles():
        print(smiles)

Get an RDKit substructure library
---------------------------------

Building on the :func:`chembl_downloader.supplier` function, the
:func:`chembl_downloader.get_substructure_library` makes the preparation of a
`substructure library
<https://www.rdkit.org/docs/cppapi/classRDKit_1_1SubstructLibrary.html>`_ automated and
reproducible. Additionally, it caches the results of the build, which takes on the order
of tens of minutes, only has to be done once and future loading from a pickle object
takes on the order of seconds.

The implementation was inspired by Greg Landrum's RDKit blog post, `Some new features in
the Substruct Library
<https://greglandrum.github.io/rdkit-blog/tutorial/substructure/2021/12/20/substructlibrary-search-order.html>`_.
The following example shows how it can be used to accomplish some of the first tasks
presented in the post:

.. code-block:: python

    from rdkit import Chem

    import chembl_downloader

    library = chembl_downloader.get_substructure_library()
    query = Chem.MolFromSmarts("[O,N]=C-c:1:c:c:n:c:c:1")
    matches = library.GetMatches(query)
