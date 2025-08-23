SQL Usage
=========

Download an extract the SQLite dump using the following:

.. code-block:: python

    import chembl_downloader

    path = chembl_downloader.download_extract_sqlite(version="28")

After it's been downloaded and extracted once, it's smart and does not need to download
again. It gets stored using :mod:`pystow` automatically in the ``~/.data/chembl``
directory.

You can modify the previous code slightly by omitting the ``version`` keyword argument
to automatically find the latest version of ChEMBL:

.. code-block:: python

    import chembl_downloader

    path = chembl_downloader.download_extract_sqlite()

The ``version`` keyword argument is available for all functions in this package (e.g.,
including :func:`chembl_downloader.connect`, :func:`chembl_downloader.cursor`, and
:func:`chembl_downloader.query`), but will be omitted below for brevity.

Automatic Connection
--------------------

Inside the archive is a single SQLite database file. Normally, people manually untar
this folder then do something with the resulting file. Don't do this, it's not
reproducible! Instead, the file can be downloaded and a connection can be opened
automatically with:

.. code-block:: python

    import chembl_downloader

    with chembl_downloader.connect() as conn:
        with conn.cursor() as cursor:
            cursor.execute(...)  # run your query string
            rows = cursor.fetchall()  # get your results

The :func:`chembl_downloader.cursor` function provides a convenient wrapper around this
operation:

.. code-block:: python

    import chembl_downloader

    with chembl_downloader.cursor() as cursor:
        cursor.execute(...)  # run your query string
        rows = cursor.fetchall()  # get your results

Run a Query and Get a Pandas DataFrame
--------------------------------------

The most powerful function is :func:`chembl_downloader.query`, which builds on the
previous :func:`chembl_downloader.connect` function in combination with
:func:`pandas.read_sql` to make a query and load the results into a pandas DataFrame for
any downstream use.

.. code-block:: python

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
    df.to_csv(..., sep="\t", index=False)

.. note::

    Suggestion 1: use :mod:`pystow` to make a reproducible file path that's portable to
    other people's machines (e.g., it doesn't have your username in the path).

    Suggestion 2: RDKit is now pip-installable with `pip install rdkit-pypi`, which
    means most users don't have to muck around with complicated conda environments and
    configurations. One of the powerful but understated tools in RDKit is the
    [rdkit.Chem.PandasTools](https://rdkit.org/docs/source/rdkit.Chem.PandasTools.html)
    module.

Querying for a Scalar
---------------------

For SQL queries that return a scalar result, you can use
:func:`chembl_downloader.query_scalar`. In the following example, this gets a summary
count over the number of activities in the database.

.. code-block:: python

    import chembl_downloader

    sql = "SELECT COUNT(activity_id) FROM activities"
    count: int = chembl_downloader.query_scalar(sql)
