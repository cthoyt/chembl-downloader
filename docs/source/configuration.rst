Configuration
=============

If you want to store the data elsewhere using :mod:`pystow`, you can use the ``prefix``
argument.

.. code-block:: python

    import chembl_downloader

    # It gets downloaded/extracted to
    # ~/.data/pyobo/raw/chembl/29/chembl_29/chembl_29_sqlite/chembl_29.db
    path = chembl_downloader.download_extract_sqlite(prefix=["pyobo", "raw", "chembl"])

See the `pystow documentation
<https://github.com/cthoyt/pystow#%EF%B8%8F-configuration>`_ on configuring the storage
location further.

The `prefix` keyword argument is available for all functions in this package (e.g.,
including `connect()`, `cursor()`, and `query()`).
