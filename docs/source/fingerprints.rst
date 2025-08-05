Fingerprints
============

ChEMBL makes a file containing pre-computed 2048 bit radius 2 Morgan fingerprints for
each molecule available. It can be downloaded using:

.. code-block:: python

    import chembl_downloader

    path = chembl_downloader.download_fps()

The `version` and other keyword arguments are also valid for this function.

Load fingerprints with :mod:`chemfp`
------------------------------------

The following wraps the :func:`chembl_downloader.download_fps` function with
:mod:`chemfp`'s fingerprint loader:

.. code-block:: python

    import chembl_downloader

    arena = chembl_downloader.chemfp_load_fps()

The ``version`` and other keyword arguments are also valid for this function. More
information on working with the :class:`chemfp.arena.FingerprintArena` object can be
found `here
<https://chemfp.readthedocs.io/en/latest/using-api.html#working-with-a-fingerprintarena>`_.
