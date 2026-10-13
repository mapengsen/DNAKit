"""Module placeholder for pybiomart: the package exposes no module-level public functions.

``pybiomart`` provides the ``Server``, ``Mart`` and ``Dataset`` classes. The BioMart queries
built on them (``query_biomart_dataset``, ``list_biomart_datasets``, ``list_biomart_attributes``)
are already wrapped in ``dnakit.databases.gene_backends``. Importing this module does not import
pybiomart, so it does not create the ``.pybiomart.sqlite`` request cache.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: tuple[str, ...] = ()

_register_public_verbose(globals(), tuple(__all__))
