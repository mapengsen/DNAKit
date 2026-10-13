"""Review module for mygene (group G18f): no further DNA-level wrappers.

Every ``other`` row of the mygene ledger is a client class or factory, a cache or endpoint
setting, or a service metadata listing (``get_fields``, ``metadata``). None of them returns gene
or sequence records, so this module exports no functions. The gene lookups are wrapped by
``gene_backends`` (``query_mygene_genes``, ``map_gene_identifiers``, ``fetch_mygene_gene``) and
``mygene_functions`` (``mygene_getgenes``).
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
