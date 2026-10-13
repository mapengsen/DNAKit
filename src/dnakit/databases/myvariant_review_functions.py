"""Review module for myvariant (group G18f): no further DNA-level wrappers.

Every ``other`` row of the myvariant ledger is a client class or factory, a cache or endpoint
setting, or a service metadata listing (``get_fields``, ``metadata``). None of them returns
variant records, so this module exports no functions. The variant lookups are wrapped by
``gene_backends`` (``fetch_variant_annotations``, ``query_variant_annotations``,
``hgvs_ids_from_vcf``, ``format_hgvs_id``) and ``myvariant_functions`` (``myvariant_querymany``).
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
