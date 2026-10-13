"""mappy functions reviewed by group G18o: no new wrapper in this batch.

The only other function in the inventory, ``mappy.verbose``, sets the minimap2 log level for the
whole process. It is a logging control, so it is internal and not wrapped. Mapping and indexing
remain available through the existing ``map_query_sequences_paf`` and ``build_minimap2_index``.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: tuple[str, ...] = ()

_register_public_verbose(globals(), tuple(__all__))
