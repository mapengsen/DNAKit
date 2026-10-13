"""Ncls interval-index entry points for DNAKit annotation.

The inventory for this package lists two index constructors and three profiling helpers:

- ``ncls.NCLS`` builds an integer interval index (``NCLS64`` or ``NCLS32``). Its DNA-level
  overlap, count, and containment queries are already exported from ``dnakit.annotation.backends``
  as ``query_ncls_overlaps``, ``count_ncls_overlaps``, and ``find_ncls_containing_intervals``.
- ``ncls.FNCLS`` builds an index over double-precision coordinates. Its intervals are not
  base-pair positions, so it is not wrapped here.
- ``ncls.src.*.profile`` are Cython profiling helpers and are not part of the public API.

No new function is needed, so this module exports nothing. It exists so that every package
in the group has a module in ``dnakit.annotation``.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
