"""Mappy (minimap2) helpers that this package does not wrap again.

The three module-level functions in the mappy inventory duplicate existing DNAKit functions or are
logging controls, so no wrapper is added here:

- ``mappy.revcomp`` duplicates ``dnakit.sequence.reverse_complement_dna``.
- ``mappy.fastx_read`` duplicates ``dnakit.bioinformatics.reads.read_fastx``.
- ``mappy.verbose`` only sets the minimap2 log level.

Mapping, indexing and region fetches through ``mappy.Aligner`` are already wrapped in backends.py as
``map_query_sequences_paf()``, ``build_minimap2_index()`` and ``fetch_minimap2_reference_region()``.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: tuple[str, ...] = ()

_register_public_verbose(globals(), tuple(__all__))
