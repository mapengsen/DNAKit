"""pyhmmer functions reviewed by group G18o: no new wrapper in this batch.

``pyhmmer.hmmpress`` writes the binary press index that ``hmmscan`` reads. ``hmmscan`` is a
protein search, and the DNA search path (``nhmmer``) reads plain HMM profiles, so the press step
produces no DNA-level result and is internal. Profile search and alignment are already wrapped in
``pyhmmer_functions.py`` and ``search_dna_hmm``.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: tuple[str, ...] = ()

_register_public_verbose(globals(), tuple(__all__))
