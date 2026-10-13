"""parasail functions reviewed by group G18o: no new wrapper in this batch.

The 13 unwrapped parasail functions are configuration setters (``set_case_sensitive``,
``set_alphabet_aliases``, ``set_encoding``), CPU probes (``can_use_*``), a timer, C build-path
helpers, the version string and ``matrix_create``. The wrappers in ``backends.py`` build their own
DNA matrices from match and mismatch values, so ``matrix_create`` is internal. The global setters
would change parasail state for every later call, so they are not wrapped.

``parasail.ssw_init`` (score_size 16) segfaults in parasail 1.3.4; it is not wrapped here, and
any future wrapper must guard its score_size before the call.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: tuple[str, ...] = ()

_register_public_verbose(globals(), tuple(__all__))
