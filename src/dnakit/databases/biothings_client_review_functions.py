"""Review of the biothings_client functions left unwrapped after the G18j pass.

No biothings_client function in this pass is DNA-level and runnable offline. ``get_client`` and
``get_async_client`` are client factories, ``cache_request`` is a request cache helper, and the
``join`` helpers combine two network query clients without any DNA-specific meaning. The ledger
records each reason, so this module exports nothing.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
