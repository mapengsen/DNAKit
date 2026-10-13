"""Review of the gget functions left unwrapped after the G18j pass.

No gget function in this pass is both DNA-level and runnable offline. ``cbio_search`` and ``gpt``
are not genomics operations, ``cellxgene`` needs the network and the ``cellxgene_census`` package
(not installed), and ``setup`` is an installer. The ledger records each reason, so this module
exports nothing.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
