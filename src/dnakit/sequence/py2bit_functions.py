"""py2bit functions with no new DNAKit wrapper.

The inventory lists only ``py2bit.open``. It is already used by
``summarize_twobit`` and ``fetch_twobit_region`` in ``dnakit.sequence.backends``,
so this module adds no public wrapper.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: tuple[str, ...] = ()

_register_public_verbose(globals(), tuple(__all__))
