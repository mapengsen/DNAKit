"""Review group G18c for pyfaidx: no new wrappers in this batch.

The two reviewed ledger rows are not DNA-level operations. ``pyfaidx.getmtime_fsspec`` reads
file modification times for local or fsspec paths, and ``pyfaidx.get_valid_filename``
sanitizes file names. Neither reads or writes sequence data. The sequence helpers of pyfaidx
are wrapped in ``dnakit.sequence.pyfaidx_functions``.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
