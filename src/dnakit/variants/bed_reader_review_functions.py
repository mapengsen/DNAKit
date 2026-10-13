"""Review wrappers for bed-reader's remaining rows (G18h).

Both rows left as ``other`` are internal helpers with no genotype operation behind them:
``bed_reader.get_num_threads`` reports the kernel thread setting and ``bed_reader.tmp_path``
returns a temporary directory for its tests. No public DNA-level function is wrapped here.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
