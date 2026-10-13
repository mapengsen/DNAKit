"""Review wrappers for bgen-reader's remaining rows (G18h).

The one row left as ``other`` is the package test runner, which is not a genotype operation.
No public function is wrapped here.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
