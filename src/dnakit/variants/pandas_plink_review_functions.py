"""Review wrappers for pandas-plink's remaining rows (G18h).

The three rows left as ``other`` are a bundled-data path helper, the package test runner and
the ``Chunk`` parameter class of ``read_plink1_bin``. None of them is a genotype operation, so
no public function is wrapped here.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
