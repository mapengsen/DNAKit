"""Optional pyvolve wrappers: none beyond the existing phylogeny simulation backends.

The public pyvolve callables in the inventory are ``read_tree`` (a Newick reader that
``simulate_nucleotide_sequences()``, ``simulate_nucleotide_alignment()`` and
``simulate_substitution_counts()`` in :mod:`dnakit.phylogeny.backends` already call) and
``print_tree`` (a text printer). Simulation itself is class-based (``Model``, ``Partition``,
``Evolver``) and is wrapped by those existing functions, so this module adds no wrappers.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
