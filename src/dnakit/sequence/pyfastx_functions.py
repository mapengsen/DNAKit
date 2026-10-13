"""pyfastx functions with no new DNAKit wrapper.

The inventory's pyfastx functions are version and gzip helpers, plus
``pyfastx.reverse_complement``, which duplicates
``dnakit.sequence.reverse_complement_dna``. FASTQ quality detection through
pyfastx is already exposed as ``detect_fastq_quality_encoding`` in
``dnakit.sequence.backends``.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: tuple[str, ...] = ()

_register_public_verbose(globals(), tuple(__all__))
