"""Review group G18c for screed: no new wrappers in this batch.

The three reviewed ledger rows are not new DNA operations. ``screed.hava.hava_iter`` and
``screed.seqparse.read_hava_sequences`` parse screed's HAVA example format, whose six
placeholder fields carry no DNA meaning. ``screed.screedRecord.write_fastx_pair`` raises
NameError upstream, and its pair output is covered by
``dnakit.sequence.screed_functions.format_fastx_records``, which writes the same records in
order.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
