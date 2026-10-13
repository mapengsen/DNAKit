"""Review group G18q for selene-sdk: no new DNA-level wrappers.

The 27 ``other`` ledger rows are logging and configuration helpers, run-directory
bookkeeping, prediction writers, the trained-model ``predict`` call, and two
``get_reverse_encoding`` aliases. ``get_reverse_encoding`` raises IndexError on one-hot
rows upstream; the operation it is meant to perform is ``reverse_complement_one_hot``.
Model inference needs trained weights and is not wrapped. The statuses are recorded in
the selene-sdk ledger, and this module exports no functions.
"""

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
