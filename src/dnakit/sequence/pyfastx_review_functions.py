"""pyfastx functions reviewed by group G18o: no new wrapper in this batch.

The one remaining pyfastx entry, ``pyfastx.version``, returns the package version string. The
wrappers already record the version in their provenance, so it is internal.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: tuple[str, ...] = ()

_register_public_verbose(globals(), tuple(__all__))
