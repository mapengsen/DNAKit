"""bgen_reader example data: fetch a sample BGEN file into the bgen_reader cache.

``fetch_bgen_example`` returns the local path of one of the BGEN sample files that bgen_reader publishes
(genotype files with their sample tables). The first call downloads the file from GitHub (network access
needed); later calls use the cache and check the file hash. Unknown file names are refused before any download.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-bgen-reader"
_SOURCE = "https://raw.githubusercontent.com/fastlmm/bgen-sample-files/main"


def fetch_bgen_example(filename: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the local path of a bgen_reader sample file, downloading it on first use (bgen_reader.example_filepath).

    ``filename`` must be one of the names bgen_reader knows (for example ``example.bgen`` or ``complex.bgen``).
    The download needs network access, the ``requests`` package, and a successful hash check.
    """
    module = require("bgen_reader._example", _EXTRA)
    require("requests", _EXTRA)
    if not isinstance(filename, str) or filename not in module._filenames:
        known = ", ".join(sorted(module._filenames))
        raise ConfigurationError(
            f"Unknown bgen_reader example file {filename!r}. Known names: {known}."
        )
    path = Path(module.example_filepath(filename))
    return {
        "status": "completed",
        "backend": "bgen_reader",
        "name": filename,
        "path": str(path),
        "size_bytes": int(path.stat().st_size),
        "source": f"{_SOURCE}/{filename}",
    }


__all__ = ["fetch_bgen_example"]
_register_public_verbose(globals())
