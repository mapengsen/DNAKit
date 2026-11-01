"""screed indexed-record access, alphabet checks and FASTX output (P41)."""
from __future__ import annotations
from typing import Any
from dnakit._verbosity import register_public as _register_public_verbose
from collections.abc import Mapping
from itertools import islice
from dnakit.bioinformatics._common import require
from dnakit.bioinformatics import _p40_helpers as h


def screed_record_at_index(database: str, index: int, *, verbose: bool = True) -> dict[str, Any]:
    """Retrieve one screed database record at its zero-based index."""
    module = require("screed", "bio-screed")
    source = h.existing(database, "database")
    db = module.ScreedDB(str(source))
    try:
        raw = db.loadRecordByIndex(h.integer(index, "index", minimum=0))
        record = {key: value if isinstance(value, (str, int, float, bool, type(None))) else str(value) for key, value in raw.items()}
    finally:
        db.close()
    return {"status": "completed", "backend": "screed", "record": record}


def screed_record_names(database: str, *, max_records: int = 100000, verbose: bool = True) -> dict[str, Any]:
    """List record keys of a screed database in upstream iteration order."""
    module = require("screed", "bio-screed")
    db = module.ScreedDB(str(h.existing(database, "database")))
    try:
        limit = h.integer(max_records, "max_records", minimum=1, maximum=100000)
        keys = list(islice(db.keys(), limit + 1))
    finally:
        db.close()
    limit = h.integer(max_records, "max_records", minimum=1)
    return {"status": "completed", "backend": "screed", "names": keys[:limit], "truncated": len(keys) > limit}


def screed_is_dna(sequence: str, *, verbose: bool = True) -> dict[str, Any]:
    """Check screed's exact uppercase ACGTN alphabet, without silently uppercasing input."""
    module = require("screed.dna", "bio-screed")
    return {"status": "completed", "backend": "screed", "is_dna": bool(module.is_DNA(sequence))}


def screed_write_record(record: Mapping[str, Any], output: str, *, file_format: str = "fasta", verbose: bool = True) -> dict[str, Any]:
    """Write a screed record as FASTA or FASTQ to a new local file."""
    from dnakit.sequence.screed_functions import format_fastx_records
    target = h.new_file(output, "output")
    result = format_fastx_records([dict(record)], output_format=file_format, verbose=False)
    # The existing formatter uses screed's FASTX writer and preserves its defaults.
    target.write_text(result["text"])
    return {"status": "completed", "backend": "screed", "outputs": [str(target)]}


__all__ = ["screed_record_at_index", "screed_record_names", "screed_is_dna", "screed_write_record"]

_register_public_verbose(globals(), tuple(__all__))
