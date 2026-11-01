"""pyfaidx index and in-region editing operations as DNAKit functions (pack P40).

The functions call the pyfaidx Python package in process. ``.fai`` index reading and writing work with pyfaidx 0.8;
the BGZF ``.gzi`` index operations need pyfaidx 0.9 or later (and Biopython for the BGZF block reader) and report
``BackendUnavailableError`` with an upgrade hint on older versions. Index files are always written to new paths and
the input FASTA is never modified: ``faidx_overwrite_region`` edits a copy.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_EXTRA = "bio-pyfaidx"
_KINDS = ("fai", "gzi")
_GZI_HINT = 'Upgrade pyfaidx: python -m pip install --upgrade "pyfaidx>=0.9" biopython'


def _pyfaidx() -> Any:
    return require("pyfaidx", _EXTRA)


def _faidx(module: Any, fasta: Path, index: Path | None, *, gzi: Path | None = None, build: bool,
           mutable: bool = False) -> Any:
    kwargs: dict[str, Any] = {"build_index": build, "rebuild": False, "mutable": mutable}
    if index is not None:
        kwargs["indexname"] = str(index)
    if gzi is not None:
        kwargs["gzi_indexname"] = str(gzi)
    try:
        return module.Faidx(str(fasta), **kwargs)
    except TypeError as exc:
        raise BackendUnavailableError(f"This pyfaidx version does not accept {sorted(kwargs)}: {exc}",
                                      hint=_GZI_HINT) from exc
    except Exception as exc:  # noqa: BLE001 - pyfaidx raises several FastaIndexingError-like classes
        raise ConfigurationError(f"pyfaidx could not open {fasta}: {exc}") from exc


def _need_gzi(faidx: Any, module: Any) -> None:
    if not hasattr(faidx, "build_gzi"):
        raise BackendUnavailableError(f"pyfaidx {module.__version__} has no BGZF .gzi support.", hint=_GZI_HINT)


def faidx_read_index(
    fasta: str | Path,
    *,
    index_path: str | Path | None = None,
    kind: str = "fai",
    max_rows: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read an existing FASTA ``.fai`` index (pyfaidx Faidx.read_fai) or BGZF ``.gzi`` index (Faidx.read_gzi).

    ``kind="fai"`` returns one record per sequence (name, length, byte offset, bases per line, bytes per
    line); ``kind="gzi"`` (pyfaidx 0.9+, bgzip-compressed FASTA) returns the compressed/uncompressed block
    offsets. ``index_path`` defaults to ``<fasta>.fai`` or ``<fasta>.gzi``; the index must already exist.
    """

    module = _pyfaidx()
    which = _h.choice(kind, _KINDS, "kind")
    limit = _h.integer(max_rows, "max_rows", minimum=1, maximum=_h.MAX_ROWS * 10)
    source = _h.existing(fasta, "fasta")
    default = Path(str(source) + "." + which)
    index = _h.existing(index_path if index_path is not None else default, "index_path")
    if which == "fai":
        faidx = _faidx(module, source, index, build=False)
        # pyfaidx reads the index when the object is created (build_index=False).
        rows = [{"name": str(name), "length": int(r.rlen), "offset": int(r.offset),
                 "line_bases": int(r.lenc), "line_bytes": int(r.lenb)} for name, r in faidx.index.items()]
        faidx.close()
        return {"status": "completed", "backend": "pyfaidx", "version": module.__version__, "kind": which,
                "index_path": str(index), "record_count": len(rows), "records": rows[:limit],
                "truncated": len(rows) > limit}
    with tempfile.TemporaryDirectory() as scratch:
        faidx = _faidx(module, source, Path(scratch) / "tmp.fai", gzi=index, build=True)
        _need_gzi(faidx, module)
        faidx.read_gzi()
        blocks = [{"compressed_start": int(b.cstart), "uncompressed_start": int(b.ustart)}
                  for b in faidx.gzi_index[1:]]
        faidx.close()
    return {"status": "completed", "backend": "pyfaidx", "version": module.__version__, "kind": which,
            "index_path": str(index), "block_count": len(blocks), "blocks": blocks[:limit],
            "truncated": len(blocks) > limit}


def faidx_write_index(
    fasta: str | Path,
    output_path: str | Path,
    *,
    kind: str = "fai",
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a FASTA index and write it to a new file.

    ``kind="fai"`` builds the samtools-compatible ``.fai`` (pyfaidx Faidx.build_index + write_fai);
    ``kind="gzi"`` builds the BGZF block index of a bgzip-compressed FASTA (Faidx.build_gzi +
    write_gzi, pyfaidx 0.9+). The FASTA is not changed.
    """

    module = _pyfaidx()
    which = _h.choice(kind, _KINDS, "kind")
    source = _h.existing(fasta, "fasta")
    target = _h.new_file(output_path, "output_path")
    if which == "fai":
        faidx = _faidx(module, source, target, build=True)
        faidx.write_fai()
        count = len(faidx.index)
        faidx.close()
        return {"status": "completed", "backend": "pyfaidx", "version": module.__version__, "kind": which,
                "outputs": [str(target)], "record_count": count}
    with tempfile.TemporaryDirectory() as scratch:
        faidx = _faidx(module, source, Path(scratch) / "tmp.fai", gzi=target, build=True)
        _need_gzi(faidx, module)
        try:
            faidx.build_gzi()
        except OSError as exc:
            raise ConfigurationError(str(exc)) from exc
        faidx.write_gzi()
        count = len(faidx.gzi_index) - 1
        faidx.close()
    return {"status": "completed", "backend": "pyfaidx", "version": module.__version__, "kind": which,
            "outputs": [str(target)], "block_count": count}


def faidx_overwrite_region(
    fasta: str | Path,
    output_fasta: str | Path,
    record: str,
    start: int,
    end: int,
    sequence: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Overwrite a FASTA region with a sequence of the same length (pyfaidx Faidx.to_file).

    The FASTA is copied to the new ``output_fasta`` and the region ``start``-``end`` (1-based, inclusive, as in
    pyfaidx and ``faidx``) of ``record`` is replaced in the copy; the original is never modified. The
    replacement must have exactly ``end - start + 1`` bases; line layout and the other records stay
    byte-identical. The ``.fai`` of the copy is written next to it.
    """

    module = _pyfaidx()
    source = _h.existing(fasta, "fasta")
    target = _h.new_file(output_fasta, "output_fasta")
    name = _h.text(record, "record")
    first = _h.integer(start, "start", minimum=1)
    last = _h.integer(end, "end", minimum=first)
    replacement = _h.text(sequence, "sequence")
    if len(replacement) != last - first + 1:
        raise ConfigurationError("sequence must have exactly end - start + 1 bases.")
    shutil.copyfile(source, target)
    try:
        faidx = _faidx(module, target, None, build=True, mutable=True)
        if name not in faidx.index:
            raise ConfigurationError(f"record {name!r} is not in the FASTA file.")
        if last > faidx.index[name].rlen:
            raise ConfigurationError(f"end {last} is beyond the record length {faidx.index[name].rlen}.")
        previous = faidx.from_file(name, first, last).replace("\n", "").replace("\r", "")
        faidx.to_file(name, first, last, replacement)
        faidx.file.flush()
        faidx.close()
    except Exception:
        target.unlink(missing_ok=True)
        Path(str(target) + ".fai").unlink(missing_ok=True)
        raise
    return {"status": "completed", "backend": "pyfaidx", "version": module.__version__, "record": name,
            "start": first, "end": last, "previous_sequence": previous, "new_sequence": replacement,
            "outputs": [str(target)]}


__all__ = ["faidx_overwrite_region", "faidx_read_index", "faidx_write_index"]

_register_public_verbose(globals(), tuple(__all__))
