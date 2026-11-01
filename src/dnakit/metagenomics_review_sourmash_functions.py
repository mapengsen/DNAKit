"""sourmash operations from the coverage review: sketch many genomes from a CSV list, grep, inflate and
migrate signatures, combine and convert SBT indexes.

The user installs sourmash (``python -m pip install sourmash``). Each function runs one sourmash subcommand with an
explicit typed signature mapped to fixed flags; no shell is used; every output is a new explicit path that must not
exist yet. Results are JSON-serialisable dictionaries with status, backend, command, outputs and a bounded log tail.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path
from typing import Any, Sequence

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_HINT = "Install sourmash 4.9 or later (python -m pip install sourmash, or conda install -c bioconda sourmash)."
_MAX_TAIL = 4000
_MOLTYPE_FLAGS = {"DNA": "--dna", "RNA": "--rna", "protein": "--protein", "dayhoff": "--dayhoff", "hp": "--hp"}
_BACKENDS = ("zip", "fs")


def _program(sourmash_path: str | None) -> str:
    if sourmash_path is not None:
        candidate = Path(sourmash_path).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"sourmash_path does not exist: {candidate}", hint=_HINT)
        return str(candidate)
    found = shutil.which("sourmash")
    if found is None:
        # sourmash installed in the same environment as the running interpreter
        beside = Path(sys.executable).with_name("sourmash")
        if beside.is_file():
            return str(beside)
        raise BackendUnavailableError("sourmash was not found on PATH or beside the Python interpreter.", hint=_HINT)
    return found


def _run(
    args: list[str], outputs: Sequence[Path], sourmash_path: str | None, verbose: bool
) -> dict[str, Any]:
    result = execute_bounded_command(
        _program(sourmash_path), args, backend_id="sourmash", cwd=Path.cwd(), verbose=verbose
    )
    status = "completed" if result.return_code == 0 and all(p.exists() for p in outputs) else "failed"
    return {
        "status": status,
        "backend": "sourmash",
        "command": ["sourmash", *args],
        "return_code": result.return_code,
        "outputs": [str(p) for p in outputs if p.exists()],
        "log_tail": (result.output or "")[-_MAX_TAIL:],
    }


def _existing(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return resolved


def _new(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _nonempty(values: Sequence[str], name: str) -> list[str]:
    items = [str(v) for v in values]
    if not items:
        raise ConfigurationError(f"{name} needs at least one entry.")
    return items


def sketch_fromfile_signatures(
    csv_paths: Sequence[str],
    output_signatures: str,
    *,
    param_strings: Sequence[str] = ("dna,k=31,scaled=1000",),
    already_done: Sequence[str] = (),
    license: str = "CC0",
    check_sequence: bool = False,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sketch every genome listed in CSV files with name, genome_filename and protein_filename columns, and save
    the signatures to one new collection (sourmash sketch fromfile). ``param_strings`` are sourmash sketch parameter
    strings such as ``dna,k=31,scaled=1000``. Signatures already present in ``already_done`` collections are not
    recalculated."""
    csvs = [str(_existing(p, "csv_path")) for p in _nonempty(csv_paths, "csv_paths")]
    out = _new(output_signatures, "output_signatures")
    args = ["sketch", "fromfile", *csvs, "--output-signatures", str(out), "--license", license]
    for param in _nonempty(param_strings, "param_strings"):
        args += ["-p", param]
    if already_done:
        args += ["--already-done", *[str(_existing(p, "already_done")) for p in already_done]]
    if check_sequence:
        args.append("--check-sequence")
    return _run(args, [out], sourmash_path, verbose)


def grep_signatures(
    pattern: str,
    signatures: Sequence[str],
    output: str,
    *,
    ignore_case: bool = False,
    invert_match: bool = False,
    ksize: int | None = None,
    moltype: str | None = None,
    csv_manifest: str | None = None,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the signatures whose name matches a regular expression and save them to a new file (sourmash sig grep).
    ``ignore_case`` matches both cases; ``invert_match`` keeps the non-matching signatures. ``ksize`` and ``moltype``
    (one of DNA, RNA, protein, dayhoff, hp) restrict the signatures; ``csv_manifest`` also writes a CSV manifest."""
    if not pattern:
        raise ConfigurationError("pattern must not be empty.")
    sigs = [str(_existing(p, "signature")) for p in _nonempty(signatures, "signatures")]
    out = _new(output, "output")
    args = ["sig", "grep", pattern, *sigs, "-o", str(out)]
    if ignore_case:
        args.append("-i")
    if invert_match:
        args.append("-v")
    if ksize is not None:
        if ksize <= 0:
            raise ConfigurationError("ksize must be positive.")
        args += ["-k", str(ksize)]
    if moltype is not None:
        if moltype not in _MOLTYPE_FLAGS:
            raise ConfigurationError(f"moltype must be one of {sorted(_MOLTYPE_FLAGS)}.")
        args.append(_MOLTYPE_FLAGS[moltype])
    outputs = [out]
    if csv_manifest is not None:
        manifest = _new(csv_manifest, "csv_manifest")
        args += ["--csv", str(manifest)]
        outputs.append(manifest)
    return _run(args, outputs, sourmash_path, verbose)


def inflate_signatures(
    signature_from: str,
    other_signatures: Sequence[str],
    output: str,
    *,
    ksize: int | None = None,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Give other signatures the abundances of a reference signature (sourmash sig inflate). The reference must track
    abundances; each other signature is inflated with ``MinHash.inflate`` against it and saved to one new file.
    ``ksize`` selects the k-mer size; when omitted the DNA molecule type is used."""
    ref = str(_existing(signature_from, "signature_from"))
    others = [str(_existing(p, "other_signatures")) for p in _nonempty(other_signatures, "other_signatures")]
    out = _new(output, "output")
    args = ["sig", "inflate", ref, *others, "-o", str(out)]
    if ksize is not None:
        if ksize <= 0:
            raise ConfigurationError("ksize must be positive.")
        args += ["-k", str(ksize)]
    else:
        args.append("--dna")
    return _run(args, [out], sourmash_path, verbose)


def combine_sbt_indexes(
    output_sbt: str,
    input_sbts: Sequence[str],
    *,
    bf_size: float = 1e5,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge several SBT indexes into one new SBT index (sourmash sbt_combine). ``bf_size`` is the Bloom filter size
    used when the index is saved."""
    inputs = [str(_existing(p, "input_sbt")) for p in _nonempty(input_sbts, "input_sbts")]
    if bf_size <= 0:
        raise ConfigurationError("bf_size must be positive.")
    out = _new(output_sbt, "output_sbt")
    args = ["sbt_combine", str(out), *inputs, "-x", repr(float(bf_size))]
    return _run(args, [out], sourmash_path, verbose)


def migrate_sbt_index(
    sbt_path: str,
    output_dir: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rewrite an SBT index in the latest format (sourmash migrate) and save the copy in ``output_dir`` under the same
    file name. The original file is not changed. The name is kept on purpose: the storage tag of a ``.sbt.zip`` index
    is derived from its file name, so a renamed copy would lose its signatures. Self-contained ``.sbt.zip`` indexes
    are supported."""
    source = _existing(sbt_path, "sbt_path")
    folder = Path(output_dir).expanduser().resolve()
    folder.mkdir(parents=True, exist_ok=True)
    out = _new(folder / source.name, "migrated SBT index")
    shutil.copyfile(source, out)
    return _run(["migrate", str(out)], [out], sourmash_path, verbose)


def convert_sbt_storage(
    sbt_path: str,
    output_sbt: str,
    backend: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert the storage of an SBT index to another backend (sourmash storage convert). ``backend`` is ``zip(PATH)``
    or ``fs(PATH)``, where PATH is a new file or folder that receives the signature data. The input is copied to
    ``output_sbt`` first and the copy is converted."""
    source = _existing(sbt_path, "sbt_path")
    if not isinstance(backend, str) or "(" not in backend or not backend.endswith(")"):
        raise ConfigurationError("backend must look like zip(PATH) or fs(PATH).")
    name = backend.split("(", 1)[0].strip().lower()
    if name not in _BACKENDS:
        raise ConfigurationError(f"backend must start with one of {list(_BACKENDS)}.")
    target = backend.split("(", 1)[1][:-1].strip()
    if not target:
        raise ConfigurationError("backend needs a target PATH inside the parentheses.")
    out = _new(output_sbt, "output_sbt")
    shutil.copyfile(source, out)
    storage = Path(target).expanduser().resolve()
    if storage.exists():
        raise ConfigurationError(f"backend target must not exist yet: {storage}")
    return _run(["storage", "convert", str(out), "-b", f"{name}({storage})"], [out, storage], sourmash_path, verbose)


_register_public_verbose(globals())
