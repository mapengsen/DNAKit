"""sourmash command-line operations as typed DNAKit facades: sketch DNA and compare signatures.

The user installs sourmash (``pip install sourmash``). Each function runs one sourmash subcommand with an
explicit typed signature mapped to fixed flags; no shell; outputs are new explicit paths; results are
JSON-serialisable dictionaries with status, backend, command, outputs and a bounded log tail.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any, Sequence

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_HINT = "Install sourmash 4.9 or later (python -m pip install sourmash, or conda install -c bioconda sourmash)."
_MAX_TAIL = 4000


def _program(sourmash_path: str | None) -> str:
    if sourmash_path is not None:
        candidate = Path(sourmash_path).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"sourmash_path does not exist: {candidate}", hint=_HINT)
        return str(candidate)
    found = shutil.which("sourmash")
    if found is None:
        raise BackendUnavailableError("sourmash was not found on PATH.", hint=_HINT)
    return found


def _run(args: list[str], cwd: Path, outputs: Sequence[Path], sourmash_path: str | None) -> dict[str, Any]:
    result = execute_bounded_command(_program(sourmash_path), args, backend_id="sourmash", cwd=cwd, verbose=False)
    status = "completed" if result.return_code == 0 and all(p.is_file() for p in outputs) else "failed"
    return {
        "status": status,
        "backend": "sourmash",
        "command": ["sourmash", *args],
        "return_code": result.return_code,
        "outputs": [str(p) for p in outputs if p.is_file()],
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


def sketch_dna(
    fasta: str,
    output_signature: str,
    *,
    ksizes: Sequence[int] = (21,),
    scaled: int = 1000,
    name: str | None = None,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a DNA MinHash signature of a FASTA file with scaled sketching (sourmash sketch dna)."""
    source = _existing(fasta, "fasta")
    out = _new(output_signature, "output_signature")
    if not ksizes or any(k < 1 for k in ksizes):
        raise ConfigurationError("ksizes must be positive integers.")
    if scaled < 1:
        raise ConfigurationError("scaled must be a positive integer.")
    params = ",".join([f"k={k}" for k in ksizes] + [f"scaled={scaled}"])
    args = ["sketch", "dna", str(source), "-o", str(out), "-p", params]
    if name is not None:
        args += ["--name", name]
    return _run(args, out.parent, [out], sourmash_path)


def compare_signatures(
    signatures: Sequence[str],
    output_matrix: str,
    *,
    ksize: int = 21,
    containment: bool = False,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare signatures pairwise and write the matrix (sourmash compare). Similarity unless containment."""
    if len(signatures) < 2:
        raise ConfigurationError("compare_signatures needs at least two signatures.")
    paths = [str(_existing(s, "signature")) for s in signatures]
    out = _new(output_matrix, "output_matrix")
    args = ["compare", *paths, "-k", str(ksize), "-o", str(out)]
    if containment:
        args.append("--containment")
    return _run(args, out.parent, [out], sourmash_path)


def search_signatures(
    query: str,
    database: str,
    output_csv: str,
    *,
    ksize: int = 21,
    threshold: float = 0.08,
    containment: bool = False,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a signature database for matches to a query signature (sourmash search)."""
    q = _existing(query, "query")
    db = _existing(database, "database")
    out = _new(output_csv, "output_csv")
    if not 0 <= threshold <= 1:
        raise ConfigurationError("threshold must be between 0 and 1.")
    args = ["search", str(q), str(db), "-k", str(ksize), "--threshold", str(threshold), "-o", str(out)]
    if containment:
        args.append("--containment")
    return _run(args, out.parent, [out], sourmash_path)


def gather_signatures(
    query: str,
    database: str,
    output_csv: str,
    *,
    ksize: int = 21,
    threshold_bp: int = 50000,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Decompose a query signature into the best database matches (sourmash gather).

    ``threshold_bp`` is the minimum overlap reported (sourmash default 50 kbp).
    """
    q = _existing(query, "query")
    db = _existing(database, "database")
    out = _new(output_csv, "output_csv")
    if threshold_bp < 0:
        raise ConfigurationError("threshold_bp must not be negative.")
    args = ["gather", str(q), str(db), "-k", str(ksize), "-o", str(out), "--threshold-bp", str(threshold_bp)]
    return _run(args, out.parent, [out], sourmash_path)


def index_signatures(
    signatures: Sequence[str],
    output_database: str,
    *,
    ksize: int = 21,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a searchable signature database from signature files (sourmash index)."""
    if not signatures:
        raise ConfigurationError("signatures must contain at least one file.")
    paths = [str(_existing(s, "signature")) for s in signatures]
    out = _new(output_database, "output_database")
    # sourmash appends ".sbt.zip" to the database name it writes; report the file that actually exists.
    written = out.parent / (out.name + ".sbt.zip")
    if written.exists():
        raise ConfigurationError(f"output_database must not exist yet: {written}")
    result = _run(["index", str(out), *paths, "-k", str(ksize)], out.parent, [written], sourmash_path)
    result["database"] = str(written)
    return result


def describe_signature(
    signature: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Describe the contents of a signature file: k-mer sizes, scaled value and name (sourmash sig describe)."""
    sig = _existing(signature, "signature")
    result = execute_bounded_command(_program(sourmash_path), ["sig", "describe", str(sig)],
                                     backend_id="sourmash", cwd=sig.parent, verbose=False)
    return {
        "status": "completed" if result.return_code == 0 else "failed",
        "backend": "sourmash",
        "return_code": result.return_code,
        "description": (result.output or "")[-_MAX_TAIL:],
    }


def merge_signatures(
    signatures: Sequence[str],
    output_signature: str,
    *,
    flatten: bool = False,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge several signatures of the same k-mer size into one (sourmash sig merge)."""
    if len(signatures) < 2:
        raise ConfigurationError("merge_signatures needs at least two signatures.")
    paths = [str(_existing(s, "signature")) for s in signatures]
    out = _new(output_signature, "output_signature")
    args = ["sig", "merge", *paths, "-o", str(out)]
    if flatten:
        args.append("--flatten")
    return _run(args, out.parent, [out], sourmash_path)


def intersect_signatures(
    signatures: Sequence[str],
    output_signature: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the k-mers shared by all signatures (sourmash sig intersect)."""
    if len(signatures) < 2:
        raise ConfigurationError("intersect_signatures needs at least two signatures.")
    paths = [str(_existing(s, "signature")) for s in signatures]
    out = _new(output_signature, "output_signature")
    return _run(["sig", "intersect", *paths, "-o", str(out)], out.parent, [out], sourmash_path)


def subtract_signatures(
    signature: str,
    others: Sequence[str],
    output_signature: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove the k-mers of other signatures from one signature (sourmash sig subtract)."""
    if not others:
        raise ConfigurationError("subtract_signatures needs at least one other signature.")
    first = str(_existing(signature, "signature"))
    rest = [str(_existing(s, "signature")) for s in others]
    out = _new(output_signature, "output_signature")
    return _run(["sig", "subtract", first, *rest, "-o", str(out)], out.parent, [out], sourmash_path)


def downsample_signature(
    signature: str,
    output_signature: str,
    *,
    scaled: int,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reduce the resolution of a signature to a larger scaled value (sourmash sig downsample)."""
    if scaled < 1:
        raise ConfigurationError("scaled must be a positive integer.")
    sig = str(_existing(signature, "signature"))
    out = _new(output_signature, "output_signature")
    return _run(["sig", "downsample", sig, "--scaled", str(scaled), "-o", str(out)], out.parent, [out], sourmash_path)


def concatenate_signatures(
    signatures: Sequence[str],
    output_signature: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Concatenate signatures into one file without merging them (sourmash sig cat)."""
    if not signatures:
        raise ConfigurationError("signatures must contain at least one file.")
    paths = [str(_existing(s, "signature")) for s in signatures]
    out = _new(output_signature, "output_signature")
    return _run(["sig", "cat", *paths, "-o", str(out)], out.parent, [out], sourmash_path)


def flatten_signature(
    signature: str,
    output_signature: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove abundance information from a signature (sourmash sig flatten)."""
    sig = str(_existing(signature, "signature"))
    out = _new(output_signature, "output_signature")
    return _run(["sig", "flatten", sig, "-o", str(out)], out.parent, [out], sourmash_path)


def rename_signature(
    signature: str,
    new_name: str,
    output_signature: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Give a signature a new name (sourmash sig rename)."""
    if not new_name.strip():
        raise ConfigurationError("new_name must be a non-empty string.")
    sig = str(_existing(signature, "signature"))
    out = _new(output_signature, "output_signature")
    return _run(["sig", "rename", sig, new_name, "-o", str(out)], out.parent, [out], sourmash_path)


def signature_file_info(
    signature: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report file-level information about a signature file (sourmash sig fileinfo)."""
    sig = _existing(signature, "signature")
    result = execute_bounded_command(_program(sourmash_path), ["sig", "fileinfo", str(sig)],
                                     backend_id="sourmash", cwd=sig.parent, verbose=False)
    return {
        "status": "completed" if result.return_code == 0 else "failed",
        "backend": "sourmash",
        "return_code": result.return_code,
        "description": (result.output or "")[-_MAX_TAIL:],
    }


def filter_signature_abundance(
    signature: str,
    output_signature: str,
    *,
    min_abundance: int | None = None,
    max_abundance: int | None = None,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep hashes whose abundance lies in an inclusive range (sourmash sig filter)."""
    if min_abundance is None and max_abundance is None:
        raise ConfigurationError("filter_signature_abundance needs min_abundance or max_abundance.")
    sig = str(_existing(signature, "signature"))
    out = _new(output_signature, "output_signature")
    args = ["sig", "filter", sig, "-o", str(out)]
    if min_abundance is not None:
        args += ["-m", str(min_abundance)]
    if max_abundance is not None:
        args += ["-M", str(max_abundance)]
    return _run(args, out.parent, [out], sourmash_path)


def extract_signatures(
    collection: Sequence[str],
    output_signature: str,
    *,
    ksize: int,
    dna: bool = True,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the signatures with one k-mer size from a collection (sourmash sig extract)."""
    if not collection:
        raise ConfigurationError("collection must contain at least one signature file.")
    paths = [str(_existing(s, "signature")) for s in collection]
    out = _new(output_signature, "output_signature")
    args = ["sig", "extract", *paths, "-k", str(ksize), "-o", str(out), "--dna" if dna else "--protein"]
    return _run(args, out.parent, [out], sourmash_path)


def split_signatures(
    signatures: Sequence[str],
    output_directory: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split signature files into one file per signature with standardised names (sourmash sig split)."""
    if not signatures:
        raise ConfigurationError("signatures must contain at least one file.")
    paths = [str(_existing(s, "signature")) for s in signatures]
    out_dir = Path(output_directory).expanduser().resolve()
    if out_dir.exists() and any(out_dir.iterdir()):
        raise ConfigurationError(f"output_directory must be empty or new: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)
    result = _run(["sig", "split", *paths, "--outdir", str(out_dir)], out_dir, [], sourmash_path)
    result["outputs"] = sorted(str(p) for p in out_dir.glob("*.sig"))
    result["status"] = "completed" if result["return_code"] == 0 and result["outputs"] else "failed"
    return result


def manifest_signature(
    signature: str,
    output_csv: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a CSV manifest describing every signature in a file (sourmash sig manifest)."""
    sig = str(_existing(signature, "signature"))
    out = _new(output_csv, "output_csv")
    return _run(["sig", "manifest", sig, "-o", str(out)], out.parent, [out], sourmash_path)


def prefetch_signatures(
    query: str,
    databases: Sequence[str],
    output_csv: str,
    *,
    ksize: int = 21,
    threshold_bp: int = 0,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select the database signatures that overlap a query signature (sourmash prefetch)."""
    q = str(_existing(query, "query"))
    if not databases:
        raise ConfigurationError("databases must list at least one signature collection.")
    dbs = [str(_existing(d, "database")) for d in databases]
    out = _new(output_csv, "output_csv")
    args = ["prefetch", q, *dbs, "-k", str(ksize), "--dna", "--threshold-bp", str(threshold_bp), "-o", str(out)]
    return _run(args, out.parent, [out], sourmash_path)


def multigather_signatures(
    queries: Sequence[str],
    databases: Sequence[str],
    output_dir: str,
    *,
    ksize: int = 21,
    threshold_bp: int = 0,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a gather search for every query signature against the databases (sourmash multigather)."""
    if not queries or not databases:
        raise ConfigurationError("queries and databases must each list at least one file.")
    qs = [str(_existing(q, "query")) for q in queries]
    dbs = [str(_existing(d, "database")) for d in databases]
    folder = Path(output_dir).expanduser().resolve()
    folder.mkdir(parents=True, exist_ok=True)
    before = {p.name for p in folder.iterdir()}
    args = ["multigather", "--query", *qs, "--db", *dbs, "-k", str(ksize), "--dna",
            "--threshold-bp", str(threshold_bp), "--output-dir", str(folder)]
    result = _run(args, folder, [], sourmash_path)
    produced = sorted(p for p in folder.iterdir() if p.name not in before)
    result["outputs"] = [str(p) for p in produced]
    result["status"] = "completed" if result["return_code"] == 0 else "failed"
    return result


def categorize_signatures(
    database: str,
    queries: Sequence[str],
    output_csv: str,
    *,
    ksize: int = 21,
    threshold: float = 0.0,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign each query signature to its best match in a database (sourmash categorize)."""
    if not queries:
        raise ConfigurationError("queries must list at least one signature file.")
    if not 0.0 <= threshold <= 1.0:
        raise ConfigurationError("threshold must be between 0 and 1.")
    db = str(_existing(database, "database"))
    qs = [str(_existing(q, "query")) for q in queries]
    out = _new(output_csv, "output_csv")
    args = ["categorize", db, *qs, "-k", str(ksize), "--dna", "--threshold", str(threshold), "--csv", str(out)]
    return _run(args, out.parent, [out], sourmash_path)


def overlap_signatures(
    signature_a: str,
    signature_b: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report Jaccard similarity, containment and shared/unique hash counts of two signatures (sourmash sig overlap)."""
    a = _existing(signature_a, "signature_a")
    b = _existing(signature_b, "signature_b")
    result = execute_bounded_command(_program(sourmash_path), ["sig", "overlap", str(a), str(b)],
                                     backend_id="sourmash", cwd=a.parent, verbose=False)
    text = result.output or ""
    values: dict[str, float] = {}
    for line in text.splitlines():
        key, sep, value = line.partition(":")
        if sep:
            try:
                values[key.strip()] = float(value.split()[0])
            except (ValueError, IndexError):
                continue
    return {
        "status": "completed" if result.return_code == 0 else "failed",
        "backend": "sourmash",
        "return_code": result.return_code,
        "values": values,
        "report": text[-_MAX_TAIL:],
    }


def extract_signature_kmers(
    signature: str,
    sequences: Sequence[str],
    output_sequences: str,
    *,
    ksize: int = 21,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the sequences that contain the k-mers hashed in a signature (sourmash sig kmers)."""
    if not sequences:
        raise ConfigurationError("sequences must list at least one FASTA/FASTQ file.")
    sig = str(_existing(signature, "signature"))
    seqs = [str(_existing(s, "sequence file")) for s in sequences]
    out = _new(output_sequences, "output_sequences")
    args = ["sig", "kmers", "--signatures", sig, "--sequences", *seqs, "-k", str(ksize), "--save-sequences", str(out)]
    return _run(args, out.parent, [out], sourmash_path)


def export_signature_to_mash(
    signature: str,
    output_json: str,
    *,
    ksize: int = 21,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export a signature to the Mash JSON dump format (sourmash sig export)."""
    sig = str(_existing(signature, "signature"))
    out = _new(output_json, "output_json")
    return _run(["sig", "export", sig, "-k", str(ksize), "-o", str(out)], out.parent, [out], sourmash_path)


def sketch_translated_dna(
    fasta: str,
    output_signature: str,
    *,
    ksize: int = 10,
    scaled: int = 200,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Translate DNA in all six frames and sketch the protein k-mers (sourmash sketch translate)."""
    source = _existing(fasta, "fasta")
    out = _new(output_signature, "output_signature")
    if ksize < 1 or scaled < 1:
        raise ConfigurationError("ksize and scaled must be positive integers.")
    args = ["sketch", "translate", str(source), "-p", f"k={ksize},scaled={scaled}", "-o", str(out)]
    return _run(args, out.parent, [out], sourmash_path)


def plot_comparison_matrix(
    comparison_matrix: str,
    output_dir: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the clustered similarity heat map of a sourmash compare matrix as a PNG (sourmash plot)."""
    matrix = _existing(comparison_matrix, "comparison_matrix")
    folder = Path(output_dir).expanduser().resolve()
    folder.mkdir(parents=True, exist_ok=True)
    before = {p.name for p in folder.iterdir()}
    result = _run(["plot", str(matrix), "--output-dir", str(folder)], folder, [], sourmash_path)
    images = sorted(p for p in folder.iterdir() if p.name not in before and p.suffix == ".png")
    result["outputs"] = [str(p) for p in sorted(p for p in folder.iterdir() if p.name not in before)]
    result["status"] = "completed" if result["return_code"] == 0 and images else "failed"
    result["image"] = str(images[0]) if images else None
    return result


def _has_data_rows(csv_path: Path) -> bool:
    if not csv_path.is_file():
        return False
    return len([ln for ln in csv_path.read_text().splitlines() if ln.strip()]) > 1


def check_signatures_against_picklist(
    signatures: Sequence[str],
    picklist: str,
    output_missing_csv: str,
    *,
    output_manifest: str | None = None,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check which picklist rows are present in signature files and list the missing ones (sourmash sig check)."""
    if not signatures:
        raise ConfigurationError("signatures must list at least one file.")
    sigs = [str(_existing(f, "signature")) for f in signatures]
    missing = _new(output_missing_csv, "output_missing_csv")
    args = ["sig", "check", *sigs, "--picklist", picklist, "-o", str(missing), "--no-require-manifest"]
    outputs = [missing]
    if output_manifest is not None:
        manifest = _new(output_manifest, "output_manifest")
        args += ["-m", str(manifest)]
        outputs.append(manifest)
    result = execute_bounded_command(_program(sourmash_path), args, backend_id="sourmash",
                                     cwd=missing.parent, verbose=False)
    return {
        "status": "completed" if result.return_code in (0, 1) else "failed",
        "backend": "sourmash",
        "command": ["sourmash", *args],
        "return_code": result.return_code,
        "all_picklist_rows_found": result.return_code == 0 and not _has_data_rows(missing),
        "outputs": [str(p) for p in outputs if p.is_file()],
        "log_tail": (result.output or "")[-_MAX_TAIL:],
    }


def collect_signature_manifests(
    signatures: Sequence[str],
    output_manifest: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect the manifests of many signature files into one standalone manifest database (sourmash sig collect)."""
    if not signatures:
        raise ConfigurationError("signatures must list at least one file.")
    sigs = [str(_existing(f, "signature")) for f in signatures]
    out = _new(output_manifest, "output_manifest")
    return _run(["sig", "collect", *sigs, "-o", str(out), "--no-require-manifest"], out.parent, [out], sourmash_path)


def ingest_mash_csv(
    csv_files: Sequence[str],
    output_signature: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert Mash CSV hash dumps into a sourmash signature file (sourmash sig ingest --csv)."""
    if not csv_files:
        raise ConfigurationError("csv_files must list at least one file.")
    files = [str(_existing(f, "csv file")) for f in csv_files]
    out = _new(output_signature, "output_signature")
    return _run(["sig", "ingest", "--csv", *files, "-o", str(out)], out.parent, [out], sourmash_path)


def index_lca_database(
    signatures: Sequence[str],
    lineages_csv: str,
    output_database: str,
    *,
    ksize: int = 21,
    scaled: int = 1000,
    split_identifiers: bool = True,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a lowest-common-ancestor taxonomy database from signatures and a lineage CSV (sourmash lca index)."""
    if not signatures:
        raise ConfigurationError("signatures must list at least one file.")
    sigs = [str(_existing(f, "signature")) for f in signatures]
    lineages = str(_existing(lineages_csv, "lineages_csv"))
    out = _new(output_database, "output_database")
    args = ["lca", "index", lineages, str(out), *sigs, "-k", str(ksize), "--scaled", str(scaled), "--dna"]
    if split_identifiers:
        args.append("--split-identifiers")
    return _run(args, out.parent, [out], sourmash_path)


def classify_with_lca(
    queries: Sequence[str],
    databases: Sequence[str],
    output_csv: str,
    *,
    threshold: int = 5,
    majority: bool = False,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign a taxonomic lineage to each query signature with lowest-common-ancestor databases (sourmash lca classify)."""
    if not queries or not databases:
        raise ConfigurationError("queries and databases must each list at least one file.")
    qs = [str(_existing(q, "query")) for q in queries]
    dbs = [str(_existing(d, "database")) for d in databases]
    out = _new(output_csv, "output_csv")
    args = ["lca", "classify", "--query", *qs, "--db", *dbs, "--threshold", str(threshold), "-o", str(out)]
    if majority:
        args.append("--majority")
    return _run(args, out.parent, [out], sourmash_path)


__all__ = [
    "check_signatures_against_picklist",
    "classify_with_lca",
    "collect_signature_manifests",
    "index_lca_database",
    "ingest_mash_csv",
    "categorize_signatures",
    "export_signature_to_mash",
    "extract_signature_kmers",
    "multigather_signatures",
    "overlap_signatures",
    "plot_comparison_matrix",
    "prefetch_signatures",
    "sketch_translated_dna",
    "compare_signatures",
    "concatenate_signatures",
    "extract_signatures",
    "filter_signature_abundance",
    "manifest_signature",
    "split_signatures",
    "downsample_signature",
    "flatten_signature",
    "rename_signature",
    "signature_file_info",
    "intersect_signatures",
    "merge_signatures",
    "subtract_signatures",
    "describe_signature",
    "gather_signatures",
    "index_signatures",
    "search_signatures",
    "sketch_dna",
]
_register_public_verbose(globals())
