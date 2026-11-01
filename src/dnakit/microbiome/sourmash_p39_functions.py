"""sourmash operations that were still missing (integration pack P39).

Two groups of functions:

* command-line operations run through the ``sourmash`` program (``watch``, ``index -F rocksdb``, ``tax genome``,
  ``tax annotate``, ``tax prepare``, ``tax grep``, ``tax summarize`` and ``sig ingest`` for Mash JSON). No shell is
  used, the arguments are a fixed whitelist, every output goes to a new path and the result is a JSON-serialisable
  dictionary (status, backend, version, command, outputs, bounded log tail).
* Python-API operations on sketches (weighted containment, intersection and union size, k-mer hashes, size accuracy and
  the MurmurHash64 of a string). Sketches are plain mappings with ``ksize``, ``scaled``, ``hashes`` and optional
  ``abundances`` (one count per hash, in the order of ``hashes``).

sourmash is imported or started only when a function runs; a missing package raises ``BackendUnavailableError`` with the
``pip install "dnakit[bio-sourmash]"`` hint. Tested with sourmash 4.9.4.
"""

from __future__ import annotations

import importlib.metadata
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.metagenomics_review_sourmash_functions import _existing, _new, _nonempty, _run

_EXTRA = "bio-sourmash"
_HASH_LIMIT = 2**64 - 1
_RANKS = ("superkingdom", "phylum", "class", "order", "family", "genus", "species")
_GENOME_FORMATS = ("csv_summary", "krona", "human", "lineage_csv")
_INDEX_TYPES = ("SBT", "rocksdb", "zip")
_MOLTYPES = ("DNA", "protein", "dayhoff", "hp")
_MAX_HASHES = 5_000_000
_MAX_SEQUENCE = 10_000_000


def _version() -> str:
    try:
        return importlib.metadata.version("sourmash")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _positive(value: object, name: str, maximum: int = 10**9) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [1, {maximum}].")
    return value


def _finish(result: dict[str, Any], operation: str, **extra: Any) -> dict[str, Any]:
    result["operation"] = operation
    result["tool_version"] = _version()
    result.update(extra)
    return result


def _new_dir(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists() and not resolved.is_dir():
        raise ConfigurationError(f"{name} must be a folder: {resolved}")
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def _listing(folder: Path) -> set[Path]:
    return {item for item in folder.iterdir()}


def _run_collect(
    args: list[str], folder: Path, sourmash_path: str | None, verbose: bool
) -> dict[str, Any]:
    """Run sourmash and report the files that the run added to ``folder``."""

    before = _listing(folder)
    result = _run(args, [], sourmash_path, verbose)
    written = sorted(_listing(folder) - before)
    result["outputs"] = [str(item) for item in written]
    if result["return_code"] == 0 and not written:
        result["status"] = "failed"
    return result


def _identifier_flags(
    keep_full_identifiers: bool, keep_identifier_versions: bool, fail_on_missing_taxonomy: bool
) -> list[str]:
    flags: list[str] = []
    if keep_full_identifiers:
        flags.append("--keep-full-identifiers")
    if keep_identifier_versions:
        flags.append("--keep-identifier-versions")
    if fail_on_missing_taxonomy:
        flags.append("--fail-on-missing-taxonomy")
    return flags


def watch_reads_against_index(
    index: str,
    reads: str,
    output_signature: str,
    *,
    ksize: int | None = None,
    num_hashes: int = 500,
    threshold: float = 0.05,
    moltype: str = "DNA",
    name: str | None = None,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Stream FASTA/FASTQ reads into a growing sketch and stop at the first index match (sourmash watch).

    ``index`` is a sourmash SBT of ``num`` sketches (``sourmash index`` of signatures made with ``num=``);
    scaled indexes are rejected by upstream. ``moltype`` is DNA, protein, dayhoff or hp (sourmash watch
    selects hp unless a molecule type is given, so DNA is passed explicitly). The sketch built from the reads
    is saved to ``output_signature``; the best match (``FOUND: name, at similarity``) is parsed from the log.
    """

    path = _existing(index, "index")
    reads_path = _existing(reads, "reads")
    out = _new(output_signature, "output_signature")
    if moltype not in _MOLTYPES:
        raise ConfigurationError(f"moltype must be one of {list(_MOLTYPES)}.")
    n_hashes = _positive(num_hashes, "num_hashes", 50_000)
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not 0 <= threshold <= 1:
        raise ConfigurationError("threshold must be a number in [0, 1].")
    flag = {"DNA": "--dna", "protein": "--protein", "dayhoff": "--dayhoff", "hp": "--hp"}[moltype]
    args = ["watch", str(path), str(reads_path), flag, "-n", str(n_hashes), "--threshold", str(threshold)]
    if ksize is not None:
        args += ["-k", str(_positive(ksize, "ksize", 100))]
    if name is not None:
        args += ["--name", str(name)]
    args += ["-o", str(out)]
    result = _run(args, [out], sourmash_path, verbose)
    found = re.search(r"FOUND: (.+?), at ([0-9.]+)", result["log_tail"])
    match = {"name": found.group(1), "similarity": float(found.group(2))} if found else None
    return _finish(result, "watch_reads_against_index", match=match, threshold=threshold)


def build_signature_index(
    signatures: Sequence[str],
    output_index: str,
    *,
    index_type: str = "rocksdb",
    ksize: int = 31,
    moltype: str = "DNA",
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build an on-disk searchable signature index of a chosen type (sourmash index -F).

    ``index_type`` is ``rocksdb`` (an inverted index folder, the sourmash recommendation), ``SBT`` (a
    sequence Bloom tree, written to ``<output>.sbt.zip``) or ``zip`` (a zip collection with a manifest).
    All signatures must be scaled and share ``ksize`` and ``moltype``. Use ``index_type="rocksdb"`` for the
    inverted-index rows of the work list.
    """

    if index_type not in _INDEX_TYPES:
        raise ConfigurationError(f"index_type must be one of {list(_INDEX_TYPES)}.")
    if moltype not in _MOLTYPES:
        raise ConfigurationError(f"moltype must be one of {list(_MOLTYPES)}.")
    sigs = [str(_existing(item, "signature")) for item in _nonempty(signatures, "signatures")]
    out = _new(output_index, "output_index")
    flag = {"DNA": "--dna", "protein": "--protein", "dayhoff": "--dayhoff", "hp": "--hp"}[moltype]
    written = out.parent / (out.name + ".sbt.zip") if index_type == "SBT" else out
    if written.exists():
        raise ConfigurationError(f"output_index must not exist yet: {written}")
    args = ["index", str(out), *sigs, "-F", index_type, "-k", str(_positive(ksize, "ksize", 100)), flag]
    result = _run(args, [written], sourmash_path, verbose)
    return _finish(result, "build_signature_index", index_type=index_type, index=str(written))


def tax_genome_classification(
    gather_csvs: Sequence[str],
    taxonomy_csvs: Sequence[str],
    output_base: str,
    output_dir: str,
    *,
    output_formats: Sequence[str] = ("csv_summary",),
    rank: str | None = None,
    containment_threshold: float = 0.1,
    ani_threshold: float | None = None,
    keep_full_identifiers: bool = False,
    keep_identifier_versions: bool = False,
    fail_on_missing_taxonomy: bool = False,
    force: bool = False,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign one lineage to each query genome from gather results (sourmash tax genome).

    ``output_formats`` may hold ``csv_summary`` (the classification table ``query_name,status,rank,fraction,
    lineage,...``), ``krona`` (a Krona table, needs ``rank``), ``human`` (a text table) and ``lineage_csv`` (a lineage
    file usable as a taxonomy); one run writes one file per format into ``output_dir``. ``rank`` forces
    classification at one rank, otherwise ``containment_threshold`` (default 0.1) selects the lowest rank with
    enough support. ``ani_threshold`` uses average-nucleotide-identity instead of containment.
    """

    gathers = [str(_existing(item, "gather_csv")) for item in _nonempty(gather_csvs, "gather_csvs")]
    taxonomies = [str(_existing(item, "taxonomy_csv")) for item in _nonempty(taxonomy_csvs, "taxonomy_csvs")]
    formats = [str(item) for item in _nonempty(output_formats, "output_formats")]
    for item in formats:
        if item not in _GENOME_FORMATS:
            raise ConfigurationError(f"output_formats entries must be among {list(_GENOME_FORMATS)}.")
    if rank is not None and rank not in _RANKS:
        raise ConfigurationError(f"rank must be one of {list(_RANKS)}.")
    if "krona" in formats and rank is None:
        raise ConfigurationError("krona output needs a rank.")
    folder = _new_dir(output_dir, "output_dir")
    base = re.sub(r"[^A-Za-z0-9_.-]", "_", str(output_base))
    if not base:
        raise ConfigurationError("output_base must not be empty.")
    args = ["tax", "genome", "-g", *gathers, "-t", *taxonomies, "-o", base, "--output-dir", str(folder)]
    args += ["-F", *formats]
    if rank is not None:
        args += ["--rank", rank]
    if ani_threshold is not None:
        args += ["--ani-threshold", str(ani_threshold)]
    else:
        args += ["--containment-threshold", str(containment_threshold)]
    args += _identifier_flags(keep_full_identifiers, keep_identifier_versions, fail_on_missing_taxonomy)
    if force:
        args.append("--force")
    result = _run_collect(args, folder, sourmash_path, verbose)
    return _finish(result, "tax_genome_classification", output_formats=formats, rank=rank)


def tax_annotate_gather_results(
    gather_csvs: Sequence[str],
    taxonomy_csvs: Sequence[str],
    output_dir: str,
    *,
    keep_full_identifiers: bool = False,
    keep_identifier_versions: bool = False,
    fail_on_missing_taxonomy: bool = False,
    force: bool = False,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add a lineage column to gather CSV files (sourmash tax annotate); writes ``<name>.with-lineages.csv``."""

    gathers = [str(_existing(item, "gather_csv")) for item in _nonempty(gather_csvs, "gather_csvs")]
    taxonomies = [str(_existing(item, "taxonomy_csv")) for item in _nonempty(taxonomy_csvs, "taxonomy_csvs")]
    folder = _new_dir(output_dir, "output_dir")
    args = ["tax", "annotate", "-g", *gathers, "-t", *taxonomies, "-o", str(folder)]
    args += _identifier_flags(keep_full_identifiers, keep_identifier_versions, fail_on_missing_taxonomy)
    if force:
        args.append("--force")
    result = _run_collect(args, folder, sourmash_path, verbose)
    return _finish(result, "tax_annotate_gather_results")


def tax_prepare_taxonomy_database(
    taxonomy_csvs: Sequence[str],
    output: str,
    *,
    database_format: str = "sql",
    keep_full_identifiers: bool = False,
    keep_identifier_versions: bool = False,
    fail_on_missing_taxonomy: bool = False,
    force: bool = False,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine taxonomy spreadsheets into one new taxonomy database, SQLite (``sql``) or CSV (sourmash tax prepare)."""

    if database_format not in ("csv", "sql"):
        raise ConfigurationError("database_format must be 'csv' or 'sql'.")
    taxonomies = [str(_existing(item, "taxonomy_csv")) for item in _nonempty(taxonomy_csvs, "taxonomy_csvs")]
    out = _new(output, "output")
    args = ["tax", "prepare", "-t", *taxonomies, "-o", str(out), "-F", database_format]
    args += _identifier_flags(keep_full_identifiers, keep_identifier_versions, fail_on_missing_taxonomy)
    if force:
        args.append("--force")
    result = _run(args, [out], sourmash_path, verbose)
    return _finish(result, "tax_prepare_taxonomy_database", database_format=database_format)


def tax_grep_lineages(
    pattern: str,
    taxonomy_csvs: Sequence[str],
    *,
    output: str | None = None,
    rank: str | None = None,
    invert_match: bool = False,
    ignore_case: bool = False,
    force: bool = False,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select the taxonomy lineages whose rank names match a string (sourmash tax grep).

    The matching lineages are saved to the new CSV ``output`` (usable as a taxonomy or picklist). Without
    ``output`` only the number of matching lineages is counted (``-c``). ``rank`` limits the search to
    one rank; ``invert_match`` selects the non-matching lineages.
    """

    if not isinstance(pattern, str) or not pattern:
        raise ConfigurationError("pattern must be a non-empty string.")
    if rank is not None and rank not in _RANKS:
        raise ConfigurationError(f"rank must be one of {list(_RANKS)}.")
    taxonomies = [str(_existing(item, "taxonomy_csv")) for item in _nonempty(taxonomy_csvs, "taxonomy_csvs")]
    args = ["tax", "grep", pattern, "-t", *taxonomies]
    out = None
    if output is None:
        args.append("--count")
    else:
        out = _new(output, "output")
        args += ["-o", str(out)]
    if rank is not None:
        args += ["--rank", rank]
    if invert_match:
        args.append("--invert-match")
    if ignore_case:
        args.append("--ignore-case")
    if force:
        args.append("--force")
    result = _run(args, [out] if out else [], sourmash_path, verbose)
    count = re.search(r"found (\d+) matches", result["log_tail"])
    return _finish(result, "tax_grep_lineages", match_count=int(count.group(1)) if count else None)


def tax_summarize_lineages(
    taxonomy_files: Sequence[str],
    *,
    output_lineage_counts: str | None = None,
    keep_full_identifiers: bool = False,
    keep_identifier_versions: bool = False,
    force: bool = False,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize taxonomy databases or lineage spreadsheets: number of distinct lineages per rank (sourmash tax summarize).

    ``summary`` holds the lines printed by sourmash, ``rank_counts`` the parsed ``{rank: distinct lineages}``. With
    ``output_lineage_counts`` a new CSV ``rank,lineage_count,lineage`` is written.
    """

    files = [str(_existing(item, "taxonomy file")) for item in _nonempty(taxonomy_files, "taxonomy_files")]
    args = ["tax", "summarize", *files]
    out = None
    if output_lineage_counts is not None:
        out = _new(output_lineage_counts, "output_lineage_counts")
        args += ["-o", str(out)]
    args += _identifier_flags(keep_full_identifiers, keep_identifier_versions, False)
    if force:
        args.append("--force")
    result = _run(args, [out] if out else [], sourmash_path, verbose)
    text = result["log_tail"]
    rank_counts = {rank: int(count) for rank, count in re.findall(r"rank (\w+):\s+(\d+) distinct", text)}
    total = re.search(r"number of distinct taxonomic lineages: (\d+)", text)
    return _finish(
        result,
        "tax_summarize_lineages",
        distinct_lineages=int(total.group(1)) if total else None,
        rank_counts=rank_counts,
    )


def ingest_mash_json_signatures(
    json_files: Sequence[str],
    output_signature: str,
    *,
    sourmash_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert JSON sketch dumps written by ``mash info -d`` into one sourmash signature file (sourmash sig ingest)."""

    files = [str(_existing(item, "json file")) for item in _nonempty(json_files, "json_files")]
    out = _new(output_signature, "output_signature")
    result = _run(["sig", "ingest", *files, "-o", str(out)], [out], sourmash_path, verbose)
    return _finish(result, "ingest_mash_json_signatures")


def _sketch(sourmash: Any, value: object, name: str) -> Any:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping with ksize, scaled and hashes.")
    ksize = _positive(value.get("ksize"), f"{name} ksize", 100)
    scaled = _positive(value.get("scaled"), f"{name} scaled")
    hashes = value.get("hashes")
    if isinstance(hashes, (str, bytes)) or not isinstance(hashes, Sequence) or not hashes:
        raise ConfigurationError(f"{name} hashes must be a non-empty list of integers.")
    if len(hashes) > _MAX_HASHES:
        raise ConfigurationError(f"{name} may hold at most {_MAX_HASHES} hashes.")
    for item in hashes:
        if isinstance(item, bool) or not isinstance(item, int) or not 0 <= item <= _HASH_LIMIT:
            raise ConfigurationError(f"{name} hashes must be unsigned 64-bit integers.")
    abundances = value.get("abundances")
    try:
        minhash = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled, track_abundance=abundances is not None)
        if abundances is None:
            minhash.add_many(list(hashes))
        else:
            if len(abundances) != len(hashes):
                raise ConfigurationError(f"{name} abundances must have one value per hash.")
            minhash.set_abundances(dict(zip(hashes, abundances, strict=True)))
    except (OverflowError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not build {name}: {exc}") from exc
    return minhash


def sketch_weighted_containment(
    query: Mapping[str, Any],
    target: Mapping[str, Any],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Abundance-weighted fraction of the query sketch found in the target sketch (MinHash.contained_by_weighted).

    The query needs ``abundances`` (``track_abundance``); the weight of every query hash is its abundance, so
    the value is the share of the query's total abundance carried by hashes that also occur in the target. Sketches are
    downsampled to the coarser scale automatically, as sourmash does. Covers the MinHash and SourmashSignature
    entry points.
    """

    sourmash = require("sourmash", _EXTRA)
    left = _sketch(sourmash, query, "query")
    right = _sketch(sourmash, target, "target")
    if not left.track_abundance:
        raise ConfigurationError("query needs abundances for weighted containment.")
    try:
        value = float(left.contained_by_weighted(right))
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ConfigurationError(f"sourmash could not compute weighted containment: {exc}") from exc
    return {
        "status": "completed",
        "operation": "sketch_weighted_containment",
        "weighted_containment": value,
        "query_hashes": len(left),
        "target_hashes": len(right),
        "tool_version": _version(),
    }


def sketch_intersection_union_size(
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Number of shared hashes and the size of the union of two sketches (MinHash.count_common and intersection_and_union_size).

    ``shared_hashes`` is ``count_common`` (both sketches need the same ksize and scaled); ``intersection`` and ``union`` come
    from ``intersection_and_union_size``; the Jaccard index is their ratio.
    """

    sourmash = require("sourmash", _EXTRA)
    first = _sketch(sourmash, left, "left")
    second = _sketch(sourmash, right, "right")
    try:
        common = int(first.count_common(second))
        intersection, union = first.intersection_and_union_size(second)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not compare the sketches: {exc}") from exc
    return {
        "status": "completed",
        "operation": "sketch_intersection_union_size",
        "shared_hashes": common,
        "intersection": int(intersection),
        "union": int(union),
        "jaccard": (int(intersection) / int(union)) if union else 0.0,
        "tool_version": _version(),
    }


def sequence_kmer_hashes(
    sequence: str,
    ksize: int,
    *,
    include_kmers: bool = False,
    force: bool = False,
    seed: int = 42,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hash every k-mer of a DNA sequence without making a sketch (MinHash.seq_to_hashes, kmers_and_hashes).

    With ``include_kmers=False`` the hash values come from ``seq_to_hashes`` (the sequence is used as given). With
    ``include_kmers=True`` the (k-mer, hash) pairs come from ``kmers_and_hashes`` (the sequence is upper-cased). Hashes are
    canonical-k-mer MurmurHash64 values, one per window; with ``force`` windows with invalid bases give 0 or None.
    """

    if not isinstance(sequence, str) or not sequence or len(sequence) > _MAX_SEQUENCE:
        raise ConfigurationError(f"sequence must be a DNA string of 1 to {_MAX_SEQUENCE} bases.")
    k = _positive(ksize, "ksize", 100)
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ConfigurationError("seed must be a non-negative integer.")
    sourmash = require("sourmash", _EXTRA)
    try:
        minhash = sourmash.MinHash(n=0, ksize=k, scaled=1, seed=seed)
        if include_kmers:
            pairs = minhash.kmers_and_hashes(sequence, force=bool(force))
            kmers = [[kmer, None if h is None else int(h)] for kmer, h in pairs]
            hashes = [item[1] for item in kmers]
        else:
            hashes = [int(h) for h in minhash.seq_to_hashes(sequence, force=bool(force), bad_kmers_as_zeroes=bool(force))]
            kmers = None
    except (OverflowError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not hash the sequence: {exc}") from exc
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "sequence_kmer_hashes",
        "ksize": k,
        "seed": seed,
        "kmer_count": len(hashes),
        "hashes": hashes,
        "tool_version": _version(),
    }
    if kmers is not None:
        result["kmers"] = kmers
    return result


def sketch_size_accuracy(
    sketch: Mapping[str, Any],
    *,
    relative_error: float = 0.2,
    confidence: float = 0.95,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tell whether a scaled sketch is large enough to estimate its k-mer set size (MinHash.size_is_accurate).

    Returns True when the probability that ``sketch size * scaled`` lies within ``relative_error`` of the true set size
    is at least ``confidence`` (sourmash defaults 0.2 and 0.95).
    """

    for value, label in ((relative_error, "relative_error"), (confidence, "confidence")):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1:
            raise ConfigurationError(f"{label} must be a number in [0, 1].")
    sourmash = require("sourmash", _EXTRA)
    minhash = _sketch(sourmash, sketch, "sketch")
    try:
        accurate = bool(minhash.size_is_accurate(relative_error=relative_error, confidence=confidence))
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not evaluate the sketch: {exc}") from exc
    return {
        "status": "completed",
        "operation": "sketch_size_accuracy",
        "is_accurate": accurate,
        "sketch_hashes": len(minhash),
        "scaled": int(minhash.scaled),
        "relative_error": float(relative_error),
        "confidence": float(confidence),
        "tool_version": _version(),
    }


def murmur_hash_string(
    text: str,
    *,
    seed: int = 42,
    verbose: bool = True,
) -> dict[str, Any]:
    """MurmurHash64 of a string exactly as given, with sourmash's seed 42 by default (sourmash.minhash.hash_murmur).

    Unlike a canonical k-mer hash, the string is not reverse-complemented or upper-cased first.
    """

    if not isinstance(text, str) or not text:
        raise ConfigurationError("text must be a non-empty string.")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed <= _HASH_LIMIT:
        raise ConfigurationError("seed must be an unsigned 64-bit integer.")
    minhash = require("sourmash.minhash", _EXTRA)
    return {
        "status": "completed",
        "operation": "murmur_hash_string",
        "text": text,
        "seed": seed,
        "hash": int(minhash.hash_murmur(text, seed)),
        "tool_version": _version(),
    }


__all__ = [
    "build_signature_index",
    "ingest_mash_json_signatures",
    "murmur_hash_string",
    "sequence_kmer_hashes",
    "sketch_intersection_union_size",
    "sketch_size_accuracy",
    "sketch_weighted_containment",
    "tax_annotate_gather_results",
    "tax_genome_classification",
    "tax_grep_lineages",
    "tax_prepare_taxonomy_database",
    "tax_summarize_lineages",
    "watch_reads_against_index",
]

_register_public_verbose(globals(), tuple(__all__))
