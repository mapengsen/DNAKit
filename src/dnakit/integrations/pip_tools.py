"""Explicit adapters for optional, permissively licensed Python DNA tools.

These adapters keep third-party packages optional and expose JSON-friendly
contracts for the same Python/MCP registry used by native DNAKit functions.
They construct only allow-listed public objects; callers cannot provide Python
callbacks or arbitrary import paths through a JSON request.
"""

from __future__ import annotations

import inspect
import math
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require, sequence_text
from dnakit.exceptions import ConfigurationError
from dnakit.tools.serialization import serialize_for_agent


def _bounded_string(value: object, name: str, *, maximum: int = 1_000_000) -> str:
    if not isinstance(value, str) or not value or len(value) > maximum:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    if "\x00" in value:
        raise ConfigurationError(f"{name} cannot contain NUL characters.")
    return value


def _bounded_limit(value: object, name: str, *, maximum: int = 100_000) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [1, {maximum}].")
    return value


def _strand(value: object) -> str:
    if value not in {"+", "-"}:
        raise ConfigurationError("strand must be '+' or '-'.")
    return value


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def fetch_pyfaidx_region(
    reference: str | Path,
    sequence_id: str,
    start: int,
    end: int,
    *,
    strand: str = "+",
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch one zero-based, half-open FASTA interval through optional pyfaidx."""

    source = input_path(reference)
    identifier = _bounded_string(sequence_id, "sequence_id", maximum=512)
    if isinstance(start, bool) or not isinstance(start, int) or start < 0:
        raise ConfigurationError("start must be a non-negative integer.")
    if isinstance(end, bool) or not isinstance(end, int) or end <= start:
        raise ConfigurationError("end must be an integer greater than start.")
    direction = _strand(strand)
    pyfaidx = require("pyfaidx", "bio-pyfaidx")
    fasta = None
    try:
        fasta = pyfaidx.Fasta(
            str(source),
            as_raw=True,
            sequence_always_upper=True,
        )
        sequence = str(fasta[identifier][start:end])
    except (KeyError, IndexError, OSError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pyfaidx could not fetch the requested interval: {exc}") from exc
    finally:
        close = getattr(fasta, "close", None)
        if callable(close):
            close()
    sequence = sequence_text(sequence)
    if direction == "-":
        sequence = _reverse_complement(sequence)
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "pyfaidx",
            "reference": str(source),
            "sequence_id": identifier,
            "start": start,
            "end": end,
            "strand": direction,
            "sequence": sequence,
            "provenance": {
                "package": "pyfaidx",
                "license": "BSD-3-Clause",
                "citation": "Shirley et al. (2015), DOI 10.7287/peerj.preprints.970v1",
                "coordinate_system": "zero-based half-open",
            },
        },
        verbose=verbose,
    )


def _vcf_value(record: Any, name: str, default: Any = None) -> Any:
    value = getattr(record, name, default)
    if value is not default:
        return value
    if isinstance(record, Mapping):
        return record.get(name, default)
    return default


def _vcf_alt_values(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        return [str(value)]
    if isinstance(value, Sequence):
        result = []
        for item in value:
            result.append(str(getattr(item, "value", item)))
        return result
    return [str(getattr(value, "value", value))]


def _vcf_record_summary(record: Any) -> tuple[str, str, list[str], str]:
    chromosome = str(_vcf_value(record, "CHROM", _vcf_value(record, "chrom", "")))
    position = str(_vcf_value(record, "POS", _vcf_value(record, "pos", "")))
    reference = str(_vcf_value(record, "REF", _vcf_value(record, "ref", "")))
    alternate = _vcf_alt_values(_vcf_value(record, "ALT", _vcf_value(record, "alts", [])))
    if not chromosome or not position or not reference or not alternate:
        raise ConfigurationError("VCF backend returned a record without CHROM/POS/REF/ALT.")
    kinds = {
        "snp"
        if len(reference) == len(alt) == 1
        else "indel"
        if len(reference) != len(alt)
        else "other"
        for alt in alternate
    }
    return chromosome, position, alternate, "+".join(sorted(kinds))


def _vcf_summary(
    reader: Any,
    *,
    package: str,
    source: Path,
    max_records: int,
    verbose: bool,
) -> dict[str, Any]:
    chromosomes: Counter[str] = Counter()
    kinds: Counter[str] = Counter()
    preview: list[dict[str, Any]] = []
    records = 0
    sample_names: object = []
    try:
        for record in reader:
            if records >= max_records:
                raise ConfigurationError(
                    f"VCF input exceeds max_records={max_records}; increase the explicit limit."
                )
            chromosome, position, alternate, kind = _vcf_record_summary(record)
            chromosomes[chromosome] += 1
            kinds[kind] += 1
            if len(preview) < 100:
                preview.append(
                    {
                        "chrom": chromosome,
                        "pos": int(position) if position.isdigit() else position,
                        "ref": str(_vcf_value(record, "REF", _vcf_value(record, "ref", ""))),
                        "alt": alternate,
                        "kind": kind,
                    }
                )
            records += 1
        samples = getattr(getattr(reader, "header", None), "samples", None)
        sample_names = getattr(samples, "names", None)
        if sample_names is None:
            sample_names = getattr(reader, "samples", [])
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"{package} VCF parsing failed: {exc}") from exc
    finally:
        close = getattr(reader, "close", None)
        if callable(close):
            close()
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": package,
            "source": str(source),
            "records": records,
            "chromosomes": dict(chromosomes),
            "variant_kinds": dict(kinds),
            "samples": [str(item) for item in sample_names],
            "preview": preview,
            "provenance": {
                "package": package,
                "license": "MIT",
                "coordinate_system": "VCF 1-based POS",
            },
        },
        verbose=verbose,
    )


def summarize_vcf_vcfpy(
    source: str | Path,
    *,
    max_records: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a VCF through the optional pure-Python vcfpy backend."""

    path = input_path(source)
    limit = _bounded_limit(max_records, "max_records")
    reader = require("vcfpy", "bio-vcfpy").Reader.from_path(str(path))
    return _vcf_summary(
        reader,
        package="vcfpy",
        source=path,
        max_records=limit,
        verbose=verbose,
    )


def summarize_vcf_cyvcf2(
    source: str | Path,
    *,
    max_records: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a VCF through the optional cyvcf2 backend."""

    path = input_path(source)
    limit = _bounded_limit(max_records, "max_records")
    reader = require("cyvcf2", "bio-cyvcf2").VCF(str(path))
    return _vcf_summary(
        reader,
        package="cyvcf2",
        source=path,
        max_records=limit,
        verbose=verbose,
    )


def align_biotite(
    query: str,
    target: str,
    *,
    mode: str = "global",
    gap_penalty: int = -10,
    max_alignments: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align DNA with Biotite's NUC matrix and a linear gap penalty.

    Input pairs are limited to four million DP cells and one million possible
    returned alignment columns. CIGAR treats target as reference and query as
    segment. Returned spans are zero-based half-open query/target coordinates.
    """

    query_text = sequence_text(_bounded_string(query, "query", maximum=100_000))
    target_text = sequence_text(_bounded_string(target, "target", maximum=100_000))
    if not isinstance(mode, str) or mode not in {"global", "local"}:
        raise ConfigurationError("mode must be 'global' or 'local'.")
    if (
        isinstance(gap_penalty, bool)
        or not isinstance(gap_penalty, int)
        or not -10_000 <= gap_penalty < 0
    ):
        raise ConfigurationError("gap_penalty must be an integer in [-10000, -1].")
    limit = _bounded_limit(max_alignments, "max_alignments", maximum=1_000)
    if (len(query_text) + 1) * (len(target_text) + 1) > 4_000_000:
        raise ConfigurationError("Biotite input exceeds the four million DP cell limit.")
    if (len(query_text) + len(target_text)) * limit > 1_000_000:
        raise ConfigurationError("Reduce max_alignments: output could exceed one million columns.")
    sequence = require("biotite.sequence", "bio-biotite")
    align = require("biotite.sequence.align", "bio-biotite")
    try:
        seq1 = sequence.NucleotideSequence(query_text)
        seq2 = sequence.NucleotideSequence(target_text)
        matrix = align.SubstitutionMatrix.std_nucleotide_matrix()
        alignments = align.align_optimal(
            seq1,
            seq2,
            matrix,
            gap_penalty=gap_penalty,
            local=mode == "local",
            max_number=limit,
        )
        results = []
        for alignment in alignments[:limit]:
            symbols = align.get_symbols(alignment)
            rendered = [
                "".join("-" if item is None else str(item) for item in row) for row in symbols
            ]
            trace = alignment.trace.tolist()
            spans = {}
            for index, name in enumerate(("query", "target")):
                positions = [row[index] for row in trace if row[index] >= 0]
                spans[name] = [min(positions), max(positions) + 1] if positions else None
            results.append(
                {
                    "query": rendered[0],
                    "target": rendered[1],
                    "score": int(alignment.score),
                    "identity": (
                        float(align.get_sequence_identity(alignment, mode="all"))
                        if rendered[0] else None
                    ),
                    "cigar": (
                        str(align.write_alignment_to_cigar(
                            alignment, reference_index=1, segment_index=0,
                            distinguish_matches=True, include_terminal_gaps=True,
                        )) if rendered[0] else None
                    ),
                    "spans": spans,
                }
            )
    except (AttributeError, TypeError, ValueError) as exc:
        raise ConfigurationError(
            f"Biotite could not align the nucleotide sequences: {exc}"
        ) from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "biotite",
            "mode": mode,
            "query_length": len(query_text),
            "target_length": len(target_text),
            "max_alignments": limit,
            "gap_penalty": gap_penalty,
            "alignments": results,
            "provenance": {
                "package": "biotite",
                "license": "BSD-3-Clause",
                "citation": "10.1186/s12859-018-2367-z",
                "matrix": "NUC",
                "identity_denominator": "all alignment columns",
                "cigar_reference": "target",
                "coordinate_system": "zero-based half-open spans",
            },
        },
        verbose=verbose,
    )


def align_edlib(
    query: str,
    target: str,
    *,
    mode: str = "NW",
    task: str = "path",
    k: int = -1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align DNA with Edlib; NW is global, HW infix and SHW prefix.

    Locations use zero-based inclusive target coordinates. IUPAC symbols are
    literal characters, not wildcards. A distance of -1 means the distance
    exceeds k, and no pretty alignment is then attempted.
    """

    query_text = sequence_text(_bounded_string(query, "query", maximum=1_000_000))
    target_text = sequence_text(_bounded_string(target, "target", maximum=1_000_000))
    if not isinstance(mode, str) or mode not in {"NW", "HW", "SHW"}:
        raise ConfigurationError("mode must be one of 'NW', 'HW', or 'SHW'.")
    if not isinstance(task, str) or task not in {"distance", "locations", "path"}:
        raise ConfigurationError("task must be one of 'distance', 'locations', or 'path'.")
    if isinstance(k, bool) or not isinstance(k, int) or not -1 <= k <= 2_147_483_647:
        raise ConfigurationError("k must be -1 or a non-negative 32-bit integer.")
    edlib = require("edlib", "bio-edlib")
    try:
        result = dict(edlib.align(query_text, target_text, mode=mode, task=task, k=k))
        if task == "path" and result.get("editDistance", -1) >= 0 and result.get("cigar"):
            result["nice_alignment"] = edlib.getNiceAlignment(
                result,
                query_text,
                target_text,
            )
    except (AttributeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"Edlib could not align the nucleotide sequences: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "edlib",
            "mode": mode,
            "task": task,
            "k": k,
            "found": result.get("editDistance", -1) >= 0,
            "query_length": len(query_text),
            "target_length": len(target_text),
            "result": result,
            "provenance": {
                "package": "edlib",
                "license": "MIT",
                "citation": "10.1093/bioinformatics/btw753",
                "coordinate_system": "zero-based inclusive target locations",
                "ambiguity": "literal IUPAC symbols; no wildcard equalities",
            },
        },
        verbose=verbose,
    )


def diff_geneblocks_sequences(
    reference: str,
    target: str,
    *,
    max_sequence_length: int = 100_000,
    max_blocks: int = 10_000,
    max_complexity: int = 4_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two bounded DNA strings through GeneBlocks ``DiffBlocks``.

    The adapter deliberately disables GeneBlocks' optional BLAST acceleration,
    so the result is deterministic and does not require an external executable.
    Inputs whose dynamic-programming complexity exceeds ``max_complexity`` are
    reported by GeneBlocks as one coarse ``change`` block.
    """

    if isinstance(max_sequence_length, bool) or not isinstance(max_sequence_length, int):
        raise ConfigurationError("max_sequence_length must be an integer.")
    if not 1 <= max_sequence_length <= 1_000_000:
        raise ConfigurationError("max_sequence_length must be in [1, 1000000].")
    reference_text = sequence_text(
        _bounded_string(reference, "reference", maximum=max_sequence_length)
    )
    target_text = sequence_text(
        _bounded_string(target, "target", maximum=max_sequence_length)
    )
    block_limit = _bounded_limit(max_blocks, "max_blocks", maximum=100_000)
    complexity_limit = _bounded_limit(
        max_complexity,
        "max_complexity",
        maximum=100_000_000,
    )
    geneblocks = require("geneblocks", "bio-geneblocks")
    diff_class = getattr(geneblocks, "DiffBlocks", None)
    from_sequences = getattr(diff_class, "from_sequences", None)
    if not callable(from_sequences):
        raise ConfigurationError(
            "The installed GeneBlocks package lacks DiffBlocks.from_sequences."
        )
    try:
        diff = from_sequences(
            reference_text,
            target_text,
            blast_over=None,
            max_complexity=complexity_limit,
        )
        blocks = list(getattr(diff, "blocks", ()))
    except (ImportError, OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"GeneBlocks sequence comparison failed: {exc}") from exc
    if len(blocks) > block_limit:
        raise ConfigurationError(
            f"GeneBlocks returned more than max_blocks={block_limit}; increase the explicit limit."
        )
    operations: Counter[str] = Counter()
    serialized_blocks: list[dict[str, Any]] = []
    for block in blocks:
        operation = str(getattr(block, "operation", ""))
        reference_location = getattr(block, "s1_location", None)
        target_location = getattr(block, "s2_location", None)
        if not operation or reference_location is None or target_location is None:
            raise ConfigurationError("GeneBlocks returned a block without operation or locations.")
        reference_start = int(getattr(reference_location, "start", 0))
        reference_end = int(getattr(reference_location, "end", 0))
        target_start = int(getattr(target_location, "start", 0))
        target_end = int(getattr(target_location, "end", 0))
        if not (0 <= reference_start <= reference_end <= len(reference_text)):
            raise ConfigurationError("GeneBlocks returned invalid reference coordinates.")
        if not (0 <= target_start <= target_end <= len(target_text)):
            raise ConfigurationError("GeneBlocks returned invalid target coordinates.")
        operations[operation] += 1
        serialized_blocks.append(
            {
                "operation": operation,
                "reference_span": [reference_start, reference_end],
                "target_span": [target_start, target_end],
                "reference": reference_text[reference_start:reference_end],
                "target": target_text[target_start:target_end],
            }
        )
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "GeneBlocks",
            "operation": "diff_sequences",
            "reference_length": len(reference_text),
            "target_length": len(target_text),
            "max_complexity": complexity_limit,
            "blocks": serialized_blocks,
            "block_count": len(serialized_blocks),
            "operations": dict(operations),
            "provenance": {
                "package": "geneblocks",
                "license": "MIT",
                "citation": (
                    "GeneBlocks documentation and Edinburgh Genome Foundry contributors; "
                    "no DOI supplied by the survey"
                ),
                "coordinate_system": "zero-based half-open spans",
                "backend_scope": "deterministic DiffBlocks comparison without BLAST acceleration",
            },
        },
        verbose=verbose,
    )


def compare_sourmash_sequences(
    reference: str,
    target: str,
    *,
    ksize: int = 31,
    scaled: int = 1_000,
    track_abundance: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two bounded DNA strings with a local sourmash sketch.

    This adapter uses only sourmash's in-memory ``MinHash`` API.  It accepts
    A/C/G/T sequences, never accesses a database or network service, and
    returns similarity and containment metrics plus bounded sketch previews.
    """

    if isinstance(ksize, bool) or not isinstance(ksize, int) or not 1 <= ksize <= 512:
        raise ConfigurationError("ksize must be an integer in [1, 512].")
    if isinstance(scaled, bool) or not isinstance(scaled, int) or not 1 <= scaled <= 1_000_000_000:
        raise ConfigurationError("scaled must be an integer in [1, 1000000000].")
    if not isinstance(track_abundance, bool):
        raise ConfigurationError("track_abundance must be a boolean.")
    reference_text = sequence_text(_bounded_string(reference, "reference", maximum=1_000_000))
    target_text = sequence_text(_bounded_string(target, "target", maximum=1_000_000))
    if set(reference_text) - set("ACGT") or set(target_text) - set("ACGT"):
        raise ConfigurationError("sourmash comparison accepts only A, C, G and T symbols.")
    if len(reference_text) < ksize or len(target_text) < ksize:
        raise ConfigurationError("Both sequences must be at least ksize bases long.")
    sourmash = require("sourmash", "bio-sourmash")
    minhash_class = getattr(sourmash, "MinHash", None)
    if not callable(minhash_class):
        raise ConfigurationError("The installed sourmash package lacks MinHash.")
    try:
        reference_sketch = minhash_class(
            n=0,
            ksize=ksize,
            scaled=scaled,
            track_abundance=track_abundance,
        )
        target_sketch = minhash_class(
            n=0,
            ksize=ksize,
            scaled=scaled,
            track_abundance=track_abundance,
        )
        reference_sketch.add_sequence(reference_text, force=True)
        target_sketch.add_sequence(target_text, force=True)
        similarity = float(
            reference_sketch.similarity(
                target_sketch,
                ignore_abundance=True,
                downsample=False,
            )
        )
        reference_containment = float(reference_sketch.contained_by(target_sketch))
        target_containment = float(target_sketch.contained_by(reference_sketch))
        reference_hashes = sorted(int(value) for value in reference_sketch.hashes)
        target_hashes = sorted(int(value) for value in target_sketch.hashes)
    except (AttributeError, ImportError, OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash DNA comparison failed: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "sourmash",
            "operation": "compare_sequences",
            "reference_length": len(reference_text),
            "target_length": len(target_text),
            "ksize": ksize,
            "scaled": scaled,
            "track_abundance": track_abundance,
            "similarity": similarity,
            "reference_containment": reference_containment,
            "target_containment": target_containment,
            "reference_sketch_size": len(reference_hashes),
            "target_sketch_size": len(target_hashes),
            "reference_hash_preview": reference_hashes[:100],
            "target_hash_preview": target_hashes[:100],
            "provenance": {
                "package": "sourmash",
                "license": "BSD-3-Clause",
                "citation": "10.21105/joss.06830",
                "coordinate_system": "k-mer sketches; no genomic coordinates",
                "backend_scope": (
                    "local in-memory MinHash comparison; no database or network access"
                ),
            },
        },
        verbose=verbose,
    )


def _pyranges_rows(value: object) -> list[dict[str, Any]]:
    frame = getattr(value, "df", value)
    to_dict = getattr(frame, "to_dict", None)
    if callable(to_dict):
        rows = to_dict(orient="records")
    elif isinstance(frame, Sequence) and not isinstance(frame, (str, bytes)):
        rows = list(frame)
    else:
        raise ConfigurationError("PyRanges backend did not return a tabular interval object.")
    if not isinstance(rows, list) or any(not isinstance(row, Mapping) for row in rows):
        raise ConfigurationError("PyRanges backend returned non-record interval rows.")
    normalized = []
    for row in rows:
        row = {
            key: None if isinstance(value, float) and math.isnan(value) else value
            for key, value in row.items()
        }
        start, end = row.get("Start"), row.get("End")
        if (
            not isinstance(row.get("Chromosome"), str)
            or not row["Chromosome"]
            or isinstance(start, bool) or not isinstance(start, int)
            or isinstance(end, bool) or not isinstance(end, int)
            or not 0 <= start < end
        ):
            raise ConfigurationError("BED rows require Chromosome and 0 <= Start < End integers.")
        normalized.append(row)
    return normalized


def _read_pyranges_bed(pyranges: Any, path: Path, limit: int) -> Any:
    try:
        ranges = pyranges.read_bed(str(path), as_df=False, nrows=limit + 1)
        rows = _pyranges_rows(ranges)
        if len(rows) > limit:
            raise ConfigurationError(
                f"BED input exceeds max_rows={limit}; no partial result returned."
            )
        return ranges
    except (AttributeError, OSError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"PyRanges could not read the BED input: {exc}") from exc


def summarize_pyranges_bed(
    source: str | Path,
    *,
    max_rows: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a bounded BED file through optional PyRanges."""

    path = input_path(source)
    limit = _bounded_limit(max_rows, "max_rows")
    pyranges = require("pyranges", "bio-pyranges")
    ranges = _read_pyranges_bed(pyranges, path, limit)
    rows = _pyranges_rows(ranges)
    chromosomes = sorted({str(row.get("Chromosome", "")) for row in rows if row.get("Chromosome")})
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "pyranges",
            "operation": "summarize_bed",
            "source": str(path),
            "row_count": len(rows),
            "columns": sorted({str(key) for row in rows for key in row}),
            "chromosomes": chromosomes,
            "preview": rows[:100],
            "provenance": {
                "package": "pyranges",
                "license": "MIT",
                "citation": "10.1093/bioinformatics/btz615",
                "coordinate_system": "zero-based half-open BED intervals",
            },
        },
        verbose=verbose,
    )


def merge_pyranges_bed(
    source: str | Path,
    *,
    slack: int = 0,
    count: bool = True,
    stranded: bool | None = None,
    max_rows: int = 20_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge BED intervals with PyRanges 0.x; reject inputs over max_rows.

    slack=0 includes bookended intervals; -1 requires overlap. stranded=None
    follows PyRanges auto-detection; True requires valid strands on all rows.
    Input metadata is discarded by merge except Strand and optional Count.
    """

    path = input_path(source)
    limit = _bounded_limit(max_rows, "max_rows", maximum=20_000)
    if isinstance(slack, bool) or not isinstance(slack, int) or not -1 <= slack <= 1_000_000_000:
        raise ConfigurationError("slack must be an integer in [-1, 1000000000].")
    if not isinstance(count, bool):
        raise ConfigurationError("count must be a boolean.")
    if stranded is not None and not isinstance(stranded, bool):
        raise ConfigurationError("stranded must be a boolean or None.")
    pyranges = require("pyranges", "bio-pyranges")
    ranges = _read_pyranges_bed(pyranges, path, limit)
    if stranded and any(row.get("Strand") not in {"+", "-"} for row in _pyranges_rows(ranges)):
        raise ConfigurationError("stranded=True requires '+' or '-' in every BED row.")
    merge = getattr(ranges, "merge", None)
    if not callable(merge):
        raise ConfigurationError("PyRanges backend does not provide merge().")
    try:
        merged = merge(slack=slack, count=count, strand=stranded)
        rows = _pyranges_rows(merged)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"PyRanges could not merge the BED intervals: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "pyranges",
            "operation": "merge_bed",
            "source": str(path),
            "input_max_rows": limit,
            "slack": slack,
            "count": count,
            "stranded": stranded,
            "row_count": len(rows),
            "intervals": rows,
            "provenance": {
                "package": "pyranges",
                "license": "MIT",
                "citation": "10.1093/bioinformatics/btz615",
                "coordinate_system": "zero-based half-open BED intervals",
            },
        },
        verbose=verbose,
    )


def summarize_pyfastx(
    source: str | Path,
    *,
    max_records: int = 100_000,
    preview_records: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a bounded FASTA/FASTQ stream through optional pyfastx.

    ``pyfastx.Fastx`` performs format detection without building an index.  The
    adapter keeps only aggregate statistics and a small record preview, and it
    rejects input beyond ``max_records`` instead of silently truncating it.
    """

    path = input_path(source)
    limit = _bounded_limit(max_records, "max_records")
    preview_limit = _bounded_limit(preview_records, "preview_records", maximum=1_000)
    pyfastx = require("pyfastx", "bio-pyfastx")
    fastx_class = getattr(pyfastx, "Fastx", None)
    if fastx_class is None:
        raise ConfigurationError("The installed pyfastx package lacks Fastx().")
    try:
        reader = fastx_class(str(path), uppercase=True)
    except (OSError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pyfastx could not open the FASTX input: {exc}") from exc
    records = 0
    total_bases = 0
    gc_bases = 0
    min_length: int | None = None
    max_length = 0
    preview: list[dict[str, Any]] = []
    fastq_records = 0
    try:
        for record in reader:
            if records >= limit:
                raise ConfigurationError(
                    f"FASTX input exceeds max_records={limit}; increase the explicit limit."
                )
            if not isinstance(record, Sequence) or isinstance(record, (str, bytes)):
                raise ConfigurationError("pyfastx returned a non-tuple FASTX record.")
            if len(record) < 2:
                raise ConfigurationError("pyfastx returned a FASTX record without a sequence.")
            name = _bounded_string(str(record[0]), "record name", maximum=512)
            sequence = str(record[1]).upper()
            if not sequence:
                raise ConfigurationError("pyfastx returned an empty sequence.")
            length = len(sequence)
            total_bases += length
            gc_bases += sequence.count("G") + sequence.count("C")
            min_length = length if min_length is None else min(min_length, length)
            max_length = max(max_length, length)
            if len(record) >= 3:
                fastq_records += 1
            if len(preview) < preview_limit:
                preview.append({"name": name, "length": length})
            records += 1
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pyfastx FASTX iteration failed: {exc}") from exc
    finally:
        close = getattr(reader, "close", None)
        if callable(close):
            close()
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "pyfastx",
            "operation": "summarize_fastx",
            "source": str(path),
            "records": records,
            "format": "fastq" if fastq_records else "fasta_or_empty",
            "fastq_records": fastq_records,
            "total_bases": total_bases,
            "min_length": min_length,
            "max_length": max_length,
            "mean_length": (total_bases / records) if records else None,
            "gc_bases": gc_bases,
            "gc_percent": (100.0 * gc_bases / total_bases) if total_bases else None,
            "preview": preview,
            "provenance": {
                "package": "pyfastx",
                "license": "MIT",
                "citation": "10.1093/bib/bbaa368",
                "coordinate_system": "record-level FASTA/FASTQ summary",
                "indexing": "Fastx iteration; no index requested",
            },
        },
        verbose=verbose,
    )


def summarize_dnaio_fastx(
    source: str | Path,
    *,
    max_records: int = 100_000,
    max_bases: int = 100_000_000,
    preview_records: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a bounded FASTA/FASTQ stream through optional ``dnaio``."""

    path = input_path(source)
    record_limit = _bounded_limit(max_records, "max_records")
    if isinstance(max_bases, bool) or not isinstance(max_bases, int) or not 1 <= max_bases <= 1_000_000_000:
        raise ConfigurationError("max_bases must be an integer in [1, 1000000000].")
    preview_limit = _bounded_limit(preview_records, "preview_records", maximum=1_000)
    dnaio = require("dnaio", "bio-dnaio")
    opener = getattr(dnaio, "open", None)
    if not callable(opener):
        raise ConfigurationError("The installed dnaio package lacks open().")
    try:
        reader = opener(str(path))
    except (OSError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"dnaio could not open the FASTX input: {exc}") from exc
    records = 0
    total_bases = 0
    gc_bases = 0
    min_length: int | None = None
    max_length = 0
    fastq_records = 0
    preview: list[dict[str, Any]] = []
    try:
        for record in reader:
            if records >= record_limit:
                raise ConfigurationError(
                    f"FASTX input exceeds max_records={record_limit}; increase the explicit limit."
                )
            sequence = getattr(record, "sequence", None)
            name = getattr(record, "name", None)
            qualities = getattr(record, "qualities", None)
            if sequence is None and isinstance(record, Sequence) and not isinstance(record, (str, bytes)):
                if len(record) < 2:
                    raise ConfigurationError("dnaio returned a record without a sequence.")
                name, sequence = record[0], record[1]
                qualities = record[2] if len(record) >= 3 else None
            if sequence is None:
                raise ConfigurationError("dnaio returned a record without a sequence.")
            sequence_text_value = str(sequence).upper()
            if not sequence_text_value:
                raise ConfigurationError("dnaio returned an empty sequence.")
            length = len(sequence_text_value)
            if total_bases + length > max_bases:
                raise ConfigurationError(
                    f"FASTX input exceeds max_bases={max_bases}; increase the explicit limit."
                )
            total_bases += length
            gc_bases += sequence_text_value.count("G") + sequence_text_value.count("C")
            min_length = length if min_length is None else min(min_length, length)
            max_length = max(max_length, length)
            if qualities is not None:
                fastq_records += 1
            if len(preview) < preview_limit:
                preview.append({"name": _bounded_string(str(name or ""), "record name", maximum=512), "length": length})
            records += 1
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"dnaio FASTX iteration failed: {exc}") from exc
    finally:
        close = getattr(reader, "close", None)
        if callable(close):
            close()
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "dnaio",
            "operation": "summarize_fastx",
            "source": str(path),
            "records": records,
            "format": "fastq" if fastq_records else "fasta_or_empty",
            "fastq_records": fastq_records,
            "total_bases": total_bases,
            "min_length": min_length,
            "max_length": max_length,
            "mean_length": (total_bases / records) if records else None,
            "gc_bases": gc_bases,
            "gc_percent": (100.0 * gc_bases / total_bases) if total_bases else None,
            "preview": preview,
            "provenance": {
                "package": "dnaio",
                "license": "MIT",
                "citation": "dnaio software archive, Zenodo DOI 10.5281/zenodo.10548864",
                "coordinate_system": "record-level FASTA/FASTQ summary",
                "backend_scope": "bounded read-only stream; no paired-end writer or mutation API",
            },
        },
        verbose=verbose,
    )


def map_mappy_sequences(
    reference: str | Path,
    reads: str | Path,
    *,
    preset: str = "map-ont",
    max_reads: int = 10_000,
    max_hits_per_read: int = 100,
    cs: bool = False,
    md: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map a bounded FASTA/FASTQ stream with optional mappy/minimap2 bindings."""

    reference_path = input_path(reference)
    reads_path = input_path(reads)
    presets = {"sr", "map-pb", "map-ont", "splice", "asm5", "asm10", "asm20"}
    if not isinstance(preset, str) or preset not in presets:
        raise ConfigurationError(f"preset must be one of {sorted(presets)}.")
    read_limit = _bounded_limit(max_reads, "max_reads", maximum=100_000)
    hit_limit = _bounded_limit(max_hits_per_read, "max_hits_per_read", maximum=1_000)
    if not isinstance(cs, bool) or not isinstance(md, bool):
        raise ConfigurationError("cs and md must be booleans.")
    mappy = require("mappy", "bio-mappy")
    aligner_class = getattr(mappy, "Aligner", None)
    fastx_read = getattr(mappy, "fastx_read", None)
    if not callable(aligner_class) or not callable(fastx_read):
        raise ConfigurationError("The installed mappy package lacks Aligner() or fastx_read().")
    try:
        aligner = aligner_class(str(reference_path), preset=preset, n_threads=1)
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"mappy could not load the reference: {exc}") from exc
    if not aligner:
        raise ConfigurationError("mappy failed to build or load the reference index.")
    mapped_reads = 0
    hits: list[dict[str, Any]] = []
    try:
        for item in fastx_read(str(reads_path)):
            if mapped_reads >= read_limit:
                raise ConfigurationError(
                    f"reads input exceeds max_reads={read_limit}; increase the explicit limit."
                )
            if not isinstance(item, Sequence) or isinstance(item, (str, bytes)) or len(item) < 2:
                raise ConfigurationError("mappy.fastx_read returned an invalid record.")
            name, sequence = str(item[0]), str(item[1])
            read_hits = list(aligner.map(sequence, cs=cs, MD=md))
            if len(read_hits) > hit_limit:
                raise ConfigurationError(
                    f"mappy returned more than max_hits_per_read={hit_limit} for {name!r}."
                )
            for hit in read_hits:
                hits.append(
                    {
                        "read": name,
                        "contig": str(getattr(hit, "ctg", "")),
                        "reference_start": int(getattr(hit, "r_st", 0)),
                        "reference_end": int(getattr(hit, "r_en", 0)),
                        "query_start": int(getattr(hit, "q_st", 0)),
                        "query_end": int(getattr(hit, "q_en", 0)),
                        "strand": int(getattr(hit, "strand", 0)),
                        "mapq": int(getattr(hit, "mapq", 0)),
                        "matching_bases": int(getattr(hit, "mlen", 0)),
                        "block_length": int(getattr(hit, "blen", 0)),
                        "cigar": str(getattr(hit, "cigar_str", "")),
                        "cs": str(getattr(hit, "cs", "")) if cs else None,
                        "md": str(getattr(hit, "MD", "")) if md else None,
                    }
                )
            mapped_reads += 1
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"mappy sequence mapping failed: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "mappy",
            "operation": "map_sequences",
            "reference": str(reference_path),
            "reads": str(reads_path),
            "preset": preset,
            "reads_mapped": mapped_reads,
            "hit_count": len(hits),
            "hits": hits,
            "provenance": {
                "package": "mappy/minimap2 Python binding",
                "license": "MIT",
                "citation": "10.1093/bioinformatics/bty191",
                "coordinate_system": "zero-based half-open reference/query spans",
                "backend_scope": "bounded local mapping; no index writing or arbitrary minimap2 flags",
            },
        },
        verbose=verbose,
    )


_PYBIGWIG_STATISTICS = frozenset({"mean", "min", "max", "coverage", "std"})


def summarize_pybigwig_region(
    source: str | Path,
    chromosome: str,
    start: int,
    end: int,
    *,
    statistic: str = "mean",
    bins: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a bounded local bigWig region through optional pyBigWig."""

    path = input_path(source)
    chrom = _bounded_string(chromosome, "chromosome", maximum=512)
    if isinstance(start, bool) or not isinstance(start, int) or start < 0:
        raise ConfigurationError("start must be a non-negative integer.")
    if isinstance(end, bool) or not isinstance(end, int) or end <= start:
        raise ConfigurationError("end must be an integer greater than start.")
    if not isinstance(statistic, str) or statistic not in _PYBIGWIG_STATISTICS:
        raise ConfigurationError(
            f"statistic must be one of {sorted(_PYBIGWIG_STATISTICS)}."
        )
    bin_count = _bounded_limit(bins, "bins", maximum=10_000)
    pybigwig = require("pyBigWig", "bio-pybigwig")
    try:
        reader = pybigwig.open(str(path))
    except (OSError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pyBigWig could not open the input: {exc}") from exc
    if reader is None:
        raise ConfigurationError("pyBigWig returned no reader for the input file.")
    try:
        chromosomes = reader.chroms()
        chromosome_length = chromosomes.get(chrom) if hasattr(chromosomes, "get") else None
        if chromosome_length is None:
            raise ConfigurationError(f"Chromosome {chrom!r} is not present in the bigWig file.")
        if end > int(chromosome_length):
            raise ConfigurationError(
                f"end={end} exceeds chromosome {chrom!r} length {int(chromosome_length)}."
            )
        values = reader.stats(
            chrom,
            start,
            end,
            type=statistic,
            nBins=bin_count,
        )
        if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
            raise ConfigurationError("pyBigWig returned a non-list statistics result.")
        normalized = [None if value is None else float(value) for value in values]
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pyBigWig region statistics failed: {exc}") from exc
    finally:
        close = getattr(reader, "close", None)
        if callable(close):
            close()
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "pyBigWig",
            "operation": "summarize_region",
            "source": str(path),
            "chromosome": chrom,
            "start": start,
            "end": end,
            "statistic": statistic,
            "bins": bin_count,
            "values": normalized,
            "provenance": {
                "package": "pyBigWig",
                "license": "MIT",
                "citation": "10.5281/zenodo.45238",
                "coordinate_system": "zero-based half-open",
                "backend_scope": "read-only local bigWig statistics",
            },
        },
        verbose=verbose,
    )


def intersect_pybedtools_bed(
    source: str | Path,
    other: str | Path,
    *,
    unique: bool = False,
    write_a: bool = False,
    write_b: bool = False,
    max_rows: int = 20_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Intersect two local BED-like files through pybedtools and BEDTools.

    ``pybedtools`` is a Python wrapper; the scientific operation still needs a
    separately installed ``bedtools`` executable.  The adapter allow-lists the
    three read-only intersection switches and streams at most ``max_rows``.
    """

    first = input_path(source)
    second = input_path(other)
    limit = _bounded_limit(max_rows, "max_rows", maximum=20_000)
    if not isinstance(unique, bool):
        raise ConfigurationError("unique must be a boolean.")
    if not isinstance(write_a, bool) or not isinstance(write_b, bool):
        raise ConfigurationError("write_a and write_b must be booleans.")
    pybedtools = require("pybedtools", "bio-pybedtools")
    bedtool_class = getattr(pybedtools, "BedTool", None)
    if bedtool_class is None:
        raise ConfigurationError("The installed pybedtools package lacks BedTool().")
    try:
        left = bedtool_class(str(first))
        right = bedtool_class(str(second))
        result = left.intersect(
            right,
            u=unique,
            wa=write_a,
            wb=write_b,
            stream=True,
        )
        rows: list[list[str]] = []
        for feature in result:
            if len(rows) >= limit:
                raise ConfigurationError(
                    f"BED intersection exceeds max_rows={limit}; increase the explicit limit."
                )
            line = str(feature).rstrip("\r\n")
            if line:
                rows.append(line.split("\t"))
    except ConfigurationError:
        raise
    except (AttributeError, OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pybedtools/BEDTools intersection failed: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "pybedtools",
            "operation": "intersect_bed",
            "source": str(first),
            "other": str(second),
            "unique": unique,
            "write_a": write_a,
            "write_b": write_b,
            "row_count": len(rows),
            "rows": rows,
            "provenance": {
                "package": "pybedtools",
                "license": "MIT",
                "citation": "10.1093/bioinformatics/btr539",
                "bedtools_citation": "10.1093/bioinformatics/btq033",
                "coordinate_system": "BED zero-based half-open when input is BED",
                "backend_scope": "pybedtools wrapper plus separately installed bedtools executable",
            },
        },
        verbose=verbose,
    )


def _gffutils_feature_summary(feature: Any) -> dict[str, Any]:
    fields = {
        name: getattr(feature, name, None)
        for name in (
            "id",
            "seqid",
            "source",
            "featuretype",
            "start",
            "end",
            "score",
            "strand",
            "frame",
        )
    }
    attributes = getattr(feature, "attributes", {})
    if isinstance(attributes, Mapping):
        fields["attributes"] = {
            str(key): (
                list(value)
                if isinstance(value, Sequence) and not isinstance(value, str)
                else value
            )
            for key, value in attributes.items()
        }
    else:
        fields["attributes"] = {}
    return fields


def query_gffutils_database(
    database: str | Path,
    *,
    feature_ids: Sequence[str] = (),
    chromosome: str | None = None,
    start: int | None = None,
    end: int | None = None,
    featuretype: str | None = None,
    max_features: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Query an existing gffutils database without creating or mutating it."""

    path = input_path(database)
    limit = _bounded_limit(max_features, "max_features", maximum=5_000)
    if isinstance(feature_ids, (str, bytes)) or not isinstance(feature_ids, Sequence):
        raise ConfigurationError("feature_ids must be a list of strings.")
    ids = [_bounded_string(value, "feature_id", maximum=512) for value in feature_ids]
    region_values = (chromosome, start, end)
    has_region = any(value is not None for value in region_values)
    if ids and has_region:
        raise ConfigurationError("Use feature_ids or a genomic region, not both.")
    if not ids and not has_region:
        raise ConfigurationError("Provide feature_ids or chromosome/start/end.")
    if len(ids) > limit:
        raise ConfigurationError(f"feature_ids exceeds max_features={limit}.")
    if has_region:
        chrom = _bounded_string(chromosome, "chromosome", maximum=512)
        if (
            isinstance(start, bool)
            or not isinstance(start, int)
            or start < 1
            or isinstance(end, bool)
            or not isinstance(end, int)
            or end < start
        ):
            raise ConfigurationError("GFF region requires 1-based inclusive start/end integers.")
    else:
        chrom = None
    if featuretype is not None:
        featuretype = _bounded_string(featuretype, "featuretype", maximum=128)
    gffutils = require("gffutils", "bio-gffutils")
    feature_db = getattr(gffutils, "FeatureDB", None)
    if feature_db is None:
        raise ConfigurationError("The installed gffutils package lacks FeatureDB().")
    try:
        db = feature_db(str(path), keep_order=True)
        missing_ids: list[str] = []
        features: list[dict[str, Any]] = []
        if ids:
            exceptions = getattr(gffutils, "exceptions", None)
            not_found = getattr(gffutils, "FeatureNotFoundError", None)
            if not_found is None:
                not_found = getattr(exceptions, "FeatureNotFoundError", KeyError)
            if not isinstance(not_found, type) or not issubclass(not_found, BaseException):
                not_found = KeyError
            for identifier in ids:
                try:
                    features.append(_gffutils_feature_summary(db[identifier]))
                except not_found:
                    missing_ids.append(identifier)
        else:
            kwargs = {"region": (chrom, start, end)}
            if featuretype is not None:
                kwargs["featuretype"] = featuretype
            for feature in db.region(**kwargs):
                if len(features) >= limit:
                    raise ConfigurationError(
                        f"GFF region exceeds max_features={limit}; increase the explicit limit."
                    )
                features.append(_gffutils_feature_summary(feature))
    except ConfigurationError:
        raise
    except (AttributeError, OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"gffutils database query failed: {exc}") from exc
    finally:
        close = locals().get("db")
        close_method = getattr(close, "close", None)
        if callable(close_method):
            close_method()
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "gffutils",
            "operation": "query_database",
            "database": str(path),
            "query": (
                {"feature_ids": ids}
                if ids
                else {
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "featuretype": featuretype,
                }
            ),
            "feature_count": len(features),
            "missing_ids": missing_ids,
            "features": features,
            "provenance": {
                "package": "gffutils",
                "license": "MIT",
                "citation": "gffutils documentation and daler contributors",
                "coordinate_system": "GFF/GTF 1-based inclusive feature coordinates",
                "backend_scope": "read-only existing FeatureDB; no database creation or mutation",
            },
        },
        verbose=verbose,
    )


def summarize_sbol3_document(
    source: str | Path,
    *,
    max_objects: int = 50_000,
    preview_objects: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize an existing SBOL3 document through the optional pySBOL3 backend.

    The adapter is intentionally read-only: it loads a local RDF document,
    inventories top-level objects, and returns bounded identity/type previews.
    It does not resolve remote RDF, mutate the document, or write a new file.
    """

    path = input_path(source)
    limit = _bounded_limit(max_objects, "max_objects")
    preview_limit = _bounded_limit(preview_objects, "preview_objects", maximum=1_000)
    sbol3 = require("sbol3", "bio-pysbol3")
    document = sbol3.Document()
    try:
        document.read(str(path))
        object_types: Counter[str] = Counter()
        preview: list[dict[str, str]] = []
        objects = 0
        for item in document:
            if objects >= limit:
                raise ConfigurationError(
                    f"SBOL3 document exceeds max_objects={limit}; increase the explicit limit."
                )
            object_types[type(item).__name__] += 1
            if len(preview) < preview_limit:
                identity = getattr(item, "identity", None)
                display_id = getattr(item, "display_id", None)
                preview.append(
                    {
                        "type": type(item).__name__,
                        "identity": str(identity) if identity is not None else "",
                        "display_id": str(display_id) if display_id is not None else "",
                    }
                )
            objects += 1
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pySBOL3 document parsing failed: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "pySBOL3",
            "source": str(path),
            "objects": objects,
            "object_types": dict(object_types),
            "preview": preview,
            "provenance": {
                "package": "sbol3",
                "license": "MIT",
                "citation": "PMID:35767721",
                "operation": "read-only SBOL3 inventory",
            },
        },
        verbose=verbose,
    )


def summarize_scadnano_design(
    source: str | Path,
    *,
    max_strands: int = 50_000,
    preview_strands: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a local scadnano ``.sc`` DNA-origami design.

    This adapter uses scadnano's documented file reader and only inspects the
    resulting design.  It does not edit designs, render a browser view, or
    download sequences or other design assets.
    """

    path = input_path(source)
    limit = _bounded_limit(max_strands, "max_strands")
    preview_limit = _bounded_limit(preview_strands, "preview_strands", maximum=1_000)
    scadnano = require("scadnano", "bio-scadnano")
    design_class = getattr(scadnano, "Design", None)
    loader = getattr(design_class, "from_scadnano_file", None)
    if not callable(loader):
        raise ConfigurationError(
            "The installed scadnano package lacks Design.from_scadnano_file."
        )
    try:
        design = loader(str(path))
        strands = list(getattr(design, "strands", ()))
        if len(strands) > limit:
            raise ConfigurationError(
                f"scadnano design exceeds max_strands={limit}; increase the explicit limit."
            )
        helices = getattr(design, "helices", {})
        groups = getattr(design, "groups", {})
        preview: list[dict[str, Any]] = []
        for index, strand in enumerate(strands[:preview_limit]):
            domains = list(getattr(strand, "domains", ()))
            sequence = getattr(strand, "sequence", None)
            item: dict[str, Any] = {
                "index": index,
                "domains": len(domains),
            }
            if sequence is not None:
                item["sequence_length"] = len(str(sequence))
            preview.append(item)
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"scadnano design parsing failed: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "scadnano",
            "source": str(path),
            "helices": len(helices),
            "groups": len(groups),
            "strands": len(strands),
            "preview": preview,
            "provenance": {
                "package": "scadnano",
                "license": "MIT",
                "citation": "arXiv:2005.11841",
                "operation": "read-only design inventory",
            },
        },
        verbose=verbose,
    )


def summarize_dna_features_viewer(
    source: str | Path,
    *,
    filetype: str | None = None,
    max_features: int = 50_000,
    preview_features: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a local GenBank or GFF record through DNA Features Viewer.

    The adapter uses the package's documented Biopython translator but keeps
    the public boundary read-only and JSON-safe.  Plot rendering and arbitrary
    translator callbacks are intentionally outside this interface.
    """

    path = input_path(source)
    limit = _bounded_limit(max_features, "max_features")
    preview_limit = _bounded_limit(preview_features, "preview_features", maximum=1_000)
    selected_filetype = filetype
    if selected_filetype is not None and selected_filetype not in {"genbank", "gff"}:
        raise ConfigurationError("filetype must be 'genbank' or 'gff'.")
    dna_features_viewer = require("dna_features_viewer", "bio-dna-features-viewer")
    translator_class = getattr(dna_features_viewer, "BiopythonTranslator", None)
    if not callable(translator_class):
        raise ConfigurationError(
            "The installed DNA Features Viewer package lacks BiopythonTranslator."
        )
    try:
        translator = translator_class()
        graphic_record = translator.translate_record(str(path), filetype=selected_filetype)
        features = list(getattr(graphic_record, "features", ()))
        if len(features) > limit:
            raise ConfigurationError(
                f"DNA Features Viewer record exceeds max_features={limit}; "
                "increase the explicit limit."
            )
        preview: list[dict[str, Any]] = []
        for feature in features[:preview_limit]:
            start = getattr(feature, "start", None)
            end = getattr(feature, "end", None)
            if start is None or end is None:
                raise ConfigurationError("DNA Features Viewer returned a feature without bounds.")
            label = getattr(feature, "label", None)
            preview.append(
                {
                    "start": int(start),
                    "end": int(end),
                    "strand": getattr(feature, "strand", None),
                    "label": "" if label is None else str(label),
                }
            )
    except ConfigurationError:
        raise
    except (ImportError, OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"DNA Features Viewer record parsing failed: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "DNA Features Viewer",
            "operation": "summarize_record",
            "source": str(path),
            "filetype": selected_filetype or "inferred",
            "sequence_length": int(getattr(graphic_record, "sequence_length", 0)),
            "features": len(features),
            "preview": preview,
            "provenance": {
                "package": "dna-features-viewer",
                "license": "MIT",
                "citation": "DNA Features Viewer documentation; no DOI supplied by the survey",
                "backend_scope": "read-only feature inventory; no plot file is written",
            },
        },
        verbose=verbose,
    )


def overlap_polars_bio_bed(
    source: str | Path,
    other: str | Path,
    *,
    distinct: bool = False,
    max_rows: int = 20_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find bounded overlaps between two local BED/CSV interval files.

    polars-bio reads the files through its interval API and returns a joined
    table.  The adapter does not accept remote URLs, arbitrary SQL, or Python
    expressions, and it collects only ``max_rows + 1`` rows to detect overflow.
    """

    source_path = input_path(source)
    other_path = input_path(other)
    limit = _bounded_limit(max_rows, "max_rows", maximum=20_000)
    if not isinstance(distinct, bool):
        raise ConfigurationError("distinct must be a boolean.")
    polars_bio = require("polars_bio", "bio-polarsbio")
    overlap = getattr(polars_bio, "overlap", None)
    if not callable(overlap):
        operations = getattr(polars_bio, "IntervalOperations", None)
        overlap = getattr(operations, "overlap", None)
    if not callable(overlap):
        raise ConfigurationError("The installed polars-bio package lacks overlap().")
    try:
        result = overlap(
            str(source_path),
            str(other_path),
            overlap_output="left",
            distinct_output=distinct,
            output_type="polars.LazyFrame",
        )
        bounded = result.limit(limit + 1) if hasattr(result, "limit") else result
        if hasattr(bounded, "collect"):
            bounded = bounded.collect()
        if hasattr(bounded, "to_dicts"):
            rows = bounded.to_dicts()
        elif hasattr(bounded, "to_dict"):
            rows = bounded.to_dict(orient="records")
        elif hasattr(bounded, "rows"):
            rows = [dict(item) for item in bounded.rows(named=True)]
        else:
            raise TypeError("polars-bio returned an unsupported table type")
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"polars-bio overlap failed: {exc}") from exc
    if len(rows) > limit:
        raise ConfigurationError(
            f"polars-bio overlap exceeds max_rows={limit}; increase the explicit limit."
        )
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "polars-bio",
            "operation": "overlap",
            "source": str(source_path),
            "other": str(other_path),
            "distinct": distinct,
            "rows": rows,
            "row_count": len(rows),
            "provenance": {
                "package": "polars-bio",
                "license": "Apache-2.0",
                "citation": "10.1093/bioinformatics/btaf640",
                "coordinate_system": "determined by polars-bio input metadata",
            },
        },
        verbose=verbose,
    )


def summarize_gtars_regions(
    source: str | Path,
    *,
    max_regions: int = 50_000,
    preview_regions: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a local BED file through the optional gtars RegionSet API.

    This adapter keeps the Rust-backed interval package optional and exposes a
    bounded, read-only inventory.  It does not accept URLs, construct binary
    caches, or expose arbitrary interval expressions through the JSON boundary.
    """

    path = input_path(source)
    limit = _bounded_limit(max_regions, "max_regions")
    preview_limit = _bounded_limit(preview_regions, "preview_regions", maximum=1_000)
    gtars_models = require("gtars.models", "bio-gtars")
    region_set_class = getattr(gtars_models, "RegionSet", None)
    if not callable(region_set_class):
        raise ConfigurationError("The installed gtars package lacks models.RegionSet.")
    try:
        region_set = region_set_class(str(path))
        regions = 0
        total_interval_bases = 0
        chromosomes: Counter[str] = Counter()
        preview: list[dict[str, Any]] = []
        for region in region_set:
            if regions >= limit:
                raise ConfigurationError(
                    f"gtars BED input exceeds max_regions={limit}; increase the explicit limit."
                )
            chromosome = str(getattr(region, "chr", getattr(region, "chrom", "")))
            start = getattr(region, "start", None)
            end = getattr(region, "end", None)
            if not chromosome or start is None or end is None:
                raise ConfigurationError("gtars returned a region without chr/start/end.")
            start_int = int(start)
            end_int = int(end)
            if end_int < start_int:
                raise ConfigurationError("gtars returned a region with end before start.")
            chromosomes[chromosome] += 1
            total_interval_bases += end_int - start_int
            if len(preview) < preview_limit:
                preview.append(
                    {
                        "chrom": chromosome,
                        "start": start_int,
                        "end": end_int,
                        "rest": getattr(region, "rest", None),
                    }
                )
            regions += 1
        covered_bases = total_interval_bases
        reduce = getattr(region_set, "reduce", None)
        if callable(reduce):
            reduced = reduce()
            covered_bases = sum(
                int(region.end) - int(region.start)
                for region in reduced
            )
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"gtars BED region parsing failed: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "gtars",
            "operation": "summarize_regions",
            "source": str(path),
            "regions": regions,
            "total_interval_bases": total_interval_bases,
            "covered_bases": covered_bases,
            "chromosomes": dict(chromosomes),
            "preview": preview,
            "provenance": {
                "package": "gtars",
                "license": "BSD-2-Clause",
                "citation": (
                    "gtars documentation and databio contributors; "
                    "no DOI supplied by the survey"
                ),
                "coordinate_system": "BED-style zero-based, half-open intervals",
                "backend_scope": "read-only local RegionSet inventory",
            },
        },
        verbose=verbose,
    )


def simulate_dnacauldron_assembly(
    parts: Mapping[str, str],
    *,
    method: str = "gibson",
    enzyme: str = "auto",
    max_constructs: int = 20,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate a bounded DNA assembly with the optional DNA Cauldron backend.

    ``parts`` maps stable part IDs to DNA strings.  Only the documented Gibson
    and Type IIS assembly classes are constructible through this JSON-safe
    boundary.  The backend predicts constructs and reports its own warnings or
    flaws; this adapter does not claim that a wet-lab assembly will succeed.
    """

    if not isinstance(parts, Mapping) or not 2 <= len(parts) <= 200:
        raise ConfigurationError("parts must contain between 2 and 200 named sequences.")
    normalized_parts: dict[str, str] = {}
    total_bases = 0
    for name, value in parts.items():
        part_name = _bounded_string(name, "part name", maximum=256)
        sequence = sequence_text(_bounded_string(value, f"parts[{part_name!r}]", maximum=200_000))
        if part_name in normalized_parts:
            raise ConfigurationError(f"Duplicate assembly part name: {part_name!r}.")
        normalized_parts[part_name] = sequence
        total_bases += len(sequence)
    if total_bases > 2_000_000:
        raise ConfigurationError("Total assembly input exceeds 2,000,000 bases.")
    if method not in {"gibson", "type2s"}:
        raise ConfigurationError("method must be 'gibson' or 'type2s'.")
    if not isinstance(enzyme, str) or not enzyme or len(enzyme) > 64:
        raise ConfigurationError("enzyme must be a non-empty string of at most 64 characters.")
    limit = _bounded_limit(max_constructs, "max_constructs", maximum=100)
    dnacauldron = require("dnacauldron", "bio-dnacauldron")
    repository = dnacauldron.SequenceRepository()
    for name, sequence in normalized_parts.items():
        repository.add_record((name, sequence))
    if method == "gibson":
        assembly_class = getattr(dnacauldron, "GibsonAssembly", None)
        kwargs = {
            "parts": list(normalized_parts),
            "max_constructs": limit,
            "expected_constructs": "any_number",
            "expect_no_unused_parts": False,
        }
    else:
        assembly_class = getattr(dnacauldron, "Type2sRestrictionAssembly", None)
        kwargs = {
            "parts": list(normalized_parts),
            "enzyme": enzyme,
            "max_constructs": limit,
            "expected_constructs": "any_number",
            "expect_no_unused_parts": False,
        }
    if not callable(assembly_class):
        raise ConfigurationError(
            f"The installed DNA Cauldron package lacks the {method} assembly class."
        )
    try:
        simulation = assembly_class(**kwargs).simulate(sequence_repository=repository)
        constructs: list[dict[str, Any]] = []
        for record in list(getattr(simulation, "construct_records", ()))[:limit]:
            sequence = str(getattr(record, "seq", record))
            constructs.append(
                {
                    "id": str(getattr(record, "id", f"construct_{len(constructs) + 1}")),
                    "length": len(sequence),
                    "sequence": sequence,
                }
            )
        warnings = [
            str(getattr(item, "message", item)) for item in getattr(simulation, "warnings", ())
        ]
        errors = [
            str(getattr(item, "message", item)) for item in getattr(simulation, "errors", ())
        ]
    except (AttributeError, OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"DNA Cauldron assembly simulation failed: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "DnaCauldron",
            "operation": "simulate_assembly",
            "method": method,
            "parts": list(normalized_parts),
            "constructs": constructs,
            "warnings": warnings,
            "errors": errors,
            "provenance": {
                "package": "dnacauldron",
                "license": "MIT",
                "citation": "DnaCauldron documentation; no DOI supplied by the survey",
                "backend_scope": "prediction only; no wet-lab result is inferred",
            },
        },
        verbose=verbose,
    )


def _sequence_groups(value: object) -> dict[str, list[str]]:
    if not isinstance(value, Mapping) or not value:
        raise ConfigurationError(
            "sequences must be a non-empty mapping of group names to sequences."
        )
    result: dict[str, list[str]] = {}
    for group, values in value.items():
        name = _bounded_string(group, "sequence group name", maximum=256)
        if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
            raise ConfigurationError(f"sequences[{name!r}] must be a non-empty sequence list.")
        records = [
            _bounded_string(item, f"sequences[{name!r}] item", maximum=100_000) for item in values
        ]
        result[name] = records
    return result


def _metric_specs(value: object) -> list[dict[str, Any]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError("metrics must be a non-empty list of metric specifications.")
    specs: list[dict[str, Any]] = []
    for item in value:
        if isinstance(item, str):
            item = {"name": item}
        if not isinstance(item, Mapping) or not isinstance(item.get("name"), str):
            raise ConfigurationError("Each metric must be a name or an object with a string name.")
        name = item["name"]
        kwargs = item.get("kwargs", {})
        if not isinstance(kwargs, Mapping) or any(not isinstance(key, str) for key in kwargs):
            raise ConfigurationError(f"Metric {name!r} kwargs must be a string-keyed mapping.")
        if any(callable(value) for value in kwargs.values()):
            raise ConfigurationError("Metric kwargs cannot contain Python callbacks.")
        specs.append({"name": name, "kwargs": dict(kwargs)})
    return specs


def evaluate_seqme(
    sequences: Mapping[str, Sequence[str]],
    metrics: Sequence[str | Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate sequence groups with allow-listed ``seqme.metrics`` classes.

    ``metrics`` entries are names such as ``"Uniqueness"`` or mappings such
    as ``{"name": "Novelty", "kwargs": {"reference": [...]}}``.  Metrics
    needing a live predictor/embedder must be supplied by a native Python
    integration instead of crossing the JSON/MCP boundary.
    """

    groups = _sequence_groups(sequences)
    specs = _metric_specs(metrics)
    seqme = require("seqme", "bio-seqme")
    metric_module = require("seqme.metrics", "bio-seqme")
    built: list[Any] = []
    names: list[str] = []
    for spec in specs:
        name = spec["name"]
        metric_class = getattr(metric_module, name, None)
        if not inspect.isclass(metric_class) or not callable(metric_class):
            raise ConfigurationError(f"Unknown or non-constructible seqme metric: {name!r}.")
        kwargs = dict(spec["kwargs"])
        if any(key in kwargs for key in ("predictor", "embedder", "model", "cache")):
            raise ConfigurationError(
                f"seqme metric {name!r} requires a live model object; use a Python-only adapter."
            )
        try:
            metric = metric_class(**kwargs)
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(f"Could not construct seqme metric {name!r}: {exc}") from exc
        built.append(metric)
        names.append(name)
    try:
        frame = seqme.evaluate(groups, built)
        table = frame.to_dict(orient="index") if hasattr(frame, "to_dict") else frame
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"seqme evaluation failed: {exc}") from exc
    return {
        "status": "completed",
        "tool": "seqme",
        "metrics": names,
        "groups": list(groups),
        "result": table,
        "provenance": {
            "package": "seqme",
            "license": "BSD-3-Clause",
            "citation": "arXiv:2511.04239",
        },
    }


def _numeric_matrix(value: object, name: str) -> list[list[float]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty rectangular numeric matrix.")
    rows: list[list[float]] = []
    width: int | None = None
    for row in value:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or not row:
            raise ConfigurationError(f"{name} must contain non-empty numeric rows.")
        converted: list[float] = []
        for item in row:
            if isinstance(item, bool) or not isinstance(item, (int, float)):
                raise ConfigurationError(f"{name} must contain only finite numbers.")
            number = float(item)
            if number != number or number in {float("inf"), float("-inf")}:
                raise ConfigurationError(f"{name} must contain only finite numbers.")
            converted.append(number)
        if width is None:
            width = len(converted)
        elif len(converted) != width:
            raise ConfigurationError(f"{name} must be rectangular.")
        rows.append(converted)
    return rows


def _ids(value: object, expected: int, name: str = "ids") -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of unique strings.")
    result = [_bounded_string(item, f"{name} item", maximum=256) for item in value]
    if len(result) != expected or len(set(result)) != len(result):
        raise ConfigurationError(f"{name} must contain exactly {expected} unique strings.")
    return result


def _labels(value: object, expected: int, name: str = "labels") -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of strings.")
    result = [_bounded_string(item, f"{name} item", maximum=256) for item in value]
    if len(result) != expected:
        raise ConfigurationError(f"{name} must contain exactly {expected} strings.")
    return result


_SKBIO_METRICS = frozenset({"braycurtis", "euclidean", "cityblock", "cosine", "jaccard"})


def beta_diversity_scikit_bio(
    counts: Sequence[Sequence[float]],
    ids: Sequence[str],
    *,
    metric: str = "braycurtis",
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a bounded scikit-bio beta-diversity distance matrix."""

    matrix = _numeric_matrix(counts, "counts")
    labels = _ids(ids, len(matrix))
    if not isinstance(metric, str) or metric not in _SKBIO_METRICS:
        raise ConfigurationError(
            f"metric must be one of {sorted(_SKBIO_METRICS)} for the JSON adapter."
        )
    if any(value < 0 for row in matrix for value in row):
        raise ConfigurationError("counts cannot contain negative values.")
    beta_diversity = require("skbio.diversity", "bio-skibio").beta_diversity
    try:
        distance = beta_diversity(metric, matrix, ids=labels)
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"scikit-bio beta-diversity failed: {exc}") from exc
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "scikit-bio",
            "operation": "beta_diversity",
            "metric": metric,
            "ids": labels,
            "distance_matrix": distance.data.tolist(),
            "provenance": {
                "package": "scikit-bio",
                "license": "BSD-3-Clause",
                "citation": "10.1038/s41592-025-02981-z",
            },
        },
        verbose=verbose,
    )


def permanova_scikit_bio(
    distance_matrix: Sequence[Sequence[float]],
    ids: Sequence[str],
    grouping: Sequence[str],
    *,
    permutations: int = 999,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run scikit-bio PERMANOVA over a bounded JSON distance matrix."""

    matrix = _numeric_matrix(distance_matrix, "distance_matrix")
    labels = _ids(ids, len(matrix))
    if any(len(row) != len(matrix) for row in matrix):
        raise ConfigurationError("distance_matrix must be square.")
    groups = _labels(grouping, len(matrix), "grouping")
    if (
        isinstance(permutations, bool)
        or not isinstance(permutations, int)
        or not 0 <= permutations <= 9_999
    ):
        raise ConfigurationError("permutations must be an integer in [0, 9999].")
    distance_module = require("skbio.stats.distance", "bio-skibio")
    try:
        distance = distance_module.DistanceMatrix(matrix, ids=labels)
        result = distance_module.permanova(distance, groups, permutations=permutations)
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"scikit-bio PERMANOVA failed: {exc}") from exc
    values = result.to_dict() if hasattr(result, "to_dict") else dict(result)
    return serialize_for_agent(
        {
            "status": "completed",
            "tool": "scikit-bio",
            "operation": "permanova",
            "permutations": permutations,
            "result": values,
            "provenance": {
                "package": "scikit-bio",
                "license": "BSD-3-Clause",
                "citation": "10.1038/s41592-025-02981-z",
            },
        },
        verbose=verbose,
    )


_DNACHISEL_SPECS = frozenset(
    {
        "AvoidPattern",
        "AvoidChanges",
        "AvoidHairpins",
        "CodonOptimize",
        "EnforceGCContent",
        "EnforceTranslation",
        "EnforceChanges",
        "UniquifyAllKmers",
    }
)


def _dnachisel_specs(value: object, field: str) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{field} must be a list of specification objects.")
    result = []
    for item in value:
        if not isinstance(item, Mapping) or not isinstance(item.get("name"), str):
            raise ConfigurationError(f"Each {field} entry needs a string name.")
        name = item["name"]
        if name not in _DNACHISEL_SPECS:
            raise ConfigurationError(
                f"Unsupported DNA Chisel specification {name!r}; use an allow-listed built-in."
            )
        args = item.get("args", [])
        kwargs = item.get("kwargs", {})
        if isinstance(args, (str, bytes)) or not isinstance(args, Sequence):
            raise ConfigurationError(f"{field}.{name}.args must be a JSON list.")
        if not isinstance(kwargs, Mapping) or any(not isinstance(key, str) for key in kwargs):
            raise ConfigurationError(f"{field}.{name}.kwargs must be a string-keyed mapping.")
        kwargs = dict(kwargs)
        if "location" in kwargs and isinstance(kwargs["location"], list):
            kwargs["location"] = tuple(kwargs["location"])
        result.append((name, list(args), kwargs))
    return result


def optimize_with_dnachisel(
    sequence: str,
    *,
    constraints: Sequence[Mapping[str, Any]] = (),
    objectives: Sequence[Mapping[str, Any]] = (),
    circular: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a bounded DNA Chisel problem through an explicit optional backend.

    Only documented built-in specification classes are constructible from the
    JSON contract.  DNA Chisel remains an optional MIT dependency; this adapter
    does not copy its source or silently download codon tables or models.
    """

    sequence = sequence_text(_bounded_string(sequence, "sequence", maximum=2_000_000))
    if not isinstance(circular, bool):
        raise ConfigurationError("circular must be a boolean.")
    dnachisel = require("dnachisel", "bio-dnachisel")
    problem_class = getattr(
        dnachisel,
        "CircularDnaOptimizationProblem" if circular else "DnaOptimizationProblem",
        None,
    )
    if problem_class is None:
        raise ConfigurationError(
            "The installed DNA Chisel package lacks its optimization problem class."
        )
    built_constraints = []
    built_objectives = []
    for name, args, kwargs in _dnachisel_specs(constraints, "constraints"):
        specification = getattr(dnachisel, name, None)
        if specification is None:
            raise ConfigurationError(f"The installed DNA Chisel package lacks {name}.")
        built_constraints.append(specification(*args, **kwargs))
    for name, args, kwargs in _dnachisel_specs(objectives, "objectives"):
        specification = getattr(dnachisel, name, None)
        if specification is None:
            raise ConfigurationError(f"The installed DNA Chisel package lacks {name}.")
        built_objectives.append(specification(*args, **kwargs))
    problem = problem_class(
        sequence=sequence,
        constraints=built_constraints,
        objectives=built_objectives,
        logger=None,
    )
    problem.resolve_constraints()
    if built_objectives:
        problem.optimize()
    return {
        "status": "completed",
        "tool": "DNA Chisel",
        "sequence": str(problem.sequence),
        "constraints": [name for name, _, _ in _dnachisel_specs(constraints, "constraints")],
        "objectives": [name for name, _, _ in _dnachisel_specs(objectives, "objectives")],
        "provenance": {
            "package": "dnachisel",
            "license": "MIT",
            "citation": "10.1093/bioinformatics/btaa558",
        },
    }


__all__ = [
    "align_biotite",
    "align_edlib",
    "beta_diversity_scikit_bio",
    "compare_sourmash_sequences",
    "diff_geneblocks_sequences",
    "evaluate_seqme",
    "fetch_pyfaidx_region",
    "intersect_pybedtools_bed",
    "map_mappy_sequences",
    "merge_pyranges_bed",
    "optimize_with_dnachisel",
    "overlap_polars_bio_bed",
    "permanova_scikit_bio",
    "query_gffutils_database",
    "simulate_dnacauldron_assembly",
    "summarize_dna_features_viewer",
    "summarize_dnaio_fastx",
    "summarize_gtars_regions",
    "summarize_pybigwig_region",
    "summarize_pyfastx",
    "summarize_pyranges_bed",
    "summarize_sbol3_document",
    "summarize_scadnano_design",
    "summarize_vcf_cyvcf2",
    "summarize_vcf_vcfpy",
]

_register_public_verbose(globals(), tuple(__all__))
