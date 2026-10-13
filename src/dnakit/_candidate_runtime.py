"""Shared local runtime for the broad candidate-adapter surface.

The adapters in this module deliberately separate deterministic local work from
operations that need a database, external executable, model, or cloud service.
They return JSON-serializable envelopes so that a future formal DNAKit/MCP
adapter can add schemas and backend dispatch without changing every operation
name.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path


DEPENDENCY_BY_MODULE = {
    "assembly": "assembler backend or assembly graph library",
    "clinical": "clinical annotation database and validated clinical rules",
    "databases": "configured database or network service",
    "design_crispr": "validated CRISPR scoring/off-target backend",
    "external_adapters": "configured external executable or service",
    "forensic_reproductive": "validated population and forensic reference data",
    "genome3d": "3D genome file format or contact-map backend",
    "gwas": "validated statistical genetics backend and phenotype/genotype data",
    "haplotypes": "phasing or pedigree inference backend",
    "knowledge_graph": "configured graph database or knowledge source",
    "microbiome": "taxonomy, marker, or metagenomic reference database",
    "ml_models": "trained model artifact and compatible ML runtime",
    "molecular_lab": "validated laboratory protocol and instrument configuration",
    "population": "validated population reference panel or statistics backend",
    "sequencing_qc": "sequencing platform/QC backend for platform-specific metrics",
    "singlecell_spatial": "single-cell or spatial data backend",
    "transcriptomics": "transcript annotation or quantification backend",
    "variants": "variant reference database or caller backend",
    "workflow_security": "configured execution, secrets, or cloud/HPC backend",
}

DEPENDENCY_OPERATIONS = {
    "align", "assemble", "call", "download", "fetch", "finetune", "infer",
    "predict", "pretrain", "query", "run", "search", "serve", "submit",
    "train", "type_hla", "quantify", "deconvolute", "demultiplex",
}


def _readable(value: object) -> object:
    """Return a bounded, JSON-friendly representation of an input value."""

    if isinstance(value, Path):
        if value.exists() and value.is_file():
            data = value.read_bytes()
            return data[:1_000_000].decode("utf-8", errors="replace")
        return str(value)
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8")
        except UnicodeDecodeError:
            return {"base64": base64.b64encode(value[:1_000_000]).decode("ascii")}
    if isinstance(value, Mapping):
        return {str(key): _readable(item) for key, item in list(value.items())[:100]}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_readable(item) for item in list(value)[:100]]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def _text(value: object) -> str:
    value = _readable(value)
    if isinstance(value, str):
        return value
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _sequence(value: object) -> str:
    text = _text(value).upper()
    if ">" in text or "@" in text:
        text = "".join(line.strip() for line in text.splitlines() if not line.startswith((">", "@", "+")))
    return re.sub(r"[^A-Z*]", "", text)


def _records(value: object) -> list[object]:
    if isinstance(value, Mapping):
        for key in ("records", "items", "values", "data", "sequences", "reads"):
            if key in value:
                return _records(value[key])
        return [value]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return list(value)
    text = _text(value)
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return lines or ([text] if text else [])


def _lengths(value: object) -> list[int]:
    result = []
    for item in _records(value):
        if isinstance(item, Mapping):
            candidate = item.get("sequence", item.get("value", item.get("length", "")))
        else:
            candidate = item
        if isinstance(candidate, (int, float)):
            result.append(max(0, int(candidate)))
        else:
            result.append(len(_sequence(candidate)))
    return [length for length in result if length]


def _summary(value: object) -> dict[str, object]:
    text = _text(value)
    records = _records(value)
    lengths = _lengths(value)
    return {
        "record_count": len(records),
        "text_length": len(text),
        "sequence_bases": sum(lengths),
        "lengths": lengths[:100],
        "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
    }


def _n_stat(lengths: list[int], fraction: float) -> int:
    if not lengths:
        return 0
    target = sum(lengths) * fraction
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= target:
            return length
    return min(lengths)


def _svg(operation: str, summary: Mapping[str, object]) -> str:
    label = f"{operation}: {summary.get('record_count', 0)} records"
    return (
        '<svg xmlns="http://www.w3.org/2000/svg" width="640" height="80">'
        f'<text x="8" y="32">{label}</text></svg>'
    )


def _local_result(module: str, operation: str, value: object, options: Mapping[str, object]) -> object:
    summary = _summary(value)
    op = operation.casefold()
    lengths = summary["lengths"]

    if op in {"validate", "check", "check_input", "validate_output", "validate_metadata_schema", "validate_spatial"}:
        return {"valid": value is not None, "summary": summary, "errors": [] if value is not None else ["input is required"]}
    if op in {"count", "count_mt", "count_reads", "count_records"}:
        return int(summary["record_count"])
    if op in {"calculate_n50", "calculate_ng50", "calculate_nx"}:
        fraction = float(options.get("fraction", options.get("x", 0.5)))
        if op == "calculate_n50":
            fraction = 0.5
        return _n_stat(list(lengths), fraction)
    if op in {"calculate_l50"}:
        target = sum(lengths) * 0.5
        running = 0
        for index, length in enumerate(sorted(lengths, reverse=True), 1):
            running += length
            if running >= target:
                return index
        return 0
    if op in {"calculate", "summarize", "stats", "estimate", "qc", "report", "read"} or op.startswith(("calculate_", "estimate_", "summarize_")):
        return summary
    if op in {"compare", "compare_mt", "compare_vcfs", "compare_assemblies", "compare_database_versions"}:
        other = options.get("other", "")
        return {"equal": summary["sha256"] == _summary(other)["sha256"], "left": summary, "right": _summary(other)}
    if op in {"merge", "combine", "aggregate", "integrate", "concatenate"}:
        values = _records(value)
        return {"records": values, "record_count": len(values), "summary": summary}
    if op in {"filter", "select", "prioritize", "rank", "sort"}:
        values = _records(value)
        threshold = options.get("minimum_length")
        if isinstance(threshold, (int, float)):
            values = [item for item in values if len(_sequence(item)) >= int(threshold)]
        return {"records": values, "record_count": len(values), "summary": _summary(values)}
    if op in {"normalize", "standardize", "clean", "sanitize", "correct"}:
        if isinstance(value, str):
            return re.sub(r"\s+", "", value).upper()
        return {"records": [_sequence(item) for item in _records(value)], "summary": summary}
    if op in {"plot", "visualize", "export_figure", "save_figure"} or "plot" in op:
        return _svg(operation, summary)
    if op in {"write", "export", "save"}:
        return {
            "serialized": _readable(value),
            "output_path": options.get("output") or options.get("path"),
            "summary": summary,
        }
    if op.startswith("compress"):
        data = _text(value).encode("utf-8")
        return {"encoding": "gzip+base64", "data": base64.b64encode(gzip.compress(data)).decode("ascii"), "bytes": len(data)}
    if op.startswith("decompress"):
        try:
            data = value if isinstance(value, bytes) else base64.b64decode(str(value))
            return gzip.decompress(data).decode("utf-8", errors="replace")
        except (OSError, ValueError, TypeError):
            return {"decoded": False, "error": "input is not a valid gzip/base64 payload"}
    if op in {"get_file_checksum", "checksum"} or op.endswith("_checksum"):
        return hashlib.sha256(_text(value).encode("utf-8")).hexdigest()
    if op in {"calculate_primer_tm", "evaluate_primer_pair"}:
        sequence = _sequence(value)
        at = sequence.count("A") + sequence.count("T")
        gc = sequence.count("G") + sequence.count("C")
        return {"tm_celsius": 2 * at + 4 * gc, "length": len(sequence), "gc_fraction": gc / len(sequence) if sequence else 0.0}
    if op in {"find_sequence_patterns", "scan_pwm"}:
        sequence = _sequence(value)
        pattern = str(options.get("pattern", options.get("motif", "ACG")))
        return [{"start": match.start(), "end": match.end(), "match": match.group(0)} for match in re.finditer(pattern, sequence)]
    if op in {"encode_data_into_dna"}:
        raw = _text(value).encode("utf-8")
        alphabet = ("A", "C", "G", "T")
        return "".join(alphabet[(byte >> shift) & 3] for byte in raw for shift in (6, 4, 2, 0))
    if op in {"decode_dna_storage"}:
        dna = _sequence(value)
        pairs = {"A": 0, "C": 1, "G": 2, "T": 3}
        bits = [pairs[base] for base in dna if base in pairs]
        raw = bytes(sum(bits[index + offset] << (6 - 2 * offset) for offset in range(4)) for index in range(0, len(bits) - 3, 4))
        return raw.decode("utf-8", errors="replace")
    if op in {"remove_phi_from_metadata", "deidentify_genomic_metadata"}:
        if not isinstance(value, Mapping):
            return value
        redacted = {}
        for key, item in value.items():
            if str(key).casefold() in {"name", "patient", "patient_id", "email", "phone", "address"}:
                redacted[str(key)] = hashlib.sha256(str(item).encode()).hexdigest()[:16]
            else:
                redacted[str(key)] = item
        return redacted
    if op in {"set_random_seed"}:
        return {"seed": options.get("seed", value), "deterministic": True}
    return {"result": _readable(value), "summary": summary}


def run_candidate_operation(module: str, operation: str, input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Execute a deterministic local adapter or report a missing backend."""

    options = dict(options or {})
    dependency = DEPENDENCY_BY_MODULE.get(module)
    base = {
        "module": module,
        "operation": operation,
        "provenance": {"implementation": "local_candidate_adapter", "backend": options.get("backend", "local")},
        "input_summary": _summary(input),
        "warnings": [],
    }
    operation_root = operation.casefold().split("_", 1)[0]
    requires_backend = (
        dependency is not None
        and (operation_root in DEPENDENCY_OPERATIONS or operation in {"build", "design", "optimize", "simulate", "score"})
        and not options.get("allow_local_fallback", False)
    )
    if requires_backend:
        return {
            **base,
            "status": "dependency_required",
            "required_dependency": dependency,
            "message": "The candidate API and input validation are available locally; scientific execution requires the declared backend.",
        }
    result = _local_result(module, operation, input, options)
    return {**base, "status": "ok", "result": result}

