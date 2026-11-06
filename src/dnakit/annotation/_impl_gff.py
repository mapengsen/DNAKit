"""Dependency-free GFF3/GTF parsing and gene-structure logic behind the ``dnakit.annotation`` functions.

Coordinates are the 1-based, closed intervals of GFF3/GTF.  Nothing here reads a reference genome or a database.
"""

from __future__ import annotations

import gzip
import itertools
import os
import re
import tempfile
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from urllib.parse import unquote

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.sequence import CODON_TABLE, normalize_dna_sequence, reverse_complement_dna

MAX_TEXT_BYTES = 50_000_000
MAX_FEATURES = 1_000_000

GENE_TYPES = frozenset(
    {"gene", "pseudogene", "ncrna_gene", "transposable_element_gene", "protein_coding_gene"}
)
EXON_TYPES = frozenset({"exon"})
CDS_TYPES = frozenset({"cds"})
FIVE_UTR_TYPES = frozenset({"five_prime_utr", "5utr", "5'utr", "utr5", "5_utr"})
THREE_UTR_TYPES = frozenset({"three_prime_utr", "3utr", "3'utr", "utr3", "3_utr"})
UTR_TYPES = frozenset({"utr"}) | FIVE_UTR_TYPES | THREE_UTR_TYPES
CODON_TYPES = frozenset({"start_codon", "stop_codon"})
MIRNA_TYPES = frozenset({"mirna", "pre_mirna", "mirna_primary_transcript", "primary_mirna"})
_CHILD_TYPES = EXON_TYPES | CDS_TYPES | UTR_TYPES | CODON_TYPES
_BIOTYPE_KEYS = ("biotype", "gene_biotype", "transcript_biotype", "gene_type", "transcript_type")
_GTF_ATTRIBUTE = re.compile(r'\s*([^\s;"]+)\s+(?:"((?:[^"\\]|\\.)*)"|([^;"]*?))\s*(?:;|$)')
_GTF_LIKE = re.compile(r'(?:^|;)\s*[A-Za-z_][\w.]*\s+"[^"]*"\s*(?:;|$)')


def is_transcript_type(feature_type: str) -> bool:
    lowered = feature_type.lower()
    return (
        lowered
        in {
            "transcript",
            "primary_transcript",
            "pre_mirna",
            "mirna",
            "lncrna",
            "lnc_rna",
            "mirna_primary_transcript",
        }
        or lowered.endswith("rna")
        or lowered.endswith("_rna")
        or lowered.endswith("_transcript")
    )


# --- reading ------------------------------------------------------------------------------------


def is_file(path: Path) -> bool:
    try:
        return path.is_file()
    except (OSError, ValueError):  # e.g. a very long sequence string used as a path
        return False


def read_text(value: object) -> str:
    """Text of an annotation given as text, ``Path``, path string (``.gz`` supported) or mapping with ``text``/``path``."""

    if isinstance(value, Mapping):
        for key in ("text", "gff", "gtf", "gff3"):
            if isinstance(value.get(key), str):
                return read_text(value[key])
        if "path" in value:
            return read_text(Path(str(value["path"])))
        raise InputFormatError(
            "The annotation mapping needs a 'text' or 'path' key.", code="INVALID_ANNOTATION_INPUT"
        )
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    is_path = isinstance(value, os.PathLike) or (
        isinstance(value, str) and "\n" not in value and "\t" not in value and value.strip() != ""
    )
    if is_path:
        path = Path(value)  # type: ignore[arg-type]
        if not is_file(path):
            raise InputFormatError(f"Annotation file not found: {path}", code="FILE_NOT_FOUND")
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rb") as stream:  # type: ignore[operator]
            data = stream.read(MAX_TEXT_BYTES + 1)
        if len(data) > MAX_TEXT_BYTES:
            raise InputFormatError(
                f"Annotation exceeds {MAX_TEXT_BYTES} bytes.", code="FILE_TOO_LARGE"
            )
        return data.decode("utf-8", errors="replace")
    if isinstance(value, str):
        if len(value) > MAX_TEXT_BYTES:
            raise InputFormatError(
                f"Annotation exceeds {MAX_TEXT_BYTES} characters.", code="FILE_TOO_LARGE"
            )
        return value
    raise InputFormatError(
        f"Unsupported annotation input type {type(value).__name__}.",
        code="INVALID_ANNOTATION_INPUT",
        hint="Pass GFF3/GTF text, a file path, or a list of feature mappings.",
    )


def detect_format(text: str) -> str:
    """``gff3`` or ``gtf``: the ``##gff-version 3`` directive wins, then the shape of the attribute column."""

    checked = 0
    for line in text.splitlines():
        if line.startswith("##gff-version"):
            return "gff3"
        if not line.strip() or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) >= 9:
            attributes = columns[8]
            if "=" in attributes and not _GTF_LIKE.search(attributes):
                return "gff3"
            if _GTF_LIKE.search(attributes):
                return "gtf"
        checked += 1
        if checked >= 50:
            break
    return "gff3"


def parse_attributes(raw: str, fmt: str) -> dict[str, object]:
    attributes: dict[str, object] = {}
    raw = raw.strip()
    if raw in {"", "."}:
        return attributes
    if fmt == "gff3":
        for part in raw.split(";"):
            part = part.strip()
            if not part:
                continue
            if "=" not in part:
                raise ValueError(f"attribute {part!r} is not key=value")
            key, _, value = part.partition("=")
            values = [unquote(item) for item in value.split(",")]
            attributes[unquote(key.strip())] = values if len(values) > 1 else values[0]
        return attributes
    position = 0
    while position < len(raw):
        match = _GTF_ATTRIBUTE.match(raw, position)
        if match is None or match.end() == position:
            if raw[position:].strip(" ;") == "":
                break
            raise ValueError(f"cannot parse GTF attributes near {raw[position : position + 30]!r}")
        key = match.group(1)
        value = match.group(2) if match.group(2) is not None else (match.group(3) or "").strip()
        value = value.replace('\\"', '"')
        if key in attributes:
            current = attributes[key]
            attributes[key] = [*current, value] if isinstance(current, list) else [current, value]
        else:
            attributes[key] = value
        position = match.end()
    return attributes


def parse_lines(text: str, fmt: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Parse feature lines; malformed lines become issues (and are skipped).  Returns ``(features, issues)``."""

    features: list[dict[str, object]] = []
    issues: list[dict[str, object]] = []

    def issue(line: int, severity: str, code: str, message: str) -> None:
        issues.append({"line": line, "severity": severity, "code": code, "message": message})

    for number, line in enumerate(text.splitlines(), start=1):
        if line.startswith("##FASTA"):
            break
        if not line.strip() or line.startswith("#"):
            continue
        columns = line.rstrip("\r").split("\t")
        if len(columns) != 9:
            issue(
                number,
                "error",
                "COLUMN_COUNT",
                f"expected 9 tab-separated columns, found {len(columns)}",
            )
            continue
        (
            seqid,
            source,
            ftype,
            start_text,
            end_text,
            score_text,
            strand,
            phase_text,
            attribute_text,
        ) = columns
        try:
            start, end = int(start_text), int(end_text)
        except ValueError:
            issue(
                number,
                "error",
                "COORDINATE_NOT_INTEGER",
                f"start/end must be integers: {start_text!r}, {end_text!r}",
            )
            continue
        bad = False
        if not seqid.strip():
            issue(number, "error", "EMPTY_SEQID", "seqid is empty")
            bad = True
        if not ftype.strip() or ftype == ".":
            issue(number, "error", "EMPTY_TYPE", "feature type is empty")
            bad = True
        if start < 1:
            issue(number, "error", "COORDINATE_RANGE", f"start must be >= 1, got {start}")
            bad = True
        if end < start:
            issue(
                number, "error", "COORDINATE_ORDER", f"end ({end}) is smaller than start ({start})"
            )
            bad = True
        if strand not in {"+", "-", ".", "?"}:
            issue(
                number, "error", "INVALID_STRAND", f"strand must be one of + - . ?, got {strand!r}"
            )
            bad = True
        score: float | None = None
        if score_text != ".":
            try:
                score = float(score_text)
            except ValueError:
                issue(
                    number,
                    "error",
                    "INVALID_SCORE",
                    f"score must be a number or '.', got {score_text!r}",
                )
                bad = True
        phase: int | None = None
        if phase_text != ".":
            if phase_text in {"0", "1", "2"}:
                phase = int(phase_text)
            else:
                issue(
                    number,
                    "error",
                    "INVALID_PHASE",
                    f"phase must be 0, 1, 2 or '.', got {phase_text!r}",
                )
                bad = True
        if ftype.lower() == "cds" and phase_text == ".":
            issue(number, "error", "CDS_WITHOUT_PHASE", "CDS features require a phase (0, 1 or 2)")
        try:
            attributes = parse_attributes(attribute_text, fmt)
        except ValueError as exc:
            issue(number, "error", "INVALID_ATTRIBUTES", str(exc))
            bad = True
            attributes = {}
        if bad:
            continue
        if len(features) >= MAX_FEATURES:
            raise InputFormatError(
                f"More than {MAX_FEATURES} features; split the annotation.",
                code="TOO_MANY_FEATURES",
            )
        features.append(
            {
                "line": number,
                "seqid": seqid,
                "source": source,
                "type": ftype,
                "start": start,
                "end": end,
                "score": score,
                "strand": strand,
                "phase": phase,
                "attributes": attributes,
            }
        )
    return features, issues


def load_features(
    value: object, options: Mapping[str, object] | None = None, *, strict: bool = True
) -> tuple[str, list[dict[str, object]]]:
    """Features from annotation text/path or from a list of feature mappings (as returned by ``parse``)."""

    settings = options or {}
    if isinstance(value, Mapping) and isinstance(value.get("features"), Sequence):
        value = value["features"]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return str(settings.get("format", "gff3")), [
            normalize_feature(item, index) for index, item in enumerate(value)
        ]
    text = read_text(value)
    fmt = str(settings.get("format", "auto"))
    if fmt == "auto":
        fmt = detect_format(text)
    if fmt not in {"gff3", "gtf", "gff"}:
        raise ConfigurationError("options.format must be auto, gff3 or gtf.", code="INVALID_OPTION")
    fmt = "gff3" if fmt == "gff" else fmt
    features, issues = parse_lines(text, fmt)
    errors = [item for item in issues if item["severity"] == "error"]
    if strict and errors:
        first = errors[0]
        raise InputFormatError(
            f"Malformed {fmt.upper()} at line {first['line']}: {first['message']}",
            code="INVALID_ANNOTATION",
            context={"line": first["line"], "error_count": len(errors)},
            hint="Run annotation.validate for the complete list of problems, or pass strict=false.",
        )
    return fmt, features


def normalize_feature(item: object, index: int) -> dict[str, object]:
    if not isinstance(item, Mapping):
        raise InputFormatError(
            f"feature {index} must be a mapping.", code="INVALID_ANNOTATION_INPUT"
        )
    seqid = item.get("seqid", item.get("chrom", item.get("chromosome")))
    start, end = item.get("start"), item.get("end")
    if (
        not isinstance(seqid, str)
        or not seqid
        or any(isinstance(v, bool) or not isinstance(v, int) for v in (start, end))
    ):
        raise InputFormatError(
            f"feature {index} needs a seqid and integer start/end.", code="INVALID_ANNOTATION_INPUT"
        )
    if start < 1 or end < start:  # type: ignore[operator]
        raise InputFormatError(
            f"feature {index} has an invalid interval {start}-{end}.",
            code="INVALID_ANNOTATION_INPUT",
        )
    attributes = item.get("attributes", {})
    if not isinstance(attributes, Mapping):
        raise InputFormatError(
            f"feature {index}: attributes must be a mapping.", code="INVALID_ANNOTATION_INPUT"
        )
    return {
        "line": item.get("line", index + 1),
        "seqid": seqid,
        "source": str(item.get("source", ".")),
        "type": str(item.get("type", ".")),
        "start": start,
        "end": end,
        "score": item.get("score"),
        "strand": str(item.get("strand", ".")),
        "phase": item.get("phase"),
        "attributes": dict(attributes),
    }


# --- hierarchy ----------------------------------------------------------------------------------


def as_list(value: object) -> list[str]:
    if value is None:
        return []
    return [str(item) for item in value] if isinstance(value, list) else [str(value)]


def attribute(feature: Mapping[str, object], *names: str) -> str | None:
    attributes = feature["attributes"]
    for name in names:
        values = as_list(attributes.get(name))  # type: ignore[union-attr]
        if values:
            return values[0]
    return None


def biotypes(feature: Mapping[str, object]) -> list[str]:
    return [
        value.lower() for key in _BIOTYPE_KEYS for value in as_list(feature["attributes"].get(key))
    ]  # type: ignore[union-attr]


def feature_id(feature: Mapping[str, object], fmt: str) -> str | None:
    ftype = str(feature["type"]).lower()
    if fmt == "gff3":
        return attribute(feature, "ID")
    if ftype in GENE_TYPES:
        return attribute(feature, "gene_id")
    if is_transcript_type(ftype):
        return attribute(feature, "transcript_id")
    return None


def parents(feature: Mapping[str, object], fmt: str) -> list[str]:
    ftype = str(feature["type"]).lower()
    if fmt == "gff3":
        return as_list(feature["attributes"].get("Parent"))  # type: ignore[union-attr]
    if ftype in GENE_TYPES:
        return []
    if is_transcript_type(ftype):
        gene = attribute(feature, "gene_id")
        return [gene] if gene else []
    transcript = attribute(feature, "transcript_id")
    return [transcript] if transcript else []


def length(feature: Mapping[str, object]) -> int:
    return int(feature["end"]) - int(feature["start"]) + 1  # type: ignore[call-overload]


def build_models(
    features: Sequence[Mapping[str, object]], fmt: str
) -> dict[str, dict[str, object]]:
    """Group features into ``{gene_id: {"feature", "transcripts": {tid: {...}}}}``.

    Transcripts are features of a transcript-like type (mRNA, transcript, *RNA, ...) or, in GTF, the distinct
    ``transcript_id`` values of exon/CDS/UTR/codon lines.  Features without a gene use the pseudo gene ``""``.
    """

    by_id: dict[str, Mapping[str, object]] = {}
    for feature in features:
        identifier = feature_id(feature, fmt)
        if identifier is not None and identifier not in by_id:
            by_id[identifier] = feature
    genes: dict[str, dict[str, object]] = {}
    transcripts: dict[str, dict[str, object]] = {}

    def gene_entry(gene_id: str) -> dict[str, object]:
        return genes.setdefault(gene_id, {"feature": None, "transcripts": {}})

    def transcript_entry(transcript_id: str, gene_id: str) -> dict[str, object]:
        if transcript_id not in transcripts:
            entry: dict[str, object] = {
                "id": transcript_id,
                "gene_id": gene_id,
                "feature": None,
                "exons": [],
                "cds": [],
                "five_prime_utr": [],
                "three_prime_utr": [],
                "utr": [],
                "start_codon": [],
                "stop_codon": [],
            }
            transcripts[transcript_id] = entry
            gene_entry(gene_id)["transcripts"][transcript_id] = entry  # type: ignore[index]
        return transcripts[transcript_id]

    for feature in features:
        ftype = str(feature["type"]).lower()
        identifier = feature_id(feature, fmt)
        if ftype in GENE_TYPES and identifier is not None:
            gene_entry(identifier)["feature"] = feature
    for feature in features:
        ftype = str(feature["type"]).lower()
        if ftype in GENE_TYPES or not is_transcript_type(ftype):
            continue
        identifier = feature_id(feature, fmt)
        if identifier is None:
            continue
        gene_id = next(
            (p for p in parents(feature, fmt) if p in genes), (parents(feature, fmt) or [""])[0]
        )
        transcript_entry(identifier, gene_id)["feature"] = feature
    for feature in features:
        ftype = str(feature["type"]).lower()
        if ftype not in _CHILD_TYPES:
            continue
        for parent in parents(feature, fmt) or [""]:
            if parent == "" and fmt == "gff3":
                continue
            entry = transcripts.get(parent)
            if entry is None:
                if fmt == "gtf":
                    gene_id = attribute(feature, "gene_id") or ""
                    entry = transcript_entry(parent, gene_id)
                elif parent in by_id and str(by_id[parent]["type"]).lower() in GENE_TYPES:
                    entry = transcript_entry(
                        parent, parent
                    )  # exon attached directly to a gene (non-standard but common)
                else:
                    continue
            bucket = (
                "exons"
                if ftype in EXON_TYPES
                else "cds"
                if ftype in CDS_TYPES
                else "five_prime_utr"
                if ftype in FIVE_UTR_TYPES
                else "three_prime_utr"
                if ftype in THREE_UTR_TYPES
                else "utr"
                if ftype in UTR_TYPES
                else ftype
            )
            entry[bucket].append(feature)  # type: ignore[union-attr]
    return genes


def intervals(items: Iterable[Mapping[str, object]]) -> list[tuple[int, int]]:
    return sorted((int(item["start"]), int(item["end"])) for item in items)  # type: ignore[call-overload]


def merge_adjacent(parts: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(parts):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def derive_introns(exons: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    """Gaps between consecutive (sorted, non-overlapping) exons: ``(end_i + 1, start_{i+1} - 1)`` when non-empty."""

    ordered = sorted(exons)
    return [
        (a_end + 1, b_start - 1)
        for (_, a_end), (b_start, _) in itertools.pairwise(ordered)
        if b_start - a_end > 1
    ]


def derive_utrs(
    exons: Sequence[tuple[int, int]], cds_span: tuple[int, int], strand: str
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """``(five_prime, three_prime)`` UTR intervals: exon parts left/right of the CDS span (swapped on the minus strand)."""

    left: list[tuple[int, int]] = []
    right: list[tuple[int, int]] = []
    low, high = cds_span
    for start, end in sorted(exons):
        if start < low:
            left.append((start, min(end, low - 1)))
        if end > high:
            right.append((max(start, high + 1), end))
    return (left, right) if strand != "-" else (right, left)


def cds_span(entry: Mapping[str, object]) -> tuple[int, int] | None:
    """Genomic span of CDS plus start/stop codon features (GTF keeps the stop codon outside CDS)."""

    parts = intervals([*entry["cds"], *entry["start_codon"], *entry["stop_codon"]])  # type: ignore[misc]
    return (min(s for s, _ in parts), max(e for _, e in parts)) if parts else None


def transcript_strand(entry: Mapping[str, object]) -> str:
    feature = entry["feature"]
    if feature is not None:
        return str(feature["strand"])  # type: ignore[index]
    for key in ("exons", "cds"):
        for item in entry[key]:  # type: ignore[union-attr]
            return str(item["strand"])
    return "."


def transcript_seqid(entry: Mapping[str, object]) -> str | None:
    if entry["feature"] is not None:
        return str(entry["feature"]["seqid"])  # type: ignore[index]
    for key in (
        "exons",
        "cds",
        "utr",
        "five_prime_utr",
        "three_prime_utr",
        "start_codon",
        "stop_codon",
    ):
        for item in entry[key]:  # type: ignore[union-attr]
            return str(item["seqid"])
    return None


def transcript_model(entry: Mapping[str, object]) -> dict[str, object]:
    exons = intervals(entry["exons"])  # type: ignore[arg-type]
    strand = transcript_strand(entry)
    cds = intervals(entry["cds"])  # type: ignore[arg-type]
    span = cds_span(entry)
    five, three = intervals(entry["five_prime_utr"]), intervals(entry["three_prime_utr"])  # type: ignore[arg-type]
    unspecified = intervals(entry["utr"])  # type: ignore[arg-type]
    derived = False
    if not five and not three and span is not None and exons:
        if unspecified:
            five, three = [], []
        else:
            five, three = derive_utrs(exons, span, strand)
            derived = True
    feature = entry["feature"]
    start = (
        int(feature["start"])
        if feature is not None
        else (min(s for s, _ in exons) if exons else None)
    )  # type: ignore[index]
    end = (
        int(feature["end"])
        if feature is not None
        else (max(e for _, e in exons) if exons else None)
    )  # type: ignore[index]
    cds_length = sum(e - s + 1 for s, e in cds)
    return {
        "transcript_id": entry["id"],
        "type": feature["type"] if feature is not None else "transcript",
        "seqid": transcript_seqid(entry),
        "start": start,
        "end": end,
        "strand": strand,
        "name": attribute(feature, "Name", "transcript_name") if feature is not None else None,
        "exons": [list(item) for item in exons],
        "introns": [list(item) for item in derive_introns(exons)],
        "cds": [list(item) for item in cds],
        "five_prime_utr": [list(item) for item in five],
        "three_prime_utr": [list(item) for item in three],
        "utr_derived": derived,
        "exon_count": len(exons),
        "spliced_length": sum(e - s + 1 for s, e in exons),
        "cds_length": cds_length,
        "protein_length_aa": cds_length // 3 if cds_length and cds_length % 3 == 0 else None,
    }


def gene_models(features: Sequence[Mapping[str, object]], fmt: str) -> list[dict[str, object]]:
    models = []
    for gene_id, gene in build_models(features, fmt).items():
        feature = gene["feature"]
        transcripts = [transcript_model(entry) for entry in gene["transcripts"].values()]  # type: ignore[union-attr]
        starts = [t["start"] for t in transcripts if t["start"] is not None]
        ends = [t["end"] for t in transcripts if t["end"] is not None]
        models.append(
            {
                "gene_id": gene_id or None,
                "name": attribute(feature, "Name", "gene_name") if feature is not None else None,
                "type": feature["type"] if feature is not None else None,
                "seqid": feature["seqid"]
                if feature is not None
                else (transcripts[0]["seqid"] if transcripts else None),
                "start": feature["start"]
                if feature is not None
                else (min(starts) if starts else None),
                "end": feature["end"] if feature is not None else (max(ends) if ends else None),
                "strand": feature["strand"]
                if feature is not None
                else (transcripts[0]["strand"] if transcripts else "."),
                "biotype": (biotypes(feature) or [None])[0] if feature is not None else None,
                "transcripts": transcripts,
            }
        )
    return models


def public_feature(feature: Mapping[str, object]) -> dict[str, object]:
    return {key: value for key, value in feature.items()}


# --- writing ------------------------------------------------------------------------------------


def write_atomic(path: Path, text: str, *, overwrite: bool) -> None:
    if not path.parent.is_dir():
        raise ConfigurationError(
            f"Output directory does not exist: {path.parent}", code="OUTPUT_DIR_MISSING"
        )
    if path.exists() and not overwrite:
        raise ConfigurationError(
            f"Refusing to overwrite existing file: {path}",
            code="OUTPUT_EXISTS",
            hint="Pass options={'overwrite': True}.",
        )
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(text)
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _escape_gff3(value: str) -> str:
    for character, code in (
        ("%", "%25"),
        (";", "%3B"),
        ("=", "%3D"),
        ("&", "%26"),
        (",", "%2C"),
        ("\t", "%09"),
        ("\n", "%0A"),
        ("\r", "%0D"),
    ):
        value = value.replace(character, code)
    return value


def format_attributes(attributes: Mapping[str, object], fmt: str) -> str:
    if not attributes:
        return "."
    if fmt == "gff3":
        return ";".join(
            f"{_escape_gff3(str(key))}={','.join(_escape_gff3(item) for item in as_list(value))}"
            for key, value in attributes.items()
        )
    parts = []
    for key, value in attributes.items():
        for item in as_list(value):
            parts.append(f'{key} "{item.replace(chr(34), chr(92) + chr(34))}";')
    return " ".join(parts)


def feature_line(feature: Mapping[str, object], attributes: Mapping[str, object], fmt: str) -> str:
    score = "." if feature["score"] is None else format(feature["score"], "g")  # type: ignore[arg-type]
    phase = "." if feature["phase"] is None else str(feature["phase"])
    return "\t".join(
        (
            str(feature["seqid"]),
            str(feature["source"]),
            str(feature["type"]),
            str(feature["start"]),
            str(feature["end"]),
            score,
            str(feature["strand"]),
            phase,
            format_attributes(attributes, fmt),
        )
    )


def to_gtf(features: Sequence[Mapping[str, object]], source_format: str) -> str:
    """GTF lines; ``gene_id``/``transcript_id`` are resolved through the hierarchy for GFF3 input."""

    by_id = {
        identifier: f for f in features if (identifier := feature_id(f, source_format)) is not None
    }
    lines = []
    for feature in features:
        ftype = str(feature["type"]).lower()
        attributes = dict(feature["attributes"])  # type: ignore[arg-type]
        if source_format == "gff3":
            own = attributes.pop("ID", None)
            parent_ids = as_list(attributes.pop("Parent", None))
            gene_id = transcript_id = None
            if ftype in GENE_TYPES:
                gene_id = str(own) if own else None
            elif is_transcript_type(ftype):
                transcript_id = str(own) if own else None
                gene_id = parent_ids[0] if parent_ids else None
            else:
                transcript_id = parent_ids[0] if parent_ids else None
                parent = by_id.get(transcript_id or "")
                gene_id = (parents(parent, "gff3") or [None])[0] if parent is not None else None  # type: ignore[assignment]
            name = attributes.pop("Name", None)
            head: dict[str, object] = {}
            if gene_id:
                head["gene_id"] = gene_id
            if transcript_id:
                head["transcript_id"] = transcript_id
            if name is not None:
                head[
                    "gene_name"
                    if ftype in GENE_TYPES
                    else "transcript_name"
                    if is_transcript_type(ftype)
                    else "Name"
                ] = name
            attributes = {**head, **attributes}
        lines.append(feature_line(feature, attributes, "gtf"))
    return "\n".join(lines) + ("\n" if lines else "")


def to_gff3(
    features: Sequence[Mapping[str, object]], source_format: str, *, synthesize_parents: bool
) -> str:
    lines = ["##gff-version 3"]
    if source_format == "gff3":
        lines += [feature_line(f, f["attributes"], "gff3") for f in features]  # type: ignore[arg-type]
        return "\n".join(lines) + "\n"
    have_gene = {attribute(f, "gene_id") for f in features if str(f["type"]).lower() in GENE_TYPES}
    have_transcript = {
        attribute(f, "transcript_id") for f in features if is_transcript_type(str(f["type"]))
    }
    groups: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    transcript_gene: dict[str, str] = {}
    for feature in features:
        ftype = str(feature["type"]).lower()
        if ftype in _CHILD_TYPES:
            transcript_id = attribute(feature, "transcript_id")
            if transcript_id:
                groups[transcript_id].append(feature)
                transcript_gene.setdefault(
                    transcript_id, attribute(feature, "gene_id") or transcript_id
                )
    synthetic: list[dict[str, object]] = []
    if synthesize_parents:
        for transcript_id, children in groups.items():
            if transcript_id not in have_transcript:
                synthetic.append(
                    {
                        **dict(children[0]),
                        "type": "mRNA"
                        if any(str(c["type"]).lower() == "cds" for c in children)
                        else "transcript",
                        "start": min(int(c["start"]) for c in children),
                        "end": max(int(c["end"]) for c in children),  # type: ignore[call-overload]
                        "score": None,
                        "phase": None,
                        "attributes": {
                            "transcript_id": transcript_id,
                            "gene_id": transcript_gene[transcript_id],
                        },
                    }
                )
        by_gene: dict[str, list[Mapping[str, object]]] = defaultdict(list)
        for item in [*synthetic, *[f for f in features if is_transcript_type(str(f["type"]))]]:
            by_gene[attribute(item, "gene_id") or ""].append(item)
        for gene_id, members in by_gene.items():
            if gene_id and gene_id not in have_gene:
                synthetic.append(
                    {
                        **dict(members[0]),
                        "type": "gene",
                        "start": min(int(m["start"]) for m in members),
                        "end": max(int(m["end"]) for m in members),  # type: ignore[call-overload]
                        "score": None,
                        "phase": None,
                        "attributes": {"gene_id": gene_id},
                    }
                )
    counters: dict[tuple[str, str], int] = defaultdict(int)
    emitted: list[tuple[int, int, Mapping[str, object], dict[str, object]]] = []
    for order, feature in enumerate([*synthetic, *features]):
        ftype = str(feature["type"]).lower()
        attributes = dict(feature["attributes"])  # type: ignore[arg-type]
        gene_id, transcript_id = (
            attributes.pop("gene_id", None),
            attributes.pop("transcript_id", None),
        )
        name = (
            attributes.pop("gene_name", None)
            if ftype in GENE_TYPES
            else attributes.pop("transcript_name", None)
            if is_transcript_type(ftype)
            else None
        )
        head: dict[str, object] = {}
        if ftype in GENE_TYPES:
            head["ID"] = gene_id
        elif is_transcript_type(ftype):
            head["ID"] = transcript_id
            if gene_id:
                head["Parent"] = gene_id
        else:
            if transcript_id:
                counters[(str(transcript_id), ftype)] += 1
                if ftype == "exon":
                    head["ID"] = f"{transcript_id}.exon{counters[(str(transcript_id), ftype)]}"
                head["Parent"] = transcript_id
        if name is not None:
            head["Name"] = name
        emitted.append(
            (
                0 if ftype in GENE_TYPES else 1 if is_transcript_type(ftype) else 2,
                order,
                feature,
                {**head, **attributes},
            )
        )
    # parents first, then children, keeping the original order within each rank
    for _, _, feature, attributes in sorted(emitted, key=lambda item: (item[0], item[1])):
        lines.append(feature_line(feature, attributes, "gff3"))
    return "\n".join(lines) + "\n"


def bed_name(feature: Mapping[str, object], fmt: str) -> str:
    return (
        feature_id(feature, fmt)
        or attribute(feature, "Name", "gene_name", "gene_id", "transcript_id", "ID")
        or "."
    )


def bed6(features: Sequence[Mapping[str, object]], fmt: str, feature_type: str) -> str:
    wanted = feature_type.lower()
    lines = []
    for feature in features:
        if str(feature["type"]).lower() == wanted:
            score = 0 if feature["score"] is None else int(float(feature["score"]))  # type: ignore[arg-type]
            lines.append(
                "\t".join(
                    (
                        str(feature["seqid"]),
                        str(int(feature["start"]) - 1),
                        str(feature["end"]),
                        bed_name(feature, fmt),  # type: ignore[call-overload]
                        str(score),
                        str(feature["strand"] if feature["strand"] in {"+", "-"} else "."),
                    )
                )
            )
    return "\n".join(lines) + ("\n" if lines else "")


def bed12(features: Sequence[Mapping[str, object]], fmt: str) -> str:
    lines = []
    for gene in build_models(features, fmt).values():
        for entry in gene["transcripts"].values():  # type: ignore[union-attr]
            exons = intervals(entry["exons"])
            if not exons:
                continue
            start, end = exons[0][0] - 1, max(e for _, e in exons)
            span = cds_span(entry)
            thick_start, thick_end = (span[0] - 1, span[1]) if span else (end, end)
            sizes = [e - s + 1 for s, e in exons]
            offsets = [s - 1 - start for s, _ in exons]
            lines.append(
                "\t".join(
                    (
                        str(transcript_seqid(entry)),
                        str(start),
                        str(end),
                        str(entry["id"]),
                        "0",
                        transcript_strand(entry) if transcript_strand(entry) in {"+", "-"} else ".",
                        str(thick_start),
                        str(thick_end),
                        "0",
                        str(len(exons)),
                        ",".join(map(str, sizes)) + ",",
                        ",".join(map(str, offsets)) + ",",
                    )
                )
            )
    return "\n".join(lines) + ("\n" if lines else "")


# --- sequences ----------------------------------------------------------------------------------


def read_sequences(value: object) -> list[tuple[str, str]]:
    """``[(id, sequence)]`` from FASTA text/path, a plain sequence, a list of sequences/(id, seq)/mappings, or ``{id: seq}``."""

    entries: list[tuple[str, object]] = []
    if isinstance(value, Mapping):
        if "sequence" in value:
            entries = [(str(value.get("id", "seq1")), value["sequence"])]
        else:
            entries = [(str(key), item) for key, item in value.items()]
    elif isinstance(value, (str, os.PathLike)):
        text = (
            read_text(value)
            if isinstance(value, os.PathLike) or ("\n" not in value and is_file(Path(value)))
            else value
        )
        if text.lstrip().startswith(">"):
            identifier, chunks = None, []
            for line in text.splitlines():
                line = line.strip()
                if line.startswith(">"):
                    if identifier is not None:
                        entries.append((identifier, "".join(chunks)))
                    identifier, chunks = (
                        (line[1:].split(None, 1) or [f"seq{len(entries) + 1}"])[0],
                        [],
                    )
                elif line and identifier is not None:
                    chunks.append(line)
            if identifier is not None:
                entries.append((identifier, "".join(chunks)))
        else:
            entries = [("seq1", text)]
    elif isinstance(value, Sequence):
        for index, item in enumerate(value, start=1):
            if isinstance(item, Mapping):
                entries.append((str(item.get("id", f"seq{index}")), item.get("sequence", "")))
            elif isinstance(item, Sequence) and not isinstance(item, str) and len(item) == 2:
                entries.append((str(item[0]), item[1]))
            else:
                entries.append((f"seq{index}", item))
    else:
        raise InputFormatError(
            "Provide FASTA text, a sequence, or a list/mapping of sequences.",
            code="INVALID_SEQUENCE_INPUT",
        )
    result = []
    for identifier, sequence in entries:
        try:
            result.append((identifier, normalize_dna_sequence(str(sequence))))
        except ValueError as exc:
            raise InputFormatError(
                f"Sequence {identifier!r}: {exc}", code="INVALID_SEQUENCE"
            ) from exc
    if not result:
        raise InputFormatError("No sequences were found.", code="INVALID_SEQUENCE_INPUT")
    return result


def scan_orfs(
    seqid: str,
    sequence: str,
    *,
    min_length: int,
    start_codons: Sequence[str],
    stop_codons: Sequence[str],
    require_stop: bool,
    both_strands: bool,
) -> list[dict[str, object]]:
    """Open reading frames on one sequence (standard genetic code).

    In every frame an ORF starts at the first start codon that follows a stop codon (or the sequence start) and ends at
    the next in-frame stop codon, which is included.  ``length`` counts nucleotides including the stop.  Coordinates
    are 1-based closed on the forward strand; ``frame`` is +1..+3 / -1..-3.
    """

    starts, stops = set(start_codons), set(stop_codons)
    found: list[dict[str, object]] = []
    strands = [("+", sequence)] + (
        [("-", reverse_complement_dna(sequence))] if both_strands else []
    )
    total = len(sequence)
    for strand, strand_sequence in strands:
        for frame in range(3):
            begin: int | None = None
            position = frame
            last_end = frame
            while position + 3 <= total:
                codon = strand_sequence[position : position + 3]
                last_end = position + 3
                if begin is None and codon in starts:
                    begin = position
                if codon in stops and begin is not None:
                    found.append(
                        _orf_record(
                            seqid, strand, frame, strand_sequence, begin, position + 3, total, True
                        )
                    )
                    begin = None
                elif codon in stops:
                    begin = None
                position += 3
            if begin is not None and not require_stop:
                found.append(
                    _orf_record(
                        seqid, strand, frame, strand_sequence, begin, last_end, total, False
                    )
                )
    result = [orf for orf in found if int(orf["length"]) >= min_length]  # type: ignore[call-overload]
    return sorted(result, key=lambda orf: (int(orf["start"]), int(orf["end"]), str(orf["strand"])))  # type: ignore[call-overload]


def _orf_record(
    seqid: str,
    strand: str,
    frame: int,
    strand_sequence: str,
    begin: int,
    end: int,
    total: int,
    has_stop: bool,
) -> dict[str, object]:
    codons = [strand_sequence[i : i + 3] for i in range(begin, end, 3)]
    protein = "".join(CODON_TABLE.get(codon, "X") for codon in codons)
    if has_stop:
        protein = protein[:-1]
    start, stop = (begin + 1, end) if strand == "+" else (total - end + 1, total - begin)
    return {
        "seqid": seqid,
        "strand": strand,
        "frame": (frame + 1) if strand == "+" else -(frame + 1),
        "start": start,
        "end": stop,
        "length": end - begin,
        "has_stop": has_stop,
        "protein": protein,
        "protein_length": len(protein),
    }
