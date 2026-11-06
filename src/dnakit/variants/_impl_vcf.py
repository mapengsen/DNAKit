"""Native VCF / variant-record engine behind the ``variants`` and ``clinical`` candidate functions.

Pure standard library. It parses VCF text or in-memory records, classifies alleles,
normalises (trims and left-aligns) indels against an optional reference, splits and merges
records, filters, summarises and compares call sets. It never calls a variant caller and
never queries a database.

Conventions: positions are 1-based VCF POS; a call is a ``dict[str, str]`` of FORMAT values
(``GT``, ``DP``, ``AD`` ...) per sample.
"""

from __future__ import annotations

import gzip
import math
import os
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError

from .statistics import parse_vcf_records

MAX_RECORDS = 500_000
MAX_SAMPLES = 20_000
MAX_FILE_BYTES = 256_000_000

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


# ---------------------------------------------------------------------------
# errors / small helpers
# ---------------------------------------------------------------------------
def fail(message: str, code: str, hint: str, *, config: bool = False) -> Any:
    """Raise a DNAKit exception with a stable code and an actionable hint."""

    raise (ConfigurationError if config else InputFormatError)(message, code=code, hint=hint)


def option_number(
    options: Mapping[str, object],
    name: str,
    default: float | None,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
) -> float | None:
    """Read a finite numeric option (``None`` allowed when default is None)."""

    value = options.get(name, default)
    if value is None:
        return None
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
    ):
        fail(
            f"options.{name} must be a finite number.",
            "VARIANTS_BAD_OPTION",
            f"Pass a number for options.{name}.",
            config=True,
        )
    value = float(value)
    if (minimum is not None and value < minimum) or (maximum is not None and value > maximum):
        fail(
            f"options.{name}={value} is outside [{minimum}, {maximum}].",
            "VARIANTS_BAD_OPTION",
            f"Pass options.{name} within the documented range.",
            config=True,
        )
    return value


def option_int(
    options: Mapping[str, object],
    name: str,
    default: int | None,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> int | None:
    """Read an integer option."""

    value = options.get(name, default)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        if isinstance(value, float) and value.is_integer():
            value = int(value)
        else:
            fail(
                f"options.{name} must be an integer.",
                "VARIANTS_BAD_OPTION",
                f"Pass an integer for options.{name}.",
                config=True,
            )
    if (minimum is not None and value < minimum) or (maximum is not None and value > maximum):
        fail(
            f"options.{name}={value} is outside [{minimum}, {maximum}].",
            "VARIANTS_BAD_OPTION",
            f"Pass options.{name} within the documented range.",
            config=True,
        )
    return int(value)


def option_bool(options: Mapping[str, object], name: str, default: bool) -> bool:
    value = options.get(name, default)
    if not isinstance(value, bool):
        fail(
            f"options.{name} must be boolean.",
            "VARIANTS_BAD_OPTION",
            f"Pass true or false for options.{name}.",
            config=True,
        )
    return value


def option_choice(
    options: Mapping[str, object], name: str, default: str, choices: Sequence[str]
) -> str:
    value = options.get(name, default)
    if not isinstance(value, str) or value not in choices:
        fail(
            f"options.{name} must be one of {list(choices)}.",
            "VARIANTS_BAD_OPTION",
            f"Pass one of {list(choices)} for options.{name}.",
            config=True,
        )
    return value


def chrom_key(chrom: str) -> tuple[int, int, str]:
    """Natural chromosome order: 1..22, X, Y, M/MT, then other contigs alphabetically."""

    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if name.upper() in special:
        return (1, special[name.upper()], "")
    return (2, 0, name)


# ---------------------------------------------------------------------------
# data model + loading
# ---------------------------------------------------------------------------
@dataclass
class VariantSet:
    """Parsed variants: ``records`` plus aligned ``samples`` names and VCF meta ``header`` lines."""

    records: list[dict[str, Any]]
    samples: list[str] = field(default_factory=list)
    header: list[str] = field(default_factory=list)


def _read_path(path: Path) -> str:
    if not path.is_file():
        fail(
            f"file not found: {path}",
            "VARIANTS_FILE_NOT_FOUND",
            "Pass an existing VCF path, VCF text or variant records.",
        )
    if path.stat().st_size > MAX_FILE_BYTES:
        fail(
            f"file is larger than {MAX_FILE_BYTES} bytes: {path}",
            "VARIANTS_FILE_TOO_LARGE",
            "Split the file or pass a region subset.",
        )
    data = path.read_bytes()
    return _decode(data)


def _decode(data: bytes) -> str:
    if data[:2] == b"\x1f\x8b":
        try:
            data = gzip.decompress(data)
        except (OSError, EOFError) as exc:
            raise InputFormatError(
                "invalid gzip/bgzip VCF data.",
                code="VARIANTS_BAD_GZIP",
                hint="Pass an uncompressed VCF or a valid .vcf.gz.",
            ) from exc
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InputFormatError(
            "VCF data is not UTF-8 text.", code="VARIANTS_BAD_ENCODING", hint="Pass UTF-8 VCF text."
        ) from exc


def _looks_like_vcf(text: str) -> bool:
    return "\n" in text or "\t" in text or text.startswith("#")


def _call_to_str(value: object) -> str:
    if value is None:
        return "."
    if isinstance(value, (list, tuple)):
        return ",".join(_call_to_str(item) for item in value)
    return str(value)


def _alt_list(value: object) -> list[str]:
    if isinstance(value, str):
        return [item for item in value.split(",") if item and item != "."]
    if isinstance(value, Sequence):
        return [str(item) for item in value if str(item) and str(item) != "."]
    return []


def _record_from_mapping(
    item: Mapping[str, Any], index: int
) -> tuple[dict[str, Any], list[str] | None, dict[str, Any]]:
    """Convert a user record (statistics-style or VCF-like) to the internal form."""

    def pick(*names: str, default: Any = None) -> Any:
        for name in names:
            if name in item:
                return item[name]
        return default

    chrom = pick("chrom", "CHROM", "chromosome")
    pos = pick("pos", "POS", "position")
    ref = pick("ref", "REF")
    if (
        chrom in (None, "")
        or isinstance(pos, bool)
        or not isinstance(pos, int)
        or pos < 1
        or not isinstance(ref, str)
        or not ref
    ):
        fail(
            f"variant record {index} needs chrom, a positive integer pos and ref.",
            "VARIANTS_BAD_RECORD",
            "Records need chrom, pos, ref and alt (VCF-style, 1-based).",
        )
    alts = _alt_list(pick("alt", "ALT", "alts", default=[]))
    if not alts:
        fail(
            f"variant record {index} has no alternate allele.",
            "VARIANTS_BAD_RECORD",
            "Provide at least one alt allele.",
        )
    qual = pick("qual", "QUAL", "quality")
    if qual in (None, ".", ""):
        qual = None
    elif isinstance(qual, bool) or not isinstance(qual, (int, float, str)):
        fail(
            f"variant record {index} has a non-numeric qual.",
            "VARIANTS_BAD_RECORD",
            "qual must be a number or null.",
        )
    else:
        try:
            qual = float(qual)
        except ValueError:
            fail(
                f"variant record {index} has a non-numeric qual.",
                "VARIANTS_BAD_RECORD",
                "qual must be a number or null.",
            )
    flt = pick("filter", "FILTER", default=".")
    if isinstance(flt, (list, tuple)):
        flt = ";".join(str(x) for x in flt) or "."
    info = pick("info", "INFO", default={})
    if isinstance(info, str):
        parsed: dict[str, Any] = {}
        for part in info.split(";"):
            if part and part != ".":
                key, _, val = part.partition("=")
                parsed[key] = val if "=" in part else True
        info = parsed
    if not isinstance(info, Mapping):
        fail(
            f"variant record {index} has an invalid info field.",
            "VARIANTS_BAD_RECORD",
            "info must be a mapping or a VCF INFO string.",
        )
    names: list[str] | None = None
    calls: list[dict[str, str]] = []
    samples = pick("samples", "sample_genotypes")
    genotypes = pick("genotypes", "GT", "genotype")
    if isinstance(samples, Mapping):
        names = [str(key) for key in samples]
        for value in samples.values():
            calls.append(
                {str(k): _call_to_str(v) for k, v in value.items()}
                if isinstance(value, Mapping)
                else {"GT": _call_to_str(value)}
            )
    elif isinstance(samples, Sequence) and not isinstance(samples, (str, bytes)):
        for value in samples:
            calls.append(
                {str(k): _call_to_str(v) for k, v in value.items()}
                if isinstance(value, Mapping)
                else {"GT": _call_to_str(value)}
            )
    elif genotypes is not None:
        values = [genotypes] if isinstance(genotypes, (str, int)) else list(genotypes)
        calls = [{"GT": _call_to_str(v)} for v in values]
    if genotypes is not None and calls and not any("GT" in c for c in calls):
        values = [genotypes] if isinstance(genotypes, (str, int)) else list(genotypes)
        for call, gt in zip(calls, values, strict=False):
            call["GT"] = _call_to_str(gt)
    fmt: list[str] = []
    for call in calls:
        for key in call:
            if key not in fmt:
                fmt.append(key)
    if "GT" in fmt:
        fmt.remove("GT")
        fmt.insert(0, "GT")
    record = {
        "chrom": str(chrom),
        "pos": pos,
        "id": None if pick("id", "ID") in (None, ".", "") else str(pick("id", "ID")),
        "ref": ref,
        "alts": alts,
        "qual": qual,
        "filter": str(flt) if flt not in (None, "") else ".",
        "info": dict(info),
        "format": fmt,
        "calls": calls,
    }
    extras = {
        k: v
        for k, v in item.items()
        if k
        not in {
            "chrom",
            "CHROM",
            "chromosome",
            "pos",
            "POS",
            "position",
            "ref",
            "REF",
            "alt",
            "ALT",
            "alts",
            "qual",
            "QUAL",
            "quality",
            "filter",
            "FILTER",
            "info",
            "INFO",
            "id",
            "ID",
            "samples",
            "sample_genotypes",
            "genotypes",
            "GT",
            "genotype",
            "format",
            "original_alt_index",
        }
    }
    return record, names, extras


def _from_records(items: Sequence[object], metadata: Mapping[str, Any] | None = None) -> VariantSet:
    if len(items) > MAX_RECORDS:
        fail(
            f"more than {MAX_RECORDS} records.",
            "VARIANTS_TOO_MANY_RECORDS",
            "Process the callset in chunks.",
        )
    metadata = metadata or {}
    records: list[dict[str, Any]] = []
    names: list[str] | None = None
    for index, item in enumerate(items):
        if not isinstance(item, Mapping):
            fail(
                f"variant record {index} must be a mapping.",
                "VARIANTS_BAD_RECORD",
                "Pass a list of variant mappings.",
            )
        record, rec_names, extras = _record_from_mapping(item, index)
        if extras:
            # keep user annotations (gene, consequence, af, ...) in a sidecar so they survive transformations
            record["extra"] = deepcopy(extras)
        if names is None and rec_names:
            names = rec_names
        records.append(record)
    declared = metadata.get("samples", metadata.get("sample_names"))
    if (
        isinstance(declared, Sequence)
        and not isinstance(declared, (str, bytes))
        and all(isinstance(x, str) for x in declared)
    ):
        names = list(declared)
    width = max((len(r["calls"]) for r in records), default=0)
    if names is None or len(names) != width:
        if names is not None and len(names) != width:
            fail(
                "declared sample names do not match the number of genotype columns.",
                "VARIANTS_SAMPLE_MISMATCH",
                "Provide one sample name per genotype column.",
            )
        names = [f"sample_{i + 1}" for i in range(width)]
    for record in records:
        if len(record["calls"]) not in (0, width):
            fail(
                "records have different numbers of sample genotypes.",
                "VARIANTS_SAMPLE_MISMATCH",
                "Every record needs one genotype per sample.",
            )
        if not record["calls"]:
            record["calls"] = [{} for _ in range(width)]
    return VariantSet(records=records, samples=names, header=[])


def _from_vcf_text(text: str) -> VariantSet:
    parsed = parse_vcf_records(text)["result"]
    records: list[dict[str, Any]] = []
    names = list(parsed["samples"])
    if len(parsed["records"]) > MAX_RECORDS:
        fail(
            f"more than {MAX_RECORDS} records.",
            "VARIANTS_TOO_MANY_RECORDS",
            "Process the callset in chunks.",
        )
    for raw in parsed["records"]:
        calls = [
            {k: ("." if v is None else str(v)) for k, v in sample.items()}
            for sample in raw.get("samples", [])
        ]
        if len(calls) < len(names):
            calls += [{} for _ in range(len(names) - len(calls))]
        records.append(
            {
                "chrom": raw["chrom"],
                "pos": raw["pos"],
                "id": raw["id"],
                "ref": raw["ref"],
                "alts": list(raw["alt"]),
                "qual": raw["qual"],
                "filter": raw["filter"],
                "info": dict(raw["info"]),
                "format": list(raw.get("format", [])),
                "calls": calls[: len(names)] if names else [],
            }
        )
    header = [line for line in parsed["headers"] if not line.startswith("#CHROM")]
    return VariantSet(records=records, samples=names, header=header)


def load_variant_set(value: object, *, what: str = "input") -> VariantSet:
    """Load VCF text/path/bytes or variant records into a :class:`VariantSet`."""

    if isinstance(value, VariantSet):
        return value
    if isinstance(value, Path):
        return _from_vcf_text(_read_path(value))
    if isinstance(value, bytes):
        return _from_vcf_text(_decode(value))
    if isinstance(value, str):
        if _looks_like_vcf(value):
            return _from_vcf_text(value)
        if value and len(value) < 4096 and os.path.isfile(value):
            return _from_vcf_text(_read_path(Path(value)))
        fail(
            f"{what} is neither VCF text nor an existing file path.",
            "VARIANTS_BAD_INPUT",
            "Pass VCF text (with tab-separated columns), a readable VCF path, or a list of variant records.",
        )
    if isinstance(value, Mapping):
        for key in ("vcf", "text", "content"):
            if key in value and isinstance(value[key], (str, bytes, Path)):
                return load_variant_set(value[key], what=what)
        for key in ("records", "variants", "data", "items"):
            if key in value:
                inner = value[key]
                if isinstance(inner, Sequence) and not isinstance(inner, (str, bytes)):
                    return _from_records(list(inner), value)
        if any(k in value for k in ("chrom", "CHROM", "chromosome")):
            return _from_records([value], None)
        fail(
            f"{what} mapping has no records/variants/vcf key.",
            "VARIANTS_BAD_INPUT",
            "Use {'records': [...]} or {'vcf': '<text>'}.",
        )
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        items = list(value)
        if (
            items
            and all(isinstance(x, str) for x in items)
            and any(_looks_like_vcf(x) for x in items)
        ):
            return _from_vcf_text("\n".join(items))
        if all(isinstance(x, Mapping) for x in items):
            return _from_records(items, None)
        fail(
            f"{what} sequence must contain only variant mappings or VCF lines.",
            "VARIANTS_BAD_INPUT",
            "Pass a list of variant records or VCF lines.",
        )
    fail(
        f"unsupported {what} type {type(value).__name__}.",
        "VARIANTS_BAD_INPUT",
        "Pass VCF text, a VCF path, or a list of variant records.",
    )
    raise AssertionError("unreachable")


# ---------------------------------------------------------------------------
# serialisation
# ---------------------------------------------------------------------------
def _info_text(info: Mapping[str, Any]) -> str:
    parts = []
    for key, value in info.items():
        if value is True:
            parts.append(str(key))
        elif value is False or value is None:
            continue
        elif isinstance(value, (list, tuple)):
            parts.append(f"{key}={','.join(str(x) for x in value)}")
        else:
            parts.append(f"{key}={value}")
    return ";".join(parts) or "."


def _qual_text(qual: float | None) -> str:
    if qual is None:
        return "."
    return str(int(qual)) if float(qual).is_integer() else repr(float(qual))


def to_vcf_text(vs: VariantSet) -> str:
    """Serialise a variant set as VCF text (meta lines are kept, ``fileformat`` is ensured)."""

    header = [line for line in vs.header if not line.startswith("#CHROM")]
    if not any(line.startswith("##fileformat") for line in header):
        header.insert(0, "##fileformat=VCFv4.2")
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if vs.samples:
        columns += ["FORMAT", *vs.samples]
    lines = [*header, "\t".join(columns)]
    for rec in vs.records:
        row = [
            rec["chrom"],
            str(rec["pos"]),
            rec["id"] or ".",
            rec["ref"],
            ",".join(rec["alts"]) or ".",
            _qual_text(rec["qual"]),
            rec["filter"] or ".",
            _info_text(rec["info"]),
        ]
        if vs.samples:
            fmt = rec["format"] or ["GT"]
            row.append(":".join(fmt))
            for call in rec["calls"]:
                values = [call.get(key, ".") for key in fmt]
                while len(values) > 1 and values[-1] == ".":
                    values.pop()
                row.append(":".join(values))
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def record_view(rec: Mapping[str, Any], samples: Sequence[str]) -> dict[str, Any]:
    """JSON-friendly view of an internal record."""

    view: dict[str, Any] = {
        "chrom": rec["chrom"],
        "pos": rec["pos"],
        "id": rec["id"],
        "ref": rec["ref"],
        "alt": list(rec["alts"]),
        "qual": rec["qual"],
        "filter": rec["filter"],
        "info": dict(rec["info"]),
    }
    if samples:
        view["genotypes"] = {
            name: dict(call) for name, call in zip(samples, rec["calls"], strict=False)
        }
    if rec.get("extra"):
        view.update({k: v for k, v in rec["extra"].items() if k not in view})
    if rec.get("hgvs_g"):
        view["hgvs_g"] = list(rec["hgvs_g"])
    return view


def write_output(options: Mapping[str, object], text: str) -> str | None:
    """Write ``text`` to ``options['output']`` when given (explicit path only). Returns the path."""

    target = options.get("output")
    if target is None:
        return None
    if not isinstance(target, (str, Path)) or not str(target):
        fail(
            "options.output must be a path.",
            "VARIANTS_BAD_OPTION",
            "Pass a file path string.",
            config=True,
        )
    path = Path(target)
    overwrite = option_bool(options, "overwrite", False)
    if path.exists() and not overwrite:
        fail(
            f"output exists: {path}",
            "VARIANTS_OUTPUT_EXISTS",
            "Pass options.overwrite=true or choose a new path.",
            config=True,
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return str(path)


# ---------------------------------------------------------------------------
# genotypes and per-call numbers
# ---------------------------------------------------------------------------
_GT_RE = re.compile(r"^(?:\.|\d+)(?:[/|](?:\.|\d+))*$")


def parse_gt(text: object) -> tuple[list[int | None], bool]:
    """Parse a GT string into (allele indices with None for missing, phased flag)."""

    if text is None:
        return [None], False
    value = str(text).strip()
    if value in {"", "."}:
        return [None], False
    if not _GT_RE.match(value):
        fail(
            f"invalid genotype {value!r}.",
            "VARIANTS_BAD_GENOTYPE",
            "Genotypes look like 0/1, 1|0, ./. or 1.",
        )
    alleles = [None if p == "." else int(p) for p in re.split(r"[/|]", value)]
    return alleles, "|" in value


def call_gt(call: Mapping[str, str]) -> tuple[list[int | None], bool]:
    return parse_gt(call.get("GT"))


def is_called(alleles: Sequence[int | None]) -> bool:
    return bool(alleles) and all(a is not None for a in alleles)


def number_or_none(value: object) -> float | None:
    if value in (None, "", "."):
        return None
    try:
        number = float(str(value))
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def call_ad(call: Mapping[str, str]) -> list[int] | None:
    """Allelic depths of a call (``AD`` = ref,alt1,...) or None."""

    raw = call.get("AD")
    if raw in (None, "", "."):
        return None
    try:
        return [int(float(part)) for part in str(raw).split(",")]
    except ValueError:
        return None


def call_depth(call: Mapping[str, str]) -> float | None:
    """Read depth of a call: ``DP`` when present, otherwise the sum of ``AD``."""

    dp = number_or_none(call.get("DP"))
    if dp is not None:
        return dp
    ad = call_ad(call)
    return float(sum(ad)) if ad else None


def site_depth(rec: Mapping[str, Any]) -> float | None:
    """Site depth: INFO/DP when present, else the sum of per-sample depths."""

    info_dp = number_or_none(rec["info"].get("DP"))
    if info_dp is not None:
        return info_dp
    values = [d for d in (call_depth(c) for c in rec["calls"]) if d is not None]
    return float(sum(values)) if values else None


def allele_counts(rec: Mapping[str, Any]) -> tuple[list[int], int]:
    """(AC per allele index 0..n_alt, AN) from the called alleles of all samples."""

    counts = [0] * (len(rec["alts"]) + 1)
    total = 0
    for call in rec["calls"]:
        alleles, _ = call_gt(call)
        for allele in alleles:
            if allele is None:
                continue
            if allele >= len(counts):
                fail(
                    f"genotype allele {allele} exceeds ALT count at {rec['chrom']}:{rec['pos']}.",
                    "VARIANTS_BAD_GENOTYPE",
                    "Genotype indices must not exceed the number of ALT alleles.",
                )
            counts[allele] += 1
            total += 1
    return counts, total


def diploid_class(alleles: Sequence[int | None]) -> str:
    """Return 'missing', 'hom_ref', 'het', 'hom_alt' or 'other' (non-diploid) for a call."""

    if len(alleles) != 2 or not is_called(alleles):
        return "missing" if not is_called(alleles) else "other"
    a, b = alleles
    if a != b:
        return "het"
    return "hom_ref" if a == 0 else "hom_alt"


def genotype_counts(rec: Mapping[str, Any]) -> dict[str, int]:
    """Counts of diploid genotype classes at a site."""

    counts = Counter({"hom_ref": 0, "het": 0, "hom_alt": 0, "missing": 0, "other": 0})
    for call in rec["calls"]:
        counts[diploid_class(call_gt(call)[0])] += 1
    return dict(counts)


# ---------------------------------------------------------------------------
# classification
# ---------------------------------------------------------------------------
def is_symbolic(alt: str) -> bool:
    return alt.startswith("<") or "[" in alt or "]" in alt or alt in {"*", "."} or alt.endswith(">")


def trim_alleles(ref: str, alt: str) -> tuple[int, str, str]:
    """Trim common suffix then prefix (may leave empty strings); returns (offset, ref, alt)."""

    offset = 0
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        offset += 1
    return offset, ref, alt


def classify_alleles(ref: str, alt: str, *, sv_min_length: int = 50) -> dict[str, Any]:
    """Classify one REF/ALT pair by allele lengths (after trimming shared flanks).

    ``type`` is one of SNV, MNV, indel, SV, ref, other; ``subtype`` is insertion/deletion/complex,
    the symbolic ALT, or ``spanning_deletion``.
    """

    ref, alt = ref.upper(), alt.upper()
    if alt == "*":
        return {"type": "other", "subtype": "spanning_deletion", "length_change": 0}
    if alt.startswith("<") or "[" in alt or "]" in alt:
        symbol = alt.strip("<>").upper() if alt.startswith("<") else "BND"
        return {"type": "SV", "subtype": symbol, "length_change": None}
    if alt == ref:
        return {"type": "ref", "subtype": None, "length_change": 0}
    _, tref, talt = trim_alleles(ref, alt)
    change = len(alt) - len(ref)
    if abs(change) >= sv_min_length:
        return {
            "type": "SV",
            "subtype": "insertion" if change > 0 else "deletion",
            "length_change": change,
        }
    if not tref and talt:
        return {"type": "indel", "subtype": "insertion", "length_change": change}
    if tref and not talt:
        return {"type": "indel", "subtype": "deletion", "length_change": change}
    if len(tref) == len(talt):
        return {"type": "SNV" if len(tref) == 1 else "MNV", "subtype": None, "length_change": 0}
    return {"type": "indel", "subtype": "complex", "length_change": change}


def is_transition(ref: str, alt: str) -> bool | None:
    """True/False for a transition/transversion SNV, None when the pair is not a plain SNV."""

    _, tref, talt = trim_alleles(ref.upper(), alt.upper())
    if len(tref) == 1 and len(talt) == 1 and tref in "ACGT" and talt in "ACGT":
        return (tref, talt) in TRANSITIONS
    return None


# ---------------------------------------------------------------------------
# reference sequence
# ---------------------------------------------------------------------------
class Reference:
    """In-memory reference: contig -> sequence, or one plain sequence applied to every contig."""

    def __init__(self, seqs: dict[str, str], start: int = 1) -> None:
        self.seqs = seqs
        self.start = start

    def fetch(self, chrom: str, pos: int, length: int = 1) -> str | None:
        """Bases at 1-based [pos, pos+length-1] or None when outside the sequence."""

        seq = self.seqs.get(chrom)
        offset = 1
        if seq is None:
            alt_name = chrom[3:] if chrom.lower().startswith("chr") else f"chr{chrom}"
            seq = self.seqs.get(alt_name)
        if seq is None and "*" in self.seqs:
            seq, offset = self.seqs["*"], self.start
        if seq is None:
            return None
        index = pos - offset
        if index < 0 or length < 0 or index + length > len(seq):
            return None
        return seq[index : index + length]


def _parse_fasta(text: str) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    current = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0] if line[1:].strip() else ""
            if not current:
                fail(
                    "FASTA record without a name.",
                    "VARIANTS_BAD_REFERENCE",
                    "Give every FASTA record a name.",
                )
            seqs.setdefault(current, [])
        elif current is None:
            fail(
                "FASTA text must start with '>'.",
                "VARIANTS_BAD_REFERENCE",
                "Pass a FASTA string, a path, a sequence or a contig mapping.",
            )
        else:
            seqs[current].append(line.upper())
    return {name: "".join(parts) for name, parts in seqs.items()}


def load_reference(value: object, start: int = 1) -> Reference | None:
    """Build a Reference from a plain sequence, FASTA text/path or contig mapping (None stays None)."""

    if value is None:
        return None
    if isinstance(value, Mapping):
        seqs = {}
        for name, seq in value.items():
            if not isinstance(seq, str):
                fail(
                    "reference mapping values must be sequences.",
                    "VARIANTS_BAD_REFERENCE",
                    "Use {'chr1': 'ACGT...'}.",
                )
            seqs[str(name)] = seq.upper()
        return Reference(seqs)
    if isinstance(value, Path) or (
        isinstance(value, str)
        and "\n" not in value
        and not value.startswith(">")
        and not re.fullmatch(r"[A-Za-z*-]+", value)
        and os.path.isfile(value)
    ):
        return Reference(_parse_fasta(_read_path(Path(value))))
    if isinstance(value, str):
        if value.lstrip().startswith(">"):
            return Reference(_parse_fasta(value))
        if re.fullmatch(r"[A-Za-z]+", value):
            return Reference({"*": value.upper()}, start=start)
    fail(
        "unsupported reference.",
        "VARIANTS_BAD_REFERENCE",
        "Pass options.reference as a plain sequence, FASTA text, a FASTA path or {contig: sequence}.",
    )
    raise AssertionError("unreachable")


def reference_from_options(options: Mapping[str, object]) -> Reference | None:
    start = option_int(options, "reference_start", 1, minimum=1) or 1
    return load_reference(options.get("reference"), start)


# ---------------------------------------------------------------------------
# normalisation (trim + left-align, Tan et al. 2015)
# ---------------------------------------------------------------------------
def normalize_alleles(
    chrom: str, pos: int, alleles: Sequence[str], ref: Reference | None
) -> tuple[int, list[str], bool]:
    """Trim and (when ``ref`` is given) left-align a set of alleles (index 0 = REF).

    Algorithm of Tan, Abecasis & Kang (2015): repeatedly drop the shared last base; whenever an allele
    becomes empty, prepend the preceding reference base to all alleles; finally drop shared leading
    bases while every allele keeps at least one base. Without a reference only trimming that keeps an
    anchor base is done. Returns (pos, alleles, left_aligned).
    """

    alleles = [a.upper() for a in alleles]
    if ref is None:
        while all(len(a) >= 2 for a in alleles) and len({a[-1] for a in alleles}) == 1:
            alleles = [a[:-1] for a in alleles]
    else:
        padded_right = False
        guard = 0
        while True:
            guard += 1
            if guard > 1_000_000:
                break
            changed = False
            if not padded_right and all(alleles) and len({a[-1] for a in alleles}) == 1:
                alleles = [a[:-1] for a in alleles]
                changed = True
            if any(not a for a in alleles):
                left = ref.fetch(chrom, pos - 1, 1) if pos > 1 else None
                if left is not None:
                    alleles = [left + a for a in alleles]
                    pos -= 1
                    changed = True
                else:
                    right = ref.fetch(chrom, pos + len(alleles[0]), 1)
                    if right is None:
                        break
                    alleles = [a + right for a in alleles]
                    padded_right = True
                    changed = True
            if not changed:
                break
    while all(len(a) >= 2 for a in alleles) and len({a[0] for a in alleles}) == 1:
        alleles = [a[1:] for a in alleles]
        pos += 1
    return pos, alleles, ref is not None


def normalize_set(
    vs: VariantSet,
    ref: Reference | None,
    *,
    check_ref: str = "warn",
    split: bool = False,
    sort: bool = True,
) -> tuple[VariantSet, dict[str, Any]]:
    """Normalise every record of ``vs``; returns the new set and a report."""

    warnings: list[str] = []
    changed = 0
    skipped = 0
    out: list[dict[str, Any]] = []
    changes: list[dict[str, Any]] = []
    for rec in vs.records:
        new = deepcopy(rec)
        alleles = [rec["ref"], *rec["alts"]]
        if any(is_symbolic(a) for a in alleles[1:]):
            out.append(new)
            skipped += 1
            continue
        if ref is not None:
            actual = ref.fetch(rec["chrom"], rec["pos"], len(rec["ref"]))
            if actual is None or actual.upper() != rec["ref"].upper():
                message = (
                    f"{rec['chrom']}:{rec['pos']} REF {rec['ref']} is not found in the reference"
                    f" (reference has {actual!r}); record left unchanged."
                )
                if check_ref == "error":
                    fail(
                        message,
                        "VARIANTS_REF_MISMATCH",
                        "Fix the REF allele or pass the matching reference.",
                    )
                warnings.append(message)
                out.append(new)
                skipped += 1
                continue
        pos, norm, _ = normalize_alleles(rec["chrom"], rec["pos"], alleles, ref)
        if (pos, norm) != (rec["pos"], [a.upper() for a in alleles]):
            changed += 1
            if len(changes) < 100:
                changes.append(
                    {
                        "chrom": rec["chrom"],
                        "from": f"{rec['pos']} {rec['ref']}>{','.join(rec['alts'])}",
                        "to": f"{pos} {norm[0]}>{','.join(norm[1:])}",
                    }
                )
        new["pos"], new["ref"], new["alts"] = pos, norm[0], norm[1:]
        out.append(new)
    result = VariantSet(records=out, samples=list(vs.samples), header=list(vs.header))
    if split:
        result = split_multiallelic_set(result)
    if sort:
        sort_records(result.records)
    return result, {
        "changed_records": changed,
        "unchanged_or_skipped": len(vs.records) - changed,
        "skipped_records": skipped,
        "left_aligned": ref is not None,
        "changes": changes,
        "warnings": warnings,
    }


def sort_records(records: list[dict[str, Any]]) -> None:
    """Sort in place by natural chromosome order, POS, REF, ALT list."""

    records.sort(key=lambda r: (chrom_key(r["chrom"]), r["pos"], r["ref"], tuple(r["alts"])))


# ---------------------------------------------------------------------------
# HGVS (genomic, g.) strings for simple variants
# ---------------------------------------------------------------------------
def build_hgvs_genomic(
    accession: str,
    pos: int,
    ref: str,
    alt: str,
    reference: Reference | None = None,
    chrom: str | None = None,
) -> str:
    """Build an HGVS ``g.`` description for one VCF-style SNV, MNV, insertion, deletion or delins.

    With a reference the indel is first 3'-shifted (HGVS rule) and an insertion that repeats the preceding
    sequence is written as ``dup``. Without a reference the description is for the allele as given (not
    3'-normalised) and duplications are written as insertions.
    """

    ref, alt = ref.upper(), alt.upper()
    if is_symbolic(alt):
        fail(
            f"cannot build HGVS for symbolic allele {alt}.",
            "VARIANTS_HGVS_UNSUPPORTED",
            "Only sequence alleles are supported.",
        )
    offset, tref, talt = trim_alleles(ref, alt)
    start = pos + offset
    prefix = f"{accession}:g."
    if len(tref) == 1 and len(talt) == 1:
        return f"{prefix}{start}{tref}>{talt}"
    if tref and talt:
        end = start + len(tref) - 1
        span = f"{start}" if start == end else f"{start}_{end}"
        return f"{prefix}{span}delins{talt}"
    if tref and not talt:  # deletion
        seq = tref
        if reference is not None and chrom is not None:
            while True:
                nxt = reference.fetch(chrom, start + len(seq), 1)
                if nxt is None or nxt.upper() != seq[0]:
                    break
                seq, start = seq[1:] + seq[0], start + 1
        end = start + len(seq) - 1
        return f"{prefix}{start}del" if start == end else f"{prefix}{start}_{end}del"
    if not tref and not talt:
        fail("alleles are identical.", "VARIANTS_HGVS_UNSUPPORTED", "REF and ALT must differ.")
    seq = talt  # insertion between start-1 and start
    after = start - 1
    if reference is not None and chrom is not None:
        while True:
            nxt = reference.fetch(chrom, after + 1, 1)
            if nxt is None or nxt.upper() != seq[0]:
                break
            seq, after = seq[1:] + seq[0], after + 1
        before = (
            reference.fetch(chrom, after - len(seq) + 1, len(seq))
            if after - len(seq) + 1 >= 1
            else None
        )
        if before is not None and before.upper() == seq:
            first = after - len(seq) + 1
            return f"{prefix}{first}dup" if first == after else f"{prefix}{first}_{after}dup"
    return f"{prefix}{after}_{after + 1}ins{seq}"


# ---------------------------------------------------------------------------
# split multi-allelic
# ---------------------------------------------------------------------------
def _remap_gt(text: str, selected: int) -> str:
    alleles, phased = parse_gt(text)
    if not any(a is not None for a in alleles) and text in {".", "./.", ".|."}:
        return text
    sep = "|" if phased else "/"
    return sep.join("." if a is None else ("1" if a == selected else "0") for a in alleles)


def split_multiallelic_set(vs: VariantSet) -> VariantSet:
    """One record per ALT allele. GT is remapped (selected ALT -> 1, every other allele -> 0), ``AD``
    (Number=R) and per-allele INFO ``AC``/``AF`` (Number=A) are subset; ``PL`` (Number=G) is dropped."""

    out: list[dict[str, Any]] = []
    for rec in vs.records:
        if len(rec["alts"]) <= 1:
            out.append(deepcopy(rec))
            continue
        for index, alt in enumerate(rec["alts"], 1):
            new = deepcopy(rec)
            new["alts"] = [alt]
            new["calls"] = []
            for call in rec["calls"]:
                fresh = dict(call)
                if "GT" in fresh:
                    fresh["GT"] = _remap_gt(fresh["GT"], index)
                ad = call_ad(call)
                if ad is not None and len(ad) == len(rec["alts"]) + 1:
                    fresh["AD"] = f"{ad[0]},{ad[index]}"
                fresh.pop("PL", None)
                new["calls"].append(fresh)
            new["format"] = [k for k in rec["format"] if k != "PL"]
            for key in ("AC", "AF", "MLEAC", "MLEAF"):
                raw = rec["info"].get(key)
                if isinstance(raw, str) and len(raw.split(",")) == len(rec["alts"]):
                    new["info"][key] = raw.split(",")[index - 1]
            new["extra"] = deepcopy(rec.get("extra", {}))
            new["split_from_alt_index"] = index
            out.append(new)
    return VariantSet(records=out, samples=list(vs.samples), header=list(vs.header))


# ---------------------------------------------------------------------------
# merge
# ---------------------------------------------------------------------------
def load_many(value: object) -> list[VariantSet]:
    """Load several call sets from a sequence of VCF texts/paths/mappings (or ``{'inputs': [...]}``)."""

    if isinstance(value, Mapping) and "inputs" in value:
        value = value["inputs"]
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        fail(
            "merge input must be a list of VCFs / call sets.",
            "VARIANTS_BAD_INPUT",
            "Pass a list with two or more VCF texts, paths or {'records': [...]} mappings.",
        )
    items = list(value)
    if items and all(
        isinstance(x, Mapping) and any(k in x for k in ("chrom", "CHROM")) for x in items
    ):
        fail(
            "got a single record list, but merging needs two or more call sets.",
            "VARIANTS_BAD_INPUT",
            "Wrap each call set: [{'records': [...]}, {'records': [...]}].",
        )
    if len(items) < 2:
        fail(
            "merging needs at least two call sets.", "VARIANTS_BAD_INPUT", "Pass two or more VCFs."
        )
    return [load_variant_set(item, what=f"call set {i}") for i, item in enumerate(items)]


def merge_sets(sets: Sequence[VariantSet], mode: str = "samples") -> VariantSet:
    """Merge call sets.

    ``mode='samples'`` (like ``bcftools merge -m none``): records with identical CHROM/POS/REF/ALT-list are
    combined, adding the samples of every file (absent samples get ``./.``); sample names must be unique.
    QUAL = max of the contributing records, FILTER = PASS only if every contributing record is PASS
    (else the union of failing labels), INFO = first record carrying the key, ID = first non-missing.
    ``mode='records'`` (like ``bcftools concat``): all files must have the same samples; identical
    sites are kept once (first occurrence).
    """

    header: list[str] = []
    for vs in sets:
        for line in vs.header:
            if line not in header and not line.startswith("#CHROM"):
                header.append(line)
    if mode == "records":
        names = list(sets[0].samples)
        if any(vs.samples != names for vs in sets):
            fail(
                "records mode requires identical sample lists in all inputs.",
                "VARIANTS_SAMPLE_MISMATCH",
                "Use mode='samples' to merge different samples.",
            )
        seen: set[tuple[Any, ...]] = set()
        records = []
        for vs in sets:
            for rec in vs.records:
                key = (rec["chrom"], rec["pos"], rec["ref"], tuple(rec["alts"]))
                if key not in seen:
                    seen.add(key)
                    records.append(deepcopy(rec))
        sort_records(records)
        return VariantSet(records, names, header)
    names = []
    for vs in sets:
        for name in vs.samples:
            if name in names:
                fail(
                    f"duplicate sample name {name!r} across inputs.",
                    "VARIANTS_DUPLICATE_SAMPLE",
                    "Rename samples so that every sample name is unique before merging.",
                )
            names.append(name)
    offsets = []
    cursor = 0
    for vs in sets:
        offsets.append(cursor)
        cursor += len(vs.samples)
    merged: dict[tuple[Any, ...], dict[str, Any]] = {}
    for set_index, vs in enumerate(sets):
        for rec in vs.records:
            key = (rec["chrom"], rec["pos"], rec["ref"], tuple(rec["alts"]))
            slot = merged.get(key)
            if slot is None:
                slot = merged[key] = {
                    "rec": deepcopy(rec),
                    "calls": [None] * len(names),
                    "filters": [],
                    "quals": [],
                }
                slot["rec"]["calls"] = []
            slot["filters"].append(rec["filter"])
            if rec["qual"] is not None:
                slot["quals"].append(rec["qual"])
            if slot["rec"]["id"] is None and rec["id"]:
                slot["rec"]["id"] = rec["id"]
            for k, v in rec["info"].items():
                slot["rec"]["info"].setdefault(k, v)
            for fk in rec["format"]:
                if fk not in slot["rec"]["format"]:
                    slot["rec"]["format"].append(fk)
            for j, call in enumerate(rec["calls"]):
                slot["calls"][offsets[set_index] + j] = dict(call)
    records = []
    for slot in merged.values():
        rec = slot["rec"]
        fmt = list(rec["format"]) or (["GT"] if names else [])
        if "GT" in fmt:
            fmt.remove("GT")
            fmt.insert(0, "GT")
        rec["format"] = fmt
        rec["calls"] = [call if call is not None else {"GT": "./."} for call in slot["calls"]]
        rec["qual"] = max(slot["quals"]) if slot["quals"] else None
        failing = []
        for flt in slot["filters"]:
            for label in flt.split(";"):
                if label not in {"PASS", "."} and label not in failing:
                    failing.append(label)
        if failing:
            rec["filter"] = ";".join(failing)
        else:
            rec["filter"] = "PASS" if any(f == "PASS" for f in slot["filters"]) else "."
        records.append(rec)
    sort_records(records)
    return VariantSet(records, names, header)


# ---------------------------------------------------------------------------
# HWE exact test (Wigginton, Cutler & Abecasis 2005)
# ---------------------------------------------------------------------------
def hwe_exact_pvalue(n_aa: int, n_ab: int, n_bb: int, *, mid_p: bool = False) -> float:
    """Two-sided exact Hardy-Weinberg p-value: the sum of the probabilities of all heterozygote
    counts (given the allele counts) that are not more likely than the observed one."""

    if min(n_aa, n_ab, n_bb) < 0:
        fail(
            "genotype counts must be non-negative.",
            "VARIANTS_BAD_COUNTS",
            "Pass non-negative integer counts.",
        )
    n = n_aa + n_ab + n_bb
    if n == 0:
        return 1.0
    hom_rare, _hom_common = sorted((n_aa, n_bb))
    rare = 2 * hom_rare + n_ab
    probs = [0.0] * (rare + 1)
    mid = rare * (2 * n - rare) // (2 * n)
    if (mid % 2) != (rare % 2):
        mid += 1
    probs[mid] = 1.0
    total = 1.0
    hets, homr = mid, (rare - mid) // 2
    homc = n - hets - homr
    while hets >= 2:
        probs[hets - 2] = probs[hets] * hets * (hets - 1.0) / (4.0 * (homr + 1.0) * (homc + 1.0))
        total += probs[hets - 2]
        hets, homr, homc = hets - 2, homr + 1, homc + 1
    hets, homr = mid, (rare - mid) // 2
    homc = n - hets - homr
    while hets <= rare - 2:
        probs[hets + 2] = probs[hets] * 4.0 * homr * homc / ((hets + 2.0) * (hets + 1.0))
        total += probs[hets + 2]
        hets, homr, homc = hets + 2, homr - 1, homc - 1
    probs = [p / total for p in probs]
    observed = probs[n_ab]
    p_value = sum(p for p in probs if p <= observed * (1 + 1e-12))
    if mid_p:
        p_value -= 0.5 * observed
    return min(1.0, max(0.0, p_value))


def site_hwe(rec: Mapping[str, Any]) -> float | None:
    """HWE exact p-value of a biallelic diploid site (None for multi-allelic or no calls)."""

    if len(rec["alts"]) != 1:
        return None
    counts = genotype_counts(rec)
    if counts["other"] or counts["hom_ref"] + counts["het"] + counts["hom_alt"] == 0:
        return None
    return hwe_exact_pvalue(counts["hom_ref"], counts["het"], counts["hom_alt"])


# ---------------------------------------------------------------------------
# filtering
# ---------------------------------------------------------------------------
def _region_list(value: object) -> list[tuple[str, int, int]]:
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, Sequence):
        fail(
            "options.regions must be a list of 'chr:start-end' strings.",
            "VARIANTS_BAD_OPTION",
            "Example: ['chr1:1000-2000'].",
            config=True,
        )
    out = []
    for item in value:
        match = re.fullmatch(r"([^:]+)(?::(\d+)(?:-(\d+))?)?", str(item).replace(",", ""))
        if not match:
            fail(
                f"invalid region {item!r}.",
                "VARIANTS_BAD_OPTION",
                "Use chrom, chrom:pos or chrom:start-end.",
                config=True,
            )
        start = int(match.group(2)) if match.group(2) else 1
        end = int(match.group(3)) if match.group(3) else (start if match.group(2) else 10**12)
        if end < start:
            fail(
                f"invalid region {item!r}: end < start.",
                "VARIANTS_BAD_OPTION",
                "Use start <= end.",
                config=True,
            )
        out.append((match.group(1), start, end))
    return out


def alt_frequency(rec: Mapping[str, Any]) -> float | None:
    """Total ALT allele frequency from called genotypes, else INFO/AF (sum of listed values)."""

    if rec["calls"]:
        counts, total = allele_counts(rec)
        if total:
            return 1.0 - counts[0] / total
    raw = rec["info"].get("AF")
    if isinstance(raw, str):
        try:
            return sum(float(x) for x in raw.split(","))
        except ValueError:
            return None
    if isinstance(raw, (int, float)) and not isinstance(raw, bool):
        return float(raw)
    return None


def filter_set(vs: VariantSet, options: Mapping[str, object]) -> tuple[VariantSet, dict[str, Any]]:
    """Filter records. See :func:`filter_variants` docstrings for the option list."""

    pass_only = option_bool(options, "pass_only", False)
    min_qual = option_number(options, "min_qual", None)
    min_depth = option_number(options, "min_depth", None, minimum=0)
    max_depth = option_number(options, "max_depth", None, minimum=0)
    min_gq = option_number(options, "min_genotype_quality", None, minimum=0)
    min_sdp = option_number(options, "min_sample_depth", None, minimum=0)
    max_missing = option_number(options, "max_missing", None, minimum=0, maximum=1)
    min_af = option_number(options, "min_af", None, minimum=0, maximum=1)
    max_af = option_number(options, "max_af", None, minimum=0, maximum=1)
    min_hwe = option_number(options, "min_hwe_pvalue", None, minimum=0, maximum=1)
    mode = option_choice(options, "mode", "remove", ("remove", "flag"))
    types = options.get("variant_types")
    if types is not None:
        if isinstance(types, str):
            types = [types]
        allowed = {"SNV", "MNV", "indel", "SV"}
        if not isinstance(types, Sequence) or any(t not in allowed for t in types):
            fail(
                f"options.variant_types must be a list from {sorted(allowed)}.",
                "VARIANTS_BAD_OPTION",
                "Use SNV, MNV, indel or SV.",
                config=True,
            )
    chroms = options.get("chroms")
    if chroms is not None:
        chroms = {str(c) for c in ([chroms] if isinstance(chroms, str) else chroms)}
    regions = _region_list(options["regions"]) if options.get("regions") is not None else None

    out: list[dict[str, Any]] = []
    reasons: Counter[str] = Counter()
    masked = 0
    for rec in vs.records:
        new = deepcopy(rec)
        if min_gq is not None or min_sdp is not None:
            for call in new["calls"]:
                gq = number_or_none(call.get("GQ"))
                dp = call_depth(call)
                if (min_gq is not None and (gq is None or gq < min_gq)) or (
                    min_sdp is not None and (dp is None or dp < min_sdp)
                ):
                    if call.get("GT") not in (None, ".", "./.", ".|."):
                        masked += 1
                    call["GT"] = "./."
        failed: list[str] = []
        if pass_only and rec["filter"] not in {"PASS", "."}:
            failed.append("not_pass")
        if min_qual is not None and (rec["qual"] is None or rec["qual"] < min_qual):
            failed.append("low_qual")
        depth = site_depth(new)
        if min_depth is not None and (depth is None or depth < min_depth):
            failed.append("low_depth")
        if max_depth is not None and depth is not None and depth > max_depth:
            failed.append("high_depth")
        if max_missing is not None and new["calls"]:
            missing = sum(1 for c in new["calls"] if not is_called(call_gt(c)[0])) / len(
                new["calls"]
            )
            if missing > max_missing:
                failed.append("high_missing")
        if min_af is not None or max_af is not None:
            freq = alt_frequency(new)
            if (
                freq is None
                or (min_af is not None and freq < min_af)
                or (max_af is not None and freq > max_af)
            ):
                failed.append("af_out_of_range")
        if min_hwe is not None:
            p_value = site_hwe(new)
            if p_value is not None and p_value < min_hwe:
                failed.append("hwe")
        if types is not None:
            kinds = {classify_alleles(rec["ref"], alt)["type"] for alt in rec["alts"]}
            if not kinds & set(types):
                failed.append("variant_type")
        if chroms is not None and rec["chrom"] not in chroms:
            failed.append("chrom")
        if regions is not None and not any(
            rec["chrom"] == c and s <= rec["pos"] <= e for c, s, e in regions
        ):
            failed.append("region")
        for reason in failed:
            reasons[reason] += 1
        if failed and mode == "remove":
            continue
        if failed:
            existing = [] if new["filter"] in {"PASS", "."} else new["filter"].split(";")
            new["filter"] = ";".join(existing + [r for r in failed if r not in existing])
        elif mode == "flag" and new["filter"] == ".":
            new["filter"] = "PASS"
        out.append(new)
    result = VariantSet(out, list(vs.samples), list(vs.header))
    return result, {
        "input_records": len(vs.records),
        "output_records": len(out),
        "failed_records": sum(1 for r in out if r["filter"] not in {"PASS", "."})
        if mode == "flag"
        else len(vs.records) - len(out),
        "reason_counts": dict(sorted(reasons.items())),
        "genotypes_masked": masked,
        "mode": mode,
    }


# ---------------------------------------------------------------------------
# summaries and QC
# ---------------------------------------------------------------------------
def sample_qc(vs: VariantSet) -> list[dict[str, Any]]:
    """Per-sample QC metrics from the diploid calls of all records."""

    rows = []
    for index, name in enumerate(vs.samples):
        counts: Counter[str] = Counter()
        ti = tv = snv = indel = 0
        depths: list[float] = []
        gqs: list[float] = []
        for rec in vs.records:
            call = rec["calls"][index]
            alleles, _ = call_gt(call)
            counts[diploid_class(alleles)] += 1
            dp = call_depth(call)
            if dp is not None:
                depths.append(dp)
            gq = number_or_none(call.get("GQ"))
            if gq is not None:
                gqs.append(gq)
            if is_called(alleles):
                for allele in {a for a in alleles if a}:
                    if allele > len(rec["alts"]):
                        continue
                    alt = rec["alts"][allele - 1]
                    kind = classify_alleles(rec["ref"], alt)["type"]
                    if kind == "SNV":
                        snv += 1
                        trans = is_transition(rec["ref"], alt)
                        if trans is True:
                            ti += 1
                        elif trans is False:
                            tv += 1
                    elif kind == "indel":
                        indel += 1
        total = len(vs.records)
        called = total - counts["missing"]
        rows.append(
            {
                "sample": name,
                "n_sites": total,
                "n_called": called,
                "n_missing": counts["missing"],
                "call_rate": called / total if total else None,
                "n_hom_ref": counts["hom_ref"],
                "n_het": counts["het"],
                "n_hom_alt": counts["hom_alt"],
                "het_hom_alt_ratio": counts["het"] / counts["hom_alt"]
                if counts["hom_alt"]
                else None,
                "n_snv": snv,
                "n_indel": indel,
                "n_transitions": ti,
                "n_transversions": tv,
                "ti_tv": ti / tv if tv else None,
                "mean_depth": sum(depths) / len(depths) if depths else None,
                "mean_genotype_quality": sum(gqs) / len(gqs) if gqs else None,
            }
        )
    return rows


def site_qc(vs: VariantSet) -> list[dict[str, Any]]:
    """Per-site QC metrics: call rate, allele counts/frequencies, HWE p-value, mean depth."""

    rows = []
    for rec in vs.records:
        n = len(rec["calls"])
        counts, an = allele_counts(rec) if n else ([0] * (len(rec["alts"]) + 1), 0)
        gcounts = genotype_counts(rec) if n else {}
        called = n - gcounts.get("missing", 0) if n else 0
        depths = [d for d in (call_depth(c) for c in rec["calls"]) if d is not None]
        rows.append(
            {
                "chrom": rec["chrom"],
                "pos": rec["pos"],
                "ref": rec["ref"],
                "alt": list(rec["alts"]),
                "n_samples": n,
                "call_rate": called / n if n else None,
                "an": an,
                "ac": counts[1:],
                "af": [c / an for c in counts[1:]] if an else None,
                "n_het": gcounts.get("het"),
                "n_hom_alt": gcounts.get("hom_alt"),
                "hwe_pvalue": site_hwe(rec) if n else None,
                "mean_depth": sum(depths) / len(depths) if depths else None,
                "qual": rec["qual"],
                "filter": rec["filter"],
            }
        )
    return rows


def summarize_set(vs: VariantSet) -> dict[str, Any]:
    """Counts by type/chromosome/filter, Ti/Tv, multi-allelic sites, and sample-level QC."""

    by_type: Counter[str] = Counter()
    by_subtype: Counter[str] = Counter()
    by_chrom: Counter[str] = Counter()
    by_filter: Counter[str] = Counter()
    ti = tv = multi = alleles = 0
    quals: list[float] = []
    indel_lengths: Counter[int] = Counter()
    for rec in vs.records:
        by_chrom[rec["chrom"]] += 1
        by_filter[rec["filter"]] += 1
        if len(rec["alts"]) > 1:
            multi += 1
        if rec["qual"] is not None:
            quals.append(rec["qual"])
        for alt in rec["alts"]:
            alleles += 1
            info = classify_alleles(rec["ref"], alt)
            by_type[info["type"]] += 1
            if info["subtype"]:
                by_subtype[f"{info['type']}:{info['subtype']}"] += 1
            if info["type"] == "indel" and info["length_change"]:
                indel_lengths[info["length_change"]] += 1
            trans = is_transition(rec["ref"], alt)
            if trans is True:
                ti += 1
            elif trans is False:
                tv += 1
    samples = sample_qc(vs) if vs.samples else []
    rates = [row["call_rate"] for row in samples if row["call_rate"] is not None]
    return {
        "n_records": len(vs.records),
        "n_alt_alleles": alleles,
        "n_multiallelic_sites": multi,
        "n_samples": len(vs.samples),
        "samples": list(vs.samples),
        "type_counts": dict(sorted(by_type.items())),
        "subtype_counts": dict(sorted(by_subtype.items())),
        "by_chromosome": dict(sorted(by_chrom.items(), key=lambda kv: chrom_key(kv[0]))),
        "filter_counts": dict(sorted(by_filter.items())),
        "transitions": ti,
        "transversions": tv,
        "ti_tv": ti / tv if tv else None,
        "mean_qual": sum(quals) / len(quals) if quals else None,
        "indel_length_counts": {str(k): v for k, v in sorted(indel_lengths.items())},
        "mean_sample_call_rate": sum(rates) / len(rates) if rates else None,
        "sample_qc": samples,
    }


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------
def _allele_sequences(
    rec: Mapping[str, Any], alleles: Sequence[int | None]
) -> tuple[str, ...] | None:
    if not is_called(alleles):
        return None
    pool = [rec["ref"], *rec["alts"]]
    out = []
    for allele in alleles:
        if allele is None or allele >= len(pool):
            return None
        out.append(pool[allele].upper())
    return tuple(sorted(out))


def compare_sets(left: VariantSet, right: VariantSet, *, match: str = "exact") -> dict[str, Any]:
    """Compare two call sets (left = truth/reference set, right = query).

    Alleles are matched on (CHROM, POS, REF, ALT) (``match='exact'``) or on (CHROM, POS) (``match='position'``).
    Reports shared / only-left / only-right alleles, precision (shared/|right|), recall (shared/|left|),
    F1, Jaccard and, for samples present in both sets, genotype concordance over the shared sites.
    """

    def keys(vs: VariantSet) -> dict[tuple[Any, ...], tuple[dict[str, Any], int]]:
        out = {}
        for rec in vs.records:
            if match == "position":
                out.setdefault((rec["chrom"], rec["pos"]), (rec, 0))
            else:
                for index, alt in enumerate(rec["alts"], 1):
                    out.setdefault(
                        (rec["chrom"], rec["pos"], rec["ref"].upper(), alt.upper()), (rec, index)
                    )
        return out

    lk, rk = keys(left), keys(right)
    shared = set(lk) & set(rk)
    only_left, only_right = set(lk) - set(rk), set(rk) - set(lk)

    def show(key: tuple[Any, ...]) -> dict[str, Any]:
        if match == "position":
            return {"chrom": key[0], "pos": key[1]}
        return {"chrom": key[0], "pos": key[1], "ref": key[2], "alt": key[3]}

    order = lambda k: (chrom_key(str(k[0])), k[1:])  # noqa: E731
    precision = len(shared) / len(rk) if rk else None
    recall = len(shared) / len(lk) if lk else None
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision and recall
        else (0.0 if precision is not None and recall is not None else None)
    )
    common = [name for name in left.samples if name in right.samples]
    concordance: dict[str, Any] = {}
    if common and match == "exact":
        left_idx = {name: left.samples.index(name) for name in common}
        right_idx = {name: right.samples.index(name) for name in common}
        right_sites: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
        for rrec in right.records:
            right_sites.setdefault((rrec["chrom"], rrec["pos"], rrec["ref"].upper()), []).append(
                rrec
            )
        for name in common:
            conc = disc = nonref_conc = nonref_total = 0
            for lrec in left.records:
                for rrec in right_sites.get((lrec["chrom"], lrec["pos"], lrec["ref"].upper()), []):
                    if not set(a.upper() for a in lrec["alts"]) & set(
                        a.upper() for a in rrec["alts"]
                    ):
                        continue
                    lg = _allele_sequences(lrec, call_gt(lrec["calls"][left_idx[name]])[0])
                    rg = _allele_sequences(rrec, call_gt(rrec["calls"][right_idx[name]])[0])
                    if lg is None or rg is None:
                        continue
                    if lg == rg:
                        conc += 1
                    else:
                        disc += 1
                    ref_pair = (lrec["ref"].upper(),) * 2
                    if lg != ref_pair or rg != ref_pair:
                        nonref_total += 1
                        nonref_conc += int(lg == rg)
            total = conc + disc
            concordance[name] = {
                "compared_genotypes": total,
                "concordant": conc,
                "discordant": disc,
                "genotype_concordance": conc / total if total else None,
                "non_reference_concordance": nonref_conc / nonref_total if nonref_total else None,
            }
    return {
        "match": match,
        "left_count": len(lk),
        "right_count": len(rk),
        "shared_count": len(shared),
        "only_left_count": len(only_left),
        "only_right_count": len(only_right),
        "precision_right_vs_left": precision,
        "recall_right_vs_left": recall,
        "f1": f1,
        "jaccard": len(shared) / len(set(lk) | set(rk)) if (lk or rk) else 1.0,
        "shared": [show(k) for k in sorted(shared, key=order)],
        "only_left": [show(k) for k in sorted(only_left, key=order)],
        "only_right": [show(k) for k in sorted(only_right, key=order)],
        "genotype_concordance": concordance,
    }
