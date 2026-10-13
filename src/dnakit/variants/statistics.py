"""Dependency-free variant summary statistics with explicit VCF-like contracts.

These functions intentionally operate on materialized, normalized variant
records. They do not call a variant caller, infer genotypes, or query a
reference database. That boundary makes the formulas reproducible in a paper
and keeps caller/database-dependent operations in ``variants.external``.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from copy import deepcopy
import math
from typing import Any

from dnakit.exceptions import ConfigurationError


def _records(value: object) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    metadata: dict[str, Any] = {}
    if isinstance(value, Mapping):
        metadata = dict(value)
        for key in ("variants", "records", "data", "items"):
            if key in value:
                value = value[key]
                break
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError("input must be a sequence of variant mappings or a mapping containing variants.")
    result: list[dict[str, Any]] = []
    for index, record in enumerate(value):
        if not isinstance(record, Mapping):
            raise ConfigurationError(f"variant record {index} must be a mapping.")
        result.append(dict(record))
    return result, metadata


def _field(record: Mapping[str, Any], *names: str, default: Any = None) -> Any:
    for name in names:
        if name in record:
            return record[name]
    return default


def _base(record: Mapping[str, Any]) -> tuple[str, int, str, list[str]]:
    chrom = str(_field(record, "chrom", "CHROM", "chromosome", default=""))
    pos = _field(record, "pos", "POS", "position")
    ref = str(_field(record, "ref", "REF", default="")).upper()
    alt_value = _field(record, "alt", "ALT", "alts", default=[])
    if isinstance(alt_value, str):
        alts = [item for item in alt_value.split(",") if item and item != "."]
    elif isinstance(alt_value, Sequence):
        alts = [str(item) for item in alt_value if str(item) and str(item) != "."]
    else:
        alts = []
    if not chrom or isinstance(pos, bool) or not isinstance(pos, int) or pos < 1:
        raise ConfigurationError("variant records require chrom and a positive integer pos.")
    if not ref or not alts:
        raise ConfigurationError("variant records require ref and at least one alt allele.")
    return chrom, pos, ref, [alt.upper() for alt in alts]


def _genotype_values(record: Mapping[str, Any]) -> list[object]:
    value = _field(record, "genotypes", "sample_genotypes", "GT", "genotype", default=None)
    if value is None:
        samples = record.get("samples")
        if isinstance(samples, Mapping):
            value = [
                sample.get("GT", sample.get("genotype"))
                for sample in samples.values()
                if isinstance(sample, Mapping)
            ]
        elif isinstance(samples, Sequence) and not isinstance(samples, (str, bytes)):
            value = [
                sample.get("GT", sample.get("genotype")) if isinstance(sample, Mapping) else sample
                for sample in samples
            ]
    if value is None:
        return []
    if isinstance(value, (str, int)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    raise ConfigurationError("genotypes must be a string, integer, or sequence.")


def _alleles(genotype: object) -> list[int | None]:
    if isinstance(genotype, bool):
        raise ConfigurationError("genotype values cannot be boolean.")
    if isinstance(genotype, int):
        return [genotype] if genotype >= 0 else [None]
    if isinstance(genotype, str):
        text = genotype.replace("|", "/")
        if text in {".", "./."}:
            return [None, None]
        values = text.split("/")
        result: list[int | None] = []
        for value in values:
            if value in {"", "."}:
                result.append(None)
            elif value.isdigit():
                result.append(int(value))
            else:
                raise ConfigurationError(f"invalid genotype allele: {value!r}.")
        return result
    if isinstance(genotype, Sequence) and not isinstance(genotype, (str, bytes)):
        result = []
        for value in genotype:
            if value is None or value == ".":
                result.append(None)
            elif isinstance(value, int) and not isinstance(value, bool) and value >= 0:
                result.append(value)
            else:
                raise ConfigurationError("genotype allele arrays must contain non-negative integers or None.")
        return result
    raise ConfigurationError("unsupported genotype representation.")


def _envelope(operation: str, result: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": operation,
        "result": dict(result),
        "provenance": {
            "implementation": "dnakit.native_variant_statistics",
            "backend": "local_dependency_free",
            "coordinate_system": "VCF-like 1-based POS",
            "scientific_boundary": "descriptive statistics only; no genotype calling or annotation",
        },
    }


def calculate_allele_frequency(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Calculate alternate-allele frequencies from explicit genotypes."""

    options = dict(options or {})
    allele_index = options.get("allele_index", 1)
    if isinstance(allele_index, bool) or not isinstance(allele_index, int) or allele_index < 1:
        raise ConfigurationError("options.allele_index must be a positive integer.")
    records, _ = _records(input)
    rows = []
    total_alt = total_called = 0
    for record in records:
        chrom, pos, ref, alts = _base(record)
        if allele_index > len(alts):
            raise ConfigurationError(f"allele_index={allele_index} is absent at {chrom}:{pos}.")
        values = _genotype_values(record)
        if not values:
            alt_count = _field(record, "alt_count", default=None)
            called_count = _field(record, "allele_count", "called_alleles", default=None)
            if not isinstance(alt_count, int) or not isinstance(called_count, int) or called_count < 1:
                raise ConfigurationError(f"{chrom}:{pos} requires explicit genotypes or alt_count/allele_count.")
        else:
            parsed = [allele for value in values for allele in _alleles(value)]
            called = [allele for allele in parsed if allele is not None]
            called_count = len(called)
            alt_count = sum(allele == allele_index for allele in called)
        frequency = alt_count / called_count if called_count else None
        total_alt += alt_count
        total_called += called_count
        rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alts, "alt_count": alt_count, "called_alleles": called_count, "frequency": frequency})
    return _envelope("calculate_allele_frequency", {"variants": rows, "total_alt_alleles": total_alt, "total_called_alleles": total_called, "aggregate_frequency": total_alt / total_called if total_called else None})


def calculate_heterozygosity(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Calculate observed heterozygosity from explicit diploid genotypes."""

    records, _ = _records(input)
    rows = []
    heterozygous = called = 0
    for record in records:
        chrom, pos, ref, alts = _base(record)
        values = _genotype_values(record)
        if not values:
            raise ConfigurationError(f"{chrom}:{pos} requires genotypes for heterozygosity.")
        site_het = site_called = 0
        for value in values:
            alleles = _alleles(value)
            if len(alleles) != 2:
                raise ConfigurationError(f"{chrom}:{pos} requires diploid genotypes.")
            if None not in alleles:
                site_called += 1
                site_het += int(alleles[0] != alleles[1])
        heterozygous += site_het
        called += site_called
        rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alts, "heterozygous_samples": site_het, "called_samples": site_called, "observed_heterozygosity": site_het / site_called if site_called else None})
    return _envelope("calculate_heterozygosity", {"variants": rows, "heterozygous_samples": heterozygous, "called_samples": called, "observed_heterozygosity": heterozygous / called if called else None})


def calculate_missingness(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Calculate missing genotype fractions without imputing or filtering data."""

    records, _ = _records(input)
    rows = []
    missing = total = 0
    for record in records:
        chrom, pos, ref, alts = _base(record)
        values = _genotype_values(record)
        if not values:
            site_missing = _field(record, "missing_count", default=None)
            site_total = _field(record, "sample_count", default=None)
            if not isinstance(site_missing, int) or not isinstance(site_total, int) or not 0 <= site_missing <= site_total:
                raise ConfigurationError(f"{chrom}:{pos} requires genotypes or missing_count/sample_count.")
        else:
            site_missing = sum(1 for value in values if any(allele is None for allele in _alleles(value)))
            site_total = len(values)
        missing += site_missing
        total += site_total
        rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alts, "missing_samples": site_missing, "sample_count": site_total, "missingness": site_missing / site_total if site_total else None})
    return _envelope("calculate_missingness", {"variants": rows, "missing_samples": missing, "sample_count": total, "missingness": missing / total if total else None})


def calculate_transition_transversion(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Count SNP transitions and transversions; exclude indels explicitly."""

    transitions = {"AG", "GA", "CT", "TC"}
    ti = tv = excluded = 0
    rows = []
    records, _ = _records(input)
    for record in records:
        chrom, pos, ref, alts = _base(record)
        for alt in alts:
            pair = ref + alt
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                excluded += 1
                classification = "excluded_non_snp"
            elif pair in transitions:
                ti += 1
                classification = "transition"
            else:
                tv += 1
                classification = "transversion"
            rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "classification": classification})
    return _envelope("calculate_transition_transversion", {"variants": rows, "transitions": ti, "transversions": tv, "excluded": excluded, "ti_tv_ratio": ti / tv if tv else None})


def calculate_variant_density(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Calculate variant density from an explicit region length in bases."""

    options = dict(options or {})
    records, metadata = _records(input)
    region_length = options.get("region_length", metadata.get("region_length"))
    if isinstance(region_length, bool) or not isinstance(region_length, int) or region_length <= 0:
        raise ConfigurationError("options.region_length or input.region_length must be a positive integer.")
    by_chromosome: Counter[str] = Counter()
    for record in records:
        chrom, _, _, _ = _base(record)
        by_chromosome[chrom] += 1
    count = len(records)
    return _envelope("calculate_variant_density", {"variant_count": count, "region_length": region_length, "density_per_base": count / region_length, "density_per_megabase": count * 1_000_000 / region_length, "by_chromosome": dict(sorted(by_chromosome.items()))})


def _text_input(value: object) -> str:
    """Convert VCF text-like input without interpreting a path as file contents."""

    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, str):
        return value
    if isinstance(value, Mapping):
        for key in ("text", "vcf", "content"):
            if key in value:
                return _text_input(value[key])
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if not all(isinstance(line, (str, bytes)) for line in value):
            raise ConfigurationError("VCF line sequences must contain only strings or bytes.")
        return "\n".join(
            line.decode("utf-8") if isinstance(line, bytes) else line for line in value
        )
    raise ConfigurationError("VCF input must be text, bytes, a line sequence, or a text mapping.")


def _parse_info(text: str) -> dict[str, object]:
    if text in {"", "."}:
        return {}
    info: dict[str, object] = {}
    for item in text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def parse_vcf_records(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Parse VCF columns into JSON-serializable records.

    The parser handles VCF text and the mandatory eight columns plus optional
    FORMAT/sample columns. It validates coordinates and column count, but does
    not normalize alleles, validate against a reference, or call variants.
    """

    settings = dict(options or {})
    strict = settings.get("strict", True)
    if not isinstance(strict, bool):
        raise ConfigurationError("options.strict must be boolean.")
    text = _text_input(input)
    headers: list[str] = []
    sample_names: list[str] = []
    records: list[dict[str, Any]] = []
    warnings: list[str] = []
    for line_number, raw_line in enumerate(text.splitlines(), 1):
        line = raw_line.rstrip("\r")
        if not line:
            continue
        if line.startswith("##"):
            headers.append(line)
            continue
        if line.startswith("#CHROM"):
            columns = line.split("\t")
            if len(columns) < 8:
                raise ConfigurationError("VCF header must contain at least eight columns.")
            sample_names = columns[9:]
            headers.append(line)
            continue
        if line.startswith("#"):
            warnings.append(f"ignored comment line {line_number}")
            continue
        columns = line.split("\t")
        if len(columns) < 8:
            message = f"line {line_number} has {len(columns)} VCF columns; at least 8 are required"
            if strict:
                raise ConfigurationError(message)
            warnings.append(message)
            continue
        chrom, pos_text, identifier, ref, alt_text, qual_text, filter_text, info_text = columns[:8]
        try:
            pos = int(pos_text)
        except ValueError as exc:
            if strict:
                raise ConfigurationError(f"line {line_number} has a non-integer POS") from exc
            warnings.append(f"line {line_number} skipped: non-integer POS")
            continue
        if not chrom or pos < 1 or not ref or not alt_text or alt_text == ".":
            message = f"line {line_number} has invalid CHROM/POS/REF/ALT"
            if strict:
                raise ConfigurationError(message)
            warnings.append(message)
            continue
        if qual_text in {"", "."}:
            qual: float | None = None
        else:
            try:
                qual = float(qual_text)
            except ValueError as exc:
                if strict:
                    raise ConfigurationError(f"line {line_number} has a non-numeric QUAL") from exc
                warnings.append(f"line {line_number} skipped: non-numeric QUAL")
                continue
            if not math.isfinite(qual):
                raise ConfigurationError(f"line {line_number} has a non-finite QUAL")
        alts = [alt for alt in alt_text.split(",") if alt and alt != "."]
        if not alts:
            raise ConfigurationError(f"line {line_number} has no alternate allele")
        record: dict[str, Any] = {
            "chrom": chrom,
            "pos": pos,
            "id": None if identifier == "." else identifier,
            "ref": ref,
            "alt": alts,
            "qual": qual,
            "filter": filter_text,
            "info": _parse_info(info_text),
        }
        if len(columns) > 8:
            format_keys = [] if columns[8] == "." else columns[8].split(":")
            samples = []
            for index, sample_text in enumerate(columns[9:]):
                fields = sample_text.split(":")
                samples.append(
                    {
                        key: fields[field_index] if field_index < len(fields) else None
                        for field_index, key in enumerate(format_keys)
                    }
                )
            record["format"] = format_keys
            record["samples"] = samples
            if format_keys and "GT" in format_keys:
                gt_index = format_keys.index("GT")
                record["genotypes"] = [
                    sample_text.split(":")[gt_index]
                    if gt_index < len(sample_text.split(":"))
                    else "."
                    for sample_text in columns[9:]
                ]
        records.append(record)
    return _envelope(
        "parse_vcf_records",
        {
            "headers": headers,
            "samples": sample_names,
            "records": records,
            "record_count": len(records),
            "warnings": warnings,
        },
    )


def _record_list(value: object) -> list[dict[str, Any]]:
    records, _ = _records(value)
    return records


def _remap_genotype(genotype: object, selected_alt: int) -> object:
    """Map a multiallelic GT to a biallelic GT for one selected ALT."""

    if not isinstance(genotype, str) or genotype in {"", "."}:
        return genotype
    separator = "|" if "|" in genotype else "/"
    values = genotype.replace("|", "/").split("/")
    remapped: list[str] = []
    for value in values:
        if value in {"", "."}:
            remapped.append(".")
        elif value.isdigit():
            remapped.append("1" if int(value) == selected_alt else "0")
        else:
            raise ConfigurationError(f"invalid genotype allele: {value!r}.")
    return separator.join(remapped)


def split_multiallelic(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Split ALT lists into biallelic records and optionally remap GT fields."""

    settings = dict(options or {})
    remap = settings.get("remap_genotypes", True)
    if not isinstance(remap, bool):
        raise ConfigurationError("options.remap_genotypes must be boolean.")
    records = _record_list(input)
    split_records: list[dict[str, Any]] = []
    for record in records:
        _, _, _, alts = _base(record)
        for selected_index, alt in enumerate(alts, 1):
            row = deepcopy(record)
            row["alt"] = alt
            row["original_alt_index"] = selected_index
            if remap and "genotypes" in row:
                row["genotypes"] = [
                    _remap_genotype(value, selected_index) for value in _genotype_values(record)
                ]
            samples = row.get("samples")
            if remap and isinstance(samples, list):
                row["samples"] = [
                    {
                        **sample,
                        "GT": _remap_genotype(sample.get("GT"), selected_index),
                    }
                    if isinstance(sample, Mapping) and "GT" in sample
                    else sample
                    for sample in samples
                ]
            split_records.append(row)
    return _envelope(
        "split_multiallelic",
        {
            "records": split_records,
            "input_record_count": len(records),
            "output_record_count": len(split_records),
            "split_record_count": max(0, len(split_records) - len(records)),
            "remap_genotypes": remap,
        },
    )


def count(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Count records and allele-level variant classes from explicit records."""

    records = _record_list(input)
    chromosomes: Counter[str] = Counter()
    classes: Counter[str] = Counter()
    alternate_alleles = 0
    for record in records:
        chrom, _, ref, alts = _base(record)
        chromosomes[chrom] += 1
        alternate_alleles += len(alts)
        for alt in alts:
            if len(ref) == len(alt) == 1:
                variant_class = "snp" if ref != alt else "reference_match"
            elif len(ref) == len(alt):
                variant_class = "mnp"
            elif len(ref) < len(alt):
                variant_class = "insertion"
            elif len(ref) > len(alt):
                variant_class = "deletion"
            else:
                variant_class = "other"
            classes[variant_class] += 1
    return _envelope(
        "count",
        {
            "record_count": len(records),
            "alternate_allele_count": alternate_alleles,
            "by_chromosome": dict(sorted(chromosomes.items())),
            "allele_type_counts": dict(sorted(classes.items())),
        },
    )


def _variant_keys(value: object) -> tuple[set[tuple[str, int, str, str]], list[dict[str, Any]]]:
    records = _record_list(value)
    keyed: set[tuple[str, int, str, str]] = set()
    for record in records:
        chrom, pos, ref, alts = _base(record)
        keyed.update((chrom, pos, ref, alt) for alt in alts)
    return keyed, records


def compare(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Compare two explicit variant sets by exact CHROM/POS/REF/ALT identity."""

    if not isinstance(input, Mapping) or "left" not in input or "right" not in input:
        raise ConfigurationError("compare input must contain left and right variant collections.")
    left, left_records = _variant_keys(input["left"])
    right, right_records = _variant_keys(input["right"])

    def serialise(keys: set[tuple[str, int, str, str]]) -> list[dict[str, Any]]:
        return [
            {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
            for chrom, pos, ref, alt in sorted(keys, key=lambda item: (item[0], item[1], item[2], item[3]))
        ]

    shared = left & right
    union = left | right
    return _envelope(
        "compare",
        {
            "left_record_count": len(left_records),
            "right_record_count": len(right_records),
            "left_allele_count": len(left),
            "right_allele_count": len(right),
            "shared_allele_count": len(shared),
            "only_left": serialise(left - right),
            "only_right": serialise(right - left),
            "shared": serialise(shared),
            "jaccard": len(shared) / len(union) if union else 1.0,
        },
    )


def calculate_minor_allele_frequency(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Calculate MAF from explicitly called genotypes, including the reference allele."""

    records = _record_list(input)
    rows: list[dict[str, Any]] = []
    for record in records:
        chrom, pos, ref, alts = _base(record)
        values = _genotype_values(record)
        if not values:
            raise ConfigurationError(f"{chrom}:{pos} requires explicit genotypes for MAF.")
        counts: Counter[int] = Counter()
        for value in values:
            for allele in _alleles(value):
                if allele is not None:
                    if allele > len(alts):
                        raise ConfigurationError(f"{chrom}:{pos} genotype allele exceeds ALT count.")
                    counts[allele] += 1
        called = sum(counts.values())
        if not called:
            rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alts, "called_alleles": 0, "allele_frequencies": {}, "minor_allele_index": None, "minor_allele": None, "maf": None})
            continue
        # Include declared ALT alleles with zero observations.  Otherwise an
        # all-reference site would incorrectly report MAF=1.0 by considering
        # only alleles seen in the genotype calls.
        frequencies = {
            str(index): counts.get(index, 0) / called
            for index in range(len(alts) + 1)
        }
        minor_index, maf = min(
            frequencies.items(), key=lambda item: (item[1], int(item[0]))
        )
        minor_index = int(minor_index)
        minor_value = ref if minor_index == 0 else alts[minor_index - 1]
        rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alts, "called_alleles": called, "allele_frequencies": frequencies, "minor_allele_index": minor_index, "minor_allele": minor_value, "maf": maf})
    valid = [row["maf"] for row in rows if row["maf"] is not None]
    return _envelope(
        "calculate_minor_allele_frequency",
        {"variants": rows, "variant_count": len(rows), "mean_maf": sum(valid) / len(valid) if valid else None},
    )


def _quality(record: Mapping[str, Any]) -> float | None:
    value = _field(record, "qual", "QUAL", "quality", default=None)
    if value is None or value == ".":
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise ConfigurationError("QUAL must be a finite number or null.")
    return float(value)


def _threshold(settings: Mapping[str, object], name: str, *, default: object = None) -> float:
    value = settings.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise ConfigurationError(f"options.{name} must be a finite number.")
    return float(value)


def _filter_envelope(operation: str, records: list[dict[str, Any]], kept: list[dict[str, Any]], metrics: list[dict[str, Any]], predicate: Mapping[str, Any]) -> dict[str, Any]:
    return _envelope(
        operation,
        {
            "variants": kept,
            "metrics": metrics,
            "input_record_count": len(records),
            "output_record_count": len(kept),
            "removed_record_count": len(records) - len(kept),
            "predicate": dict(predicate),
        },
    )


def filter_variants_by_quality(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Keep records whose explicit QUAL is at least a required threshold."""

    settings = dict(options or {})
    minimum = _threshold(settings, "min_quality")
    keep_missing = settings.get("keep_missing", False)
    if not isinstance(keep_missing, bool):
        raise ConfigurationError("options.keep_missing must be boolean.")
    records = _record_list(input)
    kept: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []
    for record in records:
        quality = _quality(record)
        keep = quality is not None and quality >= minimum or quality is None and keep_missing
        metrics.append({"chrom": str(_field(record, "chrom", "CHROM", default="")), "pos": _field(record, "pos", "POS", default=None), "quality": quality, "kept": keep})
        if keep:
            kept.append(record)
    return _filter_envelope("filter_variants_by_quality", records, kept, metrics, {"min_quality": minimum, "keep_missing": keep_missing})


def _numeric_record_value(record: Mapping[str, Any], *names: str) -> float | None:
    value = _field(record, *names, default=None)
    if value is None or value == ".":
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise ConfigurationError(f"{names[0]} must be a finite number or null.")
    return float(value)


def filter_variants_by_depth(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Keep records whose explicit read depth is within a requested range."""

    settings = dict(options or {})
    minimum = _threshold(settings, "min_depth", default=0.0)
    maximum = _threshold(settings, "max_depth") if "max_depth" in settings else float("inf")
    keep_missing = settings.get("keep_missing", False)
    if not isinstance(keep_missing, bool):
        raise ConfigurationError("options.keep_missing must be boolean.")
    if minimum < 0 or maximum < minimum:
        raise ConfigurationError("depth thresholds must satisfy 0 <= min_depth <= max_depth.")
    records = _record_list(input)
    kept: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []
    for record in records:
        depth = _numeric_record_value(record, "depth", "DP", "read_depth")
        keep = keep_missing if depth is None else minimum <= depth <= maximum
        metrics.append({"chrom": str(_field(record, "chrom", "CHROM", default="")), "pos": _field(record, "pos", "POS", default=None), "depth": depth, "kept": keep})
        if keep:
            kept.append(record)
    return _filter_envelope("filter_variants_by_depth", records, kept, metrics, {"min_depth": minimum, "max_depth": maximum, "keep_missing": keep_missing})


def _genotype_qualities(record: Mapping[str, Any]) -> list[float]:
    value = _field(record, "genotype_qualities", "sample_gq", "GQ", "genotype_quality", default=None)
    if value is None:
        samples = record.get("samples")
        if isinstance(samples, Mapping):
            value = [sample.get("GQ", sample.get("genotype_quality")) for sample in samples.values() if isinstance(sample, Mapping)]
        elif isinstance(samples, Sequence) and not isinstance(samples, (str, bytes)):
            value = [sample.get("GQ", sample.get("genotype_quality")) if isinstance(sample, Mapping) else sample for sample in samples]
    if value is None:
        return []
    values = [value] if isinstance(value, (int, float)) and not isinstance(value, bool) else value
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError("genotype quality values must be numeric or a sequence of numeric values.")
    result: list[float] = []
    for item in values:
        if item is None or item == ".":
            continue
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(float(item)):
            raise ConfigurationError("genotype quality values must be finite numbers or null.")
        result.append(float(item))
    return result


def filter_variants_by_genotype_quality(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Keep records whose explicit called-sample GQ values meet a threshold."""

    settings = dict(options or {})
    minimum = _threshold(settings, "min_genotype_quality", default=settings.get("min_gq", None))
    keep_missing = settings.get("keep_missing", False)
    if not isinstance(keep_missing, bool):
        raise ConfigurationError("options.keep_missing must be boolean.")
    records = _record_list(input)
    kept: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []
    for record in records:
        qualities = _genotype_qualities(record)
        minimum_observed = min(qualities) if qualities else None
        keep = keep_missing if minimum_observed is None else minimum_observed >= minimum
        metrics.append({"chrom": str(_field(record, "chrom", "CHROM", default="")), "pos": _field(record, "pos", "POS", default=None), "genotype_qualities": qualities, "minimum_observed": minimum_observed, "kept": keep})
        if keep:
            kept.append(record)
    return _filter_envelope("filter_variants_by_genotype_quality", records, kept, metrics, {"min_genotype_quality": minimum, "keep_missing": keep_missing})


def _missingness(record: Mapping[str, Any]) -> float:
    values = _genotype_values(record)
    if values:
        return sum(any(allele is None for allele in _alleles(value)) for value in values) / len(values)
    missing = _field(record, "missing_count", default=None)
    total = _field(record, "sample_count", default=None)
    if isinstance(missing, int) and isinstance(total, int) and total > 0 and 0 <= missing <= total:
        return missing / total
    chrom = _field(record, "chrom", "CHROM", default="?")
    pos = _field(record, "pos", "POS", default="?")
    raise ConfigurationError(f"{chrom}:{pos} requires genotypes or missing_count/sample_count.")


def filter_variants_by_missingness(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Keep records at or below an explicit missing-genotype fraction."""

    settings = dict(options or {})
    maximum = _threshold(settings, "max_missingness")
    if not 0 <= maximum <= 1:
        raise ConfigurationError("options.max_missingness must be between 0 and 1.")
    records = _record_list(input)
    kept: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []
    for record in records:
        missingness = _missingness(record)
        keep = missingness <= maximum
        metrics.append({"chrom": str(_field(record, "chrom", "CHROM", default="")), "pos": _field(record, "pos", "POS", default=None), "missingness": missingness, "kept": keep})
        if keep:
            kept.append(record)
    return _filter_envelope("filter_variants_by_missingness", records, kept, metrics, {"max_missingness": maximum})


def _alternate_frequency(record: Mapping[str, Any], allele_index: int) -> float:
    alts = _base(record)[3]
    if allele_index > len(alts):
        raise ConfigurationError(f"allele_index={allele_index} is absent at {_field(record, 'chrom', 'CHROM', default='?')}:{_field(record, 'pos', 'POS', default='?')}.")
    values = _genotype_values(record)
    if values:
        alleles = [allele for value in values for allele in _alleles(value) if allele is not None]
        if not alleles:
            raise ConfigurationError("cannot calculate allele frequency with no called alleles.")
        return sum(allele == allele_index for allele in alleles) / len(alleles)
    alt_count = _field(record, "alt_count", default=None)
    total = _field(record, "allele_count", "called_alleles", default=None)
    if isinstance(alt_count, int) and isinstance(total, int) and total > 0 and 0 <= alt_count <= total:
        return alt_count / total
    raise ConfigurationError("records require explicit genotypes or alt_count/allele_count.")


def filter_variants_by_allele_frequency(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Keep records whose explicitly calculated alternate frequency is in range."""

    settings = dict(options or {})
    minimum = _threshold(settings, "min_frequency", default=0.0)
    maximum = _threshold(settings, "max_frequency", default=1.0)
    allele_index = settings.get("allele_index", 1)
    if isinstance(allele_index, bool) or not isinstance(allele_index, int) or allele_index < 1:
        raise ConfigurationError("options.allele_index must be a positive integer.")
    if not 0 <= minimum <= maximum <= 1:
        raise ConfigurationError("frequency thresholds must satisfy 0 <= min_frequency <= max_frequency <= 1.")
    records = _record_list(input)
    kept: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []
    for record in records:
        frequency = _alternate_frequency(record, allele_index)
        keep = minimum <= frequency <= maximum
        metrics.append({"chrom": str(_field(record, "chrom", "CHROM", default="")), "pos": _field(record, "pos", "POS", default=None), "allele_index": allele_index, "frequency": frequency, "kept": keep})
        if keep:
            kept.append(record)
    return _filter_envelope("filter_variants_by_allele_frequency", records, kept, metrics, {"allele_index": allele_index, "min_frequency": minimum, "max_frequency": maximum})


__all__ = [
    "calculate_allele_frequency",
    "calculate_heterozygosity",
    "calculate_missingness",
    "calculate_transition_transversion",
    "calculate_variant_density",
    "parse_vcf_records",
    "split_multiallelic",
    "count",
    "compare",
    "calculate_minor_allele_frequency",
    "filter_variants_by_quality",
    "filter_variants_by_depth",
    "filter_variants_by_genotype_quality",
    "filter_variants_by_missingness",
    "filter_variants_by_allele_frequency",
]
