"""varcode variant-collection, effect-collection, genotype, phasing and germline operations (integration pack P32).

Every function works on a local reference (a GTF and a genome FASTA, indexed once and cached by
``varcode_gap_functions``) and on local variant files. Variants are loaded from VCF files, MAF files, CSV files
written by varcode, JSON variant collections and explicit records (``contig``, 1-based ``start``, ``ref``, ``alt``).
Contig names are kept exactly as they are in the files. Results are JSON-serialisable dictionaries;
collections are returned as lists of rows. Functions that merge several varcode methods name the values of the
selecting parameter in their docstrings. Tested with varcode 10.11.10 and pysam 0.23.

Missing packages raise ``BackendUnavailableError`` with the ``pip install "dnakit[bio-varcode]"`` hint (a local
reference also needs ``dnakit[bio-pyensembl]``, RNA read phasing needs ``dnakit[bio-pysam]``).
"""

from __future__ import annotations

import os
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab.varcode_gap_functions import (
    _VARCODE_EXTRA,
    _candidate_record,
    _document,
    _effect_record,
    _genome,
    _identifier,
    _int,
    _path,
    _text,
    _variant,
    _variant_record,
)

_MAX_ROWS = 100_000
_SOURCE_KEYS = ("vcf", "maf", "csv", "json", "variants")
_BAM_EXTRA = "bio-pysam"


# --- shared loading ----------------------------------------------------------------------------


def _vc() -> Any:
    return require("varcode", _VARCODE_EXTRA)


def _new_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, os.PathLike)) or not str(value):
        raise ConfigurationError(f"{name} must be a path.")
    target = Path(value).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"{name} must not exist yet: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {target.parent}")
    return target


def _paths(values: object, name: str) -> list[Path]:
    if values is None:
        return []
    if isinstance(values, (str, bytes, os.PathLike)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of file paths.")
    return [_path(item, f"{name} entry") for item in values]


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}.")
    return value


def _guard(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run a varcode call; upstream input errors become ConfigurationError."""
    try:
        return function(*args, **kwargs)
    except (ValueError, KeyError, TypeError, AttributeError) as exc:
        raise ConfigurationError(f"varcode could not complete the request: {exc}") from exc


def _collection(
    genome: Any,
    vcf_paths: object = (),
    maf_paths: object = (),
    csv_paths: object = (),
    json_paths: object = (),
    variants: object = (),
    *,
    only_passing: bool = True,
    required: bool = True,
) -> Any:
    """Load and union the variants of all given sources into one ``VariantCollection``."""
    varcode = _vc()
    _flag(only_passing, "only_passing")
    parts: list[Any] = []
    for path in _paths(vcf_paths, "vcf_paths"):
        parts.append(
            _guard(
                varcode.load_vcf,
                str(path),
                genome=genome,
                only_passing=only_passing,
                parse_structural_variants=True,
                convert_ucsc_contig_names=False,
            )
        )
    for path in _paths(maf_paths, "maf_paths"):
        parts.append(_guard(require("varcode.maf", _VARCODE_EXTRA).load_maf, str(path)))
    for path in _paths(csv_paths, "csv_paths"):
        parts.append(_guard(varcode.VariantCollection.from_csv, str(path), genome=genome))
    for path in _paths(json_paths, "json_paths"):
        parts.append(_guard(varcode.VariantCollection.from_json, path.read_text(encoding="utf-8")))
    if variants:
        if isinstance(variants, (str, bytes)) or not isinstance(variants, Sequence):
            raise ConfigurationError("variants must be a list of mappings.")
        parts.append(varcode.VariantCollection([_variant(genome, item) for item in variants]))
    if not parts:
        if required:
            raise ConfigurationError(
                "No variants were given: pass vcf_paths, maf_paths, csv_paths, json_paths or variants."
            )
        return None
    return parts[0] if len(parts) == 1 else varcode.VariantCollection.union(*parts)


def _sources(genome: Any, spec: object, name: str, *, only_passing: bool) -> Any:
    if not isinstance(spec, Mapping):
        raise ConfigurationError(
            f"{name} must be a mapping with any of: {', '.join(_SOURCE_KEYS)}."
        )
    unknown = sorted(set(spec) - set(_SOURCE_KEYS))
    if unknown:
        raise ConfigurationError(f"{name} has unknown keys: {', '.join(unknown)}.")
    return _collection(
        genome,
        spec.get("vcf", ()),
        spec.get("maf", ()),
        spec.get("csv", ()),
        spec.get("json", ()),
        spec.get("variants", ()),
        only_passing=only_passing,
    )


def _row(variant: Any) -> dict[str, Any]:
    return {
        "contig": str(variant.contig),
        "start": int(variant.start),
        "ref": str(variant.ref),
        "alt": str(variant.alt),
    }


def _rows(collection: Any) -> list[dict[str, Any]]:
    rows = [_row(item) for item in collection]
    if len(rows) > _MAX_ROWS:
        raise ConfigurationError(f"The result exceeds the limit of {_MAX_ROWS} variants.")
    return rows


def _table(frame: Any) -> list[dict[str, Any]]:
    records = frame.to_dict(orient="records")
    if len(records) > _MAX_ROWS:
        raise ConfigurationError(f"The table exceeds the limit of {_MAX_ROWS} rows.")
    cleaned: list[dict[str, Any]] = []
    for record in records:
        row: dict[str, Any] = {}
        for key, value in record.items():
            if hasattr(value, "item"):
                value = value.item()
            if isinstance(value, float) and value != value:
                value = None
            row[str(key)] = value
        cleaned.append(row)
    return cleaned


def _effects(collection: Any, *, skip_errors: bool = False) -> Any:
    _flag(skip_errors, "skip_errors")
    return _guard(collection.effects, raise_on_error=not skip_errors)


def _effect_row(effect: Any) -> dict[str, Any]:
    from varcode.effects import effect_priority

    row = _effect_record(effect)
    row["priority"] = int(effect_priority(effect))
    transcript = getattr(effect, "transcript", None)
    row["gene_id"] = getattr(transcript, "gene_id", None) if transcript is not None else None
    return row


def _base(gtf: object, genome_fasta: object, cache_dir: object) -> tuple[Any, dict[str, Any]]:
    return _genome(gtf, genome_fasta, cache_dir)


def _load(
    gtf: object,
    genome_fasta: object,
    cache_dir: object,
    vcf_paths: object,
    maf_paths: object,
    csv_paths: object,
    json_paths: object,
    variants: object,
    only_passing: bool,
) -> tuple[Any, Any, dict[str, Any]]:
    genome, info = _base(gtf, genome_fasta, cache_dir)
    collection = _collection(
        genome, vcf_paths, maf_paths, csv_paths, json_paths, variants, only_passing=only_passing
    )
    return genome, collection, info


# --- variant collections -----------------------------------------------------------------------


def varcode_read_variant_csv(
    gtf: str,
    genome_fasta: str,
    csv_path: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a variant collection from a CSV file written by varcode (``VariantCollection.from_csv``).

    The CSV has the columns ``chr,start,ref,alt`` (optionally gene columns) after ``#`` comment lines. The variants
    are attached to the local reference built from ``gtf`` and ``genome_fasta``.
    """

    genome, info = _base(gtf, genome_fasta, cache_dir)
    collection = _collection(genome, csv_paths=[csv_path])
    rows = _rows(collection)
    return _document(
        "read_variant_csv", variants=rows, count=len(rows), cache_directory=info["cache_directory"]
    )


def varcode_variant_table(
    gtf: str,
    genome_fasta: str,
    *,
    vcf_paths: Sequence[str] = (),
    maf_paths: Sequence[str] = (),
    csv_paths: Sequence[str] = (),
    json_paths: Sequence[str] = (),
    variants: Sequence[Mapping[str, Any]] = (),
    only_passing: bool = True,
    output_csv: str | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Table of a variant collection with the genes each variant overlaps, optionally written to a CSV file.

    Without ``output_csv`` the table comes from ``VariantCollection.to_dataframe`` (columns chr, start, ref, alt,
    gene_name, gene_id), which is also the listing of the ``varcode-genes`` command. With ``output_csv`` the collection
    is also written by ``VariantCollection.to_csv`` (an existing file is refused). Sources are VCF, MAF, varcode
    CSV and JSON files and explicit records.
    """

    target = _new_file(output_csv, "output_csv") if output_csv is not None else None
    _genome, collection, info = _load(
        gtf,
        genome_fasta,
        cache_dir,
        vcf_paths,
        maf_paths,
        csv_paths,
        json_paths,
        variants,
        only_passing,
    )
    table = _table(_guard(collection.to_dataframe))
    if target is not None:
        _guard(collection.to_csv, str(target))
    return _document(
        "variant_table",
        variants=table,
        count=len(table),
        columns=list(table[0]) if table else [],
        output_csv=str(target) if target is not None else None,
        cache_directory=info["cache_directory"],
    )


def varcode_variant_set_operation(
    gtf: str,
    genome_fasta: str,
    operation: str,
    first: Mapping[str, Any],
    second: Mapping[str, Any],
    *,
    only_passing: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine two variant call sets: ``operation`` is ``union``, ``intersection`` or ``difference`` (first minus second).

    ``first`` and ``second`` are mappings with any of the keys ``vcf``, ``maf``, ``csv``, ``json`` (lists of file
    paths) and ``variants`` (list of records). They map to ``VariantCollection.union``, ``.intersection`` and
    ``.difference``.
    """

    chosen = _choice(operation, "operation", ("union", "intersection", "difference"))
    genome, info = _base(gtf, genome_fasta, cache_dir)
    left = _sources(genome, first, "first", only_passing=only_passing)
    right = _sources(genome, second, "second", only_passing=only_passing)
    result = _guard(getattr(left, chosen), right)
    rows = _rows(result)
    return _document(
        "variant_set_operation",
        set_operation=chosen,
        first_count=len(left),
        second_count=len(right),
        variants=rows,
        count=len(rows),
        cache_directory=info["cache_directory"],
    )


def _expression_map(value: object) -> dict[str, float]:
    if not isinstance(value, Mapping) or not value:
        raise ConfigurationError(
            "expression must be a non-empty mapping from transcript or gene id to a number."
        )
    cleaned: dict[str, float] = {}
    for key, number in value.items():
        if (
            isinstance(number, bool)
            or not isinstance(number, (int, float))
            or number != number
            or number < 0
        ):
            raise ConfigurationError("expression values must be non-negative numbers.")
        cleaned[_identifier(key, "expression id")] = float(number)
    return cleaned


def varcode_filter_by_expression(
    gtf: str,
    genome_fasta: str,
    expression: Mapping[str, float],
    *,
    level: str = "transcript",
    target: str = "variants",
    min_expression_value: float = 0.0,
    vcf_paths: Sequence[str] = (),
    maf_paths: Sequence[str] = (),
    csv_paths: Sequence[str] = (),
    json_paths: Sequence[str] = (),
    variants: Sequence[Mapping[str, Any]] = (),
    only_passing: bool = True,
    skip_errors: bool = False,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the variants or effects whose transcript or gene is expressed above a threshold.

    ``expression`` maps transcript ids (``level='transcript'``, ``filter_by_transcript_expression``) or gene ids
    (``level='gene'``, ``filter_by_gene_expression``) to a level such as TPM. ``target='variants'`` filters the
    ``VariantCollection``; ``target='effects'`` annotates the variants first and filters the ``EffectCollection``.
    Items with a level at or below ``min_expression_value`` are removed.
    """

    chosen_level = _choice(level, "level", ("transcript", "gene"))
    chosen_target = _choice(target, "target", ("variants", "effects"))
    if isinstance(min_expression_value, bool) or not isinstance(min_expression_value, (int, float)):
        raise ConfigurationError("min_expression_value must be a number.")
    levels = _expression_map(expression)
    _genome, collection, info = _load(
        gtf,
        genome_fasta,
        cache_dir,
        vcf_paths,
        maf_paths,
        csv_paths,
        json_paths,
        variants,
        only_passing,
    )
    subject = (
        collection if chosen_target == "variants" else _effects(collection, skip_errors=skip_errors)
    )
    kept = _guard(
        getattr(subject, f"filter_by_{chosen_level}_expression"),
        levels,
        min_expression_value=float(min_expression_value),
    )
    if chosen_target == "variants":
        items: list[dict[str, Any]] = _rows(kept)
    else:
        items = [_effect_row(effect) for effect in kept]
    return _document(
        "filter_by_expression",
        level=chosen_level,
        target=chosen_target,
        min_expression_value=float(min_expression_value),
        input_count=len(subject),
        items=items,
        count=len(items),
        cache_directory=info["cache_directory"],
    )


_SAMPLE_MODES = ("carried", "heterozygous", "homozygous_alt")


def varcode_variants_in_sample(
    gtf: str,
    genome_fasta: str,
    vcf_paths: Sequence[str],
    sample: str,
    *,
    genotype: str = "carried",
    only_passing: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Variants of a VCF sample, optionally restricted by zygosity.

    ``genotype='carried'`` returns the variants the sample carries (``for_sample``), ``'heterozygous'`` those called
    heterozygous (``heterozygous_in``) and ``'homozygous_alt'`` those homozygous for the alternate allele
    (``homozygous_alt_in``).
    """

    mode = _choice(genotype, "genotype", _SAMPLE_MODES)
    name = _text(sample, "sample")
    _genome, collection, info = _load(
        gtf, genome_fasta, cache_dir, vcf_paths, (), (), (), (), only_passing
    )
    if name not in list(collection.samples):
        raise ConfigurationError(
            f"Sample {name} is not in the VCF files; available: {sorted(collection.samples)}."
        )
    method = {
        "carried": "for_sample",
        "heterozygous": "heterozygous_in",
        "homozygous_alt": "homozygous_alt_in",
    }[mode]
    rows = _rows(_guard(getattr(collection, method), name))
    return _document(
        "variants_in_sample",
        sample=name,
        genotype=mode,
        variants=rows,
        count=len(rows),
        cache_directory=info["cache_directory"],
    )


def _genotype_record(call: Any, alt_index: int) -> dict[str, Any]:
    zygosity = call.zygosity_for_alt(alt_index)
    return {
        "raw_gt": call.raw_gt,
        "alleles": list(call.alleles),
        "phased": bool(call.phased),
        "phase_set": call.phase_set,
        "allele_depths": list(call.allele_depths) if call.allele_depths is not None else None,
        "total_depth": call.total_depth,
        "genotype_quality": call.genotype_quality,
        "alt_index": alt_index,
        "zygosity": str(getattr(zygosity, "value", zygosity)),
        "carries_alt": bool(call.carries_alt(alt_index)),
        "copies_of_alt": int(call.copies_of_alt(alt_index)),
        "depth_for_alt": call.depth_for_alt(alt_index),
    }


def varcode_sample_genotype(
    gtf: str,
    genome_fasta: str,
    vcf_paths: Sequence[str],
    sample: str,
    variant: Mapping[str, Any],
    *,
    alt_index: int = 1,
    only_passing: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genotype and zygosity of one sample at one variant of a VCF (``VariantCollection.genotype`` and ``.zygosity``).

    The result holds the call (GT, alleles, phase, phase set, depths) and, for allele ``alt_index`` (1 = first ALT),
    the zygosity (``Genotype.zygosity_for_alt``), whether the call carries it, its copy number and its read depth.
    """

    name = _text(sample, "sample")
    index = _int(alt_index, "alt_index", minimum=1, maximum=100)
    genome, collection, info = _load(
        gtf, genome_fasta, cache_dir, vcf_paths, (), (), (), (), only_passing
    )
    item = _variant(genome, variant)
    match = next((entry for entry in collection if _row(entry) == _row(item)), None)
    if match is None:
        raise ConfigurationError("The variant is not in the VCF files.")
    if name not in list(collection.samples):
        raise ConfigurationError(
            f"Sample {name} is not in the VCF files; available: {sorted(collection.samples)}."
        )
    call = _guard(collection.genotype, match, name)
    record = _genotype_record(call, index)
    record["collection_zygosity"] = str(getattr(collection.zygosity(match, name), "value", None))
    return _document(
        "sample_genotype",
        sample=name,
        variant=_variant_record(item),
        genotype=record,
        cache_directory=info["cache_directory"],
    )


def varcode_genotype_alt_summary(
    gt: str,
    *,
    allele_depths: Sequence[int] | None = None,
    alt_index: int = 1,
    phase_set: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise one VCF genotype for an alternate allele without any reference (``varcode.Genotype``).

    ``gt`` is a VCF GT string such as ``0/1`` or ``1|1``; ``allele_depths`` is the AD list. Returns
    ``carries_alt`` (``Genotype.carries_alt``), ``copies_of_alt`` (``copies_of_alt``), ``depth_for_alt``
    (``depth_for_alt``, read depth of the alt allele from AD) and ``zygosity`` (``zygosity_for_alt``) for allele
    ``alt_index`` (1 = first ALT), plus ploidy and phase.
    """

    genotype_class = _vc().Genotype
    text = _text(gt, "gt")
    index = _int(alt_index, "alt_index", minimum=1, maximum=100)
    info: dict[str, Any] = {"GT": text}
    if allele_depths is not None:
        if isinstance(allele_depths, (str, bytes)) or not isinstance(allele_depths, Sequence):
            raise ConfigurationError("allele_depths must be a list of integers.")
        info["AD"] = [_int(item, "allele depth", maximum=10**9) for item in allele_depths]
    if phase_set is not None:
        info["PS"] = _int(phase_set, "phase_set", minimum=0, maximum=10**9)
    call = _guard(genotype_class.from_sample_info, info)
    record = _genotype_record(call, index)
    record.update({"ploidy": call.ploidy, "is_called": bool(call.is_called)})
    return _document("genotype_alt_summary", genotype=record)


# --- effect collections ------------------------------------------------------------------------


def varcode_effect_table(
    gtf: str,
    genome_fasta: str,
    *,
    vcf_paths: Sequence[str] = (),
    maf_paths: Sequence[str] = (),
    csv_paths: Sequence[str] = (),
    json_paths: Sequence[str] = (),
    variants: Sequence[Mapping[str, Any]] = (),
    only_passing: bool = True,
    skip_errors: bool = False,
    output_csv: str | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate variants and return the effects as a table (``EffectCollection.to_dataframe``), optionally also as CSV.

    Columns: variant, contig, start, ref, alt, gene and transcript ids and names, effect_type and effect description.
    With ``output_csv`` the effects are also written by ``EffectCollection.to_csv`` (an existing file is refused).
    """

    target = _new_file(output_csv, "output_csv") if output_csv is not None else None
    _genome, collection, info = _load(
        gtf,
        genome_fasta,
        cache_dir,
        vcf_paths,
        maf_paths,
        csv_paths,
        json_paths,
        variants,
        only_passing,
    )
    effects = _effects(collection, skip_errors=skip_errors)
    table = _table(_guard(effects.to_dataframe))
    if target is not None:
        _guard(effects.to_csv, str(target))
    return _document(
        "effect_table",
        effects=table,
        count=len(table),
        output_csv=str(target) if target is not None else None,
        cache_directory=info["cache_directory"],
    )


def varcode_read_effect_csv(
    gtf: str,
    genome_fasta: str,
    csv_path: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read an effect CSV written by varcode and re-annotate its variants on the local reference (``EffectCollection.from_csv``)."""

    genome, info = _base(gtf, genome_fasta, cache_dir)
    path = _path(csv_path, "csv_path")
    effects = _guard(_vc().EffectCollection.from_csv, str(path), genome=genome)
    rows = [_effect_row(effect) for effect in effects]
    return _document(
        "read_effect_csv", effects=rows, count=len(rows), cache_directory=info["cache_directory"]
    )


def varcode_effect_priority(
    effect_classes: Sequence[str] | None = None,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Consequence severity rank of varcode effect classes (``varcode.effect_priority``): a higher number is more severe.

    ``effect_classes`` lists class names such as ``Substitution`` or ``PrematureStop``; without it every ranked class
    is returned, most severe first.
    """

    ordering = require("varcode.effects.effect_ordering", _VARCODE_EXTRA)
    table = {
        cls.__name__: int(rank) for cls, rank in ordering.transcript_effect_priority_dict.items()
    }
    if effect_classes is None:
        names = sorted(table, key=lambda item: (-table[item], item))
    else:
        if isinstance(effect_classes, (str, bytes)) or not isinstance(effect_classes, Sequence):
            raise ConfigurationError("effect_classes must be a list of class names.")
        names = [_text(item, "effect class") for item in effect_classes]
        unknown = [item for item in names if item not in table]
        if unknown:
            raise ConfigurationError(f"Unknown effect class(es): {', '.join(unknown)}.")
    ranks = [{"effect_class": item, "priority": table[item]} for item in names]
    return _document("effect_priority", ranks=ranks, count=len(ranks))


_SCOPES = ("overall", "variant", "transcript", "gene")


def varcode_top_priority_effects(
    gtf: str,
    genome_fasta: str,
    *,
    scope: str = "variant",
    vcf_paths: Sequence[str] = (),
    maf_paths: Sequence[str] = (),
    csv_paths: Sequence[str] = (),
    json_paths: Sequence[str] = (),
    variants: Sequence[Mapping[str, Any]] = (),
    only_passing: bool = True,
    skip_errors: bool = False,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Highest-priority (most severe) effect of annotated variants.

    ``scope='overall'`` gives the single most severe effect of the whole collection (``top_priority_effect``; give one
    variant for "the top effect of a variant"); ``'variant'`` one effect per variant
    (``top_priority_effect_per_variant``); ``'transcript'`` one per transcript id
    (``top_priority_effect_per_transcript_id``); ``'gene'`` one per gene id (``top_priority_effect_per_gene_id``).
    """

    chosen = _choice(scope, "scope", _SCOPES)
    _genome, collection, info = _load(
        gtf,
        genome_fasta,
        cache_dir,
        vcf_paths,
        maf_paths,
        csv_paths,
        json_paths,
        variants,
        only_passing,
    )
    effects = _effects(collection, skip_errors=skip_errors)
    if chosen == "overall":
        top = _guard(effects.top_priority_effect)
        rows = [{"key": None, **_effect_row(top)}] if top is not None else []
    else:
        method = {
            "variant": "top_priority_effect_per_variant",
            "transcript": "top_priority_effect_per_transcript_id",
            "gene": "top_priority_effect_per_gene_id",
        }[chosen]
        rows = [
            {"key": _variant_label(key), **_effect_row(effect)}
            for key, effect in _guard(getattr(effects, method)).items()
        ]
    return _document(
        "top_priority_effects",
        scope=chosen,
        effects=rows,
        count=len(rows),
        cache_directory=info["cache_directory"],
    )


def _variant_label(key: Any) -> Any:
    if hasattr(key, "contig") and hasattr(key, "ref"):
        return f"{key.contig}:{key.start}:{key.ref}>{key.alt}"
    return None if key is None else str(key)


def varcode_filter_effects(
    gtf: str,
    genome_fasta: str,
    *,
    drop_silent_and_noncoding: bool = True,
    keep_unresolved: bool = True,
    min_priority_class: str | None = None,
    vcf_paths: Sequence[str] = (),
    maf_paths: Sequence[str] = (),
    csv_paths: Sequence[str] = (),
    json_paths: Sequence[str] = (),
    variants: Sequence[Mapping[str, Any]] = (),
    only_passing: bool = True,
    skip_errors: bool = False,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter annotated effects by consequence.

    ``drop_silent_and_noncoding`` removes silent and non-coding effects (``drop_silent_and_noncoding``, the
    ``varcode --only-coding`` option; ``keep_unresolved`` keeps effects whose coding status is unknown).
    ``min_priority_class`` keeps effects at least as severe as the named class, for example
    ``NonsynonymousCodingMutation`` or ``Substitution`` (``filter_by_effect_priority``). Both filters can be combined.
    """

    _flag(drop_silent_and_noncoding, "drop_silent_and_noncoding")
    _flag(keep_unresolved, "keep_unresolved")
    ordering = require("varcode.effects.effect_ordering", _VARCODE_EXTRA)
    minimum = None
    if min_priority_class is not None:
        wanted = _text(min_priority_class, "min_priority_class")
        minimum = next(
            (cls for cls in ordering.transcript_effect_priority_dict if cls.__name__ == wanted),
            None,
        )
        if minimum is None:
            raise ConfigurationError(f"Unknown effect class for min_priority_class: {wanted}.")
    _genome, collection, info = _load(
        gtf,
        genome_fasta,
        cache_dir,
        vcf_paths,
        maf_paths,
        csv_paths,
        json_paths,
        variants,
        only_passing,
    )
    effects = _effects(collection, skip_errors=skip_errors)
    total = len(effects)
    if drop_silent_and_noncoding:
        effects = _guard(effects.drop_silent_and_noncoding, keep_unresolved=keep_unresolved)
    if minimum is not None:
        effects = _guard(effects.filter_by_effect_priority, minimum)
    rows = [_effect_row(effect) for effect in effects]
    return _document(
        "filter_effects",
        input_count=total,
        effects=rows,
        count=len(rows),
        cache_directory=info["cache_directory"],
    )


def varcode_expression_effects(
    gtf: str,
    genome_fasta: str,
    expression: Mapping[str, float],
    *,
    mode: str = "annotate",
    vcf_paths: Sequence[str] = (),
    maf_paths: Sequence[str] = (),
    csv_paths: Sequence[str] = (),
    json_paths: Sequence[str] = (),
    variants: Sequence[Mapping[str, Any]] = (),
    only_passing: bool = True,
    skip_errors: bool = False,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Use transcript expression to rank the effects of annotated variants.

    ``expression`` maps transcript ids to a level such as FPKM or TPM. ``mode='annotate'`` attaches the level to every
    effect that has a transcript (``EffectCollection.effect_expression``; effects without a transcript are left
    out, transcripts not listed count as 0; splice-outcome sets are included, which varcode 10.11.10 itself cannot hash). ``mode='top'`` returns the one effect on the most highly expressed
    transcript, ties broken by effect priority (``top_expression_effect``).
    """

    chosen = _choice(mode, "mode", ("annotate", "top"))
    levels = _expression_map(expression)
    _genome, collection, info = _load(
        gtf,
        genome_fasta,
        cache_dir,
        vcf_paths,
        maf_paths,
        csv_paths,
        json_paths,
        variants,
        only_passing,
    )
    effects = _effects(collection, skip_errors=skip_errors)
    # varcode 10.11.10 hashes effects while building the mapping, which fails for splice-outcome sets; the same
    # arithmetic is applied here effect by effect (level of the effect's transcript, 0 when unlisted).
    ordering = require("varcode.effects.effect_ordering", _VARCODE_EXTRA)
    pairs = [
        (effect, float(levels.get(effect.transcript.id, 0.0)))
        for effect in effects
        if getattr(effect, "transcript", None) is not None
    ]
    if chosen == "annotate":
        rows = [{**_effect_row(effect), "expression": level} for effect, level in pairs]
    else:
        rows = []
        if pairs:
            top, level = _guard(
                max, pairs, key=lambda pair: (pair[1], ordering.multi_gene_effect_sort_key(pair[0]))
            )
            rows = [{**_effect_row(top), "expression": level}]
    return _document(
        "expression_effects",
        mode=chosen,
        effects=rows,
        count=len(rows),
        cache_directory=info["cache_directory"],
    )


# --- splice outcomes with RNA evidence ---------------------------------------------------------


def varcode_reconcile_splice_outcomes_with_rna(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    transcript_id: str,
    observed: Sequence[Mapping[str, Any]],
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reconcile the predicted splice mechanisms of a variant with RNA-observed mechanisms.

    The splice outcomes of the variant on the transcript are enumerated (``enumerate_splice_outcomes``). Each entry of
    ``observed`` names a predicted mechanism class (``mechanism``, for example ``IntronRetention``,
    ``ExonSkipping`` or ``CrypticDonor``) and may give ``read_count`` and ``source``. The set is then reconciled
    with ``SpliceOutcomeSet.with_rna_evidence``: observed mechanisms replace their broad predictions, unobserved
    predictions are listed as excluded, and ``SpliceOutcomeSet.rna_evidence_for`` gives the evidence kept for each
    remaining candidate. Variants that do not disrupt a splice signal are refused.
    """

    varcode = _vc()
    genome, info = _base(gtf, genome_fasta, cache_dir)
    item = _variant(genome, variant)
    transcript = genome.transcript_by_id(_identifier(transcript_id, "transcript_id"))
    if isinstance(observed, (str, bytes)) or not isinstance(observed, Sequence) or not observed:
        raise ConfigurationError("observed must be a non-empty list of mappings with a mechanism.")
    splice = require("varcode.splice_outcomes", _VARCODE_EXTRA)
    rna = require("varcode.rna_evidence", _VARCODE_EXTRA)
    effect = _guard(
        require("varcode.effects", _VARCODE_EXTRA).predict_variant_effect_on_transcript,
        item,
        transcript,
    )

    def provider(contig: str, start: int, end: int) -> str:
        return str(varcode.reference_range(genome, contig, start, end))

    outcome = _guard(splice.enumerate_splice_outcomes, effect, genomic_sequence=provider)
    if not isinstance(outcome, splice.SpliceOutcomeSet):
        raise ConfigurationError("The variant does not disrupt a splice signal on this transcript.")
    chosen = []
    for entry in observed:
        if not isinstance(entry, Mapping) or "mechanism" not in entry:
            raise ConfigurationError("Each observed entry needs a mechanism.")
        mechanism = _text(entry["mechanism"], "mechanism")
        match = next(
            (c for c in outcome.dna_candidates if type(c.effect).__name__ == mechanism), None
        )
        if match is None:
            available = sorted({type(c.effect).__name__ for c in outcome.dna_candidates})
            raise ConfigurationError(
                f"Mechanism {mechanism} is not predicted; predicted: {', '.join(available)}."
            )
        count = entry.get("read_count")
        if count is not None:
            count = _int(count, "read_count", maximum=10**9)
        chosen.append(
            rna.make_rna_outcome(
                match.effect, source=_text(entry.get("source", "rna"), "source"), read_count=count
            )
        )
    reconciled = _guard(outcome.with_rna_evidence, chosen)
    candidates = [
        {
            **_candidate_record(candidate),
            "rna_evidence": [_candidate_record(e) for e in reconciled.rna_evidence_for(candidate)],
        }
        for candidate in reconciled.candidates
    ]
    return _document(
        "reconcile_splice_outcomes_with_rna",
        variant={k: v for k, v in _variant_record(item).items()},
        transcript_id=str(transcript_id),
        predicted=[_candidate_record(c) for c in outcome.dna_candidates],
        candidates=candidates,
        excluded=[_candidate_record(c) for c in reconciled.excluded_candidates],
        added=[_candidate_record(c) for c in reconciled.added_candidates],
        cache_directory=info["cache_directory"],
    )


# --- phasing -----------------------------------------------------------------------------------


_PHASE_METHODS = ("vcf", "rna_reads")


def _rna_options(
    min_mapping_quality: int,
    min_base_quality: int,
    min_alt_reads: int,
    phasing_error_rate: float,
    max_p_value_for_phasing: float,
    require_proper_pair: bool,
) -> dict[str, Any]:
    for name, value in (
        ("phasing_error_rate", phasing_error_rate),
        ("max_p_value_for_phasing", max_p_value_for_phasing),
    ):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 < value < 1:
            raise ConfigurationError(f"{name} must be a number between 0 and 1.")
    return {
        "min_mapping_quality": _int(min_mapping_quality, "min_mapping_quality", maximum=255),
        "min_base_quality": _int(min_base_quality, "min_base_quality", maximum=100),
        "min_alt_reads": _int(min_alt_reads, "min_alt_reads", minimum=1, maximum=10**6),
        "phasing_error_rate": float(phasing_error_rate),
        "max_p_value_for_phasing": float(max_p_value_for_phasing),
        "require_proper_pair": _flag(require_proper_pair, "require_proper_pair"),
    }


def _rna_source(bam_path: object, options: Mapping[str, Any]) -> Any:
    require("pysam", _BAM_EXTRA)
    path = _path(bam_path, "bam_path")
    return _guard(_vc().RNAReadPhasingSource, str(path), **options)


def varcode_query_in_cis(
    gtf: str,
    genome_fasta: str,
    variant_a: Mapping[str, Any],
    variant_b: Mapping[str, Any],
    *,
    method: str = "vcf",
    vcf_path: str | None = None,
    sample: str | None = None,
    bam_path: str | None = None,
    transcript_id: str | None = None,
    min_mapping_quality: int = 20,
    min_base_quality: int = 20,
    min_alt_reads: int = 2,
    phasing_error_rate: float = 0.01,
    max_p_value_for_phasing: float = 0.05,
    require_proper_pair: bool = True,
    only_passing: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Are two variants on the same haplotype (cis), on different ones (trans), or unknown (None)?

    ``method='vcf'`` reads the genotype and phase set of ``sample`` in ``vcf_path`` (``VCFPhaseResolver.in_cis``,
    GT with ``|`` and PS). ``method='rna_reads'`` counts RNA fragments from ``bam_path`` that carry both alleles
    (``RNAReadPhasingSource.in_cis`` through ``MolecularPhaseResolver``); the read filters are the varcode defaults.
    Either resolver is queried through ``varcode.query_in_cis``. ``transcript_id`` limits the question to one transcript.
    """

    chosen = _choice(method, "method", _PHASE_METHODS)
    varcode = _vc()
    genome, info = _base(gtf, genome_fasta, cache_dir)
    first, second = (
        _variant(genome, variant_a, "variant_a"),
        _variant(genome, variant_b, "variant_b"),
    )
    transcript = None
    if transcript_id is not None:
        transcript = genome.transcript_by_id(_identifier(transcript_id, "transcript_id"))
    source = None
    if chosen == "vcf":
        if vcf_path is None or sample is None:
            raise ConfigurationError("method 'vcf' needs vcf_path and sample.")
        collection = _collection(genome, [vcf_path], only_passing=only_passing)
        if sample not in list(collection.samples):
            raise ConfigurationError(f"Sample {sample} is not in the VCF.")
        resolver = _guard(varcode.VCFPhaseResolver, collection, sample)
    else:
        if bam_path is None:
            raise ConfigurationError("method 'rna_reads' needs bam_path.")
        source = _rna_source(
            bam_path,
            _rna_options(
                min_mapping_quality,
                min_base_quality,
                min_alt_reads,
                phasing_error_rate,
                max_p_value_for_phasing,
                require_proper_pair,
            ),
        )
        resolver = varcode.MolecularPhaseResolver(source)
    try:
        answer = _guard(varcode.query_in_cis, resolver, first, second, transcript)
    finally:
        if source is not None:
            source.close()
    return _document(
        "query_in_cis",
        method=chosen,
        phase_source=str(resolver.phase_source),
        variant_a=_variant_record(first),
        variant_b=_variant_record(second),
        in_cis=answer,
        relation="cis" if answer is True else "trans" if answer is False else "unknown",
        cache_directory=info["cache_directory"],
    )


def varcode_phased_partners(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    *,
    method: str = "vcf",
    vcf_path: str | None = None,
    sample: str | None = None,
    bam_path: str | None = None,
    candidates: Sequence[Mapping[str, Any]] = (),
    transcript_id: str | None = None,
    min_mapping_quality: int = 20,
    min_base_quality: int = 20,
    min_alt_reads: int = 2,
    phasing_error_rate: float = 0.01,
    max_p_value_for_phasing: float = 0.05,
    require_proper_pair: bool = True,
    only_passing: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Variants observed in cis with a variant.

    ``method='vcf'`` lists the variants in the same VCF phase set and haplotype of ``sample``
    (``VCFPhaseResolver.phased_partners``). ``method='rna_reads'`` lists the ``candidates`` that are co-observed with
    the variant on the same RNA read or fragment of ``bam_path`` (``MolecularPhaseResolver.phased_partners`` and
    ``RNAReadPhasingSource.partners_in_cis``); an empty list means no direct evidence.
    """

    chosen = _choice(method, "method", _PHASE_METHODS)
    varcode = _vc()
    genome, info = _base(gtf, genome_fasta, cache_dir)
    item = _variant(genome, variant)
    transcript = (
        genome.transcript_by_id(_identifier(transcript_id, "transcript_id"))
        if transcript_id
        else None
    )
    source = None
    if chosen == "vcf":
        if vcf_path is None or sample is None:
            raise ConfigurationError("method 'vcf' needs vcf_path and sample.")
        collection = _collection(genome, [vcf_path], only_passing=only_passing)
        if sample not in list(collection.samples):
            raise ConfigurationError(f"Sample {sample} is not in the VCF.")
        resolver = _guard(varcode.VCFPhaseResolver, collection, sample)
        known = [entry for entry in collection if _row(entry) == _row(item)]
        item = known[0] if known else item
    else:
        if bam_path is None:
            raise ConfigurationError("method 'rna_reads' needs bam_path.")
        pool = [_variant(genome, entry, "candidate") for entry in candidates]
        source = _rna_source(
            bam_path,
            _rna_options(
                min_mapping_quality,
                min_base_quality,
                min_alt_reads,
                phasing_error_rate,
                max_p_value_for_phasing,
                require_proper_pair,
            ),
        )
        source.register_variants([item, *pool])
        resolver = varcode.MolecularPhaseResolver(source)
    try:
        partners = [_row(entry) for entry in _guard(resolver.phased_partners, item, transcript)]
    finally:
        if source is not None:
            source.close()
    return _document(
        "phased_partners",
        method=chosen,
        phase_source=str(resolver.phase_source),
        variant=_variant_record(item),
        partners=partners,
        count=len(partners),
        cache_directory=info["cache_directory"],
    )


def varcode_rna_read_support(
    gtf: str,
    genome_fasta: str,
    bam_path: str,
    variants: Sequence[Mapping[str, Any]],
    *,
    min_mapping_quality: int = 20,
    min_base_quality: int = 20,
    min_alt_reads: int = 2,
    phasing_error_rate: float = 0.01,
    max_p_value_for_phasing: float = 0.05,
    require_proper_pair: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the quality-filtered RNA reads or fragments carrying each variant's alt allele.

    For each variant returns ``alt_reads`` (``RNAReadPhasingSource.supports_variant``; None when the contig is not in the
    BAM) and ``has_evidence`` (``has_evidence``: at least ``min_alt_reads`` alt fragments). The BAM must be
    coordinate-sorted and indexed.
    """

    if isinstance(variants, (str, bytes)) or not isinstance(variants, Sequence) or not variants:
        raise ConfigurationError("variants must be a non-empty list of records.")
    genome, info = _base(gtf, genome_fasta, cache_dir)
    items = [_variant(genome, entry) for entry in variants]
    options = _rna_options(
        min_mapping_quality,
        min_base_quality,
        min_alt_reads,
        phasing_error_rate,
        max_p_value_for_phasing,
        require_proper_pair,
    )
    source = _rna_source(bam_path, options)
    try:
        rows = [
            {
                **_variant_record(item),
                "alt_reads": source.supports_variant(item),
                "has_evidence": bool(source.has_evidence(item)),
            }
            for item in items
        ]
    finally:
        source.close()
    return _document(
        "rna_read_support",
        min_alt_reads=options["min_alt_reads"],
        variants=rows,
        cache_directory=info["cache_directory"],
    )


def varcode_test_rna_haplotype(
    gtf: str,
    genome_fasta: str,
    bam_path: str,
    variants: Sequence[Mapping[str, Any]],
    *,
    flanking_bases: int = 5,
    min_mapping_quality: int = 20,
    min_base_quality: int = 20,
    min_alt_reads: int = 2,
    phasing_error_rate: float = 0.01,
    max_p_value_for_phasing: float = 0.05,
    require_proper_pair: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether RNA reads carry a known local combination of alleles on one alignment (``RNAReadPhasingSource.register_haplotype``).

    ``variants`` are non-overlapping substitutions or deletions on one contig that form one hypothesis. A read supports it
    only when its whole sequence, including ``flanking_bases`` unchanged reference bases on each side, equals the
    hypothesis. Returns per-variant support after registration and, for each pair, the cis/trans call.
    """

    if (
        isinstance(variants, (str, bytes))
        or not isinstance(variants, Sequence)
        or len(variants) < 1
    ):
        raise ConfigurationError("variants must be a non-empty list of records.")
    count = _int(flanking_bases, "flanking_bases", minimum=1, maximum=100)
    genome, info = _base(gtf, genome_fasta, cache_dir)
    items = [_variant(genome, entry) for entry in variants]
    options = _rna_options(
        min_mapping_quality,
        min_base_quality,
        min_alt_reads,
        phasing_error_rate,
        max_p_value_for_phasing,
        require_proper_pair,
    )
    source = _rna_source(bam_path, options)
    try:
        _guard(source.register_haplotype, items, flanking_bases=count)
        support = [
            {**_variant_record(item), "alt_reads": source.supports_variant(item)} for item in items
        ]
        pairs = [
            {
                "first": _variant_record(a),
                "second": _variant_record(b),
                "in_cis": source.in_cis(a, b),
            }
            for index, a in enumerate(items)
            for b in items[index + 1 :]
        ]
    finally:
        source.close()
    return _document(
        "test_rna_haplotype",
        flanking_bases=count,
        support=support,
        pairs=pairs,
        cache_directory=info["cache_directory"],
    )


# --- germline context --------------------------------------------------------------------------

_COMPLETENESS = ("complete", "sparse", "hotspots_only", "empty")


def _context(
    genome: Any,
    vcf_path: object,
    sample: object,
    completeness: object,
    only_passing: bool,
) -> Any:
    germline = require("varcode.germline", _VARCODE_EXTRA)
    level = germline.Completeness(_choice(completeness, "completeness", _COMPLETENESS))
    path = _path(vcf_path, "vcf_path")
    options = {
        "genome": genome,
        "only_passing": _flag(only_passing, "only_passing"),
        "convert_ucsc_contig_names": False,
    }
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if sample is None:
            context = _guard(
                germline.GermlineContext.from_germline_vcf, str(path), completeness=level, **options
            )
        else:
            context = _guard(
                germline.GermlineContext.from_multi_sample_vcf,
                str(path),
                _text(sample, "sample"),
                completeness=level,
                **options,
            )
    return context, [str(item.message) for item in caught]


def varcode_load_germline_context(
    gtf: str,
    genome_fasta: str,
    vcf_path: str,
    *,
    sample: str | None = None,
    completeness: str = "complete",
    only_passing: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load a germline VCF into a varcode ``GermlineContext``.

    Without ``sample`` the whole VCF is the germline call set (``GermlineContext.from_germline_vcf``). With ``sample`` that
    sample's non-reference calls are extracted from a multi-sample VCF (``from_multi_sample_vcf``). ``completeness`` states
    what a missing call means: ``complete`` (absent means reference), ``sparse``, ``hotspots_only`` or ``empty``.
    """

    genome, info = _base(gtf, genome_fasta, cache_dir)
    context, notes = _context(genome, vcf_path, sample, completeness, only_passing)
    rows = _rows(context.variants)
    return _document(
        "load_germline_context",
        sample=sample,
        completeness=context.completeness.value,
        reference_name=context.reference_name,
        variants=rows,
        count=len(rows),
        warnings=notes,
        cache_directory=info["cache_directory"],
    )


def varcode_germline_variants(
    gtf: str,
    genome_fasta: str,
    vcf_path: str,
    *,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
    homozygous_only: bool = False,
    sample: str | None = None,
    completeness: str = "complete",
    only_passing: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Germline variants of a context inside a window, or the homozygous ones.

    With ``contig``, ``start`` and ``end`` (1-based, inclusive) the variants overlapping that window are returned
    (``GermlineContext.variants_in_window``). With ``homozygous_only`` the result is limited to the variants the sample carries
    on both haplotypes (``GermlineContext.homozygous``); this needs a single-sample VCF, because varcode 10.11.10
    keeps no genotypes when it extracts one sample from a multi-sample VCF (``sample`` is refused with
    ``homozygous_only``). Without a window the whole context is considered.
    """

    _flag(homozygous_only, "homozygous_only")
    if homozygous_only and sample is not None:
        raise ConfigurationError(
            "homozygous_only needs a single-sample VCF: varcode keeps no genotypes when extracting a sample."
        )
    window = (contig, start, end)
    if any(item is not None for item in window) and any(item is None for item in window):
        raise ConfigurationError("contig, start and end must be given together.")
    genome, info = _base(gtf, genome_fasta, cache_dir)
    context, notes = _context(genome, vcf_path, sample, completeness, only_passing)
    if contig is not None:
        found = context.variants_in_window(
            _identifier(contig, "contig"),
            _int(start, "start", minimum=1),
            _int(end, "end", minimum=1),
        )
    else:
        found = tuple(context.variants)
    if homozygous_only:
        found = context.homozygous(found)
    rows = [_row(item) for item in found]
    return _document(
        "germline_variants",
        window={"contig": contig, "start": start, "end": end},
        homozygous_only=homozygous_only,
        variants=rows,
        count=len(rows),
        warnings=notes,
        cache_directory=info["cache_directory"],
    )


def varcode_validate_germline_context(
    gtf: str,
    genome_fasta: str,
    germline_vcf: str,
    somatic: Mapping[str, Any],
    *,
    sample: str | None = None,
    completeness: str = "complete",
    validate_reference: bool = True,
    only_passing: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check a germline context against somatic variants before annotation (``GermlineContext.validate_against``).

    ``somatic`` is a source mapping (keys ``vcf``, ``maf``, ``csv``, ``json``, ``variants``). The check fails on a
    reference-build mismatch (unless ``validate_reference`` is False) and warns about empty or suspicious call sets.
    Returns ``valid`` with the error text when it fails, and the warnings varcode raised.
    """

    _flag(validate_reference, "validate_reference")
    genome, info = _base(gtf, genome_fasta, cache_dir)
    context, notes = _context(genome, germline_vcf, sample, completeness, only_passing)
    collection = _sources(genome, somatic, "somatic", only_passing=only_passing)
    error = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            context.validate_against(collection, validate_reference=validate_reference)
        except Exception as exc:  # varcode raises GenomeBuildMismatchError and similar
            error = f"{type(exc).__name__}: {exc}"
    return _document(
        "validate_germline_context",
        valid=error is None,
        error=error,
        warnings=notes + [str(item.message) for item in caught],
        germline_count=len(context.variants),
        somatic_count=len(collection),
        germline_reference=context.reference_name,
        cache_directory=info["cache_directory"],
    )


# --- structural variants -----------------------------------------------------------------------


def varcode_sv_junction_inserted_sequence(
    gtf: str,
    genome_fasta: str,
    vcf_paths: Sequence[str],
    *,
    only_passing: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bases inserted at each junction of the structural variants in VCF files (``StructuralVariant.junction_inserted_sequence``).

    The sequence is oriented from the retained first breakend towards the retained second one. An empty string means no
    insertion; None means the allele is conflicting or unreadable. Small variants in the VCF are skipped.
    """

    _genome, collection, info = _load(
        gtf, genome_fasta, cache_dir, vcf_paths, (), (), (), (), only_passing
    )
    rows: list[dict[str, Any]] = []
    for item in collection:
        if not getattr(item, "is_structural", False):
            continue
        for junction in item.junctions:
            rows.append(
                {
                    "contig": str(item.contig),
                    "start": int(item.start),
                    "sv_type": str(item.sv_type),
                    "alt": str(item.alt),
                    "junction": [
                        {
                            "contig": str(end.contig),
                            "position": int(end.position),
                            "keeps": str(end.keeps),
                        }
                        for end in junction
                    ],
                    "inserted_sequence": _guard(item.junction_inserted_sequence, junction),
                }
            )
    return _document(
        "sv_junction_inserted_sequence",
        junctions=rows,
        count=len(rows),
        cache_directory=info["cache_directory"],
    )


__all__ = [
    "varcode_effect_priority",
    "varcode_effect_table",
    "varcode_expression_effects",
    "varcode_filter_by_expression",
    "varcode_filter_effects",
    "varcode_genotype_alt_summary",
    "varcode_germline_variants",
    "varcode_load_germline_context",
    "varcode_phased_partners",
    "varcode_query_in_cis",
    "varcode_read_effect_csv",
    "varcode_read_variant_csv",
    "varcode_reconcile_splice_outcomes_with_rna",
    "varcode_rna_read_support",
    "varcode_sample_genotype",
    "varcode_sv_junction_inserted_sequence",
    "varcode_test_rna_haplotype",
    "varcode_top_priority_effects",
    "varcode_validate_germline_context",
    "varcode_variant_set_operation",
    "varcode_variant_table",
    "varcode_variants_in_sample",
]

_register_public_verbose(globals(), tuple(__all__))
