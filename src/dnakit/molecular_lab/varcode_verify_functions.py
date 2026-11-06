"""varcode command-line operations that DNAKit did not wrap (Fig. 2 coverage gap, batch 0 verify).

Four public functions, each taking explicit local inputs:

- ``varcode_load_variant_sources`` loads variants from VCF files, MAF files, JSON variant collections and
  explicit variant records, and unions them, as ``varcode.cli.variant_args.variant_collection_from_args``
  does. A local reference (``gtf`` and ``genome_fasta``) or an installed pyensembl reference name
  (``genome``) is needed for explicit records and for VCF files without a reference in their header.
  The ``--download-reference-genome-data`` step of the varcode command line is a network download and is
  not part of this function; ``varcode_install_reference_data`` covers it.
- ``varcode_effects_table`` loads variants the same way, predicts the coding effects of each variant with
  ``varcode_predict_variant_effects`` and writes one CSV row per effect (the ``varcode effects`` command).
  The ``--only-coding`` and ``--one-per-variant`` filters of that command are not included.
- ``varcode_check_sample_identity_files`` runs the ``varcode check-samples`` command: it loads germline,
  somatic and mixed VCF files (and an optional manifest), compares the samples, and writes the JSON and TSV
  reports.
- ``varcode_compare_sv_calls_directory`` runs the varcode structural-variant comparison command: it reads a
  CSV of SV calls and writes ``groups.csv``, ``members.csv`` and ``comparison.json`` into a new directory.

Outputs are never overwritten: an existing output file or directory is refused. varcode is required with
the ``dnakit[bio-varcode]`` extra. A local reference also needs pyensembl (``dnakit[bio-pyensembl]``).
"""

from __future__ import annotations

import csv
import dataclasses
import hashlib
import json
import os
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab.varcode_gap_functions import (
    _VARCODE_EXTRA,
    _genome,
    varcode_predict_variant_effects,
)

_MAX_VARIANTS = 100000
_SAMPLE_KINDS = ("germline", "somatic", "mixed")
_MANIFEST_COLUMNS = {
    "path",
    "kind",
    "assembly",
    "sample",
    "role",
    "label",
    "donor_id",
    "tumor_id",
    "normal_sample",
}
_SAMPLE_TSV_FIELDS = [
    "sample_a",
    "sample_b",
    "germline_status",
    "shared_snps",
    "genotype_concordance",
    "ibs0_rate",
    "somatic_status",
    "shared_somatic",
    "somatic_overlap",
    "flags",
]
_EFFECT_FIELDS = [
    "contig",
    "start",
    "ref",
    "alt",
    "class",
    "transcript_id",
    "gene_name",
    "description",
]


def _varcode_package() -> str:
    version = require("varcode.version", _VARCODE_EXTRA).__version__
    return f"varcode {version}"


def _file(value: object, name: str) -> Path:
    if not isinstance(value, (str, os.PathLike)) or not str(value):
        raise ConfigurationError(f"{name} must be a file path.")
    target = Path(value).expanduser().resolve()
    if not target.is_file():
        raise ConfigurationError(f"{name} does not exist: {target}")
    return target


def _file_list(values: object, name: str) -> list[Path]:
    if values is None:
        return []
    if isinstance(values, (str, bytes, os.PathLike)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of file paths.")
    return [_file(item, f"{name} entry") for item in values]


def _variant_records(values: object) -> list[Mapping[str, Any]]:
    if values is None:
        return []
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(
            "variants must be a list of mappings with contig, start, ref and alt."
        )
    records: list[Mapping[str, Any]] = []
    for item in values:
        if not isinstance(item, Mapping) or not {"contig", "start", "ref", "alt"} <= set(item):
            raise ConfigurationError("Each variant needs contig, start, ref and alt fields.")
        if (
            isinstance(item["start"], bool)
            or not isinstance(item["start"], int)
            or item["start"] < 1
        ):
            raise ConfigurationError("variant start must be a 1-based integer position.")
        records.append(item)
    return records


def _new_output(path: object, name: str) -> Path:
    if not isinstance(path, (str, os.PathLike)) or not str(path):
        raise ConfigurationError(f"{name} must be a path.")
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _reference(
    gtf: object, genome_fasta: object, cache_dir: object, genome: object
) -> tuple[Any, str | None]:
    """Return the genome object for varcode and a reference label for the result."""
    if gtf is not None or genome_fasta is not None:
        if gtf is None or genome_fasta is None:
            raise ConfigurationError("gtf and genome_fasta must be given together.")
        built, info = _genome(gtf, genome_fasta, cache_dir)
        return built, f"local:{info['digest']}"
    if genome is None:
        return None, None
    if not isinstance(genome, str) or not genome.strip() or any(c.isspace() for c in genome):
        raise ConfigurationError("genome must be a reference name such as 'GRCh38'.")
    return genome, genome


def _load_collections(
    vcf_paths: list[Path],
    maf_paths: list[Path],
    json_paths: list[Path],
    records: list[Mapping[str, Any]],
    genome: Any,
    only_passing: bool,
) -> Any:
    varcode = require("varcode", _VARCODE_EXTRA)
    collections: list[Any] = []
    try:
        for path in vcf_paths:
            collections.append(
                varcode.load_vcf(
                    str(path),
                    genome=genome,
                    parse_structural_variants=True,
                    only_passing=only_passing,
                    warn_on_filtered=True,
                )
            )
        for path in maf_paths:
            collections.append(require("varcode.maf", _VARCODE_EXTRA).load_maf(str(path)))
        if records:
            if genome is None:
                raise ConfigurationError(
                    "Explicit variants need gtf and genome_fasta, or a genome name."
                )
            collections.append(
                varcode.VariantCollection(
                    [
                        varcode.Variant(
                            str(item["contig"]),
                            start=int(item["start"]),
                            ref=str(item["ref"]),
                            alt=str(item["alt"]),
                            genome=genome,
                            # UCSC names are converted only for named Ensembl references; a local
                            # reference keeps the contig names of its GTF.
                            convert_ucsc_contig_names=isinstance(genome, str),
                        )
                        for item in records
                    ]
                )
            )
        for path in json_paths:
            collections.append(
                varcode.VariantCollection.from_json(path.read_text(encoding="utf-8"))
            )
    except ValueError as exc:
        raise ConfigurationError(f"varcode could not load the variants: {exc}") from exc
    if not collections:
        raise ConfigurationError(
            "No variants were given: pass vcf_paths, maf_paths, json_paths or variants."
        )
    return varcode.VariantCollection.union(*collections)


def varcode_load_variant_sources(
    vcf_paths: Sequence[str] = (),
    maf_paths: Sequence[str] = (),
    json_paths: Sequence[str] = (),
    variants: Sequence[Mapping[str, Any]] = (),
    *,
    gtf: str | None = None,
    genome_fasta: str | None = None,
    cache_dir: str | None = None,
    genome: str | None = None,
    only_passing: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load variants from VCF, MAF and JSON files and explicit records, and union them (varcode CLI loader).

    ``variants`` are mappings with ``contig``, 1-based ``start``, ``ref`` and ``alt``. A local reference is built
    from ``gtf`` and ``genome_fasta``; ``genome`` names an installed pyensembl reference instead. Returns the
    unioned variant rows with their count and the number of records taken from each source.
    """

    vcf = _file_list(vcf_paths, "vcf_paths")
    maf = _file_list(maf_paths, "maf_paths")
    json_files = _file_list(json_paths, "json_paths")
    records = _variant_records(variants)
    if not isinstance(only_passing, bool):
        raise ConfigurationError("only_passing must be True or False.")
    reference, label = _reference(gtf, genome_fasta, cache_dir, genome)
    collection = _load_collections(vcf, maf, json_files, records, reference, only_passing)
    rows = [
        {
            "contig": str(item.contig),
            "start": int(item.start),
            "ref": str(item.ref),
            "alt": str(item.alt),
        }
        for item in collection
    ]
    if len(rows) > _MAX_VARIANTS:
        raise ConfigurationError(f"The variants exceed the limit of {_MAX_VARIANTS} records.")
    return {
        "operation": "load_variant_sources",
        "status": "completed",
        "backend": "varcode",
        "package": _varcode_package(),
        "coordinates": "1-based inclusive",
        "reference": label,
        "count": len(rows),
        "sources": {
            "vcf": len(vcf),
            "maf": len(maf),
            "json": len(json_files),
            "variant": len(records),
        },
        "variants": rows,
    }


def varcode_effects_table(
    output_csv: str,
    *,
    gtf: str,
    genome_fasta: str,
    vcf_paths: Sequence[str] = (),
    maf_paths: Sequence[str] = (),
    json_paths: Sequence[str] = (),
    variants: Sequence[Mapping[str, Any]] = (),
    cache_dir: str | None = None,
    only_passing: bool = True,
    skip_errors: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict the coding effects of loaded variants against a local reference and write them to a CSV file.

    Variants are loaded as in ``varcode_load_variant_sources``. Each variant is passed to
    ``varcode_predict_variant_effects`` with the same ``gtf`` and ``genome_fasta``. With ``skip_errors`` the
    annotation failures are kept as rows of class ``Failure`` instead of stopping the run. The CSV has one row
    per effect; a variant without any effect gets no row.
    """

    target = _new_output(output_csv, "output_csv")
    loaded = varcode_load_variant_sources(
        vcf_paths,
        maf_paths,
        json_paths,
        variants,
        gtf=gtf,
        genome_fasta=genome_fasta,
        cache_dir=cache_dir,
        only_passing=only_passing,
        verbose=False,
    )
    rows: list[dict[str, Any]] = []
    classes: Counter[str] = Counter()
    for variant in loaded["variants"]:
        result = varcode_predict_variant_effects(
            gtf,
            genome_fasta,
            variant,
            raise_on_error=not skip_errors,
            cache_dir=cache_dir,
            verbose=False,
        )
        for effect in result["effects"]:
            classes[effect["class"]] += 1
            rows.append(
                {
                    **variant,
                    "class": effect["class"],
                    "transcript_id": effect.get("transcript_id", ""),
                    "gene_name": effect.get("gene_name", ""),
                    "description": effect.get("description", ""),
                }
            )
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_EFFECT_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in _EFFECT_FIELDS})
    return {
        "operation": "effects_table",
        "status": "completed",
        "backend": "varcode",
        "package": _varcode_package(),
        "coordinates": "1-based inclusive",
        "reference": loaded["reference"],
        "variant_count": loaded["count"],
        "effect_count": len(rows),
        "class_counts": dict(sorted(classes.items())),
        "columns": _EFFECT_FIELDS,
        "output_csv": str(target),
    }


def _config(overrides: Mapping[str, Any] | None, sample_identity: Any) -> Any:
    if overrides is None:
        return sample_identity.SampleCheckConfig()
    if not isinstance(overrides, Mapping):
        raise ConfigurationError("config must be a mapping of SampleCheckConfig fields.")
    allowed = {field.name for field in dataclasses.fields(sample_identity.SampleCheckConfig)}
    unknown = sorted(set(overrides) - allowed)
    if unknown:
        raise ConfigurationError(f"Unknown config field(s): {', '.join(unknown)}.")
    try:
        return sample_identity.SampleCheckConfig(**dict(overrides))
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"Invalid config: {exc}") from exc


def _manifest_samples(
    path: Path, assembly: str | None, settings: Any, sample_identity: Any
) -> list[Any]:
    groups: dict[tuple[Path, str, str | None], list[Any]] = {}
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source, delimiter="\t" if path.suffix.lower() == ".tsv" else ",")
        if not reader.fieldnames or not {"path", "sample"} <= set(reader.fieldnames):
            raise ConfigurationError("Manifest requires path and sample columns.")
        if (
            len(set(reader.fieldnames)) != len(reader.fieldnames)
            or set(reader.fieldnames) - _MANIFEST_COLUMNS
        ):
            raise ConfigurationError("Duplicate or unknown manifest columns.")
        for number, row in enumerate(reader, 2):
            if (
                None in row
                or any(value is None for value in row.values())
                or not row["path"]
                or not row["sample"]
            ):
                raise ConfigurationError(f"Malformed manifest row {number}.")
            source_path = (path.parent / Path(row["path"]).expanduser()).resolve()
            key = (source_path, row.get("kind") or "germline", row.get("assembly") or assembly)
            spec = sample_identity.SampleSpec(
                **{k: v for k, v in row.items() if k not in {"path", "kind", "assembly"} and v}
            )
            groups.setdefault(key, []).append(spec)
    samples: list[Any] = []
    for (source_path, kind, build), specs in groups.items():
        samples.extend(
            sample_identity.load_vcf_samples(
                source_path, kind=kind, assembly=build, samples=specs, config=settings
            )
        )
    return samples


def _write_sample_tsv(path: Path, report: Mapping[str, Any]) -> None:
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=_SAMPLE_TSV_FIELDS, delimiter="\t")
        writer.writeheader()
        for pair in report["pairs"]:
            germline, somatic = pair["germline"], pair["somatic"]
            writer.writerow(
                dict(
                    sample_a=pair["sample_a"],
                    sample_b=pair["sample_b"],
                    germline_status=germline["status"],
                    shared_snps=germline.get("shared_snps"),
                    genotype_concordance=germline.get("genotype_concordance"),
                    ibs0_rate=germline.get("ibs0_rate"),
                    somatic_status=somatic["status"],
                    shared_somatic=somatic.get("shared_variants"),
                    somatic_overlap=somatic.get("overlap_coefficient"),
                    flags=";".join(pair["flags"]),
                )
            )


def varcode_check_sample_identity_files(
    *,
    germline_vcfs: Sequence[str] = (),
    somatic_vcfs: Sequence[str] = (),
    mixed_vcfs: Sequence[str] = (),
    manifest_path: str | None = None,
    assembly: str | None = None,
    config: Mapping[str, Any] | None = None,
    json_output: str | None = None,
    tsv_output: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check whether the samples in several local VCF files are the same donors (varcode check-samples).

    VCF files are grouped by role: ``germline_vcfs``, ``somatic_vcfs`` and ``mixed_vcfs``. An optional
    ``manifest_path`` (CSV or TSV with path and sample columns) selects samples and roles. ``config`` overrides
    ``SampleCheckConfig`` fields. The full report is written to ``json_output`` and a pair table to
    ``tsv_output``; the report is also returned. Input files are never overwritten.
    """

    sample_identity = require("varcode.sample_identity", _VARCODE_EXTRA)
    settings = _config(config, sample_identity)
    groups = {
        "germline": _file_list(germline_vcfs, "germline_vcfs"),
        "somatic": _file_list(somatic_vcfs, "somatic_vcfs"),
        "mixed": _file_list(mixed_vcfs, "mixed_vcfs"),
    }
    manifest = None if manifest_path is None else _file(manifest_path, "manifest_path")
    outputs = [Path(p).expanduser().resolve() for p in (json_output, tsv_output) if p is not None]
    given_inputs = {path for paths in groups.values() for path in paths} | (
        {manifest} if manifest else set()
    )
    if given_inputs & set(outputs) or len(set(outputs)) != len(outputs):
        raise ConfigurationError("Output paths must differ from the inputs and from one another.")
    json_target = _new_output(json_output, "json_output") if json_output is not None else None
    tsv_target = _new_output(tsv_output, "tsv_output") if tsv_output is not None else None
    samples: list[Any] = []
    try:
        for kind in _SAMPLE_KINDS:
            for path in groups[kind]:
                samples.extend(
                    sample_identity.load_vcf_samples(
                        path, kind=kind, assembly=assembly, config=settings
                    )
                )
        if manifest is not None:
            samples.extend(_manifest_samples(manifest, assembly, settings, sample_identity))
        if not samples:
            raise ConfigurationError(
                "No samples: pass germline, somatic or mixed VCF files, or a manifest."
            )
        report = sample_identity.check_sample_identity(samples, config=settings)
    except ValueError as exc:
        raise ConfigurationError(f"Sample identity check failed: {exc}") from exc
    if {Path(sample.path).resolve() for sample in samples} & set(outputs):
        raise ConfigurationError("Output paths must differ from the inputs and from one another.")
    if json_target is not None:
        json_target.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
        )
    if tsv_target is not None:
        _write_sample_tsv(tsv_target, report)
    return {
        "operation": "check_sample_identity_files",
        "status": "completed",
        "backend": "varcode",
        "package": _varcode_package(),
        "sample_count": len(samples),
        "pair_count": len(report["pairs"]),
        "identity_conflicts": report["identity_conflicts"],
        "json_output": None if json_target is None else str(json_target),
        "tsv_output": None if tsv_target is None else str(tsv_target),
        "report": report,
    }


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    rows = list(rows)
    with path.open("w", newline="", encoding="utf-8") as out:
        fields = sorted({key for row in rows for key in row}) or ["group_id"]
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {
                key: json.dumps(value, sort_keys=True)
                if isinstance(value, (dict, tuple, list))
                else value
                for key, value in row.items()
            }
            for row in rows
        )


def varcode_compare_sv_calls_directory(
    calls_csv: str,
    output_dir: str,
    *,
    max_distance: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group structural-variant breakpoints from a CSV of calls and write the groups to a new directory.

    ``calls_csv`` holds the documented VCF allele and provenance columns (see varcode's ``compare_sv_calls``).
    Breakpoints closer than ``max_distance`` bases are grouped. The output directory receives ``groups.csv``,
    ``members.csv`` and ``comparison.json`` (with input provenance); it must not exist yet.
    """

    sv = require("varcode.sv_comparison", _VARCODE_EXTRA)
    source = _file(calls_csv, "calls_csv")
    target = _new_output(output_dir, "output_dir")
    if isinstance(max_distance, bool) or not isinstance(max_distance, int) or max_distance < 0:
        raise ConfigurationError("max_distance must be a non-negative integer number of bases.")
    with source.open(newline="", encoding="utf-8") as handle:
        try:
            result = sv.compare_sv_calls(csv.DictReader(handle), max_distance)
        except ValueError as exc:
            raise ConfigurationError(f"SV comparison failed: {exc}") from exc
    target.mkdir()
    _write_rows(target / "groups.csv", result["groups"])
    _write_rows(target / "members.csv", result["members"])
    version = require("varcode.version", _VARCODE_EXTRA).__version__
    result["provenance"] = {
        "varcode_version": version,
        "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "input": str(source),
        "call_count": len(result["members"]),
        "group_count": len(result["groups"]),
    }
    (target / "comparison.json").write_text(
        json.dumps(result, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {
        "operation": "compare_sv_calls_directory",
        "status": "completed",
        "backend": "varcode",
        "package": _varcode_package(),
        "call_count": len(result["members"]),
        "group_count": len(result["groups"]),
        "output_dir": str(target),
        "files": ["groups.csv", "members.csv", "comparison.json"],
    }


__all__ = [
    "varcode_check_sample_identity_files",
    "varcode_compare_sv_calls_directory",
    "varcode_effects_table",
    "varcode_load_variant_sources",
]

_register_public_verbose(globals())
