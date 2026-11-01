"""Picard VCF, interval-list, FASTA and GFF tools as typed DNAKit facades (integration pack P08,
part 3).

The user installs Picard (``conda install -c bioconda picard``, or ``picard.jar`` plus a Java
runtime). Each function runs one Picard tool with an explicit, typed ``KEY=value`` argument list (no
shell). Picard is found from ``picard_path`` (a ``picard`` launcher or a ``picard.jar``), then
``picard`` on PATH, then the ``PICARD_JAR`` environment variable; a missing program raises
``BackendUnavailableError`` with an install hint. Outputs are new explicit paths: an existing output
is never overwritten. Results are JSON-serialisable dictionaries with ``status`` (``completed`` only
when the exit code is 0 and every expected output exists), ``backend``, ``tool``, ``version``,
``command``, ``return_code``, ``outputs`` and ``log_tail``.

Optional arguments left as ``None`` are not passed, so the defaults are exactly Picard's own. The
format of a VCF/BCF output follows its extension (``.vcf``, ``.vcf.gz``, ``.bcf``); Picard also
writes an index next to a compressed VCF. Not exposed on purpose: Picard's ``FilterVcf
JAVASCRIPT_FILE`` (it would run caller-supplied scripts). Tested with Picard 3.5.0.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _picard_cli as _p
from dnakit.exceptions import ConfigurationError

_INTERVAL_ACTIONS = ("CONCAT", "UNION", "INTERSECT", "SUBTRACT", "SYMDIFF", "OVERLAPS")
_N_OUTPUT_TYPES = ("N", "ACGT", "BOTH")
_COUNT_VALUES = ("BASES", "INTERVALS")


def _index_of(path: Path) -> str | None:
    for suffix in (".tbi", ".idx", ".csi"):
        candidate = Path(f"{path}{suffix}")
        if candidate.is_file():
            return str(candidate)
    return None


def _vcf_records(path: Path) -> int | None:
    """Number of data lines of a plain VCF; ``None`` for compressed or binary outputs."""
    if path.suffix in {".gz", ".bgz", ".bcf"}:
        return None
    return _p.count_lines(path, skip_prefix="#")


# ---------------------------------------------------------------------------
# VCF tools
# ---------------------------------------------------------------------------


def picard_filter_vcf(
    input_path: str,
    output_path: str,
    *,
    min_dp: int | None = None,
    min_gq: int | None = None,
    min_ab: float | None = None,
    min_qd: float | None = None,
    max_fs: float | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply hard filters to the genotypes and sites of a VCF/BCF (``FilterVcf``).

    ``min_dp`` and ``min_gq`` filter individual genotypes (Picard default 0 = off), ``min_ab`` is
    the minimum allele balance of heterozygotes, ``min_qd`` the minimum QD and ``max_fs`` the
    maximum Fisher strand value that a site may have before it is filtered. The JavaScript filter is
    not offered.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={source}", f"OUTPUT={target}"]
    args += _p.a_int("MIN_DP", min_dp, "min_dp", minimum=0)
    args += _p.a_int("MIN_GQ", min_gq, "min_gq", minimum=0)
    args += _p.a_num("MIN_AB", min_ab, "min_ab", minimum=0, maximum=1)
    args += _p.a_num("MIN_QD", min_qd, "min_qd", minimum=0)
    args += _p.a_num("MAX_FS", max_fs, "max_fs", minimum=0)
    result = _p.run(
        "FilterVcf",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    result["index"] = _index_of(target)
    result["records"] = _vcf_records(target) if target.is_file() else None
    return result


def picard_liftover_vcf(
    input_path: str,
    output_path: str,
    reject_path: str,
    chain: str,
    reference_sequence: str,
    *,
    liftover_min_match: float | None = None,
    recover_swapped_ref_alt: bool | None = None,
    allow_missing_fields_in_header: bool | None = None,
    tags_to_drop: Sequence[str] | None = None,
    tags_to_reverse: Sequence[str] | None = None,
    write_original_position: bool | None = None,
    write_original_alleles: bool | None = None,
    disable_sort: bool | None = None,
    warn_on_missing_contig: bool | None = None,
    log_failed_intervals: bool | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Lift a VCF/BCF over to another reference build with a UCSC chain file (``LiftoverVcf``).

    ``reference_sequence`` is the target FASTA (with a ``.dict`` and ``.fai`` next to it). The
    output is sorted and indexed in the target coordinates unless ``disable_sort`` is set; variants
    that cannot be lifted go to ``reject_path``. ``liftover_min_match`` defaults to 1.0. With
    ``recover_swapped_ref_alt`` a bi-allelic SNP whose alternative allele is the new reference
    allele is swapped; ``tags_to_reverse`` (default ``AF``) are then transformed to 1-x and
    ``tags_to_drop`` (default ``MAX_AF``) removed.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    rejects = _p.new_file(reject_path, "reject_path")
    chain_file = _p.existing_file(chain, "chain")
    reference = _p.existing_file(reference_sequence, "reference_sequence")
    args = [
        f"INPUT={source}",
        f"OUTPUT={target}",
        f"REJECT={rejects}",
        f"CHAIN={chain_file}",
        f"REFERENCE_SEQUENCE={reference}",
    ]
    args += _p.a_num(
        "LIFTOVER_MIN_MATCH", liftover_min_match, "liftover_min_match", minimum=0, maximum=1
    )
    args += _p.a_bool("RECOVER_SWAPPED_REF_ALT", recover_swapped_ref_alt, "recover_swapped_ref_alt")
    args += _p.a_bool(
        "ALLOW_MISSING_FIELDS_IN_HEADER",
        allow_missing_fields_in_header,
        "allow_missing_fields_in_header",
    )
    args += _p.a_texts("TAGS_TO_DROP", tags_to_drop, "tags_to_drop")
    args += _p.a_texts("TAGS_TO_REVERSE", tags_to_reverse, "tags_to_reverse")
    args += _p.a_bool("WRITE_ORIGINAL_POSITION", write_original_position, "write_original_position")
    args += _p.a_bool("WRITE_ORIGINAL_ALLELES", write_original_alleles, "write_original_alleles")
    args += _p.a_bool("DISABLE_SORT", disable_sort, "disable_sort")
    args += _p.a_bool("WARN_ON_MISSING_CONTIG", warn_on_missing_contig, "warn_on_missing_contig")
    args += _p.a_bool("LOG_FAILED_INTERVALS", log_failed_intervals, "log_failed_intervals")
    outputs = {"output": target, "reject": rejects}
    result = _p.run(
        "LiftoverVcf",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    result["index"] = _index_of(target)
    if target.is_file():
        result["lifted_records"] = _vcf_records(target)
    if rejects.is_file():
        result["rejected_records"] = _vcf_records(rejects)
    return result


def picard_make_sites_only_vcf(
    input_path: str,
    output_path: str,
    *,
    samples: Sequence[str] | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove all genotype columns from a VCF/VCF.gz/BCF, keeping the site-level information and
    annotations such as AN and AF (``MakeSitesOnlyVcf``).

    ``samples`` keeps only the named samples' contribution to the site-level annotations (Picard's
    ``SAMPLE``).
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={source}", f"OUTPUT={target}"]
    args += _p.a_texts("SAMPLE", samples, "samples")
    result = _p.run(
        "MakeSitesOnlyVcf",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    result["index"] = _index_of(target)
    result["records"] = _vcf_records(target) if target.is_file() else None
    return result


def picard_make_vcf_sample_name_map(
    input_paths: Sequence[str],
    output_path: str,
    *,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a tab-separated map from sample name to VCF/GVCF path, one line per input
    (``MakeVcfSampleNameMap``).

    Each input must be a single-sample VCF/GVCF. The map is also returned as ``entries`` (a list of
    ``sample`` and ``path`` pairs). Picard 3.5.0 writes the path in the first column and the sample
    name in the second, although its documentation says the reverse; ``entries`` tells the two apart
    by the input paths.
    """
    sources = _p.existing_files(input_paths, "input_paths")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={path}" for path in sources] + [f"OUTPUT={target}"]
    result = _p.run(
        "MakeVcfSampleNameMap",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if target.is_file():
        known = {str(path) for path in sources}
        rows = [
            line.split("\t")
            for line in target.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        result["entries"] = [
            {"sample": row[1], "path": row[0]}
            if row[0] in known
            else {"sample": row[0], "path": row[1]}
            for row in rows
            if len(row) >= 2
        ]
    return result


def picard_rename_sample_in_vcf(
    input_path: str,
    output_path: str,
    new_sample_name: str,
    *,
    old_sample_name: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rename the sample of a single-sample VCF/BCF (``RenameSampleInVcf``).

    ``old_sample_name``, when given, asserts that this is the existing name.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    args = [
        f"INPUT={source}",
        f"OUTPUT={target}",
        f"NEW_SAMPLE_NAME={_p.text(new_sample_name, 'new_sample_name')}",
    ]
    args += _p.a_text("OLD_SAMPLE_NAME", old_sample_name, "old_sample_name")
    result = _p.run(
        "RenameSampleInVcf",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    result["index"] = _index_of(target)
    return result


def picard_split_vcfs(
    input_path: str,
    snp_output: str,
    indel_output: str,
    *,
    strict: bool | None = None,
    sequence_dictionary: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a VCF/BCF into a SNP file and an indel file with identical headers (``SplitVcfs``).

    Both outputs are indexed. With Picard's default ``strict=True`` a record that is neither a SNP
    nor an indel is an error; ``strict=False`` drops such records.
    """
    source = _p.existing_file(input_path, "input_path")
    snps = _p.new_file(snp_output, "snp_output")
    indels = _p.new_file(indel_output, "indel_output")
    args = [f"INPUT={source}", f"SNP_OUTPUT={snps}", f"INDEL_OUTPUT={indels}"]
    args += _p.a_bool("STRICT", strict, "strict")
    args += _p.a_file("SEQUENCE_DICTIONARY", sequence_dictionary, "sequence_dictionary")
    result = _p.run(
        "SplitVcfs",
        args,
        outputs={"snp_output": snps, "indel_output": indels},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if snps.is_file():
        result["snp_records"] = _vcf_records(snps)
    if indels.is_file():
        result["indel_records"] = _vcf_records(indels)
    return result


# ---------------------------------------------------------------------------
# Haplotype maps and interval lists
# ---------------------------------------------------------------------------


def picard_liftover_haplotype_map(
    input_path: str,
    output_path: str,
    chain: str,
    sequence_dictionary: str,
    *,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Lift a fingerprinting haplotype database over to another reference with a UCSC chain file
    (``LiftOverHaplotypeMap``).

    ``sequence_dictionary`` (any file a dictionary can be read from) is written into the output
    header.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    chain_file = _p.existing_file(chain, "chain")
    dictionary = _p.existing_file(sequence_dictionary, "sequence_dictionary")
    args = [
        f"INPUT={source}",
        f"OUTPUT={target}",
        f"CHAIN={chain_file}",
        f"SEQUENCE_DICTIONARY={dictionary}",
    ]
    result = _p.run(
        "LiftOverHaplotypeMap",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if target.is_file():
        result["snps"] = _p.count_lines(target, skip_prefix="@") - _p.count_lines(
            target, prefix="#"
        )
    return result


def picard_interval_list_to_bed(
    input_path: str,
    output_path: str,
    *,
    score: int | None = None,
    sort: bool | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a Picard interval list into a BED file (``IntervalListToBed``).

    ``score`` (0 to 1000, Picard's default 500) is written in the BED score column; ``sort``
    (default true) sorts the intervals first. Interval lists are 1-based closed, BED is 0-based
    half-open; Picard converts.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={source}", f"OUTPUT={target}"]
    args += _p.a_int("SCORE", score, "score", minimum=0, maximum=1000)
    args += _p.a_bool("SORT", sort, "sort")
    result = _p.run(
        "IntervalListToBed",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if target.is_file():
        result["intervals"] = _p.count_lines(target)
    return result


def picard_interval_list_tools(
    input_paths: Sequence[str],
    output_path: str,
    *,
    action: str = "CONCAT",
    second_input_paths: Sequence[str] | None = None,
    padding: int | None = None,
    invert: bool | None = None,
    sort: bool | None = None,
    unique: bool | None = None,
    dont_merge_abutting: bool | None = None,
    break_bands_at_multiples_of: int | None = None,
    include_filtered: bool | None = None,
    comment: Sequence[str] | None = None,
    count_value: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine and manipulate interval lists (``IntervalListTools``).

    ``action`` picks the set operation: ``CONCAT`` (all intervals of all inputs, no sorting or
    merging; Picard's default), ``UNION`` (merged, sorted), ``INTERSECT`` (loci present in every
    input), ``SUBTRACT`` (loci of the inputs not in ``second_input_paths``), ``SYMDIFF`` (loci in
    exactly one side) and ``OVERLAPS`` (the whole intervals of the inputs that overlap any interval
    of ``second_input_paths``). Inputs may be interval lists or VCFs. ``padding`` widens (or, when
    negative, shrinks) each interval before the action; ``invert`` returns the complement;
    ``count_value`` (``BASES`` or ``INTERVALS``) adds the number of bases or intervals of the output
    as ``count``.
    """
    sources = _p.existing_files(input_paths, "input_paths")
    target = _p.new_file(output_path, "output_path")
    mode = _p.choice(action, _INTERVAL_ACTIONS, "action")
    second = (
        _p.existing_files(second_input_paths, "second_input_paths") if second_input_paths else []
    )
    if mode in ("SUBTRACT", "SYMDIFF", "OVERLAPS") and not second:
        raise ConfigurationError(f"second_input_paths is required for action={mode}.")
    counted = None if count_value is None else _p.choice(count_value, _COUNT_VALUES, "count_value")
    count_file = target.with_name(target.name + ".count") if counted is not None else None
    if count_file is not None and count_file.exists():
        raise ConfigurationError(f"count_value would overwrite an existing file: {count_file}")
    args = [f"INPUT={path}" for path in sources] + [f"OUTPUT={target}", f"ACTION={mode}"]
    args += [f"SECOND_INPUT={path}" for path in second]
    args += _p.a_int("PADDING", padding, "padding")
    args += _p.a_bool("INVERT", invert, "invert")
    args += _p.a_bool("SORT", sort, "sort")
    args += _p.a_bool("UNIQUE", unique, "unique")
    args += _p.a_bool("DONT_MERGE_ABUTTING", dont_merge_abutting, "dont_merge_abutting")
    args += _p.a_int(
        "BREAK_BANDS_AT_MULTIPLES_OF",
        break_bands_at_multiples_of,
        "break_bands_at_multiples_of",
        minimum=0,
    )
    args += _p.a_bool("INCLUDE_FILTERED", include_filtered, "include_filtered")
    args += _p.a_texts("COMMENT", comment, "comment")
    if counted is not None and count_file is not None:
        args += [f"OUTPUT_VALUE={counted}", f"COUNT_OUTPUT={count_file}"]
    outputs = {"output": target, **({"count_output": count_file} if count_file is not None else {})}
    result = _p.run(
        "IntervalListTools",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if target.is_file():
        result["intervals"] = _p.count_lines(target, skip_prefix="@")
    if count_file is not None and count_file.is_file():
        text = count_file.read_text(encoding="utf-8").strip()
        result["count"] = int(text) if text.lstrip("-").isdigit() else None
    return result


def picard_scatter_intervals_by_ns(
    reference_sequence: str,
    output_path: str,
    *,
    output_type: str | None = None,
    max_to_merge: int | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a reference into alternating ``N`` and ``ACGT`` intervals, for scattering a
    variant-calling pipeline (``ScatterIntervalsByNs``).

    The FASTA needs an index (``.fai``) and a dictionary (``.dict``). ``output_type`` is ``BOTH``
    (Picard's default), ``N`` or ``ACGT``; ``max_to_merge`` (default 1) is the number of contiguous
    N bases tolerated inside an ACGT interval.
    """
    reference = _p.existing_file(reference_sequence, "reference_sequence")
    target = _p.new_file(output_path, "output_path")
    args = [f"REFERENCE={reference}", f"OUTPUT={target}"]
    args += _p.a_choice("OUTPUT_TYPE", output_type, _N_OUTPUT_TYPES, "output_type")
    args += _p.a_int("MAX_TO_MERGE", max_to_merge, "max_to_merge", minimum=0)
    result = _p.run(
        "ScatterIntervalsByNs",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if target.is_file():
        result["intervals"] = _p.count_lines(target, skip_prefix="@")
    return result


# ---------------------------------------------------------------------------
# FASTA and GFF
# ---------------------------------------------------------------------------


def picard_normalize_fasta(
    input_path: str,
    output_path: str,
    *,
    line_length: int | None = None,
    truncate_sequence_names_at_whitespace: bool | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rewrite a FASTA so that every sequence line has the same length, except the last of each
    record (``NormalizeFasta``).

    ``line_length`` defaults to Picard's 100 bases; ``truncate_sequence_names_at_whitespace``
    shortens record names at the first whitespace. Picard writes the sequence as it reads it (it
    does not change case).
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={source}", f"OUTPUT={target}"]
    args += _p.a_int("LINE_LENGTH", line_length, "line_length", minimum=1)
    args += _p.a_bool(
        "TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE",
        truncate_sequence_names_at_whitespace,
        "truncate_sequence_names_at_whitespace",
    )
    result = _p.run(
        "NormalizeFasta",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if target.is_file():
        result["records"] = _p.count_lines(target, prefix=">")
    return result


def picard_sort_gff(
    input_path: str,
    output_path: str,
    *,
    sequence_dictionary: str | None = None,
    n_records_in_memory: int | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort a GFF3 file by coordinates so that it can be indexed, adding flush directives where
    possible (``SortGff``).

    Contigs are ordered by ``sequence_dictionary`` when given, otherwise lexicographically.
    ``n_records_in_memory`` (Picard's ``nRecordsInMemory``, default 50000) is the number of records
    held in RAM before spilling to disk.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={source}", f"OUTPUT={target}"]
    args += _p.a_file("SEQUENCE_DICTIONARY", sequence_dictionary, "sequence_dictionary")
    args += _p.a_int("nRecordsInMemory", n_records_in_memory, "n_records_in_memory", minimum=1)
    result = _p.run(
        "SortGff",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if target.is_file():
        result["features"] = _p.count_lines(target, skip_prefix="#")
    return result


__all__ = [
    "picard_filter_vcf",
    "picard_interval_list_to_bed",
    "picard_interval_list_tools",
    "picard_liftover_haplotype_map",
    "picard_liftover_vcf",
    "picard_make_sites_only_vcf",
    "picard_make_vcf_sample_name_map",
    "picard_normalize_fasta",
    "picard_rename_sample_in_vcf",
    "picard_scatter_intervals_by_ns",
    "picard_sort_gff",
    "picard_split_vcfs",
]
_register_public_verbose(globals())
