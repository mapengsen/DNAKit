"""umi_tools command-line operations as typed DNAKit facades: UMI extraction into read names, cell-barcode whitelists,
UMI grouping and deduplication of BAM files, UMI counting per gene, counting from a read-gene table, and BAM
preparation for RSEM.

The user installs umi_tools (``python -m pip install umi_tools``, or ``conda install -c bioconda umi_tools``); the
``bio-umi-tools`` extra installs it. Each function runs one umi_tools sub-command with a fixed flag list built from its
typed arguments: no shell, and no generic runner. Outputs are new files. Results are JSON-serialisable dictionaries
from ``dnakit.bioinformatics._gap_cli.run`` (status, backend, version, command, return code, outputs, log tail).
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-umi-tools"
_HINT = "Install umi_tools 1.1.6 or later: python -m pip install umi_tools (or conda install -c bioconda umi_tools)."
_UMI_METHODS = ("adjacency", "directional", "percentile", "unique", "cluster")
_EXTRACT_UMI_METHODS = ("read_id", "tag", "umis")
_EXTRACT_METHODS = ("string", "regex")
_WHITELIST_METHODS = ("reads", "umis")
_KNEE_METHODS = ("distance", "density")
_STRING_PATTERN = re.compile(r"[CNX]+")


def _token(value: object, name: str) -> str:
    """Non-empty text without NUL or newline; always passed in the --name=value form, so it may start with '-'."""
    if not isinstance(value, str) or not value or "\x00" in value or "\n" in value:
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _run(sub: str, arguments: Sequence[str], *, cwd: Path, expected: Sequence[Path], explicit: str | None) -> dict[str, Any]:
    return _gap_cli.run(
        "umi_tools",
        [sub, *arguments],
        cwd=cwd,
        expected=list(expected),
        explicit=explicit,
        hint=_HINT,
    )


def _umi_options(
    extract_umi_method: str,
    umi_tag: str | None,
    umi_separator: str,
    method: str,
    edit_distance_threshold: int,
    paired: bool,
) -> list[str]:
    """UMI source, grouping method and edit distance: the options shared by count, dedup and group."""
    source = _gap_cli.choice(extract_umi_method, _EXTRACT_UMI_METHODS, "extract_umi_method")
    options = [f"--extract-umi-method={source}"]
    if source == "tag":
        options.append(f"--umi-tag={_token(umi_tag or 'RX', 'umi_tag')}")
    elif source == "read_id":
        options.append(f"--umi-separator={_token(umi_separator, 'umi_separator')}")
    options.append(f"--method={_gap_cli.choice(method, _UMI_METHODS, 'method')}")
    options.append(
        f"--edit-distance-threshold={_gap_cli.integer(edit_distance_threshold, 'edit_distance_threshold', minimum=0)}"
    )
    if _flag(paired, "paired"):
        options.append("--paired")
    return options


def umi_extract_to_names(
    fastq: str,
    output: str,
    *,
    bc_pattern: str,
    extract_method: str = "string",
    umi_separator: str = "_",
    three_prime: bool = False,
    read2_in: str | None = None,
    read2_out: str | None = None,
    bc_pattern2: str | None = None,
    umi_tools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Move the UMI (and cell barcode) from the start of each read into its read name (umi_tools extract).

    bc_pattern describes the barcode at the start of read 1: C for cell-barcode bases, N for UMI bases, X for bases
    to discard (for example 'CCCCCCNNNN'). With extract_method='regex' bc_pattern is a regular expression with
    named groups instead. For paired reads give read2_in and read2_out (and bc_pattern2 when read 2 carries a
    barcode). three_prime says the barcode is at the 3' end. Writes new FASTQ files.
    """
    source = _gap_cli.existing_file(fastq, "fastq")
    output_path = _gap_cli.new_file(output, "output")
    method = _gap_cli.choice(extract_method, _EXTRACT_METHODS, "extract_method")
    pattern = _token(bc_pattern, "bc_pattern")
    if method == "string" and (not _STRING_PATTERN.fullmatch(pattern) or "N" not in pattern):
        raise ConfigurationError("bc_pattern for extract_method='string' uses C, N and X and must contain an N.")
    arguments = [
        f"--extract-method={method}",
        f"--bc-pattern={pattern}",
        f"--umi-separator={_token(umi_separator, 'umi_separator')}",
    ]
    if _flag(three_prime, "three_prime"):
        arguments.append("--3prime")
    expected = [output_path]
    if read2_in is not None:
        read2_source = _gap_cli.existing_file(read2_in, "read2_in")
        if read2_out is None:
            raise ConfigurationError("read2_out is required with read2_in.")
        read2_path = _gap_cli.new_file(read2_out, "read2_out")
        arguments += [f"--read2-in={read2_source}", f"--read2-out={read2_path}"]
        if bc_pattern2 is not None:
            arguments.append(f"--bc-pattern2={_token(bc_pattern2, 'bc_pattern2')}")
        expected.append(read2_path)
    elif read2_out is not None or bc_pattern2 is not None:
        raise ConfigurationError("read2_out and bc_pattern2 need read2_in.")
    arguments += [f"--stdin={source}", f"--stdout={output_path}"]
    return _run("extract", arguments, cwd=output_path.parent, expected=expected, explicit=umi_tools_path)


def umi_whitelist_cells(
    fastq: str,
    output: str,
    *,
    bc_pattern: str,
    extract_method: str = "string",
    method: str = "reads",
    knee_method: str = "distance",
    expect_cells: int | None = None,
    set_cell_number: int | None = None,
    umi_tools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the cell barcodes present in a FASTQ file (umi_tools whitelist) and write them as a whitelist.

    bc_pattern and extract_method are as in umi_extract_to_names. method counts reads ('reads') or unique UMIs
    ('umis') per barcode; knee_method picks the cut-off between real cells and background ('distance' is the
    upstream default, 'density' needs expect_cells); set_cell_number or expect_cells fix the number of cells.
    The plot options of umi_tools are not used. Writes the whitelist to a new file.
    """
    source = _gap_cli.existing_file(fastq, "fastq")
    output_path = _gap_cli.new_file(output, "output")
    method_name = _gap_cli.choice(extract_method, _EXTRACT_METHODS, "extract_method")
    pattern = _token(bc_pattern, "bc_pattern")
    if method_name == "string" and (not _STRING_PATTERN.fullmatch(pattern) or "C" not in pattern):
        raise ConfigurationError("bc_pattern for extract_method='string' uses C, N and X and must contain a C.")
    arguments = [
        f"--extract-method={method_name}",
        f"--bc-pattern={pattern}",
        f"--method={_gap_cli.choice(method, _WHITELIST_METHODS, 'method')}",
        f"--knee-method={_gap_cli.choice(knee_method, _KNEE_METHODS, 'knee_method')}",
    ]
    if expect_cells is not None:
        arguments.append(f"--expect-cells={_gap_cli.integer(expect_cells, 'expect_cells', minimum=1)}")
    if set_cell_number is not None:
        arguments.append(f"--set-cell-number={_gap_cli.integer(set_cell_number, 'set_cell_number', minimum=1)}")
    arguments += [f"--stdin={source}", f"--stdout={output_path}"]
    return _run("whitelist", arguments, cwd=output_path.parent, expected=[output_path], explicit=umi_tools_path)


def umi_group_bam(
    bam: str,
    output_bam: str,
    *,
    group_tsv: str | None = None,
    extract_umi_method: str = "read_id",
    umi_tag: str | None = None,
    umi_separator: str = "_",
    umi_group_tag: str = "BX",
    method: str = "directional",
    edit_distance_threshold: int = 1,
    paired: bool = False,
    umi_tools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group the reads of a BAM file by UMI and mapping position (umi_tools group).

    The output BAM carries each read's group in the umi_group_tag tag (BX by default). extract_umi_method is
    'read_id' (UMI after umi_separator in the read name), 'tag' (UMI in the tag umi_tag, RX by default) or 'umis'.
    method and edit_distance_threshold set the UMI clustering. group_tsv optionally receives the per-read table
    (read_id, contig, position, gene, umi, umi_count, final_umi, final_umi_count, unique_id). Writes new files.
    """
    source = _gap_cli.existing_file(bam, "bam")
    output_path = _gap_cli.new_file(output_bam, "output_bam")
    arguments = _umi_options(extract_umi_method, umi_tag, umi_separator, method, edit_distance_threshold, paired)
    arguments.append(f"--umi-group-tag={_token(umi_group_tag, 'umi_group_tag')}")
    arguments.append("--output-bam")
    expected = [output_path]
    if group_tsv is not None:
        tsv = _gap_cli.new_file(group_tsv, "group_tsv")
        arguments.append(f"--group-out={tsv}")
        expected.append(tsv)
    arguments += [f"--stdin={source}", f"--stdout={output_path}"]
    return _run("group", arguments, cwd=output_path.parent, expected=expected, explicit=umi_tools_path)


def umi_dedup_bam(
    bam: str,
    output_bam: str,
    *,
    extract_umi_method: str = "read_id",
    umi_tag: str | None = None,
    umi_separator: str = "_",
    method: str = "directional",
    edit_distance_threshold: int = 1,
    paired: bool = False,
    ignore_umi: bool = False,
    per_contig: bool = False,
    gene_tag: str | None = None,
    umi_tools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove PCR duplicates from a coordinate-sorted BAM file, keeping one read per UMI group and position.

    The UMI source and clustering are as in umi_group_bam. ignore_umi deduplicates on position only. per_contig
    or gene_tag (a BAM tag naming each read's gene) deduplicates within contigs or genes. The output is a new BAM.
    """
    source = _gap_cli.existing_file(bam, "bam")
    output_path = _gap_cli.new_file(output_bam, "output_bam")
    arguments = _umi_options(extract_umi_method, umi_tag, umi_separator, method, edit_distance_threshold, paired)
    if _flag(ignore_umi, "ignore_umi"):
        arguments.append("--ignore-umi")
    if _flag(per_contig, "per_contig"):
        arguments.append("--per-contig")
    if gene_tag is not None:
        arguments += ["--per-gene", f"--gene-tag={_token(gene_tag, 'gene_tag')}"]
    arguments += [f"--stdin={source}", f"--stdout={output_path}"]
    return _run("dedup", arguments, cwd=output_path.parent, expected=[output_path], explicit=umi_tools_path)


def umi_count_per_gene(
    bam: str,
    output: str,
    *,
    gene_tag: str | None = None,
    per_contig: bool = False,
    extract_umi_method: str = "read_id",
    umi_tag: str | None = None,
    umi_separator: str = "_",
    method: str = "directional",
    edit_distance_threshold: int = 1,
    mapping_quality: int = 0,
    paired: bool = False,
    umi_tools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count UMI groups per gene in a BAM file (umi_tools count); writes a gene TSV with columns gene and count.

    umi_tools count always counts per gene: name the gene with gene_tag (a BAM tag), or pass per_contig=True to use
    the contig as the gene. Reads are grouped by UMI as in umi_group_bam and the groups are counted. mapping_quality
    drops reads below that MAPQ. Writes a new TSV file.
    """
    source = _gap_cli.existing_file(bam, "bam")
    output_path = _gap_cli.new_file(output, "output")
    if per_contig and gene_tag is not None:
        raise ConfigurationError("Give either gene_tag or per_contig, not both.")
    arguments = _umi_options(extract_umi_method, umi_tag, umi_separator, method, edit_distance_threshold, paired)
    if _flag(per_contig, "per_contig"):
        arguments.append("--per-contig")
    elif gene_tag is not None:
        arguments.append(f"--gene-tag={_token(gene_tag, 'gene_tag')}")
    else:
        raise ConfigurationError("umi_count_per_gene needs gene_tag (the BAM tag that names each read's gene) or per_contig=True.")
    arguments.append(f"--mapping-quality={_gap_cli.integer(mapping_quality, 'mapping_quality', minimum=0)}")
    arguments += [f"--stdin={source}", f"--stdout={output_path}"]
    return _run("count", arguments, cwd=output_path.parent, expected=[output_path], explicit=umi_tools_path)


def umi_count_table(
    table: str,
    output: str,
    *,
    per_cell: bool = False,
    barcode_separator: str = "_",
    method: str = "directional",
    edit_distance_threshold: int = 1,
    umi_tools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count UMI groups per gene from a two-column tab-separated table (umi_tools count_tab).

    Each line is a read id (ending in the UMI after barcode_separator) and the gene that read was assigned to.
    The UMIs of each gene are clustered with method and edit_distance_threshold, and the number of groups is written
    as gene and count. per_cell also reads the cell barcode: the read id must end in
    cell_barcode<sep>UMI (the code takes the field before the last one as the cell; the help text of the program
    says the reverse). Writes cell, gene, count. Writes a new TSV file.
    """
    source = _gap_cli.existing_file(table, "table")
    output_path = _gap_cli.new_file(output, "output")
    arguments = [
        f"--barcode-separator={_token(barcode_separator, 'barcode_separator')}",
        f"--method={_gap_cli.choice(method, _UMI_METHODS, 'method')}",
        f"--edit-distance-threshold={_gap_cli.integer(edit_distance_threshold, 'edit_distance_threshold', minimum=0)}",
    ]
    if _flag(per_cell, "per_cell"):
        arguments.append("--per-cell")
    arguments += [f"--stdin={source}", f"--stdout={output_path}"]
    return _run("count_tab", arguments, cwd=output_path.parent, expected=[output_path], explicit=umi_tools_path)


def umi_prepare_bam_for_rsem(
    bam: str,
    output_bam: str,
    *,
    tags: str | None = None,
    sam: bool = False,
    umi_tools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Prepare a transcriptome-aligned BAM for RSEM (umi_tools prepare-for-rsem).

    The program pairs the alignments of each read name, copies the comma-separated tags given in tags from read 1
    to read 2, and writes a BAM in the form RSEM expects. Single-end input produced a BAM with no reads when run
    (checked), so the input must be paired. sam=True reads and writes SAM text. Writes a new file.
    """
    source = _gap_cli.existing_file(bam, "bam")
    output_path = _gap_cli.new_file(output_bam, "output_bam")
    arguments: list[str] = []
    if tags is not None:
        arguments.append(f"--tags={_token(tags, 'tags')}")
    if _flag(sam, "sam"):
        arguments.append("--sam")
    arguments += [f"--stdin={source}", f"--stdout={output_path}"]
    return _run("prepare-for-rsem", arguments, cwd=output_path.parent, expected=[output_path], explicit=umi_tools_path)


__all__ = [
    "umi_count_per_gene",
    "umi_count_table",
    "umi_dedup_bam",
    "umi_extract_to_names",
    "umi_group_bam",
    "umi_prepare_bam_for_rsem",
    "umi_whitelist_cells",
]
_register_public_verbose(globals())
