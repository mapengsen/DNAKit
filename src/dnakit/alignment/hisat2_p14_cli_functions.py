"""HISAT2 command-line facades (tested with HISAT2 2.2.3): splice-aware read alignment, index inspection, repeat-index
building, annotation and variant extraction helpers, read statistics and read simulation.

The user installs HISAT2 (``conda install -c bioconda hisat2``, or the release archive from
https://daehwankimlab.github.io/hisat2/) and puts it on PATH, or passes ``hisat2_path`` (the folder that holds the
HISAT2 programs and its helper scripts). Each public function runs one fixed HISAT2 program or helper script with an
argument list built from typed, validated parameters: no shell, no generic runner, outputs are new files. Only the
options that the caller sets are passed, so every other option keeps HISAT2's own default. Results are
JSON-serialisable dictionaries (status, backend, version, command, return code, outputs, log tail).
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli, _p14_cli
from dnakit.exceptions import ConfigurationError

_HINT = (
    "Install HISAT2 2.2 or later (for example `conda install -c bioconda hisat2`, or the archive from "
    "https://daehwankimlab.github.io/hisat2/), then put it on PATH or pass hisat2_path (the folder that "
    "contains the HISAT2 programs and scripts)."
)
_INPUT_FORMATS = {"fastq": "-q", "fasta": "-f", "raw": "-r"}
_PRESETS = ("fast", "sensitive", "very-sensitive")
_STRANDNESS = ("F", "R", "FR", "RF")
_INSPECT_VIEWS = {
    "snps": ["--snp"],
    "splice_sites": ["--ss"],
    "splice_sites_all": ["--ss-all"],
    "exons": ["--exon"],
    "summary": ["-s"],
    "names": ["-n"],
    "sequences": [],
}
_REFERENCE_TYPES = ("genome", "chromosome", "gene")
_FUNCTION = re.compile(r"^[A-Z],-?\d+(?:\.\d+)?,-?\d+(?:\.\d+)?$")


def _executable(name: str, hisat2_path: str | None) -> str | None:
    if hisat2_path is None:
        return None
    folder = Path(hisat2_path).expanduser()
    if not folder.is_dir():
        raise ConfigurationError(f"hisat2_path must be a folder that contains the HISAT2 programs: {folder}")
    return str(folder / name)


def _run(
    program: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    hisat2_path: str | None,
    stdout_file: Path | None = None,
) -> dict[str, Any]:
    return _p14_cli.run(
        program,
        arguments,
        cwd=cwd,
        expected=expected,
        explicit=_executable(program, hisat2_path),
        hint=_HINT,
        stdout_file=stdout_file,
    )


def _index(prefix: str) -> str:
    path = Path(prefix).expanduser().resolve()
    if next(path.parent.glob(f"{path.name}.1.ht2*"), None) is None:
        raise ConfigurationError(f"No HISAT2 index with this prefix (no {path.name}.1.ht2 or .ht2l file): {path}")
    return str(path)


def _file_list(values: Sequence[str] | None, name: str) -> str | None:
    if values is None:
        return None
    if isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of file paths.")
    paths = [str(_gap_cli.existing_file(item, name)) for item in values]
    if any("," in item for item in paths):
        raise ConfigurationError(f"{name} paths must not contain commas.")
    return ",".join(paths)


def _new_base(value: str, name: str) -> Path:
    """A path prefix whose folder exists and under which nothing exists yet (``<prefix>`` or ``<prefix>.*``)."""
    base = Path(value).expanduser().resolve()
    if not base.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {base.parent}")
    if base.exists() or any(base.parent.glob(f"{base.name}.*")) or any(base.parent.glob(f"{base.name}_*")):
        raise ConfigurationError(f"{name} must not exist yet: {base}")
    return base


def _score_function(value: str, name: str) -> str:
    if not _FUNCTION.fullmatch(value):
        raise ConfigurationError(f"{name} must be a function such as 'G,-8,1' or 'L,0,-0.5'.")
    return value


def hisat2_align_reads(
    index_prefix: str,
    output_sam: str,
    *,
    unpaired_reads: Sequence[str] | None = None,
    reads_1: Sequence[str] | None = None,
    reads_2: Sequence[str] | None = None,
    read_format: str = "fastq",
    preset: str | None = None,
    threads: int | None = None,
    no_spliced_alignment: bool = False,
    known_splice_sites: str | None = None,
    novel_splice_sites_out: str | None = None,
    novel_splice_sites_in: str | None = None,
    min_intron_length: int | None = None,
    max_intron_length: int | None = None,
    rna_strandness: str | None = None,
    downstream_transcriptome_assembly: bool = False,
    dta_cufflinks: bool = False,
    transcriptome_mapping_only: bool = False,
    max_alignments: int | None = None,
    min_fragment_length: int | None = None,
    max_fragment_length: int | None = None,
    no_mixed: bool = False,
    no_discordant: bool = False,
    trim_5_prime: int | None = None,
    trim_3_prime: int | None = None,
    skip_reads: int | None = None,
    stop_after_reads: int | None = None,
    read_group_id: str | None = None,
    no_softclip: bool = False,
    no_repeat_index: bool = False,
    summary_file: str | None = None,
    seed: int | None = None,
    hisat2_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align reads to an indexed genome with splice-aware HISAT2 (``hisat2``) and write a SAM file.

    ``index_prefix`` is the prefix of an index built with hisat2-build. Give ``unpaired_reads`` and/or paired
    ``reads_1`` and ``reads_2`` (lists of files; gzip and bzip2 are read). ``read_format`` is 'fastq', 'fasta' or
    'raw'. ``preset`` is 'fast', 'sensitive' or 'very-sensitive'. ``no_spliced_alignment=True`` aligns DNA reads
    without splicing (``min_fragment_length``/``max_fragment_length`` then bound the fragment size).
    ``known_splice_sites`` is the file written by hisat2_extract_splice_sites; ``novel_splice_sites_out`` collects the
    junctions found. ``rna_strandness`` is 'F', 'R' (single-end) or 'FR', 'RF' (paired-end).
    ``downstream_transcriptome_assembly`` (``--dta``) and ``dta_cufflinks`` tailor the output to StringTie/Cufflinks;
    ``transcriptome_mapping_only`` is ``--tmo``. ``max_alignments`` is ``-k``. Options left as None keep HISAT2's
    defaults. Writes a new SAM file (and the optional summary and splice-site files).
    """
    index = _index(index_prefix)
    output_path = _gap_cli.new_file(output_sam, "output_sam")
    if (reads_1 is None) != (reads_2 is None):
        raise ConfigurationError("reads_1 and reads_2 must be given together.")
    if unpaired_reads is None and reads_1 is None:
        raise ConfigurationError("Give unpaired_reads, or reads_1 and reads_2.")
    arguments = ["-x", index, _INPUT_FORMATS[_gap_cli.choice(read_format, tuple(_INPUT_FORMATS), "read_format")]]
    if reads_1 is not None:
        arguments += ["-1", str(_file_list(reads_1, "reads_1")), "-2", str(_file_list(reads_2, "reads_2"))]
    if unpaired_reads is not None:
        arguments += ["-U", str(_file_list(unpaired_reads, "unpaired_reads"))]
    arguments += ["-S", str(output_path)]
    if preset is not None:
        arguments.append("--" + _gap_cli.choice(preset, _PRESETS, "preset"))
    if threads is not None:
        arguments += ["-p", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
    if _p14_cli.flag(no_spliced_alignment, "no_spliced_alignment"):
        arguments.append("--no-spliced-alignment")
    expected = [output_path]
    if known_splice_sites is not None:
        arguments += ["--known-splicesite-infile", str(_gap_cli.existing_file(known_splice_sites, "known_splice_sites"))]
    if novel_splice_sites_in is not None:
        arguments += ["--novel-splicesite-infile", str(_gap_cli.existing_file(novel_splice_sites_in, "novel_splice_sites_in"))]
    if novel_splice_sites_out is not None:
        novel = _gap_cli.new_file(novel_splice_sites_out, "novel_splice_sites_out")
        arguments += ["--novel-splicesite-outfile", str(novel)]
        expected.append(novel)
    if summary_file is not None:
        summary = _gap_cli.new_file(summary_file, "summary_file")
        arguments += ["--summary-file", str(summary)]
        expected.append(summary)
    for name, value, flag_name in (
        ("min_intron_length", min_intron_length, "--min-intronlen"),
        ("max_intron_length", max_intron_length, "--max-intronlen"),
        ("max_alignments", max_alignments, "-k"),
        ("min_fragment_length", min_fragment_length, "-I"),
        ("max_fragment_length", max_fragment_length, "-X"),
        ("trim_5_prime", trim_5_prime, "-5"),
        ("trim_3_prime", trim_3_prime, "-3"),
        ("skip_reads", skip_reads, "-s"),
        ("stop_after_reads", stop_after_reads, "-u"),
        ("seed", seed, "--seed"),
    ):
        if value is not None:
            arguments += [flag_name, str(_gap_cli.integer(value, name, minimum=0 if name not in {"max_alignments", "min_intron_length", "max_intron_length"} else 1))]
    if rna_strandness is not None:
        arguments += ["--rna-strandness", _gap_cli.choice(rna_strandness, _STRANDNESS, "rna_strandness")]
    for value, flag_name, argument_name in (
        (downstream_transcriptome_assembly, "--dta", "downstream_transcriptome_assembly"),
        (dta_cufflinks, "--dta-cufflinks", "dta_cufflinks"),
        (transcriptome_mapping_only, "--tmo", "transcriptome_mapping_only"),
        (no_mixed, "--no-mixed", "no_mixed"),
        (no_discordant, "--no-discordant", "no_discordant"),
        (no_softclip, "--no-softclip", "no_softclip"),
        (no_repeat_index, "--no-repeat-index", "no_repeat_index"),
    ):
        if _p14_cli.flag(value, argument_name):
            arguments.append(flag_name)
    if read_group_id is not None:
        arguments += ["--rg-id", _gap_cli.text(read_group_id, "read_group_id")]
    return _run("hisat2", arguments, cwd=output_path.parent, expected=expected, hisat2_path=hisat2_path)


def hisat2_inspect_index(
    index_prefix: str,
    output: str,
    *,
    view: str = "snps",
    hisat2_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Print what is stored in a HISAT2 index (``hisat2-inspect``) into a new file.

    ``view`` selects the content: 'snps' (``--snp``, the SNPs of a graph index: id, type, chromosome, position,
    allele), 'splice_sites' (``--ss``), 'splice_sites_all' (``--ss-all``, including sites not in the global index),
    'exons' (``--exon``), 'summary' (reference names, lengths and index parameters), 'names' (reference names only)
    or 'sequences' (the indexed sequences as FASTA). ``rows`` holds the tab-separated records of the four annotation
    views.
    """
    index = _index(index_prefix)
    output_path = _gap_cli.new_file(output, "output")
    chosen = _gap_cli.choice(view, tuple(_INSPECT_VIEWS), "view")
    result = _run(
        "hisat2-inspect",
        [*_INSPECT_VIEWS[chosen], index],
        cwd=output_path.parent,
        expected=[],
        hisat2_path=hisat2_path,
        stdout_file=output_path,
    )
    result["view"] = chosen
    if result["status"] == "completed" and chosen in {"snps", "splice_sites", "splice_sites_all", "exons"}:
        rows = _p14_cli.tsv_rows(output_path)
        result["rows"] = rows
        result["row_count"] = len(rows)
    return result


def hisat2_build_repeat_index(
    reference: str,
    output_prefix: str,
    *,
    threads: int | None = None,
    seed_length: int | None = None,
    seed_count: int | None = None,
    min_repeat_length: int | None = None,
    max_repeat_length: int | None = None,
    repeat_count: int | None = None,
    repeat_edit_distance: int | None = None,
    forward_only: bool = False,
    hisat2_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find repeat sequences in a reference genome and write the files of a HISAT2 repeat index (``hisat2-repeat``).

    ``reference`` is a FASTA file; the outputs start with ``output_prefix`` (a new prefix). ``seed_length`` and
    ``seed_count`` control the seed k-mers (HISAT2 defaults 50 and 5), ``min_repeat_length`` and
    ``max_repeat_length`` the repeat length range (100 and 65535), ``repeat_count`` the minimum number of copies (5),
    ``repeat_edit_distance`` the largest edit distance between copies (10), and ``forward_only`` uses the forward
    strand only. The repeat information is used by ``hisat2-build`` to make a repeat-aware index.
    """
    source = _gap_cli.existing_file(reference, "reference")
    base = _new_base(output_prefix, "output_prefix")
    arguments: list[str] = []
    if threads is not None:
        arguments += ["-p", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256))]
    for name, value, flag_name in (
        ("seed_length", seed_length, "--seed-length"),
        ("seed_count", seed_count, "--seed-count"),
        ("min_repeat_length", min_repeat_length, "--min-repeat-length"),
        ("max_repeat_length", max_repeat_length, "--max-repeat-length"),
        ("repeat_count", repeat_count, "--repeat-count"),
        ("repeat_edit_distance", repeat_edit_distance, "--repeat-edit"),
    ):
        if value is not None:
            arguments += [flag_name, str(_gap_cli.integer(value, name, minimum=0 if name == "repeat_edit_distance" else 1))]
    if _p14_cli.flag(forward_only, "forward_only"):
        arguments.append("--forward-only")
    arguments += [str(source), str(base)]
    result = _run("hisat2-repeat", arguments, cwd=base.parent, expected=[], hisat2_path=hisat2_path)
    produced = sorted(str(path) for path in base.parent.glob(f"{base.name}*") if path.is_file())
    result["outputs"] = produced
    if result["status"] == "completed" and not produced:
        result["status"] = "failed"
    return result


def hisat2_extract_splice_sites(
    annotation_gtf: str,
    output: str,
    *,
    hisat2_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract splice sites from a gene annotation in GTF format (``hisat2_extract_splice_sites.py``).

    The output has one line per intron: chromosome, last base of the upstream exon, first base of the downstream
    exon (0-based coordinates) and the strand. It is the ``known_splice_sites`` input of hisat2_align_reads and
    of ``hisat2-build --ss``. Writes a new file.
    """
    source = _gap_cli.existing_file(annotation_gtf, "annotation_gtf")
    output_path = _gap_cli.new_file(output, "output")
    result = _run(
        "hisat2_extract_splice_sites.py", [str(source)], cwd=output_path.parent, expected=[],
        hisat2_path=hisat2_path, stdout_file=output_path,
    )
    return _with_rows(result, output_path)


def hisat2_extract_exons(
    annotation_gtf: str,
    output: str,
    *,
    hisat2_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract exons from a gene annotation in GTF format (``hisat2_extract_exons.py``).

    The output has one line per merged exon: chromosome, start, end (0-based) and strand; it is the input of
    ``hisat2-build --exon``. Writes a new file.
    """
    source = _gap_cli.existing_file(annotation_gtf, "annotation_gtf")
    output_path = _gap_cli.new_file(output, "output")
    result = _run(
        "hisat2_extract_exons.py", [str(source)], cwd=output_path.parent, expected=[],
        hisat2_path=hisat2_path, stdout_file=output_path,
    )
    return _with_rows(result, output_path)


def _with_rows(result: dict[str, Any], output_path: Path) -> dict[str, Any]:
    if result["status"] == "completed":
        rows = _p14_cli.tsv_rows(output_path)
        result["row_count"] = len(rows)
        result["rows"] = rows[:1000]
    return result


def hisat2_extract_snps_haplotypes_ucsc(
    genome_fasta: str,
    snp_table: str,
    output_base: str,
    *,
    inter_gap: int | None = None,
    intra_gap: int | None = None,
    hisat2_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract SNPs and haplotypes from a UCSC dbSNP table for ``hisat2-build`` (``hisat2_extract_snps_haplotypes_UCSC.py``).

    ``snp_table`` is a UCSC dbSNP table such as ``snp144Common.txt`` (plain or gzip). ``inter_gap`` is the largest
    distance between variants of one haplotype and ``intra_gap`` the distance at which a haplotype is split (script
    defaults apply when None). Writes ``<output_base>.snp`` and ``<output_base>.haplotype``.
    """
    genome = _gap_cli.existing_file(genome_fasta, "genome_fasta")
    table = _gap_cli.existing_file(snp_table, "snp_table")
    base = _new_base(output_base, "output_base")
    arguments: list[str] = []
    if inter_gap is not None:
        arguments += ["--inter-gap", str(_gap_cli.integer(inter_gap, "inter_gap", minimum=0))]
    if intra_gap is not None:
        arguments += ["--intra-gap", str(_gap_cli.integer(intra_gap, "intra_gap", minimum=0))]
    arguments += [str(genome), str(table), str(base)]
    return _glob_result(
        _run("hisat2_extract_snps_haplotypes_UCSC.py", arguments, cwd=base.parent, expected=[], hisat2_path=hisat2_path),
        base,
    )


def hisat2_extract_snps_haplotypes_vcf(
    genome_fasta: str,
    vcf_files: Sequence[str],
    output_base: str,
    *,
    reference_type: str | None = None,
    inter_gap: int | None = None,
    intra_gap: int | None = None,
    non_rs_ids: bool = False,
    genotype_vcf: str | None = None,
    genotype_gene_list: Sequence[str] | None = None,
    extra_files: bool = False,
    hisat2_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract SNPs and haplotypes from VCF files for ``hisat2-build`` (``hisat2_extract_snps_haplotypes_VCF.py``).

    ``vcf_files`` are plain or gzip VCF files (for example dbSNP). ``reference_type`` is 'genome' (script default),
    'chromosome' or 'gene'; ``inter_gap`` (script default 30) and ``intra_gap`` (50) shape the haplotypes;
    ``non_rs_ids`` allows variant ids that do not start with ``rs``; ``genotype_vcf`` and ``genotype_gene_list``
    prepare a HISAT2-genotype reference; ``extra_files`` also writes the backbone FASTA and ``.ref`` files. Writes
    ``<output_base>.snp``, ``<output_base>.haplotype`` and the optional extra files.
    """
    genome = _gap_cli.existing_file(genome_fasta, "genome_fasta")
    vcfs = _file_list(vcf_files, "vcf_files")
    base = _new_base(output_base, "output_base")
    arguments: list[str] = []
    if reference_type is not None:
        arguments += ["--reference-type", _gap_cli.choice(reference_type, _REFERENCE_TYPES, "reference_type")]
    if inter_gap is not None:
        arguments += ["--inter-gap", str(_gap_cli.integer(inter_gap, "inter_gap", minimum=0))]
    if intra_gap is not None:
        arguments += ["--intra-gap", str(_gap_cli.integer(intra_gap, "intra_gap", minimum=0))]
    if _p14_cli.flag(non_rs_ids, "non_rs_ids"):
        arguments.append("--non-rs")
    if genotype_vcf is not None:
        arguments += ["--genotype-vcf", str(_gap_cli.existing_file(genotype_vcf, "genotype_vcf"))]
    if genotype_gene_list is not None:
        if isinstance(genotype_gene_list, str) or not genotype_gene_list:
            raise ConfigurationError("genotype_gene_list must be a non-empty list of gene names.")
        arguments += ["--genotype-gene-list", ",".join(_gap_cli.text(item, "genotype_gene_list") for item in genotype_gene_list)]
    if _p14_cli.flag(extra_files, "extra_files"):
        arguments.append("--extra-files")
    arguments += [str(genome), str(vcfs), str(base)]
    return _glob_result(
        _run("hisat2_extract_snps_haplotypes_VCF.py", arguments, cwd=base.parent, expected=[], hisat2_path=hisat2_path),
        base,
    )


def _glob_result(result: dict[str, Any], base: Path) -> dict[str, Any]:
    found = [*base.parent.glob(f"{base.name}.*"), *base.parent.glob(f"{base.name}_*")]
    produced = sorted({str(path) for path in found if path.is_file()})
    result["outputs"] = produced
    if result["status"] == "completed" and not produced:
        result["status"] = "failed"
    return result


def hisat2_read_statistics(
    reads: str,
    *,
    read_count: int | None = None,
    hisat2_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count reads and report the minimum, maximum and average read length (``hisat2_read_statistics.py``).

    Only the first ``read_count`` reads are examined (script default 10000). The result holds ``reads_examined``,
    ``minimum_length``, ``maximum_length``, ``average_length`` and the sorted distinct ``read_lengths``.
    """
    source = _gap_cli.existing_file(reads, "reads")
    arguments: list[str] = []
    if read_count is not None:
        arguments += ["-n", str(_gap_cli.integer(read_count, "read_count", minimum=1))]
    arguments.append(str(source))
    import tempfile

    with tempfile.TemporaryDirectory(prefix="dnakit_hisat2_") as folder:
        work = Path(folder)
        result = _run(
            "hisat2_read_statistics.py", arguments, cwd=work, expected=[], hisat2_path=hisat2_path,
            stdout_file=work / "stats.txt",
        )
        text = (work / "stats.txt").read_text(encoding="utf-8").strip()
    result["outputs"] = []
    fields = text.split()
    if result["status"] == "completed" and len(fields) >= 4:
        result.update(
            reads_examined=int(fields[0]),
            minimum_length=int(fields[1]),
            maximum_length=int(fields[2]),
            average_length=int(fields[3]),
            read_lengths=[int(item) for item in fields[4].split(",")] if len(fields) > 4 else [],
        )
    elif result["status"] == "completed":
        result["status"] = "failed"
    return result


def hisat2_simulate_reads(
    genome_fasta: str,
    annotation_gtf: str,
    snp_file: str,
    output_base: str,
    *,
    dna: bool = False,
    single_end: bool = False,
    read_length: int | None = None,
    fragment_length: int | None = None,
    num_fragments: int | None = None,
    expression_profile: str | None = None,
    repeat_info: str | None = None,
    error_rate: float | None = None,
    max_mismatch: int | None = None,
    random_seed: int | None = None,
    snp_probability: float | None = None,
    python_path: str | None = None,
    hisat2_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate RNA-seq reads (default) or DNA-seq reads (``dna=True``) from a genome, annotation and SNPs
    (``hisat2_simulate_reads.py``).

    ``annotation_gtf`` supplies the transcripts (RNA-seq) and ``snp_file`` the variants (the ``.snp`` file written by
    hisat2_extract_snps_haplotypes_vcf or _ucsc). Paired-end reads are simulated unless ``single_end``.
    ``read_length`` (script default 100), ``fragment_length`` (250) and ``num_fragments`` (1,000,000) set the sample;
    ``expression_profile`` is 'flux' (default) or 'constant'; ``error_rate`` is the per-base sequencing error in
    percent, ``max_mismatch`` the most errors per read (3), ``snp_probability`` the chance that a read spanning a SNP
    carries it. Writes ``<output_base>_1.fa``, ``<output_base>_2.fa`` (paired) and the matching ``.sam`` truth files.

    Upstream limits: the script keeps only the chromosomes named 1-22, X and Y (human style, without ``chr``), and
    its RNA-seq mode calls ``random.shuffle(..., random=...)``, which Python 3.11 and later removed; run the RNA-seq
    mode with an older interpreter by passing ``python_path`` (for example Python 3.8). DNA mode (``dna=True``)
    works with current Python.
    """
    genome = _gap_cli.existing_file(genome_fasta, "genome_fasta")
    annotation = _gap_cli.existing_file(annotation_gtf, "annotation_gtf")
    snps = _gap_cli.existing_file(snp_file, "snp_file")
    base = _new_base(output_base, "output_base")
    arguments: list[str] = []
    if _p14_cli.flag(dna, "dna"):
        arguments.append("--dna")
    if _p14_cli.flag(single_end, "single_end"):
        arguments.append("--single-end")
    for name, value, flag_name in (
        ("read_length", read_length, "--read-length"),
        ("fragment_length", fragment_length, "--fragment-length"),
        ("num_fragments", num_fragments, "--num-fragment"),
        ("max_mismatch", max_mismatch, "--max-mismatch"),
        ("random_seed", random_seed, "--random-seed"),
    ):
        if value is not None:
            arguments += [flag_name, str(_gap_cli.integer(value, name, minimum=0 if name in {"max_mismatch", "random_seed"} else 1))]
    if expression_profile is not None:
        arguments += ["--expr-profile", _gap_cli.choice(expression_profile, ("flux", "constant"), "expression_profile")]
    if repeat_info is not None:
        arguments += ["--repeat-info", str(_gap_cli.existing_file(repeat_info, "repeat_info"))]
    if error_rate is not None:
        arguments += ["--error-rate", str(_gap_cli.number(error_rate, "error_rate", minimum=0.0, maximum=100.0))]
    if snp_probability is not None:
        arguments += ["--snp-prob", str(_gap_cli.number(snp_probability, "snp_probability", minimum=0.0, maximum=1.0))]
    arguments += [str(genome), str(annotation), str(snps), str(base)]
    if python_path is None:
        result = _run("hisat2_simulate_reads.py", arguments, cwd=base.parent, expected=[], hisat2_path=hisat2_path)
    else:
        script = _gap_cli.program("hisat2_simulate_reads.py", explicit=_executable("hisat2_simulate_reads.py", hisat2_path), hint=_HINT)
        result = _p14_cli.run(
            "python3", [script, *arguments], cwd=base.parent, expected=[], explicit=python_path,
            hint="Pass python_path as an existing Python interpreter.", program_name="hisat2_simulate_reads.py",
        )
    result = _glob_result(result, base)
    result["read_type"] = "dna" if dna else "rna"
    return result


__all__ = [
    "hisat2_align_reads",
    "hisat2_build_repeat_index",
    "hisat2_extract_exons",
    "hisat2_extract_snps_haplotypes_ucsc",
    "hisat2_extract_snps_haplotypes_vcf",
    "hisat2_extract_splice_sites",
    "hisat2_inspect_index",
    "hisat2_read_statistics",
    "hisat2_simulate_reads",
]
_register_public_verbose(globals())
