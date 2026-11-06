"""WhatsHap command-line facades (tested with WhatsHap 2.8): pedigree and polyploid phasing, phasing statistics and
comparison, un-phasing, haplotagging, read-based and pedigree genotyping, SNV candidate discovery, read splitting by
haplotype, error-profile learning and phasing a VCF from haplotagged reads.

The user installs WhatsHap (``python -m pip install whatshap``, or ``conda install -c bioconda whatshap``) and puts
``whatshap`` on PATH, or passes ``whatshap_path``. Each public function runs one fixed ``whatshap`` sub-command with an
argument list built from typed, validated parameters: no shell, no generic runner, outputs are new files. Only the
options that the caller sets are passed, so every other option keeps WhatsHap's own default. Results are
JSON-serialisable dictionaries (status, backend, version, command, return code, outputs, log tail). Plain
``whatshap phase`` without a pedigree is ``phase_variants_whatshap``.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli, _p14_cli
from dnakit.exceptions import ConfigurationError

_HINT = "Install WhatsHap 2.0 or later (python -m pip install whatshap, or conda install -c bioconda whatshap) and put it on PATH or pass whatshap_path."
_TAGS = ("PS", "HP")
_ALGORITHMS = ("whatshap", "hapchat", "heuristic")
_TECHNOLOGIES = ("illumina", "pacbio", "nanopore")
_SUPPLEMENTARY = ("skip", "copy-primary", "independent-or-skip", "independent-or-copy-primary")


def _whatshap(
    sub: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    whatshap_path: str | None,
    stdout_file: Path | None = None,
) -> dict[str, Any]:
    return _p14_cli.run(
        "whatshap",
        [sub, *arguments],
        cwd=cwd,
        expected=expected,
        explicit=whatshap_path,
        hint=_HINT,
        stdout_file=stdout_file,
    )


def _files(values: Sequence[str], name: str) -> list[str]:
    if isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of files.")
    return [str(_gap_cli.existing_file(item, name)) for item in values]


def _common(
    sample: str | None = None,
    chromosome: str | None = None,
    exclude_chromosomes: Sequence[str] | None = None,
    tag: str | None = None,
    only_snvs: bool = False,
    distrust_genotypes: bool = False,
) -> list[str]:
    arguments: list[str] = []
    if sample is not None:
        arguments += ["--sample", _gap_cli.text(sample, "sample")]
    if chromosome is not None:
        arguments += ["--chromosome", _gap_cli.text(chromosome, "chromosome")]
    for item in exclude_chromosomes or ():
        arguments += ["--exclude-chromosome", _gap_cli.text(item, "exclude_chromosomes")]
    if tag is not None:
        arguments += ["--tag", _gap_cli.choice(tag, _TAGS, "tag")]
    if _p14_cli.flag(only_snvs, "only_snvs"):
        arguments.append("--only-snvs")
    if _p14_cli.flag(distrust_genotypes, "distrust_genotypes"):
        arguments.append("--distrust-genotypes")
    return arguments


def _reference(reference: str | None) -> list[str]:
    return (
        []
        if reference is None
        else ["--reference", str(_gap_cli.existing_file(reference, "reference"))]
    )


def whatshap_phase_with_pedigree(
    vcf: str,
    alignments: Sequence[str],
    pedigree: str,
    output: str,
    *,
    reference: str | None = None,
    genetic_map: str | None = None,
    recombination_rate: float | None = None,
    recombination_list: str | None = None,
    no_genetic_haplotyping: bool = False,
    use_ped_samples: bool = False,
    algorithm: str | None = None,
    mapping_quality: int | None = None,
    sample: str | None = None,
    chromosome: str | None = None,
    tag: str | None = None,
    only_snvs: bool = False,
    distrust_genotypes: bool = False,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phase variants with reads and a pedigree (``whatshap phase --ped``).

    ``pedigree`` is a PED/FAM file; ``genetic_map`` (a genetic map file) or ``recombination_rate`` (per base, default
    1.26e-8) model recombination between the family members and ``recombination_list`` writes the recombination
    events. ``no_genetic_haplotyping`` skips the Mendelian haplotyping step and ``use_ped_samples`` phases only the
    samples named in the pedigree. ``algorithm`` is 'whatshap', 'hapchat' or 'heuristic'; ``tag`` stores PS or HP.
    ``vcf`` and ``alignments`` (BAM/CRAM) must contain all family members. Writes a new phased VCF.
    """
    vcf_path = _gap_cli.existing_file(vcf, "vcf")
    ped = _gap_cli.existing_file(pedigree, "pedigree")
    output_path = _gap_cli.new_file(output, "output")
    arguments = ["--ped", str(ped), "-o", str(output_path)] + _reference(reference)
    expected = [output_path]
    if genetic_map is not None:
        arguments += ["--genmap", str(_gap_cli.existing_file(genetic_map, "genetic_map"))]
    if recombination_rate is not None:
        arguments += [
            "--recombrate",
            str(_gap_cli.number(recombination_rate, "recombination_rate", minimum=0.0)),
        ]
    if recombination_list is not None:
        listing = _gap_cli.new_file(recombination_list, "recombination_list")
        arguments += ["--recombination-list", str(listing)]
        expected.append(listing)
    for value, flag_name, name in (
        (no_genetic_haplotyping, "--no-genetic-haplotyping", "no_genetic_haplotyping"),
        (use_ped_samples, "--use-ped-samples", "use_ped_samples"),
    ):
        if _p14_cli.flag(value, name):
            arguments.append(flag_name)
    if algorithm is not None:
        arguments += ["--algorithm", _gap_cli.choice(algorithm, _ALGORITHMS, "algorithm")]
    if mapping_quality is not None:
        arguments += [
            "--mapping-quality",
            str(_gap_cli.integer(mapping_quality, "mapping_quality", minimum=0, maximum=255)),
        ]
    arguments += _common(sample, chromosome, None, tag, only_snvs, distrust_genotypes)
    arguments += [str(vcf_path), *_files(alignments, "alignments")]
    return _whatshap(
        "phase", arguments, cwd=output_path.parent, expected=expected, whatshap_path=whatshap_path
    )


def whatshap_polyphase(
    vcf: str,
    alignments: Sequence[str],
    output: str,
    *,
    ploidy: int,
    reference: str | None = None,
    block_cut_sensitivity: int | None = None,
    min_overlap: int | None = None,
    use_prephasing: bool = False,
    include_haploid_sets: bool = False,
    no_multi_allelic: bool = False,
    mapping_quality: int | None = None,
    threads: int | None = None,
    sample: str | None = None,
    chromosome: str | None = None,
    tag: str | None = None,
    only_snvs: bool = False,
    distrust_genotypes: bool = False,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phase polyploid variants from reads (``whatshap polyphase``, clustering and threading algorithm).

    ``ploidy`` is the number of haplotypes (for example 4). ``block_cut_sensitivity`` (0-5) trades block length
    against accuracy, ``min_overlap`` is the minimum number of shared variants for linking reads, ``use_prephasing``
    keeps existing phase blocks, ``include_haploid_sets`` writes haploid sets, ``no_multi_allelic`` disables
    phasing of multi-allelic variants. Writes a new phased VCF.
    """
    vcf_path = _gap_cli.existing_file(vcf, "vcf")
    output_path = _gap_cli.new_file(output, "output")
    arguments = [
        "-o",
        str(output_path),
        "--ploidy",
        str(_gap_cli.integer(ploidy, "ploidy", minimum=2, maximum=16)),
    ] + _reference(reference)
    if block_cut_sensitivity is not None:
        arguments += [
            "--block-cut-sensitivity",
            str(
                _gap_cli.integer(
                    block_cut_sensitivity, "block_cut_sensitivity", minimum=0, maximum=5
                )
            ),
        ]
    if min_overlap is not None:
        arguments += ["--min-overlap", str(_gap_cli.integer(min_overlap, "min_overlap", minimum=1))]
    for value, flag_name, name in (
        (use_prephasing, "--use-prephasing", "use_prephasing"),
        (include_haploid_sets, "--include-haploid-sets", "include_haploid_sets"),
        (no_multi_allelic, "--no-mav", "no_multi_allelic"),
    ):
        if _p14_cli.flag(value, name):
            arguments.append(flag_name)
    if mapping_quality is not None:
        arguments += [
            "--mapping-quality",
            str(_gap_cli.integer(mapping_quality, "mapping_quality", minimum=0, maximum=255)),
        ]
    if threads is not None:
        arguments += [
            "--threads",
            str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256)),
        ]
    arguments += _common(sample, chromosome, None, tag, only_snvs, distrust_genotypes)
    arguments += [str(vcf_path), *_files(alignments, "alignments")]
    return _whatshap(
        "polyphase",
        arguments,
        cwd=output_path.parent,
        expected=[output_path],
        whatshap_path=whatshap_path,
    )


def whatshap_polyphase_with_progeny(
    vcf: str,
    pedigree: str,
    output: str,
    *,
    ploidy: int,
    progeny_file: str | None = None,
    scoring_window: int | None = None,
    complexity_support: int | None = None,
    sample: str | None = None,
    chromosome: str | None = None,
    tag: str | None = None,
    only_snvs: bool = False,
    distrust_genotypes: bool = False,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phase polyploid parents using the genotypes of their progeny (``whatshap polyphasegenetic``).

    ``pedigree`` is the pedigree file naming the parents and progeny in ``vcf``; ``progeny_file`` optionally lists the
    progeny to use. ``scoring_window`` sets the window of the scoring step and ``complexity_support`` the support
    needed for complex positions. Writes a new phased VCF.
    """
    vcf_path = _gap_cli.existing_file(vcf, "vcf")
    ped = _gap_cli.existing_file(pedigree, "pedigree")
    output_path = _gap_cli.new_file(output, "output")
    arguments = [
        "-o",
        str(output_path),
        "--ploidy",
        str(_gap_cli.integer(ploidy, "ploidy", minimum=2, maximum=16)),
    ]
    if progeny_file is not None:
        arguments += ["-P", str(_gap_cli.existing_file(progeny_file, "progeny_file"))]
    if scoring_window is not None:
        arguments += [
            "--scoring-window",
            str(_gap_cli.integer(scoring_window, "scoring_window", minimum=1)),
        ]
    if complexity_support is not None:
        arguments += [
            "--complexity-support",
            str(_gap_cli.integer(complexity_support, "complexity_support", minimum=1)),
        ]
    arguments += _common(sample, chromosome, None, tag, only_snvs, distrust_genotypes)
    arguments += [str(vcf_path), str(ped)]
    return _whatshap(
        "polyphasegenetic",
        arguments,
        cwd=output_path.parent,
        expected=[output_path],
        whatshap_path=whatshap_path,
    )


def whatshap_phasing_statistics(
    vcf: str,
    tsv: str,
    *,
    gtf: str | None = None,
    block_list: str | None = None,
    chromosome_lengths: str | None = None,
    sample: str | None = None,
    chromosome: str | None = None,
    only_snvs: bool = False,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute phasing statistics of a phased VCF (``whatshap stats``).

    ``tsv`` receives the statistics (phased variants, blocks, block N50 and so on, one row per chromosome plus an
    ``ALL`` row); ``gtf`` writes the blocks as GTF, ``block_list`` one line per block, and ``chromosome_lengths``
    overrides the lengths in the VCF header. ``statistics`` holds the parsed TSV rows. Writes new files.
    """
    vcf_path = _gap_cli.existing_file(vcf, "vcf")
    tsv_path = _gap_cli.new_file(tsv, "tsv")
    arguments = ["--tsv", str(tsv_path)]
    expected = [tsv_path]
    for value, flag_name, name in (
        (gtf, "--gtf", "gtf"),
        (block_list, "--block-list", "block_list"),
    ):
        if value is not None:
            target = _gap_cli.new_file(value, name)
            arguments += [flag_name, str(target)]
            expected.append(target)
    if chromosome_lengths is not None:
        arguments += [
            "--chr-lengths",
            str(_gap_cli.existing_file(chromosome_lengths, "chromosome_lengths")),
        ]
    arguments += _common(sample, chromosome, None, None, only_snvs, False) + [str(vcf_path)]
    result = _whatshap(
        "stats", arguments, cwd=tsv_path.parent, expected=expected, whatshap_path=whatshap_path
    )
    if result["status"] == "completed":
        lines = tsv_path.read_text(encoding="utf-8").lstrip("#").splitlines()
        if len(lines) >= 2:
            names = lines[0].split("\t")
            result["statistics"] = [
                dict(zip(names, row.split("\t"), strict=False)) for row in lines[1:] if row.strip()
            ]
    return result


def whatshap_compare_phasing(
    vcfs: Sequence[str],
    *,
    tsv_pairwise: str | None = None,
    tsv_multiway: str | None = None,
    switch_error_bed: str | None = None,
    longest_block_tsv: str | None = None,
    names: Sequence[str] | None = None,
    ignore_sample_name: bool = False,
    only_snvs: bool = False,
    sample: str | None = None,
    ploidy: int | None = None,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two or more phased VCF files (``whatshap compare``): switch errors, flips and Hamming distances.

    ``tsv_pairwise`` and ``tsv_multiway`` receive the pairwise and multiway comparison tables, ``switch_error_bed``
    the switch-error positions and ``longest_block_tsv`` the longest-block comparison; ``names`` are labels for the
    files, ``ignore_sample_name`` compares single-sample VCFs with different sample names. At least one output file
    must be requested. ``pairwise`` holds the parsed pairwise table.
    """
    files = _files(vcfs, "vcfs")
    if len(files) < 2:
        raise ConfigurationError("vcfs must list at least two VCF files.")
    if (
        tsv_pairwise is None
        and tsv_multiway is None
        and switch_error_bed is None
        and longest_block_tsv is None
    ):
        raise ConfigurationError(
            "Request at least one output file (tsv_pairwise, tsv_multiway, switch_error_bed, longest_block_tsv)."
        )
    folders = {
        Path(item).expanduser().resolve().parent
        for item in (tsv_pairwise, tsv_multiway, switch_error_bed, longest_block_tsv)
        if item is not None
    }
    if len(folders) != 1:
        raise ConfigurationError("All requested output files must be in the same folder.")
    arguments: list[str] = []
    expected: list[Path] = []
    for value, flag_name, name in (
        (tsv_pairwise, "--tsv-pairwise", "tsv_pairwise"),
        (tsv_multiway, "--tsv-multiway", "tsv_multiway"),
        (switch_error_bed, "--switch-error-bed", "switch_error_bed"),
        (longest_block_tsv, "--longest-block-tsv", "longest_block_tsv"),
    ):
        if value is not None:
            target = _gap_cli.new_file(value, name)
            arguments += [flag_name, str(target)]
            expected.append(target)
    if names is not None:
        arguments += ["--names", ",".join(_gap_cli.text(item, "names") for item in names)]
    if _p14_cli.flag(ignore_sample_name, "ignore_sample_name"):
        arguments.append("--ignore-sample-name")
    if ploidy is not None:
        arguments += ["--ploidy", str(_gap_cli.integer(ploidy, "ploidy", minimum=2, maximum=16))]
    arguments += _common(sample, None, None, None, only_snvs, False) + files
    result = _whatshap(
        "compare",
        arguments,
        cwd=next(iter(folders)),
        expected=expected,
        whatshap_path=whatshap_path,
    )
    if result["status"] == "completed" and tsv_pairwise is not None:
        lines = (
            Path(tsv_pairwise)
            .expanduser()
            .resolve()
            .read_text(encoding="utf-8")
            .lstrip("#")
            .splitlines()
        )
        if len(lines) >= 2:
            names_row = lines[0].split("\t")
            result["pairwise"] = [
                dict(zip(names_row, row.split("\t"), strict=False))
                for row in lines[1:]
                if row.strip()
            ]
    return result


def whatshap_unphase_vcf(
    vcf: str,
    output: str,
    *,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove phasing information from a VCF (``whatshap unphase``): phased genotypes ``1|0`` become ``0/1`` and the
    PS, HP and related fields are dropped. Writes a new VCF."""
    vcf_path = _gap_cli.existing_file(vcf, "vcf")
    output_path = _gap_cli.new_file(output, "output")
    return _whatshap(
        "unphase",
        [str(vcf_path)],
        cwd=output_path.parent,
        expected=[],
        whatshap_path=whatshap_path,
        stdout_file=output_path,
    )


def whatshap_haplotag_reads(
    vcf: str,
    alignments: str,
    output: str,
    *,
    reference: str | None = None,
    haplotag_list: str | None = None,
    regions: Sequence[str] | None = None,
    ploidy: int | None = None,
    sample: str | None = None,
    ignore_read_groups: bool = False,
    ignore_linked_reads: bool = False,
    linked_read_distance_cutoff: int | None = None,
    skip_missing_contigs: bool = False,
    tag_supplementary: str | None = None,
    threads: int | None = None,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tag each read with its haplotype using a phased VCF (``whatshap haplotag``): adds the HP (haplotype) and PS
    (phase set) tags to a BAM/CRAM. ``vcf`` must be bgzip-compressed and tabix-indexed (``.vcf.gz`` + ``.tbi``).

    ``output`` is the tagged alignment file (``.bam``/``.cram``); ``haplotag_list`` also writes a TSV of read name,
    haplotype, phase set and chromosome (the input of whatshap_split_reads). ``regions`` limit the tagging,
    ``ploidy`` is the number of haplotypes (default 2), ``tag_supplementary`` handles supplementary alignments
    ('skip', 'copy-primary', 'independent-or-skip', 'independent-or-copy-primary'), and ``threads`` the output
    compression threads. Writes new files.
    """
    vcf_path = _gap_cli.existing_file(vcf, "vcf")
    alignment = _gap_cli.existing_file(alignments, "alignments")
    output_path = _gap_cli.new_file(output, "output")
    arguments = ["-o", str(output_path)] + _reference(reference)
    expected = [output_path]
    if haplotag_list is not None:
        listing = _gap_cli.new_file(haplotag_list, "haplotag_list")
        arguments += ["--output-haplotag-list", str(listing)]
        expected.append(listing)
    for item in regions or ():
        arguments += ["--regions", _gap_cli.text(item, "regions")]
    if ploidy is not None:
        arguments += ["--ploidy", str(_gap_cli.integer(ploidy, "ploidy", minimum=1, maximum=16))]
    if sample is not None:
        arguments += ["--sample", _gap_cli.text(sample, "sample")]
    for value, flag_name, name in (
        (ignore_read_groups, "--ignore-read-groups", "ignore_read_groups"),
        (ignore_linked_reads, "--ignore-linked-read", "ignore_linked_reads"),
        (skip_missing_contigs, "--skip-missing-contigs", "skip_missing_contigs"),
    ):
        if _p14_cli.flag(value, name):
            arguments.append(flag_name)
    if linked_read_distance_cutoff is not None:
        arguments += [
            "--linked-read-distance-cutoff",
            str(
                _gap_cli.integer(
                    linked_read_distance_cutoff, "linked_read_distance_cutoff", minimum=1
                )
            ),
        ]
    if tag_supplementary is not None:
        arguments += [
            "--tag-supplementary",
            _gap_cli.choice(tag_supplementary, _SUPPLEMENTARY, "tag_supplementary"),
        ]
    if threads is not None:
        arguments += [
            "--output-threads",
            str(_gap_cli.integer(threads, "threads", minimum=1, maximum=256)),
        ]
    arguments += [str(vcf_path), str(alignment)]
    return _whatshap(
        "haplotag",
        arguments,
        cwd=output_path.parent,
        expected=expected,
        whatshap_path=whatshap_path,
    )


def whatshap_genotype_variants(
    vcf: str,
    alignments: Sequence[str],
    output: str,
    *,
    reference: str | None = None,
    pedigree: str | None = None,
    genetic_map: str | None = None,
    recombination_rate: float | None = None,
    use_ped_samples: bool = False,
    max_coverage: int | None = None,
    mapping_quality: int | None = None,
    genotype_quality_threshold: float | None = None,
    no_priors: bool = False,
    prior_output: str | None = None,
    sample: str | None = None,
    chromosome: str | None = None,
    only_snvs: bool = False,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genotype variants from reads (``whatshap genotype``), or from reads plus a pedigree (``--ped``).

    Without ``pedigree`` the genotypes of each sample are computed from its reads; with ``pedigree`` (PED/FAM) the
    family is genotyped jointly, with ``genetic_map`` or ``recombination_rate`` modelling recombination and
    ``use_ped_samples`` limiting the samples. ``reference`` enables re-alignment based allele detection;
    ``max_coverage`` down-samples deep regions; ``genotype_quality_threshold`` masks low-quality calls as './.';
    ``no_priors`` uses uniform priors and ``prior_output`` writes the prior genotypes. Writes a new VCF.
    """
    vcf_path = _gap_cli.existing_file(vcf, "vcf")
    output_path = _gap_cli.new_file(output, "output")
    arguments = ["-o", str(output_path)] + _reference(reference)
    expected = [output_path]
    if pedigree is not None:
        arguments += ["--ped", str(_gap_cli.existing_file(pedigree, "pedigree"))]
    elif genetic_map is not None or recombination_rate is not None or use_ped_samples:
        raise ConfigurationError(
            "genetic_map, recombination_rate and use_ped_samples need pedigree."
        )
    if genetic_map is not None:
        arguments += ["--genmap", str(_gap_cli.existing_file(genetic_map, "genetic_map"))]
    if recombination_rate is not None:
        arguments += [
            "--recombrate",
            str(_gap_cli.number(recombination_rate, "recombination_rate", minimum=0.0)),
        ]
    if _p14_cli.flag(use_ped_samples, "use_ped_samples"):
        arguments.append("--use-ped-samples")
    if max_coverage is not None:
        arguments += [
            "--max-coverage",
            str(_gap_cli.integer(max_coverage, "max_coverage", minimum=1)),
        ]
    if mapping_quality is not None:
        arguments += [
            "--mapping-quality",
            str(_gap_cli.integer(mapping_quality, "mapping_quality", minimum=0, maximum=255)),
        ]
    if genotype_quality_threshold is not None:
        arguments += [
            "--gt-qual-threshold",
            str(
                _gap_cli.number(
                    genotype_quality_threshold, "genotype_quality_threshold", minimum=0.0
                )
            ),
        ]
    if _p14_cli.flag(no_priors, "no_priors"):
        arguments.append("--no-priors")
    if prior_output is not None:
        prior = _gap_cli.new_file(prior_output, "prior_output")
        arguments += ["--prioroutput", str(prior)]
        expected.append(prior)
    arguments += _common(sample, chromosome, None, None, only_snvs, False)
    arguments += [str(vcf_path), *_files(alignments, "alignments")]
    return _whatshap(
        "genotype",
        arguments,
        cwd=output_path.parent,
        expected=expected,
        whatshap_path=whatshap_path,
    )


def whatshap_find_snv_candidates(
    reference: str,
    alignment: str,
    output: str,
    *,
    technology: str | None = None,
    min_alt_depth: int | None = None,
    min_alt_fraction: float | None = None,
    multi_allelic: bool = False,
    sample: str | None = None,
    chromosome: str | None = None,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find candidate SNV sites from read alignments (``whatshap find_snv_candidates``).

    A site is a candidate when the alternative allele has at least ``min_alt_depth`` reads (default 3) and
    ``min_alt_fraction`` of the depth. ``technology`` ('illumina', 'pacbio' or 'nanopore') sets the technology
    defaults (fraction 0.25, 0.25 and 0.4); ``multi_allelic`` also writes multi-allelic sites. Writes a new VCF.
    """
    ref = _gap_cli.existing_file(reference, "reference")
    bam = _gap_cli.existing_file(alignment, "alignment")
    output_path = _gap_cli.new_file(output, "output")
    arguments = ["-o", str(output_path)]
    if technology is not None:
        arguments.append("--" + _gap_cli.choice(technology, _TECHNOLOGIES, "technology"))
    if min_alt_depth is not None:
        arguments += ["--minabs", str(_gap_cli.integer(min_alt_depth, "min_alt_depth", minimum=1))]
    if min_alt_fraction is not None:
        arguments += [
            "--minrel",
            str(_gap_cli.number(min_alt_fraction, "min_alt_fraction", minimum=0.0, maximum=1.0)),
        ]
    if _p14_cli.flag(multi_allelic, "multi_allelic"):
        arguments.append("--multi-allelics")
    arguments += _common(sample, chromosome, None, None, False, False)
    arguments += [str(ref), str(bam)]
    return _whatshap(
        "find_snv_candidates",
        arguments,
        cwd=output_path.parent,
        expected=[output_path],
        whatshap_path=whatshap_path,
    )


def whatshap_split_reads(
    reads: str,
    haplotype_list: str,
    *,
    output_h1: str | None = None,
    output_h2: str | None = None,
    output_untagged: str | None = None,
    add_untagged: bool = False,
    only_largest_block: bool = False,
    discard_unknown_reads: bool = False,
    read_lengths_histogram: str | None = None,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the reads of each haplotype to its own file (``whatshap split``).

    ``reads`` is FASTQ (optionally gzip) or BAM; ``haplotype_list`` is a TSV of read name and haplotype (H1, H2, ...),
    for example the ``haplotag_list`` of whatshap_haplotag_reads. ``output_h1``, ``output_h2`` and
    ``output_untagged`` are the output files (``.gz`` supported); ``add_untagged`` copies untagged reads to every
    haplotype file, ``only_largest_block`` keeps tags only from the largest phase block (needs the 4-column list),
    ``discard_unknown_reads`` drops reads that are not in the list. At least one output is required.
    """
    source = _gap_cli.existing_file(reads, "reads")
    listing = _gap_cli.existing_file(haplotype_list, "haplotype_list")
    targets = [
        (output_h1, "--output-h1", "output_h1"),
        (output_h2, "--output-h2", "output_h2"),
        (output_untagged, "--output-untagged", "output_untagged"),
        (read_lengths_histogram, "--read-lengths-histogram", "read_lengths_histogram"),
    ]
    if all(item[0] is None for item in targets):
        raise ConfigurationError("Give at least one output file.")
    arguments: list[str] = []
    expected: list[Path] = []
    for value, flag_name, name in targets:
        if value is not None:
            target = _gap_cli.new_file(value, name)
            arguments += [flag_name, str(target)]
            expected.append(target)
    for value, flag_name, name in (
        (add_untagged, "--add-untagged", "add_untagged"),
        (only_largest_block, "--only-largest-block", "only_largest_block"),
        (discard_unknown_reads, "--discard-unknown-reads", "discard_unknown_reads"),
    ):
        if _p14_cli.flag(value, name):
            arguments.append(flag_name)
    arguments += [str(source), str(listing)]
    folders = {path.parent for path in expected}
    if len(folders) != 1:
        raise ConfigurationError("All output files must be in the same folder.")
    return _whatshap(
        "split", arguments, cwd=next(iter(folders)), expected=expected, whatshap_path=whatshap_path
    )


def whatshap_learn_error_profile(
    alignment: str,
    vcf: str,
    reference: str,
    output: str,
    *,
    kmer_size: int | None = None,
    window: int | None = None,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Learn a sequencing-technology specific error profile from alignments and variants (``whatshap learn``).

    The output counts k-mer pairs (read versus reference) around the variants in ``vcf``; ``kmer_size`` is k (default
    4) and ``window`` the half-width of the region around each variant. Writes a new file.
    """
    bam = _gap_cli.existing_file(alignment, "alignment")
    vcf_path = _gap_cli.existing_file(vcf, "vcf")
    ref = _gap_cli.existing_file(reference, "reference")
    output_path = _gap_cli.new_file(output, "output")
    arguments = ["--reference", str(ref), "--output", str(output_path)]
    if kmer_size is not None:
        arguments += ["-k", str(_gap_cli.integer(kmer_size, "kmer_size", minimum=1, maximum=10))]
    if window is not None:
        arguments += ["--window", str(_gap_cli.integer(window, "window", minimum=1))]
    arguments += [str(bam), str(vcf_path)]
    return _whatshap(
        "learn",
        arguments,
        cwd=output_path.parent,
        expected=[output_path],
        whatshap_path=whatshap_path,
    )


def whatshap_phase_from_haplotags(
    vcf: str,
    alignments: str,
    output: str,
    *,
    reference: str | None = None,
    gap_threshold: float | None = None,
    cut_poly: int | None = None,
    only_indels: bool = False,
    no_multi_allelic: bool = False,
    ignore_read_groups: bool = False,
    sample: str | None = None,
    chromosome: str | None = None,
    whatshap_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add phase information to a VCF from haplotagged reads (``whatshap haplotagphase``).

    ``alignments`` must carry HP/PS tags (see whatshap_haplotag_reads) and ``reference`` is required by WhatsHap; variants that the tagged reads support are
    written phased. ``gap_threshold`` is the percentage of reads with a gap that blocks phasing of an indel,
    ``cut_poly`` shortens homopolymers to the given length when comparing, and ``only_indels`` phases only indels.
    Writes a new phased VCF.
    """
    vcf_path = _gap_cli.existing_file(vcf, "vcf")
    alignment = _gap_cli.existing_file(alignments, "alignments")
    output_path = _gap_cli.new_file(output, "output")
    arguments = ["-o", str(output_path)] + _reference(reference)
    if gap_threshold is not None:
        arguments += [
            "--gap-threshold",
            str(_gap_cli.number(gap_threshold, "gap_threshold", minimum=0.0, maximum=100.0)),
        ]
    if cut_poly is not None:
        arguments += ["--cut-poly", str(_gap_cli.integer(cut_poly, "cut_poly", minimum=1))]
    for value, flag_name, name in (
        (only_indels, "--only-indels", "only_indels"),
        (no_multi_allelic, "--no-mav", "no_multi_allelic"),
        (ignore_read_groups, "--ignore-read-groups", "ignore_read_groups"),
    ):
        if _p14_cli.flag(value, name):
            arguments.append(flag_name)
    arguments += _common(sample, chromosome, None, None, False, False)
    arguments += [str(vcf_path), str(alignment)]
    return _whatshap(
        "haplotagphase",
        arguments,
        cwd=output_path.parent,
        expected=[output_path],
        whatshap_path=whatshap_path,
    )


__all__ = [
    "whatshap_compare_phasing",
    "whatshap_find_snv_candidates",
    "whatshap_genotype_variants",
    "whatshap_haplotag_reads",
    "whatshap_learn_error_profile",
    "whatshap_phase_from_haplotags",
    "whatshap_phase_with_pedigree",
    "whatshap_phasing_statistics",
    "whatshap_polyphase",
    "whatshap_polyphase_with_progeny",
    "whatshap_split_reads",
    "whatshap_unphase_vcf",
]
_register_public_verbose(globals())
