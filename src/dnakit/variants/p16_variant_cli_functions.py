"""Depth, structural-variant and BAM-annotation facades (pack P16): mosdepth, Sniffles2, cuteSV, vcfanno.

``compute_mosdepth_depth`` (per-base depth, global distribution, summary), ``compute_mosdepth_regions`` (window or BED mean
or median depth, region distributions, threshold counts), ``compute_mosdepth_quantized_coverage`` (quantised segments),
``call_structural_variants_sniffles2`` (mosaic, phasing, SNF output), ``merge_structural_variants_sniffles2_snf``
(multi-sample calling from SNF files), ``genotype_structural_variants_cutesv`` and ``annotate_vcf_from_bam_vcfanno``. Programs
run without a shell through ``dnakit.bioinformatics._p16_cli``; results are JSON-serialisable dictionaries.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p16_cli as cli
from dnakit.exceptions import ConfigurationError

_MOSDEPTH = cli.Tool(
    "mosdepth",
    "mosdepth",
    "Install mosdepth (conda install -c bioconda mosdepth) and put it on PATH or pass mosdepth_path.",
    "0.3.14",
    "https://doi.org/10.1093/bioinformatics/btx699",
    ("--version",),
    None,
)
_SNIFFLES = cli.Tool(
    "sniffles2",
    "sniffles",
    "Install Sniffles2 (pip install sniffles or conda install -c bioconda sniffles) and put it on PATH or pass sniffles_path.",
    "2.6",
    "https://doi.org/10.1038/s41592-018-0001-7",
    ("--version",),
    None,
)
_CUTESV = cli.Tool(
    "cutesv",
    "cuteSV",
    "Install cuteSV (pip install cuteSV or conda install -c bioconda cutesv) and put it on PATH or pass cutesv_path.",
    "2.1",
    "https://doi.org/10.1186/s13059-020-01990-9",
    ("--version",),
    None,
)
_VCFANNO = cli.Tool(
    "vcfanno",
    "vcfanno",
    "Install vcfanno (conda install -c bioconda vcfanno) and put it on PATH or pass vcfanno_path.",
    "0.3.7",
    "https://doi.org/10.1186/s13059-016-0973-5",
    ("-version",),
    r"(\d+(?:\.\d+)+)",
)


def _mosdepth_common(
    alignment: str,
    prefix: str,
    *,
    chromosome: str | None,
    threads: int,
    mapq: int,
    min_fragment: int | None,
    max_fragment: int | None,
    fast_mode: bool,
    fragment_mode: bool,
    exclude_flag: int | None,
    include_flag: int | None,
    reference: str | None,
    read_groups: str | None,
) -> tuple[Path, Path, list[str]]:
    bam = cli.existing_file(alignment, "alignment")
    out_prefix = Path(prefix).expanduser().resolve()
    if not out_prefix.parent.is_dir():
        raise ConfigurationError(f"The folder of prefix does not exist: {out_prefix.parent}")
    if any(out_prefix.parent.glob(out_prefix.name + ".*")):
        raise ConfigurationError(
            f"Files with the prefix {out_prefix} already exist and will not be overwritten."
        )
    arguments = [
        "--threads",
        str(cli.integer(threads, "threads", minimum=0, maximum=256)),
        "--mapq",
        str(cli.integer(mapq, "mapq", minimum=0, maximum=255)),
    ]
    if chromosome is not None:
        arguments += ["--chrom", cli.text(chromosome, "chromosome")]
    if min_fragment is not None:
        arguments += ["--min-frag-len", str(cli.integer(min_fragment, "min_fragment", minimum=0))]
    if max_fragment is not None:
        arguments += ["--max-frag-len", str(cli.integer(max_fragment, "max_fragment", minimum=0))]
    if cli.flag(fast_mode, "fast_mode"):
        arguments.append("--fast-mode")
    if cli.flag(fragment_mode, "fragment_mode"):
        arguments.append("--fragment-mode")
    if exclude_flag is not None:
        arguments += ["--flag", str(cli.integer(exclude_flag, "exclude_flag", minimum=0))]
    if include_flag is not None:
        arguments += ["--include-flag", str(cli.integer(include_flag, "include_flag", minimum=0))]
    if reference is not None:
        arguments += ["--fasta", str(cli.existing_file(reference, "reference"))]
    if read_groups is not None:
        arguments += ["--read-groups", cli.text(read_groups, "read_groups")]
    return bam, out_prefix, arguments


def _summary_rows(prefix: Path) -> list[dict[str, Any]]:
    path = Path(str(prefix) + ".mosdepth.summary.txt")
    return cli.read_table(path, text_columns=frozenset({"chrom"})) if path.is_file() else []


def compute_mosdepth_depth(
    alignment: str,
    prefix: str,
    *,
    per_base: bool = True,
    chromosome: str | None = None,
    threads: int = 0,
    mapq: int = 0,
    min_fragment: int | None = None,
    max_fragment: int | None = None,
    fast_mode: bool = False,
    fragment_mode: bool = False,
    exclude_flag: int | None = None,
    include_flag: int | None = None,
    reference: str | None = None,
    read_groups: str | None = None,
    mosdepth_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-base depth, genome-wide depth distribution and per-chromosome summary of a BAM/CRAM (``mosdepth``).

    Writes ``<prefix>.per-base.bed.gz`` (unless ``per_base=False``, ``--no-per-base``), ``<prefix>.mosdepth.global.dist.txt``
    (cumulative fraction of bases at each depth) and ``<prefix>.mosdepth.summary.txt`` (length, bases, mean, min, max per
    chromosome and total, parsed into ``summary``). The remaining parameters are mosdepth's: ``--chrom``, ``--mapq``,
    ``--min-frag-len``, ``--max-frag-len``, ``--fast-mode``, ``--fragment-mode``, ``--flag`` (default excludes 1796),
    ``--include-flag``, ``--fasta`` (CRAM) and ``--read-groups``.
    """
    bam, out_prefix, arguments = _mosdepth_common(
        alignment,
        prefix,
        chromosome=chromosome,
        threads=threads,
        mapq=mapq,
        min_fragment=min_fragment,
        max_fragment=max_fragment,
        fast_mode=fast_mode,
        fragment_mode=fragment_mode,
        exclude_flag=exclude_flag,
        include_flag=include_flag,
        reference=reference,
        read_groups=read_groups,
    )
    expected = [
        Path(str(out_prefix) + ".mosdepth.global.dist.txt"),
        Path(str(out_prefix) + ".mosdepth.summary.txt"),
    ]
    if cli.flag(per_base, "per_base"):
        expected.append(Path(str(out_prefix) + ".per-base.bed.gz"))
    else:
        arguments.append("--no-per-base")
    record = cli.invoke(
        _MOSDEPTH,
        [*arguments, str(out_prefix), str(bam)],
        cwd=out_prefix.parent,
        expected=expected,
        explicit=mosdepth_path,
    )
    distribution = cli.read_table(
        expected[0],
        header=False,
        names=["chrom", "depth", "fraction"],
        text_columns=frozenset({"chrom"}),
        limit=200,
    )
    return cli.finish(
        record,
        "compute_mosdepth_depth",
        summary=_summary_rows(out_prefix),
        distribution_head=distribution,
    )


def compute_mosdepth_regions(
    alignment: str,
    prefix: str,
    *,
    by: str,
    thresholds: Sequence[int] | None = None,
    use_median: bool = False,
    per_base: bool = False,
    chromosome: str | None = None,
    threads: int = 0,
    mapq: int = 0,
    fast_mode: bool = False,
    fragment_mode: bool = False,
    exclude_flag: int | None = None,
    include_flag: int | None = None,
    reference: str | None = None,
    read_groups: str | None = None,
    mosdepth_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mean (or median) depth, depth histogram and threshold counts per window or BED region (``mosdepth --by``).

    ``by`` is an integer window size (for example ``"1000"``) or the path of a BED file. Outputs are
    ``<prefix>.regions.bed.gz`` (mean depth per region; median with ``use_median``), ``<prefix>.mosdepth.region.dist.txt``
    (cumulative depth distribution per chromosome over the regions) and, with ``thresholds`` (``--thresholds t1,t2``),
    ``<prefix>.thresholds.bed.gz`` with the number of bases covered at least that deeply. Per-base output is off by
    default (``--no-per-base``); ``per_base=True`` also writes it. Other parameters as in :func:`compute_mosdepth_depth`.
    """
    bam, out_prefix, arguments = _mosdepth_common(
        alignment,
        prefix,
        chromosome=chromosome,
        threads=threads,
        mapq=mapq,
        min_fragment=None,
        max_fragment=None,
        fast_mode=fast_mode,
        fragment_mode=fragment_mode,
        exclude_flag=exclude_flag,
        include_flag=include_flag,
        reference=reference,
        read_groups=read_groups,
    )
    window = str(by)
    if window.isdigit():
        by_value = str(cli.integer(int(window), "by", minimum=1))
    else:
        by_value = str(cli.existing_file(window, "by"))
    arguments += ["--by", by_value]
    expected = [
        Path(str(out_prefix) + ".mosdepth.global.dist.txt"),
        Path(str(out_prefix) + ".mosdepth.summary.txt"),
        Path(str(out_prefix) + ".mosdepth.region.dist.txt"),
        Path(str(out_prefix) + ".regions.bed.gz"),
    ]
    if cli.flag(use_median, "use_median"):
        arguments.append("--use-median")
    if thresholds:
        values = [cli.integer(t, "thresholds", minimum=1) for t in thresholds]
        arguments += ["--thresholds", ",".join(map(str, values))]
        expected.append(Path(str(out_prefix) + ".thresholds.bed.gz"))
    if cli.flag(per_base, "per_base"):
        expected.append(Path(str(out_prefix) + ".per-base.bed.gz"))
    else:
        arguments.append("--no-per-base")
    record = cli.invoke(
        _MOSDEPTH,
        [*arguments, str(out_prefix), str(bam)],
        cwd=out_prefix.parent,
        expected=expected,
        explicit=mosdepth_path,
    )
    regions = cli.read_table(
        expected[3],
        header=False,
        names=["chrom", "start", "end", "depth"],
        text_columns=frozenset({"chrom"}),
        limit=1000,
    )
    return cli.finish(
        record, "compute_mosdepth_regions", summary=_summary_rows(out_prefix), regions=regions
    )


def compute_mosdepth_quantized_coverage(
    alignment: str,
    prefix: str,
    *,
    quantize: str,
    chromosome: str | None = None,
    threads: int = 0,
    mapq: int = 0,
    fast_mode: bool = False,
    reference: str | None = None,
    mosdepth_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collapse per-base depth into labelled coverage segments (``mosdepth --quantize``).

    ``quantize`` is mosdepth's bin specification, for example ``"0:1:4:"`` (bins 0, 1-3, 4+) or ``"0:1:5:150"``, and
    writes ``<prefix>.quantized.bed.gz`` with ``chrom, start, end, label`` rows (parsed into ``segments``, first 1,000).
    """
    bam, out_prefix, arguments = _mosdepth_common(
        alignment,
        prefix,
        chromosome=chromosome,
        threads=threads,
        mapq=mapq,
        min_fragment=None,
        max_fragment=None,
        fast_mode=fast_mode,
        fragment_mode=False,
        exclude_flag=None,
        include_flag=None,
        reference=reference,
        read_groups=None,
    )
    spec = cli.text(quantize, "quantize")
    if not all(part == "" or part.isdigit() for part in spec.split(":")):
        raise ConfigurationError("quantize must be colon-separated integers such as '0:1:4:'.")
    expected = [
        Path(str(out_prefix) + ".quantized.bed.gz"),
        Path(str(out_prefix) + ".mosdepth.summary.txt"),
        Path(str(out_prefix) + ".mosdepth.global.dist.txt"),
    ]
    # the program documents --quantize as a value with '=' syntax
    record = cli.invoke(
        _MOSDEPTH,
        [*arguments, "--no-per-base", f"--quantize={spec}", str(out_prefix), str(bam)],
        cwd=out_prefix.parent,
        expected=expected,
        explicit=mosdepth_path,
    )
    segments = cli.read_table(
        expected[0],
        header=False,
        names=["chrom", "start", "end", "label"],
        text_columns=frozenset({"chrom", "label"}),
        limit=1000,
    )
    return cli.finish(record, "compute_mosdepth_quantized_coverage", segments=segments)


def call_structural_variants_sniffles2(
    alignment: str,
    vcf_path: str,
    *,
    reference: str | None = None,
    snf_path: str | None = None,
    mosaic: bool = False,
    phase: bool = False,
    min_sv_length: int | None = None,
    min_support: int | str | None = None,
    sample_id: str | None = None,
    threads: int = 1,
    sniffles_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call SVs from an aligned long-read BAM/CRAM with Sniffles2, optionally mosaic, phased, and saving an SNF.

    ``mosaic=True`` (``--mosaic``) detects low-allele-fraction/somatic SVs; ``phase=True`` (``--phase``) phases calls with
    the reads' haplotype tags (HP/PS; needs a phased BAM); ``snf_path`` (``--snf``) saves the per-sample candidates in a
    new SNF file for later multi-sample calling (:func:`merge_structural_variants_sniffles2_snf`). ``min_support`` is
    ``--minsupport`` (integer or ``auto``). The VCF is written to ``vcf_path``.
    """
    bam = cli.existing_file(alignment, "alignment")
    vcf = cli.new_file(vcf_path, "vcf_path")
    arguments = [
        "--input",
        str(bam),
        "--vcf",
        str(vcf),
        "--threads",
        str(cli.threads_value(threads)),
    ]
    expected = [vcf]
    if reference is not None:
        arguments += ["--reference", str(cli.existing_file(reference, "reference"))]
    if snf_path is not None:
        snf = cli.new_file(snf_path, "snf_path")
        arguments += ["--snf", str(snf)]
        expected.append(snf)
    if cli.flag(mosaic, "mosaic"):
        arguments.append("--mosaic")
    if cli.flag(phase, "phase"):
        arguments.append("--phase")
    if min_sv_length is not None:
        arguments += ["--minsvlen", str(cli.integer(min_sv_length, "min_sv_length", minimum=1))]
    if min_support is not None:
        arguments += [
            "--minsupport",
            "auto"
            if min_support == "auto"
            else str(cli.integer(min_support, "min_support", minimum=1)),
        ]
    if sample_id is not None:
        arguments += ["--sample-id", cli.text(sample_id, "sample_id")]
    record = cli.invoke(
        _SNIFFLES, arguments, cwd=vcf.parent, expected=expected, explicit=sniffles_path
    )
    return cli.finish(
        record,
        "call_structural_variants_sniffles2",
        mosaic=mosaic,
        phase=phase,
        records=sum(1 for line in cli.open_text(vcf) if not line.startswith("#")),
    )


def merge_structural_variants_sniffles2_snf(
    snf_files: Sequence[str],
    vcf_path: str,
    *,
    reference: str | None = None,
    min_support: int | str | None = None,
    threads: int = 1,
    sniffles_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Multi-sample SV calling from per-sample SNF candidate files (``sniffles --input a.snf b.snf ...``)."""
    sources = cli.existing_files(list(snf_files), "snf_files", minimum=2)
    vcf = cli.new_file(vcf_path, "vcf_path")
    arguments = [
        "--input",
        *map(str, sources),
        "--vcf",
        str(vcf),
        "--threads",
        str(cli.threads_value(threads)),
    ]
    if reference is not None:
        arguments += ["--reference", str(cli.existing_file(reference, "reference"))]
    if min_support is not None:
        arguments += [
            "--minsupport",
            "auto"
            if min_support == "auto"
            else str(cli.integer(min_support, "min_support", minimum=1)),
        ]
    record = cli.invoke(
        _SNIFFLES, arguments, cwd=vcf.parent, expected=[vcf], explicit=sniffles_path
    )
    return cli.finish(
        record,
        "merge_structural_variants_sniffles2_snf",
        samples=len(sources),
        records=sum(1 for line in cli.open_text(vcf) if not line.startswith("#")),
    )


def genotype_structural_variants_cutesv(
    alignment: str,
    reference: str,
    vcf_path: str,
    work_directory: str,
    *,
    genotype: bool = True,
    force_calling_vcf: str | None = None,
    platform: str | None = None,
    min_support: int | None = None,
    min_size: int | None = None,
    sample: str | None = None,
    threads: int = 1,
    cutesv_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call SVs from long reads and genotype them (``cuteSV --genotype``).

    ``genotype=True`` adds genotypes and allele depths (``--genotype``). ``force_calling_vcf`` (``-Ivcf``) genotypes the SVs of
    that VCF instead of calling de novo. ``platform`` sets the preset: ``ont`` (``--max_cluster_bias_INS 100 ...``),
    ``hifi`` or ``clr`` (Pacific Biosciences) as the program's documented parameter sets. ``work_directory`` must be new.
    """
    bam = cli.existing_file(alignment, "alignment")
    ref = cli.existing_file(reference, "reference")
    vcf = cli.new_file(vcf_path, "vcf_path")
    work = cli.new_directory(work_directory, "work_directory")
    arguments = [
        str(bam),
        str(ref),
        str(vcf),
        str(work),
        "--threads",
        str(cli.threads_value(threads)),
    ]
    if cli.flag(genotype, "genotype"):
        arguments.append("--genotype")
    if force_calling_vcf is not None:
        arguments += ["-Ivcf", str(cli.existing_file(force_calling_vcf, "force_calling_vcf"))]
    presets = {
        "ont": [
            "--max_cluster_bias_INS",
            "100",
            "--diff_ratio_merging_INS",
            "0.3",
            "--max_cluster_bias_DEL",
            "200",
            "--diff_ratio_merging_DEL",
            "0.5",
        ],
        "hifi": [
            "--max_cluster_bias_INS",
            "1000",
            "--diff_ratio_merging_INS",
            "0.9",
            "--max_cluster_bias_DEL",
            "1000",
            "--diff_ratio_merging_DEL",
            "0.5",
        ],
        "clr": [
            "--max_cluster_bias_INS",
            "100",
            "--diff_ratio_merging_INS",
            "0.3",
            "--max_cluster_bias_DEL",
            "200",
            "--diff_ratio_merging_DEL",
            "0.5",
        ],
    }
    if platform is not None:
        arguments += presets[cli.choice(platform, tuple(presets), "platform")]
    if min_support is not None:
        arguments += ["--min_support", str(cli.integer(min_support, "min_support", minimum=1))]
    if min_size is not None:
        arguments += ["--min_size", str(cli.integer(min_size, "min_size", minimum=1))]
    if sample is not None:
        arguments += ["--sample", cli.text(sample, "sample")]
    record = cli.invoke(_CUTESV, arguments, cwd=vcf.parent, expected=[vcf], explicit=cutesv_path)
    return cli.finish(
        record,
        "genotype_structural_variants_cutesv",
        genotype=genotype,
        records=sum(1 for line in cli.open_text(vcf) if not line.startswith("#")),
    )


# BAM annotation in vcfanno: the *field* selects what is read (depth, mapq, seq) and the op reduces it.
_BAM_OPERATIONS = {"depth": "count", "mapq": "mean", "seq": "concat"}


def annotate_vcf_from_bam_vcfanno(
    vcf: str,
    alignment: str,
    output_path: str,
    *,
    operations: Sequence[str] = ("depth", "mapq", "seq"),
    threads: int = 1,
    vcfanno_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate VCF variants with read depth, mapping quality and read bases from a BAM (``vcfanno``).

    A vcfanno configuration is written next to the output (``<output_path>.toml``) with one ``[[annotation]]`` block
    for the BAM and the chosen ``operations`` (vcfanno BAM fields): ``depth`` (read count, INFO ``<bam>_depth``), ``mapq``
    (mean mapping quality, ``_mapq``) and ``seq`` (the read bases covering each variant, concatenated, ``_seq``). The BAM must be indexed. Existing annotate_vcf_vcfanno takes a hand-written
    configuration; this function builds it for the BAM case.
    """
    source = cli.existing_file(vcf, "vcf")
    bam = cli.existing_file(alignment, "alignment")
    target = cli.new_file(output_path, "output_path")
    chosen = list(operations)
    if not chosen or any(op not in _BAM_OPERATIONS for op in chosen):
        raise ConfigurationError(
            f"operations must be a non-empty subset of {sorted(_BAM_OPERATIONS)}."
        )
    config = target.with_name(target.name + ".toml")
    stem = bam.stem
    names = ", ".join(f'"{stem}_{field}"' for field in chosen)
    fields = ", ".join(f'"{field}"' for field in chosen)
    ops = ", ".join(f'"{_BAM_OPERATIONS[field]}"' for field in chosen)
    config.write_text(
        f'[[annotation]]\nfile="{bam}"\nnames=[{names}]\nfields=[{fields}]\nops=[{ops}]\n',
        encoding="utf-8",
    )
    record = cli.invoke(
        _VCFANNO,
        ["-p", str(cli.threads_value(threads)), str(config), str(source)],
        cwd=target.parent,
        stdout_to=target,
        explicit=vcfanno_path,
    )
    return cli.finish(
        record,
        "annotate_vcf_from_bam_vcfanno",
        configuration=str(config),
        records=sum(1 for line in cli.open_text(target) if not line.startswith("#")),
    )


__all__ = [
    "annotate_vcf_from_bam_vcfanno",
    "call_structural_variants_sniffles2",
    "compute_mosdepth_depth",
    "compute_mosdepth_quantized_coverage",
    "compute_mosdepth_regions",
    "genotype_structural_variants_cutesv",
    "merge_structural_variants_sniffles2_snf",
]

_register_public_verbose(globals(), tuple(__all__))
