"""DELLY command-line facades: genotyping, assembly-based SV calling, copy-number calling, merging, filtering and helper scripts.

DELLY (https://github.com/dellytools/delly, tested against 2.7.0) is installed by the user (``conda install -c bioconda delly`` or a
release binary) and put on PATH or passed as ``delly_path``. Each function runs it with a fixed argument whitelist and no shell, in
the folder of a new output file, and returns a JSON-serialisable dictionary (status, backend, version, command, outputs, parsed
records, log tail). Output files must not exist yet. BCF outputs are parsed with pysam when it is installed.

Short-read and long-read SV discovery is ``dnakit.variants.external.call_delly``; this module adds the other DELLY commands.
"""

from __future__ import annotations

import os
import sys
import sysconfig
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p15_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="delly",
    executable="delly",
    hint="Install DELLY 2.x (for example `conda install -c bioconda delly`, or a release binary from https://github.com/dellytools/delly/releases) and put delly on PATH or pass delly_path.",
    tested_version="2.7.0",
    citation="https://doi.org/10.1093/bioinformatics/bts378",
    version_arguments=("--version",),
    version_pattern=r"Version:?\s*v?(\d+\.\d+\.\d+)",
)
_PYTHON_TOOL = _cli.Tool(
    name="delly2bnd",
    executable="python",
    hint="A Python interpreter with cyvcf2 is required to run delly2bnd.py.",
    tested_version="2.7.0",
    citation="https://github.com/dellytools/delly",
    version_arguments=("--version",),
)
_R_TOOL = _cli.Tool(
    name="delly-rd-plot",
    executable="Rscript",
    hint="Install R with the packages DNAcopy, ggplot2, scales and gtable (for example `conda install -c bioconda -c conda-forge r-base r-ggplot2 bioconductor-dnacopy`) and put Rscript on PATH.",
    tested_version="2.7.0",
    citation="https://github.com/dellytools/delly",
    version_arguments=("--version",),
    version_pattern=r"(\d+\.\d+\.\d+)",
)
_SV_TYPES = ("DEL", "INS", "DUP", "INV", "BND", "ALL")


def _variants(path: Path, limit: int = 1000) -> list[dict[str, Any]]:
    """Records of a BCF/VCF as dictionaries (empty list when pysam is not installed)."""
    try:
        import pysam
    except ImportError:
        return []
    rows: list[dict[str, Any]] = []
    with pysam.VariantFile(str(path)) as handle:
        for record in handle:
            if len(rows) >= limit:
                break
            rows.append(
                {
                    "chrom": record.chrom,
                    "start": record.pos,
                    "end": record.stop,
                    "id": record.id,
                    "type": record.info.get("SVTYPE"),
                    "filter": list(record.filter.keys()),
                    "quality": record.qual,
                }
            )
    return rows


def _count(path: Path) -> int | None:
    try:
        import pysam
    except ImportError:
        return None
    with pysam.VariantFile(str(path)) as handle:
        return sum(1 for _ in handle)


def _bams(values: Sequence[str], name: str) -> list[str]:
    return [str(path) for path in _cli.existing_files(values, name)]


def delly_genotype_sv_sites(
    bam_files: Sequence[str],
    reference: str,
    sites_file: str,
    output_bcf: str,
    *,
    mode: str = "sr",
    technology: str = "ont",
    exclude_regions: str | None = None,
    delly_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genotype known structural-variant sites in samples from short or long reads (``delly sr -v`` / ``delly lr -v``).

    ``sites_file`` is a BCF/VCF of SV sites, for example the merged call set of ``delly merge``. ``mode`` is ``'sr'`` (short
    reads) or ``'lr'`` (long reads; ``technology`` ``'ont'`` or ``'pb'``). The coordinate-sorted, indexed ``bam_files`` give one
    sample each. ``exclude_regions`` is DELLY's tab-separated exclude file. The genotyped sites go to the new BCF ``output_bcf``
    (returned as ``variants`` and ``variant_count``).
    """
    bams = _bams(bam_files, "bam_files")
    genome = _cli.existing_file(reference, "reference")
    sites = _cli.existing_file(sites_file, "sites_file")
    target = _cli.new_file(output_bcf, "output_bcf")
    chosen = _cli.choice(mode, ("sr", "lr"), "mode")
    arguments = [chosen, "-g", str(genome), "-v", str(sites), "-o", str(target)]
    if chosen == "lr":
        arguments += ["-y", _cli.choice(technology, ("ont", "pb"), "technology")]
    if exclude_regions is not None:
        arguments += ["-x", str(_cli.existing_file(exclude_regions, "exclude_regions"))]
    record = _cli.invoke(
        _TOOL,
        [*arguments, *bams],
        cwd=target.parent,
        expected=[target],
        explicit=delly_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(
        record,
        "genotype_sv_sites",
        output=str(target),
        mode=chosen,
        variant_count=_count(target),
        variants=_variants(target),
    )


def delly_assembly_sv(
    assembly_bam_files: Sequence[str],
    reference: str,
    output_bcf: str,
    *,
    mode: str = "squashed",
    svtype: str = "ALL",
    min_mapping_quality: int = 10,
    min_clip: int = 25,
    delly_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Discover structural variants from the alignment of genome assemblies to the reference (``delly asm``).

    ``assembly_bam_files`` are assemblies aligned to ``reference`` (coordinate-sorted, indexed BAM). ``mode='squashed'`` takes
    one BAM per assembly; ``'diploid'`` takes the two haplotype BAMs of each assembly one after the other. ``svtype`` is
    ``DEL``, ``INS``, ``DUP``, ``INV``, ``BND`` or ``ALL``; ``min_mapping_quality`` and ``min_clip`` are DELLY's ``-q`` and
    ``-c``. The calls go to the new BCF ``output_bcf`` (returned as ``variants``).
    """
    bams = _bams(assembly_bam_files, "assembly_bam_files")
    genome = _cli.existing_file(reference, "reference")
    target = _cli.new_file(output_bcf, "output_bcf")
    chosen = _cli.choice(mode, ("squashed", "diploid"), "mode")
    arguments = [
        "asm",
        "-g",
        str(genome),
        "-y",
        chosen,
        "-t",
        _cli.choice(svtype, _SV_TYPES, "svtype"),
        "-q",
        str(_cli.integer(min_mapping_quality, "min_mapping_quality", minimum=0)),
        "-c",
        str(_cli.integer(min_clip, "min_clip", minimum=1)),
        "-o",
        str(target),
    ]
    record = _cli.invoke(
        _TOOL,
        [*arguments, *bams],
        cwd=target.parent,
        expected=[target],
        explicit=delly_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(
        record,
        "assembly_based_sv_discovery",
        output=str(target),
        variant_count=_count(target),
        variants=_variants(target),
    )


def delly_copy_number_variants(
    bam_file: str,
    reference: str,
    output_bcf: str,
    *,
    genotype_sites: str | None = None,
    coverage_file: str | None = None,
    segmentation_file: str | None = None,
    exclude_regions: str | None = None,
    min_cnv_size: int = 1000,
    min_mapping_quality: int = 10,
    ploidy: int = 2,
    window: int = 0,
    delly_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call copy-number variants from read depth, or re-genotype known CNV sites (``delly cnv``).

    One coordinate-sorted, indexed ``bam_file`` is analysed. DELLY counts reads in windows (``window`` bp, 0 chooses
    automatically), corrects for GC and segments the profile. With ``genotype_sites`` (a CNV BCF from an earlier run or ``delly merge -e``)
    the sites are re-genotyped instead of discovered (``-v``). ``coverage_file`` (new, ``.gz``) receives the window read depths
    and ``segmentation_file`` (new, BED) the segments; these can feed ``delly_plot_read_depth``. The calls go to the new BCF
    ``output_bcf``.
    """
    bam = _cli.existing_file(bam_file, "bam_file")
    genome = _cli.existing_file(reference, "reference")
    target = _cli.new_file(output_bcf, "output_bcf")
    arguments = [
        "cnv",
        "-g",
        str(genome),
        "-o",
        str(target),
        "-z",
        str(_cli.integer(min_cnv_size, "min_cnv_size", minimum=1)),
        "-q",
        str(_cli.integer(min_mapping_quality, "min_mapping_quality", minimum=0)),
        "-y",
        str(_cli.integer(ploidy, "ploidy", minimum=1)),
        "-w",
        str(_cli.integer(window, "window", minimum=0)),
    ]
    expected = [target]
    if genotype_sites is not None:
        arguments += ["-v", str(_cli.existing_file(genotype_sites, "genotype_sites"))]
    if exclude_regions is not None:
        arguments += ["-x", str(_cli.existing_file(exclude_regions, "exclude_regions"))]
    for flag, value, name in (
        ("-c", coverage_file, "coverage_file"),
        ("-u", segmentation_file, "segmentation_file"),
    ):
        if value is not None:
            path = _cli.new_file(value, name)
            arguments += [flag, str(path)]
            expected.append(path)
    record = _cli.invoke(
        _TOOL,
        [*arguments, str(bam)],
        cwd=target.parent,
        expected=expected,
        explicit=delly_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(
        record,
        "copy_number_variants",
        output=str(target),
        regenotyped=genotype_sites is not None,
        variant_count=_count(target),
        variants=_variants(target),
    )


def delly_merge_call_sets(
    bcf_files: Sequence[str],
    output_bcf: str,
    *,
    cnv_mode: bool = False,
    long_reads: bool = False,
    reference: str | None = None,
    min_quality: int = 200,
    min_size: int = 0,
    max_size: int = 1_000_000,
    pass_only: bool = False,
    delly_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge the SV or CNV call sets of several samples into one site list (``delly merge``, ``delly merge -e`` for CNVs).

    ``cnv_mode=True`` merges the BCF files written by ``delly_copy_number_variants`` (``-e``); ``long_reads=True`` uses the
    long-read merge mode (``-l``). ``min_quality``, ``min_size``, ``max_size`` and ``pass_only`` filter the sites; ``reference``
    is the reference FASTA (needed for some SV types). The merged site list goes to the new BCF ``output_bcf``.
    """
    files = _bams(bcf_files, "bcf_files")
    target = _cli.new_file(output_bcf, "output_bcf")
    arguments = [
        "merge",
        "-o",
        str(target),
        "-y",
        str(_cli.integer(min_quality, "min_quality", minimum=0)),
        "-m",
        str(_cli.integer(min_size, "min_size", minimum=0)),
        "-n",
        str(_cli.integer(max_size, "max_size", minimum=1)),
    ]
    if reference is not None:
        arguments += ["-g", str(_cli.existing_file(reference, "reference"))]
    for enabled, option in ((cnv_mode, "-e"), (long_reads, "-l"), (pass_only, "-p")):
        if _cli.flag(enabled, option):
            arguments.append(option)
    record = _cli.invoke(
        _TOOL,
        [*arguments, *files],
        cwd=target.parent,
        expected=[target],
        explicit=delly_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(
        record,
        "merge_call_sets",
        output=str(target),
        cnv_mode=cnv_mode,
        variant_count=_count(target),
        variants=_variants(target),
    )


def delly_filter_call_set(
    bcf_file: str,
    output_bcf: str,
    *,
    filter_mode: str = "somatic",
    samples_file: str | None = None,
    min_size: int = 0,
    max_size: int = 500_000_000,
    min_quality: int = 300,
    min_alt_fraction: float = 0.03,
    min_genotyped_fraction: float = 0.75,
    min_coverage: int = 10,
    max_control_support: float = 0.0,
    long_reads: bool = False,
    pass_only: bool = False,
    tag_only: bool = False,
    delly_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter a genotyped SV/CNV call set for somatic or germline variants (``delly filter``).

    ``filter_mode='somatic'`` (needs ``samples_file``, a tab-separated list of sample ids with ``tumor`` or ``control``) keeps
    sites with enough tumour coverage (``min_coverage``) and at most ``max_control_support`` ALT fraction in the controls;
    ``'germline'`` applies population filters (``long_reads`` for long-read call sets). ``min_size``, ``max_size``, ``min_quality``,
    ``min_alt_fraction`` and ``min_genotyped_fraction`` are DELLY's ``-m -n -y -a -r``; ``pass_only`` keeps PASS sites and
    ``tag_only`` tags filtered sites instead of removing them. The result goes to the new BCF ``output_bcf``.
    """
    source = _cli.existing_file(bcf_file, "bcf_file")
    target = _cli.new_file(output_bcf, "output_bcf")
    chosen = _cli.choice(filter_mode, ("somatic", "germline"), "filter_mode")
    arguments = [
        "filter",
        "-f",
        chosen,
        "-o",
        str(target),
        "-m",
        str(_cli.integer(min_size, "min_size", minimum=0)),
        "-n",
        str(_cli.integer(max_size, "max_size", minimum=1)),
        "-y",
        str(_cli.integer(min_quality, "min_quality", minimum=0)),
        "-a",
        format(_cli.number(min_alt_fraction, "min_alt_fraction", minimum=0.0, maximum=1.0), ".12g"),
        "-r",
        format(
            _cli.number(min_genotyped_fraction, "min_genotyped_fraction", minimum=0.0, maximum=1.0),
            ".12g",
        ),
    ]
    if chosen == "somatic":
        if samples_file is None:
            raise ConfigurationError("filter_mode='somatic' needs samples_file.")
        arguments += [
            "-s",
            str(_cli.existing_file(samples_file, "samples_file")),
            "-v",
            str(_cli.integer(min_coverage, "min_coverage", minimum=0)),
            "-c",
            format(
                _cli.number(max_control_support, "max_control_support", minimum=0.0, maximum=1.0),
                ".12g",
            ),
        ]
    elif long_reads:
        arguments.append("-l")
    for enabled, option in ((pass_only, "-p"), (tag_only, "-t")):
        if _cli.flag(enabled, option):
            arguments.append(option)
    record = _cli.invoke(
        _TOOL,
        [*arguments, str(source)],
        cwd=target.parent,
        expected=[target],
        explicit=delly_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(
        record,
        "filter_call_set",
        output=str(target),
        filter_mode=chosen,
        variant_count=_count(target),
        variants=_variants(target),
    )


def delly_translocations_to_breakends(
    vcf_file: str,
    reference: str,
    output_vcf: str,
    script_path: str,
    *,
    python_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert DELLY translocation (``TRA``/``BND``) records into pairs of breakend (BND) VCF records (``scripts/delly2bnd.py``).

    ``script_path`` is the ``delly2bnd.py`` file of the DELLY repository (it is not installed by the package managers) and
    needs ``cyvcf2`` in the Python it runs with (``python_path``, default: the running interpreter). ``vcf_file`` is the DELLY
    output, ``reference`` the reference FASTA (breakpoint bases are read from it). For every translocation two records with
    ``MATEID`` and the standard breakend ALT notation are written to the new file ``output_vcf``.
    """
    script = _cli.existing_file(script_path, "script_path")
    source = _cli.existing_file(vcf_file, "vcf_file")
    genome = _cli.existing_file(reference, "reference")
    target = _cli.new_file(output_vcf, "output_vcf")
    interpreter = python_path or sys.executable
    environment = None
    if (
        python_path is None
    ):  # the resolved interpreter path may lose the virtual environment; keep its packages visible
        environment = {
            "PYTHONPATH": os.pathsep.join(
                filter(None, [sysconfig.get_paths()["purelib"], os.environ.get("PYTHONPATH", "")])
            )
        }
    record = _cli.invoke(
        _PYTHON_TOOL,
        [str(script), "-v", str(source), "-r", str(genome), "-o", str(target)],
        cwd=target.parent,
        expected=[target],
        explicit=interpreter,
        executable_name=Path(interpreter).name,
        timeout_seconds=timeout_seconds,
        environment=environment,
    )
    breakends = 0
    with target.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("#"):
                breakends += 1
    return _cli.finish(
        record, "translocations_to_breakends", output=str(target), breakend_record_count=breakends
    )


def delly_plot_read_depth(
    coverage_file: str,
    output_directory: str,
    script_path: str,
    *,
    segmentation_file: str | None = None,
    rscript_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the read-depth copy-number profile and its segmentation (``R/rd.R`` of DELLY).

    ``coverage_file`` is the window coverage written by ``delly_copy_number_variants`` (``coverage_file``); ``script_path`` is
    DELLY's ``rd.R`` and R needs the packages DNAcopy, ggplot2, scales and gtable. Without ``segmentation_file`` the profile
    is segmented by circular binary segmentation (DNAcopy). The plots ``plot.wholegenome.png`` and ``plot.<chromosome>.png``
    (and ``segmentation.tsv`` when it segmented) are written to the new folder ``output_directory``; ``images`` lists the
    PNG paths.
    """
    script = _cli.existing_file(script_path, "script_path")
    coverage = _cli.existing_file(coverage_file, "coverage_file")
    out = _cli.new_directory(output_directory, "output_directory")
    arguments = [str(script), str(coverage)]
    if segmentation_file is not None:
        arguments.append(str(_cli.existing_file(segmentation_file, "segmentation_file")))
    record = _cli.invoke(
        _R_TOOL,
        arguments,
        cwd=out,
        expected=[out / "plot.wholegenome.png"],
        explicit=rscript_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(
        record,
        "plot_read_depth_profile",
        output_directory=str(out),
        images=sorted(str(p) for p in out.glob("plot.*.png")),
        segmentation=str(out / "segmentation.tsv")
        if (out / "segmentation.tsv").is_file()
        else None,
    )


__all__ = [
    "delly_assembly_sv",
    "delly_copy_number_variants",
    "delly_filter_call_set",
    "delly_genotype_sv_sites",
    "delly_merge_call_sets",
    "delly_plot_read_depth",
    "delly_translocations_to_breakends",
]
_register_public_verbose(globals())
