"""bcftools command-line facades (pack P06), core commands: convert to VCF, index counts, norm, polysomy.

Each public function runs one fixed ``bcftools`` command (1.24 is the tested version) without a shell and
with an explicit argument list; the user installs bcftools (conda, apt or a release) and puts it on PATH
or passes ``bcftools_path``. Outputs are written to new explicit paths (an existing file is never
overwritten); results are JSON-serialisable dictionaries with status, backend, version, command, outputs,
parsed values and a bounded log tail. A missing program raises ``BackendUnavailableError`` with an install
hint. The bcftools plugins of the same pack are in ``bcftools_plugins_p06_cli_functions`` and
``bcftools_trio_p06_cli_functions``; the older pysam-bundled ``bcftools norm`` / ``call -c`` / ``stats``
wrappers are in ``dnakit.bioinformatics.pysam_gap_functions``.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics import _htslib_p06 as h
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_CONVERT_FORMATS = ("gensample", "hapsample", "haplegendsample", "tsv", "gvcf")
_CONVERT_FILES = {"gensample": 2, "hapsample": 2, "haplegendsample": 3, "tsv": 1, "gvcf": 1}
_CONVERT_OPTION = {
    "gensample": "-G",
    "hapsample": "--hapsample2vcf",
    "haplegendsample": "-H",
    "tsv": "--tsv2vcf",
    "gvcf": "--gvcf2vcf",
}
_CHECK_REF = ("e", "w", "x", "s")
_RM_DUP = ("snps", "indels", "both", "all", "exact")
_MULTIALLELIC = re.compile(r"^[-+](?:snps|indels|both|any)$")
_SUMMARY_LINE = re.compile(r"^([A-Za-z][A-Za-z /-]*?):\s+(-?\d+)\s*$")


def _summary(log: str) -> dict[str, int]:
    """Counters bcftools prints on standard error (``Rows total:  3`` and similar)."""
    values: dict[str, int] = {}
    for line in log.splitlines():
        match = _SUMMARY_LINE.match(line.strip())
        if match:
            values[match.group(1).strip()] = int(match.group(2))
    return values


def _vcf_out(output: Path, output_type: str | None) -> tuple[list[str], str]:
    kind = h.output_type_for(output, output_type)
    return ["-O", kind, "-o", str(output)], kind


# ---------------------------------------------------------------------------
# Public operations
# ---------------------------------------------------------------------------


def bcftools_convert_to_vcf(
    input_format: str,
    input_files: Sequence[str | Path],
    output_path: str | Path,
    *,
    reference_fasta: str | Path | None = None,
    columns: str | None = None,
    sample_names: Sequence[str] | None = None,
    haploid_to_diploid: bool = False,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert genotype text formats to VCF/BCF (``bcftools convert``).

    ``input_format`` selects the source and the files that ``input_files`` must list in this order:
    ``gensample`` (IMPUTE2 ``gen`` file, ``sample`` file: ``-G``), ``hapsample`` (SHAPEIT ``hap``,
    ``sample``: ``--hapsample2vcf``), ``haplegendsample`` (``hap``, ``legend``, ``sample``: ``-H``),
    ``tsv`` (one genotype table such as a 23andMe export: ``--tsv2vcf``, needs ``reference_fasta`` and
    ``sample_names``; ``columns`` as in ``-c``, bcftools default ``ID,CHROM,POS,AA``) and ``gvcf`` (one
    gVCF whose reference blocks are expanded to one record per position: ``--gvcf2vcf``, needs
    ``reference_fasta``). ``haploid_to_diploid`` turns haploid hap genotypes into diploid homozygotes.
    The VCF/BCF is written to ``output_path``; ``output_type`` is ``v``, ``z``, ``b`` or ``u``
    (default: from the file name). Returns the output, the number of records (plain VCF only) and the
    counters bcftools reports (for example rows total, sites written, missing GTs).
    """
    kind = cli.choice(input_format, _CONVERT_FORMATS, "input_format")
    sources = h.existing_files(list(input_files), "input_files")
    if len(sources) != _CONVERT_FILES[kind]:
        raise ConfigurationError(
            f"input_format {kind!r} needs exactly {_CONVERT_FILES[kind]} input file(s)."
        )
    target = cli.new_file(output_path, "output_path")
    arguments = ["convert", _CONVERT_OPTION[kind]]
    if kind != "gvcf":
        arguments.append(",".join(str(path) for path in sources))
    if kind in {"tsv", "gvcf"} and reference_fasta is None:
        raise ConfigurationError(f"input_format {kind!r} needs reference_fasta.")
    if reference_fasta is not None:
        arguments += ["-f", str(cli.existing_file(reference_fasta, "reference_fasta"))]
    if kind == "tsv":
        if sample_names is None:
            raise ConfigurationError("input_format 'tsv' needs sample_names.")
        arguments += ["-s", ",".join(h.text_list(sample_names, "sample_names"))]
        if columns is not None:
            arguments += ["-c", cli.text(columns, "columns")]
    elif columns is not None or sample_names is not None:
        raise ConfigurationError("columns and sample_names apply only to input_format 'tsv'.")
    if h.flag(haploid_to_diploid, "haploid_to_diploid"):
        if kind not in {"hapsample", "haplegendsample"}:
            raise ConfigurationError("haploid_to_diploid applies to hapsample and haplegendsample.")
        arguments.append("--haploid2diploid")
    flags, out_kind = _vcf_out(target, output_type)
    arguments += h.selection(regions=regions, include=include, exclude=exclude)
    arguments += flags
    if kind == "gvcf":
        arguments.append(str(sources[0]))
    done = h.invoke("bcftools", bcftools_path, arguments, created=[target])
    return h.result(
        "bcftools_convert_to_vcf",
        done,
        outputs={"output": target},
        input_format=kind,
        output_type=out_kind,
        records=h.vcf_record_count(target),
        summary=_summary(done.stderr_tail),
    )


def bcftools_index_records(
    vcf_path: str | Path,
    *,
    per_contig: bool = True,
    include_empty_contigs: bool = False,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Record counts of a VCF/BCF read from its index (``bcftools index -s`` / ``-n``).

    The file needs an existing CSI or TBI index. ``per_contig`` True (``-s``) returns one row per
    contig with its name, length (None when the header has none) and number of records, plus the total;
    ``include_empty_contigs`` adds contigs without records (``-a``). ``per_contig`` False (``-n``)
    returns only the total.
    """
    source = cli.existing_file(vcf_path, "vcf_path")
    if h.flag(per_contig, "per_contig"):
        arguments = ["index", "-s"]
        if h.flag(include_empty_contigs, "include_empty_contigs"):
            arguments.append("-a")
    else:
        if include_empty_contigs:
            raise ConfigurationError("include_empty_contigs applies only with per_contig=True.")
        arguments = ["index", "-n"]
    done = h.invoke("bcftools", bcftools_path, [*arguments, str(source)])
    lines = [line for line in done.stdout.splitlines() if line.strip()]
    if per_contig:
        contigs: list[dict[str, Any]] = []
        for line in lines:
            name, length, count = [*line.split("\t"), "", ""][:3]
            contigs.append(
                {
                    "contig": name,
                    "length": None if length in {".", ""} else int(length),
                    "records": int(count),
                }
            )
        total = sum(row["records"] for row in contigs)
        return h.result("bcftools_index_records", done, contigs=contigs, total_records=total)
    return h.result("bcftools_index_records", done, total_records=int(lines[0]) if lines else 0)


def bcftools_norm_variants(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    reference_fasta: str | Path | None = None,
    multiallelics: str | None = None,
    atomize: bool = False,
    atom_overlaps: str = "*",
    check_ref: str | None = None,
    remove_duplicates: str | None = None,
    no_realign: bool = False,
    gff_annotation: str | Path | None = None,
    old_record_tag: str | None = None,
    strict_filter: bool = False,
    site_window: int = 1000,
    regions: str | None = None,
    include: str | None = None,
    exclude: str | None = None,
    output_type: str | None = None,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Normalise a VCF/BCF: left-align indels, split or join multiallelics, atomize, drop duplicates,
    check REF (``bcftools norm``).

    One function for all the modes of the program. ``reference_fasta`` (``-f``) enables left-alignment
    and the REF check; ``check_ref`` is ``e`` (error and stop, the bcftools default when -c is given),
    ``w`` (warn), ``x`` (exclude bad sites) or ``s`` (set REF from the reference). ``multiallelics`` is
    ``-snps``/``-indels``/``-both``/``-any`` to split or ``+snps``/``+indels``/``+both``/``+any`` to
    join (``-m``). ``atomize`` decomposes complex variants (``-a``; ``atom_overlaps`` ``*`` or ``.`` for
    overlapping alleles). ``remove_duplicates`` is ``snps``, ``indels``, ``both``, ``all`` or ``exact``
    (``-d``). ``no_realign`` skips realignment (``-N``). ``gff_annotation`` right-aligns in transcripts
    (HGVS 3' rule, ``-g``). ``old_record_tag`` annotates changed records with their original variant.
    At least one action must be requested. Plain VCF output is counted.
    """
    source = cli.existing_file(vcf_path, "vcf_path")
    target = cli.new_file(output_path, "output_path")
    arguments = ["norm"]
    if reference_fasta is not None:
        arguments += ["-f", str(cli.existing_file(reference_fasta, "reference_fasta"))]
    if multiallelics is not None:
        if not isinstance(multiallelics, str) or not _MULTIALLELIC.match(multiallelics):
            raise ConfigurationError(
                "multiallelics must be -/+ followed by snps, indels, both or any."
            )
        arguments += ["-m", multiallelics]
    if h.flag(atomize, "atomize"):
        arguments += [
            "-a",
            "--atom-overlaps",
            cli.choice(atom_overlaps, ("*", "."), "atom_overlaps"),
        ]
    if check_ref is not None:
        if reference_fasta is None:
            raise ConfigurationError("check_ref needs reference_fasta.")
        arguments += ["-c", cli.choice(check_ref, _CHECK_REF, "check_ref")]
    if remove_duplicates is not None:
        arguments += ["-d", cli.choice(remove_duplicates, _RM_DUP, "remove_duplicates")]
    if h.flag(no_realign, "no_realign"):
        arguments.append("-N")
    if gff_annotation is not None:
        arguments += ["-g", str(cli.existing_file(gff_annotation, "gff_annotation"))]
    if old_record_tag is not None:
        arguments += ["--old-rec-tag", cli.text(old_record_tag, "old_record_tag")]
    if h.flag(strict_filter, "strict_filter"):
        arguments.append("-s")
    if (
        reference_fasta is None
        and multiallelics is None
        and not atomize
        and remove_duplicates is None
    ):
        raise ConfigurationError(
            "Nothing to do: give reference_fasta, multiallelics, atomize or remove_duplicates."
        )
    arguments += ["-w", str(cli.integer(site_window, "site_window", minimum=1))]
    arguments += h.selection(regions=regions, include=include, exclude=exclude)
    flags, out_kind = _vcf_out(target, output_type)
    done = h.invoke("bcftools", bcftools_path, [*arguments, *flags, str(source)], created=[target])
    return h.result(
        "bcftools_norm_variants",
        done,
        outputs={"output": target},
        output_type=out_kind,
        records=h.vcf_record_count(target),
        summary=_summary(done.stderr_tail),
    )


def bcftools_polysomy(
    vcf_path: str | Path,
    output_dir: str | Path,
    *,
    sample: str | None = None,
    regions: str | None = None,
    targets: str | None = None,
    peak_size: float = 0.1,
    cn_penalty: float = 0.7,
    fit_threshold: float = 3.3,
    min_fraction: float = 0.1,
    peak_symmetry: float = 0.5,
    include_aa: bool = False,
    bcftools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Detect whole-chromosome copy-number aberrations from B-allele frequencies (``bcftools polysomy``).

    The VCF/BCF needs ``FORMAT/BAF`` (Illumina B-allele frequency). Per chromosome, peaks of the BAF
    distribution are fitted for 2, 3 and 4 copies; the best fit is the estimated copy number. The new
    folder ``output_dir`` receives ``dist.dat`` (the BAF distributions, the fitted functions and the
    ``CN`` lines) and ``dist.py`` (the matplotlib script that draws them; bcftools does not run it).
    ``sample`` selects the sample (``-s``), ``regions`` / ``targets`` restrict the sites (``-r`` / ``-t``).
    Algorithm options: ``peak_size`` minimum peak size 0-1 (``-b``), ``cn_penalty`` penalty for a higher
    copy number 0-1 (``-c``), ``fit_threshold`` goodness-of-fit threshold (``-f``, smaller is stricter),
    ``min_fraction`` minimum distinguishable fraction of aberrant cells (``-m``), ``peak_symmetry``
    peak symmetry threshold 0-1 (``-p``) and ``include_aa`` the AA peak in the 2 and 3 copy evaluation
    (``-i``). Returns ``copy_numbers`` (chromosome, estimated copy number, absolute fit deviation). Needs a
    bcftools built with GSL (``--enable-libgsl``; the conda builds are); otherwise
    ``BackendUnavailableError``.
    """
    source = cli.existing_file(vcf_path, "vcf_path")
    folder = cli.new_directory(output_dir, "output_dir")
    try:
        arguments = ["polysomy", "-o", str(folder)]
        if sample is not None:
            arguments += ["-s", cli.text(sample, "sample")]
        arguments += h.selection(regions=regions, targets=targets)
        arguments += [
            "-b",
            repr(cli.number(peak_size, "peak_size", minimum=0.0, maximum=1.0)),
            "-c",
            repr(cli.number(cn_penalty, "cn_penalty", minimum=0.0, maximum=1.0)),
            "-f",
            repr(cli.number(fit_threshold, "fit_threshold", minimum=0.0)),
            "-m",
            repr(cli.number(min_fraction, "min_fraction", minimum=0.0, maximum=1.0)),
            "-p",
            repr(cli.number(peak_symmetry, "peak_symmetry", minimum=0.0, maximum=1.0)),
        ]
        if h.flag(include_aa, "include_aa"):
            arguments.append("-i")
        try:
            done = h.invoke("bcftools", bcftools_path, [*arguments, str(source)], cwd=folder)
        except BackendExecutionError as error:
            if "unrecognized command" in str(error.context.get("log_tail", "")):
                raise BackendUnavailableError(
                    "This bcftools was built without polysomy support (GSL).",
                    hint=h.HINTS["bcftools"],
                ) from error
            raise
        data = folder / "dist.dat"
        if not data.is_file():
            raise BackendExecutionError(
                "bcftools polysomy wrote no dist.dat.", code="HTSLIB_OUTPUT_MISSING"
            )
    except BaseException:
        for leftover in folder.iterdir() if folder.is_dir() else []:
            if leftover.is_file():
                leftover.unlink()
        if folder.is_dir():
            folder.rmdir()
        raise
    copy_numbers = []
    for line in h.read_text(data).splitlines():
        fields = line.split("\t")
        if fields[0] == "CN" and len(fields) >= 3:
            copy_numbers.append(
                {
                    "chromosome": fields[1],
                    "copy_number": h.number(fields[2]),
                    "fit_deviation": h.number(fields[3]) if len(fields) > 3 else None,
                }
            )
    return h.result(
        "bcftools_polysomy",
        done,
        outputs={"output_dir": folder, "data": data},
        files=h.directory_files(folder),
        copy_numbers=copy_numbers,
    )


__all__ = [
    "bcftools_convert_to_vcf",
    "bcftools_index_records",
    "bcftools_norm_variants",
    "bcftools_polysomy",
]

_register_public_verbose(globals(), tuple(__all__))
