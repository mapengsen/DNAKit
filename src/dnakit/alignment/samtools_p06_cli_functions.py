"""samtools command-line facades (pack P06): FASTQ region extraction, coverage histograms, checksums.

Each public function runs one fixed ``samtools`` subcommand (1.24 is the tested version) without a
shell and with an explicit argument list; the user installs samtools (conda, apt or a release) and puts
it on PATH or passes ``samtools_path``. Outputs are written to new explicit paths (an existing file is
never overwritten); results are JSON-serialisable dictionaries with status, backend, version, command,
outputs, parsed rows and a bounded log tail. A missing program raises ``BackendUnavailableError`` with an
install hint. The other samtools operations that DNAKit offers live in
``dnakit.bioinformatics.pysam_gap_functions`` (consensus, flagstat, stats, ampliconstats, phase,
embedded reference).
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics import _htslib_p06 as h
from dnakit.exceptions import ConfigurationError

_STRAND_MARKS = ("rc", "no", "sign")
_FLAG_SPEC = re.compile(r"^(?:0x[0-9A-Fa-f]+|\d+|[A-Z0-9_]+(?:,[A-Z0-9_]+)*)$")
_HIST_STAT = re.compile(r"([A-Za-z][A-Za-z ]*?):\s+(\S+)")
_CONTIG = re.compile(r"^(\S+) \(([^)]*)\)$")


def _flag_spec(value: str | None, name: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, str) or not _FLAG_SPEC.match(value):
        raise ConfigurationError(f"{name} must be a number (decimal or 0x hex) or comma-joined flag names.")
    return [value]


def _fastq_records(path: Path, limit: int = 1_000_000) -> list[dict[str, Any]]:
    """Name and sequence length of every record of a (possibly line-wrapped) FASTQ file."""
    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        name: str | None = None
        length = 0
        in_quality = False
        quality = 0
        for raw in handle:
            line = raw.rstrip("\n")
            if name is None:
                if line.startswith("@"):
                    name, length, in_quality, quality = line[1:], 0, False, 0
                continue
            if not in_quality:
                if line.startswith("+"):
                    in_quality = True
                    if length == 0:  # samtools -c writes an empty record for a missing region
                        records.append({"name": name, "length": 0})
                        name = None
                else:
                    length += len(line)
                continue
            quality += len(line)
            if quality >= length:
                records.append({"name": name, "length": length})
                name = None
                if len(records) >= limit:
                    break
    return records


# ---------------------------------------------------------------------------
# Public operations
# ---------------------------------------------------------------------------


def samtools_extract_fastq_regions(
    fastq_path: str | Path,
    output_path: str | Path,
    *,
    regions: Sequence[str] = (),
    region_file: str | Path | None = None,
    line_length: int | None = None,
    reverse_complement: bool = False,
    mark_strand: str | None = None,
    continue_on_missing: bool = False,
    index_path: str | Path | None = None,
    samtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract records or sub-regions from a FASTQ file through its index (``samtools fqidx``).

    Covers the work-list rows ``samtools faidx -f/--fastq REGION`` and ``samtools fqidx REGION``: both
    run the same FASTQ indexed extraction. ``regions`` are ``name``, ``name:start-end`` (1-based,
    inclusive) strings, or ``region_file`` lists them one per line. The FASTQ index (``.fai``) is created
    next to the input unless ``index_path`` names another place (use it for a read-only folder).
    ``reverse_complement`` reverse-complements the sequence and reverses the quality string;
    ``mark_strand`` is ``rc`` (the samtools default: ``/rc`` after the name), ``no`` or ``sign``;
    ``continue_on_missing`` keeps going after a region that cannot be retrieved. ``line_length`` is the
    sequence line width (samtools default when None). Returns the written FASTQ path and the extracted
    record names with their lengths.
    """
    source = cli.existing_file(fastq_path, "fastq_path")
    target = cli.new_file(output_path, "output_path")
    names = [cli.text(region, "region") for region in regions]
    if not names and region_file is None:
        raise ConfigurationError("Give at least one region or a region_file.")
    arguments = ["fqidx"]
    if line_length is not None:
        arguments += ["-n", str(cli.integer(line_length, "line_length", minimum=0))]
    if h.flag(continue_on_missing, "continue_on_missing"):
        arguments.append("-c")
    if h.flag(reverse_complement, "reverse_complement"):
        arguments.append("-i")
    if mark_strand is not None:
        arguments += ["--mark-strand", cli.choice(mark_strand, _STRAND_MARKS, "mark_strand")]
    if region_file is not None:
        arguments += ["-r", str(cli.existing_file(region_file, "region_file"))]
    if index_path is not None:
        index = Path(index_path).expanduser().resolve()
        if not index.parent.is_dir():
            raise ConfigurationError(f"The folder for index_path does not exist: {index.parent}")
        arguments += ["--fai-idx", str(index)]
    arguments += ["-o", str(target), str(source), *names]
    done = h.invoke("samtools", samtools_path, arguments, created=[target])
    return h.result(
        "samtools_extract_fastq_regions",
        done,
        outputs={"output": target},
        records=_fastq_records(target),
    )


def samtools_coverage_histogram(
    alignment_paths: Sequence[str | Path],
    output_path: str | Path,
    *,
    region: str | None = None,
    bins: int = 40,
    ascii_only: bool = False,
    min_read_length: int = 0,
    min_mapping_quality: int = 0,
    min_base_quality: int = 0,
    required_flags: str | None = None,
    excluded_flags: str | None = None,
    max_depth: int = 1_000_000,
    min_depth: int = 1,
    reference_fasta: str | Path | None = None,
    samtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Coverage histogram drawn as text, one panel per reference sequence (``samtools coverage -m``).

    The panels are written to ``output_path`` (UTF-8 block characters, or plain ASCII when
    ``ascii_only`` is True, ``-A``). ``bins`` is the number of histogram columns (``-w``; samtools
    uses the terminal width minus 40, which depends on the output when no terminal is attached, so the
    facade fixes it to 40, the value of an 80-column terminal). ``region`` is ``chr:start-end`` (``-r``, needs an index).
    ``required_flags`` / ``excluded_flags`` are ``--rf`` / ``--ff`` (a number or flag names such as
    ``UNMAP,SECONDARY``; samtools excludes UNMAP, SECONDARY, QCFAIL and DUP when unset). Returns the
    path and, per sequence, the numbers printed beside each panel (reads, covered bases, percent
    covered, mean coverage, mean base and mapping quality).
    """
    sources = h.existing_files(list(alignment_paths), "alignment_paths")
    target = cli.new_file(output_path, "output_path")
    arguments = ["coverage", "-m"]
    if h.flag(ascii_only, "ascii_only"):
        arguments.append("-A")
    arguments += ["-w", str(cli.integer(bins, "bins", minimum=1, maximum=1000))]
    if region is not None:
        arguments += ["-r", cli.text(region, "region")]
    arguments += [
        "-l", str(cli.integer(min_read_length, "min_read_length", minimum=0)),
        "-q", str(cli.integer(min_mapping_quality, "min_mapping_quality", minimum=0)),
        "-Q", str(cli.integer(min_base_quality, "min_base_quality", minimum=0)),
        "-d", str(cli.integer(max_depth, "max_depth", minimum=0)),
        "--min-depth", str(cli.integer(min_depth, "min_depth", minimum=0)),
    ]
    for option, value in (("--rf", required_flags), ("--ff", excluded_flags)):
        for item in _flag_spec(value, option):
            arguments += [option, item]
    if reference_fasta is not None:
        arguments += ["--reference", str(cli.existing_file(reference_fasta, "reference_fasta"))]
    arguments += ["-o", str(target), *(str(path) for path in sources)]
    done = h.invoke("samtools", samtools_path, arguments, created=[target])
    contigs: list[dict[str, Any]] = []
    for line in h.read_text(target).splitlines():
        match = _CONTIG.match(line.strip())
        if match:
            contigs.append({"name": match.group(1), "length": match.group(2), "statistics": {}})
            continue
        if contigs and ("│" in line or "|" in line):
            tail = re.split(r"[│|]", line)[-1]
            for key, value in _HIST_STAT.findall(tail):
                contigs[-1]["statistics"][key.strip()] = h.number(value.rstrip("x%"))
    return h.result("samtools_coverage_histogram", done, outputs={"output": target}, contigs=contigs)


def samtools_checksum(
    input_paths: Sequence[str | Path],
    output_path: str | Path,
    *,
    merge: bool = False,
    exclude_flags: str | None = None,
    require_flags: str | None = None,
    flag_mask: str | None = None,
    no_reverse_complement: bool = False,
    tags: Sequence[str] | None = None,
    in_order: bool = False,
    check_position: bool = False,
    check_cigar: bool = False,
    check_mate: bool = False,
    max_records: int = 0,
    show_qc: bool = False,
    tab_separated: bool = False,
    bamseqchksum_format: bool = False,
    samtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Order-agnostic checksums of sequence, quality and name content (``samtools checksum``).

    With ``merge`` False (default) the SAM/BAM/CRAM files in ``input_paths`` are checksummed and the
    report is written to ``output_path``. With ``merge`` True (``samtools checksum -m``) the inputs are
    checksum reports written earlier by this function (use the same options for all of them,
    ``tab_separated`` must match how they were written) and the merged report is written. The
    checksum is independent of record order unless ``in_order`` (``-O``); ``check_position``,
    ``check_cigar`` and ``check_mate`` add CHR/POS, MAPQ/CIGAR and PNEXT/RNEXT/TLEN (``-P -C -M``);
    ``exclude_flags`` / ``require_flags`` / ``flag_mask`` are ``-F`` / ``-f`` / ``-b`` (samtools
    defaults 0x900, 0, 0x0c1); ``tags`` the auxiliary tags (``-t``, default BC,FI,QT,RT,TC);
    ``max_records`` stops after that many records (``-N``, 0 = all); ``show_qc`` adds QC pass/fail
    lines; ``bamseqchksum_format`` writes the bamseqchksum layout. Returns the report path and one row
    per group with count and checksums.
    """
    sources = h.existing_files(list(input_paths), "input_paths")
    target = cli.new_file(output_path, "output_path")
    arguments = ["checksum"]
    if h.flag(merge, "merge"):
        arguments.append("-m")
    if h.flag(tab_separated, "tab_separated"):
        arguments.append("-T")
    if h.flag(show_qc, "show_qc"):
        arguments.append("-q")
    if h.flag(bamseqchksum_format, "bamseqchksum_format"):
        arguments.append("-B")
    if not merge:
        for option, value, name in (
            ("-F", exclude_flags, "exclude_flags"),
            ("-f", require_flags, "require_flags"),
            ("-b", flag_mask, "flag_mask"),
        ):
            for item in _flag_spec(value, name):
                arguments += [option, item]
        if h.flag(no_reverse_complement, "no_reverse_complement"):
            arguments.append("-c")
        if tags is not None:
            arguments += ["-t", ",".join(h.text_list(tags, "tags"))]
        for enabled, option, name in (
            (in_order, "-O", "in_order"),
            (check_position, "-P", "check_position"),
            (check_cigar, "-C", "check_cigar"),
            (check_mate, "-M", "check_mate"),
        ):
            if h.flag(enabled, name):
                arguments.append(option)
        count = cli.integer(max_records, "max_records", minimum=0)
        if count:
            arguments += ["-N", str(count)]
    arguments += ["-o", str(target), *(str(path) for path in sources)]
    done = h.invoke("samtools", samtools_path, arguments, created=[target])
    rows: list[dict[str, Any]] = []
    columns: list[str] | None = None
    for line in h.read_text(target).splitlines():
        if line.startswith("# Group"):
            columns = re.split(r"\t+" if tab_separated else r"\s+", line[1:].strip())
            continue
        if not line.strip() or line.startswith("#") or columns is None:
            continue
        fields = re.split(r"\t" if tab_separated else r"\s+", line.strip())
        row: dict[str, Any] = {}
        for column, field in zip(columns, fields, strict=False):
            row[column] = h.number(field) if column == "count" else field
        rows.append(row)
    return h.result(
        "samtools_checksum", done, outputs={"report": target}, merged=merge, groups=rows
    )


__all__ = [
    "samtools_checksum",
    "samtools_coverage_histogram",
    "samtools_extract_fastq_regions",
]

_register_public_verbose(globals(), tuple(__all__))
