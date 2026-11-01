"""pysam-bundled samtools/bcftools facades and Phred codecs for operations DNAKit lacked.

Every samtools or bcftools call runs in a child interpreter, ``sys.executable -I -c <fixed
snippet>``. The snippet imports ``pysam.samtools`` or ``pysam.bcftools``, calls one whitelisted
subcommand with the argument list it receives as JSON, and writes the captured standard output to a
file. The child runs under ``dnakit.backends.external.execute_bounded_command`` (no shell, timeout,
output and artifact limits), so a crash or abort inside the bundled C code ends only the child and
is reported as a ``BackendExecutionError`` with a bounded log tail. Inputs are read in place; results
are written to explicit new paths and never overwrite an existing file.

The bundled tools are samtools 1.21 and bcftools 1.21 as shipped with pysam (``bio-reads`` extra).
Phred codecs are pure Python and do not need pysam.
"""

from __future__ import annotations

import importlib.metadata
import importlib.util
import json
import math
import re
import sys
import tempfile
import uuid
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_PACKAGE = "pysam"
_EXTRA = "bio-reads"
_INSTALL_HINT = 'python -m pip install "dnakit[bio-reads]" (pysam 0.22 to 0.24 bundles samtools and bcftools)'
_TIMEOUT_SECONDS = 1800.0
_MAX_OUTPUT_BYTES = 5_000_000
_MAX_ARTIFACT_BYTES = 2_000_000_000
_LOG_TAIL_CHARS = 4000
_REGION = re.compile(r"^[^\s:,]+(?::\d[\d,]*-\d[\d,]*)?$")
_MULTIALLELIC = ("-any", "+any", "-snps", "+snps", "-indels", "+indels", "-both", "+both")
_CONSENSUS_MODES = ("simple", "bayesian")
_PHASE_CODES = ("PS", "FL", "M0", "M1", "M2", "EV", "//")
_OFFSETS = (0, 33, 64)

# The only subcommands the child interpreter may call. Each public function passes one of these names.
_SUBCOMMANDS = {
    "samtools": frozenset(
        {"consensus", "flagstat", "stats", "ampliconstats", "cram_size", "depad", "phase",
         "targetcut", "reference", "rmdup"}
    ),
    "bcftools": frozenset({"mpileup", "call", "norm", "stats"}),
}

# Fixed snippet: the only Python code that runs in the child. It receives the program, the subcommand
# (checked against _SUBCOMMANDS by the parent), the arguments as a JSON list and the stdout file path.
_SNIPPET = """\
import json
import os
import sys

import pysam.bcftools
import pysam.samtools

program, subcommand, arguments_json, stdout_path = sys.argv[1:5]
module = pysam.samtools if program == "samtools" else pysam.bcftools
# Some pysam entry points write to file descriptor 1 instead of returning the text, so the output
# file becomes descriptor 1 before the call; a returned string or bytes is appended to it.
descriptor = os.open(stdout_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
os.dup2(descriptor, 1)
try:
    result = getattr(module, subcommand)(*json.loads(arguments_json))
except Exception as exc:
    sys.stderr.write("pysam {} {} failed: {}\\n".format(program, subcommand, exc))
    sys.exit(3)
if isinstance(result, str):
    result = result.encode("utf-8", "surrogateescape")
if result:
    with os.fdopen(1, "wb", closefd=False) as handle:
        handle.write(result)
"""


# ---------------------------------------------------------------------------
# Bounded child execution
# ---------------------------------------------------------------------------


def _pysam_version() -> str:
    """Read the installed version without importing pysam into this process."""
    if importlib.util.find_spec(_PACKAGE) is None:
        raise BackendUnavailableError(
            "pysam is not installed; it provides the samtools and bcftools code used here.",
            hint=f"Install with: {_INSTALL_HINT}",
        )
    return importlib.metadata.version(_PACKAGE)


def _tail(text: str) -> str:
    return text[-_LOG_TAIL_CHARS:]


def _child(
    program: str,
    subcommand: str,
    arguments: Sequence[str],
    workdir: Path,
    *,
    verbose: bool,
) -> tuple[Path, str]:
    """Run ``pysam.<program>.<subcommand>(*arguments)`` in a bounded child and keep its stdout.

    Returns the path of the file holding the captured standard output and the bounded log.
    """
    if subcommand not in _SUBCOMMANDS.get(program, frozenset()):
        raise ConfigurationError(f"pysam {program} {subcommand!r} is not allowed by this facade.")
    stdout_path = workdir / f"{program}-{subcommand}-{uuid.uuid4().hex}.stdout"
    result = execute_bounded_command(
        sys.executable,
        (
            "-I",
            "-c",
            _SNIPPET,
            program,
            subcommand,
            json.dumps([str(argument) for argument in arguments]),
            str(stdout_path),
        ),
        backend_id=f"pysam-{program}",
        cwd=workdir,
        timeout_seconds=_TIMEOUT_SECONDS,
        max_output_bytes=_MAX_OUTPUT_BYTES,
        monitored_output_paths=(stdout_path,),
        max_monitored_output_bytes=_MAX_ARTIFACT_BYTES,
        verbose=verbose,
    )
    if result.return_code != 0:
        raise BackendExecutionError(
            f"pysam {program} {subcommand} exited with code {result.return_code}.",
            code="PYSAM_EXECUTION_FAILED",
            context={"log_tail": _tail(result.output)},
        )
    if not stdout_path.is_file():
        raise BackendExecutionError(
            f"pysam {program} {subcommand} did not produce its output.",
            code="PYSAM_OUTPUT_MISSING",
            context={"log_tail": _tail(result.output)},
        )
    return stdout_path, _tail(result.output)


def _run_pysam(
    program: str,
    subcommand: str,
    arguments: Sequence[str],
    *,
    verbose: bool,
) -> tuple[bytes, str]:
    """Run one child call in a fresh temporary directory; return its stdout bytes and log tail."""
    _pysam_version()
    with tempfile.TemporaryDirectory(prefix="dnakit-pysam-") as directory:
        stdout_path, log = _child(program, subcommand, arguments, Path(directory), verbose=verbose)
        return stdout_path.read_bytes(), log


def _run_pysam_to_file(
    program: str,
    subcommand: str,
    arguments: Sequence[str],
    output: Path,
    *,
    verbose: bool,
) -> tuple[int, str]:
    """Run a child call whose stdout is the result, and write it to a new output file."""
    data, log = _run_pysam(program, subcommand, arguments, verbose=verbose)
    with output.open("xb") as handle:
        handle.write(data)
    return len(data), log


# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------


def _text(value: object, name: str) -> str:
    if isinstance(value, Path):
        value = str(value)
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n"):
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value.strip()


def _input(value: object, name: str = "source") -> Path:
    path = Path(_text(value, name)).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"Input file does not exist: {path}")
    return path


def _output(value: object, name: str = "output_path") -> Path:
    """A new file: its parent must exist and the file itself must not, so nothing is overwritten."""
    path = Path(_text(value, name)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    return path


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}.")
    return value


def _region(value: object) -> str:
    text = _text(value, "region")
    if not _REGION.match(text):
        raise ConfigurationError("region must look like CHROM or CHROM:START-END.")
    return text


def _provenance(program: str) -> dict[str, str]:
    return {"package": _PACKAGE, "version": _pysam_version(), "tool": program, "bundled_with": "pysam"}


def _result(operation: str, subcommand: str, log: str, **fields: Any) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": "completed",
        "backend": _PACKAGE,
        "operation": operation,
        "command": subcommand,
        "log_tail": log,
        "provenance": _provenance(subcommand.split()[0]),
    }
    result.update(fields)
    return result


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------


def _number(text: str) -> object:
    try:
        return int(text)
    except ValueError:
        pass
    try:
        value = float(text)
    except ValueError:
        return text
    return value if math.isfinite(value) else text


def _tab_sections(text: str) -> dict[str, list[list[object]]]:
    """Group tab-separated data lines by their first field; ``#`` lines are comments."""
    sections: dict[str, list[list[object]]] = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        sections.setdefault(fields[0], []).append([_number(field) for field in fields[1:]])
    return sections


def _flagstat_categories(text: str) -> list[dict[str, Any]]:
    categories: list[dict[str, Any]] = []
    for line in text.splitlines():
        match = re.match(r"^(\d+) \+ (\d+) (.+)$", line)
        if not match:
            continue
        label, _, detail = match.group(3).partition(" (")
        categories.append(
            {
                "label": label,
                "qc_passed": int(match.group(1)),
                "qc_failed": int(match.group(2)),
                "detail": detail.rstrip(")") if detail else None,
            }
        )
    return categories


_CRAM_BLOCK_KEYS = ("content_id", "uncompressed_size", "compressed_size", "ratio", "method", "data_series")


def _cram_blocks(text: str) -> list[dict[str, Any]]:
    blocks: list[dict[str, Any]] = []
    for line in text.splitlines():
        fields = line.split()
        if not fields or fields[0] != "BLOCK":
            continue
        values = fields[1:]
        row: dict[str, Any] = {}
        for key, value in zip(_CRAM_BLOCK_KEYS, values):
            row[key] = _number(value) if key.endswith("size") else value
        blocks.append(row)
    return blocks


def _phase_counts(text: str) -> dict[str, int]:
    counts = dict.fromkeys(_PHASE_CODES, 0)
    for line in text.splitlines():
        code = line.split("\t", 1)[0]
        if code in counts:
            counts[code] += 1
    return counts


# ---------------------------------------------------------------------------
# samtools facades
# ---------------------------------------------------------------------------


def samtools_consensus(
    source: str | Path,
    output_path: str | Path,
    *,
    mode: str = "simple",
    verbose: bool = True,
) -> dict[str, Any]:
    """Consensus FASTA from alignments (``samtools consensus``).

    ``mode`` is ``simple`` (per-column base counting, the samtools default) or ``bayesian`` (the
    Gap5-style model). The FASTA is written to ``output_path``; the result lists the sequence names
    and lengths found in it.
    """
    path = _input(source)
    model = _choice(mode, _CONSENSUS_MODES, "mode")
    output = _output(output_path)
    _, log = _run_pysam_to_file(
        "samtools", "consensus", ("--mode", model, str(path)), output, verbose=verbose
    )
    lengths = _fasta_lengths(output)
    return _result(
        "samtools_consensus",
        "samtools consensus",
        log,
        mode=model,
        output=str(output),
        sequences=lengths,
    )


def _fasta_lengths(path: Path) -> dict[str, int]:
    """Sequence name (first word of the header) to residue count, in file order."""
    lengths: dict[str, int] = {}
    name: str | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                words = line[1:].split()
                name = words[0] if words else ""
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def samtools_flagstat(source: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Counts per alignment category, split into QC-passed and QC-failed reads (``samtools flagstat``).

    Each category has ``qc_passed`` and ``qc_failed`` counts and the percentage text that samtools
    prints in parentheses (``detail``).
    """
    path = _input(source)
    data, log = _run_pysam("samtools", "flagstat", (str(path),), verbose=verbose)
    text = data.decode("utf-8", "replace")
    categories = _flagstat_categories(text)
    return _result("samtools_flagstat", "samtools flagstat", log, categories=categories)


def samtools_stats(source: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Full ``samtools stats`` report: every section (SN, FFQ, LFQ, GCF, GCL, IS, RL, ID, COV, ...).

    ``sections`` maps each record code to its rows (numbers converted where the text is numeric);
    ``summary`` maps the SN summary keys (without the trailing colon) to their values.
    """
    path = _input(source)
    data, log = _run_pysam("samtools", "stats", (str(path),), verbose=verbose)
    sections = _tab_sections(data.decode("utf-8", "replace"))
    summary = {
        str(row[0]).rstrip(":"): row[1] for row in sections.get("SN", []) if len(row) >= 2
    }
    return _result("samtools_stats", "samtools stats", log, summary=summary, sections=sections)


def samtools_ampliconstats(
    primers_bed: str | Path, source: str | Path, *, verbose: bool = True
) -> dict[str, Any]:
    """Amplicon statistics for primer-clipped alignments (``samtools ampliconstats``).

    ``primers_bed`` lists the primers in pairs (forward then reverse, per amplicon). The result has
    the SS summary, the AMPLICON layout rows and every other statistics section of samtools.
    """
    bed = _input(primers_bed, "primers_bed")
    path = _input(source)
    data, log = _run_pysam(
        "samtools", "ampliconstats", (str(bed), str(path)), verbose=verbose
    )
    sections = _tab_sections(data.decode("utf-8", "replace"))
    return _result(
        "samtools_ampliconstats", "samtools ampliconstats", log, sections=sections
    )


def samtools_cram_size(source: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Block sizes of a CRAM file: content ID, uncompressed and compressed size, method (``samtools cram_size``)."""
    path = _input(source)
    data, log = _run_pysam("samtools", "cram_size", (str(path),), verbose=verbose)
    blocks = _cram_blocks(data.decode("utf-8", "replace"))
    return _result("samtools_cram_size", "samtools cram_size", log, blocks=blocks)


def samtools_depad(
    source: str | Path,
    reference: str | Path,
    output_path: str | Path,
    *,
    output_format: str = "bam",
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert alignments made against a padded reference (``*`` pad columns) to unpadded coordinates.

    ``reference`` is the padded FASTA the reads were aligned to. The alignments must contain P
    (padding) CIGAR operations. ``output_format`` is ``bam`` (default) or ``sam``. The same
    operation is also exposed by pysam under the alias ``pad2unpad``.
    """
    path = _input(source)
    ref = _input(reference, "reference")
    fmt = _choice(output_format, ("bam", "sam"), "output_format")
    output = _output(output_path)
    options = ["-T", str(ref)]
    if fmt == "sam":
        options.append("-s")
    _, log = _run_pysam_to_file(
        "samtools", "depad", (*options, str(path)), output, verbose=verbose
    )
    return _result(
        "samtools_depad", "samtools depad", log, output=str(output), output_format=fmt
    )


def samtools_phase_report(
    source: str | Path,
    output_path: str | Path,
    *,
    block_length: int = 13,
    min_het_phred: int = 37,
    min_base_quality: int = 13,
    max_depth: int = 256,
    fix_chimeras: bool = True,
    drop_ambiguous: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call and phase heterozygous SNPs from reads (``samtools phase``, report on standard output).

    The report (phase sets ``PS``, filtered regions ``FL``, markers ``M0``/``M1``/``M2`` for
    singletons, phased and filtered sites, and supporting-read ``EV`` lines) is written to
    ``output_path``. ``counts`` gives the number of lines per code (``//`` closes a phase set).
    """
    path = _input(source)
    output = _output(output_path)
    arguments = [
        "-k", str(_integer(block_length, "block_length", minimum=1, maximum=10_000)),
        "-q", str(_integer(min_het_phred, "min_het_phred", minimum=0, maximum=1_000)),
        "-Q", str(_integer(min_base_quality, "min_base_quality", minimum=0, maximum=93)),
        "-D", str(_integer(max_depth, "max_depth", minimum=1, maximum=1_000_000)),
    ]
    if not fix_chimeras:
        arguments.append("-F")
    if drop_ambiguous:
        arguments.append("-A")
    arguments.append(str(path))
    size, log = _run_pysam_to_file("samtools", "phase", arguments, output, verbose=verbose)
    counts = _phase_counts(output.read_text(encoding="utf-8", errors="replace"))
    return _result(
        "samtools_phase_report",
        "samtools phase",
        log,
        output=str(output),
        bytes=size,
        counts=counts,
    )


def samtools_targetcut(
    source: str | Path,
    output_path: str | Path,
    *,
    reference: str | Path | None = None,
    min_base_quality: int | None = None,
    insertion_penalty: int | None = None,
    error_rate_0: int | None = None,
    error_rate_1: int | None = None,
    error_rate_2: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut target regions from pooled reads and build a haploid consensus per target (``samtools targetcut``).

    Options left as ``None`` keep the samtools defaults. The consensus records (SAM text) are written
    to ``output_path``.
    """
    path = _input(source)
    output = _output(output_path)
    arguments: list[str] = []
    if reference is not None:
        arguments += ["-f", str(_input(reference, "reference"))]
    if min_base_quality is not None:
        arguments += ["-Q", str(_integer(min_base_quality, "min_base_quality", minimum=0, maximum=93))]
    if insertion_penalty is not None:
        arguments += ["-i", str(_integer(insertion_penalty, "insertion_penalty", minimum=0, maximum=1_000))]
    for flag, value, name in (
        ("-0", error_rate_0, "error_rate_0"),
        ("-1", error_rate_1, "error_rate_1"),
        ("-2", error_rate_2, "error_rate_2"),
    ):
        if value is not None:
            arguments += [flag, str(_integer(value, name, minimum=0, maximum=1_000))]
    arguments.append(str(path))
    size, log = _run_pysam_to_file("samtools", "targetcut", arguments, output, verbose=verbose)
    return _result("samtools_targetcut", "samtools targetcut", log, output=str(output), bytes=size)


def samtools_embedded_reference(
    source: str | Path,
    output_path: str | Path,
    *,
    region: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the reference sequences embedded in a CRAM file (``samtools reference -e``).

    Only the embedded-reference mode is provided here; the MD-tag reconstruction of ``samtools
    reference`` is ``reconstruct_reference_from_alignments``. The CRAM must have been written with
    embedded references (``embed_ref=1``).
    """
    path = _input(source)
    output = _output(output_path)
    arguments = ["-e"]
    if region is not None:
        arguments += ["-r", _region(region)]
    arguments.append(str(path))
    size, log = _run_pysam_to_file("samtools", "reference", arguments, output, verbose=verbose)
    return _result(
        "samtools_embedded_reference",
        "samtools reference -e",
        log,
        output=str(output),
        sequences=_fasta_lengths(output),
    )


def samtools_rmdup(
    source: str | Path,
    output_path: str | Path,
    *,
    single_end: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove duplicate alignments from a coordinate-sorted BAM (``samtools rmdup``).

    ``single_end=True`` treats all reads as single-end (``-s``); otherwise paired-end reads are
    handled as pairs. Output is a BAM file at ``output_path``.
    """
    path = _input(source)
    output = _output(output_path)
    if output.suffix.lower() != ".bam":
        raise ConfigurationError("rmdup output must be a .bam file.")
    arguments = ["-s"] if single_end else []
    with tempfile.TemporaryDirectory(prefix="dnakit-pysam-") as directory:
        workdir = Path(directory)
        target = workdir / "rmdup-output.bam"
        _, log = _child(
            "samtools",
            "rmdup",
            (*arguments, str(path), str(target)),
            workdir,
            verbose=verbose,
        )
        if not target.is_file():
            raise BackendExecutionError(
                "samtools rmdup did not write its output file.",
                code="PYSAM_OUTPUT_MISSING",
                context={"log_tail": log},
            )
        data = target.read_bytes()
    with output.open("xb") as handle:
        handle.write(data)
    return _result(
        "samtools_rmdup",
        "samtools rmdup",
        log,
        output=str(output),
        single_end=single_end,
        bytes=len(data),
    )


# ---------------------------------------------------------------------------
# bcftools facades
# ---------------------------------------------------------------------------


def bcftools_normalize_vcf(
    vcf_path: str | Path,
    output_path: str | Path,
    *,
    reference_fasta: str | Path | None = None,
    multiallelics: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Normalise variants and write the normalised VCF (``bcftools norm``).

    ``reference_fasta`` enables left-alignment and reference checks (``-f``). ``multiallelics`` is one
    of ``-any``/``+any``/``-snps``/``+snps``/``-indels``/``+indels``/``-both``/``+both`` (``-m``) to split
    (``-``) or join (``+``) multiallelic records. The normalised VCF is written to ``output_path``.
    """
    path = _input(vcf_path, "vcf_path")
    output = _output(output_path)
    if reference_fasta is None and multiallelics is None:
        raise ConfigurationError("bcftools norm needs reference_fasta, multiallelics, or both.")
    arguments: list[str] = []
    if reference_fasta is not None:
        arguments += ["-f", str(_input(reference_fasta, "reference_fasta"))]
    if multiallelics is not None:
        arguments += ["-m", _choice(multiallelics, _MULTIALLELIC, "multiallelics")]
    arguments.append(str(path))
    size, log = _run_pysam_to_file("bcftools", "norm", arguments, output, verbose=verbose)
    records = _vcf_record_count(output)
    return _result(
        "bcftools_normalize_vcf",
        "bcftools norm",
        log,
        output=str(output),
        bytes=size,
        records=records,
    )


def _vcf_record_count(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def bcftools_call_consensus(
    source: str | Path,
    reference: str | Path,
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """SNP and indel calls with the original consensus caller (``bcftools mpileup`` then ``call -c``).

    The alignment is piled up against ``reference`` (which needs a ``.fai`` index), and the pileup
    is passed to ``bcftools call -c``. The multiallelic caller (``-m``) is ``ngs.call_variants``.
    """
    path = _input(source)
    ref = _input(reference, "reference")
    output = _output(output_path)
    with tempfile.TemporaryDirectory(prefix="dnakit-pysam-") as directory:
        workdir = Path(directory)
        pileup, pileup_log = _child(
            "bcftools", "mpileup", ("-f", str(ref), str(path)), workdir, verbose=verbose
        )
        calls, call_log = _child(
            "bcftools", "call", ("-c", "-v", str(pileup)), workdir, verbose=verbose
        )
        data = calls.read_bytes()
    with output.open("xb") as handle:
        handle.write(data)
    return _result(
        "bcftools_call_consensus",
        "bcftools mpileup | bcftools call -c -v",
        pileup_log + "\n" + call_log,
        output=str(output),
        bytes=len(data),
        records=_vcf_record_count(output),
    )


def bcftools_stats(
    vcf_path: str | Path,
    *,
    second_vcf: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Full ``bcftools stats`` report for one VCF/BCF, or the comparison of two (``second_vcf``).

    ``sections`` maps every record code (SN, TSTV, SiS, AF, QUAL, IDD, ST, DP, PSC, ...) to its rows;
    ``summary`` maps the SN keys to their values. The two-file comparison needs bgzip-compressed or
    BCF inputs (bcftools refuses plain VCF there), and its ``ID`` rows name the two sets.
    """
    path = _input(vcf_path, "vcf_path")
    arguments = [str(path)]
    if second_vcf is not None:
        second = _input(second_vcf, "second_vcf")
        for item in (path, second):
            if item.read_bytes()[:2] != b"\x1f\x8b" and item.read_bytes()[:3] != b"BCF":
                raise ConfigurationError(
                    f"Two-file bcftools stats needs bgzip-compressed input (.vcf.gz or .bcf): {item}"
                )
        arguments.append(str(second))
    data, log = _run_pysam("bcftools", "stats", tuple(arguments), verbose=verbose)
    sections = _tab_sections(data.decode("utf-8", "replace"))
    summary = {str(row[0]).rstrip(":"): row[1] for row in sections.get("SN", []) if len(row) >= 2}
    return _result(
        "bcftools_stats",
        "bcftools stats",
        log,
        summary=summary,
        sections=sections,
        files=2 if second_vcf is not None else 1,
    )


# ---------------------------------------------------------------------------
# Phred quality codecs (pysam.qualitystring_to_array and friends, in pure Python)
# ---------------------------------------------------------------------------


def _offset(value: object) -> int:
    return _integer(value, "offset", minimum=0, maximum=127)


def decode_phred_quality_string(
    quality: str,
    *,
    offset: int = 33,
    verbose: bool = True,
) -> list[int]:
    """Convert a quality string to Phred scores: each character minus ``offset`` (Sanger 33, Illumina 1.3+ 64).

    Raises ``ConfigurationError`` for a character that is below the offset or above 126.
    """
    shift = _offset(offset)
    if not isinstance(quality, str):
        raise ConfigurationError("quality must be a string.")
    scores: list[int] = []
    for character in quality:
        code = ord(character)
        if code < shift or code > 126:
            raise ConfigurationError(
                f"Quality character {character!r} is outside the range for offset {shift}."
            )
        scores.append(code - shift)
    return scores


def encode_phred_qualities(
    qualities: Sequence[int],
    *,
    offset: int = 33,
    verbose: bool = True,
) -> str:
    """Convert Phred scores to a quality string: each score plus ``offset``.

    Scores must be integers from 0 to ``126 - offset``. This also serves the pysam names
    ``array_to_qualitystring`` and ``qualities_to_qualitystring``, which take an array or a list.
    """
    shift = _offset(offset)
    characters: list[str] = []
    for score in qualities:
        if isinstance(score, bool) or not isinstance(score, int):
            raise ConfigurationError("Phred scores must be integers.")
        if not 0 <= score <= 126 - shift:
            raise ConfigurationError(f"Phred score {score} cannot be encoded with offset {shift}.")
        characters.append(chr(score + shift))
    return "".join(characters)


__all__ = (
    "samtools_consensus",
    "samtools_flagstat",
    "samtools_stats",
    "samtools_ampliconstats",
    "samtools_cram_size",
    "samtools_depad",
    "samtools_phase_report",
    "samtools_targetcut",
    "samtools_embedded_reference",
    "samtools_rmdup",
    "bcftools_normalize_vcf",
    "bcftools_call_consensus",
    "bcftools_stats",
    "decode_phred_quality_string",
    "encode_phred_qualities",
)

_register_public_verbose(globals(), tuple(__all__))
