"""pysam command-line operations on local files: padded-alignment conversion and copy-number calling.

``pysam.depad`` and ``pysam.pad2unpad`` run the samtools bundled with pysam to convert a padded alignment to
unpadded coordinates; they need the padded reference FASTA (``-T``). ``pysam.bcftools.cnv`` calls copy-number
states from a VCF that carries B-allele frequency (BAF) and log R ratio (LRR) FORMAT fields.

pysam writes command output to a temporary file unless ``catch_stdout=False`` is given, and it rewrites ``-o``
in that mode. Every call here therefore passes ``catch_stdout=False`` with an explicit output path.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-reads"
_MAX_RECORDS = 200


def _existing(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return resolved


def _new(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _new_directory(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The parent folder of {name} does not exist: {resolved.parent}")
    resolved.mkdir()
    return resolved


def _identifier(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or value.startswith("-"):
        raise ConfigurationError(f"{name} must be a non-empty name that does not start with '-'.")
    return value


def _samtools_error(exc: Exception) -> BackendExecutionError:
    return BackendExecutionError(f"pysam command failed: {exc}")


def _alignment_records(path: Path) -> tuple[int, list[dict[str, Any]]]:
    pysam = require("pysam", _EXTRA)
    count = 0
    records: list[dict[str, Any]] = []
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        for read in handle:
            count += 1
            if len(records) < _MAX_RECORDS:
                records.append(
                    {
                        "query_name": read.query_name,
                        "flag": int(read.flag),
                        "reference_name": read.reference_name,
                        "reference_start": None
                        if read.reference_start is None
                        else int(read.reference_start) + 1,
                        "cigar": read.cigarstring,
                        "sequence": read.query_sequence,
                    }
                )
    return count, records


def _depad(
    operation: str,
    alignment: str,
    padded_reference: str,
    output_alignment: str,
    *,
    output_sam: bool,
) -> dict[str, Any]:
    pysam = require("pysam", _EXTRA)
    source = _existing(alignment, "alignment")
    reference = _existing(padded_reference, "padded_reference")
    out = _new(output_alignment, "output_alignment")
    args = ["-T", str(reference), "-o", str(out)]
    if output_sam:
        args.insert(0, "-s")
    args.append(str(source))
    function = getattr(pysam, operation)
    try:
        function(*args, catch_stdout=False)
    except pysam.SamtoolsError as exc:
        raise _samtools_error(exc) from exc
    if not out.is_file():
        raise BackendExecutionError(f"samtools {operation} did not write {out}.")
    count, records = _alignment_records(out)
    return {
        "status": "completed",
        "backend": "pysam",
        "command": ["samtools", operation, *args],
        "output_format": "SAM" if output_sam else "BAM",
        "outputs": [str(out)],
        "record_count": count,
        "records": records,
    }


def depad_alignments(
    alignment: str,
    padded_reference: str,
    output_alignment: str,
    *,
    output_sam: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a padded alignment to unpadded coordinates (samtools depad through pysam.depad).

    ``padded_reference`` is the FASTA of the padded reference (pad columns written as ``*``); ``-T`` is
    passed to samtools. Reads must not use P or I CIGAR operators in the input. Output is BAM unless
    ``output_sam`` is True. The output path must not exist.
    """
    return _depad("depad", alignment, padded_reference, output_alignment, output_sam=output_sam)


def pad_to_unpad_alignments(
    alignment: str,
    padded_reference: str,
    output_alignment: str,
    *,
    output_sam: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a padded alignment to unpadded coordinates through the legacy samtools name pad2unpad (pysam.pad2unpad).

    Same inputs, output and checks as ``depad_alignments``; the samtools command ``pad2unpad`` produced
    identical records on the test fixture (see the batch notes).
    """
    return _depad("pad2unpad", alignment, padded_reference, output_alignment, output_sam=output_sam)


def _parse_tab(path: Path) -> list[list[str]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        rows.append(line.split("\t"))
    return rows


def call_copy_number_states(
    vcf: str,
    output_dir: str,
    *,
    query_sample: str,
    control_sample: str | None = None,
    allele_frequency_file: str | None = None,
    regions: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call copy-number states from BAF and LRR (bcftools cnv through pysam.bcftools.cnv).

    The VCF must contain the FORMAT fields ``BAF`` and ``LRR`` for the query sample. ``allele_frequency_file``
    (CHR, POS, REF,ALT, AF) and ``regions`` (a region string) are optional. The output directory must not exist;
    bcftools writes ``cn.<sample>.tab`` (per site), ``dat.<sample>.tab``, ``summary.<sample>.tab`` (segments) and
    ``plot.<sample>.py``. Returns the output paths, the segments (CN state per region) and the per-site states.
    """
    pysam = require("pysam", _EXTRA)
    cnv = require("pysam.bcftools", _EXTRA).cnv
    source = _existing(vcf, "vcf")
    sample = _identifier(query_sample, "query_sample")
    directory = _new_directory(output_dir, "output_dir")
    args = ["-s", sample, "-o", str(directory)]
    if control_sample is not None:
        args += ["-c", _identifier(control_sample, "control_sample")]
    if allele_frequency_file is not None:
        args += ["-f", str(_existing(allele_frequency_file, "allele_frequency_file"))]
    if regions is not None:
        if not regions.strip() or regions.startswith("-"):
            raise ConfigurationError("regions must be a region string such as chr1:1-1000.")
        args += ["-r", regions]
    args.append(str(source))
    try:
        cnv(*args, catch_stdout=False)
    except pysam.SamtoolsError as exc:
        raise _samtools_error(exc) from exc
    outputs = {
        "cn": directory / f"cn.{sample}.tab",
        "dat": directory / f"dat.{sample}.tab",
        "summary": directory / f"summary.{sample}.tab",
        "plot_script": directory / f"plot.{sample}.py",
    }
    missing = [name for name in ("cn", "summary") if not outputs[name].is_file()]
    if missing:
        raise BackendExecutionError(f"bcftools cnv did not write the expected files: {missing}")
    sites = []
    for row in _parse_tab(outputs["cn"]):
        sites.append(
            {
                "chromosome": row[0],
                "position": int(row[1]),
                "copy_number": int(row[2]),
                "probabilities": [float(value) for value in row[3:7]],
            }
        )
    segments = []
    for row in _parse_tab(outputs["summary"]):
        if row[0] != "RG":
            continue
        segments.append(
            {
                "chromosome": row[1],
                "start": int(row[2]),
                "end": int(row[3]),
                "copy_number": int(row[4]),
                "quality": float(row[5]),
                "n_sites": int(row[6]),
                "n_hets": int(row[7]),
            }
        )
    return {
        "status": "completed",
        "backend": "pysam",
        "command": ["bcftools", "cnv", *args],
        "sample": sample,
        "outputs": {name: str(path) for name, path in outputs.items() if path.is_file()},
        "segments": segments,
        "sites": sites,
    }


__all__ = [
    "call_copy_number_states",
    "depad_alignments",
    "pad_to_unpad_alignments",
]
_register_public_verbose(globals())
