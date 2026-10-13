"""Optional pysam adapters for tabix/bgzip files, genomic regions and bundled bcftools.

- ``index_tabix_file()`` and ``compress_bgzip_file()`` use pysam's htslib tabix and bgzip code.
- ``parse_genomic_region()``, ``parse_genomic_regions()`` and ``fetch_padded_sequence()`` use
  pysam's region parsers and its FASTA reader.
- The other functions call the bcftools build that pysam ships as ``pysam.bcftools``. Each takes
  explicit arguments; none accepts a free-form command line. Output files are only written to new
  paths, unless ``force=True``.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``bio-reads`` install hint.
"""

from __future__ import annotations

import importlib
import importlib.metadata
import re
import tempfile
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require, workspace
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-reads"
_TABIX_PRESETS = ("vcf", "bed", "gff")
_INDEX_SUFFIXES = (".tbi", ".csi")
_FIELD_PATTERN = re.compile(r"^(?:CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO/[A-Za-z_][A-Za-z0-9_.]*)$")
_TAG_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_.]*$")
_REGION_PATTERN = re.compile(r"^(?P<chrom>[^:\s,]+)(?::(?P<start>\d[\d,]*)-(?P<end>\d[\d,]*))?$")
_FIXED_COLUMNS = frozenset({"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"})
_DEFAULT_MAX_RECORDS = 100_000
_MAX_RECORDS_LIMIT = 1_000_000
_MAX_WINDOW = 10_000_000
_MAX_REGIONS = 10_000
_MAX_FIELDS = 64
_MAX_HEADER_LINES = 256


def _pysam() -> Any:
    return require("pysam", _EXTRA)


def _submodule(name: str) -> Any:
    require("pysam", _EXTRA)
    return importlib.import_module(name)


def _provenance(convention: str) -> dict[str, str]:
    return {
        "package": "pysam",
        "version": importlib.metadata.version("pysam"),
        "convention": convention,
    }


def _text(value: object, name: str) -> str:
    if isinstance(value, Path):
        value = str(value)
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n"):
        raise ConfigurationError(f"{name} must be a non-empty single-line string.")
    return value.strip()


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return str(value)


def _sequence(value: object, name: str, *, limit: int) -> list[object]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list.")
    if len(value) > limit:
        raise ConfigurationError(f"{name} exceeds the limit of {limit} items.")
    return list(value)


def _names(value: object, name: str) -> list[str]:
    names = [_text(item, name) for item in _sequence(value, name, limit=100_000)]
    if len(set(names)) != len(names) or any("," in item for item in names):
        raise ConfigurationError(f"{name} must be unique and must not contain commas.")
    return names


def _input(value: object, name: str) -> Path:
    return input_path(_text(value, name))


def _output(value: object, name: str, *, force: bool, suffix: str | None = None) -> Path:
    path = Path(_text(value, name)).expanduser().resolve()
    if suffix is not None and not path.name.endswith(suffix):
        raise ConfigurationError(f"{name} must end with {suffix}.")
    if not path.parent.is_dir():
        raise ConfigurationError(f"{name} directory does not exist: {path.parent}")
    _flag(force, "force")
    if path.exists() and not force:
        raise ConfigurationError(f"{name} already exists: {path}; pass force=True to replace it.")
    return path


def _indexed(path: Path) -> None:
    if not any(Path(str(path) + suffix).is_file() for suffix in _INDEX_SUFFIXES):
        raise ConfigurationError(
            f"{path.name} needs a tabix (.tbi) or CSI (.csi) index; run index_tabix_file() first."
        )


def _bgzipped(path: Path) -> None:
    if path.suffix != ".gz":
        raise ConfigurationError(
            f"{path.name} must be bgzip-compressed (.gz); run index_tabix_file() first."
        )


def _call(pysam: Any, label: str, function: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except (pysam.SamtoolsError, OSError, RuntimeError, TypeError, ValueError, IndexError) as exc:
        raise ConfigurationError(f"{label} failed: {exc}") from exc


def _write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _records(text: str) -> int:
    return sum(1 for line in text.splitlines() if line and not line.startswith("#"))


def _position(text: str) -> int:
    value = int(text.replace(",", ""))
    if value < 1:
        raise ConfigurationError("region positions are 1-based and must be at least 1.")
    return value


def _region(value: object) -> tuple[str, int | None, int | None]:
    match = _REGION_PATTERN.match(value.strip()) if isinstance(value, str) else None
    if match is None:
        raise ConfigurationError(
            "region must look like 'chr1' or 'chr1:1000-2000' (1-based, inclusive); "
            "open-ended regions such as 'chr1:1000' are not accepted."
        )
    if match["start"] is None:
        return match["chrom"], None, None
    start = _position(match["start"])
    end = _position(match["end"])
    if end < start:
        raise ConfigurationError("region end must be greater than or equal to its start.")
    return match["chrom"], start, end


def _region_text(chrom: str, start: int | None, end: int | None) -> str:
    return chrom if start is None else f"{chrom}:{start}-{end}"


def _fields(value: object) -> list[str]:
    names = [_text(item, "fields") for item in _sequence(value, "fields", limit=_MAX_FIELDS)]
    for name in names:
        if not _FIELD_PATTERN.match(name):
            raise ConfigurationError(
                f"unsupported field {name!r}; use CHROM, POS, ID, REF, ALT, QUAL, FILTER "
                "or INFO/<tag>."
            )
    if len(set(names)) != len(names):
        raise ConfigurationError("fields must not repeat.")
    return names


def _columns(value: object) -> list[str]:
    names = [_text(item, "columns") for item in _sequence(value, "columns", limit=_MAX_FIELDS)]
    for name in names:
        if name in _FIXED_COLUMNS:
            continue
        if not (name.startswith("INFO/") and _TAG_PATTERN.match(name[5:])):
            raise ConfigurationError(f"unsupported annotation column {name!r}.")
    if len(set(names)) != len(names):
        raise ConfigurationError("columns must not repeat.")
    if "CHROM" not in names or "POS" not in names:
        raise ConfigurationError("columns must include CHROM and POS.")
    if not any(name.startswith("INFO/") for name in names):
        raise ConfigurationError("columns must include at least one INFO/<tag>.")
    return names


def _header_lines(value: object) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(
            "header_lines must be a list of ##INFO, ##FILTER or ##FORMAT lines."
        )
    if len(value) > _MAX_HEADER_LINES:
        raise ConfigurationError(f"header_lines exceed the limit of {_MAX_HEADER_LINES} items.")
    lines = [_text(item, "header_lines") for item in value]
    for line in lines:
        if not line.startswith("##") or "=" not in line:
            raise ConfigurationError("each header line must start with '##' and contain '='.")
    return lines


def _header_samples(header: str) -> list[str]:
    for line in header.splitlines():
        if line.startswith("#CHROM"):
            return line.split("\t")[9:]
    raise ConfigurationError("the VCF header has no #CHROM line.")


def _fasta_lengths(text: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current: str | None = None
    for line in text.splitlines():
        if line.startswith(">"):
            current = line[1:].strip()
            lengths[current] = 0
        elif current is not None:
            lengths[current] += len(line.strip())
    return lengths


def index_tabix_file(
    source: str | Path,
    *,
    preset: str = "vcf",
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bgzip (when plain) and tabix-index a VCF, BED or GFF file for region queries.

    A plain input is compressed to ``<name>.gz`` beside it and the original is kept.
    """

    path = _input(source, "source")
    selected = _choice(preset, "preset", _TABIX_PRESETS)
    _flag(force, "force")
    pysam = _pysam()
    compressed = path if path.suffix == ".gz" else path.with_name(path.name + ".gz")
    index = Path(str(compressed) + ".tbi")
    if not force:
        if index.exists():
            raise ConfigurationError(f"index already exists: {index}; pass force=True to rebuild.")
        if compressed != path and compressed.exists():
            raise ConfigurationError(f"{compressed} already exists; pass force=True to replace.")
    _call(
        pysam,
        "tabix_index",
        pysam.tabix_index,
        str(path),
        preset=selected,
        force=force,
        keep_original=True,
    )
    return {
        "status": "completed",
        "operation": "index_tabix_file",
        "preset": selected,
        "compressed": str(compressed),
        "index": str(index),
        "provenance": _provenance(
            "tabix index; VCF POS and BED/GFF starts are 1-based in the file"
        ),
    }


def compress_bgzip_file(
    source: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a BGZF-compressed copy of a text file with pysam's htslib bgzip."""

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force, suffix=".gz")
    pysam = _pysam()
    _call(pysam, "tabix_compress", pysam.tabix_compress, str(path), str(target), force=force)
    return {
        "status": "completed",
        "operation": "compress_bgzip_file",
        "output": str(target),
        "uncompressed_bytes": path.stat().st_size,
        "compressed_bytes": target.stat().st_size,
        "provenance": _provenance("BGZF blocks; readable by htslib, tabix and bcftools"),
    }


def parse_genomic_region(region: str, *, verbose: bool = True) -> dict[str, Any]:
    """Convert ``chr1:start-end`` (1-based, inclusive) to a 0-based half-open interval."""

    chrom, start, end = _region(region)
    pysam = _pysam()
    libcutils = _submodule("pysam.libcutils")
    reference, begin, stop = _call(
        pysam,
        "parse_region",
        libcutils.parse_region,
        region=_region_text(chrom, start, end),
    )
    whole = start is None
    return {
        "status": "completed",
        "operation": "parse_genomic_region",
        "region": region.strip(),
        "reference": reference,
        "start": int(begin),
        "end": None if whole else int(stop),
        "length": None if whole else int(stop) - int(begin),
        "provenance": _provenance(
            "start is 0-based and end is exclusive; a bare contig gives end None (whole contig)"
        ),
    }


def parse_genomic_regions(regions: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse comma-separated regions such as ``chr1:10-20,chr2`` into 0-based intervals."""

    pieces = [piece.strip() for piece in _text(regions, "regions").split(",")]
    if len(pieces) > _MAX_REGIONS:
        raise ConfigurationError(f"regions exceed the limit of {_MAX_REGIONS} items.")
    parsed = [_region(piece) for piece in pieces]
    pysam = _pysam()
    libcvcf = _submodule("pysam.libcvcf")
    normalized = ",".join(_region_text(chrom, start, end) for chrom, start, end in parsed)
    upstream = _call(pysam, "parse_regions", libcvcf.parse_regions, normalized)
    intervals = []
    for (_chrom, _start, end), (reference, begin, stop) in zip(parsed, upstream, strict=True):
        intervals.append(
            {
                "reference": reference,
                "start": int(begin),
                "end": None if end is None else int(stop),
            }
        )
    return {
        "status": "completed",
        "operation": "parse_genomic_regions",
        "regions": intervals,
        "count": len(intervals),
        "provenance": _provenance("0-based start and exclusive end; a bare contig gives end None"),
    }


def fetch_padded_sequence(
    fasta: str | Path,
    reference: str,
    start: int,
    end: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return upper-case bases for a 0-based half-open window; off-contig positions become N."""

    path = _input(fasta, "fasta")
    if not Path(str(path) + ".fai").is_file():
        raise ConfigurationError("Reference index missing; call index_reference first.")
    name = _text(reference, "reference")
    begin = _integer(start, "start", minimum=-_MAX_WINDOW, maximum=_MAX_WINDOW)
    stop = _integer(end, "end", minimum=-_MAX_WINDOW, maximum=2_000_000_000)
    if stop <= begin or stop - begin > _MAX_WINDOW:
        raise ConfigurationError(f"window must be non-empty and at most {_MAX_WINDOW} bases.")
    pysam = _pysam()
    libcvcf = _submodule("pysam.libcvcf")
    handle = _call(pysam, "FastaFile", pysam.FastaFile, str(path))
    try:
        if name not in handle.references:
            raise ConfigurationError(f"reference {name!r} is not in {path.name}.")
        length = int(handle.get_reference_length(name))
        sequence = _call(pysam, "get_sequence", libcvcf.get_sequence, name, begin, stop, handle)
    finally:
        handle.close()
    return {
        "status": "completed",
        "operation": "fetch_padded_sequence",
        "reference": name,
        "start": begin,
        "end": stop,
        "reference_length": length,
        "left_padding": max(0, -begin),
        "right_padding": max(0, stop - max(length, begin, 0)),
        "sequence": sequence,
        "provenance": _provenance(
            "0-based half-open window; positions outside the contig are N; bases upper-cased"
        ),
    }


def query_vcf_fields(
    source: str | Path,
    fields: Sequence[str],
    *,
    region: str | None = None,
    genotypes: bool = False,
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract VCF columns as text rows with ``bcftools query``.

    ``fields`` are CHROM, POS, ID, REF, ALT, QUAL, FILTER or ``INFO/<tag>``. With
    ``genotypes=True`` each sample adds one ``GT:<sample>`` column. Values are returned as
    text, and missing values stay ``"."``.
    """

    path = _input(source, "source")
    names = _fields(fields)
    _flag(genotypes, "genotypes")
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    pysam = _pysam()
    bcftools = _submodule("pysam.bcftools")
    args: list[str] = []
    if region is not None:
        chrom, start, end = _region(region)
        _indexed(path)
        args = ["-r", _region_text(chrom, start, end)]
    columns = list(names)
    template = "\t".join(f"%{name}" for name in names)
    if genotypes:
        header = _call(pysam, "bcftools view -h", bcftools.view, "-h", str(path))
        columns += [f"GT:{sample}" for sample in _header_samples(header)]
        template += "[\t%GT]"
    text = _call(pysam, "bcftools query", bcftools.query, "-f", template + "\n", *args, str(path))
    rows = [line.split("\t") for line in text.splitlines() if line]
    if len(rows) > limit:
        raise ConfigurationError(f"query returned more than max_records={limit} records.")
    if any(len(row) != len(columns) for row in rows):
        raise ConfigurationError("bcftools query returned an unexpected number of columns.")
    return {
        "status": "completed",
        "operation": "query_vcf_fields",
        "columns": columns,
        "rows": rows,
        "record_count": len(rows),
        "region": region,
        "provenance": _provenance("bcftools query text; no genotype or INFO value is converted"),
    }


def select_vcf_records(
    source: str | Path,
    output_path: str | Path,
    *,
    samples: Sequence[str] | None = None,
    region: str | None = None,
    pass_only: bool = False,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a VCF subset chosen by sample names, a region and optional PASS-only filtering.

    Samples keep the order given. A region needs an indexed input.
    """

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    _flag(pass_only, "pass_only")
    names = None if samples is None else _names(samples, "samples")
    args: list[str] = []
    if names is not None:
        args += ["-s", ",".join(names)]
    if region is not None:
        chrom, start, end = _region(region)
        _indexed(path)
        args += ["-r", _region_text(chrom, start, end)]
    if pass_only:
        args += ["-f", "PASS"]
    pysam = _pysam()
    bcftools = _submodule("pysam.bcftools")
    text = _call(pysam, "bcftools view", bcftools.view, *args, str(path))
    _write(target, text)
    return {
        "status": "completed",
        "operation": "select_vcf_records",
        "output": str(target),
        "records": _records(text),
        "samples": names,
        "region": region,
        "pass_only": pass_only,
        "provenance": _provenance("uncompressed VCF output; FILTER PASS keeps records marked PASS"),
    }


def concatenate_vcf_files(
    sources: Sequence[str | Path],
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Concatenate indexed, bgzipped VCF files that share samples, in the order given."""

    paths = [_input(item, "sources") for item in _sequence(sources, "sources", limit=10_000)]
    if len(paths) < 2:
        raise ConfigurationError("sources must list at least two VCF files.")
    for path in paths:
        _bgzipped(path)
        _indexed(path)
    target = _output(output_path, "output_path", force=force)
    pysam = _pysam()
    bcftools = _submodule("pysam.bcftools")
    text = _call(pysam, "bcftools concat", bcftools.concat, *[str(path) for path in paths])
    _write(target, text)
    return {
        "status": "completed",
        "operation": "concatenate_vcf_files",
        "output": str(target),
        "inputs": len(paths),
        "records": _records(text),
        "provenance": _provenance("bcftools concat; files are joined in the order given"),
    }


def intersect_vcf_sites(
    first: str | Path,
    second: str | Path,
    output_dir: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count sites private to each of two indexed VCF files and shared by both.

    Sites match on chromosome, position and alleles. ``output_dir`` must not exist; bcftools
    writes its VCF files and ``sites.txt`` there.
    """

    left = _input(first, "first")
    right = _input(second, "second")
    for path in (left, right):
        _bgzipped(path)
        _indexed(path)
    directory = workspace(_text(output_dir, "output_dir"))
    pysam = _pysam()
    bcftools = _submodule("pysam.bcftools")
    _call(pysam, "bcftools isec", bcftools.isec, "-p", str(directory), str(left), str(right))
    counts = {"10": 0, "01": 0, "11": 0}
    for line in (directory / "sites.txt").read_text(encoding="utf-8").splitlines():
        code = line.split("\t")[-1]
        if code not in counts:
            raise ConfigurationError(f"unexpected bcftools isec site code {code!r}.")
        counts[code] += 1
    return {
        "status": "completed",
        "operation": "intersect_vcf_sites",
        "directory": str(directory),
        "first": str(left),
        "second": str(right),
        "private_first": counts["10"],
        "private_second": counts["01"],
        "shared": counts["11"],
        "files": sorted(item.name for item in directory.iterdir()),
        "provenance": _provenance("sites match on CHROM, POS and REF/ALT; each site counted once"),
    }


def sort_vcf_records(
    source: str | Path,
    output_path: str | Path,
    *,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort VCF records by chromosome order and position with ``bcftools sort``."""

    path = _input(source, "source")
    target = _output(output_path, "output_path", force=force)
    pysam = _pysam()
    bcftools = _submodule("pysam.bcftools")
    text = _call(pysam, "bcftools sort", bcftools.sort, str(path))
    _write(target, text)
    return {
        "status": "completed",
        "operation": "sort_vcf_records",
        "output": str(target),
        "records": _records(text),
        "provenance": _provenance("uncompressed VCF output; contigs follow the input header order"),
    }


def annotate_vcf_from_tabix(
    source: str | Path,
    annotations: str | Path,
    columns: Sequence[str],
    output_path: str | Path,
    *,
    header_lines: Sequence[str] = (),
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Copy INFO tags from a bgzipped, tabix-indexed annotation table onto matching records.

    ``columns`` names the annotation-file columns in order, for example
    ``("CHROM", "POS", "REF", "ALT", "INFO/NOTE")``. CHROM and POS are required. Each new tag
    needs a matching ``##INFO`` line in ``header_lines``.
    """

    path = _input(source, "source")
    annotation = _input(annotations, "annotations")
    _bgzipped(annotation)
    _indexed(annotation)
    cols = _columns(columns)
    headers = _header_lines(header_lines)
    target = _output(output_path, "output_path", force=force)
    pysam = _pysam()
    bcftools = _submodule("pysam.bcftools")
    with tempfile.TemporaryDirectory() as scratch:
        args = ["-a", str(annotation), "-c", ",".join(cols)]
        if headers:
            header_file = Path(scratch) / "header.txt"
            _write(header_file, "\n".join(headers) + "\n")
            args += ["-h", str(header_file)]
        text = _call(pysam, "bcftools annotate", bcftools.annotate, *args, str(path))
    _write(target, text)
    tags = sorted(column[len("INFO/") :] for column in cols if column.startswith("INFO/"))
    annotated = dict.fromkeys(tags, 0)
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        keys = {item.split("=", 1)[0] for item in line.split("\t")[7].split(";")}
        for tag in tags:
            if tag in keys:
                annotated[tag] += 1
    return {
        "status": "completed",
        "operation": "annotate_vcf_from_tabix",
        "output": str(target),
        "records": _records(text),
        "annotated": annotated,
        "provenance": _provenance("bcftools annotate with exact CHROM, POS, REF and ALT matches"),
    }


def build_vcf_consensus(
    source: str | Path,
    fasta: str | Path,
    sample: str,
    output_path: str | Path,
    *,
    haplotype: int = 1,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply one sample's variants on one haplotype to a FASTA reference with bcftools.

    The VCF must be bgzip-compressed and indexed, and the FASTA needs a ``.fai`` index.
    ``haplotype`` is 1 or 2; heterozygous calls use the allele at that position of the GT field.
    """

    path = _input(source, "source")
    _bgzipped(path)
    _indexed(path)
    reference = _input(fasta, "fasta")
    if not Path(str(reference) + ".fai").is_file():
        raise ConfigurationError("Reference index missing; call index_reference first.")
    name = _names([sample], "sample")[0]
    selected = _integer(haplotype, "haplotype", minimum=1, maximum=2)
    target = _output(output_path, "output_path", force=force)
    pysam = _pysam()
    bcftools = _submodule("pysam.bcftools")
    text = _call(
        pysam,
        "bcftools consensus",
        bcftools.consensus,
        "-f",
        str(reference),
        "-s",
        name,
        "-H",
        str(selected),
        str(path),
    )
    _write(target, text)
    return {
        "status": "completed",
        "operation": "build_vcf_consensus",
        "output": str(target),
        "sample": name,
        "haplotype": selected,
        "sequences": _fasta_lengths(text),
        "provenance": _provenance("bcftools consensus; REF-only sites are copied from the FASTA"),
    }


__all__ = [
    "annotate_vcf_from_tabix",
    "build_vcf_consensus",
    "compress_bgzip_file",
    "concatenate_vcf_files",
    "fetch_padded_sequence",
    "index_tabix_file",
    "intersect_vcf_sites",
    "parse_genomic_region",
    "parse_genomic_regions",
    "query_vcf_fields",
    "select_vcf_records",
    "sort_vcf_records",
]

_register_public_verbose(globals(), tuple(__all__))
