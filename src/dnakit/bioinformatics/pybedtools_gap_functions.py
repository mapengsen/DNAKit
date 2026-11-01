"""pybedtools operations as DNAKit wrappers: chromosome sizes, UCSC bigBed/bigWig facades, samtools
read counts, BAM-to-bigWig conversion, binary presence heat maps, Venn-diagram helpers and
feature colouring.

Chromosome-size wrappers call pybedtools' lookup functions (its built-in registry, genomepy or UCSC)
and record the service, the request, a UTC retrieval time and the sizes. The UCSC and samtools
programs are the user's to install; the facades take explicit file paths, fixed flags and new output
paths, and never run through a shell. Interval operations (multi-intersection, Venn set
construction, BAM coverage) run through pybedtools and therefore through the bedtools program that
pybedtools finds. Plot functions write one PNG to a new path and draw on a matplotlib Figure object,
so no pyplot global state is used.

Missing optional packages raise ``BackendUnavailableError`` with the ``pip install "dnakit[...]"``
hint. Ordinary tests never use the network.
"""

from __future__ import annotations

import datetime as dt
import math
import re
import shutil
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_PYBEDTOOLS_EXTRA = "bio-pybedtools"
_GENOMEPY_EXTRA = "bio-genomepy"
_PLOT_EXTRA = "bio-plots"
_GENOME = re.compile(r"[A-Za-z0-9_.-]{1,64}")
_HOST = re.compile(r"[A-Za-z0-9.-]{1,253}")
_BED_FORMAT = re.compile(r"bed[1-9][0-9]*(\+[0-9]*)?")
_R_NAME = re.compile(r"[A-Za-z0-9_.-][A-Za-z0-9_. -]{0,39}")
_R_ARGUMENT = re.compile(r"[A-Za-z][A-Za-z0-9_.]*\s*=[^\r\n;`]*")
_MAX_TAIL = 4000
_MAX_ROWS = 100_000
_MAX_CHROMOSOMES = 10_000
_SERVICE_REGISTRY = "pybedtools built-in genome registry (offline)"
_SERVICE_LOOKUP = "pybedtools chromsizes lookup: genomepy first, then UCSC over the network"
_SERVICE_UCSC = "UCSC Genome Browser MySQL server (through the mysql and fetchChromSizes programs)"
_SERVICE_GENOMEPY = "genomepy local genome installation"
_HINT_BEDTOBIGBED = (
    "Install the UCSC bedToBigBed program (conda install -c bioconda ucsc-bedtobigbed, "
    "or download it from "
    "https://hgdownload.soe.ucsc.edu/admin/exe/) and put it on PATH."
)
_HINT_BIGBEDTOBED = (
    "Install the UCSC bigBedToBed program (conda install -c bioconda ucsc-bigbedtobed, "
    "or download it from "
    "https://hgdownload.soe.ucsc.edu/admin/exe/) and put it on PATH."
)
_HINT_BEDGRAPHTOBIGWIG = (
    "Install the UCSC bedGraphToBigWig program "
    "(conda install -c bioconda ucsc-bedgraphtobigwig, or download it "
    "from https://hgdownload.soe.ucsc.edu/admin/exe/) and put it on PATH."
)
_HINT_BIGWIGTOBEDGRAPH = (
    "Install the UCSC bigWigToBedGraph program "
    "(conda install -c bioconda ucsc-bigwigtobedgraph, or download it "
    "from https://hgdownload.soe.ucsc.edu/admin/exe/) and put it on PATH."
)
_HINT_WIGTOBIGWIG = (
    "Install the UCSC wigToBigWig program (conda install -c bioconda ucsc-wigtobigwig, "
    "or download it from "
    "https://hgdownload.soe.ucsc.edu/admin/exe/) and put it on PATH."
)
_HINT_SAMTOOLS = (
    "Install samtools (conda install -c bioconda samtools, "
    "or apt install samtools) and put it on PATH."
)
_HINT_BEDTOOLS = (
    "Install bedtools (conda install -c bioconda bedtools, "
    "or apt install bedtools) and put it on PATH."
)
_HINT_RSCRIPT = (
    "Install R with the VennDiagram package (conda install -c conda-forge r-base, then "
    "install.packages('VennDiagram') in R) and put Rscript on PATH."
)


def _utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _pybedtools_version() -> str:
    try:
        return version("pybedtools")
    except (
        PackageNotFoundError
    ):  # pragma: no cover - pybedtools is installed whenever its wrappers run.
        return "unknown"


def _plain(value: Any) -> Any:
    """Convert numpy scalars to JSON-friendly values; NaN and infinity become None."""

    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _genome(value: object) -> str:
    if not isinstance(value, str) or not _GENOME.fullmatch(value):
        raise ConfigurationError(
            "genome must be an assembly name such as hg38 (letters, digits, '_', '.' or '-')."
        )
    return value


def _host(value: object) -> str:
    if not isinstance(value, str) or not _HOST.fullmatch(value):
        raise ConfigurationError("host_url must be a host name such as genome-mysql.cse.ucsc.edu.")
    return value


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _positive(value: object, name: str, maximum: int = 1_000_000_000) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


def _nonnegative(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError(f"{name} must be a non-negative integer.")
    return value


def _existing(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _existing_list(values: object, name: str, minimum: int, maximum: int) -> list[Path]:
    if isinstance(values, (str, bytes, Path)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of file paths.")
    if not minimum <= len(values) <= maximum:
        raise ConfigurationError(f"{name} must list {minimum} to {maximum} files.")
    return [_existing(item, name) for item in values]


def _new_path(value: object, name: str) -> Path:
    """Output path that must not exist yet; its parent directory must exist."""

    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory for {name} does not exist: {path.parent}")
    return path


def _optional_output(value: object, name: str) -> Path | None:
    return None if value is None else _new_path(value, name)


def _program(name: str, explicit: str | None, hint: str) -> str:
    if explicit is not None:
        candidate = Path(_text(explicit, "program path")).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"{name} path does not exist: {candidate}", hint=hint)
        return str(candidate)
    found = shutil.which(name)
    if found is None:
        raise BackendUnavailableError(f"{name} was not found on PATH.", hint=hint)
    return found


def _run_tool(
    name: str, args: list[str], output: Path, explicit: str | None, hint: str
) -> dict[str, Any]:
    executable = _program(name, explicit, hint)
    result = execute_bounded_command(
        executable,
        args,
        backend_id=f"ucsc-{name.lower()}",
        cwd=output.parent,
        verbose=False,
    )
    produced = output.is_file()
    return {
        "status": "completed" if result.return_code == 0 and produced else "failed",
        "backend": name,
        "executable": executable,
        "command": [name, *args],
        "return_code": result.return_code,
        "outputs": [str(output)] if produced else [],
        "log_tail": (result.output or "")[-_MAX_TAIL:],
    }


def _sizes_result(
    sizes: Any, service: str, genome: str, target: Path | None, pybedtools: Any
) -> dict[str, Any]:
    """Build the JSON result of a chromosome-size lookup and optionally write a chrom.sizes file."""

    data = {str(chrom): int(bounds[1]) for chrom, bounds in sizes.items()}
    if len(data) > _MAX_CHROMOSOMES:
        raise ConfigurationError(f"{genome} lists more than {_MAX_CHROMOSOMES} sequences.")
    if target is not None:
        pybedtools.chromsizes_to_file(
            {chrom: (0, length) for chrom, length in data.items()}, str(target)
        )
    return {
        "status": "completed",
        "service": service,
        "request": {"genome": genome},
        "retrieved_at": _utc_now(),
        "pybedtools_version": _pybedtools_version(),
        "row_count": len(data),
        "data": data,
        "output_path": None if target is None else str(target),
    }


def chromsizes(
    genome: str,
    *,
    output_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the chromosome sizes of a genome from pybedtools' registry, then genomepy or UCSC
    (pybedtools chromsizes).
    """

    name = _genome(genome)
    target = _optional_output(output_path, "output_path")
    pybedtools = require("pybedtools", _PYBEDTOOLS_EXTRA)
    offline = hasattr(pybedtools.genome_registry, name)
    try:
        sizes = pybedtools.chromsizes(name)
    except Exception as exc:  # genomepy, mysql and urllib raise unrelated exception types
        raise BackendExecutionError(
            f"pybedtools could not get chromosome sizes for {name}: {type(exc).__name__}: {exc}"
        ) from exc
    service = _SERVICE_REGISTRY if offline else _SERVICE_LOOKUP
    return _sizes_result(sizes, service, name, target, pybedtools)


def get_chromsizes_from_ucsc(
    genome: str,
    *,
    host_url: str = "genome-mysql.cse.ucsc.edu",
    timeout: float | None = None,
    output_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download the chromosome sizes of a UCSC assembly from the UCSC MySQL server (pybedtools
    get_chromsizes_from_ucsc).
    """

    name = _genome(genome)
    host = _host(host_url)
    if timeout is not None and (
        isinstance(timeout, bool)
        or not isinstance(timeout, (int, float))
        or not 0 < timeout <= 3600
    ):
        raise ConfigurationError("timeout must be a number of seconds in (0, 3600], or None.")
    target = _optional_output(output_path, "output_path")
    pybedtools = require("pybedtools", _PYBEDTOOLS_EXTRA)
    try:
        sizes = pybedtools.get_chromsizes_from_ucsc(
            name,
            saveas=None if target is None else str(target),
            timeout=timeout,
            host_url=host,
        )
    except (OSError, ValueError) as exc:
        raise BackendExecutionError(
            f"pybedtools could not download chromosome sizes for {name} from UCSC: {exc}"
        ) from exc
    result = _sizes_result(sizes, _SERVICE_UCSC, name, None, pybedtools)
    result["request"]["host_url"] = host
    result["output_path"] = None if target is None else str(target)
    return result


def get_chromsizes_from_genomepy(
    genome: str,
    *,
    output_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the chromosome sizes of a genome installed with genomepy (pybedtools
    get_chromsizes_from_genomepy).
    """

    name = _genome(genome)
    target = _optional_output(output_path, "output_path")
    require("genomepy", _GENOMEPY_EXTRA)
    pybedtools = require("pybedtools", _PYBEDTOOLS_EXTRA)
    try:
        sizes = pybedtools.helpers.get_chromsizes_from_genomepy(
            name, saveas=None if target is None else str(target)
        )
    except Exception as exc:  # genomepy raises its own exception types for unknown genomes
        raise BackendExecutionError(
            f"genomepy could not give chromosome sizes for {name}: {type(exc).__name__}: {exc}"
        ) from exc
    if sizes is None:
        return {
            "status": "not_found",
            "service": _SERVICE_GENOMEPY,
            "request": {"genome": name},
            "retrieved_at": _utc_now(),
            "pybedtools_version": _pybedtools_version(),
            "row_count": 0,
            "data": None,
            "output_path": None,
            "message": (
                "genomepy has no local sizes file for this genome; "
                "install it with genomepy install."
            ),
        }
    return _sizes_result(sizes, _SERVICE_GENOMEPY, name, None, pybedtools) | {
        "output_path": None if target is None else str(target)
    }


def _field_count(path: Path) -> int:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            return len(line.rstrip("\r\n").split("\t"))
    raise ConfigurationError(f"No BED records were found in {path}.")


def bigbed(
    bed: str,
    chromsizes: str,
    output_bigbed: str,
    *,
    bed_type: str | None = None,
    block_size: int = 256,
    items_per_slot: int = 512,
    autosql: str | None = None,
    unc: bool = False,
    tab: bool = False,
    bedtobigbed_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a BED file to an indexed bigBed file with UCSC bedToBigBed (pybedtools
    contrib.bigbed.bigbed).
    """

    bed_file = _existing(bed, "bed")
    sizes_file = _existing(chromsizes, "chromsizes")
    out = _new_path(output_bigbed, "output_bigbed")
    block = _positive(block_size, "block_size", 1_000_000)
    items = _positive(items_per_slot, "items_per_slot", 1_000_000)
    if bed_type is None:
        kind = f"bed{_field_count(bed_file)}"
    else:
        if not isinstance(bed_type, str) or not _BED_FORMAT.fullmatch(bed_type):
            raise ConfigurationError(
                "bed_type must look like bed3, bed6 or bed6+4 (bedToBigBed -type)."
            )
        kind = bed_type
    args = [
        str(bed_file),
        str(sizes_file),
        str(out),
        f"-blockSize={block}",
        f"-itemsPerSlot={items}",
        f"-type={kind}",
    ]
    if _flag(unc, "unc"):
        args.append("-unc")
    if _flag(tab, "tab"):
        args.append("-tab")
    if autosql is not None:
        args.append(f"-as={_existing(autosql, 'autosql')}")
    return _run_tool("bedToBigBed", args, out, bedtobigbed_path, _HINT_BEDTOBIGBED)


def _region_options(chrom: str | None, start: int | None, end: int | None, name: str) -> list[str]:
    if chrom is None:
        if start is not None or end is not None:
            raise ConfigurationError(f"{name} needs chrom when start or end is given.")
        return []
    options = [f"-chrom={_text(chrom, 'chrom')}"]
    if start is not None:
        options.append(f"-start={_nonnegative(start, 'start')}")
    if end is not None:
        options.append(f"-end={_nonnegative(end, 'end')}")
    return options


def bigbed_to_bed(
    bigbed: str,
    output_bed: str,
    *,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    max_items: int | None = None,
    bigbedtobed_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the intervals of an indexed bigBed file to a BED file with UCSC bigBedToBed
    (pybedtools contrib.bigbed.bigbed_to_bed).
    """

    source = _existing(bigbed, "bigbed")
    out = _new_path(output_bed, "output_bed")
    args = _region_options(chrom, start, end, "bigbed_to_bed")
    if max_items is not None:
        args.append(f"-maxItems={_positive(max_items, 'max_items')}")
    args += [str(source), str(out)]
    return _run_tool("bigBedToBed", args, out, bigbedtobed_path, _HINT_BIGBEDTOBED)


def bedgraph_to_bigwig(
    bedgraph: str,
    chromsizes: str,
    output_bigwig: str,
    *,
    bedgraphtobigwig_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a sorted bedGraph file to an indexed bigWig file with UCSC bedGraphToBigWig
    (pybedtools contrib.bigwig.bedgraph_to_bigwig).
    """

    graph_file = _existing(bedgraph, "bedgraph")
    sizes_file = _existing(chromsizes, "chromsizes")
    out = _new_path(output_bigwig, "output_bigwig")
    args = [str(graph_file), str(sizes_file), str(out)]
    return _run_tool("bedGraphToBigWig", args, out, bedgraphtobigwig_path, _HINT_BEDGRAPHTOBIGWIG)


def bigwig_to_bedgraph(
    bigwig: str,
    output_bedgraph: str,
    *,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    udc_dir: str | None = None,
    bigwigtobedgraph_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the values of an indexed bigWig file to a bedGraph file with UCSC bigWigToBedGraph
    (pybedtools contrib.bigwig.bigwig_to_bedgraph).
    """

    source = _existing(bigwig, "bigwig")
    out = _new_path(output_bedgraph, "output_bedgraph")
    args = _region_options(chrom, start, end, "bigwig_to_bedgraph")
    if udc_dir is not None:
        args.append(f"-udcDir={_text(udc_dir, 'udc_dir')}")
    args += [str(source), str(out)]
    return _run_tool("bigWigToBedGraph", args, out, bigwigtobedgraph_path, _HINT_BIGWIGTOBEDGRAPH)


def wig_to_bigwig(
    wig: str,
    chromsizes: str,
    output_bigwig: str,
    *,
    wigtobigwig_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a wiggle file to an indexed bigWig file with UCSC wigToBigWig (pybedtools
    contrib.bigwig.wig_to_bigwig).
    """

    wig_file = _existing(wig, "wig")
    sizes_file = _existing(chromsizes, "chromsizes")
    out = _new_path(output_bigwig, "output_bigwig")
    args = [str(wig_file), str(sizes_file), str(out)]
    return _run_tool("wigToBigWig", args, out, wigtobigwig_path, _HINT_WIGTOBIGWIG)


def mapped_read_count(
    bam: str,
    *,
    samtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the mapped reads of a BAM file with samtools view -c -F 0x4 (pybedtools
    contrib.bigwig.mapped_read_count).

    pybedtools caches the count in a ``<bam>.scale`` file beside the BAM; this facade does not read
    or write that cache, so the input folder is left unchanged.
    """

    source = _existing(bam, "bam")
    executable = _program("samtools", samtools_path, _HINT_SAMTOOLS)
    args = ["view", "-c", "-F", "0x4", str(source)]
    result = execute_bounded_command(
        executable, args, backend_id="samtools", cwd=source.parent, verbose=False
    )
    count: int | None = None
    if result.return_code == 0:
        lines = (result.output or "").strip().splitlines()
        try:
            count = int(lines[-1])
        except (ValueError, IndexError):
            count = None
    return {
        "status": "completed" if count is not None else "failed",
        "backend": "samtools",
        "executable": executable,
        "command": ["samtools", *args],
        "return_code": result.return_code,
        "mapped_reads": count,
        "outputs": [],
        "log_tail": (result.output or "")[-_MAX_TAIL:],
    }


def bam_to_bigwig(
    bam: str,
    chromsizes: str,
    output_bigwig: str,
    *,
    scale: bool = False,
    samtools_path: str | None = None,
    bedgraphtobigwig_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a bigWig of read coverage from a BAM file, optionally in reads per million mapped reads
    (pybedtools contrib.bigwig.bam_to_bigwig).

    Coverage comes from bedtools genomecov (-bg -split, as pybedtools does); bedGraphToBigWig then
    writes the bigWig. With scale=True the factor 1e6 / mapped reads comes from samtools. The
    intermediate bedGraph files are removed afterwards. The chromsizes argument is passed to
    bedGraphToBigWig; bedtools ignores -g for BAM input.
    """

    source = _existing(bam, "bam")
    sizes_file = _existing(chromsizes, "chromsizes")
    out = _new_path(output_bigwig, "output_bigwig")
    use_scale = _flag(scale, "scale")
    _program("bedtools", None, _HINT_BEDTOOLS)
    _program("bedGraphToBigWig", bedgraphtobigwig_path, _HINT_BEDGRAPHTOBIGWIG)
    pybedtools = require("pybedtools", _PYBEDTOOLS_EXTRA)
    coverage_options: dict[str, Any] = {"bg": True, "split": True}
    reads: int | None = None
    if use_scale:
        counted = mapped_read_count(str(source), samtools_path=samtools_path, verbose=False)
        if counted["status"] != "completed":
            raise BackendExecutionError(
                "samtools could not count the mapped reads of the BAM file."
            )
        reads = int(counted["mapped_reads"])
        if reads == 0:
            raise ConfigurationError(
                "The BAM file has no mapped reads, so no scale factor can be computed."
            )
        coverage_options["scale"] = 1 / (reads / 1e6)
    intermediate = out.parent / f"{out.name}.tmp-coverage.bedgraph"
    sorted_copy = out.parent / f"{out.name}.tmp-sorted.bedgraph"
    for extra in (intermediate, sorted_copy):
        if extra.exists():
            raise ConfigurationError(
                f"A temporary file already exists and will not be overwritten: {extra}"
            )
    try:
        try:
            pybedtools.BedTool(str(source)).genome_coverage(**coverage_options).saveas(
                str(intermediate)
            )
        except Exception as exc:  # pybedtools raises ValueError or OSError from bedtools failures
            raise BackendExecutionError(f"bedtools genomecov failed for {source}: {exc}") from exc
        result = _run_tool(
            "bedGraphToBigWig",
            [str(intermediate), str(sizes_file), str(out)],
            out,
            bedgraphtobigwig_path,
            _HINT_BEDGRAPHTOBIGWIG,
        )
        if result["status"] == "failed" and "bedSort" in result["log_tail"]:
            # The bedGraph is unsorted (BAM not coordinate-sorted): sort it with bedtools and retry,
            # as pybedtools does.
            pybedtools.BedTool(str(intermediate)).sort().saveas(str(sorted_copy))
            result = _run_tool(
                "bedGraphToBigWig",
                [str(sorted_copy), str(sizes_file), str(out)],
                out,
                bedgraphtobigwig_path,
                _HINT_BEDGRAPHTOBIGWIG,
            )
            result["sorted_before_conversion"] = True
    finally:
        for extra in (intermediate, sorted_copy):
            if extra.exists():
                extra.unlink()
    result["mapped_reads"] = reads
    result["scale_factor"] = coverage_options.get("scale")
    return result


def _figure_classes() -> tuple[Any, Any]:
    figure = require("matplotlib.figure", _PLOT_EXTRA)
    canvas = require("matplotlib.backends.backend_agg", _PLOT_EXTRA)
    return figure.Figure, canvas.FigureCanvasAgg


def _names(values: object, count: int) -> list[str]:
    if isinstance(values, str) or not isinstance(values, Sequence) or len(values) != count:
        raise ConfigurationError(f"names must list one name per file ({count} names).")
    labels: list[str] = []
    for item in values:
        if not isinstance(item, str) or not _R_NAME.fullmatch(item):
            raise ConfigurationError(
                "Each name may contain letters, digits, '_', '.', '-' and spaces only "
                "(at most 40 characters)."
            )
        labels.append(item)
    if len(set(labels)) != count:
        raise ConfigurationError("names must be unique.")
    return labels


def _binary_rows(matrix: object) -> list[list[int]]:
    if isinstance(matrix, (str, bytes)) or not isinstance(matrix, Sequence):
        raise ConfigurationError("matrix must be a sequence of rows of 0 and 1 values.")
    if len(matrix) > _MAX_ROWS:
        raise ConfigurationError(f"matrix can hold at most {_MAX_ROWS} rows.")
    rows: list[list[int]] = []
    width: int | None = None
    for row in matrix:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or not row:
            raise ConfigurationError("Every matrix row must be a non-empty sequence.")
        if width is None:
            width = len(row)
        if len(row) != width:
            raise ConfigurationError("Every matrix row must have the same length.")
        values: list[int] = []
        for value in row:
            if isinstance(value, bool) or value not in (0, 1) or not isinstance(value, int):
                raise ConfigurationError("The matrix may contain only the integers 0 and 1.")
            values.append(value)
        rows.append(values)
    return rows


def _row_order(numpy: Any, matrix: Any) -> list[int]:
    """Column-weighted row order used by pybedtools: highest weighted score first."""

    weights = [2**i for i in range(1, matrix.shape[1] + 1)[::-1]]
    score_matrix = matrix * weights
    return [int(index) for index in numpy.argsort(score_matrix.sum(axis=1))[::-1]]


def sort_binary_matrix(matrix: Sequence[Sequence[int]], *, verbose: bool = True) -> list[int]:
    """Return the row order of a 0/1 presence matrix sorted by column-weighted score, highest first
    (pybedtools sort_binary_matrix).
    """

    rows = _binary_rows(matrix)
    if not rows:
        return []
    numpy = require("numpy", _PYBEDTOOLS_EXTRA)
    return _row_order(numpy, numpy.asarray(rows, dtype=int))


def _sets_text(intervals: Any) -> list[str]:
    """Interval strings as in the Venn sets: tabs become '|' and the trailing newline is dropped."""

    return sorted({str(item).rstrip("\n").replace("\t", "|") for item in intervals})


def cleaned_intersect(beds: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Build the interval-string sets of 2 to 4 BED files whose overlapping intervals are made
    identical (pybedtools venn_maker.cleaned_intersect).
    """

    paths = _existing_list(beds, "beds", 2, 4)
    pybedtools = require("pybedtools", _PYBEDTOOLS_EXTRA)
    venn = require("pybedtools.contrib.venn_maker", _PYBEDTOOLS_EXTRA)
    cleaned = venn.cleaned_intersect([pybedtools.BedTool(str(path)) for path in paths])
    sets = [_sets_text(part) for part in cleaned]
    return {
        "status": "completed",
        "backend": "pybedtools + bedtools",
        "pybedtools_version": _pybedtools_version(),
        "files": [str(path) for path in paths],
        "sets": sets,
        "set_sizes": [len(part) for part in sets],
    }


def _region_counts(sets: Mapping[str, set[str]]) -> dict[str, int]:
    """Number of items in each exact combination of sets (for example 'a' or 'a&b'), as a Venn
    diagram shows them.
    """

    labels = list(sets)
    universe: set[str] = set().union(*sets.values())
    counts: dict[str, int] = {}
    for item in universe:
        key = "&".join(label for label in labels if item in sets[label])
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def _r_arguments(values: object) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str) or not isinstance(values, Sequence):
        raise ConfigurationError(
            "additional_args must be a sequence of R named arguments such as 'euler.d=TRUE'."
        )
    arguments: list[str] = []
    for item in values:
        if not isinstance(item, str) or not _R_ARGUMENT.fullmatch(item):
            raise ConfigurationError(
                "Each additional argument must be an R named argument such as "
                '\'cat.col=c("red","blue")\' '
                "(no line breaks or semicolons)."
            )
        arguments.append(item)
    return arguments


def venn_maker(
    beds: Sequence[str],
    output_script: str,
    *,
    names: Sequence[str] | None = None,
    figure_path: str | None = None,
    additional_args: Sequence[str] | None = None,
    render: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write an R script that draws a Venn diagram of 2 to 4 BED files with VennDiagram, and render
    it when Rscript is on PATH (pybedtools venn_maker).
    """

    paths = _existing_list(beds, "beds", 2, 4)
    labels = list("abcd"[: len(paths)]) if names is None else _names(names, len(paths))
    script = _new_path(output_script, "output_script")
    figure = None if figure_path is None else _new_path(figure_path, "figure_path")
    if figure is not None and any(c in str(figure) for c in '"\\\r\n'):
        raise ConfigurationError("figure_path may not contain quotes, backslashes or line breaks.")
    extras = _r_arguments(additional_args)
    if not isinstance(render, bool):
        raise ConfigurationError("render must be a boolean.")

    pybedtools = require("pybedtools", _PYBEDTOOLS_EXTRA)
    venn = require("pybedtools.contrib.venn_maker", _PYBEDTOOLS_EXTRA)
    cleaned = venn.cleaned_intersect([pybedtools.BedTool(str(path)) for path in paths])
    sets = {label: set(_sets_text(part)) for label, part in zip(labels, cleaned, strict=True)}
    venn.venn_maker(
        [str(path) for path in paths],
        names=labels,
        figure_filename=None if figure is None else str(figure),
        script_filename=str(script),
        additional_args=extras or None,
        run=False,
    )

    rscript = shutil.which("Rscript")
    render_result: dict[str, Any] | None = None
    if not render:
        note = "Rendering was not requested; the R script was written only."
    elif figure is None:
        note = "No figure_path was given, so the R script was written without rendering."
    elif rscript is None:
        note = (
            "Rscript was not found on PATH, so the R script was written but not run. "
            + _HINT_RSCRIPT
        )
    else:
        run = execute_bounded_command(
            rscript, [str(script)], backend_id="rscript", cwd=figure.parent, verbose=False
        )
        render_result = {
            "return_code": run.return_code,
            "log_tail": (run.output or "")[-_MAX_TAIL:],
        }
        note = (
            "Rendered with Rscript."
            if run.return_code == 0 and figure.is_file()
            else (
                "Rscript ran but did not produce the figure; see log_tail "
                "(VennDiagram must be installed in R)."
            )
        )
    return {
        "status": "completed",
        "backend": "pybedtools venn_maker + VennDiagram",
        "script_path": str(script),
        "figure_path": None if figure is None else str(figure),
        "rendered": bool(
            render_result is not None
            and render_result["return_code"] == 0
            and figure is not None
            and figure.is_file()
        ),
        "render_note": note,
        "render": render_result,
        "data": {
            "names": labels,
            "set_sizes": {label: len(sets[label]) for label in labels},
            "regions": _region_counts(sets),
        },
    }


def binary_heatmap(
    beds: Sequence[str],
    names: Sequence[str],
    output_png: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a binary heat map of which files contain each interval (multi-intersection) and return
    the class counts (pybedtools binary_heatmap).
    """

    paths = _existing_list(beds, "beds", 2, 30)
    labels = _names(names, len(paths))
    out = _new_path(output_png, "output_png")
    _program("bedtools", None, _HINT_BEDTOOLS)
    pybedtools = require("pybedtools", _PYBEDTOOLS_EXTRA)
    numpy = require("numpy", _PYBEDTOOLS_EXTRA)
    matplotlib = require("matplotlib", _PLOT_EXTRA)
    figure_class, canvas_class = _figure_classes()
    multi = pybedtools.BedTool().multi_intersect(
        i=[str(path) for path in paths], names=labels, cluster=True
    )
    class_counts: dict[str, int] = {}
    matrix: list[list[int]] = []
    for item in multi:
        if len(matrix) >= _MAX_ROWS:
            raise ConfigurationError(f"The multi-intersection has more than {_MAX_ROWS} rows.")
        cls = item[4]
        class_counts[cls] = class_counts.get(cls, 0) + 1
        matrix.append([int(value) for value in item[5:]])
    if not matrix:
        raise ConfigurationError("The input files share no intervals to draw.")
    array = numpy.asarray(matrix, dtype=int)
    order = _row_order(numpy, array)
    figure = figure_class(figsize=(max(3.0, 0.5 * len(labels) + 1.5), 6.0))
    axis = figure.add_subplot(1, 1, 1)
    axis.imshow(
        array[order],
        aspect="auto",
        interpolation="nearest",
        cmap=matplotlib.colormaps["binary"],
        origin="upper",
    )
    axis.set_xticks(list(range(len(labels))))
    axis.set_xticklabels(labels, rotation=90)
    canvas_class(figure)
    figure.savefig(str(out), format="png", dpi=120)
    return {
        "status": "completed",
        "backend": "pybedtools + bedtools multiinter + matplotlib",
        "image_path": str(out),
        "data": {
            "names": labels,
            "interval_count": len(matrix),
            "class_counts": dict(
                sorted(class_counts.items(), key=lambda pair: (-pair[1], pair[0]))
            ),
            "matrix": matrix,
            "row_order": order,
        },
    }


def add_color(
    bed: str,
    output_bed: str,
    *,
    colormap: str = "viridis",
    score_min: float | None = None,
    score_max: float | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a BED file whose itemRgb colour follows each feature's score through a matplotlib
    colormap (pybedtools featurefuncs.add_color).
    """

    source = _existing(bed, "bed")
    out = _new_path(output_bed, "output_bed")
    name = _text(colormap, "colormap")
    matplotlib = require("matplotlib", _PLOT_EXTRA)
    colors = require("matplotlib.colors", _PLOT_EXTRA)
    try:
        cmap = matplotlib.colormaps[name]
    except KeyError as exc:
        raise ConfigurationError(f"Unknown matplotlib colormap: {name}") from exc
    pybedtools = require("pybedtools", _PYBEDTOOLS_EXTRA)
    featurefuncs = require("pybedtools.featurefuncs", _PYBEDTOOLS_EXTRA)
    intervals = pybedtools.BedTool(str(source))
    scores: list[float] = []
    for feature in intervals:
        if len(feature.fields) < 5:
            raise ConfigurationError("Every feature needs a numeric score in the fifth column.")
        try:
            scores.append(float(feature.score))
        except ValueError as exc:
            raise ConfigurationError(f"Score is not a number: {feature.score!r}") from exc
        if len(scores) > _MAX_ROWS:
            raise ConfigurationError(f"The BED file has more than {_MAX_ROWS} features.")
    if not scores:
        raise ConfigurationError("The BED file has no features.")
    low = min(scores) if score_min is None else float(score_min)
    high = max(scores) if score_max is None else float(score_max)
    if not math.isfinite(low) or not math.isfinite(high) or low > high:
        raise ConfigurationError(
            "score_min and score_max must be finite numbers with score_min <= score_max."
        )
    norm = colors.Normalize(vmin=low, vmax=high)
    intervals.each(featurefuncs.add_color, cmap=cmap, norm=norm).saveas(str(out))
    rgb = [
        line.split("\t")[8]
        for line in out.read_text().splitlines()
        if line and not line.startswith(("#", "track", "browser"))
    ]
    return {
        "status": "completed",
        "backend": "pybedtools featurefuncs + matplotlib",
        "output_path": str(out),
        "data": {
            "feature_count": len(rgb),
            "colormap": name,
            "score_min": low,
            "score_max": high,
            "rgb": rgb,
        },
    }


__all__ = [
    "add_color",
    "bam_to_bigwig",
    "bedgraph_to_bigwig",
    "bigbed",
    "bigbed_to_bed",
    "bigwig_to_bedgraph",
    "binary_heatmap",
    "chromsizes",
    "cleaned_intersect",
    "get_chromsizes_from_genomepy",
    "get_chromsizes_from_ucsc",
    "mapped_read_count",
    "sort_binary_matrix",
    "venn_maker",
    "wig_to_bigwig",
]
_register_public_verbose(globals())
