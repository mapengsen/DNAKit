"""FreeBayes companion tools: region partitioning, VCF header merging and BAM left-alignment of indels.

Four of the five operations are the small helper scripts that ship with FreeBayes (1.3.x) and are re-implemented here
in Python with their exact semantics, so that no script has to be on PATH: ``fasta_generate_regions.py``,
``coverage_to_regions.py``, ``split_ref_by_bai_datasize.py`` and ``vcffirstheader``. The fifth, ``bamleftalign``, is a
compiled FreeBayes program run through the shared command-line facade. Results are JSON-serialisable dictionaries.
"""

from __future__ import annotations

import math
import struct
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p16_cli as cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_FREEBAYES = cli.Tool(
    name="freebayes",
    executable="bamleftalign",
    hint=(
        "Install FreeBayes (for example `conda install -c bioconda freebayes`), which provides bamleftalign, "
        "then put it on PATH or pass bamleftalign_path."
    ),
    tested_version="1.3.10",
    citation="https://arxiv.org/abs/1207.3907",
    version_arguments=("--help",),
    version_pattern=r"version[: ]+v?(\d+(?:\.\d+)+)",
)
_MAX_LISTED_REGIONS = 100_000


def _read_fai(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0]:
                continue
            try:
                rows.append((fields[0], int(fields[1])))
            except ValueError as exc:
                raise ConfigurationError(f"Malformed FASTA index line in {path}: {line!r}") from exc
    if not rows:
        raise ConfigurationError(f"The FASTA index has no sequences: {path}")
    return rows


def _index_path(value: str | Path, name: str) -> Path:
    """The ``.fai`` file for a reference: the file itself, or ``<fasta>.fai`` when a FASTA path is given."""
    path = Path(value).expanduser().resolve()
    if path.suffix != ".fai":
        path = path.with_name(path.name + ".fai")
    if not path.is_file():
        raise ConfigurationError(f"{name} needs an existing FASTA index (.fai); not found: {path}")
    return path


def partition_reference_into_regions(
    reference: str,
    region_size: int,
    *,
    output_path: str | None = None,
    chunks: bool = False,
    chromosomes: Sequence[str] | None = None,
    bed_prefix: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a reference into fixed-size ``chrom:start-end`` regions (FreeBayes ``fasta_generate_regions.py``).

    ``reference`` is a FASTA path or its ``.fai`` index; the index must exist. With ``chunks=False`` each region is
    ``region_size`` bases long; with ``chunks=True`` every sequence is divided into ``region_size`` equal chunks.
    ``chromosomes`` limits the sequences. Regions are 0-based half-open, exactly as the original script prints them.
    ``output_path`` (a new file) receives one region per line; ``bed_prefix`` instead writes one BED file per region
    named ``<bed_prefix>.<chrom>.region.<n>.bed`` (the Snakemake mode of the script). Returns the (bounded) region list.
    """
    index = _index_path(reference, "reference")
    size = cli.integer(region_size, "region_size", minimum=1)
    by_chunks = cli.flag(chunks, "chunks")
    wanted = (
        None if chromosomes is None else {cli.text(name, "chromosomes") for name in chromosomes}
    )
    target = None if output_path is None else cli.new_file(output_path, "output_path")
    prefix = None if bed_prefix is None else cli.text(bed_prefix, "bed_prefix")
    if target is not None and prefix is not None:
        raise ConfigurationError("Give output_path or bed_prefix, not both.")
    regions: list[str] = []
    bed_files: list[str] = []
    count = 0
    for name, length in _read_fai(index):
        if wanted is not None and name not in wanted:
            continue
        step = math.ceil(length / size) if by_chunks else size
        start = 0
        while start < length:
            end = min(start + step, length)
            count += 1
            if prefix is not None:
                number = str(math.ceil(end / step))
                bed = Path(f"{prefix}.{name}.region.{number}.bed")
                with bed.open("x", encoding="utf-8") as handle:
                    handle.write("\t".join([name, str(start), str(end)]))
                bed_files.append(str(bed.resolve()))
            elif len(regions) < _MAX_LISTED_REGIONS:
                regions.append(f"{name}:{start}-{end}")
            start = end
    if target is not None:
        target.write_text("".join(f"{region}\n" for region in regions), encoding="utf-8")
    if wanted is not None and count == 0:
        raise ConfigurationError("None of the requested chromosomes occur in the FASTA index.")
    return {
        "status": "completed",
        "backend": "freebayes-scripts",
        "operation": "partition_reference_into_regions",
        "region_count": count,
        "regions": regions,
        "regions_truncated": count > len(regions) and not bed_files,
        "outputs": [str(target)] if target is not None else bed_files,
        "provenance": {"tool": "FreeBayes fasta_generate_regions.py", "reimplemented_in": "dnakit"},
    }


def partition_reference_by_coverage(
    reference: str,
    coverage: str,
    num_regions: int,
    *,
    output_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Make regions of about equal total coverage from per-position depth (FreeBayes ``coverage_to_regions.py``).

    ``coverage`` is the three-column ``chrom<TAB>pos<TAB>depth`` output of ``bamtools coverage`` (one line per
    covered position, sorted); ``reference`` is the FASTA or its ``.fai``. A new region is closed when the running
    depth sum exceeds ``total_depth / num_regions``; the last region of each sequence ends at its length.
    Returns ``chrom:start-end`` regions as the original script prints them.
    """
    index = _index_path(reference, "reference")
    source = cli.existing_file(coverage, "coverage")
    parts = cli.integer(num_regions, "num_regions", minimum=1)
    target = None if output_path is None else cli.new_file(output_path, "output_path")
    lengths = dict(_read_fai(index))
    positions: list[tuple[str, int, int]] = []
    total = 0
    with cli.open_text(source) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) != 3:
                raise ConfigurationError(
                    f"coverage lines need 3 tab-separated columns (chrom, pos, depth): {line!r}"
                )
            chrom, pos, depth = fields[0], int(fields[1]), int(fields[2])
            if chrom not in lengths:
                raise ConfigurationError(
                    f"Sequence {chrom!r} of the coverage file is not in the FASTA index."
                )
            positions.append((chrom, pos, depth))
            total += depth
    if not positions:
        raise ConfigurationError("The coverage file has no positions.")
    per_region = total / parts
    regions: list[str] = []
    last_chrom: str | None = None
    last_pos = 0
    in_region = 0
    for chrom, pos, depth in positions:
        if last_chrom != chrom:
            if last_chrom:
                regions.append(f"{last_chrom}:{last_pos}-{lengths[last_chrom]}")
                last_pos = 0
                last_chrom = chrom
                in_region = 0
            else:
                last_chrom = chrom
        in_region += depth
        if in_region > per_region:
            regions.append(f"{chrom}:{last_pos}-{pos}")
            last_pos = pos
            in_region = 0
    assert last_chrom is not None
    regions.append(f"{last_chrom}:{last_pos}-{lengths[last_chrom]}")
    if target is not None:
        target.write_text("".join(f"{region}\n" for region in regions), encoding="utf-8")
    return {
        "status": "completed",
        "backend": "freebayes-scripts",
        "operation": "partition_reference_by_coverage",
        "region_count": len(regions),
        "regions": regions[:_MAX_LISTED_REGIONS],
        "total_depth": total,
        "outputs": [str(target)] if target is not None else [],
        "provenance": {"tool": "FreeBayes coverage_to_regions.py", "reimplemented_in": "dnakit"},
    }


def _bai_interval_offsets(path: Path, chromosome_count: int) -> list[list[int]]:
    """The linear-index (interval) virtual offsets per reference of a BAI file; the bin index is skipped."""
    offsets: list[list[int]] = []
    with path.open("rb") as handle:
        if handle.read(4) != b"BAI\x01":
            raise ConfigurationError(f"Not a BAM index (BAI): {path}")
        (n_ref,) = struct.unpack("<i", handle.read(4))
        if n_ref != chromosome_count:
            raise ConfigurationError(
                f"The FASTA index and {path.name} must have the same number of sequences ({chromosome_count} vs {n_ref})."
            )
        for _ in range(n_ref):
            (n_bin,) = struct.unpack("<i", handle.read(4))
            for _ in range(n_bin):
                handle.read(4)
                (n_chunk,) = struct.unpack("<i", handle.read(4))
                handle.seek(n_chunk * 16, 1)
            (n_intv,) = struct.unpack("<i", handle.read(4))
            offsets.append(
                list(struct.unpack(f"<{n_intv}Q", handle.read(n_intv * 8))) if n_intv > 0 else []
            )
    return offsets


def partition_reference_by_bam_size(
    reference: str,
    bam_files: Sequence[str],
    *,
    output_path: str | None = None,
    target_data_size: float = 100e6,
    bai_interval_size: int = 16384,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a reference into regions holding equal amounts of BAM data (FreeBayes ``split_ref_by_bai_datasize.py``).

    The data amount per 16 kb interval is read from the linear index of every ``<bam>.bai`` (the difference of the
    BGZF virtual file offsets, used as a size exactly as the original script does) and summed over the BAM files; each
    sequence is cut where the cumulated amount reaches multiples of ``target_data_size`` (the script's ``-s``, default
    ``100e6``, multiplied by 1,000,000 as in the script; values below 1 are refused). Sequences with fewer than two
    index intervals are returned whole. Needs numpy and scipy (``dnakit[bio-analysis]``). Rows are
    ``(chrom, start, end)``; ``output_path`` receives them tab-separated as the original script prints them, and
    ``regions`` in the result is cut at 100,000 rows.
    """
    np = require("numpy", "bio-analysis")
    interpolate = require("scipy.interpolate", "bio-analysis")
    index = _index_path(reference, "reference")
    bams = cli.existing_files(bam_files, "bam_files")
    target = None if output_path is None else cli.new_file(output_path, "output_path")
    interval = cli.integer(bai_interval_size, "bai_interval_size", minimum=1)
    size = int(float(cli.number(target_data_size, "target_data_size", minimum=1e-6))) * 1_000_000
    if size <= 0:
        raise ConfigurationError("target_data_size is too small.")
    chromosomes = _read_fai(index)
    interval_counts = [int(np.ceil(length / interval)) for _, length in chromosomes]
    cumulative = [np.zeros((count,), dtype=np.int64) for count in interval_counts]
    for bam in bams:
        bai = bam.with_name(bam.name + ".bai")
        if not bai.is_file():
            raise ConfigurationError(f"BAM index not found: {bai}")
        for ci, raw in enumerate(_bai_interval_offsets(bai, len(chromosomes))):
            if not raw:
                continue
            ioff = np.array(raw, dtype=np.int64)
            while len(ioff) < len(cumulative[ci]):
                ioff = np.append(ioff, ioff[-1] + 1)
            cumulative[ci] += ioff - ioff[0]
    regions: list[tuple[str, int, int]] = []
    for ci, (name, length) in enumerate(chromosomes):
        data = cumulative[ci]
        if len(data) < 2:
            regions.append((name, 0, length))
            continue
        positions = np.arange(0, data.shape[0]) * interval
        forward = interpolate.interp1d(positions, data, fill_value="extrapolate", kind="linear")
        data_total = forward([length])[0]
        parts = int(np.ceil(data_total / size))
        backward = interpolate.interp1d(data, positions, fill_value="extrapolate", kind="linear")
        data_points = np.linspace(0, data_total, num=parts + 1, endpoint=True, dtype=int)
        cuts = backward(data_points).astype(int)
        cuts[0] = 0
        cuts[-1] = length
        regions.extend((name, int(cuts[i]), int(cuts[i + 1])) for i in range(len(cuts) - 1))
    if target is not None:
        target.write_text("".join(f"{a}\t{b}\t{c}\n" for a, b, c in regions), encoding="utf-8")
    return {
        "status": "completed",
        "backend": "freebayes-scripts",
        "operation": "partition_reference_by_bam_size",
        "region_count": len(regions),
        "regions": [list(row) for row in regions[:_MAX_LISTED_REGIONS]],
        "outputs": [str(target)] if target is not None else [],
        "provenance": {
            "tool": "FreeBayes split_ref_by_bai_datasize.py",
            "reimplemented_in": "dnakit",
        },
    }


def keep_first_vcf_header(
    inputs: Sequence[str],
    output_path: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Concatenate VCF files keeping only the first header block (FreeBayes ``vcffirstheader``).

    The files are read in order as one stream (as ``cat a.vcf b.vcf | vcffirstheader`` would). Meta lines (``##``) and
    the ``#CHROM`` line are written until the first ``#CHROM`` line has been seen; every later header line is dropped.
    Data lines are copied with trailing white space removed. ``.gz`` inputs are read transparently.
    """
    sources = cli.existing_files(list(inputs), "inputs")
    target = cli.new_file(output_path, "output_path")
    header_open = True
    headers = records = dropped = 0
    with target.open("x", encoding="utf-8") as out:
        for source in sources:
            with cli.open_text(source) as handle:
                for line in handle:
                    if line.startswith("##"):
                        if header_open:
                            out.write(line.strip() + "\n")
                            headers += 1
                        else:
                            dropped += 1
                        continue
                    if line.startswith("#"):
                        if header_open:
                            out.write(line.strip() + "\n")
                            headers += 1
                            header_open = False
                        else:
                            dropped += 1
                        continue
                    out.write(line.strip() + "\n")
                    records += 1
    return {
        "status": "completed",
        "backend": "freebayes-scripts",
        "operation": "keep_first_vcf_header",
        "header_lines_written": headers,
        "header_lines_dropped": dropped,
        "records_written": records,
        "outputs": [str(target)],
        "provenance": {"tool": "FreeBayes vcffirstheader", "reimplemented_in": "dnakit"},
    }


def left_align_bam_indels(
    alignment: str,
    output_path: str,
    *,
    reference: str,
    compressed: bool = True,
    max_iterations: int | None = None,
    bamleftalign_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Left-align and merge indels in BAM alignments (FreeBayes ``bamleftalign``).

    ``reference`` is the FASTA the reads were mapped to (required by the program, ``--fasta-reference``). The BAM is
    streamed to the program's standard input and the realigned BAM is written to ``output_path`` (a new file).
    ``compressed=True`` writes a compressed BAM (``--compressed``; the program's own default is uncompressed).
    ``max_iterations`` limits the left-realignment passes (``--max-iterations``).
    """
    bam = cli.existing_file(alignment, "alignment")
    target = cli.new_file(output_path, "output_path")
    arguments = ["--fasta-reference", str(cli.existing_file(reference, "reference"))]
    if cli.flag(compressed, "compressed"):
        arguments.append("--compressed")
    if max_iterations is not None:
        arguments += [
            "--max-iterations",
            str(cli.integer(max_iterations, "max_iterations", minimum=1)),
        ]
    record = cli.invoke(
        _FREEBAYES,
        arguments,
        cwd=target.parent,
        stdout_to=target,
        stdin_from=bam,
        explicit=bamleftalign_path,
    )
    return cli.finish(record, "left_align_bam_indels", alignment=str(bam), output=str(target))


__all__ = [
    "keep_first_vcf_header",
    "left_align_bam_indels",
    "partition_reference_by_bam_size",
    "partition_reference_by_coverage",
    "partition_reference_into_regions",
]

_register_public_verbose(globals(), tuple(__all__))
