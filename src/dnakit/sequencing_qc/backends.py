"""Optional-backend sequencing QC helpers for cutadapt, umi-tools and pyslow5.

The functions in this module call the Python API of ``cutadapt``:

- ``find_adapter_match()`` aligns one adapter to one DNA sequence and returns the
  match coordinates and the trimmed sequence.
- ``trim_fastq_adapters()`` applies the same adapter to a FASTQ file and writes the
  reads that remain. Reads are read and written with ``dnaio``, which is installed
  together with cutadapt.
- ``find_quality_trim_bounds()`` and ``calculate_expected_errors()`` use cutadapt's
  compiled quality-trimming helpers on Phred-encoded quality strings.

The ``cluster_umi_counts()``, ``extract_read_barcodes()`` and
``map_barcodes_to_whitelist()`` functions call ``umi-tools`` (UMI clustering,
regex barcode extraction and cell-barcode whitelist correction). The
``write_slow5_signals()``, ``summarize_slow5_signals()`` and
``read_slow5_signal()`` functions call ``pyslow5`` for SLOW5/BLOW5 raw signals.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import math
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require, sequence_text
from dnakit.exceptions import ConfigurationError

_CUTADAPT_EXTRA = "bio-cutadapt"
_ADAPTER_POSITIONS = ("3prime", "5prime", "anywhere")
_DEFAULT_MAX_READS = 1_000_000
_MAX_READS_LIMIT = 10_000_000
_PHRED_OFFSETS = (33, 64)
_MAX_PHRED = 93


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _choose(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return value


def _positive_int(value: object, name: str, *, maximum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")
    if maximum is not None and value > maximum:
        raise ConfigurationError(f"{name} must not exceed {maximum}.")
    return value


def _error_rate(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError("max_errors must be a number in [0, 1).")
    number = float(value)
    if not math.isfinite(number) or not 0 <= number < 1:
        raise ConfigurationError("max_errors must be a finite error rate in [0, 1).")
    return number


def _path_text(value: object, name: str) -> str:
    text = str(value) if isinstance(value, (str, Path)) else ""
    if not text or "\0" in text:
        raise ConfigurationError(f"{name} must be a non-empty path without NUL characters.")
    return text


def _new_output(value: object, source: Path | None = None) -> Path:
    path = Path(_path_text(value, "output")).expanduser().resolve()
    if source is not None and path == source:
        raise ConfigurationError("output must be a different file from the input.")
    if path.exists():
        raise ConfigurationError(f"Output file already exists: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    return path


def _phred_text(qualities: object, phred_offset: object) -> tuple[str, int]:
    offset = _choose_offset(phred_offset)
    if not isinstance(qualities, str) or not qualities:
        raise ConfigurationError("qualities must be a non-empty Phred quality string.")
    if any(not offset <= ord(char) <= offset + _MAX_PHRED for char in qualities):
        raise ConfigurationError(f"qualities contain characters outside the Phred+{offset} range.")
    return qualities, offset


def _choose_offset(value: object) -> int:
    if isinstance(value, bool) or value not in _PHRED_OFFSETS:
        raise ConfigurationError("phred_offset must be 33 or 64.")
    return int(value)


def _cutoff(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= _MAX_PHRED:
        raise ConfigurationError(f"{name} must be an integer in [0, {_MAX_PHRED}].")
    return value


def _adapter(module: Any, sequence: str, position: str, max_errors: float, min_overlap: int) -> Any:
    options = {"max_errors": max_errors, "min_overlap": min_overlap, "name": "adapter"}
    try:
        if position == "3prime":
            return module.BackAdapter(sequence, **options)
        if position == "5prime":
            return module.FrontAdapter(sequence, **options)
        return module.AnywhereAdapter(sequence, **options)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"cutadapt rejected the adapter: {exc}") from exc


def find_adapter_match(
    sequence: str,
    adapter: str,
    *,
    position: str = "3prime",
    max_errors: float = 0.1,
    min_overlap: int = 3,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align one adapter to a DNA sequence with cutadapt and return the trimmed sequence.

    ``position`` is ``"3prime"``, ``"5prime"`` or ``"anywhere"``. ``start`` and ``stop``
    are 0-based, half-open coordinates of the aligned adapter. Only the first
    best match is reported, as in cutadapt.
    """

    cutadapt = require("cutadapt.adapters", _CUTADAPT_EXTRA)
    text = sequence_text(sequence)
    adapter_text = sequence_text(adapter)
    selected = _choose(position, _ADAPTER_POSITIONS, "position")
    matcher = _adapter(
        cutadapt,
        adapter_text,
        selected,
        _error_rate(max_errors),
        _positive_int(min_overlap, "min_overlap"),
    )
    match = matcher.match_to(text)
    if match is None:
        matched, start, stop, errors, score, trimmed = False, None, None, None, None, text
    else:
        matched = True
        start, stop = int(match.rstart), int(match.rstop)
        errors, score = int(match.errors), int(match.score)
        trimmed = str(match.trimmed(text))
    return {
        "status": "completed",
        "operation": "find_adapter_match",
        "adapter": adapter_text,
        "position": selected,
        "matched": matched,
        "start": start,
        "stop": stop,
        "errors": errors,
        "score": score,
        "trimmed_sequence": trimmed,
        "provenance": {
            "package": "cutadapt",
            "version": _version("cutadapt"),
            "coordinate_system": "0-based half-open read coordinates",
            "method": "cutadapt adapter aligner (Python API)",
        },
    }


def trim_fastq_adapters(
    source: str | Path,
    output: str | Path,
    adapter: str,
    *,
    position: str = "3prime",
    max_errors: float = 0.1,
    min_overlap: int = 3,
    minimum_length: int = 1,
    max_reads: int = _DEFAULT_MAX_READS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim one adapter from every FASTQ read with cutadapt and write the retained reads.

    Reads shorter than ``minimum_length`` after trimming are discarded and counted.
    ``output`` must not exist yet. More than ``max_reads`` input reads raises an error
    and removes the partial output file.
    """

    cutadapt_adapters = require("cutadapt.adapters", _CUTADAPT_EXTRA)
    modifiers = require("cutadapt.modifiers", _CUTADAPT_EXTRA)
    dnaio = require("dnaio", _CUTADAPT_EXTRA)
    path = input_path(_path_text(source, "source"))
    target = _new_output(output, path)
    adapter_text = sequence_text(adapter)
    selected = _choose(position, _ADAPTER_POSITIONS, "position")
    limit = _positive_int(max_reads, "max_reads", maximum=_MAX_READS_LIMIT)
    minimum = _positive_int(minimum_length, "minimum_length")
    matcher = _adapter(
        cutadapt_adapters,
        adapter_text,
        selected,
        _error_rate(max_errors),
        _positive_int(min_overlap, "min_overlap"),
    )
    cutter = modifiers.AdapterCutter([matcher], times=1, action="trim", index=False)
    counts = {"reads": 0, "reads_with_adapter": 0, "reads_written": 0, "reads_too_short": 0}
    bases_removed = 0
    try:
        with (
            dnaio.open(str(path), fileformat="fastq", mode="r") as reader,
            dnaio.open(str(target), fileformat="fastq", mode="w") as writer,
        ):
            for record in reader:
                counts["reads"] += 1
                if counts["reads"] > limit:
                    raise ConfigurationError(
                        f"FASTQ input exceeds max_reads={limit}; increase the explicit limit."
                    )
                trimmed, matches = cutter.match_and_trim(record)
                if matches:
                    counts["reads_with_adapter"] += 1
                bases_removed += len(record.sequence) - len(trimmed.sequence)
                if len(trimmed.sequence) < minimum:
                    counts["reads_too_short"] += 1
                    continue
                writer.write(trimmed)
                counts["reads_written"] += 1
    except ConfigurationError:
        target.unlink(missing_ok=True)
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"cutadapt/dnaio could not trim {path.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "trim_fastq_adapters",
        "source": str(path),
        "output": str(target),
        "adapter": adapter_text,
        "position": selected,
        **counts,
        "bases_removed": bases_removed,
        "provenance": {
            "package": "cutadapt",
            "version": _version("cutadapt"),
            "coordinate_system": "0-based half-open read coordinates",
            "method": "cutadapt AdapterCutter (Python API) with dnaio FASTQ I/O",
            "read_file_format": "fastq",
        },
    }


def find_quality_trim_bounds(
    qualities: str,
    *,
    cutoff_front: int = 0,
    cutoff_back: int = 20,
    phred_offset: int = 33,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the retained 0-based interval after cutadapt's quality trimming.

    The trimming follows the BWA algorithm: the running sum of ``cutoff - quality``
    from each end is maximized, which removes low-quality ends.
    """

    qualtrim = require("cutadapt.qualtrim", _CUTADAPT_EXTRA)
    text, offset = _phred_text(qualities, phred_offset)
    front = _cutoff(cutoff_front, "cutoff_front")
    back = _cutoff(cutoff_back, "cutoff_back")
    start, stop = qualtrim.quality_trim_index(text, front, back, offset)
    return {
        "status": "completed",
        "operation": "find_quality_trim_bounds",
        "length": len(text),
        "start": int(start),
        "stop": int(stop),
        "retained_length": int(stop - start),
        "provenance": {
            "package": "cutadapt",
            "version": _version("cutadapt"),
            "coordinate_system": "0-based half-open read coordinates",
            "method": "cutadapt quality_trim_index (BWA algorithm)",
        },
    }


def calculate_expected_errors(
    qualities: str, *, phred_offset: int = 33, verbose: bool = True
) -> dict[str, Any]:
    """Expected number of base-calling errors in a read, from its Phred qualities.

    Each base contributes ``10 ** (-Q / 10)``. The sum is computed by cutadapt.
    """

    qualtrim = require("cutadapt.qualtrim", _CUTADAPT_EXTRA)
    text, offset = _phred_text(qualities, phred_offset)
    value = float(qualtrim.expected_errors(text, offset))
    return {
        "status": "completed",
        "operation": "calculate_expected_errors",
        "length": len(text),
        "expected_errors": value if math.isfinite(value) else None,
        "provenance": {
            "package": "cutadapt",
            "version": _version("cutadapt"),
            "statistic": "sum of 10^(-Q/10) over bases (cutadapt expected_errors)",
        },
    }


_PYSLOW5_EXTRA = "bio-pyslow5"
_SLOW5_SUFFIXES = (".slow5", ".blow5")
_SLOW5_ASCII_MAGIC = b"#slow5_version"
_BLOW5_MAGIC = b"BLOW5\x01"
_DEFAULT_MAX_SIGNAL_READS = 100_000
_DEFAULT_MAX_SAMPLES = 100_000
_MAX_SAMPLES_LIMIT = 10_000_000
_INT16_MIN, _INT16_MAX = -32768, 32767


def _finite_number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    return float(value)


def _positive_number(value: object, name: str) -> float:
    number = _finite_number(value, name)
    if number <= 0:
        raise ConfigurationError(f"{name} must be a positive number.")
    return number


def _read_id(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 1024
        or any(char.isspace() or char == "\0" for char in value)
    ):
        raise ConfigurationError("read_id must be a non-empty string without whitespace.")
    return value


def _slow5_source(value: object) -> Path:
    """Check the suffix and magic bytes first; pyslow5 1.5.0 can crash on malformed headers."""

    path = input_path(_path_text(value, "source"))
    if path.suffix not in _SLOW5_SUFFIXES:
        raise ConfigurationError("source must be a .slow5 or .blow5 file.")
    with path.open("rb") as handle:
        head = handle.read(16)
    if not head.startswith((_SLOW5_ASCII_MAGIC, _BLOW5_MAGIC)):
        raise ConfigurationError(f"{path.name} does not start with a SLOW5 or BLOW5 header.")
    return path


def _slow5_output(value: object) -> Path:
    path = _new_output(value)
    if path.suffix not in _SLOW5_SUFFIXES:
        raise ConfigurationError("output must end with .slow5 or .blow5.")
    return path


def _open_slow5(pyslow5: Any, path: Path, mode: str) -> Any:
    try:
        return pyslow5.Open(str(path), mode)
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pyslow5 could not open {path.name}: {exc}") from exc


def _slow5_reads(handle: Any, *, pA: bool, limit: int) -> Iterator[dict[str, Any]]:
    for count, read in enumerate(handle.seq_reads(pA=pA), start=1):
        if count > limit:
            raise ConfigurationError(
                f"SLOW5 input exceeds max_reads={limit}; increase the explicit limit."
            )
        yield read


def _round_value(value: object) -> float | None:
    # pyslow5 converts to float32; 7 significant digits removes conversion noise.
    number = float(value)
    return float(f"{number:.7g}") if math.isfinite(number) else None


def summarize_slow5_signals(
    source: str | Path,
    *,
    max_reads: int = _DEFAULT_MAX_SIGNAL_READS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count reads and raw-signal samples in a SLOW5 or BLOW5 file, using pyslow5.

    The file is scanned once, sequentially, so no index file is written next to it.
    More than ``max_reads`` reads raises an error.
    """

    pyslow5 = require("pyslow5", _PYSLOW5_EXTRA)
    path = _slow5_source(source)
    limit = _positive_int(max_reads, "max_reads", maximum=_MAX_READS_LIMIT)
    handle = _open_slow5(pyslow5, path, "r")
    lengths: list[int] = []
    rates: set[float] = set()
    try:
        for read in _slow5_reads(handle, pA=False, limit=limit):
            lengths.append(len(read["signal"]))
            rates.add(float(read["sampling_rate"]))
        read_groups = int(handle.get_num_read_groups())
        frequency = handle.get_header_value("sample_frequency")
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pyslow5 could not read {path.name}: {exc}") from exc
    finally:
        handle.close()
    count = len(lengths)
    return {
        "status": "completed",
        "operation": "summarize_slow5_signals",
        "source": str(path),
        "reads": count,
        "read_groups": read_groups,
        "signal_samples": sum(lengths),
        "min_signal_samples": min(lengths) if lengths else None,
        "mean_signal_samples": (sum(lengths) / count) if lengths else None,
        "max_signal_samples": max(lengths) if lengths else None,
        "sampling_rates": sorted(rates),
        "header_sample_frequency": None if frequency is None else str(frequency),
        "provenance": {
            "package": "pyslow5",
            "version": _version("pyslow5"),
            "method": "sequential seq_reads scan (no index file written)",
        },
    }


def read_slow5_signal(
    source: str | Path,
    read_id: str,
    *,
    pA: bool = False,
    max_reads: int = _DEFAULT_MAX_SIGNAL_READS,
    max_samples: int = _DEFAULT_MAX_SAMPLES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the raw signal of one read from a SLOW5 or BLOW5 file, using pyslow5.

    The file is scanned sequentially until ``read_id`` is found, so no index file is
    written. With ``pA=True`` pyslow5 converts the values to picoamps; otherwise they
    are the stored integer DAC values.
    """

    pyslow5 = require("pyslow5", _PYSLOW5_EXTRA)
    path = _slow5_source(source)
    target = _read_id(read_id)
    if not isinstance(pA, bool):
        raise ConfigurationError("pA must be a boolean.")
    limit = _positive_int(max_reads, "max_reads", maximum=_MAX_READS_LIMIT)
    sample_limit = _positive_int(max_samples, "max_samples", maximum=_MAX_SAMPLES_LIMIT)
    handle = _open_slow5(pyslow5, path, "r")
    try:
        for read in _slow5_reads(handle, pA=pA, limit=limit):
            if read["read_id"] != target:
                continue
            signal = read["signal"]
            if len(signal) > sample_limit:
                raise ConfigurationError(
                    f"read {target!r} has {len(signal)} samples; max_samples={sample_limit}."
                )
            values = [_round_value(value) for value in signal] if pA else [int(v) for v in signal]
            return {
                "status": "completed",
                "operation": "read_slow5_signal",
                "source": str(path),
                "read_id": target,
                "read_group": int(read["read_group"]),
                "units": "pA" if pA else "raw DAC",
                "length": len(values),
                "digitisation": float(read["digitisation"]),
                "offset": float(read["offset"]),
                "range": float(read["range"]),
                "sampling_rate": float(read["sampling_rate"]),
                "signal": values,
                "provenance": {
                    "package": "pyslow5",
                    "version": _version("pyslow5"),
                    "method": "sequential seq_reads scan (no index file written)",
                },
            }
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pyslow5 could not read {path.name}: {exc}") from exc
    finally:
        handle.close()
    raise ConfigurationError(f"read_id {target!r} was not found in {path.name}.")


def _slow5_records(records: object, limit: int) -> list[dict[str, Any]]:
    if isinstance(records, (str, bytes, Mapping)) or not isinstance(records, Sequence):
        raise ConfigurationError("records must be a non-empty list of raw-signal mappings.")
    if not records:
        raise ConfigurationError("records must be a non-empty list of raw-signal mappings.")
    if len(records) > limit:
        raise ConfigurationError(f"records exceed max_reads={limit}; increase the explicit limit.")
    seen: set[str] = set()
    items: list[dict[str, Any]] = []
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise ConfigurationError(f"record {index} must be a mapping.")
        read_id = _read_id(record.get("read_id"))
        if read_id in seen:
            raise ConfigurationError(f"duplicate read_id {read_id!r}.")
        seen.add(read_id)
        signal = record.get("signal")
        if isinstance(signal, (str, bytes)) or not isinstance(signal, Sequence) or not signal:
            raise ConfigurationError(f"read {read_id!r} needs a non-empty list of integers.")
        values: list[int] = []
        for value in signal:
            if (
                isinstance(value, bool)
                or not isinstance(value, int)
                or not _INT16_MIN <= value <= _INT16_MAX
            ):
                raise ConfigurationError(
                    f"signal values of {read_id!r} must be integers in [-32768, 32767]."
                )
            values.append(value)
        items.append(
            {
                "read_id": read_id,
                "signal": values,
                "digitisation": _positive_number(record.get("digitisation"), "digitisation"),
                "offset": _finite_number(record.get("offset"), "offset"),
                "range": _positive_number(record.get("range"), "range"),
                "sampling_rate": _positive_number(record.get("sampling_rate"), "sampling_rate"),
            }
        )
    if len({item["sampling_rate"] for item in items}) != 1:
        raise ConfigurationError("all records must share one sampling_rate.")
    return items


def write_slow5_signals(
    records: Sequence[Mapping[str, Any]],
    output: str | Path,
    *,
    max_reads: int = _DEFAULT_MAX_SIGNAL_READS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write raw-signal records to a new SLOW5 (.slow5) or BLOW5 (.blow5) file, using pyslow5.

    Each record needs ``read_id``, ``signal`` (a list of integers in the int16 range),
    ``digitisation``, ``offset``, ``range`` and ``sampling_rate``. All records share one
    sampling rate, which is written to the ``sample_frequency`` header. The format follows
    the output suffix. ``output`` must not exist yet.
    """

    pyslow5 = require("pyslow5", _PYSLOW5_EXTRA)
    numpy = require("numpy", _PYSLOW5_EXTRA)
    target = _slow5_output(output)
    limit = _positive_int(max_reads, "max_reads", maximum=_MAX_READS_LIMIT)
    items = _slow5_records(records, limit)
    sampling_rate = items[0]["sampling_rate"]
    samples_written = 0
    try:
        handle = _open_slow5(pyslow5, target, "w")
        try:
            if handle.write_header({"sample_frequency": f"{sampling_rate:g}"}) != 0:
                raise ConfigurationError("pyslow5 could not write the SLOW5 header.")
            for item in items:
                record, aux = handle.get_empty_record(aux=True)
                record.update(
                    read_id=item["read_id"],
                    read_group=0,
                    digitisation=item["digitisation"],
                    offset=item["offset"],
                    range=item["range"],
                    sampling_rate=item["sampling_rate"],
                    len_raw_signal=len(item["signal"]),
                    signal=numpy.array(item["signal"], dtype=numpy.int16),
                )
                if handle.write_record(record, aux) != 0:
                    raise ConfigurationError(f"pyslow5 could not write read {item['read_id']!r}.")
                samples_written += len(item["signal"])
        finally:
            handle.close()
    except ConfigurationError:
        target.unlink(missing_ok=True)
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"pyslow5 could not write {target.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_slow5_signals",
        "output": str(target),
        "reads_written": len(items),
        "signal_samples_written": samples_written,
        "sampling_rate": float(sampling_rate),
        "provenance": {
            "package": "pyslow5",
            "version": _version("pyslow5"),
            "format": "BLOW5" if target.suffix == ".blow5" else "SLOW5",
            "method": "pyslow5 writer, single read group 0",
        },
    }


_UMI_TOOLS_EXTRA = "bio-umi-tools"
_CLUSTER_METHODS = ("adjacency", "cluster", "directional", "unique")
_DEFAULT_MAX_UMIS = 10_000
_MAX_UMIS_LIMIT = 100_000
_MAX_UMI_DISTANCE = 5
_MAX_PATTERN_LENGTH = 1000
_DEFAULT_MAX_BARCODES = 100_000
_MAX_BARCODES_LIMIT = 1_000_000
_MAX_WHITELIST_DISTANCE = 3


def _bounded_int(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _bases(value: object, alphabet: str, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty string of {alphabet} bases.")
    text = value.upper()
    if any(char not in alphabet for char in text):
        raise ConfigurationError(f"{name} must contain only {alphabet} bases; got {value!r}.")
    return text


def _same_length(barcodes: Sequence[str], name: str) -> int:
    lengths = {len(item) for item in barcodes}
    if len(lengths) != 1:
        raise ConfigurationError(f"all {name} must have the same length.")
    return lengths.pop()


def cluster_umi_counts(
    umi_counts: Mapping[str, int],
    *,
    method: str = "directional",
    max_distance: int = 1,
    max_umis: int = _DEFAULT_MAX_UMIS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group UMIs by Hamming distance with umi-tools' ``UMIClusterer``.

    ``method`` is ``"directional"``, ``"adjacency"``, ``"cluster"`` or ``"unique"``.
    Each cluster lists its most frequent UMI first, which is the representative.
    In the directional method a UMI joins a parent only when the parent count is at
    least twice the child count minus one. ``"unique"`` keeps every UMI separate.
    All UMIs must have the same length.
    """

    umi_tools = require("umi_tools", _UMI_TOOLS_EXTRA)
    selected = _choose(method, _CLUSTER_METHODS, "method")
    distance = _bounded_int(max_distance, "max_distance", minimum=0, maximum=_MAX_UMI_DISTANCE)
    limit = _bounded_int(max_umis, "max_umis", minimum=1, maximum=_MAX_UMIS_LIMIT)
    if not isinstance(umi_counts, Mapping) or not umi_counts:
        raise ConfigurationError("umi_counts must be a non-empty mapping from UMI to count.")
    if len(umi_counts) > limit:
        raise ConfigurationError(
            f"umi_counts exceeds max_umis={limit}; increase the explicit limit."
        )
    table: dict[str, int] = {}
    for umi, count in umi_counts.items():
        text = _bases(umi, "ACGTN", "UMI")
        if text in table:
            raise ConfigurationError(f"duplicate UMI after normalization: {text!r}.")
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ConfigurationError(f"count for UMI {text!r} must be a positive integer.")
        table[text] = count
    length = _same_length(list(table), "UMIs")
    clusterer = umi_tools.UMIClusterer(cluster_method=selected)
    raw = clusterer({umi.encode("ascii"): count for umi, count in table.items()}, distance)
    clusters = [[item.decode("ascii") for item in group] for group in raw]
    return {
        "status": "completed",
        "operation": "cluster_umi_counts",
        "method": selected,
        "max_distance": distance,
        "umi_length": length,
        "umis": len(table),
        "clusters": clusters,
        "representatives": [group[0] for group in clusters],
        "cluster_counts": [sum(table[umi] for umi in group) for group in clusters],
        "provenance": {
            "package": "umi-tools",
            "version": _version("umi-tools"),
            "distance": "Hamming distance (umi_tools._dedup_umi.edit_distance)",
            "method": f"umi_tools.UMIClusterer(cluster_method={selected!r})",
        },
    }


def _barcode_pattern(regex: Any, pattern: object) -> Any:
    if (
        not isinstance(pattern, str)
        or not pattern
        or len(pattern) > _MAX_PATTERN_LENGTH
        or "\0" in pattern
    ):
        raise ConfigurationError(
            f"pattern must be a non-empty regular expression of at most {_MAX_PATTERN_LENGTH} "
            "characters."
        )
    try:
        compiled = regex.compile(pattern)
    except regex.error as exc:
        raise ConfigurationError(f"pattern is not a valid regular expression: {exc}") from exc
    if not any(name.startswith(("cell_", "umi_")) for name in compiled.groupindex):
        raise ConfigurationError("pattern must define at least one cell_<n> or umi_<n> group.")
    return compiled


def extract_read_barcodes(
    sequence: str,
    qualities: str,
    pattern: str,
    *,
    discard: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract cell and UMI barcodes from one read with umi-tools' regex extraction.

    ``pattern`` is a regular expression matched from the start of the read. Named
    groups ``cell_<n>`` and ``umi_<n>`` mark barcode bases; ``discard_<n>`` marks bases
    to drop when ``discard`` is true. Groups must take part in every match. The
    returned read has the barcode bases removed. umi-tools does not return cell-barcode
    qualities, so only UMI qualities are reported.
    """

    extract_methods = require("umi_tools.extract_methods", _UMI_TOOLS_EXTRA)
    umi_methods = require("umi_tools.umi_methods", _UMI_TOOLS_EXTRA)
    regex = require("regex", _UMI_TOOLS_EXTRA)
    if not isinstance(discard, bool):
        raise ConfigurationError("discard must be a boolean.")
    text = sequence_text(sequence)
    quals, _offset = _phred_text(qualities, 33)
    if len(quals) != len(text):
        raise ConfigurationError("qualities must contain one character per base.")
    compiled = _barcode_pattern(regex, pattern)
    provenance = {
        "package": "umi-tools",
        "version": _version("umi-tools"),
        "method": "umi_tools.extract_methods.ExtractBarcodes (regex groups)",
        "cell_qualities": "not returned by umi-tools; omitted",
    }
    match = compiled.match(text)
    if match is None:
        return {
            "status": "completed",
            "operation": "extract_read_barcodes",
            "matched": False,
            "cell_barcode": "",
            "umi": "",
            "umi_qualities": "",
            "read_sequence": text,
            "read_qualities": quals,
            "provenance": provenance,
        }
    record = umi_methods.Record("read", text, quals)
    try:
        cell, _cell_quals, umi, umi_quals, new_seq, new_quals = extract_methods.ExtractBarcodes(
            record, match, extract_cell=True, extract_umi=True, discard=discard
        )
    except (AttributeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"barcode groups must participate in every match: {exc}") from exc
    return {
        "status": "completed",
        "operation": "extract_read_barcodes",
        "matched": True,
        "cell_barcode": cell,
        "umi": umi,
        "umi_qualities": umi_quals,
        "read_sequence": new_seq,
        "read_qualities": new_quals,
        "provenance": provenance,
    }


def _barcode_list(value: object, name: str, limit: int) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of DNA barcodes.")
    if len(value) > limit:
        raise ConfigurationError(
            f"{name} exceeds max_barcodes={limit}; increase the explicit limit."
        )
    return [_bases(item, "ACGT", f"{name} barcode") for item in value]


def map_barcodes_to_whitelist(
    observed: Sequence[str],
    whitelist: Sequence[str],
    *,
    threshold: int = 1,
    max_barcodes: int = _DEFAULT_MAX_BARCODES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Correct observed cell barcodes to a whitelist with umi-tools' Hamming mapping.

    An observed barcode that is not whitelisted is assigned only when exactly one
    whitelisted barcode lies within ``threshold`` mismatches. Barcodes with no
    candidate, or with several candidates, are reported in ``unassigned_barcodes``.
    All barcodes must have the same length.
    """

    whitelist_methods = require("umi_tools.whitelist_methods", _UMI_TOOLS_EXTRA)
    distance = _bounded_int(threshold, "threshold", minimum=1, maximum=_MAX_WHITELIST_DISTANCE)
    limit = _bounded_int(max_barcodes, "max_barcodes", minimum=1, maximum=_MAX_BARCODES_LIMIT)
    observed_items = _barcode_list(observed, "observed", limit)
    white = _barcode_list(whitelist, "whitelist", limit)
    length = _same_length(observed_items + white, "barcodes")
    unique_observed = sorted(set(observed_items))
    whitelist_set = set(white)
    mapping = whitelist_methods.getErrorCorrectMapping(
        unique_observed, sorted(whitelist_set), threshold=distance
    )
    corrected = {
        str(key): sorted(str(item) for item in values) for key, values in sorted(mapping.items())
    }
    mapped = {item for values in corrected.values() for item in values}
    return {
        "status": "completed",
        "operation": "map_barcodes_to_whitelist",
        "threshold": distance,
        "barcode_length": length,
        "observed_unique": len(unique_observed),
        "exact_matches": sum(1 for item in unique_observed if item in whitelist_set),
        "corrected": corrected,
        "corrected_barcodes": len(mapped),
        "unassigned_barcodes": [
            item for item in unique_observed if item not in whitelist_set and item not in mapped
        ],
        "provenance": {
            "package": "umi-tools",
            "version": _version("umi-tools"),
            "distance": "Hamming distance (umi_tools.whitelist_methods.getErrorCorrectMapping)",
            "method": "unique candidate within threshold; ambiguous barcodes are not assigned",
        },
    }


__all__ = [
    "calculate_expected_errors",
    "cluster_umi_counts",
    "extract_read_barcodes",
    "find_adapter_match",
    "find_quality_trim_bounds",
    "map_barcodes_to_whitelist",
    "read_slow5_signal",
    "summarize_slow5_signals",
    "trim_fastq_adapters",
    "write_slow5_signals",
]

_register_public_verbose(globals(), tuple(__all__))
