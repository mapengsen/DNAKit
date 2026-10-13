"""Optional-backend sequence records, DNA operations, FASTQ encodings and 2bit access.

The functions in this module are explicit adapters over optional packages:

- ``read_sequence_records()``, ``extract_feature_sequences()`` and
  ``convert_sequence_file()`` use Biopython's ``Bio.SeqIO`` for EMBL and GenBank.
- ``translate_dna_with_genetic_code()``, ``calculate_kmer_frequencies()`` and
  ``expand_degenerate_sequences()`` use scikit-bio DNA objects.
- ``detect_fastq_quality_encoding()`` reads FASTQ quality ranges through pyfastx.
- ``build_screed_database()`` and ``fetch_screed_records()`` write and query a
  screed SQLite record database.
- ``summarize_twobit()`` and ``fetch_twobit_region()`` read UCSC .2bit files
  through py2bit.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import sqlite3
from collections import Counter
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_BIO_EXTRA = "bio-formats"
_SKBIO_EXTRA = "bio-skibio"
_PYFASTX_EXTRA = "bio-pyfastx"
_SCREED_EXTRA = "bio-screed"
_PY2BIT_EXTRA = "bio-py2bit"

_IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_SEQUENCE_FORMATS = ("embl", "genbank")
_OUTPUT_FORMATS = ("fasta", "embl", "genbank")
_TRANSLATION_OPTIONS = ("ignore", "require", "optional")
_READING_FRAMES = (1, 2, 3, -1, -2, -3)
_SCREED_SUFFIX = "_screed"
_MAX_COORDINATE = 4_294_967_295  # 2bit offsets are 32-bit unsigned integers

_DEFAULT_MAX_RECORDS = 1_000
_MAX_RECORDS_LIMIT = 100_000
_DEFAULT_MAX_FEATURES = 10_000
_MAX_FEATURES_LIMIT = 1_000_000
_DEFAULT_MAX_BASES = 10_000_000
_MAX_BASES_LIMIT = 100_000_000
_DEFAULT_MAX_SEQUENCE = 1_000_000
_MAX_SEQUENCE_LIMIT = 10_000_000
_MAX_KMER = 32
_DEFAULT_MAX_EXPANSIONS = 1_024
_MAX_EXPANSIONS_LIMIT = 65_536
_DEFAULT_MAX_FASTQ_BYTES = 1_000_000_000
_MAX_FASTQ_BYTES_LIMIT = 10_000_000_000
_DEFAULT_MAX_CHROMOSOMES = 10_000
_MAX_CHROMOSOMES_LIMIT = 1_000_000


def _provenance(distribution: str, **notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": distribution, "package_version": installed, **notes}


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _member(value: object, allowed: Sequence[object], name: str) -> Any:
    if isinstance(value, bool) or not isinstance(value, type(allowed[0])) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _text(value: object, name: str, *, maximum: int) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > maximum:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    if "\x00" in value:
        raise ConfigurationError(f"{name} cannot contain NUL characters.")
    return value


def _dna(value: object, *, max_length: int) -> str:
    text = _text(value, "sequence", maximum=max_length).upper()
    invalid = sorted(set(text) - _IUPAC_DNA)
    if invalid:
        raise ConfigurationError(
            f"sequence contains unsupported symbols {invalid}; "
            "gaps and non-DNA characters are rejected."
        )
    return text


def _source(value: object) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError("source must be a file path.")
    _text(str(value), "source", maximum=4096)
    return input_path(value)


def _new_output(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    _text(str(value), name, maximum=4096)
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory for {name} does not exist: {path.parent}")
    return path


def _screed_target(value: object) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError("database must be a file path.")
    _text(str(value), "database", maximum=4096)
    path = Path(value).expanduser().resolve()
    if not str(path).endswith(_SCREED_SUFFIX):
        path = path.with_name(path.name + _SCREED_SUFFIX)
    return path


def _feature_type(value: object) -> str:
    text = _text(value, "feature_type", maximum=64)
    if any(character.isspace() for character in text):
        raise ConfigurationError("feature_type cannot contain whitespace.")
    return text


def _first(qualifiers: Mapping[str, Any], key: str) -> str | None:
    values = qualifiers.get(key)
    return str(values[0]) if values else None


def read_sequence_records(
    source: str | Path,
    *,
    format: str = "embl",
    max_records: int = _DEFAULT_MAX_RECORDS,
    max_bases: int = _DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse EMBL or GenBank records with Biopython, bounded by record and base counts."""

    path = _source(source)
    selected = _member(format, _SEQUENCE_FORMATS, "format")
    record_limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    base_limit = _integer(max_bases, "max_bases", minimum=1, maximum=_MAX_BASES_LIMIT)
    seqio = require("Bio.SeqIO", _BIO_EXTRA)
    records: list[dict[str, Any]] = []
    total_bases = 0
    try:
        for record in seqio.parse(str(path), selected):
            if len(records) >= record_limit:
                raise ConfigurationError(
                    f"{path.name} exceeds max_records={record_limit}; increase the explicit limit."
                )
            sequence = str(record.seq).upper()
            total_bases += len(sequence)
            if total_bases > base_limit:
                raise ConfigurationError(
                    f"{path.name} exceeds max_bases={base_limit}; increase the explicit limit."
                )
            records.append(
                {
                    "id": str(record.id),
                    "name": str(record.name),
                    "description": str(record.description),
                    "length": len(sequence),
                    "feature_counts": dict(
                        sorted(Counter(str(f.type) for f in record.features).items())
                    ),
                    "sequence": sequence,
                }
            )
    except ConfigurationError:
        raise
    except (OSError, TypeError, ValueError) as exc:
        raise ConfigurationError(
            f"Biopython could not parse {path.name} as {selected}: {exc}"
        ) from exc
    return {
        "status": "completed",
        "operation": "read_sequence_records",
        "source": str(path),
        "format": selected,
        "record_count": len(records),
        "total_bases": total_bases,
        "records": records,
        "provenance": _provenance(
            "biopython",
            parser="Bio.SeqIO.parse",
            coordinate_system="feature locations are 1-based and inclusive",
        ),
    }


def extract_feature_sequences(
    source: str | Path,
    *,
    format: str = "genbank",
    feature_type: str = "CDS",
    max_records: int = _DEFAULT_MAX_RECORDS,
    max_features: int = _DEFAULT_MAX_FEATURES,
    max_bases: int = _DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract sequences of one feature type from EMBL or GenBank records through Biopython."""

    path = _source(source)
    selected = _member(format, _SEQUENCE_FORMATS, "format")
    kind = _feature_type(feature_type)
    record_limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    feature_limit = _integer(max_features, "max_features", minimum=1, maximum=_MAX_FEATURES_LIMIT)
    base_limit = _integer(max_bases, "max_bases", minimum=1, maximum=_MAX_BASES_LIMIT)
    seqio = require("Bio.SeqIO", _BIO_EXTRA)
    features: list[dict[str, Any]] = []
    records_seen = 0
    total_bases = 0
    try:
        for record in seqio.parse(str(path), selected):
            records_seen += 1
            if records_seen > record_limit:
                raise ConfigurationError(
                    f"{path.name} exceeds max_records={record_limit}; increase the explicit limit."
                )
            for index, feature in enumerate(record.features, start=1):
                if str(feature.type) != kind:
                    continue
                if len(features) >= feature_limit:
                    raise ConfigurationError(
                        f"{path.name} exceeds max_features={feature_limit}; "
                        "increase the explicit limit."
                    )
                sequence = str(feature.extract(record.seq)).upper()
                total_bases += len(sequence)
                if total_bases > base_limit:
                    raise ConfigurationError(
                        f"{path.name} exceeds max_bases={base_limit}; increase the explicit limit."
                    )
                location = feature.location
                features.append(
                    {
                        "record_id": str(record.id),
                        "feature_index": index,
                        "type": kind,
                        "gene": _first(feature.qualifiers, "gene"),
                        "locus_tag": _first(feature.qualifiers, "locus_tag"),
                        "start": int(location.start) + 1,
                        "end": int(location.end),
                        "strand": location.strand,
                        "parts": len(location.parts),
                        "length": len(sequence),
                        "sequence": sequence,
                    }
                )
    except ConfigurationError:
        raise
    except (OSError, TypeError, ValueError) as exc:
        raise ConfigurationError(
            f"Biopython could not extract {kind} features from {path.name}: {exc}"
        ) from exc
    return {
        "status": "completed",
        "operation": "extract_feature_sequences",
        "source": str(path),
        "format": selected,
        "feature_type": kind,
        "feature_count": len(features),
        "features": features,
        "provenance": _provenance(
            "biopython",
            extraction="Bio.SeqFeature.extract; minus-strand and join() locations "
            "are reverse-complemented",
            coordinate_system="start and end span the whole feature, 1-based and inclusive",
        ),
    }


def convert_sequence_file(
    source: str | Path,
    target: str | Path,
    *,
    input_format: str = "genbank",
    output_format: str = "fasta",
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert EMBL or GenBank records to FASTA, EMBL or GenBank through Biopython."""

    path = _source(source)
    in_format = _member(input_format, _SEQUENCE_FORMATS, "input_format")
    out_format = _member(output_format, _OUTPUT_FORMATS, "output_format")
    record_limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    destination = _new_output(target, "target")
    seqio = require("Bio.SeqIO", _BIO_EXTRA)
    seen = 0

    def bounded() -> Iterator[Any]:
        nonlocal seen
        for record in seqio.parse(str(path), in_format):
            if seen >= record_limit:
                raise ConfigurationError(
                    f"{path.name} exceeds max_records={record_limit}; increase the explicit limit."
                )
            seen += 1
            yield record

    try:
        written = seqio.write(bounded(), str(destination), out_format)
    except ConfigurationError:
        destination.unlink(missing_ok=True)
        raise
    except (OSError, TypeError, ValueError) as exc:
        destination.unlink(missing_ok=True)
        raise ConfigurationError(f"Biopython could not convert {path.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "convert_sequence_file",
        "source": str(path),
        "target": str(destination),
        "input_format": in_format,
        "output_format": out_format,
        "record_count": written,
        "provenance": _provenance(
            "biopython",
            writer="Bio.SeqIO.write",
            note="FASTA output keeps sequences only; features need EMBL or GenBank",
        ),
    }


def translate_dna_with_genetic_code(
    sequence: str,
    *,
    genetic_code: int = 1,
    reading_frame: int = 1,
    start: str = "ignore",
    stop: str = "ignore",
    max_length: int = _DEFAULT_MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Translate DNA with a selected NCBI genetic code through scikit-bio."""

    length_limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_SEQUENCE_LIMIT)
    text = _dna(sequence, max_length=length_limit)
    table = _integer(genetic_code, "genetic_code", minimum=1, maximum=33)
    frame = _member(reading_frame, _READING_FRAMES, "reading_frame")
    start_policy = _member(start, _TRANSLATION_OPTIONS, "start")
    stop_policy = _member(stop, _TRANSLATION_OPTIONS, "stop")
    skbio = require("skbio", _SKBIO_EXTRA)
    try:
        protein = skbio.DNA(text).translate(
            genetic_code=table, reading_frame=frame, start=start_policy, stop=stop_policy
        )
    except ValueError as exc:
        raise ConfigurationError(f"scikit-bio could not translate the sequence: {exc}") from exc
    protein_text = str(protein)
    return {
        "status": "completed",
        "operation": "translate_dna_with_genetic_code",
        "sequence_length": len(text),
        "protein": protein_text,
        "protein_length": len(protein_text),
        "genetic_code": {
            "table_id": table,
            "name": skbio.sequence.GeneticCode.from_ncbi(table).name,
        },
        "reading_frame": frame,
        "start": start_policy,
        "stop": stop_policy,
        "provenance": _provenance(
            "scikit-bio",
            translation="skbio.DNA.translate with NCBI genetic code table IDs",
            coordinate_system="reading frames count from the 5' end; "
            "negative frames use the reverse complement",
            gaps="rejected before translation",
        ),
    }


def calculate_kmer_frequencies(
    sequence: str,
    *,
    k: int = 3,
    overlap: bool = True,
    relative: bool = False,
    max_length: int = _DEFAULT_MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count k-mers, optionally non-overlapping or as relative frequencies, through scikit-bio."""

    length_limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_SEQUENCE_LIMIT)
    text = _dna(sequence, max_length=length_limit)
    word = _integer(k, "k", minimum=1, maximum=_MAX_KMER)
    use_overlap = _flag(overlap, "overlap")
    use_relative = _flag(relative, "relative")
    skbio = require("skbio", _SKBIO_EXTRA)
    counts = skbio.DNA(text).kmer_frequencies(word, overlap=use_overlap, relative=use_relative)
    kmers = {
        str(key): (float(value) if use_relative else int(value))
        for key, value in sorted(counts.items())
    }
    return {
        "status": "completed",
        "operation": "calculate_kmer_frequencies",
        "sequence_length": len(text),
        "k": word,
        "overlap": use_overlap,
        "relative": use_relative,
        "distinct_kmers": len(kmers),
        "kmers": kmers,
        "provenance": _provenance(
            "scikit-bio",
            counting="k-mers containing N or other IUPAC symbols are counted as written",
            relative_denominator="number of k-mers counted under the chosen overlap setting",
        ),
    }


def expand_degenerate_sequences(
    sequence: str,
    *,
    max_sequences: int = _DEFAULT_MAX_EXPANSIONS,
    max_length: int = _DEFAULT_MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate every unambiguous DNA sequence an IUPAC-coded sequence represents."""

    length_limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_SEQUENCE_LIMIT)
    text = _dna(sequence, max_length=length_limit)
    limit = _integer(max_sequences, "max_sequences", minimum=1, maximum=_MAX_EXPANSIONS_LIMIT)
    skbio = require("skbio", _SKBIO_EXTRA)
    choices = skbio.DNA.degenerate_map
    expected = 1
    for symbol in text:
        expected *= len(choices.get(symbol, (symbol,)))
        if expected > limit:  # stop early so long all-degenerate inputs stay cheap
            raise ConfigurationError(
                f"expansion would produce more than {limit} sequences; "
                "raise max_sequences to allow it."
            )
    sequences = sorted(str(item) for item in skbio.DNA(text).expand_degenerates())
    return {
        "status": "completed",
        "operation": "expand_degenerate_sequences",
        "sequence_length": len(text),
        "degenerate_positions": sum(1 for symbol in text if symbol not in "ACGT"),
        "expansion_count": len(sequences),
        "sequences": sequences,
        "provenance": _provenance(
            "scikit-bio",
            method="skbio.DNA.expand_degenerates",
            ordering="lexicographic; scikit-bio does not guarantee an order",
        ),
    }


def detect_fastq_quality_encoding(
    source: str | Path,
    *,
    max_bytes: int = _DEFAULT_MAX_FASTQ_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Infer the FASTQ quality offset and quality range through pyfastx."""

    path = _source(source)
    byte_limit = _integer(max_bytes, "max_bytes", minimum=1, maximum=_MAX_FASTQ_BYTES_LIMIT)
    size = path.stat().st_size
    if size > byte_limit:
        raise ConfigurationError(
            f"{path.name} has {size} bytes, above max_bytes={byte_limit}; "
            "increase the explicit limit."
        )
    pyfastx = require("pyfastx", _PYFASTX_EXTRA)
    try:
        reader = pyfastx.Fastq(str(path), build_index=False)
        offset = int(reader.phred)
        low = int(reader.minqual)
        high = int(reader.maxqual)
        min_length = int(reader.minlen)
        max_length = int(reader.maxlen)
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pyfastx could not read {path.name} as FASTQ: {exc}") from exc
    if offset not in (0, 33, 64) or low > high:
        raise ConfigurationError(f"pyfastx reported an inconsistent quality range for {path.name}.")
    if low < 33 or (offset and low < offset):
        raise ConfigurationError(
            f"{path.name} has quality characters below the detected offset {offset}."
        )
    determined = offset != 0
    return {
        "status": "completed",
        "operation": "detect_fastq_quality_encoding",
        "source": str(path),
        "phred_offset": offset if determined else None,
        "possible_phred_offsets": [offset] if determined else [33, 64],
        "min_quality_ascii": low,
        "max_quality_ascii": high,
        "min_quality_score": low - offset if determined else None,
        "max_quality_score": high - offset if determined else None,
        "min_length": min_length,
        "max_length": max_length,
        "provenance": _provenance(
            "pyfastx",
            method="pyfastx.Fastq scan without an index",
            ambiguity="pyfastx reports offset 0 when both 33 and 64 fit the observed characters",
        ),
    }


def _close(handle: Any) -> None:
    close = getattr(handle, "close", None)
    if callable(close):
        close()


def build_screed_database(
    source: str | Path,
    *,
    database: str | Path,
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a new screed SQLite record database from one FASTA or FASTQ file."""

    path = _source(source)
    record_limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    target = _screed_target(database)
    if target.exists():
        raise ConfigurationError(f"database already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The directory for database does not exist: {target.parent}")
    screed = require("screed", _SCREED_EXTRA)
    try:
        reader = screed.open(str(path), parse_description=True)
    except (OSError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"screed could not open {path.name}: {exc}") from exc
    kind = getattr(getattr(reader, "iter_fn", None), "__name__", "")
    field_types = {
        "fasta_iter": screed.fasta.FieldTypes,
        "fastq_iter": screed.fastq.FieldTypes,
    }.get(kind)
    if field_types is None:
        _close(reader)
        raise ConfigurationError(f"screed could not identify {path.name} as FASTA or FASTQ.")
    count = 0

    def bounded() -> Iterator[Any]:
        nonlocal count
        for record in reader:
            if count >= record_limit:
                raise ConfigurationError(
                    f"{path.name} exceeds max_records={record_limit}; increase the explicit limit."
                )
            count += 1
            yield record

    try:
        screed.create_db(str(target), field_types, bounded())
    except ConfigurationError:
        target.unlink(missing_ok=True)
        raise
    except (OSError, TypeError, ValueError, sqlite3.Error) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(
            f"screed could not build a database from {path.name}: {exc}"
        ) from exc
    finally:
        _close(reader)
    return {
        "status": "completed",
        "operation": "build_screed_database",
        "source": str(path),
        "database": str(target),
        "format": "fasta" if kind == "fasta_iter" else "fastq",
        "record_count": count,
        "provenance": _provenance(
            "screed",
            storage="SQLite file written by screed.create_db",
            lookup="record name",
        ),
    }


def _screed_record(entry: Mapping[str, Any]) -> dict[str, Any]:
    sequence = str(entry["sequence"]).upper()
    record: dict[str, Any] = {
        "name": str(entry["name"]),
        "description": str(entry.get("description") or entry.get("annotations") or ""),
        "length": len(sequence),
        "sequence": sequence,
    }
    if "quality" in entry:
        record["quality"] = str(entry["quality"])
    return record


def fetch_screed_records(
    database: str | Path,
    *,
    names: Sequence[str] | None = None,
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read named records, or all records up to a limit, from a screed database."""

    target = _screed_target(database)
    if not target.is_file():
        raise ConfigurationError(f"screed database does not exist: {target}")
    record_limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS_LIMIT)
    requested: list[str] | None = None
    if names is not None:
        if isinstance(names, (str, bytes)) or not isinstance(names, Sequence) or not names:
            raise ConfigurationError("names must be a non-empty list of record names.")
        if len(names) > record_limit:
            raise ConfigurationError(f"names has more than max_records={record_limit} entries.")
        requested = [_text(name, "record name", maximum=512) for name in names]
    screed = require("screed", _SCREED_EXTRA)
    try:
        handle = screed.ScreedDB(str(target))
    except (OSError, ValueError, sqlite3.Error) as exc:
        raise ConfigurationError(f"screed could not open {target.name}: {exc}") from exc
    records: list[dict[str, Any]] = []
    missing: list[str] = []
    try:
        if requested is None:
            requested = list(map(str, handle.keys()))
            if len(requested) > record_limit:
                raise ConfigurationError(
                    f"database holds {len(requested)} records, above max_records={record_limit}; "
                    "pass names or raise the limit."
                )
        for key in requested:
            entry = handle.get(key)
            if entry is None:
                missing.append(key)
                continue
            records.append(_screed_record(entry))
    except (OSError, TypeError, ValueError, sqlite3.Error) as exc:
        raise ConfigurationError(f"screed could not read {target.name}: {exc}") from exc
    finally:
        _close(handle)
    if missing:
        raise ConfigurationError(f"screed database has no records named: {missing[:10]}")
    return {
        "status": "completed",
        "operation": "fetch_screed_records",
        "database": str(target),
        "record_count": len(records),
        "records": records,
        "provenance": _provenance("screed", lookup="SQLite key lookup by record name"),
    }


def _twobit_open(py2bit: Any, path: Path) -> Any:
    try:
        return py2bit.open(str(path))
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"py2bit could not open {path.name} as 2bit: {exc}") from exc


def summarize_twobit(
    source: str | Path,
    *,
    max_chromosomes: int = _DEFAULT_MAX_CHROMOSOMES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize chromosome lengths, hard-masked N bases and ACGT composition of a 2bit file."""

    path = _source(source)
    limit = _integer(max_chromosomes, "max_chromosomes", minimum=1, maximum=_MAX_CHROMOSOMES_LIMIT)
    py2bit = require("py2bit", _PY2BIT_EXTRA)
    handle = _twobit_open(py2bit, path)
    chromosomes: list[dict[str, Any]] = []
    try:
        lengths = handle.chroms()
        if len(lengths) > limit:
            raise ConfigurationError(
                f"{path.name} has {len(lengths)} sequences, above max_chromosomes={limit}."
            )
        for name in sorted(lengths):
            counts = handle.bases(name, fraction=False)
            acgt = sum(int(counts.get(base, 0)) for base in "ACGT")
            gc = int(counts.get("G", 0)) + int(counts.get("C", 0))
            hard = sum(int(end) - int(start) for start, end in handle.hardMaskedBlocks(name))
            chromosomes.append(
                {
                    "name": str(name),
                    "length": int(lengths[name]),
                    "acgt_bases": acgt,
                    "hard_masked_bases": hard,
                    "gc_percent_of_acgt": (100.0 * gc / acgt) if acgt else None,
                }
            )
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"py2bit could not summarize {path.name}: {exc}") from exc
    finally:
        _close(handle)
    return {
        "status": "completed",
        "operation": "summarize_twobit",
        "source": str(path),
        "chromosome_count": len(chromosomes),
        "total_length": sum(item["length"] for item in chromosomes),
        "total_hard_masked": sum(item["hard_masked_bases"] for item in chromosomes),
        "chromosomes": chromosomes,
        "provenance": _provenance(
            "py2bit",
            masking="soft-masked lowercase is ignored; N bases are hard-masked "
            "and excluded from ACGT counts",
            coordinate_system="chromosome lengths are sequence lengths",
        ),
    }


def fetch_twobit_region(
    source: str | Path,
    chromosome: str,
    start: int,
    end: int,
    *,
    max_length: int = _DEFAULT_MAX_SEQUENCE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read one 0-based, half-open interval from a 2bit file through py2bit."""

    path = _source(source)
    name = _text(chromosome, "chromosome", maximum=512)
    begin = _integer(start, "start", minimum=0, maximum=_MAX_COORDINATE)
    finish = _integer(end, "end", minimum=begin + 1, maximum=_MAX_COORDINATE)
    limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_SEQUENCE_LIMIT)
    if finish - begin > limit:
        raise ConfigurationError(
            f"requested interval has {finish - begin} bases, above max_length={limit}."
        )
    py2bit = require("py2bit", _PY2BIT_EXTRA)
    handle = _twobit_open(py2bit, path)
    try:
        length = handle.chroms(name)
        if length is None:
            raise ConfigurationError(f"chromosome {name!r} is not present in {path.name}.")
        if finish > int(length):
            raise ConfigurationError(
                f"end={finish} is beyond {name!r}, whose length is {int(length)}."
            )
        sequence = str(handle.sequence(name, begin, finish)).upper()
    except ConfigurationError:
        raise
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"py2bit could not read {name!r} from {path.name}: {exc}") from exc
    finally:
        _close(handle)
    return {
        "status": "completed",
        "operation": "fetch_twobit_region",
        "source": str(path),
        "chromosome": name,
        "start": begin,
        "end": finish,
        "length": len(sequence),
        "hard_masked_bases": sequence.count("N"),
        "sequence": sequence,
        "provenance": _provenance(
            "py2bit",
            coordinate_system="0-based start, 0-based exclusive end",
            masking="soft-masked lowercase is ignored; output is uppercase",
        ),
    }


__all__ = [
    "build_screed_database",
    "calculate_kmer_frequencies",
    "convert_sequence_file",
    "detect_fastq_quality_encoding",
    "expand_degenerate_sequences",
    "extract_feature_sequences",
    "fetch_screed_records",
    "fetch_twobit_region",
    "read_sequence_records",
    "summarize_twobit",
    "translate_dna_with_genetic_code",
]

_register_public_verbose(globals(), tuple(__all__))
