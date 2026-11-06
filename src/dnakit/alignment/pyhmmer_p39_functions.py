"""pyhmmer (HMMER 3 and Easel) operations on DNA that were still missing (integration pack P39).

The functions work on local files and return JSON-serialisable dictionaries; every written file is a new path.
They cover the HMMER pipelines ``hmmsearch`` and ``hmmscan`` run with DNA profiles, building a profile from one
sequence, emitting sequences and alignments from a profile, and the Easel alignment tools (reverse complement,
sequence weights, identity filtering, alphabet detection) and SSI sequence indexes.

pyhmmer is imported only when a function runs; a missing package raises ``BackendUnavailableError`` with the
``pip install "dnakit[bio-pyhmmer]"`` hint. Tested with pyhmmer 0.12.3.
"""

from __future__ import annotations

import importlib.metadata
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-pyhmmer"
_MAX_SEQUENCES = 1_000_000
_MAX_EMIT = 100_000
_MSA_FORMATS = (
    "stockholm",
    "pfam",
    "a2m",
    "psiblast",
    "selex",
    "afa",
    "clustal",
    "clustallike",
    "phylip",
    "phylips",
)
_WEIGHT_METHODS = ("pb", "gsc", "blosum")
_PREFERENCES = ("conscover", "origorder", "random")
_SEQUENCE_KINDS = ("sequence", "alignment")
_HEADER = re.compile(rb"^>([^\s]+)[^\n]*\n", re.MULTILINE)


def _version() -> str:
    try:
        return importlib.metadata.version("pyhmmer")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


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


def _number(value: object, name: str, low: float, high: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not low <= value <= high:
        raise ConfigurationError(f"{name} must be a number in [{low}, {high}].")
    return float(value)


def _integer(value: object, name: str, low: int, high: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ConfigurationError(f"{name} must be an integer in [{low}, {high}].")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _text(value: bytes | str | None) -> str | None:
    if value is None:
        return None
    return value.decode("ascii", "replace") if isinstance(value, bytes) else str(value)


def _dna_sequences(pyhmmer: Any, path: Path, alphabet: Any) -> Any:
    try:
        with pyhmmer.easel.SequenceFile(str(path), digital=True, alphabet=alphabet) as handle:
            block = handle.read_block(_MAX_SEQUENCES + 1)
    except (OSError, ValueError, EOFError) as exc:
        raise ConfigurationError(f"{path.name} is not a readable sequence file: {exc}") from exc
    if len(block) == 0:
        raise ConfigurationError(f"{path.name} holds no sequences.")
    if len(block) > _MAX_SEQUENCES:
        raise ConfigurationError(f"{path.name} holds more than {_MAX_SEQUENCES} sequences.")
    return block


def _read_hmms(pyhmmer: Any, path: Path) -> list[Any]:
    try:
        with pyhmmer.plan7.HMMFile(str(path)) as handle:
            hmms = list(handle)
    except (OSError, ValueError, EOFError) as exc:
        raise ConfigurationError(
            f"{path.name} is not a readable HMMER profile file: {exc}"
        ) from exc
    if not hmms:
        raise ConfigurationError(f"{path.name} holds no profile HMMs.")
    return hmms


def _read_msa(pyhmmer: Any, path: Path, alphabet: Any, file_format: str | None) -> Any:
    try:
        with pyhmmer.easel.MSAFile(
            str(path), format=file_format, digital=True, alphabet=alphabet
        ) as handle:
            return handle.read()
    except (OSError, ValueError, EOFError) as exc:
        raise ConfigurationError(f"{path.name} is not a readable alignment file: {exc}") from exc


def _write_msa(pyhmmer: Any, msa: Any, path: Path, output_format: str) -> None:
    with path.open("xb") as handle:
        msa.write(handle, output_format)


def _hit_rows(top_hits: Any, max_hits: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for hit in top_hits:
        if len(rows) >= max_hits:
            break
        domains = []
        for domain in hit.domains:
            domains.append(
                {
                    "included": bool(domain.included),
                    "score": float(domain.score),
                    "bias": float(domain.bias),
                    "c_evalue": float(domain.c_evalue),
                    "i_evalue": float(domain.i_evalue),
                    "env_from": int(domain.env_from),
                    "env_to": int(domain.env_to),
                    "ali_from": int(domain.alignment.target_from),
                    "ali_to": int(domain.alignment.target_to),
                    "hmm_from": int(domain.alignment.hmm_from),
                    "hmm_to": int(domain.alignment.hmm_to),
                }
            )
        rows.append(
            {
                "name": _text(hit.name),
                "accession": _text(hit.accession),
                "description": _text(hit.description),
                "score": float(hit.score),
                "bias": float(hit.bias),
                "evalue": float(hit.evalue),
                "included": bool(hit.included),
                "reported": bool(hit.reported),
                "domains": domains,
            }
        )
    return rows


def _thresholds(evalue: float, domain_evalue: float, bit_cutoffs: str | None) -> dict[str, Any]:
    options: dict[str, Any] = {
        "E": _number(evalue, "evalue", 1e-300, 1e9),
        "domE": _number(domain_evalue, "domain_evalue", 1e-300, 1e9),
    }
    if bit_cutoffs is not None:
        options["bit_cutoffs"] = _choice(
            bit_cutoffs, ("gathering", "noise", "trusted"), "bit_cutoffs"
        )
    return options


def pyhmmer_hmmsearch(
    profiles: str,
    sequences: str,
    *,
    evalue: float = 10.0,
    domain_evalue: float = 10.0,
    bit_cutoffs: str | None = None,
    max_hits: int = 1000,
    cpus: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a DNA sequence file with one or more DNA profile HMMs using the HMMER 3 hmmsearch pipeline (pyhmmer.hmmsearch).

    ``profiles`` is a HMMER3 ``.hmm`` file and ``sequences`` a FASTA file read in the DNA alphabet. The result has one
    entry per profile with its hits and domains (HMMER reporting thresholds ``evalue`` 10 and ``domain_evalue`` 10 by
    default; ``bit_cutoffs`` may be ``gathering``, ``noise`` or ``trusted``). This is the profile-against-database
    direction (hmmsearch, Pipeline.search_hmm); ``pyhmmer_hmmscan`` is the sequence-against-profiles direction.
    """

    options = _thresholds(evalue, domain_evalue, bit_cutoffs)
    limit = _integer(max_hits, "max_hits", 1, 1_000_000)
    workers = _integer(cpus, "cpus", 0, 256)
    profile_path = _existing(profiles, "profiles")
    sequence_path = _existing(sequences, "sequences")
    pyhmmer = require("pyhmmer", _EXTRA)
    alphabet = pyhmmer.easel.Alphabet.dna()
    hmms = _read_hmms(pyhmmer, profile_path)
    if not all(hmm.alphabet.is_dna() for hmm in hmms):
        raise ConfigurationError("pyhmmer_hmmsearch needs DNA profile HMMs.")
    targets = _dna_sequences(pyhmmer, sequence_path, alphabet)
    results = []
    try:
        for top_hits in pyhmmer.hmmsearch(hmms, targets, cpus=workers, **options):
            results.append(
                {
                    "query": _text(top_hits.query.name),
                    "targets_searched": int(top_hits.Z),
                    "hit_count": len(top_hits),
                    "included_count": int(sum(1 for hit in top_hits if hit.included)),
                    "hits": _hit_rows(top_hits, limit),
                }
            )
    except (RuntimeError, ValueError) as exc:
        raise ConfigurationError(f"pyhmmer hmmsearch failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "pyhmmer_hmmsearch",
        "tool_version": _version(),
        "target_count": len(targets),
        "results": results,
    }


def pyhmmer_hmmscan(
    sequences: str,
    profiles: str,
    *,
    evalue: float = 10.0,
    domain_evalue: float = 10.0,
    bit_cutoffs: str | None = None,
    max_hits: int = 1000,
    cpus: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan each DNA sequence of a FASTA file against a database of DNA profile HMMs with the HMMER 3 hmmscan pipeline (pyhmmer.hmmscan).

    The result has one entry per query sequence listing the profiles it matches (name, score, E-value and domains).
    Thresholds and ``bit_cutoffs`` are those of ``pyhmmer_hmmsearch``.
    """

    options = _thresholds(evalue, domain_evalue, bit_cutoffs)
    limit = _integer(max_hits, "max_hits", 1, 1_000_000)
    workers = _integer(cpus, "cpus", 0, 256)
    sequence_path = _existing(sequences, "sequences")
    profile_path = _existing(profiles, "profiles")
    pyhmmer = require("pyhmmer", _EXTRA)
    alphabet = pyhmmer.easel.Alphabet.dna()
    hmms = _read_hmms(pyhmmer, profile_path)
    if not all(hmm.alphabet.is_dna() for hmm in hmms):
        raise ConfigurationError("pyhmmer_hmmscan needs DNA profile HMMs.")
    queries = _dna_sequences(pyhmmer, sequence_path, alphabet)
    results = []
    try:
        for top_hits in pyhmmer.hmmscan(queries, hmms, cpus=workers, **options):
            results.append(
                {
                    "query": _text(top_hits.query.name),
                    "profiles_searched": int(top_hits.Z),
                    "hit_count": len(top_hits),
                    "included_count": int(sum(1 for hit in top_hits if hit.included)),
                    "hits": _hit_rows(top_hits, limit),
                }
            )
    except (RuntimeError, ValueError) as exc:
        raise ConfigurationError(f"pyhmmer hmmscan failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "pyhmmer_hmmscan",
        "tool_version": _version(),
        "query_count": len(queries),
        "results": results,
    }


def build_hmm_from_sequence(
    sequence: str,
    output_path: str,
    *,
    name: str = "query",
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a DNA profile HMM from one sequence and write it as a HMMER3 file (pyhmmer.plan7.Builder.build).

    The model has one match state per base (HMMER defaults: DNA alphabet, default score system and background).
    Use ``dnakit.alignment.build_dna_hmm`` to build a model from an alignment.
    """

    if not isinstance(sequence, str) or not sequence or len(sequence) > 1_000_000:
        raise ConfigurationError("sequence must be a DNA string of 1 to 1,000,000 bases.")
    if not isinstance(name, str) or not name or any(c.isspace() for c in name):
        raise ConfigurationError("name must be a non-empty string without spaces.")
    target = _new(output_path, "output_path")
    pyhmmer = require("pyhmmer", _EXTRA)
    alphabet = pyhmmer.easel.Alphabet.dna()
    try:
        digital = pyhmmer.easel.TextSequence(
            name=name.encode("ascii"), sequence=sequence.upper()
        ).digitize(alphabet)
        builder = pyhmmer.plan7.Builder(alphabet)
        hmm, _profile, _optimized = builder.build(digital, pyhmmer.plan7.Background(alphabet))
        with target.open("wb") as handle:
            hmm.write(handle)
    except (UnicodeEncodeError, ValueError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"pyhmmer could not build the profile: {exc}") from exc
    return {
        "status": "completed",
        "operation": "build_hmm_from_sequence",
        "tool_version": _version(),
        "name": name,
        "model_length": int(hmm.M),
        "output": str(target),
    }


def emit_hmm_sequences(
    profile: str,
    count: int,
    output_fasta: str,
    *,
    source: str = "core",
    seed: int = 42,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sample sequences from a DNA profile HMM and write them to a new FASTA file (HMM.emit_sequence, Profile.emit_sequence).

    ``source="core"`` emits from the core model (hmmemit default); ``source="profile"`` emits from the implicit search
    profile (``Profile.emit_sequence``). A fixed ``seed`` makes the sampling reproducible. Sequences are named
    ``<model>-sample<N>``.
    """

    n = _integer(count, "count", 1, _MAX_EMIT)
    seed_value = _integer(seed, "seed", 1, 2**31 - 1)
    kind = _choice(source, ("core", "profile"), "source")
    hmm_path = _existing(profile, "profile")
    target = _new(output_fasta, "output_fasta")
    pyhmmer = require("pyhmmer", _EXTRA)
    hmm = _read_hmms(pyhmmer, hmm_path)[0]
    if not hmm.alphabet.is_dna():
        raise ConfigurationError("emit_hmm_sequences needs a DNA profile HMM.")
    rng = pyhmmer.easel.Randomness(seed_value)
    emitter = hmm
    if kind == "profile":
        background = pyhmmer.plan7.Background(hmm.alphabet)
        emitter = pyhmmer.plan7.Profile(hmm.M, hmm.alphabet)
        emitter.configure(hmm, background)
    model = _text(hmm.name) or "hmm"
    rows = []
    with target.open("w", encoding="ascii") as handle:
        for index in range(1, n + 1):
            emitted = (
                emitter.emit_sequence(rng)
                if kind == "core"
                else emitter.emit_sequence(hmm, randomness=rng)
            )
            text = hmm.alphabet.decode(emitted.sequence)
            rows.append(len(text))
            handle.write(f">{model}-sample{index}\n{text}\n")
    return {
        "status": "completed",
        "operation": "emit_hmm_sequences",
        "tool_version": _version(),
        "model": model,
        "source": kind,
        "seed": seed_value,
        "sequences": n,
        "lengths": rows,
        "output": str(target),
    }


def emit_hmm_alignment(
    profile: str,
    count: int,
    output_alignment: str,
    *,
    output_format: str = "afa",
    seed: int = 42,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sample a multiple alignment of sequences emitted from the core model of a DNA profile HMM (HMM.emit_alignment).

    ``count`` rows are sampled; the alignment is written in ``output_format`` (stockholm, pfam, a2m, psiblast, selex,
    afa, clustal, clustallike, phylip or phylips) to a new file.
    """

    n = _integer(count, "count", 1, _MAX_EMIT)
    seed_value = _integer(seed, "seed", 1, 2**31 - 1)
    fmt = _choice(output_format, _MSA_FORMATS, "output_format")
    hmm_path = _existing(profile, "profile")
    target = _new(output_alignment, "output_alignment")
    pyhmmer = require("pyhmmer", _EXTRA)
    hmm = _read_hmms(pyhmmer, hmm_path)[0]
    if not hmm.alphabet.is_dna():
        raise ConfigurationError("emit_hmm_alignment needs a DNA profile HMM.")
    emitted = hmm.emit_alignment(n, pyhmmer.easel.Randomness(seed_value))
    model = _text(hmm.name) or "hmm"
    rows = [
        pyhmmer.easel.TextSequence(name=f"{model}-sample{index}".encode("ascii"), sequence=row)
        for index, row in enumerate(emitted.textize().alignment, start=1)
    ]
    msa = pyhmmer.easel.TextMSA(name=hmm.name, sequences=rows).digitize(hmm.alphabet)
    _write_msa(pyhmmer, msa, target, fmt)
    return {
        "status": "completed",
        "operation": "emit_hmm_alignment",
        "tool_version": _version(),
        "model": _text(hmm.name),
        "sequences": len(msa.sequences),
        "columns": len(msa.alignment[0]) if len(msa.alignment) else 0,
        "format": fmt,
        "output": str(target),
    }


def reverse_complement_alignment(
    alignment: str,
    output_alignment: str,
    *,
    input_format: str | None = None,
    output_format: str = "afa",
    verbose: bool = True,
) -> dict[str, Any]:
    """Reverse-complement a DNA multiple alignment, including per-column annotation, and write it to a new file (DigitalMSA.reverse_complement).

    ``input_format`` is detected when omitted. Sequence names are kept.
    """

    fmt = _choice(output_format, _MSA_FORMATS, "output_format")
    source = _existing(alignment, "alignment")
    target = _new(output_alignment, "output_alignment")
    pyhmmer = require("pyhmmer", _EXTRA)
    msa = _read_msa(pyhmmer, source, pyhmmer.easel.Alphabet.dna(), input_format)
    flipped = msa.reverse_complement()
    _write_msa(pyhmmer, flipped, target, fmt)
    return {
        "status": "completed",
        "operation": "reverse_complement_alignment",
        "tool_version": _version(),
        "sequences": len(flipped.sequences),
        "columns": len(flipped.alignment[0]) if len(flipped.alignment) else 0,
        "rows": {
            _text(name): row
            for name, row in zip(flipped.names, flipped.textize().alignment, strict=True)
        },
        "output": str(target),
    }


def guess_file_alphabet(
    path: str,
    *,
    kind: str = "sequence",
    verbose: bool = True,
) -> dict[str, Any]:
    """Detect whether a sequence file or an alignment file holds DNA, RNA or amino acids (SequenceFile.guess_alphabet, MSAFile.guess_alphabet).

    ``kind`` is ``sequence`` (FASTA, GenBank, EMBL, ...) or ``alignment`` (Stockholm, A2M, aligned FASTA, ...). The alphabet is
    inferred from the first entry; ``alphabet`` is None when it cannot be told reliably.
    """

    selected = _choice(kind, _SEQUENCE_KINDS, "kind")
    source = _existing(path, "path")
    pyhmmer = require("pyhmmer", _EXTRA)
    opener = pyhmmer.easel.SequenceFile if selected == "sequence" else pyhmmer.easel.MSAFile
    try:
        with opener(str(source)) as handle:
            guessed = handle.guess_alphabet()
    except (OSError, ValueError, EOFError) as exc:
        raise ConfigurationError(f"{source.name} could not be read: {exc}") from exc
    name = None
    if guessed is not None:
        name = "DNA" if guessed.is_dna() else "RNA" if guessed.is_rna() else "amino"
    return {
        "status": "completed",
        "operation": "guess_file_alphabet",
        "tool_version": _version(),
        "kind": selected,
        "alphabet": name,
    }


def compute_alignment_weights(
    alignment: str,
    *,
    method: str = "pb",
    max_identity: float = 0.62,
    input_format: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute sequence weights of a DNA alignment with Easel (MSA.compute_weights, esl-weight).

    ``method`` is ``pb`` (position-based, the default), ``gsc`` (Gerstein-Sonnhammer-Chothia) or ``blosum`` (BLOSUM-style
    clustering at ``max_identity``). Returns one weight per sequence, in alignment order.
    """

    way = _choice(method, _WEIGHT_METHODS, "method")
    identity = _number(max_identity, "max_identity", 0.0, 1.0)
    source = _existing(alignment, "alignment")
    pyhmmer = require("pyhmmer", _EXTRA)
    msa = _read_msa(pyhmmer, source, pyhmmer.easel.Alphabet.dna(), input_format)
    try:
        weights = msa.compute_weights(way, identity)
    except ValueError as exc:
        raise ConfigurationError(f"pyhmmer could not compute weights: {exc}") from exc
    values = [float(item) for item in weights]
    return {
        "status": "completed",
        "operation": "compute_alignment_weights",
        "tool_version": _version(),
        "method": way,
        "names": [_text(name) for name in msa.names],
        "weights": values,
        "weight_sum": float(sum(values)),
    }


def filter_alignment_by_identity(
    alignment: str,
    output_alignment: str,
    *,
    max_identity: float = 0.8,
    preference: str = "conscover",
    seed: int = 42,
    input_format: str | None = None,
    output_format: str = "afa",
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove alignment sequences whose pairwise identity reaches ``max_identity`` and write the rest to a new file (DigitalMSA.identity_filter, esl-weight -f).

    Of two similar sequences ``preference`` keeps the one covering more consensus columns (``conscover``), the first one
    (``origorder``) or a random one (``random``, uses ``seed``). The default ``max_identity`` is 0.8.
    """

    identity = _number(max_identity, "max_identity", 0.0, 1.0)
    rule = _choice(preference, _PREFERENCES, "preference")
    seed_value = _integer(seed, "seed", 0, 2**31 - 1)
    fmt = _choice(output_format, _MSA_FORMATS, "output_format")
    source = _existing(alignment, "alignment")
    target = _new(output_alignment, "output_alignment")
    pyhmmer = require("pyhmmer", _EXTRA)
    msa = _read_msa(pyhmmer, source, pyhmmer.easel.Alphabet.dna(), input_format)
    try:
        kept = msa.identity_filter(identity, preference=rule, seed=seed_value)
    except ValueError as exc:
        raise ConfigurationError(f"pyhmmer could not filter the alignment: {exc}") from exc
    _write_msa(pyhmmer, kept, target, fmt)
    kept_names = [_text(name) for name in kept.names]
    return {
        "status": "completed",
        "operation": "filter_alignment_by_identity",
        "tool_version": _version(),
        "input_sequences": len(msa.names),
        "kept_sequences": len(kept_names),
        "kept_names": kept_names,
        "output": str(target),
    }


def build_sequence_ssi_index(
    sequence_file: str,
    *,
    aliases: Mapping[str, str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Index the records of a FASTA file by name in a new SSI file next to it (SSIWriter.add_file, add_key, add_alias).

    The index is written to ``<sequence_file>.ssi`` (which must not exist yet), where pyhmmer finds it automatically.
    Each record name (the first word of the header) is a key; ``aliases`` maps extra names to existing record names.
    Use ``fetch_sequence_by_ssi`` to retrieve records by name.
    """

    source = _existing(sequence_file, "sequence_file")
    target = _new(str(source) + ".ssi", "the SSI index")
    pyhmmer = require("pyhmmer", _EXTRA)
    data = source.read_bytes()
    records: list[list[Any]] = []
    for match in _HEADER.finditer(data):
        records.append([match.group(1).decode("ascii", "replace"), match.start(), match.end()])
    if not records:
        raise ConfigurationError(f"{source.name} holds no FASTA records.")
    for position, record in enumerate(records):
        record.append(records[position + 1][1] if position + 1 < len(records) else len(data))
    names = {record[0] for record in records}
    if len(names) != len(records):
        raise ConfigurationError("record names must be unique to build an SSI index.")
    extra = dict(aliases or {})
    for alias, key in extra.items():
        if key not in names or alias in names:
            raise ConfigurationError(
                f"alias {alias!r} must point to an existing record and not repeat a record name."
            )
    try:
        with pyhmmer.easel.SSIWriter(str(target)) as writer:
            handle = writer.add_file(str(source), 1)
            for name, start, data_start, end in records:
                writer.add_key(name, handle, start, data_start, end - start)
            for alias, key in extra.items():
                writer.add_alias(alias, key)
    except (OSError, ValueError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"pyhmmer could not build the SSI index: {exc}") from exc
    return {
        "status": "completed",
        "operation": "build_sequence_ssi_index",
        "tool_version": _version(),
        "sequence_file": str(source),
        "keys": len(records),
        "aliases": len(extra),
        "index": str(target),
    }


def fetch_sequence_by_ssi(
    sequence_file: str,
    name: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch one record of a FASTA file by name or alias through its SSI index (SSIReader.find_name).

    The index ``<sequence_file>.ssi`` is made by ``build_sequence_ssi_index``. The result holds the record name,
    description, length and sequence.
    """

    if not isinstance(name, str) or not name:
        raise ConfigurationError("name must be a non-empty string.")
    source = _existing(sequence_file, "sequence_file")
    _existing(str(source) + ".ssi", "the SSI index")
    pyhmmer = require("pyhmmer", _EXTRA)
    try:
        with pyhmmer.easel.SequenceFile(str(source)) as handle:
            if handle.indexed is None:
                raise ConfigurationError("the SSI index could not be opened.")
            record = handle.indexed[name]
    except KeyError as exc:
        raise ConfigurationError(f"{name!r} is not in the SSI index.") from exc
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"pyhmmer could not read the record: {exc}") from exc
    return {
        "status": "completed",
        "operation": "fetch_sequence_by_ssi",
        "tool_version": _version(),
        "name": _text(record.name),
        "description": _text(record.description),
        "length": len(record.sequence),
        "sequence": record.sequence,
    }


def pyhmmer_read_sequences(
    path: str, *, max_records: int = 1000, verbose: bool = True
) -> dict[str, Any]:
    """Read FASTA, GenBank or EMBL sequences using Easel SequenceFile, preserving upstream text."""
    source = _existing(path, "path")
    limit = _integer(max_records, "max_records", 1, _MAX_SEQUENCES)
    module = require("pyhmmer", _EXTRA)
    records = []
    with module.easel.SequenceFile(str(source)) as handle:
        for sequence in handle:
            if len(records) == limit:
                raise ConfigurationError("Sequence file exceeds max_records.")
            records.append(
                {
                    "name": _text(sequence.name),
                    "description": _text(sequence.description),
                    "sequence": sequence.sequence,
                }
            )
    return {"status": "completed", "backend": "pyhmmer", "version": _version(), "records": records}


def pyhmmer_write_sequences(
    records: Sequence[Mapping[str, str]], output_fasta: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Write sequence records to a new FASTA using Easel TextSequence.write."""
    target = _new(output_fasta, "output_fasta")
    module = require("pyhmmer", _EXTRA)
    values = [
        module.easel.TextSequence(
            name=r["name"].encode(),
            sequence=r["sequence"],
            description=r.get("description", "").encode(),
        )
        for r in records
    ]
    with target.open("wb") as handle:
        for sequence in values:
            sequence.write(handle)
    return {
        "status": "completed",
        "backend": "pyhmmer",
        "version": _version(),
        "record_count": len(values),
        "outputs": [str(target)],
    }


def pyhmmer_read_alignment(
    path: str, *, input_format: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Read a DNA multiple alignment with Easel MSAFile and return names and gapped sequences."""
    source = _existing(path, "path")
    module = require("pyhmmer", _EXTRA)
    msa = _read_msa(module, source, module.easel.Alphabet.dna(), input_format).textize()
    return {
        "status": "completed",
        "backend": "pyhmmer",
        "version": _version(),
        "records": [
            {"name": _text(n), "sequence": str(s)}
            for n, s in zip(msa.names, msa.alignment, strict=True)
        ],
    }


def pyhmmer_write_alignment(
    records: Mapping[str, str],
    output_alignment: str,
    *,
    output_format: str = "afa",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a multiple alignment in a supported Easel format (afa, Stockholm, A2M, Clustal, PHYLIP...)."""
    target = _new(output_alignment, "output_alignment")
    fmt = _choice(output_format, _MSA_FORMATS, "output_format")
    module = require("pyhmmer", _EXTRA)
    msa = module.easel.TextMSA(
        sequences=[
            module.easel.TextSequence(name=n.encode(), sequence=s) for n, s in records.items()
        ]
    )
    _write_msa(module, msa.digitize(module.easel.Alphabet.dna()), target, fmt)
    return {
        "status": "completed",
        "backend": "pyhmmer",
        "version": _version(),
        "outputs": [str(target)],
        "record_count": len(records),
    }


__all__ = [
    "build_hmm_from_sequence",
    "build_sequence_ssi_index",
    "compute_alignment_weights",
    "emit_hmm_alignment",
    "emit_hmm_sequences",
    "fetch_sequence_by_ssi",
    "filter_alignment_by_identity",
    "guess_file_alphabet",
    "pyhmmer_hmmscan",
    "pyhmmer_hmmsearch",
    "pyhmmer_read_alignment",
    "pyhmmer_read_sequences",
    "pyhmmer_write_alignment",
    "pyhmmer_write_sequences",
    "reverse_complement_alignment",
]

_register_public_verbose(globals(), tuple(__all__))
