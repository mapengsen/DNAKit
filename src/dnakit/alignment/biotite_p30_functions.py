"""Biotite sequence and alignment functions (integration pack P30).

Plain-Python wrappers over ``biotite.sequence`` and ``biotite.sequence.align``. Sequences and
alignments are passed as plain strings (an alignment is a list of equal-length gapped strings,
gap character ``-``); results are JSON-serialisable dictionaries. ``seq_type`` selects the
alphabet (``dna`` or ``protein``). Biotite is imported only when a function runs. The tested
version is biotite 2.0.0.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-biotite"
_PACKAGE = "biotite"
_TESTED_VERSION = "2.0.0"
_SEQ_TYPES = ("dna", "protein")
_MAX_LENGTH = 200_000
_MAX_DP_LENGTH = 20_000
_MAX_SEQUENCES = 2_000
_MAX_ALIGNMENT_COLUMNS = 200_000
_MAX_KMERS_RETURNED = 100_000
_MAX_K = 12
_MODES = ("all", "not_terminal", "shortest")


def _provenance(**notes: Any) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {
        "package": _PACKAGE,
        "package_version": installed,
        "tested_version": _TESTED_VERSION,
        **notes,
    }


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _int(
    value: object, name: str, *, minimum: int | None = None, maximum: int | None = None
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum}.")
    return value


def _text(value: object, name: str, *, limit: int = _MAX_LENGTH) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    cleaned = value.strip()
    if len(cleaned) > limit:
        raise ConfigurationError(f"{name} may be at most {limit} characters long.")
    return cleaned


def _modules() -> tuple[Any, Any]:
    sequence = require("biotite.sequence", _EXTRA)
    align = require("biotite.sequence.align", _EXTRA)
    return sequence, align


def _factory(seq_type: str) -> Any:
    sequence = require("biotite.sequence", _EXTRA)
    return sequence.NucleotideSequence if seq_type == "dna" else sequence.ProteinSequence


def _make(text: object, seq_type: str, name: str, *, limit: int = _MAX_LENGTH) -> Any:
    cleaned = _text(text, name, limit=limit)
    try:
        return _factory(seq_type)(cleaned)
    except Exception as exc:
        raise ConfigurationError(f"{name} is not a valid {seq_type} sequence: {exc}") from exc


def _strings(values: object, name: str, *, minimum: int = 1) -> list[str]:
    if (
        isinstance(values, (str, bytes))
        or not isinstance(values, Sequence)
        or not minimum <= len(values) <= _MAX_SEQUENCES
    ):
        raise ConfigurationError(f"{name} must be a list of {minimum} to {_MAX_SEQUENCES} strings.")
    return [_text(value, f"{name}[{i}]") for i, value in enumerate(values)]


def _gapped(values: object, name: str, minimum: int = 2) -> list[str]:
    items = _strings(values, name, minimum=minimum)
    if len({len(item) for item in items}) != 1:
        raise ConfigurationError(f"{name} must all have the same length (gapped strings).")
    if len(items[0]) > _MAX_ALIGNMENT_COLUMNS:
        raise ConfigurationError(f"{name} may have at most {_MAX_ALIGNMENT_COLUMNS} columns.")
    return items


def _alignment(values: object, seq_type: str, name: str = "aligned", minimum: int = 2) -> Any:
    items = _gapped(values, name, minimum)
    _, align = _modules()
    try:
        return align.Alignment.from_strings(items, _factory(seq_type))
    except Exception as exc:
        raise ConfigurationError(f"{name} is not a valid {seq_type} alignment: {exc}") from exc


def _alignment_dict(alignment: Any, score: int | None = None) -> dict[str, Any]:
    return {
        "aligned": [str(item) for item in alignment.get_gapped_sequences()],
        "trace": [[int(a), int(b)] for a, b in alignment.trace.tolist()]
        if alignment.trace.shape[1] == 2
        else alignment.trace.tolist(),
        "score": score
        if score is not None
        else (None if alignment.score is None else int(alignment.score)),
        "length": int(alignment.trace.shape[0]),
    }


def _matrix(seq_type: str, matrix: str | None, match: int, mismatch: int) -> Any:
    """A (match, mismatch) pair, or a named Biotite substitution matrix such as NUC or BLOSUM62."""
    _int(match, "match")
    _int(mismatch, "mismatch")
    if matrix is None:
        return (match, mismatch)
    name = _text(matrix, "matrix", limit=40)
    sequence, align = _modules()
    alphabet = (
        sequence.NucleotideSequence.alphabet_ambiguous
        if seq_type == "dna"
        else sequence.ProteinSequence.alphabet
    )
    try:
        return align.SubstitutionMatrix(alphabet, alphabet, name)
    except Exception as exc:
        raise ConfigurationError(f"Unknown substitution matrix {name!r}: {exc}") from exc


def _gap(gap_penalty: int, gap_extension: int | None) -> Any:
    """Linear gap penalty, or (open, extend) when ``gap_extension`` is given."""
    _int(gap_penalty, "gap_penalty")
    if gap_extension is None:
        return gap_penalty
    return (gap_penalty, _int(gap_extension, "gap_extension"))


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _kmer_alphabet(seq_type: str, k: int, spacing: str | None = None) -> Any:
    _, align = _modules()
    base = (
        _factory(seq_type).alphabet
        if seq_type == "protein"
        else _factory(seq_type).alphabet_unambiguous
    )
    return align.KmerAlphabet(base, _int(k, "k", minimum=1, maximum=_MAX_K), spacing=spacing)


def _decode_kmer(alphabet: Any, code: int) -> str:
    return "".join(str(item) for item in alphabet.decode(int(code)))


def _sequence_for_kmers(text: object, seq_type: str, name: str) -> Any:
    sequence = _make(text, seq_type, name)
    if seq_type == "dna" and sequence.alphabet is not _factory("dna").alphabet_unambiguous:
        raise ConfigurationError(f"{name} may only contain A, C, G and T for k-mer functions.")
    return sequence


def codon_table_is_start_codon(
    codons: Sequence[str],
    *,
    table: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tell whether each codon is a start codon of an NCBI genetic code table.

    Wraps ``CodonTable.is_start_codon``. ``codons`` are DNA triplets (A, C, G, T); ``table`` is
    the NCBI translation table number (1 standard, 11 bacterial, ...).
    """

    items = _strings(codons, "codons")
    table_number = _int(table, "table", minimum=1)
    sequence, _ = _modules()
    numpy = require("numpy", _EXTRA)
    try:
        code_table = sequence.CodonTable.load(table_number)
    except Exception as exc:
        raise ConfigurationError(f"Unknown codon table {table_number}: {exc}") from exc
    alphabet = sequence.NucleotideSequence.alphabet_unambiguous
    rows = []
    for codon in items:
        if len(codon) != 3 or set(codon.upper()) - set("ACGT"):
            raise ConfigurationError(f"{codon!r} is not a DNA codon of three A/C/G/T bases.")
        rows.append([alphabet.encode(base) for base in codon.upper()])
    flags = code_table.is_start_codon(numpy.array(rows, dtype=int))
    return {
        "status": "completed",
        "operation": "codon_table_is_start_codon",
        "table": table_number,
        "codons": [codon.upper() for codon in items],
        "is_start": [bool(flag) for flag in flags],
        "provenance": _provenance(method="CodonTable.is_start_codon"),
    }


def _profile(aligned: Sequence[str], seq_type: str) -> Any:
    sequence, _ = _modules()
    alignment = _alignment(aligned, seq_type, minimum=1)
    return sequence.SequenceProfile.from_alignment(alignment)


def sequence_profile_from_alignment(
    aligned: Sequence[str],
    *,
    seq_type: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the symbols of every alignment column (a position-frequency matrix).

    ``aligned`` are equal-length gapped strings. Returns the symbols (including the gap column)
    and one count row per alignment column (``SequenceProfile.from_alignment``).
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    profile = _profile(aligned, kind)
    symbols = [str(item) for item in profile.alphabet.get_symbols()]
    return {
        "status": "completed",
        "operation": "sequence_profile_from_alignment",
        "seq_type": kind,
        "symbols": symbols,
        "counts": profile.symbols.tolist(),
        "gap_counts": profile.gaps.tolist(),
        "columns": int(profile.symbols.shape[0]),
        "provenance": _provenance(method="SequenceProfile.from_alignment"),
    }


def sequence_profile_probability(
    aligned: Sequence[str],
    sequence: str,
    *,
    seq_type: str = "dna",
    pseudocount: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Probability of ``sequence`` under the profile of an alignment.

    The profile is built from ``aligned`` (``SequenceProfile.from_alignment``); the sequence must
    have one symbol per alignment column. ``pseudocount`` is added to every symbol count.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    count = _int(pseudocount, "pseudocount", minimum=0)
    profile = _profile(aligned, kind)
    query = _make(sequence, kind, "sequence")
    if len(query) != profile.symbols.shape[0]:
        raise ConfigurationError(
            f"sequence has {len(query)} symbols but the alignment has "
            f"{profile.symbols.shape[0]} columns."
        )
    value = profile.sequence_probability(query, pseudocount=count)
    return {
        "status": "completed",
        "operation": "sequence_profile_probability",
        "probability": float(value),
        "pseudocount": count,
        "provenance": _provenance(method="SequenceProfile.sequence_probability"),
    }


def sequence_profile_score(
    aligned: Sequence[str],
    sequence: str,
    *,
    seq_type: str = "dna",
    background_frequencies: float | Sequence[float] | None = None,
    pseudocount: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Log-odds score of ``sequence`` against the weight matrix of an alignment profile.

    ``background_frequencies`` is one float or one frequency per alphabet symbol; None means
    uniform. ``pseudocount`` is added to every symbol count (``SequenceProfile.sequence_score``).
    """

    numpy = require("numpy", _EXTRA)
    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    count = _int(pseudocount, "pseudocount", minimum=0)
    profile = _profile(aligned, kind)
    query = _make(sequence, kind, "sequence")
    if len(query) != profile.symbols.shape[0]:
        raise ConfigurationError(
            f"sequence has {len(query)} symbols but the alignment has "
            f"{profile.symbols.shape[0]} columns."
        )
    background: Any = None
    if background_frequencies is not None:
        if isinstance(background_frequencies, (int, float)) and not isinstance(
            background_frequencies, bool
        ):
            background = float(background_frequencies)
        else:
            background = numpy.asarray([float(item) for item in background_frequencies])
            if background.shape != (profile.symbols.shape[1],):
                raise ConfigurationError(
                    "background_frequencies must list one value per alphabet symbol."
                )
    value = profile.sequence_score(query, background_frequencies=background, pseudocount=count)
    return {
        "status": "completed",
        "operation": "sequence_profile_score",
        "score": float(value),
        "pseudocount": count,
        "provenance": _provenance(method="SequenceProfile.sequence_score"),
    }


def alignment_from_strings(
    aligned: Sequence[str],
    *,
    seq_type: str = "dna",
    gap_character: str = "-",
    verbose: bool = True,
) -> dict[str, Any]:
    """Build an alignment (trace of residue indices) from gapped strings.

    Covers ``Alignment.from_strings`` and ``Alignment.trace_from_strings``. In the trace, ``-1``
    marks a gap; the ungapped sequences are returned too.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    gap = _text(gap_character, "gap_character", limit=1)
    items = _gapped(aligned, "aligned", minimum=1)
    _, align = _modules()
    if gap != "-":
        items = [item.replace(gap, "-") for item in items]
    alignment = align.Alignment.from_strings(items, _factory(kind))
    result = _alignment_dict(alignment)
    result.update(
        {
            "status": "completed",
            "operation": "alignment_from_strings",
            "sequences": [str(item) for item in alignment.sequences],
            "provenance": _provenance(method="Alignment.from_strings"),
        }
    )
    return result


def kmer_index(
    sequences: Sequence[str],
    *,
    k: int = 4,
    seq_type: str = "dna",
    spacing: str | None = None,
    bucketed: bool = False,
    n_buckets: int | None = None,
    max_kmers: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Index the positions of every k-mer of a set of sequences.

    ``bucketed=False`` uses ``KmerTable.from_sequences`` (exact k-mer to position lists);
    ``bucketed=True`` uses ``BucketKmerTable.from_sequences`` (hash buckets, ``n_buckets``).
    ``spacing`` is a spaced-seed pattern such as ``"1101"``. Each position is
    ``[sequence_index, offset]``. At most ``max_kmers`` distinct k-mers are listed.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    flag = _flag(bucketed, "bucketed")
    limit = _int(max_kmers, "max_kmers", minimum=1, maximum=_MAX_KMERS_RETURNED)
    items = [
        _sequence_for_kmers(text, kind, f"sequences[{i}]")
        for i, text in enumerate(_strings(sequences, "sequences"))
    ]
    _, align = _modules()
    alphabet = _kmer_alphabet(kind, k, spacing)
    base_alphabet = items[0].alphabet
    if flag:
        buckets = None if n_buckets is None else _int(n_buckets, "n_buckets", minimum=1)
        table = align.BucketKmerTable.from_sequences(
            k, items, alphabet=base_alphabet, spacing=spacing, n_buckets=buckets
        )
    else:
        table = align.KmerTable.from_sequences(k, items, alphabet=base_alphabet, spacing=spacing)
    entries = []
    for code in table.get_kmers() if hasattr(table, "get_kmers") else []:
        if len(entries) >= limit:
            break
        positions = table[int(code)]
        entries.append(
            {
                "kmer": _decode_kmer(alphabet, int(code)),
                "positions": [[int(a), int(b)] for a, b in positions.tolist()],
            }
        )
    return {
        "status": "completed",
        "operation": "kmer_index",
        "k": k,
        "table": "bucket" if flag else "exact",
        "sequence_count": len(items),
        "indexed_kmers": len(entries),
        "kmers": entries,
        "truncated": len(entries) >= limit,
        "provenance": _provenance(
            method="BucketKmerTable.from_sequences" if flag else "KmerTable.from_sequences"
        ),
    }


def kmer_match(
    reference_sequences: Sequence[str],
    query: str,
    *,
    k: int = 4,
    seq_type: str = "dna",
    spacing: str | None = None,
    similarity_threshold: int | None = None,
    matrix: str | None = None,
    bucketed: bool = False,
    max_matches: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find k-mer matches (alignment seeds) between a query and reference sequences.

    The references are indexed once (``KmerTable`` or, with ``bucketed=True``,
    ``BucketKmerTable``) and the query is matched against them. With ``similarity_threshold``
    similar k-mers whose substitution score reaches the threshold also match
    (``ScoreThresholdRule``; ``matrix`` names the substitution matrix, default NUC for DNA and
    BLOSUM62 for protein). Each match is ``[query_position, reference_index,
    reference_position]``.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    flag = _flag(bucketed, "bucketed")
    limit = _int(max_matches, "max_matches", minimum=1, maximum=1_000_000)
    refs = [
        _sequence_for_kmers(text, kind, f"reference_sequences[{i}]")
        for i, text in enumerate(_strings(reference_sequences, "reference_sequences"))
    ]
    query_sequence = _sequence_for_kmers(query, kind, "query")
    _, align = _modules()
    _kmer_alphabet(kind, k, spacing)
    cls = align.BucketKmerTable if flag else align.KmerTable
    table = cls.from_sequences(k, refs, alphabet=refs[0].alphabet, spacing=spacing)
    rule = None
    if similarity_threshold is not None:
        threshold = _int(similarity_threshold, "similarity_threshold")
        matrix_name = matrix or ("NUC" if kind == "dna" else "BLOSUM62")
        substitution = _matrix(kind, matrix_name, 1, -1)
        if kind == "dna":
            sequence, _ = _modules()
            alphabet = sequence.NucleotideSequence.alphabet_unambiguous
            substitution = align.SubstitutionMatrix(alphabet, alphabet, matrix_name)
        rule = align.ScoreThresholdRule(substitution, threshold)
    matches = table.match(query_sequence, similarity_rule=rule)
    rows = [[int(a), int(b), int(c)] for a, b, c in matches.tolist()]
    return {
        "status": "completed",
        "operation": "kmer_match",
        "k": k,
        "match_count": len(rows),
        "matches": rows[:limit],
        "truncated": len(rows) > limit,
        "provenance": _provenance(method=f"{cls.__name__}.match"),
    }


def syncmer_selection(
    sequence: str,
    *,
    k: int = 5,
    s: int = 3,
    seq_type: str = "dna",
    offset: Sequence[int] = (0,),
    verbose: bool = True,
) -> dict[str, Any]:
    """Select the syncmers (k-mers whose smallest s-mer is at one of the offsets) of a sequence.

    Uses ``CachedSyncmerSelector``. ``s`` is the s-mer length (s < k); ``offset`` the allowed
    positions of the minimal s-mer inside the k-mer. Returns positions and k-mer strings.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    k_value = _int(k, "k", minimum=2, maximum=_MAX_K)
    s_value = _int(s, "s", minimum=1, maximum=k_value - 1)
    offsets = tuple(_int(item, "offset", minimum=0) for item in offset)
    query = _sequence_for_kmers(sequence, kind, "sequence")
    _, align = _modules()
    selector = align.CachedSyncmerSelector(query.alphabet, k_value, s_value, offset=offsets)
    positions, codes = selector.select(query)
    alphabet = _kmer_alphabet(kind, k_value)
    return {
        "status": "completed",
        "operation": "syncmer_selection",
        "k": k_value,
        "s": s_value,
        "positions": [int(item) for item in positions],
        "kmers": [_decode_kmer(alphabet, int(code)) for code in codes],
        "provenance": _provenance(method="CachedSyncmerSelector.select"),
    }


def mincode_selection(
    sequence: str,
    *,
    k: int = 5,
    compression: float = 2.0,
    seq_type: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Select a subset of the k-mers of a sequence with the mincode method.

    Keeps the k-mers whose code is below a threshold set by ``compression`` (the desired factor
    by which the k-mer count is reduced; ``MincodeSelector``).
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    if isinstance(compression, bool) or not isinstance(compression, (int, float)):
        raise ConfigurationError("compression must be a number.")
    if compression < 1:
        raise ConfigurationError("compression must be at least 1.")
    query = _sequence_for_kmers(sequence, kind, "sequence")
    _, align = _modules()
    alphabet = align.KmerAlphabet(query.alphabet, _int(k, "k", minimum=1, maximum=_MAX_K))
    selector = align.MincodeSelector(alphabet, float(compression))
    positions, codes = selector.select(query)
    return {
        "status": "completed",
        "operation": "mincode_selection",
        "k": k,
        "compression": float(compression),
        "positions": [int(item) for item in positions],
        "kmers": [_decode_kmer(alphabet, int(code)) for code in codes],
        "provenance": _provenance(method="MincodeSelector.select"),
    }


def similar_kmers(
    kmer: str,
    *,
    threshold: int,
    seq_type: str = "dna",
    matrix: str | None = None,
    max_results: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List all k-mers whose substitution score against ``kmer`` reaches ``threshold``.

    Uses ``ScoreThresholdRule.similar_kmers``. ``matrix`` is a Biotite matrix name (default NUC
    for DNA, BLOSUM62 for protein).
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    limit = _int(max_results, "max_results", minimum=1, maximum=1_000_000)
    threshold_value = _int(threshold, "threshold")
    word = _text(kmer, "kmer", limit=_MAX_K)
    sequence, align = _modules()
    alphabet = _kmer_alphabet(kind, len(word))
    base = alphabet.base_alphabet
    name = matrix or ("NUC" if kind == "dna" else "BLOSUM62")
    try:
        substitution = align.SubstitutionMatrix(base, base, name)
        code = alphabet.encode(word)
    except Exception as exc:
        raise ConfigurationError(f"Invalid kmer or matrix: {exc}") from exc
    rule = align.ScoreThresholdRule(substitution, threshold_value)
    codes = rule.similar_kmers(alphabet, int(code))
    words = [_decode_kmer(alphabet, int(item)) for item in codes]
    del sequence
    return {
        "status": "completed",
        "operation": "similar_kmers",
        "kmer": word,
        "threshold": threshold_value,
        "matrix": name,
        "count": len(words),
        "similar_kmers": words[:limit],
        "truncated": len(words) > limit,
        "provenance": _provenance(method="ScoreThresholdRule.similar_kmers"),
    }


def _evalue_estimator(
    seq_type: str,
    matrix: str | None,
    match: int,
    mismatch: int,
    gap_penalty: int,
    gap_extension: int | None,
    frequencies: Sequence[float] | None,
    sample_length: int,
    sample_size: int,
    seed: int,
) -> Any:
    numpy = require("numpy", _EXTRA)
    sequence, align = _modules()
    alphabet = (
        sequence.NucleotideSequence.alphabet_unambiguous
        if seq_type == "dna"
        else sequence.ProteinSequence.alphabet
    )
    name = matrix or ("NUC" if seq_type == "dna" else "BLOSUM62")
    if matrix is None and seq_type == "dna" and (match, mismatch) != (1, -1):
        matrix_object = align.SubstitutionMatrix(
            alphabet,
            alphabet,
            {a: {b: match if a == b else mismatch for b in alphabet} for a in alphabet},
        )
    else:
        matrix_object = align.SubstitutionMatrix(alphabet, alphabet, name)
    size = len(alphabet)
    freq = (
        numpy.full(size, 1.0 / size)
        if frequencies is None
        else numpy.asarray([float(item) for item in frequencies])
    )
    if freq.shape != (size,):
        raise ConfigurationError(f"frequencies must list {size} values.")
    return align.EValueEstimator.from_samples(
        alphabet,
        matrix_object,
        _gap(gap_penalty, gap_extension),
        freq,
        sample_length=_int(sample_length, "sample_length", minimum=10, maximum=100_000),
        sample_size=_int(sample_size, "sample_size", minimum=10, maximum=1_000_000),
        rng=numpy.random.default_rng(_int(seed, "seed", minimum=0)),
    )


def alignment_evalue_parameters(
    *,
    seq_type: str = "dna",
    matrix: str | None = None,
    match: int = 1,
    mismatch: int = -1,
    gap_penalty: int = -5,
    gap_extension: int | None = None,
    frequencies: Sequence[float] | None = None,
    sample_length: int = 1000,
    sample_size: int = 1000,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the Gumbel parameters lambda and K of a scoring scheme by sampling.

    Random sequences are drawn with the given symbol ``frequencies`` (uniform by default) and
    aligned locally (``EValueEstimator.from_samples``). ``matrix`` names a Biotite matrix;
    otherwise the ``match``/``mismatch`` scores are used for DNA. ``seed`` makes it reproducible.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    estimator = _evalue_estimator(
        kind,
        matrix,
        match,
        mismatch,
        gap_penalty,
        gap_extension,
        frequencies,
        sample_length,
        sample_size,
        seed,
    )
    return {
        "status": "completed",
        "operation": "alignment_evalue_parameters",
        "lambda": float(estimator.lam),
        "k": float(estimator.k),
        "sample_length": sample_length,
        "sample_size": sample_size,
        "seed": seed,
        "provenance": _provenance(method="EValueEstimator.from_samples"),
    }


def alignment_log_evalue(
    score: int,
    seq1_length: int,
    seq2_length: int,
    *,
    lam: float | None = None,
    k: float | None = None,
    seq_type: str = "dna",
    matrix: str | None = None,
    match: int = 1,
    mismatch: int = -1,
    gap_penalty: int = -5,
    gap_extension: int | None = None,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Natural logarithm of the E-value of an alignment score (``EValueEstimator.log_evalue``).

    Give ``lam`` and ``k`` (from ``alignment_evalue_parameters``) or leave both None to estimate
    them by sampling for the given scoring scheme. Also returns the E-value itself.
    """

    import math

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    _int(score, "score")
    length1 = _int(seq1_length, "seq1_length", minimum=1)
    length2 = _int(seq2_length, "seq2_length", minimum=1)
    _, align = _modules()
    if (lam is None) != (k is None):
        raise ConfigurationError("Give both lam and k, or neither.")
    if lam is not None and k is not None:
        if lam <= 0 or k <= 0:
            raise ConfigurationError("lam and k must be positive.")
        estimator = align.EValueEstimator(float(lam), float(k))
    else:
        estimator = _evalue_estimator(
            kind, matrix, match, mismatch, gap_penalty, gap_extension, None, 1000, 1000, seed
        )
    value = float(estimator.log_evalue(score, length1, length2))
    return {
        "status": "completed",
        "operation": "alignment_log_evalue",
        "log_evalue": value,
        "evalue": math.exp(value) if value < 700 else None,
        "lambda": float(estimator.lam),
        "k": float(estimator.k),
        "provenance": _provenance(method="EValueEstimator.log_evalue"),
    }


def seed_extension_ungapped(
    sequence_1: str,
    sequence_2: str,
    *,
    seed: Sequence[int],
    threshold: int,
    seq_type: str = "dna",
    matrix: str | None = None,
    match: int = 1,
    mismatch: int = -1,
    direction: str = "both",
    verbose: bool = True,
) -> dict[str, Any]:
    """Extend a seed without gaps in both sequences with the X-drop rule (``SeedExtension``).

    ``seed`` is ``[position_in_sequence_1, position_in_sequence_2]``; extension stops when the
    score falls ``threshold`` below the best score seen. ``direction`` is ``both``,
    ``upstream`` or ``downstream``.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    way = _choice(direction, ("both", "upstream", "downstream"), "direction")
    first = _make(sequence_1, kind, "sequence_1")
    second = _make(sequence_2, kind, "sequence_2")
    if len(seed) != 2:
        raise ConfigurationError("seed must be [position_1, position_2].")
    seed_pair = (_int(seed[0], "seed[0]", minimum=0), _int(seed[1], "seed[1]", minimum=0))
    if seed_pair[0] >= len(first) or seed_pair[1] >= len(second):
        raise ConfigurationError("seed lies outside the sequences.")
    _, align = _modules()
    extension = align.SeedExtension(
        _matrix(kind, matrix, match, mismatch), _int(threshold, "threshold", minimum=0), way
    )
    alignment = extension.align(first, second, seed_pair)
    result = _alignment_dict(alignment)
    result.update(
        {
            "status": "completed",
            "operation": "seed_extension_ungapped",
            "seed": list(seed_pair),
            "provenance": _provenance(method="SeedExtension.align"),
        }
    )
    return result


def align_banded_pairwise(
    sequence_1: str,
    sequence_2: str,
    *,
    band: Sequence[int],
    seq_type: str = "dna",
    matrix: str | None = None,
    match: int = 1,
    mismatch: int = -1,
    gap_penalty: int = -1,
    gap_extension: int | None = None,
    local: bool = False,
    max_number: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise alignment restricted to a diagonal band (``align_banded``).

    ``band`` is ``[lower, upper]`` diagonals (``j - i``). ``local=True`` gives a local
    alignment. The default scores are biotite's: match 1, mismatch -1, linear gap -1.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    first = _make(sequence_1, kind, "sequence_1", limit=_MAX_DP_LENGTH)
    second = _make(sequence_2, kind, "sequence_2", limit=_MAX_DP_LENGTH)
    if len(band) != 2 or _int(band[0], "band[0]") > _int(band[1], "band[1]"):
        raise ConfigurationError("band must be [lower, upper] with lower <= upper.")
    _, align = _modules()
    alignments = align.align_banded(
        first,
        second,
        (band[0], band[1]),
        matrix=_matrix(kind, matrix, match, mismatch),
        gap_penalty=_gap(gap_penalty, gap_extension),
        local=_flag(local, "local"),
        max_number=_int(max_number, "max_number", minimum=1, maximum=10_000),
    )
    return {
        "status": "completed",
        "operation": "align_banded_pairwise",
        "alignment_count": len(alignments),
        "alignments": [_alignment_dict(item) for item in alignments],
        "provenance": _provenance(method="align_banded"),
    }


def align_local_gapped_from_seed(
    sequence_1: str,
    sequence_2: str,
    *,
    seed: Sequence[int],
    threshold: int,
    seq_type: str = "dna",
    matrix: str | None = None,
    match: int = 1,
    mismatch: int = -1,
    gap_penalty: int = -1,
    gap_extension: int | None = None,
    max_number: int = 1,
    direction: str = "both",
    verbose: bool = True,
) -> dict[str, Any]:
    """Gapped local alignment grown from a seed with the X-drop rule (``align_local_gapped``)."""

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    way = _choice(direction, ("both", "upstream", "downstream"), "direction")
    first = _make(sequence_1, kind, "sequence_1", limit=_MAX_DP_LENGTH)
    second = _make(sequence_2, kind, "sequence_2", limit=_MAX_DP_LENGTH)
    if len(seed) != 2:
        raise ConfigurationError("seed must be [position_1, position_2].")
    seed_pair = (_int(seed[0], "seed[0]", minimum=0), _int(seed[1], "seed[1]", minimum=0))
    if seed_pair[0] >= len(first) or seed_pair[1] >= len(second):
        raise ConfigurationError("seed lies outside the sequences.")
    _, align = _modules()
    alignments = align.align_local_gapped(
        first,
        second,
        seed_pair,
        _int(threshold, "threshold", minimum=0),
        matrix=_matrix(kind, matrix, match, mismatch),
        gap_penalty=_gap(gap_penalty, gap_extension),
        max_number=_int(max_number, "max_number", minimum=1, maximum=10_000),
        direction=way,
    )
    return {
        "status": "completed",
        "operation": "align_local_gapped_from_seed",
        "seed": list(seed_pair),
        "alignment_count": len(alignments),
        "alignments": [_alignment_dict(item) for item in alignments],
        "provenance": _provenance(method="align_local_gapped"),
    }


def align_multiple_progressive(
    sequences: Sequence[str],
    *,
    seq_type: str = "dna",
    matrix: str | None = None,
    match: int = 1,
    mismatch: int = -1,
    gap_penalty: int = -1,
    gap_extension: int | None = None,
    terminal_penalty: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Progressive multiple sequence alignment with a guide tree (``align_multiple``).

    Returns the gapped sequences (input order), the guide-tree leaf order for display, the
    guide-tree edges and the pairwise distance matrix used to build the tree.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    items = [
        _make(text, kind, f"sequences[{i}]", limit=_MAX_DP_LENGTH)
        for i, text in enumerate(_strings(sequences, "sequences", minimum=2))
    ]
    _, align = _modules()
    alignment, order, tree, distances = align.align_multiple(
        items,
        _matrix(kind, matrix, match, mismatch),
        gap_penalty=_gap(gap_penalty, gap_extension),
        terminal_penalty=_flag(terminal_penalty, "terminal_penalty"),
    )
    gapped = [str(item) for item in alignment.get_gapped_sequences()]
    return {
        "status": "completed",
        "operation": "align_multiple_progressive",
        "aligned": gapped,
        "display_order": [int(item) for item in order],
        "columns": int(alignment.trace.shape[0]),
        "guide_tree_edges": [[int(a), int(b)] for a, b in tree.edges()],
        "distances": [[float(v) for v in row] for row in distances.tolist()],
        "provenance": _provenance(method="align_multiple"),
    }


def align_ungapped_pairwise(
    sequence_1: str,
    sequence_2: str,
    *,
    seq_type: str = "dna",
    matrix: str | None = None,
    match: int = 1,
    mismatch: int = -1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align two sequences position by position without gaps (``align_ungapped``).

    The sequences must have the same length; the score is the sum of substitution scores.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    first = _make(sequence_1, kind, "sequence_1")
    second = _make(sequence_2, kind, "sequence_2")
    if len(first) != len(second):
        raise ConfigurationError("An ungapped alignment needs two sequences of equal length.")
    _, align = _modules()
    alignment = align.align_ungapped(first, second, _matrix(kind, matrix, match, mismatch))
    result = _alignment_dict(alignment)
    result.update(
        {
            "status": "completed",
            "operation": "align_ungapped_pairwise",
            "provenance": _provenance(method="align_ungapped"),
        }
    )
    return result


def alignment_remove_terminal_gaps(
    aligned: Sequence[str],
    *,
    seq_type: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim the leading and trailing columns that contain a gap (``remove_terminal_gaps``).

    Also reports the number of trimmed columns at each end (``find_terminal_gaps``).
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    alignment = _alignment(aligned, kind)
    _, align = _modules()
    start, stop = align.find_terminal_gaps(alignment)
    trimmed = align.remove_terminal_gaps(alignment)
    result = _alignment_dict(trimmed, None)
    result.update(
        {
            "status": "completed",
            "operation": "alignment_remove_terminal_gaps",
            "first_kept_column": int(start),
            "last_kept_column_exclusive": int(stop),
            "provenance": _provenance(method="remove_terminal_gaps"),
        }
    )
    return result


def alignment_remove_gaps(
    aligned: Sequence[str],
    *,
    seq_type: str = "dna",
    verbose: bool = True,
) -> dict[str, Any]:
    """Drop every alignment column that contains a gap in any sequence (``remove_gaps``)."""

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    _, align = _modules()
    trimmed = align.remove_gaps(_alignment(aligned, kind))
    result = _alignment_dict(trimmed, None)
    result.update(
        {
            "status": "completed",
            "operation": "alignment_remove_gaps",
            "provenance": _provenance(method="remove_gaps"),
        }
    )
    return result


def alignment_pairwise_identity(
    aligned: Sequence[str],
    *,
    seq_type: str = "dna",
    mode: str = "not_terminal",
    verbose: bool = True,
) -> dict[str, Any]:
    """Matrix of pairwise sequence identities of an alignment (``get_pairwise_sequence_identity``).

    ``mode``: ``all`` divides by all columns, ``not_terminal`` by the columns without terminal
    gaps of the pair (biotite default), ``shortest`` by the length of the shorter sequence.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    way = _choice(mode, _MODES, "mode")
    _, align = _modules()
    matrix = align.get_pairwise_sequence_identity(_alignment(aligned, kind), mode=way)
    return {
        "status": "completed",
        "operation": "alignment_pairwise_identity",
        "mode": way,
        "identity": [[float(v) for v in row] for row in matrix.tolist()],
        "provenance": _provenance(method="get_pairwise_sequence_identity"),
    }


def alignment_sequence_identity(
    aligned: Sequence[str],
    *,
    seq_type: str = "dna",
    mode: str = "not_terminal",
    verbose: bool = True,
) -> dict[str, Any]:
    """Overall sequence identity of an alignment, a single number (``get_sequence_identity``).

    ``mode`` is as in ``alignment_pairwise_identity``.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    way = _choice(mode, _MODES, "mode")
    _, align = _modules()
    value = align.get_sequence_identity(_alignment(aligned, kind), mode=way)
    return {
        "status": "completed",
        "operation": "alignment_sequence_identity",
        "mode": way,
        "identity": float(value),
        "provenance": _provenance(method="get_sequence_identity"),
    }


def alignment_score(
    aligned: Sequence[str],
    *,
    seq_type: str = "dna",
    matrix: str | None = None,
    match: int = 1,
    mismatch: int = -1,
    gap_penalty: int = -1,
    gap_extension: int | None = None,
    terminal_penalty: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score an existing pairwise alignment of gapped strings (``biotite...align.score``)."""

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    alignment = _alignment(aligned, kind, minimum=2)
    _, align = _modules()
    try:
        value = align.score(
            alignment,
            _matrix(kind, matrix, match, mismatch),
            _gap(gap_penalty, gap_extension),
            _flag(terminal_penalty, "terminal_penalty"),
        )
    except Exception as exc:
        raise ConfigurationError(f"The alignment could not be scored: {exc}") from exc
    return {
        "status": "completed",
        "operation": "alignment_score",
        "score": int(value),
        "terminal_penalty": terminal_penalty,
        "provenance": _provenance(method="align.score"),
    }


def plot_sequence_alignment(
    aligned: Sequence[str],
    output_path: str,
    *,
    seq_type: str = "dna",
    coloring: str = "type",
    values: Sequence[Sequence[float]] | None = None,
    labels: Sequence[str] | None = None,
    symbols_per_line: int = 50,
    show_numbers: bool = False,
    dpi: int = 150,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw an alignment to an image file and return the file path with the plotted data.

    ``coloring`` selects the symbol plotter: ``type`` colours by symbol type
    (``plot_alignment_type_based``), ``similarity`` by column similarity
    (``plot_alignment_similarity_based``), ``values`` shades each symbol by a per-position value
    (``plot_alignment_array``; ``values`` is a columns-by-sequences matrix of floats).
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    way = _choice(coloring, ("type", "similarity", "values"), "coloring")
    alignment = _alignment(aligned, kind)
    out = Path(output_path).expanduser().resolve()
    if out.exists() or not out.parent.is_dir():
        raise ConfigurationError(f"output_path must be a new file in an existing folder: {out}")
    if out.suffix.lower() not in (".png", ".svg", ".pdf"):
        raise ConfigurationError("output_path must end in .png, .svg or .pdf.")
    matplotlib = require("matplotlib", _EXTRA)
    matplotlib.use("Agg", force=True)
    pyplot = require("matplotlib.pyplot", _EXTRA)
    graphics = require("biotite.sequence.graphics", _EXTRA)
    numpy = require("numpy", _EXTRA)
    columns = int(alignment.trace.shape[0])
    row_labels = None if labels is None else [str(item) for item in labels]
    if row_labels is not None and len(row_labels) != len(alignment.sequences):
        raise ConfigurationError("labels must list one name per sequence.")
    per_line = _int(symbols_per_line, "symbols_per_line", minimum=1, maximum=500)
    common = {
        "symbols_per_line": per_line,
        "show_numbers": _flag(show_numbers, "show_numbers"),
        "labels": row_labels,
    }
    figure, axes = pyplot.subplots(
        figsize=(8.0, 1.0 + 0.4 * len(alignment.sequences) * max(1, -(-columns // per_line)))
    )
    try:
        if way == "type":
            graphics.plot_alignment_type_based(
                axes,
                alignment,
                color_scheme="clustalx" if kind == "dna" else None,
                **common,
            )
        elif way == "similarity":
            graphics.plot_alignment_similarity_based(axes, alignment, **common)
        else:
            if values is None:
                raise ConfigurationError("coloring='values' needs the values matrix.")
            array = numpy.asarray([[float(v) for v in row] for row in values])
            if array.shape != (columns, len(alignment.sequences)):
                raise ConfigurationError(
                    f"values must have shape ({columns}, {len(alignment.sequences)})."
                )
            graphics.plot_alignment_array(axes, alignment, array, **common)
        figure.savefig(out, dpi=_int(dpi, "dpi", minimum=30, maximum=600), bbox_inches="tight")
    finally:
        pyplot.close(figure)
    return {
        "status": "completed",
        "operation": "plot_sequence_alignment",
        "image": str(out),
        "coloring": way,
        "aligned": [str(item) for item in alignment.get_gapped_sequences()],
        "columns": columns,
        "provenance": _provenance(method=f"plot_alignment ({way})"),
    }


def write_a3m_alignment(
    aligned: Sequence[str],
    output_path: str,
    *,
    seq_type: str = "protein",
    query_label: str = "query",
    target_labels: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write an alignment of one query against targets as an A3M file (``set_a3m_alignments``).

    ``aligned[0]`` is the query and every other string a target aligned to it. A3M lower-cases
    target residues that are inserted relative to the query and drops the gaps of the query.
    """

    kind = _choice(seq_type, _SEQ_TYPES, "seq_type")
    items = _gapped(aligned, "aligned", minimum=2)
    labels = (
        [f"target_{i}" for i in range(1, len(items))]
        if target_labels is None
        else [str(item) for item in target_labels]
    )
    if len(labels) != len(items) - 1:
        raise ConfigurationError("target_labels must list one label per target.")
    out = Path(output_path).expanduser().resolve()
    if out.exists() or not out.parent.is_dir():
        raise ConfigurationError(f"output_path must be a new file in an existing folder: {out}")
    fasta = require("biotite.sequence.io.fasta", _EXTRA)
    _, align = _modules()
    alignments = []
    for target in items[1:]:
        alignments.append(align.Alignment.from_strings([items[0], target], _factory(kind)))
    file = fasta.FastaFile()
    fasta.set_a3m_alignments(file, alignments, _text(query_label, "query_label", limit=200), labels)
    file.write(str(out))
    return {
        "status": "completed",
        "operation": "write_a3m_alignment",
        "path": str(out),
        "targets": len(labels),
        "provenance": _provenance(method="fasta.set_a3m_alignments"),
    }


__all__ = [
    "align_banded_pairwise",
    "align_local_gapped_from_seed",
    "align_multiple_progressive",
    "align_ungapped_pairwise",
    "alignment_evalue_parameters",
    "alignment_from_strings",
    "alignment_log_evalue",
    "alignment_pairwise_identity",
    "alignment_remove_gaps",
    "alignment_remove_terminal_gaps",
    "alignment_score",
    "alignment_sequence_identity",
    "codon_table_is_start_codon",
    "kmer_index",
    "kmer_match",
    "mincode_selection",
    "plot_sequence_alignment",
    "seed_extension_ungapped",
    "sequence_profile_from_alignment",
    "sequence_profile_probability",
    "sequence_profile_score",
    "similar_kmers",
    "syncmer_selection",
    "write_a3m_alignment",
]
_register_public_verbose(globals(), tuple(__all__))
