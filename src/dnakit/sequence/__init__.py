"""Deterministic candidate implementations for sequence capabilities.

These functions accept strings or DNAKit ``DNASequence`` values and return
plain Python values that can be serialized by a future MCP adapter.  The
implementation deliberately avoids optional dependencies so that the first
batch can be validated in a clean local environment.
"""

from __future__ import annotations

import math
import random
from collections import Counter, defaultdict
from typing import Iterable, Mapping, Sequence


DNA_BASES = "ACGT"
RNA_BASES = "ACGU"
PROTEIN_BASES = set("ABCDEFGHIKLMNPQRSTVWXYZ*OU")
IUPAC_DNA = set("ACGTRYSWKMBDHVN")
IUPAC_RNA = set("ACGURYSWKMBDHVN")
DNA_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
)

CODON_TABLE = {
    codon: amino_acid
    for codon, amino_acid in (
        ("TTT", "F"), ("TTC", "F"), ("TTA", "L"), ("TTG", "L"),
        ("TCT", "S"), ("TCC", "S"), ("TCA", "S"), ("TCG", "S"),
        ("TAT", "Y"), ("TAC", "Y"), ("TAA", "*"), ("TAG", "*"),
        ("TGT", "C"), ("TGC", "C"), ("TGA", "*"), ("TGG", "W"),
        ("CTT", "L"), ("CTC", "L"), ("CTA", "L"), ("CTG", "L"),
        ("CCT", "P"), ("CCC", "P"), ("CCA", "P"), ("CCG", "P"),
        ("CAT", "H"), ("CAC", "H"), ("CAA", "Q"), ("CAG", "Q"),
        ("CGT", "R"), ("CGC", "R"), ("CGA", "R"), ("CGG", "R"),
        ("ATT", "I"), ("ATC", "I"), ("ATA", "I"), ("ATG", "M"),
        ("ACT", "T"), ("ACC", "T"), ("ACA", "T"), ("ACG", "T"),
        ("AAT", "N"), ("AAC", "N"), ("AAA", "K"), ("AAG", "K"),
        ("AGT", "S"), ("AGC", "S"), ("AGA", "R"), ("AGG", "R"),
        ("GTT", "V"), ("GTC", "V"), ("GTA", "V"), ("GTG", "V"),
        ("GCT", "A"), ("GCC", "A"), ("GCA", "A"), ("GCG", "A"),
        ("GAT", "D"), ("GAC", "D"), ("GAA", "E"), ("GAG", "E"),
        ("GGT", "G"), ("GGC", "G"), ("GGA", "G"), ("GGG", "G"),
    )
}
STOP_CODONS = {codon for codon, amino_acid in CODON_TABLE.items() if amino_acid == "*"}


def _text(sequence: object, *, uppercase: bool = True) -> str:
    """Coerce a DNAKit sequence value or text without changing whitespace."""

    if hasattr(sequence, "symbols"):
        value = str(getattr(sequence, "symbols"))
    elif isinstance(sequence, str):
        value = sequence
    else:
        raise TypeError("sequence must be a string or a DNASequence-like value")
    return value.upper() if uppercase else value


def _clean(sequence: object, *, alphabet: set[str], remove_whitespace: bool = True) -> str:
    value = _text(sequence)
    if remove_whitespace:
        value = "".join(value.split())
    invalid = sorted(set(value) - alphabet)
    if invalid:
        raise ValueError(f"invalid symbols: {', '.join(invalid)}")
    return value


def _seeded(seed: int | None) -> random.Random:
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise TypeError("seed must be an integer or None")
    return random.Random(seed)


def _kmers(sequence: str, k: int) -> Iterable[tuple[int, str]]:
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise ValueError("k must be a positive integer")
    return ((index, sequence[index : index + k]) for index in range(max(0, len(sequence) - k + 1)))


def _identity(left: str, right: str) -> float:
    if not left or not right:
        return 0.0
    return sum(a == b for a, b in zip(left, right)) / max(len(left), len(right))


def _entropy(value: str) -> float:
    if not value:
        return 0.0
    counts = Counter(value)
    total = len(value)
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def apply_reverse_complement_augmentation(sequence: str, *, include_original: bool = True) -> list[str]:
    """Return an optional original sequence and its reverse complement."""

    value = normalize_dna_sequence(sequence)
    result = [value] if include_original else []
    result.append(reverse_complement_dna(value))
    return result


def back_transcribe_rna_to_dna(sequence: str) -> str:
    """Convert normalized RNA text to DNA text."""

    return _clean(sequence, alphabet=IUPAC_RNA).replace("U", "T")


def build_minimizer_index(sequence: str, *, k: int = 15, window: int = 10) -> dict[str, list[int]]:
    """Index the leftmost lexicographic k-mer in each sliding window."""

    value = normalize_dna_sequence(sequence)
    if window < 1:
        raise ValueError("window must be positive")
    kmers = list(_kmers(value, k))
    index: dict[str, list[int]] = defaultdict(list)
    if len(kmers) < window:
        windows = [kmers] if kmers else []
    else:
        windows = (kmers[start : start + window] for start in range(len(kmers) - window + 1))
    for current in windows:
        if not current:
            continue
        position, minimizer = min(current, key=lambda item: (item[1], item[0]))
        if not index[minimizer] or index[minimizer][-1] != position:
            index[minimizer].append(position)
    return dict(index)


def calculate_sequence_metrics(sequence: str) -> dict[str, object]:
    """Calculate length, composition, GC content, ambiguity and entropy."""

    value = _clean(sequence, alphabet=IUPAC_DNA)
    counts = Counter(value)
    gc = (counts["G"] + counts["C"]) / len(value) if value else 0.0
    return {
        "length": len(value),
        "gc_fraction": gc,
        "gc_percent": gc * 100.0,
        "n_count": counts["N"],
        "ambiguous_count": sum(counts[base] for base in IUPAC_DNA - set(DNA_BASES)),
        "entropy_bits": _entropy(value),
        "composition": dict(sorted(counts.items())),
    }


def cluster_sequences_by_identity(sequences: Sequence[str], *, threshold: float = 0.9) -> list[list[int]]:
    """Greedily group sequences by identity to the first sequence in each cluster."""

    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must be between 0 and 1")
    values = [normalize_dna_sequence(value) for value in sequences]
    clusters: list[list[int]] = []
    representatives: list[str] = []
    for index, value in enumerate(values):
        for cluster_index, representative in enumerate(representatives):
            if _identity(value, representative) >= threshold:
                clusters[cluster_index].append(index)
                break
        else:
            representatives.append(value)
            clusters.append([index])
    return clusters


def compare_motif_alleles(left: str, right: str, motif: str, *, start: int = 0) -> dict[str, object]:
    """Compare a motif interval between two allele sequences."""

    motif_value = normalize_dna_sequence(motif)
    left_value = normalize_dna_sequence(left)
    right_value = normalize_dna_sequence(right)
    end = start + len(motif_value)
    left_hit = left_value[start:end]
    right_hit = right_value[start:end]
    return {
        "start": start,
        "end": end,
        "motif": motif_value,
        "left": left_hit,
        "right": right_hit,
        "left_matches": left_hit == motif_value,
        "right_matches": right_hit == motif_value,
        "changed": left_hit != right_hit,
    }


def compare_sequence_sketches(left: str, right: str, *, k: int = 7) -> dict[str, object]:
    """Compare k-mer sets with a deterministic Jaccard sketch."""

    left_set = {kmer for _, kmer in _kmers(normalize_dna_sequence(left), k)}
    right_set = {kmer for _, kmer in _kmers(normalize_dna_sequence(right), k)}
    union = left_set | right_set
    return {
        "k": k,
        "left_size": len(left_set),
        "right_size": len(right_set),
        "intersection_size": len(left_set & right_set),
        "jaccard": len(left_set & right_set) / len(union) if union else 1.0,
    }


def complement_dna(sequence: str) -> str:
    """Return the IUPAC complement without reversing the sequence."""

    return _clean(sequence, alphabet=IUPAC_DNA).translate(DNA_COMPLEMENT)


def count_canonical_kmers(sequence: str, *, k: int = 5) -> dict[str, int]:
    """Count each k-mer under lexicographic canonicalization with its reverse complement."""

    counts: Counter[str] = Counter()
    for _, kmer in _kmers(normalize_dna_sequence(sequence), k):
        counts[min(kmer, reverse_complement_dna(kmer))] += 1
    return dict(sorted(counts.items()))


def count_kmers(sequence: str, *, k: int = 5) -> dict[str, int]:
    """Count overlapping DNA k-mers."""

    return dict(sorted(Counter(kmer for _, kmer in _kmers(normalize_dna_sequence(sequence), k)).items()))


def count_kmers_batch(sequences: Sequence[str], *, k: int = 5) -> list[dict[str, int]]:
    """Count k-mers independently for each sequence."""

    return [count_kmers(sequence, k=k) for sequence in sequences]


def deduplicate_sequences(sequences: Sequence[str], *, canonical: bool = False) -> list[str]:
    """Remove duplicate sequences while preserving first occurrence order."""

    result: list[str] = []
    seen: set[str] = set()
    for sequence in sequences:
        value = normalize_dna_sequence(sequence)
        key = min(value, reverse_complement_dna(value)) if canonical else value
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def digest_sequence_in_silico(sequence: str, enzymes: Mapping[str, str]) -> dict[str, object]:
    """Find exact restriction motifs and return fragment lengths."""

    value = normalize_dna_sequence(sequence)
    sites: dict[str, list[int]] = {}
    cuts: list[int] = [0, len(value)]
    for name, motif in enzymes.items():
        motif_value = normalize_dna_sequence(motif)
        positions = [index for index, _ in _kmers(value, len(motif_value)) if value[index : index + len(motif_value)] == motif_value]
        sites[name] = positions
        cuts.extend(positions)
    ordered = sorted(set(cuts))
    return {"sites": sites, "cut_positions": ordered, "fragment_lengths": [b - a for a, b in zip(ordered, ordered[1:])]}


def extract_cds_from_transcript(transcript: str, intervals: Sequence[Sequence[int]]) -> str:
    """Extract and concatenate zero-based half-open CDS intervals."""

    value = normalize_dna_sequence(transcript)
    pieces = []
    for interval in intervals:
        if len(interval) != 2:
            raise ValueError("each interval must contain start and end")
        start, end = interval
        if start < 0 or end < start or end > len(value):
            raise ValueError("CDS interval is outside the transcript")
        pieces.append(value[start:end])
    return "".join(pieces)


def find_sequence_patterns(sequence: str, *, motif: str | None = None, min_repeat: int = 3) -> dict[str, object]:
    """Find motif hits, homopolymers and simple tandem repeats."""

    value = normalize_dna_sequence(sequence)
    result: dict[str, object] = {"motif_hits": [], "homopolymers": [], "tandem_repeats": []}
    if motif:
        motif_value = normalize_dna_sequence(motif)
        result["motif_hits"] = [index for index, _ in _kmers(value, len(motif_value)) if value[index : index + len(motif_value)] == motif_value]
    homopolymers = []
    index = 0
    while index < len(value):
        end = index + 1
        while end < len(value) and value[end] == value[index]:
            end += 1
        if end - index >= min_repeat:
            homopolymers.append({"start": index, "end": end, "base": value[index], "length": end - index})
        index = end
    result["homopolymers"] = homopolymers
    for unit_length in range(1, min(6, len(value) + 1)):
        for start in range(0, len(value) - unit_length * min_repeat + 1):
            unit = value[start : start + unit_length]
            count = 1
            while value[start + count * unit_length : start + (count + 1) * unit_length] == unit:
                count += 1
            if count >= min_repeat:
                result["tandem_repeats"].append({"start": start, "end": start + count * unit_length, "unit": unit, "count": count})
    return result


def generate_kmer_table(sequences: Sequence[str], *, k: int = 5) -> list[dict[str, object]]:
    """Return a row-oriented k-mer count table."""

    rows = []
    for sequence_id, sequence in enumerate(sequences):
        for kmer, count in count_kmers(sequence, k=k).items():
            rows.append({"sequence_index": sequence_id, "kmer": kmer, "count": count})
    return rows


def generate_markov_sequence(sequence: str, *, length: int, order: int = 1, seed: int | None = None) -> str:
    """Generate a sequence from an order-n empirical Markov model."""

    value = normalize_dna_sequence(sequence)
    if length < 0 or order < 0 or len(value) <= order:
        raise ValueError("length must be non-negative and source must be longer than order")
    rng = _seeded(seed)
    if length <= order:
        return value[:length]
    contexts: dict[str, list[str]] = defaultdict(list)
    for index in range(len(value) - order):
        contexts[value[index : index + order]].append(value[index + order])
    result = list(value[:order])
    while len(result) < length:
        context = "".join(result[-order:])
        choices = contexts.get(context)
        result.append(rng.choice(choices or list(DNA_BASES)))
    return "".join(result)


def generate_random_dna(*, length: int, seed: int | None = None, alphabet: str = DNA_BASES) -> str:
    """Generate reproducible random DNA text."""

    if length < 0 or not alphabet or any(symbol not in DNA_BASES for symbol in alphabet):
        raise ValueError("length must be non-negative and alphabet must contain A, C, G or T")
    return "".join(_seeded(seed).choice(alphabet) for _ in range(length))


def mask_low_complexity(sequence: str, *, window: int = 10, entropy_threshold: float = 1.2, mask: str = "N") -> str:
    """Mask windows whose Shannon entropy is below a threshold."""

    value = normalize_dna_sequence(sequence)
    if window < 1 or len(mask) != 1:
        raise ValueError("window must be positive and mask must be one symbol")
    result = list(value)
    for start in range(0, max(0, len(value) - window + 1)):
        if _entropy(value[start : start + window]) < entropy_threshold:
            result[start : start + window] = [mask] * window
    return "".join(result)


def merge_overlapping_windows(windows: Sequence[Sequence[int]]) -> list[list[int]]:
    """Merge sorted or unsorted zero-based half-open intervals."""

    normalized = sorted((int(start), int(end)) for start, end in windows)
    merged: list[list[int]] = []
    for start, end in normalized:
        if end < start:
            raise ValueError("window end must not be less than start")
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def mutate_sequence_at_positions(sequence: str, mutations: Mapping[str, str]) -> str:
    """Apply zero-based substitutions, accepting JSON object keys as positions.

    Python callers may still provide integer keys.  JSON object keys arrive as
    strings, so the public annotation uses ``Mapping[str, str]`` and the
    implementation converts each key to an integer before validation.
    """

    value = list(normalize_dna_sequence(sequence))
    for raw_position, replacement in mutations.items():
        try:
            position = int(raw_position)
        except (TypeError, ValueError) as exc:
            raise ValueError("mutation positions must be integers") from exc
        if position < 0 or position >= len(value) or replacement not in DNA_BASES:
            raise ValueError("mutation position or replacement is invalid")
        value[position] = replacement
    return "".join(value)


def normalize_dna_sequence(sequence: str) -> str:
    """Uppercase DNA, remove whitespace and convert U to T."""

    value = "".join(_text(sequence).split()).replace("U", "T")
    invalid = sorted(set(value) - IUPAC_DNA)
    if invalid:
        raise ValueError(f"invalid DNA symbols: {', '.join(invalid)}")
    return value


def normalize_protein_sequence(sequence: str) -> str:
    """Uppercase and whitespace-normalize a protein sequence."""

    value = "".join(_text(sequence).split())
    invalid = sorted(set(value) - PROTEIN_BASES)
    if invalid:
        raise ValueError(f"invalid protein symbols: {', '.join(invalid)}")
    return value


def normalize_rna_sequence(sequence: str) -> str:
    """Uppercase RNA, remove whitespace and convert T to U."""

    value = "".join(_text(sequence).split()).replace("T", "U")
    invalid = sorted(set(value) - IUPAC_RNA)
    if invalid:
        raise ValueError(f"invalid RNA symbols: {', '.join(invalid)}")
    return value


def one_hot_encode_dna(sequence: str, *, include_ambiguous: bool = False) -> list[list[int]]:
    """Encode DNA symbols in A,C,G,T order."""

    value = normalize_dna_sequence(sequence)
    alphabet = "ACGTN" if include_ambiguous else "ACGT"
    return [[int(symbol == base) for base in alphabet] for symbol in value]


def one_hot_encode_rna(sequence: str, *, include_ambiguous: bool = False) -> list[list[int]]:
    """Encode RNA symbols in A,C,G,U order."""

    value = normalize_rna_sequence(sequence)
    alphabet = "ACGUN" if include_ambiguous else "ACGU"
    return [[int(symbol == base) for base in alphabet] for symbol in value]


def pad_sequence(sequence: str, *, length: int, symbol: str = "N", side: str = "right") -> str:
    """Pad a sequence to an exact length."""

    value = normalize_dna_sequence(sequence)
    if length < len(value) or len(symbol) != 1 or side not in {"left", "right"}:
        raise ValueError("length, symbol or side is invalid")
    padding = symbol * (length - len(value))
    return padding + value if side == "left" else value + padding


def predict_coding_regions(sequence: str, *, min_aa: int = 1) -> list[dict[str, object]]:
    """Find simple ATG-to-stop open reading frames on the forward strand."""

    value = normalize_dna_sequence(sequence)
    result = []
    for frame in range(3):
        index = frame
        while index + 3 <= len(value):
            if value[index : index + 3] != "ATG":
                index += 3
                continue
            stop = index + 3
            while stop + 3 <= len(value) and value[stop : stop + 3] not in STOP_CODONS:
                stop += 3
            if stop + 3 <= len(value):
                protein = translate_dna_to_protein(value[index : stop + 3])
                if len(protein.rstrip("*")) >= min_aa:
                    result.append({"start": index, "end": stop + 3, "frame": frame, "protein": protein})
                index = stop + 3
            else:
                break
    return result


def remove_sequence_gaps(sequence: str) -> str:
    """Remove alignment gap characters from sequence text."""

    return _text(sequence).replace("-", "").replace(".", "")


def remove_sequence_whitespace(sequence: str) -> str:
    """Remove all Unicode whitespace from sequence text."""

    return "".join(_text(sequence, uppercase=False).split())


def replace_ambiguous_bases(sequence: str, *, replacement: str = "N") -> str:
    """Replace non-ACGT IUPAC bases with a chosen ambiguity symbol."""

    value = normalize_dna_sequence(sequence)
    if len(replacement) != 1 or replacement not in IUPAC_DNA:
        raise ValueError("replacement must be one IUPAC DNA symbol")
    return "".join(symbol if symbol in DNA_BASES else replacement for symbol in value)


def reverse_complement_dna(sequence: str) -> str:
    """Return the reverse complement of IUPAC DNA text."""

    return complement_dna(sequence)[::-1]


def reverse_sequence(sequence: str) -> str:
    """Reverse DNA text without complementing it."""

    return normalize_dna_sequence(sequence)[::-1]


def scan_pwm(sequence: str, pwm: Sequence[Mapping[str, float]], *, threshold: float = 0.0) -> list[dict[str, object]]:
    """Scan a position-weight matrix represented by base-score mappings."""

    value = normalize_dna_sequence(sequence)
    if not pwm:
        raise ValueError("pwm must not be empty")
    hits = []
    width = len(pwm)
    for start in range(len(value) - width + 1):
        window = value[start : start + width]
        score = sum(float(pwm[offset].get(base, float("-inf"))) for offset, base in enumerate(window))
        if score >= threshold:
            hits.append({"start": start, "end": start + width, "sequence": window, "score": score})
    return hits


def scan_pwm_batch(sequences: Sequence[str], pwm: Sequence[Mapping[str, float]], *, threshold: float = 0.0) -> list[list[dict[str, object]]]:
    """Scan one PWM across multiple sequences."""

    return [scan_pwm(sequence, pwm, threshold=threshold) for sequence in sequences]


def score_motif_hit(sequence: str, motif: str, *, position: int = 0, match_score: float = 1.0, mismatch_score: float = 0.0) -> dict[str, object]:
    """Score a motif at a zero-based sequence position."""

    value = normalize_dna_sequence(sequence)
    pattern = normalize_dna_sequence(motif)
    observed = value[position : position + len(pattern)]
    score = sum(match_score if left == right else mismatch_score for left, right in zip(observed, pattern))
    return {"position": position, "motif": pattern, "observed": observed, "score": score, "matches": observed == pattern}


def select_representative_sequences(sequences: Sequence[str], *, threshold: float = 0.9) -> list[int]:
    """Return indices of greedy identity-cluster representatives."""

    clusters = cluster_sequences_by_identity(sequences, threshold=threshold)
    return [cluster[0] for cluster in clusters]


def simulate_indels(sequence: str, *, count: int = 1, seed: int | None = None, max_length: int | None = None) -> str:
    """Apply reproducible single-base insertions or deletions."""

    value = list(normalize_dna_sequence(sequence))
    rng = _seeded(seed)
    for _ in range(count):
        if not value or rng.random() < 0.5:
            position = rng.randrange(len(value) + 1)
            value.insert(position, rng.choice(DNA_BASES))
        else:
            del value[rng.randrange(len(value))]
    result = "".join(value)
    if max_length is not None and len(result) > max_length:
        result = result[:max_length]
    return result


def simulate_point_mutations(sequence: str, *, rate: float = 0.01, seed: int | None = None) -> str:
    """Apply independent seeded substitutions at the requested rate."""

    if not 0.0 <= rate <= 1.0:
        raise ValueError("rate must be between 0 and 1")
    rng = _seeded(seed)
    value = []
    for symbol in normalize_dna_sequence(sequence):
        if rng.random() < rate:
            value.append(rng.choice([base for base in DNA_BASES if base != symbol]))
        else:
            value.append(symbol)
    return "".join(value)


def simulate_recombination(left: str, right: str, *, crossover: int | None = None) -> str:
    """Create a recombinant sequence from one crossover point."""

    first = normalize_dna_sequence(left)
    second = normalize_dna_sequence(right)
    point = min(len(first), len(second)) // 2 if crossover is None else crossover
    if point < 0 or point > len(first) or point > len(second):
        raise ValueError("crossover is outside one of the parent sequences")
    return first[:point] + second[point:]


def simulate_sequence_evolution(sequence: str, *, generations: int = 1, mutation_rate: float = 0.01, seed: int | None = None) -> list[str]:
    """Return the starting sequence and each seeded mutation generation."""

    rng = _seeded(seed)
    generations_out = [normalize_dna_sequence(sequence)]
    for _ in range(generations):
        generations_out.append(simulate_point_mutations(generations_out[-1], rate=mutation_rate, seed=rng.randrange(2**32)))
    return generations_out


def split_sequence_windows(sequence: str, *, window: int, step: int | None = None) -> list[dict[str, object]]:
    """Split a sequence into overlapping or adjacent windows."""

    value = normalize_dna_sequence(sequence)
    if window < 1:
        raise ValueError("window must be positive")
    step = window if step is None else step
    if step < 1:
        raise ValueError("step must be positive")
    return [{"start": start, "end": min(start + window, len(value)), "sequence": value[start : start + window]} for start in range(0, len(value), step) if start < len(value)]


def tokenize_dna_bpe(sequence: str, *, merges: Sequence[Sequence[str]] = ()) -> list[str]:
    """Tokenize DNA and apply deterministic pair merges in order."""

    tokens = list(normalize_dna_sequence(sequence))
    for merge in merges:
        if len(merge) != 2:
            raise ValueError("each merge must contain two tokens")
        index = 0
        merged = []
        while index < len(tokens):
            if index + 1 < len(tokens) and (tokens[index], tokens[index + 1]) == tuple(merge):
                merged.append(tokens[index] + tokens[index + 1])
                index += 2
            else:
                merged.append(tokens[index])
                index += 1
        tokens = merged
    return tokens


def tokenize_dna_kmer(sequence: str, *, k: int = 3, stride: int = 1) -> list[str]:
    """Return overlapping or strided DNA k-mer tokens."""

    if stride < 1:
        raise ValueError("stride must be positive")
    value = normalize_dna_sequence(sequence)
    return [value[index : index + k] for index in range(0, max(0, len(value) - k + 1), stride)]


def transcribe_dna_to_rna(sequence: str) -> str:
    """Transcribe DNA text to RNA text."""

    return normalize_dna_sequence(sequence).replace("T", "U")


def translate_codon(codon: str, *, unknown: str = "X") -> str:
    """Translate one DNA or RNA codon."""

    value = normalize_dna_sequence(codon)
    if len(value) != 3:
        raise ValueError("codon must have length three")
    if value in CODON_TABLE:
        return CODON_TABLE[value]
    return unknown


def translate_dna_to_protein(sequence: str, *, frame: int = 0, stop_at_stop: bool = False) -> str:
    """Translate DNA using frame zero, one or two and the standard code."""

    value = normalize_dna_sequence(sequence)
    if frame not in (0, 1, 2):
        raise ValueError("frame must be 0, 1 or 2")
    amino_acids = []
    for index in range(frame, len(value) - 2, 3):
        amino_acid = CODON_TABLE.get(value[index : index + 3], "X")
        amino_acids.append(amino_acid)
        if stop_at_stop and amino_acid == "*":
            break
    return "".join(amino_acids)


def translate_orf(sequence: str, *, start_codon: str = "ATG", stop_codons: Iterable[str] = STOP_CODONS) -> list[dict[str, object]]:
    """Translate all in-frame ORFs on the forward strand."""

    value = normalize_dna_sequence(sequence)
    starts = normalize_dna_sequence(start_codon)
    stops = {normalize_dna_sequence(codon) for codon in stop_codons}
    result = []
    for start in range(len(value) - 2):
        if value[start : start + 3] != starts:
            continue
        for end in range(start + 3, len(value) - 2, 3):
            if value[end : end + 3] in stops:
                result.append({"start": start, "end": end + 3, "sequence": value[start : end + 3], "protein": translate_dna_to_protein(value[start : end + 3], stop_at_stop=True)})
                break
    return result


def translate_rna_to_protein(sequence: str, *, frame: int = 0, stop_at_stop: bool = False) -> str:
    """Translate RNA text through the DNA codon table."""

    return translate_dna_to_protein(back_transcribe_rna_to_dna(sequence), frame=frame, stop_at_stop=stop_at_stop)


def translate_six_frames(sequence: str) -> dict[str, str]:
    """Translate three forward and three reverse-complement frames."""

    value = normalize_dna_sequence(sequence)
    reverse = reverse_complement_dna(value)
    result = {f"+{frame + 1}": translate_dna_to_protein(value, frame=frame) for frame in range(3)}
    result.update({f"-{frame + 1}": translate_dna_to_protein(reverse, frame=frame) for frame in range(3)})
    return result


def trim_sequence(sequence: str, *, left: int = 0, right: int = 0) -> str:
    """Trim a number of bases from the left and right ends."""

    value = normalize_dna_sequence(sequence)
    if left < 0 or right < 0 or left + right > len(value):
        raise ValueError("trim lengths are invalid")
    return value[left : len(value) - right if right else None]


def validate_dna_alphabet(sequence: str, *, allow_ambiguous: bool = True) -> dict[str, object]:
    """Validate DNA symbols and return a serializable report."""

    value = "".join(_text(sequence).split())
    alphabet = IUPAC_DNA if allow_ambiguous else set(DNA_BASES)
    invalid = sorted(set(value) - alphabet)
    return {"valid": not invalid, "length": len(value), "invalid_symbols": invalid, "alphabet": "iupac" if allow_ambiguous else "strict"}


def validate_protein_alphabet(sequence: str) -> dict[str, object]:
    """Validate a protein alphabet and return a serializable report."""

    value = "".join(_text(sequence).split())
    invalid = sorted(set(value) - PROTEIN_BASES)
    return {"valid": not invalid, "length": len(value), "invalid_symbols": invalid, "alphabet": "protein"}


def validate_rna_alphabet(sequence: str, *, allow_ambiguous: bool = True) -> dict[str, object]:
    """Validate RNA symbols and return a serializable report."""

    value = "".join(_text(sequence).split())
    alphabet = IUPAC_RNA if allow_ambiguous else set(RNA_BASES)
    invalid = sorted(set(value) - alphabet)
    return {"valid": not invalid, "length": len(value), "invalid_symbols": invalid, "alphabet": "iupac" if allow_ambiguous else "strict"}


__all__ = [
    "apply_reverse_complement_augmentation", "back_transcribe_rna_to_dna", "build_minimizer_index",
    "calculate_sequence_metrics", "cluster_sequences_by_identity", "compare_motif_alleles",
    "compare_sequence_sketches", "complement_dna", "count_canonical_kmers", "count_kmers",
    "count_kmers_batch", "deduplicate_sequences", "digest_sequence_in_silico",
    "extract_cds_from_transcript", "find_sequence_patterns", "generate_kmer_table",
    "generate_markov_sequence", "generate_random_dna", "mask_low_complexity", "merge_overlapping_windows",
    "mutate_sequence_at_positions", "normalize_dna_sequence", "normalize_protein_sequence",
    "normalize_rna_sequence", "one_hot_encode_dna", "one_hot_encode_rna", "pad_sequence",
    "predict_coding_regions", "remove_sequence_gaps", "remove_sequence_whitespace",
    "replace_ambiguous_bases", "reverse_complement_dna", "reverse_sequence", "scan_pwm",
    "scan_pwm_batch", "score_motif_hit", "select_representative_sequences", "simulate_indels",
    "simulate_point_mutations", "simulate_recombination", "simulate_sequence_evolution",
    "split_sequence_windows", "sequence_length_gc_summary", "tokenize_dna_bpe", "tokenize_dna_kmer", "transcribe_dna_to_rna",
    "translate_codon", "translate_dna_to_protein", "translate_orf", "translate_rna_to_protein",
    "translate_six_frames", "trim_sequence", "validate_dna_alphabet", "validate_protein_alphabet",
    "validate_rna_alphabet",
]


from .local import sequence_length_gc_summary

from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    build_screed_database,
    calculate_kmer_frequencies,
    convert_sequence_file,
    detect_fastq_quality_encoding,
    expand_degenerate_sequences,
    extract_feature_sequences,
    fetch_screed_records,
    fetch_twobit_region,
    read_sequence_records,
    summarize_twobit,
    translate_dna_with_genetic_code,
)

__all__ = [*__all__, "build_screed_database", "calculate_kmer_frequencies", "convert_sequence_file", "detect_fastq_quality_encoding", "expand_degenerate_sequences", "extract_feature_sequences", "fetch_screed_records", "fetch_twobit_region", "read_sequence_records", "summarize_twobit", "translate_dna_with_genetic_code"]

# Full-pass exports (see docs/api/features)
from .Bio_functions import (
    ambiguous_codon_amino_acids,
    bootstrap_consensus_tree,
    build_codon_alignment,
    convert_quality_score,
    expand_ambiguous_codons,
    gc_content_by_codon_position,
    gc_skew_windows,
    load_substitution_matrix,
    local_composition_complexity,
    mcdonald_kreitman_test,
    pairwise_dn_ds_matrix,
    parse_kegg_gene_entries,
    parse_primer3_output,
    parse_primersearch_output,
    read_ace_assembly,
    read_gaf_annotations,
    read_genepop_genotypes,
    read_phd_sequence,
    read_psl_alignments,
    sequence_checksum,
)

from .dnaio_functions import (
    check_read_mates,
    read_fastx_records,
)

from .pyfaidx_functions import (
    split_bed_entry,
    split_ucsc_region,
    translate_chromosome_names,
)

from .screed_functions import (
    export_screed_database,
    format_fastx_records,
    read_screed_fastx_records,
)

from .skbio_functions import (
    alignment_distance_matrix,
    alignment_score,
    multiple_nucleotide_alignment,
    nucleotide_distance,
    nucleotide_transition_matrix,
    sniff_sequence_format,
)

__all__ = [*__all__, "ambiguous_codon_amino_acids", "bootstrap_consensus_tree", "build_codon_alignment", "convert_quality_score", "expand_ambiguous_codons", "gc_content_by_codon_position", "gc_skew_windows", "load_substitution_matrix", "local_composition_complexity", "mcdonald_kreitman_test", "pairwise_dn_ds_matrix", "parse_kegg_gene_entries", "parse_primer3_output", "parse_primersearch_output", "read_ace_assembly", "read_gaf_annotations", "read_genepop_genotypes", "read_phd_sequence", "read_psl_alignments", "sequence_checksum", "check_read_mates", "read_fastx_records", "split_bed_entry", "split_ucsc_region", "translate_chromosome_names", "export_screed_database", "format_fastx_records", "read_screed_fastx_records", "alignment_distance_matrix", "alignment_score", "multiple_nucleotide_alignment", "nucleotide_distance", "nucleotide_transition_matrix", "sniff_sequence_format"]

# Review-pass exports (see docs/api/features)
from .Bio_review_functions import (
    adam_consensus_tree,
    bootstrap_branch_support,
    codon_back_table,
    complement_rna_sequence,
    convert_search_results,
    melting_temperature_gc,
    read_fasta_qual_pairs,
    read_mast_report,
    reverse_complement_rna_sequence,
    tree_to_graph_edges,
    write_gaf_file,
    write_motifs,
    write_pairwise_alignment_file,
    write_qual_file,
)

from .skbio_review_functions import (
    build_bme_tree,
    build_gme_tree,
    improve_tree_nni,
    tree_path_distances,
    weighted_robinson_foulds_distances,
)

__all__ = [*__all__, "adam_consensus_tree", "bootstrap_branch_support", "codon_back_table", "complement_rna_sequence", "convert_search_results", "melting_temperature_gc", "read_fasta_qual_pairs", "read_mast_report", "reverse_complement_rna_sequence", "tree_to_graph_edges", "write_gaf_file", "write_motifs", "write_pairwise_alignment_file", "write_qual_file", "build_bme_tree", "build_gme_tree", "improve_tree_nni", "tree_path_distances", "weighted_robinson_foulds_distances"]

# Family-pass exports (see docs/api/features)
from .Bio_family_functions import (
    bootstrap_replicate_trees,
    parse_primer3_records,
    read_alignment_blocks,
    read_alignment_file,
    read_blast_records,
    read_genbank_record_fields,
    read_goa_gene_product_file,
    read_phd_records,
    read_phylogeny_collection,
    read_substitution_matrix_file,
)

from .skbio_family_functions import (
    correct_jc69_distance,
)

__all__ = [*__all__, "bootstrap_replicate_trees", "parse_primer3_records", "read_alignment_blocks", "read_alignment_file", "read_blast_records", "read_genbank_record_fields", "read_goa_gene_product_file", "read_phd_records", "read_phylogeny_collection", "read_substitution_matrix_file", "correct_jc69_distance"]

_register_public_verbose(globals())
