"""Native, versioned DNA encodings; no iFeatureOmega or iLearn dependency.

All words use A,C,G,T order. Ambiguity is rejected rather than silently masked.
Physicochemical tables are supplied explicitly as property -> word -> value.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import positive, require, sequence_text
from dnakit.exceptions import ConfigurationError

DNA_ENCODINGS = (
    "Kmer",
    "RCKmer",
    "Mismatch",
    "Subsequence",
    "NAC",
    "ANF",
    "ENAC",
    "binary",
    "PS2",
    "PS3",
    "PS4",
    "CKSNAP",
    "NCP",
    "EIIP",
    "PseEIIP",
    "ASDC",
    "DBE",
    "LPDF",
    "DPCP",
    "TPCP",
    "MMI",
    "Z_curve_9bit",
    "Z_curve_12bit",
    "Z_curve_36bit",
    "Z_curve_48bit",
    "Z_curve_144bit",
    "DAC",
    "DCC",
    "DACC",
    "TAC",
    "TCC",
    "TACC",
    "NMBroto",
    "Moran",
    "Geary",
    "PseDNC",
    "PseKNC",
    "PCPseDNC",
    "PCPseTNC",
    "SCPseDNC",
    "SCPseTNC",
    "PSTNPss",
    "PSTNPds",
    "KNN",
)
_EIIP = dict(zip("ACGT", (0.1260, 0.1340, 0.0806, 0.1335), strict=True))


@dataclass(frozen=True)
class EncodingResult:
    """A rows-by-features NumPy array with stable IDs and feature names."""

    ids: tuple[str, ...]
    names: tuple[str, ...]
    values: Any
    descriptor: str
    parameters: Mapping[str, Any]
    schema_version: str = "dnakit.dna-encoding.v1"


def _words(k: int) -> tuple[str, ...]:
    return tuple("".join(word) for word in itertools.product("ACGT", repeat=k))


def _rc(text: str) -> str:
    return text.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _frequencies(text: str, k: int, normalized: bool = True) -> dict[str, float]:
    if len(text) < k:
        raise ConfigurationError(f"Sequence must contain at least {k} bases.")
    counts = Counter(text[i : i + k] for i in range(len(text) - k + 1))
    denominator = len(text) - k + 1 if normalized else 1
    return {word: counts[word] / denominator for word in _words(k)}


def _property_matrix(
    tables: Mapping[str, Mapping[str, float]], k: int
) -> tuple[tuple[str, ...], Any]:
    np = require("numpy", "bio-analysis")
    words = _words(k)
    if not tables:
        raise ConfigurationError("Provide physicochemical property_tables explicitly.")
    if any(set(table) != set(words) for table in tables.values()):
        raise ConfigurationError(f"Each property table must contain all {4**k} canonical {k}-mers.")
    values = np.asarray([[table[word] for word in words] for table in tables.values()], dtype=float)
    if not np.isfinite(values).all():
        raise ConfigurationError("Property tables must be finite.")
    sd = values.std(axis=1, keepdims=True)
    if (sd == 0).any() or not np.isfinite(sd).all():
        raise ConfigurationError("Property tables need a finite nonzero standard deviation.")
    return tuple(tables), (values - values.mean(axis=1, keepdims=True)) / sd


def _physicochemical(
    text: str,
    descriptor: str,
    parameters: Mapping[str, Any],
    tables: Mapping[str, Mapping[str, float]],
) -> dict[str, float]:
    k = 3 if descriptor in {"TPCP", "TAC", "TCC", "TACC", "PCPseTNC", "SCPseTNC"} else 2
    property_names, values = _property_matrix(tables, k)
    words = _words(k)
    lookup = {word: index for index, word in enumerate(words)}
    if len(text) < k:
        raise ConfigurationError("Sequence is shorter than the property word length.")
    profile = values[:, [lookup[text[i : i + k]] for i in range(len(text) - k + 1)]]
    if descriptor in {"DPCP", "TPCP"}:
        return {
            f"{name}:position_{pos}": float(profile[row, pos])
            for row, name in enumerate(property_names)
            for pos in range(profile.shape[1])
        }
    lag_count = positive(parameters.get("lag", 3), "lag")
    if lag_count >= profile.shape[1]:
        raise ConfigurationError("lag must be smaller than the number of property words.")
    if "Pse" in descriptor:
        weight = float(parameters.get("weight", 0.05))
        if not math.isfinite(weight) or weight < 0:
            raise ConfigurationError("Pseudo-composition weight must be finite and nonnegative.")
        composition_k = positive(parameters.get("k", 3), "k") if descriptor == "PseKNC" else k
        if composition_k > 6:
            raise ConfigurationError("Pseudo k-mer dimension is limited to k <= 6.")
        frequency = _frequencies(text, composition_k)
        correlations: dict[str, float] = {}
        for lag in range(1, lag_count + 1):
            if descriptor.startswith("SC"):
                for index, name in enumerate(property_names):
                    correlations[f"correlation:{name}:lag_{lag}"] = float(
                        (profile[index, :-lag] * profile[index, lag:]).mean()
                    )
            else:
                correlations[f"theta:lag_{lag}"] = float(
                    ((profile[:, :-lag] - profile[:, lag:]) ** 2).mean()
                )
        denominator = 1 + weight * sum(correlations.values())
        if denominator <= 0:
            raise ConfigurationError("Pseudo-composition denominator is nonpositive.")
        return {
            **{f"frequency:{key}": value / denominator for key, value in frequency.items()},
            **{key: weight * value / denominator for key, value in correlations.items()},
        }
    output: dict[str, float] = {}
    centered = profile - profile.mean(axis=1, keepdims=True)
    auto_only = descriptor in {"DAC", "TAC", "Moran", "Geary", "NMBroto"}
    cross_only = descriptor in {"DCC", "TCC"}
    if cross_only and len(property_names) < 2:
        raise ConfigurationError("Cross-covariance requires at least two property tables.")
    for lag in range(1, lag_count + 1):
        for left, name in enumerate(property_names):
            for right, other in enumerate(property_names):
                if (auto_only and left != right) or (cross_only and left == right):
                    continue
                value = float((centered[left, :-lag] * centered[right, lag:]).mean())
                variance = float((centered[left] ** 2).mean())
                if descriptor == "NMBroto":
                    value = float((profile[left, :-lag] * profile[left, lag:]).mean())
                elif descriptor in {"Moran", "Geary"}:
                    if variance == 0:
                        raise ConfigurationError(
                            f"{descriptor} is undefined for a constant profile."
                        )
                    if descriptor == "Moran":
                        value /= variance
                    else:
                        n = profile.shape[1]
                        value = float(((profile[left, :-lag] - profile[left, lag:]) ** 2).mean())
                        value *= (n - 1) / (2 * n * variance)
                output[f"{name}:{other}:lag_{lag}"] = value
    return output


def _composition(text: str, descriptor: str, parameters: Mapping[str, Any]) -> dict[str, float]:
    length = len(text)
    normalized = parameters.get("normalized", True)
    if type(normalized) is not bool:
        raise ConfigurationError("normalized must be a boolean.")
    if descriptor in {"Kmer", "RCKmer", "Mismatch", "Subsequence"}:
        k = positive(parameters.get("k", 3), "k")
        if k > 6:
            raise ConfigurationError("k is limited to 6 to bound feature dimension.")
        frequencies = _frequencies(text, k, normalized)
        if descriptor == "Kmer":
            return frequencies
        if descriptor == "RCKmer":
            result: dict[str, float] = {}
            for word, value in frequencies.items():
                canonical = min(word, _rc(word))
                result[canonical] = result.get(canonical, 0) + value
            return result
        if descriptor == "Mismatch":
            mismatches = parameters.get("mismatches", 1)
            if type(mismatches) is not int or not 0 <= mismatches < k:
                raise ConfigurationError("mismatches must be an integer in [0, k).")
            if len(frequencies) ** 2 > 2_000_000:
                raise ConfigurationError("Mismatch feature computation is limited to k <= 5.")
            return {
                word: sum(
                    value
                    for candidate, value in frequencies.items()
                    if sum(a != b for a, b in zip(word, candidate, strict=True)) <= mismatches
                )
                for word in frequencies
            }
        delta = float(parameters.get("delta", 0.5))
        if not math.isfinite(delta) or not 0 < delta <= 1:
            raise ConfigurationError("Subsequence delta must be in (0, 1].")
        if math.comb(length, k) > 1_000_000:
            raise ConfigurationError("Subsequence combinations exceed 1000000.")
        subsequence_counts = dict.fromkeys(frequencies, 0.0)
        for positions in itertools.combinations(range(length), k):
            word = "".join(text[index] for index in positions)
            subsequence_counts[word] += delta ** (positions[-1] - positions[0] + 1)
        total = sum(subsequence_counts.values()) if normalized else 1.0
        if total <= 0:
            raise ConfigurationError("Subsequence weights underflowed; increase delta.")
        return {word: value / total for word, value in subsequence_counts.items()}
    if descriptor == "NAC":
        return _frequencies(text, 1)
    if descriptor == "ANF":
        counts: Counter[str] = Counter()
        output = {}
        for index, base in enumerate(text):
            counts[base] += 1
            output[f"position_{index}"] = counts[base] / (index + 1)
        return output
    if descriptor == "ENAC":
        window = positive(parameters.get("window", 5), "window")
        if window > length:
            raise ConfigurationError("ENAC window exceeds sequence length.")
        return {
            f"window_{start}:{base}": text[start : start + window].count(base) / window
            for start in range(length - window + 1)
            for base in "ACGT"
        }
    if descriptor in {"binary", "PS2", "PS3", "PS4"}:
        k = 1 if descriptor == "binary" else int(descriptor[-1])
        if length < k:
            raise ConfigurationError("Sequence is too short for positional encoding.")
        return {
            f"position_{index}:{word}": float(text[index : index + k] == word)
            for index in range(length - k + 1)
            for word in _words(k)
        }
    if descriptor == "NCP":
        return {
            f"position_{index}:{name}": float(base in group)
            for index, base in enumerate(text)
            for name, group in (("purine", "AG"), ("amino", "AC"), ("weak", "AT"))
        }
    if descriptor == "EIIP":
        return {f"position_{index}": _EIIP[base] for index, base in enumerate(text)}
    if descriptor == "PseEIIP":
        return {
            word: value * sum(_EIIP[base] for base in word)
            for word, value in _frequencies(text, 3).items()
        }
    if descriptor in {"CKSNAP", "ASDC"}:
        gap = parameters.get("gap", 3)
        if descriptor == "CKSNAP" and (type(gap) is not int or gap < 0 or gap >= length - 1):
            raise ConfigurationError("CKSNAP gap must be in [0, sequence_length - 2].")
        result = {}
        distances = range(1, length) if descriptor == "ASDC" else range(1, gap + 2)
        totals: Counter[str] = Counter()
        for distance in distances:
            counts = Counter(text[i] + text[i + distance] for i in range(length - distance))
            if descriptor == "ASDC":
                totals.update(counts)
            else:
                for word in _words(2):
                    result[f"gap_{distance - 1}:{word}"] = counts[word] / (
                        length - distance if normalized else 1
                    )
        if descriptor == "ASDC":
            if length < 2:
                raise ConfigurationError("ASDC requires at least two bases.")
            return {word: totals[word] / (length * (length - 1) / 2) for word in _words(2)}
        return result
    if descriptor == "DBE":
        if length < 2:
            raise ConfigurationError("DBE requires at least two bases.")
        return {
            f"position_{index}:bit_{bit}": float(
                (4 * "ACGT".index(text[index]) + "ACGT".index(text[index + 1])) >> (3 - bit) & 1
            )
            for index in range(length - 1)
            for bit in range(4)
        }
    if descriptor == "LPDF":
        if length < 2:
            raise ConfigurationError("LPDF requires at least two bases.")
        counts = Counter()
        output = {}
        for index in range(length - 1):
            word = text[index : index + 2]
            counts[word] += 1
            output[f"position_{index}"] = counts[word] / (index + 1)
        return output
    if descriptor.startswith("Z_curve"):
        dimension = int(descriptor.split("_")[-1].removesuffix("bit"))
        k = 1 if dimension == 9 else (2 if dimension in {12, 36} else 3)
        phased = dimension in {9, 36, 144}
        if length < k:
            raise ConfigurationError("Sequence too short for Z-curve encoding.")
        result = {}
        for phase in range(3 if phased else 1):
            indices = range(phase, length - k + 1, 3 if phased else 1)
            counts = Counter(text[index : index + k] for index in indices)
            for prefix in _words(k - 1):
                a, c, g, t = (counts[prefix + base] / (length - k + 1) for base in "ACGT")
                for axis, value in zip(
                    "xyz", (a + g - c - t, a + c - g - t, a + t - c - g), strict=True
                ):
                    result[f"phase_{phase}:{prefix}:{axis}"] = value
        return result
    if descriptor == "MMI":
        # Ordered positional interaction-information contributions, aggregated by
        # nucleotide multiset (10 pair and 20 triplet cells). Their sums recover
        # pairwise MI and three-way co-information, respectively.
        if length < 3:
            raise ConfigurationError("MMI requires at least three bases.")
        result = {}
        for k in (2, 3):
            joint = _frequencies(text, k)
            marginals = {}
            for size in range(1, k + 1):
                for positions in itertools.combinations(range(k), size):
                    marginal: dict[str, float] = {}
                    for word, probability in joint.items():
                        projected = "".join(word[pos] for pos in positions)
                        marginal[projected] = marginal.get(projected, 0) + probability
                    marginals[positions] = marginal
            for combination in itertools.combinations_with_replacement("ACGT", k):
                result["".join(combination)] = 0.0
            for word, probability in joint.items():
                if probability:
                    contribution = sum(
                        (-1) ** len(positions)
                        * math.log2(marginal["".join(word[pos] for pos in positions)])
                        for positions, marginal in marginals.items()
                    )
                    result["".join(sorted(word))] += probability * contribution
        return result
    raise ConfigurationError(f"Unsupported encoding {descriptor!r}.")


def encode_dna(
    sequences: Mapping[str, Any],
    descriptor: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    property_tables: Mapping[str, Mapping[str, float]] | None = None,
    labels: Mapping[str, int] | None = None,
    training_ids: Iterable[str] | None = None,
    max_output_cells: int = 5_000_000,
    max_input_bases: int = 1_000_000,
    verbose: bool = True,
) -> EncodingResult:
    """Compute native DNA features. Positional encodings require equal lengths.

    PSTNP/KNN use only explicit training labels and exclude the current training
    sample from its own reference set. PSTNPds folds A/T and C/G to A/C.
    Formulas and conventions are versioned; no bitwise equivalence with other
    toolkits is implied. Undefined short/constant inputs raise ConfigurationError.
    """
    np = require("numpy", "bio-analysis")
    if descriptor not in DNA_ENCODINGS or not sequences:
        raise ConfigurationError("Provide sequences and a descriptor from DNA_ENCODINGS.")
    if any(not isinstance(key, str) or not key for key in sequences):
        raise ConfigurationError("Sequence IDs must be non-empty strings.")
    texts = {key: sequence_text(value) for key, value in sequences.items()}
    if any(set(value) - set("ACGT") for value in texts.values()):
        raise ConfigurationError("Native encodings require canonical A/C/G/T sequences.")
    options = dict(parameters or {})
    allowed: set[str] = set()
    if descriptor in {"Kmer", "RCKmer", "Mismatch", "Subsequence", "CKSNAP"}:
        allowed.add("normalized")
    if descriptor in {"Kmer", "RCKmer", "Mismatch", "Subsequence", "PseKNC"}:
        allowed.add("k")
    allowed.update(
        {
            "Mismatch": {"mismatches"},
            "Subsequence": {"delta"},
            "CKSNAP": {"gap"},
            "ENAC": {"window"},
            "KNN": {"neighbors"},
        }.get(descriptor, set())
    )
    if descriptor in {"DAC", "DCC", "DACC", "TAC", "TCC", "TACC", "Moran", "Geary", "NMBroto"}:
        allowed.add("lag")
    if "Pse" in descriptor and descriptor != "PseEIIP":
        allowed.update({"lag", "weight"})
    if unknown := options.keys() - allowed:
        raise ConfigurationError(f"Unknown {descriptor} parameters: {sorted(unknown)}.")
    positive(max_input_bases, "max_input_bases")
    if sum(map(len, texts.values())) > max_input_bases:
        raise ConfigurationError("Encoding inputs exceed max_input_bases.")
    # Bound positional allocations before constructing feature dictionaries.
    width = max(map(len, texts.values()))
    multiplier = {
        "binary": 4,
        "PS2": 16,
        "PS3": 64,
        "PS4": 256,
        "NCP": 3,
        "DBE": 4,
        "ENAC": 4,
        "DPCP": len(property_tables or {}),
        "TPCP": len(property_tables or {}),
        "ANF": 1,
        "EIIP": 1,
        "LPDF": 1,
        "PSTNPss": 1,
        "PSTNPds": 1,
    }.get(descriptor, 0)
    if width * multiplier * len(texts) > max_output_cells:
        raise ConfigurationError("Positional encoding exceeds max_output_cells.")
    if descriptor == "ASDC" and sum(len(t) ** 2 for t in texts.values()) > 10_000_000:
        raise ConfigurationError("ASDC pair count exceeds the computation limit.")
    positive(max_output_cells, "max_output_cells")
    training = set(training_ids or ())
    if not training <= texts.keys():
        raise ConfigurationError("Unknown training IDs.")
    dependent = descriptor in {"PSTNPss", "PSTNPds", "KNN"}
    if dependent:
        if not training or labels is None or not training <= labels.keys():
            raise ConfigurationError("Explicit training_ids and labels are required.")
        if any(type(labels[key]) is not int or labels[key] not in {0, 1} for key in training):
            raise ConfigurationError("Only binary training labels 0 and 1 are supported.")
        if len({len(text) for text in texts.values()}) != 1:
            raise ConfigurationError("Label-dependent positional encodings require equal lengths.")
    rows = []
    names: tuple[str, ...] | None = None
    property_family = descriptor in {
        "DPCP",
        "TPCP",
        "DAC",
        "DCC",
        "DACC",
        "TAC",
        "TCC",
        "TACC",
        "NMBroto",
        "Moran",
        "Geary",
    } or ("Pse" in descriptor and descriptor != "PseEIIP")
    for key, text in texts.items():
        if dependent:
            assert labels is not None
            references = [other for other in texts if other in training and other != key]
            if descriptor == "KNN":
                neighbors = positive(options.get("neighbors", 3), "neighbors")
                if len(references) < neighbors:
                    raise ConfigurationError("Not enough training references for KNN.")
                nearest = sorted(
                    references,
                    key=lambda other: (
                        sum(a != b for a, b in zip(text, texts[other], strict=True)),
                        other,
                    ),
                )[:neighbors]
                features = {
                    "positive_fraction": sum(labels[other] for other in nearest) / neighbors
                }
            else:
                if len(text) < 3:
                    raise ConfigurationError("PSTNP requires at least three bases.")
                groups = [
                    [other for other in references if labels[other] == label] for label in (0, 1)
                ]
                if any(not group for group in groups):
                    raise ConfigurationError("Each leave-one-out training class must be non-empty.")

                def fold(value: str) -> str:
                    return (
                        value.translate(str.maketrans("TG", "AC"))
                        if descriptor == "PSTNPds"
                        else value
                    )

                features = {}
                for index in range(len(text) - 2):
                    word = fold(text[index : index + 3])
                    rates = [
                        sum(fold(texts[other][index : index + 3]) == word for other in group)
                        / len(group)
                        for group in groups
                    ]
                    features[f"position_{index}"] = rates[1] - rates[0]
        elif property_family:
            features = _physicochemical(text, descriptor, options, property_tables or {})
        else:
            features = _composition(text, descriptor, options)
        if names is None:
            names = tuple(features)
            if not names or len(names) * len(texts) > max_output_cells:
                raise ConfigurationError("Encoding output is empty or exceeds max_output_cells.")
        elif tuple(features) != names:
            raise ConfigurationError("Positional encoding requires equal sequence lengths.")
        rows.append(tuple(features.values()))
    array = np.asarray(rows, dtype=float)
    if not np.isfinite(array).all():
        raise ConfigurationError("Encoding produced non-finite values.")
    array.setflags(write=False)
    return EncodingResult(tuple(texts), names or (), array, descriptor, options)


_register_public_verbose(globals(), ("encode_dna",))
