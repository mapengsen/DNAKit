"""Sequence encoding, tokenizer, masking and in-silico mutagenesis operations of ``ml_models``."""

from __future__ import annotations

import itertools
import json
import math
import random
from collections.abc import Callable, Mapping, Sequence

from ._impl_common import (
    Unavailable,
    bad_input,
    bad_option,
    float_list,
    get_choice,
    get_float,
    get_int,
    output_format,
    sequences_from,
    to_float,
    write_or_return,
)

BASES = "ACGT"


def _dna(sequence: str, *, name: str = "sequence") -> object:
    from dnakit.core import DNASequence

    return DNASequence(sequence.upper(), alphabet="iupac")


def encode(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Encode sequences with the existing DNAKit one-hot, integer or k-mer encoders."""

    from dnakit.fingerprints import integer_encode, kmer_fingerprint, one_hot_encode

    method = get_choice(options, "method", "one_hot", ("one_hot", "integer", "kmer"))
    pairs = sequences_from(
        value, options, limit=get_int(options, "max_sequences", 1000, low=1, high=10_000)
    )
    limit = get_int(options, "max_output_length", 100_000, low=1)
    out: list[dict[str, object]] = []
    feature_names: list[str] | None = None
    for name, sequence in pairs:
        dna = _dna(sequence)
        if method == "one_hot":
            res = one_hot_encode(
                dna,
                ambiguity_policy=str(options.get("ambiguity_policy", "error")),
                max_output_length=limit,
                verbose=False,
            )
            out.append(
                {"id": name, "length": len(res.values), "values": [list(r) for r in res.values]}
            )
            feature_names = list(res.feature_names)
        elif method == "integer":
            res = integer_encode(
                dna,
                ambiguity_policy=str(options.get("ambiguity_policy", "iupac")),
                max_output_length=limit,
                verbose=False,
            )
            out.append({"id": name, "length": len(res.values), "values": list(res.values)})
        else:
            res = kmer_fingerprint(
                dna,
                k=get_int(options, "k", 3, low=1, high=8),
                canonical=bool(options.get("canonical", False)),
                mode=str(options.get("mode", "count")),
                ambiguity_policy=str(options.get("ambiguity_policy", "ignore")),
                verbose=False,
            )
            out.append({"id": name, "length": len(res.values), "values": list(res.values)})
            feature_names = list(res.feature_names)
    result: dict[str, object] = {"method": method, "n_sequences": len(out), "encodings": out}
    if feature_names is not None:
        result["feature_names"] = feature_names
    return result


def export_onehot(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """One-hot encode sequences (A, C, G, T columns) and write them as JSON or long-format TSV."""

    from dnakit.fingerprints import one_hot_encode

    pairs = sequences_from(
        value, options, limit=get_int(options, "max_sequences", 10_000, low=1, high=100_000)
    )
    fmt = output_format(options, "json", ("json", "tsv"))
    policy = str(options.get("ambiguity_policy", "error"))
    encoded = {}
    for name, sequence in pairs:
        if name in encoded:
            raise bad_input(f"duplicate sequence id '{name}'.", code="ML_DUPLICATE_ID")
        res = one_hot_encode(
            _dna(sequence), ambiguity_policy=policy, max_output_length=1_000_000, verbose=False
        )
        encoded[name] = (res.values, res.feature_names)
    if fmt == "json":
        text = (
            json.dumps(
                {
                    "columns": list(next(iter(encoded.values()))[1]),
                    "sequences": {n: [list(r) for r in v[0]] for n, v in encoded.items()},
                }
            )
            + "\n"
        )
    else:
        lines = ["id\tposition\t" + "\t".join(next(iter(encoded.values()))[1])]
        for name, (rows, _) in encoded.items():
            lines.extend(
                f"{name}\t{pos}\t" + "\t".join(f"{x:g}" for x in row)
                for pos, row in enumerate(rows)
            )
        text = "\n".join(lines) + "\n"
    info = write_or_return(text, options)
    return {
        "format": fmt,
        "n_sequences": len(encoded),
        "total_positions": sum(len(v[0]) for v in encoded.values()),
        **info,
    }


def export_kmer(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Compute k-mer count/frequency vectors per sequence and write a matrix (TSV/JSON)."""

    from dnakit.fingerprints import kmer_fingerprint

    pairs = sequences_from(
        value, options, limit=get_int(options, "max_sequences", 10_000, low=1, high=100_000)
    )
    k = get_int(options, "k", 3, low=1, high=8)
    mode = get_choice(options, "mode", "frequency", ("frequency", "count", "binary"))
    fmt = output_format(options, "tsv", ("tsv", "json"))
    vectors = []
    names: tuple[str, ...] = ()
    for name, sequence in pairs:
        res = kmer_fingerprint(
            _dna(sequence),
            k=k,
            mode=mode,
            canonical=bool(options.get("canonical", False)),
            ambiguity_policy=str(options.get("ambiguity_policy", "ignore")),
            verbose=False,
        )
        vectors.append((name, list(res.values)))
        names = res.feature_names
    if fmt == "tsv":
        lines = ["id\t" + "\t".join(names)] + [
            f"{n}\t" + "\t".join(f"{v:.10g}" for v in vec) for n, vec in vectors
        ]
        text = "\n".join(lines) + "\n"
    else:
        text = (
            json.dumps(
                {"k": k, "mode": mode, "feature_names": list(names), "vectors": dict(vectors)}
            )
            + "\n"
        )
    info = write_or_return(text, options)
    return {
        "k": k,
        "mode": mode,
        "format": fmt,
        "dimension": len(names),
        "n_sequences": len(vectors),
        **info,
    }


def tokenizer(options: Mapping[str, object]) -> dict[str, object]:
    """Build and export a k-mer tokenizer vocabulary (special tokens first, then k-mers in ACGT order)."""

    k = get_int(options, "k", 3, low=1, high=8)
    stride = get_int(options, "stride", 1, low=1, high=k)
    specials = options.get("special_tokens", ["[PAD]", "[UNK]", "[CLS]", "[SEP]", "[MASK]"])
    if (
        not isinstance(specials, Sequence)
        or isinstance(specials, str)
        or not all(isinstance(s, str) and s for s in specials)
    ):
        raise bad_option("option 'special_tokens' must be a list of non-empty strings.")
    if len(set(specials)) != len(specials):
        raise bad_option("special_tokens must be unique.")
    words = ["".join(w) for w in itertools.product(BASES, repeat=k)]
    vocab = {token: i for i, token in enumerate([*specials, *words])}
    model = {
        "type": "dnakit.kmer_tokenizer",
        "version": 1,
        "model": {"type": "kmer", "k": k, "stride": stride, "alphabet": BASES, "vocab": vocab},
        "special_tokens": list(specials),
        "unk_token": "[UNK]" if "[UNK]" in specials else None,
    }
    info = write_or_return(json.dumps(model, indent=2) + "\n", options)
    return {
        "k": k,
        "stride": stride,
        "vocab_size": len(vocab),
        "tokenizer": None if info.get("output_path") else model,
        **info,
    }


def mask_tokens(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Mask a fraction of non-overlapping k-base tokens, as in masked-language-model pre-processing."""

    pairs = sequences_from(
        value, options, limit=get_int(options, "max_sequences", 10_000, low=1, high=100_000)
    )
    size = get_int(options, "token_size", 1, low=1, high=32)
    rate = get_float(options, "mask_rate", 0.15, low=0.0, high=1.0)
    token = options.get("mask_token", "[MASK]")
    if not isinstance(token, str) or not token:
        raise bad_option("option 'mask_token' must be a non-empty string.")
    seed = get_int(options, "seed", 42)
    out = []
    for index, (name, sequence) in enumerate(pairs):
        tokens = [sequence[i : i + size] for i in range(0, len(sequence), size)]
        count = min(len(tokens), round(rate * len(tokens))) if rate > 0 else 0
        if rate > 0 and count == 0:
            count = 1
        positions = sorted(random.Random(seed + index).sample(range(len(tokens)), count))
        masked = list(tokens)
        for p in positions:
            masked[p] = token
        out.append(
            {
                "id": name,
                "tokens": tokens,
                "masked_tokens": masked,
                "masked_positions": positions,
                "labels": {str(p): tokens[p] for p in positions},
            }
        )
    return {
        "token_size": size,
        "mask_rate": rate,
        "mask_token": token,
        "seed": seed,
        "sequences": out,
    }


# ---------------------------------------------------------------- in-silico mutagenesis


def get_scorer(options: Mapping[str, object]) -> Callable[[str], float]:
    """Return the sequence scoring callable from ``options['scorer']`` or raise ``Unavailable``."""

    scorer = options.get("scorer", options.get("model"))
    if scorer is None:
        raise Unavailable(
            "a callable scorer (options['scorer']: sequence -> float) or trained model",
            "In-silico scoring needs a model: pass options['scorer'], a Python callable mapping a DNA string to a number.",
        )
    if not callable(scorer):
        raise bad_option("options['scorer'] must be callable (sequence -> float).")

    def score(sequence: str) -> float:
        raw = scorer(sequence)
        if isinstance(raw, bool) or not isinstance(raw, (int, float)) or not math.isfinite(raw):
            raise bad_input(
                "scorer must return a finite number for every sequence.",
                code="ML_BAD_SCORER_OUTPUT",
            )
        return float(raw)

    return score


def ism(
    sequence: str, scorer: Callable[[str], float], options: Mapping[str, object]
) -> dict[str, object]:
    """Score all single-base substitutions of ``sequence`` (at ``positions``, default all)."""

    sequence = sequence.upper()
    if set(sequence) - set(BASES):
        raise bad_input(
            "mutagenesis needs an unambiguous A/C/G/T sequence.", code="ML_AMBIGUOUS_SEQUENCE"
        )
    limit = get_int(options, "max_length", 2000, low=1, high=100_000)
    if len(sequence) > limit:
        raise bad_input(f"sequence is longer than max_length={limit}.", code="ML_INPUT_TOO_LARGE")
    positions = options.get("positions")
    if positions is None:
        selected = list(range(len(sequence)))
    else:
        if (
            not isinstance(positions, Sequence)
            or isinstance(positions, str)
            or not all(
                isinstance(p, int) and not isinstance(p, bool) and 0 <= p < len(sequence)
                for p in positions
            )
        ):
            raise bad_option(
                "option 'positions' must be 0-based integer positions inside the sequence."
            )
        selected = sorted(set(positions))
    reference = scorer(sequence)
    effects = []
    per_position = []
    for p in selected:
        deltas = []
        for alt in BASES:
            if alt == sequence[p]:
                continue
            score = scorer(sequence[:p] + alt + sequence[p + 1 :])
            effects.append(
                {
                    "position": p,
                    "ref": sequence[p],
                    "alt": alt,
                    "score": score,
                    "delta": score - reference,
                }
            )
            deltas.append(score - reference)
        per_position.append(
            {
                "position": p,
                "ref": sequence[p],
                "mean_delta": math.fsum(deltas) / 3,
                "max_abs_delta": max(abs(d) for d in deltas),
            }
        )
    return {"reference_score": reference, "effects": effects, "per_position": per_position}


def one_sequence(value: object, options: Mapping[str, object]) -> str:
    """Accept one DNA sequence (string, FASTA text or single record)."""

    pairs = sequences_from(value, options, limit=1)
    return pairs[0][1]


def mutagenesis(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Saturation mutagenesis effect table."""

    sequence = one_sequence(value, options)
    result = ism(sequence, get_scorer(options), options)
    return {"sequence_length": len(sequence), "n_variants": len(result["effects"]), **result}  # type: ignore[arg-type]


def attribution(sequence: str, options: Mapping[str, object]) -> dict[str, object]:
    """Per-position ISM attribution: reference score minus the mean score of the three substitutions."""

    result = ism(sequence, get_scorer(options), options)
    ref = float(result["reference_score"])  # type: ignore[arg-type]
    per_position = result["per_position"]
    return {
        "method": "in_silico_mutagenesis",
        "definition": "attribution[i] = score(sequence) - mean(score(sequence with base i replaced by each other base))",
        "reference_score": ref,
        "positions": [p["position"] for p in per_position],  # type: ignore[union-attr]
        "attributions": [-float(p["mean_delta"]) for p in per_position],  # type: ignore[union-attr]
    }


def attribute(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """ISM attribution for one sequence."""

    sequence = one_sequence(value, options)
    return {"sequence_length": len(sequence), **attribution(sequence, options)}


def export_interpretability(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Export per-position attributions (precomputed or computed by ISM) as TSV or JSON."""

    fmt = output_format(options, "tsv", ("tsv", "json"))
    if isinstance(value, Mapping) and "attributions" in value:
        sequence = str(value.get("sequence", "")).upper()
        values = float_list(value["attributions"], what="attributions")
        positions = list(range(len(values)))
        if sequence and len(sequence) != len(values):
            raise bad_input(
                "sequence and attributions differ in length.", code="ML_LENGTH_MISMATCH"
            )
        method = str(value.get("method", "user_supplied"))
    else:
        sequence = one_sequence(value, options).upper()
        computed = attribution(sequence, options)
        values = list(computed["attributions"])  # type: ignore[arg-type]
        positions = list(computed["positions"])  # type: ignore[arg-type]
        method = str(computed["method"])
    bases = [sequence[p] if sequence else "" for p in positions]
    if fmt == "tsv":
        text = "position\tbase\tattribution\n" + "".join(
            f"{p}\t{b}\t{v:.10g}\n" for p, b, v in zip(positions, bases, values, strict=False)
        )
    else:
        text = (
            json.dumps(
                {
                    "method": method,
                    "sequence": sequence or None,
                    "positions": positions,
                    "attributions": values,
                }
            )
            + "\n"
        )
    info = write_or_return(text, options)
    return {"method": method, "format": fmt, "n_positions": len(values), **info}


def permutation_importance(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Permutation feature importance of a callable model on tabular data."""

    model = options.get("model", options.get("scorer"))
    if model is None:
        raise Unavailable(
            "a trained model passed as options['model'] (callable or object with predict)",
            "Explaining a model needs the model itself: pass options['model'], a callable mapping a list of feature rows to predictions.",
        )
    predict = model if callable(model) else getattr(model, "predict", None)
    if not callable(predict):
        raise bad_option("options['model'] must be callable or provide predict(rows).")
    if not isinstance(value, Mapping) or "X" not in value or "y" not in value:
        raise bad_input("input must be a mapping with 'X' (rows of features) and 'y' (targets).")
    x = [float_list(row, what="X") for row in value["X"]]  # type: ignore[union-attr]
    if not x or len({len(r) for r in x}) != 1 or len(x) != len(value["y"]):  # type: ignore[arg-type]
        raise bad_input(
            "X must be a non-empty rectangular matrix aligned with y.", code="ML_LENGTH_MISMATCH"
        )
    names = value.get("feature_names") or [f"x{j}" for j in range(len(x[0]))]
    if len(names) != len(x[0]):  # type: ignore[arg-type]
        raise bad_input(
            "feature_names must match the number of columns.", code="ML_LENGTH_MISMATCH"
        )
    task = get_choice(options, "task", "regression", ("regression", "classification"))
    repeats = get_int(options, "n_repeats", 5, low=1, high=100)
    seed = get_int(options, "seed", 42)
    from ._impl_metrics import label_text, regression_metrics

    y_raw = list(value["y"])  # type: ignore[arg-type]

    def score(matrix: list[list[float]]) -> float:
        predicted = list(predict(matrix))  # type: ignore[misc]
        if len(predicted) != len(y_raw):
            raise bad_input(
                "the model returned a different number of predictions than rows.",
                code="ML_LENGTH_MISMATCH",
            )
        if task == "regression":
            return float(
                regression_metrics(
                    [to_float(v, what="y") for v in y_raw],
                    [to_float(v, what="prediction") for v in predicted],
                )["r2"]
            )
        return sum(label_text(a) == label_text(b) for a, b in zip(y_raw, predicted, strict=False)) / len(y_raw)

    baseline = score(x)
    rng = random.Random(seed)
    importances = []
    for j in range(len(x[0])):
        drops = []
        for _ in range(repeats):
            column = [row[j] for row in x]
            rng.shuffle(column)
            shuffled = [[*row[:j], column[i], *row[j + 1:]] for i, row in enumerate(x)]
            drops.append(baseline - score(shuffled))
        mean = math.fsum(drops) / repeats
        std = math.sqrt(math.fsum((d - mean) ** 2 for d in drops) / repeats)
        importances.append({"feature": names[j], "importance_mean": mean, "importance_std": std})  # type: ignore[index]
    return {
        "method": "permutation_importance",
        "metric": "r2" if task == "regression" else "accuracy",
        "baseline_score": baseline,
        "n_repeats": repeats,
        "importances": importances,
    }
