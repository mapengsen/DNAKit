"""Native extensible constraint/objective system and bounded sequence search.

Search exhaustion is not a proof of infeasibility. Every result reports actual
violations, edits and evaluation counts; callbacks can add specialist constraints.
"""

from __future__ import annotations

import itertools
import json
import math
import random
from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import positive, sequence_text
from dnakit.exceptions import ConfigurationError
from dnakit.ops.translation import translate

_BASES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


@dataclass(frozen=True)
class SpecEvaluation:
    """Nonnegative constraint violation and a separate utility to maximize."""

    violation: float = 0.0
    utility: float = 0.0
    locations: tuple[tuple[int, int], ...] = ()
    message: str = ""

    def __post_init__(self) -> None:
        if (
            not math.isfinite(self.violation)
            or self.violation < 0
            or not math.isfinite(self.utility)
        ):
            raise ConfigurationError(
                "Specification evaluations must be finite with nonnegative violation."
            )


@dataclass(frozen=True)
class Specification:
    """Custom Python evaluator(sequence, circular) -> SpecEvaluation."""

    name: str
    evaluator: Callable[[str, bool], SpecEvaluation]
    weight: float = 1.0
    locked_positions: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if (
            not self.name
            or not callable(self.evaluator)
            or not math.isfinite(self.weight)
            or self.weight <= 0
        ):
            raise ConfigurationError("Specification needs a name, evaluator and positive weight.")

    def evaluate(self, sequence: str, circular: bool = False) -> SpecEvaluation:
        result = self.evaluator(sequence, circular)
        if not isinstance(result, SpecEvaluation):
            raise ConfigurationError("Custom specifications must return SpecEvaluation.")
        return result


def _rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def specification(
    kind: str,
    *,
    region: tuple[int, int] | None = None,
    weight: float = 1.0,
    verbose: bool = True,
    **parameters: Any,
) -> Specification:
    """Build a native specification. See API guide for each kind's parameters.

    Region coordinates are 0-based half-open, including circular end > length.
    Circular window scans cover the origin. Whole-sequence constraints operate
    on the declared linearization where order/reading frame is meaningful.
    """
    schema = {
        "gc": ({"minimum", "maximum", "window"}, set()),
        "avoid_pattern": ({"pattern", "both_strands"}, {"pattern"}),
        "pattern_count": ({"pattern", "count", "both_strands"}, {"pattern"}),
        "sequence": ({"pattern"}, {"pattern"}),
        "translation": ({"protein"}, {"protein"}),
        "length": ({"minimum", "maximum"}, set()),
        "avoid_changes": ({"sequence", "positions"}, {"sequence"}),
        "unique_kmers": ({"k", "both_strands"}, set()),
        "avoid_stop": ({"allow_terminal"}, set()),
        "avoid_codons": ({"codons"}, {"codons"}),
        "avoid_rare_codons": ({"usage", "threshold"}, {"usage"}),
        "hairpin": ({"match_length", "window", "minimum_loop"}, set()),
        "heterodimer": ({"match_length", "references"}, {"references"}),
        "avoid_matches": ({"match_length", "references", "both_strands"}, {"references"}),
        "tm": ({"minimum", "maximum", "conditions"}, {"minimum", "maximum"}),
        "cai": ({"usage"}, {"usage"}),
        "codon_usage": ({"usage"}, {"usage"}),
        "harmonize": ({"usage", "relative_adaptiveness"}, {"usage", "relative_adaptiveness"}),
        "regions_compatibility": ({"regions", "match_length"}, {"regions"}),
    }
    if kind not in schema:
        raise ConfigurationError(f"Unknown specification kind {kind!r}.")
    allowed, required = schema[kind]
    if parameters.keys() - allowed or required - parameters.keys():
        raise ConfigurationError(
            f"Invalid {kind} parameters: unknown={sorted(parameters.keys() - allowed)}, "
            f"missing={sorted(required - parameters.keys())}."
        )
    parameters = dict(parameters)
    locked: tuple[int, ...] = ()
    if kind == "avoid_changes":
        original = sequence_text(parameters["sequence"])
        parameters["sequence"] = original
        default_positions = range(*region) if region is not None else range(len(original))
        locked = tuple(parameters.get("positions", default_positions))
        if any(type(index) is not int or not 0 <= index < len(original) for index in locked):
            raise ConfigurationError("Locked positions must lie inside the original sequence.")

    def evaluate(sequence: str, circular: bool) -> SpecEvaluation:
        offset = 0
        text = sequence
        if region is not None:
            start, end = region
            if (
                type(start) is not int
                or type(end) is not int
                or not 0 <= start < len(sequence)
                or not start < end <= start + len(sequence)
            ):
                raise ConfigurationError("Invalid specification region.")
            if not circular and end > len(sequence):
                raise ConfigurationError("Linear specification region exceeds the sequence.")
            text = (sequence + sequence)[start:end]
            offset = start
        locations: list[tuple[int, int]] = []
        violation = utility = 0.0
        if kind == "length":
            lower, upper = parameters.get("minimum", 1), parameters.get("maximum", math.inf)
            if not 0 <= lower <= upper:
                raise ConfigurationError("Invalid length bounds.")
            violation = max(lower - len(text), len(text) - upper, 0)
        elif kind == "gc":
            lower, upper = parameters.get("minimum", 0.4), parameters.get("maximum", 0.6)
            window = positive(parameters.get("window", len(text)), "window")
            if not 0 <= lower <= upper <= 1 or not 1 <= window <= len(text):
                raise ConfigurationError("Invalid GC interval/window.")
            scan = text + text[: window - 1] if circular and region is None else text
            for start in range(len(scan) - window + 1):
                fraction = sum(base in "GC" for base in scan[start : start + window]) / window
                error = max(lower - fraction, fraction - upper, 0)
                if error:
                    locations.append((offset + start, offset + start + window))
                    violation += error
        elif kind in {"avoid_pattern", "pattern_count", "sequence"}:
            pattern = sequence_text(parameters["pattern"])
            if kind == "sequence":
                if len(pattern) != len(text):
                    violation = abs(len(pattern) - len(text))
                for index, (base, code) in enumerate(zip(text, pattern, strict=False)):
                    if base not in _BASES[code]:
                        violation += 1
                        locations.append((offset + index, offset + index + 1))
            else:
                scan = text + text[: len(pattern) - 1] if circular and region is None else text
                patterns: tuple[str, ...] = (pattern,)
                if parameters.get("both_strands", False):
                    complement = pattern.translate(
                        str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
                    )[::-1]
                    patterns = tuple(dict.fromkeys((pattern, complement)))
                for start in range(len(scan) - len(pattern) + 1):
                    if any(
                        all(
                            base in _BASES[code]
                            for base, code in zip(
                                scan[start : start + len(candidate)], candidate, strict=True
                            )
                        )
                        for candidate in patterns
                    ):
                        locations.append((offset + start, offset + start + len(pattern)))
                target = 0 if kind == "avoid_pattern" else parameters.get("count", 1)
                if type(target) is not int or target < 0:
                    raise ConfigurationError("Pattern count must be a nonnegative integer.")
                violation = abs(len(locations) - target)
        elif kind in {
            "translation",
            "avoid_stop",
            "avoid_codons",
            "avoid_rare_codons",
            "cai",
            "codon_usage",
            "harmonize",
        }:
            if len(text) % 3:
                raise ConfigurationError("Coding specifications require a complete frame-0 CDS.")
            peptide = translate(text)
            codons = [text[index : index + 3] for index in range(0, len(text), 3)]
            if kind == "translation":
                target = parameters["protein"]
                violation = abs(len(peptide) - len(target)) + sum(
                    a != b for a, b in zip(peptide, target, strict=False)
                )
            elif kind == "avoid_stop":
                positions = [
                    index
                    for index, aa in enumerate(peptide)
                    if aa == "*"
                    and not (parameters.get("allow_terminal", True) and index == len(peptide) - 1)
                ]
                violation = len(positions)
                locations = [(offset + 3 * index, offset + 3 * index + 3) for index in positions]
            elif kind == "avoid_codons":
                forbidden = set(parameters["codons"])
                locations = [
                    (offset + 3 * index, offset + 3 * index + 3)
                    for index, codon in enumerate(codons)
                    if codon in forbidden
                ]
                violation = len(locations)
            else:
                usage = parameters["usage"]
                all_codons = tuple("".join(word) for word in itertools.product("ACGT", repeat=3))
                if set(usage) != set(all_codons) or any(
                    not math.isfinite(value) or value < 0 for value in usage.values()
                ):
                    raise ConfigurationError("Codon usage needs all 64 nonnegative finite entries.")
                synonyms: dict[str, list[str]] = {}
                for codon in all_codons:
                    synonyms.setdefault(translate(codon), []).append(codon)
                relative = {
                    codon: usage[codon] / max(usage[other] for other in synonyms[translate(codon)])
                    if max(usage[other] for other in synonyms[translate(codon)])
                    else 0
                    for codon in all_codons
                }
                if kind == "avoid_rare_codons":
                    threshold = parameters.get("threshold", 0.1)
                    if not 0 <= threshold <= 1:
                        raise ConfigurationError(
                            "Rare codon relative-frequency threshold must be in [0,1]."
                        )
                    locations = [
                        (offset + 3 * i, offset + 3 * i + 3)
                        for i, codon in enumerate(codons)
                        if relative[codon] < threshold
                    ]
                    violation = len(locations)
                elif kind == "cai":
                    selected = [relative[codon] for codon in codons if translate(codon) != "*"]
                    utility = (
                        math.exp(sum(math.log(value) for value in selected) / len(selected))
                        if selected and all(selected)
                        else 0
                    )
                elif kind == "harmonize":
                    target = tuple(parameters["relative_adaptiveness"])
                    if len(target) != len(codons) or any(
                        not math.isfinite(v) or not 0 <= v <= 1 for v in target
                    ):
                        raise ConfigurationError(
                            "Harmonization target needs one relative adaptiveness per codon."
                        )
                    utility = -sum(
                        (relative[codon] - value) ** 2
                        for codon, value in zip(codons, target, strict=True)
                    ) / len(codons)
                else:
                    counts = Counter(codons)
                    for group in synonyms.values():
                        observed = sum(counts[codon] for codon in group)
                        expected = sum(usage[codon] for codon in group)
                        if observed and expected:
                            utility -= sum(
                                abs(counts[codon] / observed - usage[codon] / expected)
                                for codon in group
                            )
        elif kind == "avoid_changes":
            original = parameters["sequence"]
            if len(original) != len(sequence):
                violation = abs(len(original) - len(sequence))
            locations = [
                (index, index + 1)
                for index in locked
                if index >= len(sequence) or sequence[index] != original[index]
            ]
            violation += len(locations)
        elif kind == "unique_kmers":
            k = positive(parameters.get("k", 12), "k")
            if k > len(text):
                raise ConfigurationError("k exceeds sequence length.")
            scan = text + text[: k - 1] if circular and region is None else text
            seen = set()
            for index in range(len(scan) - k + 1):
                word = scan[index : index + k]
                if parameters.get("both_strands", True):
                    word = min(word, _rc(word))
                if word in seen:
                    locations.append((offset + index, offset + index + k))
                seen.add(word)
            violation = len(locations)
        elif kind in {"hairpin", "heterodimer", "avoid_matches"}:
            k = positive(parameters.get("match_length", 12), "match_length")
            if k > len(text):
                raise ConfigurationError("match_length exceeds sequence length.")
            references = tuple(sequence_text(value) for value in parameters.get("references", ()))
            scan = text + text[: k - 1] if circular and region is None else text
            if kind == "hairpin":
                window = positive(parameters.get("window", 200), "window")
                minimum_loop = parameters.get("minimum_loop", 3)
                if type(minimum_loop) is not int or minimum_loop < 0:
                    raise ConfigurationError("minimum_loop must be a nonnegative integer.")
                for index in range(len(scan) - k + 1):
                    complement = _rc(scan[index : index + k])
                    following = scan.find(
                        complement, index + k + minimum_loop, min(len(scan), index + window)
                    )
                    if following >= 0:
                        locations.append((offset + index, offset + following + k))
            else:
                if not references:
                    raise ConfigurationError("Provide explicit reference sequences.")
                words = {
                    ref[index : index + k]
                    for ref in references
                    for index in range(len(ref) - k + 1)
                }
                if kind == "heterodimer":
                    words = {_rc(word) for word in words}
                elif parameters.get("both_strands", True):
                    words |= {_rc(word) for word in words}
                for index in range(len(scan) - k + 1):
                    if scan[index : index + k] in words:
                        locations.append((offset + index, offset + index + k))
            violation = len(locations)
        elif kind == "regions_compatibility":
            k = positive(parameters.get("match_length", 12), "match_length")
            spans = tuple(parameters["regions"])
            if len(spans) < 2 or any(
                len(span) != 2 or not 0 <= span[0] < span[1] <= len(sequence) for span in spans
            ):
                raise ConfigurationError("Provide at least two valid linear regions.")
            for left, right in itertools.combinations(spans, 2):
                a, b = sequence[left[0] : left[1]], sequence[right[0] : right[1]]
                if min(len(a), len(b)) < k:
                    raise ConfigurationError("Compatibility match_length exceeds a region length.")
                complements = {_rc(b[i : i + k]) for i in range(len(b) - k + 1)}
                if any(a[i : i + k] in complements for i in range(len(a) - k + 1)):
                    locations.extend((tuple(left), tuple(right)))
                    violation += 1
        elif kind == "tm":
            from dnakit.core import DNASequence
            from dnakit.thermodynamics import melting_temperature

            tm = melting_temperature(
                DNASequence(text), conditions=parameters.get("conditions")
            ).tm_celsius
            lower, upper = parameters["minimum"], parameters["maximum"]
            if lower > upper:
                raise ConfigurationError("Tm minimum exceeds maximum.")
            violation = max(lower - tm, tm - upper, 0)
        return SpecEvaluation(
            float(violation),
            float(utility if kind in {"cai", "codon_usage", "harmonize"} else -violation),
            tuple(locations),
            kind,
        )

    return Specification(kind, evaluate, weight, locked)


@dataclass(frozen=True)
class OptimizationResult:
    original: str
    sequence: str
    constraints_satisfied: bool
    status: str
    constraints: tuple[SpecEvaluation, ...]
    objectives: tuple[SpecEvaluation, ...]
    edits: tuple[tuple[int, str, str], ...]
    evaluations: int
    iterations: int

    def write_report(self, path: str | Path) -> Path:
        """Write a fresh JSON report, including every failed constraint and edit."""
        from dataclasses import asdict

        target = Path(path)
        with target.open("x", encoding="utf-8") as handle:
            json.dump(asdict(self), handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
        return target


def optimize_sequence(
    sequence: Any,
    *,
    constraints: Iterable[Specification] = (),
    objectives: Iterable[Specification] = (),
    circular: bool = False,
    preserve_translation: bool = False,
    mutable_positions: Iterable[int] | None = None,
    max_iterations: int = 100,
    max_evaluations: int = 100_000,
    seed: int = 0,
    verbose: bool = True,
) -> OptimizationResult:
    """Native greedy coordinate search with hard-constraint priority and weighted utilities.

    Uses full synonymous-codon substitutions when preserve_translation=True.
    Otherwise single-base substitutions are searched. Sequence length is fixed;
    impossible length constraints are reported without pretending to insert/delete.
    A local optimum or exhausted budget is explicitly not proof of infeasibility.
    """
    text = sequence_text(sequence)
    if set(text) - set("ACGT"):
        raise ConfigurationError("Optimization requires canonical DNA.")
    positive(max_iterations, "max_iterations")
    positive(max_evaluations, "max_evaluations")
    rules, goals = tuple(constraints), tuple(objectives)
    positions = set(range(len(text))) if mutable_positions is None else set(mutable_positions)
    if any(type(pos) is not int or not 0 <= pos < len(text) for pos in positions):
        raise ConfigurationError("Mutable positions must be within the sequence.")
    positions -= {position for rule in rules for position in rule.locked_positions}
    if preserve_translation and len(text) % 3:
        raise ConfigurationError("Translation preservation requires a complete frame-0 CDS.")
    synonyms: dict[str, list[str]] = {}
    if preserve_translation:
        for word in itertools.product("ACGT", repeat=3):
            codon = "".join(word)
            synonyms.setdefault(translate(codon), []).append(codon)
    evaluations = 0

    def evaluate(
        candidate: str,
    ) -> tuple[tuple[float, float], tuple[SpecEvaluation, ...], tuple[SpecEvaluation, ...]]:
        nonlocal evaluations
        evaluations += 1
        constraint_results = tuple(rule.evaluate(candidate, circular) for rule in rules)
        goal_results = tuple(goal.evaluate(candidate, circular) for goal in goals)
        rank = (
            -sum(
                rule.weight * result.violation
                for rule, result in zip(rules, constraint_results, strict=True)
            ),
            sum(
                goal.weight * result.utility
                for goal, result in zip(goals, goal_results, strict=True)
            ),
        )
        return rank, constraint_results, goal_results

    current = text
    best_rank, rule_results, goal_results = evaluate(current)
    rng = random.Random(seed)
    iterations = 0
    for iteration in range(1, max_iterations + 1):
        iterations = iteration
        improved = False
        units = list(range(0, len(text), 3)) if preserve_translation else sorted(positions)
        rng.shuffle(units)
        for start in units:
            width = 3 if preserve_translation else 1
            segment = current[start : start + width]
            replacements = synonyms[translate(segment)] if preserve_translation else list("ACGT")
            for replacement in replacements:
                if replacement == segment or any(
                    start + index not in positions and old != new
                    for index, (old, new) in enumerate(zip(segment, replacement, strict=True))
                ):
                    continue
                if evaluations >= max_evaluations:
                    break
                candidate = current[:start] + replacement + current[start + width :]
                rank, candidate_rules, candidate_goals = evaluate(candidate)
                if rank > best_rank:
                    current, best_rank, rule_results, goal_results = (
                        candidate,
                        rank,
                        candidate_rules,
                        candidate_goals,
                    )
                    improved = True
                    break
            if evaluations >= max_evaluations:
                break
        if not improved or evaluations >= max_evaluations or (best_rank[0] == 0 and not goals):
            break
    satisfied = all(result.violation == 0 for result in rule_results)
    status = "constraints_satisfied" if satisfied else "search_exhausted"
    return OptimizationResult(
        text,
        current,
        satisfied,
        status,
        rule_results,
        goal_results,
        tuple(
            (index, old, new)
            for index, (old, new) in enumerate(zip(text, current, strict=True))
            if old != new
        ),
        evaluations,
        iterations,
    )


def specifications_from_record(
    record: Any, *, verbose: bool = True
) -> tuple[tuple[Specification, ...], tuple[Specification, ...]]:
    """Read JSON dnakit_constraint/dnakit_objective qualifiers from GenBank features.

    Example qualifier: {"kind":"gc","minimum":0.4,"maximum":0.6,"window":20}.
    Simple feature locations supply the region automatically. Compound/fuzzy
    feature constraints must be declared explicitly in Python to avoid ambiguity.
    """
    constraints: list[Specification] = []
    objectives: list[Specification] = []
    for feature in record.features:
        qualifiers = feature.qualifiers
        for key, destination in (
            ("dnakit_constraint", constraints),
            ("dnakit_objective", objectives),
        ):
            values = qualifiers.get(key, ())
            if isinstance(values, str):
                values = (values,)
            for value in values:
                try:
                    parameters = json.loads(value)
                    kind = parameters.pop("kind")
                except (ValueError, TypeError, KeyError, AttributeError) as exc:
                    raise ConfigurationError("Invalid JSON specification annotation.") from exc
                location = feature.location
                parts = getattr(location, "parts", ())
                if (
                    len(parts) > 1
                    or location is None
                    or type(location.start).__name__ != "ExactPosition"
                    or type(location.end).__name__ != "ExactPosition"
                    or location.strand == -1
                ):
                    raise ConfigurationError(
                        "Compound, fuzzy or reverse-strand feature constraints "
                        "require explicit Python declarations."
                    )
                region = (int(location.start), int(location.end))
                destination.append(specification(kind, region=region, **parameters))
    return tuple(constraints), tuple(objectives)


def blast_match_specification(
    database: str | Path,
    *,
    minimum_identity: float = 0.8,
    minimum_length: int = 20,
    evalue: float = 1e-5,
    config: Any = None,
    program: str = "blastn",
    weight: float = 1.0,
    verbose: bool = True,
) -> Specification:
    """Reject qualifying local BLAST HSPs against a separately prepared database.

    Uses blastn-short with masking disabled and both strands. A fresh temporary
    workspace is removed after each call; up to 128 sequence results are cached.
    BLAST is heuristic: an absence of HSPs is not a proof of no homology. Searches
    a doubled prefix for circular molecules and reports origin-crossing spans.
    """
    from functools import lru_cache
    from tempfile import TemporaryDirectory

    from dnakit.bioinformatics._common import ExecutionConfig, run

    positive(minimum_length, "minimum_length")
    if not 0 < minimum_identity <= 1 or not math.isfinite(evalue) or evalue <= 0:
        raise ConfigurationError("Identity must be in (0,1] and evalue positive.")
    prefix = Path(database).expanduser().resolve()
    limits = config or ExecutionConfig()

    @lru_cache(maxsize=128)
    def evaluate(sequence: str, circular: bool) -> SpecEvaluation:
        with TemporaryDirectory(prefix="dnakit-blast-spec-") as temporary:
            directory = Path(temporary)
            query, output = directory / "query.fa", directory / "hits.tsv"
            text = sequence + sequence[:-1] if circular else sequence
            query.write_text(">candidate\n" + text + "\n", encoding="ascii")
            run(
                program,
                (
                    "-task",
                    "blastn-short",
                    "-dust",
                    "no",
                    "-query",
                    str(query),
                    "-db",
                    str(prefix),
                    "-out",
                    str(output),
                    "-outfmt",
                    "6 qstart qend pident length evalue",
                    "-evalue",
                    str(evalue),
                    "-num_threads",
                    str(limits.threads),
                ),
                directory,
                limits,
                outputs=(output,),
            )
            locations = set()
            for line in output.read_text(encoding="ascii").splitlines():
                start, end, identity, length, _ = line.split()
                left, right = sorted((int(start), int(end)))
                left -= 1
                if (
                    left < len(sequence)
                    and float(identity) >= 100 * minimum_identity
                    and int(length) >= minimum_length
                ):
                    locations.add((left, min(right, left + len(sequence))))
            return SpecEvaluation(
                float(len(locations)),
                -float(len(locations)),
                tuple(sorted(locations)),
                "avoid_blast_matches",
            )

    return Specification("avoid_blast_matches", evaluate, weight)


_register_public_verbose(
    globals(),
    (
        "specification",
        "optimize_sequence",
        "specifications_from_record",
        "blast_match_specification",
    ),
)
