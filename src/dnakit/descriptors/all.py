"""Unified calculation of the fixed 240-field DNA descriptor schema."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from itertools import product
from typing import cast

from dnakit.core import (
    Citation,
    DNASequence,
    ExecutionMode,
    Gap,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
    ReferenceInfo,
    Topology,
)
from dnakit.core._json import FrozenDict, JSONScalar
from dnakit.descriptors._dinucleotide import (
    DINUCLEOTIDE_PROPERTY_SPECS,
    DINUCLEOTIDE_TABLE_SCHEMA_VERSION,
    DinucleotidePropertyTable,
)
from dnakit.descriptors._shared import (
    DescriptorAmbiguityPolicy,
    SequenceInput,
    canonical_runs,
    coerce_ambiguity_policy,
    fragments,
    iter_kmers,
    reject_ambiguity,
    sequence_and_id,
)
from dnakit.descriptors.codon import codon_statistics
from dnakit.descriptors.complexity import exact_repeat_fraction
from dnakit.descriptors.results import AllDescriptorsResult
from dnakit.descriptors.schema import (
    DESCRIPTOR_NAMES_V1,
    DESCRIPTOR_SCHEMA_VERSION,
)
from dnakit.exceptions import ConfigurationError, DNAKitError
from dnakit.thermodynamics import (
    ThermodynamicConditions,
    extinction_coefficient_260nm,
    melting_temperature,
    molecular_weight,
    nearest_neighbor,
    stacking_interactions,
)

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_STANDARD_START = "ATG"
_STANDARD_STOPS = frozenset({"TAA", "TAG", "TGA"})
_MAX_LZ76_LENGTH = 10_000


@dataclass
class _ValueBuilder:
    values: dict[str, JSONScalar]
    reasons: dict[str, str]

    def __init__(self) -> None:
        self.values = {}
        self.reasons = {}

    def add(self, name: str, value: JSONScalar, reason: str | None = None) -> None:
        if name in self.values:
            raise AssertionError(f"Descriptor {name!r} was assigned more than once.")
        if value is None:
            if reason is None or not reason.strip():
                raise AssertionError(f"Unavailable descriptor {name!r} requires a reason.")
            self.reasons[name] = reason
        elif reason is not None:
            raise AssertionError(f"Available descriptor {name!r} cannot have a reason.")
        self.values[name] = value

    def ordered(self) -> tuple[dict[str, JSONScalar], dict[str, str]]:
        missing = set(DESCRIPTOR_NAMES_V1) - set(self.values)
        extra = set(self.values) - set(DESCRIPTOR_NAMES_V1)
        if missing or extra:
            raise AssertionError(
                f"All-descriptor assignments differ from schema: missing={missing}, extra={extra}."
            )
        return (
            {name: self.values[name] for name in DESCRIPTOR_NAMES_V1},
            {name: self.reasons[name] for name in DESCRIPTOR_NAMES_V1 if name in self.reasons},
        )


def _ratio(numerator: int | float, denominator: int | float) -> float | None:
    return numerator / denominator if denominator else None


def _entropy(counts: Counter[str], total: int) -> float | None:
    if total == 0:
        return None
    return -math.fsum(
        (count / total) * math.log2(count / total) for count in counts.values() if count
    )


def _canonical_run_locations(sequence: DNASequence) -> tuple[tuple[str, int], ...]:
    located: list[tuple[str, int]] = []
    symbol_offset = 0
    for part in sequence.parts:
        if not isinstance(part, str):
            continue
        start = 0
        for index, symbol in enumerate(part):
            if symbol not in "ACGT":
                if start < index:
                    located.append((part[start:index], symbol_offset + start))
                start = index + 1
        if start < len(part):
            located.append((part[start:], symbol_offset + start))
        symbol_offset += len(part)
    return tuple(located)


def _cumulative_range(text: str, positive: str, negative: str) -> tuple[int, int, int]:
    score = 0
    minimum = 0
    maximum = 0
    for symbol in text:
        if symbol == positive:
            score += 1
        elif symbol == negative:
            score -= 1
        minimum = min(minimum, score)
        maximum = max(maximum, score)
    return maximum, minimum, maximum - minimum


def _lz76_complexity(text: str) -> int:
    """Return exhaustive-history LZ76 phrase count using the Kaspar-Schuster scan."""

    length = len(text)
    if length <= 1:
        return length
    complexity = 1
    prefix_index = 0
    phrase_start = 1
    match_length = 1
    longest_match = 1
    while True:
        if text[prefix_index + match_length - 1] == text[phrase_start + match_length - 1]:
            match_length += 1
            if phrase_start + match_length > length:
                complexity += 1
                break
            continue
        longest_match = max(longest_match, match_length)
        prefix_index += 1
        if prefix_index == phrase_start:
            complexity += 1
            phrase_start += longest_match
            if phrase_start + 1 > length:
                break
            prefix_index = 0
            match_length = 1
            longest_match = 1
        else:
            match_length = 1
    return complexity


def _union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    ordered = sorted(intervals)
    total = 0
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    return total + current_end - current_start


def _six_frame_orfs(
    runs: tuple[tuple[str, int], ...],
) -> tuple[int, int, int, int, list[tuple[int, int]]]:
    forward_count = 0
    reverse_count = 0
    longest = 0
    covered: list[tuple[int, int]] = []
    for text, symbol_offset in runs:
        run_length = len(text)
        for reverse in (False, True):
            oriented = text.translate(_COMPLEMENT)[::-1] if reverse else text
            for frame in range(3):
                open_starts: list[int] = []
                coding_end = frame + ((run_length - frame) // 3) * 3
                for position in range(frame, coding_end, 3):
                    codon = oriented[position : position + 3]
                    if codon == _STANDARD_START:
                        open_starts.append(position)
                    if codon not in _STANDARD_STOPS:
                        continue
                    if open_starts:
                        count = len(open_starts)
                        if reverse:
                            reverse_count += count
                        else:
                            forward_count += count
                        end = position + 3
                        earliest = open_starts[0]
                        longest = max(longest, end - earliest)
                        if reverse:
                            covered.append(
                                (
                                    symbol_offset + run_length - end,
                                    symbol_offset + run_length - earliest,
                                )
                            )
                        else:
                            covered.append((symbol_offset + earliest, symbol_offset + end))
                    open_starts.clear()
    return forward_count + reverse_count, forward_count, reverse_count, longest, covered


def _exception_reason(error: DNAKitError) -> str:
    return f"{error.code}: {error.message}"


def _provenance(table: DinucleotidePropertyTable | None) -> Provenance:
    return Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.REIMPLEMENTATION,
            execution_mode=ExecutionMode.INTERNAL,
            origin_class=OriginClass.PUBLISHED_ALGORITHM,
            license_expression="MIT",
            citations=(
                Citation(
                    "shannon1948",
                    title="A Mathematical Theory of Communication",
                    doi="10.1002/j.1538-7305.1948.tb01338.x",
                ),
                Citation(
                    "lempel-ziv1976",
                    title="On the Complexity of Finite Sequences",
                    doi="10.1109/TIT.1976.1055501",
                ),
                Citation(
                    "santalucia1998",
                    title="A unified view of DNA nearest-neighbor thermodynamics",
                    doi="10.1073/pnas.95.4.1460",
                ),
                Citation(
                    "warshaw-tinoco1966",
                    title="Optical properties of sixteen dinucleoside phosphates",
                    doi="10.1016/0022-2836(66)90115-X",
                ),
            ),
        ),
        reference=(
            None
            if table is None
            else ReferenceInfo(
                table.name,
                version=table.version,
                checksum=table.sha256,
                filters={
                    "source_declared_by_user": table.source,
                    "schema_version": DINUCLEOTIDE_TABLE_SCHEMA_VERSION,
                },
            )
        ),
    )


def all_descriptors(
    value: SequenceInput,
    *,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.IGNORE,
    conditions: ThermodynamicConditions | None = None,
    dinucleotide_property_table: DinucleotidePropertyTable | None = None,
) -> AllDescriptorsResult:
    """Calculate all 240 fields in ``descriptor_schema_v1``.

    Exact words never cross an ambiguity symbol or explicit gap. A field outside
    its mathematical or model domain is returned as ``None`` and receives a
    stable entry in ``unavailable_reasons``; it is never silently replaced by 0.
    The final 60 legacy ``diprodb_*`` fields require an explicitly loaded
    user-supplied table because DNAKit bundles no third-party coefficients.
    """

    sequence, sequence_id = sequence_and_id(value)
    policy = coerce_ambiguity_policy(ambiguity_policy)
    reject_ambiguity(sequence, policy)
    resolved_conditions = ThermodynamicConditions() if conditions is None else conditions
    if not isinstance(resolved_conditions, ThermodynamicConditions):
        raise ConfigurationError(
            "conditions must be ThermodynamicConditions or None.",
            code="INVALID_ALL_DESCRIPTOR_CONDITIONS",
        )
    if dinucleotide_property_table is not None and not isinstance(
        dinucleotide_property_table, DinucleotidePropertyTable
    ):
        raise ConfigurationError(
            "dinucleotide_property_table must be DinucleotidePropertyTable or None.",
            code="INVALID_DINUCLEOTIDE_PROPERTY_TABLE",
        )

    builder = _ValueBuilder()
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    run_texts = tuple(
        run
        for fragment in fragments(sequence, cross_gaps=False)
        for run in canonical_runs(fragment)
    )
    located_runs = _canonical_run_locations(sequence)
    canonical_text = "".join(run_texts)
    counts = Counter(canonical_text)
    canonical_count = len(canonical_text)
    symbol_length = sequence.symbol_length
    unknown_gap_count = sum(gap.length is None for gap in gaps)
    known_gap_nt = sum(gap.length or 0 for gap in gaps)

    builder.add("symbol_length", symbol_length)
    builder.add(
        "coordinate_span",
        sequence.coordinate_span,
        "undefined because at least one explicit gap has unknown length"
        if sequence.coordinate_span is None
        else None,
    )
    builder.add("canonical_base_count", canonical_count)
    builder.add("ambiguity_symbol_count", sequence.ambiguity_count)
    builder.add("gap_object_count", len(gaps))
    builder.add("known_gap_nt", known_gap_nt)
    builder.add("unknown_gap_count", unknown_gap_count)
    builder.add(
        "canonical_symbol_fraction",
        _ratio(canonical_count, symbol_length),
        "requires at least one nucleotide symbol" if symbol_length == 0 else None,
    )
    builder.add(
        "ambiguity_symbol_fraction",
        _ratio(sequence.ambiguity_count, symbol_length),
        "requires at least one nucleotide symbol" if symbol_length == 0 else None,
    )
    coordinate_span = sequence.coordinate_span
    builder.add(
        "known_gap_fraction",
        _ratio(known_gap_nt, coordinate_span or 0) if coordinate_span is not None else None,
        (
            "requires a positive, fully known coordinate span"
            if coordinate_span in {None, 0}
            else None
        ),
    )
    builder.add("canonical_run_count", len(run_texts))
    builder.add("longest_canonical_run_nt", max(map(len, run_texts), default=0))

    group_counts = {
        "purine": counts["A"] + counts["G"],
        "pyrimidine": counts["C"] + counts["T"],
        "amino": counts["A"] + counts["C"],
        "keto": counts["G"] + counts["T"],
        "weak": counts["A"] + counts["T"],
        "strong": counts["C"] + counts["G"],
    }
    empty_base_reason = "requires at least one canonical A/C/G/T base"
    for group in ("purine", "pyrimidine", "amino", "keto", "weak", "strong"):
        builder.add(f"{group}_count", group_counts[group])
        builder.add(
            f"{group}_fraction",
            _ratio(group_counts[group], canonical_count),
            empty_base_reason if canonical_count == 0 else None,
        )
    for name, left, right in (
        ("purine_pyrimidine_skew", "purine", "pyrimidine"),
        ("amino_keto_skew", "amino", "keto"),
        ("weak_strong_skew", "weak", "strong"),
    ):
        denominator = group_counts[left] + group_counts[right]
        builder.add(
            name,
            _ratio(group_counts[left] - group_counts[right], denominator),
            empty_base_reason if denominator == 0 else None,
        )
    builder.add(
        "gc_at_ratio",
        _ratio(group_counts["strong"], group_counts["weak"]),
        "requires at least one A or T base" if group_counts["weak"] == 0 else None,
    )

    kmer_counts: dict[int, Counter[str]] = {}
    kmer_totals: dict[int, int] = {}
    for k in range(1, 7):
        observed = Counter(iter_kmers(sequence, k=k, overlapping=True, cross_gaps=False))
        kmer_counts[k] = observed
        kmer_totals[k] = sum(observed.values())
        if k <= 3:
            for symbols in product("ACGT", repeat=k):
                word = "".join(symbols)
                builder.add(
                    f"k{k}_{word}_frequency",
                    _ratio(observed[word], kmer_totals[k]),
                    f"requires at least one valid canonical {k}-mer position"
                    if kmer_totals[k] == 0
                    else None,
                )

    gc_denominator = counts["G"] + counts["C"]
    at_denominator = counts["A"] + counts["T"]
    builder.add(
        "gc_skew",
        _ratio(counts["G"] - counts["C"], gc_denominator),
        "requires at least one G or C base" if gc_denominator == 0 else None,
    )
    builder.add(
        "at_skew",
        _ratio(counts["A"] - counts["T"], at_denominator),
        "requires at least one A or T base" if at_denominator == 0 else None,
    )
    dinucleotide_total = kmer_totals[2]
    cpg_count = kmer_counts[2]["CG"]
    gpc_count = kmer_counts[2]["GC"]
    builder.add("cpg_count", cpg_count)
    builder.add(
        "cpg_density",
        _ratio(cpg_count, dinucleotide_total),
        "requires at least one valid canonical dinucleotide position"
        if dinucleotide_total == 0
        else None,
    )
    builder.add(
        "cpg_observed_expected",
        _ratio(cpg_count * canonical_count, counts["C"] * counts["G"]),
        "requires at least one C and one G base" if counts["C"] == 0 or counts["G"] == 0 else None,
    )
    builder.add("gpc_count", gpc_count)
    builder.add(
        "gpc_density",
        _ratio(gpc_count, dinucleotide_total),
        "requires at least one valid canonical dinucleotide position"
        if dinucleotide_total == 0
        else None,
    )
    builder.add(
        "cpg_gpc_ratio",
        _ratio(cpg_count, gpc_count),
        "requires at least one GC dinucleotide" if gpc_count == 0 else None,
    )
    if canonical_count:
        gc_cumulative = _cumulative_range(canonical_text, "G", "C")
        at_cumulative = _cumulative_range(canonical_text, "A", "T")
        for name, value_item in zip(
            ("cumulative_gc_skew_max", "cumulative_gc_skew_min", "cumulative_gc_skew_range"),
            gc_cumulative,
            strict=True,
        ):
            builder.add(name, value_item)
        for name, value_item in zip(
            ("cumulative_at_skew_max", "cumulative_at_skew_min", "cumulative_at_skew_range"),
            at_cumulative,
            strict=True,
        ):
            builder.add(name, value_item)
    else:
        for name in (
            "cumulative_gc_skew_max",
            "cumulative_gc_skew_min",
            "cumulative_gc_skew_range",
            "cumulative_at_skew_max",
            "cumulative_at_skew_min",
            "cumulative_at_skew_range",
        ):
            builder.add(name, None, empty_base_reason)
    if dinucleotide_total:
        reverse_complement_frequency = 0.0
        for word in ("".join(item) for item in product("ACGT", repeat=2)):
            reverse_complement = word.translate(_COMPLEMENT)[::-1]
            reverse_complement_frequency += abs(
                kmer_counts[2][word] / dinucleotide_total
                - kmer_counts[2][reverse_complement] / dinucleotide_total
            )
        builder.add("dinucleotide_rc_total_variation", 0.5 * reverse_complement_frequency)
    else:
        builder.add(
            "dinucleotide_rc_total_variation",
            None,
            "requires at least one valid canonical dinucleotide position",
        )
    builder.add(
        "mono_chargaff_l1_distance",
        (
            abs(counts["A"] / canonical_count - counts["T"] / canonical_count)
            + abs(counts["C"] / canonical_count - counts["G"] / canonical_count)
            if canonical_count
            else None
        ),
        empty_base_reason if canonical_count == 0 else None,
    )

    entropies: dict[int, float | None] = {}
    for k in range(1, 4):
        entropy = _entropy(kmer_counts[k], kmer_totals[k])
        entropies[k] = entropy
        builder.add(
            f"shannon_entropy_k{k}_bits",
            entropy,
            f"requires at least one valid canonical {k}-mer position" if entropy is None else None,
        )
    for k in range(1, 4):
        entropy = entropies[k]
        builder.add(
            f"normalized_entropy_k{k}",
            entropy / (2 * k) if entropy is not None else None,
            f"requires at least one valid canonical {k}-mer position" if entropy is None else None,
        )
    linguistic: dict[int, float | None] = {}
    for k in range(1, 7):
        possible = min(4**k, kmer_totals[k])
        linguistic[k] = len(kmer_counts[k]) / possible if possible else None
        if k >= 2:
            builder.add(
                f"linguistic_complexity_k{k}",
                linguistic[k],
                f"requires at least one valid canonical {k}-mer position"
                if possible == 0
                else None,
            )
    defined_linguistic = tuple(item for item in linguistic.values() if item is not None)
    builder.add(
        "linguistic_complexity_product_k1_k6",
        math.prod(defined_linguistic) if defined_linguistic else None,
        "requires at least one valid canonical k-mer position" if not defined_linguistic else None,
    )

    lz_reason: str | None = None
    if sequence.topology is Topology.CIRCULAR:
        lz_reason = "LZ76 descriptor requires a linear sequence"
    elif gaps:
        lz_reason = "LZ76 descriptor requires an ungapped sequence"
    elif sequence.ambiguity_count:
        lz_reason = "LZ76 descriptor requires only canonical A/C/G/T symbols"
    elif canonical_count == 0:
        lz_reason = empty_base_reason
    elif canonical_count > _MAX_LZ76_LENGTH:
        lz_reason = f"LZ76 descriptor is limited to {_MAX_LZ76_LENGTH} nt"
    if lz_reason is None:
        lz_value = _lz76_complexity(canonical_text)
        builder.add("lz76_complexity", lz_value)
        builder.add(
            "normalized_lz76_complexity",
            lz_value * math.log(canonical_count, 4) / canonical_count
            if canonical_count > 1
            else None,
            "normalized LZ76 requires at least two canonical bases"
            if canonical_count <= 1
            else None,
        )
    else:
        builder.add("lz76_complexity", None, lz_reason)
        builder.add("normalized_lz76_complexity", None, lz_reason)

    longest_by_base = {base: 0 for base in "ACGT"}
    for run in run_texts:
        start = 0
        for index in range(1, len(run) + 1):
            if index == len(run) or run[index] != run[start]:
                longest_by_base[run[start]] = max(longest_by_base[run[start]], index - start)
                start = index
    builder.add("longest_homopolymer_nt", max(longest_by_base.values(), default=0))
    for base in "ACGT":
        builder.add(f"longest_homopolymer_{base.lower()}_nt", longest_by_base[base])
    if canonical_count == 0:
        builder.add("exact_tandem_repeat_coverage_fraction", None, empty_base_reason)
    else:
        try:
            repeat = exact_repeat_fraction(
                value,
                ambiguity_policy=policy,
                cross_gaps=False,
            )
            builder.add("exact_tandem_repeat_coverage_fraction", repeat.repeat_fraction)
        except DNAKitError as error:
            builder.add(
                "exact_tandem_repeat_coverage_fraction",
                None,
                _exception_reason(error),
            )

    codons = codon_statistics(
        value,
        frame=0,
        genetic_code=1,
        ambiguity_policy=policy,
        cross_gaps=False,
    )
    codon_counts = Counter({codon: cast(int, count) for codon, count in codons.counts.items()})
    codon_entropy = _entropy(codon_counts, codons.codon_count)
    builder.add("frame0_codon_count", codons.codon_count)
    builder.add("frame0_unique_codon_count", len(codon_counts))
    builder.add("frame0_start_codon_count", codons.start_count)
    builder.add("frame0_stop_codon_count", codons.stop_count)
    builder.add(
        "frame0_start_codon_fraction",
        codons.start_density,
        "requires at least one valid forward frame-0 codon" if codons.codon_count == 0 else None,
    )
    builder.add(
        "frame0_stop_codon_fraction",
        codons.stop_density,
        "requires at least one valid forward frame-0 codon" if codons.codon_count == 0 else None,
    )
    builder.add(
        "frame0_codon_entropy_bits",
        codon_entropy,
        "requires at least one valid forward frame-0 codon" if codon_entropy is None else None,
    )
    builder.add(
        "frame0_effective_number_of_codons",
        2**codon_entropy if codon_entropy is not None else None,
        "requires at least one valid forward frame-0 codon" if codon_entropy is None else None,
    )
    for position in range(3):
        gc_count = sum(count for codon, count in codon_counts.items() if codon[position] in "GC")
        builder.add(
            f"frame0_gc{position + 1}_fraction",
            _ratio(gc_count, codons.codon_count),
            "requires at least one valid forward frame-0 codon"
            if codons.codon_count == 0
            else None,
        )

    orf_names = (
        "six_frame_complete_orf_count",
        "six_frame_forward_complete_orf_count",
        "six_frame_reverse_complete_orf_count",
        "six_frame_longest_complete_orf_nt",
        "six_frame_complete_orf_coverage_fraction",
    )
    if sequence.topology is Topology.CIRCULAR:
        for name in orf_names:
            builder.add(name, None, "six-frame ORF descriptors require a linear sequence")
    else:
        total_orfs, forward_orfs, reverse_orfs, longest_orf, orf_intervals = _six_frame_orfs(
            located_runs
        )
        builder.add("six_frame_complete_orf_count", total_orfs)
        builder.add("six_frame_forward_complete_orf_count", forward_orfs)
        builder.add("six_frame_reverse_complete_orf_count", reverse_orfs)
        builder.add("six_frame_longest_complete_orf_nt", longest_orf)
        builder.add(
            "six_frame_complete_orf_coverage_fraction",
            _ratio(_union_length(orf_intervals), canonical_count),
            empty_base_reason if canonical_count == 0 else None,
        )

    thermo_sequence: DNASequence | None = None
    thermo_reason: str | None = None
    if sequence.topology is Topology.CIRCULAR:
        thermo_reason = "thermodynamic descriptors require a linear sequence"
    elif gaps:
        thermo_reason = "thermodynamic descriptors require an ungapped sequence"
    elif sequence.ambiguity_count:
        thermo_reason = "thermodynamic descriptors require only canonical A/C/G/T symbols"
    else:
        thermo_sequence = DNASequence(canonical_text, strandedness="single")

    mass_names = (
        "mw_ss_oh_da",
        "mw_ss_5p_phosphate_da",
        "mw_ds_oh_da",
        "mw_ds_5p_phosphate_da",
    )
    if thermo_sequence is None:
        for name in mass_names:
            builder.add(name, None, cast(str, thermo_reason))
    else:
        try:
            mass_values = (
                molecular_weight(thermo_sequence).value_dalton,
                molecular_weight(
                    thermo_sequence,
                    five_prime_phosphorylated=True,
                ).value_dalton,
                molecular_weight(thermo_sequence, strand="double").value_dalton,
                molecular_weight(
                    thermo_sequence,
                    strand="double",
                    five_prime_phosphorylated=True,
                ).value_dalton,
            )
            for name, mass in zip(mass_names, mass_values, strict=True):
                builder.add(name, mass)
        except DNAKitError as error:
            reason = _exception_reason(error)
            for name in mass_names:
                builder.add(name, None, reason)

    epsilon_names = (
        "epsilon260_ss_m_inverse_cm_inverse",
        "nmol_per_a260_1ml_1cm",
        "ug_per_a260_1ml_1cm",
    )
    if thermo_sequence is None:
        for name in epsilon_names:
            builder.add(name, None, cast(str, thermo_reason))
    else:
        try:
            extinction = extinction_coefficient_260nm(thermo_sequence)
            epsilon = extinction.value_m_inverse_cm_inverse
            mw_value = builder.values["mw_ss_oh_da"]
            builder.add("epsilon260_ss_m_inverse_cm_inverse", epsilon)
            builder.add("nmol_per_a260_1ml_1cm", 1_000_000.0 / epsilon)
            builder.add(
                "ug_per_a260_1ml_1cm",
                1000.0 * cast(float, mw_value) / epsilon if mw_value is not None else None,
                "requires an available ssDNA molecular weight" if mw_value is None else None,
            )
        except DNAKitError as error:
            reason = _exception_reason(error)
            for name in epsilon_names:
                builder.add(name, None, reason)

    if thermo_sequence is None:
        builder.add("tm_wallace_c", None, cast(str, thermo_reason))
    else:
        try:
            wallace = melting_temperature(
                thermo_sequence,
                method="wallace",
                conditions=resolved_conditions,
            )
            builder.add("tm_wallace_c", wallace.tm_celsius)
        except DNAKitError as error:
            builder.add("tm_wallace_c", None, _exception_reason(error))

    stacking_names = (
        "stacking_delta_h_kcal_per_mol",
        "stacking_delta_s_cal_per_k_mol_k",
        "stacking_delta_g37_kcal_per_mol",
    )
    if thermo_sequence is None:
        for name in stacking_names:
            builder.add(name, None, cast(str, thermo_reason))
    else:
        try:
            stacking = stacking_interactions(thermo_sequence, temperature_celsius=37.0)
            for name, item in zip(
                stacking_names,
                (
                    stacking.total_delta_h_kcal_per_mol,
                    stacking.total_delta_s_cal_per_k_mol,
                    stacking.total_delta_g_kcal_per_mol,
                ),
                strict=True,
            ):
                builder.add(name, item)
        except DNAKitError as error:
            reason = _exception_reason(error)
            for name in stacking_names:
                builder.add(name, None, reason)

    nearest_neighbor_names = (
        "nn_delta_h_kcal_per_mol",
        "nn_delta_s_cal_per_mol_k",
        "nn_delta_g37_kcal_per_mol",
        "nn_tm_c",
    )
    if thermo_sequence is None:
        for name in nearest_neighbor_names:
            builder.add(name, None, cast(str, thermo_reason))
    else:
        try:
            free_energy_conditions = ThermodynamicConditions(
                temperature_celsius=37.0,
                sodium_molar=resolved_conditions.sodium_molar,
                potassium_molar=resolved_conditions.potassium_molar,
                magnesium_molar=resolved_conditions.magnesium_molar,
                dntp_molar=resolved_conditions.dntp_molar,
                strand_concentration_molar=resolved_conditions.strand_concentration_molar,
                dmso_percent=resolved_conditions.dmso_percent,
                dmso_factor_celsius_per_percent=(
                    resolved_conditions.dmso_factor_celsius_per_percent
                ),
                formamide_molar=resolved_conditions.formamide_molar,
            )
            nearest = nearest_neighbor(thermo_sequence, conditions=free_energy_conditions)
            for name, item in zip(
                nearest_neighbor_names,
                (
                    nearest.delta_h_kcal_per_mol,
                    nearest.delta_s_cal_per_k_mol,
                    nearest.delta_g_kcal_per_mol,
                    nearest.tm_celsius,
                ),
                strict=True,
            ):
                builder.add(name, item)
        except DNAKitError as error:
            reason = _exception_reason(error)
            for name in nearest_neighbor_names:
                builder.add(name, None, reason)

    if thermo_sequence is None:
        builder.add("self_complementary", None, cast(str, thermo_reason))
    elif canonical_count == 0:
        builder.add("self_complementary", None, empty_base_reason)
    else:
        builder.add(
            "self_complementary",
            canonical_text == canonical_text.translate(_COMPLEMENT)[::-1],
        )

    for property_spec in DINUCLEOTIDE_PROPERTY_SPECS:
        names = tuple(
            f"diprodb_{property_spec.key}_{stat}" for stat in ("mean", "sd", "min", "max")
        )
        if dinucleotide_property_table is None:
            for name in names:
                builder.add(
                    name,
                    None,
                    (
                        "requires an explicit user-supplied DinucleotidePropertyTable; "
                        "DNAKit bundles no DiProDB numerical values"
                    ),
                )
            continue
        property_definition = dinucleotide_property_table.property(property_spec.key)
        weighted_values = tuple(
            (float(cast(int | float, property_definition.values[word])), count)
            for word, count in kmer_counts[2].items()
            if count
        )
        if not weighted_values:
            for name in names:
                builder.add(
                    name,
                    None,
                    "requires at least one valid canonical dinucleotide position",
                )
            continue
        property_mean = (
            math.fsum(property_value * count for property_value, count in weighted_values)
            / dinucleotide_total
        )
        statistics = (
            property_mean,
            math.sqrt(
                math.fsum(
                    count * (property_value - property_mean) ** 2
                    for property_value, count in weighted_values
                )
                / dinucleotide_total
            ),
            min(property_value for property_value, _ in weighted_values),
            max(property_value for property_value, _ in weighted_values),
        )
        for name, statistic in zip(names, statistics, strict=True):
            builder.add(name, statistic)

    ordered_values, ordered_reasons = builder.ordered()
    return AllDescriptorsResult(
        schema_version=DESCRIPTOR_SCHEMA_VERSION,
        sequence_id=sequence_id,
        values=FrozenDict(ordered_values),
        unavailable_reasons=FrozenDict(ordered_reasons),
        conditions=FrozenDict(
            {
                "ambiguity_policy": policy.value,
                "cross_gaps": False,
                "kmer_overlapping": True,
                "kmer_lengths": [1, 2, 3],
                "linguistic_k_range": [1, 6],
                "lz76_max_sequence_length": _MAX_LZ76_LENGTH,
                "orf_genetic_code": 1,
                "orf_start_codons": [_STANDARD_START],
                "orf_stop_codons": sorted(_STANDARD_STOPS),
                "orf_complete_only": True,
                "repeat_min_unit_length": 1,
                "repeat_max_unit_length": 20,
                "repeat_min_copies": 2,
                "temperature_celsius": resolved_conditions.temperature_celsius,
                "free_energy_reference_temperature_celsius": 37.0,
                "sodium_molar": resolved_conditions.sodium_molar,
                "potassium_molar": resolved_conditions.potassium_molar,
                "magnesium_molar": resolved_conditions.magnesium_molar,
                "dntp_molar": resolved_conditions.dntp_molar,
                "strand_concentration_molar": (resolved_conditions.strand_concentration_molar),
                "dmso_percent": resolved_conditions.dmso_percent,
                "dmso_factor_celsius_per_percent": (
                    resolved_conditions.dmso_factor_celsius_per_percent
                ),
                "formamide_molar": resolved_conditions.formamide_molar,
                "dinucleotide_property_table": (
                    None
                    if dinucleotide_property_table is None
                    else {
                        "name": dinucleotide_property_table.name,
                        "version": dinucleotide_property_table.version,
                        "source": dinucleotide_property_table.source,
                        "sha256": dinucleotide_property_table.sha256,
                        "schema_version": DINUCLEOTIDE_TABLE_SCHEMA_VERSION,
                    }
                ),
                "dinucleotide_sd": "population",
            }
        ),
        provenance=_provenance(dinucleotide_property_table),
    )


__all__ = ["all_descriptors"]
