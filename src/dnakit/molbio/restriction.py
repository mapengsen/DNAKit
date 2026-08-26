"""Bounded restriction digestion, end classification, and ligation."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from itertools import pairwise

from dnakit.core import DNAAlphabet, DNASequence, Issue, IssueSeverity, Topology
from dnakit.exceptions import ConfigurationError
from dnakit.patterns import RestrictionEnzyme, scan_restriction_sites

from ._shared import (
    circular_slice,
    freeze_parameters,
    materialize_bounded,
    native_provenance,
    require_sequence,
    reverse_complement_text,
    validate_positive_int,
    validate_text,
)
from .results import (
    DigestCut,
    DigestFragment,
    EndDescriptor,
    EndPolarity,
    EndSide,
    EndTypeResult,
    LigationCompatibilityResult,
    LigationResult,
    RestrictionDigestResult,
)


@dataclass(frozen=True)
class LigationFragment:
    """A linear nucleotide fragment plus its two abstract duplex ends."""

    id: str
    sequence: DNASequence
    left_end: EndDescriptor
    right_end: EndDescriptor

    def __post_init__(self) -> None:
        validate_text(self.id, "fragment id")
        require_sequence(
            self.sequence,
            operation="LigationFragment construction",
            max_length=100_000_000,
            allow_circular=False,
        )
        if self.left_end.side != "left" or self.right_end.side != "right":
            raise ConfigurationError("Fragment ends must be assigned to their matching sides.")


def _cohesive_key(overhang: str) -> str:
    reverse = reverse_complement_text(overhang)
    return min(overhang, reverse)


def _end_descriptor(
    sequence: DNASequence,
    *,
    top_cut: int,
    bottom_cut: int,
    side: EndSide,
    five_prime_phosphorylated: bool,
    cut_span: int | None,
    source: str,
) -> EndDescriptor:
    symbols = require_sequence(
        sequence,
        operation="restriction-end classification",
        max_length=100_000_000,
    )
    if not isinstance(five_prime_phosphorylated, bool):
        raise ConfigurationError("five_prime_phosphorylated must be boolean.")
    for value, name in ((top_cut, "top_cut"), (bottom_cut, "bottom_cut")):
        if isinstance(value, bool) or not isinstance(value, int):
            raise ConfigurationError(f"{name} must be an integer.")
        upper = len(symbols) - 1 if sequence.topology is Topology.CIRCULAR else len(symbols)
        if not 0 <= value <= upper:
            raise ConfigurationError(f"{name} lies outside the template boundary.")
    delta = bottom_cut - top_cut if cut_span is None else cut_span
    if isinstance(delta, bool) or not isinstance(delta, int):
        raise ConfigurationError("cut_span must be an integer or None.")
    if abs(delta) > len(symbols):
        raise ConfigurationError("The staggered cut span cannot exceed the sequence length.")
    if delta == 0:
        polarity: EndPolarity = "blunt"
        overhang = ""
    elif delta > 0:
        polarity = "5prime"
        overhang = (
            circular_slice(symbols, top_cut, delta)
            if sequence.topology is Topology.CIRCULAR
            else symbols[top_cut:bottom_cut]
        )
    else:
        polarity = "3prime"
        overhang = (
            circular_slice(symbols, bottom_cut, -delta)
            if sequence.topology is Topology.CIRCULAR
            else symbols[bottom_cut:top_cut]
        )
    if len(overhang) != abs(delta):
        raise ConfigurationError(
            "Cut coordinates and cut_span are inconsistent with the sequence boundary."
        )
    return EndDescriptor(
        polarity=polarity,
        overhang_sequence_5to3=overhang,
        cohesive_key=_cohesive_key(overhang) if overhang else "",
        side=side,
        five_prime_phosphorylated=five_prime_phosphorylated,
        top_cut=top_cut,
        bottom_cut=bottom_cut,
        source=source,
    )


def classify_restriction_end(
    sequence: DNASequence,
    top_cut: int,
    bottom_cut: int,
    *,
    side: EndSide = "right",
    five_prime_phosphorylated: bool = True,
    cut_span: int | None = None,
) -> EndTypeResult:
    """Classify one blunt, 5-prime, or 3-prime end.

    For modular circular cut coordinates, ``cut_span`` must retain the signed
    bottom-minus-top offset from the enzyme definition.
    """

    end = _end_descriptor(
        sequence,
        top_cut=top_cut,
        bottom_cut=bottom_cut,
        side=side,
        five_prime_phosphorylated=five_prime_phosphorylated,
        cut_span=cut_span,
        source="explicit-cut-coordinates",
    )
    return EndTypeResult(
        end=end,
        method="signed_staggered_cut_classification",
        algorithm_version="dnakit-end-type-v1",
        parameters=freeze_parameters(
            {
                "top_cut": top_cut,
                "bottom_cut": bottom_cut,
                "cut_span": cut_span,
                "side": side,
                "five_prime_phosphorylated": five_prime_phosphorylated,
                "coordinate_system": "0-based-interbase",
            }
        ),
        provenance=native_provenance(),
        issues=(),
    )


def _natural_end(side: EndSide) -> EndDescriptor:
    return EndDescriptor(
        polarity="blunt",
        overhang_sequence_5to3="",
        cohesive_key="",
        side=side,
        five_prime_phosphorylated=False,
        top_cut=None,
        bottom_cut=None,
        source="natural-linear-terminus",
    )


def _uncut_circular_end(side: EndSide) -> EndDescriptor:
    return EndDescriptor(
        polarity="blunt",
        overhang_sequence_5to3="",
        cohesive_key="",
        side=side,
        five_prime_phosphorylated=False,
        top_cut=None,
        bottom_cut=None,
        source="uncut-circular-no-physical-end",
    )


def _slice_product(sequence: DNASequence, start: int, length: int) -> DNASequence:
    symbols = (
        circular_slice(sequence.symbols, start, length)
        if sequence.topology is Topology.CIRCULAR
        else sequence.symbols[start : start + length]
    )
    return DNASequence(
        symbols,
        alphabet=sequence.alphabet,
        topology=Topology.LINEAR,
        strandedness=sequence.strandedness,
    )


def digest_restriction(
    sequence: DNASequence,
    enzymes: Iterable[str | RestrictionEnzyme],
    *,
    methylation_state: str = "unmethylated",
    allow_ambiguous_template: bool = False,
    five_prime_phosphorylated: bool = True,
    max_fragments: int = 100_000,
    max_enzymes: int = 10_000,
    max_sequence_length: int = 10_000_000,
) -> RestrictionDigestResult:
    """Simulate complete simultaneous digestion by an audited enzyme panel."""

    symbols = require_sequence(
        sequence,
        operation="restriction digestion",
        max_length=max_sequence_length,
        canonical=not allow_ambiguous_template,
    )
    if methylation_state != "unmethylated":
        raise ConfigurationError(
            "Native digestion models only an unmethylated template.",
            code="METHYLATION_MODEL_UNAVAILABLE",
        )
    if not isinstance(allow_ambiguous_template, bool) or not isinstance(
        five_prime_phosphorylated, bool
    ):
        raise ConfigurationError("Boolean digestion controls must be booleans.")
    validate_positive_int(max_fragments, "max_fragments", maximum=1_000_000)
    validate_positive_int(max_enzymes, "max_enzymes", maximum=100_000)
    definitions = materialize_bounded(
        enzymes,
        max_items=max_enzymes,
        name="restriction enzymes",
    )
    if not definitions:
        raise ConfigurationError("At least one restriction enzyme is required.")
    resolved_by_name: dict[str, RestrictionEnzyme] = {}
    from dnakit.patterns import BUILTIN_RESTRICTION_ENZYMES

    for definition in definitions:
        if isinstance(definition, RestrictionEnzyme):
            resolved = definition
        elif isinstance(definition, str):
            try:
                resolved = BUILTIN_RESTRICTION_ENZYMES[definition]
            except KeyError as exc:
                raise ConfigurationError(
                    "Restriction enzyme is absent from the small built-in catalog.",
                    context={"enzyme": definition},
                ) from exc
        else:
            raise ConfigurationError("Each enzyme must be a name or RestrictionEnzyme.")
        if resolved.name in resolved_by_name:
            raise ConfigurationError("Restriction enzyme names must be unique in one digest.")
        resolved_by_name[resolved.name] = resolved
    scan = scan_restriction_sites(
        sequence,
        tuple(resolved_by_name.values()),
        max_matches=max_fragments,
        max_scan_length=max_sequence_length,
        max_enzymes=max_enzymes,
    )
    if scan.truncated:
        raise ConfigurationError(
            "Restriction site count reaches max_fragments; complete digestion is not safe.",
            code="DIGEST_FRAGMENT_LIMIT_EXCEEDED",
        )

    grouped: dict[tuple[int, int], list[str]] = {}
    cut_spans: dict[tuple[int, int], int] = {}
    top_geometries: dict[int, tuple[int, int]] = {}
    for hit in scan.hits:
        if hit.top_cut is None or hit.bottom_cut is None:
            raise ConfigurationError("Restriction digestion requires resolved cut coordinates.")
        definition = resolved_by_name[hit.enzyme]
        span = definition.bottom_cut_offset - definition.top_cut_offset
        geometry = (hit.top_cut, hit.bottom_cut)
        existing = top_geometries.get(hit.top_cut)
        if existing is not None and existing != geometry:
            raise ConfigurationError(
                "Simultaneous enzymes create conflicting bottom cuts at one top-strand cut.",
                code="CONFLICTING_DIGEST_GEOMETRY",
            )
        top_geometries[hit.top_cut] = geometry
        grouped.setdefault(geometry, []).append(hit.enzyme)
        cut_spans[geometry] = span

    cuts: list[DigestCut] = []
    left_ends: dict[int, EndDescriptor] = {}
    right_ends: dict[int, EndDescriptor] = {}
    for geometry in sorted(grouped):
        top_cut, bottom_cut = geometry
        span = cut_spans[geometry]
        right = _end_descriptor(
            sequence,
            top_cut=top_cut,
            bottom_cut=bottom_cut,
            side="right",
            five_prime_phosphorylated=five_prime_phosphorylated,
            cut_span=span,
            source="restriction-digest",
        )
        left = _end_descriptor(
            sequence,
            top_cut=top_cut,
            bottom_cut=bottom_cut,
            side="left",
            five_prime_phosphorylated=five_prime_phosphorylated,
            cut_span=span,
            source="restriction-digest",
        )
        left_ends[top_cut] = left
        right_ends[top_cut] = right
        cuts.append(
            DigestCut(
                enzymes=tuple(sorted(set(grouped[geometry]))),
                top_cut=top_cut,
                bottom_cut=bottom_cut,
                polarity=right.polarity,
                overhang_sequence_5to3=right.overhang_sequence_5to3,
            )
        )

    fragments: list[DigestFragment] = []
    top_cuts = sorted(left_ends)
    if sequence.topology is Topology.LINEAR:
        boundaries = sorted({0, *top_cuts, len(symbols)})
        for index, (start, end) in enumerate(pairwise(boundaries), 1):
            fragments.append(
                DigestFragment(
                    id=f"fragment_{index}",
                    sequence=_slice_product(sequence, start, end - start),
                    source_start=start,
                    source_end=end,
                    wraps_origin=False,
                    left_end=left_ends.get(start, _natural_end("left")),
                    right_end=(right_ends.get(end, _natural_end("right"))),
                )
            )
    elif not top_cuts:
        fragments.append(
            DigestFragment(
                id="fragment_1",
                sequence=sequence,
                source_start=0,
                source_end=len(symbols),
                wraps_origin=False,
                left_end=_uncut_circular_end("left"),
                right_end=_uncut_circular_end("right"),
            )
        )
    else:
        for index, start in enumerate(top_cuts, 1):
            end = top_cuts[index % len(top_cuts)]
            length = (end - start) % len(symbols)
            if length == 0:
                length = len(symbols)
            fragments.append(
                DigestFragment(
                    id=f"fragment_{index}",
                    sequence=_slice_product(sequence, start, length),
                    source_start=start,
                    source_end=end,
                    wraps_origin=end <= start,
                    left_end=left_ends[start],
                    right_end=right_ends[end],
                )
            )
    if len(fragments) > max_fragments:
        raise ConfigurationError("Digest product exceeds max_fragments.")
    issues: list[Issue] = []
    if allow_ambiguous_template and sequence.ambiguity_count:
        issues.append(
            Issue(
                "AMBIGUOUS_DIGEST_MATCHES",
                IssueSeverity.WARNING,
                "Ambiguous template symbols were interpreted by IUPAC compatibility.",
            )
        )
    if sequence.topology is Topology.CIRCULAR and not cuts:
        issues.append(
            Issue(
                "UNCUT_CIRCULAR_MOLECULE",
                IssueSeverity.INFO,
                "No physical fragment ends exist because the circular template was not cut.",
            )
        )
    return RestrictionDigestResult(
        source_length=len(symbols),
        source_topology=sequence.topology.value,
        cuts=tuple(cuts),
        fragments=tuple(fragments),
        method="complete_simultaneous_restriction_digest",
        algorithm_version="dnakit-restriction-digest-v1",
        parameters=freeze_parameters(
            {
                "enzymes": sorted(resolved_by_name),
                "methylation_state": methylation_state,
                "allow_ambiguous_template": allow_ambiguous_template,
                "five_prime_phosphorylated": five_prime_phosphorylated,
                "partial_digest_modeled": False,
                "star_activity_modeled": False,
                "max_fragments": max_fragments,
                "max_enzymes": max_enzymes,
                "max_sequence_length": max_sequence_length,
            }
        ),
        provenance=native_provenance(
            reimplementation=True,
            reference_name="DNAKit common restriction-enzyme definitions",
            reference_version="common-enzymes-v1",
        ),
        issues=tuple(issues),
    )


def as_ligation_fragment(fragment: DigestFragment | LigationFragment) -> LigationFragment:
    if isinstance(fragment, LigationFragment):
        return fragment
    if isinstance(fragment, DigestFragment):
        return LigationFragment(
            fragment.id, fragment.sequence, fragment.left_end, fragment.right_end
        )
    raise ConfigurationError("Expected DigestFragment or LigationFragment.")


def check_end_compatibility(
    left: EndDescriptor,
    right: EndDescriptor,
    *,
    allow_blunt: bool = False,
    require_phosphorylation: bool = True,
) -> LigationCompatibilityResult:
    """Check a right-facing end against a left-facing end."""

    if not isinstance(left, EndDescriptor) or not isinstance(right, EndDescriptor):
        raise ConfigurationError("Compatibility inputs must be EndDescriptor objects.")
    if not isinstance(allow_blunt, bool) or not isinstance(require_phosphorylation, bool):
        raise ConfigurationError("Compatibility controls must be booleans.")
    if left.side != "right" or right.side != "left":
        compatible, reason = False, "ends do not face an ordered right-to-left junction"
    elif require_phosphorylation and not (
        left.five_prime_phosphorylated and right.five_prime_phosphorylated
    ):
        compatible, reason = False, "the abstract model requires both 5-prime phosphates"
    elif left.polarity != right.polarity:
        compatible, reason = False, "overhang polarities differ"
    elif left.polarity == "blunt":
        compatible, reason = (
            allow_blunt,
            ("blunt ligation is enabled" if allow_blunt else "blunt ligation is disabled"),
        )
    elif left.cohesive_key != right.cohesive_key:
        compatible, reason = False, "cohesive overhang sequences are not complementary"
    else:
        compatible, reason = True, "cohesive overhangs are complementary"
    return LigationCompatibilityResult(
        compatible=compatible,
        reason=reason,
        left=left,
        right=right,
        method="abstract_duplex_end_compatibility",
        algorithm_version="dnakit-end-compatibility-v1",
        parameters=freeze_parameters(
            {
                "allow_blunt": allow_blunt,
                "require_phosphorylation": require_phosphorylation,
                "kinetics_modeled": False,
            }
        ),
        provenance=native_provenance(),
        issues=(),
    )


def ligate_fragments(
    fragments: Iterable[DigestFragment | LigationFragment],
    *,
    circularize: bool = False,
    allow_blunt: bool = False,
    require_phosphorylation: bool = True,
    max_fragments: int = 10_000,
    max_product_length: int = 100_000_000,
) -> LigationResult:
    """Ligate an explicitly ordered list of pre-cut fragments."""

    validate_positive_int(max_fragments, "max_fragments", maximum=1_000_000)
    validate_positive_int(max_product_length, "max_product_length", maximum=100_000_000)
    if not isinstance(circularize, bool):
        raise ConfigurationError("circularize must be boolean.")
    resolved: list[LigationFragment] = []
    for fragment in fragments:
        if len(resolved) >= max_fragments:
            raise ConfigurationError("Ligation input exceeds max_fragments.")
        resolved.append(as_ligation_fragment(fragment))
    if len(resolved) < 2:
        raise ConfigurationError("Ligation requires at least two fragments.")
    junctions = list(pairwise(resolved))
    if circularize:
        junctions.append((resolved[-1], resolved[0]))
    for left_fragment, right_fragment in junctions:
        compatibility = check_end_compatibility(
            left_fragment.right_end,
            right_fragment.left_end,
            allow_blunt=allow_blunt,
            require_phosphorylation=require_phosphorylation,
        )
        if not compatibility.compatible:
            raise ConfigurationError(
                "Fragment ends are incompatible for the requested ligation.",
                code="INCOMPATIBLE_LIGATION_ENDS",
                context={
                    "left_fragment": left_fragment.id,
                    "right_fragment": right_fragment.id,
                    "reason": compatibility.reason,
                },
            )
    product_length = sum(fragment.sequence.symbol_length for fragment in resolved)
    if product_length > max_product_length:
        raise ConfigurationError("Ligated product exceeds max_product_length.")
    symbols = "".join(fragment.sequence.symbols for fragment in resolved)
    alphabet = (
        DNAAlphabet.IUPAC
        if any(fragment.sequence.alphabet is DNAAlphabet.IUPAC for fragment in resolved)
        else DNAAlphabet.STRICT
    )
    if any(
        fragment.sequence.strandedness is not resolved[0].sequence.strandedness
        for fragment in resolved[1:]
    ):
        raise ConfigurationError("All ligation fragments must use the same strandedness.")
    product = DNASequence(
        symbols,
        alphabet=alphabet,
        topology=Topology.CIRCULAR if circularize else Topology.LINEAR,
        strandedness=resolved[0].sequence.strandedness,
    )
    return LigationResult(
        product=product,
        fragment_ids=tuple(fragment.id for fragment in resolved),
        junction_count=len(junctions),
        circularized=circularize,
        method="ordered_abstract_duplex_ligation",
        algorithm_version="dnakit-ligation-v1",
        parameters=freeze_parameters(
            {
                "allow_blunt": allow_blunt,
                "require_phosphorylation": require_phosphorylation,
                "circularize": circularize,
                "kinetics_modeled": False,
                "max_fragments": max_fragments,
                "max_product_length": max_product_length,
            }
        ),
        provenance=native_provenance(),
        issues=(),
    )


__all__ = [
    "LigationFragment",
    "as_ligation_fragment",
    "check_end_compatibility",
    "classify_restriction_end",
    "digest_restriction",
    "ligate_fragments",
]
