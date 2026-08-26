"""Immutable DNA sequence value object."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass

from dnakit.core.enums import DNAAlphabet, Strandedness, Topology
from dnakit.core.gap import Gap
from dnakit.exceptions import (
    InvalidAlphabetError,
    SequenceError,
    UnknownLengthError,
    UnsupportedGapOperationError,
)

STRICT_SYMBOLS = frozenset("ACGT")
IUPAC_SYMBOLS = frozenset("ACGTRYSWKMBDHVN")
_COMPLEMENT = str.maketrans(
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


def _coerce_alphabet(value: DNAAlphabet | str) -> DNAAlphabet:
    try:
        return value if isinstance(value, DNAAlphabet) else DNAAlphabet(value)
    except (TypeError, ValueError) as exc:
        raise SequenceError(
            "Unknown DNA alphabet.",
            code="INVALID_ALPHABET_DECLARATION",
            context={"alphabet": value},
            hint=f"Choose one of: {', '.join(item.value for item in DNAAlphabet)}.",
        ) from exc


def _coerce_topology(value: Topology | str) -> Topology:
    try:
        return value if isinstance(value, Topology) else Topology(value)
    except (TypeError, ValueError) as exc:
        raise SequenceError(
            "Unknown sequence topology.",
            code="INVALID_TOPOLOGY",
            context={"topology": value},
            hint=f"Choose one of: {', '.join(item.value for item in Topology)}.",
        ) from exc


def _coerce_strandedness(value: Strandedness | str) -> Strandedness:
    try:
        return value if isinstance(value, Strandedness) else Strandedness(value)
    except (TypeError, ValueError) as exc:
        raise SequenceError(
            "Unknown sequence strandedness.",
            code="INVALID_STRANDEDNESS",
            context={"strandedness": value},
            hint=f"Choose one of: {', '.join(item.value for item in Strandedness)}.",
        ) from exc


def _decode_ascii(value: bytes) -> str:
    try:
        return value.decode("ascii")
    except UnicodeDecodeError as exc:
        raise InvalidAlphabetError(
            "DNA bytes must use ASCII encoding.",
            context={"byte_offset": exc.start},
            hint="Decode or clean the raw input with normalize() before construction.",
        ) from exc


def _normalize_parts(parts: str | bytes | Iterable[str | Gap]) -> tuple[str | Gap, ...]:
    if isinstance(parts, bytes):
        raw_parts: Iterable[str | Gap] = (_decode_ascii(parts),)
    elif isinstance(parts, str):
        raw_parts = (parts,)
    else:
        try:
            raw_parts = iter(parts)
        except TypeError as exc:
            raise SequenceError(
                "Sequence parts must be text, bytes, or an iterable of text and Gap objects.",
                code="INVALID_SEQUENCE_PARTS",
                context={"parts_type": type(parts).__name__},
            ) from exc

    normalized: list[str | Gap] = []
    text_buffer: list[str] = []

    def flush_text() -> None:
        if text_buffer:
            normalized.append("".join(text_buffer))
            text_buffer.clear()

    for index, part in enumerate(raw_parts):
        if isinstance(part, str):
            if part:
                text_buffer.append(part)
        elif isinstance(part, Gap):
            flush_text()
            normalized.append(part)
        else:
            raise SequenceError(
                "Each sequence part must be a string or Gap.",
                code="INVALID_SEQUENCE_PART",
                context={"part_index": index, "part_type": type(part).__name__},
            )
    flush_text()
    return tuple(normalized)


@dataclass(frozen=True, init=False)
class DNASequence:
    """Canonical immutable DNA symbols and explicit gap parts.

    The constructor validates already-normalized content.  It deliberately does
    not uppercase, strip whitespace, or replace uracil; callers should use the
    normalization API for those auditable transformations.
    """

    parts: tuple[str | Gap, ...]
    alphabet: DNAAlphabet
    topology: Topology
    strandedness: Strandedness

    def __init__(
        self,
        parts: str | bytes | Iterable[str | Gap],
        *,
        alphabet: DNAAlphabet | str = DNAAlphabet.STRICT,
        topology: Topology | str = Topology.LINEAR,
        strandedness: Strandedness | str = Strandedness.SINGLE,
    ) -> None:
        resolved_parts = _normalize_parts(parts)
        resolved_alphabet = _coerce_alphabet(alphabet)
        resolved_topology = _coerce_topology(topology)
        resolved_strandedness = _coerce_strandedness(strandedness)

        allowed = STRICT_SYMBOLS if resolved_alphabet is DNAAlphabet.STRICT else IUPAC_SYMBOLS
        for part_index, part in enumerate(resolved_parts):
            if isinstance(part, Gap):
                continue
            invalid = next(
                ((offset, symbol) for offset, symbol in enumerate(part) if symbol not in allowed),
                None,
            )
            if invalid is not None:
                offset, symbol = invalid
                raise InvalidAlphabetError(
                    f"Symbol {symbol!r} is not valid for the {resolved_alphabet.value} alphabet.",
                    context={
                        "alphabet": resolved_alphabet.value,
                        "part_index": part_index,
                        "part_offset": offset,
                        "symbol": symbol,
                    },
                    hint="Pass raw or mixed-case input through normalize() first.",
                )

        symbol_length = sum(len(part) for part in resolved_parts if isinstance(part, str))
        if resolved_topology is Topology.CIRCULAR and symbol_length == 0:
            raise SequenceError(
                "A circular DNA sequence must contain at least one nucleotide symbol.",
                code="EMPTY_CIRCULAR_SEQUENCE",
                hint="Use linear topology for an empty placeholder sequence.",
            )

        object.__setattr__(self, "parts", resolved_parts)
        object.__setattr__(self, "alphabet", resolved_alphabet)
        object.__setattr__(self, "topology", resolved_topology)
        object.__setattr__(self, "strandedness", resolved_strandedness)

    @classmethod
    def from_fragments(
        cls,
        fragments: Iterable[str | bytes],
        gaps: Iterable[Gap],
        *,
        alphabet: DNAAlphabet | str = DNAAlphabet.STRICT,
        topology: Topology | str = Topology.LINEAR,
        strandedness: Strandedness | str = Strandedness.SINGLE,
    ) -> DNASequence:
        """Interleave ``n`` gaps between exactly ``n + 1`` sequence fragments."""

        fragment_tuple = tuple(
            _decode_ascii(fragment) if isinstance(fragment, bytes) else fragment
            for fragment in fragments
        )
        gap_tuple = tuple(gaps)
        if len(fragment_tuple) != len(gap_tuple) + 1:
            raise SequenceError(
                "from_fragments() requires exactly one more fragment than gaps.",
                code="FRAGMENT_GAP_COUNT_MISMATCH",
                context={"fragment_count": len(fragment_tuple), "gap_count": len(gap_tuple)},
            )
        if any(not isinstance(fragment, str) for fragment in fragment_tuple):
            raise SequenceError(
                "Fragments must be strings or ASCII bytes.",
                code="INVALID_SEQUENCE_FRAGMENT",
            )
        if any(not isinstance(gap, Gap) for gap in gap_tuple):
            raise SequenceError("All gaps must be Gap objects.", code="INVALID_GAP")

        parts: list[str | Gap] = []
        for index, fragment in enumerate(fragment_tuple):
            parts.append(fragment)
            if index < len(gap_tuple):
                parts.append(gap_tuple[index])
        return cls(
            parts,
            alphabet=alphabet,
            topology=topology,
            strandedness=strandedness,
        )

    @property
    def symbol_length(self) -> int:
        """Number of nucleotide symbols, excluding explicit gap spans."""

        return sum(len(part) for part in self.parts if isinstance(part, str))

    @property
    def canonical_base_count(self) -> int:
        return sum(
            sum(symbol in STRICT_SYMBOLS for symbol in part)
            for part in self.parts
            if isinstance(part, str)
        )

    @property
    def ambiguity_count(self) -> int:
        return self.symbol_length - self.canonical_base_count

    @property
    def coordinate_span(self) -> int | None:
        span = self.symbol_length
        for part in self.parts:
            if isinstance(part, Gap):
                if part.length is None:
                    return None
                span += part.length
        return span

    @property
    def length(self) -> int | None:
        """Explicit alias for :attr:`coordinate_span`."""

        return self.coordinate_span

    @property
    def is_gapped(self) -> bool:
        return any(isinstance(part, Gap) for part in self.parts)

    @property
    def has_unknown_length(self) -> bool:
        return any(isinstance(part, Gap) and part.length is None for part in self.parts)

    @property
    def symbols(self) -> str:
        """Concatenate nucleotide symbols while explicitly omitting gaps."""

        return "".join(part for part in self.parts if isinstance(part, str))

    def to_string(self) -> str:
        """Return an ordinary string only when doing so cannot discard gaps."""

        if self.is_gapped:
            raise UnsupportedGapOperationError(
                "A gapped DNASequence cannot be represented as a plain DNA string.",
                hint="Inspect parts or use symbols explicitly if omitting gaps is intended.",
            )
        return self.symbols

    def __len__(self) -> int:
        span = self.coordinate_span
        if span is None:
            raise UnknownLengthError(
                "Exact sequence length is unknown because at least one gap has unknown length.",
                hint="Use symbol_length or resolve every Gap.length before requesting len().",
            )
        return span

    def __iter__(self) -> Iterator[str | Gap]:
        return iter(self.parts)

    def reverse(self) -> DNASequence:
        return self._transform(reverse=True, complement=False)

    def complement(self) -> DNASequence:
        return self._transform(reverse=False, complement=True)

    def reverse_complement(self) -> DNASequence:
        return self._transform(reverse=True, complement=True)

    def _transform(self, *, reverse: bool, complement: bool) -> DNASequence:
        source: Sequence[str | Gap] = self.parts[::-1] if reverse else self.parts
        transformed: list[str | Gap] = []
        for part in source:
            if isinstance(part, Gap):
                transformed.append(part)
                continue
            text = part[::-1] if reverse else part
            transformed.append(text.translate(_COMPLEMENT) if complement else text)
        return DNASequence(
            transformed,
            alphabet=self.alphabet,
            topology=self.topology,
            strandedness=self.strandedness,
        )


__all__ = ["IUPAC_SYMBOLS", "STRICT_SYMBOLS", "DNASequence"]
