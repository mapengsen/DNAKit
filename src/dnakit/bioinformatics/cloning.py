"""Native double-stranded molecules, cloning history and qualitative gel simulation."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit.bioinformatics._common import input_path, positive, require, sequence_text
from dnakit.exceptions import ConfigurationError


def _rc(text: str) -> str:
    return text.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@dataclass(frozen=True)
class Duplex:
    """Two explicit 5'->3' strands. offset locates Crick's 3' end relative to Watson's 5'.

    offset=2 means a two-base Watson 5' overhang at the left. Negative offset
    means Crick extends to the left. Only fully complementary paired regions
    are supported; gaps/nicks/mismatches are rejected.
    """

    watson: str
    crick: str
    offset: int = 0
    circular: bool = False

    def __post_init__(self) -> None:
        if (
            not self.watson
            or not self.crick
            or set(self.watson + self.crick) - set("ACGT")
            or type(self.offset) is not int
            or type(self.circular) is not bool
        ):
            raise ConfigurationError(
                "Duplex requires canonical uppercase strands and integer offset."
            )
        start, end = max(0, self.offset), min(len(self.watson), self.offset + len(self.crick))
        if (
            start >= end
            or self.watson[start:end] != _rc(self.crick)[start - self.offset : end - self.offset]
        ):
            raise ConfigurationError("Duplex strands must overlap complementarily.")
        if self.circular and (self.offset != 0 or len(self.watson) != len(self.crick)):
            raise ConfigurationError("Circular duplexes cannot contain overhangs.")

    @classmethod
    def from_sequence(cls, sequence: Any, *, circular: bool = False) -> Duplex:
        text = sequence_text(sequence)
        return cls(text, _rc(text), circular=circular)

    @property
    def span(self) -> int:
        return max(len(self.watson), self.offset + len(self.crick)) - min(0, self.offset)

    @property
    def ends(self) -> dict[str, dict[str, Any]]:
        if self.circular:
            return {}
        left = {"polarity": "blunt", "sequence": ""}
        if self.offset > 0:
            left = {"polarity": "5prime", "sequence": self.watson[: self.offset]}
        elif self.offset < 0:
            left = {"polarity": "3prime", "sequence": self.crick[self.offset :]}
        extension = self.offset + len(self.crick) - len(self.watson)
        right = {"polarity": "blunt", "sequence": ""}
        if extension > 0:
            right = {"polarity": "5prime", "sequence": self.crick[:extension]}
        elif extension < 0:
            right = {"polarity": "3prime", "sequence": self.watson[extension:]}
        return {"left": left, "right": right}

    def cut(self, watson_position: int, crick_position: int) -> tuple[Duplex, Duplex]:
        """Staggered cut; both positions use left-to-right Watson coordinates."""
        bottom = crick_position - self.offset
        if (
            self.circular
            or not 0 < watson_position < len(self.watson)
            or not 0 < bottom < len(self.crick)
        ):
            raise ConfigurationError(
                "cut currently requires interior positions on a linear duplex."
            )
        return (
            Duplex(self.watson[:watson_position], self.crick[-bottom:], self.offset),
            Duplex(
                self.watson[watson_position:],
                self.crick[:-bottom],
                crick_position - watson_position,
            ),
        )

    def reverse_complement(self) -> Duplex:
        return Duplex(
            self.crick, self.watson, self.offset + len(self.crick) - len(self.watson), self.circular
        )

    def ligate(self, other: Duplex) -> Duplex:
        """Join compatible ends, assuming ligatable 5' phosphates and no nicks remain."""
        if self.circular or other.circular:
            raise ConfigurationError("Ligation requires linear duplex fragments.")
        right, left = self.ends["right"], other.ends["left"]
        if right["polarity"] != left["polarity"] or right["sequence"] != _rc(left["sequence"]):
            raise ConfigurationError("Fragment ends are incompatible.")
        return Duplex(self.watson + other.watson, other.crick + self.crick, self.offset)


@dataclass(frozen=True)
class CloningStep:
    name: str
    operation: str
    parents: tuple[str, ...]
    product: Any
    parameters: Mapping[str, Any]


class CloningPlan:
    """Executable Python cloning history retaining every intermediate product."""

    def __init__(self) -> None:
        self._steps: dict[str, CloningStep] = {}

    @property
    def steps(self) -> tuple[CloningStep, ...]:
        return tuple(self._steps.values())

    def add(self, name: str, product: Any) -> CloningStep:
        return self.apply(name, lambda: product, operation="input")

    def apply(
        self,
        name: str,
        function: Callable[..., Any],
        *,
        parents: Iterable[str] = (),
        operation: str | None = None,
        **parameters: Any,
    ) -> CloningStep:
        parents = tuple(parents)
        if not name or name in self._steps or any(parent not in self._steps for parent in parents):
            raise ConfigurationError("Step names must be unique and parents must already exist.")
        product = function(*(self._steps[parent].product for parent in parents), **parameters)
        step = CloningStep(name, operation or function.__name__, parents, product, dict(parameters))
        self._steps[name] = step
        return step


def assemble_duplexes(
    fragments: Iterable[Duplex],
    *,
    circular: bool = False,
    minimum_overlap: int = 25,
    max_assemblies: int = 50,
) -> tuple[Duplex, ...]:
    """Specialist pydna overlap-graph assembly for Gibson/homologous recombination.

    The molecular object and ordinary cuts/ligation remain native DNAKit.
    This explicit optional backend enumerates complex assembly products.
    """
    positive(minimum_overlap, "minimum_overlap")
    dseq = require("pydna.dseq", "bio-assembly")
    record = require("pydna.dseqrecord", "bio-assembly")
    assembly = require("pydna.assembly2", "bio-assembly")
    positive(max_assemblies, "max_assemblies")
    molecules = tuple(fragments)
    if len(molecules) < 2 or any(fragment.circular for fragment in molecules):
        raise ConfigurationError("Provide at least two linear duplex fragments.")
    native = [
        record.Dseqrecord(dseq.Dseq(fragment.watson, fragment.crick, ovhg=-fragment.offset))
        for fragment in molecules
    ]
    graph = assembly.Assembly(
        native, limit=minimum_overlap, use_fragment_order=False, use_all_fragments=True
    )
    products = (
        graph.assemble_circular(max_assemblies=max_assemblies)
        if circular
        else graph.assemble_linear(max_assemblies=max_assemblies)
    )
    return tuple(
        Duplex(
            str(product.seq.watson).upper(),
            str(product.seq.crick).upper(),
            -int(product.seq.ovhg),
            bool(product.circular),
        )
        for product in products
    )


def read_sanger(path: str | Path, *, error_probability_cutoff: float = 0.05) -> dict[str, Any]:
    """ABI codec plus native Mott maximum-score trimming and chromatogram channels."""
    if not 0 < error_probability_cutoff < 1:
        raise ConfigurationError("error_probability_cutoff must lie in (0,1).")
    record = require("Bio.SeqIO", "bio-formats").read(str(input_path(path)), "abi")
    qualities = tuple(record.letter_annotations.get("phred_quality", ()))
    best_score = running = 0.0
    best_start = best_end = start = 0
    for index, quality in enumerate(qualities):
        running += error_probability_cutoff - 10 ** (-quality / 10)
        if running > best_score:
            best_score, best_start, best_end = running, start, index + 1
        if running < 0:
            running, start = 0.0, index + 1
    raw = record.annotations.get("abif_raw", {})
    order = raw.get("FWO_1", b"GATC")
    if isinstance(order, bytes):
        order = order.decode("ascii")
    channels = {base: tuple(raw.get(f"DATA{9 + index}", ())) for index, base in enumerate(order)}
    return {
        "id": record.id,
        "sequence": str(record.seq),
        "qualities": qualities,
        "trim_interval": (best_start, best_end),
        "trimmed_sequence": str(record.seq)[best_start:best_end],
        "channels": channels,
        "peak_positions": tuple(raw.get("PLOC2", ())),
    }


def plot_chromatogram(trace: Mapping[str, Any], *, axes: Any = None) -> Any:
    """Plot ABI channel intensities and called bases; return an undisplayed Axes."""
    plt = require("matplotlib.pyplot", "bio-plots")
    if axes is None:
        _, axes = plt.subplots()
    colors = {"A": "green", "C": "blue", "G": "black", "T": "red"}
    for base, values in trace["channels"].items():
        axes.plot(values, color=colors.get(base, "gray"), label=base)
    for peak, base in zip(trace["peak_positions"], trace["sequence"], strict=True):
        axes.text(peak, 1, base, transform=axes.get_xaxis_transform(), ha="center")
    axes.set(xlabel="Trace sample index", ylabel="Signal intensity")
    axes.legend()
    return axes


def sanger_consensus(
    traces: Iterable[Mapping[str, Any]],
    reference: Any,
    *,
    min_quality: int = 20,
    min_identity: float = 0.8,
    min_read_coverage: float = 0.8,
    min_depth: int = 1,
) -> dict[str, Any]:
    """Native local mapping of trimmed Sanger calls, both strands, and reference consensus.

    Reuses DNAKit's bounded pairwise aligner. Chooses the best strand (ties use
    forward), excludes weak/short mappings, then votes among bases above the
    quality cutoff. Ties/uncovered sites are N; indels remain in returned
    alignments but are not inserted into the reference-coordinate consensus.
    """
    from dnakit.alignment import AlignmentConfig, align_pairwise
    from dnakit.core import DNASequence

    np = require("numpy", "bio-analysis")
    target = sequence_text(reference)
    positive(min_depth, "min_depth")
    if not 0 <= min_quality <= 93 or not 0 < min_identity <= 1 or not 0 < min_read_coverage <= 1:
        raise ConfigurationError("Invalid Sanger quality/identity/coverage thresholds.")
    counts = np.zeros((len(target), 4), dtype=int)
    mappings = []
    for trace in traces:
        text = sequence_text(trace["sequence"])
        quality = tuple(trace["qualities"])
        if len(quality) != len(text) or any(
            type(q) is not int or not 0 <= q <= 93 for q in quality
        ):
            raise ConfigurationError("Sanger calls require one valid Phred quality per base.")
        start, end = trace.get("trim_interval", (0, len(text)))
        if not 0 <= start <= end <= len(text):
            raise ConfigurationError("Invalid Sanger trim interval.")
        if start == end:
            mappings.append({"id": trace["id"], "accepted": False, "reason": "empty trim"})
            continue
        text, quality = text[start:end], quality[start:end]
        config = AlignmentConfig(mode="local")
        forward = align_pairwise(DNASequence(text), DNASequence(target), config=config)
        reverse = align_pairwise(DNASequence(_rc(text)), DNASequence(target), config=config)
        flipped = reverse.score > forward.score
        selected = reverse if flipped else forward
        if flipped:
            quality = quality[::-1]
        accepted = (
            selected.identity is not None
            and selected.identity >= min_identity
            and selected.query_coverage >= min_read_coverage
        )
        mappings.append(
            {
                "id": trace["id"],
                "accepted": accepted,
                "strand": "-" if flipped else "+",
                "alignment": selected,
            }
        )
        if accepted:
            for column in selected.columns:
                qpos, rpos = column.query_position, column.target_position
                if (
                    qpos is not None
                    and rpos is not None
                    and quality[qpos] >= min_quality
                    and column.query_symbol in "ACGT"
                ):
                    counts[rpos, "ACGT".index(column.query_symbol)] += 1
    depth = counts.sum(axis=1)
    consensus = "".join(
        "ACGT"[int(row.argmax())] if depth[i] >= min_depth and sum(row == row.max()) == 1 else "N"
        for i, row in enumerate(counts)
    )
    return {"consensus": consensus, "counts": counts, "depth": depth, "mappings": tuple(mappings)}


def read_constructs(
    path: str | Path, *, format: str = "snapgene", permissive: bool = False
) -> tuple[Any, ...]:
    """SnapGene/GenBank codecs return Bio.SeqRecord; permissive parsing is opt-in.

    Permissive mixed GenBank/EMBL/FASTA extraction uses pydna's specialist parser.
    No silent repair of corrupt input is performed in strict mode.
    """
    source = input_path(path)
    if permissive:
        return tuple(require("pydna.parsers", "bio-assembly").parse(str(source), is_path=True))
    return tuple(require("Bio.SeqIO", "bio-formats").parse(str(source), format))


def simulate_gel(
    lanes: Iterable[Iterable[int | Duplex]],
    *,
    axes: Any = None,
    min_bp: int = 100,
    max_bp: int = 10_000,
) -> Any:
    """Qualitative log-size band plot; no claim of calibrated electrophoretic migration.

    Equal input amounts are assumed. Linear dsDNA only; topology-dependent
    migration, gel concentration and run time require experimental calibration.
    """
    positive(min_bp, "min_bp")
    if max_bp <= min_bp:
        raise ConfigurationError("max_bp must exceed min_bp.")
    plt = require("matplotlib.pyplot", "bio-plots")
    if axes is None:
        _, axes = plt.subplots()
    resolved = tuple(tuple(lane) for lane in lanes)
    for lane_index, lane in enumerate(resolved):
        for item in lane:
            if isinstance(item, Duplex):
                if item.circular:
                    raise ConfigurationError("This gel model supports linear dsDNA only.")
                length = item.span
            else:
                length = positive(item, "band length")
            axes.plot(
                [lane_index - 0.3, lane_index + 0.3], [length, length], linewidth=4, color="white"
            )
    axes.set_facecolor("#172337")
    axes.set_yscale("log")
    axes.set_ylim(min_bp, max_bp)
    axes.set_xticks(range(len(resolved)), [str(index + 1) for index in range(len(resolved))])
    axes.set(xlabel="Lane", ylabel="Fragment length (bp)", title="Qualitative linear-dsDNA gel")
    return axes
