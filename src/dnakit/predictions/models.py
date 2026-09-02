"""Registry, inputs, and configuration for pretrained property prediction."""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Literal, TypeAlias

from dnakit.exceptions import ConfigurationError, SequenceError

from .enformer_benchmarks import (
    ENFORMER_BENCHMARK_CHECKPOINTS_URL,
    ENFORMER_BENCHMARK_TASKS,
    get_enformer_benchmark_task,
    is_enformer_benchmark_task,
)

PredictionInputKind: TypeAlias = Literal["sequence", "pair", "variant"]
PredictionOutputKind: TypeAlias = Literal[
    "classification",
    "regression",
    "score",
    "segmentation",
    "tracks",
]
SequenceType: TypeAlias = Literal["gene", "protein"]
AmbiguityPolicy: TypeAlias = Literal["replace_with_n", "error"]
PredictionDType: TypeAlias = Literal["auto", "float32", "float16", "bfloat16"]

_GENE_SYMBOLS = frozenset("ACGTRYSWKMBDHVNU")
_PROTEIN_SYMBOLS = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZ*")
_CANONICAL_BASES = frozenset("ACGT")
_MAX_INPUT_SYMBOLS = 10_000_000


@dataclass(frozen=True, slots=True)
class DirectPredictionModel:
    """One model family with an official pretrained prediction path."""

    name: str
    display_name: str
    source_repository: str
    checkpoint_ids: tuple[str, ...]
    notes: str = ""


@dataclass(frozen=True, slots=True)
class DirectPredictionTask:
    """One task supported by an official pretrained checkpoint or model head."""

    model_name: str
    name: str
    display_name: str
    input_kinds: tuple[PredictionInputKind, ...]
    output_kind: PredictionOutputKind
    description: str
    requires_remote_code: bool = False
    notes: str = ""


MODEL_REGISTRY: dict[str, DirectPredictionModel] = {
    "alphagenome": DirectPredictionModel(
        "alphagenome",
        "AlphaGenome all folds",
        "https://github.com/google-deepmind/alphagenome_research",
        ("google/alphagenome-all-folds",),
        "The checkpoint is gated by Google's non-commercial model terms.",
    ),
    "enformer": DirectPredictionModel(
        "enformer",
        "Enformer",
        "https://github.com/google-deepmind/deepmind-research/tree/master/enformer",
        ("EleutherAI/enformer-official-rough", ENFORMER_BENCHMARK_CHECKPOINTS_URL),
        (
            "DNAKit uses the documented PyTorch port for the released regulatory-track "
            "weights and supports local task-specific NT Revised/GB checkpoints."
        ),
    ),
    "evo2": DirectPredictionModel(
        "evo2",
        "Evo 2",
        "https://github.com/ArcInstitute/evo2",
        (
            "arcinstitute/evo2_7b",
            "arcinstitute/evo2_7b_base",
            "schmojo/evo2-exon-classifier",
        ),
    ),
    "generator": DirectPredictionModel(
        "generator",
        "GENERator v2 eukaryote 1.2B",
        "https://github.com/GenerTeam/GENERator",
        ("GenerTeam/GENERator-v2-eukaryote-1.2b-base",),
    ),
    "lucaone": DirectPredictionModel(
        "lucaone",
        "LucaOneTasks",
        "https://github.com/LucaOne/LucaOneTasks",
        ("LucaOneTasks/DownstreamTasksTrainedModels",),
        "The official source checkout and its dedicated environment are required.",
    ),
    "segmentnt": DirectPredictionModel(
        "segmentnt",
        "SegmentNT (Nucleotide Transformer v2 head)",
        "https://github.com/instadeepai/nucleotide-transformer",
        ("InstaDeepAI/segment_nt",),
        "This is the trained SegmentNT head, not the NT-v2 base checkpoint.",
    ),
}


def _task(
    model: str,
    name: str,
    display: str,
    input_kinds: tuple[PredictionInputKind, ...],
    output_kind: PredictionOutputKind,
    description: str,
    *,
    remote_code: bool = False,
    notes: str = "",
) -> DirectPredictionTask:
    return DirectPredictionTask(
        model,
        name,
        display,
        input_kinds,
        output_kind,
        description,
        remote_code,
        notes,
    )


_ALPHAGENOME_OUTPUTS = {
    "atac": ("ATAC-seq", "chromatin accessibility"),
    "cage": ("CAGE", "transcription initiation and gene expression"),
    "dnase": ("DNase-seq", "chromatin accessibility"),
    "rna_seq": ("RNA-seq", "RNA coverage and gene expression"),
    "chip_histone": ("histone ChIP-seq", "histone modification tracks"),
    "chip_tf": ("TF ChIP-seq", "transcription-factor binding tracks"),
    "splice_sites": ("splice sites", "donor and acceptor splice-site tracks"),
    "splice_site_usage": ("splice-site usage", "splice-site usage fractions"),
    "splice_junctions": ("splice junctions", "splice-junction read-count tracks"),
    "contact_maps": ("contact maps", "3D DNA contact probabilities"),
    "procap": ("PRO-cap", "nascent transcription initiation"),
}


TASK_REGISTRY: dict[tuple[str, str], DirectPredictionTask] = {}
for _name, (_display, _description) in _ALPHAGENOME_OUTPUTS.items():
    TASK_REGISTRY[("alphagenome", _name)] = _task(
        "alphagenome",
        _name,
        _display,
        ("sequence",) if _name == "splice_junctions" else ("sequence", "variant"),
        "tracks",
        _description,
    )

TASK_REGISTRY.update(
    {
        ("enformer", "human_tracks"): _task(
            "enformer",
            "human_tracks",
            "Human regulatory tracks",
            ("sequence", "variant"),
            "tracks",
            "5,313 human CAGE, DNase/ATAC, TF and histone ChIP tracks at 128-bp bins",
        ),
        ("enformer", "mouse_tracks"): _task(
            "enformer",
            "mouse_tracks",
            "Mouse regulatory tracks",
            ("sequence", "variant"),
            "tracks",
            "1,643 mouse CAGE, DNase/ATAC, TF and histone ChIP tracks at 128-bp bins",
        ),
        ("segmentnt", "genomic_segmentation"): _task(
            "segmentnt",
            "genomic_segmentation",
            "14-class genomic segmentation",
            ("sequence",),
            "segmentation",
            "single-nucleotide probabilities for 14 human gene and regulatory annotations",
            remote_code=True,
        ),
        ("evo2", "variant_effect"): _task(
            "evo2",
            "variant_effect",
            "Zero-shot variant effect",
            ("variant",),
            "score",
            "alternate-minus-reference sequence likelihood score",
            notes="This is a continuous score, not a calibrated pathogenicity probability.",
        ),
        ("evo2", "exon_probability"): _task(
            "evo2",
            "exon_probability",
            "Exon probability",
            ("pair",),
            "classification",
            "exonic probability from forward and reverse Evo 2 context embeddings",
            remote_code=True,
        ),
        ("generator", "variant_effect"): _task(
            "generator",
            "variant_effect",
            "Zero-shot variant effect",
            ("variant",),
            "score",
            "log ratio of reference-allele and alternate-allele next-token probabilities",
            remote_code=True,
            notes="This is a continuous score, not a calibrated pathogenicity probability.",
        ),
        ("lucaone", "central_dogma"): _task(
            "lucaone",
            "central_dogma",
            "Central dogma relation",
            ("pair",),
            "classification",
            "DNA/RNA and protein relationship classification",
        ),
        ("lucaone", "supktax"): _task(
            "lucaone",
            "supktax",
            "SupKTax taxonomy",
            ("sequence",),
            "classification",
            "SupKTax checkpoint taxonomy classification for gene sequences",
        ),
        ("lucaone", "genustax"): _task(
            "lucaone",
            "genustax",
            "GenusTax taxonomy",
            ("sequence",),
            "classification",
            "GenusTax checkpoint taxonomy classification for gene sequences",
        ),
        ("lucaone", "speciestax"): _task(
            "lucaone",
            "speciestax",
            "Species taxonomy",
            ("sequence",),
            "classification",
            "species-level taxonomy classification for gene sequences",
        ),
        ("lucaone", "protein_location"): _task(
            "lucaone",
            "protein_location",
            "Protein subcellular location",
            ("sequence",),
            "classification",
            "prokaryotic protein subcellular-location classification",
        ),
        ("lucaone", "protein_stability"): _task(
            "lucaone",
            "protein_stability",
            "Protein stability",
            ("sequence",),
            "regression",
            "protein stability regression",
        ),
        ("lucaone", "ncrna_family"): _task(
            "lucaone",
            "ncrna_family",
            "ncRNA family",
            ("sequence",),
            "classification",
            "non-coding RNA family classification",
        ),
        ("lucaone", "influenza_antigenicity"): _task(
            "lucaone",
            "influenza_antigenicity",
            "Influenza A antigenic relation",
            ("pair",),
            "classification",
            "pairwise influenza A antigenic-relationship classification",
        ),
        ("lucaone", "protein_interaction"): _task(
            "lucaone",
            "protein_interaction",
            "Protein-protein interaction",
            ("pair",),
            "classification",
            "protein-protein interaction classification",
        ),
        ("lucaone", "ncrna_protein_interaction"): _task(
            "lucaone",
            "ncrna_protein_interaction",
            "ncRNA-protein interaction",
            ("pair",),
            "classification",
            "ncRNA-protein interaction classification",
        ),
    }
)

for _benchmark in ENFORMER_BENCHMARK_TASKS.values():
    TASK_REGISTRY[("enformer", _benchmark.name)] = _task(
        "enformer",
        _benchmark.name,
        _benchmark.display_name,
        ("sequence",),
        "classification",
        _benchmark.description,
        notes=(
            f"Fully fine-tuned checkpoint {_benchmark.checkpoint_filename} from "
            f"{_benchmark.dataset_name}; no fitting is performed at prediction time."
        ),
    )


_MODEL_ALIASES = {
    "alpha-genome": "alphagenome",
    "evo-2": "evo2",
    "generator-v2": "generator",
    "luca-one": "lucaone",
    "segment-nt": "segmentnt",
    "ntv2-segment": "segmentnt",
}

_TASK_ALIASES = {
    ("segmentnt", "segmentation"): "genomic_segmentation",
    ("evo2", "vep"): "variant_effect",
    ("generator", "vep"): "variant_effect",
    ("lucaone", "protloc"): "protein_location",
    ("lucaone", "protstab"): "protein_stability",
    ("lucaone", "ncrnafam"): "ncrna_family",
    ("lucaone", "infa"): "influenza_antigenicity",
    ("lucaone", "ppi"): "protein_interaction",
    ("lucaone", "ncrpi"): "ncrna_protein_interaction",
}


def _canonical_model(name: str) -> str:
    if not isinstance(name, str) or not name.strip():
        raise ConfigurationError("model must be non-empty text.", code="INVALID_PREDICTION_MODEL")
    key = name.strip().lower().replace("_", "-")
    canonical = _MODEL_ALIASES.get(key, key.replace("-", ""))
    if canonical not in MODEL_REGISTRY:
        raise ConfigurationError(
            "Unknown direct-prediction model.",
            code="INVALID_PREDICTION_MODEL",
            context={"model": name, "available": available_prediction_models()},
        )
    return canonical


def available_prediction_models() -> tuple[str, ...]:
    """Return model families with at least one official direct prediction task."""

    return tuple(sorted(MODEL_REGISTRY))


def available_prediction_tasks(model: str | None = None) -> tuple[str, ...]:
    """Return canonical task names, optionally limited to one model family."""

    if model is None:
        return tuple(sorted(f"{item.model_name}:{item.name}" for item in TASK_REGISTRY.values()))
    canonical = _canonical_model(model)
    return tuple(sorted(name for family, name in TASK_REGISTRY if family == canonical))


def get_prediction_model(name: str) -> DirectPredictionModel:
    """Resolve one canonical model name or documented alias."""

    return MODEL_REGISTRY[_canonical_model(name)]


def get_prediction_task(model: str, task: str) -> DirectPredictionTask:
    """Resolve a task under a direct-prediction model family."""

    canonical_model = _canonical_model(model)
    if not isinstance(task, str) or not task.strip():
        raise ConfigurationError("task must be non-empty text.", code="INVALID_PREDICTION_TASK")
    key = task.strip().lower().replace("-", "_").replace(" ", "_")
    if canonical_model == "enformer" and is_enformer_benchmark_task(key):
        key = get_enformer_benchmark_task(key).name
    canonical_task = _TASK_ALIASES.get((canonical_model, key), key)
    resolved = TASK_REGISTRY.get((canonical_model, canonical_task))
    if resolved is None:
        raise ConfigurationError(
            "The selected model has no registered direct prediction task with this name.",
            code="INVALID_PREDICTION_TASK",
            context={
                "model": canonical_model,
                "task": task,
                "available": available_prediction_tasks(canonical_model),
            },
        )
    return resolved


def _identifier(value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\r\n\x00"):
        raise SequenceError(
            f"{field_name} must be non-empty single-line text.",
            code="INVALID_PREDICTION_INPUT",
        )
    return value.strip()


def _sequence(value: object, sequence_type: SequenceType, field_name: str) -> str:
    if not isinstance(value, str):
        raise SequenceError(f"{field_name} must be sequence text.", code="INVALID_PREDICTION_INPUT")
    normalized = value.strip().upper()
    if not normalized:
        raise SequenceError(f"{field_name} must not be empty.", code="EMPTY_PREDICTION_INPUT")
    if len(normalized) > _MAX_INPUT_SYMBOLS:
        raise SequenceError(
            f"{field_name} exceeds the input-size limit.",
            code="PREDICTION_INPUT_LIMIT",
            context={"max_symbols": _MAX_INPUT_SYMBOLS},
        )
    allowed = _GENE_SYMBOLS if sequence_type == "gene" else _PROTEIN_SYMBOLS
    invalid = sorted(set(normalized) - allowed)
    if invalid:
        raise SequenceError(
            f"{field_name} contains symbols unsupported for {sequence_type} input.",
            code="INVALID_PREDICTION_INPUT_SYMBOL",
            context={"symbols": tuple(invalid)},
        )
    return normalized


@dataclass(frozen=True, slots=True)
class BiologicalSequence:
    """One gene (DNA/RNA) or protein sequence for a pretrained task head."""

    id: str
    sequence: str
    sequence_type: SequenceType = "gene"

    def __post_init__(self) -> None:
        if self.sequence_type not in {"gene", "protein"}:
            raise SequenceError(
                "sequence_type must be 'gene' or 'protein'.",
                code="INVALID_PREDICTION_INPUT",
            )
        object.__setattr__(self, "id", _identifier(self.id, "id"))
        object.__setattr__(
            self,
            "sequence",
            _sequence(self.sequence, self.sequence_type, "sequence"),
        )


@dataclass(frozen=True, slots=True)
class BiologicalSequencePair:
    """Two typed biological sequences used by a pairwise pretrained task."""

    id: str
    first: BiologicalSequence
    second: BiologicalSequence

    def __post_init__(self) -> None:
        object.__setattr__(self, "id", _identifier(self.id, "id"))
        if not isinstance(self.first, BiologicalSequence) or not isinstance(
            self.second, BiologicalSequence
        ):
            raise SequenceError(
                "first and second must be BiologicalSequence objects.",
                code="INVALID_PREDICTION_INPUT",
            )


@dataclass(frozen=True, slots=True)
class VariantContext:
    """Reference and alternate equal-length contexts for one SNV."""

    id: str
    reference_sequence: str
    alternate_sequence: str
    variant_index: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "id", _identifier(self.id, "id"))
        reference = _sequence(self.reference_sequence, "gene", "reference_sequence")
        alternate = _sequence(self.alternate_sequence, "gene", "alternate_sequence")
        reference = reference.replace("U", "T")
        alternate = alternate.replace("U", "T")
        if len(reference) != len(alternate):
            raise SequenceError(
                "Variant contexts must have equal lengths.",
                code="INVALID_VARIANT_CONTEXT",
            )
        differences = tuple(
            index
            for index, (ref, alt) in enumerate(zip(reference, alternate, strict=True))
            if ref != alt
        )
        if self.variant_index is None:
            if len(differences) != 1:
                raise SequenceError(
                    "A variant context must contain exactly one changed position.",
                    code="INVALID_VARIANT_CONTEXT",
                    context={"difference_count": len(differences)},
                )
            resolved_index = differences[0]
        else:
            if (
                isinstance(self.variant_index, bool)
                or not isinstance(self.variant_index, int)
                or not 0 <= self.variant_index < len(reference)
            ):
                raise SequenceError(
                    "variant_index must be a valid 0-based sequence position.",
                    code="INVALID_VARIANT_CONTEXT",
                )
            resolved_index = self.variant_index
            if differences != (resolved_index,):
                raise SequenceError(
                    "Only variant_index may differ between reference and alternate contexts.",
                    code="INVALID_VARIANT_CONTEXT",
                    context={"difference_positions": differences[:20]},
                )
        if (
            reference[resolved_index] not in _CANONICAL_BASES
            or alternate[resolved_index] not in _CANONICAL_BASES
        ):
            raise SequenceError(
                "The reference and alternate SNV alleles must be A, C, G, or T.",
                code="INVALID_VARIANT_CONTEXT",
            )
        object.__setattr__(self, "reference_sequence", reference)
        object.__setattr__(self, "alternate_sequence", alternate)
        object.__setattr__(self, "variant_index", resolved_index)

    @property
    def reference_base(self) -> str:
        """Return the reference allele at :attr:`variant_index`."""

        assert self.variant_index is not None
        return self.reference_sequence[self.variant_index]

    @property
    def alternate_base(self) -> str:
        """Return the alternate allele at :attr:`variant_index`."""

        assert self.variant_index is not None
        return self.alternate_sequence[self.variant_index]


PredictionInput: TypeAlias = BiologicalSequence | BiologicalSequencePair | VariantContext


@dataclass(frozen=True, slots=True)
class PropertyPredictionConfig:
    """Configure one official pretrained property-prediction task."""

    model: str
    task: str
    checkpoint_dir: str | os.PathLike[str] | None = None
    checkpoint_path: str | os.PathLike[str] | None = None
    model_source_path: str | os.PathLike[str] | None = None
    device: str = "auto"
    dtype: PredictionDType = "auto"
    batch_size: int = 1
    max_length: int | None = None
    max_records: int = 1_000
    ambiguity_policy: AmbiguityPolicy = "replace_with_n"
    organism: Literal["human", "mouse"] = "human"
    ontology_terms: tuple[str, ...] = ()
    threshold: float = 0.5
    top_k: int = 5
    show_progress: bool = True
    allow_remote_code: bool = False
    hf_token: str | None = field(default=None, repr=False)
    timeout_seconds: float = 86_400.0
    max_backend_output_bytes: int = 10_000_000

    def __post_init__(self) -> None:
        task_spec = get_prediction_task(self.model, self.task)
        object.__setattr__(self, "model", task_spec.model_name)
        object.__setattr__(self, "task", task_spec.name)
        for field_name in ("checkpoint_dir", "checkpoint_path", "model_source_path"):
            value = getattr(self, field_name)
            if value is None:
                continue
            try:
                resolved = os.fspath(value)
            except TypeError as exc:
                raise ConfigurationError(
                    f"{field_name} must be path-like or None.",
                    code="INVALID_PREDICTION_CONFIG",
                ) from exc
            if not resolved.strip():
                raise ConfigurationError(
                    f"{field_name} must not be empty.", code="INVALID_PREDICTION_CONFIG"
                )
            object.__setattr__(self, field_name, resolved)
        if self.checkpoint_dir is not None and self.checkpoint_path is not None:
            raise ConfigurationError(
                "checkpoint_dir and checkpoint_path are mutually exclusive.",
                code="INVALID_PREDICTION_CONFIG",
            )
        if not isinstance(self.device, str) or not self.device.strip():
            raise ConfigurationError(
                "device must be non-empty text.", code="INVALID_PREDICTION_CONFIG"
            )
        if self.dtype not in {"auto", "float32", "float16", "bfloat16"}:
            raise ConfigurationError("Unknown prediction dtype.", code="INVALID_PREDICTION_CONFIG")
        for name, value, maximum in (
            ("batch_size", self.batch_size, 1_024),
            ("max_records", self.max_records, 100_000),
            ("top_k", self.top_k, 100),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
                raise ConfigurationError(
                    f"{name} must be in [1, {maximum}].",
                    code="INVALID_PREDICTION_CONFIG",
                )
        if self.max_length is not None and (
            isinstance(self.max_length, bool)
            or not isinstance(self.max_length, int)
            or not 1 <= self.max_length <= _MAX_INPUT_SYMBOLS
        ):
            raise ConfigurationError(
                f"max_length must be None or an integer in [1, {_MAX_INPUT_SYMBOLS}].",
                code="INVALID_PREDICTION_CONFIG",
            )
        if self.ambiguity_policy not in {"replace_with_n", "error"}:
            raise ConfigurationError("Unknown ambiguity_policy.", code="INVALID_PREDICTION_CONFIG")
        if self.organism not in {"human", "mouse"}:
            raise ConfigurationError(
                "organism must be 'human' or 'mouse'.", code="INVALID_PREDICTION_CONFIG"
            )
        terms: list[str] = []
        for term in self.ontology_terms:
            if not isinstance(term, str) or not term.strip():
                raise ConfigurationError(
                    "ontology_terms must contain non-empty strings.",
                    code="INVALID_PREDICTION_CONFIG",
                )
            terms.append(term.strip())
        object.__setattr__(self, "ontology_terms", tuple(terms))
        if (
            isinstance(self.threshold, bool)
            or not isinstance(self.threshold, (int, float))
            or not math.isfinite(self.threshold)
            or not 0.0 <= self.threshold <= 1.0
        ):
            raise ConfigurationError(
                "threshold must be finite and in [0, 1].",
                code="INVALID_PREDICTION_CONFIG",
            )
        for name in ("show_progress", "allow_remote_code"):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(
                    f"{name} must be boolean.", code="INVALID_PREDICTION_CONFIG"
                )
        if self.hf_token is not None and (
            not isinstance(self.hf_token, str) or not self.hf_token.strip()
        ):
            raise ConfigurationError(
                "hf_token must be None or non-empty text.",
                code="INVALID_PREDICTION_CONFIG",
            )
        if (
            isinstance(self.timeout_seconds, bool)
            or not isinstance(self.timeout_seconds, (int, float))
            or not math.isfinite(self.timeout_seconds)
            or not 0 < self.timeout_seconds <= 86_400
        ):
            raise ConfigurationError(
                "timeout_seconds must be in (0, 86400].",
                code="INVALID_PREDICTION_CONFIG",
            )
        if (
            isinstance(self.max_backend_output_bytes, bool)
            or not isinstance(self.max_backend_output_bytes, int)
            or not 1 <= self.max_backend_output_bytes <= 100_000_000
        ):
            raise ConfigurationError(
                "max_backend_output_bytes must be in [1, 100000000].",
                code="INVALID_PREDICTION_CONFIG",
            )

    @property
    def model_spec(self) -> DirectPredictionModel:
        """Return the resolved model registry entry."""

        return get_prediction_model(self.model)

    @property
    def task_spec(self) -> DirectPredictionTask:
        """Return the resolved task registry entry."""

        return get_prediction_task(self.model, self.task)


__all__ = [
    "MODEL_REGISTRY",
    "TASK_REGISTRY",
    "AmbiguityPolicy",
    "BiologicalSequence",
    "BiologicalSequencePair",
    "DirectPredictionModel",
    "DirectPredictionTask",
    "PredictionDType",
    "PredictionInput",
    "PredictionInputKind",
    "PredictionOutputKind",
    "PropertyPredictionConfig",
    "SequenceType",
    "VariantContext",
    "available_prediction_models",
    "available_prediction_tasks",
    "get_prediction_model",
    "get_prediction_task",
]
