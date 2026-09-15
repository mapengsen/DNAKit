"""Task metadata for the released Enformer NT Revised and GB checkpoints."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError

EnformerBenchmarkFamily = Literal["nt", "gb"]

ENFORMER_BENCHMARK_CHECKPOINTS_URL = (
    "https://drive.google.com/drive/folders/1lrZXzkrgAJMqM0wAmnIeZ4DEp0XFNIRI?usp=sharing"
)


_CHECKPOINT_FILE_IDS = {
    "demo_coding_vs_intergenomic_seqs": "1deBuRxYeh8yFjNSRK5rH58IIirgkdtoa",
    "demo_human_or_worm": "1woHoMjKG4sZSqFMSwUcLLoe-w4dMc-8o",
    "drosophila_enhancers_stark": "1ZHfykdEoXFM9s3PVI_j-i3Y9GjathrFd",
    "dummy_mouse_enhancers_ensembl": "17vi4Q1NiaBaKfY_Kti1uu-Ww9rg5Ywot",
    "enhancers": "15fkN_7ZftfdHzi0In9_NcskOm2ZvjZdF",
    "enhancers_types": "129kCJHEpHqXc_TNxaOj0W73S2JUQUmgW",
    "H2AFZ": "1IX_Yqgqu-OkquLgGeGy46ftC1_kj_CVq",
    "H3K27ac": "1Mpc3CQ0D7r7QL3JXfG50dfQMdZ1DlSLY",
    "H3K27me3": "17xIy_3vRu6APrzB9nlqWziaNJ-iGu5hG",
    "H3K36me3": "1ec35jpkQWiyNpZLFBzUVKqEpjKxuTlk1",
    "H3K4me1": "1pvpjpqJ8f5DcVg6fPDaN319Nsmf-F77F",
    "H3K4me2": "171mim8ZVsaha86WRAC3dxxEu46S4olUK",
    "H3K4me3": "1PfVfHUrhkbhbzcdDu9KUtTPtsZOCrA4p",
    "H3K9ac": "1F-nXU-rwSXYBAaDslOrrdPlU0LcB7Aem",
    "H3K9me3": "1whSdl2NwgGQtW3v8Odf30FNUsrdyVtyD",
    "H4K20me1": "1ONQkaRBA4JOQkzSPXEE9wgknXWnGRy0b",
    "human_enhancers_cohn": "18yvR0d1QpUkDhZc6VSvaYYdCPMmmJtl1",
    "human_enhancers_ensembl": "1RQNkzydJZ45i0mHv-sRLy8av0oT_pY9j",
    "human_ensembl_regulatory": "11D6D_-92Jq5e5vQAEjkkUhjpz0o4LBcF",
    "human_nontata_promoters": "1YZ2VKQROLY8VkkzfCIAJzjEjfoBFK_50",
    "human_ocr_ensembl": "1HJ6uN-DeCQfF3F-IyyJehU4bT8b9LW4s",
    "promoter_all": "1fHy8I5MVzi49gpnABlILXSHdzlKnSTm_",
    "promoter_no_tata": "1fPmCojazAfQVIqHDxlVueMpcgoEhMIOH",
    "promoter_tata": "1WzJcLuUh7vnXWpSmpgSRnPVuczqtKj8P",
    "splice_sites_acceptors": "1DzJtmO9ig-weMAo6WkDU_LDbNVMELP14",
    "splice_sites_all": "16Jz8vVjOJylEn4IgIU4ohOu9xw9muR-8",
    "splice_sites_donors": "1NuiRZmf8sH1G02aanVwtC99XYKK7v2bq",
}


@dataclass(frozen=True, slots=True)
class EnformerBenchmarkTask:
    """One task-specific, fully fine-tuned Enformer checkpoint."""

    name: str
    checkpoint_task: str
    checkpoint_filename: str
    family: EnformerBenchmarkFamily
    display_name: str
    description: str
    num_classes: int
    dataset_name: str

    @property
    def checkpoint_file_id(self) -> str:
        """Google Drive file ID for this task only."""
        return _CHECKPOINT_FILE_IDS[self.checkpoint_task]

    @property
    def checkpoint_url(self) -> str:
        """Public download source for this task."""
        return f"https://drive.google.com/file/d/{self.checkpoint_file_id}/view?usp=drivesdk"


def _spec(
    checkpoint_task: str,
    family: EnformerBenchmarkFamily,
    display_name: str,
    description: str,
    *,
    num_classes: int = 2,
) -> EnformerBenchmarkTask:
    name = checkpoint_task.lower()
    dataset_name = (
        "InstaDeepAI/nucleotide_transformer_downstream_tasks_revised"
        if family == "nt"
        else f"genomic_benchmarks/{checkpoint_task}"
    )
    return EnformerBenchmarkTask(
        name=name,
        checkpoint_task=checkpoint_task,
        checkpoint_filename=f"{checkpoint_task}.ckpt",
        family=family,
        display_name=display_name,
        description=description,
        num_classes=num_classes,
        dataset_name=dataset_name,
    )


_TASKS = (
    _spec("H2AFZ", "nt", "H2A.Z occupancy", "presence of an H2A.Z ChIP-seq peak"),
    _spec("H3K27ac", "nt", "H3K27ac occupancy", "presence of an H3K27ac ChIP-seq peak"),
    _spec("H3K27me3", "nt", "H3K27me3 occupancy", "presence of an H3K27me3 ChIP-seq peak"),
    _spec("H3K36me3", "nt", "H3K36me3 occupancy", "presence of an H3K36me3 ChIP-seq peak"),
    _spec("H3K4me1", "nt", "H3K4me1 occupancy", "presence of an H3K4me1 ChIP-seq peak"),
    _spec("H3K4me2", "nt", "H3K4me2 occupancy", "presence of an H3K4me2 ChIP-seq peak"),
    _spec("H3K4me3", "nt", "H3K4me3 occupancy", "presence of an H3K4me3 ChIP-seq peak"),
    _spec("H3K9ac", "nt", "H3K9ac occupancy", "presence of an H3K9ac ChIP-seq peak"),
    _spec("H3K9me3", "nt", "H3K9me3 occupancy", "presence of an H3K9me3 ChIP-seq peak"),
    _spec("H4K20me1", "nt", "H4K20me1 occupancy", "presence of an H4K20me1 ChIP-seq peak"),
    _spec("enhancers", "nt", "Enhancer detection", "enhancer versus non-enhancer sequence"),
    _spec(
        "enhancers_types",
        "nt",
        "Enhancer type",
        "none, tissue-specific, or tissue-invariant enhancer class",
        num_classes=3,
    ),
    _spec("promoter_all", "nt", "Promoter detection", "promoter versus non-promoter sequence"),
    _spec(
        "promoter_no_tata",
        "nt",
        "Non-TATA promoter detection",
        "non-TATA promoter versus non-promoter sequence",
    ),
    _spec(
        "promoter_tata",
        "nt",
        "TATA promoter detection",
        "TATA promoter versus non-promoter sequence",
    ),
    _spec(
        "splice_sites_acceptors",
        "nt",
        "Splice acceptor detection",
        "splice acceptor versus non-acceptor sequence",
    ),
    _spec(
        "splice_sites_all",
        "nt",
        "Splice-site class",
        "no-splice, splice-acceptor, or splice-donor class",
        num_classes=3,
    ),
    _spec(
        "splice_sites_donors",
        "nt",
        "Splice donor detection",
        "splice donor versus non-donor sequence",
    ),
    _spec(
        "demo_coding_vs_intergenomic_seqs",
        "gb",
        "Coding versus intergenic sequence",
        "coding versus intergenic genomic sequence",
    ),
    _spec(
        "demo_human_or_worm",
        "gb",
        "Human versus worm sequence",
        "human versus Caenorhabditis elegans sequence",
    ),
    _spec(
        "drosophila_enhancers_stark",
        "gb",
        "Drosophila enhancer detection",
        "Drosophila enhancer versus matched negative sequence",
    ),
    _spec(
        "dummy_mouse_enhancers_ensembl",
        "gb",
        "Mouse enhancer detection",
        "mouse Ensembl enhancer versus matched negative sequence",
    ),
    _spec(
        "human_enhancers_cohn",
        "gb",
        "Human Cohn enhancer detection",
        "human Cohn enhancer versus non-enhancer sequence",
    ),
    _spec(
        "human_enhancers_ensembl",
        "gb",
        "Human Ensembl enhancer detection",
        "human Ensembl enhancer versus matched negative sequence",
    ),
    _spec(
        "human_ensembl_regulatory",
        "gb",
        "Human regulatory-element class",
        "human enhancer, promoter, or open-chromatin-region class",
        num_classes=3,
    ),
    _spec(
        "human_nontata_promoters",
        "gb",
        "Human non-TATA promoter detection",
        "human non-TATA promoter versus negative sequence",
    ),
    _spec(
        "human_ocr_ensembl",
        "gb",
        "Human open-chromatin detection",
        "human Ensembl open-chromatin region versus matched negative sequence",
    ),
)

ENFORMER_BENCHMARK_TASKS: dict[str, EnformerBenchmarkTask] = {item.name: item for item in _TASKS}

_ALIASES = {
    "enhancer_types": "enhancers_types",
    "splice_sites_acceptor": "splice_sites_acceptors",
    "splice_sites_donor": "splice_sites_donors",
    "coding_vs_intergenic": "demo_coding_vs_intergenomic_seqs",
    "human_or_worm": "demo_human_or_worm",
    "human_non_tata_promoters": "human_nontata_promoters",
}


def _normalize_task(task: str) -> str:
    if not isinstance(task, str) or not task.strip():
        raise ConfigurationError(
            "Enformer benchmark task must be non-empty text.",
            code="INVALID_PREDICTION_TASK",
        )
    normalized = task.strip().lower().replace("-", "_").replace(" ", "_")
    if normalized.endswith(".ckpt"):
        normalized = normalized[:-5]
    normalized = _ALIASES.get(normalized, normalized)
    if normalized.startswith("nt_") and normalized[3:] in ENFORMER_BENCHMARK_TASKS:
        candidate = normalized[3:]
        if ENFORMER_BENCHMARK_TASKS[candidate].family == "nt":
            return candidate
    if normalized.startswith("gb_") and normalized[3:] in ENFORMER_BENCHMARK_TASKS:
        candidate = normalized[3:]
        if ENFORMER_BENCHMARK_TASKS[candidate].family == "gb":
            return candidate
    return normalized


def get_enformer_benchmark_task(task: str, *, verbose: bool = True) -> EnformerBenchmarkTask:
    """Resolve a canonical task, checkpoint filename, or ``nt_``/``gb_`` alias."""

    canonical = _normalize_task(task)
    try:
        return ENFORMER_BENCHMARK_TASKS[canonical]
    except KeyError as exc:
        raise ConfigurationError(
            "Unknown Enformer NT Revised or Genomic Benchmarks task.",
            code="INVALID_PREDICTION_TASK",
            context={"task": task, "available": available_enformer_benchmark_tasks()},
        ) from exc


def available_enformer_benchmark_tasks(
    family: EnformerBenchmarkFamily | None = None, *, verbose: bool = True
) -> tuple[str, ...]:
    """Return the 27 task names, optionally restricted to ``nt`` or ``gb``."""

    if family not in {None, "nt", "gb"}:
        raise ConfigurationError(
            "family must be None, 'nt', or 'gb'.",
            code="INVALID_PREDICTION_CONFIG",
        )
    return tuple(item.name for item in _TASKS if family is None or item.family == family)


def is_enformer_benchmark_task(task: str, *, verbose: bool = True) -> bool:
    """Return whether a canonical Enformer task uses one of the 27 local checkpoints."""

    try:
        return _normalize_task(task) in ENFORMER_BENCHMARK_TASKS
    except ConfigurationError:
        return False


__all__ = [
    "ENFORMER_BENCHMARK_CHECKPOINTS_URL",
    "ENFORMER_BENCHMARK_TASKS",
    "EnformerBenchmarkFamily",
    "EnformerBenchmarkTask",
    "available_enformer_benchmark_tasks",
    "get_enformer_benchmark_task",
    "is_enformer_benchmark_task",
]

_register_public_verbose(globals())
