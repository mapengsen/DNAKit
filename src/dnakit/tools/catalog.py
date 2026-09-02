"""Discover stable public DNAKit functions and describe them as Agent tools."""

from __future__ import annotations

import importlib
import inspect
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

from .adapters import CallPlan, analyze_callable

ToolEffect = Literal["local", "network", "filesystem", "model", "external_process"]

PUBLIC_TOOL_MODULES: tuple[str, ...] = (
    "dnakit.alignment",
    "dnakit.annotation",
    "dnakit.backends",
    "dnakit.batch",
    "dnakit.chunking",
    "dnakit.comparative",
    "dnakit.config",
    "dnakit.core",
    "dnakit.datasets",
    "dnakit.descriptors",
    "dnakit.download",
    "dnakit.evaluation",
    "dnakit.fingerprints",
    "dnakit.io",
    "dnakit.molbio",
    "dnakit.ops",
    "dnakit.patterns",
    "dnakit.predictions",
    "dnakit.references",
    "dnakit.representations",
    "dnakit.search",
    "dnakit.secondary_structure",
    "dnakit.similarity",
    "dnakit.standardize",
    "dnakit.structure3d",
    "dnakit.thermodynamics",
    "dnakit.visualization",
    "dnakit.workflows",
)

CATEGORY_LABELS: dict[str, tuple[str, str]] = {
    "alignment": ("序列比对", "alignment"),
    "annotation": ("变异注释", "variant annotation"),
    "backends": ("计算后端", "computational backends"),
    "batch": ("批量处理", "batch processing"),
    "chunking": ("序列切分", "sequence chunking"),
    "comparative": ("比较基因组学", "comparative genomics"),
    "config": ("配置", "configuration"),
    "core": ("核心对象", "core objects"),
    "datasets": ("数据集处理", "dataset processing"),
    "descriptors": ("序列描述符", "sequence descriptors"),
    "download": ("数据下载", "data download"),
    "evaluation": ("数据评价", "dataset evaluation"),
    "fingerprints": ("DNA 指纹", "DNA fingerprints"),
    "io": ("文件读写", "file input and output"),
    "molbio": ("分子生物学", "molecular biology"),
    "ops": ("序列操作", "sequence operations"),
    "patterns": ("功能模式扫描", "functional pattern scanning"),
    "predictions": ("性质预测", "property prediction"),
    "references": ("参考基因组", "reference genomes"),
    "representations": ("模型表征", "model representations"),
    "search": ("数据库搜索", "database search"),
    "secondary_structure": ("二级结构", "secondary structure"),
    "similarity": ("相似性", "similarity"),
    "standardize": ("标准化与验证", "standardization and validation"),
    "structure3d": ("三维结构", "3D structure"),
    "thermodynamics": ("理化与热力学", "physicochemical and thermodynamic"),
    "visualization": ("可视化", "visualization"),
    "workflows": ("工作流", "workflows"),
}

_DENIED_TOOLS: dict[str, str] = {
    "dnakit.backends.execute_bounded_command": (
        "The low-level arbitrary command executor is intentionally unavailable to Agents."
    ),
}
_FILESYSTEM_NAMES = frozenset(
    {
        "dnakit.io.build_fasta_index",
        "dnakit.io.build_fastq_index",
        "dnakit.io.export_result",
        "dnakit.io.export_table",
        "dnakit.io.write",
        "dnakit.io.write_agp",
        "dnakit.io.write_bed",
        "dnakit.io.write_gff3",
        "dnakit.references.download_genome",
        "dnakit.similarity.save_sketch_index",
        "dnakit.visualization.save_html_report",
        "dnakit.visualization.save_image",
        "dnakit.visualization.save_svg",
        "dnakit.workflows.run_workflow",
        "dnakit.workflows.save_manifest",
    }
)
_MODEL_NAMES = frozenset(
    {
        "dnakit.datasets.neural_cluster_sequences",
        "dnakit.evaluation.evaluate_frechet_distance",
        "dnakit.predictions.ensure_prediction_checkpoint",
        "dnakit.predictions.predict_enformer_benchmark",
        "dnakit.predictions.predict_pair_properties",
        "dnakit.predictions.predict_properties",
        "dnakit.predictions.predict_sequence_properties",
        "dnakit.predictions.predict_variant_effects",
        "dnakit.representations.ensure_model_checkpoint",
        "dnakit.representations.extract_representations",
    }
)
_EXTERNAL_PROCESS_NAMES = frozenset(
    {
        "dnakit.download.build_index",
        "dnakit.thermodynamics.probe_primer3",
    }
)
_NETWORK_PREFIXES = (
    "dnakit.annotation.",
    "dnakit.comparative.",
    "dnakit.search.",
)
_NETWORK_NAMES = frozenset(
    {
        "dnakit.molbio.assemble_golden_gate",
        "dnakit.molbio.design_golden_gate",
        "dnakit.references.resolve_genome_assembly",
    }
)


@dataclass(frozen=True, slots=True)
class AgentToolSpec:
    """One discovered public function and its Agent-facing metadata."""

    name: str
    category: str
    function_name: str
    description: str
    full_description: str
    effect: ToolEffect
    open_world: bool
    requires_confirmation: bool
    agent_compatible: bool
    incompatibilities: tuple[str, ...]
    hidden_parameters: tuple[str, ...]
    call_plan: CallPlan
    function: Callable[..., object]

    def to_summary(self) -> dict[str, object]:
        """Return compact metadata for discovery results."""

        labels = CATEGORY_LABELS.get(self.category, (self.category, self.category))
        return {
            "name": self.name,
            "category": self.category,
            "category_zh": labels[0],
            "description": self.description,
            "effect": self.effect,
            "agent_compatible": self.agent_compatible,
            "requires_confirmation": self.requires_confirmation,
        }

    def to_dict(self) -> dict[str, object]:
        """Return the complete manifest exposed by ``describe_dnakit_tool``."""

        payload = self.to_summary()
        payload.update(
            {
                "function": self.function_name,
                "full_description": self.full_description,
                "signature": self.call_plan.signature,
                "input_schema": self.call_plan.input_schema,
                "output_schema": self.call_plan.output_schema,
                "open_world": self.open_world,
                "incompatibilities": list(self.incompatibilities),
                "python_only_optional_parameters": list(self.hidden_parameters),
            }
        )
        return payload


def _description(function: Callable[..., object]) -> tuple[str, str]:
    full = inspect.getdoc(function) or f"Call the public DNAKit function {function.__name__}."
    full = full[:8_000]
    first_paragraph = full.split("\n\n", 1)[0].replace("\n", " ").strip()
    return first_paragraph[:500], full


def _effect(name: str) -> ToolEffect:
    if name in _MODEL_NAMES:
        return "model"
    if name in _EXTERNAL_PROCESS_NAMES:
        return "external_process"
    if name.startswith("dnakit.download.") or name in _FILESYSTEM_NAMES:
        return "filesystem"
    if name.startswith(_NETWORK_PREFIXES) or name in _NETWORK_NAMES:
        return "network"
    return "local"


def _build_spec(
    module_name: str,
    export_name: str,
    function: Callable[..., object],
) -> AgentToolSpec:
    category = module_name.rsplit(".", 1)[-1]
    name = f"{module_name}.{export_name}"
    plan = analyze_callable(function)
    denied_reason = _DENIED_TOOLS.get(name)
    incompatibilities = list(plan.incompatibilities)
    if denied_reason is not None:
        incompatibilities.insert(0, denied_reason)
    effect = _effect(name)
    short_description, full_description = _description(function)
    return AgentToolSpec(
        name=name,
        category=category,
        function_name=f"{function.__module__}.{function.__qualname__}",
        description=short_description,
        full_description=full_description,
        effect=effect,
        open_world=effect in {"network", "filesystem", "model", "external_process"},
        requires_confirmation=effect in {"filesystem", "model", "external_process"},
        agent_compatible=plan.compatible and denied_reason is None,
        incompatibilities=tuple(incompatibilities),
        hidden_parameters=tuple(
            parameter.name
            for parameter in plan.parameters
            if not parameter.supported and not parameter.required
        ),
        call_plan=plan,
        function=function,
    )


def build_tool_catalog() -> tuple[AgentToolSpec, ...]:
    """Discover functions explicitly exported by stable DNAKit domain modules."""

    tools: list[AgentToolSpec] = []
    for module_name in PUBLIC_TOOL_MODULES:
        module = importlib.import_module(module_name)
        for export_name in getattr(module, "__all__", ()):
            if not isinstance(export_name, str) or export_name.startswith("_"):
                continue
            value = getattr(module, export_name, None)
            if inspect.isfunction(value):
                tools.append(_build_spec(module_name, export_name, value))
    return tuple(sorted(tools, key=lambda item: item.name))


__all__ = [
    "CATEGORY_LABELS",
    "PUBLIC_TOOL_MODULES",
    "AgentToolSpec",
    "ToolEffect",
    "build_tool_catalog",
]
