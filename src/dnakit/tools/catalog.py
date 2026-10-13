"""Discover stable public DNAKit functions and describe them as Agent tools."""

from __future__ import annotations

import importlib
import importlib.util
import inspect
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal, cast

from dnakit._verbosity import register_public as _register_public_verbose

from .adapters import CallPlan, analyze_callable

ToolEffect = Literal["local", "network", "filesystem", "model", "external_process"]

PUBLIC_TOOL_MODULES: tuple[str, ...] = (
    "dnakit.generated",
    "dnakit.alignment",
    "dnakit.annotation",
    "dnakit.assembly",
    "dnakit.backends",
    "dnakit.batch",
    "dnakit.breeding",
    "dnakit.clinical",
    "dnakit.chunking",
    "dnakit.comparative",
    "dnakit.config",
    "dnakit.core",
    "dnakit.datasets",
    "dnakit.design_crispr",
    "dnakit.databases",
    "dnakit.descriptors",
    "dnakit.download",
    "dnakit.evaluation",
    "dnakit.epigenome",
    "dnakit.encode",
    "dnakit.experimental_design",
    "dnakit.external_adapters",
    "dnakit.fasta_fastq",
    "dnakit.forensic_reproductive",
    "dnakit.fingerprints",
    "dnakit.genomics",
    "dnakit.genome3d",
    "dnakit.gwas",
    "dnakit.haplotypes",
    "dnakit.integrations",
    "dnakit.io",
    "dnakit.io_formats",
    "dnakit.knowledge_graph",
    "dnakit.microbiome",
    "dnakit.ml_models",
    "dnakit.molecular_lab",
    "dnakit.molbio",
    "dnakit.ngs",
    "dnakit.ops",
    "dnakit.patterns",
    "dnakit.phylogeny",
    "dnakit.population",
    "dnakit.predictions",
    "dnakit.references",
    "dnakit.regulation",
    "dnakit.replication_repair",
    "dnakit.representations",
    "dnakit.search",
    "dnakit.sequence",
    "dnakit.sequencing_qc",
    "dnakit.secondary_structure",
    "dnakit.similarity",
    "dnakit.singlecell_spatial",
    "dnakit.skills",
    "dnakit.standardize",
    "dnakit.structure3d",
    "dnakit.transcriptomics",
    "dnakit.thermodynamics",
    "dnakit.variants",
    "dnakit.visualization",
    "dnakit.workflow_security",
    "dnakit.workflows",
)

CATEGORY_LABELS: dict[str, tuple[str, str]] = {
    "generated": ("生成工具", "generated tools"),
    "encode": ("ENCODE 实验数据", "ENCODE experimental data"),
    "databases": ("科学数据库接口", "scientific database APIs"),
    "alignment": ("序列比对", "alignment"),
    "annotation": ("基因与变异注释", "gene and variant annotation"),
    "assembly": ("基因组组装", "genome assembly"),
    "backends": ("计算后端", "computational backends"),
    "batch": ("批量处理", "batch processing"),
    "breeding": ("育种分析", "breeding genomics"),
    "clinical": ("临床基因组学", "clinical genomics"),
    "chunking": ("序列切分", "sequence chunking"),
    "comparative": ("比较基因组学", "comparative genomics"),
    "config": ("配置", "configuration"),
    "core": ("核心对象", "core objects"),
    "datasets": ("数据集处理", "dataset processing"),
    "design_crispr": ("CRISPR 设计", "CRISPR design"),
    "descriptors": ("序列描述符", "sequence descriptors"),
    "download": ("数据下载", "data download"),
    "evaluation": ("数据评价", "dataset evaluation"),
    "epigenome": ("表观基因组学", "epigenomics"),
    "experimental_design": ("实验设计", "experimental design"),
    "external_adapters": ("外部工具适配", "external adapters"),
    "fasta_fastq": ("FASTA/FASTQ", "FASTA/FASTQ"),
    "fingerprints": ("DNA 指纹", "DNA fingerprints"),
    "forensic_reproductive": ("法医与生殖基因组学", "forensic and reproductive genomics"),
    "genomics": ("基因组与 DNA 工程", "genomics and DNA engineering"),
    "genome3d": ("三维基因组", "three-dimensional genome"),
    "gwas": ("GWAS", "GWAS"),
    "haplotypes": ("单倍型", "haplotypes"),
    "integrations": ("外部分析与 MCP", "external analyses and MCP"),
    "io": ("文件读写", "file input and output"),
    "io_formats": ("输入输出格式", "input/output formats"),
    "knowledge_graph": ("知识图谱", "knowledge graph"),
    "microbiome": ("微生物组", "microbiome"),
    "ml_models": ("机器学习模型", "machine-learning models"),
    "molecular_lab": ("分子实验", "molecular laboratory"),
    "molbio": ("分子生物学", "molecular biology"),
    "ops": ("序列操作", "sequence operations"),
    "patterns": ("功能模式扫描", "functional pattern scanning"),
    "phylogeny": ("系统发育", "phylogeny"),
    "population": ("群体遗传学", "population genetics"),
    "predictions": ("性质预测", "property prediction"),
    "references": ("参考基因组", "reference genomes"),
    "regulation": ("基因调控", "gene regulation"),
    "replication_repair": ("复制与 DNA 修复", "replication and DNA repair"),
    "representations": ("模型表征", "model representations"),
    "search": ("数据库搜索", "database search"),
    "sequence": ("序列分析", "sequence analysis"),
    "sequencing_qc": ("测序质量控制", "sequencing quality control"),
    "secondary_structure": ("二级结构", "secondary structure"),
    "similarity": ("相似性", "similarity"),
    "singlecell_spatial": ("单细胞与空间组学", "single-cell and spatial genomics"),
    "skills": ("分析工作流说明", "workflow skills"),
    "standardize": ("标准化与验证", "standardization and validation"),
    "structure3d": ("三维结构", "3D structure"),
    "transcriptomics": ("转录组学", "transcriptomics"),
    "thermodynamics": ("理化与热力学", "physicochemical and thermodynamic"),
    "variants": ("变异分析", "variant analysis"),
    "visualization": ("可视化", "visualization"),
    "workflow_security": ("工作流与数据安全", "workflow and data security"),
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
        "dnakit.datasets.diversity_sample",
        "dnakit.datasets.neural_cluster_sequences",
        "dnakit.evaluation.evaluate_frechet_distance",
        "dnakit.predictions.ensure_prediction_checkpoint",
        "dnakit.predictions.predict_enformer_benchmark",
        "dnakit.predictions.predict_custom_properties",
        "dnakit.predictions.predict_pair_properties",
        "dnakit.predictions.predict_properties",
        "dnakit.predictions.predict_sequence_properties",
        "dnakit.predictions.predict_variant_effects",
        "dnakit.predictions.train_custom_property_model",
        "dnakit.representations.ensure_model_checkpoint",
        "dnakit.representations.extract_representations",
    }
)
_EXTERNAL_PROCESS_NAMES = frozenset(
    {
        "dnakit.annotation.annotate_genes",
        "dnakit.annotation.annotate_bakta",
        "dnakit.annotation.annotate_prokka",
        "dnakit.annotation.predict_checkm2",
        "dnakit.download.build_index",
        "dnakit.thermodynamics.probe_primer3",
        "dnakit.integrations.run_analysis",
        "dnakit.integrations.list_mcp_tools",
        "dnakit.integrations.call_mcp_tool",
        "dnakit.integrations.call_mcp_workflow",
        "dnakit.sequencing_qc.run_cutadapt",
        "dnakit.sequencing_qc.run_fastp",
        "dnakit.sequencing_qc.run_fastqc",
        "dnakit.sequencing_qc.run_multiqc",
        "dnakit.sequencing_qc.filter_nanopore_reads",
        "dnakit.sequencing_qc.filter_long_reads",
        "dnakit.sequencing_qc.summarize_long_reads",
        "dnakit.sequencing_qc.trim_porechop",
        "dnakit.sequencing_qc.trim_trimmomatic",
        "dnakit.sequencing_qc.run_afterqc",
        "dnakit.assembly.assemble_shasta",
        "dnakit.assembly.evaluate_merqury",
        "dnakit.assembly.run_busco",
        "dnakit.assembly.assemble_flye",
        "dnakit.assembly.assemble_hifiasm",
        "dnakit.assembly.assemble_megahit",
        "dnakit.assembly.assemble_unicycler",
        "dnakit.ngs.map_reads",
        "dnakit.ngs.prepare_bam",
        "dnakit.ngs.mark_duplicates_picard",
        "dnakit.ngs.call_variants",
        "dnakit.ngs.compress_index_vcf_htslib",
        "dnakit.ngs.phase_variants_whatshap",
        "dnakit.ngs.sample_sequences_seqtk",
        "dnakit.ngs.sort_bam_bamtools",
        "dnakit.ngs.assemble_reads",
        "dnakit.ngs.run_dna_seq_workflow",
        "dnakit.ngs.run_rna_seq_workflow",
        "dnakit.ngs.call_chip_peaks",
        "dnakit.alignment.build_blast_database",
        "dnakit.alignment.local_blast",
        "dnakit.alignment.map_long_reads",
        "dnakit.alignment.intersect_bed",
        "dnakit.alignment.sort_bed_bedops",
        "dnakit.alignment.search_diamond_blastp",
        "dnakit.alignment.search_diamond_blastx",
        "dnakit.alignment.search_mmseqs_protein",
        "dnakit.alignment.cluster_mmseqs_protein",
        "dnakit.alignment.align_sequences",
        "dnakit.alignment.consensus_spoa",
        "dnakit.fasta_fastq.convert_seqmagick",
        "dnakit.fasta_fastq.cluster_vsearch",
        "dnakit.phylogeny.infer_ml_tree",
        "dnakit.phylogeny.infer_fast_tree",
        "dnakit.variants.call_freebayes",
        "dnakit.variants.call_cutesv",
        "dnakit.variants.call_gatk_haplotypecaller",
        "dnakit.variants.call_deepvariant",
        "dnakit.variants.call_delly",
        "dnakit.variants.call_sniffles2",
        "dnakit.variants.filter_vcf_vcftools",
        "dnakit.variants.annotate_vcf_vcfanno",
        "dnakit.variants.merge_structural_variants_survivor",
        "dnakit.variants.load_gemini_database",
        "dnakit.variants.query_gemini_database",
        "dnakit.variants.annotate_vcf_ensembl_vep",
        "dnakit.microbiome.classify_kraken2",
        "dnakit.microbiome.classify_centrifuge",
        "dnakit.microbiome.classify_kaiju",
        "dnakit.microbiome.estimate_bracken_abundance",
        "dnakit.microbiome.classify_gtdbtk",
        "dnakit.microbiome.profile_humann",
        "dnakit.epigenome.extract_methylation_methyldackel",
    }
)
_NETWORK_PREFIXES = (
    "dnakit.databases.",
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
    search_text: str | None = None

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
    if name == "dnakit.skills.install_skill":
        return "filesystem"
    if name == "dnakit.workflows.write_wes_report":
        return "filesystem"
    if name in {"dnakit.workflows.run", "dnakit.workflows.run_support"}:
        return "external_process"
    if name.startswith("dnakit.databases."):
        operation = name.rsplit(".", 1)[-1]
        if operation in {"search", "describe", "list_operations"}:
            return "local"
        if operation in {"call", "execute"}:
            return "external_process"
        from dnakit.databases._runtime import definitions

        return cast(ToolEffect, definitions().get(operation, {}).get("effect", "network"))
    if name.startswith("dnakit.encode."):
        operation = name.rsplit(".", 1)[-1]
        if operation in {
            "search_experiments",
            "get_experiment",
            "list_files",
            "search_files",
            "get_metadata",
            "get_facets",
            "get_file_info",
        }:
            return "network"
        return "filesystem"
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
    if module_name == "dnakit.generated":
        from dnakit.generated import _entries

        metadata = _entries()[export_name]
        declared_effect = metadata.get("effect", "external_process")
        effect = cast(ToolEffect, declared_effect)
        short_description = str(metadata.get("description", short_description))
        for dependency in metadata["contract"].get("dependencies", []):
            if importlib.util.find_spec(dependency) is None:
                incompatibilities.append(
                    f"Required generated-tool dependency is missing: {dependency}"
                )
    spec = AgentToolSpec(
        name=name,
        category=category,
        function_name=f"{function.__module__}.{function.__qualname__}",
        description=short_description,
        full_description=full_description,
        effect=effect,
        open_world=effect in {"network", "filesystem", "model", "external_process"},
        requires_confirmation=effect in {"filesystem", "model", "external_process"},
        agent_compatible=plan.compatible and not incompatibilities,
        incompatibilities=tuple(incompatibilities),
        hidden_parameters=tuple(
            parameter.name
            for parameter in plan.parameters
            if not parameter.supported and not parameter.required
        ),
        call_plan=plan,
        function=function,
    )
    from .descriptions import apply_native_description

    return apply_native_description(spec)


def build_tool_catalog(*, verbose: bool = True) -> tuple[AgentToolSpec, ...]:
    """Discover functions explicitly exported by stable DNAKit domain modules."""

    tools: list[AgentToolSpec] = []
    from dnakit.generated import _refresh

    _refresh()
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

_register_public_verbose(globals())
