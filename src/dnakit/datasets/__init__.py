"""Dataset deduplication, clustering, splitting, leakage, and quality evaluation."""

from dnakit._processing.dispatch import register as _register_inputs
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.datasets.advanced_deduplicate import deduplicate_approximate, deduplicate_iupac
from dnakit.datasets.advanced_split import joint_split, temporal_split
from dnakit.datasets.clustering import (
    cluster_sequences,
    hierarchical_cluster,
    select_representatives,
)
from dnakit.datasets.config import (
    ClusterConfig,
    DeduplicationConfig,
    DiversitySamplingConfig,
    HierarchicalClusteringConfig,
    IUPACDeduplicationConfig,
    JointSplitConfig,
    LeakageConfig,
    NeuralClusteringConfig,
    SplitConfig,
    TemporalSplitConfig,
)
from dnakit.datasets.deduplicate import deduplicate
from dnakit.datasets.evaluation import detect_leakage, evaluate_split_quality
from dnakit.datasets.file_results import FileRecords, FileSplitResult, FileSplitSubset
from dnakit.datasets.filter import (
    MetadataExclusionField,
    MetadataValues,
    exclude_by_metadata,
    exclude_chromosomes,
    exclude_species,
)
from dnakit.datasets.neural_clustering import neural_cluster_sequences
from dnakit.datasets.results import (
    Cluster,
    ClusteringResult,
    DedupGroup,
    DeduplicationResult,
    DiversitySamplingResult,
    HierarchicalClusteringResult,
    IUPACDeduplicationResult,
    IUPACDuplicateGroup,
    IUPACPairRelation,
    JointSplitResult,
    LeakageEvent,
    LeakageReport,
    LinkageStep,
    NeuralClusteringResult,
    RepresentativeSelectionResult,
    SamplingAllocation,
    SplitAssignment,
    SplitQualityResult,
    SplitResult,
    SplitSubset,
    TemporalSplitResult,
)
from dnakit.datasets.sampling import diversity_sample
from dnakit.datasets.split import split
from dnakit.datasets.tabular import export_grouped_split, export_kfold_table

_register_inputs(
    globals(),
    (
        "deduplicate",
        "deduplicate_iupac",
        "deduplicate_approximate",
        "cluster_sequences",
        "hierarchical_cluster",
        "neural_cluster_sequences",
        "diversity_sample",
        "select_representatives",
        "split",
        "temporal_split",
        "joint_split",
    ),
    mode="collection",
)
_register_inputs(
    globals(), ("exclude_by_metadata", "exclude_chromosomes", "exclude_species"), mode="set"
)
_register_inputs(globals(), ("detect_leakage",), mode="mapping")
_register_inputs(globals(), ("evaluate_split_quality",), mode="set")

__all__ = [
    "Cluster",
    "ClusterConfig",
    "ClusteringResult",
    "DedupGroup",
    "DeduplicationConfig",
    "DeduplicationResult",
    "DiversitySamplingConfig",
    "DiversitySamplingResult",
    "FileRecords",
    "FileSplitResult",
    "FileSplitSubset",
    "HierarchicalClusteringConfig",
    "HierarchicalClusteringResult",
    "IUPACDeduplicationConfig",
    "IUPACDeduplicationResult",
    "IUPACDuplicateGroup",
    "IUPACPairRelation",
    "JointSplitConfig",
    "JointSplitResult",
    "LeakageConfig",
    "LeakageEvent",
    "LeakageReport",
    "LinkageStep",
    "MetadataExclusionField",
    "MetadataValues",
    "NeuralClusteringConfig",
    "NeuralClusteringResult",
    "RepresentativeSelectionResult",
    "SamplingAllocation",
    "SplitAssignment",
    "SplitConfig",
    "SplitQualityResult",
    "SplitResult",
    "SplitSubset",
    "TemporalSplitConfig",
    "TemporalSplitResult",
    "cluster_sequences",
    "deduplicate",
    "deduplicate_approximate",
    "deduplicate_iupac",
    "detect_leakage",
    "diversity_sample",
    "evaluate_split_quality",
    "exclude_by_metadata",
    "exclude_chromosomes",
    "exclude_species",
    "export_grouped_split",
    "export_kfold_table",
    "hierarchical_cluster",
    "joint_split",
    "neural_cluster_sequences",
    "select_representatives",
    "split",
    "temporal_split",
]

_register_public_verbose(globals())
