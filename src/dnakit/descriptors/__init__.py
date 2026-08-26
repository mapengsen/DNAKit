"""Deterministic, dependency-free DNA sequence descriptors."""

from dnakit.descriptors._dinucleotide import (
    DEFAULT_MAX_DINUCLEOTIDE_TABLE_BYTES,
    DINUCLEOTIDE_PROPERTY_SPECS,
    DINUCLEOTIDE_TABLE_SCHEMA_VERSION,
    DINUCLEOTIDES,
    DinucleotideProperty,
    DinucleotidePropertySpec,
    DinucleotidePropertyTable,
    load_dinucleotide_property_table,
)
from dnakit.descriptors._shared import DescriptorAmbiguityPolicy, SequenceInput
from dnakit.descriptors.all import all_descriptors
from dnakit.descriptors.basic import (
    base_composition,
    base_skew,
    cpg_features,
    gc_at_content,
    length_features,
)
from dnakit.descriptors.codon import codon_statistics
from dnakit.descriptors.complexity import exact_repeat_fraction, linguistic_complexity
from dnakit.descriptors.entropy import shannon_entropy
from dnakit.descriptors.homopolymer import homopolymer_runs
from dnakit.descriptors.kmer import canonical_kmer, kmer_statistics
from dnakit.descriptors.results import (
    AllDescriptorsResult,
    CodonResult,
    ComplexityResult,
    CompositionResult,
    ContentResult,
    CpGResult,
    DescriptorResult,
    EntropyResult,
    ExactRepeatResult,
    HomopolymerResult,
    HomopolymerRun,
    KmerResult,
    LengthResult,
    RepeatRun,
    SkewResult,
    WindowDescriptorResult,
    WindowResult,
)
from dnakit.descriptors.schema import (
    DESCRIPTOR_NAMES_V1,
    DESCRIPTOR_SCHEMA_V1,
    DESCRIPTOR_SCHEMA_VERSION,
    DescriptorField,
    descriptor_schema_v1,
)
from dnakit.descriptors.window import window_descriptors

__all__ = [
    "DEFAULT_MAX_DINUCLEOTIDE_TABLE_BYTES",
    "DESCRIPTOR_NAMES_V1",
    "DESCRIPTOR_SCHEMA_V1",
    "DESCRIPTOR_SCHEMA_VERSION",
    "DINUCLEOTIDES",
    "DINUCLEOTIDE_PROPERTY_SPECS",
    "DINUCLEOTIDE_TABLE_SCHEMA_VERSION",
    "AllDescriptorsResult",
    "CodonResult",
    "ComplexityResult",
    "CompositionResult",
    "ContentResult",
    "CpGResult",
    "DescriptorAmbiguityPolicy",
    "DescriptorField",
    "DescriptorResult",
    "DinucleotideProperty",
    "DinucleotidePropertySpec",
    "DinucleotidePropertyTable",
    "EntropyResult",
    "ExactRepeatResult",
    "HomopolymerResult",
    "HomopolymerRun",
    "KmerResult",
    "LengthResult",
    "RepeatRun",
    "SequenceInput",
    "SkewResult",
    "WindowDescriptorResult",
    "WindowResult",
    "all_descriptors",
    "base_composition",
    "base_skew",
    "canonical_kmer",
    "codon_statistics",
    "cpg_features",
    "descriptor_schema_v1",
    "exact_repeat_fraction",
    "gc_at_content",
    "homopolymer_runs",
    "kmer_statistics",
    "length_features",
    "linguistic_complexity",
    "load_dinucleotide_property_table",
    "shannon_entropy",
    "window_descriptors",
]
