"""Optional genomic prediction, gene finding, methylation and DNA engineering APIs."""

from dnakit._verbosity import register_public as _register_public_verbose

from .ani import average_nucleotide_identity
from .breeding import genomic_prediction, genomic_relationship_matrix
from .genes import GenePredictionResult, PredictedGene, predict_prokaryotic_genes
from .methylation import read_bedmethyl, summarize_methylation
from .nanotechnology import design_dna_duplex, export_oxdna_design
from .population import genotype_pca
from .storage import decode_dna_storage, encode_dna_storage

__all__ = [
    "GenePredictionResult",
    "PredictedGene",
    "average_nucleotide_identity",
    "decode_dna_storage",
    "design_dna_duplex",
    "encode_dna_storage",
    "export_oxdna_design",
    "genomic_prediction",
    "genomic_relationship_matrix",
    "genotype_pca",
    "predict_prokaryotic_genes",
    "read_bedmethyl",
    "summarize_methylation",
]

_register_public_verbose(globals(), tuple(__all__))
