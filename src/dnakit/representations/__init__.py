"""DNA foundation-model checkpoint and representation APIs."""

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.representations.backends import RepresentationBackend
from dnakit.representations.checkpoints import (
    CheckpointInfo,
    default_checkpoint_root,
    ensure_model_checkpoint,
)
from dnakit.representations.extraction import extract_representations
from dnakit.representations.models import (
    DNAEmbeddingModel,
    RepresentationConfig,
    available_embedding_models,
    get_embedding_model,
)
from dnakit.representations.results import RepresentationResult

__all__ = [
    "CheckpointInfo",
    "DNAEmbeddingModel",
    "RepresentationBackend",
    "RepresentationConfig",
    "RepresentationResult",
    "available_embedding_models",
    "default_checkpoint_root",
    "ensure_model_checkpoint",
    "extract_representations",
    "get_embedding_model",
]

# Public integration of the reviewed P00–P45 operation packs.
from .lucaonetasks_p16_functions import (
    embed_dna_with_dnabert,
)

__all__ = [
    *__all__,
    "embed_dna_with_dnabert",
]

_register_public_verbose(globals())
