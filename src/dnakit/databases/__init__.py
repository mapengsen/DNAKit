"""DNAKit scientific database APIs with local, attributed provider implementations."""

from ._generated import *  # noqa: F403
from ._generated import __all__ as _functions
from ._routes import (
    cadd_get_scores,
    orphanet_get_disease_genes,
    pdb_get_structure,
    string_get_functional_enrichment,
)
from ._runtime import call, describe, execute, list_operations, search

__all__ = [
    "call",
    "describe",
    "search",
    "execute",
    "list_operations",
    *_functions,
    "pdb_get_structure",
    "orphanet_get_disease_genes",
    "cadd_get_scores",
    "string_get_functional_enrichment",
]

from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
