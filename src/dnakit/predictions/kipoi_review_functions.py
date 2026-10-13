"""Review group G18d for kipoi: no new wrappers in this batch.

Every reviewed kipoi row is plumbing, registry or CLI code, or a source listing that needs
the network. ``kipoi.list_sources`` clones each git-backed model source on first use, so it
is marked ``needs-network`` and is never run live. The DNA-free model listing and description
calls already live in ``dnakit.predictions.dna_model_backends`` as ``list_kipoi_models`` and
``describe_kipoi_model``.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
