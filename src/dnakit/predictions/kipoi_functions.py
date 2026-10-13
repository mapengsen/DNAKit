"""Kipoi wrappers for this batch: none. Every public kipoi function is recorded in the ledger.

Kipoi's public functions cover model-zoo listing and loading from git-backed sources (cloned
on first use), CLI and conda or singularity helpers, plugin registries, dataloader and tensor
utilities, and prediction writers. None of them runs on local DNA sequences or intervals
without a model zoo. The model listing and model description calls are already wrapped in
``dnakit.predictions.dna_model_backends`` as ``list_kipoi_models`` and ``describe_kipoi_model``.
"""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose

__all__: list[str] = []

_register_public_verbose(globals(), tuple(__all__))
