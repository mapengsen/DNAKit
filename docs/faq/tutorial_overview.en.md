# Feature catalog and counts

[Back to the feature tree](../api/features/function_tree.md)

> Updated: 2026-09-13. Checked against public exports in the current workspace's `src/dnakit`, model/encoding catalogs, CLI, Agent catalog, and feature guides, including uncommitted implementations. Capabilities are grouped into data preparation, processing, analysis, evaluation, visualization and large language model tools; engineering capabilities are covered in the FAQ. IDs follow feature guides or requirement tracking and are not consecutive feature counts. This overview groups user capabilities rather than enumerating configuration classes, result types, compatibility aliases, or low-level helpers. Listed capabilities have interfaces and implementations; this does not imply all optional dependencies are installed or all external backends have been tested locally.

## Audit scope and counts


This audit covered public module exports, public functions/class methods in each `bioinformatics` submodule, `tools/catalog.py`, CLI commands, and Galaxy wrappers. The feature tree groups capabilities by task without counting them again under their Python modules.

| Code catalog | Verified count | Meaning |
| --- | --- | --- |
| `descriptors.DESCRIPTOR_NAMES_V1` | 240 | Fixed-schema output fields, including statistics and property-table fields; not 240 independent APIs |
| `bioinformatics.encoding.DNA_ENCODINGS` | 44 | Encoding names; output dimensions depend on parameters/sequences |
| `representations.available_embedding_models()` | 11 | All provide full fine-tuning and linear probing; model-specific runtimes and weights are required |
| `predictions.available_machine_learning_models()` | 37 | Algorithm entries with separate classification/regression support |
| `predictions.available_prediction_tasks()` | 54 | Registered direct prediction tasks, excluding user-trained properties |

Formulas, input requirements, and limits are defined in the [bioinformatics notes](bioinformatics.md) and code docstrings. For example, gene naming uses reference alignment; NB2 differential expression omits the full DESeq2/edgeR shrinkage pipelines; constraint optimization uses a bounded search. Entries reflect these implementation boundaries.
