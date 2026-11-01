# NUPACK and DNA blocks

## Purpose

Explicit DNA operations for integration pack P39/P40, with optional upstream backends and JSON results.

## Method

Each function invokes the named upstream operation. Local paths are explicit; new outputs must not already exist. UCSC/blacklist retrieval functions use the network only when called. NUPACK and BLAST+ are not present in the validation environment; those scientific pipelines remain unvalidated. Enformer checkpoints must be local and match the supplied architecture and adapter configuration.

## API

`dnakit.secondary_structure.nupack_p39_functions`

| Function | Operation |
| --- | --- |
| `nupack_partition_function` | Complex ensemble partition function and free energy (kcal/mol); Decimal partition values are strings. |
| `nupack_structure_energy` | Free energy in kcal/mol of a specified dot-parens-plus secondary structure. |
| `nupack_sample_structures` | Draw Boltzmann secondary structures using NUPACK's native sampling (num_sample). |
| `nupack_pair_probabilities` | Equilibrium base-pair probability matrix; diagonal entries denote unpaired probabilities. |
| `nupack_mfe_structures` | MFE proxy secondary structures with proxy and stacking-state free energies. |
| `nupack_suboptimal_structures` | Proxy structures within energy_gap kcal/mol of the minimum stacking-state energy. |
| `nupack_ensemble_size` | Number of stacking states, or secondary structures for ensemble='nostacking'. |
| `nupack_design_complex` | Design one complex to a target structure, optionally with IUPAC sequence constraints. |
| `nupack_sequence_distance` | NUPACK sequence distance, including its native run-length notation and nick handling. |
| `nupack_structure_distance` | NUPACK secondary-structure distance, accepting its native structure notations. |
| `nupack_tube_analysis` | Analyze equilibrium complex concentrations of a DNA test tube; optionally save native analysis for reuse. |
| `nupack_complex_concentrations` | Recompute concentrations from a trusted native NUPACK AnalysisResult file without recomputing partition functions. |
| `nupack_design_tubes` | Design multiple target tubes from IUPAC strand constraints, target structures and molar concentrations. |
| `nupack_evaluate_design` | Evaluate fully specified designed DNA strands under the requested physical model and tube targets. |

`dnakit.molecular_lab.geneblocks_p39_functions`

| Function | Operation |
| --- | --- |
| `geneblocks_find_common` | Discover common DNA blocks shared across sequences using upstream all-versus-all BLAST+. |
| `geneblocks_unique_segments` | Coordinates of non-shared segments using CommonBlocks.compute_unique_blocks. |
| `geneblocks_export_csv` | Write the upstream semicolon-separated common-block table to a new file. |
| `geneblocks_common_records` | Export common blocks as JSON representations of upstream sequence records. |
| `geneblocks_unique_records` | Export non-shared segments as upstream sequence records including their DNA. |
| `geneblocks_annotate_common` | Annotate original DNA sequences with discovered common-block feature coordinates. |
| `geneblocks_transfer_annotations` | Transfer strand-aware features between homologous common blocks of annotated DNA records. |
| `geneblocks_plot_common` | Plot shared blocks along DNA sequences; return the image path and underlying block data. |
| `geneblocks_merge_differences` | Merge nearby edit blocks into larger changes using DiffBlocks.merged. |
| `geneblocks_difference_features` | Represent DNA differences as annotated features relative to s1 or s2. |
| `geneblocks_plot_differences` | Plot DNA edit blocks and return the image path and block coordinates. |
| `geneblocks_reconstruct_sequences` | Reconstruct both sequences from a complete user-supplied list of edit and equal blocks. |
| `geneblocks_sequence_record` | Wrap DNA as a sequence record using geneblocks.biotools.sequence_to_record. |
| `geneblocks_circular_equality` | Test whether DNA strings are rotations of one circular sequence using upstream semantics. |
| `geneblocks_extract_subsequence` | Extract a strand-aware 0-based half-open DNA subsequence with geneblocks.Location. |

## Input

DNA strings, JSON mappings/lists, local genomic files or explicit local model checkpoints. Interval coordinates and physical-model parameters follow each function docstring; the NUPACK default material is DNA. All functions accept keyword-only `verbose`.

## Output

A dictionary with completion status, backend provenance and operation-specific data. Writers return output paths. Figures return image paths and underlying data. Optional backends raise `BackendUnavailableError` when missing. P40 manifests distinguish unvalidated pipelines from tested interfaces.

**Example code:**

```python
from dnakit.secondary_structure.nupack_p39_functions import nupack_sequence_distance
from dnakit.exceptions import BackendUnavailableError
try:
    print(nupack_sequence_distance('ACGT', 'ACGA', verbose=False)['distance'])
except BackendUnavailableError:
    print('NUPACK backend unavailable')
```

**Example output:**

```text
NUPACK backend unavailable
```

The example shows the actual missing-backend outcome in this environment. Scientific NUPACK results are not claimed. Upstream API: [NUPACK utilities](https://docs.nupack.org/utilities/).

## NUPACK loop energies (corrected P39 coverage)

**Purpose:** Compute loop free energy and individual stacking-state energies.

**Method:** Call the documented public `Model.loop_energy` and `Model.stack_energies` methods: https://docs.nupack.org/model/ . The earlier internal-helper classification was incorrect.

**API:**

```text
nupack_loop_energy(loop, structure, *, model=None, verbose=True)
nupack_loop_stacking_energies(loop, structure, *, model=None, verbose=True)
```

**Input:** A loop sequence and NUPACK loop structure; `+` denotes a strand nick. `model` accepts the same physical options as the other wrappers and defaults to DNA.

**Output:** `free_energy` in kcal/mol, or `stacking_energies` mapping native state strings to kcal/mol. NUPACK is absent from this validation environment; numerical equivalence has not been tested.

**Example code:**

```python
from dnakit.secondary_structure.nupack_p39_functions import nupack_loop_energy
from dnakit.exceptions import BackendUnavailableError
try:
    print(nupack_loop_energy('AATT', '....', verbose=False)['free_energy'])
except BackendUnavailableError:
    print('NUPACK unavailable; energy calculation not run.')
```

**Example output:**

```text
NUPACK unavailable; energy calculation not run.
```
