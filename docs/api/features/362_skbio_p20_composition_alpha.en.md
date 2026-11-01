# scikit-bio Compositions and Alpha Diversity

This page covers P20 composition interfaces in `dnakit.microbiome.skbio_p20_composition_functions` and `dnakit.microbiome.skbio_p20_stats_functions`. They use the optional `bio-skibio` extra (scikit-bio 0.7.4). Composition matrices have samples in rows and features in columns; values are not silently altered.

## 1) Composition transforms

- **Purpose:** Apply closure, zero replacement, perturbation, power and log-ratio transforms, and run structural-zero or differential-abundance analyses.
- **Method:** Wraps scikit-bio composition algorithms and validates matrix dimensions, finite values and parameter ranges.
- **API:** `dnakit.microbiome.skbio_p20_composition_functions.closure_compositions(matrix, *, verbose=True)`; zero replacement uses `replace_zeros_multiplicative(matrix, *, delta=None, verbose=True)`.
- **Input:** A finite nonnegative composition matrix; some models also need metadata or grouping labels.
- **Output:** JSON-friendly transformed matrices, statistics and provenance.
- **Example code:**

```python
from dnakit.microbiome.skbio_p20_composition_functions import (
    closure_compositions,
    replace_zeros_multiplicative,
)

print(closure_compositions([[2, 3, 5], [1, 1, 2]], verbose=False)["compositions"])
print(replace_zeros_multiplicative([[0, 2, 3]], verbose=False)["compositions"])
```

- **Example output:**

```text
[[0.2, 0.3, 0.5], [0.25, 0.25, 0.5]]
[0.1111111111111111, 0.35555555555555557, 0.5333333333333333]
```

## 2) Alpha diversity

- **Purpose:** Calculate richness, evenness and phylogenetic alpha-diversity metrics.
- **Method:** `alpha_diversity_metric` selects a scikit-bio metric by `metric`; Faith's PD can use a feature tree.
- **API:** `dnakit.microbiome.skbio_p20_stats_functions.alpha_diversity_metric(counts, metric, *, sample_ids=None, feature_ids=None, tree_newick=None, verbose=True)`.
- **Input:** Nonnegative integer counts per sample and a metric name. Faith's PD also requires a Newick tree and feature IDs.
- **Output:** One diversity value per sample with provenance.
- **Example code:**

```python
from dnakit.microbiome.skbio_p20_stats_functions import alpha_diversity_metric

result = alpha_diversity_metric(
    [[1, 2, 0], [0, 1, 3]], "dominance",
    sample_ids=["s1", "s2"], verbose=False,
)
print(result["values"])
```

- **Example output:**

```text
[0.5555555555555556, 0.625]
```
