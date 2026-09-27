# Statistical and differential analysis

## 1) Test group differences using a distance matrix {#bio-35}

- **Purpose:** Compare sample groups using PERMANOVA, ANOSIM or PERMDISP.
- **Method:** Validate the distance matrix, then permute group labels to calculate a p-value.
- **API:** `dnakit.bioinformatics.statistics.distance_matrix(values)`; `dnakit.bioinformatics.statistics.distance_test(distances, groups, method='permanova', permutations=999, seed=0)`
- **Input:** A symmetric nonnegative distance matrix and matching group labels; PERMDISP requires Euclidean distances; requires `bio-analysis`.
- **Output:** `distance_matrix` returns a distance matrix; `distance_test` returns a statistic, p-value, and method.
- **Example code:**

```python
from dnakit.bioinformatics.statistics import distance_matrix, distance_test

distances = distance_matrix([[0, 1, 4, 5], [1, 0, 3, 4], [4, 3, 0, 1], [5, 4, 1, 0]])
result = distance_test(distances, ["a", "a", "b", "b"], permutations=99, seed=0)
print(result.method, result.statistic, result.pvalue)
```

- **Example output:**

```text
permanova 32.0 0.36
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#statistics).

## 2) Compare distance matrices and identify associated environmental variables {#bio-36}

- **Purpose:** Test matrix association with Mantel and select environmental-variable subsets with BIOENV.
- **Method:** Mantel permutes correlations of upper triangles; BIOENV enumerates small variable subsets.
- **API:** `dnakit.bioinformatics.statistics.mantel(left, right, method='pearson', permutations=999, seed=0)`; `dnakit.bioinformatics.statistics.bioenv(distances, environmental, max_columns=12)`
- **Input:** Distance matrices with matching sample order, or distances plus a samples×environmental-variables matrix; requires `bio-analysis`.
- **Output:** `mantel` returns correlation statistics and a p-value; variable selection returns associated variables and scores.
- **Example code:**

```python
from dnakit.bioinformatics.statistics import mantel, bioenv

distances = [[0, 1, 3], [1, 0, 2], [3, 2, 0]]
test = mantel(distances, distances, permutations=99, seed=0)
selection = bioenv(distances, [[0, 0], [1, 2], [3, 1]])
print(round(test.statistic, 6))
print(selection[0])
```

- **Example output:**

```text
1.0
{'columns': (0,), 'correlation': 1.0}
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#statistics).

## 3) Run PCoA, CA, CCA or RDA ordination {#bio-37}

- **Purpose:** Convert sample relationships into low-dimensional coordinates.
- **Method:** PCoA decomposes a distance matrix; CA/CCA/RDA perform unconstrained or constrained ordination of data.
- **API:** `dnakit.bioinformatics.statistics.pcoa(distances, dimensions=None)`; `dnakit.bioinformatics.statistics.ordination(values, method='ca', constraints=None)`
- **Input:** PCoA takes distances; CA/CCA take nonnegative data; CCA/RDA also need a constraint matrix; requires `bio-analysis`.
- **Output:** Returns `OrdinationResult` with sample coordinates, explained variance, and axis information.
- **Example code:**

```python
from dnakit.bioinformatics.statistics import pcoa, ordination

coordinates = pcoa([[0, 1, 3], [1, 0, 2], [3, 2, 0]], dimensions=1)
ca = ordination([[1, 2, 3], [3, 2, 1], [2, 4, 1]], method="ca")
print(coordinates.samples.shape, ca.samples.shape)
```

- **Example output:**

```text
(3, 1) (3, 2)
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#statistics).

## 4) Apply CLR, ILR or ALR composition transforms {#bio-38}

- **Purpose:** Transform positive compositions into log-ratio space and back.
- **Method:** CLR centers rowwise logarithms; ILR uses an orthonormal Helmert basis; ALR uses a reference column.
- **API:** `dnakit.bioinformatics.statistics.composition_transform(values, method='clr', reference=0, inverse=False)`
- **Input:** A strictly positive samples×features matrix; handle zeros explicitly. Inverse transforms return unit-sum compositions; requires `bio-analysis`.
- **Output:** Returns a transformed numeric matrix; inverse transformation restores a matrix.
- **Example code:**

```python
from dnakit.bioinformatics.statistics import composition_transform

encoded = composition_transform([[1, 2, 1]], method="clr")
recovered = composition_transform(encoded, method="clr", inverse=True)
print(encoded.round(6).tolist())
print(recovered.round(6).tolist())
```

- **Example output:**

```text
[[-0.231049, 0.462098, -0.231049]]
[[0.25, 0.5, 0.25]]
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#statistics).

## 5) Compare compositions with ANCOM {#bio-39}

- **Purpose:** Count significant pairwise log-ratio differences per taxon.
- **Method:** Test pairwise log-ratios, apply per-taxon Holm correction and return W.
- **API:** `dnakit.bioinformatics.statistics.ancom(values, groups, alpha=0.05, pseudocount=0.0)`
- **Input:** A nonnegative abundance matrix and group labels; explicitly set a pseudocount for zeros; requires `bio-analysis`. Choose your own W decision threshold.
- **Output:** Returns a mapping with ANCOM `W` statistics, p-values, and significance results.
- **Example code:**

```python
from dnakit.bioinformatics.statistics import ancom

result = ancom([[1, 2, 3], [2, 3, 4], [8, 2, 3], [9, 3, 4]],
               ["a", "a", "b", "b"])
print(result["W"].tolist())
print(result["pvalues"].shape)
```

- **Example output:**

```text
[0, 0, 0]
(3, 3)
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#statistics).

## 6) Calculate alpha and beta diversity {#bio-40}

- **Purpose:** Measure within-sample diversity or between-sample compositional distances.
- **Method:** Alpha supports observed/Shannon/Simpson/inverse Simpson; beta supports Bray–Curtis/Jaccard/Euclidean and others.
- **API:** `dnakit.bioinformatics.statistics.alpha_diversity(counts, metric='shannon')`; `dnakit.bioinformatics.statistics.beta_diversity(counts, metric='braycurtis')`
- **Input:** A nonnegative samples×taxa count matrix; requires `bio-analysis`.
- **Output:** `alpha_diversity` returns a diversity vector per sample; `beta_diversity` returns a between-sample distance matrix.
- **Example code:**

```python
from dnakit.bioinformatics.statistics import alpha_diversity, beta_diversity

counts = [[1, 1], [2, 0]]
print(alpha_diversity(counts).round(6).tolist())
print(beta_diversity(counts).round(6).tolist())
```

- **Example output:**

```text
[0.693147, -0.0]
[[0.0, 0.5], [0.5, 0.0]]
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#statistics).

## 7) Calculate phylogenetic UniFrac distance {#bio-41}

- **Purpose:** Compare community compositions using a phylogenetic tree.
- **Method:** Compute unweighted UniFrac from branch presence or weighted UniFrac from abundances.
- **API:** `dnakit.bioinformatics.phylogeny.unifrac(counts, taxa, tree, weighted=False, normalized=True)`
- **Input:** Nonnegative abundances, taxa columns exactly matching tree tips, and a rooted tree with nonnegative branch lengths; requires `bio-analysis`.
- **Output:** Returns a between-sample UniFrac distance matrix.
- **Example code:**

```python
from dnakit.bioinformatics.phylogeny import TreeNode, unifrac

tree = TreeNode(children=(TreeNode("a", 1), TreeNode("b", 1)))
print(unifrac([[1, 0], [0, 1]], ["a", "b"], tree).tolist())
```

- **Example output:**

```text
[[0.0, 1.0], [1.0, 0.0]]
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#statistics).

## 8) Estimate abundance differences with ANCOM-BC {#bio-42}

- **Purpose:** Estimate abundance effects, standard errors and p/q values for design variables.
- **Method:** Use OLS, robust sandwich variance and three-Gaussian EM sampling-bias correction.
- **API:** `dnakit.bioinformatics.counts.ancom_bc(counts, design, tolerance=1e-5, max_iterations=1000)`
- **Input:** Positive samples×taxa values with at least 10 taxa, and a full-rank design containing an intercept; requires `bio-analysis`. The example is synthetic; see [inference scope](../../faq/bioinformatics.md#statistics).
- **Output:** Returns `AbundanceResult` with coefficients, p-values, and adjusted q-values.
- **Example code:**

```python
import numpy as np
from dnakit.bioinformatics.counts import ancom_bc

rng = np.random.default_rng(13)
design = np.column_stack([np.ones(30), np.repeat([0, 1], 15)])
effects = np.r_[np.full(6, 1.2), np.full(6, -1.0), np.zeros(48)]
abundance = np.exp(rng.normal(4, 0.7, (30, 60)) + design[:, 1, None] * effects)
result = ancom_bc(abundance, design)
print(result.coefficients.shape, result.qvalues.shape)
```

- **Example output:**

```text
(60, 2) (60, 2)
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#statistics).

## 9) Analyze expression differences with NB2 GLM {#bio-43}

- **Purpose:** Estimate log2 fold change, standard errors and p/q values for a contrast.
- **Method:** Fit an unshrunk NB2 GLM, apply Wald tests and BH correction.
- **API:** `dnakit.bioinformatics.counts.differential_expression(counts, design, contrast, size_factors=None, dispersions=None)`
- **Input:** Integer samples×genes counts, a numeric design, a contrast and at least two residual degrees of freedom; real analyses need biological replicates. Requires `bio-analysis`; see [differences from DESeq2/edgeR](../../faq/bioinformatics.md#statistics).
- **Output:** Returns `DifferentialExpressionResult` with log2 fold changes, p-values, and convergence status.
- **Example code:**

```python
import numpy as np
from dnakit.bioinformatics.counts import differential_expression

rng = np.random.default_rng(18)
design = np.column_stack([np.ones(30), np.repeat([0, 1], 15)])
mean = np.exp(design @ np.array([[3, 4, 5], [0.8, -1.0, 0.0]]))
counts = rng.negative_binomial(10, 10 / (mean + 10))
result = differential_expression(counts, design, [0, 1],
                                 size_factors=np.ones(30), dispersions=0.1)
print(result.log2_fold_change.round(3).tolist())
print(result.converged.tolist())
```

- **Example output:**

```text
[1.319, -1.501, -0.011]
[True, True, True]
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#statistics).

## 10) Adjust p-values for multiple testing {#bio-44}

- **Purpose:** Convert p-values into Benjamini–Hochberg q-values.
- **Method:** Adjust sorted p-values and restore the original order.
- **API:** `dnakit.bioinformatics.counts.adjust_pvalues(values)`
- **Input:** A one-dimensional array of p-values; requires `bio-analysis`.
- **Output:** Returns an array of adjusted p-values.
- **Example code:**

```python
from dnakit.bioinformatics.counts import adjust_pvalues

print(adjust_pvalues([0.01, 0.04, 0.03]).tolist())
```

- **Example output:**

```text
[0.03, 0.04, 0.04]
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#statistics).
