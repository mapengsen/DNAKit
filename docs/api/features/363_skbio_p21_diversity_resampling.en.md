# scikit-bio diversity, resampling and cospeciation tests

These functions in `dnakit.microbiome.skbio_p21_functions` are adapters over [scikit-bio](https://scikit.bio/): 17 alpha-diversity metrics, UniFrac distances, count rarefaction, per-bin random subsampling, the Hommola cospeciation test and gradient ANOVA. scikit-bio is an optional dependency (`dnakit[bio-skibio]`) imported only when a function runs. All inputs and outputs are JSON-friendly lists and mappings.

## 1) `dnakit.microbiome.skbio_p21_functions.calculate_alpha_diversity()`

- **Purpose:** Compute one scikit-bio alpha-diversity metric for one sample count vector.
- **Method:** `metric` selects the function of `skbio.diversity.alpha`: margalef, mcintosh_d, mcintosh_e, menhinick, michaelis_menten_fit, sobs (observed features), osd, phydiv, pielou_e, renyi, robbins, shannon, simpson, simpson_e, singles, strong, tsallis. One function merges the 17 metrics; each metric uses only its own options (`base`, `exp`, `order`, `finite`, `num_repeats`, `taxa`, `tree_newick`, `rooted`, `weight`). michaelis_menten_fit subsamples randomly inside scikit-bio without a seed, so its value changes between calls.
- **API:** `calculate_alpha_diversity(counts, *, metric='shannon', base=None, exp=False, order=2.0, finite=False, num_repeats=1, taxa=None, tree_newick=None, rooted=None, weight=False)`.
- **Input:** A list of whole-number counts; phydiv also needs `taxa` (same length as counts) and a Newick tree.
- **Output:** `value` (a dict of observed/singles/doubles for osd), `n_features`, `total_count` and provenance.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import calculate_alpha_diversity
counts = [1, 2, 0, 5, 1, 7]
for metric in ("shannon", "simpson", "pielou_e", "osd", "sobs"):
    print(metric, calculate_alpha_diversity(counts, metric=metric, verbose=False)["value"])
print(calculate_alpha_diversity(counts, metric="shannon", base=2, verbose=False)["value"])
```

- **Example output:**

```text
shannon 1.331660286822432
simpson 0.6875
pielou_e 0.8274070571684498
osd {'observed': 5, 'singles': 2, 'doubles': 1}
sobs 5.0
1.9211796919474975
```

## 2) `dnakit.microbiome.skbio_p21_functions.calculate_unifrac_distance()`

- **Purpose:** Compute the pairwise UniFrac distance matrix of samples on a phylogenetic tree.
- **Method:** `weighted=False` gives unweighted UniFrac (`beta.unweighted_unifrac`); `weighted=True` gives weighted UniFrac (`beta.weighted_unifrac`) with `normalized` choosing the normalised form. All sample pairs are computed through `beta_diversity`.
- **API:** `calculate_unifrac_distance(counts, taxa, tree_newick, *, ids=None, weighted=False, normalized=False)`.
- **Input:** A samples x taxa count matrix (at least 2 samples), the tree tip names `taxa` for the columns, and a Newick string.
- **Output:** `ids` and the symmetric `distance_matrix`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import calculate_unifrac_distance
counts = [[1, 2, 0, 5, 1], [3, 0, 4, 1, 0], [0, 0, 1, 2, 9]]
tree = "(((a:1,b:1):1,c:2):1,(d:1,e:1):2);"
for weighted in (False, True):
    result = calculate_unifrac_distance(counts, list("abcde"), tree, weighted=weighted, normalized=True, verbose=False)
    print(result["metric"], [round(v, 4) for v in result["distance_matrix"][0]])
```

- **Example output:**

```text
unweighted_unifrac [0.0, 0.4, 0.5]
weighted_unifrac [0.0, 0.6157, 0.4352]
```

## 3) `dnakit.microbiome.skbio_p21_functions.rarefy_count_vector()`

- **Purpose:** Randomly subsample (rarefy) a count vector to a given depth.
- **Method:** Calls `skbio.stats.subsample_counts`. Sampling is without replacement by default; `replace=True` samples with replacement. A `seed` makes it reproducible.
- **API:** `rarefy_count_vector(counts, depth, *, replace=False, seed=None)`.
- **Input:** A list of whole-number counts and the target depth (not above the total when sampling without replacement).
- **Output:** The subsampled `counts`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import rarefy_count_vector
print(rarefy_count_vector([5, 3, 0, 2], 4, seed=1, verbose=False)["counts"])
print(rarefy_count_vector([5, 3, 0, 2], 4, seed=1, verbose=False)["counts"])
print(rarefy_count_vector([1, 1], 6, replace=True, seed=1, verbose=False)["counts"])
```

- **Example output:**

```text
[2, 1, 0, 1]
[2, 1, 0, 1]
[3, 3]
```

## 4) `dnakit.microbiome.skbio_p21_functions.subsample_items_per_bin()`

- **Purpose:** Randomly keep up to N items in each bin (streaming subsample).
- **Method:** Calls `skbio.stats.isubsample`: each bin keeps at most `maximum` items without replacement, and bins with fewer than `minimum` items are dropped. `bins` gives the bin of each item; without it all items share one bin.
- **API:** `subsample_items_per_bin(items, maximum, *, minimum=1, bins=None, seed=None)`.
- **Input:** A list of string items and an optional list of bin labels of the same length.
- **Output:** `selected`, a list of [bin, item] pairs.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import subsample_items_per_bin
result = subsample_items_per_bin(list("abcdefgh"), 2, bins=list("xxxyyyzz"), seed=1, verbose=False)
print(result["selected"])
print(subsample_items_per_bin(list("abc"), 2, minimum=2, bins=["x", "y", "y"], seed=1, verbose=False)["selected"])
```

- **Example output:**

```text
[['x', 'a'], ['x', 'b'], ['y', 'f'], ['y', 'd'], ['z', 'h'], ['z', 'g']]
[['y', 'c'], ['y', 'b']]
```

## 5) `dnakit.microbiome.skbio_p21_functions.run_hommola_cospeciation_test()`

- **Purpose:** The Hommola et al. (2009) host-parasite cospeciation permutation test.
- **Method:** Calls `skbio.stats.evolve.hommola_cospeciation`: it correlates host and parasite distances over every interaction edge and gets the p-value by permuting the interaction links.
- **API:** `run_hommola_cospeciation_test(host_distances, parasite_distances, interaction, *, host_ids=None, parasite_ids=None, permutations=999, seed=None)`.
- **Input:** A host distance matrix, a parasite distance matrix, and a 0/1 parasites x hosts interaction matrix.
- **Output:** `correlation`, `p_value`, `interaction_edges`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import run_hommola_cospeciation_test
host = [[0, 1, 2], [1, 0, 3], [2, 3, 0]]
parasite = [[0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 3], [1, 2, 3, 0]]
link = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]]
result = run_hommola_cospeciation_test(host, parasite, link, seed=1, verbose=False)
print(round(result["correlation"], 4), result["p_value"], result["interaction_edges"])
```

- **Example output:**

```text
0.2335 0.409 4
```

## 6) `dnakit.microbiome.skbio_p21_functions.calculate_gradient_anova()`

- **Purpose:** Compute gradient-ANOVA trajectories of ordination coordinates by metadata category.
- **Method:** `algorithm` selects the scikit-bio class: `average` (AverageGradientANOVA, RMS average), `trajectory` (TrajectoryGradientANOVA, RMS trajectory), `first_difference` (FirstDifferenceGradientANOVA) or `window_difference` (WindowDifferenceGradientANOVA, uses `window_size`). `axes=0` uses all axes.
- **API:** `calculate_gradient_anova(coordinates, proportion_explained, metadata, *, algorithm='average', trajectory_categories=None, sort_category=None, axes=3, weighted=False, window_size=2)`.
- **Input:** A mapping of sample to coordinates, the proportion explained by each axis, and a mapping of sample to {category: value} metadata.
- **Output:** `algorithm`, `weighted`, and `categories` (probability and, per group, trajectory, mean and extra info).
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import calculate_gradient_anova
coords = {f"s{i}": [0.1 * i, (i * 7 % 5) / 5, (i * 3 % 4) / 4] for i in range(1, 13)}
meta = {k: {"g": "A" if i % 2 else "B", "t": i // 2} for i, k in enumerate(coords)}
for algorithm in ("average", "trajectory", "first_difference", "window_difference"):
    result = calculate_gradient_anova(coords, [0.5, 0.3, 0.2], meta, algorithm=algorithm,
                                      sort_category="t", trajectory_categories=["g"], verbose=False)
    groups = result["categories"][0]["groups"]
    print(result["algorithm"], {k: round(v["mean"], 4) for k, v in groups.items()})
```

- **Example output:**

```text
avg {'A': 0.1837, 'B': 0.1878}
trajectory {'A': 0.1786, 'B': 0.1786}
diff {'A': 0.0312, 'B': 0.0}
wdiff {'A': 0.0375, 'B': 0.0}
```
