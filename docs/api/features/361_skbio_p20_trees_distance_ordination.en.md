# scikit-bio Phylogenetic Trees, Distances and Ordination

This page covers the P20 tree and distance-statistics interfaces in `dnakit.phylogeny.skbio_p20_tree_functions` and `dnakit.microbiome.skbio_p20_stats_functions`. They use the optional `bio-skibio` extra (scikit-bio 0.7.4). Inputs and outputs use JSON-friendly Newick strings, lists and dictionaries; scikit-bio is imported only when a function runs.

## 1) `reroot_tree()`

- **Purpose:** Reroot a tree at a node or on the branch above it.
- **Method:** Uses `TreeNode.root_at`, `unrooted_copy` or `unrooted_move`. Select a node by name or by a list of tip names.
- **API:** `dnakit.phylogeny.skbio_p20_tree_functions.reroot_tree(newick, node, *, above=False, method='root_at', verbose=True)`.
- **Input:** A Newick string and a node name or list of tip names.
- **Output:** The rerooted Newick string, selected method and tip count.
- **Example code:**

```python
from dnakit.phylogeny.skbio_p20_tree_functions import reroot_tree, tree_diameter

tree = "((a:1,b:2)90:1,(c:3,d:1)0.5:2,e:1);"
print(reroot_tree(tree, "c", verbose=False)["newick"])
print(tree_diameter(tree, verbose=False)["distance"])
```

- **Example output:**

```text
((d:1.0,((a:1.0,b:2.0)90:1.0,e:1.0):2.0)0.5:3.0)c;
8.0
```

## 2) Distance tests and ordination

- **Purpose:** Run ANOSIM, PERMDISP and Mantel tests, project PCoA biplots, and perform constrained ordination.
- **Method:** `anosim_test`, `permdisp_test`, `mantel_test` and `pcoa_biplot_projection` call the corresponding scikit-bio algorithms. `constrained_ordination` supports `rda` and `cca`; `plot_ordination` saves an image and returns its coordinates.
- **API:** For example, `dnakit.microbiome.skbio_p20_stats_functions.anosim_test(distance_matrix, grouping, *, ids=None, permutations=999, seed=None, verbose=True)`.
- **Input:** A symmetric distance matrix and group labels; ordination functions may also accept sample/feature tables.
- **Output:** Test statistics, p-values, ordination coordinates, or an image path with coordinate data.
- **Example code:**

```python
from dnakit.microbiome.skbio_p20_stats_functions import anosim_test

result = anosim_test(
    [[0, 1, 4, 5], [1, 0, 5, 4], [4, 5, 0, 1], [5, 4, 1, 0]],
    ["A", "A", "B", "B"],
    permutations=99,
    seed=7,
    verbose=False,
)
print(result["test_statistic_name"], result["number_of_groups"])
```

- **Example output:**

```text
R 2
```
