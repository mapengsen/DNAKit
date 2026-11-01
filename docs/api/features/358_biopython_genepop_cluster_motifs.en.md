# Biopython GenePop editing, clustering and motif scoring

This document covers three sets of wrappers over [Biopython](https://biopython.org) 1.88: GenePop file editing in `dnakit.population` (`Bio.PopGen.GenePop`), the clustering routines of `Bio.Cluster` (the C Clustering Library) in `dnakit.population`, and `Bio.motifs` scoring in `dnakit.regulation`. Matrices are lists of rows; distances are named `euclidean`, `cityblock`, `pearson`, `abs_pearson`, `uncentered`, `abs_uncentered`, `spearman` and `kendall`. The random starts of k-means, k-medoids and the self-organizing map come from the C library and cannot be seeded; `initial_clusters` makes them deterministic. Install Biopython with `python -m pip install "dnakit[bio-formats]"`; it is imported only when a function runs. The WebLogo example needs network access. The examples import from the modules `dnakit.population.biopython_p19_cluster_functions` and `dnakit.regulation.biopython_p19_motif_functions`.

## Example data

`pops.gen` is a small GenePop file with three populations and three loci (two-digit alleles). The clustering and motif examples pass their data directly: five items with three measurements, and seven DNA binding sites of length 5.

**Example file:** `pops.gen`

```text
Test genepop
Loc1
Loc2
Loc3
Pop
Ind1, 0101 0202 0303
Ind2, 0102 0203 0000
Pop
Ind3, 0201 0202 0303
Ind4, 0101 0102 0103
Pop
Ind5, 0202 0101 0101
```

## 1) `dnakit.population.biopython_p19_cluster_functions.split_genepop_by_population()`

- **Purpose:** Write one GenePop file for each population of a GenePop file.
- **Method:** Reads the file with `Bio.PopGen.GenePop.read` and calls `Record.split_in_pops`. Population names in GenePop files are not reliable, so by default the files are `pop1.gen`, `pop2.gen`, ... in file order. Every output keeps all loci. The output directory must not exist.
- **API:** `dnakit.population.biopython_p19_cluster_functions.split_genepop_by_population(genepop_path, output_dir, *, population_names=None)`.
- **Input:** A GenePop file path, a new `output_dir` and optional `population_names` (one per population).
- **Output:** `population_count`, `locus_count` and `files` (name to path).
- **Example code:**

```python
import os, tempfile
from dnakit.population.biopython_p19_cluster_functions import split_genepop_by_population

target = os.path.join(tempfile.mkdtemp(), "pops")
result = split_genepop_by_population("pops.gen", target)
print(result["population_count"], result["locus_count"], sorted(result["files"]))
print(open(result["files"]["pop2"]).read())
```

- **Example output:**

```text
3 3 ['pop1', 'pop2', 'pop3']
Test genepop
Loc1
Loc2
Loc3
Pop
Ind3, 0201 0202 0303
Ind4, 0101 0102 0103
```

## 2) `dnakit.population.biopython_p19_cluster_functions.split_genepop_by_locus()`

- **Purpose:** Write one GenePop file for each locus of a GenePop file.
- **Method:** Calls `Bio.PopGen.GenePop.Record.split_in_loci`. Each output `<locus>.gen` holds one locus and all populations. The output directory must not exist.
- **API:** `dnakit.population.biopython_p19_cluster_functions.split_genepop_by_locus(genepop_path, output_dir)`.
- **Input:** A GenePop file path and a new `output_dir`.
- **Output:** `population_count`, `locus_count` and `files` (locus name to path).
- **Example code:**

```python
import os, tempfile
from dnakit.population.biopython_p19_cluster_functions import split_genepop_by_locus

target = os.path.join(tempfile.mkdtemp(), "loci")
result = split_genepop_by_locus("pops.gen", target)
print(result["locus_count"], sorted(result["files"]))
print(open(result["files"]["Loc2"]).read())
```

- **Example output:**

```text
3 ['Loc1', 'Loc2', 'Loc3']
Test genepop
Loc2
Pop
Ind1, 0202
Ind2, 0203
Pop
Ind3, 0202
Ind4, 0102
Pop
Ind5, 0101
```

## 3) `dnakit.population.biopython_p19_cluster_functions.remove_genepop_populations()`

- **Purpose:** Write a copy of a GenePop file without some populations.
- **Method:** Calls `Bio.PopGen.GenePop.Record.remove_population` for each 0-based position (highest first) and writes the record. The file-based `FileRecord` variant is not used because it writes wrongly padded alleles in Biopython 1.88. At least one population must remain; the output path must not exist.
- **API:** `dnakit.population.biopython_p19_cluster_functions.remove_genepop_populations(genepop_path, output_path, positions)`.
- **Input:** A GenePop file path, a new `output_path` and the `positions` of the populations to remove.
- **Output:** `removed_positions`, the remaining `population_count` and the `output` path.
- **Example code:**

```python
import os, tempfile
from dnakit.population.biopython_p19_cluster_functions import remove_genepop_populations

target = os.path.join(tempfile.mkdtemp(), "kept.gen")
result = remove_genepop_populations("pops.gen", target, [0, 2])
print(result["removed_positions"], result["population_count"])
print(open(target).read())
```

- **Example output:**

```text
[0, 2] 1
Test genepop
Loc1
Loc2
Loc3
Pop
Ind3, 0201 0202 0303
Ind4, 0101 0102 0103
```

## 4) `dnakit.population.biopython_p19_cluster_functions.remove_genepop_loci()`

- **Purpose:** Write a copy of a GenePop file without some loci.
- **Method:** Calls `Record.remove_locus_by_position` or `Record.remove_locus_by_name` of `Bio.PopGen.GenePop` and writes the record. Exactly one of `positions` (0-based) and `names` is required; at least one locus must remain.
- **API:** `dnakit.population.biopython_p19_cluster_functions.remove_genepop_loci(genepop_path, output_path, *, positions=None, names=None)`.
- **Input:** A GenePop file path, a new `output_path` and either `positions` or `names`.
- **Output:** `removed_loci`, `remaining_loci` and the `output` path.
- **Example code:**

```python
import os, tempfile
from dnakit.population.biopython_p19_cluster_functions import remove_genepop_loci

target = os.path.join(tempfile.mkdtemp(), "kept.gen")
result = remove_genepop_loci("pops.gen", target, names=["Loc2"])
print(result["removed_loci"], result["remaining_loci"])
print(open(target).read())
```

- **Example output:**

```text
['Loc2'] ['Loc1', 'Loc3']
Test genepop
Loc1
Loc3
Pop
Ind1, 0101 0303
Ind2, 0102 0000
Pop
Ind3, 0201 0303
Ind4, 0101 0103
Pop
Ind5, 0202 0101
```

## 5) `dnakit.population.biopython_p19_cluster_functions.kmeans_cluster()`

- **Purpose:** Cluster items with k-means or k-medians.
- **Method:** Calls `Bio.Cluster.kcluster`. `method` is `mean` (k-means) or `median` (k-medians); `distance` is one of `euclidean`, `cityblock`, `pearson`, `abs_pearson`, `uncentered`, `abs_uncentered`, `spearman` or `kendall`. The EM algorithm runs `n_passes` times from random starts, or once, deterministically, from `initial_clusters`.
- **API:** `dnakit.population.biopython_p19_cluster_functions.kmeans_cluster(data, n_clusters, *, method='mean', distance='euclidean', n_passes=1, initial_clusters=None, mask=None, weights=None, transpose=False)`.
- **Input:** A matrix as a list of rows (items by measurements), `n_clusters`, `method`, `distance`, `n_passes`, `initial_clusters`, `mask` (0 = missing), `weights` and `transpose` (cluster columns instead of rows).
- **Output:** `cluster_ids`, the within-cluster sum of distances `error` and `times_found`.
- **Example code:**

```python
from dnakit.population.biopython_p19_cluster_functions import kmeans_cluster

data = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]
means = kmeans_cluster(data, 2, initial_clusters=[0, 0, 1, 1, 0])
print(means["cluster_ids"], round(means["error"], 4), means["times_found"])
medians = kmeans_cluster(data, 2, method="median", distance="cityblock", initial_clusters=[0, 0, 1, 1, 0])
print(medians["cluster_ids"], round(medians["error"], 4))
```

- **Example output:**

```text
[0, 0, 1, 1, 0] 0.0256 1
[0, 0, 1, 1, 0] 0.3
```

## 6) `dnakit.population.biopython_p19_cluster_functions.kmedoids_cluster()`

- **Purpose:** Cluster items with k-medoids from a distance matrix.
- **Method:** Calls `Bio.Cluster.kmedoids`; only the lower triangle of the matrix is read. The cluster number reported for an item is the index of its medoid item.
- **API:** `dnakit.population.biopython_p19_cluster_functions.kmedoids_cluster(distance_matrix, n_clusters, *, n_passes=1, initial_clusters=None)`.
- **Input:** A square distance matrix, `n_clusters`, `n_passes` and `initial_clusters` (a start cluster number below `n_clusters` for each item).
- **Output:** `cluster_ids` (medoid item indices), `error` and `times_found`.
- **Example code:**

```python
from dnakit.population.biopython_p19_cluster_functions import kmedoids_cluster, profile_distance_matrix

matrix = profile_distance_matrix([[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]])["matrix"]
result = kmedoids_cluster(matrix, 2, initial_clusters=[0, 0, 1, 1, 0])
print(result["cluster_ids"], round(result["error"], 4), result["times_found"])
```

- **Example output:**

```text
[4, 4, 2, 2, 4] 0.0467 1
```

## 7) `dnakit.population.biopython_p19_cluster_functions.hierarchical_cluster()`

- **Purpose:** Hierarchical clustering with single, complete, average or centroid linkage.
- **Method:** Calls `Bio.Cluster.treecluster`. Give either the data matrix (all four linkages) or a square distance matrix (single, complete and average only). The tree is a list of merge nodes: an item `i >= 0` is a leaf and `-k` refers to the node created at step `k` (1-based), as in Cluster 3.0.
- **API:** `dnakit.population.biopython_p19_cluster_functions.hierarchical_cluster(data=None, *, distance_matrix=None, linkage='complete', distance='euclidean', mask=None, weights=None, transpose=False)`.
- **Input:** `data` or `distance_matrix` (not both), `linkage`, `distance`, `mask`, `weights` and `transpose`.
- **Output:** `item_count` and the `nodes` (`left`, `right`, `distance`).
- **Example code:**

```python
from dnakit.population.biopython_p19_cluster_functions import hierarchical_cluster

data = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]
for linkage in ("single", "complete", "average", "centroid"):
    nodes = hierarchical_cluster(data, linkage=linkage)["nodes"]
    print(linkage, [(n["left"], n["right"]) for n in nodes], round(nodes[-1]["distance"], 3))
```

- **Example output:**

```text
single [(0, 4), (2, 3), (1, -1), (-2, -3)] 47.153
complete [(4, 0), (3, 2), (1, -1), (-2, -3)] 49.98
average [(4, 0), (3, 2), (1, -1), (-2, -3)] 48.866
centroid [(0, 4), (-1, 1), (3, 2), (-2, -3)] 48.855
```

## 8) `dnakit.population.biopython_p19_cluster_functions.cut_cluster_tree()`

- **Purpose:** Cut a hierarchical clustering tree into a given number of clusters.
- **Method:** Rebuilds the `Bio.Cluster.Tree` from the node list and calls `Tree.cut`.
- **API:** `dnakit.population.biopython_p19_cluster_functions.cut_cluster_tree(nodes, n_clusters)`.
- **Input:** The `nodes` returned by `hierarchical_cluster` and `n_clusters`.
- **Output:** The `cluster_ids` of the items.
- **Example code:**

```python
from dnakit.population.biopython_p19_cluster_functions import hierarchical_cluster, cut_cluster_tree

nodes = hierarchical_cluster([[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]], linkage="average")["nodes"]
for count in (2, 3):
    print(count, cut_cluster_tree(nodes, count)["cluster_ids"])
```

- **Example output:**

```text
2 [1, 1, 0, 0, 1]
3 [2, 1, 0, 0, 2]
```

## 9) `dnakit.population.biopython_p19_cluster_functions.sort_cluster_tree()`

- **Purpose:** Order the leaves of a hierarchical clustering tree.
- **Method:** Calls `Bio.Cluster.Tree.sort`: the left and right subnodes are swapped so that the leaf order tends to follow `order` (default: the item index).
- **API:** `dnakit.population.biopython_p19_cluster_functions.sort_cluster_tree(nodes, *, order=None)`.
- **Input:** The `nodes` of a tree and an optional `order` (one number per item).
- **Output:** `leaf_order` (item indices from left to right) and the sorted `nodes`.
- **Example code:**

```python
from dnakit.population.biopython_p19_cluster_functions import hierarchical_cluster, sort_cluster_tree

nodes = hierarchical_cluster([[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]], linkage="average")["nodes"]
print(sort_cluster_tree(nodes)["leaf_order"])
print(sort_cluster_tree(nodes, order=[4, 3, 2, 1, 0])["leaf_order"])
```

- **Example output:**

```text
[3, 2, 1, 4, 0]
[3, 2, 4, 0, 1]
```

## 10) `dnakit.population.biopython_p19_cluster_functions.som_cluster()`

- **Purpose:** Cluster items with a self-organizing map on a rectangular grid.
- **Method:** Calls `Bio.Cluster.somcluster`. The map starts from random vectors, so results vary between calls.
- **API:** `dnakit.population.biopython_p19_cluster_functions.som_cluster(data, *, nx_grid=2, ny_grid=1, init_tau=0.02, n_iterations=100, distance='euclidean', mask=None, weights=None, transpose=False)`.
- **Input:** A matrix, the grid size `nx_grid` and `ny_grid`, `init_tau`, `n_iterations`, `distance`, `mask`, `weights` and `transpose`.
- **Output:** `cell_ids` (grid cell x, y of every item) and `cell_data[x][y]` (centroid vector of every cell).
- **Example code:**

```python
from dnakit.population.biopython_p19_cluster_functions import som_cluster

result = som_cluster([[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]], nx_grid=2, ny_grid=1, n_iterations=50)
print(len(result["cell_ids"]), len(result["cell_data"]), len(result["cell_data"][0][0]))
print(all(0 <= x < 2 and y == 0 for x, y in result["cell_ids"]))
```

- **Example output:**

```text
5 2 3
True
```

## 11) `dnakit.population.biopython_p19_cluster_functions.cluster_distance()`

- **Purpose:** Compute the distance between two clusters of items.
- **Method:** Calls `Bio.Cluster.clusterdistance`. `method` is `mean` or `median` (distance between the centres), `single` (smallest pairwise distance), `complete` (largest) or `average` (mean of the pairwise distances).
- **API:** `dnakit.population.biopython_p19_cluster_functions.cluster_distance(data, index1, index2, *, method='mean', distance='euclidean', mask=None, weights=None, transpose=False)`.
- **Input:** A matrix, the member item indices `index1` and `index2`, `method`, `distance`, `mask`, `weights` and `transpose`.
- **Output:** The cluster `distance`.
- **Example code:**

```python
from dnakit.population.biopython_p19_cluster_functions import cluster_distance

data = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]
for method in ("mean", "median", "single", "complete", "average"):
    print(method, round(cluster_distance(data, [0, 1], [2, 3], method=method)["distance"], 4))
```

- **Example output:**

```text
mean 48.5417
median 48.5417
single 47.1533
complete 49.9533
average 48.5517
```

## 12) `dnakit.population.biopython_p19_cluster_functions.cluster_centroids()`

- **Purpose:** Compute the centroid of every cluster.
- **Method:** Calls `Bio.Cluster.clustercentroids`: the mean (`method="mean"`) or the median over the items of each cluster, dimension by dimension.
- **API:** `dnakit.population.biopython_p19_cluster_functions.cluster_centroids(data, cluster_ids, *, method='mean', mask=None, transpose=False)`.
- **Input:** A matrix, `cluster_ids` (non-negative cluster number per item), `method`, `mask` and `transpose`.
- **Output:** One centroid row per cluster (`centroids`) and the matching 0/1 `centroid_mask`.
- **Example code:**

```python
from dnakit.population.biopython_p19_cluster_functions import cluster_centroids

data = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]
print(cluster_centroids(data, [0, 0, 1, 1, 0])["centroids"])
print(cluster_centroids(data, [0, 0, 1, 1, 0], method="median")["centroids"])
```

- **Example output:**

```text
[[1.0, 2.033333333333333, 3.1], [8.1, 9.05, 9.95]]
[[1.0, 2.0, 3.1], [8.1, 9.05, 9.95]]
```

## 13) `dnakit.population.biopython_p19_cluster_functions.profile_distance_matrix()`

- **Purpose:** Compute the distance matrix between the profiles (rows or columns) of a data matrix.
- **Method:** Calls `Bio.Cluster.distancematrix`. `distance` is `euclidean`, `cityblock`, `pearson` (1 - Pearson r), `abs_pearson` (1 - |r|), `uncentered` (1 - uncentered correlation), `abs_uncentered`, `spearman` (1 - rank correlation) or `kendall` (1 - Kendall tau). Missing values are excluded through `mask`.
- **API:** `dnakit.population.biopython_p19_cluster_functions.profile_distance_matrix(data, *, distance='euclidean', mask=None, weights=None, transpose=False)`.
- **Input:** A matrix, `distance`, `mask`, `weights` and `transpose`.
- **Output:** `item_count` and the full symmetric `matrix`.
- **Example code:**

```python
from dnakit.population.biopython_p19_cluster_functions import profile_distance_matrix

data = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]
for distance in ("euclidean", "cityblock", "pearson", "abs_pearson", "uncentered", "abs_uncentered", "spearman", "kendall"):
    row = profile_distance_matrix(data, distance=distance)["matrix"][0]
    print(distance, [round(value, 4) for value in row])
```

- **Example output:**

```text
euclidean [0.0, 0.02, 49.0, 49.9533, 0.0067]
cityblock [0.0, 0.1333, 7.0, 7.0667, 0.0667]
pearson [0.0, 0.0004, 0.0, 0.0006, -0.0]
abs_pearson [0.0, 0.0004, 0.0, 0.0006, -0.0]
uncentered [0.0, 0.0001, 0.0438, 0.048, 0.0006]
abs_uncentered [0.0, 0.0001, 0.0438, 0.048, 0.0006]
spearman [0.0, 0.0, 0.0, 0.0, 0.0]
kendall [0.0, 0.0, 0.0, 0.0, 0.0]
```

## 14) `dnakit.population.biopython_p19_cluster_functions.cluster_pca()`

- **Purpose:** Principal component analysis of a data matrix.
- **Method:** Calls `Bio.Cluster.pca` (centred singular value decomposition). Results are sorted by decreasing eigenvalue; the column means plus `coordinates @ components` reproduce the data.
- **API:** `dnakit.population.biopython_p19_cluster_functions.cluster_pca(data)`.
- **Input:** A matrix with at least two rows and two columns.
- **Output:** `column_means`, `components`, `coordinates` and `eigenvalues`.
- **Example code:**

```python
from dnakit.population.biopython_p19_cluster_functions import cluster_pca

result = cluster_pca([[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]])
print([round(value, 4) for value in result["column_means"]])
print([round(value, 4) for value in result["eigenvalues"]])
print([round(value, 4) for value in result["coordinates"][0]])
```

- **Example output:**

```text
[3.84, 4.84, 5.84]
[13.2636, 0.181, 0.0368]
[-4.9185, -0.0734, 0.0035]
```

## 15) `dnakit.regulation.biopython_p19_motif_functions.motif_matrices()`

- **Purpose:** Build a motif from binding sites and return its count, probability and log-odds matrices.
- **Method:** `Bio.motifs.create` makes the count matrix, `Motif.pwm` normalises it into a position weight matrix (with `pseudocounts`) and `Motif.pssm` gives the log2-odds matrix against `background` (uniform by default). Biopython's default of zero pseudocounts is kept, so an unseen base scores minus infinity, shown as `null`. It also returns the consensus, degenerate consensus and anticonsensus.
- **API:** `dnakit.regulation.biopython_p19_motif_functions.motif_matrices(instances=None, *, counts=None, pseudocounts=None, background=None)`.
- **Input:** Either equal-length DNA `instances` or `counts` (base to list of counts), `pseudocounts` (number or per-base mapping) and `background` (base probabilities).
- **Output:** `counts`, `pwm` and `pssm` as `{base: [value per position]}`, `background`, `consensus`, `degenerate_consensus` and `anticonsensus`.
- **Example code:**

```python
from dnakit.regulation.biopython_p19_motif_functions import motif_matrices

sites = ["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"]
plain = motif_matrices(sites)
print(plain["consensus"], plain["degenerate_consensus"], plain["anticonsensus"])
print(plain["counts"]["A"], [round(v, 3) for v in plain["pwm"]["A"]])
print(plain["pssm"]["C"])
smooth = motif_matrices(sites, pseudocounts=0.5)
print([round(v, 3) for v in smooth["pssm"]["C"]])
```

- **Example output:**

```text
TACGC WACVC CCATG
[3.0, 7.0, 0.0, 2.0, 1.0] [0.429, 1.0, 0.0, 0.286, 0.143]
[None, None, 1.5145731728297582, 0.19264507794239583, 1.7776075786635521]
[-2.17, -2.17, 1.29, 0.152, 1.531]
```

## 16) `dnakit.regulation.biopython_p19_motif_functions.motif_anticonsensus()`

- **Purpose:** Return the anticonsensus of a motif: the least frequent base at every position.
- **Method:** Calls `Bio.motifs.Motif.anticonsensus`. It is the sequence that matches the motif worst; the consensus and degenerate consensus are returned for comparison.
- **API:** `dnakit.regulation.biopython_p19_motif_functions.motif_anticonsensus(instances=None, *, counts=None)`.
- **Input:** Either equal-length DNA `instances` or `counts`.
- **Output:** `anticonsensus`, `consensus` and `degenerate_consensus`.
- **Example code:**

```python
from dnakit.regulation.biopython_p19_motif_functions import motif_anticonsensus

result = motif_anticonsensus(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"])
print(result["anticonsensus"], result["consensus"], result["degenerate_consensus"])
```

- **Example output:**

```text
CCATG TACGC WACVC
```

## 17) `dnakit.regulation.biopython_p19_motif_functions.pssm_search()`

- **Purpose:** Scan a DNA sequence with a motif and report the hits at or above a score threshold.
- **Method:** Calls `Bio.motifs.matrix.PositionSpecificScoringMatrix.search` on the log-odds matrix (both strands by default). Each hit has Biopython's `position` (negative on the reverse strand; add the sequence length to get the forward start), the forward 0-based half-open `start` and `end`, the `strand` and the `score`. `pseudocounts` defaults to 0.5 here so that unseen bases do not score minus infinity; pass 0 for Biopython's default.
- **API:** `dnakit.regulation.biopython_p19_motif_functions.pssm_search(sequence, instances=None, *, counts=None, threshold=0.0, both_strands=True, pseudocounts=0.5, background=None, max_hits=100000)`.
- **Input:** A DNA `sequence` (A, C, G, T, N), motif `instances` or `counts`, `threshold`, `both_strands`, `pseudocounts`, `background` and `max_hits`.
- **Output:** `hit_count`, `truncated` and the `hits`.
- **Example code:**

```python
from dnakit.regulation.biopython_p19_motif_functions import pssm_search

result = pssm_search("TACACGTACGCAACGCATTTGTAAC", ["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], threshold=3.0)
print(result["hit_count"])
for hit in result["hits"]:
    print(hit["position"], hit["start"], hit["end"], hit["strand"], round(hit["score"], 3))
```

- **Example output:**

```text
5
0 0 5 + 5.709
6 6 11 + 6.194
11 11 16 + 5.832
-11 14 19 - 4.694
-7 18 23 - 3.594
```

## 18) `dnakit.regulation.biopython_p19_motif_functions.pssm_score_statistics()`

- **Purpose:** Compute the mean, standard deviation, maximum and minimum of a motif's PSSM score.
- **Method:** `mean` and `std` are the expected value and standard deviation of the score of a random sequence drawn from `background` (`PositionSpecificScoringMatrix.mean` and `.std`); `max` and `min` are the scores of the consensus and the anticonsensus.
- **API:** `dnakit.regulation.biopython_p19_motif_functions.pssm_score_statistics(instances=None, *, counts=None, pseudocounts=0.5, background=None)`.
- **Input:** Motif `instances` or `counts`, `pseudocounts` (default 0.5) and `background`.
- **Output:** `mean`, `std`, `max` and `min` (`null` for infinity).
- **Example code:**

```python
from dnakit.regulation.biopython_p19_motif_functions import pssm_score_statistics

result = pssm_score_statistics(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"])
print({key: round(result[key], 4) for key in ("mean", "std", "max", "min")})
skewed = pssm_score_statistics(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], background={"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.4})
print(round(skewed["mean"], 4))
```

- **Example output:**

```text
{'mean': 3.1603, 'std': 2.5129, 'max': 6.1944, 'min': -10.8496}
4.4366
```

## 19) `dnakit.regulation.biopython_p19_motif_functions.pssm_score_distribution()`

- **Purpose:** Approximate the score distribution of a PSSM under the motif model and under the background.
- **Method:** Calls `PositionSpecificScoringMatrix.distribution` (`Bio.motifs.thresholds.ScoreDistribution`), a dynamic programme over `precision * width` score bins. Each bin has the score `min_score + i * step` and the probability of that score under the motif and under the background.
- **API:** `dnakit.regulation.biopython_p19_motif_functions.pssm_score_distribution(instances=None, *, counts=None, pseudocounts=0.5, background=None, precision=1000)`.
- **Input:** Motif `instances` or `counts`, `pseudocounts`, `background` and `precision`.
- **Output:** `min_score`, `step`, `bin_count`, `scores`, `motif_density` and `background_density`.
- **Example code:**

```python
from dnakit.regulation.biopython_p19_motif_functions import pssm_score_distribution

result = pssm_score_distribution(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], precision=100)
print(result["bin_count"], round(result["min_score"], 3), round(result["step"], 4))
print(round(sum(result["motif_density"]), 6), round(sum(result["background_density"]), 6))
best = max(range(result["bin_count"]), key=lambda i: result["motif_density"][i])
print(round(result["scores"][best], 3))
```

- **Example output:**

```text
500 -10.85 0.0342
1.0 1.0
5.716
```

## 20) `dnakit.regulation.biopython_p19_motif_functions.pssm_threshold()`

- **Purpose:** Choose a score threshold for a motif PSSM.
- **Method:** From the score distribution (`Bio.motifs.thresholds.ScoreDistribution`): `fpr` gives the threshold for a false-positive rate `rate` against the background (`threshold_fpr`), `fnr` for a false-negative rate against the motif (`threshold_fnr`), `balanced` the threshold where the false-negative rate is `rate_proportion` times the false-positive rate (`threshold_balanced`, also returns the rate) and `patser` the threshold used by the Patser program (`threshold_patser`).
- **API:** `dnakit.regulation.biopython_p19_motif_functions.pssm_threshold(instances=None, *, counts=None, method='fpr', rate=None, rate_proportion=1.0, pseudocounts=0.5, background=None, precision=1000)`.
- **Input:** Motif `instances` or `counts`, `method`, `rate` (required for `fpr` and `fnr`), `rate_proportion`, `pseudocounts`, `background` and `precision`.
- **Output:** `threshold` (and `rate` for `balanced`).
- **Example code:**

```python
from dnakit.regulation.biopython_p19_motif_functions import pssm_threshold

sites = ["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"]
print(round(pssm_threshold(sites, method="fpr", rate=0.01)["threshold"], 4))
print(round(pssm_threshold(sites, method="fnr", rate=0.1)["threshold"], 4))
balanced = pssm_threshold(sites, method="balanced")
print(round(balanced["threshold"], 4), round(balanced["rate"], 4))
print(round(pssm_threshold(sites, method="patser")["threshold"], 4))
```

- **Example output:**

```text
4.2101
-0.4336
0.1698 0.126
0.4221
```

## 21) `dnakit.regulation.biopython_p19_motif_functions.jaspar_pseudocounts()`

- **Purpose:** Compute pseudocounts the way JASPAR does.
- **Method:** Calls `Bio.motifs.jaspar.calculate_pseudocounts`: the square root of the number of sites times the background probability of each base.
- **API:** `dnakit.regulation.biopython_p19_motif_functions.jaspar_pseudocounts(instances=None, *, counts=None, background=None)`.
- **Input:** Motif `instances` or `counts` and `background`.
- **Output:** `pseudocounts` (one number per base).
- **Example code:**

```python
from dnakit.regulation.biopython_p19_motif_functions import jaspar_pseudocounts

print(jaspar_pseudocounts(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"])["pseudocounts"])
print(jaspar_pseudocounts(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], background={"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3})["pseudocounts"])
```

- **Example output:**

```text
{'A': 0.6614378277661477, 'C': 0.6614378277661477, 'G': 0.6614378277661477, 'T': 0.6614378277661477}
{'A': 0.7937253933193772, 'C': 0.5291502622129182, 'G': 0.5291502622129182, 'T': 0.7937253933193772}
```

## 22) `dnakit.regulation.biopython_p19_motif_functions.weblogo_service_logo()`

- **Purpose:** Draw a sequence logo with the WebLogo web service.
- **Method:** Calls `Bio.motifs.Motif.weblogo`, which sends the motif to https://weblogo.threeplusone.com (WebLogo 3.3) and saves the returned image to a new file. `options` are WebLogo form fields (unknown names are refused). Needs network access; local logos without a service are made by `plot_pwm_logo`.
- **API:** `dnakit.regulation.biopython_p19_motif_functions.weblogo_service_logo(instances, output_path, *, counts=None, image_format='png', options=None)`.
- **Input:** Motif `instances` (or `counts`), a new `output_path`, `image_format` (`png`, `png_print`, `pdf`, `jpeg`, `eps`, `svg` or `txt`) and `options`.
- **Output:** `output` path, `image_format`, `bytes` and the `service` used.
- **Example code:**

```python
import os, tempfile
from dnakit.regulation.biopython_p19_motif_functions import weblogo_service_logo

target = os.path.join(tempfile.mkdtemp(), "logo.png")
result = weblogo_service_logo(["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], target, options={"logo_title": "demo"})
print(result["status"], result["image_format"], result["bytes"] > 1000, open(target, "rb").read(4))
```

- **Example output:**

```text
completed png True b'\x89PNG'
```
