# 102 Review additions (G18p): toytree, dendropy and pyranges1

This page documents the 42 wrappers added by review group G18p: toytree node queries, enumeration, tree editing, consensus support and network parsing (26, [toytree](https://github.com/eaton-lab/toytree)); DendroPy bipartition comparison and tree-shape vectors (6, [DendroPy](https://github.com/jeetsukumaran/DendroPy)); and pyranges1 readers, interval concatenation, ORFs and overlap statistics (10, [pyranges1](https://pypi.org/project/pyranges1/)). Everything runs locally: no network service or external program is called and nothing is drawn. Tree inputs are Newick text or one-tree files; pyranges1 inputs are local files. Optional extras: `bio-toytree`, `bio-dendropy`, `bio-pyranges1`; bigWig and ORF functions also need [pyBigWig](https://github.com/deeptools/pyBigWig) (`bio-pybigwig`) and [pyfaidx](https://github.com/mdshw5/pyfaidx) (`bio-pyfaidx`). The modules are `dnakit.phylogeny`, `dnakit.phylogeny` and `dnakit.annotation`; the package exports are wired by the lead.

## 1) `dnakit.phylogeny.list_tree_node_path()`

- **Purpose:** List the nodes on the path between two nodes, including both ends.
- **Method:** Calls toytree.distance.get_node_path; each node is reported with idx, name and is_tip.
- **API:** `dnakit.phylogeny.list_tree_node_path(tree, start, end, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** Newick text or a one-tree Newick file; start and end are node names or integer indexes.
- **Output:** nodes (path order), edge_count, and provenance.
- **Example code:**

```python
from dnakit.phylogeny import list_tree_node_path

result = list_tree_node_path("((A:1,B:2):1,(C:3,D:1):2);", "A", "C")["result"]
print([node["idx"] for node in result["nodes"]], result["edge_count"])
```

- **Example output:**

```text
[0, 4, 6, 5, 2] 4
```

## 2) `dnakit.phylogeny.calculate_tree_node_distance()`

- **Purpose:** Compute the branch-length distance between two nodes, or count edges only.
- **Method:** Calls toytree.distance.get_node_distance; with topology_only=True each edge counts 1.
- **API:** `dnakit.phylogeny.calculate_tree_node_distance(tree, start, end, *, topology_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree; start and end are node names or indexes; topology_only is a boolean.
- **Output:** distance (number), topology_only and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_tree_node_distance

tree = "((A:1,B:2):1,(C:3,D:1):2);"
print(calculate_tree_node_distance(tree, "A", "C")["result"]["distance"])
print(calculate_tree_node_distance(tree, "A", "C", topology_only=True)["result"]["distance"])
```

- **Example output:**

```text
7.0
4.0
```

## 3) `dnakit.phylogeny.calculate_all_node_distance_matrix()`

- **Purpose:** Compute pairwise distances between all nodes, tips and internal nodes.
- **Method:** Calls toytree.distance.get_node_distance_matrix; rows and columns follow node index, unnamed internal nodes are labelled node_<idx>.
- **API:** `dnakit.phylogeny.calculate_all_node_distance_matrix(tree, *, topology_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree; topology_only is a boolean.
- **Output:** labels (node labels), matrix (symmetric), and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_all_node_distance_matrix

result = calculate_all_node_distance_matrix("((A:1,B:2):1,(C:3,D:1):2);")["result"]
print(result["labels"])
print(result["matrix"][0])
```

- **Example output:**

```text
['A', 'B', 'C', 'D', 'node_4', 'node_5', 'node_6']
[0.0, 3.0, 7.0, 5.0, 1.0, 4.0, 2.0]
```

## 4) `dnakit.phylogeny.calculate_internal_node_distance_matrix()`

- **Purpose:** Compute pairwise distances between internal nodes, including the root.
- **Method:** Calls toytree.distance.get_internal_node_distance_matrix, ordered by node index.
- **API:** `dnakit.phylogeny.calculate_internal_node_distance_matrix(tree, *, topology_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree; topology_only is a boolean.
- **Output:** labels, matrix and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_internal_node_distance_matrix

result = calculate_internal_node_distance_matrix("((A:1,B:2):1,(C:3,D:1):2);")["result"]
print(result["labels"])
print(result["matrix"])
```

- **Example output:**

```text
['node_4', 'node_5', 'node_6']
[[0.0, 3.0, 1.0], [3.0, 0.0, 2.0], [1.0, 2.0, 0.0]]
```

## 5) `dnakit.phylogeny.calculate_descendant_distances()`

- **Purpose:** Give the distance from a node (default: the root) to each of its descendants.
- **Method:** Calls toytree.distance.get_descendant_dists; results are sorted by node index.
- **API:** `dnakit.phylogeny.calculate_descendant_distances(tree, node=None, *, topology_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree; node is a name, an index, or None for the root.
- **Output:** node, descendants (idx, name, is_tip and distance for each) and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_descendant_distances

result = calculate_descendant_distances("((A:1,B:2):1,(C:3,D:1):2);")["result"]
print([(row["idx"], row["distance"]) for row in result["descendants"]])
```

- **Example output:**

```text
[(0, 2.0), (1, 3.0), (2, 5.0), (3, 3.0), (4, 1.0), (5, 2.0), (6, 0.0)]
```

## 6) `dnakit.phylogeny.find_tree_farthest_node()`

- **Purpose:** Find the node farthest from a given node, with its distance.
- **Method:** Calls toytree.distance.get_farthest_node and get_farthest_node_distance.
- **API:** `dnakit.phylogeny.find_tree_farthest_node(tree, node, *, topology_only=False, descendants_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree; node is a name or index; topology_only and descendants_only are booleans.
- **Output:** node, farthest (the farthest node), distance and provenance.
- **Example code:**

```python
from dnakit.phylogeny import find_tree_farthest_node

result = find_tree_farthest_node("((A:1,B:2):1,(C:3,D:1):2);", "A")["result"]
print(result["farthest"]["name"], result["distance"])
```

- **Example output:**

```text
C 7.0
```

## 7) `dnakit.phylogeny.calculate_farthest_node_distance()`

- **Purpose:** Compute the distance from a node (default: the root) to the farthest node.
- **Method:** Calls toytree.distance.get_farthest_node_distance.
- **API:** `dnakit.phylogeny.calculate_farthest_node_distance(tree, node=None, *, topology_only=False, descendants_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree; node is a name, an index or None; the boolean options are as above.
- **Output:** node, distance, topology_only, descendants_only and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_farthest_node_distance

print(calculate_farthest_node_distance("((A:1,B:2):1,(C:3,D:1):2);", "A")["result"]["distance"])
```

- **Example output:**

```text
7.0
```

## 8) `dnakit.phylogeny.list_tree_quadripartitions()`

- **Purpose:** List the four leaf groups around each internal edge (quadripartitions).
- **Method:** Calls toytree.enum.iter_quadripartitions and flattens each into four sorted leaf lists.
- **API:** `dnakit.phylogeny.list_tree_quadripartitions(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree.
- **Output:** quadripartitions (four groups per item), count and provenance.
- **Example code:**

```python
from dnakit.phylogeny import list_tree_quadripartitions

result = list_tree_quadripartitions("((A,B),(C,D),E);")["result"]
print(result["count"])
for groups in result["quadripartitions"]:
    print(groups)
```

- **Example output:**

```text
2
[['A'], ['B'], ['C', 'D'], ['E']]
[['C'], ['D'], ['A', 'B'], ['E']]
```

## 9) `dnakit.phylogeny.list_tree_quartets()`

- **Purpose:** List the quartets induced by internal edges (two leaf groups on each side).
- **Method:** Calls toytree.enum.iter_quartets and sorts by the left and right groups.
- **API:** `dnakit.phylogeny.list_tree_quartets(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree.
- **Output:** quartets (left and right leaf groups), count and provenance.
- **Example code:**

```python
from dnakit.phylogeny import list_tree_quartets

result = list_tree_quartets("((A,B),(C,D),E);")["result"]
print(result["count"], result["quartets"][0])
```

- **Example output:**

```text
5 {'left': ['A', 'B'], 'right': ['C', 'D']}
```

## 10) `dnakit.phylogeny.map_consensus_tree_features()`

- **Purpose:** Map split support and branch statistics from a list of trees onto a target tree.
- **Method:** Builds a MultiTree with toytree.mtree and calls toytree.infer.get_consensus_features; support is the fraction of trees holding the split.
- **API:** `dnakit.phylogeny.map_consensus_tree_features(tree, trees, *, features=('support',), max_trees=1000, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A target Newick tree, a list of Newick trees with the same tips, and feature names.
- **Output:** newick, features, tree_count, nodes (tips and values per node) and provenance.
- **Example code:**

```python
from dnakit.phylogeny import map_consensus_tree_features

trees = ["((A,B),(C,D));", "((A,B),(C,D));", "((A,C),(B,D));"]
result = map_consensus_tree_features("((A,B),(C,D));", trees)["result"]
print(result["newick"])
node = next(row for row in result["nodes"] if row["tips"] == ["C", "D"])
print(round(node["values"]["support"], 4))
```

- **Example output:**

```text
(A:1,B:1,(C:1,D:1)0.666666666667:1);
0.6667
```

## 11) `dnakit.phylogeny.extend_tree_tips_to_align()`

- **Purpose:** Extend every tip branch so that all tips end at the same depth.
- **Method:** Calls toytree.mod.edges_extend_tips_to_align and returns the edited tree as Newick.
- **API:** `dnakit.phylogeny.extend_tree_tips_to_align(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree.
- **Output:** newick, tip_count, tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import extend_tree_tips_to_align

print(extend_tree_tips_to_align("((A:1,B:2):1,(C:3,D:1):2);")["result"]["newick"])
```

- **Example output:**

```text
((A:4,B:4):1,(C:3,D:3):2);
```

## 12) `dnakit.phylogeny.set_tree_node_heights()`

- **Purpose:** Set the height (distance from the root) of named or indexed nodes.
- **Method:** Calls toytree.mod.edges_set_node_heights; keys are node names or indexes.
- **API:** `dnakit.phylogeny.set_tree_node_heights(tree, heights, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree and a {node: height} mapping with non-negative heights.
- **Output:** newick, heights (node and height list) and provenance.
- **Example code:**

```python
from dnakit.phylogeny import set_tree_node_heights

print(set_tree_node_heights("((A:1,B:2):1,(C:3,D:1):2);", {"A": 2.0})["result"]["newick"])
```

- **Example output:**

```text
((A:2,B:2):1,(C:3,D:1):2);
```

## 13) `dnakit.phylogeny.rotate_tree_node()`

- **Purpose:** Reverse the child order of one node.
- **Method:** Calls toytree.mod.rotate_node and returns the edited tree as Newick.
- **API:** `dnakit.phylogeny.rotate_tree_node(tree, node, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree and a node name or index.
- **Output:** newick, node, tip_count, tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import rotate_tree_node

print(rotate_tree_node("((A:1,B:2):1,(C:3,D:1):2);", 4)["result"]["newick"])
```

- **Example output:**

```text
((B:2,A:1):1,(C:3,D:1):2);
```

## 14) `dnakit.phylogeny.drop_tree_tips()`

- **Purpose:** Remove the named tips and return the remaining tree.
- **Method:** Calls toytree.mod.drop_tips; the call fails if fewer than two tips would remain.
- **API:** `dnakit.phylogeny.drop_tree_tips(tree, tips, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree and a list of tip names present in the tree.
- **Output:** newick, dropped (sorted tip names), tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import drop_tree_tips

print(drop_tree_tips("((A:1,B:2):1,(C:3,D:1):2);", ["A"])["result"]["newick"])
```

- **Example output:**

```text
((C:3,D:1):2,B:3);
```

## 15) `dnakit.phylogeny.bisect_tree()`

- **Purpose:** Split a tree on the edge above a node into the clade and the remaining tree.
- **Method:** Calls toytree.mod.bisect with reroot=False; dist_partition must be positive.
- **API:** `dnakit.phylogeny.bisect_tree(tree, node, *, dist_partition=1.0, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree, a node name or index, and dist_partition.
- **Output:** node, subtree_newick and subtree_tips, remainder_newick and remainder_tips, provenance.
- **Example code:**

```python
from dnakit.phylogeny import bisect_tree

result = bisect_tree("((A:1,B:2):1,(C:3,D:1):2);", 4)["result"]
print(result["subtree_newick"], result["remainder_newick"])
```

- **Example output:**

```text
(A:1,B:2):1; ((C:3,D:1):2);
```

## 16) `dnakit.phylogeny.resolve_tree_polytomies()`

- **Purpose:** Resolve multifurcations at random with a fixed seed, giving a bifurcating tree.
- **Method:** Calls toytree.mod.resolve_polytomies with numpy default_rng(seed); no global random state is changed.
- **API:** `dnakit.phylogeny.resolve_tree_polytomies(tree, *, node=None, seed, dist=1e-06, recursive=True, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree, an optional node, a required integer seed (0 to 2**32-1), dist and recursive.
- **Output:** newick, node, seed, dist, recursive, tip_count, tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import resolve_tree_polytomies

print(resolve_tree_polytomies("(A,B,C,D);", seed=3)["result"]["newick"])
```

- **Example output:**

```text
(D:1,(A:1,(B:1,C:1):1e-06):1e-06);
```

## 17) `dnakit.phylogeny.add_tree_internal_node()`

- **Purpose:** Insert an internal node on the edge above a node.
- **Method:** Calls toytree.mod.add_internal_node; name and dist are optional.
- **API:** `dnakit.phylogeny.add_tree_internal_node(tree, node, *, name=None, dist=None, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree, a node name or index, and optional name and dist (non-negative).
- **Output:** newick, node, name, dist, tip_count, tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import add_tree_internal_node

result = add_tree_internal_node("((A:1,B:2):1,(C:3,D:1):2);", "A", name="X", dist=0.5)["result"]
print(result["newick"], result["name"])
```

- **Example output:**

```text
((B:2,(A:0.5):0.5):1,(C:3,D:1):2); X
```

## 18) `dnakit.phylogeny.add_tree_child_node()`

- **Purpose:** Add a new tip as a child of a node.
- **Method:** Calls toytree.mod.add_child_node. If the target is a tip, it becomes an unnamed internal node whose only child is the new tip.
- **API:** `dnakit.phylogeny.add_tree_child_node(tree, node, *, name=None, dist=None, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree, a node name or index, and optional name and dist.
- **Output:** newick, node, name, dist, tip_count, tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import add_tree_child_node

print(add_tree_child_node("((A:1,B:2):1,(C:3,D:1):2);", 4, name="Z", dist=0.5)["result"]["newick"])
```

- **Example output:**

```text
((A:1,B:2,Z:0.5):1,(C:3,D:1):2);
```

## 19) `dnakit.phylogeny.add_tree_sister_node()`

- **Purpose:** Add a new tip as a sister of a node.
- **Method:** Calls toytree.mod.add_sister_node.
- **API:** `dnakit.phylogeny.add_tree_sister_node(tree, node, *, name=None, dist=None, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree, a node name or index, and optional name and dist.
- **Output:** newick, node, name, dist, tip_count, tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import add_tree_sister_node

print(add_tree_sister_node("((A:1,B:2):1,(C:3,D:1):2);", "A", name="Z", dist=0.5)["result"]["newick"])
```

- **Example output:**

```text
((A:1,B:2,Z:0.5):1,(C:3,D:1):2);
```

## 20) `dnakit.phylogeny.add_tree_internal_node_and_child()`

- **Purpose:** Split the edge above a node with a new internal node that gets a new tip child.
- **Method:** Calls toytree.mod.add_internal_node_and_child; the parent node name and edge length are optional.
- **API:** `dnakit.phylogeny.add_tree_internal_node_and_child(tree, node, *, name=None, dist=None, parent_name=None, parent_dist=None, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree, a node name or index, and optional name, dist, parent_name, parent_dist.
- **Output:** newick, node, name, dist, parent_name, parent_dist, tip_count, tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import add_tree_internal_node_and_child

result = add_tree_internal_node_and_child(
    "((A:1,B:2):1,(C:3,D:1):2);", "A", name="Z", dist=0.5, parent_name="P", parent_dist=0.5
)["result"]
print(result["newick"], result["tips"])
```

- **Example output:**

```text
((B:2,(A:0.5,Z:0.5):0.5):1,(C:3,D:1):2); ['A', 'B', 'C', 'D', 'Z']
```

## 21) `dnakit.phylogeny.add_tree_internal_node_and_subtree()`

- **Purpose:** Graft a Newick subtree onto the edge above a node through a new internal node.
- **Method:** Calls toytree.mod.add_internal_node_and_subtree; subtree tips must not repeat tips of the tree.
- **API:** `dnakit.phylogeny.add_tree_internal_node_and_subtree(tree, node, subtree, *, subtree_stem_dist=None, subtree_rescale=False, parent_name=None, parent_dist=None, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree, a node, a Newick subtree, and optional subtree_stem_dist, subtree_rescale, parent_name, parent_dist.
- **Output:** newick, node, subtree_tips, subtree_stem_dist, subtree_rescale, parent_name, parent_dist, tip_count, tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import add_tree_internal_node_and_subtree

result = add_tree_internal_node_and_subtree(
    "((A:1,B:2):1,(C:3,D:1):2);", "A", "(X:1,Y:1);", subtree_stem_dist=0.5, parent_name="P", parent_dist=0.5
)["result"]
print(result["newick"], result["subtree_tips"])
```

- **Example output:**

```text
((B:2,(A:0.5,(X:1,Y:1):0.5):0.5):1,(C:3,D:1):2); ['X', 'Y']
```

## 22) `dnakit.phylogeny.remove_tree_nodes()`

- **Purpose:** Remove a list of nodes and keep their children.
- **Method:** Calls toytree.mod.remove_nodes; with preserve_dists=True branch lengths of removed nodes are added to their children. At least two tips must remain.
- **API:** `dnakit.phylogeny.remove_tree_nodes(tree, nodes, *, preserve_dists=True, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree, a list of node names or indexes, and preserve_dists.
- **Output:** newick, removed, preserve_dists, tip_count, tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import remove_tree_nodes

print(remove_tree_nodes("((A:1,B:2):1,(C:3,D:1):2);", [4])["result"]["newick"])
```

- **Example output:**

```text
((C:3,D:1):2,A:2,B:3);
```

## 23) `dnakit.phylogeny.merge_tree_nodes()`

- **Purpose:** Merge clades whose tips share one node attribute, keeping one representative per clade.
- **Method:** Calls toytree.mod.merge_nodes with merge_method set to a node attribute name, which must be an identifier.
- **API:** `dnakit.phylogeny.merge_tree_nodes(tree, feature, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree and a node attribute name such as dist.
- **Output:** newick, feature, removed_node_count, tip_count, tips and provenance.
- **Example code:**

```python
from dnakit.phylogeny import merge_tree_nodes

result = merge_tree_nodes("((A:1,B:1):1,(C:1,D:1):1);", "dist")["result"]
print(result["newick"], result["removed_node_count"])
```

- **Example output:**

```text
(A:2,C:2); 4
```

## 24) `dnakit.phylogeny.parse_tree_network()`

- **Purpose:** Parse an extended Newick network with hybrid (#H) nodes into a major tree and admixture edges.
- **Method:** Calls toytree.network.parse_network_to_tree_and_admix; the input can be network text or a file.
- **API:** `dnakit.phylogeny.parse_tree_network(network, *, max_file_bytes=10000000, verbose=True)`.
- **Input:** Network Newick text (containing ;) or a file path.
- **Output:** newick, tips, admixture_edges (destination_tips, source_tips, gamma), admixture_count and provenance.
- **Example code:**

```python
from dnakit.phylogeny import parse_tree_network

result = parse_tree_network("((A:1,(B:1,#H1:0::0.3):1):1,(#H1:0::0.7,C:2):1);")["result"]
print(result["tips"])
print(result["admixture_edges"])
```

- **Example output:**

```text
['A', 'B', 'C']
[{'destination_tips': [], 'source_tips': ['B'], 'gamma': 0.3}]
```

## 25) `dnakit.phylogeny.calculate_tip_distance_matrix_from_vcv()`

- **Purpose:** Compute tip-to-tip distances from a tip variance-covariance matrix.
- **Method:** Calls toytree.pcm.get_distance_matrix_from_vcv_matrix; each entry is V[i,i]+V[j,j]-2V[i,j], the squared Euclidean distance (upstream takes no square root).
- **API:** `dnakit.phylogeny.calculate_tip_distance_matrix_from_vcv(vcv, labels, *, verbose=True)`.
- **Input:** A symmetric square vcv matrix and a list of unique tip labels of the same length.
- **Output:** labels, matrix and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_tip_distance_matrix_from_vcv

result = calculate_tip_distance_matrix_from_vcv([[2, 1], [1, 2]], ["x", "y"])["result"]
print(result["labels"], result["matrix"])
```

- **Example output:**

```text
['x', 'y'] [[0.0, 2.0], [2.0, 0.0]]
```

## 26) `dnakit.phylogeny.calculate_tip_diversification_rates()`

- **Purpose:** Compute tip-level diversification rates (Jetz et al. 2012) for every tip.
- **Method:** Calls toytree.pcm.get_tip_level_diversification; only tip DR values are kept and the unnamed root row is dropped.
- **API:** `dnakit.phylogeny.calculate_tip_diversification_rates(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree.
- **Output:** rates (tip and diversification_rate per item), tip_count and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_tip_diversification_rates

for row in calculate_tip_diversification_rates("((A:1,B:2):1,(C:3,D:1):2);")["result"]["rates"]:
    print(row["tip"], row["diversification_rate"])
```

- **Example output:**

```text
A 2.0
B 2.0
C 1.0
D 1.0
```

## 27) `dnakit.phylogeny.list_missing_bipartitions()`

- **Purpose:** List the splits of the reference tree that the comparison tree does not contain.
- **Method:** Calls dendropy.calculate.treecompare.find_missing_bipartitions; each split is reported as its smaller side (ties: the side without the first tip name).
- **API:** `dnakit.phylogeny.list_missing_bipartitions(reference_tree, comparison_tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** Two Newick trees with the same tip labels.
- **Output:** bipartitions (tip name lists), missing_count, tip_count and provenance.
- **Example code:**

```python
from dnakit.phylogeny import list_missing_bipartitions

result = list_missing_bipartitions("((A,B),(C,D),E);", "((A,C),(B,D),E);")["result"]
print(result["bipartitions"], result["missing_count"])
```

- **Example output:**

```text
[['A', 'B'], ['C', 'D']] 2
```

## 28) `dnakit.phylogeny.calculate_n_bar_index()`

- **Purpose:** Compute the N-bar tree-shape statistic: the average number of nodes above a tip.
- **Method:** Calls dendropy.calculate.treemeasure.N_bar.
- **API:** `dnakit.phylogeny.calculate_n_bar_index(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree.
- **Output:** n_bar, tip_count and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_n_bar_index

print(calculate_n_bar_index("((A,B),C);")["result"]["n_bar"])
```

- **Example output:**

```text
1.666666666667
```

## 29) `dnakit.phylogeny.calculate_node_ages()`

- **Purpose:** Compute node ages in backward time (distance from each node to the deepest tip).
- **Method:** Calls dendropy.calculate.treemeasure.node_ages; every edge except the root edge needs a length.
- **API:** `dnakit.phylogeny.calculate_node_ages(tree, *, internal_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree with branch lengths and the internal_only boolean.
- **Output:** internal_only, values (ascending list), tip_count and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_node_ages

print(calculate_node_ages("((A:1,B:2):1,(C:3,D:1):2);")["result"]["values"])
```

- **Example output:**

```text
[0.0, 2.0, 2.0, 3.0, 3.0, 4.0, 5.0]
```

## 30) `dnakit.phylogeny.calculate_node_depths()`

- **Purpose:** Compute node depths in forward time (distance from the root to each node).
- **Method:** Calls dendropy.calculate.treemeasure.node_depths; branch lengths are required.
- **API:** `dnakit.phylogeny.calculate_node_depths(tree, *, internal_only=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree with branch lengths and the internal_only boolean.
- **Output:** internal_only, values (ascending list), tip_count and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_node_depths

print(calculate_node_depths("((A:1,B:2):1,(C:3,D:1):2);", internal_only=True)["result"]["values"])
```

- **Example output:**

```text
[0.0, 1.0, 2.0]
```

## 31) `dnakit.phylogeny.calculate_coalescence_ages()`

- **Purpose:** Compute the ages of internal nodes (coalescence events) in backward time.
- **Method:** Calls dendropy.calculate.treemeasure.coalescence_ages; branch lengths are required.
- **API:** `dnakit.phylogeny.calculate_coalescence_ages(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree with branch lengths.
- **Output:** values (ascending list), tip_count and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_coalescence_ages

print(calculate_coalescence_ages("((A:1,B:2):1,(C:3,D:1):2);")["result"]["values"])
```

- **Example output:**

```text
[3.0, 4.0, 5.0]
```

## 32) `dnakit.phylogeny.calculate_divergence_times()`

- **Purpose:** Compute the times of internal nodes (splits) in forward time from the root.
- **Method:** Calls dendropy.calculate.treemeasure.divergence_times; branch lengths are required.
- **API:** `dnakit.phylogeny.calculate_divergence_times(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick tree with branch lengths.
- **Output:** values (ascending list), tip_count and provenance.
- **Example code:**

```python
from dnakit.phylogeny import calculate_divergence_times

print(calculate_divergence_times("((A:1,B:2):1,(C:3,D:1):2);")["result"]["values"])
```

- **Example output:**

```text
[0.0, 1.0, 2.0]
```

## 33) `dnakit.annotation.read_pyranges1_gff()`

- **Purpose:** Read a GTF or GFF file with pyranges1, one record per feature with its attributes.
- **Method:** Calls pyranges1.read_gff, reading at most max_rows+1 rows; exceeding max_rows is an error. Coordinates are 0-based half-open.
- **API:** `dnakit.annotation.read_pyranges1_gff(path, *, max_rows=100000, duplicate_attr=False, ignore_bad=False, verbose=True)`.
- **Input:** A GTF or GFF file path, with optional max_rows, duplicate_attr and ignore_bad.
- **Output:** records (chrom, start, end, strand and attribute columns), record_count and provenance.

**Example file:** `a.gtf`

```text
chr1	src	gene	100	500	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	CDS	100	130	.	+	0	gene_id "g1"; transcript_id "t1";
chr1	src	CDS	200	260	.	+	0	gene_id "g1"; transcript_id "t1";
chr1	src	CDS	400	450	.	-	0	gene_id "g2"; transcript_id "t2";
```

- **Example code:**

```python
from dnakit.annotation import read_pyranges1_gff

result = read_pyranges1_gff("a.gtf")
print(result["record_count"])
first = result["records"][0]
print(first["chrom"], first["start"], first["end"], first["Feature"], first["transcript_id"])
```

- **Example output:**

```text
4
chr1 99 130 CDS t1
```

## 34) `dnakit.annotation.read_pyranges1_gff3()`

- **Purpose:** Read a GFF3 file with pyranges1; attribute keys become record fields.
- **Method:** Calls pyranges1.read_gff3; starts are converted from 1-based to 0-based.
- **API:** `dnakit.annotation.read_pyranges1_gff3(path, *, max_rows=100000, verbose=True)`.
- **Input:** A GFF3 file path and optional max_rows.
- **Output:** records, record_count and provenance.

**Example file:** `a.gff3`

```text
##gff-version 3
chr1	src	gene	101	500	.	+	.	ID=g1;Name=G1
chr1	src	CDS	101	130	.	+	0	ID=c1;Parent=g1
```

- **Example code:**

```python
from dnakit.annotation import read_pyranges1_gff3

for record in read_pyranges1_gff3("a.gff3")["records"]:
    print(record["Feature"], record["start"], record["end"], record["ID"])
```

- **Example output:**

```text
CDS 100 130 c1
gene 100 500 g1
```

## 35) `dnakit.annotation.read_pyranges1_gtf_full()`

- **Purpose:** Read a GTF file with pyranges1's full attribute parser, one record per feature.
- **Method:** Calls pyranges1.readers.read_gtf_full, reading at most max_rows+1 rows.
- **API:** `dnakit.annotation.read_pyranges1_gtf_full(path, *, max_rows=100000, duplicate_attr=False, ignore_bad=False, verbose=True)`.
- **Input:** A GTF file path with optional max_rows, duplicate_attr and ignore_bad.
- **Output:** records, record_count and provenance.
- **Example code:**

```python
from dnakit.annotation import read_pyranges1_gtf_full

result = read_pyranges1_gtf_full("a.gtf")
print(result["record_count"], sorted({record["transcript_id"] for record in result["records"]}))
```

- **Example output:**

```text
4 ['t1', 't2']
```

## 36) `dnakit.annotation.read_pyranges1_bigwig()`

- **Purpose:** Read every interval of a bigWig signal track with its value.
- **Method:** Calls pyranges1.read_bigwig; the file size is checked against max_file_bytes before the max_rows check.
- **API:** `dnakit.annotation.read_pyranges1_bigwig(path, *, max_file_bytes=50000000, max_rows=100000, verbose=True)`.
- **Input:** A bigWig file path with optional max_file_bytes and max_rows.
- **Output:** records (chrom, start, end, value), record_count and provenance.
- **Example code:**

```python
import pyBigWig

from dnakit.annotation import read_pyranges1_bigwig

handle = pyBigWig.open("track.bw", "w")
handle.addHeader([("chr1", 1000), ("chr2", 500)])
handle.addEntries(["chr1", "chr1", "chr2"], [0, 200, 50], ends=[100, 260, 80], values=[1.5, 2.0, 0.25])
handle.close()
for record in read_pyranges1_bigwig("track.bw")["records"]:
    print(record["chrom"], record["start"], record["end"], record["value"])
```

- **Example output:**

```text
chr1 0 100 1.5
chr1 200 260 2.0
chr2 50 80 0.25
```

## 37) `dnakit.annotation.concatenate_pyranges1_intervals()`

- **Purpose:** Concatenate several BED files into one coordinate-sorted list of intervals.
- **Method:** Reads each file with pyranges1.read_bed and combines them with pyranges1.concat; each file is row-checked first.
- **API:** `dnakit.annotation.concatenate_pyranges1_intervals(paths, *, max_rows=100000, verbose=True)`.
- **Input:** A list of BED file paths (at most 50) and optional max_rows.
- **Output:** records, record_count, source_count and provenance.

**Example file:** `x.bed`

```text
chr1	0	100	x1	0	+
chr1	200	300	x2	0	-
chr2	50	80	x3	0	+
```


**Example file:** `y.bed`

```text
chr1	50	150	y1	0	+
chr2	60	70	y2	0	+
```

- **Example code:**

```python
from dnakit.annotation import concatenate_pyranges1_intervals

result = concatenate_pyranges1_intervals(["x.bed", "y.bed"])
print(result["record_count"], result["source_count"])
for record in result["records"]:
    print(record["chrom"], record["start"], record["end"], record["Name"])
```

- **Example output:**

```text
5 2
chr1 0 100 x1
chr1 50 150 y1
chr1 200 300 x2
chr2 50 80 x3
chr2 60 70 y2
```

## 38) `dnakit.annotation.calculate_pyranges1_cds_frame()`

- **Purpose:** Assign the reading frame of CDS features in a GTF file, grouped by one or more attributes.
- **Method:** Reads with pyranges1.read_gff, keeps Feature equal to CDS, then calls pyranges1.orfs.calculate_frame; frame is the sum of earlier CDS lengths in the group modulo 3.
- **API:** `dnakit.annotation.calculate_pyranges1_cds_frame(path, *, feature='CDS', group_by=('transcript_id',), max_rows=100000, verbose=True)`.
- **Input:** A GTF/GFF path with optional feature (default CDS), group_by (default transcript_id) and max_rows.
- **Output:** records (chrom, start, end, strand, group columns and frame), record_count and provenance.
- **Example code:**

```python
from dnakit.annotation import calculate_pyranges1_cds_frame

for record in calculate_pyranges1_cds_frame("a.gtf")["records"]:
    print(record["transcript_id"], record["start"], record["end"], record["frame"])
```

- **Example output:**

```text
t1 99 130 0
t1 199 260 1
t2 399 450 0
```

## 39) `dnakit.annotation.extend_pyranges1_orfs()`

- **Purpose:** Extend BED6 CDS intervals to complete open reading frames using a FASTA sequence.
- **Method:** Calls pyranges1.orfs.extend_orfs (needs pyfaidx); the search goes upstream for an ATG start codon and downstream for the stop codons TAG, TGA and TAA. pyfaidx writes a .fai index next to the FASTA file.
- **API:** `dnakit.annotation.extend_pyranges1_orfs(path, fasta_path, *, group_by='Name', direction=None, keep_off_bounds=False, chunk_size=900, max_rows=100000, verbose=True)`.
- **Input:** A BED6 path and a FASTA path, with optional group_by (default Name), direction (None, up or down), keep_off_bounds and chunk_size (a multiple of 3).
- **Output:** records (extended intervals), record_count and provenance.

**Example file:** `cds.bed`

```text
chr1	9	21	g	0	+
```


**Example file:** `g.fa`

```text
>chr1
GGGATGGCTGCAAAGTGATTTCCCATGCCCGCGTAAGGGCCC
```

- **Example code:**

```python
from dnakit.annotation import extend_pyranges1_orfs

up = extend_pyranges1_orfs("cds.bed", "g.fa", direction="up")["records"][0]
down = extend_pyranges1_orfs("cds.bed", "g.fa", direction="down")["records"][0]
print(up["start"], up["end"], down["start"], down["end"])
```

- **Example output:**

```text
3 21 9 36
```

## 40) `dnakit.annotation.calculate_pyranges1_forbes()`

- **Purpose:** Compute the Forbes coefficient of two BED interval sets for a genome of given chromosome sizes.
- **Method:** Calls pyranges1.stats.forbes; strand_behavior is auto, same, opposite or ignore.
- **API:** `dnakit.annotation.calculate_pyranges1_forbes(source, other, *, chromsizes, strand_behavior='auto', max_rows=100000, verbose=True)`.
- **Input:** Two BED paths, chromsizes (chromosome to length), and optional strand_behavior and max_rows.
- **Output:** forbes, strand_behavior and provenance.
- **Example code:**

```python
from dnakit.annotation import calculate_pyranges1_forbes

result = calculate_pyranges1_forbes("x.bed", "y.bed", chromsizes={"chr1": 1000, "chr2": 500})
print(result["forbes"])
```

- **Example output:**

```text
3.5573122529644268
```

## 41) `dnakit.annotation.calculate_pyranges1_mcc()`

- **Purpose:** Compute the Matthews correlation coefficient (MCC) of base-pair overlaps between labelled BED sets.
- **Method:** Calls pyranges1.stats.mcc; each row treats truth_label as the reference and prediction_label as the prediction, with TP, FP, TN, FN in base pairs.
- **API:** `dnakit.annotation.calculate_pyranges1_mcc(paths, labels, *, use_strand=False, max_rows=100000, verbose=True)`.
- **Input:** 2 to 10 BED paths, a list of unique labels of the same length, and optional use_strand and max_rows.
- **Output:** records (truth_label, prediction_label, tp, fp, tn, fn, mcc), record_count and provenance.
- **Example code:**

```python
from dnakit.annotation import calculate_pyranges1_mcc

result = calculate_pyranges1_mcc(["x.bed", "y.bed"], ["x", "y"])
for row in result["records"]:
    print(row["truth_label"], row["prediction_label"], row["tp"], row["fp"], row["tn"], row["fn"], round(row["mcc"], 4))
```

- **Example output:**

```text
x x 230 0 150 0 1.0
x y 60 50 100 170 -0.0781
y x 60 170 100 50 -0.0781
y y 110 0 270 0 1.0
```

## 42) `dnakit.annotation.calculate_pyranges1_relative_distance()`

- **Purpose:** Compute the spatial correlation between two BED interval sets as relative-distance bins.
- **Method:** Calls pyranges1.stats.relative_distance and reports reldist, count, total and fraction unchanged.
- **API:** `dnakit.annotation.calculate_pyranges1_relative_distance(source, other, *, max_rows=100000, verbose=True)`.
- **Input:** Two BED paths and optional max_rows.
- **Output:** records (reldist, count, total, fraction), record_count and provenance.
- **Example code:**

```python
from dnakit.annotation import calculate_pyranges1_relative_distance

result = calculate_pyranges1_relative_distance("x.bed", "y.bed")
for record in result["records"]:
    print(record["reldist"], record["count"], record["total"], record["fraction"])
```

- **Example output:**

```text
0.0 1 1 1.0
```

## How these functions are verified

- Test file `others/tests/unit/test_review_G18p.py`: 31 tests, all offline, no network.
- Tree and interval results are checked against hand calculations on tiny inputs (for example node distance 7 and topological distance 4 on `((A:1,B:2):1,(C:3,D:1):2);`, node depths `[0, 1, 2, 2, 3, 3, 5]`, ORF coordinates from a 42 bp sequence) and against direct upstream calls for matrices, network parsing, the Forbes coefficient and relative distances.
- Each module has a missing-package test (upstream module set to `None`, expecting `BackendUnavailableError` with the extra name in the hint), plus validation tests (booleans, node queries, unknown tips, remaining tip count, VCV symmetry and similar).

## Current scope

- 42 functions are wrapped. Not wrapped: toytree `iter_node_path` and `iter_descendant_dists` and `io.parse_newick_string_custom` duplicate the list or read wrappers; `mtree` and `io.mtree` are constructors (internal); `get_markov_model` returns only a model object (internal); `simulate_discrete_data` simulates discrete traits with any state count (not-dna; the DNA case is `simulate_nucleotide_sequences`); `utils.show` displays canvases in a browser (plot); `utils.network.parse_network` is an older prototype that prints debug output (deprecated); `get_tree_from_vcv_matrix` calls `infer_neighbor_joining_tree`, which toytree 3.0.11 does not define (broken-upstream).
- DendroPy `tree_source_iter` and `multi_tree_source_iter` raise NotImplementedError since DendroPy 4 (deprecated); `mason_gamer_kellogg_score` always raises AttributeError (broken-upstream); the four `simulate_discrete_*` functions are generic discrete-character simulators (not-dna).
- pyranges1 `read_bam` needs `bamread`, which is not installed, and upstream calls `sys.exit` when it is missing, so it is not wrapped (needs-binary:bamread). `random` is a test-data generator and `seqs.clear_kmer_memory` a cache helper (internal); `version_info` is a version report (not-dna); the remaining aliases point to the functions above (duplicate-of).
- The pyranges1 ORF functions make pyfaidx write a `.fai` index next to the FASTA file. toytree's `write()` omits internal node names, so the `name` returned by the add functions is an echo of the input.
