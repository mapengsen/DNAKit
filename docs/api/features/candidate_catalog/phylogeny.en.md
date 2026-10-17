# Phylogeny

This category contains **8** documented capabilities; all 8 compute a result natively in pure Python from the supplied input (Newick trees, distance matrices or DNA alignments), with no external program or optional library. Only `infer` can return `status: "not_implemented"`, and only when `options.method` asks for a maximum-likelihood or Bayesian method (use `infer_fast_tree` / `infer_ml_tree`, which call external programs). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) annotate

- **Purpose:** Annotates every node of a Newick tree with structural fields (depth in edges, number of leaves, distance from the root, clade membership) and, optionally, with per-leaf metadata.

- **Method:** Parses the Newick tree in pure Python, numbers the nodes in preorder and labels unnamed internal nodes `N<preorder index>`. For each node it reports the parent id, child count, number of leaves below, depth in edges, its own branch length and the root distance (the sum of branch lengths from the root, `None` if a length on the path is missing). Internal nodes with at most `max_clade` leaves also list the sorted names of their leaves. With `metadata`, each node lists the distinct metadata values (as strings) of the leaves below it and whether they are homogeneous (exactly one distinct value). Nothing is inferred along the tree (no ancestral-state reconstruction).

- **API:** `dnakit.phylogeny.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): Newick text or a path to a Newick file. `options` (optional): `metadata` (mapping leaf name -> value, default none; a name that is not a leaf raises an error), `max_clade` (int 1-1,000,000, default 50; internal nodes with more leaves do not list their clade).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `nodes` (one row per node in preorder with `id`, `name`, `parent`, `is_leaf`, `child_count`, `leaf_count`, `depth_edges`, `branch_length`, `root_distance`, `clade` for small internal nodes, and `metadata_values` / `metadata_homogeneous` when metadata is given), `annotated_newick` (the tree with the generated labels) and `internal_nodes_labelled` (how many labels were generated).

- **Example code:**

```python
import json

from dnakit.phylogeny import annotate

result = annotate(
    "((A:1,B:2):3,C:4);",
    options={"metadata": {"A": "x", "B": "x", "C": "y"}},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "phylogeny",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 18,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a96a137a9d8a4755cd277cb9193377aab112cd1f3719f4065ede30482e9bbd5b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "nodes": [
      {
        "id": 0,
        "name": "N0",
        "parent": null,
        "is_leaf": false,
        "child_count": 2,
        "leaf_count": 3,
        "depth_edges": 0,
        "branch_length": null,
        "root_distance": 0.0,
        "clade": ["A", "B", "C"],
        "metadata_values": ["x", "y"],
        "metadata_homogeneous": false
      },
      {
        "id": 1,
        "name": "N1",
        "parent": 0,
        "is_leaf": false,
        "child_count": 2,
        "leaf_count": 2,
        "depth_edges": 1,
        "branch_length": 3.0,
        "root_distance": 3.0,
        "clade": ["A", "B"],
        "metadata_values": ["x"],
        "metadata_homogeneous": true
      },
      {
        "id": 2,
        "name": "A",
        "parent": 1,
        "is_leaf": true,
        "child_count": 0,
        "leaf_count": 1,
        "depth_edges": 2,
        "branch_length": 1.0,
        "root_distance": 4.0,
        "metadata_values": ["x"],
        "metadata_homogeneous": true
      },
      {
        "id": 3,
        "name": "B",
        "parent": 1,
        "is_leaf": true,
        "child_count": 0,
        "leaf_count": 1,
        "depth_edges": 2,
        "branch_length": 2.0,
        "root_distance": 5.0,
        "metadata_values": ["x"],
        "metadata_homogeneous": true
      },
      {
        "id": 4,
        "name": "C",
        "parent": 0,
        "is_leaf": true,
        "child_count": 0,
        "leaf_count": 1,
        "depth_edges": 1,
        "branch_length": 4.0,
        "root_distance": 4.0,
        "metadata_values": ["y"],
        "metadata_homogeneous": true
      }
    ],
    "annotated_newick": "((A:1,B:2)N1:3,C:4)N0;",
    "internal_nodes_labelled": 2
  }
}
```

## 2) build

- **Purpose:** Builds a tree from a pairwise distance matrix by neighbour-joining (Saitou and Nei 1987) or UPGMA (average linkage).

- **Method:** Neighbour-joining repeatedly joins the pair that minimises the Q criterion `(m-2)*d(i,j) - sum_k d(i,k) - sum_k d(j,k)`, computes the two branch lengths from the net divergences and updates distances to the new node as `(d(i,k) + d(j,k) - d(i,j))/2`; the last three nodes are joined at a trifurcation. The result is unrooted (written with a trifurcating root) and branch lengths can be negative; they are not clipped. UPGMA repeatedly merges the two closest clusters, sets the new node height to half their distance and averages the distances to other clusters weighted by cluster size, so the tree is rooted and ultrametric. Ties are broken by input order. The matrix must be square, symmetric, with a zero diagonal and non-negative finite entries, unique names and at most 400 taxa. Pure Python; no external program.

- **API:** `dnakit.phylogeny.build(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `{"names": [...], "matrix": [[...]]}` (without `names` the taxa are called t1, t2, ...) or a nested mapping `{a: {b: distance}}`. `options` (optional): `method` (`"nj"` default or `"upgma"`).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `newick` (the tree), `method`, `rooted` (true only for UPGMA), `leaf_count`, `total_branch_length`, `negative_branch_count` and an explanatory `note`.

- **Example code:**

```python
import json

from dnakit.phylogeny import build

matrix = {
    "names": ["a", "b", "c", "d", "e"],
    "matrix": [
        [0, 5, 9, 9, 8],
        [5, 0, 10, 10, 9],
        [9, 10, 0, 8, 7],
        [9, 10, 8, 0, 3],
        [8, 9, 7, 3, 0],
    ],
}
result = build(matrix)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "phylogeny",
  "operation": "build",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 137,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f8d804d5c145007c60bd2e00c2dd6c386da3e2bdccb29b0fececaf0e04ada64f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "newick": "(d:2,e:1,(c:4,(a:2,b:3):3):2);",
    "method": "nj",
    "rooted": false,
    "leaf_count": 5,
    "total_branch_length": 17.0,
    "negative_branch_count": 0,
    "note": "UPGMA gives a rooted ultrametric tree; neighbour-joining gives an unrooted tree drawn with a trifurcating root and may contain negative branch lengths"
  }
}
```

## 3) calculate

- **Purpose:** Calculates tree statistics or a patristic distance matrix from a Newick tree, or a pairwise distance matrix from an aligned set of DNA sequences.

- **Method:** `what="tree_statistics"` reports leaf, internal-node and node counts, the leaf names, root degree, whether the tree is rooted binary, the number of polytomies (a trifurcating root is not counted), maximum and mean leaf depth in edges, the Sackin index (sum of leaf depths in edges), total branch length, the minimum, maximum and mean root-to-tip distance, whether the tree is ultrametric (tip distances equal within a relative 1e-9), treeness (internal branch length / total branch length) and, for rooted binary trees only, the Colless index (sum over internal nodes of the absolute difference in leaf numbers of the two subtrees). `what="patristic_matrix"` returns the leaf-to-leaf path length through the lowest common ancestor (all branches need lengths; at most 1,000 leaves). `what="distance_matrix"` computes, for every pair of sequences, the p distance (mismatches / compared sites) or the Jukes-Cantor distance `d = -3/4 ln(1 - 4p/3)`; only A/C/G/T are compared, so gaps and other symbols are skipped pairwise (U is read as T, case is ignored). A saturated pair (p >= 0.75) or a pair with no comparable site raises an error unless `saturation_distance` is set.

- **API:** `dnakit.phylogeny.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): a Newick tree (text or path) for the tree results, or an alignment (FASTA text/path, `{name: sequence}` or `[{"name", "sequence"}]`; equal lengths, at most 400 sequences) for the distance matrix. `options` (optional): `what` (`"auto"` default: tree statistics for a Newick string, distance matrix for an alignment; or `"tree_statistics"`, `"patristic_matrix"`, `"distance_matrix"`), `distance` (`"jc69"` default or `"p"`), `saturation_distance` (positive number used instead of raising for saturated or undefined pairs, default none).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains Always `what`. For `tree_statistics`: the statistics listed above (`leaf_count`, `internal_node_count`, `node_count`, `leaves`, `root_degree`, `rooted_binary`, `polytomy_count`, `max_depth_edges`, `mean_leaf_depth_edges`, `sackin_index`, `has_branch_lengths`, `all_branches_have_lengths`, `negative_branch_count`, `total_branch_length`, and, when all branches have lengths, `max_root_to_tip`, `min_root_to_tip`, `mean_root_to_tip`, `ultrametric`, `treeness`, plus `colless_index` for rooted binary trees). For `patristic_matrix`: `names`, `matrix`, `distance` (`"patristic"`). For `distance_matrix`: `names`, `matrix`, `distance`, `compared_sites` (matrix of compared site counts) and `alignment_length`.

- **Example code:**

```python
import json

from dnakit.phylogeny import calculate

result = calculate("((A:1,B:2)X:3,C:4);")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "phylogeny",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 19,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "90a78a93da03a005ce42758202cb84a0e7b01040de38aa3be6ccff6b31fc785c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "what": "tree_statistics",
    "leaf_count": 3,
    "internal_node_count": 2,
    "node_count": 5,
    "leaves": ["A", "B", "C"],
    "root_degree": 2,
    "rooted_binary": true,
    "polytomy_count": 0,
    "max_depth_edges": 2,
    "mean_leaf_depth_edges": 1.6666666666666667,
    "sackin_index": 5,
    "has_branch_lengths": true,
    "all_branches_have_lengths": true,
    "negative_branch_count": 0,
    "total_branch_length": 10.0,
    "max_root_to_tip": 5.0,
    "min_root_to_tip": 4.0,
    "mean_root_to_tip": 4.333333333333333,
    "ultrametric": false,
    "treeness": 0.3,
    "colless_index": 1
  }
}
```

## 4) compare

- **Purpose:** Computes the Robinson-Foulds distance between two trees on the same leaf set.

- **Method:** Collects the non-trivial splits of each tree (a split is a bipartition of the leaves with at least two leaves on each side; with `rooted=True`, the clades with at least two leaves instead) and reports RF = size of the symmetric difference of the two split sets (Robinson and Foulds 1981). The normalised value is RF / (splits_first + splits_second), not RF / maximum RF; the maximum RF of binary trees (`2(n-3)` unrooted, `2(n-2)` rooted) is returned separately. Branch lengths and internal labels are ignored. Trees with different leaf sets raise an error.

- **API:** `dnakit.phylogeny.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): `[newick_a, newick_b]`, `{"a": ..., "b": ...}` (also `left`/`right` or `tree1`/`tree2`), or one tree plus `options["other"]`; trees are Newick text or paths. `options` (optional): `rooted` (bool, default False: unrooted splits), `other` (the second tree).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `robinson_foulds`, `normalized_rf`, `shared_splits`, `splits_only_first`, `splits_only_second`, `splits_first`, `splits_second`, `max_rf_binary`, `leaf_count`, `rooted` and a `definition` string.

- **Example code:**

```python
import json

from dnakit.phylogeny import compare

result = compare(["((A,B),(C,D),(E,F));", "(((A,B),C),(D,(E,F)));"])
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "phylogeny",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 50,
    "sequence_bases": 12,
    "lengths": [6, 6],
    "sha256": "e6b7a61850243eb35ecb4e40acdf77a9dc6d014e0a4043ba5e351a3a04fcd760"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "robinson_foulds": 2,
    "normalized_rf": 0.3333333333333333,
    "shared_splits": 2,
    "splits_only_first": 1,
    "splits_only_second": 1,
    "splits_first": 3,
    "splits_second": 3,
    "max_rf_binary": 6,
    "leaf_count": 6,
    "rooted": false,
    "definition": "RF = |S1 symmetric-difference S2| over non-trivial splits (bipartitions); normalized = RF / (|S1| + |S2|)"
  }
}
```

## 5) export

- **Purpose:** Serialises a Newick tree to normalised Newick, NEXUS, nested JSON or an edge list.

- **Method:** Parses the tree and writes it again in the requested format. Newick output re-quotes labels when needed and writes branch lengths with 12 significant digits; NEXUS contains a TAXA block (`DIMENSIONS NTAX`, `TAXLABELS`) and a TREES block with the tree (marked `[&R]`); `json` is a nested `{name, branch_length, children}` structure; `edgelist` lists `parent`, `child`, `child_name` and `length` for every edge using preorder node numbers. If `output` is given the text is also written to that file (JSON as indented text). Pure Python.

- **API:** `dnakit.phylogeny.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): Newick text or a path to a Newick file. `options` (optional): `format` (`"newick"` default, `"nexus"`, `"json"` or `"edgelist"`), `output` (optional file path to write; default none).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `format`, `content` (a string for Newick and NEXUS, an object for `json`, a list for `edgelist`) and `output_path` (the file written, or null).

- **Example code:**

```python
import json

from dnakit.phylogeny import export

result = export("((A:1,B:2):3,'C d':4);", options={"format": "nexus"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "phylogeny",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 22,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "67603f444934d92f5102461f5902da759227011cabde49307f86bde1a39967c2"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "nexus",
    "content": "#NEXUS\n\nBEGIN TAXA;\n  DIMENSIONS NTAX=3;\n  TAXLABELS A B 'C d';\nEND;\n\nBEGIN TREES;\n  TREE tree1 = [&R] ((A:1,B:2):3,'C d':4);\nEND;\n",
    "output_path": null
  }
}
```

## 6) infer

- **Purpose:** Infers a distance-based tree (neighbour-joining or UPGMA on Jukes-Cantor or p distances) from aligned DNA sequences.

- **Method:** Computes the pairwise distance matrix of the alignment as in `calculate` (p distance, or Jukes-Cantor `d = -3/4 ln(1 - 4p/3)`, comparing only A/C/G/T sites pairwise) and builds the tree with neighbour-joining or UPGMA as in `build`. Maximum-likelihood and Bayesian inference are not available here: for `options.method` other than `nj` / `upgma` the function returns `not_implemented` and points to `infer_fast_tree` / `infer_ml_tree`, which call external programs (FastTree, IQ-TREE, RAxML, MrBayes). Saturated pairs (p >= 0.75 under Jukes-Cantor) raise an error unless `saturation_distance` is set.

- **API:** `dnakit.phylogeny.infer(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): FASTA text/path, `{name: sequence}` or `[{"name", "sequence"}]`; equal lengths, at most 400 sequences. `options` (optional): `method` (`"nj"` default or `"upgma"`), `distance` (`"jc69"` default or `"p"`), `saturation_distance` (distance used for saturated or undefined pairs instead of raising, default none).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains the `build` fields (`newick`, `method`, `rooted`, `leaf_count`, `total_branch_length`, `negative_branch_count`, `note`) plus `distance`, `distance_matrix` and `names`. For a method other than `nj` / `upgma` the envelope is instead `not_implemented`.

- **Example code:**

```python
import json

from dnakit.phylogeny import infer

alignment = {
    "A": "ACGTACGTACGT",
    "B": "ACGTACGTACGA",
    "C": "TGCATGGTACGT",
    "D": "TGCATGGTACGA",
}
result = infer(alignment, options={"distance": "p"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "phylogeny",
  "operation": "infer",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 84,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "691762599c2228308bdd799aaa8d442e73fff378d08b81c0ea13b06176d43582"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "newick": "(C:0.0416666666667,D:0.0416666666667,(A:0.0416666666667,B:0.0416666666667):0.458333333333);",
    "method": "nj",
    "rooted": false,
    "leaf_count": 4,
    "total_branch_length": 0.6250000000000001,
    "negative_branch_count": 0,
    "note": "UPGMA gives a rooted ultrametric tree; neighbour-joining gives an unrooted tree drawn with a trifurcating root and may contain negative branch lengths",
    "distance": "p",
    "distance_matrix": [
      [0.0, 0.08333333333333333, 0.5, 0.5833333333333334],
      [0.08333333333333333, 0.0, 0.5833333333333334, 0.5],
      [0.5, 0.5833333333333334, 0.0, 0.08333333333333333],
      [0.5833333333333334, 0.5, 0.08333333333333333, 0.0]
    ],
    "names": ["A", "B", "C", "D"]
  }
}
```

## 7) plot

- **Purpose:** Draws a Newick tree as a left-to-right rectangular SVG phylogram or cladogram.

- **Method:** Leaves are placed in tree order on equally spaced rows, each internal node is centred on its children, and horizontal positions are the root distance (phylogram, `branch_lengths=True`) or the depth in edges (cladogram, `branch_lengths=False`), scaled to the plot width. The tree is never reordered. Leaf names are written as SVG text. A phylogram of a tree with missing branch lengths raises an error. At most 2,000 leaves; pure Python, no plotting library. If `output` is given the SVG is also written to that file.

- **API:** `dnakit.phylogeny.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): Newick text or path (at most 2,000 leaves). `options` (optional): `branch_lengths` (bool, default True), `width` (int pixels 200-4000, default 640), `output` (optional file path for the SVG; default none).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `svg` (the SVG text), `leaf_count`, `branch_lengths` (the mode used) and `output_path` (the file written, or null).

- **Example code:**

```python
import json

from dnakit.phylogeny import plot

result = plot("((A:1,B:2):3,C:4);", options={"width": 300})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "phylogeny",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 18,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a96a137a9d8a4755cd277cb9193377aab112cd1f3719f4065ede30482e9bbd5b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"300\" height=\"74\" viewBox=\"0 0 300 74\">\n<g stroke=\"#333\" stroke-width=\"1\" fill=\"none\">\n<line x1=\"10.0\" y1=\"28.0\" x2=\"10.0\" y2=\"55.0\"/>\n<line x1=\"172.0\" y1=\"19.0\" x2=\"172.0\" y2=\"37.0\"/>\n<line x1=...",
    "leaf_count": 3,
    "branch_lengths": true,
    "output_path": null
  }
}
```

(Long arrays and text in the example output are shortened with `...`.)

## 8) test

- **Purpose:** Estimates nonparametric bootstrap support for the clades of a distance-based tree built from an alignment.

- **Method:** Builds the reference tree from the original alignment (neighbour-joining or UPGMA on Jukes-Cantor or p distances). Each replicate resamples alignment columns with replacement (Felsenstein 1985), recomputes the distances, rebuilds the tree and counts which reference splits reappear (unrooted splits for `nj`, rooted clades for `upgma`). Support is the fraction of replicates in `[0, 1]`; the returned Newick carries the percentage (rounded) as the label of each internal node. A fixed `seed` gives identical results. Saturated or undefined pairs in a replicate get `saturation_distance` instead of an error. The work (pairs x distinct column patterns x replicates) is capped. This is a bootstrap of a distance tree only; it does not bootstrap maximum-likelihood trees.

- **API:** `dnakit.phylogeny.test(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` (required): FASTA text/path, `{name: sequence}` or `[{"name", "sequence"}]`; equal lengths, at most 400 sequences. `options` (optional): `replicates` (int 1-10,000, default 100), `seed` (int, default 0), `method` (`"nj"` default or `"upgma"`), `distance` (`"jc69"` default or `"p"`), `saturation_distance` (distance for saturated or undefined pairs in the replicates; default 10.0 for `jc69` and 1.0 for `p`).

- **Output:** Standard envelope: `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. `result` contains `newick_with_support` (the reference tree with percentage labels), `clades` (one entry per reference split with its leaf set `clade`, listed as the side that does not contain the alphabetically first leaf for unrooted trees, and `support`), `replicates`, `seed`, `method`, `distance` and `mean_support`.

- **Example code:**

```python
import json

from dnakit.phylogeny import test

alignment = {
    "A": "AAAAAAGGGG",
    "B": "AAAAAATTTT",
    "C": "CCCCCCGGGG",
    "D": "CCCCCCTTTT",
}
result = test(alignment, options={"replicates": 100, "seed": 1, "distance": "p"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "phylogeny",
  "operation": "test",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 76,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "79dd2c6923081105920d7bb30423f9f3c0e984f0918d4115020eb9b7ca4fbf39"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "newick_with_support": "(C:0.2,D:0.2,(A:0.2,B:0.2)82:0.4);",
    "clades": [
      {
        "clade": ["C", "D"],
        "support": 0.82
      }
    ],
    "replicates": 100,
    "seed": 1,
    "method": "nj",
    "distance": "p",
    "mean_support": 0.82
  }
}
```
