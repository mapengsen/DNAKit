# Optional phylogeny backends: dendropy, toytree, and pyvolve

This group of functions in `dnakit.phylogeny` consists of adapters over optional packages. [dendropy](https://github.com/jeetsukumaran/DendroPy) reads Newick or Nexus trees, reports tips and branch lengths, computes tip-to-tip distances, and compares branch lengths between two trees. [toytree](https://github.com/eaton-lab/toytree) computes layout coordinates (without drawing), quartet distances, and outgroup rooting. [pyvolve](https://github.com/sjspielman/pyvolve) simulates DNA sequences along a given tree with a fixed random seed.

These functions do not infer phylogenetic trees, draw figures, or access the network. The tree topology and branch lengths must come from the caller. The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command.

## Example data

The examples below use these small files. The trees have 4 tips and are only for demonstrating the interface.

**Example file:** `tree.nwk`

```text
((a:1,b:2):0.5,(c:0.25,d:3):1.5);
```

**Example file:** `trees.nex`

```text
#NEXUS
begin taxa;
  dimensions ntax=3;
  taxlabels X Y Z;
end;
begin trees;
  tree one = [&R] ((X:1,Y:2):3,Z:4);
  tree two = ((X:1,Z:2):3,Y:4);
end;
```

**Example file:** `tree_b.nwk`

```text
((a:1,b:1):0.5,(c:0.25,d:2):1.5);
```

**Example file:** `tree_c.nwk`

```text
((a:1,c:1):1,(b:1,d:1):1);
```

## 1) `dnakit.phylogeny.read_tree_summary()`

- **Purpose:** Read one tree from a Newick or Nexus file and report its tips and branch lengths.
- **Method:** Parse the file with dendropy's `TreeList.get`; `tree_index` selects which tree to read. Tips are sorted by label and carry their own branch lengths. `branch_lengths` lists every edge except the root's own edge, in preorder. Underscores in labels are kept as they are. A missing or non-finite length is reported as `None`.
- **API:** `dnakit.phylogeny.read_tree_summary(path, *, schema="newick", tree_index=0, max_tips=10000, max_file_bytes=10000000)`.
- **Input:** The path to an existing tree file; `schema` is `"newick"` or `"nexus"`; `max_tips` and `max_file_bytes` are explicit limits that raise an error when exceeded.
- **Output:** A dictionary with `tree_count`, `rooting`, `tips`, `branch_lengths`, `missing_branch_lengths`, and `total_branch_length`.
- **Example code:**

```python
from dnakit.phylogeny import read_tree_summary

result = read_tree_summary("tree.nwk")
print(result["result"]["tips"])
print(result["result"]["total_branch_length"])
nexus = read_tree_summary("trees.nex", schema="nexus", tree_index=1)
print(nexus["result"]["tree_count"], nexus["result"]["rooting"])
```

- **Example output:**

```text
[{'label': 'a', 'branch_length': 1.0}, {'label': 'b', 'branch_length': 2.0}, {'label': 'c', 'branch_length': 0.25}, {'label': 'd', 'branch_length': 3.0}]
8.25
2 undefined
```

## 2) `dnakit.phylogeny.calculate_patristic_distances()`

- **Purpose:** Compute the patristic distance matrix, which holds the tip-to-tip distances on one tree.
- **Method:** Use dendropy's `phylogenetic_distance_matrix` to sum the branch lengths along the path between each pair of tips. The root's own edge is not on any tip-to-tip path. Every non-root edge must have a branch length; a missing length raises an error rather than being treated as 0.
- **API:** `dnakit.phylogeny.calculate_patristic_distances(path, *, schema="newick", tree_index=0, max_taxa=200, max_file_bytes=10000000)`.
- **Input:** The tree file path; `max_taxa` limits the number of tips.
- **Output:** A dictionary with `taxa` (sorted by label) and `matrix` (a symmetric matrix with zeros on the diagonal).
- **Example code:**

```python
from dnakit.phylogeny import calculate_patristic_distances

result = calculate_patristic_distances("tree.nwk")
print(result["result"]["taxa"])
print(result["result"]["matrix"][0])
```

- **Example output:**

```text
['a', 'b', 'c', 'd']
[0.0, 3.0, 3.25, 6.0]
```

## 3) `dnakit.phylogeny.calculate_branch_length_distance()`

- **Purpose:** Compare the branch lengths of two trees that have identical tip labels.
- **Method:** Use dendropy's `euclidean_distance` to compute the Felsenstein (2004) branch-length distance. It is the square root of the summed squared differences in branch length across the unrooted splits. The two root-adjacent edges are merged onto one unrooted branch, so the position of the root does not change the value.
- **API:** `dnakit.phylogeny.calculate_branch_length_distance(path_a, path_b, *, schema="newick", tree_index=0, max_taxa=10000, max_file_bytes=10000000)`.
- **Input:** Two tree file paths. The two trees must have exactly the same tip labels. `tree_index` applies to both files.
- **Output:** A dictionary with `distance` and `tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_branch_length_distance

result = calculate_branch_length_distance("tree.nwk", "tree_b.nwk")
print(result["result"]["distance"])
```

- **Example output:**

```text
1.4142135623730951
```

## 4) `dnakit.phylogeny.calculate_node_coordinates()`

- **Purpose:** Compute the coordinates of every node of a tree in a given layout, for use by an external plotting tool.
- **Method:** Use toytree's `LinearLayout`, `CircularLayout`, or `UnrootedLayout`. `r`, `l`, `u`, and `d` are linear layouts, `c` is circular, and `unrooted` is unrooted. Coordinates are in toytree's layout units. The function does not draw any figure.
- **API:** `dnakit.phylogeny.calculate_node_coordinates(path, *, layout="r", use_edge_lengths=True, max_tips=10000, max_file_bytes=10000000)`.
- **Input:** A Newick tree file path. `layout` is one of `"r"`, `"l"`, `"u"`, `"d"`, `"c"`, or `"unrooted"`. With `use_edge_lengths=False`, every edge is laid out with unit length.
- **Output:** A dictionary with `layout`, `tip_count`, and `nodes`. Each node has `index`, `name` (`None` for internal nodes), `is_tip`, `parent`, `x`, and `y`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_node_coordinates

result = calculate_node_coordinates("tree.nwk", layout="r")
for node in result["result"]["nodes"][:4]:
    print(node["name"], node["x"], node["y"])
```

- **Example output:**

```text
a -3.0 0.0
b -2.0 1.0
c -2.75 2.0
d 0.0 3.0
```

## 5) `dnakit.phylogeny.calculate_quartet_distance()`

- **Purpose:** Count the 4-tip subsets (quartets) on which two trees with identical tip labels agree or disagree.
- **Method:** Use toytree's `distance.get_treedist_quartets`. For each 4-tip subset, a quartet counts as `shared` when both trees induce the same resolved topology on it, and as `different` when they do not. `total_quartets` equals `C(n, 4)`.
- **API:** `dnakit.phylogeny.calculate_quartet_distance(path_a, path_b, *, max_taxa=40, max_file_bytes=10000000)`.
- **Input:** Two Newick tree file paths with identical tip labels, and at least four tips.
- **Output:** A dictionary with `total_quartets`, `shared_quartets`, `different_quartets`, and `different_fraction`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_quartet_distance

result = calculate_quartet_distance("tree.nwk", "tree_c.nwk")
print(result["result"]["total_quartets"], result["result"]["different_quartets"])
```

- **Example output:**

```text
1 1
```

## 6) `dnakit.phylogeny.root_tree_on_outgroup()`

- **Purpose:** Root a Newick tree on an outgroup and write the rooted tree to a file.
- **Method:** Use toytree's `root` to place the new root on the edge above the most recent common ancestor (MRCA) of the outgroup tips. The length of that edge is split evenly between the two sides. The outgroup must be monophyletic, or an error is raised. The output file must be new; an existing file is never overwritten.
- **API:** `dnakit.phylogeny.root_tree_on_outgroup(path, outgroup, output_path, *, max_tips=10000, max_file_bytes=10000000)`.
- **Input:** The tree file path; `outgroup` is one tip label or a list of labels; `output_path` is a path to a file that does not yet exist.
- **Output:** A dictionary with `output_path`, `outgroup`, `tip_count`, and `newick`.
- **Example code:**

```python
from pathlib import Path
from dnakit.phylogeny import root_tree_on_outgroup

Path("rooted.nwk").unlink(missing_ok=True)
result = root_tree_on_outgroup("tree.nwk", "c", "rooted.nwk")
print(result["result"]["newick"].strip())
```

- **Example output:**

```text
(c:0.125,(d:3,(a:1,b:2):2):0.125);
```

## 7) `dnakit.phylogeny.simulate_nucleotide_alignment()`

- **Purpose:** Simulate a DNA alignment along a Newick tree and write it as a FASTA file.
- **Method:** Use pyvolve's `nucleotide` model. Transitions (A↔G and C↔T) have their rates multiplied by `kappa`, and all other exchangeabilities are 1. This gives the HKY85 rate matrix, with base frequencies set by `state_freqs`. Branch lengths are expected substitutions per site. A fixed `seed` makes the output reproducible.
- **API:** `dnakit.phylogeny.simulate_nucleotide_alignment(tree_path, output_path, *, length=1000, kappa=1.0, state_freqs=(0.25, 0.25, 0.25, 0.25), seed=0, max_length=100000, max_bases=10000000, max_file_bytes=10000000)`.
- **Input:** A Newick tree file path. Tip names must be unique and every edge must have a length. `output_path` must not exist yet. `length` is the number of sites. `state_freqs` gives the frequencies of A, C, G, and T, which must sum to 1.
- **Output:** A dictionary with `output_path`, `taxa`, `sequence_count`, `site_count`, `model`, and `seed`.
- **Example code:**

```python
from pathlib import Path
from dnakit.phylogeny import simulate_nucleotide_alignment

Path("alignment.fasta").unlink(missing_ok=True)
result = simulate_nucleotide_alignment("tree.nwk", "alignment.fasta", length=12, seed=1)
print(result["result"]["taxa"], result["result"]["site_count"])
print(Path("alignment.fasta").read_text().splitlines()[:2])
```

- **Example output:**

```text
['a', 'b', 'c', 'd'] 12
['>a', 'TTCCATCAGGTT']
```

## 8) `dnakit.phylogeny.simulate_nucleotide_sequences()`

- **Purpose:** Run the same simulation as the previous function, but return the tip sequences directly instead of writing a file.
- **Method:** Use the same HKY85 rate matrix and seeding rules, and read the tip sequences from pyvolve's `get_sequences`. An error is raised if `tips × length` exceeds `max_bases`.
- **API:** `dnakit.phylogeny.simulate_nucleotide_sequences(tree_path, *, length=1000, kappa=1.0, state_freqs=(0.25, 0.25, 0.25, 0.25), seed=0, max_length=100000, max_bases=10000000, max_file_bytes=10000000)`.
- **Input:** The same as `simulate_nucleotide_alignment()`, without an output path.
- **Output:** A dictionary with `taxa`, `site_count`, `sequences` (a mapping from tip label to sequence), `model`, and `seed`.
- **Example code:**

```python
from dnakit.phylogeny import simulate_nucleotide_sequences

result = simulate_nucleotide_sequences("tree.nwk", length=12, kappa=2.0, seed=1)
print(result["result"]["sequences"]["a"])
print(result["result"]["sequences"]["d"])
```

- **Example output:**

```text
TTCCCTCAGGTT
TCGAGTTTGCTA
```

## 9) `dnakit.phylogeny.simulate_substitution_counts()`

- **Purpose:** Simulate sequences with a fixed seed and count the substitutions that occur on each branch.
- **Method:** Use pyvolve's exact Gillespie algorithm (`algorithm=1`). Matrix exponentiation is not used, because it undercounts substitutions. Branch names are tip names for tip branches and pyvolve-generated `internalNode` names for internal branches. The root edge is not counted.
- **API:** `dnakit.phylogeny.simulate_substitution_counts(tree_path, *, length=1000, kappa=1.0, state_freqs=(0.25, 0.25, 0.25, 0.25), seed=0, max_length=100000, max_bases=10000000, max_file_bytes=10000000)`.
- **Input:** The same as `simulate_nucleotide_alignment()`.
- **Output:** A dictionary with `taxa`, `site_count`, `branch_substitution_counts`, `total_substitutions`, `model`, and `seed`.
- **Example code:**

```python
from dnakit.phylogeny import simulate_substitution_counts

result = simulate_substitution_counts("tree.nwk", length=1000, seed=1)
print(result["result"]["branch_substitution_counts"])
print(result["result"]["total_substitutions"])
```

- **Example output:**

```text
{'a': 1038, 'b': 2025, 'c': 263, 'd': 2968, 'internalNode1': 498, 'internalNode2': 1538}
8330
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-dendropy]"   # dendropy, for read_tree_summary and the tree distance functions
python -m pip install "dnakit[bio-toytree]"    # toytree, for layout coordinates, quartet distance, and rooting
python -m pip install "dnakit[bio-pyvolve]"    # pyvolve, for DNA sequence simulation
```

Project links: [dendropy](https://github.com/jeetsukumaran/DendroPy), [toytree](https://github.com/eaton-lab/toytree), [pyvolve](https://github.com/sjspielman/pyvolve).

## How these functions are verified

The test file is `others/tests/unit/test_phylogeny_backends.py`, which has 32 tests. They check the following:

- dendropy tips, branch lengths, Nexus tree selection, and the patristic matrix. The patristic matrix is also checked against the root-distance identity (distance = depth(a) + depth(b) − 2 × depth(MRCA)). The branch-length distance matches a direct call to dendropy's `euclidean_distance`, a hand-calculated value of √2, and a value of 0 when only the root placement differs.
- toytree linear layout coordinates match hand calculations and a direct `LinearLayout` call. Quartet counts match a hand count on 5 tips (3 shared, 2 different). Rooted Newick output matches a hand calculation, and the total branch length is preserved.
- pyvolve's Jukes–Cantor divergence matches the analytical formula. The transition fraction matches kappa/(kappa + 2). Per-branch substitution counts agree with the expected value of branch length × number of sites, within sampling error. Base composition matches the requested frequencies. Output is reproducible for a fixed seed, and the file and in-memory outputs are identical.
- A missing optional package raises an install hint. Missing branch lengths, duplicate tip names, invalid kappa, frequencies, or seeds, and an output file that already exists are all rejected.

## Current scope

This group wraps tree reading, distances, layout coordinates, rooting, and DNA sequence simulation. The following are not wrapped:

- dendropy tree simulation (birth-death, coalescent) and character simulation, which are outside the DNA sequence scope of this group.
- Robinson–Foulds distance in dendropy and toytree. It computes the same quantity as `dnakit.bioinformatics.phylogeny.robinson_foulds()`, so it is not wrapped again.
- toytree plotting functions and toytree tree inference.
- pyvolve codon and amino acid models, custom rate matrices, and PAML or mutation-selection models.
- Other pyvolve sequence formats (Phylip, Nexus, and so on). Only FASTA is written at present.
