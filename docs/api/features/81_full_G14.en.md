# DendroPy, toytree, and pyvolve phylogeny functions (G14)

This group of functions in `dnakit.phylogeny` consists of wrappers over optional packages, in two modules. `dnakit.phylogeny` uses DendroPy for population-genetic statistics on DNA alignments (segregating sites, mean pairwise differences, nucleotide diversity, Tajima's D, Watterson's theta, derived states and frequency spectra), and for tree comparison and tree-shape statistics. `dnakit.phylogeny` uses toytree for rooting, editing, distances, and tip matrices. It returns Newick or numbers and draws nothing.

pyvolve adds no wrappers: its `read_tree` is already used by the existing `simulate_nucleotide_*` functions, and `print_tree` only renders text.

All functions work only on sequences or trees given explicitly. None of them uses the network or calls an external program. The dependencies are optional and are imported only when the matching function is called. A missing package raises an error with the `dnakit[bio-dendropy]` or `dnakit[bio-toytree]` install command. These functions are currently imported from their module paths and are not yet exported from `dnakit.phylogeny`.

## Example data

Most examples below put the alignment or tree directly in the code. `tree_t1.nwk` is used only to show file-path input.

**Example file:** `tree_t1.nwk`

```text
((A:1,B:1):1,(C:2,D:2):0.5);
```

## 1) `dnakit.phylogeny.count_segregating_sites()`

- **Purpose:** Count polymorphic sites in a DNA alignment.
- **Method:** Calls `dendropy.calculate.popgenstat.num_segregating_sites`. With `ignore_uncertain=True`, a site that has an ambiguous or gap symbol in any sequence is skipped by dendropy's rule.
- **API:** `dnakit.phylogeny.count_segregating_sites(sequences, *, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`.
- **Input:** A dictionary from sequence ID to sequence, or a `DNA` batch. At least two sequences of equal length. `max_sequences` and `max_length` are explicit limits.
- **Output:** `result` has `segregating_sites`, `sequence_count`, `site_count`, and `ignore_uncertain`; `provenance` records the dendropy version.
- **Example code:**

```python
from dnakit.phylogeny import count_segregating_sites

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(count_segregating_sites(alignment)["result"])
```

- **Example output:**

```text
{'sequence_count': 4, 'site_count': 10, 'ignore_uncertain': True, 'segregating_sites': 3}
```

## 2) `dnakit.phylogeny.calculate_mean_pairwise_differences()`

- **Purpose:** Compute the mean number of pairwise differences (Tajima's k) of an alignment.
- **Method:** Calls `dendropy.calculate.popgenstat.average_number_of_pairwise_differences`. It counts differing sites for every pair of sequences and averages over all pairs.
- **API:** `dnakit.phylogeny.calculate_mean_pairwise_differences(sequences, *, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`.
- **Input:** The same sequence dictionary or `DNA` batch as `count_segregating_sites()`.
- **Output:** `result.mean_pairwise_differences` holds the mean; other fields match `count_segregating_sites()`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_mean_pairwise_differences

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(calculate_mean_pairwise_differences(alignment)["result"])
```

- **Example output:**

```text
{'sequence_count': 4, 'site_count': 10, 'ignore_uncertain': True, 'mean_pairwise_differences': 1.5}
```

## 3) `dnakit.phylogeny.calculate_nucleotide_diversity()`

- **Purpose:** Compute nucleotide diversity (pi) of an alignment.
- **Method:** Calls `dendropy.calculate.popgenstat.nucleotide_diversity`. The value comes from dendropy and is not recomputed by the wrapper.
- **API:** `dnakit.phylogeny.calculate_nucleotide_diversity(sequences, *, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`.
- **Input:** The same as `count_segregating_sites()`.
- **Output:** `result.nucleotide_diversity` holds pi; other fields are the same.
- **Example code:**

```python
from dnakit.phylogeny import calculate_nucleotide_diversity

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(calculate_nucleotide_diversity(alignment)["result"])
```

- **Example output:**

```text
{'sequence_count': 4, 'site_count': 10, 'ignore_uncertain': True, 'nucleotide_diversity': 0.15}
```

## 4) `dnakit.phylogeny.calculate_tajimas_d()`

- **Purpose:** Compute Tajima's D for an alignment.
- **Method:** First counts segregating sites S with dendropy. When S is 0, D is undefined and the value is `None`. Otherwise it calls `dendropy.calculate.popgenstat.tajimas_d`.
- **API:** `dnakit.phylogeny.calculate_tajimas_d(sequences, *, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`.
- **Input:** The same as `count_segregating_sites()`.
- **Output:** `result.tajimas_d` is D or `None`; `result.segregating_sites` is S.
- **Example code:**

```python
from dnakit.phylogeny import calculate_tajimas_d

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(calculate_tajimas_d(alignment)["result"])
```

- **Example output:**

```text
{'sequence_count': 4, 'site_count': 10, 'ignore_uncertain': True, 'tajimas_d': -0.754451077653, 'segregating_sites': 3}
```

## 5) `dnakit.phylogeny.calculate_wattersons_theta()`

- **Purpose:** Compute Watterson's theta for an alignment.
- **Method:** Calls `dendropy.calculate.popgenstat.wattersons_theta`. The result is S / a_n with a_n = 1 + 1/2 + ... + 1/(n-1), not divided by the number of sites.
- **API:** `dnakit.phylogeny.calculate_wattersons_theta(sequences, *, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`.
- **Input:** The same as `count_segregating_sites()`.
- **Output:** `result.wattersons_theta` is theta_W; `result.convention_note` states the denominator.
- **Example code:**

```python
from dnakit.phylogeny import calculate_wattersons_theta

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(calculate_wattersons_theta(alignment)["result"])
```

- **Example output:**

```text
{'sequence_count': 4, 'site_count': 10, 'ignore_uncertain': True, 'wattersons_theta': 1.636363636364, 'convention_note': 'segregating sites divided by a_n = sum(1/i) for i in 1..n-1'}
```

## 6) `dnakit.phylogeny.build_derived_state_matrix()`

- **Purpose:** Mark each ingroup site 0 (ancestral state) or 1 (derived state) against an ancestral sequence.
- **Method:** The sequence named by `ancestral_id` is the reference and is removed from the ingroup. Calls `dendropy.calculate.popgenstat.derived_state_matrix`; rows are sorted by ID.
- **API:** `dnakit.phylogeny.build_derived_state_matrix(sequences, *, ancestral_id, ignore_uncertain=True, max_sequences=1000, max_length=100000, verbose=True)`.
- **Input:** The same as `count_segregating_sites()`, plus `ancestral_id`, which must be one of the input IDs.
- **Output:** `result.rows` is a list of `{id, states}` with 0/1 strings; `result.ancestral_id` echoes the ancestral ID.
- **Example code:**

```python
from dnakit.phylogeny import build_derived_state_matrix

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
print(build_derived_state_matrix(alignment, ancestral_id="A")["result"]["rows"])
```

- **Example output:**

```text
[{'id': 'B', 'states': '0000000100'}, {'id': 'C', 'states': '0000100100'}, {'id': 'D', 'states': '0001000100'}]
```

## 7) `dnakit.phylogeny.calculate_unfolded_site_frequency_spectrum()`

- **Purpose:** Compute the unfolded site frequency spectrum of the ingroup relative to an ancestral sequence.
- **Method:** Calls `dendropy.calculate.popgenstat.unfolded_site_frequency_spectrum` with `pad=True`. The output has one entry per derived count from 0 to the ingroup size; entry k counts sites with exactly k derived alleles. The default is `ignore_uncertain=False`.
- **API:** `dnakit.phylogeny.calculate_unfolded_site_frequency_spectrum(sequences, *, ancestral_id, ignore_uncertain=False, max_sequences=1000, max_length=100000, verbose=True)`.
- **Input:** The same as `build_derived_state_matrix()`.
- **Output:** `result.site_frequency_spectrum` is a list; `result.ingroup_count` is the number of ingroup sequences.
- **Example code:**

```python
from dnakit.phylogeny import calculate_unfolded_site_frequency_spectrum

alignment = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "ACGTTCGAAC", "D": "ACGAACGAAC"}
result = calculate_unfolded_site_frequency_spectrum(alignment, ancestral_id="A")
print(result["result"]["site_frequency_spectrum"])
```

- **Example output:**

```text
[7, 2, 0, 1]
```

## 8) `dnakit.phylogeny.calculate_weighted_robinson_foulds_distance()`

- **Purpose:** Compute the weighted Robinson-Foulds distance between two trees on the same tips.
- **Method:** Both trees are read into one TaxonNamespace. Tip labels must match exactly, and every non-root edge needs a length. Calls `dendropy.calculate.treecompare.weighted_robinson_foulds_distance`.
- **API:** `dnakit.phylogeny.calculate_weighted_robinson_foulds_distance(tree_a, tree_b, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** Two Newick strings or paths to one-tree Newick files. `max_tips` and `max_file_bytes` are explicit limits.
- **Output:** `result.distance` is the distance; `result.tip_count` is the number of tips.
- **Example code:**

```python
from dnakit.phylogeny import calculate_weighted_robinson_foulds_distance

print(calculate_weighted_robinson_foulds_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'distance': 5.0}
```

## 9) `dnakit.phylogeny.count_false_positive_and_negative_splits()`

- **Purpose:** Count splits that appear in only one of two trees.
- **Method:** Calls `dendropy.calculate.treecompare.false_positives_and_negatives`. False positives are splits present only in the comparison tree; false negatives are splits present only in the reference tree.
- **API:** `dnakit.phylogeny.count_false_positive_and_negative_splits(reference_tree, comparison_tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A reference tree and a comparison tree in the same format as `calculate_weighted_robinson_foulds_distance()`.
- **Output:** `result.false_positives`, `result.false_negatives`, and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import count_false_positive_and_negative_splits

print(count_false_positive_and_negative_splits("((A,B),(C,D));", "((A,C),(B,D));")["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'false_positives': 1, 'false_negatives': 1}
```

## 10) `dnakit.phylogeny.calculate_colless_tree_imbalance()`

- **Purpose:** Compute the Colless imbalance of a bifurcating tree.
- **Method:** Calls `dendropy.calculate.treemeasure.colless_tree_imbalance`. A non-bifurcating tree raises `ConfigurationError`. `normalize` is `"max"` (default), `"yule"`, `"pda"`, or `None` for no normalization.
- **API:** `dnakit.phylogeny.calculate_colless_tree_imbalance(tree, *, normalize="max", max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path.
- **Output:** `result.colless_imbalance`, `result.normalize`, and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_colless_tree_imbalance

print(calculate_colless_tree_imbalance("((A,B),(C,(D,E)));", normalize=None)["result"])
```

- **Example output:**

```text
{'tip_count': 5, 'colless_imbalance': 2.0, 'normalize': 'none'}
```

## 11) `dnakit.phylogeny.calculate_sackin_index()`

- **Purpose:** Compute the Sackin index: the sum over tips of their number of ancestors.
- **Method:** Calls `dendropy.calculate.treemeasure.sackin_index`. `normalize="tips"` (default) divides by the tip count, which equals N-bar; `"yule"` and `"pda"` normalize to those models; `None` leaves it raw.
- **API:** `dnakit.phylogeny.calculate_sackin_index(tree, *, normalize="tips", max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path.
- **Output:** `result.sackin_index`, `result.normalize`, and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_sackin_index

print(calculate_sackin_index("((A,B),(C,D));", normalize=None)["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'sackin_index': 8.0, 'normalize': 'none'}
```

## 12) `dnakit.phylogeny.calculate_b1_index()`

- **Purpose:** Compute the B1 tree-shape statistic.
- **Method:** Calls `dendropy.calculate.treemeasure.B1`. For each non-root internal node, it adds 1/m, where m is one more than the largest m among its children; tips have m = 0.
- **API:** `dnakit.phylogeny.calculate_b1_index(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path.
- **Output:** `result.b1` and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_b1_index

print(calculate_b1_index("((A,B),(C,(D,E)));")["result"])
```

- **Example output:**

```text
{'tip_count': 5, 'b1': 2.5}
```

## 13) `dnakit.phylogeny.calculate_pybus_harvey_gamma()`

- **Purpose:** Compute the Pybus and Harvey gamma statistic, which tests for constant birth and death rates.
- **Method:** Calls `dendropy.calculate.treemeasure.pybus_harvey_gamma`. The tree must be ultrametric within `precision` (root-to-tip distances equal); otherwise `ConfigurationError` is raised.
- **API:** `dnakit.phylogeny.calculate_pybus_harvey_gamma(tree, *, precision=1e-05, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path. `precision` is in (0, 1].
- **Output:** `result.gamma` and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_pybus_harvey_gamma

print(calculate_pybus_harvey_gamma("((A:1,B:1):1,(C:1,D:1):1);")["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'gamma': -0.816496580928}
```

## 14) `dnakit.phylogeny.calculate_treeness()`

- **Purpose:** Compute treeness: the fraction of total branch length that lies on internal branches.
- **Method:** Calls `dendropy.calculate.treemeasure.treeness`. Every branch must have a length.
- **API:** `dnakit.phylogeny.calculate_treeness(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path.
- **Output:** `result.treeness` and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_treeness

print(calculate_treeness("((A:1,B:1):1,(C:2,D:2):0.5);")["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'treeness': 0.2}
```

## 15) `dnakit.phylogeny.calculate_parsimony_score()`

- **Purpose:** Compute the Fitch parsimony score of DNA characters on a tree.
- **Method:** Maps the sequences onto the tree's tip labels and calls `dendropy.calculate.treescore.parsimony_score`. Sequence IDs must match the tree tips exactly.
- **API:** `dnakit.phylogeny.calculate_parsimony_score(tree, sequences, *, gaps_as_missing=True, max_sequences=1000, max_length=100000, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or one-tree file path, and a sequence dictionary or `DNA` batch as above. `gaps_as_missing` decides whether gaps count as missing data.
- **Output:** `result.parsimony_score` is an integer; `result.tip_count` and `result.site_count` give the tip and site counts.
- **Example code:**

```python
from dnakit.phylogeny import calculate_parsimony_score

chars = {"A": "AAAC", "B": "AAGC", "C": "TAGC", "D": "TCGC"}
print(calculate_parsimony_score("((A,B),(C,D));", chars)["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'site_count': 4, 'parsimony_score': 3, 'gaps_as_missing': True}
```

## 16) `dnakit.phylogeny.calculate_reconciliation_discordance()`

- **Purpose:** Count the gene duplications implied by reconciling a gene tree on a species tree.
- **Method:** Both trees are read as rooted (`rooting="force-rooted"`) and must have the same tip set. Calls `dendropy.model.reconcile.reconciliation_discordance`.
- **API:** `dnakit.phylogeny.calculate_reconciliation_discordance(gene_tree, species_tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A gene tree and a species tree, each a Newick string or one-tree file path.
- **Output:** `result.duplications` is the count; `result.tip_count` is the number of tips.
- **Example code:**

```python
from dnakit.phylogeny import calculate_reconciliation_discordance

print(calculate_reconciliation_discordance("((A,C),(B,(D,E)));", "((A,B),(C,(D,E)));")["result"])
```

- **Example output:**

```text
{'tip_count': 5, 'duplications': 2}
```

## 17) `dnakit.phylogeny.calculate_monophyletic_partition_discordance()`

- **Purpose:** Compute the deep coalescences needed if each tip group were forced to be monophyletic.
- **Method:** `groups` must partition all tips of the tree, each tip exactly once. Calls `dendropy.model.reconcile.monophyletic_partition_discordance`.
- **API:** `dnakit.phylogeny.calculate_monophyletic_partition_discordance(tree, groups, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or one-tree file path, and a list of lists of tip labels.
- **Output:** `result.deep_coalescences`, `result.group_count`, and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_monophyletic_partition_discordance

groups = [["A", "B"], ["C"], ["D", "E"]]
print(calculate_monophyletic_partition_discordance("((A,C),(B,(D,E)));", groups)["result"])
```

- **Example output:**

```text
{'tip_count': 5, 'deep_coalescences': 1, 'group_count': 3}
```

## 18) `dnakit.phylogeny.root_tree_on_midpoint()`

- **Purpose:** Root a tree at the midpoint of its longest tip-to-tip path and return Newick.
- **Method:** Calls `toytree.mod.root_on_midpoint`, which finds the longest path and places the root at its midpoint.
- **API:** `dnakit.phylogeny.root_tree_on_midpoint(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path. `max_tips` and `max_file_bytes` are explicit limits.
- **Output:** `result.newick` is the rooted Newick; `result.tips` is the sorted tip labels.
- **Example code:**

```python
from dnakit.phylogeny import root_tree_on_midpoint

print(root_tree_on_midpoint("(A:1,B:1,(C:2,D:2):1.5);")["result"]["newick"])
```

- **Example output:**

```text
((C:2,D:2):0.25,(A:1,B:1):1.25);
```

## 19) `dnakit.phylogeny.root_tree_on_balanced_midpoint()`

- **Purpose:** Root a tree at the balanced midpoint so both sides have similar depth.
- **Method:** Calls `toytree.mod.root_on_balanced_midpoint`; `tolerance` controls convergence.
- **API:** `dnakit.phylogeny.root_tree_on_balanced_midpoint(tree, *, tolerance=1e-06, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path. `tolerance` is a positive number.
- **Output:** `result.newick`, `result.tolerance`, and tip information.
- **Example code:**

```python
from dnakit.phylogeny import root_tree_on_balanced_midpoint

print(root_tree_on_balanced_midpoint("(A:1,B:1,(C:2,D:2):1.5);", tolerance=1e-9)["result"]["newick"])
```

- **Example output:**

```text
((C:2,D:2):0.249999999069,(A:1,B:1):1.25000000093);
```

## 20) `dnakit.phylogeny.root_tree_on_minimal_ancestor_deviation()`

- **Purpose:** Root a tree on the edge that minimizes ancestor deviation (MAD rooting).
- **Method:** Calls `toytree.mod.root_on_minimal_ancestor_deviation` and roots on the edge that gives the smallest deviation.
- **API:** `dnakit.phylogeny.root_tree_on_minimal_ancestor_deviation(tree, *, min_dist=1e-12, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path. `min_dist` is a positive number.
- **Output:** `result.newick` and `result.min_dist`.
- **Example code:**

```python
from dnakit.phylogeny import root_tree_on_minimal_ancestor_deviation

print(root_tree_on_minimal_ancestor_deviation("(A:1,B:1,(C:2,D:2):1.5);")["result"]["newick"])
```

- **Example output:**

```text
((C:2,D:2):0.25,(A:1,B:1):1.25);
```

## 21) `dnakit.phylogeny.unroot_tree()`

- **Purpose:** Remove the root split and return the unrooted tree as Newick.
- **Method:** Calls `toytree.mod.unroot`, which merges the two root edges into one.
- **API:** `dnakit.phylogeny.unroot_tree(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path.
- **Output:** `result.newick` is the unrooted Newick.
- **Example code:**

```python
from dnakit.phylogeny import unroot_tree

print(unroot_tree("((A:1,B:1):1,(C:2,D:2):0.5);")["result"]["newick"])
```

- **Example output:**

```text
(A:1,B:1,(C:2,D:2):1.5);
```

## 22) `dnakit.phylogeny.ladderize_tree()`

- **Purpose:** Order each node's children by clade size and return Newick.
- **Method:** Calls `toytree.mod.ladderize`. With `direction=True` the order is reversed.
- **API:** `dnakit.phylogeny.ladderize_tree(tree, *, direction=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path; `direction` is a boolean.
- **Output:** `result.newick` and `result.direction`.
- **Example code:**

```python
from dnakit.phylogeny import ladderize_tree

print(ladderize_tree("((A:1,B:1):1,(C:2,(D:2,E:2):1):0.5);", direction=True)["result"]["newick"])
```

- **Example output:**

```text
(((D:2,E:2):1,C:2):0.5,(A:1,B:1):1);
```

## 23) `dnakit.phylogeny.extract_subtree()`

- **Purpose:** Extract the clade spanned by the given tips and return it as Newick.
- **Method:** Calls `toytree.mod.extract_subtree`. The clade root keeps its stem length. Every tip must exist in the tree.
- **API:** `dnakit.phylogeny.extract_subtree(tree, tips, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path; `tips` is at least two distinct tip labels.
- **Output:** `result.newick`, `result.tips`, and `result.requested_tips`.
- **Example code:**

```python
from dnakit.phylogeny import extract_subtree

print(extract_subtree("((A:1,B:1):1,(C:2,D:2):0.5);", ["C", "D"])["result"]["newick"])
```

- **Example output:**

```text
(C:2,D:2):0.5;
```

## 24) `dnakit.phylogeny.collapse_weak_nodes()`

- **Purpose:** Collapse internal nodes whose support is below a threshold and return Newick.
- **Method:** Calls `toytree.mod.collapse_nodes(min_support=...)`. Support values come from the internal-node labels in the Newick.
- **API:** `dnakit.phylogeny.collapse_weak_nodes(tree, *, min_support, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path; `min_support` is a non-negative number.
- **Output:** `result.newick` and `result.min_support`.
- **Example code:**

```python
from dnakit.phylogeny import collapse_weak_nodes

print(collapse_weak_nodes("((A:1,B:1)90:1,(C:2,D:2)30:0.5);", min_support=50)["result"]["newick"])
```

- **Example output:**

```text
((A:1,B:1)90:1,C:2.5,D:2.5);
```

## 25) `dnakit.phylogeny.remove_unary_nodes()`

- **Purpose:** Remove internal nodes that have a single child and merge their branch lengths.
- **Method:** Calls `toytree.mod.remove_unary_nodes`.
- **API:** `dnakit.phylogeny.remove_unary_nodes(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path.
- **Output:** `result.newick` is the Newick without single-child nodes.
- **Example code:**

```python
from dnakit.phylogeny import remove_unary_nodes

print(remove_unary_nodes("((A:1):1,(B:2,C:2):1);")["result"]["newick"])
```

- **Example output:**

```text
((B:2,C:2):1,A:2);
```

## 26) `dnakit.phylogeny.scale_tree_to_root_height()`

- **Purpose:** Rescale all branch lengths so the root-to-tip height equals a given value.
- **Method:** Calls `toytree.mod.edges_scale_to_root_height`. With `include_stem=True`, the root stem is also scaled.
- **API:** `dnakit.phylogeny.scale_tree_to_root_height(tree, *, height=1.0, include_stem=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path; `height` is positive; `include_stem` is a boolean.
- **Output:** `result.newick`, `result.height`, and `result.include_stem`.
- **Example code:**

```python
from dnakit.phylogeny import scale_tree_to_root_height

print(scale_tree_to_root_height("((A:1,B:1):1,(C:2,D:2):0.5);", height=2.0)["result"]["newick"])
```

- **Example output:**

```text
((A:0.8,B:0.8):0.8,(C:1.6,D:1.6):0.4);
```

## 27) `dnakit.phylogeny.calculate_tree_vcv_matrix()`

- **Purpose:** Compute the tip variance-covariance (VCV) matrix of a tree.
- **Method:** Calls `toytree.pcm.get_vcv_matrix_from_tree`. Diagonal entries are root-to-tip distances; off-diagonal entries are the shared path length from the root. Rows and columns are sorted by tip label.
- **API:** `dnakit.phylogeny.calculate_tree_vcv_matrix(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** The example file `tree_t1.nwk` or a Newick string.
- **Output:** `result.tips` and `result.matrix` (nested list).
- **Example code:**

```python
from dnakit.phylogeny import calculate_tree_vcv_matrix

print(calculate_tree_vcv_matrix("tree_t1.nwk")["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'tips': ['A', 'B', 'C', 'D'], 'matrix': [[2.0, 1.0, 0.0, 0.0], [1.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.5, 0.5], [0.0, 0.0, 0.5, 2.5]]}
```

## 28) `dnakit.phylogeny.calculate_tree_correlation_matrix()`

- **Purpose:** Compute the tip correlation matrix of a tree.
- **Method:** Calls `toytree.pcm.get_corr_matrix_from_tree`; each VCV entry is divided by the square root of the product of the two tip variances.
- **API:** `dnakit.phylogeny.calculate_tree_correlation_matrix(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path.
- **Output:** `result.tips` and `result.matrix`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_tree_correlation_matrix

print(calculate_tree_correlation_matrix("((A:1,B:1):1,(C:2,D:2):0.5);")["result"]["matrix"])
```

- **Example output:**

```text
[[1.0, 0.5, 0.0, 0.0], [0.5, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.2], [0.0, 0.0, 0.2, 1.0]]
```

## 29) `dnakit.phylogeny.calculate_information_rf_distance()`

- **Purpose:** Compute the information-corrected Robinson-Foulds distance (rfi).
- **Method:** Calls `toytree.distance.get_treedist_rfi`.
- **API:** `dnakit.phylogeny.calculate_information_rf_distance(tree_a, tree_b, *, normalize=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** Two trees (Newick strings or file paths) with identical tip labels; `normalize` defaults to `False`.
- **Output:** `result.distance`, `result.normalize`, and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_information_rf_distance

print(calculate_information_rf_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'distance': 3.169925001442, 'normalize': False}
```

## 30) `dnakit.phylogeny.calculate_matching_split_distance()`

- **Purpose:** Compute the generalized Robinson-Foulds matching split distance (rfg_ms).
- **Method:** Calls `toytree.distance.get_treedist_rfg_ms`. toytree provides no normalization for this metric.
- **API:** `dnakit.phylogeny.calculate_matching_split_distance(tree_a, tree_b, *, normalize=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** Two trees (Newick strings or file paths) with identical tip labels; `normalize` defaults to `False`.
- **Output:** `result.distance`, `result.normalize`, and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_matching_split_distance

print(calculate_matching_split_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'distance': 2.0, 'normalize': False}
```

## 31) `dnakit.phylogeny.calculate_matching_split_information_distance()`

- **Purpose:** Compute the matching split information distance (rfg_msi).
- **Method:** Calls `toytree.distance.get_treedist_rfg_msi`.
- **API:** `dnakit.phylogeny.calculate_matching_split_information_distance(tree_a, tree_b, *, normalize=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** Two trees (Newick strings or file paths) with identical tip labels; `normalize` defaults to `False`.
- **Output:** `result.distance`, `result.normalize`, and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_matching_split_information_distance

print(calculate_matching_split_information_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'distance': 2.0, 'normalize': False}
```

## 32) `dnakit.phylogeny.calculate_shared_phylogenetic_information_distance()`

- **Purpose:** Compute the generalized Robinson-Foulds shared phylogenetic information distance (rfg_spi).
- **Method:** Calls `toytree.distance.get_treedist_rfg_spi`. When both trees have identical split sets, toytree's similarity is undefined and this wrapper raises `ConfigurationError`.
- **API:** `dnakit.phylogeny.calculate_shared_phylogenetic_information_distance(tree_a, tree_b, *, normalize=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** Two trees (Newick strings or file paths) with identical tip labels; `normalize` defaults to `False`.
- **Output:** `result.distance`, `result.normalize`, and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_shared_phylogenetic_information_distance

print(calculate_shared_phylogenetic_information_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'distance': 3.169925001442, 'normalize': False}
```

## 33) `dnakit.phylogeny.calculate_mutual_clustering_information_distance()`

- **Purpose:** Compute the generalized Robinson-Foulds mutual clustering information distance (rfg_mci).
- **Method:** Calls `toytree.distance.get_treedist_rfg_mci`.
- **API:** `dnakit.phylogeny.calculate_mutual_clustering_information_distance(tree_a, tree_b, *, normalize=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** Two trees (Newick strings or file paths) with identical tip labels; `normalize` defaults to `False`.
- **Output:** `result.distance`, `result.normalize`, and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_mutual_clustering_information_distance

print(calculate_mutual_clustering_information_distance("((A:1,B:1):1,(C:2,D:2):0.5);", "((A:1,C:1):1,(B:2,D:2):0.5);")["result"])
```

- **Example output:**

```text
{'tip_count': 4, 'distance': 2.0, 'normalize': False}
```

## 34) `dnakit.phylogeny.list_tree_bipartitions()`

- **Purpose:** List the splits of a tree; each split is written as two sorted tip lists.
- **Method:** Calls `toytree.enum.iter_bipartitions`. Within a split, the side whose first tip sorts first comes first; splits are sorted. `include_singletons=True` includes splits with a single tip.
- **API:** `dnakit.phylogeny.list_tree_bipartitions(tree, *, include_singletons=False, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path; `include_singletons` is a boolean.
- **Output:** `result.splits` (list of splits), `result.split_count`, and `result.tip_count`.
- **Example code:**

```python
from dnakit.phylogeny import list_tree_bipartitions

print(list_tree_bipartitions("((A:1,B:1):1,(C:2,D:2):0.5);", include_singletons=True)["result"]["splits"])
```

- **Example output:**

```text
[[['A'], ['B', 'C', 'D']], [['A', 'B'], ['C', 'D']], [['A', 'B', 'C'], ['D']], [['A', 'B', 'D'], ['C']], [['A', 'C', 'D'], ['B']]]
```

## 35) `dnakit.phylogeny.calculate_relative_evolutionary_divergence()`

- **Purpose:** Compute the relative evolutionary divergence (RED) of every node.
- **Method:** Calls `toytree.pcm.get_relative_evolutionary_divergence` and sorts the output by node index.
- **API:** `dnakit.phylogeny.calculate_relative_evolutionary_divergence(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path.
- **Output:** `result.nodes` is a list; each item has `node_index`, `name` (null for internal nodes), and `red`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_relative_evolutionary_divergence

print(calculate_relative_evolutionary_divergence("((A:1,B:1):1,(C:2,D:2):0.5);")["result"]["nodes"])
```

- **Example output:**

```text
[{'node_index': 0, 'name': 'A', 'red': 1.0}, {'node_index': 1, 'name': 'B', 'red': 1.0}, {'node_index': 2, 'name': 'C', 'red': 1.0}, {'node_index': 3, 'name': 'D', 'red': 1.0}, {'node_index': 5, 'name': None, 'red': 0.2}, {'node_index': 6, 'name': None, 'red': 0.0}]
```

## 36) `dnakit.phylogeny.calculate_equal_splits()`

- **Purpose:** Compute the equal splits (ES) metric for every tip.
- **Method:** Calls `toytree.pcm.get_equal_splits`, keeps tips only (the root is dropped), and sorts by label.
- **API:** `dnakit.phylogeny.calculate_equal_splits(tree, *, max_tips=10000, max_file_bytes=10000000, verbose=True)`.
- **Input:** A Newick string or a one-tree file path.
- **Output:** `result.tips` is a list; each item has `tip` and `equal_splits`.
- **Example code:**

```python
from dnakit.phylogeny import calculate_equal_splits

print(calculate_equal_splits("((A:1,B:1):1,(C:2,D:2):0.5);")["result"]["tips"])
```

- **Example output:**

```text
[{'tip': 'A', 'equal_splits': 0.5}, {'tip': 'B', 'equal_splits': 0.5}, {'tip': 'C', 'equal_splits': 0.25}, {'tip': 'D', 'equal_splits': 0.25}]
```

## Current scope

- pyvolve adds no wrappers.
- `calculate_shared_phylogenetic_information_distance()` raises `ConfigurationError` when two trees have identical split sets, because toytree's similarity is undefined in that case.
- This document covers only the 36 wrapper functions above. It does not cover upstream aliases, functions that duplicate existing dnakit functions, or functions beyond the 20-per-package cap.
