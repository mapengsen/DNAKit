# Phylogenetics and Molecular Evolution

## 1) Infer a maximum-likelihood tree with IQ-TREE 2

- **Purpose:** Pass an existing multiple-sequence alignment to caller-installed IQ-TREE 2 and retain the maximum-likelihood tree and report files.
- **Method:** Run the fixed `iqtree2 -s ... -m ... -nt ... -pre ...` command; optionally add `-B` ultrafast bootstrap. DNAKit does not copy or install IQ-TREE 2 source.
- **API:** `dnakit.phylogeny.infer_ml_tree(alignment, output_directory, model="MFP", bootstrap=0, threads=1, ...)`
- **Input:** An alignment format supported by IQ-TREE 2, a new output directory, and caller-installed `iqtree2` on PATH.
- **Output:** A JSON-serializable stage envelope containing `.treefile`, `.iqtree`, `.log` and other artifact paths, plus model, bootstrap and provenance.
- **Citation and license:** Cite Minh et al. (2020), DOI `10.1093/molbev/msaa015`; the official repository records GPL-2.0. If using ModelFinder, UFBoot or another specialized feature, add the corresponding IQ-TREE 2 citations.

```python
from dnakit.phylogeny import infer_ml_tree

result = infer_ml_tree(
    "alignment.fasta",
    "iqtree-result",
    model="MFP",
    bootstrap=1000,
)
print(result["result"]["treefile"])
```

## 2) Infer an approximate maximum-likelihood tree with FastTree

- **Purpose:** Pass an existing alignment to caller-installed FastTree and capture its approximate maximum-likelihood Newick tree.
- **Method:** Run the fixed `FastTree [-nt] [-gtr] [-gamma] alignment` command; stdout becomes `tree.newick` and stderr is retained separately. DNAKit does not copy or install FastTree.
- **API:** `dnakit.phylogeny.infer_fast_tree(alignment, output_directory, nucleotide=True, gtr=False, gamma=False, ...)`
- **Input and output:** A local alignment and a new output directory; the JSON-serializable envelope contains `tree.newick`, stderr, fixed options and provenance.
- **Citation and license:** Cite Price, Dehal & Arkin (2010), [10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490); the FastTree 2 repository records GPL-3.0.

## 1) Build a tree from a distance matrix {#bio-6}

- **Purpose:** Build a tree with branch lengths using NJ or UPGMA.
- **Method:** NJ uses neighbor joining; UPGMA merges clusters using size-weighted distances.
- **API:** `dnakit.bioinformatics.phylogeny.build_tree(distances, ids, method='nj')`
- **Input:** A symmetric distance matrix with a zero diagonal and matching IDs; UPGMA assumes a molecular clock.
- **Output:** Returns a `TreeNode`; the example shows `tips`.
- **Example code:**

```python
from dnakit.bioinformatics.phylogeny import build_tree

tree = build_tree([[0, 2, 4], [2, 0, 4], [4, 4, 0]], ["a", "b", "c"], method="upgma")
print(sorted(tree.tips))
```

- **Example output:**

```text
['a', 'b', 'c']
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#alignment).

## 2) Compare, combine and prune trees {#bio-8}

- **Purpose:** Compute RF distance, form a consensus tree or retain selected tips.
- **Method:** RF compares split sets; consensus retains supported branches; pruning collapses single-child nodes.
- **API:** `dnakit.bioinformatics.phylogeny.robinson_foulds(left, right, rooted=False)`; `dnakit.bioinformatics.phylogeny.consensus_tree(trees, threshold=0.5)`; `dnakit.bioinformatics.phylogeny.prune_tree(tree, retain)`
- **Input:** RF and consensus require matching tip sets and consistent rooting; `threshold=1` requests strict consensus.
- **Output:** `robinson_foulds` returns a distance; `consensus_tree` and `prune_tree` return trees.
- **Example code:**

```python
from dnakit.bioinformatics.phylogeny import build_tree, robinson_foulds, consensus_tree, prune_tree

tree = build_tree([[0, 2, 4], [2, 0, 4], [4, 4, 0]], ["a", "b", "c"])
print(robinson_foulds(tree, tree))
print(sorted(consensus_tree([tree, tree]).tips))
print(sorted(prune_tree(tree, ["a", "c"]).tips))
```

- **Example output:**

```text
0.0
['a', 'b', 'c']
['a', 'c']
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#alignment).

## 3) dN/dS selection analysis {#12-dnds-selection-analysis}

- **Purpose:** Estimate Nei–Gojobori dN, dS, and dN/dS with Jukes–Cantor correction.
- **Method:** Count synonymous/nonsynonymous sites and differences with Nei–Gojobori and apply Jukes–Cantor correction.
- **API:** `dnakit.comparative.calculate_dn_ds(sequence_a[required], sequence_b[required])`
- **Input:** Equal-length, codon-aligned, linear `DNASequence`/`DNARecord`/`DNA` values containing canonical A/C/G/T, with lengths divisible by 3.
- **Output:** Returns `ProviderResult` with query content in `data`, plus `provider` and provenance.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.comparative import calculate_dn_ds

result = calculate_dn_ds(DNASequence("ATGGCTGAA"), DNASequence("ATGGCCGAG"))
print(result.data["dN"], result.data["dS"], result.data["dN_dS"])
```

- **Example output:** `dN`, `dS`, `dN_dS`, site/difference counts, and interpretation; the ratio is `None` when `dS=0`.

[Optional backend dependencies and return conventions](../../faq/tutorial_integrations.md#optional-bioinformatics)

<!-- BEGIN GENERATED API: selection -->

<span id="selection"></span>

## 4) Calculate dN/dS between two coding sequences {#dnakit.databases.sequence_dn_ds}

- **Purpose:** Calculate dN/dS between two coding sequences.
- **Method:** Apply the Nei–Gojobori estimator with Jukes–Cantor correction to in-frame codon-aligned sequences.
- **API：** `dnakit.databases.sequence_dn_ds(*, seq1=None, seq2=None, fasta1_path=None, fasta2_path=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `dN_dS`, `dN`, `dS`, `pN` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sequence_dn_ds(
    seq1='ATGGCCAAATTTGAAGATCTGAGCACCGTGCGTGGC',
    seq2='ATGGCCAAGTTCGAAAATCTGAGCACCGTACGTGGC',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "dN_dS": "<number>",
  "dN": "<number>",
  "dS": "<number>",
  "pN": "<number>",
  "pS": "<number>",
  "N_sites": "<number>",
  "S_sites": "<number>",
  "Nd": "<number>"
}
```

[Full parameters and usage notes](../../faq/api/30_phylogeny.md#dnakit.databases.sequence_dn_ds-parameters)

<!-- END GENERATED API: selection -->
