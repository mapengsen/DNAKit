# Optional tskit and dadi population-genetics functions (group 80)

This group wraps public tskit and dadi functions that work on DNA or gene-level data and run offline from local inputs. The tskit part is in `dnakit.population`: sequence and interval export, text-table reading and writing, PLINK `.fam` reading, topology counts, polytomy resolution, and Newick output. The dadi part is in `dnakit.population`: SNP and VCF data reading, spectrum reading, writing and pooling, bootstraps, expected spectra for demographic models, model fitting, and Fisher-information uncertainty.

The functions are exported from `dnakit.population`, and the examples import them from that package. All dependencies are optional and are imported only when the matching function is called; a missing package raises an error with the install command. Output files are never overwritten, and models are chosen by name from fixed lists, so no arbitrary code runs.

## Example data

Items 1 to 12 are the tskit examples, and items 13 to 30 are the dadi examples. The files below are written into the working directory before the examples run. `tables.json` has 4 samples, 1 tree, and 3 sites; it only demonstrates the interface, and its statistics do not describe a real population.

**Example file:** `tables.json`

```json
{
 "sequence_length": 100,
 "nodes": [
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 1
  },
  {
   "time": 2
  },
  {
   "time": 3
  }
 ],
 "edges": [
  {
   "left": 0,
   "right": 100,
   "parent": 4,
   "child": 0
  },
  {
   "left": 0,
   "right": 100,
   "parent": 4,
   "child": 1
  },
  {
   "left": 0,
   "right": 100,
   "parent": 5,
   "child": 2
  },
  {
   "left": 0,
   "right": 100,
   "parent": 5,
   "child": 3
  },
  {
   "left": 0,
   "right": 100,
   "parent": 6,
   "child": 4
  },
  {
   "left": 0,
   "right": 100,
   "parent": 6,
   "child": 5
  }
 ],
 "sites": [
  {
   "position": 10,
   "ancestral_state": "A"
  },
  {
   "position": 20,
   "ancestral_state": "A"
  },
  {
   "position": 30,
   "ancestral_state": "A"
  }
 ],
 "mutations": [
  {
   "site": 0,
   "node": 0,
   "derived_state": "G"
  },
  {
   "site": 1,
   "node": 4,
   "derived_state": "G"
  },
  {
   "site": 2,
   "node": 5,
   "derived_state": "G"
  }
 ]
}
```

**Example file:** `nodes.txt`

```text
id	is_sample	time	population	individual	metadata
0	1	0.000000	-1	-1	
1	1	0.000000	-1	-1	
2	1	0.000000	-1	-1	
3	1	0.000000	-1	-1	
4	0	1.000000	-1	-1	
5	0	2.000000	-1	-1	
6	0	3.000000	-1	-1	
```

**Example file:** `edges.txt`

```text
left	right	parent	child	metadata
0.000000	100.000000	4	0	
0.000000	100.000000	4	1	
0.000000	100.000000	5	2	
0.000000	100.000000	5	3	
0.000000	100.000000	6	4	
0.000000	100.000000	6	5	
```

**Example file:** `sites.txt`

```text
position	ancestral_state	metadata
10.000000	A	
20.000000	A	
30.000000	A	
```

**Example file:** `mutations.txt`

```text
site	node	time	derived_state	parent	metadata
0	0	unknown	G	-1	
1	4	unknown	G	-1	
2	5	unknown	G	-1	
```

**Example file:** `samples.fam`

```text
f1 i1 0 0 1
f1 i2 0 0 2
f2 i3 0 0 0
```

**Example file:** `snps.txt`

```text
Context OutgroupContext Allele1 pop1 pop2 Allele2 SNPid
ACA TCA A 2 1 G 0 1 snp1
GCT GCT C 1 2 T 1 0 snp2
TTG TAG T 2 2 G 0 0 snp3
ACG ACG A 1 1 C 1 1 snp4
```

**Example file:** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4
chr1	100	.	A	G	.	PASS	AA=A	GT	0|0	0|1	1|1	0|0
chr1	200	.	C	T	.	PASS	AA=T	GT	0|1	1|1	0|0	0|0
chr1	300	.	G	A	.	PASS	AA=G	GT	0|0	0|0	0|1	1|1
```

**Example file:** `popinfo.txt`

```text
s1 popA
s2 popA
s3 popB
s4 popB
```

**Example file:** `chunks.txt`

```text
Context OutgroupContext Allele1 pop1 pop2 Allele2 SNPid
ACG ACG A 2 1 G 0 1 chr1_100
ACG ACG A 1 0 G 1 2 chr1_200
ACG ACG A 0 2 G 2 0 chr1_1100
ACG ACG A 1 1 G 1 1 chr1_1200
ACG ACG A 2 0 G 0 2 chr1_2100
ACG ACG A 0 1 G 2 1 chr1_2200
```

**Example file:** `sample.exonic_variant_function`

```text
1	nonsynonymous SNV	GENE1	chr1 100
2	synonymous SNV	GENE2	chr1 200
```

**Example file:** `spectrum.sfs`

```text
# demo
2 2 
1.5 2 3 4
```

**Example file:** `model.yaml`

```text
description: split with migration
time_units: generations
demes:
  - name: ancestral
    epochs:
      - end_time: 2000
        start_size: 1000
  - name: A
    ancestors: [ancestral]
    epochs:
      - end_time: 0
        start_size: 1000
  - name: B
    ancestors: [ancestral]
    epochs:
      - end_time: 0
        start_size: 1000
migrations:
  - source: A
    dest: B
    start_time: 2000
    end_time: 0
    rate: 0.0001
```

## 1) `dnakit.population.export_sample_sequences()`

- **Purpose:** Write each sample's haplotype sequence to a FASTA or NEXUS file.
- **Method:** Call tskit's `write_fasta()` or `write_nexus()`. Sequences are derived from the sites and mutations; positions without a site are filled with `missing_data_character`. An existing output file raises an error.
- **API:** `dnakit.population.export_sample_sequences(tables, destination, *, file_format='fasta', wrap_width=60, missing_data_character='N', max_rows=100000)`
- **Input:** `tables` is the same as for `build_tree_sequence()`; `destination` is a new file path; `file_format` is `fasta` (default) or `nexus`; `wrap_width` is the FASTA line width, and 0 disables wrapping; `missing_data_character` is a single character.
- **Output:** A dictionary with `file_format`, `destination`, `bytes`, `num_samples`, `sequence_length`, `num_sites`, and `provenance`. With `file_format="nexus"` a NEXUS file is written instead.
- **Example code:**

```python
import json
from pathlib import Path
from dnakit.population import export_sample_sequences

Path("samples.fa").unlink(missing_ok=True)  # the function never overwrites an existing file
result = export_sample_sequences(json.load(open("tables.json")), "samples.fa", wrap_width=60)
lines = Path("samples.fa").read_text().splitlines()
print(result["num_samples"], lines[0], lines[1][10], lines[1][20], lines[1][30])
```

- **Example output:**

```text
4 >n0 G G A
```

## 2) `dnakit.population.export_ms_format()`

- **Purpose:** Write the genotypes of a tree sequence in Hudson's `ms` output format.
- **Method:** Call tskit's `write_ms()`. Each row is a sample and each column a segregating site, with 1 for the derived allele. `print_trees=True` also writes the Newick tree.
- **API:** `dnakit.population.export_ms_format(tables, destination, *, print_trees=False, precision=4, max_rows=100000)`
- **Input:** `tables` is the same as for `build_tree_sequence()`; `destination` is a new file path; `precision` is the number of digits in site positions.
- **Output:** A dictionary with `destination`, `bytes`, `num_samples`, `num_segregating_sites`, and `provenance`.
- **Example code:**

```python
import json
from pathlib import Path
from dnakit.population import export_ms_format

Path("samples.ms").unlink(missing_ok=True)
result = export_ms_format(json.load(open("tables.json")), "samples.ms", precision=4)
print(result["num_segregating_sites"], Path("samples.ms").read_text().splitlines()[6:10])
```

- **Example output:**

```text
3 ['110', '010', '001', '001']
```

## 3) `dnakit.population.complement_intervals()`

- **Purpose:** Return the half-open intervals in `[start, end)` that the given intervals do not cover.
- **Method:** Call tskit's `negate_intervals()`. Input intervals must be sorted, non-overlapping, and inside `[start, end)`. An empty list returns the whole range.
- **API:** `dnakit.population.complement_intervals(intervals, *, start=0.0, end)`
- **Input:** `intervals` is a list of `[left, right]` pairs; `start` defaults to 0; `end` is required.
- **Output:** A dictionary with `intervals` (the complement), `num_intervals`, `complement_length`, `covered_length`, and `provenance`.
- **Example code:**

```python
from dnakit.population import complement_intervals

result = complement_intervals([(10, 20), (30, 40)], end=100)
print(result["intervals"], result["complement_length"])
```

- **Example output:**

```text
[[0.0, 10.0], [20.0, 30.0], [40.0, 100.0]] 80.0
```

## 4) `dnakit.population.read_text_table()`

- **Purpose:** Read one tskit text table file and return its rows.
- **Method:** `table` selects tskit's parser: `parse_nodes`, `parse_edges`, `parse_sites`, `parse_mutations`, `parse_individuals`, `parse_populations`, or `parse_migrations`. Metadata columns are not returned. `strict=True` requires tab-delimited input.
- **API:** `dnakit.population.read_text_table(source, *, table='nodes', strict=True, max_rows=100000)`
- **Input:** `source` is the path of a text table file; `table` is the table type; `max_rows` limits the number of rows.
- **Output:** A dictionary with `table`, `num_rows`, `columns`, and `rows` (one dictionary per row).
- **Example code:**

```python
from dnakit.population import read_text_table

result = read_text_table("sites.txt", table="sites")
print(result["num_rows"], [row["position"] for row in result["rows"]])
```

- **Example output:**

```text
3 [10.0, 20.0, 30.0]
```

## 5) `dnakit.population.read_text_tree_sequence()`

- **Purpose:** Build a tree sequence from tskit text table files and summarize it.
- **Method:** Call `tskit.load_text()`. `nodes` and `edges` are required; `sites` and `mutations` are optional. With `sequence_length=0`, tskit takes the length from the edges.
- **API:** `dnakit.population.read_text_tree_sequence(nodes, edges, *, sites=None, mutations=None, sequence_length=0.0, strict=True)`
- **Input:** `nodes` and `edges` are required text files; `sites` and `mutations` are optional; `strict` is as in `read_text_table()`.
- **Output:** A dictionary with `sources`, the summary fields such as `num_samples`, `num_trees`, and `num_mutations`, and `provenance`.
- **Example code:**

```python
from dnakit.population import read_text_tree_sequence

result = read_text_tree_sequence("nodes.txt", "edges.txt", sites="sites.txt", mutations="mutations.txt")
print(result["num_samples"], result["num_trees"], result["num_mutations"])
```

- **Example output:**

```text
4 1 3
```

## 6) `dnakit.population.write_text_tree_sequence()`

- **Purpose:** Write nodes, edges, sites, and mutations as tskit text tables.
- **Method:** Call `tskit.dump_text()` to write the four tables into a new directory. An existing directory raises an error. The output can be read back with `read_text_tree_sequence()`.
- **API:** `dnakit.population.write_text_tree_sequence(tables, destination_dir, *, precision=6, max_rows=100000)`
- **Input:** `tables` is the same as for `build_tree_sequence()`; `destination_dir` must not exist; `precision` is the number of digits in times.
- **Output:** A dictionary with `destination_dir`, `files`, `num_nodes`, `num_edges`, `num_sites`, `num_mutations`, and `provenance`.
- **Example code:**

```python
import json
import shutil
from dnakit.population import write_text_tree_sequence

shutil.rmtree("text_out", ignore_errors=True)  # the directory must not exist yet
result = write_text_tree_sequence(json.load(open("tables.json")), "text_out")
print(result["files"], result["num_mutations"])
```

- **Example output:**

```text
['nodes.txt', 'edges.txt', 'sites.txt', 'mutations.txt'] 3
```

## 7) `dnakit.population.count_sample_topologies()`

- **Purpose:** Count the embedded topologies of named sample sets in every tree.
- **Method:** Call tskit's `treeseq_count_topologies()` and return the counts for each combination of sample sets, sorted by topology rank (`shape`, `label`). Labels join set names with `+`. Each set needs at least two samples.
- **API:** `dnakit.population.count_sample_topologies(tables, sample_sets, *, max_rows=100000)`
- **Input:** `tables` is the same as for `build_tree_sequence()`; `sample_sets` maps names to lists of sample node IDs.
- **Output:** A dictionary with `sample_sets`, `num_trees`, and `trees` (the counts for each tree).
- **Example code:**

```python
import json
from dnakit.population import count_sample_topologies

result = count_sample_topologies(json.load(open("tables.json")), {"A": [0, 1], "B": [2, 3]})
print(result["num_trees"], sorted(result["trees"][0]["counts"]))
```

- **Example output:**

```text
1 ['A', 'A+B', 'B']
```

## 8) `dnakit.population.resolve_polytomies()`

- **Purpose:** Resolve the multifurcations of one tree into bifurcations and return Newick.
- **Method:** Call tskit's `split_polytomies()`. `random_seed` makes the random resolution reproducible. A tree that is already binary is unchanged, with `added_nodes` equal to 0.
- **API:** `dnakit.population.resolve_polytomies(tables, *, tree_index=0, random_seed=None, max_rows=100000)`
- **Input:** `tables`; `tree_index` selects the tree; `random_seed` is optional.
- **Output:** A dictionary with `newick`, `added_nodes`, `added_edges`, and `provenance`.
- **Example code:**

```python
from dnakit.population import resolve_polytomies

star = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": 1}],
    "edges": [{"left": 0, "right": 100, "parent": 4, "child": child} for child in range(4)],
    "sites": [],
    "mutations": [],
}
result = resolve_polytomies(star, random_seed=7)
print(result["added_nodes"], result["newick"].count("("))
```

- **Example output:**

```text
2 3
```

## 9) `dnakit.population.tree_to_newick()`

- **Purpose:** Return one tree of a tree sequence as a Newick string.
- **Method:** Call tskit's `build_newick()`, with sample nodes named by `taxon_labels` (default `n<node id>`). Labels must be unique and must not contain Newick-reserved characters.
- **API:** `dnakit.population.tree_to_newick(tables, *, tree_index=0, precision=6, include_branch_lengths=True, taxon_labels=None, max_rows=100000)`
- **Input:** `tables`; `tree_index`; `precision` is the number of digits in branch lengths; `include_branch_lengths` controls whether lengths are written.
- **Output:** A dictionary with `tree_index`, `num_samples`, `newick`, and `provenance`.
- **Example code:**

```python
import json
from dnakit.population import tree_to_newick

result = tree_to_newick(json.load(open("tables.json")), taxon_labels=["a", "b", "c", "d"], precision=3)
print(result["newick"])
```

- **Example output:**

```text
((a:1.000,b:1.000):2.000,(c:2.000,d:2.000):1.000);
```

## 10) `dnakit.population.random_tree_newick()`

- **Purpose:** Generate a random or regular tree on a given number of leaves and return Newick.
- **Method:** Call tskit's combinatorics generators. `random_binary` samples uniformly over leaf-labelled binary trees, and `seed` makes it reproducible; `balanced`, `star`, and `comb` are fixed shapes. Leaves are named `n0` to `n<num_leaves - 1>`.
- **API:** `dnakit.population.random_tree_newick(num_leaves, *, shape='random_binary', seed=None, arity=2, branch_length=1.0)`
- **Input:** `num_leaves` is 2 to 1000; `shape` is one of the four shapes; `seed` applies only to `random_binary`; `arity` applies only to `balanced`.
- **Output:** A dictionary with `shape`, `num_leaves`, `seed`, `newick`, and `provenance`.
- **Example code:**

```python
from dnakit.population import random_tree_newick

result = random_tree_newick(5, seed=3)
print(result["newick"])
```

- **Example output:**

```text
(n0:4.000000,(n3:3.000000,(n4:2.000000,(n1:1.000000,n2:1.000000):1.000000):1.000000):1.000000);
```

## 11) `dnakit.population.map_allele_indices()`

- **Purpose:** Map each allele in one list to its index in another list.
- **Method:** Call tskit's `allele_remap()`. Each allele takes its first position in `alleles_to`. Alleles missing from `alleles_to` get new indices after its end and are listed in `unmapped_alleles`. `alleles_to` must not repeat an allele.
- **API:** `dnakit.population.map_allele_indices(alleles_from, alleles_to)`
- **Input:** `alleles_from` and `alleles_to` are lists of allele strings.
- **Output:** A dictionary with `indices`, `num_alleles_to`, `unmapped_alleles`, and `provenance`.
- **Example code:**

```python
from dnakit.population import map_allele_indices

result = map_allele_indices(["A", "G", "T", "X"], ["A", "G"])
print(result["indices"], result["unmapped_alleles"])
```

- **Example output:**

```text
[0, 1, 2, 3] ['T', 'X']
```

## 12) `dnakit.population.read_plink_fam()`

- **Purpose:** Read a PLINK `.fam` file into a list of individuals.
- **Method:** Call tskit's `parse_fam()`. The file has five columns: family ID, individual ID, father, mother, and sex (1 male, 2 female, other codes unknown). Individuals are returned in file order.
- **API:** `dnakit.population.read_plink_fam(path)`
- **Input:** `path` is the path of a `.fam` file.
- **Output:** A dictionary with `num_individuals`, `num_males`, `num_females`, `num_unknown_sex`, `individuals`, and `provenance`.
- **Example code:**

```python
from dnakit.population import read_plink_fam

result = read_plink_fam("samples.fam")
print(result["num_individuals"], result["num_males"], result["num_females"], result["num_unknown_sex"])
```

- **Example output:**

```text
3 1 1 1
```

## 13) `dnakit.population.read_snp_data()`

- **Purpose:** Read a dadi SNP file and summarize it.
- **Method:** Call dadi's `Misc.make_data_dict()`. The header is `Context OutgroupContext Allele1 <pops...> Allele2 <id>`, and each row is one SNP. A SNP is polarized when its outgroup allele is one of the two segregating alleles.
- **API:** `dnakit.population.read_snp_data(path)`
- **Input:** `path` is a dadi SNP file (`.gz` and `.zip` are accepted).
- **Output:** A dictionary with `populations`, `num_snps`, `num_biallelic`, `num_polarized`, and `provenance`.
- **Example code:**

```python
from dnakit.population import read_snp_data

result = read_snp_data("snps.txt")
print(result["num_snps"], result["num_polarized"], result["populations"])
```

- **Example output:**

```text
4 2 ['pop1', 'pop2']
```

## 14) `dnakit.population.snp_configuration_counts()`

- **Purpose:** Count SNP configurations for the given populations.
- **Method:** Call dadi's `Misc.count_data_dict()`. Each row gives the successful and derived calls per population, whether the SNP is polarized, and the number of SNPs with that configuration. Rows are sorted by configuration. Non-biallelic sites are skipped.
- **API:** `dnakit.population.snp_configuration_counts(path, *, pop_ids)`
- **Input:** `path` is an SNP file; `pop_ids` lists population IDs, all of which must be in the file.
- **Output:** A dictionary with `pop_ids`, `num_configurations`, `num_snps_counted`, `configurations`, and `provenance`.
- **Example code:**

```python
from dnakit.population import snp_configuration_counts

result = snp_configuration_counts("snps.txt", pop_ids=["pop1", "pop2"])
print(result["num_configurations"], result["num_snps_counted"], result["configurations"][0])
```

- **Example output:**

```text
4 4 {'successful_calls': [2, 2], 'derived_calls': [0, 0], 'polarized': False, 'count': 1}
```

## 15) `dnakit.population.read_vcf_data()`

- **Purpose:** Read a VCF file and a population table into a dadi data dictionary, and summarize it.
- **Method:** `parser="builtin"` uses dadi's own VCF reader; `parser="cyvcf2"` uses cyvcf2, which must be installed. Both give the same dictionary for the same file. `filter=True` skips sites that did not pass filtering.
- **API:** `dnakit.population.read_vcf_data(vcf_path, popinfo_path, *, parser='builtin', filter=True)`
- **Input:** `vcf_path` is a VCF file; `popinfo_path` has one `sample population` pair per line.
- **Output:** A dictionary with `parser`, `populations`, `num_snps`, and `provenance`.
- **Example code:**

```python
from dnakit.population import read_vcf_data

result = read_vcf_data("variants.vcf", "popinfo.txt")
print(result["parser"], result["num_snps"], result["populations"])
```

- **Example output:**

```text
builtin 3 ['popA', 'popB']
```

## 16) `dnakit.population.bootstrap_snp_spectra()`

- **Purpose:** Bootstrap frequency spectra from chromosome chunks of an SNP file; return the mean and standard deviation.
- **Method:** Call dadi's `Misc.fragment_data_dict()` and `Misc.bootstraps_from_dd_chunks()`. SNP IDs must end in the base-pair position (for example `chr1_1500`). The file is cut into chunks of `chunk_size` base pairs, and each bootstrap draws chunks with replacement. `seed` makes the draws reproducible and restores the global random state afterwards.
- **API:** `dnakit.population.bootstrap_snp_spectra(path, *, pop_ids, projections, chunk_size, num_boot, seed=0, polarized=True)`
- **Input:** `path` is an SNP file; `pop_ids` and `projections` (the projected sample size of each population); `chunk_size` in base pairs; `num_boot` is the number of resamples (at most 1000).
- **Output:** A dictionary with `num_chunks`, `num_boot`, `mean_spectrum`, `std_spectrum`, and `provenance`.
- **Example code:**

```python
from dnakit.population import bootstrap_snp_spectra

result = bootstrap_snp_spectra(
    "chunks.txt", pop_ids=["pop1", "pop2"], projections=[2, 2], chunk_size=1000, num_boot=20, seed=5
)
print(result["num_chunks"], result["num_boot"])
print([[None if v is None else round(v, 4) for v in row] for row in result["mean_spectrum"]])
```

- **Example output:**

```text
3 20
[[None, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, None]]
```

## 17) `dnakit.population.pool_sfs_populations()`

- **Purpose:** Pool two populations of a two- or three-population spectrum into one axis.
- **Method:** Call dadi's `Misc.combine_pops()`. The two selected axes are summed as `i + j`, and the pooled axis becomes the first axis of the result. For three populations dadi supports only (0, 1), (0, 2), and (1, 2).
- **API:** `dnakit.population.pool_sfs_populations(spectrum, *, populations=(0, 1))`
- **Input:** `spectrum` is a two- or three-dimensional list; `populations` is a pair of 0-based axis indices.
- **Output:** A dictionary with `populations`, `sample_sizes`, `spectrum`, and `provenance`.
- **Example code:**

```python
from dnakit.population import pool_sfs_populations

result = pool_sfs_populations([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
print(result["sample_sizes"], result["spectrum"])
```

- **Example output:**

```text
[4] [None, 6.0, 15.0, 14.0, None]
```

## 18) `dnakit.population.read_sfs_file()`

- **Purpose:** Read a dadi array or spectrum text file.
- **Method:** Call dadi's `Numerics.array_from_file()`. The first non-comment line gives the shape and the following lines give the values. Lines starting with `#` are returned as comments.
- **API:** `dnakit.population.read_sfs_file(path)`
- **Input:** `path` is a text file.
- **Output:** A dictionary with `shape`, `spectrum`, `comments`, and `provenance`.
- **Example code:**

```python
from dnakit.population import read_sfs_file

result = read_sfs_file("spectrum.sfs")
print(result["shape"], result["spectrum"], result["comments"])
```

- **Example output:**

```text
[2, 2] [[1.5, 2.0], [3.0, 4.0]] ['demo']
```

## 19) `dnakit.population.write_sfs_file()`

- **Purpose:** Write an array or spectrum as a dadi text file.
- **Method:** Call dadi's `Numerics.array_to_file()`. The file writes the shape first, then the values; `comment_lines` are written at the top after `# `. An existing file is never overwritten.
- **API:** `dnakit.population.write_sfs_file(spectrum, destination, *, precision=16, comment_lines=())`
- **Input:** `spectrum` is a list of numbers; `destination` is a new file path; `precision` is the number of digits (1 to 17).
- **Output:** A dictionary with `destination`, `bytes`, `shape`, `num_values`, and `provenance`.
- **Example code:**

```python
from pathlib import Path
from dnakit.population import write_sfs_file

Path("spectrum_out.sfs").unlink(missing_ok=True)  # the function never overwrites an existing file
result = write_sfs_file([[1.5, 2.0], [3.0, 4.0]], "spectrum_out.sfs", comment_lines=["demo"])
print(result["shape"], result["num_values"], Path("spectrum_out.sfs").read_text().splitlines()[0])
```

- **Example output:**

```text
[2, 2] 4 # demo
```

## 20) `dnakit.population.apply_ancestral_misidentification()`

- **Purpose:** Apply ancestral-state misidentification to a spectrum.
- **Method:** Call dadi's `Numerics.apply_anc_state_misid()`. The result is `(1 - p) * fs + p * reversed(fs)`, where `p` is the fraction of sites whose ancestral state is misidentified.
- **API:** `dnakit.population.apply_ancestral_misidentification(spectrum, p_misid)`
- **Input:** `spectrum` is a spectrum list; `p_misid` is between 0 and 1.
- **Output:** A dictionary with `p_misid`, `sample_sizes`, `spectrum`, and `provenance`.
- **Example code:**

```python
from dnakit.population import apply_ancestral_misidentification

result = apply_ancestral_misidentification([0, 1, 2, 3, 0], 0.1)
print([None if v is None else round(v, 4) for v in result["spectrum"]])
```

- **Example output:**

```text
[None, 1.2, 2.0, 2.8, None]
```

## 21) `dnakit.population.expected_sfs_1d()`

- **Purpose:** Expected site frequency spectrum for a one-population demographic model.
- **Method:** Compute the expected spectrum with dadi's one-population model family. `model` is `two_epoch` (the representative), `growth`, `bottlegrowth`, `bottlegrowth_1d`, `three_epoch`, or `three_epoch_inbreeding`. `grid_points` is the integration grid; larger values approach the exact result.
- **API:** `dnakit.population.expected_sfs_1d(model, params, sample_size, *, grid_points=1000)`
- **Input:** `model`; `params` are that model's parameters in dadi's order; `sample_size` is 2 to 2000; `grid_points` is 50 to 10000.
- **Output:** A dictionary with `model`, `params`, `sample_sizes`, `spectrum` (corners None), and `provenance`.
- **Example code:**

```python
from dnakit.population import expected_sfs_1d

result = expected_sfs_1d("two_epoch", [2.0, 0.5], 10, grid_points=50)
print(result["model"], [round(v, 4) for v in result["spectrum"][1:4]])
```

- **Example output:**

```text
two_epoch [1.6053, 0.6722, 0.3969]
```

## 22) `dnakit.population.expected_sfs_2d()`

- **Purpose:** Expected joint spectrum for a two-population demographic model.
- **Method:** Computed with dadi's two-population model family. `model` defaults to `split_mig` (the representative); other names such as `asym_mig`, `anc_sym_mig`, and `IM` also work, with the same `(params, ns, pts)` interface.
- **API:** `dnakit.population.expected_sfs_2d(model, params, sample_sizes, *, grid_points=60)`
- **Input:** `model`; `params`; `sample_sizes` are the two population sample sizes (2 to 200); `grid_points` is 20 to 400.
- **Output:** A dictionary with `model`, `params`, `sample_sizes`, `spectrum`, and `provenance`.
- **Example code:**

```python
from dnakit.population import expected_sfs_2d

result = expected_sfs_2d("split_mig", [1.0, 1.0, 0.5, 1.0], [4, 4], grid_points=30)
print(result["sample_sizes"], round(result["spectrum"][1][1], 4))
```

- **Example output:**

```text
[4, 4] 0.1778
```

## 23) `dnakit.population.expected_sfs_3d()`

- **Purpose:** Expected joint spectrum for a three-population demographic model.
- **Method:** Computed with dadi's three-population model family. `model` defaults to `split_nomig` (the representative); `out_of_africa` is the Gutenkunst et al. (2009) model. The grid is kept small because the three-population integration grows with the cube of `grid_points`.
- **API:** `dnakit.population.expected_sfs_3d(model, params, sample_sizes, *, grid_points=30)`
- **Input:** `model`; `params`; `sample_sizes` are the three population sample sizes (2 to 40); `grid_points` is 20 to 60.
- **Output:** A dictionary with `model`, `params`, `sample_sizes`, `spectrum`, and `provenance`.
- **Example code:**

```python
from dnakit.population import expected_sfs_3d

result = expected_sfs_3d("split_nomig", [1.0, 1.0, 0.5, 1.0, 0.5, 0.5], [2, 2, 2], grid_points=20)
print(result["model"], round(result["spectrum"][1][1][1], 4))
```

- **Example output:**

```text
split_nomig 0.0069
```

## 24) `dnakit.population.expected_sfs_selection()`

- **Purpose:** Expected spectrum for a demographic model with selection (dadi DFE models).
- **Method:** Call dadi's `DFE.DemogSelModels` models. `model` defaults to `split_mig_sel` (the representative); the `*_single_gamma` variants share one selection coefficient across populations. `sample_sizes` has one or two entries, depending on the model.
- **API:** `dnakit.population.expected_sfs_selection(model, params, sample_sizes, *, grid_points=40)`
- **Input:** `model`; `params` (for example `nu1, nu2, T, m, gamma1, gamma2`); `sample_sizes`; `grid_points` is 20 to 400.
- **Output:** A dictionary with `model`, `params`, `sample_sizes`, `spectrum`, and `provenance`.
- **Example code:**

```python
from dnakit.population import expected_sfs_selection

result = expected_sfs_selection("split_mig_sel", [1.0, 1.0, 0.5, 1.0, 1.0, 1.0], [4, 4], grid_points=20)
print(result["model"], round(result["spectrum"][1][1], 4))
```

- **Example output:**

```text
split_mig_sel 0.2278
```

## 25) `dnakit.population.fit_demographic_model()`

- **Purpose:** Fit a demographic model to an observed spectrum by maximum likelihood.
- **Method:** Call dadi's `Inference.optimize_log()` family. `method` picks the optimizer: `bfgs` (BFGS on log parameters), `lbfgsb`, `nelder-mead`, or `powell`. With `multinom=True` dadi rescales theta, so `p0` does not include it. The log likelihood is reported at the fitted parameters.
- **API:** `dnakit.population.fit_demographic_model(data, *, model, p0, lower_bound=None, upper_bound=None, method='bfgs', multinom=True, grid_points=40, maxiter=50)`
- **Input:** `data` is the observed spectrum (one or two populations); `model` is chosen from the family matching the dimension of `data`; `p0` are positive starting values; `lower_bound` and `upper_bound` are optional.
- **Output:** A dictionary with `model`, `method`, `params`, `log_likelihood`, `grid_points`, `maxiter`, and `provenance`.
- **Example code:**

```python
from dnakit.population import expected_sfs_1d, fit_demographic_model

expected = expected_sfs_1d("two_epoch", [2.0, 0.5], 10, grid_points=50)["spectrum"]
data = [0.0 if v is None else v * 100 for v in expected]  # counts; the corners are 0
result = fit_demographic_model(
    data, model="two_epoch", p0=[1.0, 0.3], lower_bound=[0.01, 0.001], upper_bound=[10, 5], grid_points=50
)
print([round(v, 3) for v in result["params"]], round(result["log_likelihood"], 3))
```

- **Example output:**

```text
[1.999, 0.5] -23.198
```

## 26) `dnakit.population.sfs_residuals()`

- **Purpose:** Per-entry residuals between a model spectrum and observed counts.
- **Method:** Call dadi's `Inference.Anscombe_Poisson_residual()` (default) or `linear_Poisson_residual()`, which is `(model - data) / sqrt(model)`. Both arrays must have the same shape. Corners are None, and `sum_of_squares` covers the unmasked entries.
- **API:** `dnakit.population.sfs_residuals(model, data, *, residual='anscombe')`
- **Input:** `model` is the model spectrum; `data` are observed (non-negative) counts; `residual` is `anscombe` or `poisson`.
- **Output:** A dictionary with `residual`, `residuals`, `num_entries`, `sum_of_squares`, and `provenance`.
- **Example code:**

```python
from dnakit.population import sfs_residuals

result = sfs_residuals([1.0, 4.0, 9.0, 4.0, 1.0], [1.0, 2.0, 2.0, 2.0, 1.0], residual="poisson")
print([None if v is None else round(v, 4) for v in result["residuals"]], round(result["sum_of_squares"], 4))
```

- **Example output:**

```text
[None, 1.0, 2.3333, 1.0, None] 7.4444
```

## 27) `dnakit.population.parameter_uncertainty_fim()`

- **Purpose:** Parameter standard errors from the Fisher information matrix.
- **Method:** Call dadi's `Godambe.FIM_uncert()`. The model is extrapolated over `grid_points` with `Numerics.make_extrap_func()`, and the Hessian of the log likelihood by finite differences gives the Fisher information. With `multinom=True` theta is also estimated, and its standard error is the last entry.
- **API:** `dnakit.population.parameter_uncertainty_fim(data, *, model, p0, grid_points=(20, 30, 40), log=False, multinom=True, eps=0.01)`
- **Input:** `data` are observed counts; `model` is a model name; `p0` are positive parameter values; `grid_points` lists one to five grid sizes; `eps` is the relative step of the finite differences.
- **Output:** A dictionary with `params`, `grid_points`, `uncertainties`, and `provenance`.
- **Example code:**

```python
from dnakit.population import expected_sfs_1d, parameter_uncertainty_fim

expected = expected_sfs_1d("two_epoch", [2.0, 0.5], 6, grid_points=50)["spectrum"]
data = [0.0 if v is None else v * 100 for v in expected]
result = parameter_uncertainty_fim(data, model="two_epoch", p0=[2.0, 0.5])
print([round(v, 3) for v in result["uncertainties"]])
```

- **Example output:**

```text
[0.661, 0.525, 16.189]
```

## 28) `dnakit.population.annotate_snps_from_annovar()`

- **Purpose:** Keep the SNPs of one annotated variant type, using an ANNOVAR output file.
- **Method:** Call dadi's `Misc.annotate_from_annovar()`. A site is kept when its variant-type text contains `variant_type`. The kept SNP IDs are returned in sorted order.
- **API:** `dnakit.population.annotate_snps_from_annovar(path, annovar_path, *, variant_type)`
- **Input:** `path` is an SNP file; `annovar_path` is an `.exonic_variant_function` file; `variant_type` is for example `nonsynonymous`.
- **Output:** A dictionary with `num_snps_total`, `num_snps_kept`, `kept_snp_ids`, and `provenance`.
- **Example code:**

```python
from dnakit.population import annotate_snps_from_annovar

result = annotate_snps_from_annovar(
    "chunks.txt", "sample.exonic_variant_function", variant_type="nonsynonymous"
)
print(result["num_snps_total"], result["num_snps_kept"], result["kept_snp_ids"])
```

- **Example output:**

```text
6 1 ['chr1_100']
```

## 29) `dnakit.population.trinucleotide_frequency_array()`

- **Purpose:** Convert trinucleotide frequencies into dadi's 64-entry array.
- **Method:** Call dadi's `Misc.tri_freq_dict_to_array()`. dadi orders the trinucleotides with the bases `C, G, T, A` at each position, the first position varying slowest; `trinucleotides` gives that order.
- **API:** `dnakit.population.trinucleotide_frequency_array(frequencies)`
- **Input:** `frequencies` maps each of the 64 trinucleotides over A, C, G, and T to a non-negative frequency.
- **Output:** A dictionary with `trinucleotides`, `values`, and `provenance`.
- **Example code:**

```python
from dnakit.population import trinucleotide_frequency_array

bases = "ACGT"
freqs = {a + b + c: 0.0 for a in bases for b in bases for c in bases}
freqs["GTA"] = 1.0
result = trinucleotide_frequency_array(freqs)
position = result["values"].index(1.0)
print(len(result["values"]), result["trinucleotides"][position])
```

- **Example output:**

```text
64 GTA
```

## 30) `dnakit.population.expected_sfs_from_demes()`

- **Purpose:** Expected spectrum for chosen demes of a demes demographic graph file (YAML).
- **Method:** Call dadi's `Demes.SFS()`. The graph is read with the `demes` package, which must be installed separately; it is not part of the `bio-dadi` extra. `sampled_demes` names one to three demes, and `sample_sizes` gives the sample size of each; `theta` scales the spectrum.
- **API:** `dnakit.population.expected_sfs_from_demes(graph, sampled_demes, sample_sizes, *, grid_points=40, theta=1.0)`
- **Input:** `graph` is a demes YAML file; `sampled_demes` is a list of deme names; `sample_sizes` pairs with them; `grid_points` is 20 to 400 (20 to 60 for three demes).
- **Output:** A dictionary with `demes`, `sample_sizes`, `theta`, `spectrum`, and `provenance`.
- **Example code:**

```python
from dnakit.population import expected_sfs_from_demes

result = expected_sfs_from_demes("model.yaml", ["A", "B"], [4, 4], grid_points=40)
print(result["sample_sizes"], round(result["spectrum"][1][1], 4))
```

- **Example output:**

```text
[4, 4] 0.047
```

## Install optional dependencies

```bash
python -m pip install "dnakit[bio-tskit]"    # tskit functions
python -m pip install "dnakit[bio-dadi]"     # dadi functions; the cyvcf2 reader also needs cyvcf2
```

Project pages: [tskit](https://tskit.dev/tskit), [dadi](https://github.com/RyanGutenkunst/dadi), [cyvcf2](https://github.com/brentp/cyvcf2).

## Verification

The test file is `others/tests/unit/test_full_G13.py`, with 30 tests: at least one for each wrapped function.
- tskit FASTA, NEXUS, ms, and text-table output matches direct tskit calls byte for byte; `complement_intervals()` and `map_allele_indices()` match hand calculations; `tree_to_newick()` branch lengths match those from node times; `resolve_polytomies()` resolves a four-leaf star into a binary tree with 2 added nodes.
- Topology counts, `random_tree_newick()`, and PLINK `.fam` parsing match direct tskit calls.
- dadi SNP polarization, configuration counts, and VCF reading (with both the built-in and cyvcf2 readers) match direct dadi calls.
- Bootstraps, spectrum pooling, expected model spectra, model fitting, residuals, and the Fisher information match direct dadi calls; the Poisson residuals and the misidentification formula match hand calculations; the trinucleotide order matches dadi.
- Invalid parameters, existing output files, and model names outside the fixed lists raise `ConfigurationError`.
The tests make no network calls, write their inputs under `tmp_path`, and skip with `pytest.importorskip` when a package is missing.

## Current scope

This pass wraps 12 tskit functions (six family representatives; the other members, such as `write_nexus`, `tree_count_topologies`, the `generate_*` shapes, and the other `parse_*` tables, are recorded as family variants in the ledger) and 18 dadi functions. Other models and parsers in a family are selected by an argument, for example `expected_sfs_2d(model=...)`, `read_text_table(table=...)`, and `fit_demographic_model(method=...)`.
Not wrapped: tskit plotting, HTML rendering, the command-line entry points, and the metadata and provenance helpers; dadi plotting, the Triallele and TwoLocus models (they do not run on the installed dadi 2.4.4 with NumPy 2: the first lacks `dadi.tridiag`, and the second uses the removed `np.float`), the LowPass low-coverage correction, the Godambe GIM and bootstrap variants, and VCF-subsampling bootstraps.
The tree-sequence building, reading, writing, and population statistics in `dnakit.population.backends`, and the neutral spectrum functions, are recorded in the ledger as `already-wrapped` or `duplicate-of`; this pass does not wrap them again.
