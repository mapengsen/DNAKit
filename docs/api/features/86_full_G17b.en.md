# G17b functions: selene-sdk, genomic-benchmarks and kipoi

This page documents the G17b wrappers in the `predictions` area. Each function lives in its module and can be imported directly, for example `from dnakit.predictions import decode_one_hot_sequence`.

- Local functions: all 13 wrappers read only local input, such as one-hot lists, DNA strings, a local FASTA with its VCF, or a local score table. None needs the network.
- kipoi has no DNA-level function in this group that runs offline on local input, so `kipoi_functions.py` exports nothing. The kipoi listing and download functions are recorded only in the ledger. The existing `list_kipoi_models` and `describe_kipoi_model` are tested with a stand-in module and never clone the model zoo.
- Input files are declared with the `Example file` label; the other examples run as written.

## Existing and unwrapped functions

- Already wrapped: selene-sdk `sequence_to_encoding` (called by `Genome.sequence_to_encoding`, so `encode_one_hot_selene` covers it); genomic-benchmarks `info`, `is_downloaded` and `list_datasets` (`describe_genomic_benchmark`, `check_genomic_benchmark_cache`, `list_genomic_benchmark_datasets`); kipoi `list_models` and `get_model_descr` (`list_kipoi_models`, `describe_kipoi_model`).
- Duplicates: 21 module aliases of a function already covered, 15 in selene, 3 in genomic-benchmarks and 3 in kipoi. selene's `reverse_complement_sequence` does the same DNA string reverse complement as `reverse_complement_with_grelu`, so it is marked as a duplicate.
- Not wrapped in selene-sdk (75 keys): 10 plotting functions; 7 command-line and training entry points; 27 logging, configuration, run-directory, model-inference and output-writer helpers (this includes the two keys of `get_reverse_encoding`, which is not wrapped because upstream raises IndexError on one-hot input); 13 model-weight, metric and label functions that are not DNA-level; 2 deprecated classes; 3 internal helpers.
- genomic-benchmarks: 14 dataset download, dataset class and cloud-cache functions need the network; 3 functions (`labels_in_order` and two cache deletion aliases) have side effects or depend on filesystem order.
- kipoi: 88 command-line handlers; 3 functions need a singularity binary; 37 need the network or clone the model zoo; 52 are generic array, dictionary and argument helpers; 23 are internal classes and fields; 86 are plugin, environment, configuration, file and path helpers. The full status for every inventory key is in the G17b ledger.

## 1) `decode_one_hot_sequence`

- **Purpose:** Decode an (L, 4) one-hot matrix with columns A, C, G, T back to a DNA string.
- **Method:** Each row must be one-hot (one 1, the rest 0) or uniform 0.25. Valid rows go to selene-sdk's `sequence.encoding_to_sequence`, and uniform 0.25 rows give `unk_base`. Any other row raises an error instead of being decoded silently.
- **API:** `decode_one_hot_sequence(encoding[required], unk_base[optional], verbose[optional])`
- **Input:** A list of rows with four finite numbers each. `unk_base` is A, C, G, T or N; the default is N.
- **Output:** A dictionary with `length`, `sequence`, `unk_base` and `provenance`.

**Example code:**

```python
from dnakit.predictions import decode_one_hot_sequence

rows = [[1, 0, 0, 0], [0, 0, 0, 1], [0.25, 0.25, 0.25, 0.25], [0, 1, 0, 0]]
result = decode_one_hot_sequence(rows, verbose=False)
print(result["sequence"], result["length"])
```

**Example output:**

```text
ATNC 4
```

## 2) `reverse_complement_one_hot`

- **Purpose:** Reverse-complement an (L, 4) one-hot matrix.
- **Method:** Calls selene-sdk's `model_predict.get_reverse_complement_encoding`. It swaps the A and T columns and the C and G columns, then reverses the row order. Uniform N rows stay unchanged.
- **API:** `reverse_complement_one_hot(encoding[required], verbose[optional])`
- **Input:** As for `decode_one_hot_sequence`; each row must be one-hot or uniform 0.25.
- **Output:** A dictionary with `length`, `encoding` (the reverse-complemented matrix) and `provenance`.

**Example code:**

```python
from dnakit.predictions import (
    decode_one_hot_sequence,
    reverse_complement_one_hot,
)

rows = [[1, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]
rc = reverse_complement_one_hot(rows, verbose=False)["encoding"]
print(decode_one_hot_sequence(rc, verbose=False)["sequence"])
```

**Example output:**

```text
CGTT
```

## 3) `reverse_complement_allele`

- **Purpose:** Reverse-complement one variant allele string.
- **Method:** Calls selene-sdk's `model_predict.get_reverse_complement`. A lone `-` or `*` (deletion or missing allele) returns `*`. Otherwise only A, C, G, T and N are accepted, and `-` or `*` mixed into a longer allele raises an error.
- **API:** `reverse_complement_allele(allele[required], verbose[optional])`
- **Input:** An allele string of 1 to 256 characters; case-insensitive.
- **Output:** A dictionary with `allele` (uppercase), `reverse_complement` and `provenance`.

**Example code:**

```python
from dnakit.predictions import reverse_complement_allele

for allele in ["AC", "GATTACA", "-"]:
    print(allele, reverse_complement_allele(allele, verbose=False)["reverse_complement"])
```

**Example output:**

```text
AC GT
GATTACA TGTAATC
- *
```

## 4) `apply_one_hot_mutations`

- **Purpose:** Replace the bases at given 0-based positions of a one-hot matrix.
- **Method:** Calls selene-sdk's `model_predict.mutate_sequence`, which zeroes the row at each position and sets the new base. A position may appear only once. Out-of-range positions and bases other than A, C, G or T raise an error.
- **API:** `apply_one_hot_mutations(encoding[required], mutations[required], verbose[optional])`
- **Input:** A one-hot matrix. `mutations` is a list of `[position, base]` pairs, at most 1000.
- **Output:** A dictionary with `length`, `mutation_count`, `mutations` sorted by position, the mutated `encoding` and `provenance`.

**Example code:**

```python
from dnakit.predictions import (
    apply_one_hot_mutations,
    decode_one_hot_sequence,
)

rows = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
result = apply_one_hot_mutations(rows, [[0, "T"], [2, "A"]], verbose=False)
print(result["mutation_count"], decode_one_hot_sequence(result["encoding"], verbose=False)["sequence"])
```

**Example output:**

```text
2 TCAT
```

## 5) `list_in_silico_mutants`

- **Purpose:** List every mutant that changes exactly n bases in a window of a DNA string (in silico mutagenesis).
- **Method:** Calls selene-sdk's `model_predict.in_silico_mutagenesis_sequences`. An ACGT position has 3 alternative bases and an N position has 4. The count is computed from a combinatorial formula before generation, and a count above `max_mutants` raises an error rather than truncating the list. Positions are 0-based on the full sequence.
- **API:** `list_in_silico_mutants(sequence[required], mutate_n_bases[optional], start_position[optional], end_position[optional], max_mutants[optional], verbose[optional])`
- **Input:** An uppercase DNA string of at most 10000 characters; `mutate_n_bases` is 1 to 4; `max_mutants` defaults to 100000.
- **Output:** A dictionary with `count`, `mutants` (each a list of `[position, base]` pairs) and `provenance`. Order: by position combination, then by base order.

**Example code:**

```python
from dnakit.predictions import list_in_silico_mutants

result = list_in_silico_mutants("ACG", verbose=False)
print(result["count"], result["mutants"][0], result["mutants"][-1])
```

**Example output:**

```text
9 [[0, 'C']] [[2, 'T']]
```

## 6) `read_selene_vcf`

- **Purpose:** Read VCF records and check chromosome, coordinate range and REF base against a local FASTA.
- **Method:** Opens the FASTA with selene-sdk's `Genome`; pyfaidx writes a `.fai` index next to it. Then calls `model_predict.read_vcf_file`. Upstream does not check REF or coordinates, so the wrapper adds `in_bounds` and `ref_matches_fasta` to each record. A VCF chromosome missing from the FASTA raises an error.
- **API:** `read_selene_vcf(vcf_path[required], fasta_path[required], max_variants[optional], verbose[optional])`
- **Input:** A VCF file whose first five columns are `#CHROM POS ID REF ALT`, at most 50 MB; a FASTA file; `max_variants` defaults to 100000.
- **Output:** A dictionary with `variant_count` and `variants`. Each record has `chromosome`, `position` (the 1-based VCF POS), `id`, `ref`, `alt`, `strand`, `in_bounds` and `ref_matches_fasta` (null when it cannot be compared).

**Example file:** `ref.fa`

```text
>chr1
ACGTACGTAC
>chr2
GGGGTTTTAA
```

**Example file:** `variants.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	3	rs1	G	T	.	.	.
chr1	4	.	A	G	.	.	.
chr2	5	.	T	C	.	.	.
```

**Example code:**

```python
from dnakit.predictions import read_selene_vcf

result = read_selene_vcf("variants.vcf", "ref.fa", verbose=False)
for v in result["variants"]:
    print(v["chromosome"], v["position"], v["ref"], v["alt"], v["ref_matches_fasta"])
```

**Example output:**

```text
chr1 3 G T True
chr1 4 A G False
chr2 5 T C True
```

## 7) `unpack_bit_packed_sequence`

- **Purpose:** Restore selene-sdk's bit-packed one-hot matrix (each byte stores 8 positions for one channel).
- **Method:** Calls selene-sdk's `samplers.dataloader.unpackbits_sequence`. A row with all four channels set is N and becomes uniform 0.25; other rows keep 0 or 1. The input must have at least `ceil(sequence_length / 8)` rows.
- **API:** `unpack_bit_packed_sequence(packed_rows[required], sequence_length[required], verbose[optional])`
- **Input:** A list of rows with four integers from 0 to 255 each; `sequence_length` is a positive integer.
- **Output:** A dictionary with `length`, `encoding` (the first L positions) and `provenance`.

**Example code:**

```python
from dnakit.predictions import (
    decode_one_hot_sequence,
    unpack_bit_packed_sequence,
)

result = unpack_bit_packed_sequence([[136, 72, 40, 24]], 5, verbose=False)
print(decode_one_hot_sequence(result["encoding"], verbose=False)["sequence"])
```

**Example output:**

```text
ACGTN
```

## 8) `interval_length_probabilities`

- **Purpose:** Sampling probabilities for genomic intervals, weighted by interval length.
- **Method:** Calls selene-sdk's `utils.get_indices_and_probabilities`. The weight of an interval is its length divided by the total length of the selected intervals. Zero-length intervals get weight 0, are removed by upstream, and the rest are renormalised; `dropped_indices` lists them.
- **API:** `interval_length_probabilities(interval_lengths[required], indices[optional], verbose[optional])`
- **Input:** A list of non-negative integer lengths (at most 1000000 items). `indices` is a list of distinct interval indices and defaults to all. At least one selected interval must have a positive length.
- **Output:** A dictionary with `interval_count`, `indices`, `probabilities`, `dropped_indices` and `provenance`.

**Example code:**

```python
from dnakit.predictions import interval_length_probabilities

result = interval_length_probabilities([100, 300, 0, 600], verbose=False)
print(result["indices"], [round(p, 3) for p in result["probabilities"]], result["dropped_indices"])
```

**Example output:**

```text
[0, 1, 3] [0.1, 0.3, 0.6] [2]
```

## 9) `rescale_attribution_scores`

- **Purpose:** Rescale an (L, 4) attribution score matrix base by base and position by position.
- **Method:** Calls selene-sdk's `interpret.rescale_score_matrix`. `max_effect` subtracts each position's minimum across bases, then normalises across bases and multiplies by that position's range. `probability` is the identity across bases. If a divisor is zero, the result is not finite, and the wrapper raises an error instead of returning nulls.
- **API:** `rescale_attribution_scores(score_matrix[required], base_scaling[optional], position_scaling[optional], verbose[optional])`
- **Input:** Rows of four finite numbers (rows are positions; columns are A, C, G, T). Each scaling is identity, probability or max_effect; the default is identity.
- **Output:** A dictionary with `length`, the two scaling modes, `scores` and `provenance`.

**Example code:**

```python
from dnakit.predictions import rescale_attribution_scores

scores = [[1, 2, 3, 4], [0, 1, 0, 3]]
result = rescale_attribution_scores(
    scores, base_scaling="max_effect", position_scaling="max_effect", verbose=False
)
print([[round(v, 2) for v in row] for row in result["scores"]])
```

**Example output:**

```text
[[0.0, 0.5, 1.0, 1.5], [0.0, 0.75, 0.0, 2.25]]
```

## 10) `load_variant_diff_scores`

- **Purpose:** Read a variant-by-feature score table in the format written by selene-sdk variant effect runs.
- **Method:** Calls selene-sdk's `interpret.load_variant_abs_diff_scores`. The first five columns are the variant (chromosome, position, name, REF, ALT) and the rest are features. The position must be an integer, and at least one feature column is required.
- **API:** `load_variant_diff_scores(scores_path[required], max_rows[optional], verbose[optional])`
- **Input:** A tab-separated text file of at most 100 MB; `max_rows` defaults to 100000.
- **Output:** A dictionary with `variant_count`, `feature_count`, `features`, `variants` (each with `chromosome`, `position`, `name`, `ref`, `alt`), `scores` and `provenance`.

**Example file:** `scores.tsv`

```text
chrom	pos	name	ref	alt	f1	f2
chr1	10	rs1	A	G	0.5	0.25
chr2	20	rs2	C	T	0.1	0.0
```

**Example code:**

```python
from dnakit.predictions import load_variant_diff_scores

result = load_variant_diff_scores("scores.tsv", verbose=False)
print(result["variant_count"], result["features"], result["scores"][0])
```

**Example output:**

```text
2 ['f1', 'f2'] [0.5, 0.25]
```

## 11) `order_variant_labels`

- **Purpose:** Sort variants by chromosome, position, REF and ALT, and return each variant's original index.
- **Method:** Calls selene-sdk's `interpret.ordered_variants_and_indices`. The chromosome order comes from `sort_standard_chrs`: chr1 to chr22, then chrX, chrY and chrM, then all other names. Before calling upstream, the wrapper sorts the input by the same key, then by name and input index, so ties come out the same for any input order. Upstream drops variant names; `index` keeps the link to the input.
- **API:** `order_variant_labels(labels[required], verbose[optional])`
- **Input:** A list of `[chrom, pos, name, ref, alt]` items. `chrom` must start with `chr`, `pos` is a non-negative integer, and there are at most 200000 items.
- **Output:** A dictionary with `count`, `order` (the original index of each sorted item), `variants` (each with `index`, `chromosome`, `position`, `ref`, `alt`) and `provenance`.

**Example code:**

```python
from dnakit.predictions import order_variant_labels

labels = [["chr2", 30, "r1", "A", "C"], ["chr1", 20, "r2", "G", "T"], ["chr1", 10, "r3", "C", "A"]]
result = order_variant_labels(labels, verbose=False)
print(result["order"], [v["chromosome"] + ":" + str(v["position"]) for v in result["variants"]])
```

**Example output:**

```text
[2, 1, 0] ['chr1:10', 'chr1:20', 'chr2:30']
```

## 12) `sort_standard_chromosomes`

- **Purpose:** Order a set of `chr` names in selene-sdk's chromosome order.
- **Method:** Uses selene-sdk's `interpret.sort_standard_chrs` as the key: chr1 to chr22 by number, chrX as 23, chrY as 24, chrM as 25, and all others as 26. Ties within a rank are sorted by name. Names must start with `chr` and must not repeat.
- **API:** `sort_standard_chromosomes(chromosomes[required], verbose[optional])`
- **Input:** A list of distinct `chr` names, at most 10000 items.
- **Output:** A dictionary with `count` and `chromosomes` (each with `chromosome` and `rank`).

**Example code:**

```python
from dnakit.predictions import sort_standard_chromosomes

result = sort_standard_chromosomes(["chr10", "chr2", "chrX", "chr1", "chrM", "chrUn_a"], verbose=False)
print([item["chromosome"] for item in result["chromosomes"]])
```

**Example output:**

```text
['chr1', 'chr2', 'chr10', 'chrX', 'chrM', 'chrUn_a']
```

## 13) `locate_sequences_in_reference`

- **Purpose:** Find the exact matches of each FASTA sequence on reference contigs, on the forward or the reverse-complement strand.
- **Method:** Calls genomic-benchmarks' `seq2loc.fasta2loc`. The comparison is offline and needs no network. Upstream reports one location per sequence name: the last match found scanning contigs in name order, left to right. The wrapper also counts every exact match (`occurrences`) and sets `ambiguous` when that count is above 1. It rejects two inputs that would overwrite each other upstream: identical sequences, and sequences that are reverse complements of each other. A palindromic sequence is reported on the `-` strand.
- **API:** `locate_sequences_in_reference(fasta_path[required], reference[required], max_reference_length[optional], max_sequences[optional], max_query_length[optional], verbose[optional])`
- **Input:** A FASTA file with distinct sequence names and only uppercase A, C, G, T and N. `reference` maps contig names to DNA strings; the total length defaults to 2000000 bp. `max_sequences` defaults to 100 and `max_query_length` to 100000.
- **Output:** A dictionary with `query_count`, `contig_count`, `matched_count`, `matches` (each with `sequence`, `chromosome`, `start`, `end`, `strand`, `occurrences`, `ambiguous`), `unmatched` and `messages`. Coordinates are 0-based and half-open. `-` means the reverse complement of the query is at that interval.

**Example file:** `queries.fa`

```text
>unique
GGGG
>repeat
GTAC
>absent
AAAAA
```

**Example code:**

```python
from dnakit.predictions import locate_sequences_in_reference

reference = {"chr1": "ACGTACGTAC", "chr2": "GGGGTTTTAA"}
result = locate_sequences_in_reference("queries.fa", reference, verbose=False)
for match in result["matches"]:
    print(match["sequence"], match["chromosome"], match["start"], match["end"],
          match["strand"], match["occurrences"])
print(result["unmatched"])
```

**Example output:**

```text
repeat chr1 6 10 - 2
unique chr2 0 4 + 1
['absent']
```
