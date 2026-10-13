# 85 G17a: offline DNA utilities from alphagenome, grelu and enformer-pytorch

This page describes the 32 functions wrapped for group G17a. They are the DNA-level public functions of alphagenome, grelu and enformer-pytorch that run offline on local inputs. Model inference, weight downloads, genome downloads and plotting are not wrapped. Every example output was produced by running the example. The wrappers live in the `alphagenome_functions`, `grelu_functions` and `enformer_pytorch_functions` modules of `dnakit.predictions`, and the examples import from those modules. The API lines give package-level names, which import from `dnakit.predictions` only after its `__init__` re-exports these modules.

## Coverage

| Package | Status | Count |
| --- | --- | ---: |
| alphagenome | internal | 18 |
| alphagenome | plot | 12 |
| alphagenome | wrapped | 9 |
| alphagenome | other | 8 |
| alphagenome | family-variant | 7 |
| alphagenome | not-dna | 6 |
| alphagenome | already-wrapped | 2 |
| alphagenome | needs-network | 2 |
| alphagenome | duplicate-of | 1 |
| grelu | internal | 38 |
| grelu | needs-network | 33 |
| grelu | other | 31 |
| grelu | wrapped | 20 |
| grelu | duplicate-of | 17 |
| grelu | plot | 14 |
| grelu | family-variant | 10 |
| grelu | already-wrapped | 8 |
| grelu | deprecated | 6 |
| grelu | cap | 1 |
| grelu | needs-binary | 1 |
| enformer_pytorch | internal | 23 |
| enformer_pytorch | other | 10 |
| enformer_pytorch | already-wrapped | 3 |
| enformer_pytorch | wrapped | 3 |
| enformer_pytorch | duplicate-of | 2 |
| enformer_pytorch | needs-network | 2 |
| enformer_pytorch | family-variant | 1 |

Status notes: `wrapped` marks a function documented here; `already-wrapped` and `duplicate-of` point to existing dnakit functions; `family-variant` is selected by a parameter of a wrapper; `needs-network` requires a download; `other` entries whose reason starts with needs-data need a trained model or external data.

## 1) Left-align an indel (alphagenome normalize_variant)

- **Purpose:** Left-aligns an insertion or deletion to its leftmost equivalent form, so variants from different sources compare equal.
- **Method:** Supplies the reference bases from a local window to alphagenome's normalize_variant; SNVs are returned unchanged.
- **API:** `dnakit.predictions.normalize_alphagenome_variant(variant, *, reference_sequence, reference_start)`
- **Input:** variant as `chr:pos:ref>alt` (1-based); reference_sequence covers reference_start (0-based) onward and must include the base before the variant and the reference span.
- **Output:** Returns input, normalized and changed fields; positions are 1-based.

- **Example code:**

```python
from dnakit.predictions import normalize_alphagenome_variant

result = normalize_alphagenome_variant(
    "chr1:4:AA>A", reference_sequence="ACAAAG", reference_start=0
)
print(result["normalized"])
print(result["changed"])
```

- **Example output:**

```text
{'chromosome': 'chr1', 'position': 2, 'reference': 'CA', 'alternate': 'C'}
True
```

## 2) Interval intersection and union (alphagenome interval algebra)

- **Purpose:** Intersects or unions two sets of genomic intervals.
- **Method:** Calls intersect_intervals or union_intervals; each side is first checked for overlaps on the same chromosome.
- **API:** `dnakit.predictions.combine_alphagenome_intervals(left, right, *, operation='intersect')`
- **Input:** left and right are lists of `chr:start-end` strings (0-based half-open); operation is intersect or union.
- **Output:** Returns intervals sorted by chromosome and coordinates, with counts.

- **Example code:**

```python
from dnakit.predictions import combine_alphagenome_intervals

left = ["chr1:0-10", "chr1:20-30"]
right = ["chr1:5-25"]
for mode in ("intersect", "union"):
    result = combine_alphagenome_intervals(left, right, operation=mode)
    print(mode, [(item["start"], item["end"]) for item in result["intervals"]])
```

- **Example output:**

```text
intersect [(5, 10), (20, 25)]
union [(0, 30)]
```

## 3) Extract transcription start sites from a GTF (TSS)

- **Purpose:** Gets the 5' end point of each transcript or gene from a GENCODE GTF.
- **Method:** Reads transcript or gene rows and calls alphagenome's extract_tss: the start on + strand, the exclusive end on - strand.
- **API:** `dnakit.predictions.extract_alphagenome_tss(gtf_path, *, feature='transcript', max_records=1000000)`
- **Input:** gtf_path is a local GTF (gzip allowed); feature is transcript or gene; max_records caps the rows read.
- **Output:** Returns tss records with chromosome, 0-based point position, strand, gene and transcript IDs.
- **Example file:** `small.gtf`

```text
##description: small annotation for G17a tests
chr1	HAVANA	gene	1001	2000	.	+	.	gene_id "ENSG1.1"; gene_type "protein_coding"; gene_name "GENEA";
chr1	HAVANA	transcript	1001	1800	.	+	.	gene_id "ENSG1.1"; transcript_id "ENST1.1"; gene_type "protein_coding"; gene_name "GENEA"; transcript_type "protein_coding"; transcript_support_level "1"; tag "basic MANE_Select";
chr1	HAVANA	transcript	1101	2000	.	+	.	gene_id "ENSG1.1"; transcript_id "ENST2.1"; gene_type "protein_coding"; gene_name "GENEA"; transcript_type "protein_coding"; transcript_support_level "3"; tag "basic";
chr1	HAVANA	gene	5001	5500	.	-	.	gene_id "ENSG2.1"; gene_type "lncRNA"; gene_name "GENEB";
chr1	HAVANA	transcript	5001	5500	.	-	.	gene_id "ENSG2.1"; transcript_id "ENST3.1"; gene_type "lncRNA"; gene_name "GENEB"; transcript_type "lncRNA"; transcript_support_level "NA"; tag "basic";
```

- **Example code:**

```python
from dnakit.predictions import extract_alphagenome_tss

result = extract_alphagenome_tss("small.gtf")
print([item["position"] for item in result["tss"]])
print([item["strand"] for item in result["tss"]])
```

- **Example output:**

```text
[1000, 1100, 5500]
['+', '+', '-']
```

## 4) Filter a GTF by transcript rules (protein-coding, longest, MANE, TSL, type)

- **Purpose:** Filters GENCODE annotation to protein-coding, longest-transcript, MANE Select, support-level or transcript-type rows.
- **Method:** Calls the matching alphagenome filter; method selects the rule, with protein_coding as the default.
- **API:** `dnakit.predictions.filter_alphagenome_transcripts(gtf_path, *, method='protein_coding', transcript_types=None, support_levels=None, include_gene_entries=False, max_records=1000000)`
- **Input:** gtf_path is a local GTF; method is one of protein_coding, transcript_type, longest, mane_select, support_level; transcript_types and support_levels go with their rule.
- **Output:** Returns the kept row count, feature counts, and sorted gene and transcript IDs.
- **Example file:** `small.gtf`

```text
##description: small annotation for G17a tests
chr1	HAVANA	gene	1001	2000	.	+	.	gene_id "ENSG1.1"; gene_type "protein_coding"; gene_name "GENEA";
chr1	HAVANA	transcript	1001	1800	.	+	.	gene_id "ENSG1.1"; transcript_id "ENST1.1"; gene_type "protein_coding"; gene_name "GENEA"; transcript_type "protein_coding"; transcript_support_level "1"; tag "basic MANE_Select";
chr1	HAVANA	transcript	1101	2000	.	+	.	gene_id "ENSG1.1"; transcript_id "ENST2.1"; gene_type "protein_coding"; gene_name "GENEA"; transcript_type "protein_coding"; transcript_support_level "3"; tag "basic";
chr1	HAVANA	gene	5001	5500	.	-	.	gene_id "ENSG2.1"; gene_type "lncRNA"; gene_name "GENEB";
chr1	HAVANA	transcript	5001	5500	.	-	.	gene_id "ENSG2.1"; transcript_id "ENST3.1"; gene_type "lncRNA"; gene_name "GENEB"; transcript_type "lncRNA"; transcript_support_level "NA"; tag "basic";
```

- **Example code:**

```python
from dnakit.predictions import filter_alphagenome_transcripts

for method in ("protein_coding", "longest", "mane_select"):
    result = filter_alphagenome_transcripts("small.gtf", method=method)
    print(method, result["transcript_ids"])
```

- **Example output:**

```text
protein_coding ['ENST1.1', 'ENST2.1']
longest ['ENST2.1', 'ENST3.1']
mane_select ['ENST1.1']
```

## 5) Upgrade Ensembl transcript and gene ID versions

- **Purpose:** Maps old Ensembl IDs to newer versions, for example ENST00010.1 to ENST00010.3.
- **Method:** Matches on the version-free ID and calls alphagenome's upgrade_annotation_ids; unmatched IDs map to None.
- **API:** `dnakit.predictions.upgrade_alphagenome_annotation_ids(old_ids, new_ids, *, patchless=False)`
- **Input:** old_ids and new_ids are non-empty string lists; both need versions unless patchless=True for unversioned old IDs.
- **Output:** Returns mappings in input order, plus matched and total counts.

- **Example code:**

```python
from dnakit.predictions import upgrade_alphagenome_annotation_ids

result = upgrade_alphagenome_annotation_ids(
    ["ENST00010.1", "ENST00011.1"], ["ENST00010.3", "ENST00012.1"]
)
print(result["mappings"])
print(result["matched_count"])
```

- **Example output:**

```text
[{'old_id': 'ENST00010.1', 'new_id': 'ENST00010.3'}, {'old_id': 'ENST00011.1', 'new_id': None}]
1
```

## 6) List in-silico mutagenesis (ISM) single-base variants

- **Purpose:** Enumerates every single-base substitution in an interval for in-silico mutagenesis.
- **Method:** Calls alphagenome's ism_variants; reference bases are skipped, ordered by position then vocabulary.
- **API:** `dnakit.predictions.enumerate_alphagenome_ism_variants(interval, sequence, *, vocabulary='ACGT', skip_n=False)`
- **Input:** interval as `chr:start-end`; sequence is DNA of the same width; vocabulary defaults to ACGT; skip_n controls N positions.
- **Output:** Returns variants (1-based positions), the count and interval details.

- **Example code:**

```python
from dnakit.predictions import enumerate_alphagenome_ism_variants

result = enumerate_alphagenome_ism_variants("chr1:100-103", "ACG")
print(result["count"])
print(result["variants"][0])
```

- **Example output:**

```text
9
{'chromosome': 'chr1', 'position': 101, 'reference': 'A', 'alternate': 'C'}
```

## 7) Build an ISM position-by-base matrix from scores

- **Purpose:** Arranges one score per ISM variant into a (position, base) matrix.
- **Method:** Generates variants in enumerate order and calls alphagenome's ism_matrix; each position is centred on its mean.
- **API:** `dnakit.predictions.build_alphagenome_ism_matrix(interval, sequence, scores, *, vocabulary='ACGT', multiply_by_sequence=True)`
- **Input:** interval, sequence (A, C, G, T only) and scores (one per ISM variant, 3 times the width); with multiply_by_sequence only reference-base entries stay non-zero.
- **Output:** Returns a matrix with rows as positions and columns A, C, G, T.

- **Example code:**

```python
from dnakit.predictions import build_alphagenome_ism_matrix

result = build_alphagenome_ism_matrix("chr1:100-101", "A", [0.5, 0.25, 0.75])
print(result["matrix"])
full = build_alphagenome_ism_matrix(
    "chr1:100-101", "A", [0.5, 0.25, 0.75], multiply_by_sequence=False
)
print(full["matrix"])
```

- **Example output:**

```text
[[-0.5, 0.0, 0.0, 0.0]]
[[-0.5, 0.0, -0.25, 0.25]]
```

## 8) Check a sequence length against AlphaGenome's supported lengths

- **Purpose:** Checks whether a length is one of the model's supported sizes (16 KB, 100 KB, 500 KB, 1 MB).
- **Method:** Calls alphagenome's validate_sequence_length; no model is loaded and no request is sent.
- **API:** `dnakit.predictions.validate_alphagenome_sequence_length(sequence_length)`
- **Input:** sequence_length is a positive integer.
- **Output:** Returns supported_lengths when valid; unsupported lengths raise ConfigurationError.

- **Example code:**

```python
from dnakit.predictions import validate_alphagenome_sequence_length

result = validate_alphagenome_sequence_length(16384)
print(result["supported_lengths"])
```

- **Example output:**

```text
[16384, 131072, 524288, 1048576]
```

## 9) Read AlphaGenome data-fold intervals

- **Purpose:** Returns the training, validation or test intervals for a model version from a local regions file.
- **Method:** Calls alphagenome's get_fold_intervals with a local file, so the default download is not used.
- **API:** `dnakit.predictions.get_alphagenome_fold_intervals(regions_path, *, model_version='fold_0', subset='test', organism='homo_sapiens')`
- **Input:** regions_path is a headerless tab-separated BED (chromosome, start, end, fold); model_version is all_folds or fold_0 to fold_3; subset is train, valid or test; organism defaults to homo_sapiens.
- **Output:** Returns intervals, fold_names and the count.
- **Example file:** `regions.bed`

```text
chr1	0	100	fold0
chr2	200	300	fold1
chr3	400	500	fold2
chr1	600	700	fold0
```

- **Example code:**

```python
from dnakit.predictions import get_alphagenome_fold_intervals

result = get_alphagenome_fold_intervals("regions.bed", model_version="fold_0", subset="valid")
print([(item["chromosome"], item["start"]) for item in result["intervals"]])
print(result["fold_names"])
```

- **Example output:**

```text
[('chr1', 0), ('chr1', 600)]
['fold0']
```

## 10) Convert DNA encodings (indices, one-hot, strings)

- **Purpose:** Converts a batch of sequences between gReLU's index, one-hot and string forms.
- **Method:** Calls the matching gReLU format function for source and target; indices are A=0, C=1, G=2, T=3, N=4, and N is an all-zero one-hot column.
- **API:** `dnakit.predictions.convert_dna_encoding(values, *, source, target)`
- **Input:** values is a batch list; one-hot is a 4 x L matrix (rows A, C, G, T); sequences in a batch share a length. For string input use dnakit's one_hot_encode_dna.
- **Output:** Returns values (strings, indices or one-hot lists) and count.

- **Example code:**

```python
from dnakit.predictions import convert_dna_encoding

one_hot = [[[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]]]
print(convert_dna_encoding(one_hot, source="one_hot", target="strings")["values"])
print(convert_dna_encoding([[0, 3]], source="indices", target="one_hot")["values"])
```

- **Example output:**

```text
['ACGTN']
[[[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 1.0]]]
```

## 11) Check whether DNA encodings are valid

- **Purpose:** Tells whether a batch of strings, indices or one-hot matrices is a valid DNA encoding.
- **Method:** Uses gReLU's checks and adds the rules that each one-hot column has at most one 1 and indices lie in 0 to 4; strings are uppercased first.
- **API:** `dnakit.predictions.check_dna_encoding(values, *, encoding='strings')`
- **Input:** values is the batch; encoding is strings, indices or one_hot. Malformed containers raise ConfigurationError.
- **Output:** Returns valid (boolean), encoding and count.

- **Example code:**

```python
from dnakit.predictions import check_dna_encoding

print(check_dna_encoding(["ACGTN", "acgt"], encoding="strings")["valid"])
print(check_dna_encoding(["ACGX"], encoding="strings")["valid"])
print(check_dna_encoding([[[1, 1], [1, 0], [0, 0], [0, 0]]], encoding="one_hot")["valid"])
```

- **Example output:**

```text
True
False
False
```

## 12) GC content and GC histogram

- **Purpose:** Computes each sequence's GC fraction, or the GC distribution across a set of sequences.
- **Method:** mode=fraction calls gReLU's gc; mode=distribution calls gc_distribution with bins of width binwidth starting at 0.
- **API:** `dnakit.predictions.gc_content_distribution(sequences, *, mode='distribution', binwidth=0.1, normalize=False)`
- **Input:** sequences is a DNA string or list; binwidth is in (0, 1]; with normalize the counts sum to 1.
- **Output:** fraction mode returns gc values; distribution mode returns bin_edges and counts.

- **Example code:**

```python
from dnakit.predictions import gc_content_distribution

seqs = ["ACGT", "GGGG", "AAAA"]
print(gc_content_distribution(seqs, mode="fraction")["gc"])
hist = gc_content_distribution(seqs, mode="distribution", binwidth=0.25)
print(hist["bin_edges"], hist["counts"])
```

- **Example output:**

```text
[0.5, 1.0, 0.0]
[0.0, 0.25, 0.5, 0.75] [1.0, 0.0, 1.0, 1.0]
```

## 13) Summarize sequence lengths

- **Purpose:** Gives each sequence's length and whether all sequences share one length.
- **Method:** Calls gReLU's get_lengths, check_equal_lengths and, when lengths match, get_unique_length.
- **API:** `dnakit.predictions.summarize_sequence_lengths(sequences)`
- **Input:** sequences is a DNA string or list.
- **Output:** Returns lengths, all_equal and unique_length (None when lengths differ).

- **Example code:**

```python
from dnakit.predictions import summarize_sequence_lengths

print(summarize_sequence_lengths(["ACG", "TTT"]))
print(summarize_sequence_lengths(["ACG", "TT"])["unique_length"])
```

- **Example output:**

```text
{'status': 'completed', 'operation': 'summarize_sequence_lengths', 'count': 2, 'lengths': [3, 3], 'all_equal': True, 'unique_length': 3, 'provenance': {'package': 'grelu', 'version': '1.1.0', 'length_unit': 'bases'}}
None
```

## 14) Positions where two equal-length sequences differ

- **Purpose:** Lists the positions where two equal-length DNA sequences differ, with both bases.
- **Method:** Calls gReLU's seq_differences without its printed output, then adds the reference and alternate bases.
- **API:** `dnakit.predictions.dna_sequence_differences(sequence_a, sequence_b)`
- **Input:** sequence_a and sequence_b are DNA strings of equal length.
- **Output:** Returns differences (0-based position, reference base, alternate base) and count.

- **Example code:**

```python
from dnakit.predictions import dna_sequence_differences

print(dna_sequence_differences("ACGT", "ACCT")["differences"])
```

- **Example output:**

```text
[{'position': 2, 'reference': 'G', 'alternate': 'C'}]
```

## 15) Extend intervals to a fixed length around a summit

- **Purpose:** Builds fixed-length intervals centred on a peak summit (for example narrowPeak).
- **Method:** Calls gReLU's extend_from_coord: start = start + summit - seq_len // 2, with summit an offset from start.
- **API:** `dnakit.predictions.extend_intervals_from_centers(intervals, *, seq_len, center_col='summit')`
- **Input:** intervals are mappings with chrom, start, end and summit; seq_len is a positive integer; a start below 0 raises an error.
- **Output:** Returns intervals (chrom, start, end) and count.

- **Example code:**

```python
from dnakit.predictions import extend_intervals_from_centers

rows = [{"chrom": "chr1", "start": 100, "end": 200, "summit": 50}]
print(extend_intervals_from_centers(rows, seq_len=10)["intervals"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'start': 145, 'end': 155}]
```

## 16) Clip intervals to bounds

- **Purpose:** Raises interval starts to start and lowers interval ends to end.
- **Method:** Calls gReLU's clip_intervals; rows whose start is not below their end are kept and counted in invalid_count.
- **API:** `dnakit.predictions.clip_intervals_to_bounds(intervals, *, start=None, end=None)`
- **Input:** intervals are mappings; start and end are optional.
- **Output:** Returns intervals, count and invalid_count.

- **Example code:**

```python
from dnakit.predictions import clip_intervals_to_bounds

rows = [{"chrom": "chr1", "start": 10, "end": 50}, {"chrom": "chr2", "start": 200, "end": 300}]
result = clip_intervals_to_bounds(rows, start=0, end=150)
print(result["intervals"])
print(result["invalid_count"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'start': 10, 'end': 50}, {'chrom': 'chr2', 'start': 200, 'end': 150}]
1
```

## 17) Merge intervals by group column

- **Purpose:** Merges intervals that share a group value into one span from the smallest start to the largest end.
- **Method:** Calls gReLU's merge_intervals_by_column; a group spanning two chromosomes raises an error.
- **API:** `dnakit.predictions.merge_intervals_by_group(intervals, *, group_col)`
- **Input:** intervals are mappings that include the group_col column.
- **Output:** Returns one interval per group and the count.

- **Example code:**

```python
from dnakit.predictions import merge_intervals_by_group

rows = [
    {"chrom": "chr1", "start": 0, "end": 10, "grp": "a"},
    {"chrom": "chr1", "start": 5, "end": 20, "grp": "a"},
    {"chrom": "chr2", "start": 0, "end": 5, "grp": "b"},
]
print(merge_intervals_by_group(rows, group_col="grp")["intervals"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'start': 0, 'end': 20, 'grp': 'a'}, {'chrom': 'chr2', 'start': 0, 'end': 5, 'grp': 'b'}]
```

## 18) Keep or drop intervals overlapping reference intervals

- **Purpose:** Filters intervals by their overlap with a reference set.
- **Method:** Calls gReLU's filter_overlapping; window widens each reference interval on both sides.
- **API:** `dnakit.predictions.filter_intervals_overlapping_reference(intervals, reference, *, window=0, invert=False, method='any')`
- **Input:** intervals and reference are mapping lists; method is any or all; invert keeps the non-overlapping intervals.
- **Output:** Returns the kept intervals, counts and upstream messages.

- **Example code:**

```python
from dnakit.predictions import filter_intervals_overlapping_reference

intervals = [{"chrom": "chr1", "start": 0, "end": 50}, {"chrom": "chr1", "start": 100, "end": 150}]
reference = [{"chrom": "chr1", "start": 40, "end": 60}]
print(filter_intervals_overlapping_reference(intervals, reference)["intervals"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'start': 0, 'end': 50}]
```

## 19) Split intervals into train, validation and test by chromosome

- **Purpose:** Splits intervals by chromosome lists; defaults are chr10 for validation and chr11 for test.
- **Method:** Calls gReLU's split; train_chroms=None gives training every other chromosome; the three lists must be disjoint.
- **API:** `dnakit.predictions.split_intervals_by_chromosome(intervals, *, train_chroms=None, val_chroms=('chr10',), test_chroms=('chr11',))`
- **Input:** intervals are mappings; train_chroms, val_chroms and test_chroms are chromosome name lists.
- **Output:** Returns the train, validation and test sets with counts.

- **Example code:**

```python
from dnakit.predictions import split_intervals_by_chromosome

rows = [
    {"chrom": "chr1", "start": 0, "end": 50},
    {"chrom": "chr10", "start": 0, "end": 50},
    {"chrom": "chr11", "start": 0, "end": 50},
]
result = split_intervals_by_chromosome(rows)
print(result["counts"])
print([item["chrom"] for item in result["test"]])
```

- **Example output:**

```text
{'train': 1, 'validation': 1, 'test': 1}
['chr11']
```

## 20) Remove ENCODE blacklist regions

- **Purpose:** Drops intervals overlapping the ENCODE blacklist (hg19, hg38 and mm10 ship with gReLU; no download).
- **Method:** Calls gReLU's filter_blacklist with the bundled file; window pads each blacklist region.
- **API:** `dnakit.predictions.filter_intervals_blacklist(intervals, *, genome='hg38', window=0)`
- **Input:** intervals are mappings; genome is hg19, hg38 or mm10.
- **Output:** Returns the kept intervals, input_count, count and removed_count.

- **Example code:**

```python
from dnakit.predictions import filter_intervals_blacklist

rows = [
    {"chrom": "chr1", "start": 0, "end": 500},
    {"chrom": "chr1", "start": 1000, "end": 1500},
    {"chrom": "chr2", "start": 5000, "end": 5500},
]
result = filter_intervals_blacklist(rows, genome="hg38")
print(result["intervals"])
print(result["removed_count"])
```

- **Example output:**

```text
[{'chrom': 'chr2', 'start': 5000, 'end': 5500}]
2
```

## 21) Include or exclude intervals by chromosome

- **Purpose:** Keeps intervals on the listed chromosomes, or drops the listed ones.
- **Method:** Calls gReLU's filter_chromosomes; give include, exclude, or both.
- **API:** `dnakit.predictions.filter_intervals_by_chromosome(intervals, *, include=None, exclude=None)`
- **Input:** intervals are mappings; include and exclude are chromosome name lists.
- **Output:** Returns the kept intervals and counts.

- **Example code:**

```python
from dnakit.predictions import filter_intervals_by_chromosome

rows = [{"chrom": "chr1", "start": 0, "end": 10}, {"chrom": "chr2", "start": 0, "end": 10}]
print(filter_intervals_by_chromosome(rows, include=["chr2"])["intervals"])
```

- **Example output:**

```text
[{'chrom': 'chr2', 'start': 0, 'end': 10}]
```

## 22) Consensus sequence of a PWM

- **Purpose:** Turns a 4 x L position weight matrix into a consensus DNA string.
- **Method:** Calls gReLU's motifs_to_strings: each position takes its highest base; ties go to the earlier base in A, C, G, T.
- **API:** `dnakit.predictions.pwm_to_consensus(pwm)`
- **Input:** pwm is four numeric rows (A, C, G, T).
- **Output:** Returns consensus and length.

- **Example code:**

```python
from dnakit.predictions import pwm_to_consensus

pwm = [
    [0.9, 0.05, 0.02],
    [0.05, 0.9, 0.03],
    [0.03, 0.02, 0.9],
    [0.02, 0.03, 0.05],
]
print(pwm_to_consensus(pwm)["consensus"])
```

- **Example output:**

```text
ACG
```

## 23) Trim low-information PWM flanks

- **Purpose:** Removes weak columns from both ends of a PWM.
- **Method:** Calls gReLU's trim_pwm: a column scores the sum of absolute values and is kept when its score is at least trim_threshold times the largest.
- **API:** `dnakit.predictions.trim_motif_pwm(pwm, *, trim_threshold=0.3)`
- **Input:** pwm is a 4 x L numeric matrix; use log-odds or centred values, since probability columns all score 1 and are not trimmed.
- **Output:** Returns start and end (0-based half-open), the trimmed pwm and input_length.

- **Example code:**

```python
from dnakit.predictions import trim_motif_pwm

pwm = [
    [0.0, 0.65, -0.2, -0.2, 0.0],
    [0.0, -0.2, 0.65, -0.2, 0.0],
    [0.0, -0.2, -0.2, 0.65, 0.0],
    [0.0, -0.25, -0.25, -0.25, 0.0],
]
result = trim_motif_pwm(pwm)
print(result["start"], result["end"])
```

- **Example output:**

```text
1 4
```

## 24) Scan DNA with several PWMs (FIMO-style)

- **Purpose:** Scans DNA on both strands for motif sites passing a p-value threshold.
- **Method:** Calls gReLU's scan_sequences and filters by pthresh; hits are sorted by sequence, start, strand and motif.
- **API:** `dnakit.predictions.scan_motifs_with_pwms(sequences, motifs, *, pthresh=0.001, rc=True)`
- **Input:** sequences is an ID-to-sequence mapping or a list (IDs seq0, seq1 ...); motifs maps names to 4 x L matrices; rc controls the reverse-complement scan.
- **Output:** Returns hits (motif, sequence, start, end, strand, score, p_value, matched_sequence) and the count.

- **Example code:**

```python
from dnakit.predictions import scan_motifs_with_pwms

motif = [[0.9, 0.05, 0.03], [0.05, 0.9, 0.03], [0.03, 0.03, 0.9], [0.02, 0.02, 0.04]]
result = scan_motifs_with_pwms({"s1": "TTTTACGTTT"}, {"M1": motif}, pthresh=0.01)
for hit in result["hits"]:
    print(hit["start"], hit["strand"], hit["matched_sequence"])
```

- **Example output:**

```text
4 + ACG
5 - CGT
```

## 25) Compare PWM sites between reference and alternate alleles

- **Purpose:** Shows how motif scores change between a reference and an alternate sequence.
- **Method:** Calls gReLU's compare_motifs; each site gives reference and alternate p-values and scores, and the alternate-minus-reference difference.
- **API:** `dnakit.predictions.compare_motif_pwm_alleles(reference, alternate, motifs, *, pthresh=0.001, rc=True)`
- **Input:** reference and alternate are DNA strings; motifs maps names to 4 x L matrices.
- **Output:** Returns sites and the count; missing values are None.

- **Example code:**

```python
from dnakit.predictions import compare_motif_pwm_alleles

motif = [[0.9, 0.05, 0.03], [0.05, 0.9, 0.03], [0.03, 0.03, 0.9], [0.02, 0.02, 0.04]]
result = compare_motif_pwm_alleles("TTTTACGTTT", "TTTTATGTTT", {"M1": motif}, pthresh=0.01)
for site in result["sites"]:
    print(site["start"], site["strand"], round(site["score_difference"], 4))
```

- **Example output:**

```text
4 + -5.4848
5 - -14.8755
```

## 26) Filter variants by allele type and indel length

- **Purpose:** Keeps variants with standard bases and limits insertion and deletion lengths.
- **Method:** Calls gReLU's filter_variants; each row needs chrom, pos (1-based), ref and alt.
- **API:** `dnakit.predictions.filter_variants_by_length(variants, *, standard_bases=True, max_insert_len=0, max_del_len=0)`
- **Input:** variants are mappings; with standard_bases only A, C, G, T and N alleles are kept; max_insert_len and max_del_len are non-negative integers.
- **Output:** Returns the kept variants, input_count, count and removed_count.

- **Example code:**

```python
from dnakit.predictions import filter_variants_by_length

variants = [
    {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "T"},
    {"chrom": "chr1", "pos": 20, "ref": "AC", "alt": "A"},
    {"chrom": "chr2", "pos": 30, "ref": "G", "alt": "GTT"},
]
result = filter_variants_by_length(variants)
print(result["variants"])
print(result["removed_count"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'T'}]
2
```

## 27) Fixed-length intervals around variants

- **Purpose:** Builds a window of fixed length centred on each variant, for extracting sequence context.
- **Method:** Calls gReLU's variants_to_intervals: start = pos - ceil(seq_len / 2) and end = start + seq_len.
- **API:** `dnakit.predictions.variants_to_sequence_intervals(variants, *, seq_len=1)`
- **Input:** variants are mappings (chrom, pos, ref, alt); seq_len is a positive integer, default 1.
- **Output:** Returns intervals (chrom, start, end) and count.

- **Example code:**

```python
from dnakit.predictions import variants_to_sequence_intervals

variants = [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "T"}]
print(variants_to_sequence_intervals(variants, seq_len=10)["intervals"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'start': 5, 'end': 15}]
```

## 28) Read FASTA sequences (gReLU reader)

- **Purpose:** Reads sequences from a local FASTA or gzipped FASTA in file order.
- **Method:** Calls gReLU's read_fasta; gReLU returns no headers, so the result holds sequences and lengths.
- **API:** `dnakit.predictions.read_fasta_sequences(fasta_path, *, max_sequences=1000)`
- **Input:** fasta_path is a local file; max_sequences caps the record count and raises when exceeded.
- **Output:** Returns sequences, lengths, count and source.
- **Example file:** `toy.fa`

```text
>a
ACGTTA
>b
GGCC
```

- **Example code:**

```python
from dnakit.predictions import read_fasta_sequences

result = read_fasta_sequences("toy.fa")
print(result["sequences"])
print(result["lengths"])
```

- **Example output:**

```text
['ACGTTA', 'GGCC']
[6, 4]
```

## 29) Random subset of intervals with a seed

- **Purpose:** Draws n intervals without replacement; the draw is reproducible.
- **Method:** Calls gReLU's filter_random in place on a copy (the upstream returns None otherwise); output is sorted by coordinates.
- **API:** `dnakit.predictions.sample_intervals(intervals, *, n, seed=None)`
- **Input:** intervals are mappings; n is a positive integer no larger than the interval count; seed is an optional non-negative integer.
- **Output:** Returns intervals (sorted), input_count and seed.

- **Example code:**

```python
from dnakit.predictions import sample_intervals

rows = [{"chrom": "chr1", "start": i * 10, "end": i * 10 + 5} for i in range(6)]
result = sample_intervals(rows, n=3, seed=1)
print([item["start"] for item in result["intervals"]])
```

- **Example output:**

```text
[10, 20, 40]
```

## 30) DNA to integer indices (enformer-pytorch)

- **Purpose:** Converts a DNA string to the integer indices that enformer-pytorch uses.
- **Method:** Calls enformer_pytorch.data.str_to_seq_indices; A=0, C=1, G=2, T=3, N=4.
- **API:** `dnakit.predictions.encode_enformer_indices(sequence, *, max_length=10000)`
- **Input:** sequence is an ACGTN string; max_length caps its length.
- **Output:** Returns indices, length and the index convention.

- **Example code:**

```python
from dnakit.predictions import encode_enformer_indices

print(encode_enformer_indices("ACGTN")["indices"])
```

- **Example output:**

```text
[0, 1, 2, 3, 4]
```

## 31) Integer indices to one-hot with padding

- **Purpose:** Converts integer indices to an (L, 4) one-hot matrix, including padding positions.
- **Method:** Calls enformer_pytorch.data.seq_indices_to_one_hot; N (4) is all zeros and padding (default -1) gives 0.25 in each channel.
- **API:** `dnakit.predictions.indices_to_one_hot_enformer(indices, *, padding=-1)`
- **Input:** indices are integers 0 to 4 or the padding value; padding is a negative integer, default -1.
- **Output:** Returns one_hot (one row per position), length and the padding convention.

- **Example code:**

```python
from dnakit.predictions import indices_to_one_hot_enformer

result = indices_to_one_hot_enformer([0, 1, 4, -1])
for row in result["one_hot"]:
    print(row)
```

- **Example output:**

```text
[1.0, 0.0, 0.0, 0.0]
[0.0, 1.0, 0.0, 0.0]
[0.0, 0.0, 0.0, 0.0]
[0.25, 0.25, 0.25, 0.25]
```

## 32) Reverse complement in enformer encodings

- **Purpose:** Reverse-complements DNA in integer-index or one-hot form.
- **Method:** Calls seq_indices_reverse_complement (indices) or one_hot_reverse_complement (one_hot).
- **API:** `dnakit.predictions.reverse_complement_enformer(sequence, *, representation='indices')`
- **Input:** sequence is an ACGTN string; representation is indices or one_hot.
- **Output:** Returns reverse_complement (index list or one-hot rows) and length.

- **Example code:**

```python
from dnakit.predictions import reverse_complement_enformer

print(reverse_complement_enformer("ACGTN")["reverse_complement"])
```

- **Example output:**

```text
[4, 0, 1, 2, 3]
```
