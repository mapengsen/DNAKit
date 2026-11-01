# gtars region sets, gene annotation, signal matrices, tokenizers and LOLA databases

This page covers the gtars rows of work-list pack P24 that work on BED region sets (`gtars.models`), gene annotation (GTF), signal matrices, region tokenizers (`gtars.tokenizers`) and LOLA region databases (`gtars.lola`). Coordinates are 0-based half-open. Rows that differ only by a parameter are merged: the region-set statistics (widths, mean width, nucleotide length, maximum end, per-chromosome statistics) are one function, the two neighbour-distance rows are one function with `kind`, the gene model, annotation and partition list rows are one summary function, and token/id conversion is one function with `direction`. Rows already covered by existing functions (`region_gc_content`, `partition_region_counts`, `sequence_digest_record`) are listed at the end. gtars 0.10.0 was tested.

## Example files

- **Example file:** `a.bed`

```text
chr1	100	200
chr1	150	300
chr1	500	600
chr2	10	50
```

- **Example file:** `b.bed`

```text
chr1	180	520
chr2	40	90
chr3	5	15
```

- **Example file:** `u.bed`

```text
chr1	0	1000
chr1	1000	2000
chr2	0	500
```

- **Example file:** `q.bed`

```text
chr1	900	1000
chr1	20000	20100
```

- **Example file:** `sm.tsv`

```text
label	c1	c2
chr1_100_500	1.5	2
chr1_1000_2000	0.5	3
```

- **Example file:** `g.gtf`

```text
chr1	demo	gene	1001	4000	.	+	.	gene_id "G0"; transcript_id "T0"; gene_biotype "protein_coding"; gene_name "N0";
chr1	demo	transcript	1001	4000	.	+	.	gene_id "G0"; transcript_id "T0"; gene_biotype "protein_coding"; gene_name "N0";
chr1	demo	exon	1001	1500	.	+	.	gene_id "G0"; transcript_id "T0"; gene_biotype "protein_coding"; gene_name "N0";
chr1	demo	exon	2501	4000	.	+	.	gene_id "G0"; transcript_id "T0"; gene_biotype "protein_coding"; gene_name "N0";
chr1	demo	gene	6001	9000	.	-	.	gene_id "G1"; transcript_id "T1"; gene_biotype "protein_coding"; gene_name "N1";
chr1	demo	transcript	6001	9000	.	-	.	gene_id "G1"; transcript_id "T1"; gene_biotype "protein_coding"; gene_name "N1";
chr1	demo	exon	6001	6500	.	-	.	gene_id "G1"; transcript_id "T1"; gene_biotype "protein_coding"; gene_name "N1";
chr1	demo	exon	7501	9000	.	-	.	gene_id "G1"; transcript_id "T1"; gene_biotype "protein_coding"; gene_name "N1";
chr1	demo	gene	11001	14000	.	+	.	gene_id "G2"; transcript_id "T2"; gene_biotype "protein_coding"; gene_name "N2";
chr1	demo	transcript	11001	14000	.	+	.	gene_id "G2"; transcript_id "T2"; gene_biotype "protein_coding"; gene_name "N2";
chr1	demo	exon	11001	11500	.	+	.	gene_id "G2"; transcript_id "T2"; gene_biotype "protein_coding"; gene_name "N2";
chr1	demo	exon	12501	14000	.	+	.	gene_id "G2"; transcript_id "T2"; gene_biotype "protein_coding"; gene_name "N2";
```

- **Example file:** `c.toml`

```text
universe = "u.bed"
tokenizer_type = "ailist"
```

## 1) `dnakit.annotation.region_set_summary()`

- **Purpose:** Summarise a BED region set: number of regions, widths, mean width, total nucleotides, maximum end per chromosome and per-chromosome statistics.
- **Method:** Reads the file with `gtars.models.RegionSet` and calls `widths`, `mean_region_width`, `get_nucleotide_length`, `get_max_end_per_chr` and `chromosome_statistics`. Coordinates are 0-based half-open.
- **API:** `dnakit.annotation.region_set_summary(regions, *, include_widths=False, max_rows=100000, verbose=True)`.
- **Input:** BED path, `include_widths`, `max_rows`.
- **Output:** `region_count`, `mean_region_width`, `total_nucleotides`, `max_end_per_chromosome`, `chromosomes`, optional `widths`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_set_summary

r = region_set_summary('a.bed', include_widths=True, verbose=False)
print(r['widths'], r['mean_region_width'], r['total_nucleotides'])
print(r['max_end_per_chromosome'])
print(r['chromosomes'][0]['number_of_regions'])
```

- **Example output:**

```text
[100, 150, 100, 40] 97.5 390
{'chr1': 600, 'chr2': 50}
3
```

## 2) `dnakit.annotation.region_neighbor_distances()`

- **Purpose:** Distances between neighbouring regions on each chromosome: signed gaps between consecutive regions, or each region's distance to its nearest neighbour.
- **Method:** `RegionSet.neighbor_distances` (`kind='signed_gap'`) or `RegionSet.nearest_neighbors` (`kind='nearest'`). Chromosomes with too few regions are skipped, so the result is not aligned one to one with the input; gtars merges overlapping regions before measuring gaps.
- **API:** `dnakit.annotation.region_neighbor_distances(regions, *, kind='signed_gap', max_rows=100000, verbose=True)`.
- **Input:** BED path, `kind`, `max_rows`.
- **Output:** `distances`, `distance_count`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_neighbor_distances

print(region_neighbor_distances('a.bed', verbose=False)['distances'])
print(region_neighbor_distances('a.bed', kind='nearest', verbose=False)['distances'])
```

- **Example output:**

```text
[200]
[0, 0, 200]
```

## 3) `dnakit.annotation.region_distribution_bins()`

- **Purpose:** Count region midpoints in equal-width genomic bins.
- **Method:** `RegionSet.distribution(n_bins, chrom_sizes)`. Without `chrom_sizes` the bin width comes from the largest end in the file; with it every chromosome is cut into `n_bins` bins and regions beyond the stated size are skipped.
- **API:** `dnakit.annotation.region_distribution_bins(regions, *, n_bins=250, chrom_sizes=None, verbose=True)`.
- **Input:** BED path, `n_bins`, optional `chrom_sizes`.
- **Output:** `bins` (chrom, start, end, count, bin).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_distribution_bins

r = region_distribution_bins('a.bed', n_bins=2, chrom_sizes={'chr1': 1000, 'chr2': 100}, verbose=False)
for b in r['bins']:
    print(b)
```

- **Example output:**

```text
{'chrom': 'chr1', 'start': 0, 'end': 500, 'count': 2, 'bin': 0}
{'chrom': 'chr1', 'start': 500, 'end': 1000, 'count': 1, 'bin': 1}
{'chrom': 'chr2', 'start': 0, 'end': 100, 'count': 1, 'bin': 0}
```

## 4) `dnakit.annotation.region_promoters()`

- **Purpose:** Build promoter windows `[start - upstream, start + downstream)` around every region start.
- **Method:** `RegionSet.promoters(upstream, downstream)`. gtars measures from the region start and ignores strand; windows are clamped at 0. The result can also be written to a new BED file with `output`.
- **API:** `dnakit.annotation.region_promoters(regions, upstream, downstream, *, output=None, max_rows=100000, verbose=True)`.
- **Input:** BED path, `upstream`, `downstream`, optional `output`, `max_rows`.
- **Output:** `regions`, `region_count`, `truncated`, `output`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_promoters

r = region_promoters('a.bed', 50, 10, verbose=False)
print([(x['chrom'], x['start'], x['end']) for x in r['regions']])
```

- **Example output:**

```text
[('chr1', 50, 110), ('chr1', 100, 160), ('chr1', 450, 510), ('chr2', 0, 20)]
```

## 5) `dnakit.annotation.region_pairwise_operation()`

- **Purpose:** Combine two BED region sets by index-wise intersection, concatenation, or all-vs-all intersection.
- **Method:** `operation` selects `RegionSet.pintersect` (i-th region of each set), `concat` (all regions, nothing merged) or `intersect_all` (every overlapping pair; overlapping regions of one set are merged first by gtars).
- **API:** `dnakit.annotation.region_pairwise_operation(regions_a, regions_b, *, operation='pintersect', output=None, max_rows=100000, verbose=True)`.
- **Input:** Two BED paths, `operation`, optional `output`, `max_rows`.
- **Output:** `regions`, `region_count`, `truncated`, `output`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_pairwise_operation

for op in ('pintersect', 'concat', 'intersect_all'):
    r = region_pairwise_operation('a.bed', 'b.bed', operation=op, verbose=False)
    print(op, r['region_count'], [(x['chrom'], x['start'], x['end']) for x in r['regions']][:3])
```

- **Example output:**

```text
pintersect 3 [('chr1', 180, 200), ('chr1', 150, 150), ('chr1', 500, 500)]
concat 7 [('chr1', 100, 200), ('chr1', 150, 300), ('chr1', 500, 600)]
intersect_all 3 [('chr1', 180, 300), ('chr1', 500, 520), ('chr2', 40, 50)]
```

## 6) `dnakit.annotation.region_overlap_coefficient()`

- **Purpose:** Nucleotide overlap coefficient: shared bases divided by the bases of the smaller set.
- **Method:** `RegionSet.overlap_coefficient`.
- **API:** `dnakit.annotation.region_overlap_coefficient(regions_a, regions_b, *, verbose=True)`.
- **Input:** Two BED paths.
- **Output:** `overlap_coefficient` in [0, 1].
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_overlap_coefficient

print(region_overlap_coefficient('a.bed', 'b.bed', verbose=False)['overlap_coefficient'])
```

- **Example output:**

```text
0.4411764705882353
```

## 7) `dnakit.annotation.region_subset_by_overlaps()`

- **Purpose:** Keep the regions of the first set that overlap at least one region of the second set.
- **Method:** `RegionSet.subset_by_overlaps` (an AIList index of the second set is queried by every region of the first). Regions are kept unchanged.
- **API:** `dnakit.annotation.region_subset_by_overlaps(regions_a, regions_b, *, output=None, max_rows=100000, verbose=True)`.
- **Input:** Two BED paths, optional `output`, `max_rows`.
- **Output:** `regions`, `region_count`, `input_region_count`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_subset_by_overlaps

r = region_subset_by_overlaps('a.bed', 'b.bed', verbose=False)
print(r['input_region_count'], r['region_count'])
```

- **Example output:**

```text
4 4
```

## 8) `dnakit.annotation.region_overlap_indices()`

- **Purpose:** For each region of the first set, the 0-based indices of the overlapping regions of the second set.
- **Method:** `RegionSet.find_overlaps`.
- **API:** `dnakit.annotation.region_overlap_indices(regions_a, regions_b, *, max_rows=100000, verbose=True)`.
- **Input:** Two BED paths, `max_rows`.
- **Output:** `overlaps` (list of index lists), `overlapping_region_count`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_overlap_indices

r = region_overlap_indices('a.bed', 'b.bed', verbose=False)
print(r['overlaps'], r['overlapping_region_count'])
```

- **Example output:**

```text
[[0], [0], [0], [1]] 4
```

## 9) `dnakit.annotation.region_gaps()`

- **Purpose:** Return the gaps between regions on each chromosome, from position 0 up to the chromosome size.
- **Method:** `RegionSet.gaps(chrom_sizes)`: leading, inter-region and trailing gaps, plus one full-length gap for each chromosome of `chrom_sizes` that has no region. Regions on other chromosomes are skipped.
- **API:** `dnakit.annotation.region_gaps(regions, chrom_sizes, *, output=None, max_rows=100000, verbose=True)`.
- **Input:** BED path, `chrom_sizes`, optional `output`, `max_rows`.
- **Output:** `regions`, `region_count`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_gaps

r = region_gaps('a.bed', {'chr1': 1000, 'chr2': 100}, verbose=False)
print([(x['chrom'], x['start'], x['end']) for x in r['regions']])
```

- **Example output:**

```text
[('chr1', 0, 100), ('chr1', 300, 500), ('chr1', 600, 1000), ('chr2', 0, 10), ('chr2', 50, 100)]
```

## 10) `dnakit.annotation.region_disjoin()`

- **Purpose:** Break all regions into non-overlapping pieces: every region boundary becomes a piece boundary.
- **Method:** `RegionSet.disjoin`; the pieces tile the covered positions exactly.
- **API:** `dnakit.annotation.region_disjoin(regions, *, output=None, max_rows=100000, verbose=True)`.
- **Input:** BED path, optional `output`, `max_rows`.
- **Output:** `regions`, `region_count`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_disjoin

r = region_disjoin('a.bed', verbose=False)
print([(x['chrom'], x['start'], x['end']) for x in r['regions']])
```

- **Example output:**

```text
[('chr1', 100, 150), ('chr1', 150, 200), ('chr1', 200, 300), ('chr1', 500, 600), ('chr2', 10, 50)]
```

## 11) `dnakit.annotation.gene_annotation_summary()`

- **Purpose:** Build the gtars gene model, partition list and TSS index from a GTF file (or a saved annotation binary) and report their sizes.
- **Method:** `GenomicDistAnnotation.from_gtf` (or `load_bin`), then `gene_model`, `partition_list(core_prom, prox_prom, chrom_sizes)` and `tss_index`. Defaults follow gtars: protein-coding genes only, Ensembl names converted to UCSC, core promoter 100 bp, proximal promoter 2000 bp.
- **API:** `dnakit.annotation.gene_annotation_summary(*, gtf=None, annotation_binary=None, core_promoter_bp=100, proximal_promoter_bp=2000, protein_coding_only=True, convert_ensembl_ucsc=True, chrom_sizes=None, verbose=True)`.
- **Input:** `gtf` or `annotation_binary`, promoter sizes, `protein_coding_only`, `convert_ensembl_ucsc`, optional `chrom_sizes`.
- **Output:** `n_genes`, `n_exons`, `partition_names`, `tss_count`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import gene_annotation_summary

r = gene_annotation_summary(gtf='g.gtf', verbose=False)
print(r['n_genes'], r['n_exons'], r['tss_count'])
print(r['partition_names'])
```

- **Example output:**

```text
3 6 3
['promoterCore', 'promoterProx', 'exon', 'intron']
```

## 12) `dnakit.annotation.tss_distances()`

- **Purpose:** Distance from each query region to the nearest TSS (or feature) with a gtars `TssIndex`.
- **Method:** The index is built from a BED file of TSS positions (`TssIndex.from_regionset`), from a GTF file (strand-aware TSS) or from a saved annotation binary. `signed=False` calls `calc_tss_distances`, `signed=True` calls `feature_distances` (null for chromosomes without features). The distance is measured from the region midpoint.
- **API:** `dnakit.annotation.tss_distances(regions, *, tss_regions=None, gtf=None, annotation_binary=None, signed=False, protein_coding_only=True, max_rows=100000, verbose=True)`.
- **Input:** Query BED path and exactly one of `tss_regions`, `gtf`, `annotation_binary`; `signed`.
- **Output:** `distances` (region and `distance`), `tss_count`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import tss_distances

for signed in (False, True):
    r = tss_distances('q.bed', gtf='g.gtf', signed=signed, verbose=False)
    print(signed, [x['distance'] for x in r['distances']])
```

- **Example output:**

```text
False [50, 9050]
True [50.0, -9050.0]
```

## 13) `dnakit.annotation.signal_matrix_info()`

- **Purpose:** Load a gtars signal matrix and report its number of regions, conditions and condition names.
- **Method:** `SignalMatrix.from_tsv` (first column holds labels `chrom_start_end`, further columns one condition each) or, with `binary=True`, `SignalMatrix.load_bin`.
- **API:** `dnakit.annotation.signal_matrix_info(path, *, binary=False, verbose=True)`.
- **Input:** Matrix path, `binary`.
- **Output:** `n_regions`, `n_conditions`, `condition_names`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import signal_matrix_info

print(signal_matrix_info('sm.tsv', verbose=False))
```

- **Example output:**

```text
{'status': 'completed', 'operation': 'signal_matrix_info', 'format': 'tsv', 'n_regions': 2, 'n_conditions': 2, 'condition_names': ['c1', 'c2'], 'provenance': {'package': 'gtars 0.10.0', 'convention': 'TSV first column holds labels chrom_start_end'}}
```

## 14) `dnakit.annotation.region_tokenizer_vocabulary()`

- **Purpose:** Build a gtars region tokenizer and return its vocabulary and special tokens.
- **Method:** `source_type='bed'` uses `Tokenizer.from_bed` (a BED universe), `'config'` uses `from_config` (a TOML file with `universe` and `tokenizer_type = "ailist"` or `"bits"`), `'pretrained'` uses `from_pretrained` (a Hugging Face model name, network). Tokens look like `chr1:0-1000`.
- **API:** `dnakit.annotation.region_tokenizer_vocabulary(source, *, source_type='bed', include_vocabulary=True, verbose=True)`.
- **Input:** `source`, `source_type`, `include_vocabulary`.
- **Output:** `vocab_size`, `special_tokens`, `vocabulary`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_tokenizer_vocabulary

r = region_tokenizer_vocabulary('u.bed', verbose=False)
print(r['vocab_size'])
print(r['vocabulary'])
print(r['special_tokens']['cls_token'])
c = region_tokenizer_vocabulary('c.toml', source_type='config', include_vocabulary=False, verbose=False)
print(c['vocab_size'])
```

- **Example output:**

```text
10
{'chr1:0-1000': 0, 'chr1:1000-2000': 1, 'chr2:0-500': 2, '<unk>': 3, '<pad>': 4, '<mask>': 5, '<cls>': 6, '<eos>': 7, '<bos>': 8, '<sep>': 9}
<cls>
10
```

## 15) `dnakit.annotation.region_token_ids()`

- **Purpose:** Convert between region tokens and integer ids with a gtars tokenizer.
- **Method:** `direction='to_ids'` uses `convert_tokens_to_ids` (tokens that are not in the vocabulary map to the unknown id); `direction='to_tokens'` uses `convert_ids_to_tokens`.
- **API:** `dnakit.annotation.region_token_ids(source, values, *, source_type='bed', direction='to_ids', verbose=True)`.
- **Input:** `source`, `values`, `source_type`, `direction`.
- **Output:** `outputs`, one per input.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_token_ids

print(region_token_ids('u.bed', ['chr1:1000-2000', 'unknown'], verbose=False)['outputs'])
print(region_token_ids('u.bed', [0, 6], direction='to_tokens', verbose=False)['outputs'])
```

- **Example output:**

```text
[1, 3]
['chr1:0-1000', '<cls>']
```

## 16) `dnakit.annotation.region_special_tokens_mask()`

- **Purpose:** Flag which of the given tokens are special tokens such as `<cls>` or `<pad>`.
- **Method:** `Tokenizer.get_special_tokens_mask` on a list of token strings.
- **API:** `dnakit.annotation.region_special_tokens_mask(source, tokens, *, source_type='bed', verbose=True)`.
- **Input:** `source`, `tokens`, `source_type`.
- **Output:** `mask` (true = special token).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_special_tokens_mask

print(region_special_tokens_mask('u.bed', ['<cls>', 'chr1:0-1000', '<pad>'], verbose=False)['mask'])
```

- **Example output:**

```text
[True, False, True]
```

## 17) `dnakit.annotation.region_database_list_sets()`

- **Purpose:** Build a LOLA region database from BED files (or read a LOLA folder) and list its region sets.
- **Method:** `RegionDB.from_bed_files(files, names)` or `RegionDB.from_folder(folder)`, then `list_region_sets`. Names default to the file names.
- **API:** `dnakit.annotation.region_database_list_sets(*, bed_files=None, folder=None, names=None, verbose=True)`.
- **Input:** `bed_files` or `folder`, optional `names`.
- **Output:** `set_count`, `sets` (index, name).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_database_list_sets

r = region_database_list_sets(bed_files=['a.bed', 'b.bed'], verbose=False)
print(r['set_count'], r['sets'])
```

- **Example output:**

```text
2 [{'index': 0, 'name': 'a.bed'}, {'index': 1, 'name': 'b.bed'}]
```

## 18) `dnakit.annotation.region_database_get_sets()`

- **Purpose:** Return the regions of selected region sets (by 0-based index, default all) of a LOLA region database.
- **Method:** `RegionDB.get_region_sets(indices)`.
- **API:** `dnakit.annotation.region_database_get_sets(*, bed_files=None, folder=None, names=None, indices=None, max_rows=100000, verbose=True)`.
- **Input:** `bed_files` or `folder`, `names`, `indices`, `max_rows`.
- **Output:** `sets` (index, name, region_count, regions).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_database_get_sets

r = region_database_get_sets(bed_files=['a.bed', 'b.bed'], indices=[1], verbose=False)
s = r['sets'][0]
print(s['name'], s['region_count'], [(x['chrom'], x['start'], x['end']) for x in s['regions']])
```

- **Example output:**

```text
b.bed 3 [('chr1', 180, 520), ('chr2', 40, 90), ('chr3', 5, 15)]
```

## Rows already covered

- `gtars.refget.digest_sequence` (sequence_record_from_raw_sequence): `sequence_digest_record()`.
- `gtars.genomic_distributions.calc_gc_content` and `GenomeAssembly.__init__`: `region_gc_content()` (FASTA genome).
- `gtars.genomic_distributions.calc_expected_partitions`: `partition_region_counts()` with `chrom_sizes`.

Not implemented: `BinaryGenomeAssembly` (.fab files are made only by the gtars command line `gtars prep`, not by the Python API), and `GenomicDistAnnotation.load_bin` input files likewise come from the gtars command line; the functions accept such a file when the user has one (`annotation_binary`, `binary=True`), but the tests do not build one.
