# Ensembl VEP: variant consequences, filtering, haplotypes and cache tools

This group of functions in `dnakit.variants` is a set of typed adapters over the programs of [Ensembl VEP](https://github.com/Ensembl/ensembl-vep) (tested with release 116.3): `vep`, `filter_vep`, `haplo`, `convert_cache.pl` and the input-format validator script. `vep_annotate_variants()` exposes the annotation options of `vep` as parameters (the output formats VCF, tab and JSON can also be produced by the older `annotate_vcf_ensembl_vep()`). `variant_recoder` is not wrapped: it needs a live connection to the Ensembl MySQL server, which could not be checked. VEP is licensed under Apache-2.0; DNAKit runs it as separate processes and copies none of its code.

The functions are in the module `dnakit.variants.ensembl_vep_p09_cli_functions`. Install VEP with `conda install -c bioconda ensembl-vep htslib` (or from the GitHub repository with `perl INSTALL.pl`) and put its programs on PATH, or pass `vep_path`, `filter_vep_path`, `haplo_path`, `convert_cache_path`. A missing program raises `BackendUnavailableError` with an install hint before any output directory is created. Caches, FASTA files and plugins are the caller's responsibility; nothing is downloaded. Every function writes into a new `output_directory`, validates its inputs (`ConfigurationError`) and raises `BackendExecutionError` with the log tail when VEP fails. Optional parameters are not passed when `None` or `False`, so VEP's defaults apply.

## Example data

All examples read these small files from the working directory: `genome.fa` is a 3000-base contig, `genes.gff3` has one two-exon coding gene on it, `variants.vcf` has four SNVs (two in the CDS, one in the intron, one downstream), `phased.vcf` has three phased SNVs of one sample, and `input.txt` is a three-line file in VEP's default input format with two mistakes. The examples turn the GFF3 into a bgzipped, tabix-indexed file with pysam, and assume that `vep`, `filter_vep` and `haplo` are on PATH.

**Example file:** `genome.fa`

```text
>chr1
GTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTC
AATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGG
TGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGC
TTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGG
TTAACGGGGATGGCAAGTACATTTTTTCGTAGATGTGCCTTGCTAACGAAAGTATTAAAC
ACGTCCCTCACAATAGAATCATAGTTGGACGCGCGACGGCCGTTCCAGAAAATCTTTGAA
TACTCAATCCTGCGGGTTCGGTGACCTAAAACCCATTGATTGTGTTACCCAGTTCGAGCG
CATAGGGAATTCAGGTCCACACATGGCTGGATCCCCATGATATTCAAGAACTATACATTA
AGTTGAACCTCCAGAACACAATGGCTCAAGTTCCGTCAGTATTTTCTGTTGGCTTGGACA
ACTCTCTCTATTTTATTGGCAATTTAAGGACTCCTCTGTGGATTCGAGAGACGGCAATTA
AAGTGATCCGGCACAGTTTCTGGTATCGAAGGTACTCCACAAGAGCTGCTGTGAGTTTTA
CTGACGGAAAGTGAGCAAAGGCTAACGTTATTCCGTGAGCACGGGACATCCATTCTTCGT
GAGCTACAGCTCGAGAATCAGCTTCTAACCAAGCGATGCAGAACCGGCTACTTTAAGCAT
TGATGAATGCGTCGTAAGTGATACTCGACGATTCTCATGCAACGAAGTTAACCTATAGTA
ACTTACATTTTACGCGCTAGCTTCGCTGGAACTAATATCCATGTCTCAGAACTAGCGGCC
GAGAATGGGTTCCGAATCCTAAACTCCGACATGAGTTAAGGTTGCATACTAGGTCTGATA
CTAAAAGCGGGGTCAGGAGTCCGTCCAGAATATAATATAGCTAAGAAGAAAATGGAAAAT
CTATATCAGGTGGAGGCAGCCAATAGAGTGTGGTGGAATCAGTTTTGTAATATTACGACA
TAACCCTCGGAGGTTGAGCCAAACACCCCATTCGGCGGCCAAAGCGGGTGTCCCCACACA
TCCGTCGTAATGTCTAGGGAACAGTTGGAGCTGTGCGCAATCGTGTTGGACTATTGACAT
ACCCCTATTCGACCAGTGCAGTAGACTATACCACTTTTGAATCATCAGCAGACCTAATAG
CTTCATCCCTTCTAGTCGACTTTCCGGACCGATGCGCACTAATGATCGAAGTGTGTCTTT
ACTGAATCAGAAGTCGGAGAAAATTCTGCTGTACGAGATGTACGTAGCGATGATGATACG
GGTGATTCTTACATGAGGTACATCGAAAAAGGGTCATTGCGTTTACGTGAATGGGATTTG
CCTGGCCTATGGCCTTAGTACCTCTAAGAGGGCAACTAGTCACGGCGTAAACAGACAAGG
GTCGGATCCTAGTCACTAGTAGATCAACGGCTAAGTGGGCGTCACCAGACCCTCGCCCAT
CTGGACTAGTAGACCGTGTGCCTCAATCGACGCGTGAGGACACGTTGATTTCTTATACGG
CGTCTCAGCTTTCTTTTCCGCCAATTACATCAAATTCAAGTGCCTTGAGATACCGGGCCT
ATTCGTGCGTAGACCGGCGCGATGAGACAATGCCGGCATTTAAGTAAAGGTCCGGGATTA
TTGTAGACTAACTCTGCCAAATTATGCGTGCCTTGAGAGTATACACTCCGCATCCTTCTC
GTATGATGGCAGGGAGGGCCTTGGTTCCAGAACTTATCGCTAATGCGCTGGGTGCCATGC
ATGGAGCGTCTCTTGACATCAGGAAGTACCGGCTCCTGTTGGGGAGTTCGAAAATAGCTT
ATGGGACATCAGGCCATACTACTAGAACCTTGACCTGTACTGCATATCGTTTTCGGTTAA
ACTCGGTAGGTAAAAACACTTGTCGCCAGGATTCGCCACTCAAGGTTTATAGCAAAGAAT
ACTACTGAAACTGGCGTCCCCCATTAGTCACATCCTGCGAGGCGGCTCTCGACGGCATAA
CGGGATCTCTCCGGCTAAGGTCGAATGAGTTTAGTCTGAAGGTGGCAGCGCACACGACGG
TATGTTTCACCGCCGTCTGTTCGGAGCTATAGCATCATTAGGCCTAGCTATGGCGCCCCC
GTTTCCTACCGCAGCTACTGCATGCACTCGTTACGACAAAGTTCCTTTATTGGCTACAGG
ATTCTATCGTTCCGGCCTAGTGCGGGGGCATCTGCAGGCCGTCGTCCGTCCTACGTAAGG
CTCATGCGTACGGTAACAAAAGTTGGTTAGACCTGACCACGCGTGATCCGCGCCTCCTTC
CCCGCGTCACGTTAGGTGAAACCCTGCCCCGCTTAGCACGTCTCATGCTAACCCTCTTTC
CTGAGGGTAAGGGAACGTGCAATCCCGAAAAAGGTGAATGAAAGGCTGGGACGCTCCTGG
TTTGGGGCGCTCTAGGTGTTGTGGCTTGAACGTTAGCACCTCGGCTTCCGAATTGATACT
AGGCATCCTATCGAACACTATTAGCCGCCCTTCTTCATATTACTGTGGGTACAGACAAAT
AGAGGAGGTGTACGGTGACTAGCGGGTATCTATCATTAGGAAAGGGTACCGGCATCAGGA
AGAAAGCACCGTCCAATGATCACCACGGCTCCTTGCTTGCCACCGAAGCCATAGACATAT
GAATAGCCCGCTGTGGCTCTCCGACGAAGTGCGCGGCAGATCAGTTTTCTCCCCCTGAGA
TCACCAGACGACGGAACAGAACTGCGAGGCTGATCCACTGTTTTGCGCCTCTGAAATAGA
TTAAAGGGAATATCCCAGTGCAGGCCGAGGGGGTCTGATAATATACAGTCTAAAGAGTTA
GAATATCATCGTAACAGCCCATAGAACAATTCCCGTATTTAAGCATAGCATACGGATAAT
```

**Example file:** `genes.gff3`

```text
##gff-version 3
chr1	test	gene	501	1150	.	+	.	ID=gene1;Name=GENE1;biotype=protein_coding
chr1	test	mRNA	501	1150	.	+	.	ID=tx1;Parent=gene1;Name=GENE1-201;biotype=protein_coding
chr1	test	exon	501	650	.	+	.	ID=exon1;Parent=tx1
chr1	test	CDS	501	650	.	+	0	ID=cds1;Parent=tx1
chr1	test	exon	1001	1150	.	+	.	ID=exon2;Parent=tx1
chr1	test	CDS	1001	1150	.	+	0	ID=cds2;Parent=tx1
```

**Example file:** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=3000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	521	.	A	C	.	.	.
chr1	530	.	T	A	.	.	.
chr1	700	.	C	A	.	.	.
chr1	2000	.	T	A	.	.	.
```

**Example file:** `phased.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=3000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	522	.	T	A	.	.	.	GT	0|1
chr1	531	.	G	A	.	.	.	GT	1|1
chr1	540	.	A	C	.	.	.	GT	0|1
```

**Example file:** `input.txt`

```text
5 140532 140532 T/C +
5 140540 140540 A/G f
3 100 100 T +
```

## 1) `dnakit.variants.ensembl_vep_p09_cli_functions.vep_annotate_variants()`

- **Purpose:** Predict variant consequences with Ensembl VEP, with its annotation options as typed parameters.
- **Method:** Run `vep` with a fixed argument list. One function covers `summary`, `variant_class`, `hgvs`/`hgvsg`/`spdi`, 3-prime shifting (`shift_hgvs`, `shift_3prime`, `shift_genomic`), `check_ref`/`lookup_ref`, per-sample genotypes (`individual`, `individual_zyg`), co-located variants (`check_existing`) with `af`, `af_1kg`, `af_gnomade`, `af_gnomadg`, `max_af`, frequency filtering (`check_frequency`, `freq_pop`, `freq_freq`, `freq_gt_lt`, `freq_filter`, `filter_common`), ClinVar (`clin_sig_allele`, `clinvar_somatic_classification`), `gene_phenotype`, `mirna`, `regulatory` (`cell_type`), UCSC conservation (`phylop`, `phastcons`, `ucsc_assembly`, `ucsc_data_root`), user tracks (`custom_annotations`, `gff`, `gtf`, `bigwig`), transcript correction from a `bam`, `domains`, `numbers`, `nearest` and a `stats` report. `output_format` is `vcf`, `tab` or `json`. Transcripts come from a cache (`cache_directory`) or from a bgzipped, tabix-indexed `gff`/`gtf` with the reference `fasta`. Options that need a cache, a database or the network (frequencies, ClinVar, UCSC scores, the Condel flag is not offered because current VEP ignores it) are passed through unchanged and give data only when that resource exists; their command-line mapping is tested, but they were not run against a real cache here. VEP writes an index next to the FASTA file.
- **API:** `dnakit.variants.ensembl_vep_p09_cli_functions.vep_annotate_variants(input_path: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['vcf', 'tab', 'json']" = 'vcf', input_format: "Literal['ensembl', 'vcf', 'hgvs', 'id', 'region', 'vep', 'guess'] | None" = None, species: 'str' = 'homo_sapiens', assembly: 'str | None' = None, offline: 'bool | None' = None, cache_directory: 'str | Path | None' = None, cache_version: 'int | None' = None, fasta: 'str | Path | None' = None, gff: 'str | Path | None' = None, gtf: 'str | Path | None' = None, bigwig: 'str | Path | None' = None, bam: 'str | Path | None' = None, custom_annotations: 'Sequence[Mapping[str, object]] | None' = None, summary: 'bool' = False, variant_class: 'bool' = False, hgvs: 'bool' = False, hgvsg: 'bool' = False, spdi: 'bool' = False, shift_hgvs: 'int | None' = None, shift_3prime: 'int | None' = None, shift_genomic: 'int | None' = None, check_ref: 'bool' = False, lookup_ref: 'bool' = False, individual: 'str | Sequence[str] | None' = None, individual_zyg: 'str | Sequence[str] | None' = None, check_existing: 'bool' = False, af: 'bool' = False, af_1kg: 'bool' = False, af_gnomade: 'bool' = False, af_gnomadg: 'bool' = False, max_af: 'bool' = False, check_frequency: 'bool' = False, freq_pop: 'str | None' = None, freq_freq: 'float | None' = None, freq_gt_lt: "Literal['gt', 'lt'] | None" = None, freq_filter: "Literal['exclude', 'include'] | None" = None, filter_common: 'bool' = False, clin_sig_allele: 'int | None' = None, clinvar_somatic_classification: 'bool' = False, gene_phenotype: 'bool' = False, mirna: 'bool' = False, regulatory: 'bool' = False, cell_type: 'str | None' = None, phylop: 'Sequence[str] | None' = None, phastcons: 'Sequence[str] | None' = None, ucsc_assembly: 'str | None' = None, ucsc_data_root: 'str | None' = None, domains: 'bool' = False, numbers: 'bool' = False, nearest: "Literal['transcript', 'gene', 'symbol'] | None" = None, stats: "Literal['html', 'text'] | None" = None, fork: 'int | None' = None, vep_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A variant file (VCF by default), a new output directory and the options. Without `cache_directory` and without `gff`/`gtf`, VEP runs offline and needs a cache in its default folder; `offline=False` reads the Ensembl database.
- **Output:** `outputs` (`output` = `annotations.vcf`/`.tsv`/`.json`, optional `stats`, `warnings`), `columns`, `records` (first rows; VCF rows carry the parsed `CSQ` entries) and `record_count`, plus `version` and `command`.
- **Example code:**

```python
import tempfile
import pysam

pysam.tabix_compress("genes.gff3", "genes.gff3.gz", force=True)
pysam.tabix_index("genes.gff3.gz", preset="gff", force=True)
from dnakit.variants.ensembl_vep_p09_cli_functions import vep_annotate_variants

work = tempfile.mkdtemp()
result = vep_annotate_variants("variants.vcf", work + "/vcf", gff="genes.gff3.gz", fasta="genome.fa", verbose=False)
for row in result["records"]:
    print(row["pos"], row["annotations"][0]["Consequence"])
table = vep_annotate_variants("variants.vcf", work + "/tab", output_format="tab", gff="genes.gff3.gz", fasta="genome.fa", variant_class=True, hgvs=True, hgvsg=True, spdi=True, numbers=True, verbose=False)
first = table["records"][0]
print(first["VARIANT_CLASS"], first["EXON"], first["HGVSc"], first["HGVSg"], first["SPDI"])
```

- **Example output:**

```text
521 synonymous_variant
530 synonymous_variant
700 intron_variant
2000 downstream_gene_variant
SNV 1/2 tx1:c.21A>C chr1:g.521A>C chr1:520:A:C
```

## 2) `dnakit.variants.ensembl_vep_p09_cli_functions.vep_filter_results()`

- **Purpose:** Filter VEP-annotated variant records by expression.
- **Method:** Run `filter_vep`. `filters` are expressions such as `Consequence is missense_variant and AF < 0.01`; several expressions are alternatives. `list_fields` runs `--list`, `only_matched` rewrites the VCF CSQ field with matching entries, `soft_filter` flags instead of dropping, `count_only` returns the number of matches, `input_format`, `vcf_info_field`, `start`, `limit` and `gzipped` map to `--format`, `--vcf_info_field`, `--start`, `--limit` and `--gz`. The ontology option (needs a database) is not offered.
- **API:** `dnakit.variants.ensembl_vep_p09_cli_functions.vep_filter_results(input_path: 'str | Path', output_directory: 'str | Path', *, filters: 'Sequence[str] | None' = None, list_fields: 'bool' = False, only_matched: 'bool' = False, soft_filter: 'bool' = False, count_only: 'bool' = False, input_format: "Literal['vcf', 'tab'] | None" = None, vcf_info_field: 'str | None' = None, start: 'int | None' = None, limit: 'int | None' = None, gzipped: 'bool' = False, filter_vep_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A VEP result file (VCF or tab), a new output directory and expressions (letters, digits and `_ .,:<>=!()"'+*/%-` only).
- **Output:** `outputs`, `record_count` and `records` (data lines split into columns), or `fields` with `list_fields`, or `count` with `count_only`.
- **Example code:**

```python
import tempfile
import pysam

pysam.tabix_compress("genes.gff3", "genes.gff3.gz", force=True)
pysam.tabix_index("genes.gff3.gz", preset="gff", force=True)
from dnakit.variants.ensembl_vep_p09_cli_functions import vep_annotate_variants, vep_filter_results

work = tempfile.mkdtemp()
annotated = vep_annotate_variants("variants.vcf", work + "/vcf", gff="genes.gff3.gz", fasta="genome.fa", verbose=False)
path = annotated["outputs"]["output"]
print("IMPACT" in vep_filter_results(path, work + "/fields", list_fields=True, verbose=False)["fields"])
kept = vep_filter_results(path, work + "/kept", filters=["Consequence is intron_variant"], verbose=False)
print(kept["record_count"], kept["records"][0][1])
print(vep_filter_results(path, work + "/count", filters=["IMPACT is LOW"], count_only=True, verbose=False)["count"])
```

- **Example output:**

```text
True
1 700
2
```

## 3) `dnakit.variants.ensembl_vep_p09_cli_functions.vep_predict_haplotypes()`

- **Purpose:** Predict the consequences of whole transcript haplotypes in phased genotypes.
- **Method:** Run `haplo --json`: it rebuilds the transcript sequence of each sample haplotype and reports the amino-acid changes of the whole haplotype. Transcripts come from `cache_directory` or from `gff`/`gtf` with `fasta`; `haplotype_frequencies`, `phased` and `chromosomes` map to `--haplotype_frequencies`, `--phased` and `--chr`.
- **API:** `dnakit.variants.ensembl_vep_p09_cli_functions.vep_predict_haplotypes(input_path: 'str | Path', output_directory: 'str | Path', *, input_format: "Literal['vcf', 'ensembl', 'region', 'vep', 'guess'] | None" = None, species: 'str' = 'homo_sapiens', assembly: 'str | None' = None, offline: 'bool | None' = None, cache_directory: 'str | Path | None' = None, cache_version: 'int | None' = None, fasta: 'str | Path | None' = None, gff: 'str | Path | None' = None, gtf: 'str | Path | None' = None, haplotype_frequencies: 'str | Path | None' = None, phased: 'bool' = False, chromosomes: 'str | None' = None, haplo_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A VCF with genotypes and a new output directory.
- **Output:** `outputs` (`output` = `haplotypes.json`), `records` (one object per transcript with `protein_haplotypes` and `cds_haplotypes`) and `record_count`.
- **Example code:**

```python
import tempfile
import pysam

pysam.tabix_compress("genes.gff3", "genes.gff3.gz", force=True)
pysam.tabix_index("genes.gff3.gz", preset="gff", force=True)
from dnakit.variants.ensembl_vep_p09_cli_functions import vep_predict_haplotypes

result = vep_predict_haplotypes("phased.vcf", tempfile.mkdtemp() + "/run", gff="genes.gff3.gz", fasta="genome.fa", verbose=False)
transcript = result["records"][0]
print(result["record_count"], transcript["transcript_id"])
for haplotype in transcript["protein_haplotypes"]:
    print(haplotype["name"], haplotype["frequency"])
```

- **Example output:**

```text
1 tx1
tx1.p:8F>I,11G>S,14N>H 0.5
tx1.p:11G>S 0.5
```

## 4) `dnakit.variants.ensembl_vep_p09_cli_functions.vep_validate_input_format()`

- **Purpose:** Check a file in VEP's default input format.
- **Method:** Run Ensembl's `vep_input_format_validator.py` with the current Python interpreter. It checks the column count, the allele characters and `REF/ALT` form, the strand, the coordinates against the allele length, insertions, and the order of variants. `validator_path` is the script from the `validator/` folder of the ensembl-vep repository (it is not part of the conda package). Nothing is written.
- **API:** `dnakit.variants.ensembl_vep_p09_cli_functions.vep_validate_input_format(input_path: 'str | Path', *, validator_path: 'str | Path | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** The input file and `validator_path`.
- **Output:** `valid`, `error_count` and `errors` (line number, message, line content).
- **Example code:**

```python
import os
from dnakit.variants.ensembl_vep_p09_cli_functions import vep_validate_input_format

script = os.environ.get("VEP_INPUT_VALIDATOR", "vep_input_format_validator.py")
result = vep_validate_input_format("input.txt", validator_path=script, verbose=False)
print(result["valid"], result["error_count"])
for error in result["errors"]:
    print(error["line"], error["message"])
```

- **Example output:**

```text
False 2
2 Strand is not in the expected format (+ or -)
3 Allele 'T' is not in the supported format 'REF/ALT' (e.g. A/C), except for structural variants
```

## 5) `dnakit.variants.ensembl_vep_p09_cli_functions.vep_convert_cache()`

- **Purpose:** Convert a downloaded VEP cache so that it is read through tabix.
- **Method:** Run `convert_cache.pl` (`vep_convert_cache` in the conda package) on `<cache_directory>/<species>/<version>`; `bgzip` and `tabix` must be on PATH. It modifies the cache in place; `remove_original` is `--remove`, `use_sereal` is `--sereal`, `compress_command` is `--compress` and `fork` is `--fork`. Only a run on an empty cache folder was checked here, because real caches are several gigabytes.
- **API:** `dnakit.variants.ensembl_vep_p09_cli_functions.vep_convert_cache(cache_directory: 'str | Path', *, species: 'str', version: 'str', compress_command: 'str | None' = None, use_sereal: 'bool' = False, remove_original: 'bool' = False, fork: 'int | None' = None, convert_cache_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** The folder that holds `<species>/<version>`, `species` and `version` (names, comma lists or `all`).
- **Output:** `status`, `command`, `cache_directory` and the program's `log_tail`.
- **Example code:**

```python
import os
import tempfile
from dnakit.variants.ensembl_vep_p09_cli_functions import vep_convert_cache

cache = tempfile.mkdtemp()
os.makedirs(cache + "/homo_sapiens/116_GRCh38/1")
result = vep_convert_cache(cache, species="homo_sapiens", version="116_GRCh38", verbose=False)
print(result["status"], "All done" in result["log_tail"])
```

- **Example output:**

```text
completed True
```
