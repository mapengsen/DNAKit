# GATK VCF utilities: header, dictionary, liftover, split and convert

These functions in `dnakit.variants.gatk_p02_cli_functions` run GATK (tested with 4.7.0.0) command-line tools through a fixed argument whitelist, without a shell. Install GATK and Java 17 (`conda install -c bioconda gatk4` or the release zip), put `gatk` on PATH or pass `gatk_path`. Outputs go to new paths (nothing is overwritten); GATK index files are reported in `index_files`. Each result is a dictionary with `status` (`completed` only when the exit code is 0 and the outputs exist), `version`, `command`, `outputs`, `log_tail` and tool-specific values. Unset arguments keep GATK's defaults. A missing `gatk` raises `BackendUnavailableError` with an install hint.

## Example data

The examples run in one folder, in order; later examples reuse files written by earlier ones. `gatk` must be on PATH.

**Example file:** `calls.vcf`

```text
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=QD,Number=1,Type=Float,Description="Qual by depth">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##contig=<ID=chr1,length=600>
##contig=<ID=chr2,length=300>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	50	.	G	A	50	.	QD=1.5;DP=40	GT:DP	0/1:20	0/0:20
chr1	120	.	T	A	50	.	QD=18.0;DP=40	GT:DP	1/1:20	0/1:20
chr1	200	.	CTA	C	50	.	QD=12.0;DP=40	GT:DP	0/1:20	0/1:20
chr1	310	.	A	C	50	.	QD=25.0;DP=40	GT:DP	0/1:20	1/1:20
chr1	400	.	C	CAC	50	.	QD=9.0;DP=40	GT:DP	0/1:20	0/0:20
chr2	30	.	A	C	50	.	QD=21.0;DP=40	GT:DP	1/1:20	0/0:20
chr2	150	.	A	C	50	.	QD=3.0;DP=40	GT:DP	0/1:20	0/1:20
```

**Example file:** `ref.dict`

```text
@HD	VN:1.6
@SQ	SN:chr1	LN:600
@SQ	SN:chr2	LN:300
```

**Example file:** `shift.chain`

```text
chain 1000 chr1 600 + 5 600 chr1 595 + 0 595 1
595

chain 900 chr2 300 + 0 300 chr2 300 + 0 300 2
300
```

**Example file:** `new.fa`

```text
>chr1
TTATTTCTTGCACATAATAATAAGAATCTAACTTAATAACTTGTGAACCCTATTTAGCCTGGATTTGGTC
TTTTTGAATATCATATTATGTCCCATTTACACGAAATCTGTTGATAATCAGTCTTTCGATCTCCGAGGTT
TCGATTTCCGCTATTTATCCTGATTTCCTGACCTGATTGTCGTAATTTTGTAGTCTAAAACAAGCAAGAT
ATTAAATCGTCAGCATCATCAGTGCAATAGTACATCATTATAACGTATTTAAAGTTTTAACTGTTGCAAC
TCATTTCTTTTACAAAATCGGATAACGTATTGTGTCTAAGTATCTGGTTATACAATCGTTGATGCTAATT
GCTTAAAAAAGGTTTAAATTGGACCACAGCCGGATATCGGCAAGCTGGCTTTCCCGTTCTATCCTACTTA
CTAGTATTGGTCCAGACTCTCTGTAAGCATCCGCACTTCTTTGATAACAATATATTTCTTTTTAACGGAA
ACGTAGGATTTAATCAATTGTTATTAGAGGACCCTCTACTTAAGTTATCATGACAAAAGCCCTGCCTGTA
ATCTTTATATATATGTGAATACGCTCGTATATTAT
>chr2
GTGCTGGCTTACTTGCGCATTGAGGGACAACGGTAGTACCCCCGGGGATCGTGGCATTACGGGCCGTTAG
GCTTCCTTGCTGTCATATACCGTGCTGACCTTCCAATCAGGTACAGATCTTCTAGAGGCTTTCATCACGG
CGACAACGGAGCTGAGGGGTAGCGGGCATTAACCGCGCGTCCACGGTAAGTTACAGGTTATGCGAGGATT
ATCGCTGAGTAGATTCAGAGACCCCGAAGTGAACTGAGGTCAGATAGGACTGACGCCGTATCTTGAACAA
ATGCTGGCGTCCGATTGACA
```

**Example file:** `new.fa.fai`

```text
chr1	595	6	70	71
chr2	300	616	70	71
```

**Example file:** `new.dict`

```text
@HD	VN:1.6
@SQ	SN:chr1	LN:595
@SQ	SN:chr2	LN:300
```

## 1) `fix_vcf_header_gatk()`

- **Purpose:** Replace a VCF header, or add the INFO/FORMAT/FILTER definitions the records use but the header lacks.
- **Method:** Runs `gatk FixVcfHeader` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.fix_vcf_header_gatk(vcf: 'str', output_path: 'str', *, header: 'str | None' = None, enforce_same_samples: 'bool' = True, check_first_n_records: 'int | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF and optionally a VCF whose header replaces it.
- **Output:** The fixed VCF and a `summary`.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.fix_vcf_header_gatk("calls.vcf", "fixed.vcf")
print(r["status"], r["summary"]["records"], r["summary"]["samples"])
```

- **Example output:**

```text
completed 7 ['S1', 'S2']
```

## 2) `update_vcf_sequence_dictionary_gatk()`

- **Purpose:** Update the contig lines of a VCF from a dictionary, VCF, BAM or reference.
- **Method:** Runs `gatk UpdateVCFSequenceDictionary` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.update_vcf_sequence_dictionary_gatk(vcf: 'str', output_path: 'str', *, source_dictionary: 'str | None' = None, reference: 'str | None' = None, replace: 'bool' = False, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF and a dictionary source (or a reference).
- **Output:** The updated VCF and a `summary`.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.update_vcf_sequence_dictionary_gatk("calls.vcf", "updated.vcf", source_dictionary="ref.dict", replace=True)
print(r["status"], [x for x in Path("updated.vcf").read_text().splitlines() if x.startswith("##contig")])
```

- **Example output:**

```text
completed ['##contig=<ID=chr1,length=600>', '##contig=<ID=chr2,length=300>']
```

## 3) `liftover_vcf_gatk()`

- **Purpose:** Lift a VCF over to another reference build with a chain file.
- **Method:** Runs `gatk LiftoverVcf` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.liftover_vcf_gatk(vcf: 'str', chain_file: 'str', reference: 'str', output_path: 'str', reject_path: 'str', *, min_match: 'float | None' = None, recover_swapped_ref_alt: 'bool' = False, write_original_position: 'bool' = False, write_original_alleles: 'bool' = False, tags_to_drop: 'list[str] | None' = None, tags_to_reverse: 'list[str] | None' = None, allow_missing_fields_in_header: 'bool' = False, disable_sort: 'bool' = False, warn_on_missing_contig: 'bool' = False, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF, a chain file and the target reference FASTA (with .fai and .dict).
- **Output:** Lifted and rejected VCFs with their summaries.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.liftover_vcf_gatk("calls.vcf", "shift.chain", "new.fa", "lifted.vcf", "rejected.vcf")
rows = [x.split("\t")[:2] for x in Path("lifted.vcf").read_text().splitlines() if not x.startswith("#")]
print(r["status"], r["summary"]["records"], r["reject_summary"]["records"], rows[:3])
```

- **Example output:**

```text
completed 7 0 [['chr1', '45'], ['chr1', '115'], ['chr1', '195']]
```

## 4) `make_sites_only_vcf_gatk()`

- **Purpose:** Remove genotype columns but keep site-level information.
- **Method:** Runs `gatk MakeSitesOnlyVcf` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.make_sites_only_vcf_gatk(vcf: 'str', output_path: 'str', *, samples: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF/BCF and optionally the samples to use.
- **Output:** The sites-only VCF and a `summary` (empty sample list).
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.make_sites_only_vcf_gatk("calls.vcf", "sites.vcf")
print(r["status"], r["summary"]["records"], r["summary"]["samples"])
```

- **Example output:**

```text
completed 7 []
```

## 5) `split_vcf_snps_indels_gatk()`

- **Purpose:** Split a VCF into SNP and indel files with identical headers.
- **Method:** Runs `gatk SplitVcfs` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.split_vcf_snps_indels_gatk(vcf: 'str', snp_output_path: 'str', indel_output_path: 'str', *, strict: 'bool' = True, sequence_dictionary: 'str | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF/BCF and whether other types are an error (`strict`).
- **Output:** SNP and indel VCFs with their summaries.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.split_vcf_snps_indels_gatk("calls.vcf", "snps.vcf", "indels.vcf")
print(r["status"], r["snp_summary"]["records"], r["indel_summary"]["records"])
```

- **Example output:**

```text
completed 5 2
```

## 6) `convert_vcf_format_gatk()`

- **Purpose:** Convert VCF to BCF or BCF to VCF (formats follow the file extensions).
- **Method:** Runs `gatk VcfFormatConverter` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.convert_vcf_format_gatk(input_path: 'str', output_path: 'str', *, require_index: 'bool' = True, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** An input and an output path with .vcf, .vcf.gz or .bcf extensions.
- **Output:** The converted file (a `summary` for VCF outputs).
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

a = g.convert_vcf_format_gatk("calls.vcf", "calls.bcf", require_index=False)
b = g.convert_vcf_format_gatk("calls.bcf", "back.vcf", require_index=False)
print(a["status"], b["status"], b["summary"]["records"])
```

- **Example output:**

```text
completed completed 7
```

## 7) `vcf_to_interval_list_gatk()`

- **Purpose:** Create a Picard interval list from the records of a VCF.
- **Method:** Runs `gatk VcfToIntervalList` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.vcf_to_interval_list_gatk(vcf: 'str', output_path: 'str', *, include_filtered: 'bool' = False, variant_id_method: 'str | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF/BCF, whether to include filtered records and the interval naming method.
- **Output:** The interval list file and `intervals` rows `[contig, start, end, strand, name]`.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.vcf_to_interval_list_gatk("calls.vcf", "calls.interval_list")
print(r["status"], len(r["intervals"]), r["intervals"][:2])
```

- **Example output:**

```text
completed 7 [['chr1', 50, 50, '+', 'interval-1'], ['chr1', 120, 120, '+', 'interval-2']]
```

## 8) `print_variants_spark()`

- **Purpose:** Rewrite the variants of a VCF as one file through the GATK Spark engine (sites-only input in 4.7.0.0).
- **Method:** Runs `gatk PrintVariantsSpark` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.print_variants_spark(vcf: 'str', output_path: 'str', *, reference: 'str | None' = None, variant_shard_size: 'int | None' = None, use_nio: 'bool' = False, spark_master: 'str | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A sites-only VCF, optional shard size and Spark master.
- **Output:** The printed VCF and a `summary`.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.print_variants_spark("sites.vcf", "printed.vcf", spark_master="local[1]")
print(r["status"], r["summary"]["records"])
```

- **Example output:**

```text
completed 7
```
