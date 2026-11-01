# GATK VCF validation, filtering, annotation and recalibration

These functions in `dnakit.variants.gatk_p02_cli_functions` run GATK (tested with 4.7.0.0) command-line tools through a fixed argument whitelist, without a shell. Install GATK and Java 17 (`conda install -c bioconda gatk4` or the release zip), put `gatk` on PATH or pass `gatk_path`. Outputs go to new paths (nothing is overwritten); GATK index files are reported in `index_files`. Each result is a dictionary with `status` (`completed` only when the exit code is 0 and the outputs exist), `version`, `command`, `outputs`, `log_tail` and tool-specific values. Unset arguments keep GATK's defaults. A missing `gatk` raises `BackendUnavailableError` with an install hint.

## Example data

The examples run in one folder, in order; later examples reuse files written by earlier ones. `gatk` must be on PATH.

**Example file:** `ref.fa`

```text
>chr1
TTTTTTTATTTCTTGCACATAATAATAAGAATCTAACTTAATAACTTGTGAACCCTATTTAGCCTGGATT
TGGTCTTTTTGAATATCATATTATGTCCCATTTACACGAAATCTGTTGATAATCAGTCTTTCGATCTCCG
AGGTTTCGATTTCCGCTATTTATCCTGATTTCCTGACCTGATTGTCGTAATTTTGTAGTCTAAAACAAGC
AAGATATTAAATCGTCAGCATCATCAGTGCAATAGTACATCATTATAACGTATTTAAAGTTTTAACTGTT
GCAACTCATTTCTTTTACAAAATCGGATAACGTATTGTGTCTAAGTATCTGGTTATACAATCGTTGATGC
TAATTGCTTAAAAAAGGTTTAAATTGGACCACAGCCGGATATCGGCAAGCTGGCTTTCCCGTTCTATCCT
ACTTACTAGTATTGGTCCAGACTCTCTGTAAGCATCCGCACTTCTTTGATAACAATATATTTCTTTTTAA
CGGAAACGTAGGATTTAATCAATTGTTATTAGAGGACCCTCTACTTAAGTTATCATGACAAAAGCCCTGC
CTGTAATCTTTATATATATGTGAATACGCTCGTATATTAT
>chr2
GTGCTGGCTTACTTGCGCATTGAGGGACAACGGTAGTACCCCCGGGGATCGTGGCATTACGGGCCGTTAG
GCTTCCTTGCTGTCATATACCGTGCTGACCTTCCAATCAGGTACAGATCTTCTAGAGGCTTTCATCACGG
CGACAACGGAGCTGAGGGGTAGCGGGCATTAACCGCGCGTCCACGGTAAGTTACAGGTTATGCGAGGATT
ATCGCTGAGTAGATTCAGAGACCCCGAAGTGAACTGAGGTCAGATAGGACTGACGCCGTATCTTGAACAA
ATGCTGGCGTCCGATTGACA
```

**Example file:** `ref.fa.fai`

```text
chr1	600	6	70	71
chr2	300	621	70	71
```

**Example file:** `ref.dict`

```text
@HD	VN:1.6
@SQ	SN:chr1	LN:600
@SQ	SN:chr2	LN:300
```

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

**Example file:** `bad_ref.vcf`

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
chr1	50	.	A	A	50	.	QD=1.5;DP=40	GT:DP	0/1:20	0/0:20
chr1	120	.	T	A	50	.	QD=18.0;DP=40	GT:DP	1/1:20	0/1:20
chr1	200	.	CTA	C	50	.	QD=12.0;DP=40	GT:DP	0/1:20	0/1:20
chr1	310	.	A	C	50	.	QD=25.0;DP=40	GT:DP	0/1:20	1/1:20
chr1	400	.	C	CAC	50	.	QD=9.0;DP=40	GT:DP	0/1:20	0/0:20
chr2	30	.	A	C	50	.	QD=21.0;DP=40	GT:DP	1/1:20	0/0:20
chr2	150	.	A	C	50	.	QD=3.0;DP=40	GT:DP	0/1:20	0/1:20
```

## 1) `validate_vcf_gatk()`

- **Purpose:** Validate a VCF with GATK's strict checks (reference alleles, IDs, allele and chromosome counts).
- **Method:** Runs `gatk ValidateVariants` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.validate_vcf_gatk(vcf: 'str', *, reference: 'str | None' = None, dbsnp: 'str | None' = None, validate_gvcf: 'bool' = False, fail_gvcf_on_overlap: 'bool' = False, do_not_validate_filtered_records: 'bool' = False, warn_on_errors: 'bool' = False, exclude_validation_types: 'list[str] | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF, optionally a reference FASTA (with .fai and .dict) and dbSNP.
- **Output:** `valid`, `problems` (warning/error lines), exit code and log tail.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

for name in ("calls.vcf", "bad_ref.vcf"):
    r = g.validate_vcf_gatk(name, reference="ref.fa")
    print(name, r["status"], r["valid"], r["version"])
```

- **Example output:**

```text
calls.vcf completed True 4.7.0.0
bad_ref.vcf failed False 4.7.0.0
```

## 2) `evaluate_variants_gatk()`

- **Purpose:** Evaluate call sets: counts, Ti/Tv, dbSNP overlap, genotype concordance and more.
- **Method:** Runs `gatk VariantEval` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.evaluate_variants_gatk(eval_vcfs: 'list[str]', reference: 'str', output_path: 'str', *, comparison_vcfs: 'list[str] | None' = None, dbsnp: 'str | None' = None, gold_standard: 'str | None' = None, eval_modules: 'list[str] | None' = None, stratification_modules: 'list[str] | None' = None, use_all_standard_modules: 'bool' = True, use_all_standard_stratifications: 'bool' = True, known_names: 'list[str] | None' = None, merge_evals: 'bool' = False, require_strict_allele_match: 'bool' = False, keep_ac0: 'bool' = False, samples: 'list[str] | None' = None, sample_ploidy: 'int | None' = None, pedigree: 'str | None' = None, ancestral_alignments: 'str | None' = None, known_cnvs: 'str | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** One or more VCFs to evaluate, a reference FASTA, optional comparison/dbSNP/gold-standard VCFs and module selections.
- **Output:** The GATKReport file parsed into `tables` (`columns` and `rows` per table).
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.evaluate_variants_gatk(["calls.vcf"], "ref.fa", "eval.report")
t = r["tables"]["CountVariants"]
row = dict(zip(t["columns"], t["rows"][0]))
print(r["status"], row["nSNPs"], row["nInsertions"], row["nDeletions"], row["nHets"], row["nHomVar"])
```

- **Example output:**

```text
completed 5 1 1 8 3
```

## 3) `filter_vcf_variants_gatk()`

- **Purpose:** Mark variants that fail JEXL expressions on INFO or FORMAT annotations.
- **Method:** Runs `gatk VariantFiltration` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.filter_vcf_variants_gatk(vcf: 'str', output_path: 'str', *, filter_expressions: 'list[str] | None' = None, filter_names: 'list[str] | None' = None, genotype_filter_expressions: 'list[str] | None' = None, genotype_filter_names: 'list[str] | None' = None, invert_filter_expression: 'bool' = False, invert_genotype_filter_expression: 'bool' = False, missing_values_evaluate_as_failing: 'bool' = False, set_filtered_genotype_to_no_call: 'bool' = False, invalidate_previous_filters: 'bool' = False, mask: 'str | None' = None, mask_name: 'str | None' = None, mask_extension: 'int | None' = None, mask_description: 'str | None' = None, filter_not_in_mask: 'bool' = False, cluster_size: 'int | None' = None, cluster_window_size: 'int | None' = None, apply_allele_specific_filters: 'bool' = False, reference: 'str | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF, expressions with one name each, optional mask file and SNP-cluster settings.
- **Output:** The filtered VCF and a `summary` with record count and FILTER value counts.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.filter_vcf_variants_gatk(
    "calls.vcf", "filtered.vcf",
    filter_expressions=["QD < 10.0"], filter_names=["lowQD"],
)
print(r["status"], r["summary"]["records"], r["summary"]["filters"])
```

- **Example output:**

```text
completed 7 {'PASS': 4, 'lowQD': 3}
```

## 4) `annotate_vcf_gatk()`

- **Purpose:** Add annotations (and fields from resource VCFs or dbSNP IDs) to the records of a VCF.
- **Method:** Runs `gatk VariantAnnotator` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.annotate_vcf_gatk(vcf: 'str', output_path: 'str', *, annotations: 'list[str] | None' = None, annotation_groups: 'list[str] | None' = None, annotations_to_exclude: 'list[str] | None' = None, expressions: 'list[str] | None' = None, resources: 'list[str | dict[str, Any]] | None' = None, resource_allele_concordance: 'bool' = False, dbsnp: 'str | None' = None, reads: 'list[str] | None' = None, pedigree: 'str | None' = None, founder_ids: 'list[str] | None' = None, enable_all_annotations: 'bool' = False, disable_tool_default_annotations: 'bool' = False, reference: 'str | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF, annotation or group names, optional reads, resources and dbSNP (indexed feature files).
- **Output:** The annotated VCF and a `summary`.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.annotate_vcf_gatk("calls.vcf", "annotated.vcf", annotations=["ChromosomeCounts"])
first = [x for x in Path("annotated.vcf").read_text().splitlines() if not x.startswith("#")][0]
print(r["status"], first.split("\t")[7])
```

- **Example output:**

```text
completed AC=1;AF=0.250;AN=4;DP=40;QD=1.5
```

## 5) `recalibrate_variants_vqsr()`

- **Purpose:** Build the VQSR Gaussian-mixture model and write the recalibration and tranches files.
- **Method:** Runs `gatk VariantRecalibrator` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.recalibrate_variants_vqsr(vcfs: 'list[str]', resources: 'list[dict[str, Any]]', annotations: 'list[str]', recal_path: 'str', tranches_path: 'str', *, mode: 'str' = 'SNP', truth_sensitivity_tranches: 'list[float] | None' = None, max_gaussians: 'int | None' = None, max_negative_gaussians: 'int | None' = None, max_iterations: 'int | None' = None, k_means_iterations: 'int | None' = None, minimum_bad_variants: 'int | None' = None, maximum_training_variants: 'int | None' = None, bad_lod_score_cutoff: 'float | None' = None, target_titv: 'float | None' = None, standard_deviation_threshold: 'float | None' = None, output_model: 'str | None' = None, use_allele_specific_annotations: 'bool' = False, trust_all_polymorphic: 'bool' = False, ignore_all_filters: 'bool' = False, ignore_filters: 'list[str] | None' = None, reference: 'str | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** VCF(s), tagged training/truth resources, annotation names and the mode (SNP, INDEL or BOTH).
- **Output:** `recal_file`, `tranches_file`, optional `model_file`, index files.
- **Example code:**

```python
from pathlib import Path
import random
import subprocess

rng = random.Random(7)
header = (
    "##fileformat=VCFv4.2\n"
    + "".join(f'##INFO=<ID={t},Number=1,Type=Float,Description="{t}">\n' for t in ("QD", "FS", "MQ", "SOR"))
    + "##contig=<ID=chr1,length=70000>\n"
)
call = [header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"]
truth = ["##fileformat=VCFv4.2\n##contig=<ID=chr1,length=70000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"]
for i in range(900):
    pos, good = 100 + 75 * i, i % 4 != 0
    qd = rng.gauss(18, 3) if good else rng.gauss(4, 2)
    fs = abs(rng.gauss(2, 2)) if good else abs(rng.gauss(30, 10))
    mq = rng.gauss(60, 1) if good else rng.gauss(40, 5)
    sor = abs(rng.gauss(1, 0.5)) if good else abs(rng.gauss(4, 1))
    call.append(f"chr1\t{pos}\t.\tA\tC\t50\t.\tQD={qd:.3f};FS={fs:.3f};MQ={mq:.3f};SOR={sor:.3f}\n")
    if good and i % 3 == 0:
        truth.append(f"chr1\t{pos}\t.\tA\tC\t.\t.\t.\n")
Path("scored.vcf").write_text("".join(call))
Path("truth.vcf").write_text("".join(truth))
for name in ("scored.vcf", "truth.vcf"):
    subprocess.run(["gatk", "IndexFeatureFile", "-I", name], check=True, capture_output=True)
from dnakit.variants import gatk_p02_cli_functions as g

truth = {"path": "truth.vcf", "name": "truth", "tags": {"known": False, "training": True, "truth": True, "prior": 15.0}}
r = g.recalibrate_variants_vqsr(
    ["scored.vcf"], [truth], ["QD", "FS", "MQ", "SOR"], "snp.recal", "snp.tranches",
    truth_sensitivity_tranches=[100.0, 99.0, 90.0], max_gaussians=2, max_negative_gaussians=2,
    minimum_bad_variants=100,
)
print(r["status"], Path(r["recal_file"]).name, Path(r["tranches_file"]).name)
```

- **Example output:**

```text
completed snp.recal snp.tranches
```

## 6) `apply_vqsr()`

- **Purpose:** Filter variants by a truth-sensitivity level or LOD cutoff using the recalibration table.
- **Method:** Runs `gatk ApplyVQSR` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.apply_vqsr(vcf: 'str', recal_file: 'str', output_path: 'str', *, tranches_file: 'str | None' = None, truth_sensitivity_filter_level: 'float | None' = None, lod_score_cutoff: 'float | None' = None, mode: 'str' = 'SNP', use_allele_specific_annotations: 'bool' = False, exclude_filtered: 'bool' = False, ignore_all_filters: 'bool' = False, ignore_filters: 'list[str] | None' = None, reference: 'str | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** The VCF, the recal file, the tranches file and the filter level.
- **Output:** The recalibrated VCF (VQSLOD annotated) and a `summary`.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.apply_vqsr(
    "scored.vcf", "snp.recal", "vqsr.vcf", tranches_file="snp.tranches",
    truth_sensitivity_filter_level=99.0,
)
print(r["status"], r["summary"]["records"], r["summary"]["filters"])
```

- **Example output:**

```text
completed 900 {'PASS': 832, 'VQSRTrancheSNP99.00to100.00': 68}
```

## 7) `filter_variant_tranches()`

- **Purpose:** Apply tranche filtering from a score in INFO and truth resources.
- **Method:** Runs `gatk FilterVariantTranches` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.filter_variant_tranches(vcf: 'str', resources: 'list[str | dict[str, Any]]', output_path: 'str', *, snp_tranches: 'list[float] | None' = None, indel_tranches: 'list[float] | None' = None, info_key: 'str | None' = None, invalidate_previous_filters: 'bool' = False, reference: 'str | None' = None, intervals: 'list[str] | None' = None, exclude_intervals: 'list[str] | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF with a score annotation, resource VCFs, tranche levels and the score key.
- **Output:** The filtered VCF and a `summary`.
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

r = g.filter_variant_tranches("vqsr.vcf", ["truth.vcf"], "tranched.vcf", snp_tranches=[99.9], info_key="VQSLOD")
print(r["status"], r["summary"]["records"])
```

- **Example output:**

```text
completed 900
```

## 8) `score_variants_with_cnn()`

- **Purpose:** Annotate variants with the score of GATK's PyTorch convolutional neural network (CNN_1D or CNN_2D).
- **Method:** Runs `gatk NVScoreVariants` with a fixed argument list built from the typed parameters; only arguments you set are passed.
- **API:** `dnakit.variants.gatk_p02_cli_functions.score_variants_with_cnn(vcf: 'str', reference: 'str', output_path: 'str', *, reads: 'str | None' = None, tensor_type: 'str' = 'reference', batch_size: 'int | None' = None, accelerator: 'str | None' = None, random_seed: 'int | None' = None, gatk_path: 'str | None' = None, java_max_heap_gb: 'int | None' = None, timeout_seconds: 'float' = 86400.0, verbose: bool = True) -> 'dict[str, Any]'`
- **Input:** A VCF whose header defines MQ, DP, SOR, FS, QD, MQRankSum and ReadPosRankSum, a reference FASTA, optional BAM for 2D tensors.
- **Output:** The scored VCF and a `summary`. Needs the GATK Python packages (see the example).
- **Example code:**

```python
from pathlib import Path
from dnakit.variants import gatk_p02_cli_functions as g

import os
import shutil
import sys
import zipfile

archive = Path(shutil.which("gatk")).resolve().parent / "gatkPythonPackageArchive.zip"
zipfile.ZipFile(archive).extractall("gatk_python")
os.environ["PYTHONPATH"] = str(Path("gatk_python").resolve())
os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
os.environ["PATH"] = str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]
lines = Path("calls.vcf").read_text().splitlines()
head = [x for x in lines if x.startswith("##")]
head += [f'##INFO=<ID={t},Number=1,Type=Float,Description="{t}">' for t in ("MQ", "SOR", "FS", "MQRankSum", "ReadPosRankSum")]
snps = [x for x in lines if not x.startswith("#") and len(x.split("\t")[3]) == len(x.split("\t")[4])]
chrom = [x for x in lines if x.startswith("#CHROM")]
Path("cnn_in.vcf").write_text("\n".join(head + chrom + snps) + "\n")
r = g.score_variants_with_cnn("cnn_in.vcf", "ref.fa", "cnn.vcf", accelerator="cpu")
scored = [x for x in Path("cnn.vcf").read_text().splitlines() if "CNN_1D=" in x]
print(r["status"], r["summary"]["records"], len(scored))
```

- **Example output:**

```text
completed 5 5
```
