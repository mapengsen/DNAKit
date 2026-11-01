# GATK variant, GVCF and structural-variant tools (GATK, pack P01)

This group of 30 functions in `dnakit.variants.gatk_p01_variants_cli_functions` runs tools of the Broad Institute's [GATK](https://github.com/broadinstitute/gatk) 4.7.0.0 command line (Picard tools bundled in GATK included). Each function runs one fixed tool, or a family of tools chosen by a mode argument, without a shell and with a timeout and an output size limit. Only the long options that `gatk <Tool> --help` lists are accepted through `options`; options that a parameter already sets cannot be repeated, and GATK's script-executing `JAVASCRIPT_FILE` is not reachable. Outputs are written to explicit new paths and never overwritten.

GATK is an external program: install it (for example `conda install -c bioconda gatk4` or the release zip) with Java 17, then put the `gatk` launcher on PATH, set `DNAKIT_GATK` or pass `gatk_path`. If it cannot be found, `BackendUnavailableError` (code `GATK_UNAVAILABLE`) is raised with the install hint; a failing tool raises `BackendExecutionError` (code `GATK_EXECUTION_FAILED`) with the tail of its log; a tool that exits 0 without writing a requested output raises `GATK_OUTPUT_MISSING`. DNAKit does not bundle GATK. Spark tools run in local mode (`spark-master` can be set through `options`). Reference FASTA files need `.fai` and `.dict` files next to them, BAM/VCF inputs their indexes, as for GATK itself. Results are dictionaries; the section for each function lists the extra fields.

## 1) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_reblock_gvcf()`

- **Purpose:** Condense the homozygous-reference blocks of a single-sample GVCF into fewer, wider blocks.
- **Method:** Run `gatk ReblockGVCF`.
- **API:** `gatk_reblock_gvcf(gvcf_path: 'str', reference_path: 'str', output_path: 'str', *, gvcf_gq_bands: 'Sequence[str] | None' = None, drop_low_quals: 'bool | None' = None, do_qual_score_approximation: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ReblockGVCF.
- **Input:** `gvcf_path` (required): Single-sample GVCF (indexed).; `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New reblocked GVCF.; `gvcf_gq_bands` (optional): GQ band boundaries.; `drop_low_quals` (optional): Drop low-quality variants and non-ref blocks.; `do_qual_score_approximation` (optional): Approximate QUAL from the allele-specific annotations.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the reblocked GVCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_reblock_gvcf

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    gvcf_path = folder / "gvcf_path.dat"; gvcf_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_reblock_gvcf(
            gvcf_path=str(gvcf_path),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 2) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_calibrate_dragstr_model()`

- **Purpose:** Estimate the parameters of the DRAGstr model (STR-aware indel error model) from a BAM.
- **Method:** Run `gatk CalibrateDragstrModel`.
- **API:** `gatk_calibrate_dragstr_model(input_path: 'str', reference_path: 'str', str_table_path: 'str', output_path: 'str', *, threads: 'int | None' = None, min_loci_count: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CalibrateDragstrModel.
- **Input:** `input_path` (required): Input BAM/CRAM.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `str_table_path` (required): STR table from `gatk_compose_str_table`.; `output_path` (required): New DRAGstr parameter file.; `threads` (optional): Number of threads.; `min_loci_count` (optional): Minimum loci per period/repeat combination.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the parameter file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_calibrate_dragstr_model

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    str_table_path = folder / "str_table_path.dat"; str_table_path.write_bytes(b"x")
    try:
        gatk_calibrate_dragstr_model(
            input_path=str(input_path),
            reference_path=str(reference_path),
            str_table_path=str(str_table_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 3) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_combine_gvcfs()`

- **Purpose:** Merge one or more HaplotypeCaller GVCF files into a single GVCF.
- **Method:** Run `gatk CombineGVCFs`.
- **API:** `gatk_combine_gvcfs(gvcf_paths: 'Sequence[str]', reference_path: 'str', output_path: 'str', *, convert_to_base_pair_resolution: 'bool | None' = None, break_bands_at_multiples_of: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CombineGVCFs.
- **Input:** `gvcf_paths` (required): GVCF files (indexed) to merge.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New combined GVCF.; `convert_to_base_pair_resolution` (optional): Write one record per base.; `break_bands_at_multiples_of` (optional): Break reference bands at multiples of this.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the merged GVCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_combine_gvcfs

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    gvcf_paths = folder / "gvcf_paths.dat"; gvcf_paths.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_combine_gvcfs(
            gvcf_paths=[str(gvcf_paths)],
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 4) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_genomicsdb_import()`

- **Purpose:** Import single-sample GVCFs into a GenomicsDB workspace for joint genotyping.
- **Method:** Run `gatk GenomicsDBImport` to create a workspace (`mode=create`) or to add samples to an existing one (`mode=update`). `mode`: `create` = new workspace, `update` = add samples to an existing workspace.
- **API:** `gatk_genomicsdb_import(*, mode: 'str' = 'create', workspace_path: 'str | None' = None, existing_workspace_path: 'str | None' = None, gvcf_paths: 'Sequence[str] | None' = None, sample_name_map: 'str | None' = None, intervals: 'Sequence[str] | None' = None, batch_size: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): GenomicsDBImport.
- **Input:** `workspace_path` (required in some modes): New GenomicsDB workspace directory (created by GATK).; `existing_workspace_path` (required in some modes): Existing workspace to add samples to.; `gvcf_paths` (optional): GVCF files (indexed) to import; use this or `sample_name_map`.; `sample_name_map` (optional): Two-column file `sample<TAB>gvcf path`.; `intervals` (required in some modes): Intervals to import (required for `create`).; `batch_size` (optional): Samples read per batch (0 = all).; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `mode`; `outputs.workspace_path` for a new workspace, with `files` listing it. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_genomicsdb_import

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    gvcf_paths = folder / "gvcf_paths.dat"; gvcf_paths.write_bytes(b"x")
    try:
        gatk_genomicsdb_import(
            workspace_path=str(folder / "workspace_path"),
            gvcf_paths=[str(gvcf_paths)],
            intervals=["chr1"],
            mode='create',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 5) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_genotype_gvcfs()`

- **Purpose:** Perform joint genotyping of GVCFs (standard or the fast GnarlyGenotyper model).
- **Method:** Run `gatk GenotypeGVCFs` (`engine=standard`) or the beta `gatk GnarlyGenotyper` (`engine=gnarly`, expects a GenomicsDB workspace or a combined GVCF with allele-specific annotations); give the intervals through the `intervals` option. `engine`: `standard` = GenotypeGVCFs, `gnarly` = GnarlyGenotyper.
- **API:** `gatk_genotype_gvcfs(variant_source: 'str', reference_path: 'str', output_path: 'str', *, engine: 'str' = 'standard', dbsnp: 'str | None' = None, include_non_variant_sites: 'bool | None' = None, keep_all_sites: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): GenotypeGVCFs, GnarlyGenotyper.
- **Input:** `variant_source` (required): Combined GVCF, or `gendb://<workspace>` for a GenomicsDB workspace.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New genotyped VCF.; `dbsnp` (optional): dbSNP VCF for the ID column.; `include_non_variant_sites` (optional): Also output non-variant sites.; `keep_all_sites` (optional): Keep all sites, also those without alternate alleles.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `mode`; `outputs.output_path` is the VCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_genotype_gvcfs

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_source = folder / "variant_source.dat"; variant_source.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_genotype_gvcfs(
            variant_source=str(variant_source),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            engine='standard',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 6) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_haplotype_based_variant_recaller()`

- **Purpose:** Calculate the likelihood of each allele of a VCF against the reads, limited by a set of haplotypes.
- **Method:** Run `gatk HaplotypeBasedVariantRecaller` (experimental).
- **API:** `gatk_haplotype_based_variant_recaller(input_path: 'str', reference_path: 'str', alleles_vcf: 'str', haplotypes_bam: 'str', matrix_csv: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): HaplotypeBasedVariantRecaller.
- **Input:** `input_path` (required): Reads BAM/CRAM.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `alleles_vcf` (required): VCF of the alleles to evaluate.; `haplotypes_bam` (required): BAM of the haplotypes.; `matrix_csv` (required): New CSV with the likelihood matrix.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.matrix_csv` is the matrix. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_haplotype_based_variant_recaller

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    alleles_vcf = folder / "alleles_vcf.dat"; alleles_vcf.write_bytes(b"x")
    haplotypes_bam = folder / "haplotypes_bam.dat"; haplotypes_bam.write_bytes(b"x")
    try:
        gatk_haplotype_based_variant_recaller(
            input_path=str(input_path),
            reference_path=str(reference_path),
            alleles_vcf=str(alleles_vcf),
            haplotypes_bam=str(haplotypes_bam),
            matrix_csv=str(folder / "matrix_csv.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 7) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_learn_read_orientation_model()`

- **Purpose:** Get maximum-likelihood estimates of the artifact priors of the read-orientation bias mixture model used by FilterMutectCalls.
- **Method:** Run `gatk LearnReadOrientationModel`.
- **API:** `gatk_learn_read_orientation_model(f1r2_tar_gz_files: 'Sequence[str]', output_path: 'str', *, num_em_iterations: 'int | None' = None, convergence_threshold: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): LearnReadOrientationModel.
- **Input:** `f1r2_tar_gz_files` (required): F1R2 counts archives written by Mutect2 (`f1r2-tar-gz`).; `output_path` (required): New orientation-bias model archive (.tar.gz).; `num_em_iterations` (optional): Maximum EM iterations.; `convergence_threshold` (optional): EM convergence threshold.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the model archive. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_learn_read_orientation_model

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    f1r2_tar_gz_files = folder / "f1r2_tar_gz_files.dat"; f1r2_tar_gz_files.write_bytes(b"x")
    try:
        gatk_learn_read_orientation_model(
            f1r2_tar_gz_files=[str(f1r2_tar_gz_files)],
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 8) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_mutect2()`

- **Purpose:** Call somatic SNVs and indels by local assembly of haplotypes (Mutect2).
- **Method:** Run `gatk Mutect2`; give the regions through the `intervals` option.
- **API:** `gatk_mutect2(input_paths: 'Sequence[str]', reference_path: 'str', output_path: 'str', *, tumor_sample: 'str | None' = None, normal_samples: 'Sequence[str] | None' = None, germline_resource: 'str | None' = None, panel_of_normals: 'str | None' = None, f1r2_tar_gz: 'str | None' = None, bam_output: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): Mutect2.
- **Input:** `input_paths` (required): Tumor (and normal) BAM/CRAM files.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New somatic VCF.; `tumor_sample` (optional): Tumor sample name (read-group SM).; `normal_samples` (optional): Normal sample name(s).; `germline_resource` (optional): Population allele-frequency VCF (for example gnomAD).; `panel_of_normals` (optional): Panel-of-normals VCF.; `f1r2_tar_gz` (optional): New archive of F1R2 counts for `gatk_learn_read_orientation_model`.; `bam_output` (optional): New BAM of the locally re-assembled reads.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the VCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_mutect2

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_paths = folder / "input_paths.dat"; input_paths.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_mutect2(
            input_paths=[str(input_paths)],
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 9) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_reads_pipeline_spark()`

- **Purpose:** Run BWA (optional), MarkDuplicates, BQSR and HaplotypeCaller on reads in one Spark pipeline to produce a VCF/GVCF.
- **Method:** Run `gatk ReadsPipelineSpark` (beta, Spark local mode).
- **API:** `gatk_reads_pipeline_spark(input_path: 'str', reference_path: 'str', known_sites: 'Sequence[str]', output_path: 'str', *, output_bam: 'str | None' = None, align: 'bool | None' = None, emit_ref_confidence: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ReadsPipelineSpark.
- **Input:** `input_path` (required): Unaligned or aligned BAM/CRAM.; `reference_path` (required): Reference FASTA (with .fai and .dict). A .2bit file also works for Spark.; `known_sites` (required): Known-variant VCF files (indexed) for BQSR.; `output_path` (required): New VCF/GVCF.; `output_bam` (optional): New BAM of the processed reads.; `align` (optional): Run BWA-MEM first.; `emit_ref_confidence` (optional): Reference-confidence mode of the caller.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the variant file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_reads_pipeline_spark

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    known_sites = folder / "known_sites.dat"; known_sites.write_bytes(b"x")
    try:
        gatk_reads_pipeline_spark(
            input_path=str(input_path),
            reference_path=str(reference_path),
            known_sites=[str(known_sites)],
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 10) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_collect_sv_evidence()`

- **Purpose:** Gather paired-end, split-read and depth evidence files from a BAM for the GATK-SV pipeline.
- **Method:** Run `gatk CollectSVEvidence` (beta); the evidence files that are requested are written, and at least one must be requested.
- **API:** `gatk_collect_sv_evidence(input_path: 'str', sample_name: 'str', *, pe_file: 'str | None' = None, sr_file: 'str | None' = None, sd_file: 'str | None' = None, site_depth_locs_vcf: 'str | None' = None, depth_evidence_file: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CollectSVEvidence.
- **Input:** `input_path` (required): Coordinate-sorted BAM/CRAM.; `sample_name` (required): Sample name written to the evidence files.; `pe_file` (optional): New paired-end (discordant pair) evidence file.; `sr_file` (optional): New split-read evidence file.; `sd_file` (optional): New site-depth evidence file.; `site_depth_locs_vcf` (optional): VCF of the sites for site-depth counting.; `depth_evidence_file` (optional): New read-depth evidence file.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs` lists the evidence files. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_collect_sv_evidence

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_collect_sv_evidence(
            input_path=str(input_path),
            sample_name="S1",
            pe_file=str(folder / "pe_file.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 11) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_condense_depth_evidence()`

- **Purpose:** Merge adjacent depth-evidence records of a GATK-SV depth file.
- **Method:** Run `gatk CondenseDepthEvidence` (experimental).
- **API:** `gatk_condense_depth_evidence(depth_evidence: 'str', output_path: 'str', *, max_interval_size: 'int | None' = None, min_interval_size: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CondenseDepthEvidence.
- **Input:** `depth_evidence` (required): Depth evidence file (`.rd.txt`).; `output_path` (required): New condensed depth evidence file.; `max_interval_size` (optional): Maximum merged interval size.; `min_interval_size` (optional): Minimum interval size.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the condensed file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_condense_depth_evidence

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    depth_evidence = folder / "depth_evidence.dat"; depth_evidence.write_bytes(b"x")
    try:
        gatk_condense_depth_evidence(
            depth_evidence=str(depth_evidence),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 12) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_convert_counts_to_depth_file()`

- **Purpose:** Convert a GATK read-count file (`*.counts.tsv`) to a GATK-SV depth file (`.rd.txt`).
- **Method:** Run `gatk ConvertCountsToDepthFile` (experimental).
- **API:** `gatk_convert_counts_to_depth_file(counts_path: 'str', output_path: 'str', *, sample_name: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ConvertCountsToDepthFile.
- **Input:** `counts_path` (required): Read-count file (`*.counts.tsv`).; `output_path` (required): New depth evidence file (`.rd.txt`).; `sample_name` (optional): Sample name written to the file.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the depth file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_convert_counts_to_depth_file

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    counts_path = folder / "counts_path.dat"; counts_path.write_bytes(b"x")
    try:
        gatk_convert_counts_to_depth_file(
            counts_path=str(counts_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 13) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_sv_cluster()`

- **Purpose:** Cluster structural variants by type, size and breakpoints (optionally grouped by stratum and track overlap).
- **Method:** Run `gatk SVCluster` (`mode=cluster`) or `GroupedSVCluster` (`mode=grouped`, which also needs the clustering and stratification configurations); both are beta tools. `mode`: `cluster` = SVCluster, `grouped` = GroupedSVCluster.
- **API:** `gatk_sv_cluster(variant_paths: 'Sequence[str]', reference_path: 'str', ploidy_table: 'str', output_path: 'str', *, mode: 'str' = 'cluster', clustering_config: 'str | None' = None, stratify_config: 'str | None' = None, algorithm: 'str | None' = None, variant_prefix: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): GroupedSVCluster, SVCluster.
- **Input:** `variant_paths` (required): Structural-variant VCF files.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `ploidy_table` (required): Tab-separated sample/contig ploidy table.; `output_path` (required): New clustered VCF.; `clustering_config` (required in some modes): Clustering configuration table (per variant group).; `stratify_config` (required in some modes): Stratification configuration table.; `algorithm` (optional): Clustering algorithm.; `variant_prefix` (optional): Prefix of the new variant IDs.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `mode`; `outputs.output_path` is the clustered VCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_sv_cluster

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_paths = folder / "variant_paths.dat"; variant_paths.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    ploidy_table = folder / "ploidy_table.dat"; ploidy_table.write_bytes(b"x")
    try:
        gatk_sv_cluster(
            variant_paths=[str(variant_paths)],
            reference_path=str(reference_path),
            ploidy_table=str(ploidy_table),
            output_path=str(folder / "output_path.out"),
            mode='cluster',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 14) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_joint_germline_cnv_segmentation()`

- **Purpose:** Combine segmented germline copy-number (gCNV) VCFs of many samples into one cohort VCF.
- **Method:** Run `gatk JointGermlineCNVSegmentation` (beta).
- **API:** `gatk_joint_germline_cnv_segmentation(variant_paths: 'Sequence[str]', pedigree: 'str', reference_path: 'str', output_path: 'str', *, model_call_intervals: 'str | None' = None, min_sample_set_fraction_overlap: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): JointGermlineCNVSegmentation.
- **Input:** `variant_paths` (required): Segmented gCNV VCFs of the samples.; `pedigree` (required): Pedigree (PED) file.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New combined VCF.; `model_call_intervals` (optional): Interval list of the gCNV model calls.; `min_sample_set_fraction_overlap` (optional): Minimum overlap fraction to merge segments.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the combined VCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_joint_germline_cnv_segmentation

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_paths = folder / "variant_paths.dat"; variant_paths.write_bytes(b"x")
    pedigree = folder / "pedigree.dat"; pedigree.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_joint_germline_cnv_segmentation(
            variant_paths=[str(variant_paths)],
            pedigree=str(pedigree),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 15) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_print_read_counts()`

- **Purpose:** Print per-sample count files for CNV determination from a read-count file.
- **Method:** Run `gatk PrintReadCounts` (experimental). Output files are named by GATK from `output_prefix`.
- **API:** `gatk_print_read_counts(counts_path: 'str', *, output_prefix: 'str | None' = None, output_file_list: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): PrintReadCounts.
- **Input:** `counts_path` (required): Input read-count file.; `output_prefix` (optional): Prefix of the output paths `<prefix><sample>.counts.tsv`; include a folder that exists.; `output_file_list` (optional): New two-column list `sample<TAB>output path`.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_file_list` when requested. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_print_read_counts

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    counts_path = folder / "counts_path.dat"; counts_path.write_bytes(b"x")
    try:
        gatk_print_read_counts(
            counts_path=str(counts_path),
            output_prefix="x",
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 16) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_print_sv_evidence()`

- **Purpose:** Merge or subset GATK-SV evidence records into one output file.
- **Method:** Run `gatk PrintSVEvidence` (experimental).
- **API:** `gatk_print_sv_evidence(evidence_file: 'str', output_path: 'str', *, sample_names: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): PrintSVEvidence.
- **Input:** `evidence_file` (required): SV evidence file (`.pe.txt`, `.sr.txt`, `.rd.txt`, `.sd.txt`, `.baf.txt`).; `output_path` (required): New evidence file.; `sample_names` (optional): Samples to keep.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the evidence file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_print_sv_evidence

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    evidence_file = folder / "evidence_file.dat"; evidence_file.write_bytes(b"x")
    try:
        gatk_print_sv_evidence(
            evidence_file=str(evidence_file),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 17) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_sv_annotate()`

- **Purpose:** Annotate an SV VCF with gene overlaps and predicted variant consequences.
- **Method:** Run `gatk SVAnnotate`.
- **API:** `gatk_sv_annotate(variant_path: 'str', output_path: 'str', *, protein_coding_gtf: 'str | None' = None, non_coding_bed: 'str | None' = None, promoter_window_length: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): SVAnnotate.
- **Input:** `variant_path` (required): Structural-variant VCF from the GATK-SV pipeline.; `output_path` (required): New annotated VCF.; `protein_coding_gtf` (optional): GTF of protein-coding genes.; `non_coding_bed` (optional): BED of non-coding elements.; `promoter_window_length` (optional): Promoter window upstream of a gene (bp).; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the annotated VCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_sv_annotate

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    try:
        gatk_sv_annotate(
            variant_path=str(variant_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 18) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_sv_concordance()`

- **Purpose:** Annotate structural variants with their genotype concordance against a truth set.
- **Method:** Run `gatk SVConcordance` (beta).
- **API:** `gatk_sv_concordance(evaluation: 'str', truth: 'str', output_path: 'str', *, clustering_config: 'str | None' = None, stratify_config: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): SVConcordance.
- **Input:** `evaluation` (required): Evaluation SV VCF.; `truth` (required): Truth SV VCF.; `output_path` (required): New VCF with concordance annotations.; `clustering_config` (optional): Clustering configuration table.; `stratify_config` (optional): Stratification configuration table.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the annotated VCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_sv_concordance

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    evaluation = folder / "evaluation.dat"; evaluation.write_bytes(b"x")
    truth = folder / "truth.dat"; truth.write_bytes(b"x")
    try:
        gatk_sv_concordance(
            evaluation=str(evaluation),
            truth=str(truth),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 19) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_sv_stratify()`

- **Purpose:** Annotate structural variants by SV type, size and reference-track overlap (stratification).
- **Method:** Run `gatk SVStratify` (beta).
- **API:** `gatk_sv_stratify(variant_path: 'str', stratify_config: 'str', output_path: 'str', *, track_name: 'Sequence[str] | None' = None, track_intervals: 'Sequence[str] | None' = None, split_output: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): SVStratify.
- **Input:** `variant_path` (required): SV VCF.; `stratify_config` (required): Stratification configuration table.; `output_path` (required): New stratified VCF.; `track_name` (optional): Names of the reference tracks.; `track_intervals` (optional): Interval files of the reference tracks.; `split_output` (optional): Write one file per stratum.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the stratified VCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_sv_stratify

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    stratify_config = folder / "stratify_config.dat"; stratify_config.write_bytes(b"x")
    try:
        gatk_sv_stratify(
            variant_path=str(variant_path),
            stratify_config=str(stratify_config),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 20) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_site_depth_to_baf()`

- **Purpose:** Convert site-depth evidence to B-allele-frequency (BAF) evidence for GATK-SV.
- **Method:** Run `gatk SiteDepthtoBAF` (experimental).
- **API:** `gatk_site_depth_to_baf(site_depth_paths: 'Sequence[str]', baf_sites_vcf: 'str', baf_evidence_output: 'str', *, sample_names: 'Sequence[str] | None' = None, min_het_probability: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): SiteDepthtoBAF.
- **Input:** `site_depth_paths` (required): Site-depth evidence files (`.sd.txt`).; `baf_sites_vcf` (required): VCF of the candidate heterozygous sites.; `baf_evidence_output` (required): New B-allele-frequency evidence file.; `sample_names` (optional): Sample names of the inputs.; `min_het_probability` (optional): Minimum probability of a heterozygous site.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.baf_evidence_output` is the BAF file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_site_depth_to_baf

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    site_depth_paths = folder / "site_depth_paths.dat"; site_depth_paths.write_bytes(b"x")
    baf_sites_vcf = folder / "baf_sites_vcf.dat"; baf_sites_vcf.write_bytes(b"x")
    try:
        gatk_site_depth_to_baf(
            site_depth_paths=[str(site_depth_paths)],
            baf_sites_vcf=str(baf_sites_vcf),
            baf_evidence_output=str(folder / "baf_evidence_output.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 21) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_structural_variation_discovery_pipeline_spark()`

- **Purpose:** Run the structural-variation discovery workflow (evidence, local assembly, breakpoint calling) on one sample.
- **Method:** Run `gatk StructuralVariationDiscoveryPipelineSpark` (beta, Spark local mode).
- **API:** `gatk_structural_variation_discovery_pipeline_spark(input_path: 'str', reference_path: 'str', aligner_index_image: 'str', kmers_to_ignore: 'str', contig_sam_file: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): StructuralVariationDiscoveryPipelineSpark.
- **Input:** `input_path` (required): Query-name-grouped or coordinate-sorted BAM of one sample.; `reference_path` (required): Reference FASTA (with .fai and .dict). A .2bit file also works for Spark.; `aligner_index_image` (required): BWA-MEM index image of the reference.; `kmers_to_ignore` (required): K-mer list from `gatk_find_bad_genomic_kmers_spark`.; `contig_sam_file` (required): New SAM of the assembled contig alignments.; `output_path` (required): New SV VCF.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the SV VCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_structural_variation_discovery_pipeline_spark

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    aligner_index_image = folder / "aligner_index_image.dat"; aligner_index_image.write_bytes(b"x")
    kmers_to_ignore = folder / "kmers_to_ignore.dat"; kmers_to_ignore.write_bytes(b"x")
    try:
        gatk_structural_variation_discovery_pipeline_spark(
            input_path=str(input_path),
            reference_path=str(reference_path),
            aligner_index_image=str(aligner_index_image),
            kmers_to_ignore=str(kmers_to_ignore),
            contig_sam_file=str(folder / "contig_sam_file.out"),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 22) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_calculate_genotype_posteriors()`

- **Purpose:** Calculate genotype posterior probabilities given family and/or known population genotypes.
- **Method:** Run `gatk CalculateGenotypePosteriors`.
- **API:** `gatk_calculate_genotype_posteriors(variant_path: 'str', output_path: 'str', *, pedigree: 'str | None' = None, supporting_callsets: 'Sequence[str] | None' = None, de_novo_prior: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CalculateGenotypePosteriors.
- **Input:** `variant_path` (required): Input VCF with genotype likelihoods.; `output_path` (required): New VCF with posterior genotype probabilities.; `pedigree` (optional): Pedigree (PED) file for family priors.; `supporting_callsets` (optional): Population VCFs that provide allele-count priors.; `de_novo_prior` (optional): Prior probability of a de novo mutation.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the VCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_calculate_genotype_posteriors

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    try:
        gatk_calculate_genotype_posteriors(
            variant_path=str(variant_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 23) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_concordance()`

- **Purpose:** Evaluate the concordance of a VCF against a validated truth VCF (true/false positives, false negatives).
- **Method:** Run `gatk Concordance`; the summary table is parsed into `tables.summary_path`. Give the regions through the `intervals` option.
- **API:** `gatk_concordance(evaluation: 'str', truth: 'str', summary_path: 'str', *, true_positives_and_false_negatives: 'str | None' = None, true_positives_and_false_positives: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): Concordance.
- **Input:** `evaluation` (required): VCF to evaluate.; `truth` (required): Validated truth VCF.; `summary_path` (required): New summary table (sensitivity, precision).; `true_positives_and_false_negatives` (optional): New VCF of TP and FN records.; `true_positives_and_false_positives` (optional): New VCF of TP and FP records.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `tables.summary_path` (rows per variant type). Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_concordance

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    evaluation = folder / "evaluation.dat"; evaluation.write_bytes(b"x")
    truth = folder / "truth.dat"; truth.write_bytes(b"x")
    try:
        gatk_concordance(
            evaluation=str(evaluation),
            truth=str(truth),
            summary_path=str(folder / "summary_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 24) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_count_false_positives()`

- **Purpose:** Count the PASS variants of a call set (false positives when the sample has no true variants).
- **Method:** Run `gatk CountFalsePositives` (beta).
- **API:** `gatk_count_false_positives(variant_path: 'str', intervals: 'Sequence[str]', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CountFalsePositives.
- **Input:** `variant_path` (required): VCF of calls on a truth sample (all calls are false positives).; `intervals` (required): Intervals or an interval file to restrict the count to.; `output_path` (required): New file for the count.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` holds the count. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_count_false_positives

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    try:
        gatk_count_false_positives(
            variant_path=str(variant_path),
            intervals=["chr1"],
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 25) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_evaluate_info_field_concordance()`

- **Purpose:** Evaluate the concordance of an INFO field of a VCF against the same or another field of a truth VCF.
- **Method:** Run `gatk EvaluateInfoFieldConcordance` (beta).
- **API:** `gatk_evaluate_info_field_concordance(evaluation: 'str', truth: 'str', eval_info_key: 'str', truth_info_key: 'str', summary_path: 'str', *, epsilon: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): EvaluateInfoFieldConcordance.
- **Input:** `evaluation` (required): VCF to evaluate.; `truth` (required): Truth VCF.; `eval_info_key` (required): INFO key in the evaluation VCF.; `truth_info_key` (required): INFO key in the truth VCF.; `summary_path` (required): New summary file.; `epsilon` (optional): Tolerance when comparing values.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.summary_path` is the summary. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_evaluate_info_field_concordance

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    evaluation = folder / "evaluation.dat"; evaluation.write_bytes(b"x")
    truth = folder / "truth.dat"; truth.write_bytes(b"x")
    try:
        gatk_evaluate_info_field_concordance(
            evaluation=str(evaluation),
            truth=str(truth),
            eval_info_key="AF",
            truth_info_key="AF",
            summary_path=str(folder / "summary_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 26) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_reference_block_concordance()`

- **Purpose:** Evaluate the reference-block concordance of a GVCF against a truth GVCF.
- **Method:** Run `gatk ReferenceBlockConcordance`.
- **API:** `gatk_reference_block_concordance(evaluation: 'str', truth: 'str', eval_block_histogram: 'str', truth_block_histogram: 'str', confidence_concordance_histogram: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ReferenceBlockConcordance.
- **Input:** `evaluation` (required): GVCF to evaluate.; `truth` (required): Truth GVCF.; `eval_block_histogram` (required): New histogram of the evaluation reference blocks.; `truth_block_histogram` (required): New histogram of the truth reference blocks.; `confidence_concordance_histogram` (required): New histogram of the confidence concordance.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs` lists the three histograms. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_reference_block_concordance

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    evaluation = folder / "evaluation.dat"; evaluation.write_bytes(b"x")
    truth = folder / "truth.dat"; truth.write_bytes(b"x")
    try:
        gatk_reference_block_concordance(
            evaluation=str(evaluation),
            truth=str(truth),
            eval_block_histogram=str(folder / "eval_block_histogram.out"),
            truth_block_histogram=str(folder / "truth_block_histogram.out"),
            confidence_concordance_histogram=str(folder / "confidence_concordance_histogram.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 27) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_filter_funcotations()`

- **Purpose:** Filter variants on clinically significant Funcotator annotations (ClinVar, allele frequency, ...).
- **Method:** Run `gatk FilterFuncotations` (experimental).
- **API:** `gatk_filter_funcotations(variant_path: 'str', output_path: 'str', ref_version: 'str', allele_frequency_data_source: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): FilterFuncotations.
- **Input:** `variant_path` (required): Funcotator-annotated VCF.; `output_path` (required): New filtered VCF.; `ref_version` (required): Reference build of the annotation.; `allele_frequency_data_source` (required): Population database that holds the allele frequencies.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the filtered VCF. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_filter_funcotations

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    try:
        gatk_filter_funcotations(
            variant_path=str(variant_path),
            output_path=str(folder / "output_path.out"),
            ref_version="hg38",
            allele_frequency_data_source="gnomad",
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 28) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_find_mendelian_violations()`

- **Purpose:** Find Mendelian violations of all types in a VCF with trios.
- **Method:** Run Picard `gatk FindMendelianViolations`; the metrics file is parsed into `tables.output_path`.
- **API:** `gatk_find_mendelian_violations(input_path: 'str', pedigree_path: 'str', output_path: 'str', *, min_gq: 'int | None' = None, min_dp: 'int | None' = None, skip_chroms: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): FindMendelianViolations.
- **Input:** `input_path` (required): VCF with the trio members.; `pedigree_path` (required): PED file of the trios (no genotype columns).; `output_path` (required): New metrics file.; `min_gq` (optional): Minimum genotype quality.; `min_dp` (optional): Minimum depth.; `skip_chroms` (optional): Contigs to skip.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `tables.output_path` (Picard metrics rows). Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_find_mendelian_violations

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    pedigree_path = folder / "pedigree_path.dat"; pedigree_path.write_bytes(b"x")
    try:
        gatk_find_mendelian_violations(
            input_path=str(input_path),
            pedigree_path=str(pedigree_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 29) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_funcotate()`

- **Purpose:** Annotate variants or copy-number segments functionally (genes, consequences, ClinVar, gnomAD ...) with Funcotator.
- **Method:** Run `gatk Funcotator` (`input_kind=variants`) or `gatk FuncotateSegments` (`input_kind=segments`, beta; output format not well defined upstream). `input_kind`: `variants` = Funcotator, `segments` = FuncotateSegments.
- **API:** `gatk_funcotate(input_path: 'str', reference_path: 'str', data_sources_path: 'str', ref_version: 'str', output_path: 'str', output_file_format: 'str', *, input_kind: 'str' = 'variants', transcript_selection_mode: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): FuncotateSegments, Funcotator.
- **Input:** `input_path` (required): VCF (`variants`) or segment file (`segments`) to annotate.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `data_sources_path` (required): Funcotator data-sources folder (see FuncotatorDataSourceDownloader).; `ref_version` (required): Reference build of the input, for example hg19 or hg38.; `output_path` (required): New annotated file.; `output_file_format` (required): Output format.; `transcript_selection_mode` (optional): How the transcript of a variant is chosen.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `mode`; `outputs.output_path` is the annotated file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_funcotate

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    data_sources_path = folder / "data_sources_path"; data_sources_path.mkdir()
    try:
        gatk_funcotate(
            input_path=str(input_path),
            reference_path=str(reference_path),
            data_sources_path=str(data_sources_path),
            ref_version="hg38",
            output_path=str(folder / "output_path.out"),
            output_file_format="VCF",
            input_kind='variants',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 30) `dnakit.variants.gatk_p01_variants_cli_functions.gatk_validate_basic_somatic_short_mutations()`

- **Purpose:** Check somatic short mutations against tumor-normal BAMs of the same samples sequenced differently (validation power).
- **Method:** Run `gatk ValidateBasicSomaticShortMutations` (experimental).
- **API:** `gatk_validate_basic_somatic_short_mutations(variant_path: 'str', input_path: 'str', reference_path: 'str', discovery_sample_name: 'str', val_case_sample_name: 'str', val_control_sample_name: 'str', output_path: 'str', *, summary_path: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ValidateBasicSomaticShortMutations.
- **Input:** `variant_path` (required): VCF of the discovery calls to validate.; `input_path` (required): BAM/CRAM of the validation tumor and normal samples.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `discovery_sample_name` (required): Sample name in the discovery VCF.; `val_case_sample_name` (required): Validation case (tumor) sample name in the BAM.; `val_control_sample_name` (required): Validation control (normal) sample name in the BAM.; `output_path` (required): New validation table.; `summary_path` (optional): New summary file.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the validation table. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.variants.gatk_p01_variants_cli_functions import gatk_validate_basic_somatic_short_mutations

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_validate_basic_somatic_short_mutations(
            variant_path=str(variant_path),
            input_path=str(input_path),
            reference_path=str(reference_path),
            discovery_sample_name="T",
            val_case_sample_name="T2",
            val_control_sample_name="N2",
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```
