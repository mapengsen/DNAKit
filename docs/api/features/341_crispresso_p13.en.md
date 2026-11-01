# CRISPResso2 genome-editing analysis facades and helper-script ports

This group adds [CRISPResso2](https://github.com/pinellolab/CRISPResso2) to `dnakit.design_crispr`: `crispresso_quantify_editing` runs `CRISPResso` for plain indel quantification and, through parameters, for HDR, coding-sequence frameshift, base-editing and prime-editing analysis; `crispresso_batch`, `crispresso_pooled`, `crispresso_wgs`, `crispresso_compare`, `crispresso_pooled_wgs_compare` and `crispresso_aggregate` run the other console programs; and six functions port the helper scripts of the repository (`filterReadsOnSequencePresence.py`, `countHighQualityBases.py`, `count_sgRNA_specific_edits.py`, `reconstituteReads.py`, `plotAmbiguous.py`, `plotCustomAllelePlot.py`). `CRISPRessoMeta` has no console script in 2.3.4 and is not wrapped.

Install CRISPResso2 separately: `conda install -c bioconda crispresso2` (tested with 2.3.4; it brings bowtie2, samtools and fastp). The four script ports that read a finished run or draw plots import the `CRISPResso2` Python package (`pip install "dnakit[crispresso2]"`).

## Example data

**Example file:** `reads.fq`

```text
@a
ACGTGGCACTGCGGCTGGAGGTGGAAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIII
@b
ACGTACGTACGTACGTACGTACGTAAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIII
```

## 1) `dnakit.design_crispr.crispresso_quantify_editing()`

- **Purpose:** Quantify genome-editing outcomes at one or more amplicons from sequencing reads (``CRISPResso``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_quantify_editing(fastq_r1_path, output_folder, *, amplicon_sequences, fastq_r2_path=None, amplicon_names=None, guide_sequences=None, guide_names=None, expected_hdr_amplicon_sequence=None, coding_sequences=None, base_editor_output=False, conversion_nuc_from=None, conversion_nuc_to=None, prime_editing_pegrna_spacer_seq=None, prime_editing_pegrna_extension_seq=None, prime_editing_pegrna_scaffold_seq=None, prime_editing_nicking_guide_seq=None, prime_editing_extension_quantification_window_size=None, prime_editing_scaffold_min_match_length=None, name=None, quantification_window_size=None, quantification_window_center=None, quantification_window_coordinates=None, exclude_bp_from_left=None, exclude_bp_from_right=None, min_average_read_quality=None, min_single_bp_quality=None, min_bp_quality_or_n=None, min_identity_score=None, ignore_substitutions=False, ignore_insertions=False, ignore_deletions=False, discard_indel_reads=False, trim_sequences=False, write_detailed_allele_table=False, fastq_output=False, bam_output=False, vcf_output=False, suppress_plots=False, suppress_report=False, keep_intermediate=False, n_processes=None, crispresso_path=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 2) `dnakit.design_crispr.crispresso_batch()`

- **Purpose:** Run CRISPResso on many samples listed in a tab-separated settings file (``CRISPRessoBatch``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_batch(batch_settings_path, output_folder, *, batch_name=None, amplicon_sequences=None, guide_sequences=None, quantification_window_size=None, suppress_plots=False, suppress_batch_summary_plots=False, n_processes=None, crispresso_path=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 3) `dnakit.design_crispr.crispresso_pooled()`

- **Purpose:** Demultiplex and quantify pooled amplicon sequencing (``CRISPRessoPooled``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_pooled(fastq_r1_path, output_folder, *, amplicons_file_path=None, genome_index=None, gene_annotations_path=None, fastq_r2_path=None, name=None, min_reads_to_use_region=None, demultiplex_only_at_amplicons=False, guide_sequences=None, quantification_window_size=None, suppress_plots=False, n_processes=None, crispresso_path=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 4) `dnakit.design_crispr.crispresso_wgs()`

- **Purpose:** Quantify editing at BED regions of a whole-genome BAM (``CRISPRessoWGS``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_wgs(bam_path, region_file_path, reference_fasta_path, output_folder, *, name=None, gene_annotations_path=None, min_reads_to_use_region=None, guide_sequences=None, quantification_window_size=None, suppress_plots=False, n_processes=None, crispresso_path=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 5) `dnakit.design_crispr.crispresso_compare()`

- **Purpose:** Compare the editing outcomes of two finished CRISPResso analyses (``CRISPRessoCompare``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_compare(folder1, folder2, output_folder, *, sample_names=None, name=None, reported_qvalue_cutoff=None, crispresso_path=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 6) `dnakit.design_crispr.crispresso_pooled_wgs_compare()`

- **Purpose:** Compare two CRISPRessoPooled or CRISPRessoWGS analyses region by region (``CRISPRessoPooledWGSCompare``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_pooled_wgs_compare(folder1, folder2, output_folder, *, sample_names=None, name=None, reported_qvalue_cutoff=None, crispresso_path=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 7) `dnakit.design_crispr.crispresso_aggregate()`

- **Purpose:** Aggregate the editing summaries of many CRISPResso runs into one report (``CRISPRessoAggregate``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_aggregate(search_directory, name, *, folder_prefixes=None, folder_suffix=None, min_reads_for_inclusion=None, suppress_plots=False, crispresso_path=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 8) `dnakit.design_crispr.crispresso_filter_reads_by_sequence()`

- **Purpose:** Keep FASTQ reads that contain (or do not contain) given sequences (port of ``filterReadsOnSequencePresence.py``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_filter_reads_by_sequence(fastq_r1_path, output_r1_path, *, fastq_r2_path=None, output_r2_path=None, include_sequences=(), exclude_sequences=())`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 9) `dnakit.design_crispr.crispresso_count_high_quality_bases()`

- **Purpose:** Count high- and low-quality bases at one reference position in a CRISPResso-annotated FASTQ (port of ``countHighQualityBases.py``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_count_high_quality_bases(annotated_fastq_path, position, original_base, target_base, *, reference_name='Reference', quality_cutoff=30)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 10) `dnakit.design_crispr.crispresso_count_sgrna_specific_edits()`

- **Purpose:** Count the edits that fall in each sgRNA's window in a finished CRISPResso run (port of ``count_sgRNA_specific_edits.py``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_count_sgrna_specific_edits(crispresso_folder)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 11) `dnakit.design_crispr.crispresso_reconstitute_reads()`

- **Purpose:** Regenerate a FASTQ from the allele frequency table of a finished CRISPResso run (port of ``reconstituteReads.py``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_reconstitute_reads(crispresso_folder, output_fastq_path)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 12) `dnakit.design_crispr.crispresso_plot_ambiguous_alleles()`

- **Purpose:** Allele-table plot of the ambiguous reads of a finished CRISPResso run (port of ``plotAmbiguous.py``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_plot_ambiguous_alleles(crispresso_folder, output_root, *, min_frequency=None, max_rows=None, plot_cut_point=True, save_png=False)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 13) `dnakit.design_crispr.crispresso_plot_allele_window()`

- **Purpose:** Allele-table plot over a custom window around each cut site (port of ``plotCustomAllelePlot.py``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.design_crispr.crispresso_plot_allele_window(crispresso_folder, output_root, *, plot_left=20, plot_right=20, plot_center=None, min_frequency=None, max_rows=None, plot_cut_point=True, save_png=False)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## Examples

### Example 1: Quantify editing from synthetic reads

- **Example code:**

```python
import random
import tempfile
from pathlib import Path
from dnakit.design_crispr import crispresso_p13_cli_functions as cr

rng = random.Random(11)
guide = "GGCACTGCGGCTGGAGGTGG"
amplicon = "".join(rng.choice("ACGT") for _ in range(60)) + guide + "GGG" + "".join(rng.choice("ACGT") for _ in range(70))
cut = amplicon.index(guide) + 17
work = Path(tempfile.mkdtemp())
with (work / "r1.fq").open("w") as handle:
    for i in range(200):
        seq = amplicon if rng.random() < 0.6 else amplicon[:cut] + amplicon[cut + 3:]
        handle.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
result = cr.crispresso_quantify_editing(str(work / "r1.fq"), str(work / "out"), amplicon_sequences=amplicon, guide_sequences=guide, suppress_plots=True, suppress_report=True, verbose=False)
row = result["quantification"][0]
print(result["status"], result["version"])
print(row["Reads_in_input"], row["Unmodified"], row["Modified"], row["Deletions"])
```

- **Example output:**

```text
completed 2.3.4
200 107 93 93
```

### Example 2: Keep reads that contain a guide

- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.design_crispr import crispresso_p13_cli_functions as cr

out = Path(tempfile.mkdtemp()) / "kept.fq"
result = cr.crispresso_filter_reads_by_sequence("reads.fq", str(out), include_sequences=["GGCACTGCGGCTGGAGGTGG"], verbose=False)
print(result["reads_read"], result["reads_written"])
print(out.read_text().splitlines()[0])
```

- **Example output:**

```text
2 1
@a
```
