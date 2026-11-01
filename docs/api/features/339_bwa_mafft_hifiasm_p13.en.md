# BWA, MAFFT and hifiasm command-line facades

This group adds facades for [BWA](https://github.com/lh3/bwa) beyond `bwa mem` (index build stages, `aln`, `samse`/`sampe`, `fastmap`, `pemerge`, `maxk`, and the bwa.kit scripts `bwa-postalt.js`, `typeHLA.sh`, `run-gen-ref`, plus `qualfa2fq.pl`) in `dnakit.alignment`; for [MAFFT](https://mafft.cbrc.jp/alignment/software/) (the named strategies as `mafft_align(strategy=...)`, adding sequences, seeds, profiles, partition trees, sparse core and pairwise distances) in `dnakit.alignment`; and one function `hifiasm_assemble` for [hifiasm](https://github.com/chhylp123/hifiasm) in `dnakit.assembly` whose parameters select the HiFi, ONT, ultra-long, Hi-C, trio, dual-scaffold, purging, read-output, breakpoint, misjoin, telomere and chimera options. `bwa maxk` and `bwa pemerge` are hidden/experimental BWA subcommands. BWA `fa2pac`, `pac2bwt`, `pac2bwtgen`, `bwtupdate` and `bwt2sa` are the stages of `bwa index` and are selected with `bwa_index_stage(stage=...)`.

Install the programs separately: `conda install -c bioconda bwa bwakit mafft hifiasm` (tested with BWA 0.7.19, bwa.kit 0.7.15, MAFFT 7.526 and hifiasm 0.25.0). `qinsi`/`xinsi` need MAFFT built with its extensions; `mafft-sparsecore.rb` needs ruby; `--fastaparttree` needs fasta34. hifiasm reserves 16 GB with its default Bloom filter, so pass `bloom_filter_bits=0` for small genomes.

## Example data

**Example file:** `seqs.fa`

```text
>s0
ATGGCGTACGTTAGCTAGCTAGGATCGATCGGCTAAGCTTAGCGATCGATTAGCCGATCGGGATCTAGCTAGGCTA
>s1
ATGGCGTACGTTAGCTAGCTAGGATCGATCGGCTAAGCTTAGCGATCGATTAGCCGATCGGGATCTAGCTAGGCTA
>s2
ATGGCGTACGTAGCTAGCTAGGATCGATCGGCTAAGCTTAGCGATCGATTAGCCGATCGGGATCTAGCTAGGCTA
>s3
ATGGCGTACGTTAGCTAGCTAGGATCGATCGGCAAGCTTAGCGATCGATTAGCCGATCGGGATCTAGCTAGGCTA
```

**Example file:** `ref.fa`

```text
>chr1
TTTTTTTATTTCTTGCACATAATAATAAGAATCTAACTTAATAACTTGTGAACCCTATTTAGCCTGGATTTGGTCTTTTTGAATATCATATTATGTCCCATTTACACGAAATCTGTTGATAATCAGTCTTTCGATCTCCGTGG
```

**Example file:** `query.fa`

```text
>q1
CTTGTGAACCCTATTTAGCCTGGATTTGGTCTTTTTGAATATC
>q2
ACGTACGTA
```

## 1) `dnakit.alignment.bwa_index_stage()`

- **Purpose:** Run one stage of BWA's index construction by hand (``bwa fa2pac|pac2bwt|pac2bwtgen|bwtupdate|bwt2sa``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.bwa_index_stage(stage, input_path, output_path, *, forward_only=False, sa_interval=32, bwa_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 2) `dnakit.alignment.bwa_aln()`

- **Purpose:** BWA-backtrack alignment of short reads, writing suffix-array intervals (``bwa aln``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.bwa_aln(index_prefix, reads_path, output_path, *, max_diff=0.04, max_gap_opens=1, max_gap_extensions=-1, indel_end_skip=5, max_long_deletion_occurrences=10, seed_length=32, max_seed_diff=2, max_queue_entries=2000000, threads=1, mismatch_penalty=3, gap_open_penalty=11, gap_extension_penalty=4, max_best_hits=30, trim_quality=0, barcode_length=None, log_scaled_gap_penalty=False, non_iterative=False, illumina_1_3=False, bam_input=False, bam_read=None, filter_casava=False, bwa_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 3) `dnakit.alignment.bwa_sai_to_sam()`

- **Purpose:** Generate SAM from the intervals of :func:`bwa_aln`: single-end (``bwa samse``) or paired-end (``bwa sampe``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.bwa_sai_to_sam(index_prefix, sai_paths, reads_paths, output_path, *, read_group=None, max_hits=None, max_insert_size=None, max_occurrences_per_end=None, max_discordant_hits=None, chimeric_rate=None, preload_index=False, disable_smith_waterman=False, disable_insert_size_estimate=False, bwa_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 4) `dnakit.alignment.bwa_fastmap()`

- **Purpose:** Find the super-maximal exact matches (SMEMs) of sequences against a BWA index (``bwa fastmap``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.bwa_fastmap(index_prefix, reads_path, output_path, *, min_smem_length=17, max_interval_size=20, min_interval_size=1, max_mem_length=2147483647, stop_interval_size=0, max_queries_returned=1000, bwa_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 5) `dnakit.alignment.bwa_pemerge()`

- **Purpose:** Merge overlapping paired-end reads into single reads (``bwa pemerge``, experimental in BWA).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.bwa_pemerge(read1_path, output_path, *, read2_path=None, output='all', threads=1, min_overlap=10, max_error_sum=70, bwa_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 6) `dnakit.alignment.bwa_maxk()`

- **Purpose:** Histogram of maximal-exact-match lengths of sequences against an index (``bwa maxk``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.bwa_maxk(bwt_path, sequences_path, output_path, *, self_match=False, bwa_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 7) `dnakit.alignment.bwakit_postalt()`

- **Purpose:** Post-process ``bwa mem`` alignments for ALT contigs, decoys and HLA genes (bwa.kit ``bwa-postalt.js``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.bwakit_postalt(alt_path, alignments_path, output_path, *, hla_prefix=None, min_alt_fraction_for_mapq=1.0, bwakit_directory=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 8) `dnakit.alignment.bwakit_type_hla()`

- **Purpose:** HLA typing of one gene from the reads (or contigs) that match it (bwa.kit ``typeHLA.sh``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.bwakit_type_hla(prefix, gene, *, from_contigs=False, max_genotypes=50, bwakit_directory=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 9) `dnakit.alignment.bwakit_generate_reference()`

- **Purpose:** Download and assemble a human reference analysis set with bwa.kit (``run-gen-ref``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.bwakit_generate_reference(analysis_set, output_directory, *, bwakit_directory=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 10) `dnakit.alignment.bwa_qual_fasta_to_fastq()`

- **Purpose:** Combine a FASTA file and a quality-value file into FASTQ (BWA's ``qualfa2fq.pl``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.bwa_qual_fasta_to_fastq(fasta_path, quality_path, output_path, *, script_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 11) `dnakit.alignment.mafft_align()`

- **Purpose:** Multiple sequence alignment with one of MAFFT's named strategies (``mafft-<strategy>``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.mafft_align(sequences_path, output_path, *, strategy='linsi', threads=1, sequence_type=None, gap_open_penalty=None, gap_extension_offset=None, max_iterations=None, retree=None, reorder=False, adjust_direction=None, output_format='fasta', mafft_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 12) `dnakit.alignment.mafft_add_sequences()`

- **Purpose:** Insert new sequences into an existing alignment without realigning it (``mafft --add`` family).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.mafft_add_sequences(alignment_path, new_sequences_path, output_path, *, mode='add', keep_length=False, six_mer_pair=False, threads=1, sequence_type=None, mafft_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 13) `dnakit.alignment.mafft_seed_alignment()`

- **Purpose:** Align sequences while keeping the columns of one or more fixed seed alignments (``mafft --seed``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.mafft_seed_alignment(sequences_path, seed_alignment_paths, output_path, *, threads=1, sequence_type=None, mafft_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 14) `dnakit.alignment.mafft_profile_alignment()`

- **Purpose:** Align two existing alignments to each other as profiles (``mafft-profile``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.mafft_profile_alignment(alignment1_path, alignment2_path, output_path, *, mafft_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 15) `dnakit.alignment.mafft_partition_tree_align()`

- **Purpose:** Progressive alignment guided by a partitioned tree, for tens of thousands of sequences (``mafft --parttree`` family).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.mafft_partition_tree_align(sequences_path, output_path, *, tree_method='parttree', partition_size=None, retree=2, threads=1, sequence_type=None, mafft_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 16) `dnakit.alignment.mafft_sparsecore_align()`

- **Purpose:** Large-scale alignment by aligning a sparse core first and adding the rest (``mafft-sparsecore.rb``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.mafft_sparsecore_align(sequences_path, output_path, *, core_size=500, core_candidates='50%', seed=0, core_options=None, add_options=None, keep_input_order=False, mark_core=False, mafft_path=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 17) `dnakit.alignment.mafft_pairwise_distances()`

- **Purpose:** Pairwise sequence distances computed from pairwise alignments (``mafft-distance``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.alignment.mafft_pairwise_distances(sequences_path, output_path, *, sequence_type=None, max_rows_returned=10000, mafft_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 18) `dnakit.assembly.hifiasm_assemble()`

- **Purpose:** De novo assembly with hifiasm in any of its modes (``hifiasm``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.assembly.hifiasm_assemble(reads_paths, output_prefix, *, ont=False, ultra_long_paths=None, hic_reads1_paths=None, hic_reads2_paths=None, paternal_kmers_path=None, maternal_kmers_path=None, paternal_reads_list_path=None, maternal_reads_list_path=None, dual_scaffold=False, scaffold_max_gap=None, purge_level=None, purge_similarity=None, purge_min_overlaps=None, haplotype_count=None, write_error_corrected_reads=False, write_overlaps_paf=False, low_quality_percent=None, break_low_coverage=None, break_high_coverage=None, break_exact_overlap_rate=None, misjoin_min_size=None, telomere_motif=None, chimeric_read_support=None, chimeric_flank_length=None, min_ont_read_length=None, min_ont_mean_quality=None, min_ultra_long_length=None, ultra_long_error_rate=None, primary_assembly=False, hom_coverage=None, genome_size=None, bloom_filter_bits=None, threads=1, hifiasm_path=None, timeout_seconds=86400.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## Examples

### Example 1: Align with L-INS-i

- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.alignment import mafft_p13_cli_functions as mf

out = Path(tempfile.mkdtemp()) / "aln.fa"
result = mf.mafft_align("seqs.fa", str(out), strategy="linsi", verbose=False)
print(result["status"], result["sequence_count"])
widths = {len(line) for line in out.read_text().split(">")[1:] for line in [line.split("\n", 1)[1].replace("\n", "")]}
print(len(widths))
```

- **Example output:**

```text
completed 4
1
```

### Example 2: Super-maximal exact matches with BWA fastmap

- **Example code:**

```python
import subprocess
import tempfile
from pathlib import Path
from dnakit.alignment import bwa_p13_cli_functions as bw

subprocess.run(["bwa", "index", "ref.fa"], capture_output=True, check=True)
out = Path(tempfile.mkdtemp()) / "fastmap.txt"
result = bw.bwa_fastmap("ref.fa", "query.fa", str(out), verbose=False)
for query in result["queries"]:
    print(query["query"], query["length"], query["matches"])
```

- **Example output:**

```text
q1 43 [{'start': 0, 'end': 43, 'hits': 1, 'positions': ['chr1:+45']}]
q2 9 []
```

### Example 3: Assemble synthetic HiFi reads with hifiasm

- **Example code:**

```python
import random
import tempfile
from pathlib import Path
from dnakit.assembly import hifiasm_p13_cli_functions as hf

rng = random.Random(1)
genome = "".join(rng.choice("ACGT") for _ in range(60000))
work = Path(tempfile.mkdtemp())
with (work / "hifi.fq").open("w") as handle:
    for i in range(180):
        start = rng.randrange(0, len(genome) - 8000)
        handle.write(f"@r{i}\n{genome[start:start + 8000]}\n+\n{'I' * 8000}\n")
result = hf.hifiasm_assemble([str(work / "hifi.fq")], str(work / "asm"), bloom_filter_bits=0, write_error_corrected_reads=True, verbose=False)
print(result["status"], result["version"])
print(any(item["total_bp"] > 50000 for item in result["contig_summary"]), bool(result["error_corrected_reads"]))
```

- **Example output:**

```text
completed 0.25.0
True True
```
