# Alignment and search programs: BLAST+, BamTools, HISAT2

These functions in `dnakit.alignment` run [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/), [BamTools](https://github.com/pezmaster31/bamtools) and [HISAT2](https://daehwankimlab.github.io/hisat2/) programs that the user installs (for example `conda install -c bioconda blast bamtools hisat2`). Each public function runs one fixed program or sub-command with an argument list built from typed parameters: no shell, no generic runner, outputs are new files, and only the options you set are passed, so every other option keeps the program's own default. A missing program raises `BackendUnavailableError` with an install hint; `blast_path`, `bamtools_path` and `hisat2_path` name the folder (or program) to use instead of PATH.

Results are dictionaries with `status` (`completed` only when the exit code is 0 and every expected output exists), `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and small parsed results. Tested with BLAST+ 2.17.0, BamTools 2.5.3 and HISAT2 2.2.3. Functions that merge several BLAST+ tasks or sub-commands into one parameter say which values they cover.

## Example data

The setup example below builds a small BLAST nucleotide database `genes` from six random 900-bp genes, a 200-bp query cut from `cds2`, a 22-bp primer cut from `cds1`, a 40-kb genome with a dinucleotide repeat, and a 60-read BAM (40 mapped reads on two references, 20 unmapped). The later examples use these files in the same folder.

- **Example code:**

```python
import random
import subprocess

import pysam

rng = random.Random(1)


def dna(length):
    return "".join(rng.choice("ACGT") for _ in range(length))


genes = [dna(900) for _ in range(6)]
with open("genes.fa", "w") as handle:
    for index, sequence in enumerate(genes):
        handle.write(f">cds{index}\n{sequence}\n")
subprocess.run(["makeblastdb", "-in", "genes.fa", "-dbtype", "nucl", "-out", "genes", "-parse_seqids"], check=True, capture_output=True)
open("query.fa", "w").write(f">q0\n{genes[2][100:300]}\n")
open("primer.fa", "w").write(f">primer\n{genes[1][50:72]}\n")
genome = dna(20000) + "AC" * 400 + dna(20000)
open("genome.fa", "w").write(">chr1\n" + "\n".join(genome[i:i + 70] for i in range(0, len(genome), 70)) + "\n")
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 5000}, {"SN": "chr2", "LN": 3000}]}
with pysam.AlignmentFile("raw.bam", "wb", header=header) as out:
    for index in range(60):
        read = pysam.AlignedSegment(out.header)
        read.query_name = f"r{index}"
        read.query_sequence = dna(50)
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        if index < 40:
            read.reference_id = index % 2
            read.reference_start = 100 + 10 * index
            read.cigar = [(0, 50)]
            read.flag = 16 if index % 4 == 1 else 0
            read.mapping_quality = 60 if index % 3 else 10
            read.set_tag("XX", "a" if index % 2 else "b")
        else:
            read.flag = 4
            read.reference_id = -1
        out.write(read)
pysam.sort("-o", "reads.bam", "raw.bam")
pysam.index("reads.bam")
print("ready")
```

- **Example output:**

```text
ready
```

## 1) `dnakit.alignment.blastn_search()`

- **Purpose:** Run nucleotide-nucleotide BLAST with any `blastn` task against a BLAST database or a subject FASTA file.
- **Method:** Runs `blastn` with `-task` set to one of `blastn` (classic word size 11), `blastn-short` (queries under 50 bases such as primers), `dc-megablast` (discontiguous megablast, cross-species), `megablast` (the BLAST+ default) or `rmblastn` (repeat-element scoring, needs a BLAST+ build with RMBlast). Output formats 6, 7 and 10 are parsed into `rows`.
- **API:** `dnakit.alignment.blastn_search(query, output, *, database=None, subject=None, task='megablast', outfmt=None, evalue=None, word_size=None, max_target_seqs=None, max_hsps=None, perc_identity=None, strand=None, dust=None, gapopen=None, gapextend=None, penalty=None, reward=None, template_type=None, template_length=None, num_threads=None, blast_path=None)`.
- **Input:** A query FASTA, a new output path, exactly one of `database` (prefix made by makeblastdb) or `subject` (FASTA), `task`, `outfmt` (0-18 or a string such as `6 qseqid sseqid pident`), and optional `evalue`, `word_size`, `max_target_seqs`, `max_hsps`, `perc_identity`, `strand`, `dust`, `gapopen`, `gapextend`, `penalty`, `reward`, `template_type`, `template_length`, `num_threads`.
- **Output:** The common result keys plus `task`, and `rows` / `row_count` for tabular formats.
- **Example code:**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

for task in ("blastn", "megablast", "dc-megablast"):
    result = blast.blastn_search("query.fa", f"{task}.tsv", database="genes", task=task, outfmt=6, verbose=False)
    top = result["rows"][0]
    print(task, result["status"], top["sseqid"], top["pident"], top["length"])
short = blast.blastn_search("primer.fa", "short.tsv", subject="genes.fa", task="blastn-short", outfmt="6 qseqid sseqid length", verbose=False)
print(short["rows"][0])
```

- **Example output:**

```text
blastn completed cds2 100.000 200
megablast completed cds2 100.000 200
dc-megablast completed cds2 100.000 200
{'qseqid': 'primer', 'sseqid': 'cds1', 'length': '22'}
```

## 2) `dnakit.alignment.rpstblastn_search()`

- **Purpose:** Search a DNA query, translated in six frames, against a conserved-domain PSSM database.
- **Method:** Runs `rpstblastn` with the query, the RPS-BLAST database prefix (for example NCBI's `Cdd`, made with `makeprofiledb`), the genetic code and strand options.
- **API:** `dnakit.alignment.rpstblastn_search(query, database, output, *, outfmt=None, evalue=None, max_target_seqs=None, query_gencode=None, strand=None, num_threads=None, blast_path=None)`.
- **Input:** A query FASTA, the database prefix, a new output path, and optional `outfmt`, `evalue`, `max_target_seqs`, `query_gencode`, `strand`, `num_threads`.
- **Output:** The common result keys, and `rows` / `row_count` for tabular formats.
- **Example call (not run here):** This call needs the NCBI conserved-domain database (`Cdd`, several gigabytes from `update_blastdb.pl`), so it is not run in the documentation test.

```python
from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.rpstblastn_search("contig.fa", "Cdd", "domains.tsv", outfmt=6, evalue=0.01, query_gencode=1)
print(result["status"], result["row_count"])
```

## 3) `dnakit.alignment.blastdbcmd_fetch_sequences()`

- **Purpose:** Retrieve sequences, or sub-ranges of them, from a BLAST database by identifier.
- **Method:** Runs `blastdbcmd` with `-entry` (a comma-joined list of identifiers, or `all`) plus optional `-range` and `-strand`, or with `-entry_batch` (a file with one identifier per line, optionally followed by a range).
- **API:** `dnakit.alignment.blastdbcmd_fetch_sequences(database, output, *, entries=None, entry_batch=None, database_type='nucl', sequence_range=None, strand=None, outfmt=None, target_only=False, line_length=None, long_seqids=False, blast_path=None)`.
- **Input:** A database prefix, a new output path, exactly one of `entries` or `entry_batch`, and optional `database_type`, `sequence_range` (`10-50`, `10-`), `strand`, `outfmt` (default `%f`, FASTA), `target_only`, `line_length`, `long_seqids`.
- **Output:** The common result keys; the sequences are in the output file.
- **Example code:**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

one = blast.blastdbcmd_fetch_sequences("genes", "one.fa", entries=["cds3"], sequence_range="11-40", verbose=False)
print(one["status"], open("one.fa").read().split("\n")[0])
open("ids.txt", "w").write("cds0\ncds4\n")
batch = blast.blastdbcmd_fetch_sequences("genes", "batch.txt", entry_batch="ids.txt", outfmt="%a %l", verbose=False)
print(open("batch.txt").read().split("\n")[:2])
```

- **Example output:**

```text
completed >cds3:11-40 
['cds0 900', 'cds4 900']
```

## 4) `dnakit.alignment.makembindex_build_index()`

- **Purpose:** Build a MegaBLAST index of a nucleotide BLAST database or FASTA file.
- **Method:** Runs `makembindex`. With `input_format='blastdb'` (default) the new index format is written next to the database as `<prefix>.NN.idx`; with `'fasta'` the deprecated old format is written under `output_prefix`. Use the index in blastn with `-use_index`.
- **API:** `dnakit.alignment.makembindex_build_index(input_path, *, input_format='blastdb', output_prefix=None, legacy=None, idmap=False, nmer=None, ws_hint=None, volsize=None, stride=None, blast_path=None)`.
- **Input:** A database prefix (or FASTA file), and optional `input_format`, `output_prefix` (fasta only), `legacy`, `idmap`, `nmer`, `ws_hint`, `volsize`, `stride`.
- **Output:** The common result keys; `outputs` lists the index volumes.
- **Example code:**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.makembindex_build_index("genes", verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **Example output:**

```text
completed ['genes.00.idx']
```

## 5) `dnakit.alignment.windowmasker_count_units()`

- **Purpose:** Generate the k-mer unit counts that WindowMasker needs for repeat masking.
- **Method:** Runs `windowmasker -mk_counts` over the genome sequences and writes the count file in the chosen format.
- **API:** `dnakit.alignment.windowmasker_count_units(input_path, output, *, unit_length=None, genome_size=None, input_format='fasta', input_is_list=False, check_duplicates=False, memory_mb=None, counts_format='ascii', counts_size_mb=None, metadata=None, blast_path=None)`.
- **Input:** A FASTA (or a list of FASTA files with `input_is_list`), a new output path, and optional `unit_length` (1-16), `genome_size`, `input_format`, `check_duplicates`, `memory_mb`, `counts_format` (`ascii`, `binary`, `oascii`, `obinary`), `counts_size_mb`, `metadata`.
- **Output:** The common result keys; the counts are in the output file.
- **Example code:**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.windowmasker_count_units("genome.fa", "counts.txt", verbose=False)
print(result["status"], open("counts.txt").read().splitlines()[0])
```

- **Example output:**

```text
completed ##windowmasker-statistics-format-version:ascii 1.1.0
```

## 6) `dnakit.alignment.windowmasker_convert_counts()`

- **Purpose:** Convert a WindowMasker unit-count file between the ascii, binary and optimised formats.
- **Method:** Runs `windowmasker -convert` with the target `-sformat`. The optimised formats (`oascii`, `obinary`) load faster.
- **API:** `dnakit.alignment.windowmasker_convert_counts(input_path, output, *, counts_format='ascii', counts_size_mb=None, metadata=None, blast_path=None)`.
- **Input:** A count file, a new output path, `counts_format`, and optional `counts_size_mb`, `metadata`.
- **Output:** The common result keys; the converted counts are in the output file.
- **Example code:**

```python
import os

from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.windowmasker_convert_counts("counts.txt", "counts.obin", counts_format="obinary", verbose=False)
print(result["status"], os.path.getsize("counts.obin") > 0, open("counts.obin", "rb").read(1)[0])
```

- **Example output:**

```text
completed True 3
```

## 7) `dnakit.alignment.windowmasker_mask_sequences()`

- **Purpose:** Mask repeats in genome sequences with WindowMasker.
- **Method:** Runs `windowmasker -ustat <counts>` and writes intervals (default), soft-masked FASTA or one of the maskinfo/seqloc formats. Thresholds left as None keep the values stored in the count file.
- **API:** `dnakit.alignment.windowmasker_mask_sequences(input_path, unit_counts, output, *, input_format='fasta', output_format='interval', window=None, t_thres=None, t_extend=None, t_low=None, t_high=None, set_t_low=None, set_t_high=None, dust=False, dust_level=None, parse_seqids=False, blast_path=None)`.
- **Input:** A FASTA, a count file from `windowmasker_count_units`, a new output path, `input_format`, `output_format`, and optional `window`, `t_thres`, `t_extend`, `t_low`, `t_high`, `set_t_low`, `set_t_high`, `dust`, `dust_level`, `parse_seqids`.
- **Output:** The common result keys; the masking is in the output file.
- **Example code:**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

intervals = blast.windowmasker_mask_sequences("genome.fa", "counts.obin", "mask.int", verbose=False)
soft = blast.windowmasker_mask_sequences("genome.fa", "counts.txt", "soft.fa", output_format="fasta", verbose=False)
lines = open("mask.int").read().splitlines()
lowered = sum(char.islower() for char in open("soft.fa").read().split("\n", 1)[1])
print(intervals["status"], soft["status"], lines[0], lowered > 0)
```

- **Example output:**

```text
completed completed >chr1 True
```

## 8) `dnakit.alignment.convert2blastmask_from_fasta()`

- **Purpose:** Convert lower-case masking in a FASTA file into a mask file that `makeblastdb -mask_data` accepts.
- **Method:** Runs `convert2blastmask` with the masking algorithm name and a free-text description of its options.
- **API:** `dnakit.alignment.convert2blastmask_from_fasta(input_path, output, *, masking_algorithm, masking_options, output_format='maskinfo_asn1_text', parse_seqids=False, blast_path=None)`.
- **Input:** A soft-masked FASTA, a new output path, `masking_algorithm` (`dust`, `seg`, `windowmasker`, `repeat`, `other`), `masking_options`, and optional `output_format`, `parse_seqids`.
- **Output:** The common result keys; the mask file is the output.
- **Example code:**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.convert2blastmask_from_fasta("soft.fa", "mask.asn", masking_algorithm="windowmasker", masking_options="-mk_counts", verbose=False)
print(result["status"], open("mask.asn").read().split(" ::=")[0])
```

- **Example output:**

```text
completed Blast-db-mask-info
```

## 9) `dnakit.alignment.blast_vdb_search()`

- **Purpose:** Search SRA run accessions with a nucleotide query (`blastn_vdb`) or a protein query (`tblastn_vdb`).
- **Method:** Runs `blastn_vdb` or `tblastn_vdb` with the accessions as the database. The programs read the runs from NCBI over the network.
- **API:** `dnakit.alignment.blast_vdb_search(query, accessions, output, *, program='blastn_vdb', outfmt=None, evalue=None, max_target_seqs=None, num_threads=None, task=None, genetic_code=None, sra_mode=None, include_filtered_reads=False, blast_path=None)`.
- **Input:** A query FASTA, a list of SRA run accessions, a new output path, `program`, and optional `outfmt`, `evalue`, `max_target_seqs`, `num_threads`, `task` (blastn_vdb), `genetic_code` (tblastn_vdb), `sra_mode`, `include_filtered_reads`.
- **Output:** The common result keys, and `rows` / `row_count` for tabular formats.
- **Example call (not run here):** This call downloads SRA data from NCBI and is not run in the documentation test.

```python
from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.blast_vdb_search("read.fa", ["SRR000001"], "vdb_hits.tsv", outfmt=6, max_target_seqs=5)
print(result["status"], result["row_count"])
```

## 10) `dnakit.alignment.blast_vdb_cmd_fetch_reads()`

- **Purpose:** Retrieve sequence records, run information or reference sequences from SRA run accessions.
- **Method:** Runs `blast_vdb_cmd`. `view='sequences'` writes the reads named by `entries` (or `['all']`) or `entry_batch`; `'info'` prints the number of sequences and bases of each run; `'references'` dumps reference sequences of aligned runs; `'paths'` and `'paths_all'` list the run files. The runs are read from NCBI over the network.
- **API:** `dnakit.alignment.blast_vdb_cmd_fetch_reads(accessions, output, *, view='sequences', entries=None, entry_batch=None, outfmt=None, line_length=None, include_filtered_reads=False, blast_path=None)`.
- **Input:** A list of SRA run accessions, a new output path, `view`, and for sequences exactly one of `entries` or `entry_batch`, with optional `outfmt`, `line_length`, `include_filtered_reads`.
- **Output:** The common result keys; the records are in the output file.
- **Example call (not run here):** This call contacts NCBI over the network and is not run in the documentation test. Note that `blast_vdb_cmd` is in the NCBI BLAST+ archive but not in the bioconda `blast` package.

```python
from dnakit.alignment import blast_p14_cli_functions as blast

info = blast.blast_vdb_cmd_fetch_reads(["SRR000001"], "run_info.txt", view="info")
print(info["status"])
print(open("run_info.txt").read())
```

## 11) `dnakit.alignment.bamtools_convert_alignments()`

- **Purpose:** Convert BAM alignments to JSON, YAML or another text format.
- **Method:** Runs `bamtools convert` with `-format` set to `json` (one JSON object per alignment: name, flag, reference, position, mapping quality, CIGAR, bases, qualities, tags), `yaml`, `bed`, `fasta`, `fastq`, `pileup` or `sam`.
- **API:** `dnakit.alignment.bamtools_convert_alignments(bams, output, *, output_format='json', region=None, no_header=False, reference_fasta=None, print_mapping_quality=False, bamtools_path=None)`.
- **Input:** A list of BAM files (merged), a new output path, `output_format`, and optional `region`, `no_header` (sam), `reference_fasta` and `print_mapping_quality` (pileup).
- **Output:** The common result keys; the converted records are in the output file.
- **Example code:**

```python
import json

from dnakit.alignment import bamtools_p14_cli_functions as bt

bt.bamtools_convert_alignments(["reads.bam"], "chr2.json", output_format="json", region="chr2", verbose=False)
records = [json.loads(line) for line in open("chr2.json") if line.strip()]
bt.bamtools_convert_alignments(["reads.bam"], "chr2.yaml", output_format="yaml", region="chr2", verbose=False)
print(len(records), records[0]["reference"], sorted(records[0])[:4], open("chr2.yaml").read().count("RefName: chr2"))
```

- **Example output:**

```text
20 chr2 ['alignmentFlag', 'cigar', 'filename', 'mapQuality'] 20
```

## 12) `dnakit.alignment.bamtools_count_alignments()`

- **Purpose:** Count the alignments in one or more BAM files, optionally within a region.
- **Method:** Runs `bamtools count`; the files are merged before counting and the printed number is returned.
- **API:** `dnakit.alignment.bamtools_count_alignments(bams, *, region=None, bamtools_path=None)`.
- **Input:** A list of BAM files and an optional `region`.
- **Output:** The common result keys plus `count`.
- **Example code:**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

print(bt.bamtools_count_alignments(["reads.bam"], verbose=False)["count"], bt.bamtools_count_alignments(["reads.bam"], region="chr1", verbose=False)["count"])
```

- **Example output:**

```text
60 20
```

## 13) `dnakit.alignment.bamtools_filter_alignments()`

- **Purpose:** Keep the BAM alignments that match all given property criteria.
- **Method:** Runs `bamtools filter` with the criteria built from the parameters: exact flag, insert size, length and mapping quality patterns (`50`, `>=50`, `<100`, `!0`), name (with `*`), query bases, tag value (`NM:<4`) and the twelve `is_*` flag filters, or with a filter script.
- **API:** `dnakit.alignment.bamtools_filter_alignments(bams, output, *, region=None, script=None, alignment_flag=None, insert_size=None, length=None, map_quality=None, name=None, query_bases=None, tag=None, is_duplicate=None, is_failed_qc=None, is_first_mate=None, is_mapped=None, is_mate_mapped=None, is_mate_reverse_strand=None, is_paired=None, is_primary_alignment=None, is_proper_pair=None, is_reverse_strand=None, is_second_mate=None, is_singleton=None, bamtools_path=None)`.
- **Input:** A list of BAM files, a new output path, at least one criterion or a `script`, and optional `region`.
- **Output:** The common result keys; the kept alignments are in the output BAM.
- **Example code:**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

result = bt.bamtools_filter_alignments(["reads.bam"], "filtered.bam", map_quality=">=60", is_reverse_strand=False, is_mapped=True, verbose=False)
print(result["status"], bt.bamtools_count_alignments(["filtered.bam"], verbose=False)["count"])
```

- **Example output:**

```text
completed 19
```

## 14) `dnakit.alignment.bamtools_index_bam_bti()`

- **Purpose:** Create the BamTools-specific index (`.bti`) of a coordinate-sorted BAM file.
- **Method:** Runs `bamtools index -bti`; the standard BAI index is not made. The index is written next to the BAM and must not exist yet.
- **API:** `dnakit.alignment.bamtools_index_bam_bti(bam, *, bamtools_path=None)`.
- **Input:** A sorted BAM file.
- **Output:** The common result keys; `outputs` holds the `.bti` path.
- **Example code:**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

result = bt.bamtools_index_bam_bti("reads.bam", verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **Example output:**

```text
completed ['reads.bam.bti']
```

## 15) `dnakit.alignment.bamtools_sample_alignments()`

- **Purpose:** Write a random subset of the alignments of one or more BAM files.
- **Method:** Runs `bamtools random` with `-n` and `-seed`. No duplicate check is made, so a draw can contain the same alignment twice.
- **API:** `dnakit.alignment.bamtools_sample_alignments(bams, output, *, count=None, seed=None, region=None, bamtools_path=None)`.
- **Input:** A list of BAM files, a new output path, and optional `count` (BamTools default 10000), `seed`, `region`.
- **Output:** The common result keys; the sampled alignments are in the output BAM.
- **Example code:**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

result = bt.bamtools_sample_alignments(["reads.bam"], "sample.bam", count=15, seed=1, verbose=False)
print(result["status"], bt.bamtools_count_alignments(["sample.bam"], verbose=False)["count"])
```

- **Example output:**

```text
completed 15
```

## 16) `dnakit.alignment.bamtools_resolve_pairs()`

- **Purpose:** Mark the proper-pair flag of paired-end alignments from a fragment-length model.
- **Method:** Runs `bamtools resolve`. `mode='twoPass'` learns the statistics and marks in one command; `'makeStats'` only writes the statistics file; `'markPairs'` marks from an existing statistics file (BamTools also needs its `<stats>.uniq_names.txt`).
- **API:** `dnakit.alignment.bamtools_resolve_pairs(bam, stats, *, output=None, mode='twoPass', min_mapping_quality=None, confidence_interval=None, unused_model_threshold=None, force=False, bamtools_path=None)`.
- **Input:** A BAM file, a `stats` path (new for twoPass/makeStats, existing for markPairs), an `output` BAM (not for makeStats), `mode`, and optional `min_mapping_quality`, `confidence_interval`, `unused_model_threshold`, `force`.
- **Output:** The common result keys plus `mode`.
- **Example code:**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

stats = bt.bamtools_resolve_pairs("reads.bam", "fragments.txt", mode="makeStats", verbose=False)
marked = bt.bamtools_resolve_pairs("reads.bam", "fragments.txt", output="resolved.bam", mode="markPairs", verbose=False)
print(stats["status"], marked["status"], marked["mode"])
```

- **Example output:**

```text
completed completed markPairs
```

## 17) `dnakit.alignment.bamtools_revert_alignments()`

- **Purpose:** Remove duplicate marks and restore the original base qualities (OQ tag) of BAM alignments.
- **Method:** Runs `bamtools revert`; `keep_duplicate_marks` and `keep_qualities` switch off the two actions.
- **API:** `dnakit.alignment.bamtools_revert_alignments(bams, output, *, keep_duplicate_marks=False, keep_qualities=False, bamtools_path=None)`.
- **Input:** A list of BAM files, a new output path, and optional `keep_duplicate_marks`, `keep_qualities`.
- **Output:** The common result keys; the reverted alignments are in the output BAM.
- **Example code:**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

result = bt.bamtools_revert_alignments(["reads.bam"], "reverted.bam", verbose=False)
print(result["status"], bt.bamtools_count_alignments(["reverted.bam"], verbose=False)["count"])
```

- **Example output:**

```text
completed 60
```

## 18) `dnakit.alignment.bamtools_split_alignments()`

- **Purpose:** Split a BAM file into one BAM file per value of a property.
- **Method:** Runs `bamtools split` with `-mapped` (`.MAPPED.bam`, `.UNMAPPED.bam`), `-paired` (`.PAIRED_END.bam`, `.SINGLE_END.bam`), `-reference` (`.REF_<name>.bam`) or `-tag <TAG>` (`.TAG_<tag>_<value>.bam`). Only values that occur get a file.
- **API:** `dnakit.alignment.bamtools_split_alignments(bam, output_stub, *, by='mapped', tag=None, reference_prefix=None, tag_prefix=None, tag_list_delimiter=None, bamtools_path=None)`.
- **Input:** A BAM file, an `output_stub` path prefix, `by` (`mapped`, `paired`, `reference`, `tag`), `tag` for tag splitting, and optional `reference_prefix`, `tag_prefix`, `tag_list_delimiter`.
- **Output:** The common result keys; `outputs` lists the files written.
- **Example code:**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

for by, tag in (("mapped", None), ("reference", None), ("tag", "XX")):
    result = bt.bamtools_split_alignments("reads.bam", f"split_{by}", by=by, tag=tag, verbose=False)
    print(by, [path.rsplit(".", 2)[-2] for path in result["outputs"]])
```

- **Example output:**

```text
mapped ['MAPPED', 'UNMAPPED']
reference ['REF_chr1', 'REF_chr2', 'REF_unmapped']
tag ['TAG_XX_a', 'TAG_XX_b']
```

## 19) `dnakit.alignment.hisat2_extract_splice_sites()`

- **Purpose:** Extract splice sites from a gene annotation in GTF format.
- **Method:** Runs the HISAT2 helper script `hisat2_extract_splice_sites.py`. Each output line is an intron: chromosome, last base of the upstream exon, first base of the downstream exon (0-based) and strand.
- **API:** `dnakit.alignment.hisat2_extract_splice_sites(annotation_gtf, output, *, hisat2_path=None)`.
- **Input:** A GTF file and a new output path.
- **Output:** The common result keys plus `rows` / `row_count`.
- **Example code:**

```python
import random

from dnakit.alignment import hisat2_p14_cli_functions as hs

rng = random.Random(8)
sequence = "".join(rng.choice("ACGT") for _ in range(3000))
open("ref.fa", "w").write(">chr1\n" + "\n".join(sequence[i:i + 60] for i in range(0, 3000, 60)) + "\n")
exons = [(201, 400), (601, 800), (1001, 1200)]
open("genes.gtf", "w").write("".join(f'chr1\tdemo\texon\t{a}\t{b}\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' for a, b in exons))
result = hs.hisat2_extract_splice_sites("genes.gtf", "ss.txt", verbose=False)
print(result["status"], result["rows"])
```

- **Example output:**

```text
completed [['chr1', '399', '600', '+'], ['chr1', '799', '1000', '+']]
```

## 20) `dnakit.alignment.hisat2_extract_exons()`

- **Purpose:** Extract exons from a gene annotation in GTF format.
- **Method:** Runs `hisat2_extract_exons.py`. Each output line is a merged exon: chromosome, start, end (0-based) and strand; it is the input of `hisat2-build --exon`.
- **API:** `dnakit.alignment.hisat2_extract_exons(annotation_gtf, output, *, hisat2_path=None)`.
- **Input:** A GTF file and a new output path.
- **Output:** The common result keys plus `rows` / `row_count`.
- **Example code:**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

result = hs.hisat2_extract_exons("genes.gtf", "exons.txt", verbose=False)
print(result["status"], result["row_count"], result["rows"][0])
```

- **Example output:**

```text
completed 3 ['chr1', '200', '399', '+']
```

## 21) `dnakit.alignment.hisat2_inspect_index()`

- **Purpose:** Print the SNPs, splice sites or exons stored in a HISAT2 index.
- **Method:** Runs `hisat2-inspect` with `--snp`, `--ss`, `--ss-all` or `--exon` (`view` = `snps`, `splice_sites`, `splice_sites_all`, `exons`); `summary`, `names` and `sequences` give the index summary, the reference names and the indexed FASTA.
- **API:** `dnakit.alignment.hisat2_inspect_index(index_prefix, output, *, view='snps', hisat2_path=None)`.
- **Input:** An index prefix made by `hisat2-build`, a new output path and `view`.
- **Output:** The common result keys plus `view`, and `rows` / `row_count` for the four annotation views.
- **Example code:**

```python
import subprocess

from dnakit.alignment import hisat2_p14_cli_functions as hs

subprocess.run(["hisat2-build", "--ss", "ss.txt", "--exon", "exons.txt", "ref.fa", "ref_idx"], check=True, capture_output=True)
for view in ("splice_sites", "exons"):
    result = hs.hisat2_inspect_index("ref_idx", f"inspect_{view}.txt", view=view, verbose=False)
    print(view, result["status"], result["row_count"])
```

- **Example output:**

```text
splice_sites completed 2
exons completed 3
```

## 22) `dnakit.alignment.hisat2_align_reads()`

- **Purpose:** Align reads to an indexed genome with splice-aware HISAT2.
- **Method:** Runs `hisat2 -x <index>` with unpaired (`-U`) or paired (`-1`/`-2`) reads and writes SAM. `no_spliced_alignment` aligns DNA reads; `known_splice_sites`, `novel_splice_sites_out`, `rna_strandness`, `downstream_transcriptome_assembly` (`--dta`) and the other options map to the HISAT2 flags of the same meaning.
- **API:** `dnakit.alignment.hisat2_align_reads(index_prefix, output_sam, *, unpaired_reads=None, reads_1=None, reads_2=None, read_format='fastq', preset=None, threads=None, no_spliced_alignment=False, known_splice_sites=None, novel_splice_sites_out=None, novel_splice_sites_in=None, min_intron_length=None, max_intron_length=None, rna_strandness=None, downstream_transcriptome_assembly=False, dta_cufflinks=False, transcriptome_mapping_only=False, max_alignments=None, min_fragment_length=None, max_fragment_length=None, no_mixed=False, no_discordant=False, trim_5_prime=None, trim_3_prime=None, skip_reads=None, stop_after_reads=None, read_group_id=None, no_softclip=False, no_repeat_index=False, summary_file=None, seed=None, hisat2_path=None)`.
- **Input:** An index prefix, a new SAM path, `unpaired_reads` and/or `reads_1` + `reads_2`, and optional `read_format`, `preset`, `threads`, `known_splice_sites`, `summary_file`, `max_alignments` and the other options listed in the signature.
- **Output:** The common result keys; the SAM file (and the optional summary and splice-site files) are written.
- **Example code:**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

sequence = "".join(line.strip() for line in open("ref.fa").read().split("\n")[1:])
open("reads.fq", "w").write(
    f"@spliced\n{sequence[350:400] + sequence[600:650]}\n+\n{'I' * 100}\n@plain\n{sequence[1500:1600]}\n+\n{'I' * 100}\n"
)
result = hs.hisat2_align_reads("ref_idx", "aligned.sam", unpaired_reads=["reads.fq"], known_splice_sites="ss.txt", verbose=False)
rows = {line.split("\t")[0]: line.split("\t") for line in open("aligned.sam") if not line.startswith("@")}
print(result["status"], rows["plain"][3], rows["plain"][5], rows["spliced"][5])
```

- **Example output:**

```text
completed 1501 100M 50M200N50M
```

## 23) `dnakit.alignment.hisat2_read_statistics()`

- **Purpose:** Count reads and report the minimum, maximum and average read length.
- **Method:** Runs `hisat2_read_statistics.py` over the first `read_count` reads (script default 10000).
- **API:** `dnakit.alignment.hisat2_read_statistics(reads, *, read_count=None, hisat2_path=None)`.
- **Input:** A reads file and an optional `read_count`.
- **Output:** The common result keys plus `reads_examined`, `minimum_length`, `maximum_length`, `average_length`, `read_lengths`.
- **Example code:**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

result = hs.hisat2_read_statistics("reads.fq", verbose=False)
print(result["reads_examined"], result["minimum_length"], result["maximum_length"], result["average_length"])
```

- **Example output:**

```text
2 100 100 100
```

## 24) `dnakit.alignment.hisat2_extract_snps_haplotypes_vcf()`

- **Purpose:** Extract SNPs and haplotypes from VCF files for `hisat2-build`.
- **Method:** Runs `hisat2_extract_snps_haplotypes_VCF.py` and writes `<base>.snp` and `<base>.haplotype` (plus extra files with `extra_files`).
- **API:** `dnakit.alignment.hisat2_extract_snps_haplotypes_vcf(genome_fasta, vcf_files, output_base, *, reference_type=None, inter_gap=None, intra_gap=None, non_rs_ids=False, genotype_vcf=None, genotype_gene_list=None, extra_files=False, hisat2_path=None)`.
- **Input:** A genome FASTA, a list of VCF files, an output base path, and optional `reference_type`, `inter_gap`, `intra_gap`, `non_rs_ids`, `genotype_vcf`, `genotype_gene_list`, `extra_files`.
- **Output:** The common result keys; `outputs` lists the files written.
- **Example code:**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

sequence = "".join(line.strip() for line in open("ref.fa").read().split("\n")[1:])
flip = {"A": "G", "C": "T", "G": "A", "T": "C"}
rows = "".join(f"chr1\t{pos}\trs{pos}\t{sequence[pos - 1]}\t{flip[sequence[pos - 1]]}\t.\t.\t.\n" for pos in (200, 500, 900))
open("variants.vcf", "w").write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + rows)
result = hs.hisat2_extract_snps_haplotypes_vcf("ref.fa", ["variants.vcf"], "vcf_snps", verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
print(open("vcf_snps.snp").read().splitlines()[0].split("\t")[1:4])
```

- **Example output:**

```text
completed ['vcf_snps.haplotype', 'vcf_snps.snp']
['single', 'chr1', '199']
```

## 25) `dnakit.alignment.hisat2_extract_snps_haplotypes_ucsc()`

- **Purpose:** Extract SNPs and haplotypes from a UCSC dbSNP table for `hisat2-build`.
- **Method:** Runs `hisat2_extract_snps_haplotypes_UCSC.py` on a UCSC SNP table such as `snp144Common.txt` (plain or gzip) and writes `<base>.snp` and `<base>.haplotype`.
- **API:** `dnakit.alignment.hisat2_extract_snps_haplotypes_ucsc(genome_fasta, snp_table, output_base, *, inter_gap=None, intra_gap=None, hisat2_path=None)`.
- **Input:** A genome FASTA, a UCSC SNP table, an output base path, and optional `inter_gap`, `intra_gap`.
- **Output:** The common result keys; `outputs` lists the files written.
- **Example code:**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

sequence = "".join(line.strip() for line in open("ref.fa").read().split("\n")[1:])
flip = {"A": "G", "C": "T", "G": "A", "T": "C"}
lines = []
for index, pos in enumerate((100, 300, 900)):
    ref = sequence[pos]
    fields = ["585", "chr1", str(pos), str(pos + 1), f"rs{1000 + index}", "0", "+", ref, ref, f"{ref}/{flip[ref]}", "genomic", "single",
              "unknown", "0", "0", "unknown", "exact", "1", "", "1", "", "2", f"{ref},{flip[ref]},", "10,2,", "0.8,0.2,", "0"]
    lines.append("\t".join(fields))
open("snp_table.txt", "w").write("\n".join(lines) + "\n")
result = hs.hisat2_extract_snps_haplotypes_ucsc("ref.fa", "snp_table.txt", "ucsc_snps", verbose=False)
print(result["status"], len(open("ucsc_snps.snp").read().splitlines()))
```

- **Example output:**

```text
completed 3
```

## 26) `dnakit.alignment.hisat2_simulate_reads()`

- **Purpose:** Simulate RNA-seq reads (default) or DNA-seq reads from a genome, an annotation and SNPs.
- **Method:** Runs `hisat2_simulate_reads.py`; `dna=True` adds `--dna`. The script keeps only chromosomes named 1-22, X and Y, and its RNA-seq mode needs Python 3.10 or older (it calls `random.shuffle(..., random=...)`, removed in 3.11); pass `python_path` to run it with an older interpreter.
- **API:** `dnakit.alignment.hisat2_simulate_reads(genome_fasta, annotation_gtf, snp_file, output_base, *, dna=False, single_end=False, read_length=None, fragment_length=None, num_fragments=None, expression_profile=None, repeat_info=None, error_rate=None, max_mismatch=None, random_seed=None, snp_probability=None, python_path=None, hisat2_path=None)`.
- **Input:** A genome FASTA, a GTF, a SNP file, an output base path, `dna`, `single_end`, and optional `read_length`, `fragment_length`, `num_fragments`, `expression_profile`, `repeat_info`, `error_rate`, `max_mismatch`, `random_seed`, `snp_probability`, `python_path`.
- **Output:** The common result keys plus `read_type`; `outputs` lists the read and truth files.
- **Example code:**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

open("ref_num.fa", "w").write(open("ref.fa").read().replace(">chr1", ">1"))
open("genes_num.gtf", "w").write(open("genes.gtf").read().replace("chr1", "1"))
open("snps_num.snp", "w").write(open("vcf_snps.snp").read().replace("chr1", "1"))
result = hs.hisat2_simulate_reads("ref_num.fa", "genes_num.gtf", "snps_num.snp", "sim", dna=True, single_end=True,
                                  num_fragments=50, read_length=40, random_seed=1, verbose=False)
print(result["status"], result["read_type"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **Example output:**

```text
completed dna ['sim.sam', 'sim_1.fa']
```

## 27) `dnakit.alignment.hisat2_build_repeat_index()`

- **Purpose:** Find repeat sequences in a reference genome and write the files of a HISAT2 repeat index.
- **Method:** Runs `hisat2-repeat` with the seed, repeat-length, repeat-count and edit-distance options. The program is part of the HISAT2 source build; it is not in the bioconda package.
- **API:** `dnakit.alignment.hisat2_build_repeat_index(reference, output_prefix, *, threads=None, seed_length=None, seed_count=None, min_repeat_length=None, max_repeat_length=None, repeat_count=None, repeat_edit_distance=None, forward_only=False, hisat2_path=None)`.
- **Input:** A reference FASTA, a new output prefix, and optional `threads`, `seed_length`, `seed_count`, `min_repeat_length`, `max_repeat_length`, `repeat_count`, `repeat_edit_distance`, `forward_only`.
- **Output:** The common result keys; `outputs` lists the `.rep.*` files.
- **Example call (not run here):** `hisat2-repeat` must be built from the HISAT2 source (`make hisat2-repeat`), so this call is not run in the documentation test; the unit test runs it when the program is on PATH.

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

result = hs.hisat2_build_repeat_index("genome.fa", "genome_repeats", seed_length=30)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```
