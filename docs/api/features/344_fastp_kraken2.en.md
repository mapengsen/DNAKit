# fastp read processing and Kraken 2 taxonomic classification

Two command-line facades over programs the user installs: [fastp](https://github.com/OpenGene/fastp) read trimming and filtering with the options that `run_fastp` does not expose (`dnakit.sequencing_qc`), and [Kraken 2](https://github.com/DerrickWood/kraken2) classification, inspection, merging and database construction (`dnakit.microbiome`). Install with `conda install -c bioconda fastp kraken2`.

Each function runs one fixed program with an argument list built from typed parameters: no shell, no generic runner, outputs are new files, and only the options you set are passed, so every other option keeps the program's own default. A missing program raises `BackendUnavailableError` with an install hint; `fastp_path` names the fastp executable and `kraken2_path` the folder that holds `kraken2`, `kraken2-build` and `k2`. Results are dictionaries with `status` (`completed` only when the exit code is 0 and every expected output exists), `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and small parsed results. Tested with fastp 1.4.0 and Kraken 2 2.17.2.

`run_fastp` already covers the JSON/HTML report, adapter sequences, paired-end adapter detection, poly-X trimming, sliding-window quality trimming and the unqualified-base limit; `fastp_process_reads` adds adapter FASTA files, low-complexity and index filtering, overlap correction, UMI processing and over-representation analysis (and accepts the common options too). `classify_kraken2` classifies with default options; `kraken2_classify_reads` adds the report, quick mode, several databases and the other options. The download and build functions contact NCBI or the GTDB server and can take hours, so the examples here build a toy database from local sequences instead.

## Example data

The setup example writes 300 read pairs for fastp (100 of them have short inserts, so the reads run into the Illumina adapter `AGATCGGAAGAGC`), an adapter FASTA, and, for Kraken 2, a toy taxonomy (root, Bacteria, genus Genusa with species alpha and beta), two random 4-kb genomes tagged with their taxonomy ids, and 45 reads (20 from each genome plus 5 random ones).

- **Example code:**

```python
import random

rng = random.Random(9)


def dna(length):
    return "".join(rng.choice("ACGT") for _ in range(length))


rc = str.maketrans("ACGT", "TGCA")
adapter = "AGATCGGAAGAGC"
first, second = [], []
for index in range(300):
    fragment = dna(120)
    read1, read2 = fragment[:100], fragment[::-1].translate(rc)[:100]
    if index < 100:
        insert = fragment[:70]
        read1 = (insert + adapter + "ACGTACGTAC" * 3)[:100]
        read2 = (insert[::-1].translate(rc) + adapter + "ACGTACGTAC" * 3)[:100]
    first.append(f"@p{index}/1\n{read1}\n+\n{'I' * 100}\n")
    second.append(f"@p{index}/2\n{read2}\n+\n{'I' * 100}\n")
open("r1.fq", "w").write("".join(first))
open("r2.fq", "w").write("".join(second))
open("adapters.fa", "w").write(f">illumina\n{adapter}\n")
import os

os.makedirs("db/taxonomy")
open("db/taxonomy/nodes.dmp", "w").write("1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n10\t|\t2\t|\tgenus\t|\n11\t|\t10\t|\tspecies\t|\n12\t|\t10\t|\tspecies\t|\n")
names = ((1, "root"), (2, "Bacteria"), (10, "Genusa"), (11, "Genusa alpha"), (12, "Genusa beta"))
open("db/taxonomy/names.dmp", "w").write("".join(f"{i}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for i, n in names))
genomes = {11: dna(4000), 12: dna(4000)}
for taxid, sequence in genomes.items():
    open(f"genome_{taxid}.fa", "w").write(f">seq{taxid}|kraken:taxid|{taxid}\n{sequence}\n")
reads = "".join(f">a{i}\n{genomes[11][s:s + 150]}\n" for i, s in enumerate(range(0, 3000, 150)))
reads += "".join(f">b{i}\n{genomes[12][s:s + 150]}\n" for i, s in enumerate(range(0, 3000, 150)))
reads += "".join(f">u{i}\n{dna(150)}\n" for i in range(5))
open("reads.fa", "w").write(reads)
print("ready")
```

- **Example output:**

```text
ready
```

## 1) `dnakit.sequencing_qc.fastp_process_reads()`

- **Purpose:** Trim and filter FASTQ reads with fastp, including adapter FASTA files, low-complexity and index filtering, overlap correction, UMI processing and over-representation analysis.
- **Method:** Runs `fastp -i ... -o ...` (and `-I`/`-O` for pairs) with the options built from the parameters: `--adapter_fasta`, `--adapter_sequence`, `--detect_adapter_for_pe`, `--low_complexity_filter` with `--complexity_threshold`, `--filter_by_index1/2`, `--correction`, `--umi` with `--umi_loc`, `--umi_len`, `--umi_prefix`, `--overrepresentation_analysis`, poly-X and sliding-window trimming, and the quality filters. The JSON report is parsed into `summary`.
- **API:** `dnakit.sequencing_qc.fastp_process_reads(read1, output_read1, json_report, *, read2=None, output_read2=None, html_report=None, adapter_sequence=None, adapter_sequence_read2=None, adapter_fasta=None, detect_adapter_for_pe=False, disable_adapter_trimming=False, low_complexity_filter=False, complexity_threshold=None, filter_by_index1=None, filter_by_index2=None, index_filter_threshold=None, correction=False, overlap_length_required=None, overlap_diff_limit=None, overlap_diff_percent_limit=None, umi=False, umi_location=None, umi_length=None, umi_prefix=None, umi_skip=None, umi_delimiter=None, overrepresentation_analysis=False, overrepresentation_sampling=None, trim_poly_x=False, cut_front=False, cut_tail=False, cut_right=False, cut_window_size=None, cut_mean_quality=None, qualified_quality_phred=None, unqualified_percent_limit=None, n_base_limit=None, average_quality=None, length_required=None, threads=None, fastp_path=None, timeout_seconds=86400.0)`.
- **Input:** Read 1 (and read 2 with `output_read2`), the new output paths and `json_report`, and the optional parameters named above plus `html_report`, `threads`, `timeout_seconds`.
- **Output:** The common result keys plus `summary` (before/after-filtering statistics, filtering result counts, duplication rate, adapter cutting, over-represented sequences).
- **Example code:**

```python
from dnakit.sequencing_qc import fastp_p14_cli_functions as fp

result = fp.fastp_process_reads("r1.fq", "clean_1.fq", "report.json", read2="r2.fq", output_read2="clean_2.fq", html_report="report.html",
                                adapter_fasta="adapters.fa", threads=2, verbose=False)
summary = result["summary"]
lengths = {len(line) for line in open("clean_1.fq").read().split("\n")[1::4] if line}
print(result["status"], summary["before_filtering"]["total_reads"], summary["after_filtering"]["total_reads"], summary["adapter_cutting"]["adapter_trimmed_reads"])
print(min(lengths), max(lengths))
```

- **Example output:**

```text
completed 600 600 201
70 100
```

## 2) `dnakit.microbiome.kraken2_classify_reads()`

- **Purpose:** Classify reads against one or more Kraken 2 databases and write the per-taxon clade-count report.
- **Method:** Runs `kraken2 --db <db> --output <file>` (with `--report`, `--quick`, `--confidence`, `--minimum-hit-groups`, `--paired`, `--report-minimizer-data` and the other options when set). With several databases, `merge_policy` or `interleaved`, it runs `k2 classify` instead, which resolves the calls of the databases with the lowest common ancestor by default. A protein database classifies DNA reads by six-frame translated search with the same call. The processed/classified counts are read from the log and the report is parsed into `report_rows`.
- **API:** `dnakit.microbiome.kraken2_classify_reads(reads, databases, output, *, report=None, paired=False, interleaved=False, quick=False, confidence=None, minimum_hit_groups=None, minimum_base_quality=None, merge_policy=None, threads=None, use_names=False, memory_mapping=False, report_minimizer_data=False, report_zero_counts=False, use_mpa_style=False, classified_out=None, unclassified_out=None, kraken2_path=None, timeout_seconds=86400.0)`.
- **Input:** A list of read files, a database folder (or a list of them), a new output path, and optional `report`, `paired`, `interleaved`, `quick`, `confidence`, `minimum_hit_groups`, `minimum_base_quality`, `merge_policy`, `threads`, `use_names`, `memory_mapping`, `report_minimizer_data`, `report_zero_counts`, `use_mpa_style`, `classified_out`, `unclassified_out`.
- **Output:** The common result keys plus `databases`, `sequences_processed`, `sequences_classified`, `sequences_unclassified`, `percent_classified` and `report_rows` (percent, clade_reads, taxon_reads, rank, taxid, name).
- **Example code:**

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

added = k2.kraken2_add_to_library("db", ["genome_11.fa", "genome_12.fa"], no_masking=True, verbose=False)
built = k2.kraken2_build_database("db", kmer_length=31, minimizer_length=15, minimizer_spaces=3, threads=2, verbose=False)
result = k2.kraken2_classify_reads(["reads.fa"], "db", "classified.txt", report="report.txt", verbose=False)
print(added["status"], built["status"], result["status"])
print(result["sequences_processed"], result["sequences_classified"], result["sequences_unclassified"])
for row in result["report_rows"]:
    print(row["rank"], row["taxid"], row["clade_reads"], row["taxon_reads"], row["name"])
```

- **Example output:**

```text
completed completed completed
45 40 5
U 0 5 5 unclassified
R 1 40 0 root
D 2 40 0 Bacteria
G 10 40 0 Genusa
S 11 20 20 Genusa alpha
S 12 20 20 Genusa beta
```

## 3) `dnakit.microbiome.kraken2_add_to_library()`

- **Purpose:** Add your own FASTA sequences to the library of a Kraken 2 database.
- **Method:** Runs `kraken2-build --add-to-library <file> --db <db>` once per file (`backend='kraken2-build'`), or `k2 add-to-library --file ...` once (`backend='k2'`). Every sequence needs a taxonomy id in its header (`>name|kraken:taxid|<id>`) or an accession listed in `seqid2taxid.map`. Low-complexity masking is on unless `no_masking`.
- **API:** `dnakit.microbiome.kraken2_add_to_library(database, fasta_files, *, backend='kraken2-build', protein=False, no_masking=False, threads=None, masker_threads=None, kraken2_path=None, timeout_seconds=86400.0)`.
- **Input:** A database folder (created when missing), a list of FASTA files, and optional `backend`, `protein`, `no_masking`, `threads`, `masker_threads` (k2 only).
- **Output:** The common result keys, `runs` (one result per file for `kraken2-build`), and `outputs` holding the `library` folder.
- **Example code:**

```python
import os

from dnakit.microbiome import kraken2_p14_cli_functions as k2

os.makedirs("db_k2/taxonomy")
for name in ("nodes.dmp", "names.dmp"):
    open(f"db_k2/taxonomy/{name}", "w").write(open(f"db/taxonomy/{name}").read())
result = k2.kraken2_add_to_library("db_k2", ["genome_11.fa"], backend="k2", no_masking=True, verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]], sorted(os.listdir("db_k2/library")))
```

- **Example output:**

```text
completed ['library'] ['added']
```

## 4) `dnakit.microbiome.kraken2_build_database()`

- **Purpose:** Build the minimizer hash and taxonomy index of a Kraken 2 database from its library, optionally a protein database for translated search.
- **Method:** Runs `kraken2-build --build --db <db>` (or `k2 build --db <db>`) with the k-mer, minimizer, spaces, load-factor, size, thread, `--fast-build` and `--protein` options when set. The folder must hold `taxonomy/` and `library/`. With `protein=True` the database classifies DNA reads by six-frame translated search.
- **API:** `dnakit.microbiome.kraken2_build_database(database, *, backend='kraken2-build', kmer_length=None, minimizer_length=None, minimizer_spaces=None, load_factor=None, max_db_size=None, threads=None, fast_build=False, protein=False, parts=None, kraken2_path=None, timeout_seconds=86400.0)`.
- **Input:** A database folder, and optional `backend`, `kmer_length` (default 35, 15 for protein), `minimizer_length` (31, 12), `minimizer_spaces` (7, 0), `load_factor` (0.7), `max_db_size`, `threads`, `fast_build`, `protein`, `parts` (k2 only).
- **Output:** The common result keys; `outputs` lists `hash.k2d`, `opts.k2d` and `taxo.k2d`.
- **Example code:**

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

result = k2.kraken2_build_database("db_k2", backend="k2", kmer_length=31, minimizer_length=15, minimizer_spaces=3, threads=2, verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **Example output:**

```text
completed ['hash.k2d', 'opts.k2d', 'taxo.k2d']
```

## 5) `dnakit.microbiome.kraken2_inspect_database()`

- **Purpose:** Summarise a Kraken 2 database and its per-taxon minimizer counts.
- **Method:** Runs `k2 inspect --db <db> --output <file>`, which prints the same report as `kraken2-inspect`: header lines (k-mer and minimizer length, spaced mask, taxonomy nodes, table size) and, unless `skip_counts`, one line per taxon.
- **API:** `dnakit.microbiome.kraken2_inspect_database(database, output, *, skip_counts=False, use_mpa_style=False, report_zero_counts=False, memory_mapping=False, kraken2_path=None, timeout_seconds=3600.0)`.
- **Input:** A built database folder, a new output path, and optional `skip_counts`, `use_mpa_style`, `report_zero_counts`, `memory_mapping`.
- **Output:** The common result keys plus `header` (the `#` lines) and `rows` (percent, clade_minimizers, taxon_minimizers, rank, taxid, name).
- **Example code:**

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

result = k2.kraken2_inspect_database("db", "inspect.txt", verbose=False)
print(result["status"], result["header"][1])
for row in result["rows"]:
    print(row["rank"], row["taxid"], row["name"])
```

- **Example output:**

```text
completed Database options: nucleotide db, k = 31, l = 15
R 1 root
D 2 Bacteria
G 10 Genusa
S 12 Genusa beta
S 11 Genusa alpha
```

## 6) `dnakit.microbiome.kraken2_merge_classifications()`

- **Purpose:** Merge the per-read classifications of the same reads made against different databases, resolving conflicts by the lowest common ancestor.
- **Method:** Runs `k2 merge-classifications --db <db1>,<db2> --output <file> <files...>` (and `--report` when set).
- **API:** `dnakit.microbiome.kraken2_merge_classifications(classification_files, databases, output, *, report=None, use_names=False, use_mpa_style=False, report_zero_counts=False, threads=None, kraken2_path=None, timeout_seconds=86400.0)`.
- **Input:** At least two classification files, the databases that produced them, a new output path, and optional `report`, `use_names`, `use_mpa_style`, `report_zero_counts`, `threads`.
- **Output:** The common result keys; the merged classification (and report) are written.
- **Example code:**

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

k2.kraken2_classify_reads(["reads.fa"], "db", "quick.txt", quick=True, verbose=False)
result = k2.kraken2_merge_classifications(["classified.txt", "quick.txt"], ["db", "db"], "merged.txt", report="merged_report.txt", verbose=False)
print(result["status"], len(open("merged.txt").read().splitlines()))
```

- **Example output:**

```text
completed 45
```

## 7) `dnakit.microbiome.kraken2_download_library()`

- **Purpose:** Download NCBI genomic sequences into the library folder of a Kraken 2 database.
- **Method:** Runs `kraken2-build --download-library <type> --db <db>` or `k2 download-library --library <name|taxid|project|accession>` with masking, FTP, assembly-source, assembly-level, annotation, resume and BLAST-volume options when set. The taxonomy must be downloaded separately before building.
- **API:** `dnakit.microbiome.kraken2_download_library(database, library, *, backend='kraken2-build', protein=False, no_masking=False, use_ftp=False, threads=None, assembly_source=None, assembly_levels=None, has_annotation=False, resume=False, blast_volumes=None, masker_threads=None, kraken2_path=None, timeout_seconds=86400.0)`.
- **Input:** A database folder (created when missing), a `library` (archaea, bacteria, plasmid, viral, human, fungi, plant, protozoa, nr, nt, UniVec, UniVec_Core; with `backend='k2'` also a taxid, project or accession), and optional `backend`, `protein`, `no_masking`, `use_ftp`, `threads`, and the k2-only `assembly_source`, `assembly_levels`, `has_annotation`, `resume`, `blast_volumes`, `masker_threads`.
- **Output:** The common result keys; `outputs` holds the `library` folder.
- **Example call (not run here):** This call downloads from NCBI over the network and is not run in the documentation test; it was checked by hand with `UniVec_Core` (both `kraken2-build --use-ftp` and `k2 download-library`), and the unit test runs it when `DNAKIT_LIVE_TESTS=1`.

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

result = k2.kraken2_download_library("univec_db", "UniVec_Core", use_ftp=True, no_masking=True)
print(result["status"], result["outputs"])
```

## 8) `dnakit.microbiome.kraken2_build_standard_database()`

- **Purpose:** Download and build the standard Kraken 2 database.
- **Method:** Runs `kraken2-build --standard --db <db>` or `k2 build --standard --db <db>`. The standard database holds RefSeq archaea, bacteria, viral, plasmid, human and UniVec_Core sequences with the NCBI taxonomy; the download is tens of gigabytes and the build needs tens of gigabytes of memory and hours of time.
- **API:** `dnakit.microbiome.kraken2_build_standard_database(database, *, backend='kraken2-build', no_masking=False, use_ftp=False, kmer_length=None, minimizer_length=None, minimizer_spaces=None, load_factor=None, max_db_size=None, threads=None, fast_build=False, protein=False, kraken2_path=None, timeout_seconds=86400.0)`.
- **Input:** A database folder, and optional `backend`, `no_masking`, `use_ftp`, `kmer_length`, `minimizer_length`, `minimizer_spaces`, `load_factor`, `max_db_size`, `threads`, `fast_build`, `protein`, `timeout_seconds`.
- **Output:** The common result keys; `outputs` lists `hash.k2d`, `opts.k2d`, `taxo.k2d`.
- **Example call (not run here):** This call downloads tens of gigabytes and builds for hours, so it is not run in the documentation test.

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

result = k2.kraken2_build_standard_database("standard_db", threads=16, use_ftp=True, timeout_seconds=86400)
print(result["status"])
```

## 9) `dnakit.microbiome.kraken2_build_special_database()`

- **Purpose:** Download and build a special reference database: the 16S rRNA databases Greengenes, SILVA and RDP, or the Genome Taxonomy Database (GTDB).
- **Method:** Runs `kraken2-build --special <type> --db <db>` or `k2 build --special <type> --db <db>`. `special='gtdb'` and the `gtdb_*` options (`--gtdb-files`, `--gtdb-use-ncbi-taxonomy`, `--gtdb-representatives-only`, `--gtdb-split-database`, `--gtdb-server`) need `backend='k2'`. RDP's download URLs no longer work upstream.
- **API:** `dnakit.microbiome.kraken2_build_special_database(database, special, *, backend='kraken2-build', kmer_length=None, minimizer_length=None, minimizer_spaces=None, load_factor=None, max_db_size=None, threads=None, fast_build=False, gtdb_files=None, gtdb_use_ncbi_taxonomy=False, gtdb_representatives_only=False, gtdb_split_database=False, gtdb_server=None, no_masking=False, kraken2_path=None, timeout_seconds=86400.0)`.
- **Input:** A database folder, `special` (`greengenes`, `silva`, `rdp`, `gtdb`), and optional `backend`, `kmer_length`, `minimizer_length`, `minimizer_spaces`, `load_factor`, `max_db_size`, `threads`, `fast_build`, `gtdb_files`, `gtdb_use_ncbi_taxonomy`, `gtdb_representatives_only`, `gtdb_split_database`, `gtdb_server`, `no_masking`, `timeout_seconds`.
- **Output:** The common result keys; `outputs` lists `hash.k2d`, `opts.k2d`, `taxo.k2d`.
- **Example call (not run here):** These calls download the 16S or GTDB data over the network and build for a long time, so they are not run in the documentation test.

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

silva = k2.kraken2_build_special_database("silva_db", "silva", threads=8)
gtdb = k2.kraken2_build_special_database("gtdb_db", "gtdb", backend="k2", gtdb_representatives_only=True, threads=16)
print(silva["status"], gtdb["status"])
```
