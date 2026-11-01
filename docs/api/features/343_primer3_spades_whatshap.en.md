# Primer3 design tasks, the SPAdes toolkit and WhatsHap phasing

Three groups of functions that wrap tools the user installs: [Primer3](https://primer3.org/) design tasks, alignment and repeat masking (`dnakit.molecular_lab`); the [SPAdes](https://github.com/ablab/spades) toolkit programs for read correction, assembly graphs and graph alignment (`dnakit.assembly`); and [WhatsHap](https://whatshap.readthedocs.io/) phasing, genotyping and haplotagging (`dnakit.variants`). Install with `conda install -c bioconda primer3 spades whatshap` (primer3-py: extra `bio-primer3`).

The command-line functions run one fixed program or sub-command with an argument list built from typed parameters: no shell, no generic runner, outputs are new files, and only the options you set are passed, so every other option keeps the program's own default. A missing program raises `BackendUnavailableError` with an install hint; `primer3_path`, `spades_path` and `whatshap_path` name the folder (WhatsHap: the executable) to use instead of PATH. Results are dictionaries with `status` (`completed` only when the exit code is 0 and every expected output exists), `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and small parsed results. Tested with primer3-py 2.3.1 (libprimer3 2.6.1), Primer3 2.6.1, SPAdes 4.3.0 and WhatsHap 2.8.

Oligo melting temperature, hairpin and dimer thermodynamics are already provided by `calculate_primer3_tm`, `calculate_primer3_hairpin` and `calculate_primer3_dimer`; plain `whatshap phase` is `phase_variants_whatshap`.

## Example data

The setup example writes a 600-bp random template for Primer3; an 18.8-kb genome with a 400-bp repeat, 5,000 paired 100-bp reads (30% carry one substitution), four 2-kb long reads and a 200-bp query for SPAdes; and, for WhatsHap, a 4-kb reference with 25 heterozygous SNVs, 500 reads drawn alternately from the two haplotypes (as a sorted, indexed BAM) and a VCF. It also phases the VCF once with `whatshap phase` so that the statistics, comparison and haplotagging examples have a phased input. The later examples use these files in the same folder.

- **Example code:**

```python
import random
import subprocess

import pysam

rng = random.Random(2)


def dna(length):
    return "".join(rng.choice("ACGT") for _ in range(length))


open("template.txt", "w").write(dna(600))
# SPAdes data
rc = str.maketrans("ACGT", "TGCA")
repeat = dna(400)
genome = dna(6000) + repeat + dna(6000) + repeat + dna(6000)
open("genome.fa", "w").write(">g\n" + genome + "\n")
first, second = [], []
for index in range(5000):
    start = rng.randrange(0, len(genome) - 300)
    fragment = genome[start:start + 300]
    read1, read2 = fragment[:100], fragment[::-1].translate(rc)[:100]
    if rng.random() < 0.3:
        position = rng.randrange(100)
        read1 = read1[:position] + rng.choice("ACGT") + read1[position + 1:]
    first.append(f"@p{index}/1\n{read1}\n+\n{'I' * 100}\n")
    second.append(f"@p{index}/2\n{read2}\n+\n{'I' * 100}\n")
open("r1.fq", "w").write("".join(first))
open("r2.fq", "w").write("".join(second))
open("long.fa", "w").write("".join(f">l{i}\n{genome[start:start + 2000]}\n" for i, start in enumerate((100, 5000, 12000, 16000))))
open("query.fa", "w").write(">q\n" + genome[3000:3200] + "\n")
# WhatsHap data
length = 4000
reference = dna(length)
open("ref.fa", "w").write(">chr1\n" + "\n".join(reference[i:i + 60] for i in range(0, length, 60)) + "\n")
pysam.faidx("ref.fa")
flip = {"A": "G", "C": "T", "G": "A", "T": "C"}
sites = list(range(100, length - 200, 150))
alternative = list(reference)
for position in sites:
    alternative[position] = flip[reference[position]]
haplotypes = [reference, "".join(alternative)]
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": length}], "RG": [{"ID": "rg1", "SM": "S1"}]}
with pysam.AlignmentFile("raw.bam", "wb", header=header) as out:
    for index in range(500):
        start = rng.randrange(0, length - 300)
        read = pysam.AlignedSegment(out.header)
        read.query_name = f"r{index}"
        read.query_sequence = haplotypes[index % 2][start:start + 300]
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigar = [(0, 300)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 300)
        read.set_tag("RG", "rg1")
        out.write(read)
pysam.sort("-o", "reads.bam", "raw.bam")
pysam.index("reads.bam")
with open("variants.vcf", "w") as handle:
    handle.write(f'##fileformat=VCFv4.2\n##contig=<ID=chr1,length={length}>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')
    for position in sites:
        handle.write(f"chr1\t{position + 1}\t.\t{reference[position]}\t{flip[reference[position]]}\t50\tPASS\t.\tGT\t0/1\n")
subprocess.run(["whatshap", "phase", "-o", "phased.vcf", "--reference", "ref.fa", "variants.vcf", "reads.bam"], check=True, capture_output=True)
print("ready")
```

- **Example output:**

```text
ready
```

## 1) `dnakit.molecular_lab.primer3_design_by_task()`

- **Purpose:** Run Primer3 for any task: pick primer pairs, check given primers, pick a ranked primer list, sequencing, cloning or allele-discriminative primers.
- **Method:** Calls primer3-py `design_primers` with `PRIMER_TASK` set to `generic`, `check_primers`, `pick_primer_list`, `pick_sequencing_primers`, `pick_cloning_primers` or `pick_discriminative_primers` (the `primer3_core PRIMER_TASK=...` runs). Further Primer3 tags go in `settings`; unset tags keep Primer3's defaults. Oligos are returned with 0-based start, length, Tm, GC, penalty and thermodynamic fields.
- **API:** `dnakit.molecular_lab.primer3_design_by_task(template, *, task='generic', left_primer=None, right_primer=None, internal_oligo=None, target=None, included_region=None, excluded_regions=None, product_size_range=None, num_return=None, settings=None)`.
- **Input:** A DNA template, `task`, and optional `left_primer`, `right_primer`, `internal_oligo` (check_primers), `target`, `included_region`, `excluded_regions` (0-based `[start, length]`), `product_size_range`, `num_return`, `settings` (`PRIMER_*` tags).
- **Output:** `status`, `task`, `left`, `right`, `internal` (oligo rows), `pairs`, `explain` (Primer3's counters) and `error` / `warning`.
- **Example code:**

```python
from dnakit.molecular_lab import primer3_p14_functions as p3

template = open("template.txt").read()
ranked = p3.primer3_design_by_task(template, task="pick_primer_list", num_return=3, verbose=False)
first = ranked["left"][0]
print(ranked["status"], len(ranked["left"]), len(ranked["right"]), len(ranked["pairs"]))
print(template[first["start"]:first["start"] + first["length"]] == first["sequence"], first["length"])
checked = p3.primer3_design_by_task(template, task="check_primers", left_primer=first["sequence"], verbose=False)
print(checked["left"][0]["start"] == first["start"], round(checked["left"][0]["tm"], 2) == round(first["tm"], 2), checked["explain"]["PRIMER_LEFT_EXPLAIN"])
sequencing = p3.primer3_design_by_task(template, task="pick_sequencing_primers", target=[300, 20], verbose=False)
print(sequencing["status"], len(sequencing["left"]) > 0)
```

- **Example output:**

```text
completed 3 3 0
True 20
True True considered 1, ok 1
completed True
```

## 2) `dnakit.molecular_lab.ntdpal_align_oligos()`

- **Purpose:** Align two short DNA sequences with Primer3's dynamic-programming aligner in global, global end-anchored, local or local end-anchored mode.
- **Method:** Runs `ntdpal` with mode letter `g` (global), `G` (global end-anchored, also written `e`), `l` (local) or `L` (local end-anchored), and parses the score, length and aligned position pairs.
- **API:** `dnakit.molecular_lab.ntdpal_align_oligos(sequence_1, sequence_2, *, mode='local', gap_open_penalty=None, gap_extend_penalty=None, max_gap=None, gc_weighted_scores=False, ambiguity_codes=False, primer3_path=None)`.
- **Input:** Two DNA strings, `mode`, and optional `gap_open_penalty`, `gap_extend_penalty`, `max_gap`, `gc_weighted_scores` (G/C match 3, A/T 2, mismatch -0.5), `ambiguity_codes`.
- **Output:** The common result keys plus `score`, `length`, `pairs` (0-based index pairs) and, for local modes, `end_1` / `end_2`.
- **Example code:**

```python
from dnakit.molecular_lab import primer3_p14_functions as p3

for mode in ("local", "global", "global_end_anchored", "local_end_anchored"):
    result = p3.ntdpal_align_oligos("GCTAGCTAGGATCC", "GGATCCAAAA", mode=mode, verbose=False)
    print(mode, result["score"], result["length"], result["pairs"][0], result["pairs"][-1])
weighted = p3.ntdpal_align_oligos("AAAGGCCTTT", "GGCC", mode="local", gc_weighted_scores=True, verbose=False)
print(weighted["score"], weighted["end_1"], weighted["end_2"])
```

- **Example output:**

```text
local 6.0 6 [8, 0] [13, 5]
global 6.0 6 [8, 0] [13, 5]
global_end_anchored 6.0 6 [8, 0] [13, 5]
local_end_anchored 6.0 6 [8, 0] [13, 5]
12.0 6 3
```

## 3) `dnakit.molecular_lab.primer3_mask_repeats()`

- **Purpose:** Mask repeat k-mers in FASTA sequences so that primers avoid them, by k-mer probability or by absolute k-mer count.
- **Method:** Runs `primer3_masker`. By default a position is masked when the probability that its k-mer comes from a repeat reaches `probability_cutoff` (0.1), using `<list_prefix>_11.list` and `<list_prefix>_16.list` in `kmer_lists_path`; with `absolute_cutoff` and `kmer_list` the count-based mode (`-a`) masks k-mers that occur at least that often in one list file.
- **API:** `dnakit.molecular_lab.primer3_mask_repeats(input_fasta, output, *, kmer_lists_path, list_prefix='homo_sapiens', probability_cutoff=None, absolute_cutoff=None, kmer_list=None, mask_5_prime=None, mask_3_prime=None, masking_character=None, soft_mask=False, masking_direction=None, primer3_path=None)`.
- **Input:** A FASTA, a new output path, `kmer_lists_path` (the `kmer_lists` folder of the Primer3 source: `homo_sapiens`, `mus_musculus`), and optional `list_prefix`, `probability_cutoff`, `absolute_cutoff` + `kmer_list`, `mask_5_prime`, `mask_3_prime`, `masking_character`, `soft_mask`, `masking_direction`.
- **Output:** The common result keys plus `mode` (`probability` or `absolute_count`); the masked FASTA is the output file.
- **Example call (not run here):** The k-mer list files (`homo_sapiens_11.list`, `homo_sapiens_16.list`, about 95 MB) ship with the Primer3 source and are not installed by conda, so this call is not run in the documentation test.  With the lists, k-mers of the poly-A tail are masked: masked bases become N (lower case with `soft_mask=True`).

```python
from dnakit.molecular_lab import primer3_p14_functions as p3

open("to_mask.fa", "w").write(">s\nACGTACGTAGCTAGCTAGCTAGCTAAAAAAAAAAAAAAAAAAAAAAAA\n")
result = p3.primer3_mask_repeats("to_mask.fa", "masked.fa", kmer_lists_path="kmer_lists")
print(result["status"], result["mode"], open("masked.fa").read().splitlines()[1])
```

## 4) `dnakit.assembly.spades_correct_reads()`

- **Purpose:** Correct sequencing errors in short reads without assembling: BayesHammer for Illumina-type reads, IonHammer for Ion Torrent reads.
- **Method:** Runs `spades.py --only-error-correction` (which calls `spades-hammer`); `iontorrent=True` adds `--iontorrent` (`spades-ionhammer`). The corrected FASTQ files are written to `<output_dir>/corrected`.
- **API:** `dnakit.assembly.spades_correct_reads(output_dir, *, reads_1=None, reads_2=None, unpaired_reads=None, iontorrent=False, threads=None, memory_gb=None, phred_offset=None, spades_path=None)`.
- **Input:** A new `output_dir`, `reads_1` + `reads_2` (one pair of files) and/or `unpaired_reads`, and optional `iontorrent`, `threads`, `memory_gb`, `phred_offset` (33 or 64; give it when the qualities do not allow detection).
- **Output:** The common result keys plus `corrector` and `corrected_files`.
- **Example code:**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_correct_reads("corrected_reads", reads_1=["r1.fq"], reads_2=["r2.fq"], threads=2, phred_offset=33, verbose=False)
print(result["status"], result["corrector"], sorted(path.rsplit("/", 1)[-1] for path in result["corrected_files"] if path.endswith(".gz")))
```

- **Example output:**

```text
completed bayeshammer ['r100.0_0.cor.fastq.gz', 'r200.0_0.cor.fastq.gz', 'r_unpaired00.0_0.cor.fastq.gz']
```

## 5) `dnakit.assembly.spades_filter_reads_by_kmer_coverage()`

- **Purpose:** Remove reads whose median k-mer count is low.
- **Method:** Runs `spades-read-filter -k <kmer> -c <coverage_threshold>` on a dataset description written for you. A read pair is dropped when the median k-mer count of both reads is at or below the threshold.
- **API:** `dnakit.assembly.spades_filter_reads_by_kmer_coverage(output_dir, *, kmer, coverage_threshold, reads_1=None, reads_2=None, unpaired_reads=None, threads=None, drop_names=False, drop_quality=False, spades_path=None)`.
- **Input:** A new `output_dir`, `kmer`, `coverage_threshold`, reads (`reads_1` + `reads_2`, or `unpaired_reads`), and optional `threads`, `drop_names`, `drop_quality`.
- **Output:** The common result keys; `outputs` lists the filtered FASTQ files and `dataset.yaml`.
- **Example code:**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_filter_reads_by_kmer_coverage("filtered_reads", kmer=21, coverage_threshold=3, reads_1=["r1.fq"], reads_2=["r2.fq"], verbose=False)
print(result["status"], sorted(path.rsplit("/", 1)[-1] for path in result["outputs"]))
```

- **Example output:**

```text
completed ['1.1.fastq', '1.2.fastq', 'dataset.yaml']
```

## 6) `dnakit.assembly.spades_estimate_kmer_count()`

- **Purpose:** Estimate the number of distinct k-mers in a read set.
- **Method:** Runs `spades-kmer-estimating -k <kmer>` on a dataset description written for you; k-mers from reverse-complementary reads are not counted separately.
- **API:** `dnakit.assembly.spades_estimate_kmer_count(kmer, *, reads_1=None, reads_2=None, unpaired_reads=None, threads=None, spades_path=None)`.
- **Input:** `kmer`, reads (`reads_1` + `reads_2`, or `unpaired_reads`), and optional `threads`.
- **Output:** The common result keys plus `distinct_kmers`.
- **Example code:**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_estimate_kmer_count(21, reads_1=["r1.fq"], reads_2=["r2.fq"], verbose=False)
print(result["status"], result["distinct_kmers"])
```

- **Example output:**

```text
completed 38355
```

## 7) `dnakit.assembly.spades_build_assembly_graph()`

- **Purpose:** Build a de Bruijn assembly graph from reads or from a FASTA file.
- **Method:** Runs `spades-gbuilder <dataset or FASTA> <output> -k <kmer>` with `--unitigs`, `--gfa`, `--fastg` or `--spades` for the chosen `graph_format`; `infer_coverage` adds `-c`.
- **API:** `dnakit.assembly.spades_build_assembly_graph(output, *, kmer, reads_1=None, reads_2=None, unpaired_reads=None, fasta=None, graph_format='unitigs', infer_coverage=False, threads=None, spades_path=None)`.
- **Input:** A new `output` path, `kmer`, either reads (`reads_1` + `reads_2`, or `unpaired_reads`) or a `fasta`, and optional `graph_format` (`unitigs`, `gfa`, `fastg`, `spades`), `infer_coverage`, `threads`.
- **Output:** The common result keys; the graph is the output file.
- **Example code:**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

genome_graph = sp.spades_build_assembly_graph("genome.gfa", kmer=21, fasta="genome.fa", graph_format="gfa", verbose=False)
read_graph = sp.spades_build_assembly_graph("reads.gfa", kmer=21, reads_1=["r1.fq"], reads_2=["r2.fq"], graph_format="gfa", infer_coverage=True, verbose=False)
segments = [line for line in open("genome.gfa") if line.startswith("S\t")]
print(genome_graph["status"], read_graph["status"], len(segments) >= 3)
```

- **Example output:**

```text
completed completed True
```

## 8) `dnakit.assembly.spades_simplify_assembly_graph()`

- **Purpose:** Simplify an assembly graph by removing tips, bulges and low-coverage edges.
- **Method:** Runs `spades-gsimplifier <graph.gfa> <prefix> -k <kmer> --read-length <length>` with a scratch folder; the result is `<prefix>.gfa`.
- **API:** `dnakit.assembly.spades_simplify_assembly_graph(graph, output_prefix, *, kmer, read_length, coverage=None, use_coverage_ratios=False, threads=None, spades_path=None)`.
- **Input:** A GFA graph, a new `output_prefix`, `kmer`, `read_length`, and optional `coverage` (number or `auto`), `use_coverage_ratios`, `threads`.
- **Output:** The common result keys; `outputs` holds `<prefix>.gfa`.
- **Example code:**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_simplify_assembly_graph("genome.gfa", "simplified", kmer=21, read_length=100, verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **Example output:**

```text
completed ['simplified.gfa']
```

## 9) `dnakit.assembly.spades_split_assembly_graph()`

- **Purpose:** Split an assembly graph into one GFA file per connected component.
- **Method:** Runs `spades-gfa-split <graph.gfa> <output folder>`; the program creates the folder and writes `subgraph_<n>.gfa`.
- **API:** `dnakit.assembly.spades_split_assembly_graph(graph, output_dir, *, spades_path=None)`.
- **Input:** A GFA graph and a new `output_dir`.
- **Output:** The common result keys plus `components`; `outputs` lists the component files.
- **Example code:**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_split_assembly_graph("genome.gfa", "components", verbose=False)
print(result["status"], result["components"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **Example output:**

```text
completed 1 ['subgraph_0.gfa']
```

## 10) `dnakit.assembly.spades_map_reads_to_graph()`

- **Purpose:** Align reads to an assembly graph and write the paired-end distance histogram.
- **Method:** Runs `spades-gmapper <dataset> <graph.gfa> <output> -k <kmer>` with the library index, inner alignment mode (`-X`), retained-alignment policy (`-R`) and `--hic` when set.
- **API:** `dnakit.assembly.spades_map_reads_to_graph(graph, output, *, kmer, reads_1=None, reads_2=None, unpaired_reads=None, library_index=None, mode=None, retain=None, hic=False, threads=None, spades_path=None)`.
- **Input:** A GFA graph, a new `output` path, `kmer`, reads (`reads_1` + `reads_2`, or `unpaired_reads`), and optional `library_index`, `mode` (`illumina`, `pacbio`, `intractg`, `hic`), `retain` (`default`, `all`, `primary`, `quality`), `hic`, `threads`.
- **Output:** The common result keys; the distance histogram is the output file.
- **Example code:**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_map_reads_to_graph("genome.gfa", "distances.txt", kmer=21, reads_1=["r1.fq"], reads_2=["r2.fq"], verbose=False)
print(result["status"], "weight" in open("distances.txt").read())
```

- **Example output:**

```text
completed True
```

## 11) `dnakit.assembly.spaligner_align_long_reads()`

- **Purpose:** Align long reads or contigs to an assembly graph with SPAligner.
- **Method:** Runs `spaligner <config.yaml> -d <datatype> -s <sequences> -g <graph> -k <kmer> -o <output_dir>`. Without `config` the `spaligner_config.yaml` installed with SPAligner is used.
- **API:** `dnakit.assembly.spaligner_align_long_reads(sequences, graph, output_dir, *, datatype, kmer, config=None, threads=None, spades_path=None)`.
- **Input:** A FASTA/FASTQ of sequences, a GFA graph (or SPAdes saves folder), a new `output_dir`, `datatype` (`nanopore`, `pacbio`), `kmer` (odd), and optional `config`, `threads`.
- **Output:** The common result keys; `outputs` holds `alignment.tsv` (one row per aligned read with its graph path).
- **Example code:**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spaligner_align_long_reads("long.fa", "genome.gfa", "long_alignments", datatype="nanopore", kmer=21, verbose=False)
rows = open("long_alignments/alignment.tsv").read().splitlines()
print(result["status"], len(rows), [row.split("\t")[0] for row in rows])
```

- **Example output:**

```text
completed 4 ['l0', 'l1', 'l2', 'l3']
```

## 12) `dnakit.assembly.pathracer_search_graph()`

- **Purpose:** Find paths in an assembly graph that match a profile HMM or a sequence (PathRacer).
- **Method:** Runs `pathracer <query> <graph> --output <dir>` with `--hmm`, `--nt` or `--aa` for the query type and `--global` or `--local` for the matching mode; `minimum_length`, `indel_rate`, `top`, `seeding`, `evalue`, `threads` and `memory_gb` map to the PathRacer options.
- **API:** `dnakit.assembly.pathracer_search_graph(query, graph, output_dir, *, query_type='hmm', mode='global', minimum_length=None, indel_rate=None, top=None, threads=None, memory_gb=None, seeding=None, evalue=None, spades_path=None)`.
- **Input:** A query file (HMMER3 profile, or FASTA of nucleotide/amino-acid strings), a GFA graph or saves folder, a new `output_dir`, `query_type`, `mode`, and the optional parameters above.
- **Output:** The common result keys; `outputs` lists the matched sequences (`*.seqs.fa`), edges, annotated graph and log.
- **Example code:**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.pathracer_search_graph("query.fa", "genome.gfa", "pathracer_out", query_type="nt", threads=2, verbose=False)
print(result["status"], sorted(path.rsplit("/", 1)[-1] for path in result["outputs"]))
```

- **Example output:**

```text
completed ['graph_with_hmm_paths.gfa', 'pathracer.log', 'q.edges.fa', 'q.seqs.fa']
```

## 13) `dnakit.assembly.pathracer_align_hmm_to_sequences()`

- **Purpose:** Align an amino-acid profile HMM to nucleotide sequences while allowing frame shifts.
- **Method:** Runs `pathracer-seq-fs <hmm> <sequences> --output <dir>` with `--global` or `--local`, `--max-fs`, `--indel-rate`, `--top`, `--cutoff` and `--exhaustive` when set.
- **API:** `dnakit.assembly.pathracer_align_hmm_to_sequences(hmm, sequences, output_dir, *, mode='global', indel_rate=None, top=None, max_frameshifts=None, cutoff=None, exhaustive=False, threads=None, memory_gb=None, evalue=None, spades_path=None)`.
- **Input:** A HMMER3 amino-acid profile, a nucleotide FASTA, a new `output_dir`, and optional `mode`, `indel_rate`, `top`, `max_frameshifts`, `cutoff`, `exhaustive`, `threads`, `memory_gb`, `evalue`.
- **Output:** The common result keys; `outputs` lists the result files.
- **Example call (not run here):** This call needs a HMMER3 amino-acid profile of a real protein family, which cannot be generated here without HMMER, so it is not run in the documentation test; the unit test checks the argument list with a stand-in program.

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.pathracer_align_hmm_to_sequences("protein.hmm", "contigs.fa", "frameshift_out", max_frameshifts=4)
print(result["status"], result["outputs"])
```

## 14) `dnakit.assembly.splitter_resolve_repeats()`

- **Purpose:** Resolve repeats in an assembly graph with synthetic long reads or linked reads (SPlitteR).
- **Method:** Runs `splitter <graph> <library.yaml> <output_dir>` with the graph type (`-G`), resolution mode (`-M`), score thresholds, `--scaffold-links`, `--mapping-k` and barcode down-sampling options when set.
- **API:** `dnakit.assembly.splitter_resolve_repeats(graph, library_yaml, output_dir, *, library_index=None, graph_type=None, resolution_mode=None, score_threshold=None, relative_threshold=None, scaffold_links=False, mapping_kmer=None, sampling_factor=None, seed=None, assembly_info=None, threads=None, spades_path=None)`.
- **Input:** A GFA (or binary) graph, a SPAdes YAML description of the barcoded-read library, a new `output_dir`, and optional `library_index`, `graph_type`, `resolution_mode`, `score_threshold`, `relative_threshold`, `scaffold_links`, `mapping_kmer`, `sampling_factor`, `seed`, `assembly_info`, `threads`.
- **Output:** The common result keys; `outputs` lists the files of the resolved graph.
- **Example call (not run here):** This call needs synthetic-long-read (barcoded) data, so it is not run in the documentation test; the unit test checks the argument list with a stand-in program.

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.splitter_resolve_repeats("assembly_graph.gfa", "slr_library.yaml", "resolved", resolution_mode="meta", graph_type="mdbg")
print(result["status"])
```

## 15) `dnakit.assembly.binspreader_refine_bins()`

- **Purpose:** Refine and propagate metagenome bins over the assembly graph (BinSPreader).
- **Method:** Runs `binspreader <graph> <binning.tsv> <output.tsv>` with the contig paths, assignment strategy (`-Smax`/`-Smle`), refiner (`-Rcorr`/`-Rprop`), multiple assignment (`-m`), CAMI format, sparse propagation and read splitting options when set.
- **API:** `dnakit.assembly.binspreader_refine_bins(graph, binning, output, *, contig_paths=None, dataset_yaml=None, multiple_assignment=False, assignment_strategy=None, refiner=None, cami_format=False, zero_bin=False, labeled_alpha=None, sparse_propagation=False, split_reads=False, read_bin_weight_threshold=None, threads=None, spades_path=None)`.
- **Input:** A GFA (or binary) graph, a binner's contig-to-bin TSV, a new output path, and optional `contig_paths`, `dataset_yaml`, `multiple_assignment`, `assignment_strategy`, `refiner`, `cami_format`, `zero_bin`, `labeled_alpha`, `sparse_propagation`, `split_reads` (needs `dataset_yaml`), `read_bin_weight_threshold`, `threads`.
- **Output:** The common result keys; the refined binning is the output file.
- **Example call (not run here):** This call needs a metagenome assembly graph and a binner's output, so it is not run in the documentation test; the unit test checks the argument list with a stand-in program.

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.binspreader_refine_bins("assembly_graph.gfa", "bins.tsv", "refined.tsv", contig_paths="contigs.paths", assignment_strategy="mle", refiner="prop")
print(result["status"])
```

## 16) `dnakit.variants.whatshap_phasing_statistics()`

- **Purpose:** Compute phasing statistics of a phased VCF.
- **Method:** Runs `whatshap stats --tsv <file>` (and `--gtf`, `--block-list`, `--chr-lengths` when set) and parses the TSV.
- **API:** `dnakit.variants.whatshap_phasing_statistics(vcf, tsv, *, gtf=None, block_list=None, chromosome_lengths=None, sample=None, chromosome=None, only_snvs=False, whatshap_path=None)`.
- **Input:** A phased VCF, a new `tsv` path, and optional `gtf`, `block_list`, `chromosome_lengths`, `sample`, `chromosome`, `only_snvs`.
- **Output:** The common result keys plus `statistics` (one row per chromosome plus `ALL`: variants, phased, blocks, block N50 and so on).
- **Example code:**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_phasing_statistics("phased.vcf", "stats.tsv", gtf="blocks.gtf", block_list="blocks.txt", verbose=False)
row = result["statistics"][0]
print(result["status"], row["chromosome"], row["variants"], row["phased"], row["blocks"], row["block_n50"])
```

- **Example output:**

```text
completed chr1 25 25 1 3600
```

## 17) `dnakit.variants.whatshap_unphase_vcf()`

- **Purpose:** Remove phasing information from a VCF.
- **Method:** Runs `whatshap unphase`: phased genotypes such as `1|0` become `0/1` and the PS and HP fields are dropped.
- **API:** `dnakit.variants.whatshap_unphase_vcf(vcf, output, *, whatshap_path=None)`.
- **Input:** A VCF and a new output path.
- **Output:** The common result keys; the unphased VCF is the output file.
- **Example code:**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_unphase_vcf("phased.vcf", "unphased.vcf", verbose=False)
genotypes = {line.split("\t")[9].split(":")[0].strip() for line in open("unphased.vcf") if not line.startswith("#")}
print(result["status"], genotypes)
```

- **Example output:**

```text
completed {'0/1'}
```

## 18) `dnakit.variants.whatshap_compare_phasing()`

- **Purpose:** Compare two or more phased VCF files: switch errors, flips and Hamming distance.
- **Method:** Runs `whatshap compare` with `--tsv-pairwise`, `--tsv-multiway`, `--switch-error-bed` and/or `--longest-block-tsv`; at least one output is required and all outputs must be in one folder.
- **API:** `dnakit.variants.whatshap_compare_phasing(vcfs, *, tsv_pairwise=None, tsv_multiway=None, switch_error_bed=None, longest_block_tsv=None, names=None, ignore_sample_name=False, only_snvs=False, sample=None, ploidy=None, whatshap_path=None)`.
- **Input:** A list of at least two phased VCFs, output paths, and optional `names`, `ignore_sample_name`, `only_snvs`, `sample`, `ploidy`.
- **Output:** The common result keys plus `pairwise` (rows of the pairwise table) when `tsv_pairwise` is given.
- **Example code:**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_compare_phasing(["phased.vcf", "phased.vcf"], tsv_pairwise="compare.tsv", names=["a", "b"], verbose=False)
row = result["pairwise"][0]
print(result["status"], row["all_assessed_pairs"], row["all_switches"], row["blockwise_hamming"])
```

- **Example output:**

```text
completed 24 0 0
```

## 19) `dnakit.variants.whatshap_haplotag_reads()`

- **Purpose:** Tag each read with its haplotype (HP) and phase set (PS) using a phased VCF.
- **Method:** Runs `whatshap haplotag -o <bam>` and optionally `--output-haplotag-list`. The VCF must be bgzip-compressed and tabix-indexed.
- **API:** `dnakit.variants.whatshap_haplotag_reads(vcf, alignments, output, *, reference=None, haplotag_list=None, regions=None, ploidy=None, sample=None, ignore_read_groups=False, ignore_linked_reads=False, linked_read_distance_cutoff=None, skip_missing_contigs=False, tag_supplementary=None, threads=None, whatshap_path=None)`.
- **Input:** A bgzipped, indexed phased VCF, a BAM/CRAM, a new output path, and optional `reference`, `haplotag_list`, `regions`, `ploidy`, `sample`, `ignore_read_groups`, `ignore_linked_reads`, `linked_read_distance_cutoff`, `skip_missing_contigs`, `tag_supplementary`, `threads`.
- **Output:** The common result keys; the tagged BAM and the optional TSV (read, haplotype, phase set, chromosome) are written.
- **Example code:**

```python
import pysam

from dnakit.variants import whatshap_p14_cli_functions as wh

pysam.tabix_index("phased.vcf", preset="vcf", force=True, keep_original=True)
result = wh.whatshap_haplotag_reads("phased.vcf.gz", "reads.bam", "tagged.bam", reference="ref.fa", haplotag_list="haplotags.tsv", verbose=False)
tags = {}
for line in open("haplotags.tsv").read().splitlines()[1:]:
    name, haplotype = line.split("\t")[:2]
    tags[name] = haplotype
even = {tags[name] for name in tags if int(name[1:]) % 2 == 0}
odd = {tags[name] for name in tags if int(name[1:]) % 2 == 1}
print(result["status"], len(tags), len(even), len(odd), even != odd)
```

- **Example output:**

```text
completed 500 1 1 True
```

## 20) `dnakit.variants.whatshap_split_reads()`

- **Purpose:** Write the reads of each haplotype to its own file.
- **Method:** Runs `whatshap split <reads> <list>` with `--output-h1`, `--output-h2`, `--output-untagged` and `--read-lengths-histogram` as requested.
- **API:** `dnakit.variants.whatshap_split_reads(reads, haplotype_list, *, output_h1=None, output_h2=None, output_untagged=None, add_untagged=False, only_largest_block=False, discard_unknown_reads=False, read_lengths_histogram=None, whatshap_path=None)`.
- **Input:** A FASTQ/BAM of reads, a haplotype list (read name and H1/H2, for example the output of `whatshap_haplotag_reads`), at least one output path, and optional `add_untagged`, `only_largest_block`, `discard_unknown_reads`.
- **Output:** The common result keys; the per-haplotype files are written.
- **Example code:**

```python
import pysam

from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_split_reads("reads.bam", "haplotags.tsv", output_h1="h1.bam", output_h2="h2.bam", verbose=False)
counts = [sum(1 for _ in pysam.AlignmentFile(name)) for name in ("h1.bam", "h2.bam")]
print(result["status"], counts, sum(counts))
```

- **Example output:**

```text
completed [250, 250] 500
```

## 21) `dnakit.variants.whatshap_phase_from_haplotags()`

- **Purpose:** Add phase information to a VCF from haplotagged reads.
- **Method:** Runs `whatshap haplotagphase -o <vcf>`; the alignments must carry HP/PS tags and be indexed, and WhatsHap requires `reference`.
- **API:** `dnakit.variants.whatshap_phase_from_haplotags(vcf, alignments, output, *, reference=None, gap_threshold=None, cut_poly=None, only_indels=False, no_multi_allelic=False, ignore_read_groups=False, sample=None, chromosome=None, whatshap_path=None)`.
- **Input:** A VCF, an indexed haplotagged BAM/CRAM, a new output path, `reference`, and optional `gap_threshold`, `cut_poly`, `only_indels`, `no_multi_allelic`, `ignore_read_groups`, `sample`, `chromosome`.
- **Output:** The common result keys; the phased VCF is the output file.
- **Example code:**

```python
import pysam

from dnakit.variants import whatshap_p14_cli_functions as wh

pysam.index("tagged.bam")
result = wh.whatshap_phase_from_haplotags("unphased.vcf", "tagged.bam", "rephased.vcf", reference="ref.fa", verbose=False)
phased = sum("|" in line.split("\t")[9] for line in open("rephased.vcf") if not line.startswith("#"))
print(result["status"], phased)
```

- **Example output:**

```text
completed 25
```

## 22) `dnakit.variants.whatshap_genotype_variants()`

- **Purpose:** Genotype variants from reads, or from reads and a pedigree.
- **Method:** Runs `whatshap genotype`; with `pedigree` it adds `--ped` (with `--genmap` or `--recombrate`, and `--use-ped-samples`) so that the family is genotyped jointly.
- **API:** `dnakit.variants.whatshap_genotype_variants(vcf, alignments, output, *, reference=None, pedigree=None, genetic_map=None, recombination_rate=None, use_ped_samples=False, max_coverage=None, mapping_quality=None, genotype_quality_threshold=None, no_priors=False, prior_output=None, sample=None, chromosome=None, only_snvs=False, whatshap_path=None)`.
- **Input:** A VCF, a list of BAM/CRAM files, a new output path, and optional `reference`, `pedigree`, `genetic_map`, `recombination_rate`, `use_ped_samples`, `max_coverage`, `mapping_quality`, `genotype_quality_threshold`, `no_priors`, `prior_output`, `sample`, `chromosome`, `only_snvs`.
- **Output:** The common result keys; the genotyped VCF is the output file.
- **Example code:**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_genotype_variants("variants.vcf", ["reads.bam"], "genotyped.vcf", reference="ref.fa", verbose=False)
calls = {line.split("\t")[9].split(":")[0] for line in open("genotyped.vcf") if not line.startswith("#")}
print(result["status"], sorted(calls))
```

- **Example output:**

```text
completed ['1/0']
```

## 23) `dnakit.variants.whatshap_find_snv_candidates()`

- **Purpose:** Find candidate SNV sites from read alignments.
- **Method:** Runs `whatshap find_snv_candidates`; a site is a candidate when the alternative allele has at least `min_alt_depth` reads (default 3) and `min_alt_fraction` of the depth. `technology` applies the Illumina, PacBio or Nanopore defaults.
- **API:** `dnakit.variants.whatshap_find_snv_candidates(reference, alignment, output, *, technology=None, min_alt_depth=None, min_alt_fraction=None, multi_allelic=False, sample=None, chromosome=None, whatshap_path=None)`.
- **Input:** A reference FASTA, an indexed BAM, a new output path, and optional `technology`, `min_alt_depth`, `min_alt_fraction`, `multi_allelic`, `sample`, `chromosome`.
- **Output:** The common result keys; the candidate VCF is the output file.
- **Example code:**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_find_snv_candidates("ref.fa", "reads.bam", "candidates.vcf", technology="illumina", verbose=False)
found = {int(line.split("\t")[1]) for line in open("candidates.vcf") if not line.startswith("#")}
truth = {int(line.split("\t")[1]) for line in open("variants.vcf") if not line.startswith("#")}
print(result["status"], len(found), truth <= found)
```

- **Example output:**

```text
completed 25 True
```

## 24) `dnakit.variants.whatshap_learn_error_profile()`

- **Purpose:** Learn a sequencing-technology specific error profile from alignments and variants.
- **Method:** Runs `whatshap learn --reference <fasta> --output <file> <bam> <vcf>`; the output counts k-mer pairs (read versus reference) around the variants.
- **API:** `dnakit.variants.whatshap_learn_error_profile(alignment, vcf, reference, output, *, kmer_size=None, window=None, whatshap_path=None)`.
- **Input:** An indexed BAM, a VCF of variants, a reference FASTA, a new output path, and optional `kmer_size` (default 4), `window`.
- **Output:** The common result keys; the profile is the output file.
- **Example code:**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_learn_error_profile("reads.bam", "phased.vcf", "ref.fa", "error_profile.txt", verbose=False)
print(result["status"], len(open("error_profile.txt").read().splitlines()) > 0)
```

- **Example output:**

```text
completed True
```

## 25) `dnakit.variants.whatshap_polyphase()`

- **Purpose:** Phase polyploid variants from reads with the clustering and threading algorithm.
- **Method:** Runs `whatshap polyphase --ploidy <n>` with block-cut sensitivity, minimum overlap, prephasing and haploid-set options when set. A ploidy of 2 phases a diploid sample.
- **API:** `dnakit.variants.whatshap_polyphase(vcf, alignments, output, *, ploidy, reference=None, block_cut_sensitivity=None, min_overlap=None, use_prephasing=False, include_haploid_sets=False, no_multi_allelic=False, mapping_quality=None, threads=None, sample=None, chromosome=None, tag=None, only_snvs=False, distrust_genotypes=False, whatshap_path=None)`.
- **Input:** A VCF, a list of BAM/CRAM files, a new output path, `ploidy`, and optional `reference`, `block_cut_sensitivity` (0-5), `min_overlap`, `use_prephasing`, `include_haploid_sets`, `no_multi_allelic`, `mapping_quality`, `threads`, `sample`, `chromosome`, `tag`, `only_snvs`, `distrust_genotypes`.
- **Output:** The common result keys; the phased VCF is the output file.
- **Example code:**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_polyphase("variants.vcf", ["reads.bam"], "polyphased.vcf", ploidy=2, reference="ref.fa", verbose=False)
phased = sum("|" in line.split("\t")[9] for line in open("polyphased.vcf") if not line.startswith("#"))
print(result["status"], phased)
```

- **Example output:**

```text
completed 25
```

## 26) `dnakit.variants.whatshap_phase_with_pedigree()`

- **Purpose:** Phase variants with reads and a pedigree.
- **Method:** Runs `whatshap phase --ped <ped> -o <vcf>` with a genetic map or recombination rate, the recombination list, `--no-genetic-haplotyping` and `--use-ped-samples` when set.
- **API:** `dnakit.variants.whatshap_phase_with_pedigree(vcf, alignments, pedigree, output, *, reference=None, genetic_map=None, recombination_rate=None, recombination_list=None, no_genetic_haplotyping=False, use_ped_samples=False, algorithm=None, mapping_quality=None, sample=None, chromosome=None, tag=None, only_snvs=False, distrust_genotypes=False, whatshap_path=None)`.
- **Input:** A VCF and BAM/CRAM files that contain all family members, a PED/FAM file, a new output path, and optional `reference`, `genetic_map`, `recombination_rate`, `recombination_list`, `no_genetic_haplotyping`, `use_ped_samples`, `algorithm`, `mapping_quality`, `sample`, `chromosome`, `tag`, `only_snvs`, `distrust_genotypes`.
- **Output:** The common result keys; the phased VCF is the output file.
- **Example call (not run here):** This call needs family data (a trio VCF, one BAM per member and a PED file), so it is not run in the documentation test; the unit test checks the argument list with a stand-in program.

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_phase_with_pedigree("trio.vcf", ["child.bam", "father.bam", "mother.bam"], "trio.ped", "trio_phased.vcf", reference="ref.fa")
print(result["status"])
```

## 27) `dnakit.variants.whatshap_polyphase_with_progeny()`

- **Purpose:** Phase polyploid parents using the genotypes of their progeny.
- **Method:** Runs `whatshap polyphasegenetic --ploidy <n> <vcf> <pedigree>` with the optional progeny file, scoring window and complexity support.
- **API:** `dnakit.variants.whatshap_polyphase_with_progeny(vcf, pedigree, output, *, ploidy, progeny_file=None, scoring_window=None, complexity_support=None, sample=None, chromosome=None, tag=None, only_snvs=False, distrust_genotypes=False, whatshap_path=None)`.
- **Input:** A VCF, a pedigree file, a new output path, `ploidy`, and optional `progeny_file`, `scoring_window`, `complexity_support`, `sample`, `chromosome`, `tag`, `only_snvs`, `distrust_genotypes`.
- **Output:** The common result keys; the phased VCF is the output file.
- **Example call (not run here):** This call needs a polyploid family VCF and pedigree, so it is not run in the documentation test; the unit test checks the argument list with a stand-in program.

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_polyphase_with_progeny("tetraploid_family.vcf", "family.ped", "parents_phased.vcf", ploidy=4)
print(result["status"])
```
