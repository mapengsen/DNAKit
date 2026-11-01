# MEGAHIT, GTDB-Tk and SPOA command-line facades

This group of functions runs three caller-installed programs through fixed argument whitelists, without a shell: [MEGAHIT](https://github.com/voutcn/megahit)'s `megahit_core` helper (read statistics, FASTG conversion, local assembly of contigs with reads), [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) (marker identification, alignment, placement-based classification, de novo tree, ANI to representatives, MSA trimming and export) and [SPOA](https://github.com/rvaser/spoa) (partial order alignment as an MSA or a graph). A missing program raises an error with the install command; output files and folders must not exist yet.

Tested versions: MEGAHIT 1.2.9, SPOA 4.1.6. GTDB-Tk and its reference package (tens of gigabytes) were not available in the test environment, so its seven functions are tested with a stand-in program that follows the GTDB-Tk command line and output file names; their example calls are not executed. The other `megahit_core` sub-programs (`count`, `read2sdbg`, `seq2sdbg`, `iterate`, `assemble`) exchange MEGAHIT-private binary files and are not exposed; `dnakit.assembly.external.assemble_megahit` runs the whole pipeline.

## Example data

The setup writes a 6 kb random genome, a 2 kb contig file in MEGAHIT's format, two split contigs and 3000 read pairs; `s.fa` is a four-sequence FASTA for SPOA.

**Example file:** `s.fa`

```text
>a
ACGTACGTACGTAAC
>b
ACGTACCTACGTAAC
>c
ACGTACGTACGAAC
>d
GTTACGTACGTACGT
```

- **Example code:**

```python
import random
from pathlib import Path

random.seed(5)
genome = "".join(random.choice("ACGT") for _ in range(6000))
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
Path("k21.contigs.fa").write_text(f">k21_0 flag=1 multi=10.0 len=2000\n{genome[:2000]}\n")
Path("split.fa").write_text(f">a\n{genome[:2500]}\n>b\n{genome[3000:5500]}\n")
with open("p1.fq", "w") as f1, open("p2.fq", "w") as f2:
    for i in range(3000):
        insert = random.randint(300, 400)
        start = random.randrange(0, len(genome) - insert)
        frag = genome[start:start + insert]
        f1.write(f"@r{i}/1\n{frag[:100]}\n+\n{'I' * 100}\n")
        f2.write(f"@r{i}/2\n{rc(frag)[:100]}\n+\n{'I' * 100}\n")
print(len(genome), Path("p1.fq").read_text().count("@r"))
```

- **Example output:**

```text
6000 3000
```

## 1) `dnakit.assembly.megahit_p15_cli_functions.megahit_read_statistics()`

- **Purpose:** Read count, total bases, longest, shortest and mean read length of a FASTA/FASTQ file.
- **Method:** Runs `megahit_core readstat` on the file (gzip accepted).
- **API:** `dnakit.assembly.megahit_p15_cli_functions.megahit_read_statistics(reads, stats_path, *, megahit_core_path=None)`.
- **Input:** A reads file and a new path for the table.
- **Output:** `read_count`, `total_bases`, `longest`, `shortest`, `mean_length`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.assembly.megahit_p15_cli_functions import megahit_read_statistics

result = megahit_read_statistics("p1.fq", str(out / "stats.txt"))
print(result["read_count"], result["total_bases"], result["longest"], result["shortest"], result["mean_length"])
```

- **Example output:**

```text
3000 300000 100 100 100.0
```

## 2) `dnakit.assembly.megahit_p15_cli_functions.megahit_contigs_to_fastg()`

- **Purpose:** Convert MEGAHIT contigs into a FASTG assembly graph.
- **Method:** Runs `megahit_core contig2fastg <k> <k*.contigs.fa>`; each contig yields a forward and a reverse node.
- **API:** `dnakit.assembly.megahit_p15_cli_functions.megahit_contigs_to_fastg(contigs, k, fastg_path, *, megahit_core_path=None)`.
- **Input:** An intermediate contig file with `flag=` and `multi=` headers, `k` and a new FASTG path.
- **Output:** `fastg_file`, `node_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.assembly.megahit_p15_cli_functions import megahit_contigs_to_fastg

result = megahit_contigs_to_fastg("k21.contigs.fa", 21, str(out / "graph.fastg"))
print(result["node_count"], (out / "graph.fastg").read_text().splitlines()[0])
```

- **Example output:**

```text
2 >NODE_1_length_2000_cov_10.0000_ID_1;
```

## 3) `dnakit.assembly.megahit_p15_cli_functions.megahit_local_assembly()`

- **Purpose:** Extend contigs or fill gaps by local assembly of the reads that map to their ends.
- **Method:** Converts the reads to MEGAHIT's binary library (`megahit_core buildlib`), then runs `megahit_core local`. Paired reads let the assembler bridge gaps.
- **API:** `dnakit.assembly.megahit_p15_cli_functions.megahit_local_assembly(contigs, output_directory, *, single_reads=(), paired_reads=(), interleaved_reads=(), kmin=11, kmax=41, step=6, seed_kmer=31, min_contig_len=200, min_mapping_len=75, sparsity=8, similarity=0.8, threads=1, megahit_core_path=None)`.
- **Input:** A contig file, a new output folder and `single_reads`, `paired_reads` or `interleaved_reads`.
- **Output:** `local_contig_file`, `local_contig_count`, `local_contig_bases`, `read_library`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.assembly.megahit_p15_cli_functions import megahit_local_assembly

result = megahit_local_assembly("split.fa", str(out / "local"), paired_reads=[["p1.fq", "p2.fq"]], min_contig_len=100)
print(result["input_contig_count"], result["local_contig_count"] >= 1, result["local_contig_bases"] > 0)
```

- **Example output:**

```text
2 True True
```

## 4) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_identify_markers()`

- **Purpose:** Identify the GTDB archaeal and bacterial marker genes in genomes.
- **Method:** Runs `gtdbtk identify` (Prodigal and HMMER). The reference package folder is given as `database` and set as `GTDBTK_DATA_PATH` for that process.
- **API:** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_identify_markers(out_dir, *, genome_dir=None, batchfile=None, database=None, extension='fna', prefix='gtdbtk', genes=False, write_single_copy_genes=False, cpus=1, force=False, gtdbtk_path=None)`.
- **Input:** `genome_dir` or `batchfile`, a new `out_dir`, `extension`.
- **Output:** `markers_summary` tables, `translation_table_summary`, `failed_genomes`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_identify_markers

result = gtdbtk_identify_markers("gtdb_identify", genome_dir="genomes", database="/data/gtdbtk/release226", cpus=8)
print(list(result["markers_summary"]))
```

## 5) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_align_markers()`

- **Purpose:** Align identified markers into a masked multiple sequence alignment.
- **Method:** Runs `gtdbtk align`; `skip_trimming`, `custom_msa_filters`, `taxa_filter` and the consensus options are passed through.
- **API:** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_align_markers(out_dir, identify_dir, *, database=None, prefix='gtdbtk', skip_gtdb_refs=False, taxa_filter=None, min_perc_aa=10.0, custom_msa_filters=False, skip_trimming=False, cols_per_gene=42, min_consensus=25.0, max_consensus=95.0, min_perc_taxa=50.0, rnd_seed=None, cpus=1, gtdbtk_path=None)`.
- **Input:** The `identify` output folder and a new `out_dir`.
- **Output:** `files` (alignment files), `filtered_genomes`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_align_markers

result = gtdbtk_align_markers("gtdb_align", "gtdb_identify", database="/data/gtdbtk/release226")
print(result["files"]["*.user_msa.fasta*"])
```

## 6) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_classify_placement()`

- **Purpose:** Classify genomes by pplacer placement of a precomputed marker alignment.
- **Method:** Runs `gtdbtk classify` with the `align` output; taxonomy comes from placement, relative evolutionary divergence and ANI.
- **API:** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_classify_placement(out_dir, align_dir, *, genome_dir=None, batchfile=None, database=None, extension='fna', prefix='gtdbtk', full_tree=False, place_species=False, min_af=None, skip_pplacer=False, pplacer_cpus=None, cpus=1, gtdbtk_path=None)`.
- **Input:** The `align` folder, the genomes and a new `out_dir`.
- **Output:** `summary` (user genome, classification, closest reference, ANI, warnings), `trees`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_classify_placement

result = gtdbtk_classify_placement("gtdb_classify", "gtdb_align", genome_dir="genomes", database="/data/gtdbtk/release226")
print(result["summary"])
```

## 7) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_de_novo_workflow()`

- **Purpose:** Infer a de novo genome tree rooted and decorated with GTDB taxonomy.
- **Method:** Runs `gtdbtk de_novo_wf` (identify, align, infer, root, decorate) for `bacteria` or `archaea` with an outgroup taxon.
- **API:** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_de_novo_workflow(out_dir, outgroup_taxon, *, domain='bacteria', genome_dir=None, batchfile=None, database=None, extension='fna', prefix='gtdbtk', skip_gtdb_refs=False, taxa_filter=None, min_perc_aa=10.0, custom_msa_filters=False, cols_per_gene=42, min_consensus=25.0, max_consensus=95.0, min_perc_taxa=50.0, prot_model='WAG', no_support=False, gamma=False, rnd_seed=None, genes=False, write_single_copy_genes=False, keep_intermediates=False, force=False, cpus=1, gtdbtk_path=None)`.
- **Input:** Genomes, `outgroup_taxon`, `domain` and a new `out_dir`.
- **Output:** `trees` (decorated and rooted Newick files), `files`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_de_novo_workflow

result = gtdbtk_de_novo_workflow("gtdb_tree", "p__Patescibacteria", genome_dir="genomes", database="/data/gtdbtk/release226", prot_model="LG")
print(result["trees"]["*.decorated.tree"])
```

## 8) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_ani_to_representatives()`

- **Purpose:** Average nucleotide identity of genomes to the GTDB representative genomes.
- **Method:** Runs `gtdbtk ani_rep` (skani against all species representatives).
- **API:** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_ani_to_representatives(out_dir, *, genome_dir=None, batchfile=None, database=None, extension='fna', prefix='gtdbtk', min_af=None, cpus=1, gtdbtk_path=None)`.
- **Input:** Genomes and a new `out_dir`; optional `min_af`.
- **Output:** `ani_summary` rows: closest representative, ANI, alignment fraction, taxonomy.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_ani_to_representatives

result = gtdbtk_ani_to_representatives("gtdb_ani", genome_dir="genomes", database="/data/gtdbtk/release226")
print(result["ani_summary"][0])
```

## 9) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_trim_msa()`

- **Purpose:** Trim an untrimmed marker alignment with a column mask.
- **Method:** Runs `gtdbtk trim_msa`; columns where the mask has `1` are kept. `reference_mask` (`bac` or `arc`) uses the reference package.
- **API:** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_trim_msa(untrimmed_msa, output_path, *, mask_file=None, reference_mask=None, database=None, gtdbtk_path=None)`.
- **Input:** An alignment, a new output path and `mask_file` or `reference_mask`.
- **Output:** `output`, `sequence_count`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_trim_msa

result = gtdbtk_trim_msa("untrimmed.fasta", "trimmed.fasta", reference_mask="bac", database="/data/gtdbtk/release226")
print(result["sequence_count"])
```

## 10) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_export_reference_msa()`

- **Purpose:** Export the untrimmed GTDB reference marker alignment of one domain.
- **Method:** Runs `gtdbtk export_msa --domain bac|arc`.
- **API:** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_export_reference_msa(output_path, *, domain='bac', database=None, gtdbtk_path=None)`.
- **Input:** `domain` and a new output path.
- **Output:** `output`, `sequence_count`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_export_reference_msa

result = gtdbtk_export_reference_msa("bac120_untrimmed.fasta", domain="bac", database="/data/gtdbtk/release226")
print(result["sequence_count"])
```

## 11) `dnakit.alignment.spoa_p15_cli_functions.spoa_partial_order_alignment()`

- **Purpose:** Align sequences into a partial order graph and write the MSA, the consensus, or the graph in GFA or DOT format, optionally picking the better strand.
- **Method:** Runs `spoa` with `-r` (`result`: `msa`, `consensus_and_msa`, `consensus`, `graph`, `consensus_and_graph`), `-d` (`dot_path`), `-s` (`strand_ambiguous`), `-l` (`algorithm`) and the score options. Defaults are SPOA's, so the default alignment is local.
- **API:** `dnakit.alignment.spoa_p15_cli_functions.spoa_partial_order_alignment(sequences, output_path, *, result='msa', algorithm='local', match=5, mismatch=-4, gap_open=-8, gap_extend=-6, gap_open_second=-10, gap_extend_second=-4, strand_ambiguous=False, min_coverage=None, dot_path=None, spoa_path=None)`.
- **Input:** A FASTA/FASTQ file, a new output path and the options above.
- **Output:** `alignment` records, `consensus`, `alignment_length`; for graphs `node_count`, `edge_count`, `path_count`; `dot_file`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.spoa_p15_cli_functions import spoa_partial_order_alignment

msa = spoa_partial_order_alignment("s.fa", str(out / "msa.fa"), result="consensus_and_msa", algorithm="global", strand_ambiguous=True, dot_path=str(out / "g.dot"))
print(msa["consensus"], msa["alignment_length"])
for record in msa["alignment"]:
    print(record["name"], record["sequence"])
graph = spoa_partial_order_alignment("s.fa", str(out / "g.gfa"), result="consensus_and_graph", algorithm="global")
print(graph["node_count"], graph["edge_count"], graph["path_count"])
```

- **Example output:**

```text
ACGTACGTACGTAAC 15
a ACGTACGTACGTAAC
b ACGTACCTACGTAAC
c ACGTACGTACG-AAC
d ACGTACGTACGTAAC
19 20 5
```
