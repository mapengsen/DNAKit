# Metagenome classification, AMR screening and microbial typing facades

This group of functions in `dnakit.microbiome` runs [Kaiju](https://github.com/bioinformatics-centre/kaiju), [Bracken](https://github.com/jenniferlu717/Bracken), [geNomad](https://github.com/apcamargo/genomad), [AMRFinderPlus](https://github.com/ncbi/amr), [ABRicate](https://github.com/tseemann/abricate), [mlst](https://github.com/tseemann/mlst), [fastANI](https://github.com/ParBLiSS/FastANI), [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper), [DIAMOND](https://github.com/bbuchfink/diamond) and two Prokka database converters. Each function runs one program with a fixed argument whitelist, without a shell, and returns a JSON-serialisable dictionary. The Bracken table merge and the fastANI plot are implemented in Python and need no program.

The programs are not installed by DNAKit. Put them on `PATH` or pass the `*_path` argument; a missing program raises `BackendUnavailableError` with an install hint. Outputs are written to new paths and never overwrite existing files. The examples below were run with Kaiju 1.10.3, DIAMOND 2.2.8, fastANI 1.34 and Bracken 3.1 on `PATH`. geNomad, eggNOG-mapper and ABRicate database download need large databases or the network and are not run in the examples; their command assembly is covered by unit tests with a recording stand-in program.

## Example data

The declared input files are written before the examples run; the other examples create their data in code.

**Example file:** `s1.bracken`

```text
name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
Species_one	100	S	30	0	30	0.5
Species_two	200	S	30	0	30	0.5
```

**Example file:** `s2.bracken`

```text
name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
Species_one	100	S	10	0	10	1.0
```

**Example file:** `nodes.dmp`

```text
1	|	1	|	no rank	|
10	|	1	|	superkingdom	|
60	|	10	|	genus	|
100	|	60	|	species	|
200	|	60	|	species	|
```

**Example file:** `names.dmp`

```text
1	|	root	|		|	scientific name	|
10	|	Bacteria	|		|	scientific name	|
60	|	Genus_A	|		|	scientific name	|
100	|	Species_one	|		|	scientific name	|
200	|	Species_two	|		|	scientific name	|
```

## 1) `dnakit.microbiome.p16_microbial_cli_functions.add_kaiju_taxon_names()`

- **Purpose:** Append taxon names to Kaiju per-read output (`kaiju-addTaxonNames`).
- **Method:** `drop_unclassified` (`-u`) removes unclassified reads, `full_path` (`-p`) prints the full lineage and `path_ranks` (`-r`) prints only the listed ranks.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.add_kaiju_taxon_names(kaiju_output, nodes, names, output_path, *, drop_unclassified=False, full_path=False, path_ranks=None, kaiju_names_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `line_count`.
- **See the example under:** `summarize_kaiju_table()`

## 2) `dnakit.microbiome.p16_microbial_cli_functions.annotate_dna_eggnog_mapper()`

- **Purpose:** Orthology annotation of DNA with eggNOG-mapper (`emapper.py`).
- **Method:** `input_type` is `CDS` (coding sequences, optionally `translate`), `genome` or `metagenome`; for the last two genes are predicted first with `gene_prediction` = `prodigal` (ab initio) or `search` (homology-based blastx search), the program default when unset. `decorate_gff` is `yes` (create a GFF from the predictions) or the path of an existing GFF to decorate. `data_directory` is the eggNOG data folder (`--data_dir`). Needs the eggNOG database (about 50 GB with DIAMOND data); outputs are `<prefix>.emapper.annotations` (parsed into `rows`).
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.annotate_dna_eggnog_mapper(fasta, output_directory, *, prefix='emapper', input_type='CDS', gene_prediction=None, translate=False, decorate_gff=None, data_directory=None, threads=1, emapper_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `input_type`, `files`, `rows`.

## 3) `dnakit.microbiome.p16_microbial_cli_functions.build_bracken_kmer_distribution()`

- **Purpose:** Build the read-length k-mer distribution of a Kraken database (`bracken-build`).
- **Method:** Runs the two Bracken stages: simulate read-length k-mer assignments from the library (`kmer2read_distr` -> `database<L>mers.kraken`) and aggregate them (`generate_kmer_distribution.py` -> `database<L>mers.kmer_distrib`). Both files are written inside the existing Kraken `database` folder, as the program does; existing files are not overwritten (the call is refused). `kraken_type` is kraken2, kraken or krakenuniq.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.build_bracken_kmer_distribution(database, *, read_length=100, kmer_length=35, threads=1, kraken_type='kraken2', bracken_build_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `read_length`; outputs `database<L>mers.kraken` and `.kmer_distrib`.

## 4) `dnakit.microbiome.p16_microbial_cli_functions.check_gff_annotation_amrfinder()`

- **Purpose:** Validate a gene-annotation GFF against nucleotide and protein FASTA (AMRFinderPlus `gff_check`).
- **Method:** `annotation_format` is one of bakta, genbank, microscope, patric, pgap, prodigal, prokka, pseudomonasdb, rast, standard. Exit status 0 means the three files are consistent.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.check_gff_annotation_amrfinder(gff, protein_fasta, nucleotide_fasta, *, annotation_format='genbank', gff_check_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `gff`.

## 5) `dnakit.microbiome.p16_microbial_cli_functions.check_nucleotide_fasta_amrfinder()`

- **Purpose:** Validate a nucleotide FASTA and report the sequence count and lengths (AMRFinderPlus `fasta_check`).
- **Method:** `fasta_check` prints the number of sequences followed by the cumulative sequence length after each record; these numbers are parsed into `sequence_count`, `cumulative_lengths` and the per-sequence `lengths`. A malformed FASTA makes the program exit with an error, which is raised.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.check_nucleotide_fasta_amrfinder(fasta, *, fasta_check_path=None, work_directory=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `sequence_count`, `lengths`, `cumulative_lengths`.

## 6) `dnakit.microbiome.p16_microbial_cli_functions.combine_bracken_tables()`

- **Purpose:** Merge per-sample Bracken output tables into one multi-sample table (Bracken `combine_bracken_outputs.py`).
- **Method:** Implemented in Python with the script's semantics. Inputs are Bracken *output* tables (not reports) estimated at the same level. Columns: `name`, `taxonomy_id`, `taxonomy_lvl` and for each sample `<sample>_num` and `<sample>_frac` (fraction of the sample's re-estimated reads, five decimals). `names` give the sample names (default: the file base names).
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.combine_bracken_tables(tables, output_path, *, names=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `samples`, `taxa`.
- **Example code:**

```python
from dnakit.microbiome.p16_microbial_cli_functions import combine_bracken_tables

result = combine_bracken_tables(["s1.bracken", "s2.bracken"], "merged.tsv", names=["A", "B"], verbose=False)
print(result["samples"], result["taxa"])
print(open("merged.tsv").read())
```

- **Example output:**

```text
['A', 'B'] 2
name	taxonomy_id	taxonomy_lvl	A_num	A_frac	B_num	B_frac
Species_one	100	S	30	0.50000	10	1.00000
Species_two	200	S	30	0.50000	0	0.00000
```

## 7) `dnakit.microbiome.p16_microbial_cli_functions.convert_kaiju_to_krona()`

- **Purpose:** Convert Kaiju output to a Krona text file (`kaiju2krona`).
- **Method:** `nodes`/`names` are the taxonomy `nodes.dmp`/`names.dmp`; `ranks` limits the taxon path to the listed ranks (`-l`); `include_unclassified` adds the count of unclassified reads (`-u`).
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.convert_kaiju_to_krona(kaiju_output, nodes, names, output_path, *, ranks=None, include_unclassified=False, kaiju2krona_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `line_count`.
- **See the example under:** `summarize_kaiju_table()`

## 8) `dnakit.microbiome.p16_microbial_cli_functions.detect_amr_genes_nucleotide()`

- **Purpose:** AMR, stress and virulence genes and point mutations in nucleotide contigs (`amrfinder --nucleotide`).
- **Method:** Without `organism` genes are reported (`amr_report`); with `organism` (see `amrfinder --list_organisms`) organism-specific point mutations are also detected (`dna_mutation`). `plus` adds the plus genes; `min_identity` and `min_coverage` are 0-1 (program defaults: curated threshold and 0.5). `database` is an AMRFinderPlus database folder (default: the installed one). The tab-separated report is written to `output_path` and parsed into `rows`.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.detect_amr_genes_nucleotide(contigs, output_path, *, organism=None, database=None, plus=False, min_identity=None, min_coverage=None, translation_table=11, name=None, threads=1, amrfinder_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `organism`, `rows` (report table).

## 9) `dnakit.microbiome.p16_microbial_cli_functions.download_abricate_database()`

- **Purpose:** Download and format an ABRicate reference database (`abricate-get_db`; network operation).
- **Method:** `database` is argannot, card, ncbi, ncbibetalactamase, plasmidfinder, resfinder or vfdb. The files are written under `database_directory` (default: the installation's `db` folder). Needs network access, BLAST+ and the program.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.download_abricate_database(database, *, database_directory=None, force=False, abricate_get_db_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `database`.

## 10) `dnakit.microbiome.p16_microbial_cli_functions.genbank_to_protein_fasta_prokka()`

- **Purpose:** GenBank CDS features to a protein FASTA reference (Prokka `prokka-genbank_to_fasta_db`).
- **Method:** Each CDS is translated (`/transl_table` or `genetic_code`) and written as `>id EC~~~gene~~~product`. `id_tag` selects the qualifier used as FASTA id (default: first of protein_id, locus_tag, db_xref).
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.genbank_to_protein_fasta_prokka(genbank_files, output_path, *, id_tag=None, include_hypothetical=False, include_pseudo=False, min_length=0, genetic_code=None, separator='~~~', blank='', prokka_script_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `proteins`.

## 11) `dnakit.microbiome.p16_microbial_cli_functions.map_fastani_conserved_regions()`

- **Purpose:** ANI of one genome pair plus the reciprocal conserved-region mapping records (`fastANI --visualize`).
- **Method:** The summary line is returned (`ani`, `matched_fragments`, `total_fragments`); the per-fragment mapping file is `<output_path>.visual` and holds query/reference names, positions, identity and mapping coordinates.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.map_fastani_conserved_regions(query, reference, output_path, *, kmer=16, fragment_length=3000, min_fraction=0.2, threads=1, fastani_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `ani`, `matched_fragments`, `total_fragments`, `visual_file`.
- **Example code:**

```python
import random
from dnakit.microbiome.p16_microbial_cli_functions import map_fastani_conserved_regions, plot_fastani_conserved_regions

rng = random.Random(3)
genome = "".join(rng.choice("ACGT") for _ in range(12000))
open("a.fa", "w").write(">a\n" + genome + "\n")
mutated = list(genome)
for position in range(500, 12000, 400):
    mutated[position] = "A" if mutated[position] != "A" else "C"
open("b.fa", "w").write(">b\n" + "".join(mutated) + "\n")
result = map_fastani_conserved_regions("b.fa", "a.fa", "ani.txt", verbose=False)
print(round(result["ani"], 2), result["matched_fragments"], result["total_fragments"])
plot = plot_fastani_conserved_regions(result["visual_file"], "ani.png", verbose=False)
print(plot["fragments"], open("ani.png", "rb").read(4) == b"\x89PNG")
```

- **Example output:**

```text
99.66 4 4
4 True
```

## 12) `dnakit.microbiome.p16_microbial_cli_functions.merge_kaiju_outputs()`

- **Purpose:** Merge the per-read classifications of two Kaiju/Kraken outputs (`kaiju-mergeOutputs`).
- **Method:** Both inputs must be sorted by read name (column 2). `conflict` resolves disagreements: `1` or `2` (take that file), `lca` (default) or `lowest`; `lca` and `lowest` need `nodes` (`nodes.dmp`). `use_score` prefers the better 4th-column score (`-s`).
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.merge_kaiju_outputs(first, second, output_path, *, conflict='lca', nodes=None, use_score=False, kaiju_merge_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `conflict`, `rows`.
- **See the example under:** `summarize_kaiju_table()`

## 13) `dnakit.microbiome.p16_microbial_cli_functions.plot_fastani_conserved_regions()`

- **Purpose:** Dot plot of conserved regions between two genomes from a fastANI `.visual` file.
- **Method:** The original R script (`scripts/visualize.R`) is not run; the same data (reference start/end against query start/end of every mapped fragment, coloured by identity) is drawn with matplotlib (`dnakit[bio-plots]`). Returns the PNG path and the fragment table (first 1,000 rows).
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.plot_fastani_conserved_regions(visual_file, image_path, *, min_identity=0.0, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `image`, `fragments`, `rows`.
- **See the example under:** `map_fastani_conserved_regions()`

## 14) `dnakit.microbiome.p16_microbial_cli_functions.reestimate_bracken_report()`

- **Purpose:** Bracken abundance table plus a Kraken-style report with the re-estimated read counts (`bracken -w`).
- **Method:** `kraken_report` is a Kraken 2 report, `database` the Kraken database holding `database<read_length>mers.kmer_distrib` (see :func:`build_bracken_kmer_distribution`). `level` is D, P, C, O, F, G or S; `threshold` is the minimum read count before re-estimation (Bracken default 0).
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.reestimate_bracken_report(kraken_report, database, output_table, output_report, *, read_length=100, level='S', threshold=0, bracken_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `rows` (name, taxonomy_id, new_est_reads, fraction_total_reads …).

## 15) `dnakit.microbiome.p16_microbial_cli_functions.run_genomad()`

- **Purpose:** Classify virus and plasmid contigs with geNomad, whole or one stage at a time.
- **Method:** `step` selects the sub-command: `end-to-end` (aggregated classification of the full pipeline), `annotate` (gene prediction and marker annotation), `find-proviruses`, `marker-classification`, `nn-classification`, `aggregated-classification`, `score-calibration` (batch-corrected scores) and `summary` (filtering and reports). Stages after `annotate` read the intermediate files that earlier stages wrote into the same `output_directory`; steps that need the geNomad database require `database` (the folder from `genomad download-database`). `output_directory` is created if it does not exist and reused across stages. `min_score` is `--min-score` (summary/end-to-end), `enable_score_calibration` is `--enable-score-calibration`, `splits` is `--splits`.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.run_genomad(contigs, output_directory, *, step='end-to-end', database=None, threads=1, cleanup=True, splits=None, min_score=None, enable_score_calibration=False, genomad_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `step`, `output_directory`, `files`.

## 16) `dnakit.microbiome.p16_microbial_cli_functions.screen_contigs_abricate()`

- **Purpose:** Screen contigs for AMR, virulence or plasmid genes by BLASTn (`abricate`).
- **Method:** `database` is an installed ABRicate database (argannot, card, ncbi, ncbibetalactamase, plasmidfinder, resfinder, vfdb or a custom one); `min_identity` and `min_coverage` are percentages (ABRicate defaults 75 and 0). The tab-separated report is written to `output_path` and parsed into `rows`.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.screen_contigs_abricate(contigs, output_path, *, database='resfinder', min_identity=75.0, min_coverage=0.0, data_directory=None, abricate_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `database`, `rows`.

## 17) `dnakit.microbiome.p16_microbial_cli_functions.search_dna_translated_diamond()`

- **Purpose:** Translated search of DNA queries against a protein database (`diamond blastx`).
- **Method:** `database` is a DIAMOND `.dmnd` file. `taxonomic_classification=True` writes per-read LCA taxonomy (output format 102; the database needs taxonomy data); otherwise `output_format` is a DIAMOND format number or name (`6` tabular default, `0`, `5`, `100`...). `frameshift_penalty` (`--frameshift`) and `long_reads` (`--long-reads`, exclusive with `frameshift_penalty`) make the alignment frameshift-aware; `exhaustive` (`--swipe`) searches all database sequences. `sensitivity` is one of fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensitive. The first 1,000 rows of tabular output are parsed.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.search_dna_translated_diamond(queries, database, output_path, *, output_format='6', taxonomic_classification=False, frameshift_penalty=None, long_reads=False, exhaustive=False, max_target_seqs=None, evalue=0.001, sensitivity=None, threads=1, diamond_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `output_format`, `rows`.
- **Example code:**

```python
import random, subprocess
codon = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
         "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}
rng = random.Random(5)
proteins = {f"p{i}_{taxon}": "".join(rng.choice(list(codon)) for _ in range(300)) for i, taxon in enumerate((100, 200))}
open("db.faa", "w").write("".join(f">{k}\n{v}\n" for k, v in proteins.items()))
reads = {f"read_{k}_{j}": "".join(codon[a] for a in v[j * 30:j * 30 + 50]) for k, v in proteins.items() for j in range(3)}
reads["read_random"] = "".join(rng.choice("ACGT") for _ in range(150))
open("reads.fa", "w").write("".join(f">{k}\n{v}\n" for k, v in reads.items()))
import subprocess
from dnakit.microbiome.p16_microbial_cli_functions import search_dna_translated_diamond

subprocess.run(["diamond", "makedb", "--in", "db.faa", "-d", "proteins", "--quiet"], check=True, capture_output=True)
result = search_dna_translated_diamond("reads.fa", "proteins.dmnd", "hits.tsv", exhaustive=True, verbose=False)
for row in result["rows"][:3]:
    print(row["column_1"], row["column_2"], row["column_3"])
```

- **Example output:**

```text
read_p0_100_0 p0_100 100
read_p0_100_1 p0_100 100
read_p0_100_2 p0_100 100
```

## 18) `dnakit.microbiome.p16_microbial_cli_functions.search_kaiju_best_match()`

- **Purpose:** Best-matching database protein for DNA reads without taxonomy (`kaijux`).
- **Method:** `database` is a Kaiju `.fmi` index of proteins; `reads` (and `reads_2` for pairs) are FASTA/FASTQ DNA reads. `mode` is `greedy` (uses `mismatches`, `min_score`, `max_evalue`) or `mem`. `seg_filter` is the SEG low-complexity filter (`-x`/`-X`). The output has classified (C/U), read name, score and the matched database identifiers per read; the first 1,000 rows are parsed into `rows`.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.search_kaiju_best_match(reads, database, output_path, *, reads_2=None, mode='greedy', mismatches=3, min_match_length=11, min_score=65, max_evalue=0.01, seg_filter=True, threads=1, kaijux_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `rows` (classified, read, score, matches), `classified`.
- **Example code:**

```python
import random, subprocess
codon = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
         "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}
rng = random.Random(5)
proteins = {f"p{i}_{taxon}": "".join(rng.choice(list(codon)) for _ in range(300)) for i, taxon in enumerate((100, 200))}
open("db.faa", "w").write("".join(f">{k}\n{v}\n" for k, v in proteins.items()))
reads = {f"read_{k}_{j}": "".join(codon[a] for a in v[j * 30:j * 30 + 50]) for k, v in proteins.items() for j in range(3)}
reads["read_random"] = "".join(rng.choice("ACGT") for _ in range(150))
open("reads.fa", "w").write("".join(f">{k}\n{v}\n" for k, v in reads.items()))
subprocess.run(["kaiju-mkbwt", "-n", "2", "-l", "1", "-a", "ACDEFGHIKLMNPQRSTVWY", "-o", "db", "db.faa"], check=True, capture_output=True)
subprocess.run(["kaiju-mkfmi", "db"], check=True, capture_output=True)
from dnakit.microbiome.p16_microbial_cli_functions import search_kaiju_best_match

result = search_kaiju_best_match("reads.fa", "db.fmi", "best.out", verbose=False)
print(result["classified"])
for row in result["rows"][:3] + result["rows"][-1:]:
    print(f'{row["classified"]} {row["read"]} {row.get("matches", "")}'.rstrip())
```

- **Example output:**

```text
6
C read_p0_100_0 p0_100,
C read_p0_100_1 p0_100,
C read_p0_100_2 p0_100,
U read_random
```

## 19) `dnakit.microbiome.p16_microbial_cli_functions.summarize_abricate_reports()`

- **Purpose:** Gene presence/absence matrix from several ABRicate reports (`abricate --summary`).
- **Method:** 
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.summarize_abricate_reports(reports, output_path, *, abricate_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `rows`.

## 20) `dnakit.microbiome.p16_microbial_cli_functions.summarize_kaiju_table()`

- **Purpose:** Per-sample classification summary at a taxonomic rank (`kaiju2table`).
- **Method:** `rank` is one of phylum, class, order, family, genus, species. Only one of `min_percent` (`-m`) and `min_reads` (`-c`) may be given. `expand_viruses` (`-e`), `exclude_unclassified` (`-u`), `full_path` (`-p`) and `path_ranks` (`-l`) follow the program. One row per taxon is parsed into `rows`.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.summarize_kaiju_table(kaiju_outputs, nodes, names, output_path, *, rank='species', min_percent=None, min_reads=None, expand_viruses=False, exclude_unclassified=False, full_path=False, path_ranks=None, kaiju2table_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `rank`, `rows` (file, percent, reads, taxon_id, taxon_name).
- **Example code:**

```python
import subprocess
from dnakit.microbiome.p16_microbial_cli_functions import add_kaiju_taxon_names, merge_kaiju_outputs, summarize_kaiju_table

subprocess.run(["kaiju", "-t", "nodes.dmp", "-f", "db.fmi", "-i", "reads.fa", "-o", "kaiju.out"], check=True, capture_output=True)
table = summarize_kaiju_table(["kaiju.out"], "nodes.dmp", "names.dmp", "species.tsv", rank="species", verbose=False)
for row in table["rows"]:
    if row["reads"]:
        print(row["taxon_name"], row["reads"])
add_kaiju_taxon_names("kaiju.out", "nodes.dmp", "names.dmp", "named.out", verbose=False)
print(open("named.out").read().splitlines()[0].split("\t")[-1])
merged = merge_kaiju_outputs("kaiju.out", "kaiju.out", "merged.out", nodes="nodes.dmp", verbose=False)
print(open("merged.out").read() == open("kaiju.out").read())
```

- **Example output:**

```text
Species_one 3
Species_two 3
unclassified 1
Species_one
True
```

## 21) `dnakit.microbiome.p16_microbial_cli_functions.translate_abricate_genes_prokka()`

- **Purpose:** Translate an ABRicate nucleotide gene FASTA into a Prokka protein reference (`prokka-abricate_to_fasta_db`).
- **Method:** Input headers must be ABRicate's `db~~~gene~~~accession~~~description`. The program forces a Met start, removes the stop codon and refuses sequences with internal stops or duplicate proteins.
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.translate_abricate_genes_prokka(abricate_sequences, output_path, *, prokka_script_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `proteins`.

## 22) `dnakit.microbiome.p16_microbial_cli_functions.type_sequence_type_mlst()`

- **Purpose:** Multi-locus sequence type of assembled contigs (`mlst`).
- **Method:** The scheme is detected automatically unless `scheme` forces one. `novel_alleles_path` (a new FASTA file) receives novel alleles (`--novel`). Output columns: file, scheme, ST, then the allele calls (`rows` holds the parsed first row as `fields`). Scoring defaults follow mlst (95 / 10 / 50).
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.type_sequence_type_mlst(contigs, output_path, *, scheme=None, novel_alleles_path=None, min_identity=95.0, min_coverage=10.0, min_score=50.0, threads=1, mlst_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `fields` (file, scheme, ST, alleles).

## 23) `dnakit.microbiome.p16_microbial_cli_functions.type_shiga_toxin_operons()`

- **Purpose:** Shiga-toxin (stx) operon type of nucleotide assemblies (`stxtyper`).
- **Method:** 
- **API:** `dnakit.microbiome.p16_microbial_cli_functions.type_shiga_toxin_operons(contigs, output_path, *, amrfinder_format=False, name=None, threads=1, stxtyper_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `rows` (stx_type, operon, identities).
