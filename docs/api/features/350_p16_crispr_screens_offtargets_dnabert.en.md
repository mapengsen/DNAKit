# CRISPR screen analysis, off-target search and DNABERT embeddings

`dnakit.design_crispr` gets facades for [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) (count, RRA test, MLE, pathway, plot), [DrugZ](https://github.com/hart-lab/drugz) and [Cas-OFFinder](https://github.com/snugel/cas-offinder); `dnakit.representations` gets `embed_dna_with_dnabert`, which reproduces the LucaOneTasks DNABERT-S / DNABERT-2 inference scripts (token matrix or `[CLS]` vector). The programs are not installed by DNAKit; put them on `PATH` or pass the `*_path` argument. DrugZ is a script: download `drugz.py` and pass it as `drugz_script`. Tested versions: MAGeCK 0.5.9.5, Cas-OFFinder 2.4.1 (CPU through pocl), DrugZ 1.1.0.2. The DNABERT function loads checkpoint code, so it needs `allow_remote_code=True` and a local checkpoint; it was compared with the LucaOneTasks scripts and gives identical matrices (maximum difference 0).

## 1) `dnakit.design_crispr.p16_crispr_screen_cli_functions.count_sgrnas_mageck()`

- **Purpose:** Count sgRNA reads from FASTQ (or aligned reads) against an sgRNA library (`mageck count`).
- **Method:** `library` is the sgRNA table (id, sequence, gene). `fastq` is one FASTQ file per sample (gzip allowed), in the order of `sample_labels`; `alignment_files` can be given instead of `fastq` for SAM/BAM input. `normalization` is none, median, total or control. Writes `<prefix>.count.txt`, `<prefix>.count_normalized.txt` and `<prefix>.countsummary.txt` (parsed into `summary`).
- **API:** `dnakit.design_crispr.p16_crispr_screen_cli_functions.count_sgrnas_mageck(library, output_prefix, *, fastq=None, sample_labels=None, alignment_files=None, normalization='median', trim_5=None, sgrna_length=None, reverse_complement=False, count_n=False, mageck_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `summary`, `count_table`.
- **Example code:**

```python
import random
from dnakit.design_crispr.p16_crispr_screen_cli_functions import count_sgrnas_mageck, rank_genes_mageck

rng = random.Random(2)
genes = [f"G{i}" for i in range(40)]
library = [(f"{g}_sg{j}", "".join(rng.choice("ACGT") for _ in range(20)), g) for g in genes for j in range(4)]
open("library.txt", "w").write("".join(f"{a}\t{b}\t{c}\n" for a, b, c in library))
for name, factor in (("c1", 1.0), ("c2", 1.0), ("t1", 0.1), ("t2", 0.1)):
    with open(name + ".fq", "w") as handle:
        n = 0
        for _sid, seq, gene in library:
            count = rng.randint(80, 150)
            for _ in range(int(count * factor) if gene in ("G1", "G2") else count):
                handle.write(f"@r{n}\nTT{seq}AA\n+\n{'I' * 24}\n")
                n += 1
labels = ["c1", "c2", "t1", "t2"]
counted = count_sgrnas_mageck("library.txt", "screen", fastq=[n + ".fq" for n in labels], sample_labels=labels, trim_5=2, verbose=False)
print([(row["Label"], row["Mapped"] == row["Reads"]) for row in counted["summary"]])
ranked = rank_genes_mageck(counted["count_table"], "rra", treatment=["t1", "t2"], control=["c1", "c2"], verbose=False)
for row in ranked["genes"][:2]:
    print(row["id"], row["neg|rank"], row["neg|fdr"])
```

- **Example output:**

```text
[('c1', True), ('c2', True), ('t1', True), ('t2', True)]
G2 1 0.0001
G1 2 0.0001
```

## 2) `dnakit.design_crispr.p16_crispr_screen_cli_functions.enrich_pathways_mageck()`

- **Purpose:** Pathway/gene-set enrichment of a ranked gene list (`mageck pathway`) with GSEA or RRA.
- **Method:** `gene_ranking` is a `mageck test`/`mle` gene summary (or, with `single_ranking`, a two-column gene ranking); `gene_sets` is a GMT file. Writes `<prefix>.pathway_summary.txt` (top rows in `pathways`).
- **API:** `dnakit.design_crispr.p16_crispr_screen_cli_functions.enrich_pathways_mageck(gene_ranking, gene_sets, output_prefix, *, method='gsea', single_ranking=False, ranking_column=None, sort_criteria=None, permutations=None, mageck_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `pathways`.

## 3) `dnakit.design_crispr.p16_crispr_screen_cli_functions.estimate_gene_essentiality_mageck()`

- **Purpose:** Maximum-likelihood gene essentiality with a design matrix (`mageck mle`, EM with permutation p-values).
- **Method:** `design_matrix` is the sample-by-condition table (first column sample label, `baseline` column of 1). The `cnv_*` parameters correct (`--cnv-norm` with `--cell-line`) or estimate (`--cnv-est`) copy-number effects as in :func:`rank_genes_mageck`. Writes `<prefix>.gene_summary.txt` and `<prefix>.sgrna_summary.txt`.
- **API:** `dnakit.design_crispr.p16_crispr_screen_cli_functions.estimate_gene_essentiality_mageck(count_table, design_matrix, output_prefix, *, beta_labels=None, include_samples=None, normalization='median', cnv_normalization=None, cell_line=None, cnv_estimation=None, control_sgrna=None, permutation_rounds=None, remove_outliers=False, threads=1, mageck_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `genes`, `gene_summary`.

## 4) `dnakit.design_crispr.p16_crispr_screen_cli_functions.plot_screen_results_mageck()`

- **Purpose:** Per-gene sgRNA count plots for the listed genes (`mageck plot`; MAGeCK runs R, so `Rscript` is required).
- **Method:** Writes `<prefix>_sgrnacounts.pdf` and the data table `<prefix>_summary.txt` when the program produces them; the files actually written are listed in `outputs`.
- **API:** `dnakit.design_crispr.p16_crispr_screen_cli_functions.plot_screen_results_mageck(count_table, gene_summary, output_prefix, *, genes, samples=None, normalization='median', mageck_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `genes`.

## 5) `dnakit.design_crispr.p16_crispr_screen_cli_functions.score_drug_gene_interactions_drugz()`

- **Purpose:** Score drug-gene interactions in a paired sgRNA count screen with DrugZ (`drugz.py`).
- **Method:** `count_table` is the tab-separated sgRNA read count table (first column `GENE_CLONE`, second the gene). With `unpaired=False` control and drug samples are compared pairwise in the order given; `unpaired=True` (`-unpaired`) compares the mean fold change of the drug against the control group. `drugz_script` is the downloaded `drugz.py` and `drugz_python` the interpreter that has pandas, scipy and six (default: the current one). The gene table is written to `output_path` (columns normZ, pval/rank/fdr for synergy and suppression), the guide fold changes to `fold_change_path` when given.
- **API:** `dnakit.design_crispr.p16_crispr_screen_cli_functions.score_drug_gene_interactions_drugz(count_table, output_path, *, control_samples, drug_samples, drugz_script, unpaired=False, fold_change_path=None, remove_genes=None, pseudocount=5, min_observations=1, half_window_size=500, drugz_python=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `unpaired`, `genes`.

## 6) `dnakit.design_crispr.p16_crispr_screen_cli_functions.search_offtargets_cas_offinder()`

- **Purpose:** Search a genome for CRISPR off-target sites tolerant to mismatches and the PAM (`cas-offinder`).
- **Method:** `genome` is a FASTA file or a folder of FASTA/2bit files. `pattern` is the IUPAC guide-plus-PAM pattern (default 20 nt + NRG, i.e. SpCas9 NGG/NAG); `guides` maps each guide (same length as `pattern`, PAM positions as N) to the maximum number of mismatches. `device` is `C` (CPU), `G` (GPU) or `A` (accelerators), optionally followed by device ids. Output rows: guide, chromosome, position (0-based), site, strand, mismatches (parsed into `sites`, first 1,000). The Cas-OFFinder input file is written next to the output (`<output_path>.input.txt`).
- **API:** `dnakit.design_crispr.p16_crispr_screen_cli_functions.search_offtargets_cas_offinder(genome, guides, output_path, *, pattern='NNNNNNNNNNNNNNNNNNNNNRG', device='C', cas_offinder_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `site_count`, `sites`, `input_file`.
- **Example code:**

```python
import random
from dnakit.design_crispr.p16_crispr_screen_cli_functions import search_offtargets_cas_offinder

rng = random.Random(6)
chrom = "".join(rng.choice("ACGT") for _ in range(6000))
guide = chrom[1000:1020]
chrom = chrom[:1020] + "TGG" + chrom[1023:]
open("target.fa", "w").write(">chrT\n" + chrom + "\n")
result = search_offtargets_cas_offinder("target.fa", {guide + "NNN": 2}, "hits.txt", pattern="N" * 20 + "NGG", verbose=False)
exact = [site for site in result["sites"] if site["mismatches"] == 0]
print(exact[0]["chrom"], exact[0]["position"], exact[0]["strand"], exact[0]["site"][-3:])
```

- **Example output:**

```text
chrT 1000 + TGG
```

## 7) `dnakit.design_crispr.p16_crispr_screen_cli_functions.rank_genes_mageck()`

- **Purpose:** Gene ranking of a count table by robust rank aggregation (`mageck test`), negative and positive selection.
- **Method:** `treatment`/`control` are sample labels (or indices) of the count table; `day0_label` compares to day 0 instead of `control`. `cnv_normalization` is a copy-number matrix file (tab-separated, first column `SYMBOL`, one column per cell line) and `cell_line` the column to use (`--cnv-norm`/`--cell-line`: copy-number-corrected sgRNA scores); `cnv_estimation` is a BED file of genome regions for `--cnv-est` (estimate copy-number profiles from the screen). Writes `<prefix>.gene_summary.txt` and `<prefix>.sgrna_summary.txt`; the top 20 gene rows are in `genes`.
- **API:** `dnakit.design_crispr.p16_crispr_screen_cli_functions.rank_genes_mageck(count_table, output_prefix, *, treatment, control=None, day0_label=None, paired=False, normalization='median', adjust_method='fdr', gene_lfc_method=None, cnv_normalization=None, cell_line=None, cnv_estimation=None, control_sgrna=None, remove_zero=None, write_normalized_counts=False, mageck_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `genes`, `gene_summary`, `sgrna_summary`.
- **See the example under:** `count_sgrnas_mageck()`

## 8) `dnakit.representations.lucaonetasks_p16_functions.embed_dna_with_dnabert()`

- **Purpose:** Embed DNA sequences with DNABERT-S or DNABERT-2 (LucaOneTasks `dnaberts` / `dnabert2` inference scripts).
- **Method:** `model` is `"dnabert-s"` (checkpoint `zhihan1996/DNABERT-S`) or `"dnabert2"` (`zhihan1996/DNABERT-2-117M`); `checkpoint_path` is a local checkpoint directory (a Hugging Face cache snapshot or a download); without it the model is looked up in the local Hugging Face cache only, nothing is downloaded. `embedding_type` is `"matrix"` (token matrix, `(tokens, 768)`; `include_special_tokens` keeps the `[CLS]` and `[SEP]` rows), `"vector"` (the `[CLS]` vector of the last layer) or `"mean"` (mean of every token row including `[CLS]` and `[SEP]`, the DNABERT-S paper protocol). `truncation_length` nucleotides are kept from the `truncate_side` (right keeps the start). The checkpoint ships Python code, so `allow_remote_code=True` is required. With `output_path` (a new `.npz` file) the arrays are saved under the sequence ids; `matrix` results are otherwise only summarised. Needs torch and transformers (`dnakit[neural]`).
- **API:** `dnakit.representations.lucaonetasks_p16_functions.embed_dna_with_dnabert(sequences, *, model='dnabert-s', checkpoint_path=None, embedding_type='matrix', include_special_tokens=False, truncation_length=4094, truncate_side='right', output_path=None, device='auto', allow_remote_code=False, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status` `status`, `model`, `embedding_dimension`, `sequences`, `embeddings`, `outputs`, `provenance`.
