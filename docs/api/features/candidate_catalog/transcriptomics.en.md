# Transcriptomics

This category contains **18** documented functions for bulk RNA-seq. Eleven compute a real result from the data passed in (10 are pure Python; `deconvolute` imports numpy and scipy when called and returns `status: "not_implemented"` if either is missing). Five wrap external programs (`align` runs STAR or HISAT2, `quantify` dispatches to the next three, `run_featurecounts`, `run_kallisto_quantification` and `run_salmon_quantification` run featureCounts, kallisto and salmon): they validate their inputs and return `status: "not_implemented"` while the needed program is not on PATH. Two functions, `predict_protein_function` and `predict_protein_localization`, always return `status: "not_implemented"` because no trained model is bundled. Count matrices are genes x samples. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError`. The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) align

- **Purpose:** Aligns RNA-seq reads to a genome by running STAR or HISAT2 (an external program).

- **Method:** The function does not align reads itself. It checks the inputs and then runs `STAR` (`--genomeDir <index> --readFilesIn ... --outSAMtype BAM SortedByCoordinate`, with `zcat` for `.gz` input) or `hisat2` (`-x <index> -U` or `-1/-2 ... -S`), whichever is found on PATH (`aligner: "auto"` prefers STAR). When neither program is on PATH it returns `status: "not_implemented"` instead of failing. The STAR/HISAT2 results were not run in this documentation (neither program is installed here).

- **API:** `dnakit.transcriptomics.align(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A FASTQ path (single-end) or a list of two FASTQ paths (paired-end); `.gz` is accepted for STAR. Options: `index` (required; STAR genome directory or HISAT2 index prefix), `output_prefix` (required; the file must not exist and its directory must exist), `aligner` (`"auto"` default, `"star"` or `"hisat2"`), `threads` (1), `timeout_seconds` (3600).

- **Output:** With STAR or HISAT2 on PATH: the envelope with `status: "ok"` and `result` holding `program`, `aligner`, `command`, `alignment_file` (STAR: `<prefix>Aligned.sortedByCoord.out.bam`; HISAT2: `<prefix>.sam`), `alignment_bytes` and `stderr_tail`. A tool that exits with an error raises `BackendExecutionError`, a timeout raises `BackendTimeoutError`. Paths must exist and must not start with `-`; output paths must not already exist. Without either program: `status: "not_implemented"` with `required_dependency` (`"STAR or hisat2 on PATH"`) and `message` (the program was not found on PATH), as in the example.

- **Example code:**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import align

Path("reads.fastq").write_text("@r1\nACGTACGT\n+\nIIIIIIII\n")
Path("genome_index").mkdir()
Path("out").mkdir()
result = align("reads.fastq", options={"index": "genome_index", "output_prefix": "out/sample_"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "align",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 11,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d9d1590a8b2978eb1aa73ef39296ead72a0a304c3b4f3d46e696cca37a3baae9"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "STAR or hisat2 on PATH",
  "message": "transcriptomics.align runs the external program STAR / hisat2, which was not found on PATH."
}
```

## 2) annotate

- **Purpose:** Annotates gene identifiers (Ensembl IDs or symbols) with coordinates, strand and biotype from a GTF file.

- **Method:** The GTF (plain or `.gz`, at most 2 GiB) is scanned for `gene` features only. An identifier matches `gene_id` (with or without the `.version` suffix) or `gene_name`; the first matching gene feature is used. Coordinates are 1-based inclusive, as in GTF; `biotype` is taken from `gene_biotype` or `gene_type`. Identifiers without a match are listed in `not_found` and get `null`. Pure Python.

- **API:** `dnakit.transcriptomics.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A list of gene IDs or symbols, text with one per line, a file path, or a mapping or count table keyed by gene. Options: `gtf` (required; path to an existing `.gtf` or `.gtf.gz` that contains `gene` features).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `gtf`, `n_requested`, `n_found`, `not_found` and `annotations` ({identifier: {`gene_id`, `gene_name`, `chromosome`, `start`, `end`, `strand`, `biotype`, `source`} or `null`}).

- **Example code:**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import annotate

gtf = """chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG1.5"; gene_name "TP53"; gene_biotype "protein_coding";
chr2\tsrc\tgene\t1000\t2000\t.\t-\t.\tgene_id "ENSG2.1"; gene_name "MYC"; gene_biotype "lncRNA";
"""
Path("genes.gtf").write_text(gtf)

result = annotate(["TP53", "ENSG2", "UNKNOWN"], options={"gtf": "genes.gtf"})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "transcriptomics",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 28,
    "sequence_bases": 13,
    "lengths": [
      2,
      4,
      7
    ],
    "sha256": "d35e660823da88a58853e330de3187d8719d657015cb617d8028c633735ff4fa"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "lookup of each identifier (gene_id with or without version suffix, or gene_name) among the 'gene' fe...",
    "gtf": "genes.gtf",
    "n_requested": 3,
    "n_found": 2,
    "not_found": [
      "UNKNOWN"
    ],
    "annotations": {
      "TP53": {
        "gene_id": "ENSG1.5",
        "gene_name": "TP53",
        "chromosome": "chr1",
        "start": 100,
        "end": 500,
        "strand": "+",
        "biotype": "protein_coding",
        "source": "src"
      },
      "ENSG2": {
        "gene_id": "ENSG2.1",
        "gene_name": "MYC",
        "chromosome": "chr2",
        "start": 1000,
        "end": 2000,
        "strand": "-",
        "biotype": "lncRNA",
        "source": "src"
      },
      "UNKNOWN": null
    }
  }
}
```

## 3) calculate_differential_expression

- **Purpose:** Tests differential gene expression between two sample groups (Welch t-test or Wilcoxon rank-sum test on log2 CPM).

- **Method:** Per gene: `norm = count / library size * target_sum` (with `normalization: "none"` the supplied values are used), the test runs on `log2(norm + pseudocount)`, `log2fc = log2(mean(norm, first group) + pseudocount) - log2(mean(norm, second group) + pseudocount)`, `pct_*` is the fraction of cells with raw count > 0, `p_value` comes from Welch's unequal-variance t-test (default) or a Mann-Whitney U / Wilcoxon rank-sum test (normal approximation with tie and continuity correction), and `q_value` is the Benjamini-Hochberg adjustment over all genes tested. Genes are listed by increasing p-value (ties: larger |log2fc| first). `group1` is the first group and `group2` the second; `target_sum` defaults to 1000000, so the values are CPM. The test is a simple per-gene t-test or rank-sum test, not a negative-binomial model such as DESeq2 or edgeR; Welch's test needs at least 2 samples in each group (otherwise `InputFormatError` with code `GROUP_TOO_SMALL` is raised; use `wilcoxon`). Unlike `singlecell_spatial.calculate_markers`, which skips too-small groups with a warning, this function raises. The comparison is `group1` versus `group2`. Pure Python.

- **API:** `dnakit.transcriptomics.calculate_differential_expression(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-sample count matrix (columns are samples): a mapping `{gene: [counts per sample]}` or `{gene: {sample: count}}`, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory. Options: `groups` (required; `{sample: condition}` or a list in column order), `group1` and `group2` (labels; optional when there are exactly two, then the sorted labels are used), `method` (`"welch"` default or `"wilcoxon"`), `normalization` (`"cpm"` default = scale each sample to `target_sum`, or `"none"` = test the supplied values), `target_sum` (1000000), `pseudocount` (1.0), `top_n` (100), `min_pct` (0.0; minimum fraction of samples with count > 0 in either group), `orientation` (`"genes_by_cells"` default = rows are genes, or `"cells_by_genes"` when rows are samples), `genes` and `cells` (name lists for unnamed forms; `cells` names the samples), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `group1`, `group2`, `n_group1`, `n_group2`, `method`, `normalization`, `target_sum`, `pseudocount`, `n_genes_tested`, `definition` and `genes` (top_n genes by increasing p-value; each gene row has `gene`, `mean_group1`, `mean_group2` (means of the normalised values), `log2fc`, `pct_group1`, `pct_group2`, `statistic`, `p_value` and `q_value`, where `pct_*` is the fraction of samples with count > 0).

- **Example code:**

```python
import json

from dnakit.transcriptomics import calculate_differential_expression

counts = """gene,ctrl1,ctrl2,ctrl3,treat1,treat2,treat3
G1,100,120,110,400,380,420
G2,50,55,45,52,48,50
G3,200,210,190,60,70,65
G4,10,12,8,11,9,10
G5,500,520,480,510,490,505
"""
groups = ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"]

result = calculate_differential_expression(counts, options={"groups": groups, "top_n": 3})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "transcriptomics",
  "operation": "calculate_differential_expression",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 162,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6bd60df05ff99d244ba3316a53a33fb300eb356cb150f138b797452f18d5f8ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "group1": "ctrl",
    "group2": "treat",
    "n_group1": 3,
    "n_group2": 3,
    "method": "welch",
    "normalization": "cpm",
    "target_sum": 1000000.0,
    "pseudocount": 1.0,
    "n_genes_tested": 5,
    "definition": "value = count / library size * target_sum; test on log2(value + pseudocount); log2fc = log2(mean(gro...",
    "genes": [
      {
        "gene": "G5",
        "mean_group1": 574897.4540996943,
        "mean_group2": 488711.4839501525,
        "log2fc": 0.23432119327432943,
        "pct_group1": 1.0,
        "pct_group2": 1.0,
        "statistic": 14.920024280439547,
        "p_value": 0.000125519935865829,
        "q_value": 0.000627599679329145
      },
      {
        "gene": "G1",
        "mean_group1": 126397.79860106656,
        "mean_group2": 389455.03823506663,
        "log2fc": -1.6234777351714555,
        "pct_group1": 1.0,
        "pct_group2": 1.0,
        "statistic": -25.979707994055687,
        "p_value": 0.000443874559514615,
        "q_value": 0.0011096863987865376
      },
      {
        "gene": "G3",
        "mean_group1": 229885.66987242262,
        "mean_group2": 63399.54882086602,
        "log2fc": 1.8583555787894248,
        "pct_group1": 1.0,
        "pct_group2": 1.0,
        "statistic": 23.047641606359313,
        "p_value": 0.0016789159197535598,
        "q_value": 0.002798193199589266
      }
    ]
  }
}
```

## 4) calculate_differential_splicing

- **Purpose:** Tests differential exon inclusion (percent spliced in, PSI) between two groups from junction read counts.

- **Method:** For each event, `PSI = (inc / inclusion_length) / (inc / inclusion_length + exc / exclusion_length)` per sample (samples without reads are ignored); `psi_group` is the mean over the samples of the group and `delta_psi = psi_group1 - psi_group2`. `p_value` is a two-sided Fisher exact test on the pooled inclusion/exclusion read counts of the two groups; it ignores variation between replicates and is therefore anti-conservative. `q_value` is the Benjamini-Hochberg adjustment over the tested events. Events with fewer than `min_reads` pooled reads in either group are listed but not tested (`delta_psi`, `p_value`, `q_value` are `null`). It is not rMATS or a beta-binomial model. Pure Python.

- **API:** `dnakit.transcriptomics.calculate_differential_splicing(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Inclusion-junction counts, events x samples: an inclusion-junction count matrix with one row per splicing event, in any of the matrix forms (a mapping `{event: [counts per sample]}`, a list of rows, CSV/TSV text or a file path, a numpy array or DataFrame, or a 10x directory); counts must be integers. Options: `exclusion` (required; exclusion-junction counts with identical events and samples in the same order), `groups` (required; `{sample: condition}` or a list), `group1` and `group2` (optional when there are exactly two labels), `inclusion_length` (1.0) and `exclusion_length` (1.0; effective lengths used to normalise the read counts), `min_reads` (10; pooled reads required in each group), `top_n` (100), `orientation` (`"genes_by_cells"` default = rows are genes, or `"cells_by_genes"` when rows are samples), `genes` and `cells` (name lists for unnamed forms; `cells` names the samples), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `group1`, `group2`, `n_events`, `n_tested`, `definition`, `min_reads` and `events` (top_n events by increasing p-value, untested ones last), each with `event`, `psi_group1`, `psi_group2`, `inclusion_group1`, `exclusion_group1`, `inclusion_group2`, `exclusion_group2` (pooled read counts), `delta_psi`, `p_value` and `q_value`.

- **Example code:**

```python
import json

from dnakit.transcriptomics import calculate_differential_splicing

inclusion = """event,a1,a2,b1,b2
exon1,80,90,20,25
exon2,50,55,52,48
"""
exclusion = """event,a1,a2,b1,b2
exon1,20,10,80,75
exon2,50,45,48,52
"""

result = calculate_differential_splicing(
    inclusion, options={"exclusion": exclusion, "groups": ["A", "A", "B", "B"]}
)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "transcriptomics",
  "operation": "calculate_differential_splicing",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 54,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "611a697ca02148bba4442fbcd806ec930e8f3c16da9b6d0a0aa3311d22f8b3c3"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "group1": "A",
    "group2": "B",
    "n_events": 2,
    "n_tested": 2,
    "definition": "PSI per sample = (inc/inclusion_length) / (inc/inclusion_length + exc/exclusion_length); psi_group =...",
    "min_reads": 10,
    "events": [
      {
        "event": "exon1",
        "psi_group1": 0.8500000000000001,
        "psi_group2": 0.225,
        "inclusion_group1": 170,
        "exclusion_group1": 30,
        "inclusion_group2": 45,
        "exclusion_group2": 155,
        "delta_psi": 0.6250000000000001,
        "p_value": 3.5477193042748607e-38,
        "q_value": 7.095438608549721e-38
      },
      {
        "event": "exon2",
        "psi_group1": 0.525,
        "psi_group2": 0.5,
        "inclusion_group1": 105,
        "exclusion_group1": 95,
        "inclusion_group2": 100,
        "exclusion_group2": 100,
        "delta_psi": 0.025000000000000022,
        "p_value": 0.6891172486827962,
        "q_value": 0.6891172486827962
      }
    ]
  }
}
```

## 5) calculate_translation_efficiency

- **Purpose:** Calculates translation efficiency (ribosome-profiling abundance divided by RNA-seq abundance) per gene.

- **Method:** For every shared gene and sample, `TE = ribo abundance / RNA abundance`, each abundance being the counts per million of its own library (RPKM when `gene_lengths` is given; the same length divides the ribo and the RNA abundance, so it cancels in the ratio: `gene_lengths` has no effect on the TE values with the default `pseudocount` 0 and only changes the unit named in the result `definition`). With `pseudocount` both abundances get it added after the length division (so with a pseudocount above 0 the lengths do slightly change TE), and with `log2` `log2(TE)` is reported. `TE` is `null` when the RNA abundance is zero. `mean_te` is the mean over the samples with a defined TE. Ribo-seq and RNA-seq samples are paired by identical names, otherwise by position when both have the same number of samples. Pure Python.

- **API:** `dnakit.transcriptomics.calculate_translation_efficiency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Ribosome-footprint counts, genes x samples: a ribosome-footprint count matrix in any of the matrix forms (a mapping `{gene: [counts per sample]}`, a list of rows, CSV/TSV text or a file path, a numpy array or DataFrame, or a 10x directory). Options: `rna` (required; RNA-seq count matrix of the same genes; samples are paired by identical names, otherwise by position and reported under the ribo-seq names), `gene_lengths` (`{gene: bp}`; only relabels the abundances RPKM instead of CPM and does not change TE, see Method), `log2` (false), `pseudocount` (0.0), `orientation` (`"genes_by_cells"` default = rows are genes, or `"cells_by_genes"` when rows are samples), `genes` and `cells` (name lists for unnamed forms; `cells` names the samples), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_genes`, `samples` and `genes` ({gene: {`te`: {sample: TE or `null`}, `mean_te`}}).

- **Example code:**

```python
import json

from dnakit.transcriptomics import calculate_translation_efficiency

ribo = """gene,s1,s2
g1,100,200
g2,300,400
"""
rna = """gene,s1,s2
g1,50,100
g2,600,800
"""

result = calculate_translation_efficiency(ribo, options={"rna": rna})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "calculate_translation_efficiency",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 33,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "18571d7b56efc4b7e5fbfa750d75cb89f9eec543617624660f0c7d8783b9dffd"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "TE = ribo abundance / RNA abundance, each as CPM (counts per million) of its own library; null when the RNA abundance is zero",
    "n_genes": 2,
    "samples": [
      "s1",
      "s2"
    ],
    "genes": {
      "g1": {
        "te": {
          "s1": 3.25,
          "s2": 3.0
        },
        "mean_te": 3.125
      },
      "g2": {
        "te": {
          "s1": 0.8125,
          "s2": 0.75
        },
        "mean_te": 0.78125
      }
    }
  }
}
```

## 6) compare

- **Purpose:** Compares expression profiles by Pearson or Spearman correlation across shared genes.

- **Method:** The correlation is computed between every sample of the input and every sample of `other` (or between all sample pairs of the input when `other` is not given) over the genes present in both matrices (at least 3). With `log2` (default true) the values are `log2(value + 1)` first; `normalization: "cpm"` scales each sample to 1000000 beforehand. A correlation is `null` when one of the profiles is constant. Pure Python.

- **API:** `dnakit.transcriptomics.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-sample count matrix (columns are samples): a mapping `{gene: [counts per sample]}` or `{gene: {sample: count}}`, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory. Options: `other` (optional second matrix; without it all sample pairs of the input are compared), `method` (`"pearson"` default or `"spearman"`), `log2` (true; correlate `log2(value + 1)`), `normalization` (`"none"` default, or `"cpm"`), `orientation` (`"genes_by_cells"` default = rows are genes, or `"cells_by_genes"` when rows are samples), `genes` and `cells` (name lists for unnamed forms; `cells` names the samples), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `method`, `n_shared_genes`, `mode` (`"all pairs within the input"` or `"input versus other"`) and `correlation` ({input sample: {other sample: r}}).

- **Example code:**

```python
import json

from dnakit.transcriptomics import compare

counts = """gene,s1,s2,s3
g1,1,2,10
g2,2,4,5
g3,3,6,1
g4,4,8,0
"""

result = compare(counts, options={"log2": False})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 51,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "0b35d186e0b0240699683fa5bfc3679b87c7da2fc9ae6d575703917756be5849"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "pearson correlation across 4 shared genes of the values",
    "method": "pearson",
    "n_shared_genes": 4,
    "mode": "all pairs within the input",
    "correlation": {
      "s1": {
        "s1": 1.0,
        "s2": 1.0,
        "s3": -0.9655351182200101
      },
      "s2": {
        "s1": 1.0,
        "s2": 1.0,
        "s3": -0.9655351182200101
      },
      "s3": {
        "s1": -0.9655351182200101,
        "s2": -0.9655351182200101,
        "s3": 1.0
      }
    }
  }
}
```

## 7) deconvolute

- **Purpose:** Estimates cell-type proportions of bulk RNA-seq samples by non-negative least squares on a signature matrix.

- **Method:** Each sample's profile is fitted as a non-negative combination of the signature columns with non-negative least squares (`scipy.optimize.nnls`) over the genes shared with the signature; the coefficients divided by their sum give the proportions, and `rmse` is the fit residual. With `normalize` (default true) every column is first scaled to sum 1. If `reference` and `reference_labels` are given instead of `signature`, the signature is the mean library-size-normalised profile per cell type. Requires numpy and scipy (returns `not_implemented` if either is missing). It is not CIBERSORT (no nu-SVR) or a model-based method.

- **API:** `dnakit.transcriptomics.deconvolute(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-sample count matrix (columns are samples): a mapping `{gene: [counts per sample]}` or `{gene: {sample: count}}`, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory. Options: `signature` (genes x cell types matrix in any matrix form) OR `reference` (count matrix of annotated cells) plus `reference_labels` (`{cell: type}` or a list), `normalize` (true), `orientation` (`"genes_by_cells"` default = rows are genes, or `"cells_by_genes"` when rows are samples), `genes` and `cells` (name lists for unnamed forms; `cells` names the samples), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `cell_types`, `n_shared_genes`, `fractions` ({sample: {cell type: proportion}}) and `rmse` ({sample: residual}).

- **Example code:**

```python
import json

from dnakit.transcriptomics import deconvolute

bulk = """gene,mix1,mix2
A,100,10
B,10,100
C,50,50
"""
signature = """gene,typeX,typeY
A,10,1
B,1,10
C,5,5
"""

result = deconvolute(bulk, options={"signature": signature})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "transcriptomics",
  "operation": "deconvolute",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 41,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "de454e1e3f3d2c243a5fd4e7192b152fcd71e2561fa478617bee4a934cc5c305"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "non-negative least squares of each mixture column on the signature columns (scipy.optimize.nnls); fr...",
    "cell_types": [
      "typeX",
      "typeY"
    ],
    "n_shared_genes": 3,
    "fractions": {
      "mix1": {
        "typeX": 1.0,
        "typeY": 4.636616585542447e-17
      },
      "mix2": {
        "typeX": 0.0,
        "typeY": 1.0
      }
    },
    "rmse": {
      "mix1": 3.0913372807800604e-17,
      "mix2": 2.1817011766850644e-17
    }
  }
}
```

## 8) detect

- **Purpose:** Detects expressed genes in RNA-seq samples with a read-count (or CPM) threshold.

- **Method:** A gene is detected in a sample when `count >= min_count`, or, when `min_cpm` is given, when `count / library size * 1e6 >= min_cpm` (which replaces the count rule). It counts as expressed overall when it is detected in at least `min_samples` samples. This is a plain threshold rule, not a statistical test. Pure Python.

- **API:** `dnakit.transcriptomics.detect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-sample count matrix (columns are samples): a mapping `{gene: [counts per sample]}` or `{gene: {sample: count}}`, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory. Options: `min_count` (1), `min_cpm` (`null`; if set, CPM >= `min_cpm` is used instead), `min_samples` (1; samples that must detect the gene), `orientation` (`"genes_by_cells"` default = rows are genes, or `"cells_by_genes"` when rows are samples), `genes` and `cells` (name lists for unnamed forms; `cells` names the samples), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `rule` (the rule applied, as text), `n_genes`, `n_detected`, `detected_genes` (at most 5000 listed) with `detected_truncated`, and `per_sample_detected` ({sample: number of genes detected}).

- **Example code:**

```python
import json

from dnakit.transcriptomics import detect

counts = """gene,s1,s2,s3
g1,0,0,0
g2,5,0,0
g3,5,7,9
"""

result = detect(counts, options={"min_count": 5, "min_samples": 2})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "detect",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 41,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e0dd490a7d446a9a0f50bf0dd5e73f4ae92c1b6c496f7e711eb6246ee72eb249"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "rule": "a gene is detected in a sample when count >= 5.0; it is expressed overall when detected in >= 2 samples",
    "n_genes": 3,
    "n_detected": 1,
    "detected_genes": [
      "g3"
    ],
    "detected_truncated": false,
    "per_sample_detected": {
      "s1": 2,
      "s2": 1,
      "s3": 1
    }
  }
}
```

## 9) normalize

- **Purpose:** Normalises an RNA-seq count matrix to CPM, TPM, FPKM, median-of-ratios or upper-quartile scale.

- **Method:** `cpm`: `count / library size * 1e6`. `fpkm`/`rpkm`: `count / (gene length in kb * library size in millions)`. `tpm`: `(count / length in kb) / sum over genes of (count / length in kb) * 1e6`. `median_of_ratios`: the DESeq2 size factor of a sample is the median, over genes with a positive count in every sample, of `count / geometric mean across samples`, and values are `count / size factor`. `upper_quartile`: the size factor is the 75th percentile of the sample's non-zero counts divided by the mean of these percentiles over samples, and values are `count / size factor`. TPM and FPKM need `gene_lengths` for every gene. With `log2` the result is `log2(value + pseudocount)`. A sample with zero total counts raises an error for cpm/fpkm/tpm. Pure Python.

- **API:** `dnakit.transcriptomics.normalize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-sample count matrix (columns are samples): a mapping `{gene: [counts per sample]}` or `{gene: {sample: count}}`, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory. Options: `method` (`"cpm"` default; `"tpm"`, `"fpkm"`/`"rpkm"`, `"median_of_ratios"`, `"upper_quartile"`), `gene_lengths` (`{gene: bp}`; required for tpm/fpkm), `log2` (false), `pseudocount` (1.0), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default = rows are genes, or `"cells_by_genes"` when rows are samples), `genes` and `cells` (name lists for unnamed forms; `cells` names the samples), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `method`, `definition`, `size_factors` ({sample: factor}; only for `median_of_ratios` and `upper_quartile`), `pseudocount` (only with `log2`), `n_genes`, `n_cells` (the number of samples), `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the normalised matrix).

- **Example code:**

```python
import json

from dnakit.transcriptomics import normalize

counts = """gene,s1,s2
g1,10,20
g2,90,180
"""

result = normalize(counts, options={"method": "cpm"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "normalize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 30,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7346eedbe68369b8e087ddc644b6a50133cb1b06e9fd39344d46c93c7998ca09"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "cpm",
    "definition": "CPM = count / library size * 1e6",
    "n_genes": 2,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "s1",
        "s2"
      ],
      "values": [
        [
          100000.0,
          100000.0
        ],
        [
          900000.0,
          900000.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 10) plot

- **Purpose:** Draws an RNA-seq expression heatmap or a sample-versus-sample scatter plot as SVG.

- **Method:** `kind: "heatmap"` draws the `max_genes` genes with the highest variance across samples; `kind: "scatter"` plots one point per gene for the samples `x` and `y` and reports their Pearson correlation `pearson_r`. With `log2` (default true) the plotted values are `log2(value + 1)`. The SVG is built as text without any plotting library; it is also written to `output_path` when given (the name must end in `.svg`; an existing file raises `ConfigurationError`).

- **API:** `dnakit.transcriptomics.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-sample count matrix (columns are samples): a mapping `{gene: [counts per sample]}` or `{gene: {sample: count}}`, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory. Options: `kind` (`"heatmap"` default, or `"scatter"`), `x` and `y` (sample names, required for `"scatter"`; points are genes), `log2` (true), `max_genes` (50; the most variable genes are shown in the heatmap), `output_path` (optional `.svg` file to write; never overwrites), `title`, `orientation` (`"genes_by_cells"` default = rows are genes, or `"cells_by_genes"` when rows are samples), `genes` and `cells` (name lists for unnamed forms; `cells` names the samples), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `format` (`"svg"`), `svg` (the SVG text), `output_path` (the `.svg` file written, else `null`), `kind`; for a heatmap also `genes_shown`; for a scatter plot also `x`, `y` and `pearson_r`.

- **Example code:**

```python
import json

from dnakit.transcriptomics import plot

counts = """gene,s1,s2,s3
g1,1,2,10
g2,2,4,5
g3,3,6,1
"""

result = plot(counts, options={"kind": "scatter", "x": "s1", "y": "s2"})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "transcriptomics",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "eab73a065bf739aaf7971512f189b5aca923aad5062f3d120371b0fc77bb69ae"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "svg",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"480\" viewBox=\"0 0 720 480\"><rect width=\"...",
    "kind": "scatter",
    "output_path": null,
    "x": "s1",
    "y": "s2",
    "pearson_r": 0.9997907600601286
  }
}
```

## 11) predict_protein_function

- **Purpose:** Intended to predict protein function; it is not implemented and always returns `not_implemented`.

- **Method:** No computation is performed. Protein function prediction needs a trained model or a functional annotation database (for example InterPro or DeepGOPlus weights) that is not bundled, so the function exists only to keep its name and call signature stable.

- **API:** `dnakit.transcriptomics.predict_protein_function(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any input is accepted (it is only summarised in `input_summary`); a protein sequence is the intended input. `options` is optional and unused.

- **Output:** Always the envelope with `status: "not_implemented"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`, `required_dependency` and `message`; there is no `result`.

- **Example code:**

```python
import json

from dnakit.transcriptomics import predict_protein_function

result = predict_protein_function("MKTAYIAKQR")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "predict_protein_function",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 10,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "0b2b6b8411d220f8648b67b5f839a4336470f4fc93d9604009a1cb3a093ae29e"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "trained protein-function model or annotation database",
  "message": "Protein function prediction needs a trained model or functional annotation database (for example InterPro/DeepGOPlus weights) that is not bundled."
}
```

## 12) predict_protein_localization

- **Purpose:** Intended to predict protein subcellular localization; it is not implemented and always returns `not_implemented`.

- **Method:** No computation is performed. Localisation prediction needs a trained predictor (for example DeepLoc or SignalP weights) that is not bundled, so the function exists only to keep its name and call signature stable.

- **API:** `dnakit.transcriptomics.predict_protein_localization(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any input is accepted (it is only summarised in `input_summary`); a protein sequence is the intended input. `options` is optional and unused.

- **Output:** Always the envelope with `status: "not_implemented"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`, `required_dependency` and `message`; there is no `result`.

- **Example code:**

```python
import json

from dnakit.transcriptomics import predict_protein_localization

result = predict_protein_localization("MKTAYIAKQR")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "predict_protein_localization",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 10,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "0b2b6b8411d220f8648b67b5f839a4336470f4fc93d9604009a1cb3a093ae29e"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "trained subcellular-localisation model",
  "message": "Protein localisation prediction needs a trained predictor (for example DeepLoc or SignalP weights) that is not bundled."
}
```

## 13) qc

- **Purpose:** Computes per-sample RNA-seq quality metrics: library size, detected genes, top-gene fractions and mitochondrial percentage.

- **Method:** Per sample: `library_size` (total counts), `genes_detected` (genes with count >= `min_count`), `fraction_zero_genes`, `top_gene_fraction` (share of the most abundant gene), `top10_genes_fraction` (share of the 10 most abundant genes) and `pct_mito` (counts of genes whose name starts with a `mito_prefix`, divided by the total, times 100). Samples failing a given threshold are listed with the names of the failed checks; thresholds that are not given are only reported. Pure Python.

- **API:** `dnakit.transcriptomics.qc(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-sample count matrix (columns are samples): a mapping `{gene: [counts per sample]}` or `{gene: {sample: count}}`, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory. Options: `min_count` (1; count threshold for "detected"), `min_library_size`, `min_genes`, `max_pct_mito` (all default `null` = report only), `mito_prefix` (`"MT-"` and `"mt-"`), `orientation` (`"genes_by_cells"` default = rows are genes, or `"cells_by_genes"` when rows are samples), `genes` and `cells` (name lists for unnamed forms; `cells` names the samples), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `n_samples`, `n_genes`, `min_count`, `n_mito_genes`, `median_library_size`, `samples` ({sample: metrics above plus `fails`}), `failing_samples` and `thresholds`.

- **Example code:**

```python
import json

from dnakit.transcriptomics import qc

counts = """gene,s1,s2
MT-ND1,5,0
A,45,50
B,50,50
"""

result = qc(counts, options={"max_pct_mito": 3})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "qc",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 38,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c3047d0b68648cb9963418c78c82ee1c83131f4a520b723a7153e31e7acf07bc"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_samples": 2,
    "n_genes": 3,
    "min_count": 1.0,
    "n_mito_genes": 1,
    "median_library_size": 100.0,
    "samples": {
      "s1": {
        "library_size": 100.0,
        "genes_detected": 3,
        "fraction_zero_genes": 0.0,
        "top_gene_fraction": 0.5,
        "top10_genes_fraction": 1.0,
        "pct_mito": 5.0,
        "fails": [
          "max_pct_mito"
        ]
      },
      "s2": {
        "library_size": 100.0,
        "genes_detected": 2,
        "fraction_zero_genes": 0.33333333333333337,
        "top_gene_fraction": 0.5,
        "top10_genes_fraction": 1.0,
        "pct_mito": 0.0,
        "fails": []
      }
    },
    "failing_samples": [
      "s1"
    ],
    "thresholds": {
      "min_library_size": null,
      "min_genes": null,
      "max_pct_mito": 3.0
    }
  }
}
```

## 14) quantify

- **Purpose:** Quantifies transcript or gene expression by running salmon or kallisto (FASTQ input) or featureCounts (BAM/SAM input).

- **Method:** A dispatcher over `run_salmon_quantification`, `run_kallisto_quantification` and `run_featurecounts`; it does not quantify reads itself. Inputs ending in `.bam` or `.sam` are sent to featureCounts, all other inputs are treated as FASTQ and sent to the first of salmon, kallisto found on PATH (`method: "auto"`). When no suitable program is on PATH it returns `status: "not_implemented"`. The external runs were not executed in this documentation (none of the programs is installed here).

- **API:** `dnakit.transcriptomics.quantify(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** FASTQ path(s) (one single-end, two paired-end) for salmon/kallisto, or BAM/SAM path(s) for featureCounts. Options: `method` (`"auto"` default; or `"salmon"`, `"kallisto"`, `"featurecounts"`, which must fit the input type), plus the options of the chosen runner (`run_salmon_quantification`, `run_kallisto_quantification` or `run_featurecounts`).

- **Output:** With a suitable program on PATH: `status: "ok"` and the `result` of the chosen runner plus `method`. Without one: `status: "not_implemented"` with `required_dependency` (for example `"salmon or kallisto on PATH"`) and `message`, as in the example. A `method` that does not fit the input type raises `ConfigurationError`.

- **Example code:**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import quantify

Path("reads.fastq").write_text("@r1\nACGTACGT\n+\nIIIIIIII\n")
result = quantify("reads.fastq", options={"index": "transcripts_index", "output_dir": "quant_out"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "quantify",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 11,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d9d1590a8b2978eb1aa73ef39296ead72a0a304c3b4f3d46e696cca37a3baae9"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "salmon or kallisto on PATH",
  "message": "transcriptomics.quantify runs the external program salmon / kallisto, which was not found on PATH."
}
```

## 15) read

- **Purpose:** Reads an expression count matrix (mapping, rows, CSV/TSV, array or 10x directory) and summarises genes, samples and sparsity.

- **Method:** The shared matrix parser accepts all listed input forms and gives the same genes x samples matrix for each. `sparsity = 1 - nonzero entries / (genes * samples)`. The matrix values are returned when the matrix has at most `max_return_elements` elements, otherwise they are omitted with a warning. Pure Python.

- **API:** `dnakit.transcriptomics.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-sample count matrix (columns are samples): a mapping `{gene: [counts per sample]}` or `{gene: {sample: count}}`, a list of rows, CSV/TSV text or a file path (header = sample names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory. Options: `max_return_elements` (10000), `orientation` (`"genes_by_cells"` default = rows are genes, or `"cells_by_genes"` when rows are samples), `genes` and `cells` (name lists for unnamed forms; `cells` names the samples), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `n_genes`, `n_cells` (the number of samples), `total_counts`, `nonzero_entries`, `sparsity`, `genes_preview`, `cells_preview` (the first 10 names), `orientation` and, when small enough, `matrix` (`orientation`, `genes`, `cells`, `values`).

- **Example code:**

```python
import json

from dnakit.transcriptomics import read

counts = """gene,s1,s2
g1,5,0
g2,3,8
"""

result = read(counts)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 25,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e9a9a4e08438f684e4fe883bbb5c626311d31f99954fd46f604ece73b69a4b4d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_genes": 2,
    "n_cells": 2,
    "total_counts": 16.0,
    "nonzero_entries": 3,
    "sparsity": 0.25,
    "genes_preview": [
      "g1",
      "g2"
    ],
    "cells_preview": [
      "s1",
      "s2"
    ],
    "orientation": "genes_by_cells",
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "s1",
        "s2"
      ],
      "values": [
        [
          5.0,
          0.0
        ],
        [
          3.0,
          8.0
        ]
      ]
    }
  }
}
```

## 16) run_featurecounts

- **Purpose:** Counts reads per gene by running Subread featureCounts on alignment files (an external program).

- **Method:** The function does not count reads itself. It checks the inputs and then runs `featureCounts -a <annotation> -o <output> -T <threads> -s <strand> -t <feature_type> -g <attribute>` (plus `-p --countReadPairs` when `paired`) on the alignment files, and summarises the output table. When featureCounts is not on PATH it returns `status: "not_implemented"`. The featureCounts run was not executed in this documentation (the program is not installed here).

- **API:** `dnakit.transcriptomics.run_featurecounts(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A BAM/SAM path or a list of paths. Options: `annotation` (required; existing GTF/GFF path), `output` (required; counts file to create, must not exist), `strand` (0 unstranded, 1 stranded, 2 reversely stranded; default 0), `feature_type` (`"exon"`), `attribute` (`"gene_id"`), `paired` (false; count read pairs), `threads` (1), `timeout_seconds` (3600).

- **Output:** With featureCounts on PATH: `status: "ok"` and `result` holding `program`, `command`, `counts_file`, `summary_file`, `n_genes`, `samples`, `assigned_counts` ({sample: summed assigned counts}) and `stderr_tail`. A tool that exits with an error raises `BackendExecutionError`, a timeout raises `BackendTimeoutError`. Paths must exist and must not start with `-`; output paths must not already exist. Without featureCounts: `status: "not_implemented"` with `required_dependency` (`"featureCounts (subread) on PATH"`) and `message`, as in the example.

- **Example code:**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import run_featurecounts

Path("sample.bam").write_bytes(b"BAM\x01")
Path("genes.gtf").write_text("")
result = run_featurecounts("sample.bam", options={"annotation": "genes.gtf", "output": "counts.txt"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "run_featurecounts",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 10,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "46082487f59dd44b1aabcca96ce194fd906f7cd4b90ecb2c50b79aa9e1cd0637"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "featureCounts (subread) on PATH",
  "message": "transcriptomics.run_featurecounts runs the external program featureCounts (subread), which was not found on PATH."
}
```

## 17) run_kallisto_quantification

- **Purpose:** Quantifies transcript abundance by running `kallisto quant` on FASTQ reads (an external program).

- **Method:** The function does not quantify reads itself. It checks the inputs and then runs `kallisto quant -i <index> -o <output_dir> -t <threads>` (single-end reads add `--single -l <fragment_length> -s <fragment_sd>`) and summarises `abundance.tsv`. When kallisto is not on PATH it returns `status: "not_implemented"`. The kallisto run was not executed in this documentation (the program is not installed here).

- **API:** `dnakit.transcriptomics.run_kallisto_quantification(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A FASTQ path (single-end) or a list of two FASTQ paths (paired-end). Options: `index` (required; existing kallisto index file), `output_dir` (required; must not exist or be empty, and its parent must exist), `fragment_length` and `fragment_sd` (required for single-end reads), `threads` (1), `timeout_seconds` (3600).

- **Output:** With kallisto on PATH: `status: "ok"` and `result` holding `program`, `command`, `abundance_file`, `n_transcripts`, `total_est_counts`, `total_tpm`, `top_transcripts` (the 20 transcripts with the highest TPM) and `stderr_tail`. A tool that exits with an error raises `BackendExecutionError`, a timeout raises `BackendTimeoutError`. Paths must exist and must not start with `-`; output paths must not already exist. Without kallisto: `status: "not_implemented"` with `required_dependency` (`"kallisto on PATH"`) and `message`, as in the example.

- **Example code:**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import run_kallisto_quantification

Path("reads.fastq").write_text("@r1\nACGTACGT\n+\nIIIIIIII\n")
Path("transcripts.idx").write_bytes(b"")
result = run_kallisto_quantification("reads.fastq", options={"index": "transcripts.idx", "output_dir": "kallisto_out", "fragment_length": 200, "fragment_sd": 20})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "run_kallisto_quantification",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 11,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d9d1590a8b2978eb1aa73ef39296ead72a0a304c3b4f3d46e696cca37a3baae9"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "kallisto on PATH",
  "message": "transcriptomics.run_kallisto_quantification runs the external program kallisto, which was not found on PATH."
}
```

## 18) run_salmon_quantification

- **Purpose:** Quantifies transcript abundance by running `salmon quant` on FASTQ reads (an external program).

- **Method:** The function does not quantify reads itself. It checks the inputs and then runs `salmon quant -i <index> -l <library_type> -o <output_dir> -p <threads> --validateMappings` with `-r` (single-end) or `-1/-2` (paired-end), and summarises `quant.sf`. When salmon is not on PATH it returns `status: "not_implemented"`. The salmon run was not executed in this documentation (the program is not installed here).

- **API:** `dnakit.transcriptomics.run_salmon_quantification(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A FASTQ path (single-end) or a list of two FASTQ paths (paired-end). Options: `index` (required; existing salmon index directory), `output_dir` (required; must not exist or be empty, and its parent must exist), `library_type` (`"A"` = automatic), `threads` (1), `timeout_seconds` (3600).

- **Output:** With salmon on PATH: `status: "ok"` and `result` holding `program`, `command`, `quant_file`, `n_transcripts`, `total_reads_assigned`, `total_tpm`, `top_transcripts` (the 20 transcripts with the highest TPM) and `stderr_tail`. A tool that exits with an error raises `BackendExecutionError`, a timeout raises `BackendTimeoutError`. Paths must exist and must not start with `-`; output paths must not already exist. Without salmon: `status: "not_implemented"` with `required_dependency` (`"salmon on PATH"`) and `message`, as in the example.

- **Example code:**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import run_salmon_quantification

Path("reads.fastq").write_text("@r1\nACGTACGT\n+\nIIIIIIII\n")
Path("salmon_index").mkdir()
result = run_salmon_quantification("reads.fastq", options={"index": "salmon_index", "output_dir": "salmon_out"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "transcriptomics",
  "operation": "run_salmon_quantification",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 11,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d9d1590a8b2978eb1aa73ef39296ead72a0a304c3b4f3d46e696cca37a3baae9"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "salmon on PATH",
  "message": "transcriptomics.run_salmon_quantification runs the external program salmon, which was not found on PATH."
}
```
