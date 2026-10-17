# Single-cell and spatial genomics

This category contains **54** documented functions for count matrices of single cells and spatial spots. Fifty-two compute a real result from the data passed in: 43 are pure Python, 7 import numpy when called (`cluster`, `cluster_spatial`, `infer_trajectory`, `infer_trajectory_spatial`, `map_cells`, `map_cells_spatial`, `map_single_cell_to_reference`) and 2 import numpy and scipy (`deconvolute_spots`, `deconvolute_spots_spatial`); if a needed package is missing they return `status: "not_implemented"` instead. Two functions, `remove_doublets` and `remove_doublets_spatial`, always return `status: "not_implemented"` because no doublet-detection model is bundled. Count matrices are genes x cells (columns are cells or spots); the `_spatial` variants treat spots as cells and take `coordinates` where the method needs them. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError`. The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) analyze_crispr_perturbations

- **Purpose:** Summarises a CRISPR perturbation screen: the number of cells per guide and the knockdown of each guide's target gene against control cells.

- **Method:** Expression is library-size normalised (`count / library size * target_sum`, default 10000). For every label other than the control label, the target gene (given by `targets`, by default the guide label itself) is compared between the guide's cells and the control cells: `log2fc = log2(mean_guide + 1) - log2(mean_control + 1)` and `percent_knockdown = 100 * (1 - mean_guide / mean_control)` (`null` when the control mean is 0). There is no significance test. Pure Python.

- **API:** `dnakit.singlecell_spatial.analyze_crispr_perturbations(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `perturbations` (required; `{cell: guide label}` or a list in column order), `control` (`"control"`; label of control cells), `targets` (`{guide: target gene}`; default the guide label is the gene), `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `control`, `n_control_cells`, `cells_per_label`, `target_sum`, `knockdown_definition` and `guides`, a mapping guide -> `n_cells`, `target_gene`, `target_in_matrix` and (when the gene is in the matrix) `mean_normalized_guide`, `mean_normalized_control`, `log2fc`, `percent_knockdown`. A target gene missing from the matrix gives `target_in_matrix: false` and a warning.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import analyze_crispr_perturbations

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
perturbations = {"c1": "control", "c2": "control", "c3": "control", "c4": "CD3D", "c5": "CD3D", "c6": "CD3D"}

result = analyze_crispr_perturbations(counts, options={"perturbations": perturbations})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "analyze_crispr_perturbations",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "control": "control",
    "n_control_cells": 3,
    "cells_per_label": {
      "CD3D": 3,
      "control": 3
    },
    "target_sum": 10000.0,
    "knockdown_definition": "100 * (1 - mean CP10k(guide cells) / mean CP10k(control cells)) of the target gene",
    "guides": {
      "CD3D": {
        "n_cells": 3,
        "target_gene": "CD3D",
        "target_in_matrix": true,
        "mean_normalized_guide": 61.7283950617284,
        "mean_normalized_control": 1363.1378477316855,
        "log2fc": -4.442726976505825,
        "percent_knockdown": 95.47159554226693
      }
    }
  }
}
```

## 2) annotate_cells

- **Purpose:** Annotates cells with cell-type labels by scoring user-supplied marker genes (no trained model).

- **Method:** Values are ln(1 + CP10k). For each marker gene the values are z-scored across all cells (population standard deviation; a constant gene scores 0). The score of a cell type in a cell is the mean z-score of that type's marker genes, and the cell gets the type with the highest score, or `"Unassigned"` when the best score is not above `min_score`. Marker genes missing from the matrix give a warning (an error if none of a type's markers is present). There is no trained model. Pure Python.

- **API:** `dnakit.singlecell_spatial.annotate_cells(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `markers` (required; `{cell type: [marker genes]}`), `min_score` (0.0), `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `assignments` ({cell: type}), `best_score` ({cell: score}), `counts` (cells per type) and `scores` ({type: {cell: score}}).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import annotate_cells

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
markers = {"T cell": ["CD3D", "CD3E"], "B cell": ["CD79A", "MS4A1"]}

result = annotate_cells(counts, options={"markers": markers})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "annotate_cells",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score(type, cell) = mean over the type's marker genes of the per-gene z-score of ln(1 + CP10k); the highest score wins if > min_score",
    "assignments": {
      "c1": "T cell",
      "c2": "T cell",
      "c3": "T cell",
      "c4": "B cell",
      "c5": "B cell",
      "c6": "B cell"
    },
    "best_score": {
      "c1": 0.8356997436691734,
      "c2": 0.8567110675970135,
      "c3": 0.8259936626183231,
      "c4": 0.8350643975343444,
      "c5": 0.8410178115070588,
      "c6": 0.8413837216201503
    },
    "counts": {
      "B cell": 3,
      "T cell": 3
    },
    "scores": {
      "T cell": {
        "c1": 0.8356997436691734,
        "c2": 0.8567110675970135,
        "c3": 0.8259936626183231,
        "c4": -0.5599901482596032,
        "c5": -0.5763799279190517,
        "c6": -1.3820343977058553
      },
      "B cell": {
        "c1": -0.5632664728908073,
        "c2": -0.5716923966489433,
        "c3": -1.3825070611218027,
        "c4": 0.8350643975343444,
        "c5": 0.8410178115070588,
        "c6": 0.8413837216201503
      }
    }
  }
}
```

## 3) annotate_cells_spatial

- **Purpose:** Annotates spatial spots with cell-type labels by scoring user-supplied marker genes (no trained model).

- **Method:** Values are ln(1 + CP10k). For each marker gene the values are z-scored across all cells (population standard deviation; a constant gene scores 0). The score of a cell type in a cell is the mean z-score of that type's marker genes, and the cell gets the type with the highest score, or `"Unassigned"` when the best score is not above `min_score`. Marker genes missing from the matrix give a warning (an error if none of a type's markers is present). There is no trained model. Pure Python. Same definition as `annotate_cells`; spots play the role of cells.

- **API:** `dnakit.singlecell_spatial.annotate_cells_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `markers` (required; `{cell type: [marker genes]}`), `min_score` (0.0), `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `assignments` ({cell: type}), `best_score` ({cell: score}), `counts` (cells per type) and `scores` ({type: {cell: score}}).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import annotate_cells_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
markers = {"T cell": ["CD3D", "CD3E"], "B cell": ["CD79A", "MS4A1"]}

result = annotate_cells_spatial(counts, options={"markers": markers})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "annotate_cells_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score(type, cell) = mean over the type's marker genes of the per-gene z-score of ln(1 + CP10k); the highest score wins if > min_score",
    "assignments": {
      "c1": "T cell",
      "c2": "T cell",
      "c3": "T cell",
      "c4": "B cell",
      "c5": "B cell",
      "c6": "B cell"
    },
    "best_score": {
      "c1": 0.8356997436691734,
      "c2": 0.8567110675970135,
      "c3": 0.8259936626183231,
      "c4": 0.8350643975343444,
      "c5": 0.8410178115070588,
      "c6": 0.8413837216201503
    },
    "counts": {
      "B cell": 3,
      "T cell": 3
    },
    "scores": {
      "T cell": {
        "c1": 0.8356997436691734,
        "c2": 0.8567110675970135,
        "c3": 0.8259936626183231,
        "c4": -0.5599901482596032,
        "c5": -0.5763799279190517,
        "c6": -1.3820343977058553
      },
      "B cell": {
        "c1": -0.5632664728908073,
        "c2": -0.5716923966489433,
        "c3": -1.3825070611218027,
        "c4": 0.8350643975343444,
        "c5": 0.8410178115070588,
        "c6": 0.8413837216201503
      }
    }
  }
}
```

## 4) build_cell_gene_matrix

- **Purpose:** Builds a gene-by-cell count matrix from per-read or per-UMI (cell, gene) assignment records.

- **Method:** When the records carry a UMI, the count of a (cell, gene) is the number of distinct UMIs (run `correct_umis` first to merge sequencing errors); otherwise it is the number of records, or the sum of the optional `count` column. Genes and cells are sorted alphabetically. Pure Python.

- **API:** `dnakit.singlecell_spatial.build_cell_gene_matrix(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Records as TSV/CSV text or a file path (one record per line: cell, gene[, umi]; an optional header may name the columns `cell`/`barcode`, `gene`/`feature`, `umi` and `count`), a list of `[cell, gene[, umi]]` lists, or a list of mappings with those keys. Options: `max_matrix_elements` (5000000), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_records`, `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import build_cell_gene_matrix

records = """cell\tgene\tumi
c1\tCD3D\tAAAA
c1\tCD3D\tCCCC
c1\tCD3D\tCCCC
c1\tACTB\tGGGG
c2\tCD3D\tAAAA
c2\tCD79A\tTTTT
"""

result = build_cell_gene_matrix(records)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "build_cell_gene_matrix",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 7,
    "text_length": 93,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1a65abffca3e87366af0ad395ad70b1ccd7cd20630b9717e39096fb6da41429d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "counts = number of distinct UMIs per (cell, gene)",
    "n_records": 6,
    "n_genes": 3,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "ACTB",
        "CD3D",
        "CD79A"
      ],
      "cells": [
        "c1",
        "c2"
      ],
      "values": [
        [
          1.0,
          0.0
        ],
        [
          2.0,
          1.0
        ],
        [
          0.0,
          1.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 5) build_cell_peak_matrix

- **Purpose:** Builds a peak-by-cell count matrix from scATAC fragments and a list of peaks.

- **Method:** Fragments are BED-like intervals (0-based, half-open) with a cell barcode. With `count: "fragments"` a fragment is counted once for every peak it overlaps by at least 1 bp (weighted by the optional fifth column); with `"cut_sites"` the first and the last base of each fragment are counted when they fall inside a peak. Peaks are the matrix rows, named `chrom:start-end`. Fragments on chromosomes without peaks only contribute their barcode to the cell list. Pure Python.

- **API:** `dnakit.singlecell_spatial.build_cell_peak_matrix(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Fragments as TSV text or a file path (`chrom, start, end, barcode[, count]`; `#` lines are ignored) or a list of such lists. Options: `peaks` (required; `"chr1:100-600"` strings, `[chrom, start, end]` lists or BED text), `count` (`"fragments"` default or `"cut_sites"`), `max_matrix_elements` (5000000), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_fragments`, `mode`, `n_peaks` (the number of peaks, i.e. matrix rows; this field was called `n_genes` before), `n_cells`, `matrix` (`orientation`, `genes` (holds the peak names), `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import build_cell_peak_matrix

fragments = """chr1\t100\t300\tAAAC\t1
chr1\t150\t250\tAAAC\t1
chr1\t900\t1000\tAAAC\t1
chr1\t120\t220\tTTTG\t1
chr2\t50\t150\tTTTG\t1
"""
peaks = ["chr1:100-600", "chr1:800-1200", "chr2:0-100"]

result = build_cell_peak_matrix(fragments, options={"peaks": peaks})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "build_cell_peak_matrix",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 100,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4f988916a04d2fe97bb9a93fe5194995673409cb28d6eef8bb12999615188d02"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "counts = fragments (weighted by the optional 5th column) overlapping each peak by >= 1 bp; each fragment counts once per peak",
    "n_fragments": 5,
    "mode": "fragments",
    "n_peaks": 3,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "chr1:100-600",
        "chr1:800-1200",
        "chr2:0-100"
      ],
      "cells": [
        "AAAC",
        "TTTG"
      ],
      "values": [
        [
          2.0,
          1.0
        ],
        [
          1.0,
          0.0
        ],
        [
          0.0,
          1.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 6) calculate_cell_cycle_score

- **Purpose:** Scores the S and G2/M cell-cycle programmes per cell from user-supplied gene sets and assigns G1, S or G2M.

- **Method:** Values are ln(1 + CP10k). The score of a gene set is the mean expression of its genes minus the mean expression of expression-matched control genes: control genes are taken from the same mean-expression rank bin (`n_bins`, nearest rank first, up to `ctrl_size` per set gene). The control choice is deterministic, so values differ slightly from scanpy's random controls. Phase is `G1` when both scores are below 0, otherwise `S` or `G2M`, whichever score is larger. Pure Python.

- **API:** `dnakit.singlecell_spatial.calculate_cell_cycle_score(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `s_genes` (required), `g2m_genes` (required), `n_bins` (25), `ctrl_size` (50), `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `s_score` and `g2m_score` ({cell: score}), `phase` ({cell: `G1`/`S`/`G2M`}) and `phase_counts`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import calculate_cell_cycle_score

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
options = {"s_genes": ["CD3D", "CD3E"], "g2m_genes": ["CD79A", "MS4A1"], "n_bins": 3, "ctrl_size": 2}

result = calculate_cell_cycle_score(counts, options=options)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_cell_cycle_score",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score = mean ln(1+CP10k) of the gene set minus mean of expression-matched control genes (genes in th...",
    "s_score": {
      "c1": 4.309330241732457,
      "c2": 4.2525698094105415,
      "c3": 4.305041081459723,
      "c4": -2.6471780733568737,
      "c5": -3.891133860838857,
      "c6": -5.096677528135468
    },
    "g2m_score": {
      "c1": -2.4136427067829653,
      "c2": -2.537043500514193,
      "c3": -6.331498665182803,
      "c4": 4.257662662904112,
      "c5": 3.040942474661665,
      "c6": 5.755455155631239
    },
    "phase": {
      "c1": "S",
      "c2": "S",
      "c3": "S",
      "c4": "G2M",
      "c5": "G2M",
      "c6": "G2M"
    },
    "phase_counts": {
      "G1": 0,
      "S": 3,
      "G2M": 3
    }
  }
}
```

## 7) calculate_cell_interactions

- **Purpose:** Scores ligand-receptor interactions between cell types from mean expression (no significance test).

- **Method:** Values are ln(1 + CP10k). For every ordered pair of cell types (sender, receiver) and every supplied ligand-receptor pair, `score = mean(ligand in sender cells) * mean(receptor in receiver cells)`; `min_pct` drops interactions where fewer than that fraction of sender cells express the ligand or receiver cells the receptor. A pair with a gene missing from the matrix is skipped with a warning. No ligand-receptor database is bundled, so the pairs must be supplied, and no permutation or other significance test is performed. Pure Python.

- **API:** `dnakit.singlecell_spatial.calculate_cell_interactions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `cell_types` (required; `{cell: type}` or a list in column order), `pairs` (required; `[[ligand, receptor], ...]`), `min_pct` (0.0), `top_n` (100), `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `no_significance_test` (true), `n_interactions` (before the `top_n` cut) and `interactions`, a list sorted by decreasing `score` (at most `top_n`) whose items have `ligand`, `receptor`, `sender`, `receiver`, `ligand_mean`, `receptor_mean`, `score`, `pct_ligand`, `pct_receptor`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import calculate_cell_interactions

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
cell_types = ["T", "T", "T", "B", "B", "B"]
pairs = [["CD3D", "MS4A1"], ["CD79A", "CD3E"]]

result = calculate_cell_interactions(counts, options={"cell_types": cell_types, "pairs": pairs})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_cell_interactions",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score = mean ln(1+CP10k) of the ligand in sender cells x mean of the receptor in receiver cells",
    "no_significance_test": true,
    "n_interactions": 8,
    "interactions": [
      {
        "ligand": "CD3D",
        "receptor": "MS4A1",
        "sender": "T",
        "receiver": "B",
        "ligand_mean": 7.212575686641653,
        "receptor_mean": 7.161175211635025,
        "score": 51.65051821921968,
        "pct_ligand": 1.0,
        "pct_receptor": 1.0
      },
      {
        "ligand": "CD79A",
        "receptor": "CD3E",
        "sender": "B",
        "receiver": "T",
        "ligand_mean": 7.296230322614124,
        "receptor_mean": 7.056039166585788,
        "score": 51.48248692479612,
        "pct_ligand": 1.0,
        "pct_receptor": 1.0
      },
      {
        "ligand": "CD79A",
        "receptor": "CD3E",
        "sender": "B",
        "receiver": "B",
        "ligand_mean": 7.296230322614124,
        "receptor_mean": 1.7422472658960821,
        "score": 12.711837330922547,
        "pct_ligand": 1.0,
        "pct_receptor": 0.3333333333333333
      },
      {
        "ligand": "CD3D",
        "receptor": "MS4A1",
        "sender": "T",
        "receiver": "T",
        "ligand_mean": 7.212575686641653,
        "receptor_mean": 1.7612005927668049,
        "score": 12.702792574688724,
        "pct_ligand": 1.0,
        "pct_receptor": 0.3333333333333333
      },
      "..."
    ]
  }
}
```

## 8) calculate_cell_interactions_spatial

- **Purpose:** Scores ligand-receptor interactions between neighbouring spots or cells that lie within a distance radius.

- **Method:** Values are ln(1 + CP10k). For every ordered pair (a, b) of different spots within `radius` of each other, the product `ligand in a * receptor in b` is formed; `score` is the mean of these products over the pairs grouped by (sender type, receiver type), and `n_cell_pairs` is the number of such pairs. A pair with a gene missing from the matrix is skipped with a warning, and type pairs without any neighbouring spot pair give no row. No ligand-receptor database is bundled and no significance test is performed. Pure Python; at most 5000 spots.

- **API:** `dnakit.singlecell_spatial.calculate_cell_interactions_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `cell_types` (required), `pairs` (required; `[[ligand, receptor], ...]`), `radius` (required; maximum sender-receiver distance in coordinate units), `coordinates` (required; `{spot: [x, y]}`, a list of `[x, y]` in column order, or a CSV/TSV/Visium `tissue_positions` file), `min_pct` (0.0), `top_n` (100), `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `no_significance_test` (true), `n_interactions` (before the `top_n` cut) and `interactions`, a list sorted by decreasing `score` (at most `top_n`) whose items have `ligand`, `receptor`, `sender`, `receiver`, `ligand_mean`, `receptor_mean`, `score`, `pct_ligand`, `pct_receptor` and `n_cell_pairs`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import calculate_cell_interactions_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
coordinates = {"c1": [0, 0], "c2": [1, 0], "c3": [2, 0], "c4": [0, 1], "c5": [1, 1], "c6": [2, 1]}
cell_types = ["T", "T", "T", "B", "B", "B"]
pairs = [["CD3D", "MS4A1"], ["CD79A", "CD3E"]]

result = calculate_cell_interactions_spatial(
    counts,
    options={"cell_types": cell_types, "pairs": pairs, "radius": 1.5, "coordinates": coordinates},
)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_cell_interactions_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score = mean over sender-receiver cell pairs within `radius` of ln(1+CP10k)(ligand in sender) x ln(1+CP10k)(receptor in receiver)",
    "no_significance_test": true,
    "n_interactions": 8,
    "interactions": [
      {
        "ligand": "CD79A",
        "receptor": "CD3E",
        "sender": "B",
        "receiver": "T",
        "ligand_mean": 7.296979225009493,
        "receptor_mean": 7.085909282507124,
        "score": 51.70338903539261,
        "pct_ligand": 1.0,
        "pct_receptor": 1.0,
        "n_cell_pairs": 7
      },
      {
        "ligand": "CD3D",
        "receptor": "MS4A1",
        "sender": "T",
        "receiver": "B",
        "ligand_mean": 7.198078570902774,
        "receptor_mean": 7.1621555536358175,
        "score": 51.54990974557287,
        "pct_ligand": 1.0,
        "pct_receptor": 1.0,
        "n_cell_pairs": 7
      },
      {
        "ligand": "CD3D",
        "receptor": "MS4A1",
        "sender": "B",
        "receiver": "B",
        "ligand_mean": 2.6133708988441233,
        "receptor_mean": 7.16289081013641,
        "score": 18.705839909215634,
        "pct_ligand": 0.5,
        "pct_receptor": 1.0,
        "n_cell_pairs": 4
      },
      {
        "ligand": "CD79A",
        "receptor": "CD3E",
        "sender": "T",
        "receiver": "T",
        "ligand_mean": 2.661704053969166,
        "receptor_mean": 7.108311869448126,
        "score": 18.502819124425383,
        "pct_ligand": 0.5,
        "pct_receptor": 1.0,
        "n_cell_pairs": 4
      },
      "..."
    ]
  }
}
```

## 9) calculate_markers

- **Purpose:** Finds marker genes of each group (cluster or cell type) by testing it against all other cells.

- **Method:** Each group is the first group and all remaining cells the second group. Per gene: `norm = count / library size * target_sum` (with `normalization: "none"` the supplied values are used), the test runs on `log2(norm + pseudocount)`, `log2fc = log2(mean(norm, first group) + pseudocount) - log2(mean(norm, second group) + pseudocount)`, `pct_*` is the fraction of cells with raw count > 0, `p_value` comes from Welch's unequal-variance t-test (default) or a Mann-Whitney U / Wilcoxon rank-sum test (normal approximation with tie and continuity correction), and `q_value` is the Benjamini-Hochberg adjustment over all genes tested. Genes are listed by increasing p-value (ties: larger |log2fc| first). With `only_positive` (default true) only genes with `log2fc > 0` are kept. Under `welch`, groups (or remainders) with fewer than 2 cells are skipped with a warning (with the message naming the group); an error (`GROUPS_TOO_SMALL`) is raised only when no group could be tested, or `GROUPS_TOO_FEW` when there is a single label. `find_differential_features` and `transcriptomics.calculate_differential_expression` behave differently: they raise `InputFormatError` (`GROUP_TOO_SMALL`) for a too-small group. In each gene row `mean_group1`/`pct_group1` describe the tested group and `mean_group2`/`pct_group2` all other cells. Pure Python.

- **API:** `dnakit.singlecell_spatial.calculate_markers(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `groups` (required; `{cell: label}` or a list in column order), `only_positive` (true), `method` (`"welch"` default or `"wilcoxon"`), `normalization` (`"cpm"` default = scale each cell to `target_sum`, or `"none"`), `target_sum` (10000), `pseudocount` (1.0), `top_n` (100), `min_pct` (0.0; keep genes whose pct is at least this in one of the two groups), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `method`, `only_positive`, `comparison` and `markers`, a mapping group -> `n_cells` and `genes` (top_n ranked genes; each gene row has `gene`, `mean_group1`, `mean_group2` (means of the normalised values), `log2fc`, `pct_group1`, `pct_group2`, `statistic`, `p_value` and `q_value`).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import calculate_markers

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
groups = ["T", "T", "T", "B", "B", "B"]

result = calculate_markers(counts, options={"groups": groups, "top_n": 2})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_markers",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "welch",
    "only_positive": true,
    "comparison": "each group versus all other cells",
    "markers": {
      "B": {
        "n_cells": 3,
        "genes": [
          {
            "gene": "CD79A",
            "mean_group1": 1481.4814814814815,
            "mean_group2": 68.02721088435375,
            "log2fc": 4.424705081255279,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.1091841711508663,
            "p_value": 0.0893622027560642,
            "q_value": 0.18438518211850805
          },
          {
            "gene": "MS4A1",
            "mean_group1": 1296.2962962962963,
            "mean_group2": 65.35947712418302,
            "log2fc": 4.289061692972825,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.062672412366318,
            "p_value": 0.09161271184455913,
            "q_value": 0.18438518211850805
          }
        ]
      },
      "T": {
        "n_cells": 3,
        "genes": [
          {
            "gene": "CD3D",
            "mean_group1": 1363.1378477316855,
            "mean_group2": 61.7283950617284,
            "log2fc": 4.442726976505825,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.136919192646547,
            "p_value": 0.08794875669662708,
            "q_value": 0.18438518211850805
          },
          {
            "gene": "CD3E",
            "mean_group1": 1173.3582321817617,
            "mean_group2": 61.7283950617284,
            "log2fc": 4.226610101334473,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.0438371186837414,
            "p_value": 0.09219259105925402,
            "q_value": 0.18438518211850805
          }
        ]
      }
    }
  }
}
```

## 10) calculate_markers_spatial

- **Purpose:** Finds marker genes of each spatial domain or spot cluster by testing it against all other spots.

- **Method:** Same method as `calculate_markers`, with spots in place of cells. Per gene: `norm = count / library size * target_sum` (with `normalization: "none"` the supplied values are used), the test runs on `log2(norm + pseudocount)`, `log2fc = log2(mean(norm, first group) + pseudocount) - log2(mean(norm, second group) + pseudocount)`, `pct_*` is the fraction of cells with raw count > 0, `p_value` comes from Welch's unequal-variance t-test (default) or a Mann-Whitney U / Wilcoxon rank-sum test (normal approximation with tie and continuity correction), and `q_value` is the Benjamini-Hochberg adjustment over all genes tested. Genes are listed by increasing p-value (ties: larger |log2fc| first). With `only_positive` (default true) only genes with `log2fc > 0` are kept. As in `calculate_markers`, under `welch` a domain (or remainder) with fewer than 2 spots is skipped with a warning and an error (`GROUPS_TOO_SMALL`) is raised only when no domain could be tested, whereas `find_differential_features_spatial` raises `GROUP_TOO_SMALL`. `mean_group1`/`pct_group1` describe the tested domain and `mean_group2`/`pct_group2` all other spots. Pure Python.

- **API:** `dnakit.singlecell_spatial.calculate_markers_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `groups` (required; `{spot: domain label}` or a list in column order), `only_positive` (true), `method` (`"welch"` default or `"wilcoxon"`), `normalization` (`"cpm"` default = scale each cell to `target_sum`, or `"none"`), `target_sum` (10000), `pseudocount` (1.0), `top_n` (100), `min_pct` (0.0; keep genes whose pct is at least this in one of the two groups), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `method`, `only_positive`, `comparison` and `markers`, a mapping domain -> `n_cells` and `genes` (top_n ranked genes; each gene row has `gene`, `mean_group1`, `mean_group2` (means of the normalised values), `log2fc`, `pct_group1`, `pct_group2`, `statistic`, `p_value` and `q_value`).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import calculate_markers_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
domains = ["tissue1", "tissue1", "tissue1", "tissue2", "tissue2", "tissue2"]

result = calculate_markers_spatial(counts, options={"groups": domains, "top_n": 2})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_markers_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "welch",
    "only_positive": true,
    "comparison": "each group versus all other cells",
    "markers": {
      "tissue1": {
        "n_cells": 3,
        "genes": [
          {
            "gene": "CD3D",
            "mean_group1": 1363.1378477316855,
            "mean_group2": 61.7283950617284,
            "log2fc": 4.442726976505825,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.136919192646547,
            "p_value": 0.08794875669662708,
            "q_value": 0.18438518211850805
          },
          {
            "gene": "CD3E",
            "mean_group1": 1173.3582321817617,
            "mean_group2": 61.7283950617284,
            "log2fc": 4.226610101334473,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.0438371186837414,
            "p_value": 0.09219259105925402,
            "q_value": 0.18438518211850805
          }
        ]
      },
      "tissue2": {
        "n_cells": 3,
        "genes": [
          {
            "gene": "CD79A",
            "mean_group1": 1481.4814814814815,
            "mean_group2": 68.02721088435375,
            "log2fc": 4.424705081255279,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.1091841711508663,
            "p_value": 0.0893622027560642,
            "q_value": 0.18438518211850805
          },
          {
            "gene": "MS4A1",
            "mean_group1": 1296.2962962962963,
            "mean_group2": 65.35947712418302,
            "log2fc": 4.289061692972825,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.062672412366318,
            "p_value": 0.09161271184455913,
            "q_value": 0.18438518211850805
          }
        ]
      }
    }
  }
}
```

## 11) calculate_perturbation_effect

- **Purpose:** Tests the transcriptome-wide effect of each perturbation against the control cells.

- **Method:** Each non-control label is compared with the control cells gene by gene with the two-group test used by `find_differential_features` (the perturbation is the first group, the control the second). Per gene: `norm = count / library size * target_sum` (with `normalization: "none"` the supplied values are used), the test runs on `log2(norm + pseudocount)`, `log2fc = log2(mean(norm, first group) + pseudocount) - log2(mean(norm, second group) + pseudocount)`, `pct_*` is the fraction of cells with raw count > 0, `p_value` comes from Welch's unequal-variance t-test (default) or a Mann-Whitney U / Wilcoxon rank-sum test (normal approximation with tie and continuity correction), and `q_value` is the Benjamini-Hochberg adjustment over all genes tested. Genes are listed by increasing p-value (ties: larger |log2fc| first). A perturbation with fewer than 2 cells (or fewer than 2 control cells) is skipped with a warning under `welch`; if no perturbation can be tested an error is raised. In every gene row the `_group1` fields (`mean_group1`, `pct_group1`) describe the perturbed cells and the `_group2` fields the control cells; positive `log2fc` means higher in the perturbed cells (perturbation versus control, not group versus all other cells). Pure Python.

- **API:** `dnakit.singlecell_spatial.calculate_perturbation_effect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `perturbations` (required; `{cell: perturbation label}` or a list), `control` (`"control"`), `method` (`"welch"` default or `"wilcoxon"`), `normalization` (`"cpm"` default = scale each cell to `target_sum`, or `"none"`), `target_sum` (10000), `pseudocount` (1.0), `top_n` (100), `min_pct` (0.0; keep genes whose pct is at least this in one of the two groups), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `control`, `n_control_cells`, `comparison`, `method` and `effects`, a mapping perturbation -> `n_cells` and `genes` (top_n ranked genes; each gene row has `gene`, `mean_group1`, `mean_group2` (means of the normalised values), `log2fc`, `pct_group1`, `pct_group2`, `statistic`, `p_value` and `q_value`).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import calculate_perturbation_effect

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
perturbations = ["KO", "KO", "KO", "control", "control", "control"]

result = calculate_perturbation_effect(counts, options={"perturbations": perturbations, "top_n": 2})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_perturbation_effect",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "control": "control",
    "n_control_cells": 3,
    "comparison": "each perturbation versus control cells",
    "method": "welch",
    "effects": {
      "KO": {
        "n_cells": 3,
        "genes": [
          {
            "gene": "CD3D",
            "mean_group1": 1363.1378477316855,
            "mean_group2": 61.7283950617284,
            "log2fc": 4.442726976505825,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.136919192646547,
            "p_value": 0.08794875669662708,
            "q_value": 0.18438518211850805
          },
          {
            "gene": "CD79A",
            "mean_group1": 68.02721088435375,
            "mean_group2": 1481.4814814814815,
            "log2fc": -4.424705081255279,
            "pct_group1": 0.3333333333333333,
            "pct_group2": 1.0,
            "statistic": -3.1091841711508663,
            "p_value": 0.0893622027560642,
            "q_value": 0.18438518211850805
          }
        ]
      }
    }
  }
}
```

## 12) cluster

- **Purpose:** Clusters cells with k-means on the top principal components of log-normalised expression.

- **Method:** Expression is ln(1 + CP10k) (each cell scaled to `target_sum` = 10000 counts, then natural log of 1 + x) (values are ln(1 + CP10k), centred per gene (optionally scaled to unit variance and restricted to the `n_top_genes` most variable genes); the principal components come from the SVD of that matrix, each component's sign fixed so its largest-magnitude loading is positive). k-means uses a deterministic farthest-point start and Lloyd iterations; labels are numbered by first appearance. With `n_clusters: "auto"` k is chosen from 2 to `max_clusters` by the best mean silhouette coefficient. Requires numpy (returns `not_implemented` if it is missing); at most 20000 cells (5000 with `"auto"`).

- **API:** `dnakit.singlecell_spatial.cluster(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `n_clusters` (`"auto"` default, or an integer), `max_clusters` (10), `n_pcs` (10; 0 = no PCA), `n_top_genes` (all; otherwise the highest-variance genes), `scale` (false), `target_sum` (10000), `return_embedding` (false), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_clusters`, `labels` ({cell: cluster number as a string, numbered by first appearance}), `cluster_sizes`, `inertia`, `explained_variance_ratio` (`null` when `n_pcs` is 0), `silhouette_by_k` (mean silhouette per tried k when `n_clusters` is `"auto"`, else `null`) and, with `return_embedding`, `embedding` ({cell: coordinates}).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import cluster

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""

result = cluster(counts, options={"n_clusters": 2, "n_pcs": 2})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "cluster",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "k-means (deterministic farthest-point start, Lloyd iterations) on the top PCs of ln(1+CP10k) values; k from options or the best mean silhouette",
    "n_clusters": 2,
    "labels": {
      "c1": "0",
      "c2": "0",
      "c3": "0",
      "c4": "1",
      "c5": "1",
      "c6": "1"
    },
    "cluster_sizes": {
      "0": 3,
      "1": 3
    },
    "inertia": 61.986183555554646,
    "explained_variance_ratio": [
      0.650159103044713,
      0.17804888640013142
    ],
    "silhouette_by_k": null
  }
}
```

## 13) cluster_spatial

- **Purpose:** Clusters spots into spatial domains with k-means on spatially smoothed principal components.

- **Method:** Heuristic, not a published method: the principal components are computed as in `cluster`, then each spot's PC vector is replaced by `(1 - alpha) * own + alpha * mean of its n_neighbors nearest spatial neighbours` (nearest by Euclidean distance between coordinates), and k-means (with the `n_clusters`/silhouette rule of `cluster`) is run on the smoothed vectors. Requires numpy (returns `not_implemented` if it is missing); at most 5000 spots.

- **API:** `dnakit.singlecell_spatial.cluster_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `coordinates` (required; `{spot: [x, y]}`, a list of `[x, y]` in column order, or a CSV/TSV/Visium `tissue_positions` file), `n_neighbors` (6; must be smaller than the number of spots), `alpha` (0.5), `n_clusters` (`"auto"` or an integer), `max_clusters` (10), `n_pcs` (10; 0 = no PCA), `n_top_genes` (all; otherwise the highest-variance genes), `scale` (false), `target_sum` (10000), `return_embedding` (false), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_clusters`, `labels` ({cell: cluster number as a string, numbered by first appearance}), `cluster_sizes`, `inertia`, `explained_variance_ratio` (`null` when `n_pcs` is 0), `silhouette_by_k` (mean silhouette per tried k when `n_clusters` is `"auto"`, else `null`) and, with `return_embedding`, `embedding` ({cell: coordinates}) (`inertia` and `embedding` refer to the smoothed vectors).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import cluster_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
coordinates = {"c1": [0, 0], "c2": [1, 0], "c3": [2, 0], "c4": [0, 1], "c5": [1, 1], "c6": [2, 1]}

result = cluster_spatial(counts, options={"coordinates": coordinates, "n_clusters": 2, "n_pcs": 2, "n_neighbors": 2})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "cluster_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "spatial-domain heuristic: each spot's PC vector is replaced by (1-alpha)*own + alpha*mean of its 2 n...",
    "n_clusters": 2,
    "labels": {
      "c1": "0",
      "c2": "0",
      "c3": "0",
      "c4": "1",
      "c5": "1",
      "c6": "1"
    },
    "cluster_sizes": {
      "0": 3,
      "1": 3
    },
    "inertia": 11.41227429529395,
    "explained_variance_ratio": [
      0.650159103044713,
      0.17804888640013142
    ],
    "silhouette_by_k": null
  }
}
```

## 14) correct_umis

- **Purpose:** Corrects sequencing errors in UMIs by collapsing groups around an abundant UMI (UMI-tools directional, adjacency and cluster algorithms).

- **Method:** Implements the UMI-tools network methods on Hamming distance (at most `threshold` mismatches). `directional`: starting from the most abundant UMI, absorb neighbours `b` with `count_a >= 2 * count_b - 1`, transitively; `adjacency`: the most abundant unassigned UMI absorbs all its unassigned neighbours in one step; `cluster`: connected components, represented by the most abundant UMI; `unique`: no merging. All UMIs must have the same length and use only A, C, G, T, N. Pure Python.

- **API:** `dnakit.singlecell_spatial.correct_umis(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A list of UMI strings (each occurrence counts as one read), `{UMI: read count}`, or text or a file with one `UMI` or `UMI<TAB>count` per line. Options: `method` (`"directional"` default, `"adjacency"`, `"cluster"` or `"unique"`), `threshold` (1; maximum Hamming distance, 1 to 4).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `method`, `threshold`, `n_reads`, `n_unique_before`, `n_unique_after`, `corrected` ({merged UMI: representative UMI}) and `groups` (representative, members and read count of each group) with `groups_truncated`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import correct_umis

umis = ["AAAA", "AAAA", "AAAA", "AAAA", "AAAT", "CCCC", "CCCC", "GGGG"]

result = correct_umis(umis)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "correct_umis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 8,
    "text_length": 64,
    "sequence_bases": 32,
    "lengths": [4, 4, 4, 4, 4, 4, 4, 4],
    "sha256": "7051d7f210a596683b53aab27d96e3402faa226fa42320cbb7d67085c757c7ec"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "UMI-tools directional: starting from the most abundant UMI, absorb neighbours within the Hamming thr...",
    "method": "directional",
    "threshold": 1,
    "n_reads": 8,
    "n_unique_before": 4,
    "n_unique_after": 3,
    "corrected": {
      "AAAT": "AAAA"
    },
    "groups": [
      {
        "representative": "AAAA",
        "members": [
          "AAAA",
          "AAAT"
        ],
        "reads": 5
      },
      {
        "representative": "CCCC",
        "members": [
          "CCCC"
        ],
        "reads": 2
      },
      {
        "representative": "GGGG",
        "members": [
          "GGGG"
        ],
        "reads": 1
      }
    ],
    "groups_truncated": false
  }
}
```

## 15) correct_umis_spatial

- **Purpose:** Corrects sequencing errors in spatial-barcoded UMIs by collapsing groups around an abundant UMI (UMI-tools algorithms).

- **Method:** Implements the UMI-tools network methods on Hamming distance (at most `threshold` mismatches). `directional`: starting from the most abundant UMI, absorb neighbours `b` with `count_a >= 2 * count_b - 1`, transitively; `adjacency`: the most abundant unassigned UMI absorbs all its unassigned neighbours in one step; `cluster`: connected components, represented by the most abundant UMI; `unique`: no merging. All UMIs must have the same length and use only A, C, G, T, N. Pure Python. The same algorithms as `correct_umis`; call it on the UMIs of one spot and gene.

- **API:** `dnakit.singlecell_spatial.correct_umis_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A list of UMI strings, `{UMI: read count}`, or text or a file with one `UMI` or `UMI<TAB>count` per line (the UMIs of one spot and gene). Options: `method` (`"directional"` default, `"adjacency"`, `"cluster"` or `"unique"`), `threshold` (1; maximum Hamming distance, 1 to 4).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `method`, `threshold`, `n_reads`, `n_unique_before`, `n_unique_after`, `corrected` ({merged UMI: representative UMI}) and `groups` (representative, members and read count of each group) with `groups_truncated`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import correct_umis_spatial

umi_counts = {"ACGTAC": 10, "ACGTAA": 1, "TTTTTT": 3}

result = correct_umis_spatial(umi_counts)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "correct_umis_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 40,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5b6db5fa8fe58914529df0bf23f8e0e2578f99c5ee8538e1b252c9973fa8c775"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "UMI-tools directional: starting from the most abundant UMI, absorb neighbours within the Hamming thr...",
    "method": "directional",
    "threshold": 1,
    "n_reads": 14,
    "n_unique_before": 3,
    "n_unique_after": 2,
    "corrected": {
      "ACGTAA": "ACGTAC"
    },
    "groups": [
      {
        "representative": "ACGTAC",
        "members": [
          "ACGTAC",
          "ACGTAA"
        ],
        "reads": 11
      },
      {
        "representative": "TTTTTT",
        "members": [
          "TTTTTT"
        ],
        "reads": 3
      }
    ],
    "groups_truncated": false
  }
}
```

## 16) deconvolute_spots

- **Purpose:** Estimates the cell-type fractions of each spot by non-negative least squares on a signature matrix.

- **Method:** Each spot's profile is fitted as a non-negative combination of the signature columns with non-negative least squares (`scipy.optimize.nnls`) over the genes shared with the signature; the coefficients divided by their sum give the cell-type fractions, and `rmse` is the fit residual. With `normalize` (default true) every column is first scaled to sum 1. If `reference` and `reference_labels` are given instead of `signature`, the signature is the mean library-size-normalised profile per cell type. Requires numpy and scipy (returns `not_implemented` if either is missing). It is not a model-based method such as RCTD or cell2location.

- **API:** `dnakit.singlecell_spatial.deconvolute_spots(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Spots to deconvolve: a gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `signature` (genes x cell types matrix in any matrix form) OR `reference` (count matrix of annotated cells) plus `reference_labels` (`{cell: type}` or a list), `normalize` (true), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `cell_types`, `n_shared_genes`, `fractions` ({spot: {cell type: fraction}}) and `rmse` ({spot: residual}).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import deconvolute_spots

spots = """gene,s1,s2
A,100,10
B,10,100
C,50,50
"""
signature = """gene,typeX,typeY
A,10,1
B,1,10
C,5,5
"""

result = deconvolute_spots(spots, options={"signature": signature})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "deconvolute_spots",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 37,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d0743d02cb300ebb3b2966930e2cef0eac47e89410bcaa86785506703f7c27f3"
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
      "s1": {
        "typeX": 1.0,
        "typeY": 4.636616585542447e-17
      },
      "s2": {
        "typeX": 0.0,
        "typeY": 1.0
      }
    },
    "rmse": {
      "s1": 3.0913372807800604e-17,
      "s2": 2.1817011766850644e-17
    }
  }
}
```

## 17) deconvolute_spots_spatial

- **Purpose:** Estimates the cell-type fractions of each spatial spot by non-negative least squares on a signature matrix.

- **Method:** Each spot's profile is fitted as a non-negative combination of the signature columns with non-negative least squares (`scipy.optimize.nnls`) over the genes shared with the signature; the coefficients divided by their sum give the cell-type fractions, and `rmse` is the fit residual. With `normalize` (default true) every column is first scaled to sum 1. If `reference` and `reference_labels` are given instead of `signature`, the signature is the mean library-size-normalised profile per cell type. Requires numpy and scipy (returns `not_implemented` if either is missing). It is not a model-based method such as RCTD or cell2location. Same method as `deconvolute_spots`.

- **API:** `dnakit.singlecell_spatial.deconvolute_spots_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Spots to deconvolve: a gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `signature` (genes x cell types matrix in any matrix form) OR `reference` (count matrix of annotated cells) plus `reference_labels` (`{cell: type}` or a list), `normalize` (true), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `cell_types`, `n_shared_genes`, `fractions` ({spot: {cell type: fraction}}) and `rmse` ({spot: residual}).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import deconvolute_spots_spatial

spots = """gene,s1,s2
A,100,10
B,10,100
C,50,50
"""
signature = """gene,typeX,typeY
A,10,1
B,1,10
C,5,5
"""

result = deconvolute_spots_spatial(spots, options={"signature": signature})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "deconvolute_spots_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 37,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d0743d02cb300ebb3b2966930e2cef0eac47e89410bcaa86785506703f7c27f3"
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
      "s1": {
        "typeX": 1.0,
        "typeY": 4.636616585542447e-17
      },
      "s2": {
        "typeX": 0.0,
        "typeY": 1.0
      }
    },
    "rmse": {
      "s1": 3.0913372807800604e-17,
      "s2": 2.1817011766850644e-17
    }
  }
}
```

## 18) demultiplex

- **Purpose:** Assigns cells to samples from hashtag-oligo (HTO) counts with a simple count-threshold rule.

- **Method:** A threshold heuristic on hashtag-oligo (HTO) counts, not the mixture-model fit of HTODemux. With the highest hashtag count `top` and the second highest `second`: `Negative` if `top < min_counts`; `Doublet` if `second >= min_counts` and `top < ratio * second`; otherwise `Singlet`, assigned the name of the top hashtag. Pure Python.

- **API:** `dnakit.singlecell_spatial.demultiplex(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A hashtag-by-cell count matrix (rows are hashtags, columns are cells) in any matrix form: a mapping `{hashtag: [counts per cell]}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = hashtag names), a numpy array or pandas DataFrame, or a 10x directory. Options: `min_counts` (5.0), `ratio` (3.0), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `min_counts`, `ratio`, `assignments` ({cell: hashtag name, `Doublet` or `Negative`}) and `counts` (number of cells per assignment).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import demultiplex

hashtags = """hashtag,c1,c2,c3,c4
HTO1,50,2,30,1
HTO2,3,40,28,0
"""

result = demultiplex(hashtags)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "demultiplex",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 50,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "209261238cd9bd6c030720432ba3df609b7a02316b34df54fa4d4bc42bd778f7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "Negative if top hashtag count < 5.0; Doublet if second count >= 5.0 and top < 3.0 x second; otherwis...",
    "min_counts": 5.0,
    "ratio": 3.0,
    "assignments": {
      "c1": "HTO1",
      "c2": "HTO2",
      "c3": "Doublet",
      "c4": "Negative"
    },
    "counts": {
      "Doublet": 1,
      "HTO1": 1,
      "HTO2": 1,
      "Negative": 1
    }
  }
}
```

## 19) demultiplex_spatial

- **Purpose:** Assigns spots to samples from hashtag-oligo counts with a simple count-threshold rule.

- **Method:** A threshold heuristic on hashtag-oligo (HTO) counts, not the mixture-model fit of HTODemux. With the highest hashtag count `top` and the second highest `second`: `Negative` if `top < min_counts`; `Doublet` if `second >= min_counts` and `top < ratio * second`; otherwise `Singlet`, assigned the name of the top hashtag. Pure Python. Same rule as `demultiplex`, with spots in place of cells.

- **API:** `dnakit.singlecell_spatial.demultiplex_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A hashtag-by-spot count matrix (rows are hashtags, columns are spots) in any matrix form: a mapping `{hashtag: [counts per spot]}`, a list of rows, CSV/TSV text or a file path (header = spot names, first column = hashtag names), a numpy array or pandas DataFrame, or a 10x directory. Options: `min_counts` (5.0), `ratio` (3.0), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `min_counts`, `ratio`, `assignments` ({cell: hashtag name, `Doublet` or `Negative`}) and `counts` (number of cells per assignment).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import demultiplex_spatial

hashtags = """hashtag,s1,s2,s3,s4
HTO1,50,2,30,1
HTO2,3,40,28,0
"""

result = demultiplex_spatial(hashtags)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "demultiplex_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 50,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3613a5996b191d009f5788ae9d7e173227292c4dae93e4de467f7299851a7c9a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "Negative if top hashtag count < 5.0; Doublet if second count >= 5.0 and top < 3.0 x second; otherwis...",
    "min_counts": 5.0,
    "ratio": 3.0,
    "assignments": {
      "s1": "HTO1",
      "s2": "HTO2",
      "s3": "Doublet",
      "s4": "Negative"
    },
    "counts": {
      "Doublet": 1,
      "HTO1": 1,
      "HTO2": 1,
      "Negative": 1
    }
  }
}
```

## 20) extract_barcodes

- **Purpose:** Extracts the cell barcode (a fixed read window) from reads and counts barcodes, optionally against a whitelist.

- **Method:** Barcode = `read[barcode_start : barcode_start + barcode_length]` (0-based). Reads shorter than the window are counted in `n_too_short` and skipped (`n_reads = n_extracted + n_unmatched + n_too_short`); barcodes containing N are counted in `n_barcode_with_N` and cannot be corrected. With a whitelist a barcode is accepted when it matches exactly or lies at Hamming distance 1 from exactly one whitelist barcode (then it is corrected to it; whitelist entries are upper-cased and a `-1` style suffix is removed). Each listed read gets a `corrected` field (`null` when unmatched), unmatched reads are not counted in `barcode_counts`, and `whitelist` reports `size`, `exact`, `corrected` and `unmatched`. Pure Python.

- **API:** `dnakit.singlecell_spatial.extract_barcodes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** FASTQ/FASTA text or a file path, a list of sequences, or `{read name: sequence}` (read 1 of a droplet protocol). Options: `barcode_start` (0), `barcode_length` (16), `whitelist` (list, text or path of valid barcodes), `top_n` (100).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_reads` (all reads), `n_extracted` (reads whose barcode was counted: without a whitelist every read long enough, with a whitelist only exact or corrected matches), `n_unmatched` (reads whose barcode is in no whitelist entry within one mismatch; always 0 without a whitelist), `n_too_short`, `n_barcode_with_N`, `n_unique_barcodes`, `barcode_counts` ({barcode: reads}, most frequent first, at most `top_n`) and `reads` (the barcode of each read, at most 1000 listed) with `reads_truncated`; with a whitelist also `whitelist`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import extract_barcodes

reads = """@r1
ACGTTTGGAACC
+
IIIIIIIIIIII
@r2
ACGATTGGAACC
+
IIIIIIIIIIII
@r3
TTTTTTGGAACC
+
IIIIIIIIIIII
@r4
AC
+
II
@r5
GGGGTTGGAACC
+
IIIIIIIIIIII
"""

result = extract_barcodes(reads, options={"barcode_start": 0, "barcode_length": 4, "whitelist": ["ACGT", "TTTT"]})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "extract_barcodes",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 140,
    "sequence_bases": 50,
    "lengths": [12, 12, 12, 2, 12],
    "sha256": "c0075b27af6dfe28d96714bc947d410ffc32d803fd54970c92ac95498e32e0a3"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "barcode = read[0:4] (0-based, 4 bases) matched to the whitelist exactly or by a unique 1-mismatch correction",
    "n_reads": 5,
    "n_extracted": 3,
    "n_unmatched": 1,
    "n_too_short": 1,
    "n_barcode_with_N": 0,
    "n_unique_barcodes": 2,
    "barcode_counts": {
      "ACGT": 2,
      "TTTT": 1
    },
    "reads": [
      {
        "read": "r1",
        "barcode": "ACGT",
        "corrected": "ACGT"
      },
      {
        "read": "r2",
        "barcode": "ACGA",
        "corrected": "ACGT"
      },
      {
        "read": "r3",
        "barcode": "TTTT",
        "corrected": "TTTT"
      },
      {
        "read": "r5",
        "barcode": "GGGG",
        "corrected": null
      }
    ],
    "reads_truncated": false,
    "whitelist": {
      "size": 2,
      "exact": 2,
      "corrected": 1,
      "unmatched": 1
    }
  }
}
```

## 21) extract_barcodes_spatial

- **Purpose:** Extracts the spatial (spot) barcode from reads and counts barcodes, optionally against a whitelist.

- **Method:** Barcode = `read[barcode_start : barcode_start + barcode_length]` (0-based), as in `extract_barcodes`. Reads shorter than the window are counted in `n_too_short` and skipped; `n_reads = n_extracted + n_unmatched + n_too_short`. With a whitelist a barcode is accepted when it matches exactly or lies at Hamming distance 1 from exactly one whitelist barcode (then it is corrected to it; whitelist entries are upper-cased and a `-1` style suffix is removed). Each listed read gets a `corrected` field (`null` when unmatched), unmatched reads are not counted in `barcode_counts`, and `whitelist` reports `size`, `exact`, `corrected` and `unmatched`. Pure Python.

- **API:** `dnakit.singlecell_spatial.extract_barcodes_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** FASTQ/FASTA text or a file path, a list of sequences, or `{read name: sequence}` (for example Visium read 1). Options: `barcode_start` (0), `barcode_length` (16), `whitelist` (list, text or path of valid barcodes), `top_n` (100).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_reads` (all reads), `n_extracted` (reads whose barcode was counted: without a whitelist every read long enough, with a whitelist only exact or corrected matches), `n_unmatched` (reads whose barcode is in no whitelist entry within one mismatch; always 0 without a whitelist), `n_too_short`, `n_barcode_with_N`, `n_unique_barcodes`, `barcode_counts` ({barcode: reads}, most frequent first, at most `top_n`) and `reads` (the barcode of each read, at most 1000 listed) with `reads_truncated`; with a whitelist also `whitelist`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import extract_barcodes_spatial

reads = {"r1": "ACGTTTGGAACC", "r2": "ACGTTTGGAAGG", "r3": "TTTTCCGGAACC"}

result = extract_barcodes_spatial(reads, options={"barcode_length": 4})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "extract_barcodes_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 66,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d8f2418172f6179af1f8f4228df8b4c23252a51e8cb59cf96a82640b296f3947"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "barcode = read[0:4] (0-based, 4 bases)",
    "n_reads": 3,
    "n_extracted": 3,
    "n_unmatched": 0,
    "n_too_short": 0,
    "n_barcode_with_N": 0,
    "n_unique_barcodes": 2,
    "barcode_counts": {
      "ACGT": 2,
      "TTTT": 1
    },
    "reads": [
      {
        "read": "r1",
        "barcode": "ACGT"
      },
      {
        "read": "r2",
        "barcode": "ACGT"
      },
      {
        "read": "r3",
        "barcode": "TTTT"
      }
    ],
    "reads_truncated": false
  }
}
```

## 22) extract_umis

- **Purpose:** Extracts the UMI (a fixed read window) from reads and counts the UMIs.

- **Method:** UMI = `read[umi_start : umi_start + umi_length]` (0-based); the defaults are the 10x Chromium v3 layout (16 bp barcode followed by a 12 bp UMI). Reads shorter than the window are counted in `n_too_short` and skipped; UMIs containing N are counted in `n_umi_with_N`. Pure Python.

- **API:** `dnakit.singlecell_spatial.extract_umis(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** FASTQ/FASTA text or a file path, a list of sequences, or `{read name: sequence}` (read 1 of a droplet protocol). Options: `umi_start` (16), `umi_length` (12), `top_n` (100).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_reads`, `n_extracted`, `n_too_short`, `n_umi_with_N`, `n_unique_umis`, `umi_counts` ({UMI: reads}, most frequent first, at most `top_n`) and `reads` (the UMI of each read, at most 1000 listed) with `reads_truncated`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import extract_umis

reads = """@r1
ACGTTTGGAACC
+
IIIIIIIIIIII
@r2
ACGTTTGGAACC
+
IIIIIIIIIIII
@r3
ACGTCCCCAAAA
+
IIIIIIIIIIII
"""

result = extract_umis(reads, options={"umi_start": 4, "umi_length": 4})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "extract_umis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 96,
    "sequence_bases": 36,
    "lengths": [12, 12, 12],
    "sha256": "872d4f5e9ac65abcc42c51edf93e26d961e0279f6e8b0fa5ebab99dd354d59ed"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "umi = read[4:8] (0-based, 4 bases)",
    "n_reads": 3,
    "n_extracted": 3,
    "n_too_short": 0,
    "n_umi_with_N": 0,
    "n_unique_umis": 2,
    "umi_counts": {
      "TTGG": 2,
      "CCCC": 1
    },
    "reads": [
      {
        "read": "r1",
        "umi": "TTGG"
      },
      {
        "read": "r2",
        "umi": "TTGG"
      },
      {
        "read": "r3",
        "umi": "CCCC"
      }
    ],
    "reads_truncated": false
  }
}
```

## 23) extract_umis_spatial

- **Purpose:** Extracts the UMI (a fixed read window) from spatial reads and counts the UMIs.

- **Method:** UMI = `read[umi_start : umi_start + umi_length]` (0-based), as in `extract_umis`. Reads shorter than the window are counted in `n_too_short` and skipped; UMIs containing N are counted in `n_umi_with_N`. Pure Python.

- **API:** `dnakit.singlecell_spatial.extract_umis_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** FASTQ/FASTA text or a file path, a list of sequences, or `{read name: sequence}` (for example Visium read 1). Options: `umi_start` (16), `umi_length` (12), `top_n` (100).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_reads`, `n_extracted`, `n_too_short`, `n_umi_with_N`, `n_unique_umis`, `umi_counts` ({UMI: reads}, most frequent first, at most `top_n`) and `reads` (the UMI of each read, at most 1000 listed) with `reads_truncated`. Reads given as a plain list are named `read_1`, `read_2`, ....

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import extract_umis_spatial

reads = ["ACGTTTGGAACC", "ACGTTTGGAACC", "ACGTCCCCAAAA"]

result = extract_umis_spatial(reads, options={"umi_start": 4, "umi_length": 4})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "extract_umis_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 48,
    "sequence_bases": 36,
    "lengths": [
      12,
      12,
      12
    ],
    "sha256": "ea8969ac1fdce4d906a276683a791ef5ae593aed7e80670ff58a88fe43c9eaa4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "umi = read[4:8] (0-based, 4 bases)",
    "n_reads": 3,
    "n_extracted": 3,
    "n_too_short": 0,
    "n_umi_with_N": 0,
    "n_unique_umis": 2,
    "umi_counts": {
      "TTGG": 2,
      "CCCC": 1
    },
    "reads": [
      {
        "read": "read_1",
        "umi": "TTGG"
      },
      {
        "read": "read_2",
        "umi": "TTGG"
      },
      {
        "read": "read_3",
        "umi": "CCCC"
      }
    ],
    "reads_truncated": false
  }
}
```

## 24) find_differential_features

- **Purpose:** Tests differential expression between two groups of cells (Welch t-test or Wilcoxon rank-sum).

- **Method:** Per gene: `norm = count / library size * target_sum` (with `normalization: "none"` the supplied values are used), the test runs on `log2(norm + pseudocount)`, `log2fc = log2(mean(norm, first group) + pseudocount) - log2(mean(norm, second group) + pseudocount)`, `pct_*` is the fraction of cells with raw count > 0, `p_value` comes from Welch's unequal-variance t-test (default) or a Mann-Whitney U / Wilcoxon rank-sum test (normal approximation with tie and continuity correction), and `q_value` is the Benjamini-Hochberg adjustment over all genes tested. Genes are listed by increasing p-value (ties: larger |log2fc| first). `group1` is the first group and `group2` the second. Welch's test needs at least 2 cells in each group (otherwise `InputFormatError` with code `GROUP_TOO_SMALL` is raised; use `wilcoxon`). This differs from `calculate_markers`, which skips too-small groups with a warning. The comparison is `group1` versus `group2`, not a group against all other cells; `mean_group1`/`pct_group1` describe `group1` and positive `log2fc` means higher in `group1`. Pure Python.

- **API:** `dnakit.singlecell_spatial.find_differential_features(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `groups` (required; `{cell: label}` or a list), `group1` and `group2` (labels to compare; optional when there are exactly two labels, then the sorted labels are used), `method` (`"welch"` default or `"wilcoxon"`), `normalization` (`"cpm"` default = scale each cell to `target_sum`, or `"none"`), `target_sum` (10000), `pseudocount` (1.0), `top_n` (100), `min_pct` (0.0; keep genes whose pct is at least this in one of the two groups), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `group1`, `group2`, `n_group1`, `n_group2`, `method`, `normalization`, `target_sum`, `pseudocount`, `n_genes_tested` and `genes` (top_n ranked genes; each gene row has `gene`, `mean_group1`, `mean_group2` (means of the normalised values), `log2fc`, `pct_group1`, `pct_group2`, `statistic`, `p_value` and `q_value`).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import find_differential_features

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
groups = ["T", "T", "T", "B", "B", "B"]

result = find_differential_features(counts, options={"groups": groups, "top_n": 2})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "find_differential_features",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "group1": "B",
    "group2": "T",
    "n_group1": 3,
    "n_group2": 3,
    "method": "welch",
    "normalization": "cpm",
    "target_sum": 10000.0,
    "pseudocount": 1.0,
    "n_genes_tested": 8,
    "genes": [
      {
        "gene": "CD3D",
        "mean_group1": 61.7283950617284,
        "mean_group2": 1363.1378477316855,
        "log2fc": -4.442726976505825,
        "pct_group1": 0.3333333333333333,
        "pct_group2": 1.0,
        "statistic": -3.136919192646547,
        "p_value": 0.08794875669662708,
        "q_value": 0.18438518211850805
      },
      {
        "gene": "CD79A",
        "mean_group1": 1481.4814814814815,
        "mean_group2": 68.02721088435375,
        "log2fc": 4.424705081255279,
        "pct_group1": 1.0,
        "pct_group2": 0.3333333333333333,
        "statistic": 3.1091841711508663,
        "p_value": 0.0893622027560642,
        "q_value": 0.18438518211850805
      }
    ]
  }
}
```

## 25) find_differential_features_spatial

- **Purpose:** Tests differential expression between two groups of spots or spatial domains.

- **Method:** Per gene: `norm = count / library size * target_sum` (with `normalization: "none"` the supplied values are used), the test runs on `log2(norm + pseudocount)`, `log2fc = log2(mean(norm, first group) + pseudocount) - log2(mean(norm, second group) + pseudocount)`, `pct_*` is the fraction of cells with raw count > 0, `p_value` comes from Welch's unequal-variance t-test (default) or a Mann-Whitney U / Wilcoxon rank-sum test (normal approximation with tie and continuity correction), and `q_value` is the Benjamini-Hochberg adjustment over all genes tested. Genes are listed by increasing p-value (ties: larger |log2fc| first). `group1` is the first group and `group2` the second, spots playing the role of cells. Welch's test needs at least 2 spots in each group (otherwise `InputFormatError` with code `GROUP_TOO_SMALL` is raised; use `wilcoxon`). This differs from `calculate_markers_spatial`, which skips too-small domains with a warning. The comparison is `group1` versus `group2`, not a group against all other spots; positive `log2fc` means higher in `group1`. Pure Python.

- **API:** `dnakit.singlecell_spatial.find_differential_features_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `groups` (required), `group1` and `group2` (optional when there are exactly two labels), `method` (`"welch"` default or `"wilcoxon"`), `normalization` (`"cpm"` default = scale each cell to `target_sum`, or `"none"`), `target_sum` (10000), `pseudocount` (1.0), `top_n` (100), `min_pct` (0.0; keep genes whose pct is at least this in one of the two groups), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `group1`, `group2`, `n_group1`, `n_group2`, `method`, `normalization`, `target_sum`, `pseudocount`, `n_genes_tested` and `genes` (top_n ranked genes; each gene row has `gene`, `mean_group1`, `mean_group2` (means of the normalised values), `log2fc`, `pct_group1`, `pct_group2`, `statistic`, `p_value` and `q_value`).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import find_differential_features_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
domains = ["tissue1", "tissue1", "tissue1", "tissue2", "tissue2", "tissue2"]

result = find_differential_features_spatial(counts, options={"groups": domains, "top_n": 2})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "find_differential_features_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "group1": "tissue1",
    "group2": "tissue2",
    "n_group1": 3,
    "n_group2": 3,
    "method": "welch",
    "normalization": "cpm",
    "target_sum": 10000.0,
    "pseudocount": 1.0,
    "n_genes_tested": 8,
    "genes": [
      {
        "gene": "CD3D",
        "mean_group1": 1363.1378477316855,
        "mean_group2": 61.7283950617284,
        "log2fc": 4.442726976505825,
        "pct_group1": 1.0,
        "pct_group2": 0.3333333333333333,
        "statistic": 3.136919192646547,
        "p_value": 0.08794875669662708,
        "q_value": 0.18438518211850805
      },
      {
        "gene": "CD79A",
        "mean_group1": 68.02721088435375,
        "mean_group2": 1481.4814814814815,
        "log2fc": -4.424705081255279,
        "pct_group1": 0.3333333333333333,
        "pct_group2": 1.0,
        "statistic": -3.1091841711508663,
        "p_value": 0.0893622027560642,
        "q_value": 0.18438518211850805
      }
    ]
  }
}
```

## 26) infer_cell_type_from_markers

- **Purpose:** Infers the cell type of each cell from user-supplied marker genes by marker-score ranking.

- **Method:** Values are ln(1 + CP10k). For each marker gene the values are z-scored across all cells (population standard deviation; a constant gene scores 0). The score of a cell type in a cell is the mean z-score of that type's marker genes, and the cell gets the type with the highest score, or `"Unassigned"` when the best score is not above `min_score`. Marker genes missing from the matrix give a warning (an error if none of a type's markers is present). There is no trained model. Pure Python. This is the same computation as `annotate_cells`.

- **API:** `dnakit.singlecell_spatial.infer_cell_type_from_markers(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `markers` (required; `{cell type: [marker genes]}`), `min_score` (0.0), `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `assignments` ({cell: type}), `best_score` ({cell: score}), `counts` (cells per type) and `scores` ({type: {cell: score}}).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import infer_cell_type_from_markers

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
markers = {"T cell": ["CD3D", "CD3E"], "B cell": ["CD79A", "MS4A1"]}

result = infer_cell_type_from_markers(counts, options={"markers": markers})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_cell_type_from_markers",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score(type, cell) = mean over the type's marker genes of the per-gene z-score of ln(1 + CP10k); the highest score wins if > min_score",
    "assignments": {
      "c1": "T cell",
      "c2": "T cell",
      "c3": "T cell",
      "c4": "B cell",
      "c5": "B cell",
      "c6": "B cell"
    },
    "best_score": {
      "c1": 0.8356997436691734,
      "c2": 0.8567110675970135,
      "c3": 0.8259936626183231,
      "c4": 0.8350643975343444,
      "c5": 0.8410178115070588,
      "c6": 0.8413837216201503
    },
    "counts": {
      "B cell": 3,
      "T cell": 3
    },
    "scores": {
      "T cell": {
        "c1": 0.8356997436691734,
        "c2": 0.8567110675970135,
        "c3": 0.8259936626183231,
        "c4": -0.5599901482596032,
        "c5": -0.5763799279190517,
        "c6": -1.3820343977058553
      },
      "B cell": {
        "c1": -0.5632664728908073,
        "c2": -0.5716923966489433,
        "c3": -1.3825070611218027,
        "c4": 0.8350643975343444,
        "c5": 0.8410178115070588,
        "c6": 0.8413837216201503
      }
    }
  }
}
```

## 27) infer_cnv

- **Purpose:** Infers copy-number-like expression profiles along chromosomes relative to reference cells (inferCNV-style, no HMM).

- **Method:** inferCNV-style expression profile without an HMM: counts are converted to `log2(CPM + 1)` (`target_sum` 1000000), the mean of the reference cells is subtracted per gene, the result is smoothed with a moving average over `window` genes ordered by position within each chromosome, and finally centred on its median per cell. Genes without an entry in `gene_positions` or with a mean raw count below `min_mean` are dropped; chromosomes are ordered naturally and genes by start position, and the moving-average window shrinks at chromosome ends. `cnv_score` of a cell is the mean absolute value of its modified expression. No CNV segments or gains/losses are called. Pure Python.

- **API:** `dnakit.singlecell_spatial.infer_cnv(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `gene_positions` (required; `{gene: [chromosome, start]}`), `groups` (required; `{cell: label}` or a list), `reference_group` (required; label of the normal cells), `window` (101 genes), `min_mean` (0.1), `target_sum` (1000000), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `window`, `reference_group`, `n_reference_cells`, `chromosome_of_gene`, `cnv_score` ({cell: score}), `cnv_score_definition`, `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the smoothed, centred profile).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import infer_cnv

counts = """gene,ref1,ref2,t1,t2
g1,10,12,10,11
g2,10,11,10,12
g3,10,10,20,22
g4,10,12,21,20
g5,10,11,10,12
g6,10,10,11,10
"""
options = {
    "gene_positions": {"g1": ["chr1", 100], "g2": ["chr1", 200], "g3": ["chr1", 300],
                       "g4": ["chr1", 400], "g5": ["chr2", 100], "g6": ["chr2", 200]},
    "groups": ["ref", "ref", "tumor", "tumor"],
    "reference_group": "ref",
    "window": 2,
}

result = infer_cnv(counts, options=options)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_cnv",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 7,
    "text_length": 111,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fde96bb232f2a9eec82e6b8896860c6131be1a2adbd848214e2b2f93e05d43ad"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "inferCNV-style expression profile: log2(CPM+1) minus the mean of the reference cells per gene, movin...",
    "window": 2,
    "reference_group": "ref",
    "n_reference_cells": 2,
    "chromosome_of_gene": {
      "g1": "chr1",
      "g2": "chr1",
      "g3": "chr1",
      "g4": "chr1",
      "g5": "chr2",
      "g6": "chr2"
    },
    "cnv_score": {
      "ref1": 0.016522693522969274,
      "ref2": 0.01652269352296878,
      "t1": 0.31512640468826447,
      "t2": 0.31465201165528606
    },
    "cnv_score_definition": "mean absolute value of the modified expression of the cell",
    "n_genes": 6,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2",
        "g3",
        "g4",
        "g5",
        "g6"
      ],
      "cells": [
        "ref1",
        "ref2",
        "t1",
        "t2"
      ],
      "values": [
        [
          -0.03387681047634672,
          0.03387681047634524,
          -0.2506093734737208,
          -0.17289386641218188
        ],
        [
          -0.0004988540582449019,
          0.0004988540582446057,
          0.1160999446066177,
          0.1728938664121819
        ],
        [
          -0.0004988540582449019,
          0.0004988540582446057,
          0.47289432167668544,
          0.4603909802784264
        ],
        [
          0.0004988540582449019,
          -0.0004988540582446062,
          0.8189548991593275,
          0.7299596507252003
        ],
        [
          0.03188139424336711,
          -0.03188139424336681,
          -0.11609994460661768,
          -0.17588685305186283
        ],
        [
          0.03188139424336711,
          -0.03188139424336681,
          -0.11609994460661768,
          -0.17588685305186283
        ]
      ]
    },
    "output_path": null
  }
}
```

## 28) infer_cnv_spatial

- **Purpose:** Infers copy-number-like expression profiles along chromosomes for spots relative to reference spots (inferCNV-style).

- **Method:** inferCNV-style expression profile without an HMM: counts are converted to `log2(CPM + 1)` (`target_sum` 1000000), the mean of the reference cells is subtracted per gene, the result is smoothed with a moving average over `window` genes ordered by position within each chromosome, and finally centred on its median per cell. Genes without an entry in `gene_positions` or with a mean raw count below `min_mean` are dropped; chromosomes are ordered naturally and genes by start position, and the moving-average window shrinks at chromosome ends. `cnv_score` of a cell is the mean absolute value of its modified expression. No CNV segments or gains/losses are called. Pure Python. Same definition as `infer_cnv`, with spots in place of cells.

- **API:** `dnakit.singlecell_spatial.infer_cnv_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `gene_positions` (required; `{gene: [chromosome, start]}`), `groups` (required; `{spot: label}` or a list), `reference_group` (required; label of the normal cells), `window` (101 genes), `min_mean` (0.1), `target_sum` (1000000), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `window`, `reference_group`, `n_reference_cells`, `chromosome_of_gene`, `cnv_score` ({cell: score}), `cnv_score_definition`, `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the smoothed, centred profile).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import infer_cnv_spatial

counts = """gene,ref1,ref2,t1,t2
g1,10,12,10,11
g2,10,11,10,12
g3,10,10,20,22
g4,10,12,21,20
g5,10,11,10,12
g6,10,10,11,10
"""
options = {
    "gene_positions": {"g1": ["chr1", 100], "g2": ["chr1", 200], "g3": ["chr1", 300],
                       "g4": ["chr1", 400], "g5": ["chr2", 100], "g6": ["chr2", 200]},
    "groups": ["ref", "ref", "tumor", "tumor"],
    "reference_group": "ref",
    "window": 2,
}

result = infer_cnv_spatial(counts, options=options)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_cnv_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 7,
    "text_length": 111,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fde96bb232f2a9eec82e6b8896860c6131be1a2adbd848214e2b2f93e05d43ad"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "inferCNV-style expression profile: log2(CPM+1) minus the mean of the reference cells per gene, movin...",
    "window": 2,
    "reference_group": "ref",
    "n_reference_cells": 2,
    "chromosome_of_gene": {
      "g1": "chr1",
      "g2": "chr1",
      "g3": "chr1",
      "g4": "chr1",
      "g5": "chr2",
      "g6": "chr2"
    },
    "cnv_score": {
      "ref1": 0.016522693522969274,
      "ref2": 0.01652269352296878,
      "t1": 0.31512640468826447,
      "t2": 0.31465201165528606
    },
    "cnv_score_definition": "mean absolute value of the modified expression of the cell",
    "n_genes": 6,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2",
        "g3",
        "g4",
        "g5",
        "g6"
      ],
      "cells": [
        "ref1",
        "ref2",
        "t1",
        "t2"
      ],
      "values": [
        [
          -0.03387681047634672,
          0.03387681047634524,
          -0.2506093734737208,
          -0.17289386641218188
        ],
        [
          -0.0004988540582449019,
          0.0004988540582446057,
          0.1160999446066177,
          0.1728938664121819
        ],
        [
          -0.0004988540582449019,
          0.0004988540582446057,
          0.47289432167668544,
          0.4603909802784264
        ],
        [
          0.0004988540582449019,
          -0.0004988540582446062,
          0.8189548991593275,
          0.7299596507252003
        ],
        [
          0.03188139424336711,
          -0.03188139424336681,
          -0.11609994460661768,
          -0.17588685305186283
        ],
        [
          0.03188139424336711,
          -0.03188139424336681,
          -0.11609994460661768,
          -0.17588685305186283
        ]
      ]
    },
    "output_path": null
  }
}
```

## 29) infer_trajectory

- **Purpose:** Computes geodesic pseudotime from a root cell through a k-nearest-neighbour graph of expression space.

- **Method:** Cells are represented by their centred ln(1 + CP10k) vectors (or by principal components when `n_pcs` > 0). A symmetric k-nearest-neighbour graph is built with Euclidean edge weights, and `pseudotime` is the shortest-path distance (Dijkstra) from the root divided by the largest finite distance. Cells not connected to the root get `null` and a warning. It is a graph-distance pseudotime, not a trajectory-inference method with branches or lineages. Requires numpy (returns `not_implemented` if it is missing); at most 5000 cells.

- **API:** `dnakit.singlecell_spatial.infer_trajectory(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `root` (required; cell name), `n_neighbors` (15), `n_pcs` (0 = use the centred vectors; > 0 = PCA first), `n_top_genes`, `scale`, `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `root`, `n_neighbors` and `pseudotime` ({cell: value from 0 to 1, `null` when not connected to the root}).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import infer_trajectory

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""

result = infer_trajectory(counts, options={"root": "c1", "n_neighbors": 3})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_trajectory",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "shortest-path (geodesic) distance from the root cell in a symmetric kNN graph of the cells' expressi...",
    "root": "c1",
    "n_neighbors": 3,
    "pseudotime": {
      "c1": 0.0,
      "c2": 0.600610963142177,
      "c3": 0.5933128032311472,
      "c4": 0.8570303182661878,
      "c5": 0.931692072807665,
      "c6": 1.0
    }
  }
}
```

## 30) infer_trajectory_spatial

- **Purpose:** Computes pseudotime of spots along the tissue: the shortest path from a root spot over a spatial neighbour graph.

- **Method:** A symmetric k-nearest-neighbour graph is built on the spot coordinates (`n_neighbors`, default 6); the weight of an edge is the Euclidean distance between the two spots' centred ln(1 + CP10k) expression vectors (or principal components when `n_pcs` > 0). `pseudotime` is the shortest-path length from the root spot divided by the maximum. Spots not connected to the root get `null` and a warning. Requires numpy (returns `not_implemented` if it is missing); at most 5000 spots.

- **API:** `dnakit.singlecell_spatial.infer_trajectory_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `root` (required; spot name), `coordinates` (required; `{spot: [x, y]}`, a list of `[x, y]` in column order, or a CSV/TSV/Visium `tissue_positions` file), `n_neighbors` (6; spatial neighbours), `n_pcs` (0), `n_top_genes`, `scale`, `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `root`, `n_neighbors` and `pseudotime` ({spot: value from 0 to 1, `null` when not connected to the root}).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import infer_trajectory_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
coordinates = {"c1": [0, 0], "c2": [1, 0], "c3": [2, 0], "c4": [0, 1], "c5": [1, 1], "c6": [2, 1]}

result = infer_trajectory_spatial(counts, options={"root": "c1", "coordinates": coordinates, "n_neighbors": 2})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_trajectory_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "shortest-path (geodesic) distance from the root cell in a symmetric kNN graph built on spatial coord...",
    "root": "c1",
    "n_neighbors": 2,
    "pseudotime": {
      "c1": 0.0,
      "c2": 0.28163935878865004,
      "c3": 0.5616005800468636,
      "c4": 0.40187989249504213,
      "c5": 0.7227424793087153,
      "c6": 1.0
    }
  }
}
```

## 31) infer_velocity

- **Purpose:** Estimates RNA velocity per gene and cell with the steady-state model from spliced and unspliced counts.

- **Method:** Steady-state model of La Manno et al. (2018): for every gene the degradation ratio `gamma = sum(u * s) / sum(s^2)` is fitted by least squares through the origin on the cells in the top `quantile` of spliced (`s`) or unspliced (`u`) expression of that gene, and `velocity = u - gamma * s`. Pass normalised layers (or `normalize: true`, which scales each layer to 10000 per cell). Velocities are not projected onto an embedding. Pure Python.

- **API:** `dnakit.singlecell_spatial.infer_velocity(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Spliced counts: a gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `unspliced` (required; matrix with the same genes and cells in the same order), `quantile` (0.05), `normalize` (false), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `quantile`, `gamma` ({gene: gamma}), `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the velocity matrix).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import infer_velocity

spliced = """gene,c1,c2,c3,c4
g1,1,2,3,4
g2,4,3,2,1
"""
unspliced = """gene,c1,c2,c3,c4
g1,2,3,5,6
g2,1,1,1,1
"""

result = infer_velocity(spliced, options={"unspliced": unspliced, "quantile": 0.5})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_velocity",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 39,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fad0b78cf99c22d8ad13bf8a0c0cc01340461e0a4940257b6b5569cf483f5dba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "steady-state model (La Manno et al. 2018): gamma_g = sum(u*s)/sum(s^2) over the cells whose s or u i...",
    "quantile": 0.5,
    "gamma": {
      "g1": 1.56,
      "g2": 0.3333333333333333
    },
    "n_genes": 2,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "c1",
        "c2",
        "c3",
        "c4"
      ],
      "values": [
        [
          0.43999999999999995,
          -0.1200000000000001,
          0.3200000000000003,
          -0.2400000000000002
        ],
        [
          -0.33333333333333326,
          0.0,
          0.33333333333333337,
          0.6666666666666667
        ]
      ]
    },
    "output_path": null
  }
}
```

## 32) infer_velocity_spatial

- **Purpose:** Estimates RNA velocity per gene and spot with the steady-state model from spliced and unspliced counts.

- **Method:** Steady-state model of La Manno et al. (2018): for every gene the degradation ratio `gamma = sum(u * s) / sum(s^2)` is fitted by least squares through the origin on the cells in the top `quantile` of spliced (`s`) or unspliced (`u`) expression of that gene, and `velocity = u - gamma * s`. Pass normalised layers (or `normalize: true`, which scales each layer to 10000 per cell). Velocities are not projected onto an embedding. Pure Python. Same model as `infer_velocity`.

- **API:** `dnakit.singlecell_spatial.infer_velocity_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Spliced counts: a gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `unspliced` (required; matrix with the same genes and cells in the same order), `quantile` (0.05), `normalize` (false), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `quantile`, `gamma` ({gene: gamma}), `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the velocity matrix).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import infer_velocity_spatial

spliced = """gene,s1,s2,s3,s4
g1,1,2,3,4
g2,4,3,2,1
"""
unspliced = """gene,s1,s2,s3,s4
g1,2,3,5,6
g2,1,1,1,1
"""

result = infer_velocity_spatial(spliced, options={"unspliced": unspliced, "quantile": 0.5})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_velocity_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 39,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "41d0822c1a43c1e406482fe3ce6fa42ad899944a4ee67f6be5a69ff110be1048"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "steady-state model (La Manno et al. 2018): gamma_g = sum(u*s)/sum(s^2) over the cells whose s or u i...",
    "quantile": 0.5,
    "gamma": {
      "g1": 1.56,
      "g2": 0.3333333333333333
    },
    "n_genes": 2,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "s1",
        "s2",
        "s3",
        "s4"
      ],
      "values": [
        [
          0.43999999999999995,
          -0.1200000000000001,
          0.3200000000000003,
          -0.2400000000000002
        ],
        [
          -0.33333333333333326,
          0.0,
          0.33333333333333337,
          0.6666666666666667
        ]
      ]
    },
    "output_path": null
  }
}
```

## 33) integrate

- **Purpose:** Integrates batches of cells by per-gene batch mean-centring of log-normalised expression.

- **Method:** Per-gene batch mean-centring: `x' = x - mean(gene in its batch) + mean(gene in all cells)`, computed on ln(1 + CP10k) values (`log_normalize` true, default) or on the supplied values. It aligns batch means only; it does not correct batch-specific variance or cell-type composition and is not Harmony, scVI or MNN. Pure Python.

- **API:** `dnakit.singlecell_spatial.integrate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `batches` (required; `{cell: batch label}` or a list in column order), `log_normalize` (true), `target_sum` (10000), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `batches` ({batch: number of cells}), `scale`, `note`, `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the corrected matrix).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import integrate

counts = """gene,a1,a2,b1,b2
g1,10,12,20,22
g2,5,6,5,7
"""

result = integrate(counts, options={"batches": ["A", "A", "B", "B"], "log_normalize": False})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "integrate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 43,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5f82895f5b4d8d8706f997e655451c6413349069dbf3487fe01fb0ba40faf310"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "per-gene batch mean-centring on the supplied values",
    "batches": {
      "A": 2,
      "B": 2
    },
    "scale": "as supplied",
    "note": "Aligns batch means only; it does not correct batch-specific variance or cell-type composition differences (not Harmony/scVI/MNN).",
    "n_genes": 2,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "a1",
        "a2",
        "b1",
        "b2"
      ],
      "values": [
        [
          15.0,
          17.0,
          15.0,
          17.0
        ],
        [
          5.25,
          6.25,
          4.75,
          6.75
        ]
      ]
    },
    "output_path": null
  }
}
```

## 34) integrate_single_cell_batches

- **Purpose:** Integrates single-cell batches by per-gene batch mean-centring of log-normalised expression.

- **Method:** Per-gene batch mean-centring: `x' = x - mean(gene in its batch) + mean(gene in all cells)`, computed on ln(1 + CP10k) values (`log_normalize` true, default) or on the supplied values. It aligns batch means only; it does not correct batch-specific variance or cell-type composition and is not Harmony, scVI or MNN. Pure Python.

- **API:** `dnakit.singlecell_spatial.integrate_single_cell_batches(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `batches` (required; `{cell: batch label}` or a list in column order), `log_normalize` (true), `target_sum` (10000), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `batches` ({batch: number of cells}), `scale`, `note`, `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the corrected matrix).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import integrate_single_cell_batches

counts = """gene,a1,a2,b1,b2
g1,10,12,20,22
g2,5,6,5,7
"""

result = integrate_single_cell_batches(counts, options={"batches": ["A", "A", "B", "B"], "log_normalize": False})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "integrate_single_cell_batches",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 43,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5f82895f5b4d8d8706f997e655451c6413349069dbf3487fe01fb0ba40faf310"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "per-gene batch mean-centring on the supplied values",
    "batches": {
      "A": 2,
      "B": 2
    },
    "scale": "as supplied",
    "note": "Aligns batch means only; it does not correct batch-specific variance or cell-type composition differences (not Harmony/scVI/MNN).",
    "n_genes": 2,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "a1",
        "a2",
        "b1",
        "b2"
      ],
      "values": [
        [
          15.0,
          17.0,
          15.0,
          17.0
        ],
        [
          5.25,
          6.25,
          4.75,
          6.75
        ]
      ]
    },
    "output_path": null
  }
}
```

## 35) integrate_spatial

- **Purpose:** Integrates spatial slices by per-gene batch mean-centring of log-normalised expression.

- **Method:** Per-gene batch mean-centring: `x' = x - mean(gene in its batch) + mean(gene in all cells)`, computed on ln(1 + CP10k) values (`log_normalize` true, default) or on the supplied values. It aligns batch means only; it does not correct batch-specific variance or cell-type composition and is not Harmony, scVI or MNN. Pure Python.

- **API:** `dnakit.singlecell_spatial.integrate_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `batches` (required; `{spot: batch label}` or a list in column order), `log_normalize` (true), `target_sum` (10000), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `batches` ({batch: number of cells}), `scale`, `note`, `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the corrected matrix).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import integrate_spatial

counts = """gene,a1,a2,b1,b2
g1,10,12,20,22
g2,5,6,5,7
"""

result = integrate_spatial(counts, options={"batches": ["slice1", "slice1", "slice2", "slice2"], "log_normalize": False})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "integrate_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 43,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5f82895f5b4d8d8706f997e655451c6413349069dbf3487fe01fb0ba40faf310"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "per-gene batch mean-centring on the supplied values",
    "batches": {
      "slice1": 2,
      "slice2": 2
    },
    "scale": "as supplied",
    "note": "Aligns batch means only; it does not correct batch-specific variance or cell-type composition differences (not Harmony/scVI/MNN).",
    "n_genes": 2,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "a1",
        "a2",
        "b1",
        "b2"
      ],
      "values": [
        [
          15.0,
          17.0,
          15.0,
          17.0
        ],
        [
          5.25,
          6.25,
          4.75,
          6.75
        ]
      ]
    },
    "output_path": null
  }
}
```

## 36) map_cells

- **Purpose:** Maps query cells onto an annotated reference by k-nearest-neighbour label transfer.

- **Method:** The query and the reference are restricted to their shared genes (at least 2 are required) and both are put on the ln(1 + CP10k) scale. Each query cell takes the majority label among its `n_neighbors` nearest reference cells by Euclidean distance (ties go to the label of the nearest tied neighbour); `confidence = winning votes / k`. This is plain k-nearest-neighbour label transfer, not Seurat anchors or scANVI. Requires numpy (returns `not_implemented` if it is missing); at most 5000 query and 5000 reference cells.

- **API:** `dnakit.singlecell_spatial.map_cells(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Query: a gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `reference` (required; annotated reference count matrix in any matrix form), `reference_labels` (required; `{cell: label}` or a list in the reference's column order), `n_neighbors` (5), `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_shared_genes` and `assignments`, a mapping query cell -> `label`, `confidence` and `nearest_reference_cell`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import map_cells

reference = """gene,r1,r2,r3,r4
CD3D,9,8,0,1
CD79A,0,1,9,8
ACTB,20,19,21,20
"""
query = """gene,q1,q2
CD3D,7,0
CD79A,1,10
ACTB,20,22
"""

result = map_cells(query, options={"reference": reference, "reference_labels": ["T", "T", "B", "B"], "n_neighbors": 3})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "map_cells",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1c42bc93471c1856987fb453c0089d4e8835c4e8d79ab188f49b6b39dd729a51"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "k-nearest-neighbour label transfer (k=3): Euclidean distance on ln(1+CP10k) values of the 3 shared g...",
    "n_shared_genes": 3,
    "assignments": {
      "q1": {
        "label": "T",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r2"
      },
      "q2": {
        "label": "B",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r3"
      }
    }
  }
}
```

## 37) map_cells_spatial

- **Purpose:** Maps spatial spots onto an annotated single-cell reference by k-nearest-neighbour label transfer.

- **Method:** The query and the reference are restricted to their shared genes (at least 2 are required) and both are put on the ln(1 + CP10k) scale. Each query cell takes the majority label among its `n_neighbors` nearest reference cells by Euclidean distance (ties go to the label of the nearest tied neighbour); `confidence = winning votes / k`. This is plain k-nearest-neighbour label transfer, not Seurat anchors or scANVI. Requires numpy (returns `not_implemented` if it is missing); at most 5000 query and 5000 reference cells. Spots of the query play the role of query cells.

- **API:** `dnakit.singlecell_spatial.map_cells_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Query spots: a gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `reference` (required; annotated reference count matrix in any matrix form), `reference_labels` (required; `{cell: label}` or a list in the reference's column order), `n_neighbors` (5), `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_shared_genes` and `assignments`, a mapping query cell -> `label`, `confidence` and `nearest_reference_cell`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import map_cells_spatial

reference = """gene,r1,r2,r3,r4
CD3D,9,8,0,1
CD79A,0,1,9,8
ACTB,20,19,21,20
"""
spots = """gene,s1,s2
CD3D,7,0
CD79A,1,10
ACTB,20,22
"""

result = map_cells_spatial(spots, options={"reference": reference, "reference_labels": ["T", "T", "B", "B"], "n_neighbors": 3})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "map_cells_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "27a116a639b79546f99b5ef15311b99dbe8acc233d94e46701791788b26a3da5"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "k-nearest-neighbour label transfer (k=3): Euclidean distance on ln(1+CP10k) values of the 3 shared g...",
    "n_shared_genes": 3,
    "assignments": {
      "s1": {
        "label": "T",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r2"
      },
      "s2": {
        "label": "B",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r3"
      }
    }
  }
}
```

## 38) map_single_cell_to_reference

- **Purpose:** Maps single cells onto an annotated reference atlas by k-nearest-neighbour label transfer.

- **Method:** The query and the reference are restricted to their shared genes (at least 2 are required) and both are put on the ln(1 + CP10k) scale. Each query cell takes the majority label among its `n_neighbors` nearest reference cells by Euclidean distance (ties go to the label of the nearest tied neighbour); `confidence = winning votes / k`. This is plain k-nearest-neighbour label transfer, not Seurat anchors or scANVI. Requires numpy (returns `not_implemented` if it is missing); at most 5000 query and 5000 reference cells.

- **API:** `dnakit.singlecell_spatial.map_single_cell_to_reference(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Query: a gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `reference` (required; annotated reference count matrix in any matrix form), `reference_labels` (required; `{cell: label}` or a list in the reference's column order), `n_neighbors` (5), `target_sum` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `n_shared_genes` and `assignments`, a mapping query cell -> `label`, `confidence` and `nearest_reference_cell`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import map_single_cell_to_reference

reference = """gene,r1,r2,r3,r4
CD3D,9,8,0,1
CD79A,0,1,9,8
ACTB,20,19,21,20
"""
query = """gene,q1,q2
CD3D,7,0
CD79A,1,10
ACTB,20,22
"""

result = map_single_cell_to_reference(
    query, options={"reference": reference, "reference_labels": ["T", "T", "B", "B"], "n_neighbors": 3}
)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "map_single_cell_to_reference",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1c42bc93471c1856987fb453c0089d4e8835c4e8d79ab188f49b6b39dd729a51"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "k-nearest-neighbour label transfer (k=3): Euclidean distance on ln(1+CP10k) values of the 3 shared g...",
    "n_shared_genes": 3,
    "assignments": {
      "q1": {
        "label": "T",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r2"
      },
      "q2": {
        "label": "B",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r3"
      }
    }
  }
}
```

## 39) normalize

- **Purpose:** Normalises a count matrix by library size, optionally followed by ln(1 + x).

- **Method:** Each cell (column) is scaled by its library size to `target_sum`; `method: "log1p"` (default) then applies ln(1 + x), `"total"` stops after the scaling and `"cpm"` scales to 1000000 instead. Cells with zero total counts stay all-zero. Pure Python.

- **API:** `dnakit.singlecell_spatial.normalize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `method` (`"log1p"` default, `"total"` or `"cpm"`), `target_sum` (10000), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `method`, `target_sum`, `library_sizes` ({cell: total counts}), `log_base` (`null` here), `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the normalised matrix).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import normalize

counts = """gene,c1,c2
g1,1,10
g2,3,30
"""

result = normalize(counts, options={"method": "total", "target_sum": 100})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "normalize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 27,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "94e87f829256b4e9de053af8b28b46acfc7b746db0aa5b9f2626d508d2123e26"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "total",
    "target_sum": 100.0,
    "library_sizes": {
      "c1": 4.0,
      "c2": 40.0
    },
    "log_base": null,
    "n_genes": 2,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "c1",
        "c2"
      ],
      "values": [
        [
          25.0,
          25.0
        ],
        [
          75.0,
          75.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 40) normalize_spatial

- **Purpose:** Normalises a spot-by-gene count matrix by library size, optionally followed by ln(1 + x).

- **Method:** Each cell (column) is scaled by its library size to `target_sum`; `method: "log1p"` (default) then applies ln(1 + x), `"total"` stops after the scaling and `"cpm"` scales to 1000000 instead. Cells with zero total counts stay all-zero. Pure Python. Same definition as `normalize`, with spots in place of cells.

- **API:** `dnakit.singlecell_spatial.normalize_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `method` (`"log1p"` default, `"total"` or `"cpm"`), `target_sum` (10000), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `method`, `target_sum`, `library_sizes` ({cell: total counts}), `log_base` (`null` here), `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the normalised matrix).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import normalize_spatial

counts = """gene,s1,s2
g1,1,10
g2,3,30
"""

result = normalize_spatial(counts, options={"method": "total", "target_sum": 100})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "normalize_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 27,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ee2c9bba2f5b8bbaa511f940886f514dbe9fa98bbca3baef34e40ab8a12dce54"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "total",
    "target_sum": 100.0,
    "library_sizes": {
      "s1": 4.0,
      "s2": 40.0
    },
    "log_base": null,
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
          25.0,
          25.0
        ],
        [
          75.0,
          75.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 41) plot_embedding

- **Purpose:** Draws a 2-D embedding (UMAP, t-SNE or PCA coordinates) as an SVG scatter plot.

- **Method:** One point per cell at the given coordinates. With `labels` the points are coloured by category with a legend (at most 20 labels are listed; `labels` take precedence over `values`); with `values` a continuous colour ramp is used; without either all points share one colour. The SVG is built as text without any plotting library; it is also written to `output_path` when given (the name must end in `.svg`; an existing file raises `ConfigurationError`).

- **API:** `dnakit.singlecell_spatial.plot_embedding(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Coordinates with one row per cell: `{cell: [x, y]}`, a list of `[x, y]`, CSV/TSV text or a file (header = dimension names, first column = cell), or a numpy array. Options: `labels` (`{cell: label}` or a list; categorical colours and legend), `values` (`{cell: number}` or a list; continuous colour ramp), `output_path` (optional `.svg` file to write; never overwrites), `title`.

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `format` (`"svg"`), `svg` (the SVG text), `output_path` (the `.svg` file written, else `null`), `n_points`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import plot_embedding

embedding = {"c1": [0.0, 0.1], "c2": [0.2, 0.0], "c3": [1.0, 1.1], "c4": [1.2, 0.9]}

result = plot_embedding(embedding, options={"labels": ["T", "T", "B", "B"], "title": "UMAP"})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "plot_embedding",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 72,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "29c10aeb08c17f75aa107412cf896fb1568539d195e831ede3be015149e5ab12"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "svg",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"480\" viewBox=\"0 0 720 480\"><rect width=\"...",
    "n_points": 4,
    "output_path": null
  }
}
```

## 42) plot_embedding_spatial

- **Purpose:** Draws the spot coordinates of a spatial dataset as an SVG scatter plot coloured by label or value.

- **Method:** One point per spot at its tissue coordinates (the y axis points down, like an image). With `labels` the points are coloured by category with a legend (`labels` take precedence over `values`); with `values` a continuous colour ramp is used. The SVG is built as text without any plotting library; it is also written to `output_path` when given (the name must end in `.svg`; an existing file raises `ConfigurationError`).

- **API:** `dnakit.singlecell_spatial.plot_embedding_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Spot coordinates with one row per spot: `{spot: [x, y]}`, a list of `[x, y]`, CSV/TSV text or a file, or a numpy array. Options: `labels`, `values`, `output_path` (optional `.svg` file to write; never overwrites), `title`.

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `format` (`"svg"`), `svg` (the SVG text), `output_path` (the `.svg` file written, else `null`), `n_points`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import plot_embedding_spatial

spots = {"s1": [0, 0], "s2": [1, 0], "s3": [0, 1], "s4": [1, 1]}

result = plot_embedding_spatial(spots, options={"values": [1.0, 2.0, 3.0, 4.0], "title": "Spots"})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "plot_embedding_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "314985b9143128b56c67c5110e2e2e99c929a010e0138a149feb01de78c4bcd6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "svg",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"480\" viewBox=\"0 0 720 480\"><rect width=\"...",
    "n_points": 4,
    "output_path": null
  }
}
```

## 43) plot_heatmap

- **Purpose:** Draws a gene-by-cell expression heatmap as an SVG.

- **Method:** The `max_genes` genes with the highest variance across the shown cells and the first `max_cells` cells are drawn; values are `log1p(x)` by default (`log1p`), optionally z-scored per gene (`row_zscore`). The SVG is built as text without any plotting library; it is also written to `output_path` when given (the name must end in `.svg`; an existing file raises `ConfigurationError`).

- **API:** `dnakit.singlecell_spatial.plot_heatmap(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `max_genes` (50), `max_cells` (200), `log1p` (true), `row_zscore` (false), `output_path` (optional `.svg` file to write; never overwrites), `title`, `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `format` (`"svg"`), `svg` (the SVG text), `output_path` (the `.svg` file written, else `null`), `genes_shown`, `n_cells_shown` and `selection`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import plot_heatmap

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""

result = plot_heatmap(counts, options={"max_genes": 4, "title": "Expression"})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "plot_heatmap",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "svg",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"230\" height=\"140\" viewBox=\"0 0 230 140\"><rect width=\"...",
    "genes_shown": [
      "CD79A",
      "CD3D",
      "MS4A1",
      "CD3E"
    ],
    "n_cells_shown": 6,
    "selection": "genes with the highest variance across the shown cells",
    "output_path": null
  }
}
```

## 44) plot_heatmap_spatial

- **Purpose:** Draws the spatial expression map of one gene (spots coloured by ln(1 + CP10k)) as an SVG.

- **Method:** Each spot is drawn at its coordinates and coloured by the gene's ln(1 + CP10k) value (`target_sum` 10000) on a continuous colour ramp. The SVG is built as text without any plotting library; it is also written to `output_path` when given (the name must end in `.svg`; an existing file raises `ConfigurationError`).

- **API:** `dnakit.singlecell_spatial.plot_heatmap_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `gene` (required), `coordinates` (required; `{spot: [x, y]}`, a list of `[x, y]` in column order, or a CSV/TSV/Visium `tissue_positions` file), `target_sum` (10000), `output_path` (optional `.svg` file to write; never overwrites), `title`, `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `format` (`"svg"`), `svg` (the SVG text), `output_path` (the `.svg` file written, else `null`), `gene`, `n_spots`, `min_value` and `max_value`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import plot_heatmap_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
coordinates = {"c1": [0, 0], "c2": [1, 0], "c3": [2, 0], "c4": [0, 1], "c5": [1, 1], "c6": [2, 1]}

result = plot_heatmap_spatial(counts, options={"gene": "CD3D", "coordinates": coordinates})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "plot_heatmap_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "svg",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"480\" viewBox=\"0 0 720 480\"><rect width=\"...",
    "gene": "CD3D",
    "n_spots": 6,
    "min_value": 0.0,
    "max_value": 7.358593577814888,
    "output_path": null
  }
}
```

## 45) qc

- **Purpose:** Computes per-cell and per-gene quality-control metrics of a count matrix and flags cells that fail thresholds.

- **Method:** Per cell: total counts, number of genes detected (count > 0) and percent mitochondrial counts (counts of genes whose name starts with a `mito_prefix` divided by the total). Cells failing any given threshold are listed with the names of the failed checks, and genes detected in fewer than `min_cells` cells are listed. Thresholds that are not given are only reported. Pure Python.

- **API:** `dnakit.singlecell_spatial.qc(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `min_genes`, `min_counts`, `max_counts`, `max_pct_mito`, `min_cells` (all default `null` = report only), `mito_prefix` (`"MT-"` and `"mt-"`), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `n_cells`, `n_genes`, `n_mito_genes`, `median_total_counts`, `median_genes_detected`, `median_pct_mito`, `genes_detected_in_any_cell`, `thresholds`, `n_cells_failing`, `failing_cells`, `n_genes_below_min_cells`, `genes_below_min_cells` and `per_cell` (`cell`, `total_counts`, `n_genes_detected`, `pct_mito`, `fails`) with `per_cell_truncated`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import qc

matrix = {"MT-A": [2, 0], "B": [4, 5], "C": [4, 5]}

result = qc(matrix, options={"cells": ["x", "y"], "max_pct_mito": 10})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "qc",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "edd1c156fa59a3ff619c27e7b9283d63d5a49d9dfb445ea4617ec6867f80b124"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_cells": 2,
    "n_genes": 3,
    "n_mito_genes": 1,
    "median_total_counts": 10.0,
    "median_genes_detected": 2.5,
    "median_pct_mito": 10.0,
    "genes_detected_in_any_cell": 3,
    "thresholds": {
      "min_genes": null,
      "min_counts": null,
      "max_counts": null,
      "max_pct_mito": 10.0,
      "min_cells": null
    },
    "n_cells_failing": 1,
    "failing_cells": [
      "x"
    ],
    "n_genes_below_min_cells": 0,
    "genes_below_min_cells": [],
    "per_cell": [
      {
        "cell": "x",
        "total_counts": 10.0,
        "n_genes_detected": 3,
        "pct_mito": 20.0,
        "fails": [
          "max_pct_mito"
        ]
      },
      {
        "cell": "y",
        "total_counts": 10.0,
        "n_genes_detected": 2,
        "pct_mito": 0.0,
        "fails": []
      }
    ],
    "per_cell_truncated": false
  }
}
```

## 46) qc_spatial

- **Purpose:** Computes per-spot and per-gene quality-control metrics of a spatial count matrix and flags spots that fail thresholds.

- **Method:** Per cell: total counts, number of genes detected (count > 0) and percent mitochondrial counts (counts of genes whose name starts with a `mito_prefix` divided by the total). Cells failing any given threshold are listed with the names of the failed checks, and genes detected in fewer than `min_cells` cells are listed. Thresholds that are not given are only reported. Pure Python. Optional `coordinates` add the number of spots with coordinates and the coordinate ranges.

- **API:** `dnakit.singlecell_spatial.qc_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `min_genes`, `min_counts`, `max_counts`, `max_pct_mito`, `min_cells` (all default `null` = report only), `mito_prefix` (`"MT-"` and `"mt-"`), `coordinates` (optional; `{spot: [x, y]}`, a list of `[x, y]` in column order, or a CSV/TSV/Visium `tissue_positions` file), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `n_cells`, `n_genes`, `n_mito_genes`, `median_total_counts`, `median_genes_detected`, `median_pct_mito`, `genes_detected_in_any_cell`, `thresholds`, `n_cells_failing`, `failing_cells`, `n_genes_below_min_cells`, `genes_below_min_cells` and `per_cell` (`cell`, `total_counts`, `n_genes_detected`, `pct_mito`, `fails`) with `per_cell_truncated`; with `coordinates` also `n_spots_with_coordinates`, `x_range` and `y_range`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import qc_spatial

matrix = {"MT-A": [2, 0], "B": [4, 5], "C": [4, 5]}
spots = {"x": [0, 0], "y": [3, 4]}

result = qc_spatial(matrix, options={"cells": ["x", "y"], "max_pct_mito": 10, "coordinates": spots})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "qc_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "edd1c156fa59a3ff619c27e7b9283d63d5a49d9dfb445ea4617ec6867f80b124"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_cells": 2,
    "n_genes": 3,
    "n_mito_genes": 1,
    "median_total_counts": 10.0,
    "median_genes_detected": 2.5,
    "median_pct_mito": 10.0,
    "genes_detected_in_any_cell": 3,
    "thresholds": {
      "min_genes": null,
      "min_counts": null,
      "max_counts": null,
      "max_pct_mito": 10.0,
      "min_cells": null
    },
    "n_cells_failing": 1,
    "failing_cells": [
      "x"
    ],
    "n_genes_below_min_cells": 0,
    "genes_below_min_cells": [],
    "per_cell": [
      {
        "cell": "x",
        "total_counts": 10.0,
        "n_genes_detected": 3,
        "pct_mito": 20.0,
        "fails": [
          "max_pct_mito"
        ]
      },
      {
        "cell": "y",
        "total_counts": 10.0,
        "n_genes_detected": 2,
        "pct_mito": 0.0,
        "fails": []
      }
    ],
    "per_cell_truncated": false,
    "n_spots_with_coordinates": 2,
    "x_range": [
      0.0,
      3.0
    ],
    "y_range": [
      0.0,
      4.0
    ]
  }
}
```

## 47) read

- **Purpose:** Reads a count matrix (mapping, rows, CSV/TSV, array or 10x directory) and summarises its shape and sparsity.

- **Method:** The shared matrix parser accepts all listed input forms and gives the same genes x cells matrix for each. `sparsity = 1 - nonzero entries / (genes * cells)`. The matrix values are returned when the matrix has at most `max_return_elements` elements, otherwise they are omitted with a warning. Pure Python.

- **API:** `dnakit.singlecell_spatial.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `max_return_elements` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `n_genes`, `n_cells`, `total_counts`, `nonzero_entries`, `sparsity`, `genes_preview`, `cells_preview` (the first 10 names), `orientation` and, when small enough, `matrix` (`orientation`, `genes`, `cells`, `values`).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import read

matrix = {"g1": [1, 0, 3], "g2": [0, 0, 2]}

result = read(matrix, options={"cells": ["c1", "c2", "c3"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 34,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ac109c8082501eb269579a666d0f9049f1507060d522b056e4cd9c762cdf5f4e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_genes": 2,
    "n_cells": 3,
    "total_counts": 6.0,
    "nonzero_entries": 3,
    "sparsity": 0.5,
    "genes_preview": [
      "g1",
      "g2"
    ],
    "cells_preview": [
      "c1",
      "c2",
      "c3"
    ],
    "orientation": "genes_by_cells",
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "c1",
        "c2",
        "c3"
      ],
      "values": [
        [
          1.0,
          0.0,
          3.0
        ],
        [
          0.0,
          0.0,
          2.0
        ]
      ]
    }
  }
}
```

## 48) read_spatial

- **Purpose:** Reads a spatial count matrix together with spot coordinates and summarises it.

- **Method:** The count matrix is parsed as in `read`. Coordinates are matched to the spots by name (or taken in column order for a list); a headerless Visium v1 `tissue_positions_list` file is recognised (x = pixel column, y = pixel row). Without coordinates only the count matrix is read and a warning is given. Pure Python.

- **API:** `dnakit.singlecell_spatial.read_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `coordinates` (optional; `{spot: [x, y]}`, a list of `[x, y]`, or a CSV/TSV/Visium `tissue_positions` file), `x_column` and `y_column` (header names), `max_return_elements` (10000), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: the fields of `read` plus, with coordinates, `coordinates` ({spot: [x, y]}), `n_spots_with_coordinates`, `x_range` and `y_range`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import read_spatial

matrix = {"g1": [1, 0, 3], "g2": [0, 0, 2]}
spots = {"c1": [0, 0], "c2": [10, 0], "c3": [10, 5]}

result = read_spatial(matrix, options={"cells": ["c1", "c2", "c3"], "coordinates": spots})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "read_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 34,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ac109c8082501eb269579a666d0f9049f1507060d522b056e4cd9c762cdf5f4e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_genes": 2,
    "n_cells": 3,
    "total_counts": 6.0,
    "nonzero_entries": 3,
    "sparsity": 0.5,
    "genes_preview": [
      "g1",
      "g2"
    ],
    "cells_preview": [
      "c1",
      "c2",
      "c3"
    ],
    "orientation": "genes_by_cells",
    "coordinates": {
      "c1": [
        0.0,
        0.0
      ],
      "c2": [
        10.0,
        0.0
      ],
      "c3": [
        10.0,
        5.0
      ]
    },
    "n_spots_with_coordinates": 3,
    "x_range": [
      0.0,
      10.0
    ],
    "y_range": [
      0.0,
      5.0
    ],
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "c1",
        "c2",
        "c3"
      ],
      "values": [
        [
          1.0,
          0.0,
          3.0
        ],
        [
          0.0,
          0.0,
          2.0
        ]
      ]
    }
  }
}
```

## 49) remove_ambient

- **Purpose:** Subtracts ambient-RNA counts from a droplet count matrix using an ambient profile and a contamination fraction.

- **Method:** Simplified subtractive, SoupX-style correction with a user-given contamination fraction: `corrected = max(0, count - contamination * cell_total * ambient_fraction(gene))`. The ambient profile is given directly (`ambient_profile`, normalised to sum 1), or estimated as the summed counts of the `empty_cells` columns, or of the columns whose total is at most `empty_max_counts`; the profile is used in this priority order, and with the two empty-droplet options those empty droplets are removed from the output. The contamination fraction is not estimated, and the removed counts are not redistributed as integers (`round` merely rounds the corrected values). Pure Python.

- **API:** `dnakit.singlecell_spatial.remove_ambient(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `contamination` (required; fraction in [0, 1) of each cell's counts that is ambient), `ambient_profile` (`{gene: weight}`) OR `empty_cells` (list of empty-droplet columns) OR `empty_max_counts` (columns with total <= this are empty droplets and are removed from the output), `round` (false; round the corrected values), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `contamination`, `ambient_source`, `ambient_profile` ({gene: fraction}), `counts_removed`, `excluded_empty_droplets`, `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the corrected matrix).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import remove_ambient

counts = """gene,c1,c2
g1,10,0
g2,90,100
"""

result = remove_ambient(counts, options={"contamination": 0.1, "ambient_profile": {"g1": 1, "g2": 9}})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "remove_ambient",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 29,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fb38f8f48b60feca759bbab925123eb87142fcf8fa7b1b5f22563bed973635e4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "corrected = max(0, count - contamination * cell_total * ambient_fraction(gene)); simplified subtract...",
    "contamination": 0.1,
    "ambient_source": "ambient_profile",
    "ambient_profile": {
      "g1": 0.1,
      "g2": 0.9
    },
    "counts_removed": 19.0,
    "excluded_empty_droplets": [],
    "n_genes": 2,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "c1",
        "c2"
      ],
      "values": [
        [
          9.0,
          0.0
        ],
        [
          81.0,
          91.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 50) remove_ambient_spatial

- **Purpose:** Subtracts ambient-RNA counts from a spatial count matrix using an ambient profile and a contamination fraction.

- **Method:** Simplified subtractive, SoupX-style correction with a user-given contamination fraction: `corrected = max(0, count - contamination * cell_total * ambient_fraction(gene))`. The ambient profile is given directly (`ambient_profile`, normalised to sum 1), or estimated as the summed counts of the `empty_cells` columns, or of the columns whose total is at most `empty_max_counts`; the profile is used in this priority order, and with the two empty-droplet options those empty droplets are removed from the output. The contamination fraction is not estimated, and the removed counts are not redistributed as integers (`round` merely rounds the corrected values). Pure Python. Same correction as `remove_ambient`, with spots in place of cells.

- **API:** `dnakit.singlecell_spatial.remove_ambient_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `contamination` (required; fraction in [0, 1) of each cell's counts that is ambient), `ambient_profile` (`{gene: weight}`) OR `empty_cells` (list of empty-droplet columns) OR `empty_max_counts` (columns with total <= this are empty droplets and are removed from the output), `round` (false; round the corrected values), `return_matrix` (true), `output_path` (optional TSV file to write; never overwrites), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `definition`, `contamination`, `ambient_source`, `ambient_profile` ({gene: fraction}), `counts_removed`, `excluded_empty_droplets`, `n_genes`, `n_cells`, `matrix` (`orientation`, `genes`, `cells`, `values`; omitted when `return_matrix` is false) and `output_path` (the TSV file written, else `null`) (the corrected matrix).

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import remove_ambient_spatial

counts = """gene,s1,s2
g1,10,0
g2,90,100
"""

result = remove_ambient_spatial(counts, options={"contamination": 0.1, "ambient_profile": {"g1": 1, "g2": 9}})
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "remove_ambient_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 29,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8bd3aebc702036c3c3c2c8bce05e801f0123aa42609892abb4c79e24c6682b4d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "corrected = max(0, count - contamination * cell_total * ambient_fraction(gene)); simplified subtract...",
    "contamination": 0.1,
    "ambient_source": "ambient_profile",
    "ambient_profile": {
      "g1": 0.1,
      "g2": 0.9
    },
    "counts_removed": 19.0,
    "excluded_empty_droplets": [],
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
          9.0,
          0.0
        ],
        [
          81.0,
          91.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 51) remove_doublets

- **Purpose:** Intended to remove doublets from a count matrix; it is not implemented and always returns `not_implemented`.

- **Method:** No computation is performed. Doublet detection needs a simulation-and-classifier implementation (for example Scrublet or DoubletFinder): artificial doublets are simulated and a trained classifier or neighbour model scores each cell. DNAKit does not call any such package, so installing one changes nothing, and no deterministic reference definition is bundled; so the function exists only to keep its name and call signature stable.

- **API:** `dnakit.singlecell_spatial.remove_doublets(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any input is accepted (it is only summarised in `input_summary`); a count matrix is the intended input. `options` is optional and unused.

- **Output:** Always the envelope `status: "not_implemented"` with `module`, `operation`, `provenance`, `input_summary`, `warnings`, `required_dependency` (a description of what is missing, not an install command: a doublet-detection implementation based on simulation and a classifier, which DNAKit does not call) and `message` (why it is not implemented). There is no `result`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import remove_doublets

counts = """gene,c1,c2
g1,10,0
g2,90,100
"""

result = remove_doublets(counts)
print(json.dumps(result, indent=2))
```

- **Example output:** (long lists and strings are shortened with `...`)

```json
{
  "module": "singlecell_spatial",
  "operation": "remove_doublets",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 29,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fb38f8f48b60feca759bbab925123eb87142fcf8fa7b1b5f22563bed973635e4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a doublet-detection implementation based on simulation and a classifier (for example Scrublet or DoubletFinder); DNAKit does not call any such package",
  "message": "Doublet detection needs a simulation-and-classifier implementation (for example Scrublet or DoubletFinder) which DNAKit does not call; no deterministic reference definition is bundled."
}
```

## 52) remove_doublets_spatial

- **Purpose:** Intended to remove multiplets from spatial spot data; it is not implemented and always returns `not_implemented`.

- **Method:** No computation is performed. Spot-level doublet detection needs a simulation-and-classifier implementation (for example Scrublet or DoubletFinder) which DNAKit does not call, so installing one changes nothing; none is bundled, so the function exists only to keep its name and call signature stable.

- **API:** `dnakit.singlecell_spatial.remove_doublets_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Any input is accepted (it is only summarised in `input_summary`); a spot count matrix is the intended input. `options` is optional and unused.

- **Output:** Always the envelope `status: "not_implemented"` with `module`, `operation`, `provenance`, `input_summary`, `warnings`, `required_dependency` (a description of what is missing, not an install command: a doublet-detection implementation based on simulation and a classifier, which DNAKit does not call) and `message` (why it is not implemented). There is no `result`.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import remove_doublets_spatial

counts = """gene,s1,s2
g1,10,0
g2,90,100
"""

result = remove_doublets_spatial(counts)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "remove_doublets_spatial",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 29,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8bd3aebc702036c3c3c2c8bce05e801f0123aa42609892abb4c79e24c6682b4d"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a doublet-detection implementation based on simulation and a classifier (for example Scrublet or DoubletFinder); DNAKit does not call any such package",
  "message": "Spot-level doublet detection needs a simulation-and-classifier implementation (for example Scrublet or DoubletFinder) which DNAKit does not call; none is bundled."
}
```

## 53) validate

- **Purpose:** Validates a count matrix: shape, numeric values, non-negativity, integer counts, unique names and empty cells/genes.

- **Method:** Errors (make `valid` false): the matrix cannot be parsed, duplicate gene or cell names, negative values, and (with `require_integer`, default true) non-integer values. Warnings: cells with zero total counts and genes that are zero in every cell. Pure Python.

- **API:** `dnakit.singlecell_spatial.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`). Options: `require_integer` (true; non-integer values are errors), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: `n_genes`, `n_cells`, `empty_cells`, `empty_genes`, `valid`, `errors` and `warnings` (lists of messages); a matrix that cannot be parsed gives `valid: false` with the parse error instead of raising.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import validate

result = validate("gene,c1,c2\na,1,2\na,3,4\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 23,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f6eb24e753294ea677e9b11c564784e44b668c8f1c28050fd650d961bb38203e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_genes": 2,
    "n_cells": 2,
    "empty_cells": 0,
    "empty_genes": 0,
    "valid": false,
    "errors": [
      "duplicate gene names: ['a']"
    ],
    "warnings": []
  }
}
```

## 54) validate_spatial

- **Purpose:** Validates a spatial count matrix and its spot coordinates (every spot needs finite coordinates).

- **Method:** The checks of `validate`, plus: `coordinates` must be given and every spot must have finite coordinates, otherwise an error is added. Pure Python.

- **API:** `dnakit.singlecell_spatial.validate_spatial(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** A gene-by-cell count matrix (columns are cells): a mapping `{gene: [counts per cell]}` or `{gene: {cell: count}}`, a list of rows, CSV/TSV text or a file path (header = cell names, first column = gene names), a numpy array or pandas DataFrame, or a 10x directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`); spots are treated as cells, so the columns are spots. Options: `coordinates` (required; `{spot: [x, y]}`, a list of `[x, y]` in column order, or a CSV/TSV/Visium `tissue_positions` file), `require_integer` (true), `orientation` (`"genes_by_cells"` default, or `"cells_by_genes"` when rows are cells), `genes` and `cells` (name lists for unnamed forms), `max_matrix_elements` (5000000).

- **Output:** The result envelope has `status` (`"ok"`), `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`. Malformed input or a missing required option raises `InputFormatError` / `ConfigurationError` instead of returning an envelope. `result`: the fields of `validate` plus, when the coordinates are valid, `n_spots_with_coordinates`, `x_range` and `y_range`; missing coordinates are reported as errors.

- **Example code:**

```python
import json

from dnakit.singlecell_spatial import validate_spatial

matrix = {"g1": [1, 0, 3], "g2": [0, 0, 2]}
spots = {"c1": [0, 0], "c2": [1, 1]}

result = validate_spatial(matrix, options={"cells": ["c1", "c2", "c3"], "coordinates": spots})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "singlecell_spatial",
  "operation": "validate_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 34,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ac109c8082501eb269579a666d0f9049f1507060d522b056e4cd9c762cdf5f4e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_genes": 2,
    "n_cells": 3,
    "empty_cells": 1,
    "empty_genes": 0,
    "valid": false,
    "errors": [
      "coordinates: [MISSING_COORDINATES] 1 cells have no coordinates (first: 'c3'). Hint: Provide coordinates for every cell/spot in the matrix."
    ],
    "warnings": [
      "1 cells have zero total counts (first: 'c2')"
    ]
  }
}
```
