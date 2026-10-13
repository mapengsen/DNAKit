# Transcriptomics single-cell backends: AnnData, scanpy, muon, and scVI

This group of functions in `dnakit.transcriptomics` is a set of gene-level single-cell adapters over optional packages. [anndata](https://github.com/scverse/anndata) reads and writes `.h5ad`, [scanpy](https://github.com/scverse/scanpy) computes cell and gene QC, highly variable genes, marker genes, and pseudobulk sums, [muon](https://github.com/scverse/muon) summarizes multi-modal `.h5mu` files, and [scvi-tools](https://github.com/scverse/scvi-tools) sets up and trains an scVI model. None of these functions performs clustering, dimensionality reduction, cell-type annotation, or trajectory analysis, and none downloads public data. Each one works only on an explicit file.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command.

## Example data

The examples below run in order and share one working directory. The first example writes `counts.csv` as `tiny.h5ad`, and later examples reuse it. Any other input file is created by the example code that needs it. The data only demonstrates the interface, and its statistics do not describe real cells.

**Example file:** `counts.csv`

```text
cell,group,sample,MT-A,GENE1,GENE2,RPL3
c0,x,s1,3,0,1,0
c1,x,s2,0,2,0,0
c2,y,s1,1,1,1,0
c3,y,s2,0,0,5,2
```

## 1) `dnakit.transcriptomics.read_anndata_summary()`

- **Purpose:** Read the dimensions, cell and gene metadata column names, the first few genes' metadata, and the names of attached data, without loading the count matrix.
- **Method:** Open the file in backed mode and read only the structure of `obs`, `var`, `layers`, `obsm`, `varm`, and `uns`. The matrix contributes only its dtype. A file larger than `max_file_bytes` raises an error.
- **API:** `dnakit.transcriptomics.read_anndata_summary(source, *, max_file_bytes=1073741824)`.
- **Input:** The path to an `.h5ad` file. `max_file_bytes` is the largest file size in bytes that will be read.
- **Output:** A dictionary with `n_obs`, `n_vars`, `x_dtype`, `obs_columns`, `var_columns`, `var_preview`, `layers`, `obsm`, `varm`, `uns`, and `provenance`. `var_preview` lists at most five genes with their metadata.
- **Example code:**

```python
import anndata as ad
import pandas as pd
from dnakit.transcriptomics import read_anndata_summary

table = pd.read_csv("counts.csv", index_col="cell")
var = pd.DataFrame({"chromosome": ["chrM", "chr1", "chr2", "chr3"]}, index=table.columns[2:])
adata = ad.AnnData(X=table.iloc[:, 2:].to_numpy(dtype="float32"), obs=table[["group", "sample"]], var=var)
adata.write_h5ad("tiny.h5ad")

result = read_anndata_summary("tiny.h5ad")
print(result["n_obs"], result["n_vars"], result["obs_columns"])
print(result["var_preview"][0])
```

- **Example output:**

```text
4 4 ['group', 'sample']
{'gene': 'MT-A', 'chromosome': 'chrM'}
```

## 2) `dnakit.transcriptomics.calculate_gene_qc_metrics()`

- **Purpose:** Compute, for each gene, the number of detecting cells, total counts, mean counts, and dropout rate, and list the genes with the largest totals.
- **Method:** Call `scanpy.pp.calculate_qc_metrics` and use its gene table. `n_cells_by_counts` is the number of cells with a count above 0, and `pct_dropout_by_counts` is the percentage of cells in which the gene is not detected. Genes with equal totals are sorted by name.
- **API:** `dnakit.transcriptomics.calculate_gene_qc_metrics(source, *, top_n=10, max_file_bytes=1073741824, max_matrix_elements=50000000)`.
- **Input:** The path to an `.h5ad` file. `top_n` is the number of genes to list (1 to 1000). `max_matrix_elements` limits cells times genes, and exceeding it raises an error.
- **Output:** `n_genes`, `genes_detected`, `median_cells_per_gene`, `top_genes` (for each gene: total counts, detecting cells, mean counts, and dropout rate), and `provenance`.
- **Example code:**

```python
from dnakit.transcriptomics import calculate_gene_qc_metrics

result = calculate_gene_qc_metrics("tiny.h5ad", top_n=2)
print(result["genes_detected"], result["median_cells_per_gene"])
for gene in result["top_genes"]:
    print(gene["gene"], gene["total_counts"], gene["n_cells_by_counts"])
```

- **Example output:**

```text
4 2.0
GENE2 7.0 3
MT-A 4.0 2
```

## 3) `dnakit.transcriptomics.calculate_cell_qc_metrics()`

- **Purpose:** Summarize each cell's gene count, total counts, and mitochondrial count fraction, and count the cells above a threshold.
- **Method:** Genes whose names start with `mito_prefix` are treated as mitochondrial (`"MT-"` is common for human data and `"mt-"` for mouse). Call `scanpy.pp.calculate_qc_metrics(qc_vars=["mt"])` to get `pct_counts_mt` for each cell, which is the mitochondrial count as a percentage of that cell's total counts. The function then reports the median, the maximum, and the number of cells above the threshold. The result contains only these summary values, not a per-cell table.
- **API:** `dnakit.transcriptomics.calculate_cell_qc_metrics(source, *, mito_prefix="MT-", mito_threshold_percent=20.0, max_file_bytes=1073741824, max_matrix_elements=50000000)`.
- **Input:** The path to an `.h5ad` file. `mito_threshold_percent` is a number in (0, 100], expressed as a percentage.
- **Output:** `n_cells`, `n_genes`, `mito_genes`, `median_genes_per_cell`, `median_total_counts`, `median_pct_counts_mt`, `max_pct_counts_mt`, `cells_above_mito_threshold`, and `provenance`.
- **Example code:**

```python
from dnakit.transcriptomics import calculate_cell_qc_metrics

result = calculate_cell_qc_metrics("tiny.h5ad", mito_threshold_percent=20)
print(result["mito_genes"], result["median_total_counts"], result["cells_above_mito_threshold"])
print(result["max_pct_counts_mt"])
```

- **Example output:**

```text
1 3.5 2
75.0
```

## 4) `dnakit.transcriptomics.filter_genes_by_min_cells()`

- **Purpose:** Remove genes detected in fewer than `min_cells` cells, and optionally write the result to a new `.h5ad`.
- **Method:** Call `scanpy.pp.filter_genes(min_cells=...)`. A gene is "detected" in a cell when its count is above 0. If the file at `output_path` already exists, the function raises an error and does not overwrite it.
- **API:** `dnakit.transcriptomics.filter_genes_by_min_cells(source, *, min_cells=3, output_path=None, max_file_bytes=1073741824, max_matrix_elements=50000000)`.
- **Input:** The path to an `.h5ad` file. `min_cells` is a positive integer. `output_path` is optional and must end with `.h5ad`.
- **Output:** `n_cells`, `n_genes_before`, `n_genes_after`, `n_genes_removed`, `output_path`, and `provenance`.
- **Example code:**

```python
from pathlib import Path
import anndata as ad
from dnakit.transcriptomics import filter_genes_by_min_cells

Path("filtered.h5ad").unlink(missing_ok=True)  # the function never overwrites an existing file
result = filter_genes_by_min_cells("tiny.h5ad", min_cells=2, output_path="filtered.h5ad")
print(result["n_genes_before"], result["n_genes_after"], result["n_genes_removed"])
print(ad.read_h5ad("filtered.h5ad").var_names.tolist())
```

- **Example output:**

```text
4 3 1
['MT-A', 'GENE1', 'GENE2']
```

## 5) `dnakit.transcriptomics.select_highly_variable_genes()`

- **Purpose:** Select highly variable genes.
- **Method:** When `log_normalize=True`, first run `scanpy.pp.normalize_total(target_sum)` and `scanpy.pp.log1p` on a copy, so the input file is unchanged. Then call `scanpy.pp.highly_variable_genes` with `flavor` set to `"seurat"` or `"cell_ranger"`. Results are sorted by normalized dispersion in descending order, and genes with equal dispersion are sorted by name. Genes with tied dispersion are flagged together, so the count can be larger than `n_top_genes`. `seurat_v3` needs raw counts and scikit-misc, and it is not wrapped here.
- **API:** `dnakit.transcriptomics.select_highly_variable_genes(source, *, n_top_genes=2000, flavor="seurat", log_normalize=True, target_sum=10000.0, max_file_bytes=1073741824, max_matrix_elements=50000000)`.
- **Input:** The path to an `.h5ad` file. `n_top_genes` cannot exceed the number of genes. `flavor` is `"seurat"` or `"cell_ranger"`. With `log_normalize=False`, the input must already be log-normalized.
- **Output:** `flavor`, `n_top_genes`, `n_genes_input`, `n_highly_variable`, `highly_variable_genes` (gene names in dispersion order), `log_normalized`, `target_sum`, and `provenance`.
- **Example code:**

```python
import anndata as ad
import numpy as np
import pandas as pd
from dnakit.transcriptomics import select_highly_variable_genes

rng = np.random.default_rng(0)
counts = rng.poisson(3, size=(30, 20)).astype("float32")
counts[:15, :4] += 8
obs = pd.DataFrame({"group": ["A"] * 15 + ["B"] * 15}, index=[f"c{i}" for i in range(30)])
ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(20)])).write_h5ad("cells.h5ad")

result = select_highly_variable_genes("cells.h5ad", n_top_genes=4)
print(result["n_highly_variable"], result["highly_variable_genes"])
```

- **Example output:**

```text
7 ['g5', 'g19', 'g10', 'g16', 'g18', 'g2', 'g3']
```

## 6) `dnakit.transcriptomics.rank_marker_genes()`

- **Purpose:** For each group, compare its cells with the rest and list the top marker genes.
- **Method:** Call `scanpy.tl.rank_genes_groups` with `use_raw=False` and Benjamini–Hochberg correction of the p-values. `method="wilcoxon"` uses the rank-sum test, and `method="t-test"` uses the t-test. When `log_normalize=True`, normalize and log-transform a copy first. The `logfoldchange` is computed by scanpy as the log2 ratio of the expm1 means of X.
- **API:** `dnakit.transcriptomics.rank_marker_genes(source, *, groupby, method="wilcoxon", n_genes=10, log_normalize=True, target_sum=10000.0, max_file_bytes=1073741824, max_matrix_elements=50000000)`.
- **Input:** The path to an `.h5ad` file. `groupby` is a column of `obs`, with no missing values and at least two groups. `n_genes` is the number of genes listed per group.
- **Output:** `groupby`, `method`, `n_groups`, and `groups` (a mapping from group name to a list of `{gene, score, logfoldchange, pval_adj}`, with group names sorted), plus `provenance`.
- **Example code:**

```python
from dnakit.transcriptomics import rank_marker_genes

result = rank_marker_genes("cells.h5ad", groupby="group", n_genes=2)
for group, rows in result["groups"].items():
    print(group, [row["gene"] for row in rows])
```

- **Example output:**

```text
A ['g0', 'g1']
B ['g13', 'g15']
```

## 7) `dnakit.transcriptomics.aggregate_pseudobulk()`

- **Purpose:** Sum the raw counts for each (group, sample) pair, giving pseudobulk counts.
- **Method:** The sums come from `scanpy.get.aggregate(func="sum")` and are not normalized. The number of cells in each pair is counted with pandas, grouping `obs`. Pairs with fewer than `min_cells_per_sample` cells are dropped and counted in `dropped_samples`. Output is sorted by (group, sample).
- **API:** `dnakit.transcriptomics.aggregate_pseudobulk(source, *, group_key, sample_key, min_cells_per_sample=1, output_path=None, max_file_bytes=1073741824, max_matrix_elements=50000000)`.
- **Input:** The path to an `.h5ad` file. `group_key` and `sample_key` are two different `obs` columns. `output_path` is optional. The written `.h5ad` stores the sums in `layers["sum"]` and the cell counts in `obs["n_cells"]`.
- **Output:** `n_pseudobulk_samples`, `dropped_samples`, `n_genes`, `samples` (for each kept pair: `group`, `sample`, `n_cells`, and `total_counts`), `output_path`, and `provenance`.
- **Example code:**

```python
from dnakit.transcriptomics import aggregate_pseudobulk

result = aggregate_pseudobulk("tiny.h5ad", group_key="group", sample_key="sample")
print(result["n_pseudobulk_samples"], result["dropped_samples"])
for sample in result["samples"]:
    print(sample["group"], sample["sample"], sample["n_cells"], sample["total_counts"])
```

- **Example output:**

```text
4 0
x s1 1 4.0
x s2 1 2.0
y s1 1 3.0
y s2 1 7.0
```

## 8) `dnakit.transcriptomics.summarize_mudata()`

- **Purpose:** Summarize the modalities, the dimensions of each modality, and the cell counts in a multi-modal `.h5mu` file.
- **Method:** Read the whole file with `muon.read_h5mu`, subject to `max_file_bytes`. `n_obs` is the union of cells across modalities, and `n_obs_shared` is the number of cells present in every modality.
- **API:** `dnakit.transcriptomics.summarize_mudata(source, *, max_file_bytes=1073741824)`.
- **Input:** The path to an `.h5mu` file.
- **Output:** `n_obs`, `n_obs_shared`, `modalities` (for each modality: `n_obs` and `n_vars`), `obs_columns`, and `provenance`.
- **Example code:**

```python
import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
from dnakit.transcriptomics import summarize_mudata

rng = np.random.default_rng(0)
rna = ad.AnnData(
    X=rng.poisson(1, size=(5, 3)).astype("float32"),
    obs=pd.DataFrame(index=[f"c{i}" for i in range(5)]),
    var=pd.DataFrame(index=["G1", "G2", "G3"]),
)
atac = ad.AnnData(
    X=rng.poisson(1, size=(4, 6)).astype("float32"),
    obs=pd.DataFrame(index=["c0", "c1", "c2", "c9"]),
    var=pd.DataFrame(index=[f"p{i}" for i in range(6)]),
)
mu.write_h5mu("multi.h5mu", mu.MuData({"rna": rna, "atac": atac}))

result = summarize_mudata("multi.h5mu")
print(result["n_obs"], result["n_obs_shared"])
print(result["modalities"])
```

- **Example output:**

```text
6 3
{'atac': {'n_obs': 4, 'n_vars': 6}, 'rna': {'n_obs': 5, 'n_vars': 3}}
```

## 9) `dnakit.transcriptomics.train_scvi_model()`

- **Purpose:** Set up and train an scVI model on raw integer counts, and return the shape of the latent space and a training summary.
- **Method:** First check that the matrix contains finite, non-negative integers, and raise an error otherwise. Then call `SCVI.setup_anndata`, `SCVI(...)`, and `train(...)` in turn, and take the posterior mean with `get_latent_representation()`. `seed` sets `scvi.settings.seed`, which is a process-wide setting. `accelerator` defaults to `"cpu"`. `"gpu"` and `"auto"` leave device choice to Lightning, and `"gpu"` raises an error when no GPU is available.
- **API:** `dnakit.transcriptomics.train_scvi_model(source, *, layer=None, batch_key=None, n_latent=10, n_hidden=128, max_epochs=1, seed=0, accelerator="cpu", max_file_bytes=1073741824, max_matrix_elements=50000000)`.
- **Input:** The path to an `.h5ad` file. `layer` names the layer that holds the counts. `batch_key` is a batch column of `obs`. `n_latent` is 1 to 256, `n_hidden` is 1 to 4096, and `max_epochs` is 1 to 1000.
- **Output:** `n_cells`, `n_genes`, `latent_shape`, `latent_finite`, `latent_dimension_means`, `epochs_trained`, `final_train_elbo`, `seed`, `accelerator`, and `provenance`.
- **Example code:**

```python
import anndata as ad
import numpy as np
import pandas as pd
from dnakit.transcriptomics import train_scvi_model

rng = np.random.default_rng(0)
counts = rng.poisson(1.5, size=(60, 40)).astype("float32")
obs = pd.DataFrame({"batch": ["b1", "b2"] * 30}, index=[f"c{i}" for i in range(60)])
ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(40)])).write_h5ad(
    "scvi.h5ad"
)

result = train_scvi_model(
    "scvi.h5ad", batch_key="batch", n_latent=3, n_hidden=16, verbose=False
)
print(result["latent_shape"], result["epochs_trained"], result["latent_finite"])
```

- **Example output:**

```text
[60, 3] 1 True
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-scanpy]"   # anndata and scanpy: summary, QC, HVG, marker genes, pseudobulk, filter
python -m pip install "dnakit[bio-muon]"     # muon: summarize_mudata
python -m pip install "dnakit[bio-scvi]"     # scvi-tools and anndata: train_scvi_model, requires PyTorch
```

Project pages: [anndata](https://github.com/scverse/anndata), [scanpy](https://github.com/scverse/scanpy), [muon](https://github.com/scverse/muon), [scvi-tools](https://github.com/scverse/scvi-tools).

## How these functions are verified

The test file is `others/tests/unit/test_transcriptomics_backends.py`, with 18 tests. All of them run in an environment with scanpy, muon, and scvi-tools installed. When a package is missing, the tests that need it are skipped. The tests check the following:

- Gene and cell QC values match hand calculations, including detection rates, medians, and mitochondrial fractions.
- `filter_genes_by_min_cells` keeps the expected genes, writes the expected file, and refuses to overwrite an existing file.
- Pseudobulk sums match hand-added totals, and pairs below the minimum cell count are dropped.
- Marker-gene Wilcoxon z-scores, Benjamini–Hochberg adjusted p-values, and log2 fold changes match hand formulas. Results after normalization match a direct scanpy call.
- Highly variable genes match direct scanpy calls for both `seurat` and `cell_ranger`.
- `summarize_mudata` returns the union and shared cell counts that were counted by hand.
- scVI trains for one epoch on CPU, returns a latent space of the right shape with finite values, and gives identical results for two runs with the same seed. Non-integer counts and invalid options are rejected.
- A missing `scanpy`, `anndata`, `scvi`, or `muon` raises `BackendUnavailableError` with the matching extra name.

This batch needs no network access, so it has no live test.

## Current scope

This batch wraps AnnData summaries and `.h5ad` output, scanpy gene and cell QC, highly variable genes (`seurat`, `cell_ranger`), marker genes (`wilcoxon`, `t-test`), pseudobulk sums, the muon `.h5mu` summary, and scVI setup and training.

Not wrapped: scanpy clustering, dimensionality reduction, and neighbor graphs (for example `pca`, `neighbors`, `umap`, `leiden`); regression correction and cell-cycle scoring; the `seurat_v3` highly variable gene flavor; scVI models other than SCVI (for example totalVI and scANVI) and scVI differential expression; muon multi-modal preprocessing; and readers for 10x, loom, and other formats. These need a separate batch, or they can be called directly from the upstream package.
