# Biopython tree files, PHYLIP distances, motif files and gene-data files (GEO, CEL, UniGene, Cluster)

This document covers two modules that wrap [Biopython](https://biopython.org/). `dnakit.phylogeny.biopython_p18_tree_functions` writes trees to Newick, NEXUS, phyloXML, NeXML and CDAO files and formats PHYLIP distance matrices. `dnakit.regulation.biopython_p18_motif_data_functions` reads motif-finder output (AlignACE, MEME XML, Cluster Buster, PFM tables, XMS), computes motif GC content, reads GEO SOFT, Affymetrix CEL, UniGene and Cluster/TreeView expression files, and saves clustering results as CDT, GTR, ATR and KGG files.

Readers read local files. Writers never overwrite. Biopython is an optional dependency (`pip install "dnakit[bio-formats]"`) imported only when a function is called. Position-specific scoring matrices have no GC content in Biopython (it raises), so `motif_gc_content()` takes sites or count matrices. The examples import from the two modules by their full path; the tree functions come from `dnakit.phylogeny.biopython_p18_tree_functions`.

## Example data

Two declared inputs: an Affymetrix CEL version 3 file (5 x 5 cells) and a UniGene `.data` file with three clusters.

**Example file:** `demo.CEL`

```text
[CEL]
Version=3

[HEADER]
Cols=5
Rows=5
TotalX=5
TotalY=5
OffsetX=0
OffsetY=0
GridCornerUL=206 129
GridCornerUR=3570 107
GridCornerLR=3597 3470
GridCornerLL=234 3492
Axis-invertX=0
AxisInvertY=0
swapXY=0
DatHeader=[11..65533]  1g_A9AF:CLS=3684 RWS=3684 XIN=1  YIN=1  VE=30        2.0 08/23/07 11:23:24 50205880  M10      Tgondii_SNP1.1sq          570  25356.509766  3.500000  1.5600  6
Algorithm=Percentile
AlgorithmParameters=Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:TRUE;FullFeatureWidth:7;FullFeatureHeight:7;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:TRUE;PoolWidthExtenstion:2;PoolHeightExtension:2;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000

[INTENSITY]
NumberCells=25
CellHeader=X	Y	MEAN	STDV	NPIXELS
  0	  0	234.0	24.0	 25
  1	  0	22428.0	2378.1	 25
  2	  0	249.0	27.7	 25
  3	  0	22339.0	2264.8	 25
  4	  0	156.0	27.2	 25
  1	  0	170.0	34.5	 25
  2	  0	22177.0	2669.0	 25
  3	  0	164.0	19.7	 25
  4	  0	22104.0	3661.2	 25
  0   1 188.0 29.8   25
  1	  1	188.0	29.8	 25
  2	  1	21871.0	2795.9	 25
  3	  1	168.0	67.9	 25
  4	  1	21883.0	2792.4	 25
  0   2 188.0 29.8   25
  1	  2	193.0	88.7	 25
  2	  2	21455.0	2976.5	 25
  3	  2	198.0	62.0	 25
  4	  2	21300.0	2914.5	 25
  0   3 188.0 29.8   25
  1	  3	182.0	76.2	 25
  2	  3	21438.0	2759.5	 25
  3	  3	188.0	49.2	 25
  4	  3	20945.0	2762.0	 25
  0	  4	193.0	38.8	 25
  1	  4	20370.0	2611.8	 25
  2	  4	174.0	26.6	 25
  3	  4	20605.0	2810.7	 25
  4	  4	168.0	24.1	 25

[MASKS]
NumberCells=3
CellHeader=X  Y
4 1
4 2
3 1

[OUTLIERS]
NumberCells=3
CellHeader=X  Y
1 1
2 1
1 3

[MODIFIED]
NumberCells=3
CellHeader=X  Y ORIGMEAN
3 1 189.0
3 2 21775.0
4 1 220.0
```

**Example file:** `demo.unigene`

```text
ID          Mdm.1
TITLE       Similar to cytochrome P450 CYP2C75
EXPRESS     liver
CHROMOSOME  *
SCOUNT      8
SEQUENCE    ACC=XM_001377126.1; NID=g126272339; PID=g126272340; SEQTYPE=Model
SEQUENCE    ACC=XM_001377113.1; NID=g126272337; PID=g126272338; SEQTYPE=Model
SEQUENCE    ACC=XM_001377092.1; NID=g126272333; PID=g126272334; SEQTYPE=Model
SEQUENCE    ACC=XM_001377142.1; NID=g126272341; PID=g126272342; SEQTYPE=Model
SEQUENCE    ACC=XM_001377105.1; NID=g126272335; PID=g126272336; SEQTYPE=Model
SEQUENCE    ACC=DR038220.1; NID=g66840115; CLONE=U41; LID=17438; SEQTYPE=EST
SEQUENCE    ACC=DR038229.1; NID=g66840124; CLONE=P9; LID=17438; SEQTYPE=EST
SEQUENCE    ACC=DR038222.1; NID=g66840117; CLONE=U84; LID=17438; SEQTYPE=EST
//
ID          Mdm.11
TITLE       Transcribed locus, weakly similar to XP_576039.2  similar to Glutathione S-transferase alpha-4 (Glutathione S-transferase Yk) (GST Yk) (GST 8-8) (GST K) (GST A4-4) [Rattus norvegicus]
EXPRESS     small intestine; liver
CHROMOSOME  *
SCOUNT      5
SEQUENCE    ACC=EC091552.1; NID=g104746629; CLONE=AD76; LID=19764; SEQTYPE=EST
SEQUENCE    ACC=EC091553.1; NID=g104746633; CLONE=AG12; LID=19764; SEQTYPE=EST
SEQUENCE    ACC=EC091554.1; NID=g104746636; CLONE=AG87; LID=19764; SEQTYPE=EST
SEQUENCE    ACC=EC091555.1; NID=g104746639; CLONE=AD12; LID=19764; SEQTYPE=EST
SEQUENCE    ACC=DR038233.1; NID=g66840128; CLONE=U61; LID=17438; SEQTYPE=EST
//
ID          Mdm.14
TITLE       Transcribed locus, strongly similar to XP_001364858.1  similar to albumin [Monodelphis domestica]
EXPRESS     liver
SCOUNT      1
SEQUENCE    ACC=DR038236.1; NID=g66840131; CLONE=P1; LID=17438; SEQTYPE=EST
//
```

## 1) `dnakit.phylogeny.biopython_p18_tree_functions.write_tree_file()`

- **Purpose:** Write one or more trees, given as Newick strings, to a new tree file.
- **Method:** Parse each string with `Bio.Phylo.read()` and write all trees with `Bio.Phylo.write()` as `newick`, `nexus`, `phyloxml`, `nexml` or `cdao` (RDF). Numeric internal labels become branch confidences, as Biopython's Newick reader does.
- **API:** `dnakit.phylogeny.biopython_p18_tree_functions.write_tree_file(newick_trees, output_path, *, file_format='nexus')`.
- **Input:** A list of Newick strings, a new output path and `file_format`.
- **Output:** `tree_count`, `bytes_written`, `file_format` and `output_path`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny.biopython_p18_tree_functions import write_tree_file

out = Path(tempfile.mkdtemp()) / "trees.nex"
result = write_tree_file(["((A:0.1,B:0.2)90:0.3,(C:0.4,D:0.5):0.1);"], out, file_format="nexus")
print(result["tree_count"], result["file_format"])
print(out.read_text(encoding="utf-8"))
```

- **Example output:**

```text
1 nexus
#NEXUS
Begin Taxa;
 Dimensions NTax=4;
 TaxLabels A B C D;
End;
Begin Trees;
 Tree tree1=((A:0.1,B:0.2)90.00:0.3,(C:0.4,D:0.5):0.1):0;
End;
```

## 2) `dnakit.phylogeny.biopython_p18_tree_functions.format_phylip_distance_matrix()`

- **Purpose:** Format a symmetric distance matrix as PHYLIP text.
- **Method:** Check that the matrix is square, symmetric, non-negative with a zero diagonal, then print it with `DistanceMatrix.format_phylip()` (four decimals). When `output_path` is given the text is also written to that new file.
- **API:** `dnakit.phylogeny.biopython_p18_tree_functions.format_phylip_distance_matrix(names, distances, *, output_path=None)`.
- **Input:** Taxon `names`, the full square `distances` matrix and an optional new `output_path`.
- **Output:** `taxon_count`, `text` and `output_path`.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p18_tree_functions import format_phylip_distance_matrix

result = format_phylip_distance_matrix(
    ["A", "B", "C"], [[0, 0.2, 0.3], [0.2, 0, 0.4], [0.3, 0.4, 0]]
)
print(result["text"].lstrip())
```

- **Example output:**

```text
3
A           0.0000  0.2000  0.3000
B           0.2000  0.0000  0.4000
C           0.3000  0.4000  0.0000
```

## 3) `dnakit.regulation.biopython_p18_motif_data_functions.read_motif_collection()`

- **Purpose:** Read every motif of an AlignACE, MEME, Cluster Buster, PFM or XMS file.
- **Method:** Call `Bio.motifs.parse()` for `alignace`, `meme` (MEME XML output), `clusterbuster`, `pfm-four-columns`, `pfm-four-rows` or `xms`. With `include_instances` the binding sites (and for MEME their sequence name, start, strand and p-value) are listed, at most 1000 per motif.
- **API:** `dnakit.regulation.biopython_p18_motif_data_functions.read_motif_collection(path, *, file_format='meme', include_instances=False, max_motifs=1000)`.
- **Input:** A motif file path, `file_format`, `include_instances` and `max_motifs`.
- **Output:** `motif_count` and `motifs` with name, length, consensus, base `counts`, `gc_content`, and where present `id`, `evalue`, `num_occurrences`, `score` and `mask`.
- **Example code:**

```python
from pathlib import Path

from dnakit.regulation.biopython_p18_motif_data_functions import read_motif_collection

Path("motifs.cb").write_text(
    ">M1\n3\t0\t0\t0\n0\t3\t0\t0\n0\t0\t3\t0\n1\t0\t0\t2\n>M2\n0\t0\t2\t0\n0\t0\t2\t0\n",
    encoding="utf-8",
)
result = read_motif_collection("motifs.cb", file_format="clusterbuster")
print(result["motif_count"])
for motif in result["motifs"]:
    print(motif["name"], motif["length"], motif["consensus"], motif["gc_content"])
```

- **Example output:**

```text
2
M1 4 ACGT 0.5
M2 2 GG 1.0
```

## 4) `dnakit.regulation.biopython_p18_motif_data_functions.motif_gc_content()`

- **Purpose:** Compute the GC fraction of motifs from their sites or base-count matrices.
- **Method:** Build the count matrix and read `GenericPositionMatrix.gc_content`: the total count of C and G over the total count of all bases across all positions. Position-specific scoring matrices are not accepted because Biopython cannot compute their GC content.
- **API:** `dnakit.regulation.biopython_p18_motif_data_functions.motif_gc_content(motifs)`.
- **Input:** A mapping of motif name to a list of equal-length DNA sites or to a count matrix `{A: [...], C: [...], G: [...], T: [...]}`.
- **Output:** `motif_count` and `motifs` with `name`, `input`, `length`, `gc_content`.
- **Example code:**

```python
from dnakit.regulation.biopython_p18_motif_data_functions import motif_gc_content

result = motif_gc_content(
    {"sites": ["ACGT", "ACGA"], "counts": {"A": [1, 0], "C": [0, 3], "G": [0, 0], "T": [2, 0]}}
)
for motif in result["motifs"]:
    print(motif["name"], motif["input"], motif["length"], motif["gc_content"])
```

- **Example output:**

```text
sites sites 4 0.5
counts counts 2 0.5
```

## 5) `dnakit.regulation.biopython_p18_motif_data_functions.read_geo_soft()`

- **Purpose:** Read the entities (platform, sample, series, dataset) of a GEO SOFT file.
- **Method:** Call `Bio.Geo.parse()`. Each record has its entity type, identifier, attributes (without the leading `!` marker), column definitions and table: the header row, up to `max_table_rows` data rows and the total row count.
- **API:** `dnakit.regulation.biopython_p18_motif_data_functions.read_geo_soft(path, *, max_records=1000, max_table_rows=1000)`.
- **Input:** A SOFT file path, `max_records` and `max_table_rows`.
- **Output:** `record_count` and `records`.
- **Example code:**

```python
from pathlib import Path

from dnakit.regulation.biopython_p18_motif_data_functions import read_geo_soft

Path("demo.soft").write_text(
    "^SAMPLE = GSM1\n!Sample_title = first sample\n#ID_REF = probe id\n#VALUE = signal\n"
    "!sample_table_begin\nID_REF\tVALUE\np1\t1.5\np2\t2.5\n!sample_table_end\n",
    encoding="utf-8",
)
result = read_geo_soft("demo.soft")
record = result["records"][0]
print(record["entity_type"], record["entity_id"], record["entity_attributes"])
print(record["table_header"], record["table_row_count"], record["table_rows"])
```

- **Example output:**

```text
SAMPLE GSM1 {'Sample_title': 'first sample', 'sample_table_begin': '', 'sample_table_end': ''}
['ID_REF', 'VALUE'] 2 [['p1', '1.5'], ['p2', '2.5']]
```

## 6) `dnakit.regulation.biopython_p18_motif_data_functions.read_affymetrix_cel()`

- **Purpose:** Read an Affymetrix CEL file (version 3 text or version 4 binary) of probe intensities.
- **Method:** Detect the version and call `Bio.Affy.CelFile.read()`. Returned are the array size, algorithm and parameters, grid corners, DAT header, masked and outlier counts and summary statistics of the intensities. With `include_values` the intensity, standard-deviation and pixel-count matrices are returned as lists of rows (at most `max_cells` cells).
- **API:** `dnakit.regulation.biopython_p18_motif_data_functions.read_affymetrix_cel(path, *, include_values=False, max_cells=1000000)`.
- **Input:** A CEL file path, `include_values` and `max_cells`.
- **Output:** `version`, `rows`, `columns`, `cell_count`, `algorithm`, intensity minimum, maximum and mean, and optionally `intensities`, `standard_deviations`, `pixel_counts`.
- **Example code:**

```python
from dnakit.regulation.biopython_p18_motif_data_functions import read_affymetrix_cel

result = read_affymetrix_cel("demo.CEL", include_values=True)
print(result["version"], result["rows"], result["columns"], result["cell_count"], result["algorithm"])
print(result["intensity_minimum"], result["intensity_maximum"], round(result["intensity_mean"], 2))
print(result["intensities"][0])
```

- **Example output:**

```text
3 5 5 25 Percentile
164.0 22177.0 8677.28
[234.0, 170.0, 22177.0, 164.0, 22104.0]
```

## 7) `dnakit.regulation.biopython_p18_motif_data_functions.read_unigene()`

- **Purpose:** Read the clusters of a UniGene `.data` flat file.
- **Method:** Call `Bio.UniGene.parse()`. Each record has the cluster identifier (for example `Hs.2`), species, title, symbol, cytoband, chromosome, expression sources, gene and LocusLink identifiers, and the parsed protein-similarity, sequence, STS and transcript-map lines.
- **API:** `dnakit.regulation.biopython_p18_motif_data_functions.read_unigene(path, *, max_records=1000)`.
- **Input:** A `.data` file path and `max_records`.
- **Output:** `record_count` and `records`.
- **Example code:**

```python
from dnakit.regulation.biopython_p18_motif_data_functions import read_unigene

result = read_unigene("demo.unigene")
print(result["record_count"])
for record in result["records"]:
    print(record["ID"], len(record["sequence"]), record["sequence"][0]["acc"], record["title"][:30])
```

- **Example output:**

```text
3
Mdm.1 8 XM_001377126.1 Similar to cytochrome P450 CYP
Mdm.11 5 EC091552.1 Transcribed locus, weakly simi
Mdm.14 1 DR038236.1 Transcribed locus, strongly si
```

## 8) `dnakit.regulation.biopython_p18_motif_data_functions.read_cluster_expression_data()`

- **Purpose:** Read a Cluster/TreeView tab-delimited expression table (genes in rows, samples in columns).
- **Method:** Call `Bio.Cluster.read()`. The first column holds gene identifiers, optionally followed by a `NAME` column, a `GWEIGHT` column, the sample columns and optional `EWEIGHT`, `GORDER` and `EORDER` rows. With `include_data` the numeric matrix and the mask (1 = value present) are listed too.
- **API:** `dnakit.regulation.biopython_p18_motif_data_functions.read_cluster_expression_data(path, *, include_data=False)`.
- **Input:** A table path and `include_data`.
- **Output:** `gene_count`, `experiment_count`, identifiers, names, weights, orders, `has_missing_values` and optionally `data` and `mask`.
- **Example code:**

```python
from pathlib import Path

rows = ["YORF\tNAME\tGWEIGHT\ts1\ts2\ts3"]
for i in range(12):
    rows.append(f"g{i}\tname{i}\t1\t{i % 3}\t{(i * 2) % 5}\t{(i * 3) % 7}")
Path("expression.txt").write_text("\n".join(rows) + "\n", encoding="utf-8")

from dnakit.regulation.biopython_p18_motif_data_functions import read_cluster_expression_data

result = read_cluster_expression_data("expression.txt", include_data=True)
print(result["gene_count"], result["experiment_count"], result["experiment_ids"])
print(result["gene_ids"][:3], result["gene_names"][:3], result["data"][1])
```

- **Example output:**

```text
12 3 ['s1', 's2', 's3']
['g0', 'g1', 'g2'] ['name0', 'name1', 'name2'] [1.0, 2.0, 3.0]
```

## 9) `dnakit.regulation.biopython_p18_motif_data_functions.cluster_expression_data()`

- **Purpose:** Cluster an expression table and save the result for Java TreeView (CDT, GTR, ATR, KGG).
- **Method:** `algorithm="hierarchical"` clusters the genes with `linkage` single, complete (the Biopython default), centroid or average, and the samples too with `cluster_samples`; it writes `<prefix>.cdt`, `.gtr` and `.atr`. `algorithm="kmeans"` runs k-means (`center="mean"`) or k-medians with `cluster_count` clusters and `passes` random restarts and writes `<prefix>_K_G<k>.cdt` and `.kgg`; give `initial_clusters` for a deterministic run. `distance` is euclidean, city_block, pearson, absolute_pearson, uncentered, absolute_uncentered, spearman or kendall. The function stops if any file starting with the prefix exists.
- **API:** `dnakit.regulation.biopython_p18_motif_data_functions.cluster_expression_data(path, output_prefix, *, algorithm='hierarchical', linkage='complete', distance='euclidean', cluster_samples=False, cluster_count=2, passes=1, center='mean', initial_clusters=None)`.
- **Input:** A table path, an output prefix and the options above.
- **Output:** `output_files`, `gene_count` and, for hierarchical clustering, the tree sizes or, for k-means, `cluster_labels`, `within_cluster_error` and `solutions_found`.
- **Example code:**

```python
from pathlib import Path

rows = ["YORF\tNAME\tGWEIGHT\ts1\ts2\ts3"]
for i in range(12):
    rows.append(f"g{i}\tname{i}\t1\t{i % 3}\t{(i * 2) % 5}\t{(i * 3) % 7}")
Path("expression.txt").write_text("\n".join(rows) + "\n", encoding="utf-8")
import tempfile

from dnakit.regulation.biopython_p18_motif_data_functions import cluster_expression_data

work = Path(tempfile.mkdtemp())
hier = cluster_expression_data("expression.txt", work / "hier", cluster_samples=True)
print(sorted(Path(p).name for p in hier["output_files"]), hier["gene_tree_nodes"], hier["sample_tree_nodes"])
kmeans = cluster_expression_data(
    "expression.txt", work / "km", algorithm="kmeans", cluster_count=3,
    initial_clusters=[i % 3 for i in range(12)],
)
print(sorted(Path(p).name for p in kmeans["output_files"]))
print(kmeans["cluster_labels"], round(kmeans["within_cluster_error"], 3))
```

- **Example output:**

```text
['hier.atr', 'hier.cdt', 'hier.gtr'] 11 2
['km_K_G3.cdt', 'km_K_G3.kgg']
[0, 1, 2, 0, 2, 0, 2, 1, 0, 2, 0, 2] 8.833
```
