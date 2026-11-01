# Biopython 树文件、PHYLIP 距离、基序文件与基因数据文件（GEO、CEL、UniGene、Cluster）

本文档涵盖两个封装 [Biopython](https://biopython.org/) 的模块。`dnakit.phylogeny.biopython_p18_tree_functions` 把树写成 Newick、NEXUS、phyloXML、NeXML 和 CDAO 文件，并格式化 PHYLIP 距离矩阵。`dnakit.regulation.biopython_p18_motif_data_functions` 读取基序程序输出（AlignACE、MEME XML、Cluster Buster、PFM 表、XMS），计算基序 GC 含量，读取 GEO SOFT、Affymetrix CEL、UniGene 和 Cluster/TreeView 表达文件，并把聚类结果保存为 CDT、GTR、ATR 和 KGG 文件。

读取函数只读本地文件，写入函数从不覆盖。Biopython 是可选依赖（`pip install "dnakit[bio-formats]"`），只有调用函数时才导入。Biopython 的位置特异性打分矩阵没有 GC 含量（会抛错），所以 `motif_gc_content()` 接受位点或计数矩阵。示例使用模块的完整路径导入；树函数来自 `dnakit.phylogeny.biopython_p18_tree_functions`。

## 示例数据

声明了两个输入：一个 Affymetrix CEL 第 3 版文件（5 x 5 个单元）和一个含三个簇的 UniGene `.data` 文件。

**示例文件：** `demo.CEL`

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

**示例文件：** `demo.unigene`

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

- **作用：** 把一棵或多棵以 Newick 字符串给出的树写成新的树文件。
- **处理方法：** 用 `Bio.Phylo.read()` 解析每个字符串，再用 `Bio.Phylo.write()` 把所有树写成 `newick`、`nexus`、`phyloxml`、`nexml` 或 `cdao`（RDF）。数值型内部节点标签会成为分支置信度，与 Biopython 的 Newick 读取器一致。
- **API：** `dnakit.phylogeny.biopython_p18_tree_functions.write_tree_file(newick_trees, output_path, *, file_format='nexus')`。
- **输入：** Newick 字符串列表、新的输出路径和 `file_format`。
- **输出：** `tree_count`、`bytes_written`、`file_format` 和 `output_path`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny.biopython_p18_tree_functions import write_tree_file

out = Path(tempfile.mkdtemp()) / "trees.nex"
result = write_tree_file(["((A:0.1,B:0.2)90:0.3,(C:0.4,D:0.5):0.1);"], out, file_format="nexus")
print(result["tree_count"], result["file_format"])
print(out.read_text(encoding="utf-8"))
```

- **示例输出：**

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

- **作用：** 把对称距离矩阵格式化为 PHYLIP 文本。
- **处理方法：** 检查矩阵为方阵、对称、非负且对角线为零，然后用 `DistanceMatrix.format_phylip()` 输出（四位小数）。给出 `output_path` 时文本同时写入该新文件。
- **API：** `dnakit.phylogeny.biopython_p18_tree_functions.format_phylip_distance_matrix(names, distances, *, output_path=None)`。
- **输入：** 类群 `names`、完整方阵 `distances` 和可选的新 `output_path`。
- **输出：** `taxon_count`、`text` 和 `output_path`。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p18_tree_functions import format_phylip_distance_matrix

result = format_phylip_distance_matrix(
    ["A", "B", "C"], [[0, 0.2, 0.3], [0.2, 0, 0.4], [0.3, 0.4, 0]]
)
print(result["text"].lstrip())
```

- **示例输出：**

```text
3
A           0.0000  0.2000  0.3000
B           0.2000  0.0000  0.4000
C           0.3000  0.4000  0.0000
```

## 3) `dnakit.regulation.biopython_p18_motif_data_functions.read_motif_collection()`

- **作用：** 读取 AlignACE、MEME、Cluster Buster、PFM 或 XMS 文件中的全部基序。
- **处理方法：** 对 `alignace`、`meme`（MEME XML 输出）、`clusterbuster`、`pfm-four-columns`、`pfm-four-rows` 或 `xms` 调用 `Bio.motifs.parse()`。启用 `include_instances` 时列出结合位点（MEME 还含序列名、起点、链方向和 p 值），每个基序最多 1000 个。
- **API：** `dnakit.regulation.biopython_p18_motif_data_functions.read_motif_collection(path, *, file_format='meme', include_instances=False, max_motifs=1000)`。
- **输入：** 基序文件路径、`file_format`、`include_instances` 与 `max_motifs`。
- **输出：** `motif_count` 与 `motifs`，含名称、长度、共有序列、碱基 `counts`、`gc_content`，以及存在时的 `id`、`evalue`、`num_occurrences`、`score`、`mask`。
- **示例代码：**

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

- **示例输出：**

```text
2
M1 4 ACGT 0.5
M2 2 GG 1.0
```

## 4) `dnakit.regulation.biopython_p18_motif_data_functions.motif_gc_content()`

- **作用：** 由位点或碱基计数矩阵计算基序的 GC 比例。
- **处理方法：** 构建计数矩阵并读取 `GenericPositionMatrix.gc_content`：所有位置上 C 与 G 的总计数除以全部碱基总计数。不接受位置特异性打分矩阵，因为 Biopython 无法计算其 GC 含量。
- **API：** `dnakit.regulation.biopython_p18_motif_data_functions.motif_gc_content(motifs)`。
- **输入：** 基序名到等长 DNA 位点列表或计数矩阵 `{A: [...], C: [...], G: [...], T: [...]}` 的映射。
- **输出：** `motif_count` 与 `motifs`，含 `name`、`input`、`length`、`gc_content`。
- **示例代码：**

```python
from dnakit.regulation.biopython_p18_motif_data_functions import motif_gc_content

result = motif_gc_content(
    {"sites": ["ACGT", "ACGA"], "counts": {"A": [1, 0], "C": [0, 3], "G": [0, 0], "T": [2, 0]}}
)
for motif in result["motifs"]:
    print(motif["name"], motif["input"], motif["length"], motif["gc_content"])
```

- **示例输出：**

```text
sites sites 4 0.5
counts counts 2 0.5
```

## 5) `dnakit.regulation.biopython_p18_motif_data_functions.read_geo_soft()`

- **作用：** 读取 GEO SOFT 文件中的实体（平台、样本、系列、数据集）。
- **处理方法：** 调用 `Bio.Geo.parse()`。每条记录含实体类型、标识符、属性（不含开头的 `!`）、列定义和表格：表头行、最多 `max_table_rows` 行数据以及总行数。
- **API：** `dnakit.regulation.biopython_p18_motif_data_functions.read_geo_soft(path, *, max_records=1000, max_table_rows=1000)`。
- **输入：** SOFT 文件路径、`max_records` 与 `max_table_rows`。
- **输出：** `record_count` 与 `records`。
- **示例代码：**

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

- **示例输出：**

```text
SAMPLE GSM1 {'Sample_title': 'first sample', 'sample_table_begin': '', 'sample_table_end': ''}
['ID_REF', 'VALUE'] 2 [['p1', '1.5'], ['p2', '2.5']]
```

## 6) `dnakit.regulation.biopython_p18_motif_data_functions.read_affymetrix_cel()`

- **作用：** 读取 Affymetrix CEL 探针强度文件（第 3 版文本或第 4 版二进制）。
- **处理方法：** 自动识别版本并调用 `Bio.Affy.CelFile.read()`。返回阵列尺寸、算法及参数、网格角点、DAT 头、掩蔽与离群计数以及强度汇总统计。启用 `include_values` 时，强度、标准差和像素数矩阵以行列表返回（最多 `max_cells` 个单元）。
- **API：** `dnakit.regulation.biopython_p18_motif_data_functions.read_affymetrix_cel(path, *, include_values=False, max_cells=1000000)`。
- **输入：** CEL 文件路径、`include_values` 与 `max_cells`。
- **输出：** `version`、`rows`、`columns`、`cell_count`、`algorithm`、强度最小值/最大值/均值，以及可选的 `intensities`、`standard_deviations`、`pixel_counts`。
- **示例代码：**

```python
from dnakit.regulation.biopython_p18_motif_data_functions import read_affymetrix_cel

result = read_affymetrix_cel("demo.CEL", include_values=True)
print(result["version"], result["rows"], result["columns"], result["cell_count"], result["algorithm"])
print(result["intensity_minimum"], result["intensity_maximum"], round(result["intensity_mean"], 2))
print(result["intensities"][0])
```

- **示例输出：**

```text
3 5 5 25 Percentile
164.0 22177.0 8677.28
[234.0, 170.0, 22177.0, 164.0, 22104.0]
```

## 7) `dnakit.regulation.biopython_p18_motif_data_functions.read_unigene()`

- **作用：** 读取 UniGene `.data` 平面文件中的簇。
- **处理方法：** 调用 `Bio.UniGene.parse()`。每条记录含簇标识符（如 `Hs.2`）、物种、标题、基因符号、细胞带、染色体、表达来源、基因与 LocusLink 标识符，以及解析后的蛋白相似性、序列、STS 和转录图谱行。
- **API：** `dnakit.regulation.biopython_p18_motif_data_functions.read_unigene(path, *, max_records=1000)`。
- **输入：** `.data` 文件路径与 `max_records`。
- **输出：** `record_count` 与 `records`。
- **示例代码：**

```python
from dnakit.regulation.biopython_p18_motif_data_functions import read_unigene

result = read_unigene("demo.unigene")
print(result["record_count"])
for record in result["records"]:
    print(record["ID"], len(record["sequence"]), record["sequence"][0]["acc"], record["title"][:30])
```

- **示例输出：**

```text
3
Mdm.1 8 XM_001377126.1 Similar to cytochrome P450 CYP
Mdm.11 5 EC091552.1 Transcribed locus, weakly simi
Mdm.14 1 DR038236.1 Transcribed locus, strongly si
```

## 8) `dnakit.regulation.biopython_p18_motif_data_functions.read_cluster_expression_data()`

- **作用：** 读取 Cluster/TreeView 制表符分隔的表达表（行为基因，列为样本）。
- **处理方法：** 调用 `Bio.Cluster.read()`。第一列是基因标识符，可接 `NAME` 列、`GWEIGHT` 列、样本列，以及可选的 `EWEIGHT`、`GORDER`、`EORDER` 行。启用 `include_data` 时同时列出数值矩阵和掩码（1 表示有值）。
- **API：** `dnakit.regulation.biopython_p18_motif_data_functions.read_cluster_expression_data(path, *, include_data=False)`。
- **输入：** 表格路径与 `include_data`。
- **输出：** `gene_count`、`experiment_count`、标识符、名称、权重、顺序、`has_missing_values`，以及可选的 `data` 和 `mask`。
- **示例代码：**

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

- **示例输出：**

```text
12 3 ['s1', 's2', 's3']
['g0', 'g1', 'g2'] ['name0', 'name1', 'name2'] [1.0, 2.0, 3.0]
```

## 9) `dnakit.regulation.biopython_p18_motif_data_functions.cluster_expression_data()`

- **作用：** 对表达表聚类，并保存为 Java TreeView 可读的结果（CDT、GTR、ATR、KGG）。
- **处理方法：** `algorithm="hierarchical"` 用 `linkage`（single、complete（Biopython 默认）、centroid 或 average）对基因聚类，启用 `cluster_samples` 时也对样本聚类，写出 `<prefix>.cdt`、`.gtr` 和 `.atr`。`algorithm="kmeans"` 运行 k-means（`center="mean"`）或 k-medians，聚成 `cluster_count` 类、随机重启 `passes` 次，写出 `<prefix>_K_G<k>.cdt` 和 `.kgg`；给出 `initial_clusters` 可得到确定性结果。`distance` 可取 euclidean、city_block、pearson、absolute_pearson、uncentered、absolute_uncentered、spearman 或 kendall。若已有以该前缀开头的文件则报错。
- **API：** `dnakit.regulation.biopython_p18_motif_data_functions.cluster_expression_data(path, output_prefix, *, algorithm='hierarchical', linkage='complete', distance='euclidean', cluster_samples=False, cluster_count=2, passes=1, center='mean', initial_clusters=None)`。
- **输入：** 表格路径、输出前缀以及上述选项。
- **输出：** `output_files`、`gene_count`，层次聚类还有树的规模，k-means 还有 `cluster_labels`、`within_cluster_error` 和 `solutions_found`。
- **示例代码：**

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

- **示例输出：**

```text
['hier.atr', 'hier.cdt', 'hier.gtr'] 11 2
['km_K_G3.cdt', 'km_K_G3.kgg']
[0, 1, 2, 0, 2, 0, 2, 1, 0, 2, 0, 2] 8.833
```
