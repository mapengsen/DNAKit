# AUGUSTUS：基因模型、序列与屏蔽转换

`dnakit.annotation` 中的这组函数是 AUGUSTUS 3.5.0 程序与辅助脚本（[AUGUSTUS](https://github.com/Gaius-Augustus/Augustus)）的类型化封装。每个函数用显式参数列表运行一个固定程序，不经过 shell。请先安装 AUGUSTUS（`conda install -c bioconda augustus`）并把其 `bin` 目录放入 PATH，或传入 `augustus_bin_dir`；缺少程序时抛出带安装提示的 `BackendUnavailableError`。物种参数集从 `AUGUSTUS_CONFIG_PATH`（或 `config_path`）读取，训练类函数会写入该目录，因此需可写。输出路径必须是新路径，不会覆盖已有文件。每个结果是包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及少量解析计数的字典。

## 1) `dnakit.annotation.augustus_add_names_to_gff3()`

- **作用：** 为 GFF3 的基因/mRNA 行添加 Name 属性（add_name_to_gff3.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`add_name_to_gff3.pl`、`only_gene_structure`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_add_names_to_gff3(gff3_path, output_path, *, only_gene_structure=False, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`gff3_path`、`output_path`；可选关键字参数：`only_gene_structure`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
g = m.augustus_gtf_to_gff(str(w / "genes.gtf"), str(w / "genes.gff3"), gff3=True, print_exons=True, verbose=False)
r = m.augustus_add_names_to_gff3(g["outputs"]["output"], str(w / "named.gff3"), verbose=False)
print(r["preview"][1])
```

- **示例输出：**

```text
chr1	AUGUSTUS	mRNA	101	400	.	+	.	ID=g1.t1;Parent=g1;Name=g1.t1
```

## 2) `dnakit.annotation.augustus_cegma_to_gff()`

- **作用：** 把 CEGMA GFF 转换为 AUGUSTUS 训练用 GTF（cegma2gff.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`cegma2gff.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_cegma_to_gff(cegma_gff_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`cegma_gff_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "cegma.gff").write_text("CHR\tcegma\tFirst\t100\t200\t5.0\t+\t0\tK1.1\nCHR\tcegma\tExon\t100\t200\t5.0\t+\t.\tK1.1\n")
r = m.augustus_cegma_to_gff(str(w / "cegma.gff"), str(w / "cds.gtf"), verbose=False)
print(r["preview"][0])
```

- **示例输出：**

```text
CHR	cegma	CDS	100	200	5.0	+	0	transcript_id "g1.K1.1"
```

## 3) `dnakit.annotation.augustus_check_gff_names_in_fasta()`

- **作用：** 检查 GFF 的序列名是否都出现在基因组 FASTA 头部（findGffNamesInFasta.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`findGffNamesInFasta.pl`、`names_missing`、`output_path`、`wrong_column_count`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_check_gff_names_in_fasta(gff_path, genome_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`gff_path`、`genome_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "x.gff").write_text("chr1\t.\tgene\t1\t2\t.\t+\t.\tID=a\nchrX\t.\tgene\t1\t2\t.\t+\t.\tID=b\n")
r = m.augustus_check_gff_names_in_fasta(str(w / "x.gff"), str(w / "genome.fa"), str(w / "chk.txt"), verbose=False)
print(r["names_missing"], r["wrong_column_count"])
```

- **示例输出：**

```text
True False
```

## 4) `dnakit.annotation.augustus_cluster_gene_dependencies()`

- **作用：** 由基因依赖关系文件构建无簇间依赖的基因簇（clusterAndSplitGenes.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`clusterAndSplitGenes.pl`、`clusters_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_cluster_gene_dependencies(dependency_path, clusters_path, *, threshold=None, max_cluster_size=None, transcript_format=False, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`dependency_path`、`clusters_path`；可选关键字参数：`threshold`、`max_cluster_size`、`transcript_format`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 5) `dnakit.annotation.augustus_compare_masking()`

- **作用：** 比较同一组装两种软屏蔽注释（compare_masking.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`compare_masking.pl`、`output_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_compare_masking(softmasked_path, second_softmasked_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`softmasked_path`、`second_softmasked_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
r = m.augustus_compare_masking(str(w / "genome.fa"), str(w / "genome.fa"), str(w / "cmp.txt"), verbose=False)
print(r["preview"])
```

- **示例输出：**

```text
['Masked in both files: 0', 'Masked in File1, only: 0', 'Masked in File2, only: 0']
```

## 6) `dnakit.annotation.augustus_compute_flanking_region()`

- **作用：** 由 GTF 计算平均基因长度与侧翼区长度（computeFlankingRegion.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`computeFlankingRegion.pl`、`numbers`、`output_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_compute_flanking_region(gtf_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`gtf_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
r = m.augustus_compute_flanking_region(str(w / "genes.gtf"), str(w / "flank.txt"), verbose=False)
print(r["preview"][1])
```

- **示例输出：**

```text
The flanking_DNA value is: 150 (the Minimum of 10 000 and 150)
```

## 7) `dnakit.annotation.augustus_convert_genepred()`

- **作用：** 转换 UCSC genePredExt 到 bigGenePred 就绪布局（gp2othergp.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`genePredExt`、`gp2othergp.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_convert_genepred(genepred_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genepred_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "gp.txt").write_text("g1\tchr1\t+\t100\t900\t0\t1\t2\t100,500\t300,900\n")
r = m.augustus_convert_genepred(str(w / "gp.txt"), str(w / "big.gp"), verbose=False)
print(r["preview"][0])
```

- **示例输出：**

```text
chr1	100	900	g1	0	+	0	1	255,0,0	2	200,400	0,400	g1				coding	?	?	?
```

## 8) `dnakit.annotation.augustus_extract_annotation_fasta()`

- **作用：** 由 AUGUSTUS 输出注释块提取蛋白与编码序列（getAnnoFasta.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`.cdsexons`、`.codingseq`、`.mrna`、`augustus.aa`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_extract_annotation_fasta(prediction_path, output_dir, *, genome_path=None, chop_cds=False, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`prediction_path`、`output_dir`；可选关键字参数：`genome_path`、`chop_cds`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 9) `dnakit.annotation.augustus_extract_cds_sequences()`

- **作用：** 由 GTF/GFF3 与基因组提取 CDS 与蛋白序列（getAnnoFastaFromJoingenes.py、stringtie2fa.py 或 gtf2aa.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`drop_invalid_stops`、`getAnnoFastaFromJoingenes.py`、`gtf2aa.pl`、`stringtie2fa.py`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_extract_cds_sequences(genome_path, annotation_path, output_stem, *, method='getAnnoFastaFromJoingenes', annotation_format='gtf', translation_table=None, drop_invalid_stops=False, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genome_path`、`annotation_path`、`output_stem`；可选关键字参数：`method`、`annotation_format`、`translation_table`、`drop_invalid_stops`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
r = m.augustus_extract_cds_sequences(str(w / "genome.fa"), str(w / "genes.gtf"), str(w / "prot"), method="gtf2aa", verbose=False)
print(r["records"])
```

- **示例输出：**

```text
{'protein': 1}
```

## 10) `dnakit.annotation.augustus_extract_mrna_sequences()`

- **作用：** 由 GFF 基因模型从基因组剪接出 mRNA 序列（gffGetmRNA.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`gffGetmRNA.pl`、`output_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_extract_mrna_sequences(gff_path, genome_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`gff_path`、`genome_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 11) `dnakit.annotation.augustus_filter_genbank_genes()`

- **作用：** 按基因名单筛选 GenBank 中的基因模型（filterGenes*.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`exclude_loci`、`exclude_mrna`、`filterGenes.pl`、`filterGenesIn.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_filter_genbank_genes(genbank_path, names_path, output_path, *, mode='exclude_loci', augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genbank_path`、`names_path`、`output_path`；可选关键字参数：`mode`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 12) `dnakit.annotation.augustus_fix_joingenes_gtf()`

- **作用：** 为 joingenes 的 GTF 补上缺失的 gene 行（fix_joingenes_gtf.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`fix_joingenes_gtf.pl`、`gene`、`joingenes`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_fix_joingenes_gtf(gtf_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`gtf_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
j = __import__("dnakit.annotation.augustus_p04_prediction_functions", fromlist=["x"]).augustus_join_gene_sets([str(w / "genes.gtf")], str(w / "joined.gtf"), verbose=False)
r = m.augustus_fix_joingenes_gtf(j["outputs"]["merged"], str(w / "fixed.gtf"), verbose=False)
print(r["preview"][0].split("\t")[2])
```

- **示例输出：**

```text
gene
```

## 13) `dnakit.annotation.augustus_gbrowse_old_to_gff3()`

- **作用：** 把旧版 GBrowse 格式预测转换为 GFF3（gbrowseold2gff3.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`gbrowseold2gff3.pl`、`hints`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_gbrowse_old_to_gff3(prediction_path, output_path, *, hints=False, source=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`prediction_path`、`output_path`；可选关键字参数：`hints`、`source`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 14) `dnakit.annotation.augustus_genbank_loci()`

- **作用：** 提取 GenBank 训练文件的转录本-基因座列表（get_loci_from_gb.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`get_loci_from_gb.pl`、`output_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_genbank_loci(genbank_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genbank_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
gb = m.augustus_gff_to_genbank(str(w / "genes.gtf"), str(w / "genome.fa"), 50, str(w / "train.gb"), verbose=False)
r = m.augustus_genbank_loci(gb["outputs"]["output"], str(w / "loci.tsv"), verbose=False)
print(r["preview"])
```

- **示例输出：**

```text
['g1.t1\tchr1_51-450']
```

## 15) `dnakit.annotation.augustus_genbank_to_gff()`

- **作用：** 把 GenBank 训练基因转换为基因组坐标 GFF（gbSmallDNA2gff.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`gbSmallDNA2gff.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_genbank_to_gff(genbank_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genbank_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
gb = m.augustus_gff_to_genbank(str(w / "genes.gtf"), str(w / "genome.fa"), 50, str(w / "train.gb"), verbose=False)
r = m.augustus_genbank_to_gff(gb["outputs"]["output"], str(w / "back.gff"), verbose=False)
print(r["lines"], r["preview"][0])
```

- **示例输出：**

```text
2 chr1	someSource	CDS	101	200	.	+	.	gene_id "g1.t1"; transcript_id "g1.t1";
```

## 16) `dnakit.annotation.augustus_gff_to_genbank()`

- **作用：** 把 GFF 基因模型与基因组转换为带侧翼 DNA 的 GenBank 训练文件（gff2gbSmallDNA.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`gff2gbSmallDNA.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_gff_to_genbank(gff_path, genome_path, flanking_dna, output_path, *, exclude_genes_path=None, include_genes_path=None, allow_overlaps=False, connected=False, keep_softmasking=False, hardmask=False, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`gff_path`、`genome_path`、`flanking_dna`、`output_path`；可选关键字参数：`exclude_genes_path`、`include_genes_path`、`allow_overlaps`、`connected`、`keep_softmasking`、`hardmask`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
r = m.augustus_gff_to_genbank(str(w / "genes.gtf"), str(w / "genome.fa"), 50, str(w / "train.gb"), verbose=False)
print(r["loci"], r["preview"][0].split()[1], r["preview"][3].strip())
```

- **示例输出：**

```text
1 chr1_51-450 CDS             join(51..150,251..350)
```

## 17) `dnakit.annotation.augustus_gtf_to_gff()`

- **作用：** 把 AUGUSTUS GTF 基因模型转换为 GFF/GFF3（gtf2gff.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`ID`、`Parent`、`gtf2gff.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_gtf_to_gff(gtf_path, output_path, *, gff3=False, print_exons=False, print_utrs=False, print_introns=False, include_stop_in_cds=False, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`gtf_path`、`output_path`；可选关键字参数：`gff3`、`print_exons`、`print_utrs`、`print_introns`、`include_stop_in_cds`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
r = m.augustus_gtf_to_gff(str(w / "genes.gtf"), str(w / "genes.gff3"), gff3=True, print_exons=True, verbose=False)
print(r["preview"][0])
print(r["preview"][1])
```

- **示例输出：**

```text
chr1	AUGUSTUS	gene	101	400	.	+	.	ID=g1;
chr1	AUGUSTUS	mRNA	101	400	.	+	.	ID=g1.t1;Parent=g1;
```

## 18) `dnakit.annotation.augustus_gth_to_gtf()`

- **作用：** 把 GenomeThreader GFF3 比对转换为 GTF（gth2gtf.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`gth2gtf.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_gth_to_gtf(gth_gff3_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`gth_gff3_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 19) `dnakit.annotation.augustus_mask_n_regions()`

- **作用：** 列出基因组 FASTA 中 N（空位）区域（maskNregions.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`maskNregions.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_mask_n_regions(genome_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genome_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
(w / "n.fa").write_text(">s1\nACGTNNNNACGT\n")
r = m.augustus_mask_n_regions(str(w / "n.fa"), str(w / "n.gff"), verbose=False)
print(r["preview"])
```

- **示例输出：**

```text
['s1\tmaskN\tN-region\t5\t8\t.\t.\t.\tN-region']
```

## 20) `dnakit.annotation.augustus_merge_masking()`

- **作用：** 合并同一组装两种软屏蔽注释（merge_masking.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`merge_masking.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_merge_masking(softmasked_path, second_softmasked_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`softmasked_path`、`second_softmasked_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
r = m.augustus_merge_masking(str(w / "genome.fa"), str(w / "genome.fa"), str(w / "merged.fa"), verbose=False)
print(r["records"])
```

- **示例输出：**

```text
1
```

## 21) `dnakit.annotation.augustus_partition_gene_models()`

- **作用：** 把单外显子基因模型划分为基因簇并写成 GenBank（partition_gtf2gb.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`partition_gtf2gb.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_partition_gene_models(genome_path, gtf_path, output_path, *, min_cluster_size=None, max_flanking_region_size=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genome_path`、`gtf_path`、`output_path`；可选关键字参数：`min_cluster_size`、`max_flanking_region_size`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 22) `dnakit.annotation.augustus_random_split()`

- **作用：** 把 GenBank 基因集随机拆分为训练集与测试集（randomSplit.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`randomSplit.pl`、`test_size`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_random_split(genbank_path, test_size, train_path, test_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genbank_path`、`test_size`、`train_path`、`test_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "genes.gtf").write_text(
    'chr1\tAUGUSTUS\tgene\t101\t400\t.\t+\t.\tg1\n'
    'chr1\tAUGUSTUS\ttranscript\t101\t400\t.\t+\t.\tg1.t1\n'
    'chr1\tAUGUSTUS\tstart_codon\t101\t103\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t101\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    'chr1\tAUGUSTUS\tCDS\t301\t400\t.\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";\n'
)
gb = m.augustus_gff_to_genbank(str(w / "genes.gtf"), str(w / "genome.fa"), 50, str(w / "train.gb"), verbose=False)
r = m.augustus_random_split(gb["outputs"]["output"], 1, str(w / "train.split.gb"), str(w / "test.split.gb"), verbose=False)
print(r["train_loci"] + r["test_loci"])
```

- **示例输出：**

```text
1
```

## 23) `dnakit.annotation.augustus_scipio_gff_to_gff()`

- **作用：** 把 Scipio GFF 转换为 AUGUSTUS GFF（scipiogff2gff.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`augustus_gff_to_genbank`、`scipiogff2gff.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_scipio_gff_to_gff(scipio_gff_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`scipio_gff_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 24) `dnakit.annotation.augustus_scipio_yaml_to_gff()`

- **作用：** 把 Scipio YAML 比对转换为 GFF（yaml2gff.1.4.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`filter_status`、`incomplete`、`yaml2gff.1.4.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_scipio_yaml_to_gff(scipio_yaml_path, output_path, *, filter_status=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`scipio_yaml_path`、`output_path`；可选关键字参数：`filter_status`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 25) `dnakit.annotation.augustus_summarize_base_content()`

- **作用：** 汇总每条序列的 A/C/G/T/N 碱基组成（summarizeACGTcontent.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`summarizeACGTcontent.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_summarize_base_content(fasta_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`fasta_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_conversion_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
r = m.augustus_summarize_base_content(str(w / "genome.fa"), str(w / "acgt.txt"), verbose=False)
print(r["preview"][-2])
```

- **示例输出：**

```text
total 4000bp
```

