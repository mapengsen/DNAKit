# 基因级数据库可选后端

这一组函数是 `dnakit.databases` 中基于可选依赖的基因级适配器，覆盖 [gget](https://github.com/pachterlab/gget)、[pyensembl](https://github.com/openvax/pyensembl)、[genomepy](https://github.com/vanheeringen-lab/genomepy)、[biothings_client](https://github.com/biothings/biothings_client.py)、[mygene](https://github.com/biothings/mygene.py)、[myvariant](https://github.com/biothings/myvariant.py)、[pybiomart](https://github.com/jrderuiter/pybiomart) 和 [ncbi-datasets-pylib](https://github.com/ncbi/datasets)。本地函数只读取显式给出的 GTF、FASTA、BED 或 VCF 文件，计算基因位点、序列切片、突变窗口和 HGVS 标识；网络函数访问 Ensembl、MyGene.info、MyVariant.info、BioMart 或 NCBI 在线服务，不做本地计算。所有依赖都是可选的，只有调用对应函数时才会导入，缺少包时会报错并给出安装提示。

本地示例只使用下面四个小文件，它们是为演示接口而写的人工数据，不代表真实基因组。

## 示例数据

**示例文件：** `genes.gtf`

```text
##description: tiny Ensembl-style annotation for DNAKit examples
chr1	test	gene	11	200	.	+	.	gene_id "ENSG00000000001"; gene_name "GENEA"; gene_biotype "protein_coding";
chr1	test	transcript	11	200	.	+	.	gene_id "ENSG00000000001"; gene_name "GENEA"; transcript_id "ENST00000000011"; transcript_name "GENEA-201"; transcript_biotype "protein_coding";
chr1	test	exon	11	60	.	+	.	gene_id "ENSG00000000001"; gene_name "GENEA"; transcript_id "ENST00000000011"; transcript_name "GENEA-201"; exon_number "1"; exon_id "ENSE00000000001";
chr1	test	exon	121	200	.	+	.	gene_id "ENSG00000000001"; gene_name "GENEA"; transcript_id "ENST00000000011"; transcript_name "GENEA-201"; exon_number "2"; exon_id "ENSE00000000002";
chr1	test	gene	301	380	.	-	.	gene_id "ENSG00000000002"; gene_name "GENEB"; gene_biotype "lncRNA";
chr1	test	transcript	301	380	.	-	.	gene_id "ENSG00000000002"; gene_name "GENEB"; transcript_id "ENST00000000022"; transcript_name "GENEB-201"; transcript_biotype "lncRNA";
chr1	test	exon	301	380	.	-	.	gene_id "ENSG00000000002"; gene_name "GENEB"; transcript_id "ENST00000000022"; transcript_name "GENEB-201"; exon_number "1"; exon_id "ENSE00000000003";
```

**示例文件：** `genome.fa`

```text
>chr1 example contig
ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT
TGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT
TGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT
TGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
```

**示例文件：** `regions.bed`

```text
chr1	10	20	r1	0	+
chr1	300	310	r2	0	-
```

**示例文件：** `sample.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	10	rs1	A	G,T	.	PASS	.
chr1	20	.	C	.	.	.	.
chrUn	5	.	C	T	.	.	.
chr2	30	.	GA	G	.	.	.
```

基因模型是两个人工基因：GENEA 在 `chr1:11-200` 的正链上，外显子为 11-60 和 121-200；GENEB 在 `chr1:301-380` 的负链上。基因组序列是一条 400 bp 的周期性序列，便于人工核对切片。

## 1) `dnakit.databases.search_ensembl_genes()`

- **作用：** 用自由文本词在 Ensembl 中搜索基因。
- **处理方法：** 调用 gget 的 Ensembl 搜索接口，以 JSON 形式返回结果，并限制返回条数。
- **API：** `dnakit.databases.search_ensembl_genes(searchwords, *, species="homo_sapiens", release=None, andor="or", max_results=100)`。
- **输入：** 一个词或词列表；`species` 为 `genus_species` 形式；`andor` 取 `"or"` 或 `"and"`。
- **输出：** 包含 `count`、`results` 和 `provenance` 的字典，`results` 中每条记录含 Ensembl ID、基因名、生物类型和描述。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import search_ensembl_genes

result = search_ensembl_genes(["CDK2"], max_results=5)
print(result["count"], [row["gene_name"] for row in result["results"]])
```

- **示例输出（未本地运行）：** 需要访问 Ensembl 在线服务，结果随 Ensembl 版本变化，因此不附固定输出。返回字典中的 `results` 为列表，每项包含 `ensembl_id` 和 `gene_name`。

## 2) `dnakit.databases.fetch_ensembl_gene_info()`

- **作用：** 按 Ensembl ID 获取基因元数据，包括坐标、链、规范转录本和各转录本列表。
- **处理方法：** 调用 gget 的 `info` 接口。默认不查询 NCBI 和 UniProt，只保留 Ensembl 的信息，以减少请求数量。
- **API：** `dnakit.databases.fetch_ensembl_gene_info(ensembl_ids, *, external_sources=False)`。
- **输入：** Ensembl ID 列表，最多 200 个；`external_sources=True` 时同时查询 NCBI 和 UniProt。
- **输出：** `records` 字典，键为请求的 ID，值为 gget 返回的元数据或 `None`（未找到）；另有 `requested` 和 `found` 计数。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import fetch_ensembl_gene_info

result = fetch_ensembl_gene_info(["ENSG00000123374"])
print(result["records"]["ENSG00000123374"]["ensembl_gene_name"])
```

- **示例输出（未本地运行）：** 需要访问 Ensembl 在线服务。成功时该基因的 `ensembl_gene_name` 为 `CDK2`，但本文不附固定输出。

## 3) `dnakit.databases.fetch_ensembl_sequences()`

- **作用：** 获取基因或转录本的核苷酸序列，或其蛋白质翻译，并解析为记录。
- **处理方法：** 调用 gget 的 `seq` 接口。gget 返回的是 FASTA 行的交替列表（标题、序列），本函数将它们配对，并检查格式是否正确。
- **API：** `dnakit.databases.fetch_ensembl_sequences(ensembl_ids, *, translate=False, isoforms=False, max_records=100)`。
- **输入：** Ensembl ID 列表；`translate=True` 返回蛋白质序列；`isoforms=True` 返回所有转录本。
- **输出：** 含 `count` 和 `records` 的字典，每条记录包含 `id`、`description` 和 `sequence`。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import fetch_ensembl_sequences

result = fetch_ensembl_sequences(["ENSG00000123374"], translate=True)
print(result["records"][0]["id"], result["records"][0]["sequence"][:10])
```

- **示例输出（未本地运行）：** 需要访问 Ensembl 或 UniProt 在线服务，因此不附固定输出。成功时 `records` 中每项的 `id` 为 Ensembl ID，`sequence` 为蛋白质序列。

## 4) `dnakit.databases.apply_cdna_mutations()`

- **作用：** 对每条序列应用一个 `c.` 表示法的突变，并返回突变位点两侧的序列窗口。
- **处理方法：** 逐条调用 gget 的 `mutate`，每次只处理一个序列和一个突变，这样可以准确对应输入和输出。突变位置与参考碱基不符时，gget 会丢弃该突变，本函数将其标记为 `dropped`。
- **API：** `dnakit.databases.apply_cdna_mutations(sequences, mutations, *, flank=30)`。
- **输入：** `sequences` 是序列 ID 到 DNA 字符串的映射；`mutations` 是序列 ID 到一个 `c.` 突变的映射（如 `"c.2C>T"`，位置从 1 开始）；`flank` 为突变两侧保留的碱基数。
- **输出：** `mutations` 列表，每项包含 `seq_id`、`mutation`、`status`（`applied` 或 `dropped`）和 `mutant_window`；另有 `applied` 计数。
- **示例代码：**

```python
from dnakit.databases import apply_cdna_mutations

result = apply_cdna_mutations({"s1": "ACGTTGCA" * 5}, {"s1": "c.2C>T"}, flank=3)
print(result["mutations"][0]["status"], result["mutations"][0]["mutant_window"])
```

- **示例输出：**

```text
applied ATGTT
```

## 5) `dnakit.databases.index_gene_models()`

- **作用：** 用 pyensembl 为本地 Ensembl 风格的 GTF 建立索引，并统计染色体、基因、转录本和外显子数量。
- **处理方法：** pyensembl 把 GTF 解析为缓存目录中的 SQLite 索引。索引名包含 GTF 内容的 SHA-256 摘要前 16 位，所以修改过的 GTF 会得到新索引，不会误用旧索引。
- **API：** `dnakit.databases.index_gene_models(gtf, *, cache_dir=None)`。
- **输入：** GTF 文件路径；`cache_dir` 是索引缓存目录，省略时使用 pyensembl 的默认缓存目录。
- **输出：** 含 `contigs`、`genes`、`transcripts`、`exons` 计数和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.databases import index_gene_models

summary = index_gene_models("genes.gtf", cache_dir="ensembl_cache")
print(summary["contigs"], summary["genes"], summary["transcripts"], summary["exons"])
```

- **示例输出：**

```text
1 2 2 3
```

## 6) `dnakit.databases.query_genes_by_name()`

- **作用：** 按精确基因名查找本地 GTF 中的基因。
- **处理方法：** 调用 pyensembl 的 `gene_ids_of_gene_name`。pyensembl 对没有匹配的名称抛出 `ValueError`，本函数将其转换为空列表。
- **API：** `dnakit.databases.query_genes_by_name(gtf, gene_names, *, cache_dir=None, max_results=100)`。
- **输入：** GTF 路径、基因名列表和 `max_results`（匹配总数上限）。
- **输出：** `results` 字典，键为查询的基因名，值为基因记录列表（含 `gene_id`、`contig`、`start`、`end`、`strand`、`biotype` 和 `transcript_ids`）；另有 `matches` 计数。
- **示例代码：**

```python
from dnakit.databases import query_genes_by_name

result = query_genes_by_name("genes.gtf", ["GENEA", "NOPE"], cache_dir="ensembl_cache")
print(result["results"]["GENEA"][0]["gene_id"], result["results"]["NOPE"])
```

- **示例输出：**

```text
ENSG00000000001 []
```

## 7) `dnakit.databases.query_genes_at_locus()`

- **作用：** 列出与一个基因组位点重叠的基因。
- **处理方法：** 调用 pyensembl 的 `genes_at_locus`。位点坐标为 1-based 闭区间；给出 `end` 时，返回与整个区间重叠的基因。`strand` 可限制为 `"+"` 或 `"-"`。
- **API：** `dnakit.databases.query_genes_at_locus(gtf, contig, start, end=None, *, strand=None, cache_dir=None, max_results=100)`。
- **输入：** GTF 路径、染色体名、起始位置、可选终止位置和可选链方向。染色体名不在注释中时报错。
- **输出：** 含 `count` 和 `genes` 的字典，`genes` 中每项与 `query_genes_by_name()` 的基因记录格式相同。
- **示例代码：**

```python
from dnakit.databases import query_genes_at_locus

result = query_genes_at_locus("genes.gtf", "chr1", 190, 310, cache_dir="ensembl_cache")
print([gene["gene_name"] for gene in result["genes"]])
```

- **示例输出：**

```text
['GENEA', 'GENEB']
```

## 8) `dnakit.databases.fetch_gene_sequence()`

- **作用：** 从本地基因组 FASTA 中取出一个基因的基因组位点序列，按基因链方向排列。
- **处理方法：** 先用 pyensembl 查到基因坐标，再调用 `Genome.sequence()` 取序列。负链基因返回反向互补序列。序列长度超过 `max_length` 时报错。
- **API：** `dnakit.databases.fetch_gene_sequence(gtf, genome_fasta, gene_id, *, cache_dir=None, max_length=1000000)`。
- **输入：** GTF 路径、基因组 FASTA 路径、Ensembl 基因 ID 和长度上限。
- **输出：** 含 `gene_id`、`gene_name`、`contig`、`start`、`end`、`strand`、`length` 和 `sequence` 的字典。
- **示例代码：**

```python
from dnakit.databases import fetch_gene_sequence

result = fetch_gene_sequence("genes.gtf", "genome.fa", "ENSG00000000002", cache_dir="ensembl_cache")
print(result["length"], result["sequence"][:20])
```

- **示例输出：**

```text
80 ACGTTGCAACGTTGCAACGT
```

## 9) `dnakit.databases.read_genome_summary()`

- **作用：** 统计本地基因组 FASTA 的染色体、长度和 N 碱基数。
- **处理方法：** 调用 genomepy 的 `Genome`。它读取 pyfaidx 索引，并在 FASTA 旁边写入 `.fai`、`.fa.sizes` 和 `.gaps.bed` 辅助文件。FASTA 本身不会被修改。README 中的 `na` 元数据转换为 `None`。
- **API：** `dnakit.databases.read_genome_summary(fasta, *, max_contigs=100000)`。
- **输入：** FASTA 路径；`max_contigs` 是记录数上限。
- **输出：** 含 `name`、`tax_id`、`assembly_accession`、`contig_count`、`total_length` 和 `contigs` 的字典。
- **示例代码：**

```python
from dnakit.databases import read_genome_summary

summary = read_genome_summary("genome.fa")
print(summary["name"], summary["contig_count"], summary["total_length"])
```

- **示例输出：**

```text
genome 1 400
```

## 10) `dnakit.databases.fetch_spliced_sequence()`

- **作用：** 把一条染色体上的若干区间按顺序拼接为一条序列，常用于外显子拼接。
- **处理方法：** 调用 genomepy（继承自 pyfaidx）的 `get_spliced_seq`。区间为 1-based 闭区间，因此 `[11, 60]` 包含第 11 到第 60 个碱基。`reverse_complement=True` 时，先拼接再取反向互补。
- **API：** `dnakit.databases.fetch_spliced_sequence(fasta, contig, intervals, *, reverse_complement=False, max_length=1000000)`。
- **输入：** FASTA 路径、染色体名、`[start, end]` 区间列表和长度上限。区间超出染色体长度时报错。
- **输出：** 含 `contig`、`intervals`、`length` 和 `sequence` 的字典。
- **示例代码：**

```python
from dnakit.databases import fetch_spliced_sequence

result = fetch_spliced_sequence("genome.fa", "chr1", [[11, 60], [121, 200]])
print(result["length"], result["sequence"][:20])
```

- **示例输出：**

```text
130 GTTGCAACGTTGCAACGTTG
```

## 11) `dnakit.databases.write_bed_sequences()`

- **作用：** 按 BED 区间从基因组 FASTA 中取出序列，并写入一个新的 FASTA 文件。
- **处理方法：** 调用 genomepy 的 `track2fasta`。BED 为 0-based 半开区间；`stranded=True` 时负链区间返回反向互补序列。`output` 已存在时报错，不会覆盖。
- **API：** `dnakit.databases.write_bed_sequences(fasta, bed, output, *, stranded=False, extend_up=0, extend_down=0, max_regions=100000)`。
- **输入：** FASTA 路径、BED 路径、输出路径；`extend_up` 和 `extend_down` 按链方向向两端延伸的碱基数。
- **输出：** 含 `regions`、`records`、`skipped`、`output` 和 `provenance` 的字典。
- **示例代码：**

```python
from pathlib import Path
from dnakit.databases import write_bed_sequences

output = Path("regions.fa")
output.unlink(missing_ok=True)  # the example is re-runnable; the function itself never overwrites
result = write_bed_sequences("genome.fa", "regions.bed", output, stranded=True)
print(result["records"], output.read_text().splitlines()[1])
```

- **示例输出：**

```text
2 GTTGCAACGT
```

## 12) `dnakit.databases.install_genome_assembly()`

- **作用：** 通过 genomepy 下载并安装一个基因组装配，可选包括注释文件。
- **处理方法：** 调用 `genomepy.install_genome`。它会从所选提供方下载大型文件，耗时和磁盘占用都较大。
- **API：** `dnakit.databases.install_genome_assembly(name, *, provider=None, genomes_dir=None, annotation=False, threads=1)`。
- **输入：** 装配名（如 `hg38`）、提供方（如 `ucsc`、`ensembl`）、安装目录、是否下载注释和线程数（1 到 64）。
- **输出：** 含 `name`、`provider`、`fasta` 和 `annotation_gtf` 的字典。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import install_genome_assembly

result = install_genome_assembly("hg38", provider="ucsc", annotation=True)
print(result["fasta"])
```

- **示例输出（未本地运行）：** 需要下载数 GB 的基因组文件，本文不在本地运行。成功时返回的 `fasta` 指向安装目录中的 FASTA 文件。

## 13) `dnakit.databases.format_hgvs_id()`

- **作用：** 把一条 VCF 风格的变异格式化为基因组 HGVS 标识。
- **处理方法：** 调用 biothings_client 变异客户端的 `format_hgvs`。该方法纯粹做字符串格式化，不访问网络。HGVS 标识不包含基因组版本，需与后续 MyVariant 查询的 `assembly` 保持一致。
- **API：** `dnakit.databases.format_hgvs_id(chrom, pos, ref, alt)`。
- **输入：** 人类染色体名（可带 `chr` 前缀）、1-based 位置、参考碱基和替代碱基，只允许 A、C、G、T、N。
- **输出：** 含 `hgvs_id`、`chrom`、`pos`、`ref` 和 `alt` 的字典。
- **示例代码：**

```python
from dnakit.databases import format_hgvs_id

print(format_hgvs_id("chr9", 107620835, "G", "A")["hgvs_id"])
```

- **示例输出：**

```text
chr9:g.107620835G>A
```

## 14) `dnakit.databases.hgvs_ids_from_vcf()`

- **作用：** 为 VCF 中每个序列替代等位基因生成基因组 HGVS 标识。
- **处理方法：** 逐行解析 VCF，把每个替代等位基因交给 `format_hgvs`。跳过的记录会计入 `skipped`：染色体不是人类染色体、替代等位基因为 `.` 或结构变异，或与参考相同。上游的 VCF 迭代器不做这些检查，所以这里不直接使用它。
- **API：** `dnakit.databases.hgvs_ids_from_vcf(path, *, max_records=100000)`。
- **输入：** VCF 路径和记录数上限。
- **输出：** 含 `records`、`hgvs_ids` 和 `skipped` 的字典。
- **示例代码：**

```python
from dnakit.databases import hgvs_ids_from_vcf

result = hgvs_ids_from_vcf("sample.vcf")
print(result["hgvs_ids"])
print(result["skipped"])
```

- **示例输出：**

```text
['chr1:g.10A>G', 'chr1:g.10A>T', 'chr2:g.31del']
{'unsupported_contig': 1, 'non_sequence_allele': 1, 'same_as_ref': 0}
```

## 15) `dnakit.databases.query_mygene_genes()`

- **作用：** 用查询字符串在 MyGene.info 中搜索基因记录。
- **处理方法：** 调用 mygene 的 `query` 接口，返回总数和命中记录，并限定返回条数。
- **API：** `dnakit.databases.query_mygene_genes(query, *, species="human", fields="symbol,name,entrezgene,ensembl.gene,taxid", size=10)`。
- **输入：** MyGene 查询语法字符串，如 `"symbol:CDK2"`；`species`、字段列表和返回条数（1 到 1000）。
- **输出：** 含 `total` 和 `hits` 的字典。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import query_mygene_genes

result = query_mygene_genes("symbol:CDK2", size=1)
print(result["hits"][0]["symbol"])
```

- **示例输出（未本地运行）：** 需要访问 MyGene.info 在线服务。成功时 `hits` 的第一项 `symbol` 为 `CDK2`。

## 16) `dnakit.databases.fetch_mygene_gene()`

- **作用：** 按 Entrez 或 Ensembl ID 获取一条 MyGene.info 基因记录。
- **处理方法：** 调用 mygene 的 `getgene`。未找到时返回 `found=False`，而不是报错。
- **API：** `dnakit.databases.fetch_mygene_gene(gene_id, *, species="human", fields="symbol,name,entrezgene,ensembl.gene,taxid")`。
- **输入：** 基因 ID（如 `"1017"`）和物种。
- **输出：** 含 `gene_id`、`found` 和 `gene` 的字典。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import fetch_mygene_gene

result = fetch_mygene_gene("1017")
print(result["found"], result["gene"]["symbol"])
```

- **示例输出（未本地运行）：** 需要访问 MyGene.info 在线服务。成功时输出 `True` 和 `CDK2`。

## 17) `dnakit.databases.map_gene_identifiers()`

- **作用：** 批量把一种标识符（如 Entrez ID）映射为 MyGene.info 基因记录。
- **处理方法：** 调用 mygene 的 `querymany`，每个输入保留 `query` 字段，便于对应。
- **API：** `dnakit.databases.map_gene_identifiers(identifiers, *, scopes="entrezgene", species="human", fields="symbol,name,taxid")`。
- **输入：** 标识符列表（最多 1000 个）、`scopes`（标识符类型）和物种。
- **输出：** 含 `requested` 和 `results` 的字典。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import map_gene_identifiers

result = map_gene_identifiers(["1017", "7157"])
print([row["symbol"] for row in result["results"]])
```

- **示例输出（未本地运行）：** 需要访问 MyGene.info 在线服务。成功时输出的符号列表包含 `CDK2` 和 `TP53`。

## 18) `dnakit.databases.fetch_variant_annotations()`

- **作用：** 获取一组基因组 HGVS 标识的 MyVariant.info 注释。
- **处理方法：** 调用 myvariant 的 `getvariants`，把 `assembly` 作为请求参数传给服务。
- **API：** `dnakit.databases.fetch_variant_annotations(hgvs_ids, *, assembly="hg38", fields="dbsnp.rsid,clinvar.rcv.clinical_significance,cadd.phred,dbnsfp.genename")`。
- **输入：** HGVS 标识列表（最多 1000 个），格式如 `"chr9:g.107620835G>A"`；`assembly` 取 `"hg19"` 或 `"hg38"`。
- **输出：** 含 `requested`、`results` 和 `provenance` 的字典。`results` 中未找到的标识带有 `notfound` 标记。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import fetch_variant_annotations

result = fetch_variant_annotations(["chr9:g.107620835G>A"], assembly="hg19")
print(result["results"][0]["_id"])
```

- **示例输出（未本地运行）：** 需要访问 MyVariant.info 在线服务。成功时输出 `chr9:g.107620835G>A`。

## 19) `dnakit.databases.query_variant_annotations()`

- **作用：** 用查询字符串在 MyVariant.info 中搜索变异记录。
- **处理方法：** 调用 myvariant 的 `query` 接口，返回总数和命中记录，并限定返回条数。
- **API：** `dnakit.databases.query_variant_annotations(query, *, assembly="hg38", fields=..., size=10)`。
- **输入：** 查询字符串，如 `"dbsnp.rsid:rs58991260"`；`assembly` 和返回条数（1 到 1000）。
- **输出：** 含 `total` 和 `hits` 的字典。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import query_variant_annotations

result = query_variant_annotations("dbsnp.rsid:rs58991260", size=1)
print(result["total"])
```

- **示例输出（未本地运行）：** 需要访问 MyVariant.info 在线服务，结果数随数据库更新而变化，因此不附固定输出。

## 20) `dnakit.databases.query_biomart_dataset()`

- **作用：** 查询一个 BioMart 数据集，返回指定属性和过滤条件下的表格。
- **处理方法：** 调用 pybiomart 的 `Dataset.query`，它把请求组装为 BioMart XML 并解析 TSV 结果。空单元格记为 `None`。pybiomart 会在导入时在当前工作目录创建 `.pybiomart.sqlite` 请求缓存。
- **API：** `dnakit.databases.query_biomart_dataset(dataset, attributes, *, filters=None, host="http://www.ensembl.org", only_unique=True, use_attr_names=False, max_rows=10000)`。
- **输入：** 数据集名（如 `hsapiens_gene_ensembl`）、属性列表、可选过滤器字典、服务地址和行数上限。服务器返回的完整结果会先被读取，`max_rows` 在读取后检查。
- **输出：** 含 `columns`、`row_count` 和 `rows` 的字典。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import query_biomart_dataset

result = query_biomart_dataset(
    "hsapiens_gene_ensembl",
    ["ensembl_gene_id", "external_gene_name"],
    filters={"chromosome_name": "21"},
    max_rows=5,
)
print(result["rows"][0])
```

- **示例输出（未本地运行）：** 需要访问 Ensembl BioMart 在线服务。成功时 `rows` 的每项是一条包含两列的字典。

## 21) `dnakit.databases.list_biomart_attributes()`

- **作用：** 列出一个 BioMart 数据集可用的属性。
- **处理方法：** 调用 pybiomart 的 `Dataset.list_attributes`，返回属性名、显示名和描述。
- **API：** `dnakit.databases.list_biomart_attributes(dataset, *, host="http://www.ensembl.org", max_rows=5000)`。
- **输入：** 数据集名和服务地址。
- **输出：** 含 `columns`、`count` 和 `attributes` 的字典。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import list_biomart_attributes

result = list_biomart_attributes("hsapiens_gene_ensembl")
print(result["count"])
```

- **示例输出（未本地运行）：** 需要访问 Ensembl BioMart 在线服务。属性数量随 Ensembl 版本变化，因此不附固定输出。

## 22) `dnakit.databases.list_biomart_datasets()`

- **作用：** 列出一个 BioMart mart 中的数据集。
- **处理方法：** 调用 pybiomart 的 `Server.marts` 找到 mart，再调用 `Mart.list_datasets`。mart 名不在服务列表中时报错。
- **API：** `dnakit.databases.list_biomart_datasets(mart="ENSEMBL_MART_ENSEMBL", *, host="http://www.ensembl.org", max_rows=1000)`。
- **输入：** mart 名、服务地址和行数上限。
- **输出：** 含 `mart`、`columns`、`count` 和 `datasets` 的字典。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import list_biomart_datasets

result = list_biomart_datasets(max_rows=5)
print(result["datasets"][0]["name"])
```

- **示例输出（未本地运行）：** 需要访问 Ensembl BioMart 在线服务，结果随版本变化，因此不附固定输出。

## 23) `dnakit.databases.fetch_ncbi_gene_by_id()`

- **作用：** 按 NCBI Gene ID 获取基因元数据。
- **处理方法：** 调用 ncbi-datasets-pylib 的 `GeneApi.gene_metadata_by_id`，并把服务端返回条数限制为 `max_genes`。返回的 openapi 模型通过 `to_dict()` 转换为普通字典。
- **API：** `dnakit.databases.fetch_ncbi_gene_by_id(gene_ids, *, max_genes=100)`。
- **输入：** 正整数 Gene ID 列表，最多 200 个；`max_genes` 为返回条数上限（1 到 1000）。
- **输出：** 含 `requested`、`total_count` 和 `genes` 的字典，`genes` 中每项为 NCBI 的 gene match 记录。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import fetch_ncbi_gene_by_id

result = fetch_ncbi_gene_by_id([1017], max_genes=1)
print(result["genes"][0]["gene"]["symbol"])
```

- **示例输出（未本地运行）：** 需要访问 NCBI Datasets 在线服务，并且需要安装 `ncbi-datasets-pylib`。成功时输出基因符号 `CDK2`。

## 24) `dnakit.databases.fetch_ncbi_gene_by_symbol()`

- **作用：** 按基因符号和物种获取 NCBI 基因元数据。
- **处理方法：** 调用 `GeneApi.gene_metadata_by_tax_and_symbol`，符号和物种名先经过格式检查。
- **API：** `dnakit.databases.fetch_ncbi_gene_by_symbol(symbols, taxon, *, max_genes=100)`。
- **输入：** 基因符号列表、物种名（如 `"human"` 或 `"9606"`）和返回条数上限。
- **输出：** 含 `symbols`、`taxon`、`total_count` 和 `genes` 的字典。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import fetch_ncbi_gene_by_symbol

result = fetch_ncbi_gene_by_symbol(["CDK2"], "human", max_genes=1)
print(result["genes"][0]["gene"]["gene_id"])
```

- **示例输出（未本地运行）：** 需要访问 NCBI Datasets 在线服务，并且需要安装 `ncbi-datasets-pylib`。成功时输出 `CDK2` 的 Gene ID。

## 25) `dnakit.databases.fetch_ncbi_gene_orthologs()`

- **作用：** 获取一个 NCBI Gene ID 的直系同源基因集。
- **处理方法：** 调用 `GeneApi.gene_orthologs_by_id`。`taxon_filter` 省略时不做物种过滤。
- **API：** `dnakit.databases.fetch_ncbi_gene_orthologs(gene_id, *, taxon_filter=None)`。
- **输入：** 一个正整数 Gene ID，以及可选的物种 ID 列表（如 `["9606"]`）。
- **输出：** 含 `gene_id`、`taxon_filter`、`ortholog_set_id`、`method` 和 `genes` 的字典。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import fetch_ncbi_gene_orthologs

result = fetch_ncbi_gene_orthologs(1017, taxon_filter=["9606"])
print(result["method"], len(result["genes"]))
```

- **示例输出（未本地运行）：** 需要访问 NCBI Datasets 在线服务，并且需要安装 `ncbi-datasets-pylib`。输出为同源方法名和同源基因数量，本文不附固定输出。

## 验证方式

- 测试文件是 `others/tests/unit/test_database_gene_backends.py`，共 37 个测试：33 个通过，4 个跳过。跳过的测试是两个 pybiomart 测试（依赖冲突，见下文）、一个 ncbi-datasets-pylib 签名测试（包未安装）和一个需要 `DNAKIT_LIVE_TESTS=1` 的在线测试。
- 本地测试直接运行上游包：pyensembl 基于临时目录中的 GTF 建立索引，并与直接调用 `pyensembl.Genome` 的计数比较；genomepy 的拼接结果与直接调用 `get_seq` 的结果比较；gget.mutate 的输出与直接调用结果比较；biothings_client 的 HGVS 格式与上游 `format_hgvs` 和 `get_hgvs_from_vcf` 比较。
- 另外的基因序列、外显子拼接、突变窗口和 HGVS 结果都与人工计算的切片比较，例如 GENEA 为 190 bp，GENEB 为反向互补序列。
- 网络函数替换上游客户端的调用，测试解析、缺失值处理和参数校验，不访问网络。ncbi 函数用一个假的 `ncbi.datasets` 模块测试调用参数、`limit` 限制和错误转换，真实包的方法签名测试在包可导入时运行。
- 缺少任一可选包时，函数抛出 `BackendUnavailableError`，并给出对应的 `bio-*` 额外依赖名。
- 设置 `DNAKIT_LIVE_TESTS=1` 时运行一个真实在线测试，检查 CDK2 的 Ensembl 和 MyGene 结果。
- 参数校验覆盖布尔值、越界数字、非法 ID、过长序列、已存在的输出文件和未知染色体等情况。

## 当前范围

- 已覆盖：gget 的 Ensembl 搜索、基因信息、序列和 cDNA 突变窗口；pyensembl 的本地 GTF 索引、基因名与位点查询、基因组位点序列；genomepy 的 FASTA 摘要、外显子拼接、BED 序列输出和装配安装；biothings_client 的 HGVS 格式化；mygene 与 myvariant 的查询；pybiomart 的属性、数据集和查询。
- 未包装：gget 的 `ref`（只返回 FTP 链接）、`blat`、`blast`、`cbio`、`cellxgene` 等在线分析接口；pyensembl 的 `EnsemblRelease` 下载流程，本模块只使用显式给出的本地 GTF 和 FASTA；biothings_client 的基因和变异客户端与 mygene、myvariant 是同一服务，因此只包装其独有的 HGVS 工具。
- 依赖（extras）：`bio-gget` = `gget>=0.30,<1`；`bio-pyensembl` = `pyensembl>=2.27,<3`；`bio-genomepy` = `genomepy>=0.16,<1`；`bio-biothings` = `biothings_client>=0.5,<1`；`bio-mygene` = `mygene>=3.2,<4`；`bio-myvariant` = `myvariant>=1.0,<2`；`bio-pybiomart` = `pybiomart>=0.2,<1`；`bio-ncbi-datasets` = `ncbi-datasets-pylib>=16.6,<17`。
- ncbi-datasets-pylib（16.6.1）在当前 venv 中无法安装：PyPI 只发布 sdist，构建时缺少 `pkg_resources`。因此第 23 至 25 个函数的接口依据该 sdist 的源码编写（`GeneApi` 的方法签名和 `limit` 参数），没有对真实包执行。
- pybiomart 在当前 venv 中不可导入：`cattrs` 26.2.1 需要更新的 `attrs`，而当前安装的是 attrs 21.4.0。这三个 BioMart 函数的测试在该环境中跳过。
- 在线服务的结果随时间变化，因此相关示例未在本地运行，文档中也不附输出。
