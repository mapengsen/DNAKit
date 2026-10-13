# pyensembl、genomepy 与 NCBI Datasets 的基因组函数

这一组函数包装 [pyensembl](https://github.com/openvax/pyensembl)、[genomepy](https://github.com/vanheeringen-lab/genomepy) 和 [ncbi-datasets-pylib](https://github.com/ncbi/datasets) 中的 DNA 序列、基因组区间、注释文件与基因标识符操作；多数函数只读取调用者给出的本地文件或记录，不下载数据，也不做比对或基因型分析。NCBI Datasets 的装配体查询需要网络，本地测试使用模拟客户端。[pybiomart](https://github.com/jrowen/pybiomart) 没有模块级公开函数，因此本组不包含它。依赖都是可选的，额外安装项分别为 `bio-pyensembl`、`bio-genomepy` 和 `bio-ncbi-datasets`。

## 示例数据

**示例文件：** `toy.fa`

```text
>chr1 test
ACGT
GG
>chr2
TTTT
```

小型 FASTA，两条序列。

**示例文件：** `toy.gtf`

```text
chr1	src	gene	10	100	.	+	.	gene_id "G1";
chr2	src	gene	5	50	.	-	.	gene_id "G2";
```

两行 GTF 注释（制表符分隔）。

**示例文件：** `toy_n.fa`

```text
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNACGTACGTACGTACGTACGT
>chr2 desc
GGGGCCCCAAAATTTT
>scaf_9
NNNNACGTAC
```

含 N 段的 FASTA：chr1 中部有 4 个 N，scaf_9 开头有 4 个 N。

## 1) `dnakit.databases.read_fasta_records()`

- **作用：** 读取 FASTA 或 gzip 压缩的 FASTA，返回按标识符排序的序列记录。
- **处理方法：** 调用 `pyensembl.fasta.parse_fasta_dictionary`。读取前检查文件大小不超过 `max_bytes`（默认 64 MB）；不计算序列统计。
- **API：** `dnakit.databases.read_fasta_records(fasta_path, *, max_bytes=64000000, verbose=True)`。
- **输入：** FASTA 路径（可为 `.gz`）；`max_bytes` 为文件大小上限。
- **输出：** 含 `records`（`id`、`length`、`sequence`）、`record_count` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.databases import read_fasta_records

result = read_fasta_records("toy.fa", verbose=False)
print(result["record_count"])
for record in result["records"]:
    print(record["id"], record["length"], record["sequence"])
```

- **示例输出：**

```text
2
chr1 6 ACGTGG
chr2 4 TTTT
```

## 2) `dnakit.databases.read_first_fasta_id()`

- **作用：** 返回 FASTA 第一条记录的标识符；没有记录时返回 `None`。
- **处理方法：** 调用 `pyensembl.fasta.first_fasta_id`，读到第一行 `>` 标题为止。
- **API：** `dnakit.databases.read_first_fasta_id(fasta_path, *, verbose=True)`。
- **输入：** FASTA 路径（可为 `.gz`）。
- **输出：** 含 `identifier`（字符串或 `None`）的字典。
- **示例代码：**

```python
from dnakit.databases import read_first_fasta_id

print(read_first_fasta_id("toy.fa", verbose=False)["identifier"])
```

- **示例输出：**

```text
chr1
```

## 3) `dnakit.databases.merge_genomic_intervals()`

- **作用：** 合并重叠或首尾相接的 1-based 闭区间。
- **处理方法：** 调用 `pyensembl.common.merge_intervals`。输入先排序；`(a, b)` 与 `(b + 1, c)` 视为相接并合并。
- **API：** `dnakit.databases.merge_genomic_intervals(intervals, *, max_intervals=1000000, verbose=True)`。
- **输入：** `(start, end)` 二元组列表，坐标从 1 开始，`start` 不大于 `end`；`max_intervals` 限制数量。
- **输出：** `intervals`（按起点升序的 `start`、`end` 字典）、`input_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.databases import merge_genomic_intervals

result = merge_genomic_intervals([(5, 9), (1, 3), (4, 4), (12, 15)], verbose=False)
print(result["input_count"])
for item in result["intervals"]:
    print(item["start"], item["end"])
```

- **示例输出：**

```text
4
1 9
12 15
```

## 4) `dnakit.databases.find_nearest_genomic_locus()`

- **作用：** 找出距离给定区间最近的位点。
- **处理方法：** 用 `pyensembl.locus.Locus` 构造位点，再调用 `pyensembl.search.find_nearest_locus` 计算距离。距离相同时取输入中靠前的位点。
- **API：** `dnakit.databases.find_nearest_genomic_locus(start, end, loci, *, max_loci=100000, verbose=True)`。
- **输入：** 查询区间 `start`、`end`（1-based 闭区间）；`loci` 为含 `contig`、`start`、`end`、`strand` 的字典列表。
- **输出：** `distance`（没有位点时为 `None`）、`nearest_index`、`nearest` 和 `provenance`。
- **示例代码：**

```python
from dnakit.databases import find_nearest_genomic_locus

loci = [
    {"contig": "1", "start": 100, "end": 200, "strand": "+"},
    {"contig": "2", "start": 500, "end": 600, "strand": "-"},
]
result = find_nearest_genomic_locus(300, 310, loci, verbose=False)
print(result["distance"], result["nearest_index"], result["nearest"]["contig"])
```

- **示例输出：**

```text
100 0 1
```

## 5) `dnakit.databases.stable_ensembl_id()`

- **作用：** 去掉 Ensembl 标识符末尾的版本号，例如 `.15`。
- **处理方法：** 调用 `pyensembl.gene_name_aliases.stable_ensembl_gene_id`。不符合 Ensembl 格式的标识符保持原样。
- **API：** `dnakit.databases.stable_ensembl_id(identifier, *, verbose=True)`。
- **输入：** 标识符字符串。
- **输出：** `identifier` 与 `stable_id`。
- **示例代码：**

```python
from dnakit.databases import stable_ensembl_id

print(stable_ensembl_id("ENSG00000139618.15", verbose=False)["stable_id"])
print(stable_ensembl_id("AT1G01010.1", verbose=False)["stable_id"])
```

- **示例输出：**

```text
ENSG00000139618
AT1G01010.1
```

## 6) `dnakit.databases.normalize_chromosome_name()`

- **作用：** 规范化染色体或序列名，例如 `chrx` 变为 `chrX`，整数 `7` 变为 `"7"`。
- **处理方法：** 调用 `pyensembl.normalization.normalize_chromosome`。空字符串和 `"0"` 报 `ConfigurationError`。
- **API：** `dnakit.databases.normalize_chromosome_name(name, *, verbose=True)`。
- **输入：** 字符串或整数（不接受布尔值）。
- **输出：** `input` 与 `chromosome`。
- **示例代码：**

```python
from dnakit.databases import normalize_chromosome_name

for name in ["chrx", "x", 7, "chrUn_gl000212"]:
    print(normalize_chromosome_name(name, verbose=False)["chromosome"])
```

- **示例输出：**

```text
chrX
X
7
chrUn_gl000212
```

## 7) `dnakit.databases.normalize_strand_symbol()`

- **作用：** 把链方向统一写成 `+` 或 `-`。
- **处理方法：** 调用 `pyensembl.normalization.normalize_strand`，接受 `+`、`-`、`1`、`-1`、`1` 和 `-1`（整数）。
- **API：** `dnakit.databases.normalize_strand_symbol(strand, *, verbose=True)`。
- **输入：** 字符串或整数（不接受布尔值）。
- **输出：** `input` 与 `strand`。
- **示例代码：**

```python
from dnakit.databases import normalize_strand_symbol

print(normalize_strand_symbol("1", verbose=False)["strand"])
print(normalize_strand_symbol(-1, verbose=False)["strand"])
```

- **示例输出：**

```text
+
-
```

## 8) `dnakit.databases.normalize_gene_aliases()`

- **作用：** 规范化基因名到基因 ID 的映射。
- **处理方法：** 调用 `pyensembl.gene_name_aliases.normalize_aliases`。单个字符串 ID 视为一个值；每个基因名下的 ID 去重并排序。
- **API：** `dnakit.databases.normalize_gene_aliases(aliases, *, max_names=100000, verbose=True)`。
- **输入：** 映射：键为基因名，值为 ID 字符串或 ID 列表；`max_names` 限制基因名数量。
- **输出：** `records`（`name`、`gene_ids`，按基因名排序）、`name_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.databases import normalize_gene_aliases

result = normalize_gene_aliases({"BRCA2": "ENSG1.2", "X": ["b", "a", "a"]}, verbose=False)
for record in result["records"]:
    print(record["name"], record["gene_ids"])
```

- **示例输出：**

```text
BRCA2 ['ENSG1.2']
X ['a', 'b']
```

## 9) `dnakit.databases.match_identifier_version()`

- **作用：** 在已安装的 Ensembl 标识符中找出与给定标识符对应的形式。
- **处理方法：** 先只比较稳定 ID（去掉版本号），再调用 `pyensembl.versioned_ids.match_version`。稳定 ID 未安装时返回 `None`，不会匹配到其他基因。
- **API：** `dnakit.databases.match_identifier_version(identifier, installed, *, verbose=True)`。
- **输入：** 标识符；`installed` 为 `{已安装标识符: 记录的版本或 None}`。
- **输出：** `matched`（匹配到的已安装形式或 `None`）、`installed_count`、`compared_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.databases import match_identifier_version

print(match_identifier_version("ENSG1.2", {"ENSG1": 2}, verbose=False)["matched"])
print(match_identifier_version("ENSG9", {"ENSG1": 2}, verbose=False)["matched"])
```

- **示例输出：**

```text
ENSG1
None
```

## 10) `dnakit.databases.lookup_sequence_version_fallback()`

- **作用：** 按标识符查找序列，并容忍 `.N` 版本号的差异。
- **处理方法：** 调用 `pyensembl.sequence_data.lookup_sequence_with_version_fallback`。普通字典没有 pyensembl 为 FASTA 建立的版本索引，所以裸 ID 不能匹配带版本号的键。
- **API：** `dnakit.databases.lookup_sequence_version_fallback(sequences, identifier, *, version=None, verbose=True)`。
- **输入：** `{标识符: 序列}` 字典、标识符，以及可选的 `version`。
- **输出：** `found`、`sequence`、`identifier`、`version` 和 `provenance`。
- **示例代码：**

```python
from dnakit.databases import lookup_sequence_version_fallback

sequences = {"ENSP1.3": "MKV", "ENST2": "ACGT"}
print(lookup_sequence_version_fallback(sequences, "ENST2.5", verbose=False)["sequence"])
print(lookup_sequence_version_fallback(sequences, "ENSP1.3", verbose=False)["sequence"])
```

- **示例输出：**

```text
ACGT
MKV
```

## 11) `dnakit.databases.normalize_genome_reference_name()`

- **作用：** 把参考基因组名称规范为标准写法，例如 `grch38` 变为 `GRCh38`。
- **处理方法：** 调用 `pyensembl.reference_name.normalize_reference_name`，比较时忽略大小写；不在 pyensembl 登记表中的名称报 `ConfigurationError`。
- **API：** `dnakit.databases.normalize_genome_reference_name(name, *, verbose=True)`。
- **输入：** 参考基因组名称字符串。
- **输出：** `input` 与 `reference`。
- **示例代码：**

```python
from dnakit.databases import normalize_genome_reference_name

print(normalize_genome_reference_name("grch38", verbose=False)["reference"])
```

- **示例输出：**

```text
GRCh38
```

## 12) `dnakit.databases.ensembl_reference_name()`

- **作用：** 返回某物种在某个 Ensembl 版本中使用的参考基因组名称。
- **处理方法：** 调用 `pyensembl.reference_name.which_reference`，只读取 pyensembl 内置的登记表，不访问网络。
- **API：** `dnakit.databases.ensembl_reference_name(species_name, ensembl_release, *, verbose=True)`。
- **输入：** 物种名（如 `human`）；Ensembl 版本号（整数）或注释日期（如 `2026_04`）。
- **输出：** `species`、`ensembl_release`、`reference` 和 `provenance`。
- **示例代码：**

```python
from dnakit.databases import ensembl_reference_name

print(ensembl_reference_name("human", 75, verbose=False)["reference"])
print(ensembl_reference_name("human", 110, verbose=False)["reference"])
```

- **示例输出：**

```text
GRCh37
GRCh38
```

## 13) `dnakit.databases.ensembl_gtf_url()`

- **作用：** 生成某物种某个 Ensembl 版本的基因注释 GTF 下载地址。
- **处理方法：** 调用 `pyensembl.ensembl_url_templates.make_gtf_url`，只拼接字符串，不发起请求。
- **API：** `dnakit.databases.ensembl_gtf_url(ensembl_release, species, *, verbose=True)`。
- **输入：** Ensembl 版本号（整数）；物种名。
- **输出：** `url`、`species` 和 `ensembl_release`。
- **示例代码：**

```python
from dnakit.databases import ensembl_gtf_url

print(ensembl_gtf_url(110, "human", verbose=False)["url"])
```

- **示例输出：**

```text
https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.chr_patch_hapl_scaff.gtf.gz
```

## 14) `dnakit.databases.ensembl_genome_fasta_url()`

- **作用：** 生成参考 DNA FASTA 的下载地址。
- **处理方法：** 调用 `pyensembl.ensembl_url_templates.make_genome_fasta_url`。`fasta_type` 取 `toplevel` 或 `primary_assembly`；`mask` 取 `none`、`soft` 或 `hard`。只拼接字符串。
- **API：** `dnakit.databases.ensembl_genome_fasta_url(ensembl_release, species, *, fasta_type="toplevel", mask="none", verbose=True)`。
- **输入：** Ensembl 版本号、物种名、`fasta_type` 和 `mask`。
- **输出：** `url` 以及所选的 `fasta_type`、`mask`。
- **示例代码：**

```python
from dnakit.databases import ensembl_genome_fasta_url

result = ensembl_genome_fasta_url(110, "human", fasta_type="primary_assembly", mask="soft", verbose=False)
print(result["url"])
```

- **示例输出：**

```text
https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```

## 15) `dnakit.databases.ensembl_sequence_fasta_url()`

- **作用：** 生成 cDNA、ncRNA 或蛋白质 FASTA 的下载地址。
- **处理方法：** 调用 `pyensembl.ensembl_url_templates.make_fasta_url`。`sequence_type` 只接受 `cdna`、`ncrna` 和 `pep`。只拼接字符串。
- **API：** `dnakit.databases.ensembl_sequence_fasta_url(ensembl_release, species, sequence_type, *, verbose=True)`。
- **输入：** Ensembl 版本号、物种名和 `sequence_type`。
- **输出：** `url`、`sequence_type`、`species` 和 `ensembl_release`。
- **示例代码：**

```python
from dnakit.databases import ensembl_sequence_fasta_url

print(ensembl_sequence_fasta_url(110, "mouse", "pep", verbose=False)["url"])
```

- **示例输出：**

```text
https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz
```

## 16) `dnakit.databases.read_annotation_records()`

- **作用：** 读取 GTF（9 列）或 BED12（12 列）注释文件，每行返回一条记录。
- **处理方法：** 调用 `genomepy.annotation.utils.read_annot`，按列数判断格式；列数不是 9 或 12 时报 `ConfigurationError`。文件大小受 `max_bytes` 限制（默认 64 MB）。
- **API：** `dnakit.databases.read_annotation_records(annotation_path, *, max_bytes=64000000, max_records=1000000, verbose=True)`。
- **输入：** 注释文件路径；`max_bytes` 与 `max_records` 为上限。
- **输出：** `file_format`（`gtf` 或 `bed12`）、`columns`、`record_count`、`records` 和 `provenance`。
- **示例代码：**

```python
from dnakit.databases import read_annotation_records

result = read_annotation_records("toy.gtf", verbose=False)
print(result["file_format"], result["record_count"])
for record in result["records"]:
    print(record["seqname"], record["start"], record["end"], record["strand"])
```

- **示例输出：**

```text
gtf 2
chr1 10 100 +
chr2 5 50 -
```

## 17) `dnakit.databases.write_annotation_records()`

- **作用：** 把 GTF 或 BED12 记录写成新的制表符分隔文件（无表头）。
- **处理方法：** 检查每条记录的键与所选格式的列完全一致，然后调用 `genomepy.annotation.utils.write_annot`。输出文件必须不存在，不会覆盖已有文件。
- **API：** `dnakit.databases.write_annotation_records(records, output_path, *, file_format="gtf", max_records=1000000, verbose=True)`。
- **输入：** 记录列表；`output_path`；`file_format` 取 `gtf` 或 `bed12`。
- **输出：** `output_path`、`file_format`、`record_count` 和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.databases import write_annotation_records

Path("written.gtf").unlink(missing_ok=True)
record = {
    "seqname": "chr1",
    "source": "s",
    "feature": "gene",
    "start": 1,
    "end": 5,
    "score": ".",
    "strand": "+",
    "frame": ".",
    "attribute": 'gene_id "X";',
}
result = write_annotation_records([record], "written.gtf", verbose=False)
print(result["record_count"])
print(Path("written.gtf").read_text(encoding="utf-8"), end="")
```

- **示例输出：**

```text
1
chr1	s	gene	1	5	.	+	.	gene_id "X";
```

## 18) `dnakit.databases.filter_annotation_contigs()`

- **作用：** 按染色体名（contig）的正则表达式筛选注释记录。
- **处理方法：** 调用 `genomepy.annotation.filter_regex`，用 `re.match` 匹配第一列。`invert_match` 为真时保留不匹配的记录。
- **API：** `dnakit.databases.filter_annotation_contigs(records, pattern, *, file_format="gtf", invert_match=False, max_records=1000000, verbose=True)`。
- **输入：** 注释记录列表；`pattern`（正则表达式）；`file_format`；`invert_match`。
- **输出：** `records`（保留的记录，顺序不变）、`kept_count`、`removed_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.databases import filter_annotation_contigs, read_annotation_records

records = read_annotation_records("toy.gtf", verbose=False)["records"]
kept = filter_annotation_contigs(records, "chr2", verbose=False)
print(kept["kept_count"], kept["removed_count"], kept["records"][0]["seqname"])
```

- **示例输出：**

```text
1 1 chr2
```

## 19) `dnakit.databases.write_gap_bed()`

- **作用：** 把 FASTA 中连续的 N 段（gap）写成 BED 文件。
- **处理方法：** 调用 `genomepy.genome.generate_gap_bed`，经 pyfaidx 读取。坐标为 0-based 半开区间。pyfaidx 会在输入 FASTA 旁写入 `.fai` 索引。
- **API：** `dnakit.databases.write_gap_bed(fasta_path, output_path, *, verbose=True)`。
- **输入：** FASTA 路径；新的输出路径。
- **输出：** `output_path`、`gap_count` 和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.databases import write_gap_bed

Path("gaps.bed").unlink(missing_ok=True)
result = write_gap_bed("toy_n.fa", "gaps.bed", verbose=False)
print(result["gap_count"])
print(Path("gaps.bed").read_text(encoding="utf-8"), end="")
```

- **示例输出：**

```text
2
chr1	40	44
scaf_9	0	4
```

## 20) `dnakit.databases.write_fasta_sizes()`

- **作用：** 写出 FASTA 中每条序列的名称和长度。
- **处理方法：** 调用 `genomepy.genome.generate_fa_sizes`，输出两列、制表符分隔的 sizes 文件。
- **API：** `dnakit.databases.write_fasta_sizes(fasta_path, output_path, *, verbose=True)`。
- **输入：** FASTA 路径；新的输出路径。
- **输出：** `output_path`、`sequence_count` 和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.databases import write_fasta_sizes

Path("sizes.txt").unlink(missing_ok=True)
result = write_fasta_sizes("toy_n.fa", "sizes.txt", verbose=False)
print(result["sequence_count"])
print(Path("sizes.txt").read_text(encoding="utf-8"), end="")
```

- **示例输出：**

```text
3
chr1	64
chr2	16
scaf_9	10
```

## 21) `dnakit.databases.filter_fasta_records()`

- **作用：** 按标题行的正则表达式筛选 FASTA 记录，写入新文件。
- **处理方法：** 调用 `genomepy.files.filter_fasta`，返回被删除的序列名（已排序）。输出文件必须与输入不同且不存在。
- **API：** `dnakit.databases.filter_fasta_records(fasta_path, output_path, *, pattern=".*", invert_match=False, verbose=True)`。
- **输入：** FASTA 路径；新的输出路径；`pattern`；`invert_match`。
- **输出：** `output_path`、`pattern`、`invert_match`、`kept_count`、`removed`（已排序）和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.databases import filter_fasta_records

Path("kept.fa").unlink(missing_ok=True)
result = filter_fasta_records("toy_n.fa", "kept.fa", pattern="chr", verbose=False)
print(result["removed"], result["kept_count"])
```

- **示例输出：**

```text
['scaf_9'] 2
```

## 22) `dnakit.databases.sample_random_regions()`

- **作用：** 从 FASTA 中按随机种子抽取固定长度的区间。
- **处理方法：** 调用 `genomepy.Genome` 的 `get_random_sequences`。调用前保存 Python 全局随机状态，设置 `seed`，结束后恢复，不影响调用者的随机数。N 比例超过 `max_n` 的区间不会被选中。
- **API：** `dnakit.databases.sample_random_regions(fasta_path, *, n=10, length=200, chroms=None, max_n=0.1, seed=0, verbose=True)`。
- **输入：** FASTA 路径；`n`、`length`、可选的 `chroms`、`max_n`（0 到 1）和 `seed`。
- **输出：** `regions`（`chrom`、`start`、`end`，0-based 半开区间）、`region_count`、`seed` 和 `provenance`。
- **示例代码：**

```python
import random
from pathlib import Path

from dnakit.databases import sample_random_regions

rng = random.Random(5)
chr1 = "".join(rng.choice("ACGT") for _ in range(3000))
Path("random.fa").write_text(f">chr1\n{chr1}\n>scaf_9\nNNNNACGTAC\n", encoding="utf-8")
result = sample_random_regions("random.fa", n=3, length=20, seed=7, verbose=False)
print(result["region_count"], result["seed"])
for region in result["regions"]:
    print(region["chrom"], region["start"], region["end"])
```

- **示例输出：**

```text
3 7
chr1 215 235
chr1 1596 1616
chr1 1089 1109
```

## 23) `dnakit.databases.fetch_ncbi_assembly_metadata()`

- **作用：** 按 GCF 或 GCA 登录号查询 NCBI Datasets 装配体元数据。
- **处理方法：** 调用 `ncbi.datasets.metadata.genome.get_assembly_metadata_by_asm_accessions`，需要访问网络。上游打印的警告和错误收集到 `messages`；超过 `max_records` 时报 `ConfigurationError`。
- **API：** `dnakit.databases.fetch_ncbi_assembly_metadata(accessions, *, max_records=100, verbose=True)`。
- **输入：** 登录号列表（1 到 100 个，形如 `GCF_000001405.40`）；`max_records`。
- **输出：** `records`、`record_count`、`messages`、`requested` 和 `provenance`。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import fetch_ncbi_assembly_metadata

result = fetch_ncbi_assembly_metadata(["GCF_000001405.40"], max_records=5)
print(result["record_count"])
```

- **示例输出（未本地运行）：**

```text
未本地运行：该函数需要访问 NCBI Datasets 网络服务；测试使用模拟客户端。
```

## 24) `dnakit.databases.fetch_ncbi_assemblies_by_taxon()`

- **作用：** 按分类 ID 或物种名列出 NCBI Datasets 装配体元数据。
- **处理方法：** 调用 `ncbi.datasets.metadata.genome.get_assembly_metadata_by_taxon`，需要访问网络。输出结构与 `fetch_ncbi_assembly_metadata` 相同。
- **API：** `dnakit.databases.fetch_ncbi_assemblies_by_taxon(taxon, *, max_records=100, verbose=True)`。
- **输入：** 整数分类 ID（如 `9606`）或物种名（如 `Homo sapiens`）；`max_records`。
- **输出：** `records`、`record_count`、`messages`、`requested` 和 `provenance`。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import fetch_ncbi_assemblies_by_taxon

result = fetch_ncbi_assemblies_by_taxon(9606, max_records=5)
print(result["record_count"])
```

- **示例输出（未本地运行）：**

```text
未本地运行：该函数需要访问 NCBI Datasets 网络服务；测试使用模拟客户端。
```

## 25) `dnakit.databases.fetch_ncbi_assemblies_by_bioproject()`

- **作用：** 按 BioProject 登录号列出 NCBI Datasets 装配体元数据。
- **处理方法：** 调用 `ncbi.datasets.metadata.genome.get_assembly_metadata_by_bioproject_accessions`，需要访问网络。
- **API：** `dnakit.databases.fetch_ncbi_assemblies_by_bioproject(bioproject_accessions, *, max_records=100, verbose=True)`。
- **输入：** 以 `PRJ` 开头的 BioProject 登录号列表（如 `PRJNA168`）；`max_records`。
- **输出：** `records`、`record_count`、`messages`、`requested` 和 `provenance`。
- **示例代码（未本地运行）：**

```python
from dnakit.databases import fetch_ncbi_assemblies_by_bioproject

result = fetch_ncbi_assemblies_by_bioproject(["PRJNA168"], max_records=5)
print(result["record_count"])
```

- **示例输出（未本地运行）：**

```text
未本地运行：该函数需要访问 NCBI Datasets 网络服务；测试使用模拟客户端。
```

## 26) `dnakit.databases.select_ncbi_gene_fields()`

- **作用：** 从一条 NCBI 基因记录中取出 `gene_id` 和 `symbol`。
- **处理方法：** 用 openapi 的 `V1GeneDescriptor` 构造记录，再调用 `ncbi.datasets.metadata.gene.gene_values_by_fields`。只支持这两个字段，其他字段名报错。不访问网络。
- **API：** `dnakit.databases.select_ncbi_gene_fields(gene, fields, *, verbose=True)`。
- **输入：** 含 `gene_id` 与 `symbol` 的映射；`fields` 是 `gene_id`、`symbol` 的非空子集。
- **输出：** `fields`、`values` 和 `provenance`。
- **示例代码：**

```python
from dnakit.databases import select_ncbi_gene_fields

gene = {"gene_id": "7157", "symbol": "TP53", "description": "tumor protein p53"}
print(select_ncbi_gene_fields(gene, ["gene_id", "symbol"], verbose=False)["values"])
```

- **示例输出：**

```text
{'gene_id': '7157', 'symbol': 'TP53'}
```

## 27) `dnakit.databases.select_ncbi_assembly_fields()`

- **作用：** 从一条装配体记录中取出顶层字段。
- **处理方法：** 用 openapi 的 `V1AssemblyDatasetDescriptor` 构造记录，再调用 `ncbi.datasets.metadata.genome.assembly_values_by_fields`。`fields` 为空时返回全部顶层字段。不访问网络。
- **API：** `dnakit.databases.select_ncbi_assembly_fields(assembly, fields=(), *, verbose=True)`。
- **输入：** 装配体记录映射（字符串键）；`fields` 为顶层字段名列表。
- **输出：** `fields`、`values`、`missing_fields`（请求但记录中不存在的字段）和 `provenance`。
- **示例代码：**

```python
from dnakit.databases import select_ncbi_assembly_fields

assembly = {"accession": "GCF_000001405.40", "organism": {"tax_id": 9606}}
print(select_ncbi_assembly_fields(assembly, ["accession", "organism"], verbose=False)["values"])
```

- **示例输出：**

```text
{'accession': 'GCF_000001405.40', 'organism': {'tax_id': 9606}}
```

## 验证方式

- 测试文件：`others/tests/unit/test_full_G8.py`，共 36 个测试（含参数化的关键字 `verbose` 签名检查）。
- 每个封装函数至少有一个测试直接调用真实上游包，并与手工计算的结果比较（FASTA 长度、区间合并、GTF 列值、BED gap 坐标、Ensembl 下载路径），或与直接调用上游的结果比较。
- NCBI Datasets 的装配体函数使用模拟 API 客户端运行真实的上游分页生成器；测试不访问网络。
- 每个模块都有缺包测试：把对应包的 `sys.modules` 项置为 `None`，检查抛出 `BackendUnavailableError` 且消息包含安装额外项名称。

## 当前范围

- 已封装 27 个函数：pyensembl 15 个（FASTA 读取、区间合并与最近位点、标识符与别名规范化、版本匹配、参考基因组名与 Ensembl 下载地址）；genomepy 7 个（注释读写与过滤、gap 与 sizes 文件、FASTA 过滤、随机区间采样）；NCBI Datasets 5 个（3 个装配体元数据查询和 2 个字段选择）。
- 未封装：pyensembl 的 `Genome` 与 `EnsemblRelease` 类、缓存与日期版本函数，以及 CLI 函数（它们依赖下载或本地缓存目录）；genomepy 的 `Genome` 方法族（`track2fasta` 已由 `write_bed_sequences` 封装）、安装与 provider 查询（`install_genome` 已由 `install_genome_assembly` 封装）、需要 UCSC 工具的函数、配置与插件函数；NCBI Datasets 的 openapi 模型工具、gRPC 注册和 `print_*` 输出函数。
- `dnakit.databases` 之外的导出尚未接入，因此示例直接从模块文件导入。
