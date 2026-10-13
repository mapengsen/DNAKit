# pysam、cyvcf2 与 vcfpy 的变异与区域辅助函数

这一组函数是 `dnakit.variants` 中基于可选依赖的适配器，覆盖 [pysam](https://github.com/pysam-developers/pysam)、[cyvcf2](https://github.com/brentp/cyvcf2) 和 [vcfpy](https://github.com/dahlo/vcfpy) 中可以离线运行的 VCF、基因组区间、参考序列和群体基因型功能。pysam 的 bcftools 功能使用它随包附带的 htslib 与 bcftools，无需另外安装 bcftools 或 tabix。这些函数只处理显式给出的文件和序列，不访问网络，不调用变异检测程序，也不查询数据库。依赖都是可选的，只在调用时导入；缺少包时会报错，并给出安装命令。

## 示例数据

示例按顺序运行，后面的示例使用前面生成的文件。示例从 `dnakit.variants` 导入。下面的 `sample.vcf` 有 4 个位点和 3 个样本，带有 AD 与 DP 字段；`ref.fa` 是 50 bp 的参考序列；`annotation.tsv` 是两行注释表。这些文件只用于演示接口，统计值不代表真实群体。

**示例文件：** `sample.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=50>
##FILTER=<ID=q10,Description="Quality below 10">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	3	rs1	G	T	50	PASS	DP=20;AF=0.5	GT:AD:DP	0/1:10,10:20	1/1:0,12:12	0/0:15,0:15
chr1	10	.	C	CA	40	PASS	DP=10	GT:AD:DP	0/1:5,5:10	0/0:10,0:10	0/1:4,6:10
chr1	22	rs3	C	G,A	60	q10	DP=30;AF=0.3,0.1	GT:AD:DP	1|2:0,10,5:15	0/0:20,0,0:20	0/2:10,0,10:20
chr1	35	rs4	G	C	70	PASS	DP=12;AF=0.4	GT:AD:DP	0/1:6,6:12	0/1:5,5:10	./.:0,0:0
```

**示例文件：** `ref.fa`

```text
>chr1
ACGTACGTAC
GTACGTTTGA
CCAAGGTTAC
GATCGATCGA
ACGTACGTAC
```

**示例文件：** `annotation.tsv`

```text
chr1	3	G	T	common
chr1	22	C	G	rare
```

## 1) `dnakit.variants.index_tabix_file()`

- **作用：** 对 VCF、BED 或 GFF 文件建立 tabix 索引，之后可以按区域读取。
- **处理方法：** 调用 `pysam.tabix_index(..., keep_original=True)`。未压缩的输入会先压缩为同目录下的 `<文件名>.gz`，原文件保留。`preset` 只接受 `vcf`、`bed` 或 `gff`。已有索引时拒绝覆盖，除非 `force=True`。
- **API：** `dnakit.variants.index_tabix_file(source, *, preset="vcf", force=False)`。
- **输入：** 文本或 bgzip 压缩的 VCF、BED 或 GFF 路径。
- **输出：** 字典，包含 `compressed`（`.gz` 路径）、`index`（`.tbi` 路径）、`preset` 和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import index_tabix_file

result = index_tabix_file("sample.vcf", force=True, verbose=False)
print(Path(result["compressed"]).name, Path(result["index"]).name)
print(result["preset"])
```

- **示例输出：**

```text
sample.vcf.gz sample.vcf.gz.tbi
vcf
```

## 2) `dnakit.variants.compress_bgzip_file()`

- **作用：** 把文本文件写成 BGZF 压缩文件，供 tabix 与 htslib 读取。
- **处理方法：** 调用 `pysam.tabix_compress`。输出名必须以 `.gz` 结尾，已有文件不会被覆盖，除非 `force=True`。
- **API：** `dnakit.variants.compress_bgzip_file(source, output_path, *, force=False)`。
- **输入：** 源文本文件和新的 `.gz` 输出路径。
- **输出：** 字典，包含 `output`、两个文件大小和 `provenance`。
- **示例代码：**

```python
from dnakit.variants import compress_bgzip_file

result = compress_bgzip_file("annotation.tsv", "annotation.tsv.gz", force=True, verbose=False)
print(result["operation"], result["uncompressed_bytes"])
```

- **示例输出：**

```text
compress_bgzip_file 35
```

## 3) `dnakit.variants.parse_genomic_region()`

- **作用：** 把 `chr1:起点-终点`（1-based，闭区间）转换为 pysam 使用的 0-based 半开区间。
- **处理方法：** 先严格检查格式，再调用 `pysam.libcutils.parse_region`。只给染色体名时返回整条染色体，`end` 为 `None`。不接受 `chr1:1000` 这种没有终点的写法。
- **API：** `dnakit.variants.parse_genomic_region(region)`。
- **输入：** 区域字符串，可以在数字中使用逗号，例如 `chr1:1,001-1,020`。
- **输出：** 字典，包含 `reference`、`start`（0-based）、`end`（不含）和 `length`。
- **示例代码：**

```python
from dnakit.variants import parse_genomic_region

result = parse_genomic_region("chr1:11-20", verbose=False)
print(result["reference"], result["start"], result["end"], result["length"])
```

- **示例输出：**

```text
chr1 10 20 10
```

## 4) `dnakit.variants.parse_genomic_regions()`

- **作用：** 解析逗号分隔的区域列表，例如 `chr1:10-20,chr2`。
- **处理方法：** 每一段都按严格格式检查，然后调用 `pysam.libcvcf.parse_regions`。上游会把 `chr1:10-20-30` 悄悄读成 `chr1:10-20`，本函数先拒绝这类输入。
- **API：** `dnakit.variants.parse_genomic_regions(regions)`。
- **输入：** 逗号分隔的区域字符串，最多 10000 段。
- **输出：** 字典，`regions` 中每项包含 `reference`、`start`（0-based）、`end`（不含，整条染色体为 `None`），以及 `count`。
- **示例代码：**

```python
from dnakit.variants import parse_genomic_regions

result = parse_genomic_regions("chr1:10-20,chr2", verbose=False)
for item in result["regions"]:
    print(item["reference"], item["start"], item["end"])
```

- **示例输出：**

```text
chr1 9 20
chr2 0 None
```

## 5) `dnakit.variants.fetch_padded_sequence()`

- **作用：** 读取参考序列的一个窗口，窗口超出染色体的部分补 N。
- **处理方法：** 用 `pysam.FastaFile` 打开 FASTA，再调用 `pysam.libcvcf.get_sequence`。窗口为 0-based 半开区间，`start` 可以为负数。输入必须已有 `.fai` 索引。
- **API：** `dnakit.variants.fetch_padded_sequence(fasta, reference, start, end)`。
- **输入：** FASTA 路径（已建立 `.fai`）、染色体名、`start` 与 `end`。窗口最长 10000000 bp。
- **输出：** 字典，包含 `sequence`（大写）、`reference_length`、`left_padding` 和 `right_padding`。
- **示例代码：**

```python
import pysam

from dnakit.variants import fetch_padded_sequence

pysam.faidx("ref.fa")
right = fetch_padded_sequence("ref.fa", "chr1", 45, 60, verbose=False)
left = fetch_padded_sequence("ref.fa", "chr1", -3, 4, verbose=False)
print(right["sequence"], right["right_padding"])
print(left["sequence"], left["left_padding"])
```

- **示例输出：**

```text
CGTACNNNNNNNNNN 10
NNNACGT 3
```

## 6) `dnakit.variants.query_vcf_fields()`

- **作用：** 用 `bcftools query` 把 VCF 的指定列输出为文本行。
- **处理方法：** `fields` 可以是 `CHROM`、`POS`、`ID`、`REF`、`ALT`、`QUAL`、`FILTER` 或 `INFO/<标签>`。`genotypes=True` 时，每个样本增加一列 `GT:<样本名>`。数值保持文本，缺失值保留为 `.`。给出 `region` 时需要索引。
- **API：** `dnakit.variants.query_vcf_fields(source, fields, *, region=None, genotypes=False, max_records=100000)`。
- **输入：** VCF 路径、字段列表，可选的区域与 `max_records`（1 到 1000000）。
- **输出：** 字典，包含 `columns`、`rows`（字符串列表）、`record_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.variants import query_vcf_fields

result = query_vcf_fields("sample.vcf", ["POS", "REF", "ALT", "FILTER", "INFO/DP"], verbose=False)
print(result["columns"])
for row in result["rows"]:
    print(row)
```

- **示例输出：**

```text
['POS', 'REF', 'ALT', 'FILTER', 'INFO/DP']
['3', 'G', 'T', 'PASS', '20']
['10', 'C', 'CA', 'PASS', '10']
['22', 'C', 'G,A', 'q10', '30']
['35', 'G', 'C', 'PASS', '12']
```

## 7) `dnakit.variants.select_vcf_records()`

- **作用：** 按样本、区域和 PASS 过滤写出一个新的未压缩 VCF。
- **处理方法：** 调用 `bcftools view -s -r -f PASS`。样本顺序与给出的顺序一致。区域需要索引。输出文件默认不覆盖。
- **API：** `dnakit.variants.select_vcf_records(source, output_path, *, samples=None, region=None, pass_only=False, force=False)`。
- **输入：** VCF 路径、新输出路径，以及可选的样本名列表、区域和 `pass_only`。
- **输出：** 字典，包含 `output`、保留的记录数 `records`、`samples` 和 `provenance`。
- **示例代码：**

```python
from dnakit.variants import select_vcf_records

result = select_vcf_records(
    "sample.vcf.gz",
    "subset.vcf",
    samples=["S2", "S1"],
    region="chr1:9-22",
    pass_only=True,
    force=True,
    verbose=False,
)
print(result["records"], result["samples"])
```

- **示例输出：**

```text
1 ['S2', 'S1']
```

## 8) `dnakit.variants.concatenate_vcf_files()`

- **作用：** 按给出的顺序拼接多个 VCF，常用于合并同一样本集的不同区段。
- **处理方法：** 调用 `bcftools concat`。每个输入都必须是已建立索引的 bgzip 文件。输出为未压缩 VCF。
- **API：** `dnakit.variants.concatenate_vcf_files(sources, output_path, *, force=False)`。
- **输入：** 至少两个已索引 VCF 路径和新的输出路径。
- **输出：** 字典，包含 `output`、输入文件数 `inputs`、记录数 `records` 和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

import pysam

from dnakit.variants import concatenate_vcf_files

lines = Path("sample.vcf").read_text().splitlines(keepends=True)
header = [line for line in lines if line.startswith("#")]
body = [line for line in lines if not line.startswith("#")]
Path("part1.vcf").write_text("".join(header + body[:2]))
Path("part2.vcf").write_text("".join(header + body[2:]))
for name in ("part1.vcf", "part2.vcf"):
    pysam.tabix_index(name, preset="vcf", force=True, keep_original=True)
result = concatenate_vcf_files(
    ["part1.vcf.gz", "part2.vcf.gz"], "joined.vcf", force=True, verbose=False
)
print(result["inputs"], result["records"])
```

- **示例输出：**

```text
2 4
```

## 9) `dnakit.variants.intersect_vcf_sites()`

- **作用：** 统计两个 VCF 中各自独有的位点数和共享位点数。
- **处理方法：** 调用 `bcftools isec -p`。位点按染色体、位置和等位基因匹配，每个位点只计一次。`output_dir` 必须不存在，bcftools 会在其中写入 VCF 与 `sites.txt`。
- **API：** `dnakit.variants.intersect_vcf_sites(first, second, output_dir)`。
- **输入：** 两个已索引 VCF 路径和一个新目录路径。
- **输出：** 字典，包含 `private_first`、`private_second`、`shared`、`directory`、`files` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.variants import intersect_vcf_sites

extra = Path("sample.vcf").read_text() + "chr1\t40\t.\tA\tT\t30\tPASS\t.\tGT\t0/1\t0/1\t0/1\n"
Path("more.vcf").write_text(extra)
pysam.tabix_index("more.vcf", preset="vcf", force=True, keep_original=True)
target = Path(tempfile.mkdtemp(), "isec")
result = intersect_vcf_sites("sample.vcf.gz", "more.vcf.gz", target, verbose=False)
print(result["private_first"], result["private_second"], result["shared"])
```

- **示例输出：**

```text
0 1 4
```

## 10) `dnakit.variants.sort_vcf_records()`

- **作用：** 按染色体顺序和位置重新排列 VCF 记录。
- **处理方法：** 调用 `bcftools sort`，输出未压缩 VCF。染色体顺序沿用输入头部的 contig 顺序。
- **API：** `dnakit.variants.sort_vcf_records(source, output_path, *, force=False)`。
- **输入：** VCF 路径（可以未压缩）和新的输出路径。
- **输出：** 字典，包含 `output`、记录数 `records` 和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import sort_vcf_records

lines = Path("sample.vcf").read_text().splitlines(keepends=True)
header = [line for line in lines if line.startswith("#")]
body = [line for line in lines if not line.startswith("#")]
Path("shuffled.vcf").write_text("".join(header + body[::-1]))
result = sort_vcf_records("shuffled.vcf", "sorted.vcf", force=True, verbose=False)
positions = [line.split("\t")[1] for line in Path("sorted.vcf").read_text().splitlines()
             if not line.startswith("#")]
print(result["records"], positions)
```

- **示例输出：**

```text
4 ['3', '10', '22', '35']
```

## 11) `dnakit.variants.annotate_vcf_from_tabix()`

- **作用：** 把注释表中的 INFO 标签复制到 VCF 中位点相同的记录上。
- **处理方法：** 调用 `bcftools annotate -a -c -h`。注释表必须是 bgzip 压缩并已建立 tabix 索引的文件。`columns` 按注释表列的顺序给出，必须包含 `CHROM` 和 `POS`，并至少含一个 `INFO/<标签>`。每个新标签都需要在 `header_lines` 中有对应的 `##INFO` 行。
- **API：** `dnakit.variants.annotate_vcf_from_tabix(source, annotations, columns, output_path, *, header_lines=(), force=False)`。
- **输入：** VCF 路径、注释表路径、列名列表、新输出路径，以及可选的头部行。
- **输出：** 字典，包含 `output`、记录数 `records`、每个标签的注释数 `annotated` 和 `provenance`。
- **示例代码：**

```python
import pysam

from dnakit.variants import annotate_vcf_from_tabix

pysam.tabix_compress("annotation.tsv", "annotation.tsv.gz", force=True)
pysam.tabix_index("annotation.tsv.gz", seq_col=0, start_col=1, end_col=1, force=True)
result = annotate_vcf_from_tabix(
    "sample.vcf",
    "annotation.tsv.gz",
    ["CHROM", "POS", "REF", "ALT", "INFO/NOTE"],
    "annotated.vcf",
    header_lines=['##INFO=<ID=NOTE,Number=1,Type=String,Description="note">'],
    force=True,
    verbose=False,
)
print(result["annotated"], result["records"])
```

- **示例输出：**

```text
{'NOTE': 2} 4
```

## 12) `dnakit.variants.build_vcf_consensus()`

- **作用：** 把一个样本某一条单倍型的变异应用到参考序列上，生成一条共识序列。
- **处理方法：** 调用 `bcftools consensus -f -s -H`。VCF 必须是 bgzip 压缩并已建立索引；FASTA 需要 `.fai` 索引。杂合位点取 GT 中该单倍型对应的等位基因。
- **API：** `dnakit.variants.build_vcf_consensus(source, fasta, sample, output_path, *, haplotype=1, force=False)`。
- **输入：** 已索引的 bgzip VCF、FASTA、样本名、单倍型（1 或 2）和新的输出路径。
- **输出：** 字典，包含 `output`、`sample`、`haplotype`、每条序列的长度 `sequences` 和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import build_vcf_consensus

result = build_vcf_consensus(
    "sample.vcf.gz", "ref.fa", "S2", "s2_hap1.fa", haplotype=1, force=True, verbose=False
)
print(result["sample"], result["haplotype"], result["sequences"])
print(Path("s2_hap1.fa").read_text().splitlines()[1][:20])
```

- **示例输出：**

```text
S2 1 {'chr1': 50}
ACTTACGTACGTACGTTTGA
```

## 13) `dnakit.variants.summarize_sample_heterozygosity()`

- **作用：** 在指定位点上计算每个样本的杂合比例，以及杂合位点等位基因比例的分位数。
- **处理方法：** 调用 `cyvcf2.par_het`，位点由调用者明确给出，不使用上游的稀疏抽样规则。深度来自 FORMAT 的 AD；只有深度大于 `min_depth` 的杂合位点进入分位数。上游会跳过等位基因频率低于 0.01 或调用率低于 0.8 的位点，结果中的 `sites_used` 列出实际使用的位点。
- **API：** `dnakit.variants.summarize_sample_heterozygosity(source, samples, sites, *, ncpus=1, min_depth=8, percentiles=(10, 90))`。
- **输入：** 已索引的 bgzip VCF、样本名、位点列表（`(染色体, 位置)` 或 `(染色体, 位置, REF, ALT)`）。
- **输出：** 字典，`samples` 中每个样本包含 `het_count`、`het_ratio`、`p10`、`p90`、`mean_depth` 等值；另有 `sites_used` 和 `warnings`。
- **示例代码：**

```python
from dnakit.variants import summarize_sample_heterozygosity

result = summarize_sample_heterozygosity(
    "sample.vcf.gz", ["S1", "S2", "S3"], [("chr1", 3), ("chr1", 10)], verbose=False
)
for name, values in result["samples"].items():
    print(name, values["het_count"], values["het_ratio"], values["p10"])
```

- **示例输出：**

```text
S1 2 1.0 0.5
S2 0 0.0 0.0
S3 1 0.5 0.6
```

## 14) `dnakit.variants.summarize_pairwise_relatedness()`

- **作用：** 在指定位点上计算样本两两之间的 IBS 计数和亲缘度系数。
- **处理方法：** 调用 `cyvcf2.par_relatedness`，与 `VCF.site_relatedness` 的结果一致。每对样本只出现一次，按样本名排序。没有杂合位点的样本会产生警告，此时系数不可靠。
- **API：** `dnakit.variants.summarize_pairwise_relatedness(source, samples, sites, *, ncpus=1, min_depth=5)`。
- **输入：** 已索引的 bgzip VCF、样本名和位点列表。
- **输出：** 字典，`pairs` 中每项包含 `sample_a`、`sample_b`、`relatedness`、`ibs0`、`ibs2`、`informative_sites` 等；另有 `warnings`。
- **示例代码：**

```python
from dnakit.variants import summarize_pairwise_relatedness

result = summarize_pairwise_relatedness(
    "sample.vcf.gz",
    ["S1", "S2", "S3"],
    [("chr1", 3, "G", "T"), ("chr1", 10, "C", "CA")],
    verbose=False,
)
for pair in result["pairs"]:
    print(pair["sample_a"], pair["sample_b"], pair["relatedness"], pair["ibs2"])
```

- **示例输出：**

```text
S1 S2 0.0 0
S1 S3 1.0 1
S2 S3 -2.0 0
```

## 15) `dnakit.variants.pearson_genotype_correlation()`

- **作用：** 计算两个基因型剂量向量（0、1、2）之间的 Pearson 相关系数。
- **处理方法：** 调用 `cyvcf2.r_unphased`。两侧任一缺失的位置被跳过。可用位置少于两个，或某一侧在可用位置上恒定时，返回 `None`。上游的频率参数不被使用，包装函数传入 0.0。
- **API：** `dnakit.variants.pearson_genotype_correlation(first, second)`。
- **输入：** 两个等长列表，元素为 0、1、2 或 `None`（缺失）。
- **输出：** 字典，包含 `correlation`、`usable_pairs`、`missing_pairs` 和 `provenance`。
- **示例代码：**

```python
from dnakit.variants import pearson_genotype_correlation

print(pearson_genotype_correlation([0, 1, 2, 1], [0, 2, 1, 1], verbose=False)["correlation"])
print(pearson_genotype_correlation([0, 1, None, 1], [0, 2, 1, 1], verbose=False)["correlation"])
```

- **示例输出：**

```text
0.5
0.8660253882408142
```

## 16) `dnakit.variants.classify_vcf_alt()`

- **作用：** 判断一个 ALT 等位基因的类别，例如替换、符号等位基因或断点。
- **处理方法：** 调用 `vcfpy.parser.process_alt`。替换类等位基因还会调用 `process_sub`，给出 SNV、MNV、INS、DEL 或 INDEL。断点以 `breakend` 字段返回。
- **API：** `dnakit.variants.classify_vcf_alt(ref, alt)`。
- **输入：** REF 碱基（A、C、G、T、N）和 ALT 字符串。ALT 不能是 `.`。
- **输出：** 字典，包含 `kind`、`subtype`、`value`、`breakend` 和 `provenance`。
- **示例代码：**

```python
from dnakit.variants import classify_vcf_alt

pairs = [("G", "T"), ("CAT", "C"), ("C", "CA"), ("C", "<DEL>"), ("C", "G]17:198982]")]
for ref, alt in pairs:
    result = classify_vcf_alt(ref, alt, verbose=False)
    print(ref, alt, result["kind"], result["subtype"], result["value"])
```

- **示例输出：**

```text
G T substitution SNV T
CAT C substitution DEL C
C CA substitution INS CA
C <DEL> symbolic None DEL
C G]17:198982] breakend None None
```

## 17) `dnakit.variants.classify_vcf_substitution()`

- **作用：** 根据 REF 与 ALT 的长度和碱基，判断替换的类型。
- **处理方法：** 调用 `vcfpy.parser.process_sub`。长度相等时为 SNV 或 MNV，ALT 更长为 INS 或 INDEL，ALT 更短为 DEL 或 INDEL。`process_sub_grow` 与 `process_sub_shrink` 由它统一调用。
- **API：** `dnakit.variants.classify_vcf_substitution(ref, alt)`。
- **输入：** REF 与 ALT 两段非空碱基序列。
- **输出：** 字典，包含 `subtype`、`value` 和 `provenance`。
- **示例代码：**

```python
from dnakit.variants import classify_vcf_substitution

for ref, alt in [("AT", "A"), ("AC", "G"), ("ACG", "TTT")]:
    print(ref, alt, classify_vcf_substitution(ref, alt, verbose=False)["subtype"])
```

- **示例输出：**

```text
AT A DEL
AC G INDEL
ACG TTT MNV
```

## 18) `dnakit.variants.parse_vcf_breakend()`

- **作用：** 解析方括号形式的断点 ALT，例如 `G]17:198982]`，得到配对位置和方向。
- **处理方法：** 调用 `vcfpy.parser.parse_breakend`。非括号形式会被拒绝。`orientation` 与 `mate_orientation` 为 `+` 或 `-`；尖括号中的 contig 会把 `within_main_assembly` 设为假。
- **API：** `dnakit.variants.parse_vcf_breakend(alt)`。
- **输入：** 断点 ALT 字符串。
- **输出：** 字典，包含 `mate_chrom`、`mate_pos`、`orientation`、`mate_orientation`、`sequence`、`within_main_assembly` 和 `provenance`。
- **示例代码：**

```python
from dnakit.variants import parse_vcf_breakend

result = parse_vcf_breakend("]13:123456]T", verbose=False)
print(result["mate_chrom"], result["mate_pos"], result["orientation"],
      result["mate_orientation"], result["sequence"])
```

- **示例输出：**

```text
13 123456 + - T
```

## 19) `dnakit.variants.parse_vcf_field_value()`

- **作用：** 按声明的 Number 与 Type，把一个 INFO 或 FORMAT 值转换为 Python 对象。
- **处理方法：** 调用 `vcfpy.parser.parse_field_value`。`number` 为整数，或 `A`、`R`、`G`、`.`。Number 为 1 时返回单个值，其他情况返回列表；`.` 表示缺失。无法转换的值会报错，而不是保留为字符串。
- **API：** `dnakit.variants.parse_vcf_field_value(text, *, value_type, number, field_id='')`。
- **输入：** 字段文本、类型（`Integer`、`Float`、`Flag`、`Character`、`String`）和 Number。
- **输出：** 字典，包含 `value`、`value_type`、`number` 和 `provenance`。
- **示例代码：**

```python
from dnakit.variants import parse_vcf_field_value

print(parse_vcf_field_value("0.5,.", value_type="Float", number="A", verbose=False)["value"])
print(parse_vcf_field_value("12", value_type="Integer", number=1, verbose=False)["value"])
```

- **示例输出：**

```text
[0.5, None]
12
```

## 20) `dnakit.variants.format_vcf_field_value()`

- **作用：** 把一个 INFO 或 FORMAT 值格式化为 VCF 文本，是解析函数的逆过程。
- **处理方法：** 调用 `vcfpy.writer.format_value`。类型不匹配会被拒绝。`Flag` 类型不在这里格式化，应直接写出标签名。`None` 写为 `.`。
- **API：** `dnakit.variants.format_vcf_field_value(value, *, value_type, number, section='INFO', field_id='')`。
- **输入：** 值、类型、Number、`section`（`INFO` 或 `FORMAT`）和可选的字段 ID。
- **输出：** 字典，包含 `text`、`value_type`、`number`、`section` 和 `provenance`。
- **示例代码：**

```python
from dnakit.variants import format_vcf_field_value

print(format_vcf_field_value([0.5, 0.25], value_type="Float", number="A", verbose=False)["text"])
print(format_vcf_field_value(None, value_type="Integer", number=1, verbose=False)["text"])
```

- **示例输出：**

```text
0.5,0.25
.
```

## 21) `dnakit.variants.remove_vcf_header_lines()`

- **作用：** 读取 VCF 头部，并删去指定的 `(键, ID)` 行。
- **处理方法：** 用 `vcfpy.Reader` 读取头部，再调用 `vcfpy.header_without_lines`。没有 ID 的行可以把整个值放在 ID 的位置。`#CHROM` 行不在结果中。
- **API：** `dnakit.variants.remove_vcf_header_lines(source, remove)`。
- **输入：** VCF 路径（可以是压缩文件）和 `(键, ID)` 列表。
- **输出：** 字典，包含 `header_lines`（`##` 开头的行）、`samples`、`removed_lines` 和 `provenance`。
- **示例代码：**

```python
from dnakit.variants import remove_vcf_header_lines

result = remove_vcf_header_lines("sample.vcf", [("FILTER", "q10")], verbose=False)
print(result["removed_lines"], len(result["header_lines"]), result["samples"])
```

- **示例输出：**

```text
1 7 ['S1', 'S2', 'S3']
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-reads]"        # pysam (bcftools and htslib bundled)
python -m pip install "dnakit[bio-cyvcf2]"       # cyvcf2
python -m pip install "dnakit[bio-vcfpy]"        # vcfpy
python -m pip install "dnakit[bio-analysis]"     # numpy, used by pearson_genotype_correlation
```

项目页面：[pysam](https://github.com/pysam-developers/pysam)、[cyvcf2](https://github.com/brentp/cyvcf2)、[vcfpy](https://github.com/dahlo/vcfpy)。

## 验证方式

- 测试文件是 `others/tests/unit/test_full_G5.py`，共 24 个测试，运行时间约 1 秒。
- pysam 的 tabix 索引、bgzip、区域解析和 FASTA 窗口结果按手算值检查；bcftools 的 query、view、concat、isec、sort、annotate 结果按手算值检查，consensus 同时与 `pysam.bcftools.consensus` 直接调用的结果对照。
- cyvcf2 的杂合比例与 `VCF.het_check` 直接对照，亲缘度与 `VCF.site_relatedness` 直接对照；`pearson_genotype_correlation` 的取值按手算值检查。
- vcfpy 的 ALT 与替换类型对照 `vcfpy.parser.process_alt` 和 `process_sub`；断点解析、字段值解析与格式化按手算值检查；header 删除结果与 `vcfpy.header_without_lines` 对照。
- 每个模块都检查缺少依赖时会抛出 `BackendUnavailableError`，并带有安装额外项的提示。非法输入（例如未压缩或未索引的 VCF、未知样本、过长窗口、不支持的字段）会被拒绝。

## 当前范围

- 本批包装的上游操作：pysam 的 tabix 与 bgzip、区域解析、FASTA 窗口读取，以及 bcftools 的 query、view、concat、isec、sort、annotate 和 consensus；cyvcf2 的 `par_het`、`par_relatedness` 和 `r_unphased`；vcfpy 的 ALT 与替换类型分类、断点解析、字段值解析与格式化，以及 header 删除。
- 未包装的 pysam 功能：`bcftools csq` 需要 GFF3 基因模型以及杂合调用的相位策略；`roh`、`gtcheck`、`merge`、`cnv`、`convert`、`filter` 等不在本批；`pysam.Pileup` 是过时的 samtools `pileup -c` 格式接口，`tabix_iterator` 需要 pysam 解析器对象。BAM 与 SAM 相关的 samtools 命令不属于变异区域，未包装。
- 与已有函数重复、未包装：`bcftools norm` 与 `bcftools stats` 由 `dnakit.databases.vcf_normalize` 和 `vcf_summary_stats` 覆盖；`bcftools call` 由 `dnakit.ngs.call_variants` 覆盖；`samtools faidx` 由 `dnakit.bioinformatics.reads.index_reference` 覆盖。
- cyvcf2 的 `cli` 命令行打印函数和 bcftools 的 `plugin` 命令不包装。
