# 72 变异：varcode、hgvs 与 scikit-allel 离线函数

本文档记录 G6a 批次包装的 varcode、hgvs 与 scikit-allel 函数。所有函数只读取本地输入，不需要参考基因组、pyensembl 缓存、seqrepo、UTA、PostgreSQL 或网络。

包装模块为 `dnakit.variants`、`dnakit.variants` 与 `dnakit.variants`，示例直接从这些模块导入。

依赖为可选包：`bio-varcode`、`bio-hgvs`、`bio-scikit-allel` 与 `bio-analysis`（numpy）。缺少依赖时抛出 `BackendUnavailableError` 并给出安装提示。

## 概览

| 包 | 上游公开函数条目 | 本批包装 | 已由 dnakit 包装 | 族变体 | 与现有 dnakit 函数重复 |
| --- | ---: | ---: | ---: | ---: | ---: |
| varcode | 333 | 18 | 1 | 2 | 6 |
| hgvs | 31 | 6 | 0 | 2 | 1 |
| scikit-allel | 634 | 20 | 11 | 24 | 6 |

其余条目在账本中记录原因：需要 pyensembl 或 Ensembl 数据、需要 seqrepo 或 UTA、已有 dnakit 函数、已弃用、命令行入口、或超出包装上限。

## 1) 枚举隐蔽外显子候选

- **作用：** 在一段本地序列中找出可能组成新内含子外显子的受体-供体位点对。
- **计算方法：** 用 varcode 的 9 bp 供体与 4 bp 受体共识模型逐窗口打分，保留分数不低于阈值的位点，再按内含子长度范围配对；结果按分数和、起点、终点排序。
- **API：** `dnakit.variants.enumerate_cryptic_exon_candidates(contig: str, sequence: str, sequence_start: int = 1, *, donor_threshold: float = 0.7, acceptor_threshold: float = 0.7, min_intron_length: int = 30, max_intron_length: int = 10000)`
- **输入：** contig、序列（仅含 A/C/G/T/N）、序列起点（1-based）；可设阈值与内含子长度范围。
- **输出：** 返回 `candidates` 列表，每项含 contig、interval_start、interval_end、donor_score、acceptor_score。
- **示例代码：**

```python
from dnakit.variants import enumerate_cryptic_exon_candidates

sequence = "TTTT" + "CAGG" + "C" * 30 + "AAG" + "GTAAGT" + "TTTT"
result = enumerate_cryptic_exon_candidates("chr1", sequence, 1000, min_intron_length=10)
print(result["count"], result["candidates"][0])
```

- **示例输出：**

```text
1 {'contig': 'chr1', 'interval_start': 1007, 'interval_end': 1041, 'donor_score': 1.0, 'acceptor_score': 1.0}
```

## 2) 给剪接位点窗口打分

- **作用：** 对一个 9 bp 供体窗口或 4 bp 受体窗口给出共识匹配分数。
- **计算方法：** 检查 GT（供体）或 AG（受体）锚定碱基；锚定不符时得 0 分，否则为与共识匹配的碱基比例。
- **API：** `dnakit.variants.score_splice_site(window: str, *, site: str = 'donor')`
- **输入：** 窗口序列；`site` 为 `donor`（9 bp）或 `acceptor`（4 bp）。
- **输出：** 返回 `score`（0 到 1 之间）和窗口长度。
- **示例代码：**

```python
from dnakit.variants import score_splice_site

print(score_splice_site("CAGGTAAGT", site="donor")["score"])
print(score_splice_site("CAGG", site="acceptor")["score"])
print(score_splice_site("CAGATAAGT", site="donor")["score"])
```

- **示例输出：**

```text
1.0
1.0
0.0
```

## 3) 解析符号型结构变异 ALT

- **作用：** 把 `<DEL>`、`<DUP>`、`<INV>`、`<CNV>` 或断点 ALT 解析为结构变异字段。
- **计算方法：** 按 VCF 约定读取 END、SVTYPE、MATEID 等 INFO 字段，END 不大于 POS 时报错。
- **API：** `dnakit.variants.parse_symbolic_sv_alt(contig: str, start: int, ref: str, alt: str, *, info: Mapping[str, Any] | None = None)`
- **输入：** contig、POS、REF、ALT 和可选的 INFO 映射。
- **输出：** 返回 `sv_type`、`start`、`end`、`affected_start`、`affected_end`；断点记录还含 `mate_contig` 与 `mate_start`。
- **示例代码：**

```python
from dnakit.variants import parse_symbolic_sv_alt

record = parse_symbolic_sv_alt("chr1", 100, "N", "<DEL>", info={"END": 500})
print(record["sv_type"], record["start"], record["end"], record["affected_start"], record["affected_end"])
```

- **示例输出：**

```text
DEL 100 500 101 500
```

## 4) 读取断点 ALT 保留的一侧

- **作用：** 判断 VCF 断点 ALT 保留断点的左侧还是右侧。
- **计算方法：** 匹配括号与前后缀形式；`t]p]` 等形式分别映射为 left 或 right。
- **API：** `dnakit.variants.breakend_alt_sides(alt: str)`
- **输入：** VCF 断点 ALT 字符串。
- **输出：** 返回 `valid`、`this_side` 与 `mate_side`；非断点 ALT 返回 `valid=False`。
- **示例代码：**

```python
from dnakit.variants import breakend_alt_sides

sides = breakend_alt_sides("A]chr1:500]")
print(sides["valid"], sides["this_side"], sides["mate_side"])
```

- **示例输出：**

```text
True left left
```

## 5) 读取断点 ALT 的配对位置

- **作用：** 从断点 ALT 中取出配对染色体、位置与括号方向。
- **计算方法：** 用 VCF 断点正则解析 ALT；无效 ALT 返回 `valid=False`。
- **API：** `dnakit.variants.breakend_mate_position(alt: str)`
- **输入：** VCF 断点 ALT 字符串。
- **输出：** 返回 `mate_contig`、`mate_position` 与 `mate_orientation`。
- **示例代码：**

```python
from dnakit.variants import breakend_mate_position

mate = breakend_mate_position("A]chr1:500]")
print(mate["mate_contig"], mate["mate_position"], mate["mate_orientation"])
```

- **示例输出：**

```text
chr1 500 ]]
```

## 6) 读取断点 ALT 插入的碱基

- **作用：** 给出断点 ALT 相对一个碱基 REF 锚点插入的碱基。
- **计算方法：** 前缀形式取锚点之后、后缀形式取锚点之前的碱基；锚点不是单碱基或不匹配时视为未知。
- **API：** `dnakit.variants.breakend_insert_bases(ref: str, alt: str)`
- **输入：** 单碱基 REF 与断点 ALT。
- **输出：** 返回 `known` 与 `inserted_bases`；未知时 `inserted_bases` 为 None。
- **示例代码：**

```python
from dnakit.variants import breakend_insert_bases

inserted = breakend_insert_bases("A", "AGT[chr1:500[")
print(inserted["known"], inserted["inserted_bases"])
```

- **示例输出：**

```text
True GT
```

## 7) 判断成对断点的事件类型

- **作用：** 把两条配对断点记录判定为 DEL、DUP 或 INV。
- **计算方法：** 两条记录需指向彼此、位于同一染色体，并携带一致的 SVTYPE；调用 varcode 的成对断点判定。
- **API：** `dnakit.variants.classify_breakend_pair(first: Mapping[str, Any], second: Mapping[str, Any])`
- **输入：** 两条断点记录，每条含 contig、start、ref、alt 与可选 info。
- **输出：** 返回 `typed`、`sv_type`、`start` 与 `end`（`start` 为事件前的碱基）。
- **示例代码：**

```python
from dnakit.variants import classify_breakend_pair

first = {"contig": "chr1", "start": 100, "ref": "A", "alt": "A]chr1:500]",
         "info": {"SVTYPE": "DEL", "MATEID": "b2"}}
second = {"contig": "chr1", "start": 500, "ref": "T", "alt": "[chr1:100[T",
          "info": {"SVTYPE": "DEL", "MATEID": "b1"}}
typed = classify_breakend_pair(first, second)
print(typed["sv_type"], typed["start"], typed["end"])
```

- **示例输出：**

```text
DEL 100 499
```

## 8) 比较结构变异调用的断点

- **作用：** 把断点相近的结构变异调用归为一组，并保留每条输入记录。
- **计算方法：** 按规范化后的断点拓扑分桶，同一桶内用完全链接的贪心分组，断点距离超过 `max_distance` 不合并。
- **API：** `dnakit.variants.compare_sv_calls(rows: Sequence[Mapping[str, Any]], *, max_distance: int = 100)`
- **输入：** 调用记录列表；每条需含 call_id、caller、sample、build、chrom、pos、ref、alt，build 为 GRCh37 或 GRCh38。
- **输出：** 返回 `result['groups']`，每组含 call_ids、callers、status 等字段。
- **示例代码：**

```python
from dnakit.variants import compare_sv_calls

rows = [
    {"call_id": "a", "caller": "x", "sample": "s", "build": "GRCh38", "chrom": "chr1",
     "pos": 100, "ref": "A", "alt": "<DEL>", "info": {"END": 500}},
    {"call_id": "b", "caller": "y", "sample": "s", "build": "GRCh38", "chrom": "chr1",
     "pos": 100, "ref": "A", "alt": "<DEL>", "info": {"END": 500}},
    {"call_id": "c", "caller": "x", "sample": "s", "build": "GRCh38", "chrom": "chr1",
     "pos": 5000, "ref": "A", "alt": "<DEL>", "info": {"END": 5400}},
]
result = compare_sv_calls(rows, max_distance=100)
print(sorted(group["call_ids"] for group in result["result"]["groups"]))
```

- **示例输出：**

```text
[['a', 'b'], ['c']]
```

## 9) 去除等位基因共享的两端

- **作用：** 去掉 REF 与 ALT 共有的前缀和后缀，得到最小等位基因表示。
- **计算方法：** 先去共同前缀，再在剩余部分去共同后缀，与 varcode 分类前的修剪一致。
- **API：** `dnakit.variants.trim_shared_allele_flanks(ref: str, alt: str)`
- **输入：** REF 与 ALT 字符串（可为空）。
- **输出：** 返回修剪后的 `trimmed_ref`、`trimmed_alt`，以及 `prefix` 与 `suffix`。
- **示例代码：**

```python
from dnakit.variants import trim_shared_allele_flanks

result = trim_shared_allele_flanks("ACGT", "ACTT")
print(result["trimmed_ref"], result["trimmed_alt"], result["prefix"], result["suffix"])
```

- **示例输出：**

```text
G T AC T
```

## 10) 解析 VCF 基因型 GT

- **作用：** 把 `0|1`、`./.` 等 GT 字符串解析为等位基因编号与相位状态。
- **计算方法：** 先用 VCF GT 正则检查格式，再调用 varcode 解析；缺失等位基因为 None。
- **API：** `dnakit.variants.parse_vcf_genotype(genotype: str)`
- **输入：** GT 字符串。
- **输出：** 返回 `alleles` 列表、`phased` 与 `ploidy`。
- **示例代码：**

```python
from dnakit.variants import parse_vcf_genotype

genotype = parse_vcf_genotype("1|.")
print(genotype["alleles"], genotype["phased"], genotype["ploidy"])
```

- **示例输出：**

```text
[1, None] True 2
```

## 11) 规范化等位基因字符串

- **作用：** 把等位基因转为大写，并把间隔字符（如 `.`、`-`）记为空等位基因。
- **计算方法：** 调用 varcode 的核苷酸规范化函数；默认只接受标准碱基与间隔字符。
- **API：** `dnakit.variants.normalize_nucleotide_allele(allele: str, *, allow_extended_nucleotides: bool = False, empty_chars: str = '.-', treat_nan_as_empty: bool = True)`
- **输入：** 等位基因字符串；可选 `allow_extended_nucleotides` 与 `empty_chars`。
- **输出：** 返回 `normalized` 与 `is_empty`。
- **示例代码：**

```python
from dnakit.variants import normalize_nucleotide_allele

print(repr(normalize_nucleotide_allele("acgt")["normalized"]))
print(normalize_nucleotide_allele(".")["is_empty"])
```

- **示例输出：**

```text
'ACGT'
True
```

## 12) 判断变异是否与区间重叠

- **作用：** 检查长度为 n 的变异是否与 1-based 闭区间相交。
- **计算方法：** 变异占据 start 到 start+n-1；插入（n=0）只在区间包含插入点时重叠。
- **API：** `dnakit.variants.variant_overlaps_interval(variant_start: int, n_ref_bases: int, interval_start: int, interval_end: int)`
- **输入：** 变异起点、参考碱基数、区间起点与终点（均为 1-based 闭区间）。
- **输出：** 返回 `overlaps`，以及变异终点 `variant_end`。
- **示例代码：**

```python
from dnakit.variants import variant_overlaps_interval

print(variant_overlaps_interval(10, 3, 12, 20)["overlaps"])
print(variant_overlaps_interval(10, 3, 13, 20)["overlaps"])
```

- **示例输出：**

```text
True
False
```

## 13) 四配子检验判断顺反式

- **作用：** 根据两个变异的四类读段片段计数判断它们位于同一条染色体（cis）还是不同染色体（trans）。
- **计算方法：** 调用 varcode 的四配子检验；计数对应 AB、Ab、aB、ab 四种组合，显著性由 `max_p_value` 控制，证据不足时返回 unknown。
- **API：** `dnakit.variants.four_gamete_phase_test(both: int, first: int, second: int, neither: int, *, min_fragments: int = 2, error_rate: float = 0.01, max_p_value: float = 0.05)`
- **输入：** 四个非负整数计数（both、first、second、neither）与检验参数。
- **输出：** 返回 `relation`（cis、trans 或 unknown）与 `in_cis`。
- **示例代码：**

```python
from dnakit.variants import four_gamete_phase_test

cis = four_gamete_phase_test(10, 1, 2, 10)
print(cis["relation"], cis["in_cis"])
```

- **示例输出：**

```text
cis True
```

## 14) 检测体细胞变异与种系变异的等位基因重叠

- **作用：** 找出与体细胞变异具有相同染色体、位置、REF 与 ALT 的种系变异。
- **计算方法：** 对每条种系变异调用 varcode 的重叠判断，不考虑基因型、等位基因比例或 LOH。
- **API：** `dnakit.variants.detect_germline_allele_overlap(somatic: Mapping[str, Any], germline: Sequence[Mapping[str, Any]])`
- **输入：** 体细胞变异与种系变异列表；每条含 contig、start、ref、alt。
- **输出：** 返回 `overlap`、`matched_indices` 与 `germline_count`。
- **示例代码：**

```python
from dnakit.variants import detect_germline_allele_overlap

somatic = {"contig": "chr1", "start": 100, "ref": "A", "alt": "G"}
germline = [{"contig": "chr1", "start": 100, "ref": "A", "alt": "G"},
            {"contig": "chr1", "start": 200, "ref": "C", "alt": "T"}]
result = detect_germline_allele_overlap(somatic, germline)
print(result["overlap"], result["matched_indices"])
```

- **示例输出：**

```text
True [0]
```

## 15) 按相位划分种系变异

- **作用：** 把种系变异划分为与体细胞变异同向（cis）、异向（trans）、相连或未定相的组。
- **计算方法：** 纯合变异恒为 cis；其余变异按外部给出的同向证据做并查集传递，缺证据的留在 `unphased`。
- **API：** `dnakit.variants.partition_germline_by_phase(somatic: Mapping[str, Any], germline: Sequence[Mapping[str, Any]], *, homozygous: Sequence[Mapping[str, Any]] | None = None, cis_variants: Sequence[Mapping[str, Any]] | None = None, in_cis_evidence: Sequence[Mapping[str, Any]] | None = None)`
- **输入：** 体细胞变异、种系变异，以及可选的纯合变异与 `in_cis_evidence`（每条含 first、second 与 in_cis）。
- **输出：** 返回 `cis`、`trans`、`linked` 与 `unphased`。
- **示例代码：**

```python
from dnakit.variants import partition_germline_by_phase

somatic = {"contig": "chr1", "start": 100, "ref": "A", "alt": "G"}
germline = {"contig": "chr1", "start": 200, "ref": "C", "alt": "T"}
evidence = [{"first": somatic, "second": germline, "in_cis": False}]
partition = partition_germline_by_phase(somatic, [germline], in_cis_evidence=evidence)
print(partition["cis"], partition["trans"], partition["unphased"])
```

- **示例输出：**

```text
[] [{'contig': 'chr1', 'start': 200, 'ref': 'C', 'alt': 'T'}] []
```

## 16) 列出相位假设

- **作用：** 列出体细胞变异相对种系变异可能的相位构型。
- **计算方法：** 先做相位划分，再按已知同向证据生成假设；超过 `max_hypotheses` 时报错。
- **API：** `dnakit.variants.enumerate_phase_hypotheses(somatic: Mapping[str, Any], germline: Sequence[Mapping[str, Any]], *, in_cis_evidence: Sequence[Mapping[str, Any]] | None = None, max_hypotheses: int = 8)`
- **输入：** 体细胞变异、种系变异、可选的 `in_cis_evidence` 与 `max_hypotheses`（1 到 64）。
- **输出：** 返回 `count` 与 `hypotheses`，每项含 phase_state、haplotype、cis 与 trans。
- **示例代码：**

```python
from dnakit.variants import enumerate_phase_hypotheses

somatic = {"contig": "chr1", "start": 100, "ref": "A", "alt": "G"}
germline = {"contig": "chr1", "start": 200, "ref": "C", "alt": "T"}
evidence = [{"first": somatic, "second": germline, "in_cis": True}]
result = enumerate_phase_hypotheses(somatic, [germline], in_cis_evidence=evidence)
print(result["count"], result["hypotheses"][0]["phase_state"])
```

- **示例输出：**

```text
1 phased
```

## 17) 读取本地 MAF 文件

- **作用：** 读取 TCGA 风格的 MAF 文件，得到每条变异的染色体、位置、REF 与 ALT。
- **计算方法：** 调用 varcode 的 MAF 读取器，按染色体与位置排序，并保留 reference_name。
- **API：** `dnakit.variants.load_maf_variants(path: str)`
- **输入：** 本地 MAF 文件路径，需含 varcode 要求的 19 个标准列。
- **输出：** 返回 `reference`、`count` 与 `variants`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import load_maf_variants

COLUMNS = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position",
           "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
           "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status",
           "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1",
           "Match_Norm_Seq_Allele2"]
ROWS = [
    {"Hugo_Symbol": "TP53", "NCBI_Build": "GRCh38", "Chromosome": "17", "Start_Position": "7676154",
     "End_Position": "7676154", "Strand": "+", "Variant_Type": "SNP", "Reference_Allele": "C",
     "Tumor_Seq_Allele1": "C", "Tumor_Seq_Allele2": "T", "Tumor_Sample_Barcode": "S1"},
    {"Hugo_Symbol": "KRAS", "NCBI_Build": "GRCh38", "Chromosome": "12", "Start_Position": "25245350",
     "End_Position": "25245350", "Strand": "+", "Variant_Type": "SNP", "Reference_Allele": "C",
     "Tumor_Seq_Allele1": "C", "Tumor_Seq_Allele2": "A", "Tumor_Sample_Barcode": "S1"},
]
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "sample.maf"
    lines = ["\t".join(COLUMNS)]
    lines += ["\t".join(row.get(column, "") for column in COLUMNS) for row in ROWS]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    result = load_maf_variants(str(path))
print(result["reference"], result["count"], result["variants"][0]["contig"], result["variants"][0]["start"])
```

- **示例输出：**

```text
GRCh38 2 12 25245350
```

## 18) 检查 VCF 样本身份

- **作用：** 读取本地 VCF 中的样本指纹，并报告样本两两之间的种系一致性。
- **计算方法：** 调用 varcode 的样本指纹读取器与身份检查；结论为启发式筛查，不是身份概率。
- **API：** `dnakit.variants.check_vcf_sample_identity(path: str, *, kind: str = 'germline', assembly: str | None = None, config: Mapping[str, Any] | None = None)`
- **输入：** 本地 VCF 路径；可选 `kind`（germline、somatic 或 mixed）、`assembly` 与 `config` 覆盖阈值。
- **输出：** 返回 `sample_count` 与 `report`，报告中含各样本对的种系统计。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import check_vcf_sample_identity

VCF = (
    "##fileformat=VCFv4.2\n"
    "##reference=GRCh38\n"
    "##contig=<ID=chr1,length=248956422>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30\t0/1:30\n"
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:DP\t1/1:30\t1/1:30\n"
)
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "pair.vcf"
    path.write_text(VCF, encoding="utf-8")
    result = check_vcf_sample_identity(str(path))
germline = result["report"]["pairs"][0]["germline"]
print(result["sample_count"], germline["shared_snps"], germline["exact_matches"])
```

- **示例输出：**

```text
2 2 2
```

## 19) 分窗口计算 Tajima's D

- **作用：** 按碱基对窗口计算 Tajima's D，用于检测同一条染色体上的选择信号。
- **计算方法：** 调用 scikit-allel 的 windowed_tajima_d，窗口按位置划分，每窗口至少含 `min_sites` 个位点。
- **API：** `dnakit.variants.calculate_windowed_tajimas_d(positions: Sequence[int], allele_counts: Sequence[Sequence[int]], *, window_size: int, step: int | None = None, min_sites: int = 3)`
- **输入：** 按位置升序排列的 1-based 位点，以及每个位点的等位基因计数 (variants x alleles)。
- **输出：** 返回 `windows` 列表，每项含 start、stop、sites 与 tajima_d。
- **示例代码：**

```python
from dnakit.variants import calculate_windowed_tajimas_d

positions = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
counts = [[3, 1], [2, 2], [4, 0], [1, 3], [0, 4], [2, 2], [3, 1], [1, 3], [2, 2], [4, 0]]
result = calculate_windowed_tajimas_d(positions, counts, window_size=500)
for window in result["windows"]:
    print(window["start"], window["stop"], window["sites"], round(window["tajima_d"], 6))
```

- **示例输出：**

```text
100 599 5 0.167656
600 1000 5 0.650102
```

## 20) 分窗口计算核苷酸多样性

- **作用：** 按窗口计算核苷酸多样性 pi，窗口内的可调用碱基数也一并输出。
- **计算方法：** 调用 scikit-allel 的 windowed_diversity；pi 为窗口内各位点平均两两差异之和除以可调用碱基数。
- **API：** `dnakit.variants.calculate_windowed_nucleotide_diversity(positions: Sequence[int], allele_counts: Sequence[Sequence[int]], *, window_size: int, step: int | None = None)`
- **输入：** 位点位置、等位基因计数与窗口大小（碱基对）。
- **输出：** 返回 `windows` 列表，每项含 start、stop、callable_bases、sites 与 pi。
- **示例代码：**

```python
from dnakit.variants import calculate_windowed_nucleotide_diversity

positions = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
counts = [[3, 1], [2, 2], [4, 0], [1, 3], [0, 4], [2, 2], [3, 1], [1, 3], [2, 2], [4, 0]]
result = calculate_windowed_nucleotide_diversity(positions, counts, window_size=500)
for window in result["windows"]:
    print(window["start"], window["stop"], window["callable_bases"], round(window["pi"], 8))
```

- **示例输出：**

```text
100 599 500 0.00333333
600 1000 401 0.00581879
```

## 21) 计算每个位点的平均两两差异

- **作用：** 对每个位点计算等位基因两两比较的平均差异数。
- **计算方法：** 由等位基因计数计算 2·n_a·n_b / (n·(n−1))，不同等位基因计数组合的结果由 scikit-allel 给出。
- **API：** `dnakit.variants.calculate_mean_pairwise_difference(allele_counts: Sequence[Sequence[int]])`
- **输入：** 每个位点的等位基因计数 (variants x alleles)。
- **输出：** 返回 `values` 列表，与输入位点一一对应。
- **示例代码：**

```python
from dnakit.variants import calculate_mean_pairwise_difference

result = calculate_mean_pairwise_difference([[3, 1], [2, 2], [4, 0]])
print([round(value, 6) for value in result["values"]])
```

- **示例输出：**

```text
[0.5, 0.666667, 0.0]
```

## 22) 按位点计算 Hudson Fst

- **作用：** 按位点给出 Hudson Fst 的分子与分母，并汇总为总体 Fst。
- **计算方法：** 调用 scikit-allel 的 hudson_fst；总体值为分子之和除以分母之和。
- **API：** `dnakit.variants.calculate_hudson_fst_per_site(first_counts: Sequence[Sequence[int]], second_counts: Sequence[Sequence[int]])`
- **输入：** 两个群体的等位基因计数，位点顺序一致。
- **输出：** 返回 `numerator`、`denominator` 与 `fst`。
- **示例代码：**

```python
from dnakit.variants import calculate_hudson_fst_per_site

first = [[4, 0], [3, 1], [2, 2], [0, 4]]
second = [[0, 4], [1, 3], [2, 2], [4, 0]]
result = calculate_hudson_fst_per_site(first, second)
print([round(value, 6) for value in result["numerator"]], round(result["fst"], 6))
```

- **示例输出：**

```text
[1.0, 0.125, -0.166667, 1.0] 0.626667
```

## 23) 计算 Patterson's D 统计量

- **作用：** 用四个群体的等位基因计数计算 ABBA-BABA 检验的 D 统计量。
- **计算方法：** 调用 scikit-allel 的 patterson_d，分子为 (a−b)(c−d)，分母为对应的方差项；D 为分子之和除以分母之和。
- **API：** `dnakit.variants.calculate_patterson_d(first_counts: Sequence[Sequence[int]], second_counts: Sequence[Sequence[int]], third_counts: Sequence[Sequence[int]], fourth_counts: Sequence[Sequence[int]])`
- **输入：** 四个群体 A、B、C、D 的等位基因计数，位点顺序一致。
- **输出：** 返回 `numerator`、`denominator` 与 `d`。
- **示例代码：**

```python
from dnakit.variants import calculate_patterson_d

result = calculate_patterson_d([[4, 0]], [[0, 4]], [[4, 0]], [[0, 4]])
print(result["numerator"], result["denominator"], result["d"])
```

- **示例输出：**

```text
[1.0] [1.0] 1.0
```

## 24) 计算期望杂合度

- **作用：** 由每个位点的等位基因频率计算期望杂合度 1 − Σp²。
- **计算方法：** 对二倍体按 Hardy-Weinberg 期望计算；缺失位点按 scikit-allel 的规则给出 None。
- **API：** `dnakit.variants.calculate_expected_heterozygosity(allele_frequencies: Sequence[Sequence[float]])`
- **输入：** 每个位点的等位基因频率行向量，每行和为 1。
- **输出：** 返回 `values` 与 `ploidy`（固定为 2）。
- **示例代码：**

```python
from dnakit.variants import calculate_expected_heterozygosity

result = calculate_expected_heterozygosity([[0.5, 0.5], [0.75, 0.25], [1.0, 0.0]])
print(result["values"])
```

- **示例输出：**

```text
[0.5, 0.375, 0.0]
```

## 25) 计算近交系数

- **作用：** 按位点计算近交系数 F = 1 − Ho / He。
- **计算方法：** 由二倍体基因型计算观测杂合率 Ho 与期望杂合率 He；He 为零时返回 None。
- **API：** `dnakit.variants.calculate_inbreeding_coefficient(genotypes: Sequence[Sequence[Sequence[int]]])`
- **输入：** 二倍体基因型数组 (variants x samples x 2)，缺失为 −1。
- **输出：** 返回每个位点的 `values`。
- **示例代码：**

```python
from dnakit.variants import calculate_inbreeding_coefficient

genotypes = [[[0, 0], [0, 1]], [[1, 1], [0, 1]], [[0, 1], [1, 1]], [[0, 0], [0, 0]]]
print([None if value is None else round(value, 6) for value in calculate_inbreeding_coefficient(genotypes)["values"]])
```

- **示例输出：**

```text
[-0.333333, -0.333333, -0.333333, None]
```

## 26) 计算 Rogers-Huff r

- **作用：** 计算每对位点之间的 Rogers-Huff r，即基因型剂量的相关系数。
- **计算方法：** 对每个位点的剂量向量做 Pearson 相关；任一位点方差为零时返回 None。
- **API：** `dnakit.variants.calculate_rogers_huff_r(alt_counts: Sequence[Sequence[int]])`
- **输入：** 等位基因剂量矩阵 (variants x samples)，取值 0、1、2 或缺失 −1，位点数不超过 200。
- **输出：** 返回 `pairs` 列表，每项含 variant_a、variant_b 与 r。
- **示例代码：**

```python
from dnakit.variants import calculate_rogers_huff_r

result = calculate_rogers_huff_r([[0, 1], [2, 1], [1, 2], [0, 0]])
for pair in result["pairs"]:
    print(pair["variant_a"], pair["variant_b"], pair["r"])
```

- **示例输出：**

```text
0 1 -1.0
0 2 1.0
0 3 None
1 2 -1.0
1 3 None
2 3 None
```

## 27) 选取近似连锁平衡的位点

- **作用：** 用滑动窗口的 r² 剪枝，保留近似连锁平衡的位点。
- **计算方法：** 调用 scikit-allel 的 locate_unlinked；窗口大小与步长按位点数计算，r² 超过阈值的位点被剔除。
- **API：** `dnakit.variants.select_unlinked_variants(alt_counts: Sequence[Sequence[int]], *, size: int = 100, step: int = 20, threshold: float = 0.1)`
- **输入：** 剂量矩阵 (variants x samples)，以及窗口大小、步长与阈值。
- **输出：** 返回 `kept_indices`（保留位点的下标，升序）与 `kept_count`。
- **示例代码：**

```python
import numpy as np

from dnakit.variants import select_unlinked_variants

rng = np.random.default_rng(1)
dosage = rng.integers(0, 3, size=(30, 8))
result = select_unlinked_variants(dosage.tolist(), size=10, step=5, threshold=0.1)
print(result["kept_indices"])
```

- **示例输出：**

```text
[0, 1, 5, 12, 14, 15, 20, 27, 28]
```

## 28) 计算样本两两距离

- **作用：** 用位点剂量计算样本之间的距离矩阵。
- **计算方法：** 调用 scikit-allel 的 pairwise_distance 得到压缩距离，再展开为对称方阵；不允许缺失值。
- **API：** `dnakit.variants.calculate_pairwise_distance(alt_counts: Sequence[Sequence[int]], *, metric: str = 'euclidean')`
- **输入：** 剂量矩阵 (variants x samples)，取值非负；可选距离度量。
- **输出：** 返回 `matrix`（样本 x 样本）与 `metric`。
- **示例代码：**

```python
from dnakit.variants import calculate_pairwise_distance

result = calculate_pairwise_distance([[0, 1, 2], [1, 1, 2], [0, 0, 1]])
print([[round(value, 6) for value in row] for row in result["matrix"]])
```

- **示例输出：**

```text
[[0.0, 1.0, 2.44949], [1.0, 0.0, 1.732051], [2.44949, 1.732051, 0.0]]
```

## 29) 计算未标准化的 iHS

- **作用：** 对每个位点计算单群体的未标准化综合单倍型评分（iHS）。
- **计算方法：** 调用 scikit-allel 的 ihs，按 EHH 衰减积分；等位基因频率低于 `min_maf` 的位点不计算。
- **API：** `dnakit.variants.calculate_ihs(haplotypes: Sequence[Sequence[int]], positions: Sequence[int], *, min_ehh: float = 0.05, min_maf: float = 0.05, gap_scale: int = 20000, max_gap: int = 200000, include_edges: bool = False)`
- **输入：** 单倍型矩阵 (variants x haplotypes，取值 0 或 1)、位置与 EHH 参数。
- **输出：** 返回 `ihs` 列表，与位点一一对应。
- **示例代码：**

```python
import numpy as np

from dnakit.variants import calculate_ihs

rng = np.random.default_rng(3)
haplotypes = rng.integers(0, 2, size=(40, 12))
positions = list(range(1000, 1000 + 100 * 40, 100))
values = calculate_ihs(haplotypes.tolist(), positions)["ihs"]
print(len(values), [round(value, 6) for value in values[5:8]])
```

- **示例输出：**

```text
40 [0.143101, -0.169076, 0.077962]
```

## 30) 计算未标准化的 XP-EHH

- **作用：** 比较两个群体的单倍型，计算每个位点的跨群体 EHH 评分。
- **计算方法：** 调用 scikit-allel 的 xpehh；两组单倍型的位点数必须一致。
- **API：** `dnakit.variants.calculate_xpehh(haplotypes_a: Sequence[Sequence[int]], haplotypes_b: Sequence[Sequence[int]], positions: Sequence[int], *, min_ehh: float = 0.05, gap_scale: int = 20000, max_gap: int = 200000, include_edges: bool = False)`
- **输入：** 两组单倍型矩阵、位置与 EHH 参数。
- **输出：** 返回 `xpehh` 列表。
- **示例代码：**

```python
import numpy as np

from dnakit.variants import calculate_xpehh

rng = np.random.default_rng(4)
first = rng.integers(0, 2, size=(40, 12))
second = rng.integers(0, 2, size=(40, 12))
positions = list(range(1000, 1000 + 100 * 40, 100))
values = calculate_xpehh(first.tolist(), second.tolist(), positions)["xpehh"]
print(len(values), [round(value, 6) for value in values[5:8]])
```

- **示例输出：**

```text
40 [0.157004, 0.055152, 0.196028]
```

## 31) 计算未标准化的 nSL

- **作用：** 按位点计算按长度加权的分离位点数（nSL），用于检测选择信号。
- **计算方法：** 调用 scikit-allel 的 nsl，逐位点比较单倍型之间共享的前缀长度。
- **API：** `dnakit.variants.calculate_nsl(haplotypes: Sequence[Sequence[int]])`
- **输入：** 单倍型矩阵 (variants x haplotypes，取值 0 或 1)。
- **输出：** 返回 `nsl` 列表。
- **示例代码：**

```python
import numpy as np

from dnakit.variants import calculate_nsl

rng = np.random.default_rng(6)
haplotypes = rng.integers(0, 2, size=(40, 12))
values = calculate_nsl(haplotypes.tolist())["nsl"]
print(len(values), [round(value, 6) for value in values[:3]])
```

- **示例输出：**

```text
40 [-0.02353, 0.092288, 0.059423]
```

## 32) 计算 Garud H 统计量

- **作用：** 从单倍型频率计算 H1、H12、H123 与 H2/H1，用于检测软扫荡。
- **计算方法：** 调用 scikit-allel 的 garud_h：H1 为单倍型频率平方和，H12 与 H123 分别合并最频繁的前两、三个单倍型。
- **API：** `dnakit.variants.calculate_garud_h(haplotypes: Sequence[Sequence[int]])`
- **输入：** 单倍型矩阵 (variants x haplotypes，取值 0 或 1)。
- **输出：** 返回 `h1`、`h12`、`h123` 与 `h2_h1`。
- **示例代码：**

```python
import numpy as np

from dnakit.variants import calculate_garud_h

rng = np.random.default_rng(7)
haplotypes = rng.integers(0, 2, size=(40, 12))
result = calculate_garud_h(haplotypes.tolist())
print(round(result["h1"], 6), round(result["h12"], 6), round(result["h2_h1"], 6))
```

- **示例输出：**

```text
0.083333 0.097222 0.916667
```

## 33) 计算 EHH 衰减曲线

- **作用：** 计算从第一个位点向外扩展时，单倍型纯合度（EHH）的衰减。
- **计算方法：** 调用 scikit-allel 的 ehh_decay；`truncate=True` 时去掉末尾的零值。
- **API：** `dnakit.variants.calculate_ehh_decay(haplotypes: Sequence[Sequence[int]], *, truncate: bool = False)`
- **输入：** 单倍型矩阵 (variants x haplotypes，取值 0 或 1)。
- **输出：** 返回 `ehh` 列表，按距离递减排列。
- **示例代码：**

```python
import numpy as np

from dnakit.variants import calculate_ehh_decay

rng = np.random.default_rng(8)
haplotypes = rng.integers(0, 2, size=(40, 12))
values = calculate_ehh_decay(haplotypes.tolist())["ehh"]
print([round(value, 6) for value in values[:4]])
```

- **示例输出：**

```text
[0.454545, 0.242424, 0.121212, 0.030303]
```

## 34) 查找孟德尔遗传错误

- **作用：** 统计三代家系中与孟德尔遗传不一致的基因型数。
- **计算方法：** 调用 scikit-allel 的 mendel_errors；父母为两列（母亲在前），子代按列给出。
- **API：** `dnakit.variants.find_mendel_errors(parents: Sequence[Sequence[Sequence[int]]], progeny: Sequence[Sequence[Sequence[int]]])`
- **输入：** 父母基因型 (variants x 2 x 2) 与子代基因型 (variants x n x 2)。
- **输出：** 返回每个位点、每个子代的 `errors`，以及 `total_errors`。
- **示例代码：**

```python
from dnakit.variants import find_mendel_errors

parents = [[[0, 0], [1, 1]], [[0, 1], [0, 1]], [[1, 1], [0, 0]]]
progeny = [[[0, 1]], [[1, 1]], [[0, 0]]]
result = find_mendel_errors(parents, progeny)
print(result["errors"], result["total_errors"])
```

- **示例输出：**

```text
[[0], [0], [1]] 1
```

## 35) 按遗传方式为三代家系定相

- **作用：** 利用孟德尔遗传把母亲、父亲与子代的基因型定相。
- **计算方法：** 调用 scikit-allel 的 phase_by_transmission，在 `window_size` 个位点的窗口内求解可确定的相位。
- **API：** `dnakit.variants.phase_trio_by_transmission(genotypes: Sequence[Sequence[Sequence[int]]], window_size: int)`
- **输入：** 三代基因型 (variants x 3 x 2)，样本顺序为母亲、父亲、子代；窗口大小按位点数。
- **输出：** 返回 `genotypes`（等位基因编号）与 `phased` 标记矩阵。
- **示例代码：**

```python
from dnakit.variants import phase_trio_by_transmission

genotypes = [[[0, 1], [0, 0], [0, 1]], [[1, 1], [0, 1], [1, 1]], [[0, 0], [1, 0], [0, 0]], [[0, 1], [0, 1], [1, 1]]]
result = phase_trio_by_transmission(genotypes, 2)
print(result["genotypes"][0], result["phased"][0])
```

- **示例输出：**

```text
[[0, 1], [0, 0], [1, 0]] [True, True, True]
```

## 36) 用隐马尔可夫模型识别纯合片段

- **作用：** 从单个二倍体个体的基因型识别纯合片段（ROH）。
- **计算方法：** 调用 scikit-allel 的 roh_mhmm；多项隐马尔可夫模型用杂合概率与转移概率区分 ROH 与非 ROH 状态。
- **API：** `dnakit.variants.call_roh_hmm(genotypes: Sequence[Sequence[int]], positions: Sequence[int], *, phet_roh: float = 0.001, transition: float = 1e-06, min_roh: int = 0, contig_size: int | None = None)`
- **输入：** 单个个体的基因型 (variants x 2)、位置、杂合概率与转移概率；可选 `contig_size`。
- **输出：** 返回 `runs`（start、stop、length、is_marginal）与 `froh`。
- **示例代码：**

```python
import numpy as np

from dnakit.variants import call_roh_hmm

rng = np.random.default_rng(11)
calls = rng.choice(np.array([[0, 0], [0, 1]]), size=3000, p=[0.5, 0.5])
calls[1000:2000] = 0
positions = (np.arange(1, 3001) * 10).tolist()
result = call_roh_hmm(calls.tolist(), positions)
longest = max(result["runs"], key=lambda run: run["length"])
print(len(result["runs"]) > 0, longest["start"] <= 10100, longest["stop"] >= 19900)
```

- **示例输出：**

```text
True True True
```

## 37) 计算区域核苷酸分歧 dxy

- **作用：** 计算两个群体在一段区域内的平均核苷酸分歧 dxy。
- **计算方法：** 对区域内每个位点求群体间两两差异的平均值之和，再除以区域碱基数（默认为首末位点之间）。
- **API：** `dnakit.variants.calculate_sequence_divergence(positions: Sequence[int], first_counts: Sequence[Sequence[int]], second_counts: Sequence[Sequence[int]], *, start: int | None = None, stop: int | None = None)`
- **输入：** 位点位置、两个群体的等位基因计数；可选区域起点与终点。
- **输出：** 返回 `dxy`、区域 `start` 与 `stop`。
- **示例代码：**

```python
from dnakit.variants import calculate_sequence_divergence

positions = [100, 200, 300, 400]
first = [[4, 0], [3, 1], [2, 2], [0, 4]]
second = [[0, 4], [1, 3], [2, 2], [4, 0]]
result = calculate_sequence_divergence(positions, first, second)
print(round(result["dxy"], 8), result["start"], result["stop"])
```

- **示例输出：**

```text
0.01038206 100 400
```

## 38) 计算群体分支统计量 PBS

- **作用：** 用三个群体在滑动窗口内计算 PBS，用于定位某一群体特有的分化信号。
- **计算方法：** 调用 scikit-allel 的 pbs，由三对两两 Fst 经对数变换组合而成；窗口大小按位点数计算。
- **API：** `dnakit.variants.calculate_population_branching_statistic(first_counts: Sequence[Sequence[int]], second_counts: Sequence[Sequence[int]], third_counts: Sequence[Sequence[int]], *, window_size: int, window_step: int | None = None, normed: bool = True)`
- **输入：** 三个群体的等位基因计数、窗口大小（位点数）、步长与是否标准化。
- **输出：** 返回每个窗口的 `values`。
- **示例代码：**

```python
from dnakit.variants import calculate_population_branching_statistic

counts = [[3, 1], [2, 2], [4, 0], [1, 3], [0, 4], [2, 2], [3, 1], [1, 3], [2, 2], [4, 0]]
result = calculate_population_branching_statistic(counts, counts, counts, window_size=2)
print(result["values"])
```

- **示例输出：**

```text
[0.0, 0.0, 0.0, 0.0, 0.0]
```

## 39) 规范化 indel 等位基因位置

- **作用：** 把 indel 平移到参考序列中最左（或最右）的等价位置。
- **计算方法：** 调用 hgvs 的等位基因规范化：先去共同后缀与前缀，再沿参考序列逐步平移直到不能再平移。坐标为 0-based 半开区间。
- **API：** `dnakit.variants.normalize_indel_alleles(reference: str, start: int, stop: int, alleles: Sequence[str], *, left: bool = True, bound: int | None = None, ref_step: int = 1, shuffle: bool = True)`
- **输入：** 参考序列、起止坐标、两个等位基因字符串；可设方向 `left`、边界 `bound` 与步长 `ref_step`。
- **输出：** 返回规范化后的 `start`、`stop`、`alleles` 与是否变化 `changed`。
- **示例代码：**

```python
from dnakit.variants import normalize_indel_alleles

left = normalize_indel_alleles("ACGTTTTA", 5, 6, ("T", ""))
right = normalize_indel_alleles("ACGTTTTA", 3, 4, ("T", ""), left=False, bound=8)
print(left["start"], left["stop"], right["start"], right["stop"])
```

- **示例输出：**

```text
3 4 6 7
```

## 40) 去除等位基因字符串的共同前缀

- **作用：** 去掉一组等位基因字符串的共同前缀，并报告去掉的长度。
- **计算方法：** 调用 hgvs 的前缀修剪，用字符串的最小与最大值比较前缀。
- **API：** `dnakit.variants.trim_common_prefix_strings(strings: Sequence[str], *, min_length: int = 0)`
- **输入：** 至少两个等位基因字符串；可选 `min_length` 保留的最短长度。
- **输出：** 返回 `trimmed_length` 与修剪后的 `strings`。
- **示例代码：**

```python
from dnakit.variants import trim_common_prefix_strings

result = trim_common_prefix_strings(["ACGT", "ACTT"])
print(result["trimmed_length"], result["strings"])
```

- **示例输出：**

```text
2 ['GT', 'TT']
```

## 41) 拼接转录本 CIGAR

- **作用：** 把各外显子的 CIGAR 拼接为整条转录本的 CIGAR，内含子以 N 表示。
- **计算方法：** 调用 hgvs 的 build_tx_cigar；负链时逐外显子反转 CIGAR 顺序，内含子长度为相邻外显子的间隔。
- **API：** `dnakit.variants.build_transcript_cigar(exons: Sequence[Mapping[str, Any]], strand: int)`
- **输入：** 按转录本顺序排列的外显子列表（含 cigar、alt_start_i、alt_end_i）与链方向（1 或 −1）。
- **输出：** 返回 `cigar` 字符串与 `exon_count`。
- **示例代码：**

```python
from dnakit.variants import build_transcript_cigar

exons = [
    {"cigar": "10=", "alt_start_i": 0, "alt_end_i": 10},
    {"cigar": "5=", "alt_start_i": 20, "alt_end_i": 25},
]
print(build_transcript_cigar(exons, 1)["cigar"])
```

- **示例输出：**

```text
10=10N5=
```

## 42) 检查 HGVS 类型与登录号是否匹配

- **作用：** 判断 NM_、NP_、NC_ 等登录号是否适用于某个 HGVS 变异类型（c、g、m、n、p、r）。
- **计算方法：** 调用 hgvs 的登录号校验，按正则结果返回 VALID、WARNING 或 ERROR。
- **API：** `dnakit.variants.validate_hgvs_type_accession(variant_type: str, accession: str)`
- **输入：** 变异类型字母与登录号字符串。
- **输出：** 返回 `level`（VALID、WARNING 或 ERROR）与 `message`。
- **示例代码：**

```python
from dnakit.variants import validate_hgvs_type_accession

print(validate_hgvs_type_accession("c", "NM_000088.3")["level"])
print(validate_hgvs_type_accession("g", "NM_000088.3")["level"])
```

- **示例输出：**

```text
VALID
ERROR
```

## 43) 检测串联重复块

- **作用：** 在序列中找出串联重复单元及其重复次数。
- **计算方法：** 调用 hgvs 的重复块检测；可从左到右或从右到左扫描，并输出每个块的单元、次数与长度。
- **API：** `dnakit.variants.detect_repeat_blocks(sequence: str, *, reverse: bool = False)`
- **输入：** 序列（仅含 A/C/G/T/N）与扫描方向 `reverse`。
- **输出：** 返回 `blocks` 列表，每项含 repeat_unit、repeat_count、block_size 与 block。
- **示例代码：**

```python
from dnakit.variants import detect_repeat_blocks

for block in detect_repeat_blocks("ACACACGTTTT")["blocks"]:
    print(block["repeat_unit"], block["repeat_count"], block["block_size"])
```

- **示例输出：**

```text
AC 3 6
T 4 4
```

## 44) 组装重复块为 HGVS 写法

- **作用：** 把检测到的重复块组装为 HGVS 风格的重复表示，例如 AC[3]G[1]T[4]。
- **计算方法：** 调用 hgvs 的重复字符串组装；未被重复块覆盖的碱基按单碱基计数。
- **API：** `dnakit.variants.assemble_repeat_blocks(sequence: str, blocks: Sequence[Mapping[str, Any]], *, reverse: bool = False)`
- **输入：** 原序列、由 detect_repeat_blocks 得到的块列表与方向 `reverse`。
- **输出：** 返回 `formatted` 字符串。
- **示例代码：**

```python
from dnakit.variants import assemble_repeat_blocks, detect_repeat_blocks

sequence = "ACACACGTTTT"
blocks = detect_repeat_blocks(sequence)["blocks"]
print(assemble_repeat_blocks(sequence, blocks)["formatted"])
```

- **示例输出：**

```text
AC[3]G[1]T[4]
```

## 限制与说明

- varcode 的 `variant_overlaps_interval` 上游把变异终点记为 `start + n_ref_bases`，比闭区间终点多一个碱基；包装层传入 `n_ref_bases - 1`，使判断符合闭区间约定。
- scikit-allel 的 `heterozygosity_expected` 把 ploidy 当作指数使用，只在二倍体下给出 Hardy-Weinberg 期望；包装层固定 ploidy 为 2。
- `enumerate_phase_hypotheses` 与 `partition_germline_by_phase` 的相位证据只能通过 `in_cis_evidence` 提供，没有 pyensembl 转录本数据时不做转录本感知的相位推断。
- hgvs 包装不调用 seqrepo 或 UTA。参考依赖的投影与比对函数（`full_house`、`variant_context_w_alignment`）未包装。
- 缺失基因型用 −1 表示；缺失统计量在输出中为 None。
