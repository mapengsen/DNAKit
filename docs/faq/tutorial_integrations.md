# 教程补充说明与集成配置

功能教程的安装、参数和运行说明。

## 反向互补：文件格式、进度与大文件处理 {#reverse-complement}

包含点号或路径分隔符的字符串按路径处理；无后缀文件请显式传 `Path`。对象输入保留原返回类型，并同步映射 feature、方向和逐碱基注释。


大型未压缩 FASTA/FNA（至少 16 MiB）在默认读取配置下会使用文件区间处理：
各工作进程直接读取自己的区间、校验并转换碱基，按原记录顺序写入临时输出的独立区间。
进程间只传递路径和偏移，不传递整条染色体；校验成功后才安装输出文件。
`workers=1` 使用相同的分块算法串行执行，`workers=-1` 在可用 CPU 和内存范围内并行。
使用多进程的脚本应保留 `if __name__ == "__main__":` 入口保护。
不规则换行、压缩文件、其他格式、显式 `read_config` 或输出字节限制继续使用通用处理路径。
读取配置和错误检查不会因启用多核而被放宽。


也可以直接处理文本或文件，不必自行循环读取、转换和写出：

```python
from dnakit.ops import reverse_complement

print(reverse_complement("AATCGC"))  # GCGATT；接受小写和空白，输出大写 IUPAC DNA
result = reverse_complement("hg38.fna")  # 自动生成 hg38_reverse_complement.fna
result = reverse_complement("sequences.txt", "rc.txt")
result = reverse_complement("sequences.csv", "rc.csv", sequence_column="sequence")
```

文件处理逐条读取记录，默认显示开始提示、已完成记录数和完成路径，可用 `show_progress=False` 关闭。
返回 `WriteResult`，其中 `record_count` 是记录数，`target_artifact.relative_path` 是相对于当前工作目录的输出路径。
未指定读取配置时，未压缩文件按实际文件大小设置读取上限，完整 hg38 无需手写三个上限参数；
gzip 默认限制解压后读取量为 20 GB，可传 `read_config=ReadConfig(...)` 显式调整。
原始输入不能作为输出，已有输出默认不覆盖；重新生成时传 `write_config=WriteConfig(overwrite=True)`。

TXT 默认每个非空行是一条 DNA 序列，输出每条序列一行，不包含记录名称。
若整个 TXT 是一条换行排版的序列，传 `read_config=ReadConfig(text_mode="sequence")`。
CSV 默认读取 `sequence` 列，可通过 `sequence_column` 指定其他列名。CSV 输出保留原列名、顺序及标签等附加列；记录注释按对应的反向互补规则映射。

## 描述符：240 项字段、计算条件与来源 {#descriptors}

结果对象可通过 `to_dict()` 转为可序列化字典。凡接口提供 `cross_gaps`，`True` 也只允许跨越可跨 Gap；显式设置 `Gap(crossable=False)` 的位置始终是硬边界。接受环状序列的描述符仍以当前 origin 为扫描边界，不自动补算末端到开头的相邻对或窗口；需要跨原点语义时应先旋转或显式线性化序列。

以下先列出模糊碱基统计、长度、组成等共 13 个可单独调用的基础功能；需要一次输出更完整的固定特征向量时，请继续查看本页后面的[“全部描述符计算”](../api/features/05_all_descriptors.md#all-descriptors)。


### 14) 全部描述符计算（240项）

**完整 240 项描述符**

一次性计算 DNA 序列的组成、统计、理化和二核苷酸等描述符，从而得到固定顺序的 240 项特征结果。

`all_descriptors()` 使用固定版本 `descriptor_schema_v1`，并同时返回不可计算原因、计算条件和来源信息。字段顺序不会随输入改变；前 180 项按输入及各自适用域计算，DNAKit 不内置二核苷酸数值表，因此后 60 项默认均为 `None`。

**最短用法**

```python
from dnakit import DNASequence
from dnakit.descriptors import all_descriptors

result = all_descriptors(DNASequence("ACGT"))
print(len(result.values))
print(result.values["epsilon260_ss_m_inverse_cm_inverse"])
print(sum(name.startswith("diprodb_") and value is None for name, value in result.values.items()))
print(result.unavailable_reasons["diprodb_twist_mean"])
```

```text
240
40300.0
60
requires an explicit user-supplied DinucleotidePropertyTable; DNAKit bundles no DiProDB numerical values
```

CLI 默认输出同一套完整字段；如需旧版四组精简结果，可显式使用 `--compact`：

```bash
dnakit describe ACGT
dnakit describe ACGT --compact
```

**API**

- `dnakit.descriptors.all_descriptors(value[必须], ambiguity_policy[可选], conditions[可选], dinucleotide_property_table[可选])`：返回全部字段；只有显式提供表时才计算后 60 项。
- `dnakit.descriptors.load_dinucleotide_property_table(path[必须])`：有界读取并严格验证用户 JSON 表，记录文件 SHA-256。
- `dnakit.descriptors.descriptor_schema_v1()`：返回 240 个不可变 `DescriptorField`，每项包含 `index`、`name`、`category`、`unit`、`formula` 和 `source`。
- `dnakit.descriptors.DESCRIPTOR_SCHEMA_V1`：同一固定 schema 的常量形式。
- `dnakit.descriptors.DESCRIPTOR_NAMES_V1`：仅包含有序字段名。

结果对象结构：

| 属性                    | 含义                                                              |
| ----------------------- | ----------------------------------------------------------------- |
| `schema_version`      | 固定为`descriptor_schema_v1`                                    |
| `sequence_id`         | 输入为`DNARecord` 时保留记录 ID，否则为 `None`                |
| `values`              | 严格按 schema 排列的 240 个值                                     |
| `unavailable_reasons` | 仅包含值为`None` 的字段，且每个字段恰有一个原因                 |
| `conditions`          | IUPAC、Gap、k-mer、ORF、重复、热力学和参数表条件                  |
| `provenance`          | DNAKit 版本、实现类型，以及用户表声明的名称、版本、来源和 SHA-256 |

**统一计算口径**

- 仅 A/C/G/T 计入 canonical 分母；默认 `ambiguity_policy="ignore"`，也可设为 `"error"`。
- k-mer、CpG/GpC、用户二核苷酸表和其他相邻词统计均允许重叠，但不会跨越 IUPAC 模糊符号或显式 Gap。
- 所有比例在分母为 0 时返回 `None`，不会用伪造的 `0` 代替；原因写入 `unavailable_reasons`。
- 用户表统计的 `sd` 使用总体标准差；表必须提供规定的 15 组属性和每组 16 个 DNA 二核苷酸有限数值。
- frame 0 使用 NCBI standard genetic code 1；完整 ORF 使用 `ATG` 起始和 `TAA/TAG/TGA` 终止，并扫描正反链共六个阅读框。
- 热力学默认条件为 37 °C、Na⁺ 0.05 M、链浓度 250 nM；Wallace 仅限 2–13 nt，SantaLucia NN 仅限 2–60 nt。
- 分子量和 ε260 只用于线性、无 Gap、无修饰、A/C/G/T DNA；ε260 是理论系数，不是实验 A260。
- 当前 LZ76 对线性、无 Gap、无模糊符号且不超过 10,000 nt 的序列计算；更长序列返回 `None`，避免二次时间复杂度失控。

**类别与字段范围**

| 序号     | 类别                                     | 字段数 |
| -------- | ---------------------------------------- | -----: |
| 1–12    | 基础长度与数据质量                       |     12 |
| 13–28   | 碱基化学分组组成                         |     16 |
| 29–112  | 1/2/3-mer 频率                           |     84 |
| 113–128 | skew、CpG/GpC 与 Chargaff 指标           |     16 |
| 129–148 | 熵、复杂度、homopolymer 与重复           |     20 |
| 149–164 | 密码子与六阅读框 ORF                     |     16 |
| 165–180 | 分子量、ε260、Tm 与 NN 热力学           |     16 |
| 181–240 | 15 组用户二核苷酸参数 × mean/sd/min/max |     60 |

**240 项完整字段表**

|   # | 字段                                         | 类别                      | 单位            | 公式                                                                                         | 来源                                                               |
| --: | -------------------------------------------- | ------------------------- | --------------- | -------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
|   1 | `symbol_length`                            | `basic`                 | `nt`          | `number of nucleotide symbols; explicit gaps excluded`                                     | DNAKit descriptor_schema_v1                                        |
|   2 | `coordinate_span`                          | `basic`                 | `nt`          | `symbol_length + sum(known gap lengths); undefined with an unknown gap`                    | DNAKit descriptor_schema_v1                                        |
|   3 | `canonical_base_count`                     | `basic`                 | `count`       | `count(A,C,G,T)`                                                                           | DNAKit descriptor_schema_v1                                        |
|   4 | `ambiguity_symbol_count`                   | `basic`                 | `count`       | `symbol_length - canonical_base_count`                                                     | DNAKit descriptor_schema_v1                                        |
|   5 | `gap_object_count`                         | `basic`                 | `count`       | `number of explicit Gap objects`                                                           | DNAKit descriptor_schema_v1                                        |
|   6 | `known_gap_nt`                             | `basic`                 | `nt`          | `sum(length of known explicit gaps)`                                                       | DNAKit descriptor_schema_v1                                        |
|   7 | `unknown_gap_count`                        | `basic`                 | `count`       | `number of explicit gaps with unknown length`                                              | DNAKit descriptor_schema_v1                                        |
|   8 | `canonical_symbol_fraction`                | `basic`                 | `fraction`    | `canonical_base_count / symbol_length`                                                     | DNAKit descriptor_schema_v1                                        |
|   9 | `ambiguity_symbol_fraction`                | `basic`                 | `fraction`    | `ambiguity_symbol_count / symbol_length`                                                   | DNAKit descriptor_schema_v1                                        |
|  10 | `known_gap_fraction`                       | `basic`                 | `fraction`    | `known_gap_nt / coordinate_span`                                                           | DNAKit descriptor_schema_v1                                        |
|  11 | `canonical_run_count`                      | `basic`                 | `count`       | `number of uninterrupted A/C/G/T runs split by ambiguity or gaps`                          | DNAKit descriptor_schema_v1                                        |
|  12 | `longest_canonical_run_nt`                 | `basic`                 | `nt`          | `maximum uninterrupted A/C/G/T run length`                                                 | DNAKit descriptor_schema_v1                                        |
|  13 | `purine_count`                             | `composition`           | `count`       | `count(A)+count(G)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  14 | `purine_fraction`                          | `composition`           | `fraction`    | `purine_count / canonical_base_count`                                                      | DNAKit descriptor_schema_v1                                        |
|  15 | `pyrimidine_count`                         | `composition`           | `count`       | `count(C)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  16 | `pyrimidine_fraction`                      | `composition`           | `fraction`    | `pyrimidine_count / canonical_base_count`                                                  | DNAKit descriptor_schema_v1                                        |
|  17 | `amino_count`                              | `composition`           | `count`       | `count(A)+count(C)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  18 | `amino_fraction`                           | `composition`           | `fraction`    | `amino_count / canonical_base_count`                                                       | DNAKit descriptor_schema_v1                                        |
|  19 | `keto_count`                               | `composition`           | `count`       | `count(G)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  20 | `keto_fraction`                            | `composition`           | `fraction`    | `keto_count / canonical_base_count`                                                        | DNAKit descriptor_schema_v1                                        |
|  21 | `weak_count`                               | `composition`           | `count`       | `count(A)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  22 | `weak_fraction`                            | `composition`           | `fraction`    | `weak_count / canonical_base_count`                                                        | DNAKit descriptor_schema_v1                                        |
|  23 | `strong_count`                             | `composition`           | `count`       | `count(C)+count(G)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  24 | `strong_fraction`                          | `composition`           | `fraction`    | `strong_count / canonical_base_count`                                                      | DNAKit descriptor_schema_v1                                        |
|  25 | `purine_pyrimidine_skew`                   | `composition`           | `ratio`       | `(purine_count-pyrimidine_count)/(purine_count+pyrimidine_count)`                          | DNAKit descriptor_schema_v1                                        |
|  26 | `amino_keto_skew`                          | `composition`           | `ratio`       | `(amino_count-keto_count)/(amino_count+keto_count)`                                        | DNAKit descriptor_schema_v1                                        |
|  27 | `weak_strong_skew`                         | `composition`           | `ratio`       | `(weak_count-strong_count)/(weak_count+strong_count)`                                      | DNAKit descriptor_schema_v1                                        |
|  28 | `gc_at_ratio`                              | `composition`           | `ratio`       | `strong_count / weak_count`                                                                | DNAKit descriptor_schema_v1                                        |
|  29 | `k1_A_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(A) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  30 | `k1_C_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(C) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  31 | `k1_G_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(G) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  32 | `k1_T_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(T) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  33 | `k2_AA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  34 | `k2_AC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  35 | `k2_AG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  36 | `k2_AT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  37 | `k2_CA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  38 | `k2_CC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  39 | `k2_CG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  40 | `k2_CT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  41 | `k2_GA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  42 | `k2_GC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  43 | `k2_GG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  44 | `k2_GT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  45 | `k2_TA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  46 | `k2_TC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  47 | `k2_TG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  48 | `k2_TT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  49 | `k3_AAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  50 | `k3_AAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  51 | `k3_AAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  52 | `k3_AAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  53 | `k3_ACA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  54 | `k3_ACC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  55 | `k3_ACG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  56 | `k3_ACT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  57 | `k3_AGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  58 | `k3_AGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  59 | `k3_AGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  60 | `k3_AGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  61 | `k3_ATA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  62 | `k3_ATC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  63 | `k3_ATG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  64 | `k3_ATT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  65 | `k3_CAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  66 | `k3_CAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  67 | `k3_CAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  68 | `k3_CAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  69 | `k3_CCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  70 | `k3_CCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  71 | `k3_CCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  72 | `k3_CCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  73 | `k3_CGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  74 | `k3_CGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  75 | `k3_CGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  76 | `k3_CGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  77 | `k3_CTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  78 | `k3_CTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  79 | `k3_CTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  80 | `k3_CTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  81 | `k3_GAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  82 | `k3_GAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  83 | `k3_GAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  84 | `k3_GAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  85 | `k3_GCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  86 | `k3_GCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  87 | `k3_GCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  88 | `k3_GCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  89 | `k3_GGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  90 | `k3_GGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  91 | `k3_GGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  92 | `k3_GGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  93 | `k3_GTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  94 | `k3_GTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  95 | `k3_GTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  96 | `k3_GTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  97 | `k3_TAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  98 | `k3_TAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  99 | `k3_TAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 100 | `k3_TAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 101 | `k3_TCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 102 | `k3_TCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 103 | `k3_TCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 104 | `k3_TCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 105 | `k3_TGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 106 | `k3_TGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 107 | `k3_TGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 108 | `k3_TGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 109 | `k3_TTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 110 | `k3_TTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 111 | `k3_TTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 112 | `k3_TTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 113 | `gc_skew`                                  | `skew_cpg`              | `ratio`       | `(count(G)-count(C))/(count(G)+count(C))`                                                  | DNAKit descriptor_schema_v1                                        |
| 114 | `at_skew`                                  | `skew_cpg`              | `ratio`       | `(count(A)-count(T))/(count(A)+count(T))`                                                  | DNAKit descriptor_schema_v1                                        |
| 115 | `cpg_count`                                | `skew_cpg`              | `count`       | `overlapping count(CG)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 116 | `cpg_density`                              | `skew_cpg`              | `fraction`    | `count(CG) / valid canonical dinucleotide positions`                                       | DNAKit descriptor_schema_v1                                        |
| 117 | `cpg_observed_expected`                    | `skew_cpg`              | `ratio`       | `count(CG)*canonical_base_count/(count(C)*count(G))`                                       | DNAKit descriptor_schema_v1                                        |
| 118 | `gpc_count`                                | `skew_cpg`              | `count`       | `overlapping count(GC)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 119 | `gpc_density`                              | `skew_cpg`              | `fraction`    | `count(GC) / valid canonical dinucleotide positions`                                       | DNAKit descriptor_schema_v1                                        |
| 120 | `cpg_gpc_ratio`                            | `skew_cpg`              | `ratio`       | `count(CG) / count(GC)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 121 | `cumulative_gc_skew_max`                   | `skew_cpg`              | `count`       | `max prefix cumulative score where G=+1,C=-1,A/T=0`                                        | DNAKit descriptor_schema_v1                                        |
| 122 | `cumulative_gc_skew_min`                   | `skew_cpg`              | `count`       | `min prefix cumulative score where G=+1,C=-1,A/T=0`                                        | DNAKit descriptor_schema_v1                                        |
| 123 | `cumulative_gc_skew_range`                 | `skew_cpg`              | `count`       | `cumulative_gc_skew_max - cumulative_gc_skew_min`                                          | DNAKit descriptor_schema_v1                                        |
| 124 | `cumulative_at_skew_max`                   | `skew_cpg`              | `count`       | `max prefix cumulative score where A=+1,T=-1,C/G=0`                                        | DNAKit descriptor_schema_v1                                        |
| 125 | `cumulative_at_skew_min`                   | `skew_cpg`              | `count`       | `min prefix cumulative score where A=+1,T=-1,C/G=0`                                        | DNAKit descriptor_schema_v1                                        |
| 126 | `cumulative_at_skew_range`                 | `skew_cpg`              | `count`       | `cumulative_at_skew_max - cumulative_at_skew_min`                                          | DNAKit descriptor_schema_v1                                        |
| 127 | `dinucleotide_rc_total_variation`          | `skew_cpg`              | `fraction`    | `0.5*sum_xy(abs(f_xy-f_reverse_complement(xy)))`                                           | DNAKit descriptor_schema_v1                                        |
| 128 | `mono_chargaff_l1_distance`                | `skew_cpg`              | `fraction`    | `abs(f_A-f_T)+abs(f_C-f_G)`                                                                | DNAKit descriptor_schema_v1                                        |
| 129 | `shannon_entropy_k1_bits`                  | `complexity`            | `bits`        | `-sum(p(1-mer)*log2(p(1-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 130 | `shannon_entropy_k2_bits`                  | `complexity`            | `bits`        | `-sum(p(2-mer)*log2(p(2-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 131 | `shannon_entropy_k3_bits`                  | `complexity`            | `bits`        | `-sum(p(3-mer)*log2(p(3-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 132 | `normalized_entropy_k1`                    | `complexity`            | `fraction`    | `shannon_entropy_k1_bits / log2(4**1)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 133 | `normalized_entropy_k2`                    | `complexity`            | `fraction`    | `shannon_entropy_k2_bits / log2(4**2)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 134 | `normalized_entropy_k3`                    | `complexity`            | `fraction`    | `shannon_entropy_k3_bits / log2(4**3)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 135 | `linguistic_complexity_k2`                 | `complexity`            | `fraction`    | `unique 2-mers / min(4**2, valid 2-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 136 | `linguistic_complexity_k3`                 | `complexity`            | `fraction`    | `unique 3-mers / min(4**3, valid 3-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 137 | `linguistic_complexity_k4`                 | `complexity`            | `fraction`    | `unique 4-mers / min(4**4, valid 4-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 138 | `linguistic_complexity_k5`                 | `complexity`            | `fraction`    | `unique 5-mers / min(4**5, valid 5-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 139 | `linguistic_complexity_k6`                 | `complexity`            | `fraction`    | `unique 6-mers / min(4**6, valid 6-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 140 | `linguistic_complexity_product_k1_k6`      | `complexity`            | `fraction`    | `product of defined linguistic_complexity_k values for k=1..6`                             | Observed/possible k-word linguistic complexity                     |
| 141 | `lz76_complexity`                          | `complexity`            | `count`       | `number of phrases in exhaustive Lempel-Ziv 1976 parsing`                                  | Lempel and Ziv 1976; DOI 10.1109/TIT.1976.1055501                  |
| 142 | `normalized_lz76_complexity`               | `complexity`            | `ratio`       | `lz76_complexity*log_base4(canonical_base_count)/canonical_base_count`                     | Lempel and Ziv 1976; DOI 10.1109/TIT.1976.1055501                  |
| 143 | `longest_homopolymer_nt`                   | `complexity`            | `nt`          | `max canonical homopolymer run`                                                            | DNAKit descriptor_schema_v1                                        |
| 144 | `longest_homopolymer_a_nt`                 | `complexity`            | `nt`          | `max homopolymer run of A`                                                                 | DNAKit descriptor_schema_v1                                        |
| 145 | `longest_homopolymer_c_nt`                 | `complexity`            | `nt`          | `max homopolymer run of C`                                                                 | DNAKit descriptor_schema_v1                                        |
| 146 | `longest_homopolymer_g_nt`                 | `complexity`            | `nt`          | `max homopolymer run of G`                                                                 | DNAKit descriptor_schema_v1                                        |
| 147 | `longest_homopolymer_t_nt`                 | `complexity`            | `nt`          | `max homopolymer run of T`                                                                 | DNAKit descriptor_schema_v1                                        |
| 148 | `exact_tandem_repeat_coverage_fraction`    | `complexity`            | `fraction`    | `union bases covered by exact tandem repeats / canonical_base_count`                       | DNAKit exact tandem repeat scanner; units 1..20; minimum repeats 2 |
| 149 | `frame0_codon_count`                       | `coding`                | `count`       | `number of valid forward frame-0 codons`                                                   | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 150 | `frame0_unique_codon_count`                | `coding`                | `count`       | `number of distinct forward frame-0 codons`                                                | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 151 | `frame0_start_codon_count`                 | `coding`                | `count`       | `count(ATG) in forward frame 0`                                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 152 | `frame0_stop_codon_count`                  | `coding`                | `count`       | `count(TAA,TAG,TGA) in forward frame 0`                                                    | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 153 | `frame0_start_codon_fraction`              | `coding`                | `fraction`    | `frame0_start_codon_count / frame0_codon_count`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 154 | `frame0_stop_codon_fraction`               | `coding`                | `fraction`    | `frame0_stop_codon_count / frame0_codon_count`                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 155 | `frame0_codon_entropy_bits`                | `coding`                | `bits`        | `-sum(frame0 codon frequency*log2(frequency))`                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 156 | `frame0_effective_number_of_codons`        | `coding`                | `count`       | `2**frame0_codon_entropy_bits`                                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 157 | `frame0_gc1_fraction`                      | `coding`                | `fraction`    | `GC bases at position 1 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 158 | `frame0_gc2_fraction`                      | `coding`                | `fraction`    | `GC bases at position 2 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 159 | `frame0_gc3_fraction`                      | `coding`                | `fraction`    | `GC bases at position 3 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 160 | `six_frame_complete_orf_count`             | `coding`                | `count`       | `complete start-to-next-stop ORFs across three frames on both strands`                     | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 161 | `six_frame_forward_complete_orf_count`     | `coding`                | `count`       | `complete ORFs across three forward frames`                                                | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 162 | `six_frame_reverse_complete_orf_count`     | `coding`                | `count`       | `complete ORFs across three reverse-complement frames`                                     | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 163 | `six_frame_longest_complete_orf_nt`        | `coding`                | `nt`          | `maximum complete six-frame ORF length including terminal stop`                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 164 | `six_frame_complete_orf_coverage_fraction` | `coding`                | `fraction`    | `union of symbol positions covered by complete six-frame ORFs / canonical_base_count`      | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 165 | `mw_ss_oh_da`                              | `physicochemical`       | `Da`          | `anhydrous mass of one ssDNA strand with 5-prime OH`                                       | DNAKit anhydrous DNA residue mass table v1                         |
| 166 | `mw_ss_5p_phosphate_da`                    | `physicochemical`       | `Da`          | `anhydrous mass of one ssDNA strand with 5-prime phosphate`                                | DNAKit anhydrous DNA residue mass table v1                         |
| 167 | `mw_ds_oh_da`                              | `physicochemical`       | `Da`          | `anhydrous mass of sequence plus complete reverse complement with 5-prime OH`              | DNAKit anhydrous DNA residue mass table v1                         |
| 168 | `mw_ds_5p_phosphate_da`                    | `physicochemical`       | `Da`          | `anhydrous mass of sequence plus complete reverse complement; both 5-prime phosphorylated` | DNAKit anhydrous DNA residue mass table v1                         |
| 169 | `epsilon260_ss_m_inverse_cm_inverse`       | `physicochemical`       | `M^-1 cm^-1`  | `nearest-neighbor epsilon260 pair sum minus internal-base sum`                             | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 170 | `nmol_per_a260_1ml_1cm`                    | `physicochemical`       | `nmol`        | `1e6 / epsilon260 for A260=1, volume=1 mL, path=1 cm`                                      | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 171 | `ug_per_a260_1ml_1cm`                      | `physicochemical`       | `ug`          | `1000*mw_ss_oh_da/epsilon260 for A260=1, volume=1 mL, path=1 cm`                           | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 172 | `tm_wallace_c`                             | `physicochemical`       | `degree C`    | `2*(A+T)+4*(G+C)`                                                                          | Wallace short-oligo 2AT+4GC rule                                   |
| 173 | `stacking_delta_h_kcal_per_mol`            | `physicochemical`       | `kcal/mol`    | `sum SantaLucia nearest-neighbor stacking delta H`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 174 | `stacking_delta_s_cal_per_k_mol_k`         | `physicochemical`       | `cal/(K mol)` | `sum SantaLucia nearest-neighbor stacking delta S`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 175 | `stacking_delta_g37_kcal_per_mol`          | `physicochemical`       | `kcal/mol`    | `stacking delta H - 310.15*stacking delta S/1000`                                          | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 176 | `nn_delta_h_kcal_per_mol`                  | `physicochemical`       | `kcal/mol`    | `SantaLucia complete-duplex delta H with initiation and symmetry`                          | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 177 | `nn_delta_s_cal_per_mol_k`                 | `physicochemical`       | `cal/(K mol)` | `SantaLucia complete-duplex delta S with initiation, symmetry, and salt`                   | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 178 | `nn_delta_g37_kcal_per_mol`                | `physicochemical`       | `kcal/mol`    | `complete-duplex delta H - 310.15*delta S/1000`                                            | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 179 | `nn_tm_c`                                  | `physicochemical`       | `degree C`    | `SantaLucia concentration- and sodium-adjusted Tm`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 180 | `self_complementary`                       | `physicochemical`       | `boolean`     | `sequence == reverse_complement(sequence)`                                                 | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 181 | `diprodb_twist_mean`                       | `dinucleotide_property` | `degree`      | `population mean of Twist values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 182 | `diprodb_twist_sd`                         | `dinucleotide_property` | `degree`      | `population sd of Twist values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 183 | `diprodb_twist_min`                        | `dinucleotide_property` | `degree`      | `population min of Twist values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 184 | `diprodb_twist_max`                        | `dinucleotide_property` | `degree`      | `population max of Twist values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 185 | `diprodb_tilt_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Tilt values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 186 | `diprodb_tilt_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Tilt values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 187 | `diprodb_tilt_min`                         | `dinucleotide_property` | `degree`      | `population min of Tilt values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 188 | `diprodb_tilt_max`                         | `dinucleotide_property` | `degree`      | `population max of Tilt values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 189 | `diprodb_roll_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Roll values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 190 | `diprodb_roll_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Roll values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 191 | `diprodb_roll_min`                         | `dinucleotide_property` | `degree`      | `population min of Roll values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 192 | `diprodb_roll_max`                         | `dinucleotide_property` | `degree`      | `population max of Roll values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 193 | `diprodb_shift_mean`                       | `dinucleotide_property` | `angstrom`    | `population mean of Shift values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 194 | `diprodb_shift_sd`                         | `dinucleotide_property` | `angstrom`    | `population sd of Shift values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 195 | `diprodb_shift_min`                        | `dinucleotide_property` | `angstrom`    | `population min of Shift values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 196 | `diprodb_shift_max`                        | `dinucleotide_property` | `angstrom`    | `population max of Shift values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 197 | `diprodb_slide_mean`                       | `dinucleotide_property` | `angstrom`    | `population mean of Slide values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 198 | `diprodb_slide_sd`                         | `dinucleotide_property` | `angstrom`    | `population sd of Slide values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 199 | `diprodb_slide_min`                        | `dinucleotide_property` | `angstrom`    | `population min of Slide values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 200 | `diprodb_slide_max`                        | `dinucleotide_property` | `angstrom`    | `population max of Slide values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 201 | `diprodb_rise_mean`                        | `dinucleotide_property` | `angstrom`    | `population mean of Rise values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 202 | `diprodb_rise_sd`                          | `dinucleotide_property` | `angstrom`    | `population sd of Rise values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 203 | `diprodb_rise_min`                         | `dinucleotide_property` | `angstrom`    | `population min of Rise values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 204 | `diprodb_rise_max`                         | `dinucleotide_property` | `angstrom`    | `population max of Rise values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 205 | `diprodb_bend_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Bend values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 206 | `diprodb_bend_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Bend values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 207 | `diprodb_bend_min`                         | `dinucleotide_property` | `degree`      | `population min of Bend values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 208 | `diprodb_bend_max`                         | `dinucleotide_property` | `degree`      | `population max of Bend values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 209 | `diprodb_inclination_mean`                 | `dinucleotide_property` | `degree`      | `population mean of Inclination values over valid overlapping dinucleotides`               | Caller-supplied table; DNAKit bundles no numerical values          |
| 210 | `diprodb_inclination_sd`                   | `dinucleotide_property` | `degree`      | `population sd of Inclination values over valid overlapping dinucleotides`                 | Caller-supplied table; DNAKit bundles no numerical values          |
| 211 | `diprodb_inclination_min`                  | `dinucleotide_property` | `degree`      | `population min of Inclination values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |
| 212 | `diprodb_inclination_max`                  | `dinucleotide_property` | `degree`      | `population max of Inclination values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |
| 213 | `diprodb_direction_mean`                   | `dinucleotide_property` | `degree`      | `population mean of Direction values over valid overlapping dinucleotides`                 | Caller-supplied table; DNAKit bundles no numerical values          |
| 214 | `diprodb_direction_sd`                     | `dinucleotide_property` | `degree`      | `population sd of Direction values over valid overlapping dinucleotides`                   | Caller-supplied table; DNAKit bundles no numerical values          |
| 215 | `diprodb_direction_min`                    | `dinucleotide_property` | `degree`      | `population min of Direction values over valid overlapping dinucleotides`                  | Caller-supplied table; DNAKit bundles no numerical values          |
| 216 | `diprodb_direction_max`                    | `dinucleotide_property` | `degree`      | `population max of Direction values over valid overlapping dinucleotides`                  | Caller-supplied table; DNAKit bundles no numerical values          |
| 217 | `diprodb_propeller_twist_mean`             | `dinucleotide_property` | `degree`      | `population mean of Propeller twist values over valid overlapping dinucleotides`           | Caller-supplied table; DNAKit bundles no numerical values          |
| 218 | `diprodb_propeller_twist_sd`               | `dinucleotide_property` | `degree`      | `population sd of Propeller twist values over valid overlapping dinucleotides`             | Caller-supplied table; DNAKit bundles no numerical values          |
| 219 | `diprodb_propeller_twist_min`              | `dinucleotide_property` | `degree`      | `population min of Propeller twist values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 220 | `diprodb_propeller_twist_max`              | `dinucleotide_property` | `degree`      | `population max of Propeller twist values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 221 | `diprodb_major_groove_width_mean`          | `dinucleotide_property` | `angstrom`    | `population mean of Major groove width values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 222 | `diprodb_major_groove_width_sd`            | `dinucleotide_property` | `angstrom`    | `population sd of Major groove width values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 223 | `diprodb_major_groove_width_min`           | `dinucleotide_property` | `angstrom`    | `population min of Major groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 224 | `diprodb_major_groove_width_max`           | `dinucleotide_property` | `angstrom`    | `population max of Major groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 225 | `diprodb_minor_groove_width_mean`          | `dinucleotide_property` | `angstrom`    | `population mean of Minor groove width values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 226 | `diprodb_minor_groove_width_sd`            | `dinucleotide_property` | `angstrom`    | `population sd of Minor groove width values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 227 | `diprodb_minor_groove_width_min`           | `dinucleotide_property` | `angstrom`    | `population min of Minor groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 228 | `diprodb_minor_groove_width_max`           | `dinucleotide_property` | `angstrom`    | `population max of Minor groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 229 | `diprodb_persistence_length_mean`          | `dinucleotide_property` | `nanometer`   | `population mean of Persistence length values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 230 | `diprodb_persistence_length_sd`            | `dinucleotide_property` | `nanometer`   | `population sd of Persistence length values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 231 | `diprodb_persistence_length_min`           | `dinucleotide_property` | `nanometer`   | `population min of Persistence length values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 232 | `diprodb_persistence_length_max`           | `dinucleotide_property` | `nanometer`   | `population max of Persistence length values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 233 | `diprodb_stacking_energy_mean`             | `dinucleotide_property` | `kcal/mol`    | `population mean of Stacking energy values over valid overlapping dinucleotides`           | Caller-supplied table; DNAKit bundles no numerical values          |
| 234 | `diprodb_stacking_energy_sd`               | `dinucleotide_property` | `kcal/mol`    | `population sd of Stacking energy values over valid overlapping dinucleotides`             | Caller-supplied table; DNAKit bundles no numerical values          |
| 235 | `diprodb_stacking_energy_min`              | `dinucleotide_property` | `kcal/mol`    | `population min of Stacking energy values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 236 | `diprodb_stacking_energy_max`              | `dinucleotide_property` | `kcal/mol`    | `population max of Stacking energy values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 237 | `diprodb_free_energy_mean`                 | `dinucleotide_property` | `kcal/mol`    | `population mean of Free energy values over valid overlapping dinucleotides`               | Caller-supplied table; DNAKit bundles no numerical values          |
| 238 | `diprodb_free_energy_sd`                   | `dinucleotide_property` | `kcal/mol`    | `population sd of Free energy values over valid overlapping dinucleotides`                 | Caller-supplied table; DNAKit bundles no numerical values          |
| 239 | `diprodb_free_energy_min`                  | `dinucleotide_property` | `kcal/mol`    | `population min of Free energy values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |
| 240 | `diprodb_free_energy_max`                  | `dinucleotide_property` | `kcal/mol`    | `population max of Free energy values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |

方法、论文、数据库与网址已统一移至[致谢与主要来源](../acknowledgements.md#methods-and-references)；许可和用户责任见[第三方声明](../acknowledgements.md#third-party-notices)。

## 数据划分：稳定哈希、K-fold 与复现 {#dataset-splits}

### 1) 随机与稳定哈希划分 {#1-data-012}

**作用：** 按目标比例把记录分配到 `train / valid / test` 等命名子集，或用固定 seed 生成 K-fold 划分并选择一折作为 `valid`、其余折作为 `train`。返回各子集和分配清单。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项支持 `method="random"`、`method="hash"` 和 `method="kfold"`。

**输入：** 必填记录；比例划分通过 `ratios` 指定总和为 1 的比例，默认 `train=0.8, valid=0.1, test=0.1`。K-fold 通过 `n_splits` 指定折数、`fold_id` 指定验证折。默认 `seed=42`、`shuffle=True`。`hash` 模式要求每条记录具有稳定且唯一的 `record.id`。

`SplitConfig` 和 `dnakit split` 的默认 seed 已由 `0` 调整为 `42`；需要复现此前省略 seed 的比例划分时，请显式传入 `seed=0`。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("AC"), f"r{i}") for i in range(10)]
result = split(
    records,
    config=SplitConfig(
        method="random", ratios={"train": 0.8, "test": 0.2}, seed=17
    ),
)
print(dict(result.counts))
```

**示例结果：**

```text
{'train': 8, 'test': 2}
```

**顺序无关的 `hash` 划分**

`hash` 模式在 `shuffle=True` 时使用版本化 SHA-256 计算 `seed + record.id` 的稳定排序键，在 `shuffle=False` 时按 `record.id` 的稳定字节序排列，再按比例分配记录。因此，同一批记录即使输入顺序变化，按 `record.id` 查看时分类仍完全一致；它不使用 Python 内置的进程随机化 `hash()`。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("AC"), f"r{i}") for i in range(10)]
config = SplitConfig(
    method="hash",
    ratios={"train": 0.6, "valid": 0.2, "test": 0.2},
    seed=17,
    preserve_order=False,
)
first = split(records, config=config)
second = split(list(reversed(records)), config=config)
first_by_id = {item.record_id: item.split for item in first.assignments}
second_by_id = {item.record_id: item.split for item in second.assignments}
print(first_by_id == second_by_id)
print(first.get("train").ids == second.get("train").ids)
```

**示例结果：**

```text
True
True
```

#### 固定 K-fold 划分：train / valid

设置 `method="kfold"`，行为对应 [scikit-learn KFold](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.KFold.html)：

```python
KFold(n_splits=10, shuffle=True, random_state=42)
```

**同一随机种子 → 同一套固定的 K-fold 划分 → 依次选择 K 个验证折。** 保持输入记录及其顺序、`seed`、`shuffle`、`n_splits` 不变，只改变 `fold_id`。例如 10 折使用 `fold_id=0..9`，每个样本恰好作为验证样本一次；不要使用 `seed + fold_id`。每次调用会重建相同的划分并返回所选折。

| `SplitConfig` 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `method` | `"random"` | 使用 `"kfold"` 启用 K-fold；`"random"` / `"hash"` 保留比例划分 |
| `n_splits` | K-fold 模式下为 `10` | 总折数，可改为 5、10 等；要求 `2 <= n_splits <= 样本数` |
| `fold_id` | K-fold 模式下为 `0` | 当前验证折，范围为 `0` 到 `n_splits - 1` |
| `seed` | `42` | 生成整套划分的固定随机种子，K-fold 要求 `0 <= seed < 2**32` |
| `shuffle` | `True` | 划分前打乱；`False` 时按输入顺序生成连续的验证折，seed 不影响分配 |
| `ratios` | 比例模式为 `0.8 / 0.1 / 0.1` | K-fold 模式请省略；`train / valid` 目标比例由折数推导 |

K-fold 返回的子集名称固定为 `train` 和 `valid`，每个子集内部保持原始输入顺序，与 scikit-learn 一致；`preserve_order=False` 也不会额外打乱子集。样本数不能整除 K 时，前 `样本数 % K` 个验证折多一个样本，各验证折大小相差至多 1。结果和导出的分配清单记录 `seed`、`shuffle`、`n_splits`、`fold_id` 及每条记录的归属。

此模式需要 scikit-learn；可安装 `python -m pip install -e ".[bio-analysis]"`。已安装 `.[ml]` 或 `.[training]` 的环境也包含所需依赖。普通比例划分不需要这项依赖。

**示例：先留出固定 test，再对 train 做 10-fold。**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("ACGT"), f"r{i}") for i in range(100)]
seed = 42
n_splits = 10  # 改为 5 即可运行 5 折

## 第一步只做一次，留出固定测试集。
outer = split(
    records,
    config=SplitConfig(method="random", ratios={"train": 0.8, "test": 0.2}, seed=seed),
)
training_pool = outer.get("train")
test = outer.get("test")

## 每轮只改变 fold_id；测试集不传入 K-fold。
for fold_id in range(n_splits):
    fold = split(
        training_pool,
        config=SplitConfig(
            method="kfold", n_splits=n_splits, fold_id=fold_id,
            shuffle=True, seed=seed,
        ),
    )
    train = fold.get("train")
    valid = fold.get("valid")
    print(f"fold={fold_id}: train={len(train)}, valid={len(valid)}, test={len(test)}")
    # 在这里创建并训练一个新的模型，再用 valid 评估。
```

上述示例的 `fold=0` 到 `fold=9` 均为 `train=72, valid=8, test=20`，验证样本逐折轮换，20 条测试样本保持固定。若已准备好训练数据，直接将其作为 `training_pool` 即可。`n_splits` 控制总折数，`fold_id` 选择某一折；遍历所有折才形成完整 K-fold 流程，`split()` 本身只划分数据，不会自动训练模型或汇总指标。

接入 `train_custom_property_model()` 时，需要在每一轮把当前 `train`、`valid` 和固定 `test` 写成训练接口支持的 split 列值，并用 `DataMapping.split_field` 映射；不能直接把 `fold_0` 等名称当作训练接口的 split 值。相关说明见[自定义模型训练](../api/features/23_prediction_custom_training.md)。

**CLI：** 输入已留出的训练池，输出当前验证折对应的两个子集和分配清单：

```bash
dnakit split train.fasta runs/fold_0 --method kfold --n-splits 10 --fold-id 0 --seed 42
```

运行其余折时保持输入和 seed 不变，修改 `--fold-id` 和输出目录；可用 `--no-shuffle` 关闭打乱。Agent/MCP 的 `dnakit.datasets.split` 使用相同的 `config` 参数；JSON/YAML workflow 的 `split.params` 也接受 `method`、`n_splits`、`fold_id`，随机种子使用 workflow 的全局 `seed`，请设为 `42`。

## 下载：格式、进度、重试与物种名称 {#downloads}

- **Agent 参数：** `query` 是组装标识或物种名，不是工具搜索词。请求下载人类 hg38 时使用 `query="hg38"`，不要传入 `"human hg38 reference genome FASTA"`。工具 Schema 包含此约束；下载接口的 `api_key` 指 NCBI 密钥，不是 DeepSeek 密钥。


根据链接、物种名称、数据库编号或查询结果下载文件、参考基因组、序列和注释。

NCBI／ENA 序列下载、通用链接下载和 ENCODE 批量下载的参数及示例均在本页。实验查询见[数据查询](../api/features/22_database_query_encode.md#encode)，追踪、引用与来源记录见[数据集整理](../api/features/10_datasets.md#encode-collection)。

### 按 accession 下载序列

`dnakit.download.ncbi.sequence(accession, output_path, *, database="nuccore", format="fasta", seq_start=None, seq_stop=None, strand=1, config=None, search_config=None)`
直接访问 NCBI EFetch；`dnakit.download.ena.sequence(accession, output_path, *, format="fasta", config=None, search_config=None)` 访问 ENA Browser。这两个接口无需 ToolUniverse。

```python
import dnakit

result = dnakit.download.ncbi.sequence(
    "NC_000913.3", "ecoli_k12.fasta", format="fasta",
)
print(result.files[0].path)
print(result.metadata["resolved_accession"], result.metadata["sequence_length"])
print(result.manifest_path)
```

| 提供方                     | 基础包直接获取                    | 需要`dnakit[bio-formats]` 的格式转换 |
| -------------------------- | --------------------------------- | -------------------------------------- |
| NCBI`nuccore`（DNA/RNA） | `fasta`、`genbank`            | `embl`：通过 Biopython 转换 GenBank  |
| NCBI`protein`            | `fasta`、`genbank`（GenPept） | 不支持 EMBL 蛋白输出                   |
| ENA INSDC 核酸记录         | `fasta`、`embl`               | `genbank`：通过 Biopython 转换 EMBL  |

NCBI 支持 RefSeq（包括 `NZ_`）和 INSDC 序列 accession；ENA 需要 INSDC accession，例如
`U00096.3`。格式转换不代表 accession 已转成另一个数据库的编号。`GCA_/GCF_` 是组装编号，
请用已有的 `genome_package()`。

返回 `DatasetDownloadResult`，保存序列和相邻的 `.manifest.json`，记录请求及解析出的
accession、序列长度、来源 URL、下载/输出格式、转换标记和校验值。正式写入前校验记录的
accession/version、格式和序列长度；异常响应不会作为序列文件保存。这里的校验值用于本地完整性
记录，并非提供方独立签名。已有文件默认受保护；需要覆盖时传入 `DownloadConfig(overwrite=True)`。

可选 `seq_start/seq_stop` 是 **1-based 闭区间**，必须同时指定；`strand=2` 获取反向互补序列，
只适用于 NCBI 核酸序列。蛋白下载目前获取完整记录。此处与 Ensembl 坐标下载的 0-based 半开区间不同。

`DownloadConfig` 管理保存、覆盖和文件大小上限；`SearchConfig` 管理请求超时、响应大小
（默认 20 MB、最高 100 MB）及 NCBI `api_key/email`。显式提供 `search_config` 时以它的
凭据为准，否则使用 `DownloadConfig.api_key`；持久化 URL 会隐藏密钥和邮箱。
较大的组装请使用分包下载。`format="embl"` 的可选依赖可通过
`python -m pip install 'dnakit[bio-formats]'` 一次安装，数据本身仍从上游数据库获取。

MCP 对应名称为 `dnakit.download.ncbi_sequence` 和 `dnakit.download.ena_sequence`；
Python 也支持这两个扁平别名。原有 `dnakit.download.sequence` 继续用于 Ensembl 坐标序列。

#### 下载进度

下载默认向终端标准错误流显示开始提示、进度条和完成路径。参考基因组还显示解压进度；
服务器未提供文件大小时显示已下载字节数，不估算百分比。只有下载、校验和文件写入成功后才提示完成。
`download_genome(..., show_progress=False)` 可关闭显示；`dnakit.download` 中接受 `config` 的网络下载接口
使用 `DownloadConfig(show_progress=False)`。显式传入 `progress`（序列下载为 `query_progress`）回调时，
由回调负责显示，避免重复进度条。NCBI 数据包、ENA、ENCODE、UCSC 和通用文件下载共享这一行为。

#### 下载失败与重试

基因组包下载、通用 HTTPS 文件下载及共享远程 GET 查询遇到超时、连接中断或 HTTP
408/429/500/502/503/504 时，最多尝试 3 次，重试前分别等待 0.25、0.5 秒。
404、参数错误、校验和不匹配和无效 ZIP 不自动重试；POST 查询也不自动重发。
每次失败的下载从暂存文件重新开始，尚不支持程序退出后的 HTTP Range 断点续传。

服务器提供 `Content-Length` 时会核对实际字节数；基因组包还核对 FASTA 校验和及
报告中的 assembly accession，不能把其他版本误报为当前请求成功。
只有完整下载和校验通过后才安装正式文件。覆盖冲突保留 `FileExistsError`，目录权限错误保留
相应文件系统异常；损坏输入不会替换已有完整结果。瞬时故障重试耗尽后，异常的
`context["attempts"]` 为 3，并提示检查网络后重新运行。

**常见物种名称**

本节列出可以直接复制到 `download_genome(query, ...)` 或
`dnakit download-genome` 命令中的常见 NCBI 物种名称。

- **API**：`dnakit.references.download_genome(query[必须], output_dir[必须], overwrite[可选], keep_package[可选], timeout[可选], api_key[可选], api_base_url[可选], chunk_size[可选], max_download_bytes[可选], progress[可选])`、`dnakit download-genome query[必须] output_dir[必须] --overwrite[可选] --keep-package[可选] --api-key[可选] --progress/--no-progress[可选]`。

**人和常用模式动物**

| 中文名称 | 可直接复制的`query`      | 当前参考组装 accession |
| -------- | -------------------------- | ---------------------- |
| 人类     | `Homo sapiens`           | `GCF_000001405.40`   |
| 小鼠     | `Mus musculus`           | `GCF_000001635.27`   |
| 大鼠     | `Rattus norvegicus`      | `GCF_036323735.1`    |
| 猪       | `Sus scrofa`             | `GCF_054392235.1`    |
| 牛       | `Bos taurus`             | `GCF_002263795.3`    |
| 羊       | `Ovis aries`             | `GCF_016772045.2`    |
| 山羊     | `Capra hircus`           | `GCF_001704415.2`    |
| 马       | `Equus caballus`         | `GCF_041296265.1`    |
| 狗       | `Canis lupus familiaris` | `GCF_011100685.1`    |
| 猕猴     | `Macaca mulatta`         | `GCF_049350105.2`    |
| 鸡       | `Gallus gallus`          | `GCF_016699485.2`    |
| 斑马鱼   | `Danio rerio`            | `GCF_049306965.2`    |
| 爪蟾     | `Xenopus tropicalis`     | `GCF_000004195.4`    |

**常见农作物和植物**

| 中文名称     | 可直接复制的`query`           | 当前参考组装 accession |
| ------------ | ------------------------------- | ---------------------- |
| 拟南芥       | `Arabidopsis thaliana`        | `GCF_000001735.4`    |
| 水稻         | `Oryza sativa Japonica Group` | `GCF_034140825.1`    |
| 玉米         | `Zea mays`                    | `GCF_902167145.1`    |
| 小麦         | `Triticum aestivum`           | `GCF_018294505.1`    |
| 大豆         | `Glycine max`                 | `GCF_000004515.6`    |
| 番茄         | `Solanum lycopersicum`        | `GCF_036512215.1`    |
| 马铃薯       | `Solanum tuberosum`           | `GCF_000226075.1`    |
| 陆地棉       | `Gossypium hirsutum`          | `GCF_007990345.1`    |
| 烟草         | `Nicotiana tabacum`           | `GCF_000715075.1`    |
| 毛果杨       | `Populus trichocarpa`         | `GCF_000002775.5`    |
| 阿拉比卡咖啡 | `Coffea arabica`              | `GCF_036785885.1`    |

**其他常用模式生物**

| 中文名称     | 可直接复制的`query`              | 当前参考组装 accession |
| ------------ | ---------------------------------- | ---------------------- |
| 果蝇         | `Drosophila melanogaster`        | `GCF_000001215.4`    |
| 秀丽隐杆线虫 | `Caenorhabditis elegans`         | `GCF_000002985.6`    |
| 家蚕         | `Bombyx mori`                    | `GCF_030269925.1`    |
| 蜜蜂         | `Apis mellifera`                 | `GCF_003254395.2`    |
| 酿酒酵母     | `Saccharomyces cerevisiae S288C` | `GCF_000146045.2`    |


### 22) 通用链接下载 {#22-dbd-022}

`download_url(url, output_dir)` 只需提供 HTTPS 文件直链和目标文件夹，目录不存在时会自动创建。
默认使用链接路径末尾的文件名，自动解码 URL 转义字符，不把 `?` 后的查询参数放入文件名。

```python
from dnakit.download import download_url

result = download_url(
    "https://example.org/files/data.fa.gz",  # 替换为需要下载的文件直链
    "downloads/my_data",
)
print(result.path)       # 下载文件的绝对路径
print(result.byte_size)  # 文件字节数
```

可选参数 `filename="reference.fa.gz"` 可以指定保存名称；链接没有文件名时必须提供该参数。
默认同名文件已存在就抛出 `FileExistsError`；需要覆盖时使用
`config=DownloadConfig(overwrite=True)`（从 `dnakit.download` 导入 `DownloadConfig`）。
超时、单文件大小上限和进度显示也通过 `DownloadConfig` 设置。
可传入 `expected_md5` / `expected_sha256` 验证下载内容，或用 `headers` 提供请求头。

此接口不限文件类型，按响应原样保存；例如 `.gz` 下载后仍是压缩文件。
它接收文件直链，不解析网页或网盘分享页面，不自动解压。
Agent/MCP 可通过 `dnakit.download.download_url` 调用，必填参数同样只有 `url` 和 `output_dir`。
需要批量下载并生成来源清单时，可以使用 `dataset([RemoteFile(url1), RemoteFile(url2)], output_dir)`。


### 23) 从检索结果下载序列 {#downloads-从检索结果下载序列}

本例依次使用原有且已扩展的 `entrez`、新增的 `accessions` 和新增的 NCBI 序列下载函数。查询会联网，下载会写入当前目录；实际命中可能变化。

```python
import dnakit

found = dnakit.search.ncbi.entrez(
    "nuccore", organism="Escherichia coli", strain="K-12",
    seq_type="complete_genome", retmax=3,
)
accessions = dnakit.search.ncbi.accessions(found.metadata["uids"])
if accessions:
    downloaded = dnakit.download.ncbi.sequence(
        accessions[0], "ecoli_k12.fasta", format="fasta",
    )
    print(downloaded.files[0].path)
    print(downloaded.metadata["sequence_length"])
    print(downloaded.manifest_path)
else:
    print("没有匹配序列，请核对查询条件。")
```

ENA 使用 INSDC 编号，例如 `dnakit.download.ena.sequence("U00096.3", "ecoli_ena.fasta")`。这与数据库模块的 `ncbi_get_sequence` 不同：下载命名空间把文件和来源清单保存到磁盘，数据库包装按对应服务返回记录。格式转换依赖与坐标说明见[下载教程](../api/features/15_download.md#accession-sequences)。

## 基因注释：安装、坐标、参考组装与可靠性 {#gene-annotation}

`dnakit.annotation.annotate_genes()` 输入一段基因组 DNA 和物种，返回具体基因名称、ID、方向与输入序列中的区间。`assembly` 可选：首次未填写时，查询该物种当前的 NCBI RefSeq 参考组装；成功后把准确的组装和注释文件固定在项目缓存中。

这是基于参考序列比对的已有注释映射。所选组装必须有基因注释；无法获得正式名称时保留基因 ID 或 locus tag，不推测名称。本接口面向基因组 DNA，不用于拼接后的 cDNA、转录本或从头基因预测。

### 安装和命令行

需要安装 DNAKit 和外部程序 minimap2。导入模块不会下载参考文件或执行外部程序。

安装当前工作区的 DNAKit：

```bash
python -m pip install -e .
```

通过 Conda 安装 minimap2：

```bash
conda install -c conda-forge -c bioconda minimap2
```

仅提供物种，自动选择并缓存默认组装；`query.fasta` 应含一条基因组 DNA 序列：

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json
```

明确指定组装版本：

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --assembly GCF_000001405.40 \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json
```

使用已准备的同版本 FASTA 和 GFF3，完全离线运行：

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --reference-fasta reference.fna \
  --annotation-gff annotation.gff3 \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json \
  --no-download
```

可通过 `--minimap2 /path/to/minimap2` 指定可执行文件，`--threads 4` 设置线程数。已有输出默认不覆盖，需要覆盖时显式增加 `--overwrite`。输入 FASTA 和本地参考文件支持 gzip；一次命令处理一条序列，多个样本可通过 Python 循环调用并复用同一缓存。

### Python 接口

```python
from dnakit.annotation import GeneAnnotationConfig, annotate_genes
from dnakit.io import read

with read("query.fasta", format="fasta") as source:
    record = next(iter(source))

result = annotate_genes(
    record.sequence.symbols,
    species="Homo sapiens",  # 也接受 NCBI TaxID，例如 9606
    assembly=None,
    sequence_id=record.id,
    config=GeneAnnotationConfig(
        cache_dir=".dnakit/gene_annotation",
        threads=4,
    ),
)
print(result.reference["assembly"])
for gene in result.gene_records:
    print(gene.to_dict())
for segment in result.segments:
    print(segment.to_dict())
```

`assembly` 接受版本化 `GCA_/GCF_` 编号、现有固定别名，以及能在指定物种下唯一解析的组装名称。如果名称同时匹配 RefSeq/GenBank 两套记录，会提示改用明确编号。显式指定旧版时不自动替换为新版。

### 输出和坐标

- `reference`：实际组装编号、组装名、注释元数据、文件 SHA256、来源和选择方式。自动下载同时验证物种与组装的归属；菌株可匹配其所属物种。
- `gene_records`：一个基因实例一条记录，允许重叠。`gene_name` 缺失时为 `null`，`display_name` 回退到 locus tag 或 ID。重复位点的映射身份由 `mapping_id` 保留。
- `segments`：完整覆盖输入序列、彼此不重叠的分段；类型为 `gene`、`intergenic` 或 `unknown`。重叠区间保留多个基因。
- `alignments`：比对位置、方向、MAPQ、相似度和是否达到接受阈值。

Python 对象中的 `Interval` 为 **0-based 半开区间**；`to_dict()` 和 CLI 的 JSON 均为 **1-based 闭区间**。例如输入位置 1–12000，在 Python 内部是 `Interval(0, 12000)`。

对示意基因 A:1–12000、B:18001–45000、C:43000–52000，`gene_records` 保留三条记录；`segments` 为：

```text
1–12000      gene         A
12001–18000  intergenic
18001–42999  gene         B
43000–45000  gene         B + C
45001–52000  gene         C
```

`gene_records.start/end` 是输入中的外包区间，`parts` 才是可靠映射到该基因的具体片段；存在插入或歧义时，中间可能有空洞。正反链和插入/缺失通过 CIGAR 逐块映射，不进行简单的整体坐标平移。GFF3 明确声明环状染色体时，支持跨原点的基因区间。

### 默认组装与可靠性

未提供 `assembly` 时，查询该物种及其菌株下的当前 RefSeq 参考组装。如果找不到参考组装、存在多个候选、物种不匹配或没有 GFF3 基因注释，会返回明确错误；不会随意挑选其他组装。尚未收录的物种或组装可以提供本地文件。

同一缓存中的相同请求会复用第一次成功获取的组装和注释，不再联网解析默认值。参考 FASTA、GFF3 与 minimap2 索引都校验 SHA256；索引按参考内容、preset 和 minimap2 版本隔离。要更新默认组装或注释，请使用新的 `cache_dir`。`--no-download` 可复用已固定的缓存，或者使用一对本地参考文件。

本地文件检查序列 ID、长度、基因边界及可获得的组装声明。离线时物种由调用者声明，结果会记录 `species_validated=false`；FASTA/GFF3 的物种和版本配套关系仍由提供者负责。GFF3 必须包含 gene、pseudogene 或免疫基因片段这类基因级行，不能只含 CDS/转录本。染色体 ID 必须完全匹配，不自动猜测 `chr1` 与 accession 的对应关系。

默认 minimap2 preset 为 `asm5`，适合较长且准确的基因组 DNA；可根据数据修改 `preset`。默认只接受主比对、MAPQ ≥ 20、相似度 ≥ 0.90、实际比对长度 ≥ 100 bp；MAPQ=255 视为未知。阈值是可配置筛选条件，不代表经过校准的基因识别准确率。

**`intergenic` 仅表示可靠映射位置不被本次所选 GFF3 的基因特征覆盖。** 未映射、低置信度、多条可信主比对重叠、输入中的插入，以及参考 contig 没有注释覆盖信息时均输出 `unknown`。缺少注释不等于生物学上不存在基因。

首次下载大型基因组和建立索引可能需要较多磁盘、内存和时间。配置可限制参考长度、下载量、GFF 行数、索引/比对输出大小和运行时间；单批索引避免多批索引导致 MAPQ 不可靠。

### 验证与来源

运行坐标映射、默认组装解析、缓存完整性和 CLI 测试；安装 minimap2 后还会运行真实比对测试：

```bash
PYTHONPATH=src python -m pytest -q others/tests/unit/test_gene_annotation.py
```

参考：[NCBI 参考组装选择规则](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genome-processing/reference-selection/)、[NCBI 组装数据包](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/)、[NCBI GFF3 说明](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/)、[minimap2 手册](https://lh3.github.io/minimap2/minimap2.html)。

本版新增：多物种、可选 assembly 的参考基因注释流程、重叠区间输出及可复现缓存。

## 可选后端的依赖、返回值与使用边界 {#optional-bioinformatics}

`annotation`、`comparative`、`molbio` 入口返回 `ProviderResult`；[数据库变异注释](../api/features/35_variant_annotation.md#variant-annotation)、[选择压力计算](../api/features/30_phylogeny.md#selection)和[富集接口](../api/features/37_functional_enrichment.md#enrichment)使用 `dnakit.databases` 时，保留各自的返回结构。依赖和坐标约定以具体函数说明为准，数据库调用方式见[查询教程](../api/features/26_gene_annotation.md#dnakit.databases.ensembl_lookup_gene)。

以下说明适用于通过可选科学后端调用的功能。公开 API 按功能放在
`annotation`、`comparative` 和 `molbio` 中，
不需要使用 `integrations.tooluniverse` 一类入口。

安装可选依赖：

```bash
python -m pip install "dnakit[external-tools]"
```

这些方法不会训练或微调模型。Ensembl、ClinVar、dbSNP 和 gnomAD 需要联网；dN/dS 与
Golden Gate 在本地执行。每次调用只加载对应的一个
后端工具。统一返回 `ProviderResult`：主要结果位于 `result.data`，调用参数、后端版本和
来源分别位于 `result.parameters`、`result.metadata` 和 `result.provenance`。


### 使用边界

- 未安装可选依赖时会抛出 `BackendUnavailableError`，并提示安装 `dnakit[external-tools]`。
- 远程提供方可能限流、维护或改变数据库版本；请保存 `result.provenance`。
- ClinVar、VEP、dbSNP 和 gnomAD 输出仅用于研究注释，不能直接作为临床诊断。
- 本页列出的功能使用 DNAKit 随包提供的封装，不再导入 ToolUniverse。更多接口见[统一 Skill 接口](../llm_tools/skills.md#unified-skill-api)。

## 性质预测：安装、权重缓存、训练参数与数据分享 {#property-prediction}

本页包含自定义训练与预测流程，以及 AlphaGenome、ESM、Evo 2、NVIDIA NIM 等已训练模型接口的参数、返回值和凭据要求。调用示例分别见[自定义模型训练](../api/features/23_prediction_custom_training.md)、[调控轨道预测](../api/features/23_prediction_regulatory_tracks.md)、[变异效应](../api/features/23_prediction_variant_effects.md)和[蛋白质特征与结构预测](../api/features/23_prediction_protein_features.md)。CADD、AlphaMissense 等已收录的分数见[变异效应评分与实验数据](../api/features/36_variant_effect_scores.md)；SpliceAI 与 Pangolin 见[变异剪接效应](../api/features/23_prediction_variant_effects.md#variant-splicing)。

缺少或无法导入预测依赖时，DNAKit 会先向终端标准错误输出安装提示，再抛出原有异常；异常中也保留安装命令。`verbose=False` 关闭运行进度，依赖错误提示仍保留。DNAKit 不会自动执行安装。

在 **DNAKit 项目根目录**，按所用功能安装：

| 功能 | 安装命令 |
| --- | --- |
| Enformer（含 27 个分类任务） | `python -m pip install -e ".[neural,neural-enformer]"` |
| 通用 PyTorch / Transformers 主干 | `python -m pip install -e ".[neural]"` |
| 自定义深度学习训练 | `python -m pip install -e ".[training]"` |
| 基础机器学习 | `python -m pip install -e ".[ml]"` |
| XGBoost / LightGBM / CatBoost | 分别使用 `.[xgboost]` / `.[lightgbm]` / `.[catboost]` |

通过发布包安装的用户，可将 `-e ".[选项]"` 改成 `"dnakit[选项]"`，无需源码目录。Evo 2、Caduceus 及需要官方源码的模型会显示各自的安装提示；Evo 2 的 extra 要求 Python 3.11–3.12，AlphaGenome、LucaOneTasks、JanusDNA 仍需满足官方环境要求。

多个模型可以共用满足各自版本约束的 PyTorch 或 TensorFlow；二者不能互相替代。上述命令不强制升级或重新安装已满足要求的包，只补齐缺失或不兼容的依赖。当前 `neural` / `training` 的 PyTorch 范围是 `>=2.5,<3`，没有为每个模型固定不同版本；第三方模型仍可能有更严格的约束。参见 [pip 安装行为](https://pip.pypa.io/en/stable/cli/pip_install/)。

本页包含三条路径：

1. **用户自定义性质：** 从用户文件训练分类或回归模型，支持全量微调和线性探测；
2. **官方现成性质：** 加载官方已经发布的任务头直接推理，不在本地重新训练。
3. **机器学习性质：** 使用 scikit-learn 或可选的 XGBoost、LightGBM、CatBoost，训练自己的分类或回归模型。

原有的[神经网络表征](../api/features/08_fingerprints.md)仍只负责提取 embedding。下面先介绍自定义训练；随后编号 1–54 的内容仍表示 54 个官方直接推理功能。

本页当前共 54 个可直接推理的功能。下面每个功能都单独说明作用、预测方法、API、输入与输出。示例结果中的占位符只表示输出结构；具体数值、标签数和轨道数取决于输入、物种、组织筛选条件及实际 checkpoint。

本章涉及的模型论文、官方仓库和权重来源统一收录在 [FAQ：深度学习性质预测的参考文献](../faq.md#deep-learning-property-prediction-references)。

### 输入与输出的异常处理

直接预测拒绝空序列和非法字符，并在序列错误中提供 `record_id`。Enformer 任务允许
1、199、200、201 bp 等长度；默认上限为 196608 bp，显式 `max_length` 可以收紧限制，
不能放宽模型上限。不会静默截断输入；IUPAC 符号按 `ambiguity_policy` 处理。

`PredictionOutput` 拒绝空、非有限、复数或无法转换的数值，以及重复输出列名。
直接分类／分割任务要求输出为 `[0,1]` 概率；Enformer 分类还要求输出维度等于类别数、
概率和在 `1e-5` 容差内等于 1。连续分数、回归和轨道输出不套用概率范围。
自定义后端也必须遵守所选任务的输出约定。

直接预测及 embedding 公共入口将内存不足（包括后端异常链中的 OOM）报告为
`MODEL_OUT_OF_MEMORY`，提示减小批次、释放内存或使用支持的其他设备；不自动输出替代预测。
显式指定不存在的 CUDA 设备或编号会在模型加载前报错。

### 权重缓存与自动下载 {#property-prediction-checkpoint-cache}

预测前会检查本地预训练权重：缺失时只下载当前模型或任务所需的文件，并提示开始、完成及本地绝对路径；已有权重直接复用。`verbose=False` 可关闭提示。Enformer 任务默认存放在 `./ckpt/enformer-benchmarks/`，其他主干存放在 `./ckpt/<模型目录>/`。机器学习使用的 embedding 主干同样按需下载。用户训练的分类器（`estimator.joblib`）和神经网络任务头（`head.npz`）没有登记公共下载来源，需保留 `model_dir` 中的训练输出。LucaOneTasks 使用所提供源码目录中的官方下载流程。

Enformer 的 27 个分类任务无需指定权重路径即可预测。使用 `checkpoint_dir` 可指定缓存目录：已有文件也会在其 `enformer-benchmarks/` 及对应任务类别的 `nt-revised/` 或 `genomic-benchmarks/` 子目录中查找；缺失时直接下载到指定目录。也可使用 `checkpoint_path` 指定精确文件路径，所选 Enformer 任务的文件缺失时会下载到该路径。两个参数不能同时使用。

安装命令为 `python -m pip install "dnakit[neural,neural-enformer]"`，其中已包含 Google Drive 下载所需的 `gdown`。下载先写入临时文件，检查后再放入正式路径；失败或中断的文件不会被当作完整权重。预测时继续校验任务和模型参数，类别含义以 `output.output_names` 中的标签为准。

### 机器学习分类与回归

与深度学习共用 `CustomPropertyTrainingConfig`、`train_custom_property_model()` 和 `predict_custom_properties()`。设置 `finetune_method="machine_learning"` 后，`model` 表示机器学习算法；embedding 主干另用 `embedding_model` 指定。

安装基础机器学习及表格读取依赖（不需要 PyTorch）：

```bash
python -m pip install -e ".[ml]"
```

安装三个可选提升树库：

```bash
python -m pip install -e ".[xgboost,lightgbm,catboost]"
```

如需 DNA 模型 embedding，再安装 `.[neural]` 及所选主干要求的附加依赖。缺少某个可选算法库时只会影响该算法，调用时会给出对应安装命令。

#### 模型覆盖

目前登记 37 类算法。以下是完整模型名称；`available_machine_learning_models()` 返回每类支持的任务和依赖，不会导入或安装模型库。

| 任务范围 | `model` |
| --- | --- |
| 分类和回归 | `ridge`、`sgd`、`svm`、`linear_svm`、`nu_svm`、`knn`、`decision_tree`、`extra_tree`、`random_forest`、`extra_trees`、`gradient_boosting`、`hist_gradient_boosting`、`adaboost`、`bagging`、`gaussian_process`、`mlp`、`dummy`、`xgboost`、`lightgbm`、`catboost` |
| 仅分类 | `logistic_regression`、`perceptron`、`gaussian_nb`、`multinomial_nb`、`bernoulli_nb`、`complement_nb`、`lda`、`qda` |
| 仅回归 | `linear_regression`、`lasso`、`elastic_net`、`bayesian_ridge`、`ard`、`huber`、`ransac`、`theil_sen`、`kernel_ridge` |

`mlp` 是 scikit-learn 的浅层神经网络，`dummy` 是基线。这里的范围是监督性质预测，不包括聚类、异常检测及自动模型搜索。各算法仍有自己的样本量、参数和计算成本约束；例如高斯过程不适合直接用于超大样本集。

#### 三种输入特征

| `feature_source` | 行为 |
| --- | --- |
| `sequence`（默认） | 序列长度、GC 比例、模糊碱基比例，加固定顺序的 k-mer 频率；`kmer_k=3`，允许 1–6，维度为 `3 + 4**kmer_k`。含模糊碱基的窗口跳过，GC 分母只含 A/C/G/T；无有效窗口或碱基时相应频率为 0。无需深度学习权重。 |
| `embedding` | 冻结 `embedding_model`（默认 `dnabert2`）提取向量，再训练机器学习模型。沿用 checkpoint、pooling、batch_size、device 等表征参数。 |
| `numeric` | 使用显式 `DataMapping.feature_fields`，按给出的列顺序读取数值；可设置 `sequence_field=None`，输入文件不需要序列。缺失、无穷、非数值特征会报错，当前不自动填补。 |

数值特征列必须显式指定，不能与标签、ID、split 或序列字段重叠。推理时必须保留训练时的特征列名称与顺序。不要把目标派生列或其他泄漏信息选为特征。

下面使用自己的序列数据训练随机森林回归：

```python
from dnakit.predictions import (
    CustomPropertyTrainingConfig, DataMapping,
    train_custom_property_model, predict_custom_properties,
)

config = CustomPropertyTrainingConfig(
    task_type="regression",  # 分类改为 classification
    finetune_method="machine_learning",
    model="random_forest",
    model_params={"n_estimators": 300, "max_depth": 12},
    feature_source="sequence",
    kmer_k=3,
    mapping=DataMapping(sequence_field="sequence", label_fields=("target",)),
    output_dir="runs/random-forest",
    workers=6,
    seed=42,
)
result = train_custom_property_model("my_labels.csv", config=config)
print(result.metrics["test"])

predictions = predict_custom_properties(
    "new_sequences.csv", model_dir=result.output_dir,
    trust_model=True,  # 仅加载自己训练或来源可信的模型
)
print(predictions.predictions)
```

使用数值特征时改为 `feature_source="numeric"`，并传入例如 `DataMapping(sequence_field=None, feature_fields=("gc", "length"), label_fields=("target",))`；使用 embedding 时改为 `feature_source="embedding", embedding_model="dnabert2"`，并按主干要求设置 `allow_remote_code` 等参数。

#### 训练、评估与保存

默认仍为种子 42、80%/10%/10% 划分；已映射的 split 列优先。标签模式、测试指标和深度学习路径一致。多标签、多输出分类及多目标回归按目标分别训练估计器；训练中恒定的二值标签使用恒定分类器。分类预测取最大概率类别，多标签使用 `threshold`。

`model_params` 直接传给所选估计器构造函数，接受可序列化的参数，不接受任意 Python 估计器对象。`workers` 控制默认 `n_jobs` 和数值库线程数（0 按 1 处理），CatBoost 使用 `thread_count`；模型参数可覆盖估计器的并行数和随机种子。`standardize_features=True` 默认在训练集上拟合标准化器，验证和测试只做变换。`multinomial_nb`、`complement_nb` 要求非负输入，使用保持非负性的 MaxAbsScaler；不会把负 embedding 静默平移成非负数据。

无 `predict_proba` 的分类器（如 RidgeClassifier、LinearSVC）使用训练集内部最多 3 折 sigmoid 概率校准，每类至少需要 2 个训练样本。标准化也在每个校准训练折内拟合，验证集和测试集不参与。校准行为记录在结果 warnings 中。参见 [scikit-learn 概率校准 API](https://scikit-learn.org/stable/modules/generated/sklearn.calibration.CalibratedClassifierCV.html)。

机器学习路径不使用深度学习的 `epochs`、AdamW 参数或梯度累积；迭代次数、树数量等通过 `model_params` 设置。验证集只用于评估，不会自动调参或选择最优模型。

回归标签沿用现有接口的训练集均值/标准差归一化，报告和预测恢复到原始单位。分类 `class_weight` / `class_weights` 字典使用从 0 开始的编码类别索引；类别顺序保存在 `model.json` 的 `target_schema.class_names` 中，也可使用算法支持的 `"balanced"` 等选项。

输出包括 `model.json`、`estimator.joblib`、`metrics.json`、`data_inspection.json`、`splits.csv`、`test_predictions.jsonl`。估计器文件包含拟合好的特征变换；manifest 记录特征顺序、配置、标签语义和 Python/依赖版本。机器学习没有神经网络训练 history，结果的 `history` 为空。

joblib 加载可能执行 Python 代码，因此机器学习推理必须显式设置 `trust_model=True`；仅用于自己训练或来源可信的工件。重新加载会核对 manifest 记录的 Python 和依赖版本，版本不一致时要求使用原训练环境。该约束依据 [scikit-learn 模型持久化说明](https://scikit-learn.org/stable/model_persistence.html)。

可直接编辑并运行 `others/examples/train_machine_learning_properties.py`；`--demo` 使用人工数据验证完整训练和重新加载流程：

```bash
python others/examples/train_machine_learning_properties.py --demo
```

### 用户自定义深度学习分类与回归

先安装训练与常见表格格式依赖：`python -m pip install "dnakit[training]"`。

`train_custom_property_model()` 接受一个本地数据文件，自动或按用户指定的映射找出 ID、DNA/RNA 序列、一个或多个标签以及可选 split。默认随机种子为 42，训练/验证/测试比例为 0.8/0.1/0.1；若映射了已有 split 列，则使用文件中的划分。`task_type` 必须由用户明确指定为 `classification` 或 `regression`，不会仅凭标签外观替用户决定科学任务。

内置读取器支持 CSV、TSV、其他单字符分隔文本、JSON、JSONL、YAML、XML、Excel、Parquet、SQLite、FASTA、FASTQ 和 GenBank。嵌套 JSON/YAML/XML 可使用点路径；FASTA/FASTQ 描述中的 `key=value` 会映射为 `header.key`。其他格式可实现 `DataReader` 协议。大模型只生成受约束的声明式 `DataMapping`，不会生成或执行 Python、SQL、正则表达式或 shell 代码。

#### 标签语义

- 一个分类标签列：二分类或多分类；
- 多个二值标签列：自动识别 one-hot 多分类或多标签分类，也可用 `classification_mode` 明确指定；
- 多个非二值分类列：多输出分类，每列保留独立类别表；
- 一个或多个数值标签列：单目标或多目标回归。

分类报告包括 Accuracy、Balanced Accuracy、Precision、Recall、F1、Jaccard、MCC、Cohen's kappa、Log loss、AUROC、Average Precision、混淆矩阵及逐类别/逐标签结果；多标签任务另含 subset accuracy、Hamming loss、LRAP 和 coverage error。回归报告包括 MAE、MSE、RMSE、R²、explained variance、median absolute error、max error、bias、MAPE、SMAPE、Pearson 和 Spearman，并同时给出汇总与逐目标结果。F1、AUROC 等分类指标不会伪装成回归指标；样本不足或数学上未定义的指标返回 `null`。

#### 模型与微调方式

默认 `finetune_method="full"`。当前 11 个模型均提供全量微调和冻结骨干的 `linear_probe`：DNABERT-2、NT-v2、HyenaDNA、Caduceus、GROVER、LucaOne、GENERator、Enformer、Evo 2、JanusDNA、AlphaGenome。全量微调更新骨干与新任务头，线性探测只更新线性任务头。请用 `available_custom_training_models()` 获取支持矩阵；能力登记不表示本机已安装模型依赖或已完成全部官方大权重训练。

需要 checkpoint 自带代码的模型仍须设置 `allow_remote_code=True`。Enformer 需要 `enformer-pytorch`；Evo 2 需要官方 Evo2/Vortex CUDA 环境；JanusDNA 需要官方源码环境；AlphaGenome 需要本地研究版源码、JAX/Haiku/Orbax 和 PyTorch，通过 JAX 梯度更新骨干，不调用远程推理 API。可通过 `model_source_path`（CLI：`--model-source-path`）指定源码目录；迁移模型后，预测接口也接受该参数覆盖原目录。加载原始 `.bin` 权重须使用 PyTorch 2.6 或更新版本以满足 Transformers 的安全加载要求。

旧版 AlphaGenome 源码仅提供推理归一化统计，因此全量微调保留这些运行统计，但仍更新可学习的 scale/offset 和其他骨干权重；支持训练模式的新版源码会同步更新运行统计，实际策略会写入骨干清单。

Enformer 全量微调默认以 4096 bp 分块，关闭原轨道输出的中心裁剪并使用可变长度位置编码，按有效位置池化后平均各块；`max_length` 可调整分块长度。Enformer 线性探测显式设置较短 `max_length` 时采用同样的短上下文处理。Evo 2 全量微调使用最后一层归一化表征，让全部骨干层参与训练；线性探测仍使用登记的中间层。两种方法的表征定义应随实验记录。Evo 2 官方构造器会使用所有可见 GPU，应先设置 `CUDA_VISIBLE_DEVICES` 仅暴露训练用的一张卡，再使用 `device="cuda:0"`；当前训练器仍为单卡。

所有全量模型均保存骨干，按最佳验证轮次恢复后评估；重新加载不需要原始预训练权重，但仍需要模型运行环境。冻结探测保存任务头及原骨干位置，推理时需要相同预训练权重。训练和推理均显示进度，可用 `show_progress=False` 关闭。当前验证包含官方源码缩小架构的梯度、冻结、保存重载和端到端测试，不代表 11 个官方大权重均已完成训练实验。

AlphaGenome 另外通过了本地官方预训练权重的单卡 GPU 验证：2048 bp 输入完成一次前向、反向和 AdamW 参数更新；该检查不代表完整训练、多轮收敛或完整上下文长度验证。

#### 训练参数（Python API）

直接在 `CustomPropertyTrainingConfig` 中修改参数，无需命令行参数。全量微调与线性探测使用相同默认值：

全量微调使用 `dtype="auto"` 或 `"float32"`；支持的 CUDA 后端也可使用 `"bfloat16"`。纯 `"float16"` 优化器参数会被拒绝，以避免 AdamW 产生非有限权重。AlphaGenome 使用 `"auto"` 或 `"float32"`。

| Python 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `workers` | `6` | CPU 数据加载进程数；`0` 表示在主进程加载 |
| `batch_size` | `32` | 每个训练设备的单批样本数，不是梯度累积后的批大小 |
| `learning_rate` | `1e-4` | AdamW 学习率，即 lr |
| `beta1` | `0.9` | AdamW 一阶矩衰减系数 |
| `beta2` | `0.95` | AdamW 二阶矩衰减系数 |
| `weight_decay` | `0.1` | AdamW 权重衰减 |

`workers` 用于训练、验证和测试的数据取样及 CPU 批组装，编码器和 GPU 运算仍在主进程；它不控制初次文件解析或大模型请求的并发数。固定种子决定划分和每轮样本顺序，调整 workers 不会改变顺序。当前训练器在 `device` 指定的单个设备运行（例如 `"cuda:0"`），尚未实现多卡 DDP；不会根据显卡数量自动增大批大小。`gradient_accumulation_steps` 默认 `1`，可通过 Python 配置。

使用 `workers>0` 的 Python 脚本需要下例的 `if __name__ == "__main__":` 保护；Notebook 或不支持子进程的环境建议设置 `workers=0`。数据加载进程在训练完成、提前停止或异常时释放。

```python
import os

from dnakit.predictions import (
    CustomPropertyTrainingConfig,
    LLMSchemaConfig,
    train_custom_property_model,
)

def main():
    result = train_custom_property_model(
        "my_labels.xlsx",
        config=CustomPropertyTrainingConfig(
            task_type="regression",  # 或 "classification"
            model="dnabert2",
            finetune_method="full",  # 或 "linear_probe"
            output_dir="runs/my-property",
            workers=6,
            batch_size=32,
            learning_rate=1e-4,
            beta1=0.9,
            beta2=0.95,
            weight_decay=0.1,
            seed=42,
            train_ratio=0.8,
            validation_ratio=0.1,
            test_ratio=0.1,
            allow_remote_code=True,
            llm=LLMSchemaConfig(
                provider="deepseek",
                model="deepseek-v4-flash",
                api_key=os.environ["DEEPSEEK_API_KEY"],
            ),  # 不需要大模型辅助时省略 llm
        ),
    )
    print(result.metrics["test"])


if __name__ == "__main__":
    main()
```

训练结果目录包含 `model.json`、`head.npz`、`metrics.json`、`history.json`、`data_inspection.json`、`splits.csv` 和 `test_predictions.jsonl`；全量微调还包含 `backbone/`。`model.json` 的 `training` 字段记录实际 workers、批大小、学习率、两个 beta 和 weight decay。`predict_custom_properties()` 可重新加载这些产物进行推理。

仓库提供可直接编辑参数的脚本 `others/examples/train_custom_properties.py`：修改输入文件与配置后，在 IDE 中运行即可；脚本默认不启用大模型辅助。

#### 大模型服务商与数据分享范围

不仅支持 DeepSeek。`LLMSchemaConfig(provider=..., model=...)` 内置以下连接预设，统一调用 OpenAI 兼容的 `POST /chat/completions`，不需要安装各厂商 SDK。`available_schema_llm_providers()` 可离线查询这些预设。除兼容旧用法的 DeepSeek 默认模型外，必须自行填写账号可用的模型 ID；配置不会自动切换模型或服务商。

| `provider` | 默认 `base_url` | 默认密钥环境变量 |
| --- | --- | --- |
| `deepseek` | `https://api.deepseek.com` | `DEEPSEEK_API_KEY` |
| `gemini` | `https://generativelanguage.googleapis.com/v1beta/openai` | `GEMINI_API_KEY` |
| `qwen` | `https://dashscope.aliyuncs.com/compatible-mode/v1` | `DASHSCOPE_API_KEY` |
| `kimi` | `https://api.moonshot.cn/v1` | `MOONSHOT_API_KEY` |
| `openrouter` | `https://openrouter.ai/api/v1` | `OPENROUTER_API_KEY` |
| `openai` | `https://api.openai.com/v1` | `OPENAI_API_KEY` |
| `custom` | 必须显式指定 | `LLM_API_KEY` |

例如，仅分享字段与类型结构，使用 Google Gemini 辅助分析：

```python
from dnakit.predictions import LLMSchemaConfig, inspect_property_data

llm = LLMSchemaConfig(
    provider="gemini",
    model="gemini-3.8-flash",  # 可替换为账号可用的文本/JSON 模型。
    include_sample_values=False,
)
inspection = inspect_property_data("my_labels.xlsx", llm=llm)
print(inspection.mapping.to_dict())
print(inspection.recommendations)
## 训练时将同一个 llm 传给 CustomPropertyTrainingConfig(llm=llm, ...)。
```

`api_key` 仍可在 Python 中传入，或用 `api_key_env` 指定自己的环境变量名；不要将真实密钥写进示例或版本库。`base_url` 可覆盖预设，密钥只发送到你配置的地址，HTTP 重定向被拒绝。Gemini 此处使用其[官方 OpenAI 兼容接口](https://ai.google.dev/gemini-api/docs/openai)，不是原生 `generateContent` 或 Vertex AI 接口。

千问默认地址为北京地域，官方仍支持该域名；也可以填控制台提供的业务空间专属地址，例如 `https://<业务空间ID>.cn-beijing.maas.aliyuncs.com/compatible-mode/v1`。地域地址和 API Key 必须匹配，见[千问接口说明](https://help.aliyun.com/zh/model-studio/compatibility-of-openai-with-dashscope)。Kimi 国际站可使用 `base_url="https://api.moonshot.ai/v1"`，参见 [Kimi 快速开始](https://platform.kimi.ai/docs/overview)。OpenRouter 的模型名须使用其平台的完整 ID，且应检查目标模型支持的参数，见 [OpenRouter 快速开始](https://openrouter.ai/docs/quickstart)和[结构化输出说明](https://openrouter.ai/docs/guides/features/structured-outputs)；数据将由其平台路由到所选模型的提供方。

兼容参数可以直接在 Python 中调整：

- `temperature=None`（默认）：不发送温度参数，适用于不接受自定义温度的模型；需要时可显式设置。
- `token_limit_parameter="auto"`（默认）：OpenAI 使用 `max_completion_tokens`，其他兼容服务使用 `max_tokens`；也可显式选择这两个字段之一。依据 [OpenAI Docs 的 Chat Completions 参数说明](https://developers.openai.com/api/reference/python/resources/chat/subresources/completions/methods/create)，旧 `max_tokens` 不兼容部分推理模型。
- `json_mode=True`（默认）：要求 JSON Object 输出。若模型不支持该接口参数，可设为 `False`；提示词仍要求 JSON，返回结果仍须经过严格映射校验。
- `extra_body`：传递模型专用参数，例如支持关闭思考的千问模型可设 `{"enable_thinking": False}`；支持该选项的 Kimi 模型可设 `{"thinking": {"type": "disabled"}}`。不会统一强行关闭所有模型的思考模式，需按[千问](https://help.aliyun.com/zh/model-studio/qwen-structured-output)或 [Kimi 模型参数](https://platform.kimi.ai/docs/api/models-overview)确认支持情况。请求当前为非流式；只支持流式思考的模型需关闭思考或改用支持非流式的模型。

`extra_body` 不允许覆盖消息、鉴权信息、模型、输出 token 上限等核心字段。输出被截断时会明确报错，应调整 `max_output_tokens` 或模型专用思考参数，不会自动扩大额度重试。服务商预设表示接口接入能力，不保证该平台的每个模型都支持此任务；当前已做模拟接口测试，未用各平台付费密钥逐一联网验证。

默认最多分享 8 条样本、12,000 个字符及结构摘要，**不会整表上传**；小文件也至少保留一条完整记录不分享。`include_sample_values=False` 只发送字段/类型结构。`DataInspection` 记录分享的样本数和字符数；密钥不进入请求正文、结果、异常展示或训练产物。

## Enformer 各任务权重路径 {#task-checkpoints}

- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`H2AFZ.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H2AFZ.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`H3K27ac.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K27ac.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`H3K27me3.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K27me3.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`H3K36me3.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K36me3.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`H3K4me1.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K4me1.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`H3K4me2.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K4me2.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`H3K4me3.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K4me3.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`H3K9ac.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K9ac.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`H3K9me3.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K9me3.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`H4K20me1.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H4K20me1.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`enhancers.ckpt`；默认放在 `./ckpt/enformer-benchmarks/enhancers.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`enhancers_types.ckpt`；默认放在 `./ckpt/enformer-benchmarks/enhancers_types.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`promoter_all.ckpt`；默认放在 `./ckpt/enformer-benchmarks/promoter_all.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`promoter_no_tata.ckpt`；默认放在 `./ckpt/enformer-benchmarks/promoter_no_tata.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`promoter_tata.ckpt`；默认放在 `./ckpt/enformer-benchmarks/promoter_tata.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`splice_sites_acceptors.ckpt`；默认放在 `./ckpt/enformer-benchmarks/splice_sites_acceptors.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`splice_sites_all.ckpt`；默认放在 `./ckpt/enformer-benchmarks/splice_sites_all.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`splice_sites_donors.ckpt`；默认放在 `./ckpt/enformer-benchmarks/splice_sites_donors.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`demo_coding_vs_intergenomic_seqs.ckpt`；默认放在 `./ckpt/enformer-benchmarks/demo_coding_vs_intergenomic_seqs.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`demo_human_or_worm.ckpt`；默认放在 `./ckpt/enformer-benchmarks/demo_human_or_worm.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`drosophila_enhancers_stark.ckpt`；默认放在 `./ckpt/enformer-benchmarks/drosophila_enhancers_stark.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`dummy_mouse_enhancers_ensembl.ckpt`；默认放在 `./ckpt/enformer-benchmarks/dummy_mouse_enhancers_ensembl.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`human_enhancers_cohn.ckpt`；默认放在 `./ckpt/enformer-benchmarks/human_enhancers_cohn.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`human_enhancers_ensembl.ckpt`；默认放在 `./ckpt/enformer-benchmarks/human_enhancers_ensembl.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`human_ensembl_regulatory.ckpt`；默认放在 `./ckpt/enformer-benchmarks/human_ensembl_regulatory.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`human_nontata_promoters.ckpt`；默认放在 `./ckpt/enformer-benchmarks/human_nontata_promoters.ckpt`，也可使用自定义目录或精确路径。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#property-prediction-checkpoint-cache)）`human_ocr_ensembl.ckpt`；默认放在 `./ckpt/enformer-benchmarks/human_ocr_ensembl.ckpt`，也可使用自定义目录或精确路径。
## 工程化：脚本执行、报告渲染与验证范围 {#engineering}

### 执行随包坐标转换脚本 {#engineering-bundled-script-example}

此例在当前 Python 环境运行随包脚本，把 BED 0-based 半开坐标转换成 GFF 1-based 闭区间；它不是基因组版本转换。

```python
import dnakit

print(dnakit.workflows.scripts("genomic-coordinates", verbose=False))
path = dnakit.workflows.resource_path(
    "genomic-coordinates", "scripts/convert_coords.py", verbose=False,
)
print(path)
result = dnakit.workflows.run(
    "genomic-coordinates", "scripts/convert_coords.py",
    arguments=["--from", "bed", "--to", "gff", "chr1", "999", "1000", "--format", "json"],
    timeout_seconds=60, verbose=False,
)
if result["returncode"] != 0:
    raise RuntimeError(result["output"])
print(result["output"])
```

`run_support()` 运行已打包的支持流程，允许的名称为 `clinical-variant-reporter`、`clinical-trial-finder`、`omics-target-evidence-mapper`、`gwas-lookup`、`fine-mapping`。先通过 `arguments=["--help"]` 查看各自参数；脚本导入的科学依赖仍需准备。参数、资源校验和返回值详见[工作流接口](../api/features/14_engineering.md#execution)。

### 从已有报告生成 PDF {#engineering-report-example}

安装报告渲染依赖：

```bash
python -m pip install -e ".[reports]"
```

下面的示例使用临时目录，展示真实的输入、输出和来源记录：

```python
from pathlib import Path
from tempfile import TemporaryDirectory
import dnakit

with TemporaryDirectory() as directory:
    root = Path(directory)
    report = root / "report.md"
    report.write_text("# Results\nUser-supplied report content.\n", encoding="utf-8")
    result = dnakit.workflows.write_wes_report(report, root / "pdf", verbose=False)
    print(result["output"])
    print(result["manifest"])
    print(Path(result["output"]).exists())  # True
```

实际保存时使用自己的输出目录。该函数渲染已有内容，不增加医学解释；`language` 只切换英文或西班牙文标签，不翻译正文。标题、Logo、字体和覆盖选项见[完整参数](../api/features/14_engineering.md#dnakit.workflows.write_wes_report)。

### MCP 与验证范围 {#engineering-mcp-validation}

同一函数也可由 MCP 的工具搜索、参数描述和调用入口使用，详见[MCP 教程](../llm_tools/mcp.md)。Python、CLI、MCP 是同一能力的不同入口，统计时不重复计数。

文档校验检查全部 484 个函数是否有唯一条目、签名及参数是否与源码一致，并构建中英文网页。文档覆盖与离线示例验证不等于所有远程服务、模型和完整科学流程都已端到端实测。

## ENCODE 查询与下载预览 {#encode-download}

```python
import json
import dnakit

experiments = json.loads(dnakit.encode.search_experiments(
    assay_title="ATAC-seq", organism="Homo sapiens", limit=3,
))
print(experiments)
preview = json.loads(dnakit.encode.batch_download(
    download_dir="downloads/encode", assay_title="ATAC-seq",
    assembly="GRCh38", file_format="bed", limit=3, dry_run=True,
))
print(preview)
```

`dry_run=True` 查询文件列表；设为 `False` 后下载文件。ENCODE 接口返回文本：多数为 JSON，引用与表格导出也可返回 BibTeX、RIS、CSV 和 TSV；这些格式不使用 `json.loads()`。追踪、凭据管理与下载可能写入本地状态。相关入口见[实验查询](../api/features/22_database_query_encode.md#encode)、[集合管理](../api/features/10_datasets.md#encode-collection)和[凭据管理](../api/features/14_engineering.md#encode-credentials)。
