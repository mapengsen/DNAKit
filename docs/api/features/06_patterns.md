# 序列模式搜索

在 DNA 序列中搜索 motif/PWM、功能位点、回文、倒置重复、串联重复和微卫星等模式。

命中结果使用 0-based 半开区间坐标；这些 API 只报告序列模式命中，不替代实验活性、结合强度或编辑效率预测。

## 1) PAT-001 / SIM-003 通用motif搜索

- **作用：** 按 exact、IUPAC 或受限正则规则扫描指定 motif，返回每次命中的坐标、匹配文本和链方向，用于定位已知序列模式。
- **API：** `dnakit.patterns.scan_motif(value[必须], motif[必须], mode[可选], name[可选], strand[可选], overlapping[可选], merge_strands[可选], max_matches[可选], max_scan_length[可选], max_pattern_length[可选], max_scan_cells[可选])`、`dnakit.patterns.scan_pwm(value[必须], pwm[必须], threshold[必须], background[可选], pseudocount[可选], strand[可选], max_matches[可选], max_scan_length[可选], max_pwm_length[可选], max_score_cells[可选])`、`dnakit.patterns.PWM(name[必须], matrix[必须])`
- **输入：** 必填序列和 motif/正则/PWM；可选模式、strand、重叠、背景、阈值和命中上限。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_motif

result = scan_motif(DNASequence("AAAA"), "AA", strand="forward")
print([(hit.symbol_location.start, hit.symbol_location.end) for hit in result.hits])
```

- **示例结果：**

```text
[(0, 2), (1, 3), (2, 4)]
```

- **限制：** 正则仅支持安全的 DNA 子集并标记为内部重实现；PWM 由调用方提供。固定长度 motif 可跨环状原点，正则不可跨原点。

## 2) PAT-003 起止密码子

- **作用：** 在正反两条链的六个阅读框中查找起始和终止密码子，返回类型、阅读框和位置，为 ORF 与编码区域分析提供候选位点。
- **API：** `dnakit.patterns.scan_codon_sites(value[必须], genetic_code[可选], start_codons[可选], stop_codons[可选], strand[可选], max_matches[可选], max_codon_checks[可选])`
- **输入：** 必填序列；可选遗传密码表、自定义 start/stop 集合和 strand。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_codon_sites

result = scan_codon_sites(DNASequence("ATGAAATAA"), strand="forward")
print([(hit.kind, hit.codon, hit.frame) for hit in result.hits])
```

- **示例结果：**

```text
[('start', 'ATG', 1), ('stop', 'TAA', 1), ('stop', 'TGA', 2)]
```

- **限制：** 返回候选密码子位点而非完整 ORF；Gap 分段扫描，工作量和命中数受硬上限约束。

## 3) PAT-004 启动子motif

- **作用：** 扫描内置共识序列或用户提供的启动子 motif，返回候选位置和链方向，用于规则型启动子区域筛查，但不预测启动子活性。
- **API：** `dnakit.patterns.scan_promoter_motifs(value[必须], motifs[可选], strand[可选], max_matches[可选], max_scan_length[可选], max_scan_cells[可选], max_motifs[可选])`
- **输入：** 必填序列；可选 motif 映射、strand 和资源上限。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_promoter_motifs

result = scan_promoter_motifs(DNASequence("GGTATAATCC"), strand="forward")
print([hit.motif_name for hit in result.hits])
```

- **示例结果：**

```text
['bacterial_minus_10_consensus']
```

- **限制：** 内置目录很小且版本固定；只报告模式命中，不预测真实启动子或转录活性。

## 4) PAT-005 TF motif

- **作用：** 使用用户提供的位置权重矩阵逐窗打分，返回超过阈值的候选转录因子结合位点及分数，用于 motif 候选筛查而非实际结合强度预测。
- **API：** `dnakit.patterns.PWM(name[必须], matrix[必须])`、`dnakit.patterns.scan_tf_pwm(value[必须], tf_name[必须], pwm[必须], threshold[必须], background[可选], pseudocount[可选], strand[可选], max_matches[可选], max_scan_length[可选], max_pwm_length[可选], max_score_cells[可选])`
- **输入：** 必填序列、TF 名称、PWM 和阈值；可选背景、pseudocount 与 strand。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.patterns import PWM, scan_tf_pwm

pwm = PWM("input", {"A": [4], "C": [0], "G": [0], "T": [0]})
result = scan_tf_pwm(DNASequence("AA"), "TF-X", pwm, threshold=1.0, strand="forward")
print([hit.motif_name for hit in result.hits])
```

- **示例结果：**

```text
['TF-X', 'TF-X']
```

- **限制：** 不内置 JASPAR 库，也不预测真实结合强度；PWM、来源版本和阈值由调用方负责。

## 5) PAT-006 限制酶位点

- **作用：** 按内置或自定义限制酶定义查找识别序列，返回正反链命中及切割坐标，用于规划酶切实验和判断片段边界。
- **API：** `dnakit.patterns.scan_restriction_sites(value[必须], enzymes[必须], max_matches[可选], max_scan_length[可选], max_scan_cells[可选], max_enzymes[可选])`、`dnakit.patterns.RestrictionEnzyme(name[必须], recognition_sequence[必须], top_cut_offset[必须], bottom_cut_offset[必须], source[可选])`
- **输入：** 必填序列和酶名称/定义列表；可选扫描与命中上限。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_restriction_sites

result = scan_restriction_sites(DNASequence("AGAATTCCCGGG"), ["EcoRI"])
hit = result.hits[0]
print(hit.enzyme, hit.top_cut, hit.bottom_cut)
```

- **示例结果：**

```text
EcoRI 2 6
```

- **限制：** 不分发完整 REBASE，也不建模甲基化状态；未知长度 Gap 后的绝对切点可能为 `None`。

## 6) PAT-007 CRISPR PAM

- **作用：** 按指定 PAM 规则在正反链查找 CRISPR guide 候选，返回 spacer、PAM、坐标和方向，用于候选枚举，不直接预测编辑效率或脱靶风险。
- **API：** `dnakit.patterns.scan_pam_candidates(value[必须], rule[必须], guide_length[可选], strand[可选], min_gc[可选], max_gc[可选], exclude_motifs[可选], allow_ambiguous_guides[可选], max_matches[可选], max_scan_length[可选], max_pam_length[可选], max_scan_cells[可选], max_exclude_motifs[可选], max_filter_cells[可选])`、`dnakit.patterns.PAMRule(name[必须], pam[必须], pam_side[必须], guide_length[必须], source[可选])`
- **输入：** 必填序列和核酸酶名称/PAMRule；可选 guide 长度、strand、GC 范围和排除 motif。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.patterns import scan_pam_candidates

result = scan_pam_candidates(DNASequence("A" * 20 + "TGG"), "SpCas9", strand="forward")
print(result.hits[0].guide_sequence, result.hits[0].pam_sequence)
```

- **示例结果：**

```text
AAAAAAAAAAAAAAAAAAAA TGG
```

- **限制：** 只生成序列候选，不预测编辑效率或 off-target 风险；候选不能跨 Gap。

## 7) PAT-009 回文序列

- **作用：** 查找一段序列与其反向互补完全一致的回文区域，返回位置和长度，用于识别可能的限制酶位点或对称结构模式。
- **API：** `dnakit.patterns.find_reverse_complement_palindromes(value[必须], min_length[可选], max_length[可选], maximal_per_start[可选], max_comparisons[可选], max_comparison_cells[可选], max_matches[可选])`
- **输入：** 必填序列；可选最小/最大长度、每起点仅保留最大命中及资源上限。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.patterns import find_reverse_complement_palindromes

result = find_reverse_complement_palindromes(
    DNASequence("GAATTC"), min_length=4, max_length=6
)
print(any(hit.sequence == "GAATTC" for hit in result.hits))
```

- **示例结果：**

```text
True
```

- **限制：** 使用符号相等或 IUPAC 集合相交的兼容规则做反向互补匹配；不允许超出 IUPAC 兼容范围的 mismatch 或中心 Gap，也不支持环状序列。

## 8) PAT-010 倒置重复

- **作用：** 查找由两段反向互补序列和中间 loop 构成的倒置重复，返回两臂及 loop 的坐标，用于筛查潜在发卡相关序列模式。
- **API：** `dnakit.patterns.find_inverted_repeats(value[必须], min_arm_length[可选], max_arm_length[可选], min_loop_length[可选], max_loop_length[可选], max_comparisons[可选], max_comparison_cells[可选], max_matches[可选])`
- **输入：** 必填序列；可选 arm/loop 长度范围及比较、命中上限。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.patterns import find_inverted_repeats

result = find_inverted_repeats(
    DNASequence("ACGTAAACGT"),
    min_arm_length=4,
    max_arm_length=4,
    min_loop_length=2,
    max_loop_length=2,
)
print(result.hits[0].left_arm, result.hits[0].loop_length)
```

- **示例结果：**

```text
ACGT 2
```

- **限制：** arm 使用符号相等或 IUPAC 集合相交的兼容规则；不支持超出该范围的 mismatch/indel。它是潜在发卡的序列 proxy，不是折叠或能量预测。

## 9) PAT-011 串联重复

- **作用：** 查找连续排列的短重复单元，返回重复单元、次数、区间和覆盖长度，用于识别串联重复及量化其规模。
- **API：** `dnakit.patterns.find_tandem_repeats(value[必须], min_unit_length[可选], max_unit_length[可选], min_repeats[可选], min_repeats_by_unit[可选], overlapping[可选], max_comparisons[可选], max_comparison_cells[可选], max_matches[可选])`
- **输入：** 必填序列；可选单元长度、最少重复次数、重叠策略和资源上限。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.patterns import find_tandem_repeats

result = find_tandem_repeats(
    DNASequence("ATATAT"), min_unit_length=1, max_unit_length=3, min_repeats=2
)
print(result.hits[0].unit, result.hits[0].repeat_count)
```

- **示例结果：**

```text
AT 3
```

- **限制：** 不允许 mismatch/indel/中断，不等同于 Tandem Repeats Finder；不跨 Gap，不支持环状序列。

## 10) PAT-012 微卫星

- **作用：** 专门查找重复单元长度为 1–6 bp 的微卫星，返回 motif、重复次数和坐标，用于标记 STR 候选区域。
- **API：** `dnakit.patterns.find_microsatellites(value[必须], min_repeats_by_unit[可选], overlapping[可选], max_comparisons[可选], max_comparison_cells[可选], max_matches[可选])`
- **输入：** 必填序列；可选每种单元长度的最少重复次数、重叠策略和资源上限。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.patterns import find_microsatellites

result = find_microsatellites(DNASequence("AAAAAACACACA"))
print([(hit.unit, hit.repeat_count) for hit in result.hits])
```

- **示例结果：**

```text
[('A', 6), ('CA', 3)]
```

- **限制：** 当前仅检测不中断的 exact STR，不允许 mismatch/indel；自定义阈值必须覆盖所有 1–6 bp 单元长度。
