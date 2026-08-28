# 序列功能搜索

在 DNA 序列中搜索 motif/PWM、功能位点、回文、倒置重复、串联重复和微卫星等模式。

## 1) PAT-001 / SIM-003 motif搜索

- **含义：** motif 是 DNA 中具有某种特征的短序列模式，可以理解为 DNA 序列里的“关键词”。
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

## 2) PAT-003 起止密码子

- **含义：** 起始密码子表示蛋白质翻译可能从这里开始，终止密码子表示蛋白质翻译到这里结束。
- **计算规则：** 确定性规则（[FAQ 详细解释](../../faq.md#pattern-matching-strategy)）。
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

## 3) PAT-004 启动子motif

- **含义：** 启动子是帮助基因开始转录的调控区域；启动子 motif 是其中常见的短序列特征，可以理解为基因转录的“启动标记”。
- **计算规则：** 确定性规则（[FAQ 详细解释](../../faq.md#pattern-matching-strategy)）。
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

## 4) PAT-005 TF motif

- **含义：** TF motif 是转录因子偏好识别和结合的 DNA 序列模式，可能参与调节基因是否表达以及表达强弱。
- **计算规则：** 确定性规则（[FAQ 详细解释](../../faq.md#pattern-matching-strategy)）。
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

## 5) PAT-006 限制酶位点

- **含义：** 限制酶位点是限制酶能够识别并切割的特定 DNA 序列，可以理解为“分子剪刀”的下刀位置。
- **计算规则：** 确定性规则（[FAQ 详细解释](../../faq.md#pattern-matching-strategy)）。
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

## 6) PAT-007 CRISPR PAM

- **含义：** PAM 是紧邻 CRISPR 目标序列的短标记，Cas 蛋白通常需要先识别 PAM，才能结合并切割附近的 DNA。
- **计算规则：** 确定性规则（[FAQ 详细解释](../../faq.md#pattern-matching-strategy)）。
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

## 7) PAT-009 回文序列

- **含义：** DNA 回文序列是与自身反向互补序列相同的片段，例如 `GAATTC`。
- **计算规则：** 确定性规则（[FAQ 详细解释](../../faq.md#pattern-matching-strategy)）。
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

## 8) PAT-010 倒置重复

- **含义：** 倒置重复是两段彼此反向互补、且中间可以隔着一段序列的 DNA 片段，可能折叠形成发卡结构。
- **计算规则：** 确定性规则（[FAQ 详细解释](../../faq.md#pattern-matching-strategy)）。
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

## 9) PAT-011 串联重复

- **含义：** 串联重复是同一个短序列单元首尾相接、连续出现多次，例如 `ATATAT` 是 `AT` 连续重复 3 次。
- **计算规则：** 确定性规则（[FAQ 详细解释](../../faq.md#pattern-matching-strategy)）。
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

## 10) PAT-012 微卫星

- **含义：** 微卫星又称 STR，是重复单元长度为 1–6 bp 的短串联重复，例如 `CACACACA`。
- **计算规则：** 确定性规则（[FAQ 详细解释](../../faq.md#pattern-matching-strategy)）。
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
