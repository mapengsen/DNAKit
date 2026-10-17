# 单倍型、家系与亲缘

本类收录 **17** 项有文档的功能；其中 16 项直接由输入的基因型或单倍型用纯 Python 在本地计算出结果（三联体与系谱逻辑、单倍型统计、连锁不平衡区块与标签 SNP、EM 单倍型重建、模拟、SVG 绘图）；`count` 是通用的记录计数函数，只统计非空行或列表元素的个数，不把它们解释成单倍型。只有 `impute` 始终返回 `status: "not_implemented"`，因为它需要输入中并不包含的已定相参考面板和填充引擎。不需要外部程序或可选库。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) analyze_trio_inheritance

- **作用：** 检查孩子的基因型与双亲的孟德尔一致性，并报告各位点传递的等位基因。

- **处理方法：** 每个位点上，孩子的两个等位基因必须能由父亲的一个等位基因加母亲的一个等位基因组成，否则该位点为 `inconsistent`（孟德尔错误）。任一调用缺失的位点记为 `missing` 并跳过。若恰有一种（父源，母源）等位基因分配能解释孩子的基因型，就报告该分配；若有多种（例如三人都是杂合）则不确定来源。孟德尔错误率 = 不一致位点数 / 有信息位点数（一致 + 不一致）。假定为常染色体二倍体基因型；没有基因分型错误模型，也不处理 X 染色体。

- **API：** `dnakit.haplotypes.analyze_trio_inheritance(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`{"child": [...], "father": [...], "mother": [...]}`，每个位点一个基因型调用（三个列表等长），可选 `"sites"` 标签列表。基因型调用可写作 `"a/b"` 或 `"a|b"`（等位基因为整数）、数对 `(a, b)`、alt 等位基因剂量 0/1/2，缺失用 `None` 或 `"./."`。 `options`：无。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`site_count`、`informative_sites`、`consistent`、`inconsistent`、`missing`、`mendelian_error_rate`（没有有信息位点时为 null）、`inconsistent_indices`、`sites`（每个位点一项，含 `index`、可选的 `site` 标签、`status`，可确定时含 `paternal_allele` / `maternal_allele`）以及说明字符串 `rule`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import analyze_trio_inheritance

trio = {
    "sites": ["rs1", "rs2", "rs3", "rs4"],
    "child": ["0/1", "1/1", "0/0", "./."],
    "father": ["0/0", "0/1", "0/1", "0/0"],
    "mother": ["0/1", "0/1", "1/1", "0/0"],
}
result = analyze_trio_inheritance(trio)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "analyze_trio_inheritance",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 158,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "dec8cb6e39f7932d76a0417af1afaff5105efb3538d6d64a2e3718957d0cb659"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "site_count": 4,
    "informative_sites": 3,
    "consistent": 2,
    "inconsistent": 1,
    "missing": 1,
    "mendelian_error_rate": 0.3333333333333333,
    "inconsistent_indices": [2],
    "sites": [
      {
        "index": 0,
        "site": "rs1",
        "status": "consistent",
        "paternal_allele": 0,
        "maternal_allele": 1
      },
      {
        "index": 1,
        "site": "rs2",
        "status": "consistent",
        "paternal_allele": 1,
        "maternal_allele": 1
      },
      {
        "index": 2,
        "site": "rs3",
        "status": "inconsistent"
      },
      {
        "index": 3,
        "site": "rs4",
        "status": "missing"
      }
    ],
    "rule": "a site is consistent if the child's two alleles can be formed by one allele from each parent"
  }
}
```

## 2) calculate_pedigree_likelihood

- **作用：** 对小型系谱上观察到的双等位基因型，通过精确枚举计算其似然，不考虑基因分型错误。

- **处理方法：** 奠基者（founder，系谱中没有父母记录的个体）服从 Hardy-Weinberg 比例，alt 等位基因频率为 p（`(1-p)^2`、`2p(1-p)`、`p^2`）；子代基因型服从孟德尔传递，每位亲本以 剂量/2 的概率传递 alt 等位基因；位点之间相互独立。`genotypes` 中没有的成员或调用缺失的成员，通过枚举其可能基因型被精确地求和消去（每个位点最多 8 个未分型成员）。似然为各位点概率之积；任一位点概率为 0 时数据不可能出现，`log10_likelihood` 为 null。只支持双等位基因型（0/1）；不考虑位点间连锁、突变和基因分型错误。

- **API：** `dnakit.haplotypes.calculate_pedigree_likelihood(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`{"pedigree": [{"id", "father", "mother"}, ...], "genotypes": {id: [每个位点的调用]}}`。父母必须同时给出或同时省略（省略即奠基者）。调用可为 0/1 等位基因的 `"0/1"` 形式字符串、数对或剂量 0/1/2；所有列表等长。`options`（可选）：`allele_frequency`（(0, 1) 内的数，或每个位点一个数，默认 0.5）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`site_likelihoods`、`likelihood`、`log10_likelihood`（不可能时为 null）、`impossible`、`individuals`、`sites` 以及说明字符串 `model`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import calculate_pedigree_likelihood

pedigree = [
    {"id": "f"},
    {"id": "m"},
    {"id": "c", "father": "f", "mother": "m"},
]
genotypes = {
    "f": ["0/0", "0/1"],
    "m": ["1/1", "0/1"],
    "c": ["0/1", "0/1"],
}
result = calculate_pedigree_likelihood(
    {"pedigree": pedigree, "genotypes": genotypes},
    options={"allele_frequency": 0.2},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "calculate_pedigree_likelihood",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 161,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "62ffa5c42689f1bfee4521017ab68ee6092ee94bf57e8516c57147ec51c9f566"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "site_likelihoods": [0.02560000000000001, 0.05120000000000002],
    "likelihood": 0.0013107200000000011,
    "log10_likelihood": -2.882490073712319,
    "impossible": false,
    "individuals": 3,
    "sites": 2,
    "model": "founders in Hardy-Weinberg proportions with the given alt-allele frequency, Mendelian transmission, no genotyping error, sites independent; ungenotyped members summed out exactly"
  }
}
```

## 3) cluster

- **作用：** 按汉明距离对单倍型做单连接聚类（相互距离不超过 `max_distance` 个错配的单倍型构成的连通分量）。

- **处理方法：** 计算每对单倍型的汉明距离（相异位点数；`ignore_missing` 为 true 时，一对中出现 `?`、`N`、`.`、`-` 的位点被跳过），并用并查集连接所有距离 <= `max_distance` 的单倍型对；每个连通分量即一个簇。单连接意味着通过一串相近的邻居，簇内单倍型之间的距离可能大于 `max_distance`。簇按各自第一个成员的顺序列出。最多 5,000 条单倍型。

- **API：** `dnakit.haplotypes.cluster(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：等长的单倍型字符串列表，或 `{名称: 单倍型 | [单倍型]}`。`options`（可选）：`max_distance`（整数 >= 0，默认 1）、`ignore_missing`（布尔，默认 True）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`cluster_count`、`clusters`（每个簇含 `size`、以名称表示的 `members`、不重复的 `haplotypes` 和输入单倍型的 `indices`）、`max_distance` 以及 `linkage` 说明。

- **示例代码：**

```python
import json

from dnakit.haplotypes import cluster

result = cluster(["0000", "0001", "1111", "1110"], options={"max_distance": 1})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "cluster",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 32,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "35c2ec48596e2606e14353d91ff34b1fa688576c2075122e1582d191f9a9069c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "cluster_count": 2,
    "clusters": [
      {
        "size": 2,
        "members": ["hap_1", "hap_2"],
        "haplotypes": ["0000", "0001"],
        "indices": [0, 1]
      },
      {
        "size": 2,
        "members": ["hap_3", "hap_4"],
        "haplotypes": ["1110", "1111"],
        "indices": [2, 3]
      }
    ],
    "max_distance": 1,
    "linkage": "single (connected components of the graph joining haplotypes within max_distance mismatches)"
  }
}
```

## 4) compare

- **作用：** 按成对汉明距离比较单倍型，或按切换错误率（switch error rate）比较两组已定相的单倍型对。

- **处理方法：** `mode="hamming"`（默认）：给出成对相异位点数矩阵（除非 `ignore_missing` 为 false，一对中出现 `?`、`N`、`.`、`-` 的位点被跳过）；恰好两条单倍型时还列出相异位置、距离和参与比较的位点数。含多条单倍型的样本命名为 `name_1`、`name_2`……。`mode="switch_error"`：在真值和推断结果中都是杂合、且两者等位基因相同的位点上，记录推断的两条单倍型相对真值的方向；相邻的此类位点之间方向发生变化即为一次切换，`switch_error_rate = 切换次数 / (位点数 - 1)`（此类位点少于两个时为 null）。汉明模式最多 2,000 条单倍型。

- **API：** `dnakit.haplotypes.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`hamming` 模式为至少两条等长单倍型字符串的列表，或 `{名称: 单倍型 | [单倍型]}`；`switch_error` 模式为 `{"truth": [h1, h2], "inferred": [h1, h2]}`。`options`（可选）：`mode`（默认 `"hamming"`，或 `"switch_error"`）、`ignore_missing`（布尔，默认 True；仅汉明模式）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：汉明模式：`mode`、`names`、`distance_matrix`、`ignore_missing`，两条单倍型时另有 `differing_positions`、`distance`、`compared_sites`。切换错误模式：`mode`、`heterozygous_sites`、`switches`、`switch_error_rate` 以及说明字符串 `definition`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import compare

result = compare({"truth": ["0000", "1111"], "inferred": ["0011", "1100"]}, options={"mode": "switch_error"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 57,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3268637c7959c932534b3c1dee8c926c2330c85dcc4a5353cd07a2a734bda330"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "mode": "switch_error",
    "heterozygous_sites": 4,
    "switches": 1,
    "switch_error_rate": 0.3333333333333333,
    "definition": "among heterozygous sites called identically in truth and inferred, a switch is a change of relative phase between consecutive sites; rate = switches / (sites - 1)"
  }
}
```

## 5) construct

- **作用：** 用“全部成对连锁不平衡都达标”的规则，从已定相的 0/1 单倍型构建单倍型区块。

- **处理方法：** 单态位点和次要等位基因频率低于 `min_maf` 的位点被排除。两个位点间的连锁不平衡由所给染色体的单倍型频率计算：D = p(AB) - p(A)p(B)，r^2 = D^2 / (p(A)(1-p(A))p(B)(1-p(B)))，|D'| = |D| / Dmax。从左到右贪心地扩展区块，只要区块内每一对位点的所选 LD 统计量都 >= `threshold`，区块含至少 `min_snps` 个位点则保留，下一个区块从其后开始。这不是 Gabriel 等人的置信区间方法，样本很小时 D' 会被高估。

- **API：** `dnakit.haplotypes.construct(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：已定相单倍型列表（等长的字符串或 0/1 等位基因列表；至少两条单倍型，最多 2,000 个位点），或 `{名称: 单倍型}`。`options`（可选）：`statistic`（默认 `"dprime"`，或 `"r2"`）、`threshold`（默认 0.8）、`min_snps`（整数 >= 2，默认 2）、`min_maf`（默认 0.0）、`positions`（每个位点一个坐标，会作为各区块的 `start_position` / `end_position` 返回）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`blocks`（每个区块含 `first_site`、`last_site`、`sites`、`snp_count`，给出 `positions` 时还有 `start_position` / `end_position`）、`block_count`、`statistic`、`threshold`、`sites_total`、`sites_analysed`、`sites_excluded_monomorphic_or_low_maf` 以及说明字符串 `rule`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import construct

haplotypes = ["0000", "0011", "1100", "1111"]
result = construct(haplotypes, options={"positions": [100, 200, 300, 400]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "construct",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 32,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "095e1dbb5d3cd7bea5f10e204f15a9aeaff2bbb73b167c7c21679c30b16c01b7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "blocks": [
      {
        "first_site": 0,
        "last_site": 1,
        "sites": [0, 1],
        "snp_count": 2,
        "start_position": 100,
        "end_position": 200
      },
      {
        "first_site": 2,
        "last_site": 3,
        "sites": [2, 3],
        "snp_count": 2,
        "start_position": 300,
        "end_position": 400
      }
    ],
    "block_count": 2,
    "statistic": "dprime",
    "threshold": 0.8,
    "sites_total": 4,
    "sites_analysed": 4,
    "sites_excluded_monomorphic_or_low_maf": [],
    "rule": "greedy left-to-right: a block is a maximal run of >= min_snps polymorphic sites in which every pair has LD (|D'| or r^2) >= threshold; pairwise LD from the haplotype frequencies of the given chromosomes"
  }
}
```

## 6) count

- **作用：** 统计输入的记录数（文本的非空行数或列表的元素个数）；不会把记录解释成单倍型。

- **处理方法：** 在共享本地运行时（`src/dnakit/_candidate_runtime.py` 中的 `run_candidate_operation` 与辅助函数 `_local_result`）中统计记录数：文本（或字节，或 `Path` 的文件内容）按非空行拆分，每行一条记录；列表按元素个数计；含有键 `records`、`items`、`values`、`data`、`sequences`、`reads` 之一的映射按该键的值同样方式计数，其他映射算 1 条记录。不做任何解析或校验：含两条序列的 FASTA 会得到 4（标题行和序列行都算），一行以制表符分隔的内容无论是什么都算 1 条记录；`None` 算 1 条记录，空字符串为 0。

- **API：** `dnakit.haplotypes.count(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本、`bytes`、已存在文件的 `Path`（读取前 1,000,000 字节）、列表或映射（见处理方法）。`options` 可以传入但不使用。

- **输出：** 标准封装：`status`（`"ok"`）、`module`（`"haplotypes"`）、`operation`（`"count"`）、`provenance`（`implementation` 为 `"local_candidate_adapter"`，另有 `backend`）、`input_summary`（通用文本摘要：`record_count` 为非空行数，另有 `text_length`、`sequence_bases`（字母个数）、`lengths`、`sha256`）、`warnings`（空）和 `result`。`result` 是整数，即记录数。

- **示例代码：**

```python
import json
from dnakit.haplotypes import count

result = count("H1\tACGT\nH2\tACGA\nH3\tTCGA\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "count",
  "provenance": {
    "implementation": "local_candidate_adapter",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 24,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "27781814e3c6429e6ce9b19cf327192ea6a3cb3e4364cf11bac5c2b35c79e069"
  },
  "warnings": [],
  "status": "ok",
  "result": 3
}
```

## 7) detect_compound_heterozygotes

- **作用：** 找出在两条亲本单倍型上都带有杂合变异的基因（复合杂合子，变异处于反式 trans）。

- **处理方法：** 只使用孩子为杂合的位点。三联体模式（给出父母）：每个杂合的孩子等位基因的来源，由孩子等位基因到（父，母）的唯一孟德尔分配决定；父源一侧的非参考等位基因为父源变异，母源一侧为母源变异；没有唯一分配的位点（例如双亲都杂合）列为 `undetermined`。定相模式（不给父母）：由孩子已定相的调用（如 `0|1`）区分单倍型 1 与 2；未定相的杂合调用会报错。某基因在两条单倍型上各至少有 `min_variants_per_haplotype` 个非参考变异时被报告。不按变异后果、频率或质量过滤。

- **API：** `dnakit.haplotypes.detect_compound_heterozygotes(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`{"genes": [每个位点所属基因], "child": [...]}`，再加 `"father"` 和 `"mother"`（三联体模式），或不给父母但孩子的调用已定相如 `"0|1"`（定相模式）。基因型调用可写作 `"a/b"` 或 `"a|b"`（等位基因为整数）、数对 `(a, b)`、alt 等位基因剂量 0/1/2，缺失用 `None` 或 `"./."`。 `options`（可选）：`min_variants_per_haplotype`（整数，默认 1）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`mode`（`"trio"` 或 `"phased_child"`）、`haplotype_labels`、`genes_with_het_variants`、`compound_heterozygous_genes`（每项含 `gene` 以及位点序号列表 `paternal`、`maternal`、`undetermined`；定相模式下 `paternal` / `maternal` 分别表示单倍型 1 / 2）、`compound_heterozygous_count`、`undetermined_origin_sites` 以及说明字符串 `rule`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import detect_compound_heterozygotes

trio = {
    "genes": ["G1", "G1", "G2"],
    "child": ["0/1", "0/1", "0/1"],
    "father": ["0/1", "0/0", "0/1"],
    "mother": ["0/0", "0/1", "0/1"],
}
result = detect_compound_heterozygotes(trio)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "detect_compound_heterozygotes",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 127,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "0d3afcdd898ed1cdbda697dc7b7022a7b3492b561ceb8e9561c8618bd1e12284"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "mode": "trio",
    "haplotype_labels": "paternal/maternal",
    "genes_with_het_variants": 2,
    "compound_heterozygous_genes": [
      {
        "gene": "G1",
        "paternal": [0],
        "maternal": [1],
        "undetermined": []
      }
    ],
    "compound_heterozygous_count": 1,
    "undetermined_origin_sites": 1,
    "rule": "compound heterozygote = a gene with a heterozygous variant on each of the two parental haplotypes (in trans); trio mode assigns origin when only one parent carries the alt allele, phased mode uses the phase of the child's calls"
  }
}
```

## 8) detect_denovo_variants

- **作用：** 在三联体中检测候选新生（de novo）变异：孩子带有、而双亲基因型中都没有的等位基因。

- **处理方法：** 当孩子带有一个既不在父亲基因型中、也不在母亲基因型中的等位基因（例如双亲均为 0/0，孩子为 0/1）时，该位点为候选新生变异。其他孟德尔违例（例如孩子基因型中缺少应有的亲本等位基因）单独列在 `other_mendelian_inconsistent_indices` 中。任一调用缺失的位点被跳过。没有测序深度、基因型质量和等位基因平衡等信息，因此不做这类过滤；这只是基因型层面的筛查，而不是新生变异检出工具。

- **API：** `dnakit.haplotypes.detect_denovo_variants(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`{"child": [...], "father": [...], "mother": [...]}`，每个位点一个基因型调用（三个列表等长），可选 `"sites"` 标签列表。基因型调用可写作 `"a/b"` 或 `"a|b"`（等位基因为整数）、数对 `(a, b)`、alt 等位基因剂量 0/1/2，缺失用 `None` 或 `"./."`。 `options`：无。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`evaluated_sites`、`denovo_count`、`denovo`（每项含 `index`、可选的 `site` 标签、`child_genotype` 和 `novel_alleles`）、`other_mendelian_inconsistent_indices` 以及说明字符串 `rule`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import detect_denovo_variants

trio = {
    "child": ["0/1", "0/1", "1/1", "0/0"],
    "father": ["0/0", "0/1", "0/0", "0/0"],
    "mother": ["0/0", "0/0", "0/1", "0/0"],
}
result = detect_denovo_variants(trio)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "detect_denovo_variants",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 119,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "60c25594b0884721304b995108bca1551648ed3b3b4291db0bd284e5c721c7e7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "evaluated_sites": 4,
    "denovo_count": 1,
    "denovo": [
      {
        "index": 0,
        "child_genotype": "0/1",
        "novel_alleles": [1]
      }
    ],
    "other_mendelian_inconsistent_indices": [2],
    "rule": "de novo = the child carries an allele that neither parent's genotype contains; other Mendelian violations (e.g. missing parental allele transmitted) are listed separately"
  }
}
```

## 9) frequency

- **作用：** 统计单倍型的个数，计算频率和 Nei 单倍型多样性。

- **处理方法：** 统计不同单倍型字符串的个数。频率 = 个数 / 染色体总数；单倍型多样性（Nei）= n/(n-1) * (1 - sum p_i^2)，n < 2 时为 null。各行按个数降序、再按单倍型字符串排序。值全为整数的映射被当作 `{单倍型: 个数}`。

- **API：** `dnakit.haplotypes.frequency(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：单倍型字符串列表、`{样本: [h1, h2]}`、`{名称: 单倍型}` 或 `{单倍型: 个数}`。`options`：无。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`total`（染色体数）、`distinct`、`haplotypes`（每行含 `haplotype`、`count`、`frequency`）、`haplotype_diversity` 以及说明字符串 `definition`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import frequency

result = frequency({"s1": ["01", "10"], "s2": ["01", "01"], "s3": ["11", "01"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "frequency",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 60,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d5e0fa610438f0870bd6a844d434431de3eb2c92f255ac64c3012dffc59bd5d1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "total": 6,
    "distinct": 3,
    "haplotypes": [
      {
        "haplotype": "01",
        "count": 4,
        "frequency": 0.6666666666666666
      },
      {
        "haplotype": "10",
        "count": 1,
        "frequency": 0.16666666666666666
      },
      {
        "haplotype": "11",
        "count": 1,
        "frequency": 0.16666666666666666
      }
    ],
    "haplotype_diversity": 0.6,
    "definition": "frequency = count / total chromosomes; haplotype diversity (Nei) = n/(n-1) * (1 - sum p_i^2)"
  }
}
```

## 10) impute

- **作用：** 基因型填充（imputation）尚未实现：它需要已定相的参考单倍型面板和基于隐马尔可夫模型的填充引擎。

- **处理方法：** 填充是用统计（隐马尔可夫）模型根据参考面板推断缺失基因型。输入中既没有面板也没有这样的模型，也不会用启发式方法代替，因此函数不会编造基因型。对任何输入它都返回结构化的 `not_implemented` 封装。请在 DNAKit 之外使用 Beagle、IMPUTE5 或 Minimac4。

- **API：** `dnakit.haplotypes.impute(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：任意内容（只在 `input_summary` 中被概括）。`options`（可选）：被忽略。

- **输出：** 始终返回 `status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（缺少的东西：已定相的参考单倍型面板，以及 Beagle、IMPUTE5 或 Minimac4 等填充后端）和 `message`（说明为何无法仅凭输入计算）。不产生 `result`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import impute

result = impute({"sample1": ["0/1", "./.", "1/1"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "impute",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 34,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "89df9963733624e389004dc0717e7ddaf1bb4794e04960464957ab1395e8e1c2"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "phased reference haplotype panel and an imputation backend (Beagle, IMPUTE5, Minimac4)",
  "message": "Imputation infers missing genotypes from a reference panel with a statistical model; neither is part of the input."
}
```

## 11) infer_parentage

- **作用：** 用孟德尔排除法筛查孩子的候选亲本（亲本与孩子没有共享等位基因的位点）。

- **处理方法：** 对每个候选亲本，比较孩子与候选亲本都有调用的位点，统计二者没有共享等位基因的位点数（例如相反的纯合子 0/0 与 1/1），因为真正的亲本必须向孩子传递一个等位基因。该计数超过 `max_mismatches` 的候选被排除。未被排除的候选中错配最少者为最佳候选。这只是排除法，不是亲子关系或亲缘关系似然，也不使用等位基因频率或基因分型错误模型。没有共同调用位点的候选无法评估：它既不算被排除，也不会列入 `not_excluded` 或 `best_candidates`，而是列入 `not_assessed`，并在警告中点名。

- **API：** `dnakit.haplotypes.infer_parentage(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`{"child": [调用], "candidates": {名称: [调用]}}`，每个位点一个调用，所有列表等长。基因型调用可写作 `"a/b"` 或 `"a|b"`（等位基因为整数）、数对 `(a, b)`、alt 等位基因剂量 0/1/2，缺失用 `None` 或 `"./."`。 `options`（可选）：`max_mismatches`（整数 >= 0，默认 0）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`candidates`（每个候选含 `compared_sites`、`mismatches`、`mismatch_rate`、`excluded`）、`not_excluded`、`not_assessed`（没有共同调用位点的候选名称，已排序）、`best_candidates`、`max_mismatches` 以及说明字符串 `rule`；封装的 `warnings` 会列出无法评估的候选。

- **示例代码：**

```python
import json

from dnakit.haplotypes import infer_parentage

data = {
    "child": ["0/0", "1/1", "0/1"],
    "candidates": {
        "A": ["0/1", "1/1", "0/0"],
        "B": ["1/1", "0/0", "0/1"],
    },
}
result = infer_parentage(data)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "infer_parentage",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 104,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "501c4881df823b8b5ba158820fe055215be69bcb17b1bcb9ba0d68f4d3085860"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "candidates": {
      "A": {
        "compared_sites": 3,
        "mismatches": 0,
        "mismatch_rate": 0.0,
        "excluded": false
      },
      "B": {
        "compared_sites": 3,
        "mismatches": 2,
        "mismatch_rate": 0.6666666666666666,
        "excluded": true
      }
    },
    "not_excluded": ["A"],
    "not_assessed": [],
    "best_candidates": ["A"],
    "max_mismatches": 0,
    "rule": "exclusion method: a candidate parent is excluded when the number of sites where parent and child share no allele (e.g. opposite homozygotes) exceeds max_mismatches; a candidate with no jointly called site cannot be assessed and is listed in not_assessed (never in not_excluded or best_candidates)"
  }
}
```

## 12) phase

- **作用：** 用孟德尔三联体逻辑为孩子的基因型定相（顺序为 父源|母源）。

- **处理方法：** 当孩子的两个等位基因到（父，母）的分配恰有一种与双亲基因型相容时，杂合位点被定相，写作 `父源|母源`。若双亲都杂合（或有多种分配相容），位点保持未定相（`a/b`，状态 `ambiguous`）。孟德尔违例被标出（`mendelian_inconsistent`）；孩子的纯合位点写作 `a|a`；孩子调用缺失则保持 `None`；孩子杂合而某位亲本缺失时状态为 `unresolved_parent_missing`。不做统计（群体）定相和基于测序读段的定相，因此不会在位点之间传播相位。

- **API：** `dnakit.haplotypes.phase(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`{"child": [...], "father": [...], "mother": [...]}`，每个位点一个基因型调用（三个列表等长），可选 `"sites"` 标签列表。基因型调用可写作 `"a/b"` 或 `"a|b"`（等位基因为整数）、数对 `(a, b)`、alt 等位基因剂量 0/1/2，缺失用 `None` 或 `"./."`。 `options`：无。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`phased_child`（每个位点一个调用：已定相 `"a|b"`、未定相 `"a/b"` 或 null）、`status`（每个位点：`phased`、`ambiguous`、`homozygous`、`missing`、`mendelian_inconsistent` 或 `unresolved_parent_missing`）、`phased_het_sites`、`ambiguous_het_sites`、`inconsistent_sites`、`fraction_het_sites_phased`、`convention`（`"paternal|maternal"`）以及说明字符串 `rule`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import phase

trio = {
    "child": ["0/1", "0/1", "0/1", "1/1", "./.", "0/1"],
    "father": ["0/0", "0/1", "0/0", "1/1", "0/0", "0/1"],
    "mother": ["0/1", "0/1", "1/1", "0/1", "0/0", "0/0"],
}
result = phase(trio)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "phase",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 161,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "81ac1fa8bb2f1d1a55174a56e17a30af2f26a94ebd9aa8ead62f94a35869ad8f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "phased_child": ["0|1", "0/1", "0|1", "1|1", null, "1|0"],
    "status": ["phased", "ambiguous", "phased", "homozygous", "missing", "phased"],
    "phased_het_sites": 3,
    "ambiguous_het_sites": 1,
    "inconsistent_sites": 0,
    "fraction_het_sites_phased": 0.75,
    "convention": "paternal|maternal",
    "rule": "Mendelian trio phasing: a heterozygous child site is phased when exactly one assignment of its two alleles to (father, mother) is compatible with the parents' genotypes"
  }
}
```

## 13) reconstruct

- **作用：** 用 EM 算法，从未定相的双等位基因型估计单倍型频率和每个个体的单倍型对。

- **处理方法：** 在 Hardy-Weinberg 平衡下的 Excoffier-Slatkin (1995) EM 算法。对每个个体枚举与其基因型相容的所有单倍型对（每个个体最多 12 个杂合位点，候选单倍型总数最多 4,096 个）。从均匀的单倍型频率出发，E 步按 `freq(h1) * freq(h2)`（杂合对再乘 2）给每个单倍型对加权，M 步根据期望的染色体计数重新估计频率，直到某个频率的最大变化小于 `tolerance` 或达到 `max_iterations`。每个个体取后验概率最大的单倍型对。报告的 `log_likelihood` 是数据在最终单倍型频率（最后一次 M 步之后）下的对数似然，因此与返回的 `frequencies` 一致。不支持缺失基因型和多等位位点，也没有重组或连锁模型。

- **API：** `dnakit.haplotypes.reconstruct(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`{样本: [每个位点的调用]}`，调用为 alt 剂量 0/1/2 或 `"0/1"` 形式的字符串，不允许缺失，所有样本的位点数相同。`options`（可选）：`max_iterations`（整数，默认 1000）、`tolerance`（停止所需的最大频率变化，默认 1e-8）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`frequencies`（估计频率高于 1e-12 的单倍型，按频率降序）、`individuals`（每个样本：最可能的单倍型对 `haplotypes`、其后验概率 `posterior` 以及候选数 `alternatives`）、`iterations`、`converged`、`log_likelihood` 以及说明字符串 `model`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import reconstruct

genotypes = {
    "a": [0, 0],
    "b": [0, 0],
    "c": [2, 2],
    "d": [2, 2],
    "e": [1, 1],
}
result = reconstruct(genotypes)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "reconstruct",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 65,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "0f79d304afa9d889f84bd34456e47910635e840ace689d24058c1d3a3d2b5d85"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "frequencies": [
      {
        "haplotype": "00",
        "frequency": 0.5
      },
      {
        "haplotype": "11",
        "frequency": 0.5
      }
    ],
    "individuals": {
      "a": {
        "haplotypes": ["00", "00"],
        "posterior": 1.0,
        "alternatives": 1
      },
      "b": {
        "haplotypes": ["00", "00"],
        "posterior": 1.0,
        "alternatives": 1
      },
      "c": {
        "haplotypes": ["11", "11"],
        "posterior": 1.0,
        "alternatives": 1
      },
      "d": {
        "haplotypes": ["11", "11"],
        "posterior": 1.0,
        "alternatives": 1
      },
      "e": {
        "haplotypes": ["00", "11"],
        "posterior": 1.0,
        "alternatives": 2
      }
    },
    "iterations": 5,
    "converged": true,
    "log_likelihood": -6.238324625039508,
    "model": "Excoffier-Slatkin (1995) EM for haplotype frequencies from unphased biallelic genotypes under Hardy-Weinberg, uniform initial frequencies; individual haplotypes = maximum posterior pair"
  }
}
```

## 14) select

- **作用：** 按最低频率、最低个数和/或前 N 名，从单倍型样本中选取单倍型。

- **处理方法：** 按 `frequency` 的方式建立频率表（频率 = 个数 / 染色体总数），保留频率 >= `min_frequency` 且个数 >= `min_count` 的单倍型，再取按频率降序排列后的前 `top_n` 行。`covered_frequency` 是保留单倍型的频率之和。不使用连锁不平衡或功能方面的标准。

- **API：** `dnakit.haplotypes.select(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：同 `frequency`（单倍型字符串列表、`{样本: [h1, h2]}` 或 `{单倍型: 个数}`）。`options`（可选）：`min_frequency`（默认 0.05）、`min_count`（整数，默认 1）、`top_n`（整数，默认无）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`selected`（每行含 `haplotype`、`count`、`frequency`，按频率降序）、`selected_count`、`total_distinct`、`covered_frequency`、`min_frequency`、`min_count` 和 `top_n`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import select

sample = ["01"] * 5 + ["10"] * 4 + ["11"]
result = select(sample, options={"min_frequency": 0.15})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "select",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 10,
    "text_length": 60,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e80a1b0d3639709d7057db6f0431259c8919ad8b80b9dce7ee4d93f1f54f67fb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "selected": [
      {
        "haplotype": "01",
        "count": 5,
        "frequency": 0.5
      },
      {
        "haplotype": "10",
        "count": 4,
        "frequency": 0.4
      }
    ],
    "selected_count": 2,
    "total_distinct": 3,
    "covered_frequency": 0.9,
    "min_frequency": 0.15,
    "min_count": 1,
    "top_n": null
  }
}
```

## 15) simulate

- **作用：** 按给定的单倍型频率为每个二倍体个体抽取两条单倍型（随机交配）进行模拟。

- **处理方法：** 每个个体的两条染色体各自独立、有放回地从按给定频率或个数加权（归一化为和为 1）的单倍型中抽取，使用 Python 的 `random.Random(seed)`；相同的 `seed` 得到相同的结果。不模拟突变、重组和群体历史；这不是溯祖模拟器。单倍型字符串只是标签，不会被解释。

- **API：** `dnakit.haplotypes.simulate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`{单倍型: 频率或个数}`（非负数，不能全为 0）。`options`（可选）：`n_individuals`（整数 1-100,000，默认 10）、`seed`（整数，默认 0）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`individuals`（每项含 `id`（`sim_<k>`）和两条 `haplotypes`）、`observed_haplotype_counts`、`input_frequencies`（已归一化）、`n_individuals`、`seed` 以及说明字符串 `model`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import simulate

result = simulate({"01": 0.3, "10": 0.7}, options={"n_individuals": 3, "seed": 7})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "simulate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 22,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f1a0a6b72feeb29090802690be6ea8477c7f021b831c1a860459360b5dd23bfd"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "individuals": [
      {
        "id": "sim_1",
        "haplotypes": ["10", "01"]
      },
      {
        "id": "sim_2",
        "haplotypes": ["10", "01"]
      },
      {
        "id": "sim_3",
        "haplotypes": ["10", "10"]
      }
    ],
    "observed_haplotype_counts": {
      "01": 2,
      "10": 4
    },
    "input_frequencies": {
      "01": 0.3,
      "10": 0.7
    },
    "n_individuals": 3,
    "seed": 7,
    "model": "random mating: each of the 2 chromosomes of an individual is drawn independently from the haplotype frequencies (python random.Random(seed)); no mutation or recombination"
  }
}
```

## 16) tag

- **作用：** 用贪心的 r^2 分箱（Carlson 等 2004）从已定相的 0/1 单倍型中选取标签 SNP。

- **处理方法：** 单态位点和次要等位基因频率低于 `min_maf` 的位点被排除。成对 r^2 由所给染色体的单倍型频率计算。反复选取这样的位点作为标签：它与最多个尚未被标记的位点（含自身）的 r^2 >= `r2_threshold`（并列时取序号最小者），这些位点构成它的箱并被标记，直到所有位点都被标记。箱内每个位点与标签的 r^2 都 >= 阈值，但箱内其他位点之间不一定满足。

- **API：** `dnakit.haplotypes.tag(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：已定相单倍型列表（等长的字符串或 0/1 等位基因列表；至少两条单倍型，最多 2,000 个位点），或 `{名称: 单倍型}`。`options`（可选）：`r2_threshold`（默认 0.8）、`min_maf`（默认 0.0）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`tag_sites`（位点序号）、`bins`（标签位点序号 -> 它所标记的位点）、`tag_count`、`sites_analysed`、`r2_threshold`、`excluded_monomorphic_or_low_maf` 以及说明字符串 `rule`。

- **示例代码：**

```python
import json

from dnakit.haplotypes import tag

result = tag(["0000", "0011", "1100", "1111"], options={"r2_threshold": 0.8})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "tag",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 32,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "095e1dbb5d3cd7bea5f10e204f15a9aeaff2bbb73b167c7c21679c30b16c01b7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "tag_sites": [0, 2],
    "bins": {
      "0": [0, 1],
      "2": [2, 3]
    },
    "tag_count": 2,
    "sites_analysed": 4,
    "r2_threshold": 0.8,
    "excluded_monomorphic_or_low_maf": [],
    "rule": "greedy tagging (Carlson et al. 2004): repeatedly choose the site with r^2 >= threshold to the most still-untagged sites (ties: lowest index), assign those sites to its bin"
  }
}
```

## 17) visualize

- **作用：** 把单倍型画成 SVG 等位基因网格（行为单倍型，列为位点）。

- **处理方法：** 每条单倍型一行并带标签，每个位点一个方块；等位基因 `1` 为深红色，`0` 为浅蓝灰色，其他字符（如 `N`）为灰色。每个方块带有 `<title>` 提示，包含单倍型、位点序号和等位基因。不做重新排序或聚类。最多 200,000 个格子；纯 Python 实现，不依赖绘图库。给出 `output` 时，SVG 同时写入该文件。

- **API：** `dnakit.haplotypes.visualize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：等长的单倍型字符串列表，或 `{名称: 单倍型 | [单倍型]}`。`options`（可选）：`output`（SVG 的输出文件路径，默认不写）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`svg`（SVG 文本）、`haplotype_count`、`site_count` 和 `output_path`（写入的文件，未写入时为 null）。

- **示例代码：**

```python
import json

from dnakit.haplotypes import visualize

result = visualize({"h1": "01", "h2": "10", "h3": "1N"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "haplotypes",
  "operation": "visualize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 36,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f161d5549dd126e6deff5d0160767f75f7ee3854cbc850ab2e9865de4fa0089c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"64\" height=\"50\" viewBox=\"0 0 64 50\">\n<text x=\"4\" y=\"15\" font-size=\"11\">h1</text>\n<rect x=\"28\" y=\"4\" width=\"13\" height=\"13\" fill=\"#dfe7ef\"><title>h1 site 0: 0</title></rect>\n<rect x=\"42\" y=\"4\" w...",
    "haplotype_count": 3,
    "site_count": 2,
    "output_path": null
  }
}
```

（示例输出中的长数组和长文本已用 `...` 缩短。）
