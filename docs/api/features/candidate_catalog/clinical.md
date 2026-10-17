# 临床基因组学

本类收录 **18** 项有文档的函数。其中 12 项根据传入的证据（三体基因型、血浆 read 计数、突变矩阵、ACMG/AMP 标准、VCF 或 ClinVar 风格记录）计算结果；6 项始终返回 `status: "not_implemented"`，因为它们需要参考数据库、训练好的模型、外部检出程序，或本身没有确定含义（`annotate`、`calculate`、`call`、`compare`、`predict`、`type_hla`）。所有可计算的函数都只基于调用者提供的证据，没有一个是经过验证的诊断工具。尤其是：`classify` 把 ACMG/AMP 2015 的组合规则应用于调用者已判定满足的标准（标准本身不会根据数据评估，致病与良性证据相互矛盾时得到意义不明）；`detect_circulating_tumor_dna` 与 `detect_minimal_residual_disease` 默认的背景错误率 0.001 只是占位值，必须针对具体检测实测；`generate_clinical_variant_report` 与 `report` 只是汇总所提供的分类，不能代替临床遗传学家的审核。输入格式错误或缺少必需选项时，函数抛出 `InputFormatError` 或 `ConfigurationError`（带 `code` 与 `hint`），不返回封装。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) analyze_trio_variants

- **作用：** 按遗传模式对父-母-子三体的变异分类：新发突变、孟德尔错误、来自单亲或双亲、双亲携带的纯合，以及复合杂合候选。

- **处理方法：** 多等位位点先拆分；对每个 ALT，把子代二倍体基因型与亲代比较。若子代基因型不能由双亲各提供一个等位基因构成，则该位点为 `de_novo`（子代携带 ALT 而双亲都是不含它的纯合）或 `mendelian_error`（包括子代为参考纯合而一方亲本为 ALT 纯合）。携带 ALT 且一致的位点分为 `homozygous_both_parents_carriers`、`inherited_biparental`、`inherited_paternal` 或 `inherited_maternal`；某亲本基因型缺失而子代携带 ALT 的位点为 `undetermined_missing_parent`。复合杂合候选是同一基因内，一个仅由父亲携带的子代杂合变异与一个仅由母亲携带的变异所成的配对（基因符号取自 `gene` 键或 INFO 的 `GENE`/`SYMBOL`/`GENEINFO`；没有基因的变异不参与）。`require_quality_pass` 为真时，三个样本中任一个不满足 `min_depth`/`min_genotype_quality` 的位点被排除。只做常染色体二倍体逻辑：不模拟 X 连锁和线粒体遗传，也不使用定相信息，因此复合杂合只是候选列表，而非已确认的反式（in trans）结论。

- **API：** `dnakit.clinical.analyze_trio_variants(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）；必须包含三个三体样本。`options`：`child`、`father`、`mother`（样本名；缺省且输入恰有三个样本时按列顺序 child、father、mother 处理并给出警告）、`min_depth`（0）、`min_genotype_quality`（0）、`require_quality_pass`（True）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`category_counts`、`variants`（每项含 `site`、`chrom`、`pos`、`ref`、`alt`、`gene`、`category`、`quality_pass`、`child_gt`、`father_gt`、`mother_gt`）、`compound_heterozygous_candidates`（`gene`、`paternal_variant`、`maternal_variant`）和 `n_filtered_low_quality`。与亲代一致的子代纯参考位点不列出。样本名错误或缺失时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.clinical import analyze_trio_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother",
    "chr1\t100\t.\tA\tG\t50\tPASS\tGENE=GENE1\tGT:DP:GQ\t0/1:30:60\t0/0:30:60\t0/0:30:60",
    "chr1\t200\t.\tC\tT\t50\tPASS\tGENE=GENE2\tGT:DP:GQ\t0/1:30:60\t0/1:30:60\t0/0:30:60",
    "chr1\t300\t.\tG\tA\t50\tPASS\tGENE=GENE2\tGT:DP:GQ\t0/1:30:60\t0/0:30:60\t0/1:30:60",
    "chr1\t400\t.\tT\tC\t50\tPASS\tGENE=GENE3\tGT:DP:GQ\t1/1:30:60\t0/1:30:60\t0/1:30:60",
    "chr1\t500\t.\tA\tC\t50\tPASS\tGENE=GENE4\tGT:DP:GQ\t0/0:30:60\t1/1:30:60\t0/0:30:60",
]) + "\n"

result = analyze_trio_variants(vcf, options={"child": "child", "father": "father", "mother": "mother"})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "analyze_trio_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 7,
    "text_length": 452,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "eb0bb0a2db5f091abacfdcce919bbe3dd13373aaf43148678619612ad63bf6df"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "per-variant inheritance category from child/parent genotypes (after splitting multi-alleli...",
    "category_counts": {
      "de_novo": 1,
      "homozygous_both_parents_carriers": 1,
      "inherited_maternal": 1,
      "inherited_paternal": 1,
      "mendelian_error": 1
    },
    "variants": [
      {
        "site": "chr1:100:A:G",
        "chrom": "chr1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
        "gene": "GENE1",
        "category": "de_novo",
        "quality_pass": true,
        "child_gt": "0/1",
        "father_gt": "0/0",
        "mother_gt": "0/0"
      },
      {
        "site": "chr1:200:C:T",
        "chrom": "chr1",
        "pos": 200,
        "ref": "C",
        "alt": "T",
        "gene": "GENE2",
        "category": "inherited_paternal",
        "quality_pass": true,
        "child_gt": "0/1",
        "father_gt": "0/1",
        "mother_gt": "0/0"
      },
      "..."
    ],
    "compound_heterozygous_candidates": [{"gene": "GENE2", "paternal_variant": "chr1:200:C:T", "maternal_variant": "chr1:300:G:A"}],
    "n_filtered_low_quality": 0
  }
}
```

## 2) annotate

- **作用：** 用 ClinVar、gnomAD 和基因-疾病信息注释临床变异。尚未实现：函数只报告缺少什么（ClinVar、gnomAD 与基因-疾病参考数据库）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。单个变异请使用 `dnakit.databases` / `dnakit.annotation` 的查询函数；`read` 可提取 VCF 中已有的 ClinVar 字段。

- **API：** `dnakit.clinical.annotate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"ClinVar / gnomAD / gene-disease databases (local files or network access)"）与 `message`（"clinical annotation needs reference databases that are not bundled."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.clinical import annotate

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "clinical",
  "operation": "annotate",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "ClinVar / gnomAD / gene-disease databases (local files or network access)",
  "message": "clinical annotation needs reference databases that are not bundled."
}
```

## 3) build_cancer_evolution_tree

- **作用：** 由二值的突变 x 样本（或细胞）矩阵构建体细胞突变的有根完美系统发育（演化）树。

- **处理方法：** 每个突变用携带它的样本/细胞集合表示。集合相同的突变共享一个节点；节点的父节点是样本集合严格更大且最小的节点，根为种系。只有当各样本集合两两嵌套或不相交时才存在这样的树；若两个突变的集合部分重叠（违反四配子条件），结果报告 `perfect_phylogeny: false` 并列出冲突对，不构建树。因此它假设无限位点模型（无回复突变或重复突变）和无错误的调用；没有似然模型、容错或克隆频率估计。在所有样本中都不出现的突变会被忽略并给出警告。

- **API：** `dnakit.clinical.build_cancer_evolution_tree(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'mutations': {'TP53': [1, 1, 0], ...}, 'samples': [...]}`（每个样本 0/1 表示有无；也可为 `突变 -> {样本: 0/1}`；至多 5000 个突变，所有向量等长；`samples` 缺省为 `sample_1`……）。`options`：无。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`perfect_phylogeny`、`violations`（`mutation_1`、`mutation_2`、`reason`）、`n_mutations`，存在树时还有 `nodes`（`id`、`label`、`mutations`、`samples`、`parent`、`children`、`depth`；第一个是 `root` = `germline`）、`newick` 和 `sample_to_node`（每个样本所属的最深节点）。存在冲突时 `tree` 为 `null`，没有 `nodes`。

- **示例代码：**

```python
import json

from dnakit.clinical import build_cancer_evolution_tree

result = build_cancer_evolution_tree({
    "samples": ["S1", "S2", "S3", "S4"],
    "mutations": {
        "TP53": [1, 1, 1, 1],
        "KRAS": [1, 1, 0, 0],
        "PIK3CA": [1, 0, 0, 0],
        "APC": [0, 0, 1, 1],
    },
})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "build_cancer_evolution_tree",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 141,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d564184bd9487dc45696fc6e0f66a5f2e7510cc777268e6351a60433e7f8d170"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "perfect-phylogeny (infinite-sites) tree from a binary mutation x sample matrix: mutations...",
    "perfect_phylogeny": true,
    "violations": [],
    "n_mutations": 4,
    "nodes": [
      {
        "id": "root",
        "label": "germline",
        "mutations": [],
        "samples": ["S1", "S2", "..."],
        "parent": null,
        "children": ["n1"],
        "depth": 0
      },
      {
        "id": "n1",
        "label": "TP53",
        "mutations": ["TP53"],
        "samples": ["S1", "S2", "..."],
        "parent": "root",
        "children": ["n2", "n3"],
        "depth": 1
      },
      "..."
    ],
    "newick": "(((PIK3CA)KRAS,APC)TP53)germline;",
    "sample_to_node": {"S1": "n4", "S2": "n2", "S3": "n3", "S4": "n3"}
  }
}
```

## 4) calculate

- **作用：** 计算一个未指明的临床量。尚未实现：函数只报告缺少什么（具体的临床计算；泛称 `calculate` 没有确定含义）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。请改用具体函数（`classify`、`qc`、`detect_de_novo_variants` 等）。

- **API：** `dnakit.clinical.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"a specific clinical calculation (this generic name has no defined meaning)"）与 `message`（"'calculate' is too unspecific to have a defined meaning; call a specific clinical function."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.clinical import calculate

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = calculate(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "clinical",
  "operation": "calculate",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a specific clinical calculation (this generic name has no defined meaning)",
  "message": "'calculate' is too unspecific to have a defined meaning; call a specific clinical function."
}
```

## 5) call

- **作用：** 从测序 reads 中检出临床变异。尚未实现：函数只报告缺少什么（外部变异检出程序与比对后的 reads）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。请使用 `dnakit.variants.call_freebayes`、`call_gatk_haplotypecaller` 或 `call_deepvariant`，这些函数在程序已安装时调用外部检出程序。

- **API：** `dnakit.clinical.call(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"variant caller binary and aligned reads"）与 `message`（"variant calling needs an external caller and aligned reads."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.clinical import call

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = call(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "clinical",
  "operation": "call",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "variant caller binary and aligned reads",
  "message": "variant calling needs an external caller and aligned reads."
}
```

## 6) classify

- **作用：** 按 ACMG/AMP 2015 规则把证据标准组合成五级变异分类（致病、可能致病、意义不明、可能良性、良性）。

- **处理方法：** 适用于该变异的标准由调用者提供，函数不会根据数据去评估这些标准。每个代码（PVS1；PS1-4；PM1-6；PP1-5；BA1；BS1-4；BP1-7）都有默认强度（极强、强、中等、支持、独立），可用后缀修改（如 `PM2_Supporting`、`PVS1_Strong`）；同一代码以不同强度重复给出、未知代码或该方向上未定义的强度都会报错。组合采用 Richards 等 2015 年表 5 的规则：致病（PVS1 + >=1 个 PS；PVS1 + >=2 个 PM；PVS1 + PM + PP；PVS1 + >=2 个 PP；>=2 个 PS；PS + >=3 个 PM；PS + 2 个 PM + >=2 个 PP；PS + PM + >=4 个 PP）、可能致病（PVS1 + PM；PS + 1-2 个 PM；PS + >=2 个 PP；>=3 个 PM；2 个 PM + >=2 个 PP；PM + >=4 个 PP）、良性（BA1；>=2 个 BS）、可能良性（BS + BP；>=2 个 BP）。致病一侧与良性一侧分别判定：只有一侧得出分类时即为结果；两侧都没有得出分类时为意义不明；两侧都得出分类（证据互相矛盾，例如 PVS1 + PS1 同时有两个 BS）时，结果为意义不明（VUS），`rule` 中同时写出冲突的两条规则（这是 2015 指南对证据冲突的默认处理）。只实现了这套离散的计数规则：后来的贝叶斯点数框架（Tavtigian 等 2018/2020）、ClinGen 的基因或疾病特异规范以及标准的自动评估均未实现，结果也不是临床解读。

- **API：** `dnakit.clinical.classify(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：已满足的标准，可为代码列表（`['PVS1', 'PM2_Supporting']`，元素也可为 `{'code', 'strength'}`）、`{'criteria': [...], 'id': ...}`，或 `{'variants': [{'id': ..., 'criteria': [...]}, ...]}`（多个变异）。`options`：无。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `framework`；单个变异时还有 `classification`、`rule`（所适用的规则，或冲突的说明，或 "criteria do not meet any pathogenic or benign combination"）、`evidence_counts`（按方向和强度的计数）、`criteria`（规范化后的 `code`、`strength`、`direction`）以及给定时的 `id`；多个变异时 `results` 为每个变异一项。标准列表为空时给出意义不明并附警告。无效代码会抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.clinical import classify

result = classify({
    "variants": [
        {"id": "v1", "criteria": ["PVS1", "PS1"]},
        {"id": "v2", "criteria": ["PM2_Supporting", "PP3"]},
        {"id": "v3", "criteria": ["PVS1", "PS1", "BS1", "BS2"]},  # pathogenic and benign rules both met
    ]
})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "classify",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 167,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2e728bced9ea44abed7aabb390830829cd366055831cd1a3492f622e97162582"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "framework": "ACMG/AMP 2015 (Richards et al., Genet Med 17:405-424) Table 5 combining rules applied to c...",
    "results": [
      {
        "id": "v1",
        "classification": "Pathogenic",
        "rule": "Pathogenic (i)(a): 1 Very strong + >=1 Strong",
        "evidence_counts": {
          "pathogenic_very_strong": 1,
          "pathogenic_strong": 1,
          "pathogenic_moderate": 0,
          "pathogenic_supporting": 0,
          "benign_stand_alone": 0,
          "benign_strong": 0,
          "benign_supporting": 0
        },
        "criteria": [
          {"code": "PVS1", "strength": "very_strong", "direction": "pathogenic"},
          {"code": "PS1", "strength": "strong", "direction": "pathogenic"}
        ]
      },
      {
        "id": "v2",
        "classification": "Uncertain significance",
        "rule": "criteria do not meet any pathogenic or benign combination",
        "evidence_counts": {
          "pathogenic_very_strong": 0,
          "pathogenic_strong": 0,
          "pathogenic_moderate": 0,
          "pathogenic_supporting": 2,
          "benign_stand_alone": 0,
          "benign_strong": 0,
          "benign_supporting": 0
        },
        "criteria": [
          {"code": "PM2", "strength": "supporting", "direction": "pathogenic"},
          {"code": "PP3", "strength": "supporting", "direction": "pathogenic"}
        ]
      },
      {
        "id": "v3",
        "classification": "Uncertain significance",
        "rule": "contradictory evidence: Pathogenic (i)(a): 1 Very strong + >=1 Strong vs Benign (ii): >=2 Strong",
        "evidence_counts": {
          "pathogenic_very_strong": 1,
          "pathogenic_strong": 1,
          "pathogenic_moderate": 0,
          "pathogenic_supporting": 0,
          "benign_stand_alone": 0,
          "benign_strong": 2,
          "benign_supporting": 0
        },
        "criteria": [
          {"code": "PVS1", "strength": "very_strong", "direction": "pathogenic"},
          {"code": "PS1", "strength": "strong", "direction": "pathogenic"},
          "..."
        ]
      }
    ]
  }
}
```

## 7) compare

- **作用：** 比较临床数据集。尚未实现：函数只报告缺少什么（具体的比较定义；泛称 `compare` 在此没有确定含义）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。比较变异调用集请使用 `dnakit.variants.compare_vcfs`。

- **API：** `dnakit.clinical.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"a specific comparison definition (this generic name has no defined meaning)"）与 `message`（"'compare' is too unspecific here; use dnakit.variants.compare_vcfs for variant call sets."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.clinical import compare

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = compare(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "clinical",
  "operation": "compare",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a specific comparison definition (this generic name has no defined meaning)",
  "message": "'compare' is too unspecific here; use dnakit.variants.compare_vcfs for variant call sets."
}
```

## 8) detect_circulating_tumor_dna

- **作用：** 在追踪的肿瘤变异位点上检验血浆中是否存在循环肿瘤 DNA（ctDNA；基于肿瘤信息、使用合并二项检验）。

- **处理方法：** 把追踪的肿瘤变异位点上血浆的 ALT reads 数与深度合并，用单侧二项检验 P(X >= alt | depth, error) 检验合并 ALT 数相对背景错误率是否过高。检出 = 合并 p <= `alpha` 且至少有 `min_positive_sites` 个位点单独为阳性（ALT >= `min_alt_reads` 且 p <= `site_alpha`）。肿瘤分数估计 = `vaf_to_fraction_factor` * max(0, 合并 VAF - 背景错误率)，上限 1（系数 2 对应克隆性杂合、拷贝数中性的 SNV）。默认的 `background_error_rate` = 0.001 只是占位值：必须针对具体检测实测（例如用配对的健康人血浆，或追踪位置上的非参考碱基）后再传入，因为检出结论和 p 值直接取决于它。各 read 假定独立且共用一个错误率（可按位点提供 `error_rate`），没有 UMI/duplex 纠错、链或上下文模型，因此这只是统计示意，不是经过验证的临床检测。

- **API：** `dnakit.clinical.detect_circulating_tumor_dna(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：追踪位点 `[{'alt_count', 'depth'[, 'id', 'error_rate']}, ...]` 或 `{'sites': [...]}`（血浆 read 计数，整数且 0 <= alt_count <= depth），或带血浆样本 `FORMAT/AD` 的 VCF/记录（`options.sample`，默认第一个样本；多等位位点先拆分）。`options`：`background_error_rate`（0.001；需实测的占位值）、`alpha`（0.01）、`site_alpha`（0.01）、`min_alt_reads`（2）、`min_positive_sites`（1）、`vaf_to_fraction_factor`（2.0）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`parameters`、`detected`、`n_sites`、`n_positive_sites`、`pooled_alt_reads`、`pooled_depth`、`pooled_vaf`、`pooled_p_value`、`estimated_tumor_fraction` 和 `sites`（每个位点：`id`、`alt_count`、`depth`、`error_rate`、`vaf`、`p_value`、`positive`）。

- **示例代码：**

```python
import json

from dnakit.clinical import detect_circulating_tumor_dna

# plasma reads at 4 tumour-informed sites: alt_count out of depth
sites = [
    {"id": "TP53:R175H", "alt_count": 12, "depth": 5000},
    {"id": "KRAS:G12D", "alt_count": 9, "depth": 4800},
    {"id": "APC:R1450*", "alt_count": 1, "depth": 5200},
    {"id": "PIK3CA:E545K", "alt_count": 10, "depth": 4900},
]

result = detect_circulating_tumor_dna({"sites": sites}, options={"background_error_rate": 0.0005})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "detect_circulating_tumor_dna",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 226,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d40e8a937f2b6f0485bbf2cd42cad1a6a7f68bf85efc96365e6faad9baa1ebba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "ctDNA detection: tumour-informed detection. The plasma ALT reads and depths at the tracked...",
    "parameters": {"background_error_rate": 0.0005, "alpha": 0.01, "site_alpha": 0.01, "min_alt_reads": 2, "min_positive_sites": 1},
    "detected": true,
    "n_sites": 4,
    "n_positive_sites": 3,
    "pooled_alt_reads": 32,
    "pooled_depth": 19900,
    "pooled_vaf": 0.00160804,
    "pooled_p_value": 2.17455e-08,
    "estimated_tumor_fraction": 0.00221608,
    "sites": [
      {
        "id": "TP53:R175H",
        "alt_count": 12,
        "depth": 5000,
        "error_rate": 0.0005,
        "vaf": 0.0024,
        "p_value": 1.24946e-05,
        "positive": true
      },
      {
        "id": "KRAS:G12D",
        "alt_count": 9,
        "depth": 4800,
        "error_rate": 0.0005,
        "vaf": 0.001875,
        "p_value": 0.000858515,
        "positive": true
      },
      "..."
    ]
  }
}
```

## 9) detect_de_novo_variants

- **作用：** 根据基因型、深度和等位基因深度检测父-母-子三体中的新发变异候选。

- **处理方法：** 候选位点是（多等位位点拆分后）子代携带 ALT 等位基因、而双亲均被调用为不含它的纯合的位点。随后依次过滤：子代深度 >= `min_child_depth`，亲代深度 >= `min_parent_depth`，三个样本的 GQ 都 >= `min_genotype_quality`，并且在存在 `FORMAT/AD` 时，子代 ALT 比例 >= `min_child_vaf`、亲代 ALT reads <= `max_parent_alt_reads`。每个被剔除的候选只在 `rejected_by` 中计一次，计入它在固定顺序中第一个未通过的检查（先是定义候选的基因型检查，然后依次为：子代、父亲、母亲的深度；子代、父亲、母亲的 GQ；子代 ALT 比例；父亲、母亲的 ALT reads），因此原因是确定的（如 `child_low_depth`、`child_low_vaf`、`father_alt_reads`）。这些固定阈值是通用默认值，函数不了解低可比对区域、亲代嵌合或测序错误谱，因此输出只是需要正交验证的候选列表，而不是经过验证的新发突变判定。

- **API：** `dnakit.clinical.detect_de_novo_variants(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）；必须包含三个三体样本。`options`：`child`、`father`、`mother`（样本名；恰有三个样本时按 child、father、mother 顺序处理并给出警告）、`min_child_depth`（10）、`min_parent_depth`（10）、`min_genotype_quality`（20）、`min_child_vaf`（0.25，存在 AD 时应用）、`max_parent_alt_reads`（0，存在 AD 时应用）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`parameters`、`n_candidates_before_filters`、`rejected_by`、`n_de_novo` 和 `variants`（每项含 `site`、`chrom`、`pos`、`ref`、`alt`、`gene`、`category` = `de_novo`、`quality_pass`、`child_gt`、`father_gt`、`mother_gt`、`child_ad`）。

- **示例代码：**

```python
import json

from dnakit.clinical import detect_de_novo_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:35:80:20,15\t0/0:30:60:30,0\t0/0:32:60:32,0",
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:8:80:5,3\t0/0:30:60:30,0\t0/0:32:60:32,0",
    "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:40:80:36,4\t0/0:30:60:30,0\t0/0:32:60:32,0",
    "chr1\t400\t.\tT\tC\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:35:80:20,15\t0/0:30:60:28,2\t0/0:32:60:32,0",
]) + "\n"

result = detect_de_novo_variants(vcf, options={"child": "child", "father": "father", "mother": "mother"})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "detect_de_novo_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 415,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c5fffaae813fd43f7d1b4c5850ada93d33e140429b1303180591f6048a810982"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "de novo candidates = child carries the ALT allele while both parents are called homozygous...",
    "parameters": {
      "min_child_depth": 10.0,
      "min_parent_depth": 10.0,
      "min_genotype_quality": 20.0,
      "min_child_vaf": 0.25,
      "max_parent_alt_reads": 0
    },
    "n_candidates_before_filters": 4,
    "rejected_by": {"child_low_depth": 1, "child_low_vaf": 1, "father_alt_reads": 1},
    "n_de_novo": 1,
    "variants": [
      {
        "site": "chr1:100:A:G",
        "chrom": "chr1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
        "gene": null,
        "category": "de_novo",
        "quality_pass": true,
        "child_gt": "0/1",
        "father_gt": "0/0",
        "mother_gt": "0/0",
        "child_ad": [20, 15]
      }
    ]
  }
}
```

## 10) detect_minimal_residual_disease

- **作用：** 根据追踪的肿瘤变异位点上的血浆 read 计数，用更严格的肿瘤信息引导合并二项检验判断微小残留病灶（MRD）。

- **处理方法：** 统计方法与 `detect_circulating_tumor_dna` 相同，但默认值更严格（`alpha` 为 0.001，且至少有两个位点单独为阳性）。把追踪的肿瘤变异位点上血浆的 ALT reads 数与深度合并，用单侧二项检验 P(X >= alt | depth, error) 检验合并 ALT 数相对背景错误率是否过高。检出 = 合并 p <= `alpha` 且至少有 `min_positive_sites` 个位点单独为阳性（ALT >= `min_alt_reads` 且 p <= `site_alpha`）。肿瘤分数估计 = `vaf_to_fraction_factor` * max(0, 合并 VAF - 背景错误率)，上限 1（系数 2 对应克隆性杂合、拷贝数中性的 SNV）。默认的 `background_error_rate` = 0.001 只是占位值：必须针对具体检测实测（例如用配对的健康人血浆，或追踪位置上的非参考碱基）后再传入，因为检出结论和 p 值直接取决于它。各 read 假定独立且共用一个错误率（可按位点提供 `error_rate`），没有 UMI/duplex 纠错、链或上下文模型，因此这只是统计示意，不是经过验证的临床检测。

- **API：** `dnakit.clinical.detect_minimal_residual_disease(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：追踪位点 `[{'alt_count', 'depth'[, 'id', 'error_rate']}, ...]` 或 `{'sites': [...]}`（血浆 read 计数，整数且 0 <= alt_count <= depth），或带血浆样本 `FORMAT/AD` 的 VCF/记录（`options.sample`，默认第一个样本；多等位位点先拆分）。`options`：`background_error_rate`（0.001；需实测的占位值）、`alpha`（0.001）、`site_alpha`（0.01）、`min_alt_reads`（2）、`min_positive_sites`（2）、`vaf_to_fraction_factor`（2.0）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`parameters`、`detected`、`n_sites`、`n_positive_sites`、`pooled_alt_reads`、`pooled_depth`、`pooled_vaf`、`pooled_p_value`、`estimated_tumor_fraction` 和 `sites`（每个位点：`id`、`alt_count`、`depth`、`error_rate`、`vaf`、`p_value`、`positive`）。

- **示例代码：**

```python
import json

from dnakit.clinical import detect_minimal_residual_disease

# plasma reads at 4 tumour-informed sites: alt_count out of depth
sites = [
    {"id": "TP53:R175H", "alt_count": 12, "depth": 5000},
    {"id": "KRAS:G12D", "alt_count": 9, "depth": 4800},
    {"id": "APC:R1450*", "alt_count": 1, "depth": 5200},
    {"id": "PIK3CA:E545K", "alt_count": 10, "depth": 4900},
]

result = detect_minimal_residual_disease({"sites": sites}, options={"background_error_rate": 0.0005})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "detect_minimal_residual_disease",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 226,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d40e8a937f2b6f0485bbf2cd42cad1a6a7f68bf85efc96365e6faad9baa1ebba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "MRD detection: tumour-informed detection. The plasma ALT reads and depths at the tracked t...",
    "parameters": {"background_error_rate": 0.0005, "alpha": 0.001, "site_alpha": 0.01, "min_alt_reads": 2, "min_positive_sites": 2},
    "detected": true,
    "n_sites": 4,
    "n_positive_sites": 3,
    "pooled_alt_reads": 32,
    "pooled_depth": 19900,
    "pooled_vaf": 0.00160804,
    "pooled_p_value": 2.17455e-08,
    "estimated_tumor_fraction": 0.00221608,
    "sites": [
      {
        "id": "TP53:R175H",
        "alt_count": 12,
        "depth": 5000,
        "error_rate": 0.0005,
        "vaf": 0.0024,
        "p_value": 1.24946e-05,
        "positive": true
      },
      {
        "id": "KRAS:G12D",
        "alt_count": 9,
        "depth": 4800,
        "error_rate": 0.0005,
        "vaf": 0.001875,
        "p_value": 0.000858515,
        "positive": true
      },
      "..."
    ]
  }
}
```

## 11) generate_clinical_variant_report

- **作用：** 根据已分类的变异生成临床变异报告（Markdown 或纯文本）。

- **处理方法：** 每个变异被赋予一个层级：给出的 `classification`（`Pathogenic`、`Likely pathogenic`、`Uncertain significance`、`Likely benign`、`Benign`，或 `P`/`LP`/`VUS`/`LB`/`B`）原样使用；只给 `criteria` 时按 `classify` 的 ACMG/AMP 2015 规则组合（因此互相矛盾的标准得到意义不明）；否则为 `Not classified`。VCF 输入时，基因、ClinVar 术语和后果从记录中读取并作为证据显示，但 ClinVar 意义不会被转换成分类，所以除非另行提供分类，这些变异为 `Not classified`；`hgvs_accession`（及 `reference`）只用于添加 HGVS `g.` 字符串。变异按层级和基因排序，报告列出可报告结果（致病与可能致病）、意义不明或未分类的变异（除非 `include_vus` 为 False）、未列出的良性/可能良性变异数量，以及固定的方法与局限说明。它只汇总所提供的证据：不查询数据库，也不重新评估标准，因此文本不是经过验证的诊断报告，必须由合格的临床遗传学家审核。

- **API：** `dnakit.clinical.generate_clinical_variant_report(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'variants': [{'gene', 'hgvs', 'zygosity', 'classification' 或 'criteria', 'clinvar', ...}], 'sample': ..., 'indication': ...}`（`sample`、`patient`、`patient_id`、`indication`、`date` 会印在报告头部）、这类行的列表，或 VCF/记录（基因、ClinVar 与后果字段从记录中读取）。`options`：`title`、`format`（默认 `'markdown'`，另有 `'text'`）、`include_vus`（True）、用于为 VCF 输入生成 HGVS `g.` 字符串的 `hgvs_accession`（字符串或 `{contig: accession}`）与 `reference`、`output`（写出报告的路径）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `format`、`report`（完整报告文本）、`n_variants`、`classification_counts`、`n_reportable`、`variants`（带已确定 `classification` 的各行，由标准推出时还有 `classification_basis`）和 `output_path`。分类标签未知或变异列表为空时抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.clinical import generate_clinical_variant_report

result = generate_clinical_variant_report(
    {
        "sample": "P001",
        "indication": "hereditary cancer panel",
        "variants": [
            {"gene": "BRCA1", "hgvs": "NM_007294.4:c.68_69del", "zygosity": "heterozygous", "criteria": ["PVS1", "PS3"]},
            {"gene": "MSH6", "hgvs": "NM_000179.3:c.1A>G", "zygosity": "heterozygous", "classification": "VUS"},
            {"gene": "TTN", "hgvs": "NM_001267550.2:c.100C>T", "zygosity": "heterozygous", "classification": "LB"},
        ],
    },
    options={"title": "Hereditary cancer panel"},
)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "generate_clinical_variant_report",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 388,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3a1fb883b0a897f14eb0628ee788f198be4b6e220d8bea57bdbd15627135a355"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "markdown",
    "report": "# Hereditary cancer panel\n\nSample: P001\nIndication: hereditary cancer panel\nVariants evalu...",
    "n_variants": 3,
    "classification_counts": {"Likely benign": 1, "Pathogenic": 1, "Uncertain significance": 1},
    "n_reportable": 1,
    "variants": [
      {
        "gene": "BRCA1",
        "hgvs": "NM_007294.4:c.68_69del",
        "zygosity": "heterozygous",
        "criteria": ["PVS1", "PS3"],
        "classification_basis": "ACMG/AMP 2015 combination of supplied criteria",
        "classification": "Pathogenic"
      },
      {
        "gene": "MSH6",
        "hgvs": "NM_000179.3:c.1A>G",
        "zygosity": "heterozygous",
        "classification": "Uncertain significance"
      },
      {"gene": "TTN", "hgvs": "NM_001267550.2:c.100C>T", "zygosity": "heterozygous", "classification": "Likely benign"}
    ],
    "output_path": null
  }
}
```

完整的 `result["report"]` 字符串如下：

    # Hereditary cancer panel

    Sample: P001
    Indication: hereditary cancer panel
    Variants evaluated: 3 (1 Likely benign, 1 Pathogenic, 1 Uncertain significance)

    ## Reportable findings (Pathogenic / Likely pathogenic)

    | Gene | Variant | Zygosity | Classification | Evidence |
    |---|---|---|---|---|
    | BRCA1 | NM_007294.4:c.68_69del | heterozygous | Pathogenic | PVS1, PS3 |

    ## Variants of uncertain significance / not classified

    | Gene | Variant | Zygosity | Classification | Evidence |
    |---|---|---|---|---|
    | MSH6 | NM_000179.3:c.1A>G | heterozygous | Uncertain significance | - |

    ## Not reported

    1 benign / likely benign variants and 0 uncertain variants are not listed.

    ## Method and limitations

    Classifications are taken from the input or computed from the supplied ACMG/AMP criteria with the 2015 combining rules (Richards et al., Genet Med 2015). The criteria were not re-evaluated against databases. This automatically generated summary is not a validated diagnostic report and must be reviewed by a qualified clinical geneticist.

## 12) infer_tumor_clones

- **作用：** 用二项混合模型（EM，克隆数由 BIC 选择）对突变的 ALT 比例聚类，推断肿瘤克隆。

- **处理方法：** 把（alt_count, depth）数据建模为 k 个二项分布的混合，用 EM 从确定性的分位数初值开始拟合；比较 k = 1..`max_clones`（不超过突变数）的 BIC = -2 log L + (2k - 1) log n，取 BIC 最小者。突变被分配到后验责任最大的成分。各簇按 ALT 比例排序，最高者标为克隆性（clonal），其余为亚克隆（subclonal）。每个突变的癌细胞比例为 CCF = VAF * (purity * CN + 2 * (1 - purity)) / (purity * multiplicity)，CN 默认为 2、multiplicity 默认为 1（可逐突变给出），每个簇报告平均 CCF。不推断拷贝数变化（须由用户提供），不模拟克隆系统发育/先后顺序；突变数少或深度低时 BIC 的选择不可靠；它不能替代 PyClone、SciClone 等工具。

- **API：** `dnakit.clinical.infer_tumor_clones(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：突变列表 `[{'alt_count', 'depth'[, 'id', 'copy_number', 'multiplicity']}, ...]` 或 `{'sites': [...]}`，或带 `FORMAT/AD` 的 VCF/记录（`options.sample`）；至多 20,000 个突变。`options`：`purity`（1.0，范围 0.01-1）、`max_clones`（4，至多 8）、`copy_number`（默认 CN 2）、`sample`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`purity`、`n_mutations`、`n_clones`、`model_selection`（每个尝试的 k 的 `k`、`log_likelihood`、`bic`）、`clones`（每个簇：`cluster`、`role`、`cluster_vaf`、`pooled_vaf`、`weight`、`n_mutations`、`mean_ccf`、`mutations`）和 `assignments`（`id`、`cluster`）。

- **示例代码：**

```python
import json

from dnakit.clinical import infer_tumor_clones

# (alt_count, depth) per mutation: a clonal group near VAF 0.40 and a subclonal group near 0.15
data = [(40, 100), (38, 100), (42, 100), (41, 100), (15, 100), (14, 100), (16, 100)]
sites = [{"id": f"m{i + 1}", "alt_count": a, "depth": d} for i, (a, d) in enumerate(data)]

result = infer_tumor_clones({"sites": sites}, options={"purity": 0.8})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "infer_tumor_clones",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 326,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b9fd33439bd546877ffc7c62a5baa0e90ee2b2d603c91c538205fc2ebcd69939"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "binomial mixture model of ALT/depth counts fitted by EM (deterministic quantile start), nu...",
    "purity": 0.8,
    "n_mutations": 7,
    "n_clones": 2,
    "model_selection": [{"k": 1, "log_likelihood": -424.161, "bic": 850.268}, {"k": 2, "log_likelihood": -401.198, "bic": 808.234}, "..."],
    "clones": [
      {
        "cluster_vaf": 0.4025,
        "pooled_vaf": 0.4025,
        "weight": 0.571429,
        "n_mutations": 4,
        "mean_ccf": 1.00625,
        "mutations": ["m1", "m2", "..."],
        "cluster": 1,
        "role": "clonal (highest ALT fraction)"
      },
      {
        "cluster_vaf": 0.15,
        "pooled_vaf": 0.15,
        "weight": 0.428571,
        "n_mutations": 3,
        "mean_ccf": 0.375,
        "mutations": ["m5", "m6", "m7"],
        "cluster": 2,
        "role": "subclonal"
      }
    ],
    "assignments": [{"id": "m1", "cluster": 1}, {"id": "m2", "cluster": 1}, "..."]
  }
}
```

## 13) predict

- **作用：** 用训练好的模型预测临床结局或变异效应。尚未实现：函数只报告缺少什么（训练好的临床/致病性预测器及其权重）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。函数不会加载或运行任何模型。

- **API：** `dnakit.clinical.predict(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"trained clinical/pathogenicity prediction model"）与 `message`（"prediction needs trained model weights that are not available."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.clinical import predict

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = predict(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "clinical",
  "operation": "predict",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "trained clinical/pathogenicity prediction model",
  "message": "prediction needs trained model weights that are not available."
}
```

## 14) prioritize

- **作用：** 用透明的加和得分（后果、ClinVar、稀有度、CADD/REVEL、基因面板）对临床候选变异排序。

- **处理方法：** 评分与 `dnakit.variants.prioritize_candidate_variants` 相同。读取记录上已有的注释（后果术语、`clinvar`/`CLNSIG`、`gnomad_af`、`cadd`、`revel`、`gene`），不做任何查询。得分 = 后果严重度（固定的 Sequence Ontology 表，0-1）+ ClinVar（致病 1.0，可能致病 0.8，意义不明 0.1，良性 -1.0）+ 稀有度（人群 AF 缺失或 <= 1e-4 为 0.5，至 1e-3 为 0.3，至 1e-2 为 0.1）+ 0.5 * max(CADD / 40, REVEL) + 面板基因 0.5。人群 AF 高于 `max_af` 的变异，以及 `exclude_benign` 为真时 ClinVar 良性/可能良性的变异，列入排除项。权重是人为设定的，所以这是启发式排序而不是经过验证的分类器；没有人群频率注释的变异保留完整的稀有度得分（0.5），但会标记 `frequency_unknown: true`，并计入一条汇总警告，从而能与确知稀有的变异区分开。每个得分分量都会返回。

- **API：** `dnakit.clinical.prioritize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。注释从记录的键或 INFO 读取（`consequence`、`clinvar`/`CLNSIG`、`gnomad_af`、`cadd`、`revel`、`gene`）。`options`：`max_af`（0.05；高于它的变异被排除）、`exclude_benign`（True）、`genes`（面板列表，+0.5）、`top_n`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`parameters`、`ranked`（每项含 `chrom`、`pos`、`ref`、`alt`、`id`、`score`、`components`、`gene`、`rank`）、`excluded`（带 `excluded_reason`）、`n_ranked`、`n_excluded` 和 `n_frequency_unknown`（缺少人群频率注释的变异数，含已排除的；此时封装的 `warnings` 中有一条带该数量的消息）。每个 `ranked` 和 `excluded` 行都有布尔字段 `frequency_unknown`。

- **示例代码：**

```python
import json

from dnakit.clinical import prioritize

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\tGENE=BRCA1;consequence=stop_gained;clinvar=Pathogenic;gnomad_af=0.00001;cadd=35",
    "chr1\t200\t.\tC\tT\t50\tPASS\tGENE=TTN;consequence=missense_variant;gnomad_af=0.002;revel=0.6",
    "chr1\t300\t.\tG\tA\t50\tPASS\tGENE=OBSCN;consequence=synonymous_variant;gnomad_af=0.2",
    "chr1\t400\t.\tT\tC\t50\tPASS\tGENE=APOE;consequence=missense_variant;clinvar=Benign",
]) + "\n"

result = prioritize(vcf, options={"genes": ["BRCA1"], "top_n": 3})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "prioritize",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 406,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1462d42c295ee82b6444284c676a34284ef9d2be57066750c3aee1651621c7f2"
  },
  "warnings": ["1 of 4 variants have no population-frequency annotation and received the full rarity score..."],
  "status": "ok",
  "result": {
    "method": "additive heuristic score: consequence severity (Sequence Ontology table, 0-1) + ClinVar (p...",
    "parameters": {"max_af": 0.05, "exclude_benign": true, "genes": ["BRCA1"]},
    "ranked": [
      {
        "chrom": "chr1",
        "pos": 100,
        "ref": "A",
        "alt": ["G"],
        "id": null,
        "score": 3.4375,
        "components": {"consequence": 1.0, "clinvar": 1.0, "rarity": 0.5, "predictors": 0.4375, "panel_gene": 0.5},
        "excluded_reason": null,
        "gene": "BRCA1",
        "frequency_unknown": false,
        "rank": 1
      },
      {
        "chrom": "chr1",
        "pos": 200,
        "ref": "C",
        "alt": ["T"],
        "id": null,
        "score": 0.9,
        "components": {"consequence": 0.5, "rarity": 0.1, "predictors": 0.3},
        "excluded_reason": null,
        "gene": "TTN",
        "frequency_unknown": false,
        "rank": 2
      }
    ],
    "excluded": [
      {
        "chrom": "chr1",
        "pos": 300,
        "ref": "G",
        "alt": ["A"],
        "id": null,
        "score": 0.05,
        "components": {"consequence": 0.05, "rarity": 0.0},
        "excluded_reason": "population AF 0.2 > 0.05",
        "gene": "OBSCN",
        "frequency_unknown": false
      },
      {
        "chrom": "chr1",
        "pos": 400,
        "ref": "T",
        "alt": ["C"],
        "id": null,
        "score": 0.0,
        "components": {"consequence": 0.5, "clinvar": -1.0, "rarity": 0.5},
        "excluded_reason": "ClinVar benign / likely benign",
        "gene": "APOE",
        "frequency_unknown": true
      }
    ],
    "n_ranked": 2,
    "n_excluded": 2,
    "n_frequency_unknown": 1
  }
}
```

## 15) qc

- **作用：** 对临床调用集做质控：逐样本调用率、het/hom 比、Ti/Tv 与深度，逐位点调用率与精确 HWE 检验，并给出阈值标记。

- **处理方法：** 逐样本指标与 `dnakit.variants.summarize` 相同（调用率、hom-ref/het/hom-alt 计数、het/hom-alt 比、Ti/Tv、平均深度、平均 GQ）；逐位点指标为调用率、等位基因计数/频率以及 `calculate_hwe_pvalue` 的精确 HWE p 值。样本调用率低于 `min_sample_call_rate` 时标记 `low_call_rate`，仅在给出相应选项时才标记 `low_mean_depth`、`high_het_hom_ratio` 或 `ti_tv_out_of_range`；位点调用率低于 `min_site_call_rate` 标记 `low_call_rate`，HWE p 值低于 `min_hwe_pvalue` 标记 `hwe_violation`。没有样本被标记时 `passed` 为真（位点标记不影响它）。它只涉及已调用的基因型，不含原始测序运行指标（碱基质量、重复率、污染），且阈值是通用默认值，并非实验室验证过的标准。

- **API：** `dnakit.clinical.qc(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`min_sample_call_rate`（0.95）、`min_site_call_rate`（0.9）、`min_hwe_pvalue`（1e-6）、`min_mean_depth`、`max_het_hom_alt_ratio`、`ti_tv_range`（`[low, high]`）；后三项不给则不启用。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`thresholds`、`n_sites`、`n_samples`、`samples`（逐样本指标加 `flags`）、`samples_failing`、`n_sites_flagged`、`flagged_sites`（至多 1000 个，带逐位点指标与 `flags`）和 `passed`。

- **示例代码：**

```python
import json

from dnakit.clinical import qc

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/0:10\t0/1:12\t1/1:11",
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:DP\t0/1:9\t0/1:10\t./.:0",
    "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT:DP\t0/0:10\t./.:0\t./.:0",
    "chr1\t400\t.\tT\tC\t50\tPASS\t.\tGT:DP\t0/1:11\t1/1:12\t0/0:10",
]) + "\n"

result = qc(vcf, options={"min_sample_call_rate": 0.8})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "qc",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 280,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9a7e1e3b47679b6be206402306d35eae65b6c2089bdcde71517c09d310d2e168"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "call-set QC: per-sample call rate, het/hom-alt ratio, Ti/Tv and depth; per-site call rate...",
    "thresholds": {
      "min_sample_call_rate": 0.8,
      "min_site_call_rate": 0.9,
      "min_hwe_pvalue": 1e-06,
      "min_mean_depth": null,
      "max_het_hom_alt_ratio": null,
      "ti_tv_range": null
    },
    "n_sites": 4,
    "n_samples": 3,
    "samples": [
      {
        "sample": "S1",
        "n_sites": 4,
        "n_called": 4,
        "n_missing": 0,
        "call_rate": 1.0,
        "n_hom_ref": 2,
        "n_het": 2,
        "n_hom_alt": 0,
        "het_hom_alt_ratio": null,
        "n_snv": 2,
        "n_indel": 0,
        "n_transitions": 2,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 10.0,
        "mean_genotype_quality": null,
        "flags": []
      },
      {
        "sample": "S2",
        "n_sites": 4,
        "n_called": 3,
        "n_missing": 1,
        "call_rate": 0.75,
        "n_hom_ref": 0,
        "n_het": 2,
        "n_hom_alt": 1,
        "het_hom_alt_ratio": 2.0,
        "n_snv": 3,
        "n_indel": 0,
        "n_transitions": 3,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 8.5,
        "mean_genotype_quality": null,
        "flags": ["low_call_rate"]
      },
      {
        "sample": "S3",
        "n_sites": 4,
        "n_called": 2,
        "n_missing": 2,
        "call_rate": 0.5,
        "n_hom_ref": 1,
        "n_het": 0,
        "n_hom_alt": 1,
        "het_hom_alt_ratio": 0.0,
        "n_snv": 1,
        "n_indel": 0,
        "n_transitions": 1,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 5.25,
        "mean_genotype_quality": null,
        "flags": ["low_call_rate"]
      }
    ],
    "samples_failing": ["S2", "S3"],
    "n_sites_flagged": 2,
    "flagged_sites": [
      {
        "chrom": "chr1",
        "pos": 200,
        "ref": "C",
        "alt": ["T"],
        "n_samples": 3,
        "call_rate": 0.666667,
        "an": 4,
        "ac": [2],
        "af": [0.5],
        "n_het": 2,
        "n_hom_alt": 0,
        "hwe_pvalue": 1.0,
        "mean_depth": 6.33333,
        "qual": 50.0,
        "filter": "PASS",
        "flags": ["low_call_rate"]
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "ref": "G",
        "alt": ["A"],
        "n_samples": 3,
        "call_rate": 0.333333,
        "an": 2,
        "ac": [0],
        "af": [0.0],
        "n_het": 0,
        "n_hom_alt": 0,
        "hwe_pvalue": 1.0,
        "mean_depth": 3.33333,
        "qual": 50.0,
        "filter": "PASS",
        "flags": ["low_call_rate"]
      }
    ],
    "passed": false
  }
}
```

## 16) read

- **作用：** 把（ClinVar 风格的）VCF 或变异记录读成临床字段：基因、ClinVar 意义与审核状态、疾病、后果和基因型。

- **处理方法：** 这是解析器而不是解读器。对每条记录读取 INFO 的 `GENEINFO`/`GENE`（基因符号）、`CLNSIG`（原样保存在 `clinvar_significance_raw`，并规范化为小写下划线术语放入 `clinvar_terms`，如 `likely_pathogenic`）、`CLNREVSTAT`、`CLNDN`（按 `|` 拆分为疾病列表）、`MC`（Sequence Ontology 术语，或 `consequence` 键），以及每个样本的基因型。映射列表形式的输入可使用同名记录键（`gene`、`clinvar`、`consequence` 等）。不做外部查询，也不推导分类。

- **API：** `dnakit.clinical.read(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：无。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `samples`、`n_variants`、`method` 和 `variants`（每条记录：`chrom`、`pos`、`id`、`ref`、`alt`、`gene`、`clinvar_significance_raw`、`clinvar_terms`、`review_status`、`disease`、`consequence`、`filter`，VCF 含样本时还有 `genotypes`）。

- **示例代码：**

```python
import json

from dnakit.clinical import read

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr17\t43094464\trs80357906\tC\tCT\t.\t.\tGENEINFO=BRCA1:672;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Hereditary_breast_ovarian_cancer_syndrome;MC=SO:0001589|frameshift_variant",
    "chr2\t47806520\trs1\tA\tG\t.\t.\tGENEINFO=MSH6:2956;CLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=Lynch_syndrome|not_provided;MC=SO:0001583|missense_variant",
]) + "\n"

result = read(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "read",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 4,
    "text_length": 436,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b9a5b4d38401d514fb565af6b5e4ea890fab99311f3ea90124d1512a3de8c6d0"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "samples": [],
    "n_variants": 2,
    "variants": [
      {
        "chrom": "chr17",
        "pos": 43094464,
        "id": "rs80357906",
        "ref": "C",
        "alt": ["CT"],
        "gene": "BRCA1",
        "clinvar_significance_raw": "Pathogenic",
        "clinvar_terms": ["pathogenic"],
        "review_status": "reviewed_by_expert_panel",
        "disease": ["Hereditary_breast_ovarian_cancer_syndrome"],
        "consequence": ["frameshift_variant"],
        "filter": "."
      },
      {
        "chrom": "chr2",
        "pos": 47806520,
        "id": "rs1",
        "ref": "A",
        "alt": ["G"],
        "gene": "MSH6",
        "clinvar_significance_raw": "Uncertain_significance",
        "clinvar_terms": ["uncertain_significance"],
        "review_status": "criteria_provided,_single_submitter",
        "disease": ["Lynch_syndrome", "not_provided"],
        "consequence": ["missense_variant"],
        "filter": "."
      }
    ],
    "method": "VCF/records parsed into clinical fields: ClinVar CLNSIG (normalised to lower-case terms),..."
  }
}
```

## 17) report

- **作用：** 根据已分类的变异生成临床变异报告（Markdown 或纯文本）。

- **处理方法：** 每个变异被赋予一个层级：给出的 `classification`（`Pathogenic`、`Likely pathogenic`、`Uncertain significance`、`Likely benign`、`Benign`，或 `P`/`LP`/`VUS`/`LB`/`B`）原样使用；只给 `criteria` 时按 `classify` 的 ACMG/AMP 2015 规则组合（因此互相矛盾的标准得到意义不明）；否则为 `Not classified`。VCF 输入时，基因、ClinVar 术语和后果从记录中读取并作为证据显示，但 ClinVar 意义不会被转换成分类，所以除非另行提供分类，这些变异为 `Not classified`；`hgvs_accession`（及 `reference`）只用于添加 HGVS `g.` 字符串。变异按层级和基因排序，报告列出可报告结果（致病与可能致病）、意义不明或未分类的变异（除非 `include_vus` 为 False）、未列出的良性/可能良性变异数量，以及固定的方法与局限说明。它只汇总所提供的证据：不查询数据库，也不重新评估标准，因此文本不是经过验证的诊断报告，必须由合格的临床遗传学家审核。

- **API：** `dnakit.clinical.report(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'variants': [{'gene', 'hgvs', 'zygosity', 'classification' 或 'criteria', 'clinvar', ...}], 'sample': ..., 'indication': ...}`（`sample`、`patient`、`patient_id`、`indication`、`date` 会印在报告头部）、这类行的列表，或 VCF/记录（基因、ClinVar 与后果字段从记录中读取）。`options`：`title`、`format`（默认 `'markdown'`，另有 `'text'`）、`include_vus`（True）、用于为 VCF 输入生成 HGVS `g.` 字符串的 `hgvs_accession`（字符串或 `{contig: accession}`）与 `reference`、`output`（写出报告的路径）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `format`、`report`（完整报告文本）、`n_variants`、`classification_counts`、`n_reportable`、`variants`（带已确定 `classification` 的各行，由标准推出时还有 `classification_basis`）和 `output_path`。分类标签未知或变异列表为空时抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.clinical import report

rows = [
    {"gene": "CFTR", "hgvs": "NM_000492.4:c.1521_1523del", "zygosity": "homozygous", "classification": "P"},
    {"gene": "SCN5A", "hgvs": "NM_198056.3:c.5A>G", "zygosity": "heterozygous", "criteria": ["PM2", "PP3", "PP2"]},
    {"gene": "TTN", "hgvs": "NM_001267550.2:c.100C>T", "zygosity": "heterozygous", "classification": "B"},
]

result = report(rows, options={"format": "text", "include_vus": False, "title": "Cardiac panel"})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "clinical",
  "operation": "report",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 320,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "083c03cb00fc96f648d730193daf78d3db4b45b82aff7d4b708e50d754ecbe62"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "text",
    "report": "Cardiac panel\n\nVariants evaluated: 3 (1 Benign, 1 Pathogenic, 1 Uncertain significance)\n\nR...",
    "n_variants": 3,
    "classification_counts": {"Benign": 1, "Pathogenic": 1, "Uncertain significance": 1},
    "n_reportable": 1,
    "variants": [
      {"gene": "CFTR", "hgvs": "NM_000492.4:c.1521_1523del", "zygosity": "homozygous", "classification": "Pathogenic"},
      {
        "gene": "SCN5A",
        "hgvs": "NM_198056.3:c.5A>G",
        "zygosity": "heterozygous",
        "criteria": ["PM2", "PP3", "PP2"],
        "classification_basis": "ACMG/AMP 2015 combination of supplied criteria",
        "classification": "Uncertain significance"
      },
      {"gene": "TTN", "hgvs": "NM_001267550.2:c.100C>T", "zygosity": "heterozygous", "classification": "Benign"}
    ],
    "output_path": null
  }
}
```

完整的 `result["report"]` 字符串如下：

```text
Cardiac panel

Variants evaluated: 3 (1 Benign, 1 Pathogenic, 1 Uncertain significance)

Reportable findings (Pathogenic / Likely pathogenic)

Gene=CFTR; Variant=NM_000492.4:c.1521_1523del; Zygosity=homozygous; Classification=Pathogenic; Evidence=-

Not reported

1 benign / likely benign variants and 1 uncertain variants are not listed.

Method and limitations

Classifications are taken from the input or computed from the supplied ACMG/AMP criteria with the 2015 combining rules (Richards et al., Genet Med 2015). The criteria were not re-evaluated against databases. This automatically generated summary is not a validated diagnostic report and must be reviewed by a qualified clinical geneticist.
```

## 18) type_hla

- **作用：** 从测序 reads 中进行 HLA 分型。尚未实现：函数只报告缺少什么（HLA 分型程序（OptiType、HLA-LA、arcasHLA、T1K）、IPD-IMGT/HLA 参考库与 reads）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。函数不会运行任何外部程序。

- **API：** `dnakit.clinical.type_hla(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"HLA typing program (OptiType/HLA-LA/arcasHLA) with IPD-IMGT/HLA reference and reads"）与 `message`（"HLA typing needs an external program and the IMGT/HLA reference."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.clinical import type_hla

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = type_hla(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "clinical",
  "operation": "type_hla",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 85,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d5f263542fc785f67a77f789993c401e91b0654ff0a57995475ebe6e41f3f69"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "HLA typing program (OptiType/HLA-LA/arcasHLA) with IPD-IMGT/HLA reference and reads",
  "message": "HLA typing needs an external program and the IMGT/HLA reference."
}
```
