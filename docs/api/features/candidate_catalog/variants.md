# 变异过滤、规范化、注释与统计

本类收录 **63** 项有文档的函数。其中 51 项按各节所述方法，根据传入的数据（VCF 文本、VCF 路径或变异记录、基因型计数、read 计数、DNA 序列）计算结果；12 项始终返回 `status: "not_implemented"`，因为它们需要 DNAKit 未内置的注释数据库、训练好的预测器、外部检出程序或比对后的 reads（`annotate`、`annotate_conservation_score`、`annotate_mt`、`annotate_pathogenicity_score`、`annotate_population_frequency`、`annotate_regulatory_effect`、`annotate_splice_effect`、`annotate_variant_consequence`、`call`、`call_mt`、`detect_breakpoints`、`genotype_structural_variants`）。可计算的函数中有五个需要某个选项，缺少时抛出 `ConfigurationError`：`left_align`、`left_align_mt`、`left_align_variants` 需要 `options.reference`，`detect_repeat_expansions` 需要 `options.motif`，`phase_genotypes` 需要 `options.trio`。有些函数是有意的近似，并已在各自的“处理方法”中写明：`fill_missing_genotypes_with_mode` 是朴素的逐位点众数基因型填补，不是 Beagle 式的 imputation；`calculate_kinship_king` 是 KING-robust 亲缘系数，不是 IBD 片段检测；启发式方法 `detect_low_frequency_variants`、`detect_mosaic_variants`、`calculate_runs_of_homozygosity`、`prioritize_candidate_variants` 不是经过验证的检出工具或分类器。`_mt` 系列函数与对应的普通函数一样接受 VCF 文本、路径或记录（`count_mt` 除外，它只统计非空行或列表元素的个数），不读取 Hail MatrixTable 对象。输入格式错误或缺少必需选项时，函数抛出 `InputFormatError` 或 `ConfigurationError`（带 `code` 与 `hint`），不返回封装。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) annotate

- **作用：** 为变异添加功能、人群和临床注释。尚未实现：函数只报告缺少什么（转录本模型、人群频率与临床数据库）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。请改用 `annotate_coding_effect`（基于用户提供的 CDS）或外部程序封装 `annotate_vcf_ensembl_vep`、`annotate_vcf_vcfanno`。

- **API：** `dnakit.variants.annotate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"annotation databases (Ensembl/VEP cache, gnomAD, ClinVar) or the vep / vcfanno programs"）与 `message`（"annotation needs transcript, population-frequency and clinical reference data that DNAKit does not bundle."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import annotate

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
  "module": "variants",
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
  "required_dependency": "annotation databases (Ensembl/VEP cache, gnomAD, ClinVar) or the vep / vcfanno programs",
  "message": "annotation needs transcript, population-frequency and clinical reference data that DNAKit does not bundle."
}
```

## 2) annotate_coding_effect

- **作用：** 在用户提供的编码序列上预测变异的蛋白水平效应（同义、错义、无义、起始/终止丢失、移码、框内插入缺失）。

- **处理方法：** 用标准遗传密码翻译参考 CDS 与突变后的 CDS 并比较（不使用转录本模型，也不使用基因组）。位置为 1 起始的 CDS 坐标，REF/ALT 为编码链上的等位基因；先裁掉 REF/ALT 共有的两侧碱基，REF 必须与 CDS 一致。等长变化逐密码子判定为 `start_lost`、`stop_gained`、`stop_lost`、`missense_variant` 或 `synonymous_variant`（跨多个密码子时取最严重者）；长度变化不是 3 的倍数时为 `frameshift_variant`（适用时为 `stop_gained`/`start_lost`）；其余长度变化为 `inframe_insertion`、`inframe_deletion` 或 `protein_altering_variant`。蛋白 HGVS 字符串（`p.Ala2Val`、`p.Lys3Ter`、`p.Ala2=`、`p.Met1?`、`fs`、`ins`、`del`、`delins`）由第一个改变的残基生成，未做 3' 规范化。它不能替代 VEP 等基于转录本的注释工具：基因组坐标变异、UTR、内含子和剪接均不在范围内。

- **API：** `dnakit.variants.annotate_coding_effect(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'cds': 'ATG...', 'variants': [...]}`，每个变异是 `{'cds_position': 5, 'ref': 'C', 'alt': 'T'}` 或 `'c.5C>T'` 字符串（也接受带 `cds_position` 的单个字典；`'-'` 或 `''` 表示空等位基因）。`options.cds`：当 CDS 不在 `input` 中时在此提供（仅 ACGT，最长 1,000,000 个碱基）。最多 10,000 个变异。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`cds_length`、`protein_length` 和 `effects`；`effects` 每个变异一项，含 `cds_position`、`ref`、`alt`（裁剪后）、`consequence`、`protein_position`、`hgvs_p`，碱基替换另有 `codons`（`codon_index`、`ref_codon`、`alt_codon`、`ref_aa`、`alt_aa`）。CDS 长度不是 3 的倍数或不以 ATG 开头时 `warnings` 会提示。REF 与 CDS 不符或位置超出 CDS 时抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.variants import annotate_coding_effect

cds = "ATGGCCAAATGGTAG"  # M A K W stop

result = annotate_coding_effect({
    "cds": cds,
    "variants": [
        "c.5C>T",                                   # GCC -> GTC  Ala2Val
        {"cds_position": 7, "ref": "A", "alt": "T"},  # AAA -> TAA  Lys3Ter
        {"cds_position": 6, "ref": "C", "alt": "A"},  # GCC -> GCA  Ala2Ala
    ],
})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "annotate_coding_effect",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 140,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8a61ab07dc9daab0ac0e6b74d5f6a734627313a7cf8477a2dc29a261291ad048"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "translation of the reference and mutated CDS with the standard genetic code; REF/ALT are o...",
    "cds_length": 15,
    "protein_length": 4,
    "effects": [
      {
        "cds_position": 5,
        "ref": "C",
        "alt": "T",
        "consequence": "missense_variant",
        "codons": [{"codon_index": 2, "ref_codon": "GCC", "alt_codon": "GTC", "ref_aa": "A", "alt_aa": "V"}],
        "protein_position": 2,
        "hgvs_p": "p.Ala2Val"
      },
      {
        "cds_position": 7,
        "ref": "A",
        "alt": "T",
        "consequence": "stop_gained",
        "codons": [{"codon_index": 3, "ref_codon": "AAA", "alt_codon": "TAA", "ref_aa": "K", "alt_aa": "*"}],
        "protein_position": 3,
        "hgvs_p": "p.Lys3Ter"
      },
      {
        "cds_position": 6,
        "ref": "C",
        "alt": "A",
        "consequence": "synonymous_variant",
        "codons": [{"codon_index": 2, "ref_codon": "GCC", "alt_codon": "GCA", "ref_aa": "A", "alt_aa": "A"}],
        "protein_position": 2,
        "hgvs_p": "p.Ala2="
      }
    ]
  }
}
```

## 3) annotate_conservation_score

- **作用：** 为变异添加 phyloP、phastCons 或 GERP 保守性得分。尚未实现：函数只报告缺少什么（全基因组逐碱基保守性轨道（bigWig/bedGraph））。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。函数不会读取任何保守性轨道。

- **API：** `dnakit.variants.annotate_conservation_score(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"per-base conservation track (phyloP, phastCons or GERP)"）与 `message`（"conservation scoring needs a genome-wide conservation track that is not available."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import annotate_conservation_score

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_conservation_score(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "annotate_conservation_score",
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
  "required_dependency": "per-base conservation track (phyloP, phastCons or GERP)",
  "message": "conservation scoring needs a genome-wide conservation track that is not available."
}
```

## 4) annotate_mt

- **作用：** 注释线粒体变异。尚未实现：函数只报告缺少什么（MITOMAP、gnomAD-mt、HmtVar 等线粒体注释资源）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。函数不会查询任何线粒体数据库。

- **API：** `dnakit.variants.annotate_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"mitochondrial annotation resources (MITOMAP, gnomAD-mt, HmtVar)"）与 `message`（"mitochondrial annotation needs reference annotation databases that are not available."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import annotate_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_mt(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "annotate_mt",
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
  "required_dependency": "mitochondrial annotation resources (MITOMAP, gnomAD-mt, HmtVar)",
  "message": "mitochondrial annotation needs reference annotation databases that are not available."
}
```

## 5) annotate_pathogenicity_score

- **作用：** 为变异添加致病性预测得分（CADD、REVEL、AlphaMissense 等）。尚未实现：函数只报告缺少什么（训练好的预测器或预计算得分表）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。函数不会估算任何得分；已经写在记录上的得分可由 `prioritize_candidate_variants` 综合使用。

- **API：** `dnakit.variants.annotate_pathogenicity_score(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"trained pathogenicity predictor or precomputed score tables (CADD/REVEL/AlphaMissense)"）与 `message`（"pathogenicity scores come from trained models or precomputed score files that are not available."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import annotate_pathogenicity_score

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_pathogenicity_score(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "annotate_pathogenicity_score",
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
  "required_dependency": "trained pathogenicity predictor or precomputed score tables (CADD/REVEL/AlphaMissense)",
  "message": "pathogenicity scores come from trained models or precomputed score files that are not available."
}
```

## 6) annotate_population_frequency

- **作用：** 为变异添加人群等位基因频率（gnomAD、1000 Genomes、dbSNP）。尚未实现：函数只报告缺少什么（频率数据库或网络访问，二者此处均不使用）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。函数不会发起查询，也不读取频率文件。

- **API：** `dnakit.variants.annotate_population_frequency(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"population frequency database (gnomAD, 1000 Genomes, dbSNP)"）与 `message`（"population frequencies need a frequency database or network access, neither of which is used here."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import annotate_population_frequency

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_population_frequency(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "annotate_population_frequency",
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
  "required_dependency": "population frequency database (gnomAD, 1000 Genomes, dbSNP)",
  "message": "population frequencies need a frequency database or network access, neither of which is used here."
}
```

## 7) annotate_regulatory_effect

- **作用：** 为变异添加调控元件效应注释（Ensembl Regulatory Build、ENCODE cCRE、序列到活性模型）。尚未实现：函数只报告缺少什么（调控区域轨道或训练好的模型）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。函数不做任何调控注释。

- **API：** `dnakit.variants.annotate_regulatory_effect(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"regulatory annotation track or trained regulatory model"）与 `message`（"regulatory annotation needs regulatory-region tracks or trained models that are not available."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import annotate_regulatory_effect

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_regulatory_effect(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "annotate_regulatory_effect",
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
  "required_dependency": "regulatory annotation track or trained regulatory model",
  "message": "regulatory annotation needs regulatory-region tracks or trained models that are not available."
}
```

## 8) annotate_splice_effect

- **作用：** 为变异添加剪接效应预测（SpliceAI、MaxEntScan）。尚未实现：函数只报告缺少什么（转录本模型与训练好的剪接预测器）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。函数不做任何剪接预测。

- **API：** `dnakit.variants.annotate_splice_effect(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"transcript annotation and a splice-effect predictor (SpliceAI)"）与 `message`（"splice-effect prediction needs transcript models and a trained predictor that are not available."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import annotate_splice_effect

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_splice_effect(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "annotate_splice_effect",
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
  "required_dependency": "transcript annotation and a splice-effect predictor (SpliceAI)",
  "message": "splice-effect prediction needs transcript models and a trained predictor that are not available."
}
```

## 9) annotate_variant_consequence

- **作用：** 按转录本注释变异的 Sequence Ontology 后果。尚未实现：函数只报告缺少什么（转录本模型（GTF/GFF）与参考基因组）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。只有 `annotate_coding_effect` 可用，且仅基于用户提供的 CDS；基因组坐标变异请用 `annotate_vcf_ensembl_vep`。

- **API：** `dnakit.variants.annotate_variant_consequence(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"transcript annotation (GTF/GFF) and genome FASTA, e.g. via Ensembl VEP"）与 `message`（"consequence calling needs transcript models and a reference genome; only annotate_coding_effect works from a user-supplied CDS."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import annotate_variant_consequence

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = annotate_variant_consequence(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "annotate_variant_consequence",
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
  "required_dependency": "transcript annotation (GTF/GFF) and genome FASTA, e.g. via Ensembl VEP",
  "message": "consequence calling needs transcript models and a reference genome; only annotate_coding_e..."
}
```

## 10) calculate_allele_frequency

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `calculate_allele_frequency` 功能。候选编号 CAP-0768，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_allele_frequency` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.calculate_allele_frequency(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref` 和 `alt`，基于基因型的统计还需显式提供 `genotypes`。

- **输出：** 本地函数 `calculate_allele_frequency` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import calculate_allele_frequency

result = calculate_allele_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0768",
  "operation": "calculate_allele_frequency",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_allele_frequency', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'alt_count': 3, 'called_alleles': 6, 'frequency': 0.5}], 'total_alt_alleles': 3, 'total_called_alleles': 6, 'aggregate_frequency': 0.5}, 'provenance': ...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 11) calculate_copy_number_profile

- **作用：** 把逐 bin 测序深度转换为拷贝数谱：相对中位数的深度比、拷贝数以及合并后的区段。

- **处理方法：** 比值 = 深度 / 深度中位数（给定配对对照时再除以对照的同一比值）；拷贝数 = (比值 * (纯度 * 倍性 + 2 * (1 - 纯度)) - 2 * (1 - 纯度)) / 纯度，下限为 0。中位数取自 `median_exclude_chroms` 之外且深度为正的 bin。同一染色体上取整后拷贝数相同的相邻 bin 合并为一个区段，标注为 `homozygous_deletion`、`loss`、`neutral`、`gain`（至倍性 + 2）或 `amplification`。这只是深度比的换算：不做统计分段（CBS、HMM），不做 GC 或可比对性校正，也不使用 B 等位基因信息。

- **API：** `dnakit.variants.calculate_copy_number_profile(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：bin 列表 `{'chrom', 'start', 'end', 'depth'[, 'control_depth']}`（或 `[chrom, start, end, depth]`），或 `{'bins': [...]}`；`control_depth` 要么每个 bin 都给，要么都不给。`options`：`purity`（默认 1.0，范围 0.01-1）、`ploidy`（默认 2.0）、`median_exclude_chroms`（默认 X、Y、M、MT，带不带 `chr` 均可）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`purity`、`ploidy`、`median_depth`、`control_median_depth`、`bins`（每个 bin：`chrom`、`start`、`end`、`depth`、`ratio`、`log2_ratio`、`copy_number`、`rounded_copy_number`）和 `segments`（`chrom`、`start`、`end`、`copy_number`、`n_bins`、`mean_ratio`、`mean_log2_ratio`、`state`）。

- **示例代码：**

```python
import json

from dnakit.variants import calculate_copy_number_profile

bins = [
    ("chr1", 1, 1000, 30), ("chr1", 1001, 2000, 31), ("chr1", 2001, 3000, 60),
    ("chr1", 3001, 4000, 62), ("chr2", 1, 1000, 29), ("chr2", 1001, 2000, 15),
]

result = calculate_copy_number_profile(bins)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "calculate_copy_number_profile",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 150,
    "sequence_bases": 18,
    "lengths": [3, 3, 3, 3, 3, 3],
    "sha256": "711b9ccb7ed0139eca19642183f362a0f51ae054fb02c72c67491cc8e2648b6f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "per-bin depth ratio to the median (or to the control depth ratio), copy number = (ratio*(p...",
    "purity": 1.0,
    "ploidy": 2.0,
    "median_depth": 30.5,
    "control_median_depth": null,
    "bins": [
      {
        "chrom": "chr1",
        "start": 1,
        "end": 1000,
        "depth": 30.0,
        "ratio": 0.983607,
        "log2_ratio": -0.0238467,
        "copy_number": 1.96721,
        "rounded_copy_number": 2
      },
      {
        "chrom": "chr1",
        "start": 1001,
        "end": 2000,
        "depth": 31.0,
        "ratio": 1.01639,
        "log2_ratio": 0.023459,
        "copy_number": 2.03279,
        "rounded_copy_number": 2
      },
      "..."
    ],
    "segments": [
      {
        "chrom": "chr1",
        "start": 1,
        "end": 2000,
        "copy_number": 2,
        "n_bins": 2,
        "mean_ratio": 1.0,
        "mean_log2_ratio": 0.0,
        "state": "neutral"
      },
      {
        "chrom": "chr1",
        "start": 2001,
        "end": 4000,
        "copy_number": 4,
        "n_bins": 2,
        "mean_ratio": 2.0,
        "mean_log2_ratio": 1.0,
        "state": "gain"
      },
      "..."
    ]
  }
}
```

## 12) calculate_genotype_likelihood

- **作用：** 在二项错误模型下，根据 reads 证据计算双等位位点的基因型似然（log10）、Phred 标度 PL 与 GQ。

- **处理方法：** 对 m 倍体中含 k 个 ALT 拷贝的基因型，每条 read 独立地以概率 (k/m)(1-e) + (1-k/m)e 支持 ALT；e 为每条 read 的错误概率（由碱基质量按 10^(-q/10) 得到，上限 0.5；未给质量时用 `error_rate`）。log10 L 为各 read 之和；PL = -10 * (log10 L - 最大 log10 L)，取整；GQ 为第二小的 PL，上限 99。各 read 视为相互独立，没有先验、比对质量或链偏倚模型，因此这是教科书模型，而不是 GATK 或 bcftools 的实现。

- **API：** `dnakit.variants.calculate_genotype_likelihood(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'ref_count': int, 'alt_count': int}`、`{'reads': [{'allele': 'ref'|'alt', 'base_quality': q}, ...]}`、这类字典的列表，或带逐样本 `FORMAT/AD` 的 VCF/记录（多等位位点先拆分）。`options`：`error_rate`（默认 0.01，未给碱基质量时使用）、`ploidy`（默认 2，范围 1-6）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `model`、`error_rate`、`ploidy`，以及 `result`（单个计数或 reads 输入）或 `results`（列表输入或 VCF 输入，此时每项另有 `site` 与 `sample`）。每个似然条目含 `genotypes`、`log10_likelihoods`、`phred_scaled`、`best_genotype`、`genotype_quality`、`ref_reads`、`alt_reads`。没有可用 AD 的调用会被跳过并给出警告。

- **示例代码：**

```python
import json

from dnakit.variants import calculate_genotype_likelihood

result = calculate_genotype_likelihood({"ref_count": 10, "alt_count": 0})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "calculate_genotype_likelihood",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 33,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b06a2868736086281799ef5b0bab69005c6fccac827322fd56bb70f49b9e82b6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": "each read independently supports ALT with probability (k/m)(1-e)+(1-k/m)e for a genotype w...",
    "error_rate": 0.01,
    "ploidy": 2,
    "result": {
      "genotypes": ["0/0", "0/1", "1/1"],
      "log10_likelihoods": [-0.0436481, -3.0103, -20.0],
      "phred_scaled": [0, 30, 200],
      "best_genotype": "0/0",
      "genotype_quality": 30,
      "ref_reads": 10,
      "alt_reads": 0
    }
  }
}
```

## 13) calculate_haplotype_frequency

- **作用：** 估计两个或更多位点的单倍型频率：对已定相单倍型直接计数，或对未定相基因型用 EM 算法估计。

- **处理方法：** 已定相的单倍型直接计数。未定相的二倍体基因型按 Excoffier 与 Slatkin（1995）的 EM 算法解析：枚举每个个体基因型可能的所有单倍型配对，初始频率在候选单倍型上均匀分布，E/M 步骤假定 Hardy-Weinberg 比例，直到最大变化小于 `tolerance`。含缺失基因型的个体被排除；EM 支持 2-12 个位点，且每个个体至多 10 个杂合位点。等位基因用整数编码（0 为 REF，1 为第一个 ALT）。

- **API：** `dnakit.variants.calculate_haplotype_frequency(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'haplotypes': ['010', '011', ...]}`（已定相）、`{'genotypes': [['0/1', '1/1'], ...]}` 或逐个体的 GT 列表（每个位点一个 GT），或带样本基因型的 VCF/记录（位点为双等位记录）。`options`：`use_phase`（默认 True；VCF 输入且所有基因型均已定相时直接计数）、`max_iterations`（默认 1000）、`tolerance`（默认 1e-10）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`n_chromosomes`、`loci` 和 `haplotypes`（每个单倍型：形如 `0-1` 的 `haplotype`、`alleles`、`count`、`frequency`，按频率降序）；EM 分支另有 `n_individuals` 与 `iterations`。`warnings` 报告被排除的个体或跳过的多等位位点。EM 的计数与频率是浮点估计值。

- **示例代码：**

```python
import json

from dnakit.variants import calculate_haplotype_frequency

# unphased genotypes of 4 individuals at 2 loci (EM estimation)
result = calculate_haplotype_frequency({"genotypes": [["0/0", "0/0"], ["0/1", "0/1"], ["1/1", "1/1"], ["0/0", "0/1"]]})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "calculate_haplotype_frequency",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 79,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "17383568bffa3f8dde9185cc6b71f19f7de028428802b7bc3f05c24d779147b1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "EM algorithm for haplotype frequencies from unphased genotypes (Excoffier & Slatkin 1995;...",
    "n_individuals": 4,
    "n_chromosomes": 8,
    "iterations": 11,
    "loci": ["locus_1", "locus_2"],
    "haplotypes": [
      {"haplotype": "0-0", "alleles": [0, 0], "count": 4.0, "frequency": 0.5},
      {"haplotype": "1-1", "alleles": [1, 1], "count": 3.0, "frequency": 0.375},
      "..."
    ]
  }
}
```

## 14) calculate_heterozygosity

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `calculate_heterozygosity` 功能。候选编号 CAP-0772，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_heterozygosity` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.calculate_heterozygosity(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref` 和 `alt`，基于基因型的统计还需显式提供 `genotypes`。

- **输出：** 本地函数 `calculate_heterozygosity` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import calculate_heterozygosity

result = calculate_heterozygosity(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0772",
  "operation": "calculate_heterozygosity",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_heterozygosity', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'heterozygous_samples': 1, 'called_samples': 3, 'observed_heterozygosity': 0.3333333333333333}], 'heterozygous_samples': 1, 'called_samples': 3, 'observe...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 15) calculate_hwe_pvalue

- **作用：** 用 Wigginton、Cutler 与 Abecasis（2005）的精确检验检验 Hardy-Weinberg 平衡。

- **处理方法：** 给定基因型计数，用精确递推计算在已观察到的等位基因数下所有杂合子数的概率；双侧 p 值是所有不比观察到的表更可能的表的概率之和。使用 `mid_p` 时再减去观察表概率的一半。VCF 输入时对每个双等位二倍体位点各算一个 p 值（多等位及非二倍体调用会被跳过并给出警告）。

- **API：** `dnakit.variants.calculate_hwe_pvalue(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：基因型计数 `[hom_ref, het, hom_alt]`（或 `{'counts': [...]}`、`{'hom_ref', 'het', 'hom_alt'}`），或带样本基因型的 VCF/记录。`options`：`mid_p`（默认 False）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。 输入为计数时，`result` 含 `method`、`mid_p`、`n_hom_ref`、`n_het`、`n_hom_alt`、`n_samples`、`observed_het_fraction`、`expected_het_fraction`、`ref_allele_frequency` 和 `p_value`。VCF 输入时，`result.sites` 为每个位点给出同样字段并附 `site`。

- **示例代码：**

```python
import json

from dnakit.variants import calculate_hwe_pvalue

result = calculate_hwe_pvalue([10, 20, 30])  # hom_ref, het, hom_alt counts

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "calculate_hwe_pvalue",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 12,
    "sequence_bases": 60,
    "lengths": [10, 20, 30],
    "sha256": "15d98dad2ba94449f66d415d63b1e633ed896a72c15c6b91407410b8e06f5e20"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "exact test of Hardy-Weinberg equilibrium (Wigginton, Cutler & Abecasis 2005), two-sided",
    "mid_p": false,
    "n_hom_ref": 10,
    "n_het": 20,
    "n_hom_alt": 30,
    "n_samples": 60,
    "observed_het_fraction": 0.333333,
    "expected_het_fraction": 0.444444,
    "ref_allele_frequency": 0.333333,
    "p_value": 0.0776747
  }
}
```

## 16) calculate_kinship_king

- **作用：** 用 KING-robust 估计器估计样本两两之间的亲缘系数。这是基于基因型的亲缘系数，而不是 IBD 片段检测。

- **处理方法：** 对每一对样本，在两者都有调用的双等位位点上：kinship = (N_AaAa - 2 N_AAaa) / (N_Aa(i) + N_Aa(j))（Manichaikul 等 2010），其中 N_AaAa 为两者均杂合的位点数，N_AAaa 为相反纯合子位点数，N_Aa(i) 为样本 i 的杂合子数。`pi_hat_approx` 就是 2 * kinship，`ibs0_fraction` 为相反纯合子比例，`relationship` 依据 0.3535（重复样本或同卵双胞胎）、0.1768（一级）、0.0884（二级）、0.0442（三级）的阈值给出标签，否则为 `unrelated`（恰好等于阈值的值归入较低一级；阈值即 2^-1.5、2^-2.5、2^-3.5、2^-4.5 取四位小数）。函数不推断 IBD 片段，不使用连锁或重组模型，只使用基因型计数（不使用等位基因频率）；估计需要大量独立的常见 SNP，位点很少时（如示例）数值噪声大，甚至可能为负。最多 150 个样本。

- **API：** `dnakit.variants.calculate_kinship_king(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：多样本 VCF（文本或路径）或带二倍体基因型的记录；仅使用双等位位点。`options`：无。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method` 与 `pairs`，每对样本一项，含 `sample_1`、`sample_2`、`n_sites`、`n_het_het`、`n_opposite_homozygous`、`ibs0_fraction`、`kinship`、`pi_hat_approx` 和 `relationship`（两个样本都没有杂合位点时 `kinship` 为 `null`）。样本少于 2 个、多于 150 个或没有双等位位点时抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.variants import calculate_kinship_king

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t1/1\t0/1",
    "chr1\t200\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/0\t0/1",
    "chr1\t300\t.\tA\tG\t.\t.\t.\tGT\t1/1\t1/1\t0/1\t0/0",
    "chr1\t400\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/1\t1/1",
    "chr1\t500\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t0/1\t0/1",
    "chr1\t600\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t1/1\t0/0",
    "chr1\t700\t.\tA\tG\t.\t.\t.\tGT\t1/1\t1/1\t0/0\t0/1",
    "chr1\t800\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/1\t0/1",
]) + "\n"

result = calculate_kinship_king(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "calculate_kinship_king",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 10,
    "text_length": 399,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7258ad9577e33ae92458e95ff141dfec5403d4cf8329ffb6e15f1b20c0d8b391"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "KING-robust kinship estimator (Manichaikul et al. 2010): (N_AaAa - 2 N_AAaa) / (N_Aa(i) +...",
    "pairs": [
      {
        "sample_1": "S1",
        "sample_2": "S2",
        "n_sites": 8,
        "n_het_het": 4,
        "n_opposite_homozygous": 0,
        "ibs0_fraction": 0.0,
        "kinship": 0.5,
        "pi_hat_approx": 1.0,
        "relationship": "duplicate_or_monozygotic_twin"
      },
      {
        "sample_1": "S1",
        "sample_2": "S3",
        "n_sites": 8,
        "n_het_het": 2,
        "n_opposite_homozygous": 2,
        "ibs0_fraction": 0.25,
        "kinship": -0.25,
        "pi_hat_approx": -0.5,
        "relationship": "unrelated"
      },
      "..."
    ]
  }
}
```

## 17) calculate_ld_matrix

- **作用：** 计算双等位变异之间的两两连锁不平衡矩阵（r2、r 或 D'）。

- **处理方法：** 已定相时，用每对变异的单倍型频率计算 D = p_AB - p_A p_B、r = D / sqrt(p_A(1-p_A) p_B(1-p_B))、r2 = r^2 与 D' = |D| / D_max。未定相时使用 ALT 剂量（0/1/2）的 Pearson 相关（Rogers-Huff r）；此时 D' 不可用，会报错。每一对变异中缺失任一基因型的样本被丢弃；可用样本少于 2 个或变异单态时结果为 `null`。`phase: 'auto'` 仅在所有已调用杂合基因型都已定相时使用单倍型。多等位位点被跳过并给出警告。这是未经小样本校正的样本统计量。

- **API：** `dnakit.variants.calculate_ld_matrix(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：带样本基因型的 VCF 文本/路径或记录（至少 2 个样本和 2 个双等位变异）。`options`：`statistic`（默认 `'r2'`，另有 `'r'`、`'d_prime'`）、`phase`（默认 `'auto'`，另有 `'phased'`、`'unphased'`）、`max_variants`（默认 200，至多 1000）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `statistic`、`method`（`haplotype-based (phased)` 或剂量相关的说明）、`variants`（形如 `chr1:100:A:G` 的位点名）和 `matrix`（对称方阵，无定义处为 `null`）。

- **示例代码：**

```python
import json

from dnakit.variants import calculate_ld_matrix

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t0/1",
    "chr1\t110\t.\tC\tT\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t0/1",
    "chr1\t120\t.\tG\tA\t.\t.\t.\tGT\t1/1\t0/1\t0/0\t0/1",
    "chr1\t130\t.\tT\tC\t.\t.\t.\tGT\t0/0\t0/0\t0/1\t0/1",
]) + "\n"

result = calculate_ld_matrix(vcf)  # unphased genotypes: r^2 of ALT dosages

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "calculate_ld_matrix",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 239,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ce2023f5d41f80b38428eb4e237b65db249da7073fb2981324bde0150924df8c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "statistic": "r2",
    "method": "squared/plain Pearson correlation of ALT-allele dosages (unphased, Rogers-Huff style r)",
    "variants": ["chr1:100:A:G", "chr1:110:C:T", "..."],
    "matrix": [[1.0, 1.0, "..."], [1.0, 1.0, "..."], "..."]
  }
}
```

## 18) calculate_minor_allele_frequency

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `calculate_minor_allele_frequency` 功能。候选编号 CAP-0776，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_minor_allele_frequency` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.calculate_minor_allele_frequency(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref` 和 `alt`，基于基因型的统计还需显式提供 `genotypes`。必须显式提供 `genotypes`；不进行缺失基因型填补。

- **输出：** 本地函数 `calculate_minor_allele_frequency` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import calculate_minor_allele_frequency

result = calculate_minor_allele_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0776",
  "operation": "calculate_minor_allele_frequency",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_minor_allele_frequency', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'called_alleles': 6, 'allele_frequencies': {'0': 0.5, '1': 0.5}, 'minor_allele_index': 0, 'minor_allele': 'A', 'maf': 0.5}], 'variant_count': 1, ...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 19) calculate_missingness

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `calculate_missingness` 功能。候选编号 CAP-0777，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_missingness` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.calculate_missingness(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref` 和 `alt`，基于基因型的统计还需显式提供 `genotypes`。

- **输出：** 本地函数 `calculate_missingness` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import calculate_missingness

result = calculate_missingness(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0777",
  "operation": "calculate_missingness",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_missingness', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'missing_samples': 0, 'sample_count': 3, 'missingness': 0.0}], 'missing_samples': 0, 'sample_count': 3, 'missingness': 0.0}, 'provenance': {'implementation'...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 20) calculate_runs_of_homozygosity

- **作用：** 根据连续的纯合基因型在单个样本中查找纯合子区段（ROH）。

- **处理方法：** 对每条染色体从左到右贪心扫描：从一个纯合调用开始，只要杂合/缺失调用数不超过 `max_heterozygotes`、`max_missing` 且相邻变异间距不超过 `max_gap_bp` 就继续延伸；末尾的杂合/缺失调用会被裁掉。少于 `min_snps` 个调用或短于 `min_length_bp` 的区段被丢弃。这是类似 PLINK `--homozyg` 的连续 SNP 启发式方法，不是 HMM，也不使用遗传图谱、等位基因频率或基因型错误模型。默认值（25 个 SNP、1 Mb）面向稠密的全基因组数据。

- **API：** `dnakit.variants.calculate_runs_of_homozygosity(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：带基因型的 VCF 文本/路径或记录。`options`：`sample`（默认第一个样本）、`min_snps`（25）、`min_length_bp`（1000000）、`max_heterozygotes`（1）、`max_missing`（2）、`max_gap_bp`（1000000）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `sample`、`method`、`parameters`、`runs`（每项含 `chrom`、`start`、`end`、`length_bp`、`n_snps`、`n_heterozygous`、`n_missing`）、`n_runs` 和 `total_roh_bp`。

- **示例代码：**

```python
import json

from dnakit.variants import calculate_runs_of_homozygosity

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t1000\t.\tA\tG\t.\t.\t.\tGT\t0/0",
    "chr1\t2000\t.\tA\tG\t.\t.\t.\tGT\t1/1",
    "chr1\t3000\t.\tA\tG\t.\t.\t.\tGT\t0/0",
    "chr1\t4000\t.\tA\tG\t.\t.\t.\tGT\t1/1",
    "chr1\t5000\t.\tA\tG\t.\t.\t.\tGT\t0/0",
    "chr1\t6000\t.\tA\tG\t.\t.\t.\tGT\t1/1",
    "chr1\t7000\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    "chr1\t8000\t.\tA\tG\t.\t.\t.\tGT\t0/0",
    "chr1\t9000\t.\tA\tG\t.\t.\t.\tGT\t1/1",
    "chr1\t10000\t.\tA\tG\t.\t.\t.\tGT\t0/0",
]) + "\n"

# the defaults (25 SNPs, 1 Mb) are for genome-scale data; relaxed here for the toy example
result = calculate_runs_of_homozygosity(vcf, options={"min_snps": 5, "min_length_bp": 3000})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "calculate_runs_of_homozygosity",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 12,
    "text_length": 361,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "13481682ccc52cdb8c6279ba8202488902c4a2975beda2550094f35953e2d450"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sample": "S1",
    "method": "greedy left-to-right runs of homozygous genotypes tolerating a bounded number of heterozyg...",
    "parameters": {"min_snps": 5, "min_length_bp": 3000, "max_heterozygotes": 1, "max_missing": 2, "max_gap_bp": 1000000},
    "runs": [
      {
        "chrom": "chr1",
        "start": 1000,
        "end": 10000,
        "length_bp": 9001,
        "n_snps": 10,
        "n_heterozygous": 1,
        "n_missing": 0
      }
    ],
    "n_runs": 1,
    "total_roh_bp": 9001
  }
}
```

## 21) calculate_transition_transversion

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `calculate_transition_transversion` 功能。候选编号 CAP-0779，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_transition_transversion` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.calculate_transition_transversion(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref` 和 `alt`，基于基因型的统计还需显式提供 `genotypes`。

- **输出：** 本地函数 `calculate_transition_transversion` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import calculate_transition_transversion

result = calculate_transition_transversion(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0779",
  "operation": "calculate_transition_transversion",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_transition_transversion', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'classification': 'transition'}], 'transitions': 1, 'transversions': 0, 'excluded': 0, 'ti_tv_ratio': None}, 'provenance': {'implementation': 'dna...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 22) calculate_variant_density

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `calculate_variant_density` 功能。候选编号 CAP-0780，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_variant_density` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.calculate_variant_density(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref` 和 `alt`，基于基因型的统计还需显式提供 `genotypes`。计算变异密度时必须提供 `options.region_length`。

- **输出：** 本地函数 `calculate_variant_density` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import calculate_variant_density

result = calculate_variant_density(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]},
    options={"region_length": 100}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0780",
  "operation": "calculate_variant_density",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_variant_density', 'result': {'variant_count': 1, 'region_length': 100, 'density_per_base': 0.01, 'density_per_megabase': 10000.0, 'by_chromosome': {'chr1': 1}}, 'provenance': {'implementation': 'dnakit.native_variant_statistics', 'backend': 'local_d...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 23) call

- **作用：** 从测序 reads 中检出变异。尚未实现：函数只报告缺少什么（外部变异检出程序与比对后的 reads）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。请使用 `call_freebayes`、`call_gatk_haplotypecaller`、`call_deepvariant`、`call_delly`、`call_cutesv` 或 `call_sniffles2`，这些函数在对应程序已安装时调用外部程序。

- **API：** `dnakit.variants.call(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"variant caller binary (bcftools/FreeBayes/GATK/DeepVariant) and an aligned BAM/CRAM with a reference"）与 `message`（"variant calling needs an external caller and aligned reads; use the call_<caller> functions."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import call

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
  "module": "variants",
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
  "required_dependency": "variant caller binary (bcftools/FreeBayes/GATK/DeepVariant) and an aligned BAM/CRAM with a reference",
  "message": "variant calling needs an external caller and aligned reads; use the call_<caller> functions."
}
```

## 24) call_mt

- **作用：** 从 reads 中检出线粒体变异。尚未实现：函数只报告缺少什么（支持线粒体的变异检出程序（Mutect2 线粒体模式、mutserve）与比对后的 reads）。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。函数不会运行任何检出程序。

- **API：** `dnakit.variants.call_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"mitochondrial variant caller (Mutect2 --mitochondria-mode, mutserve) and aligned reads"）与 `message`（"mitochondrial variant calling needs an external caller and aligned reads."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import call_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = call_mt(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "call_mt",
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
  "required_dependency": "mitochondrial variant caller (Mutect2 --mitochondria-mode, mutserve) and aligned reads",
  "message": "mitochondrial variant calling needs an external caller and aligned reads."
}
```

## 25) compare

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `compare` 功能。候选编号 CAP-0783，来自 14 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `compare` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `variant_type` 用于选择一种支持的变体（取值见源目录记录）。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`generic_operation_name`。已有同名函数证据：`compare|dnakit.similarity.compare`。

- **API：** `dnakit.variants.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含 `left` 和 `right` 两组变异的映射；按 CHROM/POS/REF/ALT 完全一致进行匹配，不做规范化。

- **输出：** 本地函数 `compare` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import compare

result = compare(
    input={"left": {"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}]}, "right": {"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}, {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T"}]}}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0783",
  "operation": "compare",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'compare', 'result': {'left_record_count': 1, 'right_record_count': 2, 'left_allele_count': 1, 'right_allele_count': 2, 'shared_allele_count': 1, 'only_left': [], 'only_right': [{'chrom': 'chr1', 'pos': 20, 'ref': 'C', 'alt': 'T'}], 'shared': [{'chrom': 'chr1'...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 26) compare_mt

- **作用：** 比较两个变异集合：共有和私有等位基因、精确度、召回率、F1、Jaccard 与基因型一致性。

- **处理方法：** `compare_mt` 是 `compare_vcfs` 的“变异 x 样本表（matrix table 风格）”别名：它对 VCF 文本/路径/记录运行同一套实现，不读取 Hail MatrixTable 对象。左侧集合为真值/参考集合，右侧为待测集合。等位基因按每个 ALT 在 (CHROM, POS, REF, ALT)（`match: 'exact'`）或 (CHROM, POS)（`'position'`）上匹配。precision = 共有 / |right|，recall = 共有 / |left|，F1 为二者调和平均，Jaccard = 共有 / |并集|。对两个集合中都有的样本，在共有位点上计算基因型一致性（仅 `exact` 匹配），等位基因按序列比较且与顺序无关；`non_reference_concordance` 只统计任一基因型含非参考等位基因的位点。未给参考序列时，位置不同的等价插入缺失表示被视为不同等位基因；给出 `options.reference` 时先对两个集合做规范化（裁剪并左对齐）。这是集合比较，不是 vcfeval、hap.py 那样考虑单倍型的工具。

- **API：** `dnakit.variants.compare_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'left': <VCF|records>, 'right': <VCF|records>}`，或单个集合并把另一个放在 `options.other`。`options`：`other`、`match`（默认 `'exact'`，另有 `'position'`）、`reference`（纯序列、FASTA 文本/路径或 `{contig: sequence}`，配合 `reference_start`，默认 1）用于先规范化两个集合，`check_ref`（默认 `'warn'`，另有 `'error'`）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `match`、`left_count`、`right_count`、`shared_count`、`only_left_count`、`only_right_count`、`precision_right_vs_left`、`recall_right_vs_left`、`f1`、`jaccard`，列表 `shared`、`only_left`、`only_right`，`genotype_concordance`（每个共有样本：`compared_genotypes`、`concordant`、`discordant`、`genotype_concordance`、`non_reference_concordance`）以及 `normalized_with_reference`。只给一个集合又没有 `options.other` 时抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.variants import compare_mt

left = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    "chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/1",
    "chr1\t300\t.\tAT\tA\t.\t.\t.\tGT\t0/1",
]) + "\n"
right = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    "chr1\t300\t.\tAT\tA\t.\t.\t.\tGT\t1/1",
    "chr1\t400\t.\tG\tA\t.\t.\t.\tGT\t0/1",
]) + "\n"

result = compare_mt({"left": left, "right": right})  # left = truth set

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "compare_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 417,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d82b05108a1b36c57a3f1642b1b8a7f45cb9b0ba3b49730adc578adc27802a35"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "match": "exact",
    "left_count": 3,
    "right_count": 3,
    "shared_count": 2,
    "only_left_count": 1,
    "only_right_count": 1,
    "precision_right_vs_left": 0.666667,
    "recall_right_vs_left": 0.666667,
    "f1": 0.666667,
    "jaccard": 0.5,
    "shared": [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G"}, {"chrom": "chr1", "pos": 300, "ref": "AT", "alt": "A"}],
    "only_left": [{"chrom": "chr1", "pos": 200, "ref": "C", "alt": "T"}],
    "only_right": [{"chrom": "chr1", "pos": 400, "ref": "G", "alt": "A"}],
    "genotype_concordance": {
      "S1": {
        "compared_genotypes": 2,
        "concordant": 1,
        "discordant": 1,
        "genotype_concordance": 0.5,
        "non_reference_concordance": 0.5
      }
    },
    "normalized_with_reference": false
  }
}
```

## 27) compare_vcfs

- **作用：** 比较两个变异集合：共有和私有等位基因、精确度、召回率、F1、Jaccard 与基因型一致性。

- **处理方法：** 左侧集合为真值/参考集合，右侧为待测集合。等位基因按每个 ALT 在 (CHROM, POS, REF, ALT)（`match: 'exact'`）或 (CHROM, POS)（`'position'`）上匹配。precision = 共有 / |right|，recall = 共有 / |left|，F1 为二者调和平均，Jaccard = 共有 / |并集|。对两个集合中都有的样本，在共有位点上计算基因型一致性（仅 `exact` 匹配），等位基因按序列比较且与顺序无关；`non_reference_concordance` 只统计任一基因型含非参考等位基因的位点。未给参考序列时，位置不同的等价插入缺失表示被视为不同等位基因；给出 `options.reference` 时先对两个集合做规范化（裁剪并左对齐）。这是集合比较，不是 vcfeval、hap.py 那样考虑单倍型的工具。

- **API：** `dnakit.variants.compare_vcfs(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'left': <VCF|records>, 'right': <VCF|records>}`，或单个集合并把另一个放在 `options.other`。`options`：`other`、`match`（默认 `'exact'`，另有 `'position'`）、`reference`（纯序列、FASTA 文本/路径或 `{contig: sequence}`，配合 `reference_start`，默认 1）用于先规范化两个集合，`check_ref`（默认 `'warn'`，另有 `'error'`）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `match`、`left_count`、`right_count`、`shared_count`、`only_left_count`、`only_right_count`、`precision_right_vs_left`、`recall_right_vs_left`、`f1`、`jaccard`，列表 `shared`、`only_left`、`only_right`，`genotype_concordance`（每个共有样本：`compared_genotypes`、`concordant`、`discordant`、`genotype_concordance`、`non_reference_concordance`）以及 `normalized_with_reference`。只给一个集合又没有 `options.other` 时抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.variants import compare_vcfs

left = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    "chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/1",
    "chr1\t300\t.\tAT\tA\t.\t.\t.\tGT\t0/1",
]) + "\n"
right = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    "chr1\t300\t.\tAT\tA\t.\t.\t.\tGT\t1/1",
    "chr1\t400\t.\tG\tA\t.\t.\t.\tGT\t0/1",
]) + "\n"

result = compare_vcfs({"left": left, "right": right})  # left = truth set

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "compare_vcfs",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 417,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d82b05108a1b36c57a3f1642b1b8a7f45cb9b0ba3b49730adc578adc27802a35"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "match": "exact",
    "left_count": 3,
    "right_count": 3,
    "shared_count": 2,
    "only_left_count": 1,
    "only_right_count": 1,
    "precision_right_vs_left": 0.666667,
    "recall_right_vs_left": 0.666667,
    "f1": 0.666667,
    "jaccard": 0.5,
    "shared": [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G"}, {"chrom": "chr1", "pos": 300, "ref": "AT", "alt": "A"}],
    "only_left": [{"chrom": "chr1", "pos": 200, "ref": "C", "alt": "T"}],
    "only_right": [{"chrom": "chr1", "pos": 400, "ref": "G", "alt": "A"}],
    "genotype_concordance": {
      "S1": {
        "compared_genotypes": 2,
        "concordant": 1,
        "discordant": 1,
        "genotype_concordance": 0.5,
        "non_reference_concordance": 0.5
      }
    },
    "normalized_with_reference": false
  }
}
```

## 28) construct_haplotypes

- **作用：** 由已定相的基因型构建每个样本的单倍型（等位基因列表，给出参考序列时还有序列）。

- **处理方法：** 每个已定相 GT（`0|1`）的第一个等位基因为单倍型 1，第二个为单倍型 2。对每个样本和染色体，等位基因按位置顺序列出。给出 `options.reference` 时，窗口内的 REF 碱基被每个单倍型的等位基因替换，得到单倍型序列（REF 必须与参考一致，重叠变异会报错，符号等位基因被忽略并给出警告）。函数本身不做定相（见 `phase_genotypes`）：未定相的杂合基因型报错，或在 `on_unphased: 'skip'` 时跳过；缺失和非二倍体调用被跳过并给出警告。

- **API：** `dnakit.variants.construct_haplotypes(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）；基因型应已定相。`options`：`sample`/`samples`（默认全部样本）、`on_unphased`（默认 `'error'`，另有 `'skip'`）、`reference`（纯序列、FASTA 文本/路径或 `{contig: sequence}`）配合 `reference_start` 与 `window_start`/`window_end`（默认取变异跨度）以同时构建序列。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`haplotypes`（`{sample: {chrom: {'1': {...}, '2': {...}}}}`，每个单倍型有 `alleles`，为 `{pos, ref, allele}` 列表，给出参考序列时还有 `sequence`）和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import construct_haplotypes

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\t1|1",
    "chr1\t5\t.\tA\tG\t.\tPASS\t.\tGT\t1|0\t0|0",
    "chr1\t8\t.\tT\tC\t.\tPASS\t.\tGT\t0|0\t0|1",
]) + "\n"

# reference positions: A1 C2 G3 T4 A5 C6 G7 T8 A9 C10; sequences are built over 2..8
result = construct_haplotypes(vcf, options={"reference": "ACGTACGTAC", "sample": "S1"})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "construct_haplotypes",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 5,
    "text_length": 172,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "94af7e21244e9ed9a82bd69757fe1112631c1f2ae81aadf937d83f5c4e4f8377"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "per-sample haplotype alleles taken from the phased GT (first allele = haplotype 1); option...",
    "haplotypes": {
      "S1": {
        "chr1": {
          "1": {
            "alleles": [
              {"pos": 2, "ref": "C", "allele": "C"},
              {"pos": 5, "ref": "A", "allele": "G"},
              {"pos": 8, "ref": "T", "allele": "T"}
            ],
            "sequence": "CGTGCGT"
          },
          "2": {
            "alleles": [
              {"pos": 2, "ref": "C", "allele": "T"},
              {"pos": 5, "ref": "A", "allele": "A"},
              {"pos": 8, "ref": "T", "allele": "T"}
            ],
            "sequence": "TGTACGT"
          }
        }
      }
    }
  }
}
```

## 29) count

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `count` 功能。候选编号 CAP-0787，来自 14 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `count` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `variant_type` 用于选择一种支持的变体（取值见源目录记录）。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.count(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF 记录的序列或映射；计数基于显式的 CHROM/POS/REF/ALT 字段。

- **输出：** 本地函数 `count` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import count

result = count(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0787",
  "operation": "count",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'count', 'result': {'record_count': 1, 'alternate_allele_count': 1, 'by_chromosome': {'chr1': 1}, 'allele_type_counts': {'snp': 1}}, 'provenance': {'implementation': 'dnakit.native_variant_statistics', 'backend': 'local_dependency_free', 'coordinate_system': '...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 30) count_mt

- **作用：** 统计输入的记录数（文本的非空行数或列表的元素个数）。尽管名称带 `mt`，它不会把计数限制在线粒体变异上，也不读取 VCF 的结构。

- **处理方法：** 在共享本地运行时（`src/dnakit/_candidate_runtime.py` 中的 `run_candidate_operation` 与辅助函数 `_local_result`）中统计记录数：文本（或字节，或 `Path` 的文件内容）按非空行拆分，每行一条记录；列表按元素个数计；含有键 `records`、`items`、`values`、`data`、`sequences`、`reads` 之一的映射按该键的值同样方式计数，其他映射算 1 条记录。不做任何解析或校验：VCF 的头部行（`##...`、`#CHROM`）和其他染色体上的记录与任何其他行一样被计数，所以头部两行加一个变异的 VCF 得到 3；`None` 算 1 条记录，空字符串为 0。若需要按 CHROM/POS/REF/ALT 记录区分 VCF 结构、给出按染色体和变异类型分类的计数，请用 `dnakit.variants.count`。

- **API：** `dnakit.variants.count_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本、`bytes`、已存在文件的 `Path`（读取前 1,000,000 字节）、记录列表或映射（见处理方法）。`options` 可以传入但不使用。

- **输出：** 标准封装：`status`（`"ok"`）、`module`（`"variants"`）、`operation`（`"count_mt"`）、`provenance`（`implementation` 为 `"local_candidate_adapter"`，另有 `backend`）、`input_summary`（通用文本摘要：`record_count` 为非空行数，另有 `text_length`、`sequence_bases`（字母个数）、`lengths`、`sha256`）、`warnings`（空）和 `result`。`result` 是整数，即记录数。

- **示例代码：**

```python
import json
from dnakit.variants import count_mt

result = count_mt("chrM\t73\t.\tA\tG\nchrM\t263\t.\tA\tG\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "count_mt",
  "provenance": {
    "implementation": "local_candidate_adapter",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 29,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "af190e6e8fb98b032e563306294f3551e543d377ca9a2b055fa3af98f76693e4"
  },
  "warnings": [],
  "status": "ok",
  "result": 2
}
```

## 31) detect_breakpoints

- **作用：** 从比对结果中检测结构变异断点。尚未实现：函数只报告缺少什么（分裂/不一致比对 reads 与 SV 检出程序（DELLY、Manta、GRIDSS））。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。`call_delly`、`call_cutesv`、`call_sniffles2` 在程序已安装时会运行此类工具。

- **API：** `dnakit.variants.detect_breakpoints(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"SV caller (DELLY/Manta/GRIDSS) and aligned reads"）与 `message`（"breakpoint detection needs read alignments and an SV caller; call_delly, call_cutesv and call_sniffles2 run such tools."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import detect_breakpoints

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = detect_breakpoints(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "detect_breakpoints",
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
  "required_dependency": "SV caller (DELLY/Manta/GRIDSS) and aligned reads",
  "message": "breakpoint detection needs read alignments and an SV caller; call_delly, call_cutesv and c..."
}
```

## 32) detect_low_frequency_variants

- **作用：** 利用逐样本等位基因深度和二项错误率检验检测低频变异。

- **处理方法：** 对每个带 `FORMAT/AD` 的样本和 ALT 等位基因（多等位位点先拆分），当 `min_vaf <= alt/depth <= max_vaf`、`depth >= min_depth`、`alt >= min_alt_count` 且单侧二项概率 P(X >= alt | depth, error_rate) <= alpha 时报告该调用。整次运行使用同一个错误率，没有链偏倚、序列上下文、碱基质量或 UMI 模型，因此只是启发式筛查而非经过验证的低 VAF 检出工具；请按检测的实测背景设置 `error_rate`。

- **API：** `dnakit.variants.detect_low_frequency_variants(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：带逐样本 `FORMAT/AD`（ref,alt）的 VCF 文本/路径或记录。`options`：`min_vaf`（0.01）、`max_vaf`（0.2）、`min_depth`（100）、`min_alt_count`（3）、`error_rate`（0.001）、`alpha`（0.01）、`samples`/`sample`（默认全部）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`parameters`、`n_tested`（深度足够的样本-位点对数）、`n_detected` 和 `variants`（每项含 `site`、`chrom`、`pos`、`ref`、`alt`、`sample`、`alt_count`、`depth`、`vaf`、`p_value_vs_error`）。没有任何调用带 AD 时 `warnings` 会提示。

- **示例代码：**

```python
import json

from dnakit.variants import detect_low_frequency_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/0:995,5\t0/1:50,50",
    "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/0:1000,0\t0/0:998,2",
    "chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT:AD\t0/0:300,12\t0/0:40,1",
]) + "\n"

result = detect_low_frequency_variants(vcf)  # defaults: depth >= 100, 1% <= VAF <= 20%, >= 3 ALT reads

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "detect_low_frequency_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 5,
    "text_length": 224,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4fcd1dcadbc0bb575ea2cfc614e6d321fe5cadfa4d7763786d68656459193fe9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "calls where min_vaf <= alt/depth <= max_vaf, alt count >= min_alt_count, depth >= min_dept...",
    "parameters": {"min_vaf": 0.01, "max_vaf": 0.2, "min_depth": 100, "min_alt_count": 3, "error_rate": 0.001, "alpha": 0.01},
    "n_tested": 5,
    "n_detected": 1,
    "variants": [
      {
        "site": "chr1:300:G:A",
        "chrom": "chr1",
        "pos": 300,
        "ref": "G",
        "alt": "A",
        "sample": "S1",
        "alt_count": 12,
        "depth": 312,
        "vaf": 0.0384615,
        "p_value_vs_error": 1.08712e-15
      }
    ]
  }
}
```

## 33) detect_mobile_element_insertions

- **作用：** 从结构变异 VCF 中挑出移动元件插入记录（Alu、LINE-1、SVA、HERV）。

- **处理方法：** 当 ALT 是 `<INS:ME:ALU>` 之类的符号等位基因、`INFO/SVTYPE` 是已知元件（`ALU`、`L1`、`LINE1`、`SVA`、`HERV`、`HERVK`）或存在 `INFO/MEINFO` 时，该记录被判为移动元件插入；家族名会被规范化（如 `L1` 变为 `LINE-1`）。函数只是挑选 SV/MEI 检出程序（MELT、Scramble 等）已经给出的调用，不检查比对，也不自行搜索插入。

- **API：** `dnakit.variants.detect_mobile_element_insertions(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：带符号 ALT 或 `INFO/SVTYPE`/`MEINFO` 的 SV VCF 文本/路径或记录。`options`：`families`（可选列表，如 `['Alu', 'LINE-1']`；默认全部）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`n_records`、`n_insertions` 和 `insertions`（每项含 `chrom`、`pos`、`id`、`alt`、`family`、`svlen`、`filter`、`qual`）。

- **示例代码：**

```python
import json

from dnakit.variants import detect_mobile_element_insertions

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t1000\tsv1\tN\t<INS:ME:ALU>\t60\tPASS\tSVTYPE=INS;SVLEN=300",
    "chr1\t2000\tsv2\tN\t<INS:ME:L1>\t50\tPASS\tSVTYPE=INS;SVLEN=1500",
    "chr1\t3000\tsv3\tN\t<DEL>\t40\tPASS\tSVTYPE=DEL;SVLEN=-500",
]) + "\n"

result = detect_mobile_element_insertions(vcf, options={"families": ["Alu"]})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "detect_mobile_element_insertions",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 5,
    "text_length": 228,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "23103667ee2bb74d5020ced8f2e18261c7bf4e483df3da110656b25803022793"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "selection of mobile-element-insertion records from an SV VCF: ALT <INS:ME:FAMILY>, INFO/SV...",
    "n_records": 3,
    "n_insertions": 1,
    "insertions": [
      {
        "chrom": "chr1",
        "pos": 1000,
        "id": "sv1",
        "alt": "<INS:ME:ALU>",
        "family": "Alu",
        "svlen": "300",
        "filter": "PASS",
        "qual": 60.0
      }
    ]
  }
}
```

## 34) detect_mosaic_variants

- **作用：** 检测候选嵌合变异：ALT 比例显著高于测序错误且显著低于种系的 0.5。

- **处理方法：** 对每个带 `FORMAT/AD` 的样本和 ALT 等位基因，当 `min_vaf <= alt/depth <= max_vaf`、`depth >= min_depth`、`alt >= min_alt_count`、单侧上尾二项概率 P(X >= alt | depth, error_rate) <= alpha 且下尾概率 P(X <= alt | depth, 0.5) <= alpha 时报告。这是启发式方法：不应用人群频率、链偏倚、正常样本池或亲代基因型过滤，且单一的 `error_rate` 并非由数据实测。

- **API：** `dnakit.variants.detect_mosaic_variants(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：带逐样本 `FORMAT/AD`（ref,alt）的 VCF 文本/路径或记录。`options`：`min_vaf`（0.03）、`max_vaf`（0.4）、`min_depth`（20）、`min_alt_count`（3）、`error_rate`（0.002）、`alpha`（0.01）、`samples`/`sample`（默认全部）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`parameters`、`n_tested`、`n_detected` 和 `variants`（每项含 `site`、`chrom`、`pos`、`ref`、`alt`、`sample`、`alt_count`、`depth`、`vaf`、`p_value_vs_error`、`p_value_vs_heterozygous`）。

- **示例代码：**

```python
import json

from dnakit.variants import detect_mosaic_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:80,20\t0/1:50,50",
    "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/0:99,1\t0/0:98,2",
    "chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT:AD\t0/1:45,55\t0/1:60,40",
]) + "\n"

result = detect_mosaic_variants(vcf)  # defaults: depth >= 20, 3% <= VAF <= 40%, >= 3 ALT reads

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "detect_mosaic_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 5,
    "text_length": 221,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "62da7281e20b1bd52f80c4f9a982472fe23027aa997553fe9c0cc010bfbff1e7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "mosaic candidates: ALT fraction in [min_vaf, max_vaf], enough ALT reads, and the ALT count...",
    "parameters": {"min_vaf": 0.03, "max_vaf": 0.4, "min_depth": 20, "min_alt_count": 3, "error_rate": 0.002, "alpha": 0.01},
    "n_tested": 6,
    "n_detected": 1,
    "variants": [
      {
        "site": "chr1:100:A:G",
        "chrom": "chr1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
        "sample": "S1",
        "alt_count": 20,
        "depth": 100,
        "vaf": 0.2,
        "p_value_vs_error": 4.82527e-34,
        "p_value_vs_heterozygous": 5.57954e-10
      }
    ]
  }
}
```

## 35) detect_repeat_expansions

- **作用：** 在 DNA 序列中查找指定重复单元最长的连续串联重复，并标记超过阈值的扩增。

- **处理方法：** 对每条序列和每个 motif，用正则表达式找出 motif 精确、不间断拷贝的最长连续区段（任意相位；出现一次中断即结束）。当该区段至少有 `expansion_threshold` 个单元时 `expansion` 为 true。启用 `both_strands` 时同时搜索 motif 的反向互补。函数读取用户提供的位点序列，不像 ExpansionHunter 那样从测序 reads 中对重复分型，也不容忍不纯的重复。

- **API：** `dnakit.variants.detect_repeat_expansions(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：DNA 字符串、FASTA 文本、字符串列表或 `{name: sequence}`（每条序列至多 50 Mb）。`options`：`motif`（必需；1-12 个 ACGT 碱基，或列表）、`expansion_threshold`（单元数；可选）、`both_strands`（默认 False）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`expansion_threshold` 和 `results`；每条序列、每个 motif、每条链一项，含 `sequence`、`motif`、`strand`、`longest_pure_run_units`、`start`、`end`（1 起始；motif 从未出现时为 `null`）和 `expansion`（未给阈值时为 `null`）。缺少或无效的 `motif` 会抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.variants import detect_repeat_expansions

seq = "GGC" + "CAG" * 12 + "CAA" + "CAG" * 5 + "TTT"

result = detect_repeat_expansions(seq, options={"motif": "CAG", "expansion_threshold": 10})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "detect_repeat_expansions",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 60,
    "sequence_bases": 60,
    "lengths": [60],
    "sha256": "ffe9dfbc3a65a2a5108399d186ed70cdd7013fe8582ec07b76f6992253373167"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "longest uninterrupted tandem run of the repeat unit in each sequence (exact matches, every...",
    "expansion_threshold": 10,
    "results": [
      {
        "sequence": "sequence_1",
        "motif": "CAG",
        "strand": "+",
        "longest_pure_run_units": 12,
        "start": 4,
        "end": 39,
        "expansion": true
      }
    ]
  }
}
```

## 36) filter

- **作用：** 按质量、深度、缺失率、等位基因频率、HWE、变异类型和区域等标准过滤变异记录（位点级和基因型级）。

- **处理方法：** 每条记录依次接受已设置的各项标准的检验；任一检验失败则记录被移除（或被标记），失败原因分类计数（`not_pass`、`low_qual`、`low_depth`、`high_depth`、`high_missing`、`af_out_of_range`、`hwe`、`variant_type`、`chrom`、`region`）。缺少 `QUAL` 或深度时对应阈值判为不通过，`pass_only` 接受 `PASS` 或 `.`。基因型级阈值不移除记录，而是把不合格的调用置为缺失（缺少 `GQ`/`DP` 的调用同样不合格），之后缺失率、等位基因频率和 HWE 检验看到的是掩盖后的基因型。HWE 使用 `calculate_hwe_pvalue` 的精确检验。结果可通过 `output` 写成 VCF。

- **API：** `dnakit.variants.filter(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`（均可选）：`pass_only`、`min_qual`、`min_depth`、`max_depth`（位点 `DP`）、`min_genotype_quality` 与 `min_sample_depth`（低于阈值的调用被置为 `./.`）、`max_missing`（比例）、`min_af`、`max_af`（由基因型计算的 ALT 频率，否则取 `INFO/AF`）、`min_hwe_pvalue`、`variant_types`（`SNV`、`MNV`、`indel`、`SV`）、`chroms`、`regions`（`'chr:start-end'`）、`mode`（默认 `'remove'`，或 `'flag'`：不丢弃记录，而把未通过的原因写入 `FILTER`）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`（过滤后的 VCF 文本）、`output_path`、`records`（`include_records` 为 False 时省略；每条记录含 `chrom`、`pos`、`id`、`ref`、`alt`、`qual`、`filter`、`info`、`genotypes`）、`input_records`、`output_records`、`failed_records`、`reason_counts`、`genotypes_masked` 和 `mode`。

- **示例代码：**

```python
import json

from dnakit.variants import filter

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP:GQ\t0/0:10:50\t0/1:10:40\t1/1:10:30",
    "chr1\t200\t.\tC\tT\t20\tq10\t.\tGT\t0/1\t0/1\t./.",
    "chr1\t300\t.\tAT\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t0/0",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\t.\tGT\t1/2\t0/1\t0/0",
]) + "\n"

result = filter(vcf, options={"pass_only": True, "min_qual": 30, "max_missing": 0.2})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "filter",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 266,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4d0c6896b7e941e2a0cee1d3f019a90a0bddf5eb0e55588f80561ddb9152227c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 3,
    "samples": ["S1", "S2", "S3"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tG...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {"DP": "30"},
        "genotypes": {
          "S1": {"GT": "0/0", "DP": "10", "GQ": "50"},
          "S2": {"GT": "0/1", "DP": "10", "GQ": "40"},
          "S3": {"GT": "1/1", "DP": "10", "GQ": "30"}
        }
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "AT",
        "alt": ["A"],
        "qual": 60.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/0"}, "S2": {"GT": "0/1"}, "S3": {"GT": "0/0"}}
      },
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["A", "C"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "1/2"}, "S2": {"GT": "0/1"}, "S3": {"GT": "0/0"}}
      }
    ],
    "input_records": 4,
    "output_records": 3,
    "failed_records": 1,
    "reason_counts": {"high_missing": 1, "low_qual": 1, "not_pass": 1},
    "genotypes_masked": 0,
    "mode": "remove"
  }
}
```

## 37) filter_mt

- **作用：** 按质量、深度、缺失率、等位基因频率、HWE、类型和区域等标准过滤变异 x 样本表（matrix table 风格）的调用集。

- **处理方法：** `filter_mt` 是 `filter` 的“变异 x 样本表（matrix table 风格）”别名：它对 VCF 文本/路径/记录运行同一套实现，不读取 Hail MatrixTable 对象。每条记录依次接受已设置的各项标准的检验；任一检验失败则记录被移除（或被标记），失败原因分类计数（`not_pass`、`low_qual`、`low_depth`、`high_depth`、`high_missing`、`af_out_of_range`、`hwe`、`variant_type`、`chrom`、`region`）。缺少 `QUAL` 或深度时对应阈值判为不通过，`pass_only` 接受 `PASS` 或 `.`。基因型级阈值不移除记录，而是把不合格的调用置为缺失（缺少 `GQ`/`DP` 的调用同样不合格），之后缺失率、等位基因频率和 HWE 检验看到的是掩盖后的基因型。HWE 使用 `calculate_hwe_pvalue` 的精确检验。结果可通过 `output` 写成 VCF。

- **API：** `dnakit.variants.filter_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`（均可选）：`pass_only`、`min_qual`、`min_depth`、`max_depth`（位点 `DP`）、`min_genotype_quality` 与 `min_sample_depth`（低于阈值的调用被置为 `./.`）、`max_missing`（比例）、`min_af`、`max_af`（由基因型计算的 ALT 频率，否则取 `INFO/AF`）、`min_hwe_pvalue`、`variant_types`（`SNV`、`MNV`、`indel`、`SV`）、`chroms`、`regions`（`'chr:start-end'`）、`mode`（默认 `'remove'`，或 `'flag'`：不丢弃记录，而把未通过的原因写入 `FILTER`）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`（过滤后的 VCF 文本）、`output_path`、`records`（`include_records` 为 False 时省略；每条记录含 `chrom`、`pos`、`id`、`ref`、`alt`、`qual`、`filter`、`info`、`genotypes`）、`input_records`、`output_records`、`failed_records`、`reason_counts`、`genotypes_masked` 和 `mode`。

- **示例代码：**

```python
import json

from dnakit.variants import filter_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP:GQ\t0/0:10:50\t0/1:10:40\t1/1:10:30",
    "chr1\t200\t.\tC\tT\t20\tq10\t.\tGT\t0/1\t0/1\t./.",
    "chr1\t300\t.\tAT\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t0/0",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\t.\tGT\t1/2\t0/1\t0/0",
]) + "\n"

result = filter_mt(vcf, options={"pass_only": True, "min_qual": 30, "max_missing": 0.2})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "filter_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 266,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4d0c6896b7e941e2a0cee1d3f019a90a0bddf5eb0e55588f80561ddb9152227c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 3,
    "samples": ["S1", "S2", "S3"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tG...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {"DP": "30"},
        "genotypes": {
          "S1": {"GT": "0/0", "DP": "10", "GQ": "50"},
          "S2": {"GT": "0/1", "DP": "10", "GQ": "40"},
          "S3": {"GT": "1/1", "DP": "10", "GQ": "30"}
        }
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "AT",
        "alt": ["A"],
        "qual": 60.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/0"}, "S2": {"GT": "0/1"}, "S3": {"GT": "0/0"}}
      },
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["A", "C"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "1/2"}, "S2": {"GT": "0/1"}, "S3": {"GT": "0/0"}}
      }
    ],
    "input_records": 4,
    "output_records": 3,
    "failed_records": 1,
    "reason_counts": {"high_missing": 1, "low_qual": 1, "not_pass": 1},
    "genotypes_masked": 0,
    "mode": "remove"
  }
}
```

## 38) filter_variants_by_allele_frequency

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `filter_variants_by_allele_frequency` 功能。候选编号 CAP-0796，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `filter_variants_by_allele_frequency` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.filter_variants_by_allele_frequency(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref` 和 `alt`，基于基因型的统计还需显式提供 `genotypes`。可选的 `options.min_frequency`、`max_frequency` 和 `allele_index` 用于指定频率范围。

- **输出：** 本地函数 `filter_variants_by_allele_frequency` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_allele_frequency

result = filter_variants_by_allele_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]},
    options={"min_frequency": 0.1, "max_frequency": 0.9}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0796",
  "operation": "filter_variants_by_allele_frequency",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_allele_frequency', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'genotypes': ['0/1', '1/1', '0/0']}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'allele_index': 1, 'frequency': 0.5, 'kept': True}], 'input_recor...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 39) filter_variants_by_depth

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `filter_variants_by_depth` 功能。候选编号 CAP-0797，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `filter_variants_by_depth` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.filter_variants_by_depth(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref`、`alt` 以及显式的 `DP`/`depth` 值。必须提供 `options.min_depth`；可选的 `max_depth` 和 `keep_missing` 用于限定 DP 值。

- **输出：** 本地函数 `filter_variants_by_depth` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_depth

result = filter_variants_by_depth(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "DP": 30}]},
    options={"min_depth": 10, "max_depth": 100}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0797",
  "operation": "filter_variants_by_depth",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_depth', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'DP': 30}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'depth': 30.0, 'kept': True}], 'input_record_count': 1, 'output_record_count': 1, 'removed_record_coun...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 40) filter_variants_by_genotype_quality

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `filter_variants_by_genotype_quality` 功能。候选编号 CAP-0798，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `filter_variants_by_genotype_quality` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.filter_variants_by_genotype_quality(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要显式的 `GQ`、`genotype_quality`、`genotype_qualities` 或样本 GQ 值。必须提供 `options.min_genotype_quality` 或 `min_gq`；GQ 值取自显式记录或样本。

- **输出：** 本地函数 `filter_variants_by_genotype_quality` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_genotype_quality

result = filter_variants_by_genotype_quality(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "GQ": 45}]},
    options={"min_genotype_quality": 20}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0798",
  "operation": "filter_variants_by_genotype_quality",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_genotype_quality', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'GQ': 45}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'genotype_qualities': [45.0], 'minimum_observed': 45.0, 'kept': True}], 'input_record_count...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 41) filter_variants_by_missingness

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `filter_variants_by_missingness` 功能。候选编号 CAP-0799，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `filter_variants_by_missingness` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.filter_variants_by_missingness(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref`、`alt`，以及显式的 `genotypes` 或 `missing_count/sample_count`。必须提供 `options.max_missingness`，取值在 0 到 1 之间。

- **输出：** 本地函数 `filter_variants_by_missingness` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_missingness

result = filter_variants_by_missingness(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "./.", "0/0"]}]},
    options={"max_missingness": 0.5}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0799",
  "operation": "filter_variants_by_missingness",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_missingness', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'genotypes': ['0/1', './.', '0/0']}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'missingness': 0.3333333333333333, 'kept': True}], 'input_record_count...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 42) filter_variants_by_quality

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `filter_variants_by_quality` 功能。候选编号 CAP-0800，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `filter_variants_by_quality` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.filter_variants_by_quality(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF `variants` 的映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref`、`alt` 以及有限的 `qual`/`QUAL` 值。必须提供 `options.min_quality`；除非设置 `keep_missing=True`，否则排除缺失的 QUAL 值。

- **输出：** 本地函数 `filter_variants_by_quality` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_quality

result = filter_variants_by_quality(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "qual": 50}]},
    options={"min_quality": 20}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0800",
  "operation": "filter_variants_by_quality",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_quality', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'qual': 50}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'quality': 50.0, 'kept': True}], 'input_record_count': 1, 'output_record_count': 1, 'removed_recor...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 43) genotype_structural_variants

- **作用：** 利用 reads 对样本中已知的结构变异分型。尚未实现：函数只报告缺少什么（比对 reads 与 SV 分型程序（SVTyper、Paragraph、GraphTyper））。

- **处理方法：** 不做任何计算，也不编造结果；函数不检查数据内容，直接返回如实的 `not_implemented` 状态。函数不会运行任何分型程序。

- **API：** `dnakit.variants.genotype_structural_variants(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 任意值，仅在 `input_summary` 中做摘要；`options` 被接受但不使用。

- **输出：** 始终返回 `status: "not_implemented"`，并带 `required_dependency`（"SV genotyper (SVTyper/Paragraph) and aligned reads"）与 `message`（"SV genotyping needs read alignments and an external genotyper."），没有 `result`。封装中同样包含 `module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary` 和 `warnings`。

- **示例代码：**

```python
import json

from dnakit.variants import genotype_structural_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = genotype_structural_variants(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "genotype_structural_variants",
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
  "required_dependency": "SV genotyper (SVTyper/Paragraph) and aligned reads",
  "message": "SV genotyping needs read alignments and an external genotyper."
}
```

## 44) fill_missing_genotypes_with_mode

- **作用：** 用同一位点上最常见的已调用基因型填补缺失的基因型调用（朴素的逐位点众数填补）。

- **处理方法：** 对每条记录统计已调用基因型（忽略等位基因顺序，`1/0` 计为 `0/1`）；每个缺失调用（`./.`）都被赋予该位点最常见的基因型，计数相同时按基因型字符串排序取舍，并标记 `FORMAT/IMPUTED=1`。没有任何已调用基因型的位点保持缺失，并计入 `sites_without_any_call`。这不是 Beagle 或 IMPUTE5 意义上的基因型填充（imputation）：没有参考面板、没有连锁不平衡、没有单倍型模型，也不考虑亲缘关系，因此稀有变异携带者会被系统性地替换为多数基因型，填补的调用也不带不确定性。它只是为无法处理缺失值的工具快速去除缺失；正式分析请使用真正的填充程序。

- **API：** `dnakit.variants.fill_missing_genotypes_with_mode(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）；需要样本基因型。`options`：`output`（可选，写出 VCF 的路径）、`overwrite`（默认 False）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`imputed_genotypes`（被填补的调用数）、`sites_without_any_call`、`records`、`vcf`（填补后的 VCF 文本）和 `output_path`。没有样本列的输入会抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.variants import fill_missing_genotypes_with_mode

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t./.\t0/1",
    "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT\t./.\t./.\t0/0",
]) + "\n"

result = fill_missing_genotypes_with_mode(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "fill_missing_genotypes_with_mode",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 5,
    "text_length": 196,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "09e61722a5752bcc121f14eeac1a1019fd4d1f2abe1186409b50bdb06f7a06fb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "naive per-site modal-genotype imputation: each missing call receives the most frequent cal...",
    "imputed_genotypes": 3,
    "sites_without_any_call": 0,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/0"}, "S2": {"GT": "0/1"}, "S3": {"GT": "1/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 200,
        "id": null,
        "ref": "C",
        "alt": ["T"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/1"}, "S2": {"GT": "0/1", "IMPUTED": "1"}, "S3": {"GT": "0/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "G",
        "alt": ["A"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/0", "IMPUTED": "1"}, "S2": {"GT": "0/0", "IMPUTED": "1"}, "S3": {"GT": "0/0"}}
      }
    ],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tG...",
    "output_path": null
  }
}
```

## 45) left_align

- **作用：** 依据参考序列对变异记录中的插入缺失做左对齐与裁剪。

- **处理方法：** 采用 Tan、Abecasis 与 Kang（2015）的算法：只要各等位基因共有末位碱基就去掉它；某个等位基因变空时，所有等位基因向左补一个参考碱基；重复到不再变化，最后去掉共有的前导碱基。多等位记录联合处理，所有等位基因一起移动。必须提供参考序列：每条记录都会在 `options.reference` 中查找，REF 与参考不符的记录给出警告并保持原样（`check_ref: 'error'` 时抛出异常）。左移到所给参考序列的起点为止。`hgvs_accession` 还会写出 HGVS 基因组字符串，其按 HGVS 标准做 3' 移位，串联插入使用 `dup`。

- **API：** `dnakit.variants.left_align(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`reference`（纯序列、FASTA 文本/路径或 `{contig: sequence}`；左对齐必需）、`reference_start`（纯序列第一个碱基的 1 起始位置，默认 1）、`check_ref`（默认 `'warn'`：REF 与参考不符的记录保持不变；或 `'error'`）、`split_multiallelic`（默认 False）、`sort`（默认 True）、`hgvs_accession`（登录号字符串或 `{contig: accession}`，为每个 ALT 添加 HGVS `g.` 字符串）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`、`output_path`、`records`（给出 `hgvs_accession` 时带 `hgvs_g`）、`changed_records`、`unchanged_or_skipped`、`skipped_records`、`left_aligned`（未使用参考序列时为 False）和 `changes`（`chrom`、`from`、`to`）。`warnings` 列出与参考不符的记录。缺少 `options.reference` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.variants import left_align

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},    # deletion of CA, written at the right end
    {"chrom": "chr1", "pos": 7, "ref": "A", "alt": "ACA"},    # insertion of CA, written at the right end
]

result = left_align(records, options={"reference": reference})

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "left_align",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "96ebf655d007619126a3fc9b6e332bee1d199e16bf4cd031730ea9f980c84ada"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tG\tGCA\t.\t.\t.\nchr1\t1\t.\tGCA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "G", "alt": ["GCA"], "qual": null, "filter": ".", "info": {}},
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "GCA", "alt": ["G"], "qual": null, "filter": ".", "info": {}}
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 A>ACA", "to": "1 G>GCA"}]
  }
}
```

## 46) left_align_mt

- **作用：** 依据参考序列对变异 x 样本表（matrix table 风格）调用集中的插入缺失做左对齐与裁剪。

- **处理方法：** `left_align_mt` 是 `left_align` 的“变异 x 样本表（matrix table 风格）”别名：它对 VCF 文本/路径/记录运行同一套实现，不读取 Hail MatrixTable 对象。采用 Tan、Abecasis 与 Kang（2015）的算法：只要各等位基因共有末位碱基就去掉它；某个等位基因变空时，所有等位基因向左补一个参考碱基；重复到不再变化，最后去掉共有的前导碱基。多等位记录联合处理，所有等位基因一起移动。必须提供参考序列：每条记录都会在 `options.reference` 中查找，REF 与参考不符的记录给出警告并保持原样（`check_ref: 'error'` 时抛出异常）。左移到所给参考序列的起点为止。`hgvs_accession` 还会写出 HGVS 基因组字符串，其按 HGVS 标准做 3' 移位，串联插入使用 `dup`。

- **API：** `dnakit.variants.left_align_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`reference`（纯序列、FASTA 文本/路径或 `{contig: sequence}`；左对齐必需）、`reference_start`（纯序列第一个碱基的 1 起始位置，默认 1）、`check_ref`（默认 `'warn'`：REF 与参考不符的记录保持不变；或 `'error'`）、`split_multiallelic`（默认 False）、`sort`（默认 True）、`hgvs_accession`（登录号字符串或 `{contig: accession}`，为每个 ALT 添加 HGVS `g.` 字符串）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`、`output_path`、`records`（给出 `hgvs_accession` 时带 `hgvs_g`）、`changed_records`、`unchanged_or_skipped`、`skipped_records`、`left_aligned`（未使用参考序列时为 False）和 `changes`（`chrom`、`from`、`to`）。`warnings` 列出与参考不符的记录。缺少 `options.reference` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.variants import left_align_mt

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},    # deletion of CA, written at the right end
    {"chrom": "chr1", "pos": 7, "ref": "A", "alt": "ACA"},    # insertion of CA, written at the right end
]

result = left_align_mt(records, options={"reference": reference})

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "left_align_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "96ebf655d007619126a3fc9b6e332bee1d199e16bf4cd031730ea9f980c84ada"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tG\tGCA\t.\t.\t.\nchr1\t1\t.\tGCA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "G", "alt": ["GCA"], "qual": null, "filter": ".", "info": {}},
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "GCA", "alt": ["G"], "qual": null, "filter": ".", "info": {}}
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 A>ACA", "to": "1 G>GCA"}]
  }
}
```

## 47) left_align_variants

- **作用：** 依据参考序列对变异记录中的插入缺失做左对齐与裁剪。

- **处理方法：** 采用 Tan、Abecasis 与 Kang（2015）的算法：只要各等位基因共有末位碱基就去掉它；某个等位基因变空时，所有等位基因向左补一个参考碱基；重复到不再变化，最后去掉共有的前导碱基。多等位记录联合处理，所有等位基因一起移动。必须提供参考序列：每条记录都会在 `options.reference` 中查找，REF 与参考不符的记录给出警告并保持原样（`check_ref: 'error'` 时抛出异常）。左移到所给参考序列的起点为止。`hgvs_accession` 还会写出 HGVS 基因组字符串，其按 HGVS 标准做 3' 移位，串联插入使用 `dup`。

- **API：** `dnakit.variants.left_align_variants(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`reference`（纯序列、FASTA 文本/路径或 `{contig: sequence}`；左对齐必需）、`reference_start`（纯序列第一个碱基的 1 起始位置，默认 1）、`check_ref`（默认 `'warn'`：REF 与参考不符的记录保持不变；或 `'error'`）、`split_multiallelic`（默认 False）、`sort`（默认 True）、`hgvs_accession`（登录号字符串或 `{contig: accession}`，为每个 ALT 添加 HGVS `g.` 字符串）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`、`output_path`、`records`（给出 `hgvs_accession` 时带 `hgvs_g`）、`changed_records`、`unchanged_or_skipped`、`skipped_records`、`left_aligned`（未使用参考序列时为 False）和 `changes`（`chrom`、`from`、`to`）。`warnings` 列出与参考不符的记录。缺少 `options.reference` 时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.variants import left_align_variants

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},    # deletion of CA, written at the right end
    {"chrom": "chr1", "pos": 7, "ref": "A", "alt": "ACA"},    # insertion of CA, written at the right end
]

result = left_align_variants(records, options={"reference": reference})

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "left_align_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "96ebf655d007619126a3fc9b6e332bee1d199e16bf4cd031730ea9f980c84ada"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tG\tGCA\t.\t.\t.\nchr1\t1\t.\tGCA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "G", "alt": ["GCA"], "qual": null, "filter": ".", "info": {}},
      {"chrom": "chr1", "pos": 1, "id": null, "ref": "GCA", "alt": ["G"], "qual": null, "filter": ".", "info": {}}
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 A>ACA", "to": "1 G>GCA"}]
  }
}
```

## 48) merge

- **作用：** 合并多个 VCF/调用集，默认合并样本列，也可只串联记录。

- **处理方法：** `mode: 'samples'`（默认，类似 `bcftools merge -m none`）把 CHROM/POS/REF/ALT 完全相同的记录合并为一条，各输入的样本并排排列；样本名必须唯一，某输入中不存在的样本记为 `./.`。`mode: 'records'`（类似 `bcftools concat`）要求样本列表一致，重复位点只保留一条。合并后的 QUAL 取最大值；只有所有参与记录均为 `PASS` 时 FILTER 才是 `PASS`，否则取未通过标签的并集；INFO 取第一条含该键的记录；输出已排序。GT 之外的 FORMAT 字段按样本保留；不做等位基因统一或重新规范化，写法不同的等位基因仍是不同记录。

- **API：** `dnakit.variants.merge(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：至少两个调用集的列表（VCF 文本、VCF 路径或 `{'records': [...]}`），或 `{'inputs': [...]}`。`options`：`mode`（默认 `'samples'`，另有 `'records'`）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`（合并后的 VCF 文本）、`output_path`、`records` 和 `mode`。调用集少于两个、`samples` 模式下样本名重复或 `records` 模式下样本列表不一致时抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.variants import merge

vcf_a = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
    "chr1\t300\t.\tT\tC\t40\tPASS\t.\tGT\t1/1",
]) + "\n"
vcf_b = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2",
    "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
    "chr1\t200\t.\tC\tT\t10\tPASS\t.\tGT\t0/1",
]) + "\n"

result = merge([vcf_a, vcf_b])  # default mode "samples": one column per input sample

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "merge",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 338,
    "sequence_bases": 146,
    "lengths": [73, 73],
    "sha256": "ca994842191b4f8afcf8dcbeb6d46b23e2b17b99c6fa7d6b175fab7f3018422a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 3,
    "samples": ["S1", "S2"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t100\t.\tA\tG\t50...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/1"}, "S2": {"GT": "1/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 200,
        "id": null,
        "ref": "C",
        "alt": ["T"],
        "qual": 10.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "./."}, "S2": {"GT": "0/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "T",
        "alt": ["C"],
        "qual": 40.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "1/1"}, "S2": {"GT": "./."}}
      }
    ],
    "mode": "samples"
  }
}
```

## 49) merge_mt

- **作用：** 合并多个变异 x 样本表（matrix table 风格）调用集，合并样本或串联记录。

- **处理方法：** `merge_mt` 是 `merge` 的“变异 x 样本表（matrix table 风格）”别名：它对 VCF 文本/路径/记录运行同一套实现，不读取 Hail MatrixTable 对象。`mode: 'samples'`（默认，类似 `bcftools merge -m none`）把 CHROM/POS/REF/ALT 完全相同的记录合并为一条，各输入的样本并排排列；样本名必须唯一，某输入中不存在的样本记为 `./.`。`mode: 'records'`（类似 `bcftools concat`）要求样本列表一致，重复位点只保留一条。合并后的 QUAL 取最大值；只有所有参与记录均为 `PASS` 时 FILTER 才是 `PASS`，否则取未通过标签的并集；INFO 取第一条含该键的记录；输出已排序。GT 之外的 FORMAT 字段按样本保留；不做等位基因统一或重新规范化，写法不同的等位基因仍是不同记录。

- **API：** `dnakit.variants.merge_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：至少两个调用集的列表（VCF 文本、VCF 路径或 `{'records': [...]}`），或 `{'inputs': [...]}`。`options`：`mode`（默认 `'samples'`，另有 `'records'`）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`（合并后的 VCF 文本）、`output_path`、`records` 和 `mode`。调用集少于两个、`samples` 模式下样本名重复或 `records` 模式下样本列表不一致时抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.variants import merge_mt

vcf_a = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
    "chr1\t300\t.\tT\tC\t40\tPASS\t.\tGT\t1/1",
]) + "\n"
vcf_b = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2",
    "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
    "chr1\t200\t.\tC\tT\t10\tPASS\t.\tGT\t0/1",
]) + "\n"

result = merge_mt([vcf_a, vcf_b])  # default mode "samples": one column per input sample

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "merge_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 338,
    "sequence_bases": 146,
    "lengths": [73, 73],
    "sha256": "ca994842191b4f8afcf8dcbeb6d46b23e2b17b99c6fa7d6b175fab7f3018422a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 3,
    "samples": ["S1", "S2"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t100\t.\tA\tG\t50...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/1"}, "S2": {"GT": "1/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 200,
        "id": null,
        "ref": "C",
        "alt": ["T"],
        "qual": 10.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "./."}, "S2": {"GT": "0/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "T",
        "alt": ["C"],
        "qual": 40.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "1/1"}, "S2": {"GT": "./."}}
      }
    ],
    "mode": "samples"
  }
}
```

## 50) merge_vcf_files

- **作用：** 合并多个 VCF 文件，默认合并样本列，也可只串联记录。

- **处理方法：** `mode: 'samples'`（默认，类似 `bcftools merge -m none`）把 CHROM/POS/REF/ALT 完全相同的记录合并为一条，各输入的样本并排排列；样本名必须唯一，某输入中不存在的样本记为 `./.`。`mode: 'records'`（类似 `bcftools concat`）要求样本列表一致，重复位点只保留一条。合并后的 QUAL 取最大值；只有所有参与记录均为 `PASS` 时 FILTER 才是 `PASS`，否则取未通过标签的并集；INFO 取第一条含该键的记录；输出已排序。GT 之外的 FORMAT 字段按样本保留；不做等位基因统一或重新规范化，写法不同的等位基因仍是不同记录。

- **API：** `dnakit.variants.merge_vcf_files(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：至少两个调用集的列表（VCF 文本、VCF 路径或 `{'records': [...]}`），或 `{'inputs': [...]}`。`options`：`mode`（默认 `'samples'`，另有 `'records'`）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`（合并后的 VCF 文本）、`output_path`、`records` 和 `mode`。调用集少于两个、`samples` 模式下样本名重复或 `records` 模式下样本列表不一致时抛出 `InputFormatError`。

- **示例代码：**

```python
import json

from dnakit.variants import merge_vcf_files

vcf_a = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
    "chr1\t300\t.\tT\tC\t40\tPASS\t.\tGT\t1/1",
]) + "\n"
vcf_b = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2",
    "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
    "chr1\t200\t.\tC\tT\t10\tPASS\t.\tGT\t0/1",
]) + "\n"

result = merge_vcf_files([vcf_a, vcf_b])  # default mode "samples": one column per input sample

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "merge_vcf_files",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 338,
    "sequence_bases": 146,
    "lengths": [73, 73],
    "sha256": "ca994842191b4f8afcf8dcbeb6d46b23e2b17b99c6fa7d6b175fab7f3018422a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 3,
    "samples": ["S1", "S2"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t100\t.\tA\tG\t50...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 100,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": 50.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "0/1"}, "S2": {"GT": "1/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 200,
        "id": null,
        "ref": "C",
        "alt": ["T"],
        "qual": 10.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "./."}, "S2": {"GT": "0/1"}}
      },
      {
        "chrom": "chr1",
        "pos": 300,
        "id": null,
        "ref": "T",
        "alt": ["C"],
        "qual": 40.0,
        "filter": "PASS",
        "info": {},
        "genotypes": {"S1": {"GT": "1/1"}, "S2": {"GT": "./."}}
      }
    ],
    "mode": "samples"
  }
}
```

## 51) normalize

- **作用：** 规范化变异记录：裁掉共有碱基，给出参考序列时再对插入缺失左对齐（类似 bcftools norm）。

- **处理方法：** 采用 Tan、Abecasis 与 Kang（2015）的算法：只要各等位基因共有末位碱基就去掉它；某个等位基因变空时，所有等位基因向左补一个参考碱基；重复到不再变化，最后去掉共有的前导碱基。多等位记录联合处理，所有等位基因一起移动。未给 `options.reference` 时只做裁剪（插入缺失保留一个锚定碱基），并报告 `left_aligned` 为 false。给出参考序列时，REF 与参考不符的记录给出警告并保持不变（`check_ref`），`split_multiallelic` 会在规范化之后再为每个 ALT 写一条记录，`hgvs_accession` 则添加做过 3' 移位的 HGVS 基因组字符串。除非 `sort` 为 False，输出按位置排序。只使用用户提供的参考序列，不会下载基因组。

- **API：** `dnakit.variants.normalize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`reference`（纯序列、FASTA 文本/路径或 `{contig: sequence}`；此处可选）、`reference_start`（纯序列第一个碱基的 1 起始位置，默认 1）、`check_ref`（默认 `'warn'`：REF 与参考不符的记录保持不变；或 `'error'`）、`split_multiallelic`（默认 False）、`sort`（默认 True）、`hgvs_accession`（登录号字符串或 `{contig: accession}`，为每个 ALT 添加 HGVS `g.` 字符串）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`、`output_path`、`records`（给出 `hgvs_accession` 时带 `hgvs_g`）、`changed_records`、`unchanged_or_skipped`、`skipped_records`、`left_aligned`（未使用参考序列时为 False）和 `changes`（`chrom`、`from`、`to`）。`warnings` 列出与参考不符的记录。未给参考序列时只做裁剪。

- **示例代码：**

```python
import json

from dnakit.variants import normalize

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},
    {"chrom": "chr1", "pos": 7, "ref": "AT", "alt": "GT"},    # trailing T is shared and is trimmed
]

result = normalize(records, options={"reference": reference, "hgvs_accession": "NC_000001.11"})

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "normalize",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2a10b4bc4a63171a786de098b3cf4d8b1edaa19ac0703ce86b5f47c2cd6b8fab"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tGCA\tG\t.\t.\t.\nchr1\t7\t.\tA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 1,
        "id": null,
        "ref": "GCA",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.6_7del"]
      },
      {
        "chrom": "chr1",
        "pos": 7,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.7A>G"]
      }
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 AT>GT", "to": "7 A>G"}]
  }
}
```

## 52) normalize_mt

- **作用：** 规范化变异 x 样本表（matrix table 风格）调用集：裁掉共有碱基，给出参考序列时再对插入缺失左对齐。

- **处理方法：** `normalize_mt` 是 `normalize` 的“变异 x 样本表（matrix table 风格）”别名：它对 VCF 文本/路径/记录运行同一套实现，不读取 Hail MatrixTable 对象。采用 Tan、Abecasis 与 Kang（2015）的算法：只要各等位基因共有末位碱基就去掉它；某个等位基因变空时，所有等位基因向左补一个参考碱基；重复到不再变化，最后去掉共有的前导碱基。多等位记录联合处理，所有等位基因一起移动。未给 `options.reference` 时只做裁剪（插入缺失保留一个锚定碱基），并报告 `left_aligned` 为 false。给出参考序列时，REF 与参考不符的记录给出警告并保持不变（`check_ref`），`split_multiallelic` 会在规范化之后再为每个 ALT 写一条记录，`hgvs_accession` 则添加做过 3' 移位的 HGVS 基因组字符串。除非 `sort` 为 False，输出按位置排序。只使用用户提供的参考序列，不会下载基因组。

- **API：** `dnakit.variants.normalize_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`reference`（纯序列、FASTA 文本/路径或 `{contig: sequence}`；此处可选）、`reference_start`（纯序列第一个碱基的 1 起始位置，默认 1）、`check_ref`（默认 `'warn'`：REF 与参考不符的记录保持不变；或 `'error'`）、`split_multiallelic`（默认 False）、`sort`（默认 True）、`hgvs_accession`（登录号字符串或 `{contig: accession}`，为每个 ALT 添加 HGVS `g.` 字符串）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`、`output_path`、`records`（给出 `hgvs_accession` 时带 `hgvs_g`）、`changed_records`、`unchanged_or_skipped`、`skipped_records`、`left_aligned`（未使用参考序列时为 False）和 `changes`（`chrom`、`from`、`to`）。`warnings` 列出与参考不符的记录。未给参考序列时只做裁剪。

- **示例代码：**

```python
import json

from dnakit.variants import normalize_mt

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},
    {"chrom": "chr1", "pos": 7, "ref": "AT", "alt": "GT"},    # trailing T is shared and is trimmed
]

result = normalize_mt(records, options={"reference": reference, "hgvs_accession": "NC_000001.11"})

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "normalize_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2a10b4bc4a63171a786de098b3cf4d8b1edaa19ac0703ce86b5f47c2cd6b8fab"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tGCA\tG\t.\t.\t.\nchr1\t7\t.\tA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 1,
        "id": null,
        "ref": "GCA",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.6_7del"]
      },
      {
        "chrom": "chr1",
        "pos": 7,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.7A>G"]
      }
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 AT>GT", "to": "7 A>G"}]
  }
}
```

## 53) normalize_vcf

- **作用：** 规范化变异记录：裁掉共有碱基，给出参考序列时再对插入缺失左对齐（类似 bcftools norm）。

- **处理方法：** 采用 Tan、Abecasis 与 Kang（2015）的算法：只要各等位基因共有末位碱基就去掉它；某个等位基因变空时，所有等位基因向左补一个参考碱基；重复到不再变化，最后去掉共有的前导碱基。多等位记录联合处理，所有等位基因一起移动。未给 `options.reference` 时只做裁剪（插入缺失保留一个锚定碱基），并报告 `left_aligned` 为 false。给出参考序列时，REF 与参考不符的记录给出警告并保持不变（`check_ref`），`split_multiallelic` 会在规范化之后再为每个 ALT 写一条记录，`hgvs_accession` 则添加做过 3' 移位的 HGVS 基因组字符串。除非 `sort` 为 False，输出按位置排序。只使用用户提供的参考序列，不会下载基因组。

- **API：** `dnakit.variants.normalize_vcf(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`reference`（纯序列、FASTA 文本/路径或 `{contig: sequence}`；此处可选）、`reference_start`（纯序列第一个碱基的 1 起始位置，默认 1）、`check_ref`（默认 `'warn'`：REF 与参考不符的记录保持不变；或 `'error'`）、`split_multiallelic`（默认 False）、`sort`（默认 True）、`hgvs_accession`（登录号字符串或 `{contig: accession}`，为每个 ALT 添加 HGVS `g.` 字符串）、`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`samples`、`vcf`、`output_path`、`records`（给出 `hgvs_accession` 时带 `hgvs_g`）、`changed_records`、`unchanged_or_skipped`、`skipped_records`、`left_aligned`（未使用参考序列时为 False）和 `changes`（`chrom`、`from`、`to`）。`warnings` 列出与参考不符的记录。未给参考序列时只做裁剪。

- **示例代码：**

```python
import json

from dnakit.variants import normalize_vcf

# reference positions: G1 C2 A3 C4 A5 C6 A7 T8
reference = "GCACACAT"
records = [
    {"chrom": "chr1", "pos": 5, "ref": "ACA", "alt": "A"},
    {"chrom": "chr1", "pos": 7, "ref": "AT", "alt": "GT"},    # trailing T is shared and is trimmed
]

result = normalize_vcf(records, options={"reference": reference, "hgvs_accession": "NC_000001.11"})

print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "variants",
  "operation": "normalize_vcf",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 110,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2a10b4bc4a63171a786de098b3cf4d8b1edaa19ac0703ce86b5f47c2cd6b8fab"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": [],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tGCA\tG\t.\t.\t.\nchr1\t7\t.\tA\tG\t.\t.\t.\n",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 1,
        "id": null,
        "ref": "GCA",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.6_7del"]
      },
      {
        "chrom": "chr1",
        "pos": 7,
        "id": null,
        "ref": "A",
        "alt": ["G"],
        "qual": null,
        "filter": ".",
        "info": {},
        "hgvs_g": ["NC_000001.11:g.7A>G"]
      }
    ],
    "changed_records": 2,
    "unchanged_or_skipped": 0,
    "skipped_records": 0,
    "left_aligned": true,
    "changes": [{"chrom": "chr1", "from": "5 ACA>A", "to": "1 GCA>G"}, {"chrom": "chr1", "from": "7 AT>GT", "to": "7 A>G"}]
  }
}
```

## 54) parse_vcf_records

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `parse_vcf_records` 功能。候选编号 CAP-0812，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `parse_vcf_records` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.parse_vcf_records(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是 VCF 文本、UTF-8 字节、VCF 行序列或包含 `text` 的映射；会校验 VCF 必需列。

- **输出：** 本地函数 `parse_vcf_records` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import parse_vcf_records

result = parse_vcf_records(
    input='##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\nchr1\t10\trs1\tA\tG,T\t50\tPASS\tDP=12\tGT\t0/1',
    options={"strict": True}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0812",
  "operation": "parse_vcf_records",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'parse_vcf_records', 'result': {'headers': ['##fileformat=VCFv4.3', '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\ts1'], 'samples': ['s1'], 'records': [{'chrom': 'chr1', 'pos': 10, 'id': 'rs1', 'ref': 'A', 'alt': ['G', 'T'], 'qual': 50.0, 'filter': 'P...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 55) phase_genotypes

- **作用：** 利用亲代基因型按传递关系对子代的杂合基因型定相（三体定相）。

- **处理方法：** 对子代每个杂合基因型 x/y，检查两种指派（x 来自父亲、y 来自母亲，或相反）是否与亲代等位基因相容；恰有一种指派相容时，把子代 GT 改写为 父源|母源。纯合的子代基因型只标记为已定相；两种指派都相容（如双亲均为杂合）、都不相容（孟德尔错误）以及子代基因型缺失的位点保持不变。亲代基因型缺失时对该亲本的等位基因不施加约束。这只是依据传递关系的确定性定相，不做统计定相（Beagle、SHAPEIT），因此在无亲缘关系的样本中大部分位点仍不会被定相。

- **API：** `dnakit.variants.phase_genotypes(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）；必须包含子代、父亲、母亲三个样本。`options`：`trio` = `{'child': ..., 'father': ..., 'mother': ...}`（也可直接给 `child`/`father`/`mother` 键；必需）、`output`（可选的 VCF 路径）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`trio`、`counts`（`phased`、`ambiguous`、`homozygous`、`mendelian_error`、`child_missing` 各类位点数；未出现的类别不列出）、`vcf`（含子代已定相基因型的 VCF 文本）和 `output_path`。缺少 `trio` 或样本名不存在时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json

from dnakit.variants import phase_genotypes

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother",
    "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/1\t0/0",
    "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0",
    "chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1",
    "chr1\t400\t.\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t0/1",
]) + "\n"

result = phase_genotypes(vcf, options={"trio": {"child": "child", "father": "father", "mother": "mother"}})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "phase_genotypes",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 243,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "267f2e692be1c3ef4bfe87ba4f02cbf34e14b5201bf1f1c67af2e98903ab95bb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "trio phasing by transmission: the child's heterozygous genotype is written paternal|matern...",
    "trio": {"child": "child", "father": "father", "mother": "mother"},
    "counts": {"phased": 2, "ambiguous": 1, "homozygous": 1},
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother\nchr...",
    "output_path": null
  }
}
```

## 56) plot_density

- **作用：** 把变异沿染色体的密度画成分箱 SVG 柱状图，并返回各 bin 的计数。

- **处理方法：** 按每条染色体将变异计入大小为 `bin_size` 的固定 bin（bin 序号 = (POS - 1) // bin_size）；染色体按自然顺序排列，给出 `chrom_lengths` 时会包含直到染色体末端的空 bin。每 Mb 密度 = 计数 * 1,000,000 / bin_size。SVG 是以文本构建的简单柱状图（柱高与计数成正比，按最大 bin 缩放）并带染色体标签，不需要绘图库。bin 数超过 20,000 会报错（请增大 `bin_size`）。

- **API：** `dnakit.variants.plot_density(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`bin_size`（默认 1000000）、`chroms`（要包含的 contig 列表）、`chrom_lengths`（`{contig: length}`）、`svg_output`（可选，写出 SVG 的路径）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `bin_size`、`bins`（每项含 `chrom`、`start`、`end`、`count`、`density_per_mb`）、`svg`（SVG 文本）、`output_path` 和 `max_count`。

- **示例代码：**

```python
import json

from dnakit.variants import plot_density

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
    "chr1\t500\t.\tA\tG\t50\tPASS\t.",
    "chr1\t1500\t.\tA\tG\t50\tPASS\t.",
    "chr1\t2100\t.\tA\tG\t50\tPASS\t.",
    "chr1\t2200\t.\tA\tG\t50\tPASS\t.",
    "chr2\t300\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = plot_density(vcf, options={"bin_size": 1000})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "plot_density",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 8,
    "text_length": 213,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1b43dee4028171f8b561a43f9b1b12264f39c2cfb625d601dd8e188346be3938"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "bin_size": 1000,
    "bins": [
      {"chrom": "chr1", "start": 1, "end": 1000, "count": 2, "density_per_mb": 2000.0},
      {"chrom": "chr1", "start": 1001, "end": 2000, "count": 1, "density_per_mb": 1000.0},
      "..."
    ],
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"320\" height=\"220\" viewBox=\"0 0 320 220\"><te...",
    "output_path": null,
    "max_count": 2
  }
}
```

## 57) plot_density_mt

- **作用：** 把变异 x 样本表（matrix table 风格）调用集的变异密度画成分箱 SVG 柱状图。

- **处理方法：** `plot_density_mt` 是 `plot_density` 的“变异 x 样本表（matrix table 风格）”别名：它对 VCF 文本/路径/记录运行同一套实现，不读取 Hail MatrixTable 对象。按每条染色体将变异计入大小为 `bin_size` 的固定 bin（bin 序号 = (POS - 1) // bin_size）；染色体按自然顺序排列，给出 `chrom_lengths` 时会包含直到染色体末端的空 bin。每 Mb 密度 = 计数 * 1,000,000 / bin_size。SVG 是以文本构建的简单柱状图（柱高与计数成正比，按最大 bin 缩放）并带染色体标签，不需要绘图库。bin 数超过 20,000 会报错（请增大 `bin_size`）。

- **API：** `dnakit.variants.plot_density_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`bin_size`（默认 1000000）、`chroms`（要包含的 contig 列表）、`chrom_lengths`（`{contig: length}`）、`svg_output`（可选，写出 SVG 的路径）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `bin_size`、`bins`（每项含 `chrom`、`start`、`end`、`count`、`density_per_mb`）、`svg`（SVG 文本）、`output_path` 和 `max_count`。

- **示例代码：**

```python
import json

from dnakit.variants import plot_density_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
    "chr1\t500\t.\tA\tG\t50\tPASS\t.",
    "chr1\t1500\t.\tA\tG\t50\tPASS\t.",
    "chr1\t2100\t.\tA\tG\t50\tPASS\t.",
    "chr1\t2200\t.\tA\tG\t50\tPASS\t.",
    "chr2\t300\t.\tA\tG\t50\tPASS\t.",
]) + "\n"

result = plot_density_mt(vcf, options={"bin_size": 1000})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "plot_density_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 8,
    "text_length": 213,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1b43dee4028171f8b561a43f9b1b12264f39c2cfb625d601dd8e188346be3938"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "bin_size": 1000,
    "bins": [
      {"chrom": "chr1", "start": 1, "end": 1000, "count": 2, "density_per_mb": 2000.0},
      {"chrom": "chr1", "start": 1001, "end": 2000, "count": 1, "density_per_mb": 1000.0},
      "..."
    ],
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"320\" height=\"220\" viewBox=\"0 0 320 220\"><te...",
    "output_path": null,
    "max_count": 2
  }
}
```

## 58) prioritize_candidate_variants

- **作用：** 利用记录上已有的注释，用透明的加和得分对候选变异排序。

- **处理方法：** 注释从记录的键或 INFO 中读取（`consequence`/`most_severe_consequence`/`effect` 或 INFO `MC`、`clinvar`/`CLNSIG`、`gnomad_af`、`cadd`、`revel`、`gene`/`GENEINFO`）；函数不做任何查询。得分是以下各项之和：后果严重度，取自固定的 Sequence Ontology 表（0-1：功能丧失类术语 1.0，`stop_lost` 0.8，框内变化 0.6，错义 0.5，剪接区域 0.3，同义 0.05 等；取最严重的术语）；ClinVar（致病 1.0，可能致病 0.8，意义不明 0.1，良性 -1.0）；稀有度（人群 AF 缺失或 <= 1e-4 为 0.5，至 1e-3 为 0.3，至 1e-2 为 0.1，其余为 0）；预测器，0.5 * max(CADD / 40, REVEL)；以及 `options.genes` 中的基因 +0.5。人群 AF 高于 `max_af` 的变异，以及 `exclude_benign` 为真时 ClinVar 良性/可能良性的变异，会单独列入排除项。这是权重人为设定的启发式排序，不是经过验证的分类器；没有任何人群频率注释的变异保留完整的稀有度得分（0.5），但会标记 `frequency_unknown: true`，并计入 `warnings` 中的一条汇总警告，从而能与确知稀有的变异区分开。每个得分分量都会返回，以便审查排序。

- **API：** `dnakit.variants.prioritize_candidate_variants(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`max_af`（0.05；高于它的变异被排除）、`exclude_benign`（True）、`genes`（基因面板列表，+0.5）、`top_n`（只保留得分最高的 n 个）。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `method`、`parameters`、`ranked`（每项含 `chrom`、`pos`、`ref`、`alt`、`id`、`score`、`components`、`excluded_reason`（为 `null`）、`gene`、`rank`）、`excluded`（字段相同，带 `excluded_reason`）、`n_ranked`（`top_n` 截断之前）、`n_excluded` 和 `n_frequency_unknown`（缺少人群频率注释的变异数，含已排除的；此时封装的 `warnings` 中有一条带该数量的消息）。每个 `ranked` 和 `excluded` 行都有布尔字段 `frequency_unknown`。

- **示例代码：**

```python
import json

from dnakit.variants import prioritize_candidate_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\tGENE=BRCA1;consequence=stop_gained;clinvar=Pathogenic;gnomad_af=0.00001;cadd=35",
    "chr1\t200\t.\tC\tT\t50\tPASS\tGENE=TTN;consequence=missense_variant;gnomad_af=0.002;revel=0.6",
    "chr1\t300\t.\tG\tA\t50\tPASS\tGENE=OBSCN;consequence=synonymous_variant;gnomad_af=0.2",
    "chr1\t400\t.\tT\tC\t50\tPASS\tGENE=APOE;consequence=missense_variant;clinvar=Benign",
]) + "\n"

result = prioritize_candidate_variants(vcf, options={"genes": ["BRCA1"]})

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "prioritize_candidate_variants",
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

## 59) split_multiallelic

- **作用：** 为“变异过滤、规范化、注释与统计”类分析提供可复用的 `split_multiallelic` 功能。候选编号 CAP-0817，来自 14 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `split_multiallelic` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `variant_type` 用于选择一种支持的变体（取值见源目录记录）。无外部依赖的实现位于 `src/dnakit/variants/statistics.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的类 VCF 记录计算描述性统计，不调用变异检测工具、注释数据库、参考基因组或基因型填补后端。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.variants.split_multiallelic(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 必须是包含类 VCF 记录的序列或映射；每条记录需要 `chrom`、从 1 开始的正整数 `pos`、`ref` 和 `alt`。

- **输出：** 本地函数 `split_multiallelic` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的来源信息。结果仅为描述性统计：不推断基因型、注释、参考等位基因或数据库信息。

- **示例代码：**

```python
from pathlib import Path

from dnakit.variants import split_multiallelic

result = split_multiallelic(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": ["G", "T"], "genotypes": ["0/1", "1/2"]}]},
    options={"remap_genotypes": True}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0817",
  "operation": "split_multiallelic",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'split_multiallelic', 'result': {'records': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'genotypes': ['0/1', '1/0'], 'original_alt_index': 1}, {'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'T', 'genotypes': ['0/0', '0/1'], 'original_alt_index': 2}]...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 60) split_multiallelic_mt

- **作用：** 把变异 x 样本表（matrix table 风格）调用集中的多等位位点拆分为每个 ALT 等位基因一条记录。

- **处理方法：** `split_multiallelic_mt` 是 `split_multiallelic_variants` 的“变异 x 样本表（matrix table 风格）”别名：它对 VCF 文本/路径/记录运行同一套实现，不读取 Hail MatrixTable 对象。多等位记录的每个 ALT 成为一条独立记录，CHROM/POS/REF 不变。GT 会被重新映射（选中的 ALT 记为 1，其余等位基因包括 REF 一律记为 0，因此 `1/2` 在第一条记录中变为 `1/0`，在第二条中变为 `0/1`）；`FORMAT/AD` 与 `INFO` 中的 `AC`/`AF` 只保留所选等位基因，`FORMAT/PL` 因无法子集化而被丢弃。不会裁剪或左对齐等位基因（需要时请用带 `split_multiallelic: True` 的 `normalize`）。双等位记录原样复制。

- **API：** `dnakit.variants.split_multiallelic_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`（拆分后）、`samples`、`vcf`、`output_path`、`records`、`input_records` 和 `split_records`（新增的记录数）。

- **示例代码：**

```python
import json

from dnakit.variants import split_multiallelic_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\tAC=2,1\tGT:AD\t1/2:5,3,2\t0/1:6,4,0",
]) + "\n"

result = split_multiallelic_mt(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "split_multiallelic_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 131,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2cdf597390051f5fe42bb240739a02d0f21405e197041be6028c4f9d58996f97"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": ["S1", "S2"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t400\t.\tG\tA\t30...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["A"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {"AC": "2"},
        "genotypes": {"S1": {"GT": "1/0", "AD": "5,3"}, "S2": {"GT": "0/1", "AD": "6,4"}}
      },
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["C"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {"AC": "1"},
        "genotypes": {"S1": {"GT": "0/1", "AD": "5,2"}, "S2": {"GT": "0/0", "AD": "6,0"}}
      }
    ],
    "input_records": 1,
    "split_records": 1
  }
}
```

## 61) split_multiallelic_variants

- **作用：** 把多等位变异记录拆分为每个 ALT 等位基因一条记录。

- **处理方法：** 多等位记录的每个 ALT 成为一条独立记录，CHROM/POS/REF 不变。GT 会被重新映射（选中的 ALT 记为 1，其余等位基因包括 REF 一律记为 0，因此 `1/2` 在第一条记录中变为 `1/0`，在第二条中变为 `0/1`）；`FORMAT/AD` 与 `INFO` 中的 `AC`/`AF` 只保留所选等位基因，`FORMAT/PL` 因无法子集化而被丢弃。不会裁剪或左对齐等位基因（需要时请用带 `split_multiallelic: True` 的 `normalize`）。双等位记录原样复制。

- **API：** `dnakit.variants.split_multiallelic_variants(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：`output`（VCF 路径）、`include_records`（默认 True）、`overwrite`。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`（拆分后）、`samples`、`vcf`、`output_path`、`records`、`input_records` 和 `split_records`（新增的记录数）。

- **示例代码：**

```python
import json

from dnakit.variants import split_multiallelic_variants

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\tAC=2,1\tGT:AD\t1/2:5,3,2\t0/1:6,4,0",
]) + "\n"

result = split_multiallelic_variants(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "split_multiallelic_variants",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 131,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2cdf597390051f5fe42bb240739a02d0f21405e197041be6028c4f9d58996f97"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 2,
    "samples": ["S1", "S2"],
    "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t400\t.\tG\tA\t30...",
    "output_path": null,
    "records": [
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["A"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {"AC": "2"},
        "genotypes": {"S1": {"GT": "1/0", "AD": "5,3"}, "S2": {"GT": "0/1", "AD": "6,4"}}
      },
      {
        "chrom": "chr1",
        "pos": 400,
        "id": null,
        "ref": "G",
        "alt": ["C"],
        "qual": 30.0,
        "filter": "PASS",
        "info": {"AC": "1"},
        "genotypes": {"S1": {"GT": "0/1", "AD": "5,2"}, "S2": {"GT": "0/0", "AD": "6,0"}}
      }
    ],
    "input_records": 1,
    "split_records": 1
  }
}
```

## 62) summarize

- **作用：** 汇总一个调用集：按变异类型和染色体的计数、Ti/Tv、FILTER 计数、多等位位点以及逐样本质控。

- **处理方法：** 按裁掉共有碱基后的等位基因长度，对每个 ALT 分类：`SNV`、`MNV`、`indel`（子类型为插入、缺失或复杂）和 `SV`（符号等位基因或长度变化至少 50 个碱基）。Ti/Tv 只统计 SNV（转换 A<->G、C<->T）。逐样本质控给出调用率、hom-ref/het/hom-alt 计数、het/hom-alt 比、SNV 与 indel 数、样本级 Ti/Tv 以及 `DP`、`GQ` 的平均值。这些都是文件内容的朴素计数，不做过滤，也不使用参考序列。

- **API：** `dnakit.variants.summarize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：无。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`n_alt_alleles`、`n_multiallelic_sites`、`n_samples`、`samples`、`type_counts`、`subtype_counts`、`by_chromosome`、`filter_counts`、`transitions`、`transversions`、`ti_tv`、`mean_qual`、`indel_length_counts`、`mean_sample_call_rate` 和 `sample_qc`（每个样本：`call_rate`、`n_hom_ref`、`n_het`、`n_hom_alt`、`het_hom_alt_ratio`、`n_snv`、`n_indel`、`ti_tv`、`mean_depth`、`mean_genotype_quality` 等）。分母为 0 时比值为 `null`。

- **示例代码：**

```python
import json

from dnakit.variants import summarize

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP:GQ\t0/0:10:50\t0/1:10:40\t1/1:10:30",
    "chr1\t200\t.\tC\tT\t20\tq10\t.\tGT\t0/1\t0/1\t./.",
    "chr1\t300\t.\tAT\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t0/0",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\t.\tGT\t1/2\t0/1\t0/0",
]) + "\n"

result = summarize(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "summarize",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 266,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4d0c6896b7e941e2a0cee1d3f019a90a0bddf5eb0e55588f80561ddb9152227c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 4,
    "n_alt_alleles": 5,
    "n_multiallelic_sites": 1,
    "n_samples": 3,
    "samples": ["S1", "S2", "S3"],
    "type_counts": {"SNV": 4, "indel": 1},
    "subtype_counts": {"indel:deletion": 1},
    "by_chromosome": {"chr1": 4},
    "filter_counts": {"PASS": 3, "q10": 1},
    "transitions": 3,
    "transversions": 1,
    "ti_tv": 3.0,
    "mean_qual": 40.0,
    "indel_length_counts": {"-1": 1},
    "mean_sample_call_rate": 0.916667,
    "sample_qc": [
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
        "n_snv": 3,
        "n_indel": 0,
        "n_transitions": 2,
        "n_transversions": 1,
        "ti_tv": 2.0,
        "mean_depth": 10.0,
        "mean_genotype_quality": 50.0
      },
      {
        "sample": "S2",
        "n_sites": 4,
        "n_called": 4,
        "n_missing": 0,
        "call_rate": 1.0,
        "n_hom_ref": 0,
        "n_het": 4,
        "n_hom_alt": 0,
        "het_hom_alt_ratio": null,
        "n_snv": 3,
        "n_indel": 1,
        "n_transitions": 3,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 10.0,
        "mean_genotype_quality": 40.0
      },
      {
        "sample": "S3",
        "n_sites": 4,
        "n_called": 3,
        "n_missing": 1,
        "call_rate": 0.75,
        "n_hom_ref": 2,
        "n_het": 0,
        "n_hom_alt": 1,
        "het_hom_alt_ratio": 0.0,
        "n_snv": 1,
        "n_indel": 0,
        "n_transitions": 1,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 10.0,
        "mean_genotype_quality": 30.0
      }
    ]
  }
}
```

## 63) summarize_mt

- **作用：** 汇总变异 x 样本表（matrix table 风格）调用集：类型计数、Ti/Tv、FILTER 计数及逐样本质控。

- **处理方法：** `summarize_mt` 是 `summarize` 的“变异 x 样本表（matrix table 风格）”别名：它对 VCF 文本/路径/记录运行同一套实现，不读取 Hail MatrixTable 对象。按裁掉共有碱基后的等位基因长度，对每个 ALT 分类：`SNV`、`MNV`、`indel`（子类型为插入、缺失或复杂）和 `SV`（符号等位基因或长度变化至少 50 个碱基）。Ti/Tv 只统计 SNV（转换 A<->G、C<->T）。逐样本质控给出调用率、hom-ref/het/hom-alt 计数、het/hom-alt 比、SNV 与 indel 数、样本级 Ti/Tv 以及 `DP`、`GQ` 的平均值。这些都是文件内容的朴素计数，不做过滤，也不使用参考序列。

- **API：** `dnakit.variants.summarize_mt(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：VCF 文本（制表符分隔，也可以是 gzip 字节）、VCF/VCF.gz 路径，或变异字典列表（`chrom`、`pos`、`ref`、`alt`，可选 `qual`、`filter`、`info` 以及 `genotypes`/`samples`）。`options`：无。

- **输出：** 返回标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 含 `n_records`、`n_alt_alleles`、`n_multiallelic_sites`、`n_samples`、`samples`、`type_counts`、`subtype_counts`、`by_chromosome`、`filter_counts`、`transitions`、`transversions`、`ti_tv`、`mean_qual`、`indel_length_counts`、`mean_sample_call_rate` 和 `sample_qc`（每个样本：`call_rate`、`n_hom_ref`、`n_het`、`n_hom_alt`、`het_hom_alt_ratio`、`n_snv`、`n_indel`、`ti_tv`、`mean_depth`、`mean_genotype_quality` 等）。分母为 0 时比值为 `null`。

- **示例代码：**

```python
import json

from dnakit.variants import summarize_mt

vcf = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:DP:GQ\t0/0:10:50\t0/1:10:40\t1/1:10:30",
    "chr1\t200\t.\tC\tT\t20\tq10\t.\tGT\t0/1\t0/1\t./.",
    "chr1\t300\t.\tAT\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t0/0",
    "chr1\t400\t.\tG\tA,C\t30\tPASS\t.\tGT\t1/2\t0/1\t0/0",
]) + "\n"

result = summarize_mt(vcf)

print(json.dumps(result, indent=2))
```

- **示例输出：** 较长的数组和字符串用 `...` 缩略，数值保留 6 位有效数字；输出由运行上面的代码得到。

```json
{
  "module": "variants",
  "operation": "summarize_mt",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 6,
    "text_length": 266,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4d0c6896b7e941e2a0cee1d3f019a90a0bddf5eb0e55588f80561ddb9152227c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_records": 4,
    "n_alt_alleles": 5,
    "n_multiallelic_sites": 1,
    "n_samples": 3,
    "samples": ["S1", "S2", "S3"],
    "type_counts": {"SNV": 4, "indel": 1},
    "subtype_counts": {"indel:deletion": 1},
    "by_chromosome": {"chr1": 4},
    "filter_counts": {"PASS": 3, "q10": 1},
    "transitions": 3,
    "transversions": 1,
    "ti_tv": 3.0,
    "mean_qual": 40.0,
    "indel_length_counts": {"-1": 1},
    "mean_sample_call_rate": 0.916667,
    "sample_qc": [
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
        "n_snv": 3,
        "n_indel": 0,
        "n_transitions": 2,
        "n_transversions": 1,
        "ti_tv": 2.0,
        "mean_depth": 10.0,
        "mean_genotype_quality": 50.0
      },
      {
        "sample": "S2",
        "n_sites": 4,
        "n_called": 4,
        "n_missing": 0,
        "call_rate": 1.0,
        "n_hom_ref": 0,
        "n_het": 4,
        "n_hom_alt": 0,
        "het_hom_alt_ratio": null,
        "n_snv": 3,
        "n_indel": 1,
        "n_transitions": 3,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 10.0,
        "mean_genotype_quality": 40.0
      },
      {
        "sample": "S3",
        "n_sites": 4,
        "n_called": 3,
        "n_missing": 1,
        "call_rate": 0.75,
        "n_hom_ref": 2,
        "n_het": 0,
        "n_hom_alt": 1,
        "het_hom_alt_ratio": 0.0,
        "n_snv": 1,
        "n_indel": 0,
        "n_transitions": 1,
        "n_transversions": 0,
        "ti_tv": null,
        "mean_depth": 10.0,
        "mean_genotype_quality": 30.0
      }
    ]
  }
}
```
