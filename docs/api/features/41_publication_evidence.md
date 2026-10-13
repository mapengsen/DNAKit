# 论文引用与变异统计证据

这一页规定 DNAKit 在论文中如何记录第三方软件，并展示首批无外部依赖的变异统计函数。软件可以被论文引用，并不等于其源码、模型、数据库或远程服务可以随 DNAKit 重新分发。

## 1) 查询软件引用记录

- **作用：** 查询 Biopython、SeqKit、scikit-bio、pydna、DNA Chisel、seqme、Edlib、Cutadapt 等工具的论文、许可证和集成边界。
- **处理方法：** 从随包固定的离线证据登记表读取 DOI、项目地址、许可证状态、引用状态和建议接入方式。
- **API：** `dnakit.integrations.describe_publication_record(name[必须])`。
- **输入：** 一个已登记的软件名或别名，例如 `"Biopython"`、`"SeqKit"` 或 `"minimap2"`。
- **输出：** 返回软件引用、许可证、来源链接、论文引用状态和 `reuse_route`；它不表示该软件的全部 API 已由 DNAKit 重导出。
- **示例代码：**

```python
from dnakit.integrations import describe_publication_record

record = describe_publication_record("Biopython", verbose=False)
print(record["doi"])
print(record["reuse_route"])
```

- **示例结果：**

```text
10.1093/bioinformatics/btp163
optional_pip_adapter
```

## 2) 检查论文证据是否完整

- **作用：** 检查软件版本、后端版本、数据库版本、访问日期、引用和许可证证据是否足够写入方法部分。
- **处理方法：** 分开检查论文引用证据、许可证证据和实验运行版本；许可证限制不会被 DOI 检查掩盖。
- **API：** `dnakit.integrations.publication_readiness(name[必须], version[可选], backend_version[可选], database_version[可选], accessed[可选])`。
- **输入：** 软件名和实际使用的版本；使用外部 CLI、模型或数据库时，同时填写相应版本。
- **输出：** 返回 `status`、逐项 `checks`、引用文本、许可证、来源链接和复用边界。`ready_for_manuscript_record` 只表示记录证据齐全，不表示科学结果已经验证。
- **示例代码：**

```python
from dnakit.integrations import publication_readiness

result = publication_readiness(
    "SeqKit",
    version="2.10.0",
    accessed="2026-10-12",
    verbose=False,
)
print(result["status"])
print(result["checks"]["citation_evidence"])
```

- **示例结果：**

```text
ready_for_manuscript_record
True
```

## 3) 生成软件引用记录

- **作用：** 生成可复制到论文方法或补充材料的软件引用记录。
- **处理方法：** 把软件名、版本、DOI 或项目地址、访问日期、许可证和外部后端版本组合成结构化记录。
- **API：** `dnakit.integrations.build_software_citation(name[必须], version[必须], accessed[可选], backend_version[可选], database_version[可选])`。
- **输入：** 实际运行的 DNAKit 或第三方软件版本；不要用“最新版”代替版本号。
- **输出：** 返回引用文本、DOI、许可证、来源 URL、复用边界和证据检查结果。
- **示例代码：**

```python
from dnakit.integrations import build_software_citation

record = build_software_citation(
    "scikit-bio",
    "0.7.1",
    accessed="2026-10-12",
    verbose=False,
)
print(record["citation"])
```

- **示例结果：**

```text
scikit-bio version 0.7.1; Aton et al. (2026), Scikit-bio: a fundamental Python library for biological omic data analysis. DOI: 10.1038/s41592-025-02981-z.
```

## 4) 计算等位基因频率

- **作用：** 从明确给出的 VCF 风格基因型计算指定替代等位基因频率。
- **处理方法：** 解析 `0/0`、`0/1`、`1/1`、`./.` 或等价整数数组；只统计已调用等位基因，不进行基因型推断或填补。
- **API：** `dnakit.variants.calculate_allele_frequency(input[必须], options[可选])`。
- **输入：** `input` 为包含 `variants` 的映射或变异映射序列；每条记录至少包含 `chrom`、`pos`、`ref`、`alt` 和 `genotypes`。
- **输出：** 返回每个位点的替代等位基因计数、已调用等位基因数和频率，并给出合计频率。
- **示例代码：**

```python
from dnakit.variants import calculate_allele_frequency

result = calculate_allele_frequency({
    "variants": [{
        "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G",
        "genotypes": ["0/1", "1/1", "0/0"],
    }]
}, options={})
print(result["result"]["aggregate_frequency"])
```

- **示例结果：**

```text
0.5
```

## 5) 计算观察杂合度

- **作用：** 计算明确二倍体基因型中的观察杂合度。
- **处理方法：** 对每个位点统计非缺失样本中两个等位基因不同的比例。
- **API：** `dnakit.variants.calculate_heterozygosity(input[必须], options[可选])`。
- **输入：** 与等位基因频率相同，但每个位点必须提供二倍体 `genotypes`。
- **输出：** 返回每个位点和全体位点的观察杂合度。
- **示例代码：**

```python
from dnakit.variants import calculate_heterozygosity

result = calculate_heterozygosity({
    "variants": [{
        "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G",
        "genotypes": ["0/1", "1/1", "0/0"],
    }]
})
print(result["result"]["observed_heterozygosity"])
```

- **示例结果：**

```text
0.3333333333333333
```

## 6) 计算缺失率

- **作用：** 计算每个位点的缺失样本比例。
- **处理方法：** 把含有 `.` 的基因型视为缺失，不进行插补。
- **API：** `dnakit.variants.calculate_missingness(input[必须], options[可选])`。
- **输入：** VCF 风格变异记录和样本基因型。
- **输出：** 返回缺失样本数、样本总数和缺失率。
- **示例代码：**

```python
from dnakit.variants import calculate_missingness

result = calculate_missingness({
    "variants": [{
        "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G",
        "genotypes": ["0/1", "./.", "0/0"],
    }]
})
print(result["result"]["missingness"])
```

- **示例结果：**

```text
0.3333333333333333
```

## 7) 计算转换/颠换比

- **作用：** 对 SNP 统计 transitions、transversions 和 Ti/Tv 比值。
- **处理方法：** A↔G、C↔T 计为 transition；其他单碱基替换计为 transversion；插入缺失被明确排除。
- **API：** `dnakit.variants.calculate_transition_transversion(input[必须], options[可选])`。
- **输入：** 变异记录中的 `ref` 和 `alt` 必须能区分单碱基替换与 indel。
- **输出：** 返回两类计数、排除数和 `ti_tv_ratio`。
- **示例代码：**

```python
from dnakit.variants import calculate_transition_transversion

result = calculate_transition_transversion({
    "variants": [
        {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"},
        {"chrom": "chr1", "pos": 20, "ref": "A", "alt": "C"},
    ]
})
print(result["result"]["ti_tv_ratio"])
```

- **示例结果：**

```text
1.0
```

## 8) 计算变异密度

- **作用：** 在给定区域长度下计算变异数和每 Mb 变异密度。
- **处理方法：** 以输入记录数除以显式提供的区域长度；不把参考基因组长度默认为已知。
- **API：** `dnakit.variants.calculate_variant_density(input[必须], options[必须])`。
- **输入：** 变异记录和 `options["region_length"]`（单位为碱基）。
- **输出：** 返回区域长度、变异数、每碱基密度、每 Mb 密度和按染色体分组的计数。
- **示例代码：**

```python
from dnakit.variants import calculate_variant_density

result = calculate_variant_density(
    {"variants": [
        {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"},
        {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T"},
    ]},
    options={"region_length": 1_000_000},
)
print(result["result"]["density_per_megabase"])
```

- **示例结果：**

```text
2.0
```

论文使用这些函数时，应同时记录 DNAKit 版本、输入数据版本、坐标约定、过滤规则和本页所列公式。它们是描述性统计，不替代变异检测、注释、群体模型或临床解释流程。

完整的 320 项来源逐条期刊资格审计见 [43) 320 项 DNA/基因工具的论文引用审计](43_publication_audit.md)。

## 9) 论文引用边界

- **作用：** 把“可以在论文中引用”与“可以随 DNAKit 重新分发”分开记录。
- **处理方法：** 优先使用登记表中的 DOI；同时记录实际软件版本、Python/CLI 后端版本、模型或数据库版本及访问日期。
- **输入：** 论文方法中实际运行的软件和资源，不要用项目名代替版本号。
- **输出：** `verified_doi` 表示已有 DOI 证据；`preprint_or_project` 表示必须明确写成预印本或项目引用；许可证状态独立判断源码、模型、数据库和远程服务的复用方式。
- **示例：** SeqKit 2.x 优先使用 SeqKit2 DOI `10.1002/imt2.191`；pyfaidx 使用 `10.7287/peerj.preprints.970v1` 并注明预印本；dnaio 使用软件归档 DOI `10.5281/zenodo.10548864`，不要把它写成同行评审论文。
