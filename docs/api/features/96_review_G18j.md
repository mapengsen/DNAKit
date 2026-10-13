# 96 G18j 复审功能：pyensembl、gget、biothings_client 与 genomepy

本页记录 G18j 组对四个包的复审结果。复审只处理台账中状态为 `other` 的 88 行：包装 5 个函数，其余 83 行各记录一个状态与一行理由。包装函数位于 `dnakit.databases` 的新模块中，可直接从模块导入，例如 `from dnakit.databases import ensembl_species_record`。

- 离线运行：5 个包装函数都不访问网络。URL 类函数只拼接或匹配字符串；物种函数读取 pyensembl 内置的物种注册表；genomepy 函数只读取基因组目录中已保存的 `README.txt` 与 `assembly_report.txt`，并在调用上游前检查两个文件是否存在，因此不会触发下载。
- gget 与 biothings_client 在本组没有可离线运行的 DNA 级函数，`gget_review_functions.py` 与 `biothings_client_review_functions.py` 的 `__all__` 为空，理由记录在台账中。
- 示例中的输入文件以 `示例文件` 标签给出，其余示例直接运行。

## 复审结果

| 包 | 复审前 other/cap | 新包装 | 其余状态 |
| --- | ---: | ---: | --- |
| pyensembl | 56 | 4 | family-variant 15，needs-data 4，not-dna 3，internal 30 |
| gget | 4 | 0 | not-dna 2，needs-network 1，cli 1 |
| biothings_client | 5 | 0 | internal 3，not-dna 2 |
| genomepy | 23 | 1 | family-variant 1，internal 18，not-dna 2，cli 1 |
| 合计 | 88 | 5 | family-variant 16，internal 51，not-dna 9，needs-data 4，needs-network 1，cli 2 |

`family-variant` 行是同一上游函数的再导出（模块别名），台账的 `dnakit_name` 列记录对应的包装名。

## 1) `ensembl_dated_release_urls`

- **作用：** 给出一个 dated 注释版本（Ensembl 新平台，如 `2026_04`）的 GTF、cDNA、蛋白与参考 DNA 下载地址。
- **处理方法：** 检查 GCA 或 GCF 登录号、来源名与注释日期（YYYY_MM），再调用 pyensembl 的 `make_dated_release_urls` 拼接地址。`include_alt` 只改变 GTF 文件名；`genome_fasta_mask` 决定参考 DNA 文件名（none、soft 或 hard）。不发起网络请求。
- **API：** `ensembl_dated_release_urls(assembly_accession[必须], provider[必须], annotation_date[必须], *, include_alt[可选], genome_fasta_mask[可选], verbose[可选])`
- **输入：** 登录号如 `GCA_000001405.29`；来源名如 `ensembl`；注释日期如 `2026_04`；`include_alt` 为布尔值，默认 False；`genome_fasta_mask` 取 none、soft 或 hard，默认 none。
- **输出：** 字典，含 `directory`、`gtf`、`cdna`、`pep`、`genome_fasta` 五个地址，以及输入回显与 `provenance`。

**示例代码：**

```python
from dnakit.databases import ensembl_dated_release_urls

result = ensembl_dated_release_urls(
    "GCA_000001405.29", "ensembl", "2026_04", genome_fasta_mask="hard", verbose=False
)
print(result["directory"])
print(result["gtf"])
print(result["genome_fasta"])
```

**示例输出：**

```text
https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04
https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/geneset/genes.gtf.gz
https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/genome/hardmasked.fa.bgz
```

## 2) `ensembl_dated_genome_fasta_mask`

- **作用：** 从 dated 版本的参考 DNA 地址读出屏蔽类型（none、soft 或 hard）。
- **处理方法：** 按地址末尾的 `unmasked.fa.bgz`、`softmasked.fa.bgz` 或 `hardmasked.fa.bgz` 判断；不匹配时 `genome_fasta_mask` 为 None，`matched` 为 False。
- **API：** `ensembl_dated_genome_fasta_mask(url[必须], *, verbose[可选])`
- **输入：** 地址字符串，至多 2048 个字符。
- **输出：** `url`、`genome_fasta_mask` 与 `matched`。

**示例代码：**

```python
from dnakit.databases import ensembl_dated_genome_fasta_mask

base = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/genome"
print(ensembl_dated_genome_fasta_mask(f"{base}/hardmasked.fa.bgz", verbose=False)["genome_fasta_mask"])
print(ensembl_dated_genome_fasta_mask(f"{base}/cdna.fa.bgz", verbose=False)["matched"])
```

**示例输出：**

```text
hard
False
```

## 3) `ensembl_genome_fasta_source_options`

- **作用：** 读出官方 Ensembl 参考 DNA 地址的 FASTA 类型（toplevel 或 primary_assembly）与屏蔽类型。
- **处理方法：** 地址须为 `ftp.ensembl.org` 或 `ftp.ensemblgenomes.ebi.ac.uk` 的 https 地址，且文件名符合 pyensembl 的官方格式；`dna_sm` 为 soft，`dna_rm` 为 hard，`dna` 为 none。其他地址返回 `canonical` 为 False，类型与屏蔽为 None。
- **API：** `ensembl_genome_fasta_source_options(source[必须], *, verbose[可选])`
- **输入：** 地址字符串，至多 2048 个字符。
- **输出：** `canonical`、`genome_fasta_type` 与 `genome_fasta_mask`。

**示例代码：**

```python
from dnakit.databases import ensembl_genome_fasta_source_options

url = (
    "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/"
    "Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz"
)
result = ensembl_genome_fasta_source_options(url, verbose=False)
print(result["canonical"], result["genome_fasta_type"], result["genome_fasta_mask"])
print(ensembl_genome_fasta_source_options("/data/genome.fa.gz", verbose=False)["canonical"])
```

**示例输出：**

```text
True toplevel soft
False
```

## 4) `ensembl_species_record`

- **作用：** 返回物种的拉丁名、同义名、各参考组装对应的 Ensembl 发布版本范围，以及 dated 注释的登录号与来源名。
- **处理方法：** `by="name"` 通过 pyensembl 物种注册表按拉丁名或常用名（如 `human`）查找；`by="reference"` 按参考组装名（如 `GRCh38`，大小写不敏感）查找。未知名称抛出 `ConfigurationError`。
- **API：** `ensembl_species_record(query[必须], *, by[可选], verbose[可选])`
- **输入：** 物种名或参考组装名；`by` 取 `name` 或 `reference`，默认 `name`。
- **输出：** `latin_name`、`synonyms`、`division`、`ensembl_genomes`、`reference_assemblies`（每项含 `reference`、`first_release`、`last_release`）与 `dated_releases`。

**示例代码：**

```python
from dnakit.databases import ensembl_species_record

record = ensembl_species_record("human", verbose=False)
print(record["latin_name"], record["division"])
for item in record["reference_assemblies"]:
    print(item["reference"], item["first_release"], item["last_release"])
print(ensembl_species_record("GRCh38", by="reference", verbose=False)["latin_name"])
```

**示例输出：**

```text
homo_sapiens vertebrates
GRCh37 55 75
GRCh38 76 116
NCBI36 54 54
homo_sapiens
```

## 5) `map_genome_contig_names`

- **作用：** 把本地已安装基因组的染色体或 contig 名映射到 UCSC、NCBI 或 Ensembl 的命名方式。
- **处理方法：** 读取基因组文件夹中的 `README.txt`（取来源 provider）与 `assembly_report.txt`（genomepy 安装时保存的报告）。两个文件缺失时直接报错，不调用上游下载。来源与目标相同，或报告中没有目标命名时，`available` 为 False。目标名为 `na` 的行返回 None。
- **API：** `map_genome_contig_names(genome_name[必须], target_provider[必须], *, genomes_dir[必须], verbose[可选])`
- **输入：** 基因组文件夹名（如 `demo`）、目标 provider（`ucsc`、`ncbi` 或 `ensembl`），以及包含该文件夹的 `genomes_dir`。下面的示例文件放在一个新建的 `genomes/demo` 目录中。
- **输出：** `genome`、`source_provider`、`target_provider`、`available`、`reason`、`record_count`，以及按来源名排序的 `records`（每项含 `source_name` 与 `target_name`）。

**示例文件：** `demo_README.txt`

```text
name: demo
provider: NCBI
assembly_accession: GCF_000001405.40
tax_id: 9606
species: Homo sapiens
genome_file: demo.fa
```

**示例文件：** `demo_assembly_report.txt`

```text
Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
1	assembled-molecule	1	Chromosome	CM000663.2	=	NC_000001.11	Primary Assembly	248956422	chr1
2	assembled-molecule	2	Chromosome	CM000664.2	=	NC_000002.12	Primary Assembly	242193529	chr2
MT	assembled-molecule	MT	Mitochondrion	J01415.2	=	NC_012920.1	non-nuclear	16569	chrM
KI270706.1	unplaced-scaffold	na	Unplaced-Scaffold	KI270706.1	=	na	Unlocalized	175055	na
```

**示例代码：**

```python
import os
import shutil

from dnakit.databases import map_genome_contig_names

os.makedirs("genomes/demo", exist_ok=True)
shutil.copy("demo_README.txt", "genomes/demo/README.txt")
shutil.copy("demo_assembly_report.txt", "genomes/demo/assembly_report.txt")
result = map_genome_contig_names("demo", "UCSC", genomes_dir="genomes", verbose=False)
print(result["source_provider"], result["target_provider"], result["record_count"])
for record in result["records"]:
    print(record["source_name"], record["target_name"])
```

**示例输出：**

```text
ncbi ucsc 4
1 chr1
2 chr2
KI270706.1 None
MT chrM
```

## 未包装的函数

- `needs-data` 4 个：`cached_release` 与 `genome_for_reference_name` 各有两个模块别名，返回依赖本地注释缓存的 EnsemblRelease 对象。
- `needs-network` 1 个：gget 的 `cellxgene`，需要 CELLxGENE 网络查询与未安装的 `cellxgene_census`。
- `not-dna` 9 个：物种名称规范化 2 个、Ensembl 注释日期列表解析 1 个、provider 列表 2 个、gget 的 cBioPortal 检索与 OpenAI 文本生成 2 个、biothings_client 的两个通用 join（连接两个网络查询客户端）2 个。
- `cli` 2 个：gget 的 `setup` 与 genomepy 的 `generate_env`，分别是安装器与 shell 导出文件写入器。
- `internal` 51 个：缓存路径与缓存检查、销毁性清理（`prune_genome_fastas`、`clean`）、发布版本校验函数、genomepy 的 README 读写与更新、客户端工厂、请求缓存与配置写入等辅助函数。
- `family-variant` 16 个：上面各包装函数的模块别名，见各包装节的台账说明。

## 测试

测试文件为 `others/tests/unit/test_review_G18j.py`，共 12 个用例，全部离线运行。每个包装函数至少有一个用例，将输出与手工计算的字符串、上游的直接调用或注册表中的已知值比较。genomepy 的缺失报告用例替换了上游的下载函数，并断言它不会被调用。
