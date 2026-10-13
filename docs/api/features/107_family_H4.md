# H4 家族变体复核：polars-bio、pyensembl、gtars、bioframe 等 9 个包

本页记录第 H4 组对台账中 `family-variant` 行的复核结果。每一行原先都被标为“已包装函数的变体”（别名、再导出，或由参数选择的分支）。这里逐行核对：别名与再导出改为 `duplicate-of`；确实能由已有包装函数的参数到达的保留为 `family-variant`；没有任何已有包装函数能到达的新增包装。新包装位于 `src/dnakit/annotation/gtars_family_functions.py`、`src/dnakit/annotation/bioframe_family_functions.py` 与 `src/dnakit/databases/pyensembl_family_functions.py`；这些函数已从各区域包导出，示例从区域包导入。

- 共复核 132 行，涉及 polars_bio、bioframe、gtars、pybedtools、gffutils、gtfparse、genomicranges（annotation 区）与 pyensembl、genomepy（databases 区）。
- 72 行改为 `duplicate-of`，48 行经参数验证后保留为 `family-variant`，8 行新增包装，2 行为 `internal`，2 行为 `broken-upstream`。
- 测试位于 `others/tests/unit/test_family_H4.py`。保留的 `family-variant` 行用已有包装函数加上所声明的参数调用，并与上游直接调用比较；`duplicate-of` 行用对象同一性或上游源码核对。
- 示例中的输出由实际运行得到，文档检查器应报告 0 个不一致。

## 复核结果

| 包 | 复核前 family-variant | duplicate-of | family-variant（已验证） | wrapped（新增） | 其他状态 |
| --- | --- | --- | --- | --- | --- |
| polars_bio | 52 | 22 | 28 | 0 | internal 2 |
| pyensembl | 42 | 38 | 1 | 3 | 无 |
| gffutils | 10 | 0 | 10 | 0 | 无 |
| genomepy | 9 | 5 | 4 | 0 | 无 |
| bioframe | 6 | 3 | 0 | 3 | 无 |
| gtars | 6 | 1 | 3 | 2 | 无 |
| pybedtools | 4 | 0 | 2 | 0 | broken-upstream 2 |
| gtfparse | 2 | 2 | 0 | 0 | 无 |
| genomicranges | 1 | 1 | 0 | 0 | 无 |

主要发现：

- pyensembl 的 38 行是再导出：`pyensembl.genome.normalize_chromosome` 等路径与定义模块里的函数是同一个对象（用 `is` 核对），不是变体。只有定义在 `pyensembl.reference_name` 中的 `find_species_by_reference` 保留为 `ensembl_species_record(by='reference')`。
- `make_gtf_filename` 与 `make_fasta_filename` 只返回文件名，现有 URL 包装函数没有参数能只返回文件名，因此新增包装。
- `is_canonical_source` 的声明被推翻：它对官方 `dna/` 目录下的任何文件名都为真，而 `ensembl_genome_fasta_source_options` 的 `canonical` 字段只对 toplevel 或 primary_assembly 的 DNA 文件名为真。例如 `.../dna/Homo_sapiens.GRCh38.cdna.all.fa.gz` 前者为 True，后者为 False。
- gtars 的 `digest_sequence` 与 `load_fasta` 返回检测到的字母表（`dna2bit`、`dna3bit`、`dnaio` 等），已有的 `sequence_digest` 与 `digest_fasta_sequences` 不返回，因此新增包装。
- bioframe 的 `select_mask`、`select_indices`、`select_labels` 返回行位置而不是区间记录，`select_bioframe_region` 不含这些信息，因此新增一个包装；`bioframe.ops` 下的三个别名与顶层函数是同一个对象。
- pybedtools 的 `random_intersection` 与 `random_intersection_bp` 不是变体，而且在当前环境中不可用：上游 `close_or_delete` 对已关闭的生成器抛出 `StopIteration`，每次调用都失败，标为 `broken-upstream`。
- polars_bio 的 `scan_*` 与 `sink_*` 不能由参数到达：`read_*` 的源码就是 `scan_*(...).collect()`，`write_*` 与 `sink_*` 调用同一个写入函数，改为 `duplicate-of`。`print_metadata_json` 与 `print_metadata_summary` 只向控制台打印 `get_metadata` 的内容并返回 None，标为 `internal`。
- gtfparse 的 `parse_gtf_pandas` 并不是“返回 pandas 的变体”：`parse_gtf` 本来就返回 pandas，前者只是加了进度回调，改为 `duplicate-of`。genomicranges 的 `parse_gtf` 是 `read_gtf` 的第一步，同样改为 `duplicate-of`。
- gffutils 的 10 个 `merge_criteria` 函数都能由 `merge_gffutils_features` 的 `criteria` 与 `threshold` 参数到达，测试逐个与 `db.merge` 的直接调用比较。
- genomepy 的 `region_to_seq`、`regions_to_seqs`、`bed_to_seq`、`bed_to_seqs` 由 `Genome.track2fasta` 调用，`write_bed_sequences` 的 `bed` 参数传入区域文件或 BED 文件即可到达。
- polars_bio 的 `write_bam` 与 `write_cram` 保留为 `write_alignment_file(output_format=...)`；上游的 `sort_on_write` 与 `tag_type_overrides` 参数没有被暴露。

## 通用约定

- 每个公开函数的 `verbose` 都是仅关键字参数，默认 `True`。
- 返回值是可 JSON 序列化的字典，包含 `status`、`operation`、结果字段与 `provenance`。
- `pyensembl`、`gtars` 与 `bioframe` 只在函数运行时导入；缺少依赖时抛出 `BackendUnavailableError`。
- 输入不合法时抛出 `ConfigurationError`。
- 三个 pyensembl 函数只处理字符串，不访问网络。

## 1) `ensembl_gtf_filename`

- **作用：** 返回某物种与 Ensembl 版本的基因注释 GTF 文件名。
- **处理方法：** 调用上游 `ensembl_url_templates.make_gtf_filename`。文件名与 `ensembl_gtf_url` 返回的 URL 末段相同，但 URL 包装函数没有参数只返回文件名。
- **API：** `ensembl_gtf_filename(ensembl_release[必须], species[必须], *, verbose[可选])`
- **输入：** `ensembl_release` 是 1 到 1000 的整数；`species` 是物种名，如 `homo_sapiens`。
- **输出：** 字典，含 `species`、`ensembl_release`、`filename` 与 `provenance`。

**示例代码：**

```python
from dnakit.databases import ensembl_gtf_filename

print(ensembl_gtf_filename(75, "homo_sapiens", verbose=False)["filename"])
print(ensembl_gtf_filename(109, "homo_sapiens", verbose=False)["filename"])
```

**示例输出：**

```text
Homo_sapiens.GRCh37.75.gtf.gz
Homo_sapiens.GRCh38.109.chr_patch_hapl_scaff.gtf.gz
```

## 2) `ensembl_sequence_fasta_filename`

- **作用：** 返回某物种与版本的 cDNA、ncRNA 或蛋白 FASTA 文件名。
- **处理方法：** 调用上游 `ensembl_url_templates.make_fasta_filename`；新版目录布局（版本大于 75）与旧版布局的文件名不同，由上游决定。
- **API：** `ensembl_sequence_fasta_filename(ensembl_release[必须], species[必须], sequence_type[必须], *, verbose[可选])`
- **输入：** `sequence_type` 是 `cdna`、`ncrna` 或 `pep`，其余同上。
- **输出：** 字典，含 `species`、`ensembl_release`、`sequence_type`、`filename` 与 `provenance`。

**示例代码：**

```python
from dnakit.databases import ensembl_sequence_fasta_filename

for kind in ("cdna", "ncrna", "pep"):
    print(kind, ensembl_sequence_fasta_filename(109, "homo_sapiens", kind, verbose=False)["filename"])
```

**示例输出：**

```text
cdna Homo_sapiens.GRCh38.cdna.all.fa.gz
ncrna Homo_sapiens.GRCh38.ncrna.fa.gz
pep Homo_sapiens.GRCh38.pep.all.fa.gz
```

## 3) `check_ensembl_canonical_source`

- **作用：** 判断一个 URL 是否位于官方 Ensembl 参考 FASTA 的 `.../fasta/<物种>/dna/<文件>` 路径。
- **处理方法：** 调用上游 `genome_fasta_cache.is_canonical_source`，只检查协议（https）、主机与路径布局，不检查文件名。要读取 DNA 类型与掩码，用 `ensembl_genome_fasta_source_options`。
- **API：** `check_ensembl_canonical_source(source[必须], *, verbose[可选])`
- **输入：** `source` 是 URL 字符串。
- **输出：** 字典，含 `source`、`official_source`（布尔）与 `provenance`。

**示例代码：**

```python
from dnakit.databases import ensembl_genome_fasta_source_options
from dnakit.databases import check_ensembl_canonical_source

base = "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/"
for name in ("Homo_sapiens.GRCh38.dna.toplevel.fa.gz", "Homo_sapiens.GRCh38.cdna.all.fa.gz"):
    url = base + name
    official = check_ensembl_canonical_source(url, verbose=False)["official_source"]
    dna = ensembl_genome_fasta_source_options(url, verbose=False)["canonical"]
    print(name, official, dna)
```

**示例输出：**

```text
Homo_sapiens.GRCh38.dna.toplevel.fa.gz True True
Homo_sapiens.GRCh38.cdna.all.fa.gz True False
```

## 4) `sequence_digest_record`

- **作用：** 返回一条序列的 GA4GH sha512t24u 与 MD5 摘要、长度和 gtars 检测到的字母表。
- **处理方法：** 序列先转大写并检查 DNA 符号，再调用上游 `gtars.refget.digest_sequence`。`name` 默认为空字符串。
- **API：** `sequence_digest_record(sequence[必须], *, name[可选], description[可选], verbose[可选])`
- **输入：** `sequence` 是非空 DNA 字符串（允许 IUPAC 符号）；`name` 与 `description` 是可选文本。
- **输出：** 字典，含 `name`、`description`、`length`、`alphabet`、`sha512t24u`、`md5` 与 `provenance`。

**示例代码：**

```python
from dnakit.annotation import sequence_digest_record

result = sequence_digest_record("acgtacgtnn", name="s1", verbose=False)
print(result["length"], result["alphabet"], result["sha512t24u"], result["md5"])
print(sequence_digest_record("ACGT", verbose=False)["alphabet"])
```

**示例输出：**

```text
10 dna3bit 5rfV7JOOAEVoyCNHpPwNcJJwGyVAjmH5 a0e5e77e1f26e1b3d8a27baa3db84c40
dna2bit
```

## 5) `load_fasta_sequence_collection`

- **作用：** 用 gtars 载入 FASTA 文件，返回每条序列的摘要、字母表与（可选的）残基。
- **处理方法：** 调用上游 `gtars.refget.load_fasta`。集合摘要与 `digest_fasta_sequences` 的一致；结果按序列名排序。总碱基数超过 `max_bases` 时抛出 ConfigurationError。
- **API：** `load_fasta_sequence_collection(fasta[必须], *, include_sequences[可选], max_rows[可选], max_bases[可选], verbose[可选])`
- **输入：** `fasta` 是本地 FASTA 路径；`include_sequences` 默认 `True`；`max_rows` 默认 100000；`max_bases` 默认 100000000。
- **输出：** 字典，含 `sequence_count`、`total_bases`、`collection_digest`、`sequences` 与 `provenance`。每条序列含 `name`、`description`、`length`、`alphabet`、`sha512t24u`、`md5` 与 `sequence`（大写，无换行）。

**示例文件：** `t.fa`

```text
>a first record
ACGTAC
GGTT
>b
NNACGU
```

**示例代码：**

```python
from dnakit.annotation import load_fasta_sequence_collection

result = load_fasta_sequence_collection("t.fa", verbose=False)
print(result["sequence_count"], result["total_bases"], result["collection_digest"])
for row in result["sequences"]:
    print(row["name"], row["length"], row["alphabet"], row["sequence"])
```

**示例输出：**

```text
2 16 KvfiMp1LuIdD2NYCZudr72JOMy2jnJNr
a 10 dna2bit ACGTACGGTT
b 6 dnaio NNACGU
```

## 6) `select_bioframe_overlap_rows`

- **作用：** 返回 BED 文件中哪些行与一个区域重叠，给出布尔掩码、行位置与索引标签。
- **处理方法：** 读取 BED3 到 BED6 文件，依次调用 bioframe 的 `select_mask`、`select_indices` 与 `select_labels`。区间为 0 起点、半开；只写染色体名时选中该染色体上的全部区间。BED 文件读入时使用默认整数索引，因此 `labels` 等于 `indices`。
- **API：** `select_bioframe_overlap_rows(source[必须], region[必须], *, max_rows[可选], verbose[可选])`
- **输入：** `source` 是本地 BED 路径；`region` 是 UCSC 风格区域字符串；`max_rows` 默认 100000。
- **输出：** 字典，含 `region`、`row_count`、`selected_count`、`mask`（每个数据行一个布尔值，按文件顺序）、`indices`、`labels` 与 `provenance`。

**示例文件：** `a.bed`

```text
chr1	10	30
chr1	50	80
chr2	5	15
chr1	100	120
```

**示例代码：**

```python
from dnakit.annotation import select_bioframe_overlap_rows

result = select_bioframe_overlap_rows("a.bed", "chr1:25-55", verbose=False)
print(result["mask"], result["indices"], result["labels"])
print(select_bioframe_overlap_rows("a.bed", "chr1", verbose=False)["indices"])
```

**示例输出：**

```text
[True, True, False, False] [0, 1] [0, 1]
[0, 1, 3]
```
