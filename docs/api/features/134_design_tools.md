# 设计工具：pydna、DNA Chisel 与 DnaCauldron

这一组 24 个函数位于 `dnakit.molecular_lab`，用来补全最初那批适配器在三个构建体设计包上的缺口：[pydna](https://github.com/pydna-group/pydna)（凝胶迁移预测与绘图、密码子适应指数、SnapGene 文件、NCBI BLAST 与 NEB Tm 服务）、[DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel)（序列组成、模式注释、约束违反表与图、优化报告、记录读取、BLAST 与 Bowtie 检索）和 [DnaCauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron)（装配方案模拟、分层方案图与 PDF 摘要）。导入方式为 `from dnakit.molecular_lab import ...`。这些函数除上游例程本身计算的内容外，不设计也不预测任何湿实验结果。

上游包都是可选依赖，只有调用对应函数时才会导入；缺少包时会抛出带安装命令的 `BackendUnavailableError`。根据已安装包的元数据，许可证为：pydna 为 BSD，DnaChisel 为 MIT，DnaCauldron 为 MIT，不涉及 GPL 包。有三个函数通过 DNAKit 的受限命令运行器调用外部程序（[BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) 的 `blastn`，以及 [Bowtie 1](https://github.com/BenLangmead/bowtie) 的 `bowtie` 与 `bowtie-build`），两个函数封装网络服务（NCBI BLAST 和 NEB Tm 计算器）。报告与绘图类函数还可能需要 `pdf_reports` 与 `sequenticon`、`geneblocks` 或 WeasyPrint，具体见各节。所有写文件的函数都不会覆盖已存在的路径。

## 示例数据

有两个示例读取声明的输入文件：一个单条记录的 FASTA 文件，其中 `GGTCTC` 位点位于第 3 至 9 位（0 起始）；以及一个把约束写在特征标签里的 GenBank 记录。其余示例直接传入输入，或把文件写入临时目录。

**示例文件：** `part.fa`

```text
>part1
ATGGGTCTCATGCCCAAAGCTTTAA
```

**示例文件：** `annotated.gb`

```text
LOCUS       annotated                 25 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  .
FEATURES             Location/Qualifiers
     misc_feature    4..9
                     /label="@AvoidPattern(GGTCTC)"
ORIGIN
        1 atgggtctca tgcccaaagc tttaa
//
```

## 1) `dnakit.molecular_lab.pydna_blast_best_nt_hit()`

- **作用：** 用 BLAST 为一条 DNA 序列查找 NCBI nt 中的最佳命中，并返回匹配区域。
- **处理方法：** 先校验序列，再调用 pydna 的 `sequence_picker.genbank_accession`：它向 NCBI 的 nt 数据库提交 megablast 检索，并把最佳命中区域作为记录返回。尽管名字如此，它并不查询登录号。`name` 只是原样回传为 `query_name`，不会发送给 BLAST。
- **API：** `pydna_blast_best_nt_hit(sequence: 'str', *, name: 'str' = 'query', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pydna（extra `dnakit[bio-assembly]`）和访问 NCBI BLAST 服务的网络（无需密钥）。离线无法运行，所以示例展示的是在任何网络调用之前就会执行的输入检查。装好依赖后的用法见示例之后。
- **输入：** `sequence` 是由 A、C、G、T 组成的字符串（忽略空白和大小写），最长 20000 个碱基。`name` 是长度不超过 128 的非空字符串。
- **输出：** `service`、`source`（BLAST 网址）和 `result`，其中含 `query_name`、`sequence`（匹配区域，大写）、`length`、`circular` 和 `description`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import pydna_blast_best_nt_hit

try:
    pydna_blast_best_nt_hit("ATGCNNGTAC", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "only A, C, G and T" in str(error))
```

- **示例输出：**

```text
ConfigurationError True
```

安装依赖后的用法：

```python
from dnakit.molecular_lab import pydna_blast_best_nt_hit

hit = pydna_blast_best_nt_hit("ATGGCTAAAGCT" * 6, name="insert", verbose=False)
print(hit["result"]["length"], hit["result"]["description"])
```

## 2) `dnakit.molecular_lab.pydna_codon_adaptation_index()`

- **作用：** 计算一条编码序列的密码子适应指数（CAI）。
- **处理方法：** 给出显式 `weights` 时，指数在 DNAKit 内部计算：对被计分密码子的权重取几何平均（Sharp 和 Li，1987），并像 pydna 一样四舍五入到 3 位小数。ATG、TGG 以及终止密码子 TAA、TAG、TGA 不计分。不给 `weights` 时，调用 pydna 的 `utils.cai` 并使用物种密码子表，这需要 cai2 包。
- **API：** `pydna_codon_adaptation_index(seq: 'str', *, weights: 'dict[str, float] | None' = None, organism: 'str' = 'sce', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 给出 `weights` 时安装基础包即可运行（不导入任何可选包）。不给 `weights` 时需要 pydna（extra `dnakit[bio-assembly]`）和 cai2（`pip install "pydna[express]"`）。
- **输入：** `seq` 是由 A、C、G、T 组成的字符串，长度为 3 的整数倍，最长 300000 个碱基。`weights` 把三字母密码子映射到正的有限数；`seq` 中每个被计分的密码子都必须有权重，否则 `ConfigurationError` 会列出缺失的密码子。文档字符串说权重位于 (0, 1]，但代码只要求为正数。`organism`（默认 `sce`）必须是 pydna 密码子表中的键，且仅在不给 `weights` 时使用。
- **输出：** `status`、`operation`、`method`、`cai`、`codon_count`、`scored_codon_count`（仅显式权重路径）和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import pydna_codon_adaptation_index

weights = {"GCT": 1.0, "GCC": 0.5, "AAA": 0.8, "AAG": 0.4}
result = pydna_codon_adaptation_index("ATGGCTGCCAAAAAGTAA", weights=weights, verbose=False)
print(result["cai"], result["codon_count"], result["scored_codon_count"])
```

- **示例输出：**

```text
0.632 6 4
```

## 3) `dnakit.molecular_lab.pydna_draw_gel()`

- **作用：** 根据片段长度绘制琼脂糖凝胶图，并写成一个新的 PNG。
- **处理方法：** 建立 pydna 的凝胶尺寸模型（对所选 ladder 的迁移率做三次样条），把每个片段长度变成 `FakeSeq`，调用 `pydna.gel.gel`，再把 Pillow 图像保存为 PNG。每个条带长度都必须落在 ladder 的范围内。
- **API：** `pydna_draw_gel(output_path: 'str', *, lanes: 'list[list[int]] | None' = None, ladder: 'str' = 'GeneRuler_1kb_plus', gel_length: 'int' = 600, margin: 'int' = 50, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pydna（extra `dnakit[bio-assembly]`）及其 gel extra，后者提供 Pillow 和 SciPy（`pip install "pydna[gel]"`）。
- **输入：** `output_path` 以 `.png` 结尾且不存在。`lanes` 是 1 到 20 条泳道的列表，每条泳道是 1 到 50 个片段长度（bp）的列表；为 `None` 时把 ladder 本身画成一条泳道。`ladder` 为 `GeneRuler_1kb_plus` 或 `HI_LO_DNA_MARKER`。`gel_length` 为 100 到 5000，`margin` 为 0 到 `gel_length - 10`。
- **输出：** `path` 和 `data`，其中含 `ladder`、`lanes`、`lane_count`、`band_count`、`gel_length`、`margin`、`image_width`、`image_height` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import pydna_draw_gel

with tempfile.TemporaryDirectory() as folder:
    result = pydna_draw_gel(
        str(Path(folder) / "gel.png"),
        lanes=[[500, 1000], [2000, 3000, 5000]],
        verbose=False,
    )
    data = result["data"]
    print(data["lane_count"], data["band_count"], data["image_width"], data["image_height"])
    print(Path(result["path"]).read_bytes()[:4])
```

- **示例输出：**

```text
2 5 180 600
b'\x89PNG'
```

## 4) `dnakit.molecular_lab.pydna_neb_melting_temperature()`

- **作用：** 用 NEB Tm 计算器返回引物的解链温度（摄氏度）。
- **处理方法：** 先校验引物、浓度和产品代码，再调用 pydna 的 `tm.tm_neb`，它把引物发送到 NEB 服务 `https://tmapi.neb.com/tm`。
- **API：** `pydna_neb_melting_temperature(primer: 'str', *, conc: 'float' = 0.5, prodcode: 'str' = 'q5-0', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pydna（extra `dnakit[bio-assembly]`）和访问 NEB Tm API 的网络（无需密钥）。离线无法运行，所以示例展示的是在任何网络调用之前就会执行的输入检查。装好依赖后的用法见示例之后。
- **输入：** `primer` 是由 A、C、G、T 组成的字符串，最长 200 个碱基。`conc` 是引物浓度，单位微摩尔，范围 0.001 到 100。`prodcode` 是 NEB 产品代码，由 1 到 64 个字母、数字、`-` 或 `_` 组成（例如 Q5 聚合酶用 `q5-0`）。
- **输出：** `service`、`source` 和 `result`，其中含 `primer`、`conc_uM`、`prodcode` 和 `tm_celsius`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import pydna_neb_melting_temperature

try:
    pydna_neb_melting_temperature("ATGCGTACGTTAGC", prodcode="q5 0", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "prodcode must be" in str(error))
```

- **示例输出：**

```text
ConfigurationError True
```

安装依赖后的用法：

```python
from dnakit.molecular_lab import pydna_neb_melting_temperature

result = pydna_neb_melting_temperature("ATGCGTACGTTAGC", conc=0.5, prodcode="q5-0", verbose=False)
print(result["result"]["tm_celsius"])
```

## 5) `dnakit.molecular_lab.pydna_parse_snapgene_history()`

- **作用：** 用 pydna 解析本地的 SnapGene `.dna` 文件，文件带克隆历史时一并解析出来。
- **处理方法：** 对文件调用 pydna 的 `snapgene_history_parser.parse_snapgene_history`，汇总根记录及其特征，并用 `sgffp.SgffReader` 判断文件是否带有历史。没有历史的文件只返回根序列。
- **API：** `pydna_parse_snapgene_history(path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pydna（extra `dnakit[bio-assembly]`），它依赖 sgffp。
- **输入：** `path` 是以 `.dna` 结尾的本地路径；网址会被拒绝。
- **输出：** `name`、`length`、`circular`、`sequence`（大写）、`feature_count`、`features`（每项含 `type`、`start`、`end`、`strand`、`label`）、`has_history`、`source_type`（记录来源对象的类名，或 `None`）和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import sgffp

from dnakit.molecular_lab import pydna_parse_snapgene_history

sequence = "ATGCGTACGTTAGCATGCGTACGTTAGC"
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "construct.dna"
    obj = sgffp.SgffObject.new(sequence=sequence, topology="circular")
    path.write_bytes(sgffp.SgffWriter.to_bytes(obj))
    result = pydna_parse_snapgene_history(str(path), verbose=False)
print(result["name"], result["length"], result["circular"])
print(result["has_history"], result["feature_count"], result["sequence"] == sequence)
```

- **示例输出：**

```text
construct 28 True
False 0 True
```

## 6) `dnakit.molecular_lab.pydna_predict_gel_relative_front()`

- **作用：** 预测 DNA 片段在琼脂糖凝胶上的相对迁移前沿（Rf）。
- **处理方法：** 建立 pydna 的凝胶尺寸模型，即对内置 ladder 的 Rf 值做自然三次样条（`pydna.gel.interpolator`），并在每个片段长度处取值。落在 ladder 范围之外的长度没有预测值（`None`）。不会绘图；要出图请用 `pydna_draw_gel`。
- **API：** `pydna_predict_gel_relative_front(lengths: 'list[int]', *, ladder: 'str' = 'GeneRuler_1kb_plus', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pydna（extra `dnakit[bio-assembly]`）以及 SciPy（`pip install "pydna[gel]"`）。
- **输入：** `lengths` 是 1 到 500 个整数的列表，取值 1 到 10000000（bp）。`ladder` 为 `GeneRuler_1kb_plus` 或 `HI_LO_DNA_MARKER`。
- **输出：** `ladder`、`ladder_sizes`、`ladder_relative_fronts`、`lengths`、`relative_front`（每个长度对应一个浮点数或 `None`）和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import pydna_predict_gel_relative_front

result = pydna_predict_gel_relative_front([500, 1000, 3000, 50000], verbose=False)
print(result["ladder"], len(result["ladder_sizes"]))
print([None if rf is None else round(rf, 3) for rf in result["relative_front"]])
```

- **示例输出：**

```text
GeneRuler_1kb_plus 15
[0.755, 0.562, 0.261, None]
```

## 7) `dnakit.molecular_lab.annotate_pattern_occurrences_features()`

- **作用：** 在两条链上查找 IUPAC DNA 模式的每一次出现并加以注释。
- **处理方法：** 把模式中每个 IUPAC 符号展开成它代表的全部碱基，再用前瞻式正则扫描序列，因此会保留相互重叠的命中。除非模式是回文，否则还会扫描模式的反向互补以得到负链（strand -1）的命中；回文模式只搜索正链，与 DnaChisel 一致。可选地用 Biopython 把特征写入新的 GenBank 文件。函数说明 DnaChisel 自带的模式正则把 N、R、Y 映射到了错误的碱基集合，所以含这些符号的模式，其注释与 DnaChisel 不同。
- **API：** `annotate_pattern_occurrences_features(sequence: 'str', pattern: 'str', *, feature_type: 'str' = 'misc_feature', prefix: 'str' = '!', output_path: 'str | Path | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 安装基础包即可运行。写出 `output_path` 需要 Biopython（extra `dnakit[bio-formats]`）。
- **输入：** `sequence` 和 `pattern` 是 IUPAC DNA 字符串（序列可以为空；模式最长 1000 个碱基）。`feature_type` 是 GenBank 特征键；`prefix`（最长 20 个字符）加在模式前面构成标签。`output_path` 必须不存在，并以 `.gb`、`.gbk` 或 `.genbank` 结尾。
- **输出：** `sequence_length`、`pattern`、`palindromic`、`occurrence_count`、`features`（每项含 `type`、含起不含止的 `start` 与 `end`、取值 +1 或 -1 的 `strand`、`label`）、`output_path` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import annotate_pattern_occurrences_features

result = annotate_pattern_occurrences_features("AAGGTCTCTTGAGACCAA", "GGTCTC", verbose=False)
print(result["palindromic"], result["occurrence_count"])
for feature in result["features"]:
    print(feature["start"], feature["end"], feature["strand"], feature["label"])
```

- **示例输出：**

```text
False 2
2 8 1 !GGTCTC
10 16 -1 !GGTCTC
```

## 8) `dnakit.molecular_lab.blast_sequence_hits()`

- **作用：** 用 `blastn` 把一条序列与 BLAST 数据库或一个目标 FASTA 文件做比对。
- **处理方法：** 把查询序列写入临时 FASTA 文件，通过受限运行器（不经过 shell）以 DnaChisel 的默认值（字长 4、一致性 80 %、最多 1000 条比对、关闭 dust）运行 `blastn -outfmt 5`，再用 Biopython 的 `Bio.Blast.NCBIXML` 把 XML 解析成比对和高分片段对。先用 `-version` 读取 `blastn` 的版本。
- **API：** `blast_sequence_hits(sequence: 'str', *, blast_db: 'str | None' = None, subject_fasta: 'str | Path | None' = None, word_size: 'int' = 4, perc_identity: 'float' = 80.0, num_alignments: 'int' = 1000, ungapped: 'bool' = False, num_threads: 'int' = 3, culling_limit: 'int | None' = None, e_value: 'float | None' = None, task: 'str | None' = None, dust: 'str' = 'no', blastn_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 Biopython（extra `dnakit[bio-formats]`）和外部程序 `blastn`（NCBI BLAST+），它须在 PATH 中，或通过 `blastn_path` 传入。编写本页所用的环境没有安装 `blastn`，所以示例把 `blastn_path` 指向一个不存在的文件并展示相应的错误；这样在任何机器上输出都相同。装好程序后的用法见示例之后。
- **输入：** `sequence` 是 IUPAC DNA 字符串。`blast_db`（用 `makeblastdb` 建立的数据库前缀，不能以 `-` 开头）与 `subject_fasta`（已存在的文件）必须恰好给出一个。`word_size` 为 4 到 28，`perc_identity` 为 0 到 100，`num_alignments` 为 1 到 10000，`num_threads` 为 1 到 64，`task` 取 `blastn`、`blastn-short`、`dc-megablast`、`megablast` 之一，`dust` 取 `no`、`yes` 或形如 `20 64 1` 的三个整数。
- **输出：** `status`、`operation`、`backend`、`version`、`blast_db`、`subject_fasta`、`queries`（每项含 `query_id`、`query_length` 和 `alignments`；每个比对含 `title`、`length` 和 `hsps`，后者含 `identities`、`align_length`、`gaps`、`query_start`、`query_end`、`subject_start`、`subject_end`、`score`、`evalue`）、`hit_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import blast_sequence_hits

try:
    blast_sequence_hits(
        "ATGGGTCTCATGCCCAAA",
        blast_db="nt_database",
        blastn_path="/no/such/dir/blastn",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "missing or not executable" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

安装依赖后的用法：

```python
from dnakit.molecular_lab import blast_sequence_hits

result = blast_sequence_hits("ATGGGTCTCATGCCCAAA", subject_fasta="genes.fa", verbose=False)
print(result["hit_count"])
```

## 9) `dnakit.molecular_lab.constraints_breaches_table()`

- **作用：** 返回多条序列的约束违反表，每条序列一行。
- **处理方法：** 按 `evaluate_dnachisel_specs` 使用的允许列表构建 DnaChisel 约束，并调用 DnaChisel 的 `constraints_breaches_dataframe`。每一列对应一个约束，列名采用 DnaChisel 的违反标签形式；每个单元格以 `start-end(strand)` 的形式列出违反的位置，用逗号分隔，约束通过时为空。
- **API：** `constraints_breaches_table(sequences: 'object', *, constraints: 'Sequence[Mapping[str, Any]]', display_constraints_locations: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaChisel（extra `dnakit[bio-dnachisel]`）。
- **输入：** `sequences` 是映射 `{name: sequence}`，或 `{"name", "sequence"}` 映射的列表（1 到 1000 个互不重复的名字，IUPAC DNA，总计最多 2000000 个碱基）。`constraints` 是非空的映射列表，其 `name` 取自允许列表（`AvoidPattern`、`AvoidChanges`、`AvoidHairpins`、`AvoidStopCodons`、`CodonOptimize`、`EnforceChanges`、`EnforceGCContent`、`EnforceMeltingTemperature`、`EnforceTranslation`、`MaximizeCAI`），其余键是关键字参数。`display_constraints_locations` 会传给 DnaChisel，并改变列标签。
- **输出：** `sequence_count`、`constraint_labels`、`rows`（每条序列一个字典，含序列名 `sequence` 以及每个标签对应的键）和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import constraints_breaches_table

constraints = [
    {"name": "AvoidPattern", "pattern": "GGTCTC"},
    {"name": "EnforceGCContent", "mini": 0.3, "maxi": 0.7, "window": 18},
]
result = constraints_breaches_table(
    {"a": "ATGGGTCTCATGCCCAAA", "b": "ATGGAGACCAAGCCCAAA", "c": "A" * 22},
    constraints=constraints,
    verbose=False,
)
print(result["constraint_labels"])
for row in result["rows"]:
    print(row)
```

- **示例输出：**

```text
['GGTCTC', 'GC outside 30-70%/18bp']
{'sequence': 'a', 'GGTCTC': '3-9(+)', 'GC outside 30-70%/18bp': ''}
{'sequence': 'b', 'GGTCTC': '3-9(-)', 'GC outside 30-70%/18bp': ''}
{'sequence': 'c', 'GGTCTC': '', 'GC outside 30-70%/18bp': '0-22'}
```

## 10) `dnakit.molecular_lab.create_bowtie_index()`

- **作用：** 用若干序列建立 Bowtie 1 索引；序列以其位置（0、1、...）命名。
- **处理方法：** 创建新目录 `index_dir`，把序列写入其中的 `sequences.fa`，再通过受限运行器（不经过 shell）运行 `bowtie-build -f sequences.fa bowtie --quiet`。索引文件为 `bowtie.*.ebwt`；返回的 `index_prefix` 就是 `find_bowtie_matches` 所需的参数。
- **API：** `create_bowtie_index(sequences: 'Sequence[str]', *, index_dir: 'str | Path', bowtie_build_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `bowtie-build`（Bowtie 1），它须在 PATH 中，或通过 `bowtie_build_path` 传入。编写本页所用的环境没有安装它，所以示例把 `bowtie_build_path` 指向一个不存在的文件并展示相应的错误。装好程序后的用法见示例之后。
- **输入：** `sequences` 是非空的 IUPAC DNA 字符串列表，最多 1000 条（总计最多 10000000 个碱基）。`index_dir` 必须尚不存在。
- **输出：** `status`、`operation`、`backend`、`version`、`index_prefix`、`index_files`（相对文件名）、`sequence_count` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import create_bowtie_index

with tempfile.TemporaryDirectory() as folder:
    try:
        create_bowtie_index(
            ["ATGGGTCTCATGCCCAAAGCTTTAA"],
            index_dir=Path(folder) / "index",
            bowtie_build_path="/no/such/dir/bowtie-build",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, "missing or not executable" in str(error))
    print((Path(folder) / "index").exists())
```

- **示例输出：**

```text
BackendUnavailableError True
False
```

安装依赖后的用法：

```python
from dnakit.molecular_lab import create_bowtie_index

index = create_bowtie_index(["ATGGGTCTCATGCCCAAAGCTTTAA"], index_dir="bowtie_index", verbose=False)
print(index["index_prefix"], index["sequence_count"])
```

## 11) `dnakit.molecular_lab.find_bowtie_matches()`

- **作用：** 对序列的每个 k-mer，在索引中查找最佳的 Bowtie 1 命中（DnaChisel 的匹配器）。
- **处理方法：** 把序列切成长度为 `match_length` 的 k-mer，写入临时 FASTA 文件，并通过受限运行器运行 `bowtie --best -k 1 -v <max_mismatches> <index> -f ...`。每条报告的比对变成 `{"start", "end", "mismatches"}`，其中 `start` 与 `end` 是 k-mer 在 `sequence` 中的位置，`mismatches` 统计 Bowtie 描述符中的错配数。函数先检查 `<index_prefix>.1.ebwt` 是否存在，再去查找程序。
- **API：** `find_bowtie_matches(sequence: 'str', *, index_prefix: 'str | Path', match_length: 'int', max_mismatches: 'int' = 0, bowtie_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `bowtie`（Bowtie 1），它须在 PATH 中，或通过 `bowtie_path` 传入，并需要由 `create_bowtie_index` 建立的索引。编写本页所用的环境没有安装它，所以示例使用一个占位索引文件和一个不存在的 `bowtie_path`，并展示相应的错误。装好程序后的用法见示例之后。
- **输入：** `sequence` 是 IUPAC DNA 字符串。`index_prefix` 是已存在索引的前缀。`match_length` 为 1 到 `sequence` 的长度。`max_mismatches` 为 0 到 3。
- **输出：** `status`、`operation`、`backend`、`version`、`kmer_length`、`match_count`、`matches`（按 `start` 排序）和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import find_bowtie_matches

with tempfile.TemporaryDirectory() as folder:
    prefix = Path(folder) / "bowtie"
    Path(f"{prefix}.1.ebwt").write_bytes(b"")  # placeholder; only its existence is checked first
    try:
        find_bowtie_matches(
            "ATGGGTCTCATGCCCAAA",
            index_prefix=prefix,
            match_length=8,
            bowtie_path="/no/such/dir/bowtie",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, "missing or not executable" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

安装依赖后的用法：

```python
from dnakit.molecular_lab import create_bowtie_index, find_bowtie_matches

index = create_bowtie_index(["ATGGGTCTCATGCCCAAAGCTTTAA"], index_dir="bowtie_index", verbose=False)
found = find_bowtie_matches(
    "ATGGGTCTCATGCCC", index_prefix=index["index_prefix"], match_length=8, verbose=False
)
print(found["match_count"], found["matches"][:2])
```

## 12) `dnakit.molecular_lab.gc_fraction_or_windows()`

- **作用：** 返回整体 GC 比例，或每个滑动窗口的 GC 比例。
- **处理方法：** 先把序列转为大写，再调用 DnaChisel 的 `gc_content`。整体值是 G 与 C 的个数除以序列长度，包括 N 在内的每个符号都计入分母。给出 `window_size` 时，每个窗口起点 `i`（`0 <= i <= len - window_size`）得到一个值，覆盖 `[i, i + window_size)`。
- **API：** `gc_fraction_or_windows(sequence: 'str', *, window_size: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaChisel（extra `dnakit[bio-dnachisel]`）。
- **输入：** `sequence` 是 IUPAC DNA 字符串；只有给出 `window_size` 时才允许为空。`window_size` 是 1 到 200000 的整数；窗口比序列长时没有窗口。
- **输出：** `sequence_length`、`window_size`、`gc_fraction`（窗口模式下为 `None`）、`window_count`、`window_gc_fractions` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import gc_fraction_or_windows

whole = gc_fraction_or_windows("ATGCGCNNAT", verbose=False)
print(whole["sequence_length"], whole["gc_fraction"])
windows = gc_fraction_or_windows("ATGCGCAT", window_size=4, verbose=False)
print(windows["window_count"], windows["window_gc_fractions"])
```

- **示例输出：**

```text
10 0.4
5 [0.5, 0.75, 1.0, 0.75, 0.5]
```

## 13) `dnakit.molecular_lab.generate_random_dna_with_composition()`

- **作用：** 按给定的 GC 比例或各碱基概率生成随机 DNA 序列。
- **处理方法：** 调用 DnaChisel 的 `random_dna_sequence`。`gc_share` 是 G+C 的总概率，在 G 与 C 之间、A 与 T 之间平均分配。`probas` 把碱基映射到概率，且总和必须为 1。DnaChisel 使用 numpy 的全局随机数生成器，所以 `seed` 会给这个生成器设种子（这是原始例程的副作用）。
- **API：** `generate_random_dna_with_composition(length: 'int', *, gc_share: 'float | None' = None, probas: 'Mapping[str, float] | None' = None, seed: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaChisel（extra `dnakit[bio-dnachisel]`）。
- **输入：** `length` 为 0 到 200000。`gc_share` 为 0 到 1。`probas` 的键是 A、C、G、T 中互不相同的碱基。同时给出 `gc_share` 和 `probas` 会报错；两者都不给时四种碱基等概率。`seed` 为 0 到 4294967295。
- **输出：** `sequence`、`length`、`composition_source`（`gc_share`、`probas` 或 `uniform`）、`gc_share`、`probas`、`gc_fraction_observed`（长度为 0 时为 `None`）、`seed` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import generate_random_dna_with_composition

first = generate_random_dna_with_composition(20, gc_share=0.5, seed=1, verbose=False)
second = generate_random_dna_with_composition(20, gc_share=0.5, seed=1, verbose=False)
print(first["sequence"], first["composition_source"], first["gc_fraction_observed"])
print(first["sequence"] == second["sequence"])
```

- **示例输出：**

```text
CAGCGGGCCACAGTGACAGG gc_share 0.7
True
```

## 14) `dnakit.molecular_lab.list_common_restriction_enzymes_filtered()`

- **作用：** 按位点长度、最适温度、供应商数量和互斥位点筛选限制性内切酶。
- **处理方法：** 用完整参数集调用 DnaChisel 的 `list_common_enzymes`（DNAKit 的 `list_common_restriction_enzymes` 只保留单一位点长度和供应商数量），并对名称排序。`site_unlike` 给出若干酶名，其识别位点（所有 IUPAC 变体）不得出现在返回酶的位点之中。酶表是 DnaChisel 自带的表（Biopython Restriction），所以具体列表取决于已安装的版本。
- **API：** `list_common_restriction_enzymes_filtered(*, site_length: 'Sequence[int]' = (6,), opt_temp: 'Sequence[int]' = (37,), min_suppliers: 'int' = 1, site_unlike: 'Sequence[str]' = (), verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaChisel（extra `dnakit[bio-dnachisel]`）。
- **输入：** `site_length` 是 1 到 20 个 4 到 12 的整数。`opt_temp` 是 1 到 20 个 0 到 120 的整数（摄氏度）。`min_suppliers` 为 0 到 20。`site_unlike` 是酶名列表。
- **输出：** `site_length`、`opt_temp`、`min_suppliers`、`site_unlike`、`enzyme_count`、`enzymes`（排序后的名称）和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import list_common_restriction_enzymes_filtered

six = list_common_restriction_enzymes_filtered(verbose=False)
eight = list_common_restriction_enzymes_filtered(site_length=[8], verbose=False)
print("EcoRI" in six["enzymes"], "EcoRI" in eight["enzymes"], "NotI" in eight["enzymes"])
print(six["enzymes"] == sorted(six["enzymes"]), six["enzyme_count"] == len(six["enzymes"]))
```

- **示例输出：**

```text
True False True
True True
```

## 15) `dnakit.molecular_lab.load_dna_record()`

- **作用：** 从 FASTA、GenBank、EMBL 或 SnapGene（`.dna`）读取一条序列记录并返回。
- **处理方法：** `file_format="auto"` 时按扩展名选择读取器：`.gb`、`.gbk`、`.genbank` 是 GenBank，`.fa`、`.fasta`、`.fna` 是 FASTA，`.dna` 是 SnapGene。GenBank、EMBL 和 FASTA 经由 `Bio.SeqIO.read` 读取，因此含多条记录的文件会报错；SnapGene 文件经由 `snapgene_reader` 读取。序列保持文件中的大小写。
- **API：** `load_dna_record(source: 'str | Path', *, file_format: 'str' = 'auto', name: 'str | None' = None, linear: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 读取 FASTA、GenBank 和 EMBL 需要 Biopython（extra `dnakit[bio-formats]`），读取 `.dna` 文件需要 snapgene_reader（extra `dnakit[bio-snapgene]`）。
- **输入：** `source` 是已存在的文件，不超过 50 MB。`file_format` 为 `auto`、`genbank`、`embl`、`fasta` 或 `snapgene`。给出 `name` 时，它会替换记录的 id 和名称。`linear` 只决定返回的拓扑标记；不会去检查文件本身的拓扑。
- **输出：** `format`、`id`、`name`、`description`、`sequence`、`length`、`topology`、`feature_count`、`source`（原样给出的路径）、`status`、`operation` 和 `provenance`。
- **示例代码：**

```python
from dnakit.molecular_lab import load_dna_record

record = load_dna_record("part.fa", verbose=False)
print(record["format"], record["id"], record["length"], record["topology"])
print(record["sequence"])
```

- **示例输出：**

```text
fasta part1 25 linear
ATGGGTCTCATGCCCAAAGCTTTAA
```

## 16) `dnakit.molecular_lab.plot_constraint_breaches_png()`

- **作用：** 绘制一条序列并高亮其约束违反位置，保存为 PNG。
- **处理方法：** 为这条序列建立约束违反表，用 DnaChisel 的 `records_from_breaches_dataframe` 转换，再用 `plot_breaches_record`（matplotlib 与 DNA Features Viewer）绘制。违反位置从表中解析出来，与图像路径一起返回。
- **API：** `plot_constraint_breaches_png(sequence: 'str', *, output_path: 'str | Path', constraints: 'Sequence[Mapping[str, Any]]', name: 'str' = 'sequence', figure_width: 'float' = 10.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaChisel（extra `dnakit[bio-dnachisel]`）、matplotlib（extra `dnakit[bio-plots]`）和 DNA Features Viewer（extra `dnakit[bio-dna-features-viewer]`）。
- **输入：** `sequence` 是 IUPAC DNA 字符串。`output_path` 必须不存在且以 `.png` 结尾。`constraints` 是非空的规格映射列表，写法同第 9 节。`figure_width` 为 2 到 40。
- **输出：** `image_path`、`sequence_length`、`breaches`（对每个约束标签给出 `[start, end, strand]` 的列表，其中 strand 为 `1`、`-1` 或 `None`）、`constraint_labels`、`status`、`operation` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import plot_constraint_breaches_png

constraints = [
    {"name": "AvoidPattern", "pattern": "GGTCTC"},
    {"name": "EnforceGCContent", "mini": 0.3, "maxi": 0.7, "window": 18},
]
with tempfile.TemporaryDirectory() as folder:
    result = plot_constraint_breaches_png(
        "ATGGGTCTCATGCCCAAAGAGACCAA",
        output_path=Path(folder) / "breaches.png",
        constraints=constraints,
        verbose=False,
    )
    print(Path(result["image_path"]).read_bytes()[:4])
print(result["breaches"])
```

- **示例输出：**

```text
b'\x89PNG'
{'GGTCTC': [[3, 9, 1], [18, 24, -1]], 'GC outside 30-70%/18bp': []}
```

## 17) `dnakit.molecular_lab.plot_optimization_changes_png()`

- **作用：** 优化一条序列，并保存显示变化区域的 PNG（DnaChisel 的 `plot_optimization_changes`）。
- **处理方法：** 建立问题，解决约束（有目标时再优化目标），然后调用 DnaChisel 的 `plot_optimization_changes` 并保存图。变化的片段通过比较优化前后的序列得到。
- **API：** `plot_optimization_changes_png(sequence: 'str', *, output_path: 'str | Path', constraints: 'Sequence[Mapping[str, Any]]' = (), objectives: 'Sequence[Mapping[str, Any]]' = (), circular: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaChisel（extra `dnakit[bio-dnachisel]`）、geneblocks（extra `dnakit[bio-geneblocks]`，同时会安装 DNA Features Viewer）和 matplotlib（extra `dnakit[bio-plots]`）。编写本页所用的环境没有安装 geneblocks，所以示例展示相应的 `BackendUnavailableError`。装好依赖后的用法见示例之后。
- **输入：** `sequence` 是 IUPAC DNA 字符串。`output_path` 必须不存在且以 `.png` 结尾。`constraints` 与 `objectives` 是规格映射列表，写法同第 9 节；两者至少有一个非空。`circular=True` 把序列当作环状。约束无法满足时抛出 `ConfigurationError`。
- **输出：** `image_path`、`sequence_length`、`changed_segments`（左闭右开的 `[start, end)`）、`edits`、`all_constraints_pass`、`status`、`operation` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import plot_optimization_changes_png

with tempfile.TemporaryDirectory() as folder:
    try:
        plot_optimization_changes_png(
            "ATGGGTCTCATGCCCAAA",
            output_path=Path(folder) / "changes.png",
            constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, "geneblocks could not be imported" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

安装依赖后的用法：

```python
from dnakit.molecular_lab import plot_optimization_changes_png

result = plot_optimization_changes_png(
    "ATGGGTCTCATGCCCAAA",
    output_path="changes.png",
    constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
    verbose=False,
)
print(result["changed_segments"], result["edits"])
```

## 18) `dnakit.molecular_lab.write_constraint_breaches_pdf()`

- **作用：** 为每条序列写一页约束违反图到一个新的 PDF，并同时返回违反表。
- **处理方法：** 为所有序列建立违反表，用 `records_from_breaches_dataframe` 转换，再用 DnaChisel 的 `breaches_records_to_pdf`（matplotlib 与 DNA Features Viewer）写出各页。
- **API：** `write_constraint_breaches_pdf(sequences: 'object', *, output_path: 'str | Path', constraints: 'Sequence[Mapping[str, Any]]', figure_width: 'float' = 10.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaChisel（extra `dnakit[bio-dnachisel]`）、matplotlib（extra `dnakit[bio-plots]`）和 DNA Features Viewer（extra `dnakit[bio-dna-features-viewer]`）。
- **输入：** `sequences` 与 `constraints` 同第 9 节。`output_path` 必须不存在且以 `.pdf` 结尾。`figure_width` 为 2 到 40。
- **输出：** `pdf_path`、`page_count`、`sequence_count`、`constraint_labels`、`rows`（违反表）、`status`、`operation` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import write_constraint_breaches_pdf

constraints = [{"name": "AvoidPattern", "pattern": "GGTCTC"}]
with tempfile.TemporaryDirectory() as folder:
    result = write_constraint_breaches_pdf(
        {"a": "ATGGGTCTCATGCCCAAA", "b": "ATGGAGACCAAGCCCAAA"},
        output_path=Path(folder) / "breaches.pdf",
        constraints=constraints,
        verbose=False,
    )
    print(Path(result["pdf_path"]).read_bytes()[:4])
print(result["page_count"], result["rows"])
```

- **示例输出：**

```text
b'%PDF'
2 [{'sequence': 'a', 'GGTCTC': '3-9(+)'}, {'sequence': 'b', 'GGTCTC': '3-9(-)'}]
```

## 19) `dnakit.molecular_lab.write_dnachisel_no_solution_report()`

- **作用：** 报告约束在何处无法满足，或报告它们可以被满足。
- **处理方法：** 建立问题并运行 `resolve_constraints()`（给出目标时再运行 `optimize()`）。如果 DnaChisel 抛出带位置的 `NoSolutionError`，就调用 DnaChisel 的 `write_no_solution_report`，把 `plots.pdf`、局部违反 GenBank 文件和日志写入新的 `output_dir`，并返回违反位置。如果错误没有位置，状态为 `no_solution_without_location`，且不写任何文件。如果约束可以满足，状态为 `solved`，同样不写任何文件。
- **API：** `write_dnachisel_no_solution_report(sequence: 'str', *, output_dir: 'str | Path', constraints: 'Sequence[Mapping[str, Any]]' = (), objectives: 'Sequence[Mapping[str, Any]]' = (), circular: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaChisel（extra `dnakit[bio-dnachisel]`）、matplotlib（extra `dnakit[bio-plots]`）和 DNA Features Viewer（extra `dnakit[bio-dna-features-viewer]`）。不需要报告包 `pdf_reports` 和 `sequenticon`。
- **输入：** `sequence` 是 IUPAC DNA 字符串。`output_dir` 必须不存在（其父目录会被创建）。`constraints` 与 `objectives` 是规格映射列表，写法同第 9 节；两者至少有一个非空。`circular=True` 把序列当作环状。
- **输出：** `status`（`solved`、`no_solution` 或 `no_solution_without_location`）、`report_path`（只有写了报告时才不是 `None`）、`report_files`、`message`、`location`（`[start, end, strand]` 或 `None`，按 DnaChisel 报告的值给出；下面示例中它的终点 11 超出了 6 个碱基的序列）、`circular`、`operation` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import write_dnachisel_no_solution_report

with tempfile.TemporaryDirectory() as folder:
    solvable = write_dnachisel_no_solution_report(
        "ATGGGTCTCATGCCCAAA",
        output_dir=Path(folder) / "solvable",
        constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
        verbose=False,
    )
    print(solvable["status"], solvable["report_files"])
    stuck = write_dnachisel_no_solution_report(
        "GGTCTC",
        output_dir=Path(folder) / "stuck",
        constraints=[
            {"name": "AvoidPattern", "pattern": "GGTCTC"},
            {"name": "AvoidChanges", "location": [0, 6]},
        ],
        verbose=False,
    )
    print(stuck["status"], stuck["report_files"], stuck["location"])
```

- **示例输出：**

```text
solved []
no_solution ['local_constraints_breaches.gb', 'logs.txt', 'plots.pdf'] [0, 11, 1]
```

## 20) `dnakit.molecular_lab.write_dnachisel_optimization_report()`

- **作用：** 优化一条序列，并写出 DnaChisel 的优化报告文件夹。
- **处理方法：** 建立问题，解决约束，有目标时再优化目标，然后调用 DnaChisel 的 `write_optimization_report`。文件夹中包含 PDF 摘要、优化前后的约束表与目标表、GenBank 文件和日志。两张 CSV 表会被读回并放入结果。
- **API：** `write_dnachisel_optimization_report(sequence: 'str', *, output_dir: 'str | Path', constraints: 'Sequence[Mapping[str, Any]]' = (), objectives: 'Sequence[Mapping[str, Any]]' = (), circular: 'bool' = False, project_name: 'str' = 'DnaChisel optimization', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaChisel，以及 DnaChisel 的报告包 `pdf_reports` 和 `sequenticon`（extra `dnakit[bio-dnachisel-reports]`）。编写本页所用的环境没有安装它们，所以示例展示相应的 `BackendUnavailableError`。装好依赖后的用法见示例之后。
- **输入：** `sequence` 是 IUPAC DNA 字符串。`output_dir` 必须不存在。至少需要一个约束或目标，写法为规格映射，同第 9 节。`circular` 把序列当作环状。`project_name` 是长度不超过 120 的非空字符串。约束无法满足时抛出 `ConfigurationError`；这种情况请用 `write_dnachisel_no_solution_report`。
- **输出：** `report_path`、`report_files`、`sequence_length`、`edits`、`all_constraints_pass`、`constraints_before_after` 与 `objectives_before_after`（两张 CSV 表的行）、`circular`、`status`、`operation` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import write_dnachisel_optimization_report

with tempfile.TemporaryDirectory() as folder:
    try:
        write_dnachisel_optimization_report(
            "ATGGGTCTCATGCCCAAA",
            output_dir=Path(folder) / "report",
            constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, "pdf_reports could not be imported" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

安装依赖后的用法：

```python
from dnakit.molecular_lab import write_dnachisel_optimization_report

result = write_dnachisel_optimization_report(
    "ATGGGTCTCATGCCCAAA",
    output_dir="dnachisel_report",
    constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
    verbose=False,
)
print(result["edits"], result["all_constraints_pass"], result["report_files"])
```

## 21) `dnakit.molecular_lab.optimize_dnachisel_source_record()`

- **作用：** 按记录自身注释的约束和目标进行优化，并写出结果。
- **处理方法：** 用 `DnaOptimizationProblem.from_record`（`circular=True` 时用 `CircularDnaOptimizationProblem`）载入记录，它从记录的特征标签中读取约束和目标，与 DnaChisel 命令行一致。目标路径以 `.gb` 结尾时得到带注释的 GenBank 记录：先解决约束，再优化目标，然后用 `to_record` 写文件。其他目标路径则是新的文件夹或 zip 文件，用 `optimize_with_report` 写出。优化前后各取一次约束与目标的文本摘要。
- **API：** `optimize_dnachisel_source_record(source_path: 'str | Path', target_path: 'str | Path', *, circular: 'bool' = False, with_sequence_edits: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaChisel（extra `dnakit[bio-dnachisel]`）。非 `.gb` 的目标（文件夹或 zip 报告）还需要报告包（extra `dnakit[bio-dnachisel-reports]`）。
- **输入：** `source_path` 是已存在的 GenBank（或 DnaChisel 能读取的其他）记录，其特征按 DnaChisel 规格命名，例如 `@AvoidPattern(GGTCTC)`。`target_path` 必须不存在，其父目录必须存在。`with_sequence_edits=True` 会在 GenBank 输出中标注被修改的位置。
- **输出：** `source`（文件名）、`output_path`、`output_kind`（`genbank` 或 `report`）、`circular`、`before` 与 `after`（各含约束和目标的文本摘要 `constraints`、`objectives`）、`status`、`operation` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import optimize_dnachisel_source_record

with tempfile.TemporaryDirectory() as folder:
    target = Path(folder) / "optimized.gb"
    result = optimize_dnachisel_source_record("annotated.gb", target, verbose=False)
    lines = target.read_text().splitlines()
print(result["output_kind"], result["source"])
print(result["before"]["constraints"].splitlines()[0])
print(result["after"]["constraints"].splitlines()[0])
print([line for line in lines if line.startswith("        1")][0])
```

- **示例输出：**

```text
genbank annotated.gb
===> FAILURE: 1 constraints evaluations failed
===> SUCCESS - all constraints evaluations pass
        1 atgggtctaa tgcccaaagc tttaa
```

## 22) `dnakit.molecular_lab.dnacauldron_plot_leveled_graph()`

- **作用：** 把分层的装配图（节点按层成列，层间有边）画成 PNG。
- **处理方法：** 调用 DnaCauldron 的 `AssemblyPlan.plot_leveled_graph` 布局，并传入一个把每个节点画成带标签圆角框的回调，在 Agg 后端的 matplotlib 图上绘制并保存为 PNG（150 dpi，紧凑边框）。布局给出的节点位置以占绘图区比例的形式返回。
- **API：** `dnacauldron_plot_leveled_graph(output_path: 'str', levels: 'list[list[str]]', edges: 'list[list[str]]', *, width_factor: 'float' = 2.5, height_factor: 'float' = 2.0, scale: 'float' = 1.0, interlevel_shift: 'float' = 0.0, margin: 'float | None' = None, font_size: 'float' = 9.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaCauldron（extra `dnakit[bio-dnacauldron]`）和 matplotlib（extra `dnakit[bio-plots]`）。
- **输入：** `output_path` 以 `.png` 结尾且不存在。`levels` 是 1 到 20 个非空列表，由互不重复的节点名组成（每个名字 1 到 64 个字符，总计最多 500 个节点）。`edges` 是 `[source, target]` 对的列表，两端必须是 `levels` 中的节点。`width_factor` 与 `height_factor` 为 0.1 到 50，`scale` 为 0.01 到 20，`interlevel_shift` 为 -1 到 1，`margin` 为 0 到 1 或 `None`，`font_size` 为 1 到 72。
- **输出：** `path` 和 `data`，其中含 `level_count`、`node_count`、`edge_count`、`levels`、`edges`、`positions`（节点到 `[x, y]`）和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import dnacauldron_plot_leveled_graph

levels = [["partA", "partB", "partC"], ["construct_1"]]
edges = [["partA", "construct_1"], ["partB", "construct_1"], ["partC", "construct_1"]]
with tempfile.TemporaryDirectory() as folder:
    result = dnacauldron_plot_leveled_graph(
        str(Path(folder) / "plan.png"), levels, edges, verbose=False
    )
    print(Path(result["path"]).read_bytes()[:4])
data = result["data"]
print(data["level_count"], data["node_count"], data["edge_count"])
print({node: [round(value, 3) for value in xy] for node, xy in data["positions"].items()})
```

- **示例输出：**

```text
b'\x89PNG'
2 4 3
{'partA': [0.333, 0.25], 'partB': [0.333, 0.5], 'partC': [0.333, 0.75], 'construct_1': [0.667, 0.5]}
```

## 23) `dnakit.molecular_lab.dnacauldron_write_simulation_pdf_report()`

- **作用：** 把 DnaCauldron 的模拟摘要写成 PDF 报告。
- **处理方法：** 检查输出路径和各行数据，导入 `pdf_reports`，构建含 `Outcome` 与 `Number of assemblies` 两列的 pandas 数据框，再调用 DnaCauldron 的 `reports.write_simulation_pdf_report`。这两列正是 DnaCauldron 的 `AssemblyPlanSimulation` 所建表的列。
- **API：** `dnacauldron_write_simulation_pdf_report(output_path: 'str', simulation_info: 'list[dict[str, Any]]', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaCauldron（extra `dnakit[bio-dnacauldron]`）及其 reports extra（`pip install "dnacauldron[reports]"`），后者带来 `pdf_reports`，并通过它带来 WeasyPrint。编写本页所用的环境没有安装它们，所以示例展示相应的 `BackendUnavailableError`。装好依赖后的用法见示例之后。
- **输入：** `output_path` 以 `.pdf` 结尾且不存在。`simulation_info` 是 1 到 50 个字典的列表，键恰为 `Outcome`（不超过 64 个字符的字符串，例如 `Valid`、`Cancelled` 或 `Errored`）和 `Number of assemblies`（非负整数）。
- **输出：** `path` 和 `data`，其中含 `rows`、`total_assemblies` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.molecular_lab import dnacauldron_write_simulation_pdf_report

rows = [
    {"Outcome": "Valid", "Number of assemblies": 3},
    {"Outcome": "Errored", "Number of assemblies": 1},
]
with tempfile.TemporaryDirectory() as folder:
    try:
        dnacauldron_write_simulation_pdf_report(
            str(Path(folder) / "summary.pdf"), rows, verbose=False
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, "pdf_reports could not be imported" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

安装依赖后的用法：

```python
from dnakit.molecular_lab import dnacauldron_write_simulation_pdf_report

result = dnacauldron_write_simulation_pdf_report(
    "summary.pdf",
    [{"Outcome": "Valid", "Number of assemblies": 3}],
    verbose=False,
)
print(result["data"]["total_assemblies"])
```

## 24) `dnakit.molecular_lab.simulate_dnacauldron_plan()`

- **作用：** 模拟 DnaCauldron 方案中的每一个装配，并写出模拟报告。
- **处理方法：** 把 `sequences_path` 中的部件序列导入 `SequenceRepository`，把 `connectors_path` 的每个子文件夹作为一个命名集合导入，用 `AssemblyPlan.from_spreadsheet` 读取方案，运行 `simulate`，取 `compute_stats()`，并用 `write_report` 写出报告，与 DnaCauldron 命令行的做法一致。
- **API：** `simulate_dnacauldron_plan(assembly_plan_path: 'str | Path', sequences_path: 'str | Path', target_path: 'str | Path', *, connectors_path: 'str | Path | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 DnaCauldron（extra `dnakit[bio-dnacauldron]`）。在本环境中，写出报告文件夹不需要 `pdf_reports`。
- **输入：** `assembly_plan_path` 是已存在的 CSV 或电子表格文件（不超过 50 MB），每行一个装配，其中以 `class: ...` 单元格指定装配类。`sequences_path` 是已存在的部件序列文件夹。`connectors_path` 可选：其每个子文件夹存放一个连接子集合。`target_path` 必须不存在，但其父目录必须存在；要生成 zip 报告请使用 `.zip` 名称。
- **输出：** `assembly_plan`（文件名）、`assemblies`（装配数）、`outcomes`（`valid_assemblies`、`cancelled_assemblies`、`errored_assemblies`）、`report_path`、`status`、`operation` 和 `provenance`。
- **示例代码：**

```python
import random
import tempfile
from pathlib import Path

from dnakit.molecular_lab import simulate_dnacauldron_plan

rng = random.Random(3)


def dna(length):
    return "".join(rng.choice("ACGT") for _ in range(length))


o0, o1, o2 = dna(40), dna(40), dna(40)  # shared overlaps between neighbouring parts
parts = {"partA": o2 + dna(160) + o0, "partB": o0 + dna(160) + o1, "partC": o1 + dna(160) + o2}
with tempfile.TemporaryDirectory() as folder:
    root = Path(folder)
    (root / "parts").mkdir()
    for name, sequence in parts.items():
        (root / "parts" / f"{name}.fa").write_text(f">{name}\n{sequence}\n")
    (root / "plan.csv").write_text("class: gibson_assembly,construct_1,partA,partB,partC\n")
    result = simulate_dnacauldron_plan(
        root / "plan.csv", root / "parts", root / "report", verbose=False
    )
    summary = root / "report" / "plan_simulation" / "plan_summary.csv"
    print(summary.is_file())
print(result["assembly_plan"], result["assemblies"])
print(result["outcomes"])
```

- **示例输出：**

```text
True
plan.csv 1
{'cancelled_assemblies': 0, 'errored_assemblies': 0, 'valid_assemblies': 1}
```
