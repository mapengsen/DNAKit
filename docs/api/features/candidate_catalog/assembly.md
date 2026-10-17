# 基因组组装与泛基因组

本类收录 **37** 项有文档的功能。其中 27 项直接由输入用纯 Python 在本地计算出结果（5 个连续性统计函数 `calculate_l50`、`calculate_n50`、`calculate_ng50`、`calculate_nx`、`calculate_stats`，其小节保持不变；另有 22 个函数，涉及泛基因组矩阵、k-mer 模型、组装检查与比较、重复序列与端粒、GFA 图、解析、按显式布局的脚手架构建以及污染 contig 的移除）。9 项始终返回 `status: "not_implemented"`，因为它们需要输入中并不包含的测序读段、比对、参考基因组或数据库（`build_chromosome_level_assembly`、`build_haplotype_assembly`、`detect_assembly_contamination`、`detect_chimeric_contigs`、`detect_misassemblies`、`identify_centromeres`、`phase_assembly`、`polish`、`purge_haplotigs`）。`run` 会执行调用者自行安装的组装程序（Flye、hifiasm、MEGAHIT、Shasta 或 Unicycler），在该程序不在 PATH 中时返回 `not_implemented`。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) add_genome_to_pangenome

- **作用：** 把一个基因组的基因集合作为新的一列加入基因有/无矩阵。

- **处理方法：** 读取已有的泛基因组矩阵，为新基因组追加一列（基因在 `genes` 中则为 1，否则为 0），并为此前未见过的每个基因追加一行（对旧基因组均为缺失；新基因按字母顺序排列）。不做同源搜索：基因名即其所属簇，因此以不同名称出现的同一直系同源基因会被当作新基因。新基因组的名称不能已存在。

- **API：** `dnakit.assembly.add_genome_to_pangenome(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：基因有/无数据：`{基因组: [基因]}`、`{基因: {基因组: 0/1}}` 或 `{"genomes": [...], "genes": [...], "matrix": [[...]]}`（行为基因，列为基因组；取值为布尔值或非负数，> 0 表示存在）。基因名被视为已聚类好的同源组；不做序列聚类或同源搜索。`options`（必需的键）：`genome`（字符串，新基因组的唯一名称）、`genes`（字符串列表，新基因组中存在的基因）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：更新后泛基因组的 `genomes`、`genes`、`matrix`（行为基因，列为基因组，取值 0/1），以及 `new_gene_count`、`new_genes`、`pangenome_size_before`、`pangenome_size_after`。

- **示例代码：**

```python
import json

from dnakit.assembly import add_genome_to_pangenome

pangenome = {"A": ["g1", "g2"], "B": ["g2", "g3"]}
result = add_genome_to_pangenome(pangenome, options={"genome": "C", "genes": ["g2", "g4"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "add_genome_to_pangenome",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 38,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6d290492cc572135bfe94b6f4b9597f82691599f0ace9fd6682f9798aab95803"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genomes": ["A", "B", "C"],
    "genes": ["g1", "g2", "g3", "g4"],
    "matrix": [
      [1, 0, 0],
      [1, 1, 1],
      [0, 1, 0],
      [0, 0, 1]
    ],
    "new_gene_count": 1,
    "new_genes": ["g4"],
    "pangenome_size_before": 3,
    "pangenome_size_after": 4
  }
}
```

## 2) build_chromosome_level_assembly

- **作用：** 染色体水平组装尚未实现：把 contig 排序、定向成染色体需要 Hi-C、光学图谱或共线性证据。

- **处理方法：** 尚未实现。把 contig 连成染色体需要连接证据（Hi-C 接触、光学图谱或参考基因组）以及 YaHS、SALSA 之类的脚手架工具，这些都不在输入中；如果顺序和方向已知，请使用带显式布局的 `scaffold`。不会用启发式方法代替，因此函数不会编造结果；除了在 `input_summary` 中概括输入外，不对输入做任何校验，只返回结构化的 `not_implemented` 封装。

- **API：** `dnakit.assembly.build_chromosome_level_assembly(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：任意内容（只在 `input_summary` 中被概括）。`options`（可选）：被忽略。

- **输出：** 始终返回 `status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（缺少的东西：脚手架证据（Hi-C 接触、光学图谱或参考基因组）以及 YaHS/SALSA 之类的脚手架工具）和 `message`（说明为何无法仅凭输入计算）。不产生 `result`。

- **示例代码：**

```python
import json

from dnakit.assembly import build_chromosome_level_assembly

result = build_chromosome_level_assembly(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "build_chromosome_level_assembly",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "scaffolding evidence (Hi-C contacts, optical maps or a reference) and a scaffolder such as YaHS/SALSA",
  "message": "Ordering and orienting contigs into chromosomes needs linking evidence that is not part of the input; use scaffold() with an explicit layout instead."
}
```

## 3) build_command

- **作用：** 构建（但不运行）某个受支持组装程序（Flye、hifiasm、MEGAHIT、Shasta 或 Unicycler）的命令行参数列表。

- **处理方法：** 按照 `assemble_flye`、`assemble_hifiasm`、`assemble_megahit`、`assemble_shasta`、`assemble_unicycler` 这些封装所用的规则拼装参数向量，并校验各工具必需的键（例如 Flye 需要 `reads` 和 `genome_size`，Unicycler 需要 `read1` 或 `long_reads`）及取值范围。不执行任何程序，不读取也不创建文件，程序也不需要已安装。

- **API：** `dnakit.assembly.build_command(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：设置映射；`options` 会合并到其上并优先生效。键：`tool`（必填：`flye`、`hifiasm`、`megahit`、`shasta` 或 `unicycler`）、`output_directory`（必填）、`program`（可执行文件名，默认等于工具名）、`threads`（1-256，默认 1）、`reads`（flye、hifiasm、shasta）、`read1` / `read2`（megahit、unicycler）、`long_reads`（unicycler）、`genome_size`（flye，例如 `"5m"`）、`platform`（flye，默认 `"nano-hq"`）、`config`（shasta，默认 `"Nanopore-May2022"`）、`mode`（unicycler，默认 `"normal"`）、`prefix`（hifiasm，默认 `"assembly"`）、`min_contig_length`（megahit）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`tool`、`argv`（参数列表，第一项为程序名）和 `command_line`（经过 shell 引号处理的命令字符串）。

- **示例代码：**

```python
import json

from dnakit.assembly import build_command

settings = {
    "tool": "flye",
    "reads": "reads.fastq",
    "genome_size": "5m",
    "output_directory": "flye_out",
    "threads": 4,
}
result = build_command(settings)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "build_command",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 107,
    "sequence_bases": 10,
    "lengths": [10],
    "sha256": "b9a386c642fbd16778c7d5400bc52c7bba3407abbe68afa7505c06e2602b5108"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "tool": "flye",
    "argv": ["flye", "--nano-hq", "reads.fastq", "--genome-size", "5m", "--out-dir", "flye_out", "--threads", "4"],
    "command_line": "flye --nano-hq reads.fastq --genome-size 5m --out-dir flye_out --threads 4"
  }
}
```

## 4) build_graph_genome

- **作用：** 由输入序列构建压缩的正向链 de Bruijn 图（unitig，以 GFA1 格式输出）。

- **处理方法：** 每条序列的每个 k-mer（在 A/C/G/T 以外的字符处切开）是一个节点；连续观察到的两个 k-mer（即有 k-1 个碱基重叠）之间连一条边。unitig 是该图中最大的无分支路径（孤立的环也会被报告）；每个 unitig 记录它出现在哪些输入序列中（`sources`）。反向互补链不会被合并（该图不是双向图，一条序列与其反向互补会产生不同的节点），也不做比对或变异图构建。输入最多 2,000,000 个碱基。

- **API：** `dnakit.assembly.build_graph_genome(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：FASTA 文本或路径、`{名称: 序列}` 或序列列表（每条是一个输入序列，无名称时命名为 `contig_<i>`）。`options`（可选）：`k`（整数 3-127，默认 21）、`output`（写入 GFA1 文本的可选文件路径）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`k`、`kmer_count`、`unitig_count`、`link_count`、`unitigs`（每项含 `name`（`utg<i>`）、`sequence`、`length`、`kmers`、`sources`）、`links`（`from`、`to`、`overlap`）、`gfa`（GFA1 文本：带 `LN` 标签的 `S` 行和带 `(k-1)M` 重叠的 `L` 行，方向均为 `+`）、`note` 和 `output_path`。

- **示例代码：**

```python
import json

from dnakit.assembly import build_graph_genome

result = build_graph_genome({"s1": "ACGTT", "s2": "ACGTA"}, options={"k": 3})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "build_graph_genome",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 30,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "db2285ea0f9ff1176712fd147fdb9ef5c1c1bf968cc8ca6a4a4ecdbe64e768ad"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "k": 3,
    "kmer_count": 4,
    "unitig_count": 3,
    "link_count": 2,
    "unitigs": [
      {
        "name": "utg1",
        "sequence": "ACGT",
        "length": 4,
        "kmers": 2,
        "sources": ["s1", "s2"]
      },
      {
        "name": "utg2",
        "sequence": "GTA",
        "length": 3,
        "kmers": 1,
        "sources": ["s2"]
      },
      {
        "name": "utg3",
        "sequence": "GTT",
        "length": 3,
        "kmers": 1,
        "sources": ["s1"]
      }
    ],
    "links": [
      {
        "from": "utg1",
        "to": "utg2",
        "overlap": 2
      },
      {
        "from": "utg1",
        "to": "utg3",
        "overlap": 2
      }
    ],
    "gfa": "H\tVN:Z:1.0\nS\tutg1\tACGT\tLN:i:4\nS\tutg2\tGTA\tLN:i:3\nS\tutg3\tGTT\tLN:i:3\nL\tutg1\t+\tutg2\t+\t2M\nL\tutg1\t+\tutg3\t+\t2M\n",
    "note": "forward-strand (non-canonical) de Bruijn graph; nodes are k-mers, edges are observed (k-1)-overlaps; unitigs are maximal non-branching paths",
    "output_path": null
  }
}
```

## 5) build_haplotype_assembly

- **作用：** 单倍型分辨组装尚未实现：它需要测序读段（HiFi/ONT、Hi-C 或亲本 k-mer）和支持定相的组装程序。

- **处理方法：** 尚未实现。单倍型分辨的组装不能仅由 contig 得出；它需要测序读段和 hifiasm 这样支持定相的组装程序（见 `build_command` / `run`）。不会用启发式方法代替，因此函数不会编造结果；除了在 `input_summary` 中概括输入外，不对输入做任何校验，只返回结构化的 `not_implemented` 封装。

- **API：** `dnakit.assembly.build_haplotype_assembly(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：任意内容（只在 `input_summary` 中被概括）。`options`（可选）：被忽略。

- **输出：** 始终返回 `status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（缺少的东西：支持定相的组装程序（例如 hifiasm）和测序读段）和 `message`（说明为何无法仅凭输入计算）。不产生 `result`。

- **示例代码：**

```python
import json

from dnakit.assembly import build_haplotype_assembly

result = build_haplotype_assembly(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "build_haplotype_assembly",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "an assembler with phasing support (for example hifiasm) and sequencing reads",
  "message": "Haplotype-resolved assembly cannot be derived from contigs alone; run assemble_hifiasm on reads instead."
}
```

## 6) build_pangenome

- **作用：** 由各基因组的基因集合构建基因有/无泛基因组矩阵，并给出核心（core）、壳（shell）、云（cloud）基因的数目。

- **处理方法：** 把基因取并集，得到 基因 x 基因组 的 0/1 矩阵，并按携带该基因的基因组比例划分：核心 = 至少存在于 `ceil(core_threshold * 基因组数)` 个基因组；其余的辅助基因中，存在比例低于 `cloud_threshold` 的为云基因，其余为壳基因。同时统计每个基因组的基因数和仅它独有的基因数。基因名被视为已聚类好的同源组；不做序列聚类。

- **API：** `dnakit.assembly.build_pangenome(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：基因有/无数据：`{基因组: [基因]}`、`{基因: {基因组: 0/1}}` 或 `{"genomes": [...], "genes": [...], "matrix": [[...]]}`（行为基因，列为基因组；取值为布尔值或非负数，> 0 表示存在）。基因名被视为已聚类好的同源组；不做序列聚类或同源搜索。`options`（可选）：`core_threshold`（(0, 1] 内的比例，默认 0.99）、`cloud_threshold`（[0, 1] 内的比例，默认 0.15）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`genomes`、`genes`、`matrix`、`pangenome_size`、`core_count`、`shell_count`、`cloud_count`、`core_threshold`、`genes_per_genome` 和 `unique_genes_per_genome`。

- **示例代码：**

```python
import json

from dnakit.assembly import build_pangenome

genomes = {"G1": ["a", "b", "c"], "G2": ["a", "b", "d"], "G3": ["a", "e"]}
result = build_pangenome(genomes, options={"cloud_threshold": 0.5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "build_pangenome",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 64,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e4372679c8d0e6836c36cd361cef3caf9656e4da871974b5c7cf4f7878ed448b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genomes": ["G1", "G2", "G3"],
    "genes": ["a", "b", "c", "d", "e"],
    "matrix": [
      [1, 1, 1],
      [1, 1, 0],
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ],
    "pangenome_size": 5,
    "core_count": 1,
    "shell_count": 1,
    "cloud_count": 3,
    "core_threshold": 0.99,
    "genes_per_genome": {
      "G1": 3,
      "G2": 3,
      "G3": 2
    },
    "unique_genes_per_genome": {
      "G1": 1,
      "G2": 1,
      "G3": 1
    }
  }
}
```

## 7) calculate_assembly_completeness

- **作用：** 以长度比例和/或参考序列 k-mer 的回收率来估计组装的完整度。

- **处理方法：** 长度模式：`length_ratio = 组装总长度 / reference_length`，`length_completeness = min(ratio, 1)`。k-mer 模式：`kmer_completeness` 为参考序列中不同的规范（canonical）k-mer（k-mer 与其反向互补视为同一个，只使用 A/C/G/T 片段）在组装中也出现的比例。这是长度 / k-mer 层面的度量，而不是基因内容层面的度量（见 `calculate_busco_completeness`）。`reference_length` 与 `reference` 至少要给一个；k-mer 模式需要序列，并且总长度限制为 5,000,000 个碱基。

- **API：** `dnakit.assembly.calculate_assembly_completeness(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：组装结果，可以是 FASTA 文本或 FASTA 文件路径、`{名称: 序列}`，或序列 / `{"name", "sequence"}` 记录的列表（也接受 contig 长度，形式为 `{名称: 长度}` 或整数列表）。`options`（`reference_length` / `reference` 至少给一个）：`reference_length`（正数，预期基因组大小）、`reference`（含序列的组装，形式同上）、`k`（整数 3-63，默认 21；k-mer 模式）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`assembly_length`、`contig_count`、`methods`（`"length"` 和/或 `"kmer"`），以及各模式的 `reference_length`、`length_ratio`、`length_completeness`；`k`、`reference_distinct_kmers`、`shared_kmers`、`kmer_completeness`。

- **示例代码：**

```python
import json

from dnakit.assembly import calculate_assembly_completeness

result = calculate_assembly_completeness(
    {"c1": "ACGTA"},
    options={"reference_length": 8, "reference": "ACGTACGT", "k": 4},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "calculate_assembly_completeness",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 15,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6fbe1ad637a90ef6b20a533c11a2dfc6669bbeaa5ffcabfcfda1d2880ded589e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "assembly_length": 5,
    "contig_count": 1,
    "reference_length": 8,
    "length_ratio": 0.625,
    "length_completeness": 0.625,
    "k": 4,
    "reference_distinct_kmers": 3,
    "shared_kmers": 2,
    "kmer_completeness": 0.6666666666666666,
    "methods": ["length", "kmer"]
  }
}
```

## 8) calculate_busco_completeness

- **作用：** 根据 BUSCO 计数或 BUSCO 简要摘要行计算 BUSCO 完整度百分比（不会运行 BUSCO）。

- **处理方法：** complete = 单拷贝（S）+ 重复（D）；各百分比为计数除以 n，其中 n 为摘要文本中给出的总数，或为 S + D + F + M（用计数输入时）。对 `C:98.5%[S:97.0%,D:1.5%],F:0.5%,M:1.0%,n:1614` 这样的摘要行，计数由 `round(百分比 * n / 100)` 还原。它只汇总已有的 BUSCO 结果；要运行 BUSCO 本身请使用 `run_busco`。

- **API：** `dnakit.assembly.calculate_busco_completeness(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：含单拷贝（`single` 或 `S`）、重复（`duplicated` / `D`）、碎片化（`fragmented` / `F`）和缺失（`missing` / `M`）整数计数的映射，或包含 `C:..%[S:..%,D:..%],F:..%,M:..%,n:..` 行的 BUSCO 简要摘要文本或路径。`options`：无。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`n`、`single`、`duplicated`、`complete`、`fragmented`、`missing`、`complete_percent`、`single_percent`、`duplicated_percent`、`fragmented_percent`、`missing_percent` 和 `summary`（BUSCO 风格的摘要字符串）。

- **示例代码：**

```python
import json

from dnakit.assembly import calculate_busco_completeness

summary = "# BUSCO\n\tC:98.5%[S:97.0%,D:1.5%],F:0.5%,M:1.0%,n:1614\n"
result = calculate_busco_completeness(summary)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "calculate_busco_completeness",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 54,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b543a46533d7392dd244a0fae5395eb128a9c9e83904ce3338824e9c42dda654"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n": 1614,
    "single": 1566,
    "duplicated": 24,
    "complete": 1590,
    "fragmented": 8,
    "missing": 16,
    "complete_percent": 98.51301115241635,
    "single_percent": 97.02602230483271,
    "duplicated_percent": 1.486988847583643,
    "fragmented_percent": 0.49566294919454773,
    "missing_percent": 0.9913258983890955,
    "summary": "C:98.5%[S:97.0%,D:1.5%],F:0.5%,M:1.0%,n:1614"
  }
}
```

## 9) calculate_l50

- **作用：** 为“基因组组装与泛基因组”类分析提供可复用的 `calculate_l50` 功能。候选编号 CAP-0061，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_l50` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/assembly/local.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的 contig 长度计算确定性的连续性统计，不调用组装器、组装图、参考基因组或 BUSCO 数据库。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.assembly.calculate_l50(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`options`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `calculate_l50` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的本地来源信息。该指标只根据正的 contig 长度计算，不调用组装器、组装图、参考基因组或 BUSCO 数据库。

- **示例代码：**

```python
from pathlib import Path

from dnakit.assembly import calculate_l50

result = calculate_l50(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0061",
  "operation": "calculate_l50",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_l50', 'result': {'l50': 1, 'n_at_percentile': 100, 'percentile': 50.0, 'contig_count': 3, 'total_length': 200}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scientific_boundary': 'contiguity...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 10) calculate_n50

- **作用：** 为“基因组组装与泛基因组”类分析提供可复用的 `calculate_n50` 功能。候选编号 CAP-0062，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_n50` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/assembly/local.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的 contig 长度计算确定性的连续性统计，不调用组装器、组装图、参考基因组或 BUSCO 数据库。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.assembly.calculate_n50(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`options`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `calculate_n50` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的本地来源信息。该指标只根据正的 contig 长度计算，不调用组装器、组装图、参考基因组或 BUSCO 数据库。

- **示例代码：**

```python
from pathlib import Path

from dnakit.assembly import calculate_n50

result = calculate_n50(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0062",
  "operation": "calculate_n50",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_n50', 'result': {'n50': 100, 'l50': 1, 'contig_count': 3, 'total_length': 200, 'percentile': 50.0, 'options_used': []}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scientific_boundary': 'co...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 11) calculate_ng50

- **作用：** 为“基因组组装与泛基因组”类分析提供可复用的 `calculate_ng50` 功能。候选编号 CAP-0063，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_ng50` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/assembly/local.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的 contig 长度计算确定性的连续性统计，不调用组装器、组装图、参考基因组或 BUSCO 数据库。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.assembly.calculate_ng50(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`options`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `calculate_ng50` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的本地来源信息。该指标只根据正的 contig 长度计算，不调用组装器、组装图、参考基因组或 BUSCO 数据库。

- **示例代码：**

```python
from pathlib import Path

from dnakit.assembly import calculate_ng50

result = calculate_ng50(
    input={"lengths": [100, 60, 40]},
    options={"genome_size": 250}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0063",
  "operation": "calculate_ng50",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_ng50', 'result': {'ng': 60, 'ng50': 60, 'contigs_to_threshold': 2, 'percentile': 50.0, 'expected_genome_size': 250, 'assembly_length': 200}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scie...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 12) calculate_nx

- **作用：** 为“基因组组装与泛基因组”类分析提供可复用的 `calculate_nx` 功能。候选编号 CAP-0064，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_nx` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/assembly/local.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的 contig 长度计算确定性的连续性统计，不调用组装器、组装图、参考基因组或 BUSCO 数据库。审核依据：`no_exact_public_match`。已有同名函数证据：`calculate_dn_ds`。

- **API：** `dnakit.assembly.calculate_nx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`options`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `calculate_nx` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的本地来源信息。该指标只根据正的 contig 长度计算，不调用组装器、组装图、参考基因组或 BUSCO 数据库。

- **示例代码：**

```python
from pathlib import Path

from dnakit.assembly import calculate_nx

result = calculate_nx(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0064",
  "operation": "calculate_nx",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_nx', 'result': {'percentile': 50.0, 'nx': 100, 'contigs_to_threshold': 1, 'contig_count': 3, 'total_length': 200}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scientific_boundary': 'contigu...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 13) calculate_pangenome_openness

- **作用：** 用 Heaps 定律拟合每加入一个基因组带来的新基因数，据此判断泛基因组是开放还是封闭。

- **处理方法：** 模型（Tettelin 等 2008）：加入第 N 个基因组时发现的新基因数的平均值为 `kappa * N^(-alpha)`。基因组以 `n_permutations` 种随机顺序依次加入（固定 `seed` 即可复现；`0` 表示只按给定顺序加入一次），对每个 N 求新基因数的平均，再在对数-对数坐标上用普通最小二乘对 N = 2..G 拟合 `alpha` 和 `kappa`（新基因数为 0 的 N 会被舍弃）。`alpha < 1` 判为开放，`alpha >= 1` 判为封闭，结果中的 `model` 字符串也给出同样的规则。至少需要 4 个基因组，且至少 3 个 N 值有新基因。第一个基因组（N = 1）不参与拟合。

- **API：** `dnakit.assembly.calculate_pangenome_openness(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：基因有/无数据：`{基因组: [基因]}`、`{基因: {基因组: 0/1}}` 或 `{"genomes": [...], "genes": [...], "matrix": [[...]]}`（行为基因，列为基因组；取值为布尔值或非负数，> 0 表示存在）。基因名被视为已聚类好的同源组；不做序列聚类或同源搜索（至少 4 个基因组）。`options`（可选）：`n_permutations`（整数 0-10,000，默认 100）、`seed`（整数，默认 0）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`alpha`、`kappa`、`points_used`、`points_dropped`、`open`（`alpha < 1` 时为 true）、`genome_count`、`permutations`、`seed`（`n_permutations` 为 0 时为 null）、`mean_new_genes_by_genome_number`、`mean_pangenome_size_by_genome_number` 以及说明字符串 `model`。

- **示例代码：**

```python
import json

from dnakit.assembly import calculate_pangenome_openness

genomes = {
    "g1": list("abcdefgh"),
    "g2": list("abcd") + list("ijkl"),
    "g3": list("ab") + list("mn"),
    "g4": ["a", "o"],
    "g5": ["a", "p"],
}
result = calculate_pangenome_openness(genomes, options={"n_permutations": 0})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "calculate_pangenome_openness",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 160,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e6415dfbb17e1fa60f462f0dea884a8f08171b89837cce9dfad1af2832d0e517"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "alpha": 1.6324718703854002,
    "kappa": 11.866529443446138,
    "points_used": 4,
    "points_dropped": 0,
    "open": false,
    "genome_count": 5,
    "permutations": 1,
    "seed": null,
    "mean_new_genes_by_genome_number": [8.0, 4.0, 2.0, 1.0, 1.0],
    "mean_pangenome_size_by_genome_number": [8.0, 12.0, 14.0, 15.0, 16.0],
    "model": "Heaps' law: mean number of new genes found when the N-th genome is added = kappa * N^(-alpha), fitted by OLS on log-log for N=2..G; open if alpha < 1, closed if alpha >= 1 (Tettelin et al. 2008; Rasko et al. 2008)"
  }
}
```

## 14) calculate_stats

- **作用：** 为“基因组组装与泛基因组”类分析提供可复用的 `calculate_stats` 功能。候选编号 CAP-0066，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_stats` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `assembler` 用于选择一种支持的变体（取值见源目录记录）。无外部依赖的实现位于 `src/dnakit/assembly/local.py`，并通过正式的 DNAKit 命名空间导出。它根据显式给出的 contig 长度计算确定性的连续性统计，不调用组装器、组装图、参考基因组或 BUSCO 数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.assembly.calculate_stats(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`options`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 本地函数 `calculate_stats` 返回已完成的可 JSON 序列化封装，包含 `status`、`operation`、`result` 和明确的本地来源信息。该指标只根据正的 contig 长度计算，不调用组装器、组装图、参考基因组或 BUSCO 数据库。

- **示例代码：**

```python
from pathlib import Path

from dnakit.assembly import calculate_stats

result = calculate_stats(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "capability_id": "CAP-0066",
  "operation": "calculate_stats",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_stats', 'result': {'contig_count': 3, 'total_length': 200, 'min_length': 40, 'max_length': 100, 'mean_length': 66.66666666666667, 'n50': 100, 'l50': 1, 'n90': 40, 'l90': 3, 'sorted_lengths': [100, 60, 40], 'options_used': []}, 'provenance': {'implem...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 15) check_input

- **作用：** 校验组装输入，并报告检测到的格式、contig 数、总长度、重复的 contig 名称以及非 IUPAC 符号。

- **处理方法：** 先把输入解析为组装结果（序列或仅有长度），再检查 contig 名称是否唯一，以及每条序列是否只含 IUPAC 核苷酸符号 `ACGTNRYSWKMBDHV`（不区分大小写；`-` 和 `*` 会被报告为无效）。无法使用的输入（空输入、类型错误、不是 FASTA）不会抛出异常，而是报告 `valid: false` 并在 `errors` 中给出原因。

- **API：** `dnakit.assembly.check_input(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：组装结果，可以是 FASTA 文本或 FASTA 文件路径、`{名称: 序列}`，或序列 / `{"name", "sequence"}` 记录的列表（也接受 contig 长度，形式为 `{名称: 长度}` 或整数列表）。`options`：无。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`valid`、`format`（`"sequences"`、`"lengths"`，无法使用时为 null）、`contig_count`、`total_length`、`duplicate_names`、`invalid_symbols`（contig -> 出现的无效符号）和 `errors`（消息列表）；输入无法使用时只有 `valid`、`format`、`contig_count` 和 `errors`。

- **示例代码：**

```python
import json

from dnakit.assembly import check_input

result = check_input(">c1\nACGT\n>c1\nAXGT\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "check_input",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 18,
    "sequence_bases": 8,
    "lengths": [4, 4],
    "sha256": "0681afaaa49db2047869730361f15d87a0c5b07c6027752a89404eba33894188"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "format": "sequences",
    "contig_count": 2,
    "total_length": 8,
    "duplicate_names": ["c1"],
    "invalid_symbols": {
      "c1": ["X"]
    },
    "errors": ["duplicate contig names: c1", "non-IUPAC symbols in 1 contig(s) (c1: X)"]
  }
}
```

## 16) compare_assemblies

- **作用：** 通过连续性统计量（数量、长度、N50、auN、GC、gap）以及共有的 contig 序列来比较两个组装。

- **处理方法：** 对每个组装计算 contig 数、总长度/最短/最长/平均长度、N50/L50 和 N90/L90（按长度从长到短累加，达到组装总长 50% / 90% 所需的长度）、auN = sum(L_i^2) / sum(L_i)（Nx 曲线下面积），给出 `genome_size` 时还有 NG50/LG50。有序列时再加上 GC 比例（G+C 占 A+C+G+T 的比例，排除 N 和其他符号）、N 碱基数、gap 数（N 的最长连续段）和其他符号数。然后报告 右减左 的差值、N50 之比，并且在两个组装都有序列时，报告在任一条链上序列完全相同的 contig（按规范序列集合比较）及共有碱基数。不做比对，因此相似但不完全相同的 contig 不会被匹配。

- **API：** `dnakit.assembly.compare_assemblies(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`[assembly_a, assembly_b]`、`{"left": a, "right": b}`（或 `{"a": ..., "b": ...}`），或单个组装加 `options["other"]`；每个组装是组装结果，可以是 FASTA 文本或 FASTA 文件路径、`{名称: 序列}`，或序列 / `{"name", "sequence"}` 记录的列表（也接受 contig 长度，形式为 `{名称: 长度}` 或整数列表）。`options`（可选）：`other`（第二个组装）、`genome_size`（正整数；加入 NG50 / LG50）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`left` 和 `right`（统计量：`contig_count`、`total_length`、`min_length`、`max_length`、`mean_length`、`n50`、`l50`、`n90`、`l90`、`aun`，可选的 `ng50` / `lg50`，有序列时还有 `gc_fraction`、`n_bases`、`gap_count`、`other_symbol_bases`）、`difference_right_minus_left`、`n50_ratio_right_to_left`，以及两侧都有序列时的 `shared_contig_sequences`、`shared_bases`、`left_only_contigs`、`right_only_contigs`。

- **示例代码：**

```python
import json

from dnakit.assembly import compare_assemblies

left = {"c1": "AAAA", "c2": "CCCCCC"}
right = {"x": "TTTT", "y": "GGGGGGG"}
result = compare_assemblies([left, right])
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "compare_assemblies",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 63,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7b5c153e3e35be332a4aa54b930ef40bbb2ada115430172b19f7f47de8429f13"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "left": {
      "contig_count": 2,
      "total_length": 10,
      "min_length": 4,
      "max_length": 6,
      "mean_length": 5.0,
      "n50": 6,
      "l50": 1,
      "n90": 4,
      "l90": 2,
      "aun": 5.2,
      "gc_fraction": 0.6,
      "n_bases": 0,
      "gap_count": 0,
      "other_symbol_bases": 0
    },
    "right": {
      "contig_count": 2,
      "total_length": 11,
      "min_length": 4,
      "max_length": 7,
      "mean_length": 5.5,
      "n50": 7,
      "l50": 1,
      "n90": 4,
      "l90": 2,
      "aun": 5.909090909090909,
      "gc_fraction": 0.6363636363636364,
      "n_bases": 0,
      "gap_count": 0,
      "other_symbol_bases": 0
    },
    "difference_right_minus_left": {
      "contig_count": 0,
      "total_length": 1,
      "n50": 1,
      "l50": 0,
      "max_length": 1
    },
    "n50_ratio_right_to_left": 1.1666666666666667,
    "shared_contig_sequences": 1,
    "shared_bases": 4,
    "left_only_contigs": 1,
    "right_only_contigs": 1
  }
}
```

## 17) compare_assembly_graphs

- **作用：** 按片段（segment）和连接（link）的内容比较两个 GFA1 组装图。

- **处理方法：** 读取两个 GFA1 图的 `S`（片段）行和 `L`（连接）行。当每个片段都有序列时，片段按规范序列（序列与其反向互补中较小者）匹配，否则按名称匹配。连接通过已匹配的片段来匹配，一个连接与其反向互补的对应连接视为同一个连接；连接的重叠（CIGAR）被忽略。报告各项计数、共有和独有的片段与连接、Jaccard 指数以及两个图是否相同。不解析 GFA2 以及路径/walk 行。

- **API：** `dnakit.assembly.compare_assembly_graphs(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：`[gfa_a, gfa_b]`（GFA1 文本或路径），或单个 GFA 加 `options["other"]`。`options`（可选）：`other`（第二个图）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`compared_by`（`"canonical segment sequence"` 或 `"segment name"`）、`left` 和 `right`（`segments`、`links`、`total_segment_length`）、`shared_segments`、`segment_jaccard`、`left_only_segments`、`right_only_segments`、`shared_links`、`link_jaccard` 和 `identical`。

- **示例代码：**

```python
import json

from dnakit.assembly import compare_assembly_graphs

gfa_a = "S\ta\tACGT\nS\tb\tAAAC\nL\ta\t+\tb\t+\t0M\n"
gfa_b = "S\tp\tACGT\nS\tq\tGTTT\nS\tr\tCCCC\nL\tp\t+\tq\t-\t0M\nL\tq\t+\tr\t+\t0M\n"
result = compare_assembly_graphs([gfa_a, gfa_b])
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "compare_assembly_graphs",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 125,
    "sequence_bases": 42,
    "lengths": [16, 26],
    "sha256": "e3c4b1e9641ea453179fce08375aa304f932df788456ec3bee4896872a1f7e64"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "compared_by": "canonical segment sequence",
    "left": {
      "segments": 2,
      "links": 1,
      "total_segment_length": 8
    },
    "right": {
      "segments": 3,
      "links": 2,
      "total_segment_length": 12
    },
    "shared_segments": 2,
    "segment_jaccard": 0.6666666666666666,
    "left_only_segments": 0,
    "right_only_segments": 1,
    "shared_links": 1,
    "link_jaccard": 0.5,
    "identical": false
  }
}
```

## 18) detect_assembly_contamination

- **作用：** 污染筛查尚未实现：它需要分类学、载体或接头数据库，或读段覆盖度。

- **处理方法：** 尚未实现。判断某个 contig 是否为外源序列，需要与参考数据库（例如 BLAST nt、Kraken2 或 FCS-GX）比较；GC 或长度异常并不等于污染判定。如果已知污染 contig 的名称，可用 `remove_contamination` 按名称移除。不会用启发式方法代替，因此函数不会编造结果；除了在 `input_summary` 中概括输入外，不对输入做任何校验，只返回结构化的 `not_implemented` 封装。

- **API：** `dnakit.assembly.detect_assembly_contamination(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：任意内容（只在 `input_summary` 中被概括）。`options`（可选）：被忽略。

- **输出：** 始终返回 `status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（缺少的东西：分类学或污染物参考数据库（例如 BLAST nt、Kraken2、FCS-GX））和 `message`（说明为何无法仅凭输入计算）。不产生 `result`。

- **示例代码：**

```python
import json

from dnakit.assembly import detect_assembly_contamination

result = detect_assembly_contamination(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "detect_assembly_contamination",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a taxonomic or contaminant reference database (for example BLAST nt, Kraken2, FCS-GX)",
  "message": "Deciding whether a contig is foreign requires comparison against reference databases; a GC or length outlier is not a contamination call."
}
```

## 19) detect_chimeric_contigs

- **作用：** 嵌合 contig 检测尚未实现：它需要读段或参考序列的比对。

- **处理方法：** 尚未实现。嵌合体是通过分裂的或不一致的读段比对，或与参考组装的比对来检测的，这些都不在输入中。不会用启发式方法代替，因此函数不会编造结果；除了在 `input_summary` 中概括输入外，不对输入做任何校验，只返回结构化的 `not_implemented` 封装。

- **API：** `dnakit.assembly.detect_chimeric_contigs(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：任意内容（只在 `input_summary` 中被概括）。`options`（可选）：被忽略。

- **输出：** 始终返回 `status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（缺少的东西：读段对 contig 的比对（BAM/PAF）或与参考组装的比对）和 `message`（说明为何无法仅凭输入计算）。不产生 `result`。

- **示例代码：**

```python
import json

from dnakit.assembly import detect_chimeric_contigs

result = detect_chimeric_contigs(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "detect_chimeric_contigs",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "read alignments (BAM/PAF) to the contigs or a reference assembly",
  "message": "Chimeras are detected from split or discordant alignments, which are not part of the input."
}
```

## 20) detect_misassemblies

- **作用：** 错误组装检测尚未实现：它需要与参考序列的比对或读段比对证据。

- **处理方法：** 尚未实现。错误组装是相对于参考基因组或读段证据定义的，二者都不在输入中。不会用启发式方法代替，因此函数不会编造结果；除了在 `input_summary` 中概括输入外，不对输入做任何校验，只返回结构化的 `not_implemented` 封装。

- **API：** `dnakit.assembly.detect_misassemblies(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：任意内容（只在 `input_summary` 中被概括）。`options`（可选）：被忽略。

- **输出：** 始终返回 `status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（缺少的东西：与参考基因组的比对（QUAST/MUMmer）或读段比对（CRAQ/Inspector））和 `message`（说明为何无法仅凭输入计算）。不产生 `result`。

- **示例代码：**

```python
import json

from dnakit.assembly import detect_misassemblies

result = detect_misassemblies(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "detect_misassemblies",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a reference genome alignment (QUAST/MUMmer) or read mappings (CRAQ/Inspector)",
  "message": "Misassemblies are defined relative to a reference or to read evidence, neither of which is part of the input."
}
```

## 21) estimate_heterozygosity

- **作用：** 用简单的双峰模型，从 k-mer 谱估计二倍体的 SNP 杂合率。

- **处理方法：** 启发式方法。设 `peak` 为纯合覆盖度峰，覆盖度在 `min_coverage` 与 `het_upper * peak` 之间的 k-mer 记为杂合（H，位于约一半覆盖度处），在 (`het_upper * peak`, `hom_upper * peak`] 之间的记为纯合（O）。对间隔大于 k 的孤立 SNP，H = 2*k*r*L，O = L*(1 - k*r)，令 R = H/O，则杂合率 `r = R / (k * (2 + R))`。该模型忽略重复序列、截断值以上的测序错误和成簇的 SNP，结果的好坏取决于截断值和峰位置；返回结果时会附带说明这一点的警告。如果同型窗口内没有 k-mer，会报错。

- **API：** `dnakit.assembly.estimate_heterozygosity(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：k-mer 直方图（`{覆盖度: 不同 k-mer 的数量}`、`[[覆盖度, 数量], ...]`，或 `jellyfish histo` 输出的 `覆盖度 数量` 文本/路径），或序列 / 读段（FASTA 文本、`{名称: 序列}` 或序列列表；统计其规范 k-mer，限制为 5,000,000 个碱基）。`options`（可选）：`k`（整数 3-127，默认 21；解释直方图需要它，对序列计数时也使用它）、`min_coverage`（错误截断值；默认取第一个“该覆盖度的数量不大于下一覆盖度的数量”的位置，即第一个局部谷）和 `peak_coverage`（默认取截断值之上最高的一格，若给定了 `min_coverage` 则包含截断值本身；并列时取覆盖度最低者）、`het_upper`（默认 0.75）、`hom_upper`（默认 1.5）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`heterozygosity`（每个碱基）、`het_distinct_kmers`、`hom_distinct_kmers`、`het_to_hom_ratio`、`k`、`peak_coverage`、`min_coverage` 以及说明字符串 `model`；封装的 `warnings` 中含有“heuristic two-peak model; see docstring for assumptions”。

- **示例代码：**

```python
import json

from dnakit.assembly import estimate_heterozygosity

histogram = {1: 500, 10: 20000, 20: 90000}
result = estimate_heterozygosity(histogram, options={"k": 10})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "estimate_heterozygosity",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 36,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4c5a162b6bc453f6cdb977dc5af239f06552601c4aadcdd2717baea800cbe65b"
  },
  "warnings": ["heuristic two-peak model; see docstring for assumptions"],
  "status": "ok",
  "result": {
    "heterozygosity": 0.01,
    "het_distinct_kmers": 20000,
    "hom_distinct_kmers": 90000,
    "het_to_hom_ratio": 0.2222222222222222,
    "k": 10,
    "peak_coverage": 20,
    "min_coverage": 2,
    "model": "diploid, isolated SNPs spaced > k apart: het k-mers = 2*k*r*L, hom k-mers = L*(1-k*r); r = R / (k*(2+R)) with R = het/hom"
  }
}
```

## 22) estimate_repeats

- **作用：** 根据 k-mer 重数估计重复序列含量：出现次数超过唯一序列覆盖度峰 `repeat_factor` 倍的 k-mer 视为重复。

- **处理方法：** 假设覆盖度均匀的启发式方法。设 T = 对所有 x >= `min_coverage` 求 x*n(x) 之和（n(x) 为出现 x 次的不同 k-mer 数），R = 对 x > `repeat_factor` * peak 求 x*n(x) 之和：`repeat_fraction = R/T`（处于多拷贝 k-mer 中的基因组位置所占比例），`repeat_positions_estimate = R/peak`，`genome_size_estimate = T/peak`。重复 k-mer 的估计拷贝数为 `max(2, round(x/peak))`，汇总在 `copy_number_classes` 中。不对覆盖度波动、测序偏差和截断值以上的错误建模；结果附带说明这一点的警告。

- **API：** `dnakit.assembly.estimate_repeats(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：k-mer 直方图（`{覆盖度: 不同 k-mer 的数量}`、`[[覆盖度, 数量], ...]`，或 `jellyfish histo` 输出的 `覆盖度 数量` 文本/路径），或序列 / 读段（FASTA 文本、`{名称: 序列}` 或序列列表；统计其规范 k-mer，限制为 5,000,000 个碱基）。`options`（可选）：`min_coverage`（错误截断值；默认取第一个“该覆盖度的数量不大于下一覆盖度的数量”的位置，即第一个局部谷）和 `peak_coverage`（默认取截断值之上最高的一格，若给定了 `min_coverage` 则包含截断值本身；并列时取覆盖度最低者）、`repeat_factor`（> 1，默认 1.5）、`k`（整数 3-63，默认 21；仅用于序列输入）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`repeat_fraction`、`repeat_positions_estimate`、`genome_size_estimate`、`repeat_distinct_kmers`、`copy_number_classes`（估计拷贝数 -> 不同 k-mer 数）、`peak_coverage`、`min_coverage`、`repeat_factor` 以及说明字符串 `definition`；`warnings` 中含有“k-mer multiplicity estimate; assumes uniform coverage”。

- **示例代码：**

```python
import json

from dnakit.assembly import estimate_repeats

histogram = {1: 100, 10: 1000, 20: 50, 30: 10}
result = estimate_repeats(histogram)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "estimate_repeats",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d923a1562f061521c7247589eee78b5f69b2a419717251a52cc75a7348a63cfc"
  },
  "warnings": ["k-mer multiplicity estimate; assumes uniform coverage"],
  "status": "ok",
  "result": {
    "repeat_fraction": 0.11504424778761062,
    "repeat_positions_estimate": 130.0,
    "genome_size_estimate": 1130.0,
    "repeat_distinct_kmers": 60,
    "copy_number_classes": {
      "2": 50,
      "3": 10
    },
    "peak_coverage": 10,
    "min_coverage": 2,
    "repeat_factor": 1.5,
    "definition": "repeat_fraction = sum_{x > repeat_factor*peak} x*n(x) / sum_{x >= min_coverage} x*n(x)"
  }
}
```

## 23) estimate_size

- **作用：** 根据 k-mer 直方图，以 k-mer 总个数除以主覆盖度峰来估计基因组大小。

- **处理方法：** Liu 等 2013 的模型：`G = 对所有 x >= min_coverage 求 x*n(x) 之和 / peak_coverage`，其中 n(x) 为出现 x 次的不同 k-mer 数（可用 `max_coverage` 限定上限）。低于截断值的 k-mer 视为测序错误。结果以碱基计，其好坏取决于截断值和峰位置，它们默认取第一个局部谷和其上最高的一格。

- **API：** `dnakit.assembly.estimate_size(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：k-mer 直方图（`{覆盖度: 不同 k-mer 的数量}`、`[[覆盖度, 数量], ...]`，或 `jellyfish histo` 输出的 `覆盖度 数量` 文本/路径），或序列 / 读段（FASTA 文本、`{名称: 序列}` 或序列列表；统计其规范 k-mer，限制为 5,000,000 个碱基）。`options`（可选）：`min_coverage`（错误截断值；默认取第一个“该覆盖度的数量不大于下一覆盖度的数量”的位置，即第一个局部谷）和 `peak_coverage`（默认取截断值之上最高的一格，若给定了 `min_coverage` 则包含截断值本身；并列时取覆盖度最低者）、`max_coverage`（可选的上限）、`k`（整数 3-63，默认 21；仅用于序列输入）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`genome_size`（碱基数）、`peak_coverage`、`min_coverage`、`max_coverage`、`total_kmer_instances`、`distinct_kmers_used` 以及说明字符串 `model`。

- **示例代码：**

```python
import json

from dnakit.assembly import estimate_size

histogram = {1: 500, 10: 20000, 20: 90000}
result = estimate_size(histogram)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "estimate_size",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 36,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4c5a162b6bc453f6cdb977dc5af239f06552601c4aadcdd2717baea800cbe65b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genome_size": 100000.0,
    "peak_coverage": 20,
    "min_coverage": 2,
    "max_coverage": null,
    "total_kmer_instances": 2000000,
    "distinct_kmers_used": 110000,
    "model": "G = sum_{x>=min_coverage} x * n(x) / peak_coverage"
  }
}
```

## 24) extract_accessory_pangenome

- **作用：** 列出泛基因组中的辅助基因，即存在于少于 `core_threshold` 比例基因组中的基因，并分为壳基因和云基因。

- **处理方法：** 基因至少存在于 `ceil(core_threshold * 基因组数)` 个基因组时为核心基因，其余基因均为辅助基因。存在比例低于 `cloud_threshold` 的辅助基因是云基因，其余是壳基因。基因名被视为已聚类好的同源组。

- **API：** `dnakit.assembly.extract_accessory_pangenome(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：基因有/无数据：`{基因组: [基因]}`、`{基因: {基因组: 0/1}}` 或 `{"genomes": [...], "genes": [...], "matrix": [[...]]}`（行为基因，列为基因组；取值为布尔值或非负数，> 0 表示存在）。基因名被视为已聚类好的同源组；不做序列聚类或同源搜索。`options`（可选）：`core_threshold`（(0, 1] 内的比例，默认 0.99）、`cloud_threshold`（[0, 1] 内的比例，默认 0.15）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`genome_count`、`total_genes`、`core_threshold`、`accessory_genes`、`accessory_count`、`shell_genes`、`cloud_genes`、`cloud_threshold` 以及 `gene_fractions`（携带每个辅助基因的基因组比例）。

- **示例代码：**

```python
import json

from dnakit.assembly import extract_accessory_pangenome

genomes = {"G1": ["a", "b", "c"], "G2": ["a", "b", "d"], "G3": ["a", "e"]}
result = extract_accessory_pangenome(genomes, options={"core_threshold": 0.66, "cloud_threshold": 0.5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "extract_accessory_pangenome",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 64,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e4372679c8d0e6836c36cd361cef3caf9656e4da871974b5c7cf4f7878ed448b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genome_count": 3,
    "total_genes": 5,
    "core_threshold": 0.66,
    "accessory_genes": ["c", "d", "e"],
    "accessory_count": 3,
    "shell_genes": [],
    "cloud_genes": ["c", "d", "e"],
    "cloud_threshold": 0.5,
    "gene_fractions": {
      "c": 0.3333333333333333,
      "d": 0.3333333333333333,
      "e": 0.3333333333333333
    }
  }
}
```

## 25) extract_core_pangenome

- **作用：** 列出泛基因组中的核心基因，即存在于至少 `core_threshold` 比例基因组中的基因。

- **处理方法：** 基因至少存在于 `ceil(core_threshold * 基因组数)` 个基因组时为核心基因（向上取整前使用 1e-9 的容差）。默认值 0.99 要求基因几乎存在于所有基因组中；严格的核心基因请用 1.0。基因名被视为已聚类好的同源组。

- **API：** `dnakit.assembly.extract_core_pangenome(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：基因有/无数据：`{基因组: [基因]}`、`{基因: {基因组: 0/1}}` 或 `{"genomes": [...], "genes": [...], "matrix": [[...]]}`（行为基因，列为基因组；取值为布尔值或非负数，> 0 表示存在）。基因名被视为已聚类好的同源组；不做序列聚类或同源搜索。`options`（可选）：`core_threshold`（(0, 1] 内的比例，默认 0.99）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`genome_count`、`total_genes`、`core_threshold`、`core_min_genomes`（所需的最少基因组数）、`core_genes`、`core_count` 以及 `gene_fractions`（携带每个核心基因的基因组比例）。

- **示例代码：**

```python
import json

from dnakit.assembly import extract_core_pangenome

genomes = {"G1": ["a", "b", "c"], "G2": ["a", "b", "d"], "G3": ["a", "e"]}
result = extract_core_pangenome(genomes, options={"core_threshold": 0.66})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "extract_core_pangenome",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 64,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e4372679c8d0e6836c36cd361cef3caf9656e4da871974b5c7cf4f7878ed448b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genome_count": 3,
    "total_genes": 5,
    "core_threshold": 0.66,
    "core_min_genomes": 2,
    "core_genes": ["a", "b"],
    "core_count": 2,
    "gene_fractions": {
      "a": 1.0,
      "b": 0.6666666666666666
    }
  }
}
```

## 26) identify_centromeres

- **作用：** 着丝粒识别尚未实现：不存在仅凭序列即可使用的规则，它需要卫星序列 / 重复序列注释。

- **处理方法：** 尚未实现。着丝粒位置因物种而异（卫星单体、CENH3 ChIP、Hi-C），不存在可直接用于裸 contig 的通用公认规则；它需要串联重复 / 卫星序列注释（TRF、SRF、RepeatMasker）以及物种特异的知识。不会用启发式方法代替，因此函数不会编造结果；除了在 `input_summary` 中概括输入外，不对输入做任何校验，只返回结构化的 `not_implemented` 封装。

- **API：** `dnakit.assembly.identify_centromeres(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：任意内容（只在 `input_summary` 中被概括）。`options`（可选）：被忽略。

- **输出：** 始终返回 `status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（缺少的东西：串联重复 / 卫星序列注释（TRF、SRF、RepeatMasker）以及物种特异的着丝粒重复序列知识）和 `message`（说明为何无法仅凭输入计算）。不产生 `result`。

- **示例代码：**

```python
import json

from dnakit.assembly import identify_centromeres

result = identify_centromeres(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "identify_centromeres",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "tandem-repeat/satellite annotation (TRF, SRF, RepeatMasker) and species-specific centromere repeat knowledge",
  "message": "Centromere positions are species-specific (satellite monomers, CENH3 ChIP, Hi-C); no universal documented rule exists to apply to bare contigs."
}
```

## 27) identify_repeats_in_assembly

- **作用：** 把组装中被至少出现 `min_copies` 次的 k-mer 覆盖的区间，作为重复区域找出来。

- **处理方法：** 在所有 contig 上合并统计规范 k-mer（k-mer 与其反向互补视为同一个）；含非 ACGT 碱基的 k-mer 被忽略。被计数 >= `min_copies` 的 k-mer 覆盖的每个位置都是重复位置；相互重叠或相邻的 k-mer 合并为 0 起始、半开区间，短于 `min_length` 的区间被丢弃。不使用重复序列库，也不做重复类型分类（转座子、卫星序列等），因此这纯粹是组装内部重复的度量；同一 contig 内和 contig 之间的拷贝都计入。上限为 5,000,000 个碱基。

- **API：** `dnakit.assembly.identify_repeats_in_assembly(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：含序列的组装（组装结果，可以是 FASTA 文本或 FASTA 文件路径、`{名称: 序列}`，或序列 / `{"name", "sequence"}` 记录的列表；只有长度不够）。`options`（可选）：`k`（整数 3-63，默认 21）、`min_copies`（整数 >= 2，默认 2）、`min_length`（整数，默认等于 `k`）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`k`、`min_copies`、`min_length`、`repeat_bases`、`total_bases`、`repeat_fraction`、`repeats`（contig -> `{start, end, length}` 列表）以及说明字符串 `definition`。

- **示例代码：**

```python
import json

from dnakit.assembly import identify_repeats_in_assembly

result = identify_repeats_in_assembly({"c": "GATTACACCGATTACA"}, options={"k": 5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "identify_repeats_in_assembly",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 25,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f3afd198c51e3955cb13540d57c5586de1dc62cd41af6c7fb5cc2900f62127a8"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "k": 5,
    "min_copies": 2,
    "min_length": 5,
    "repeat_bases": 14,
    "total_bases": 16,
    "repeat_fraction": 0.875,
    "repeats": {
      "c": [
        {
          "start": 0,
          "end": 7,
          "length": 7
        },
        {
          "start": 9,
          "end": 16,
          "length": 7
        }
      ]
    },
    "definition": "positions covered by canonical k-mers that occur >= min_copies times in the assembly; overlapping/adjacent k-mers merged; 0-based half-open intervals"
  }
}
```

## 28) identify_telomeres

- **作用：** 把 contig 末端由端粒基序完美串联重复构成的阵列识别为端粒重复阵列。

- **处理方法：** 端粒末端指由基序或其反向互补（二者的任意旋转，因此 `TTAGGG` 也能匹配 `CCCTAA` 和 `TAGGGT`）至少 `min_repeats` 个拷贝完美串联而成的一段，且它的起点距 contig 起点（左端）不超过 `max_distance` 个碱基，或终点距 contig 终点（右端）不超过 `max_distance` 个碱基；一段重复归属于较近的一端（并列时归左端），每端报告最长的合格重复段，坐标为 0 起始的半开区间。只能找到完美重复：简并的或较短的阵列会被漏掉，也不检查染色体内部的类端粒重复。

- **API：** `dnakit.assembly.identify_telomeres(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：含序列的组装（组装结果，可以是 FASTA 文本或 FASTA 文件路径、`{名称: 序列}`，或序列 / `{"name", "sequence"}` 记录的列表）。`options`（可选）：`motif`（2-12 个 ACGT 碱基，默认 `"TTAGGG"`，即脊椎动物重复；植物为 `"TTTAGGG"`）、`min_repeats`（整数，默认 5）、`max_distance`（整数，默认 100）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`motif`、`min_repeats`、`max_distance`、`contigs`（每个 contig 含 `length`、`left` 和 `right`（`{start, end, copies}` 或 null）以及 `telomeres_at_ends`）、`contigs_with_telomere`、`contigs_with_both_ends` 以及说明字符串 `definition`。

- **示例代码：**

```python
import json

from dnakit.assembly import identify_telomeres

contigs = {
    "left_tel": "TTAGGG" * 6 + "ACGTACGTACGTACGGATC",
    "right_tel": "GATCGATCGATCGATC" + "CCCTAA" * 5,
    "none": "ACGTACGTACGTACGTACGTACGT",
}
result = identify_telomeres(contigs)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "identify_telomeres",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 170,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f308faab612a4f4db560c3649114db4f013646b30b0661f9caf32b501fa001b9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "motif": "TTAGGG",
    "min_repeats": 5,
    "max_distance": 100,
    "contigs": {
      "left_tel": {
        "length": 55,
        "left": {
          "start": 0,
          "end": 36,
          "copies": 6
        },
        "right": null,
        "telomeres_at_ends": 1
      },
      "right_tel": {
        "length": 46,
        "left": null,
        "right": {
          "start": 16,
          "end": 46,
          "copies": 5
        },
        "telomeres_at_ends": 1
      },
      "none": {
        "length": 24,
        "left": null,
        "right": null,
        "telomeres_at_ends": 0
      }
    },
    "contigs_with_telomere": 2,
    "contigs_with_both_ends": 0,
    "definition": "a telomere end is a perfect tandem run of >= min_repeats copies of the motif or its reverse complement (any rotation) lying within max_distance bases of the nearer contig end (ties go to the left end)"
  }
}
```

## 29) parse_log

- **作用：** 按严重级别汇总组装程序的日志：ERROR / WARNING / INFO / DEBUG 行的数量，以及前几条错误和警告行。

- **处理方法：** 每一行取其中出现的第一个级别词（整词，不区分大小写）作为级别：CRITICAL 或 FATAL（计为 CRITICAL）、ERROR、WARNING 或 WARN（计为 WARNING）、INFO、DEBUG；不含这类词的行只计入 `line_count`。返回前 `max_lines` 条 ERROR / CRITICAL 行和前 `max_lines` 条 WARNING 行及其行号（文本截断为 500 个字符）。不解析特定工具（SPAdes、Flye 等）的字段。

- **API：** `dnakit.assembly.parse_log(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：日志文本或日志文件路径。`options`（可选）：`max_lines`（整数 1-1000，默认 20；返回多少条错误行和多少条警告行）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`line_count`、`level_counts`（`CRITICAL`、`ERROR`、`WARNING`、`INFO`、`DEBUG`）、`has_errors`（存在 ERROR 或 CRITICAL 行时为 true）、`errors` 和 `warnings`（`{line, text}` 列表）以及说明字符串 `definition`。

- **示例代码：**

```python
import json

from dnakit.assembly import parse_log

log = "INFO started\nwarning: low coverage\nERROR failed to open file\nINFO done\n[debug] details\n"
result = parse_log(log)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "parse_log",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 87,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8e4b46eec2308c01ef1b5cd36cb9c0a55cf35e4c545db1835b4313f5c8a4cd47"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "line_count": 5,
    "level_counts": {
      "CRITICAL": 0,
      "ERROR": 1,
      "WARNING": 1,
      "INFO": 2,
      "DEBUG": 1
    },
    "has_errors": true,
    "errors": [
      {
        "line": 3,
        "text": "ERROR failed to open file"
      }
    ],
    "warnings": [
      {
        "line": 2,
        "text": "warning: low coverage"
      }
    ],
    "definition": "a line's level is the first word ERROR/CRITICAL/FATAL/WARNING/WARN/INFO/DEBUG found in it (case-insensitive)"
  }
}
```

## 30) parse_output

- **作用：** 把组装程序的输出（FASTA contig 或 GFA1）解析为每条记录的长度、GC、覆盖度字段以及组装统计量。

- **处理方法：** 第一个非空行以 `>` 开头时格式为 FASTA，否则按 GFA1 处理。FASTA 记录给出长度和 GC 比例（G+C 占 A+C+G+T 的比例）；识别 SPAdes 的标题（`NODE_1_length_100_cov_5.2`）和 MEGAHIT 的标题（`flag=1 multi=2.0 len=300`），并返回其中声明的长度和覆盖度 / 重数。GFA 片段给出长度（来自序列或 `LN` 标签）和深度（来自首个出现的 `DP`、`dp` 或 `KC` 标签）。连续性统计量（N50 等）由这些长度计算。不解释其他组装程序的标题约定。

- **API：** `dnakit.assembly.parse_output(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：FASTA 或 GFA1 文本，或路径。`options`：无。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`format`（`"fasta"` 或 `"gfa"`）、`header_convention`（`"spades"`、`"megahit"` 或 null）、`records`（FASTA：`name`、`length`、`gc_fraction`，识别到标题约定时另有 `declared_length` 及 `coverage` 或 `multiplicity`；GFA：`name`、`length`、`depth`）、`stats`（与 `compare_assemblies` 相同的连续性统计量，由记录长度计算），GFA 另有 `link_count`。

- **示例代码：**

```python
import json

from dnakit.assembly import parse_output

fasta = ">NODE_1_length_8_cov_5.5\nACGTACGT\n>NODE_2_length_4_cov_2.0\nGGCC\n"
result = parse_output(fasta)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "parse_output",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 64,
    "sequence_bases": 12,
    "lengths": [8, 4],
    "sha256": "3d35bb543426e9cf6848f211f438b43042146cba719a22a02ae49eb460d310b5"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "fasta",
    "header_convention": "spades",
    "records": [
      {
        "name": "NODE_1_length_8_cov_5.5",
        "length": 8,
        "gc_fraction": 0.5,
        "declared_length": 8,
        "coverage": 5.5
      },
      {
        "name": "NODE_2_length_4_cov_2.0",
        "length": 4,
        "gc_fraction": 1.0,
        "declared_length": 4,
        "coverage": 2.0
      }
    ],
    "stats": {
      "contig_count": 2,
      "total_length": 12,
      "min_length": 4,
      "max_length": 8,
      "mean_length": 6.0,
      "n50": 8,
      "l50": 1,
      "n90": 4,
      "l90": 2,
      "aun": 6.666666666666667
    }
  }
}
```

## 31) phase_assembly

- **作用：** 对组装进行定相尚未实现：它需要读段、Hi-C 或亲本数据。

- **处理方法：** 尚未实现。把 contig 或碱基归入单倍型需要连接证据（长读段、Hi-C 或三联体 k-mer）和 hifiasm、HapCUT2、WhatsHap 之类的定相工具，这些都不在输入中。不会用启发式方法代替，因此函数不会编造结果；除了在 `input_summary` 中概括输入外，不对输入做任何校验，只返回结构化的 `not_implemented` 封装。

- **API：** `dnakit.assembly.phase_assembly(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：任意内容（只在 `input_summary` 中被概括）。`options`（可选）：被忽略。

- **输出：** 始终返回 `status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（缺少的东西：长读段、Hi-C 或三联体 k-mer，以及定相工具（hifiasm、HapCUT2、WhatsHap））和 `message`（说明为何无法仅凭输入计算）。不产生 `result`。

- **示例代码：**

```python
import json

from dnakit.assembly import phase_assembly

result = phase_assembly(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "phase_assembly",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "long reads, Hi-C or trio k-mers and a phasing assembler (hifiasm, HapCUT2, WhatsHap)",
  "message": "Assigning contigs or bases to haplotypes needs linking evidence that is not part of the input."
}
```

## 32) plot_graph

- **作用：** 把 GFA1 组装图画成 SVG：片段排布在圆周上，连接画成直线。

- **处理方法：** 片段按名称排序后，以相等的角度排布在圆周上（640 x 640 像素）；节点半径为 `4 + 2*log10(长度 + 1)` 像素；每个连接是其两个片段之间的直线（不画方向）；仅当片段不超过 60 个时才在节点旁写出名称，每个节点带有包含名称和长度的提示。布局是确定性的，没有生物学含义（不是 Bandage 式的布局）。最多 2,000 个片段；纯 Python 实现。给出 `output` 时，SVG 同时写入该文件。

- **API：** `dnakit.assembly.plot_graph(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：GFA1 文本或路径（最多 2,000 个片段）。`options`（可选）：`output`（写入 SVG 的文件路径，默认不写）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`svg`（SVG 文本）、`node_count`、`edge_count` 和 `output_path`（写入的文件，未写入时为 null）。

- **示例代码：**

```python
import json

from dnakit.assembly import plot_graph

gfa = "S\tp\tACGT\nS\tq\tGTTT\nS\tr\tCCCC\nL\tp\t+\tq\t-\t0M\nL\tq\t+\tr\t+\t0M\n"
result = plot_graph(gfa)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "plot_graph",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 53,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "32ffc4c7e7ccc5721dfef1ca795811be545710992f65d6e40a4a80807358f2e0"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"640\" height=\"640\" viewBox=\"0 0 640 640\">\n<line x1=\"600.0\" y1=\"320.0\" x2=\"180.0\" y2=\"562.5\" stroke=\"#888\" stroke-width=\"1\"/>\n<line x1=\"180.0\" y1=\"562.5\" x2=\"180.0\" y2=\"77.5\" stroke=\"#888\" stroke...",
    "node_count": 3,
    "edge_count": 2,
    "output_path": null
  }
}
```

（示例输出中的长数组和长文本已用 `...` 缩短。）

## 33) polish

- **作用：** 抛光（polish）尚未实现：它需要比对到组装上的读段，以及 Racon、Medaka、Pilon 之类的抛光工具。

- **处理方法：** 尚未实现。共识序列校正是由读段在组装上的堆叠（pileup）计算得到的，而读段不在输入中。不会用启发式方法代替，因此函数不会编造结果；除了在 `input_summary` 中概括输入外，不对输入做任何校验，只返回结构化的 `not_implemented` 封装。

- **API：** `dnakit.assembly.polish(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：任意内容（只在 `input_summary` 中被概括）。`options`（可选）：被忽略。

- **输出：** 始终返回 `status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（缺少的东西：测序读段、比对工具和抛光工具（Racon、Medaka、Pilon、NextPolish））和 `message`（说明为何无法仅凭输入计算）。不产生 `result`。

- **示例代码：**

```python
import json

from dnakit.assembly import polish

result = polish(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "polish",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "sequencing reads, an aligner and a polishing tool (Racon, Medaka, Pilon, NextPolish)",
  "message": "Consensus correction is computed from read pileups, which are not part of the input."
}
```

## 34) purge_haplotigs

- **作用：** 清除单倍型重复 contig（haplotig）尚未实现：它需要读段深度和自身比对证据（purge_dups、Purge_Haplotigs）。

- **处理方法：** 尚未实现。haplotig 是依据半覆盖度和等位基因间的自身比对证据识别的，这些证据不在输入中。不会用启发式方法代替，因此函数不会编造结果；除了在 `input_summary` 中概括输入外，不对输入做任何校验，只返回结构化的 `not_implemented` 封装。

- **API：** `dnakit.assembly.purge_haplotigs(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：任意内容（只在 `input_summary` 中被概括）。`options`（可选）：被忽略。

- **输出：** 始终返回 `status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（缺少的东西：读段覆盖度和 contig 自身比对（purge_dups、Purge_Haplotigs））和 `message`（说明为何无法仅凭输入计算）。不产生 `result`。

- **示例代码：**

```python
import json

from dnakit.assembly import purge_haplotigs

result = purge_haplotigs(">c1\nACGTACGT\n>c2\nGGCCGGCC\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "purge_haplotigs",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 26,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "eaa8633c7014be534dfd6571d98374f8f6d24defb95b394e7f9ad5bbe9ac56f4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "read coverage and contig self-alignments (purge_dups, Purge_Haplotigs)",
  "message": "Haplotigs are identified from half-coverage and allelic alignment evidence that is not part of the input."
}
```

## 35) remove_contamination

- **作用：** 从组装中移除 `options.contaminants` 所列名称的 contig（本函数不做污染检测）。

- **处理方法：** 丢弃名称在列表中的 contig，其余保持原顺序；组装中不存在的名称列在 `not_found` 中，不会报错。报告移除前后的连续性统计量，并在有序列时给出保留 contig 的 FASTA 文本。若会移除全部 contig 则报错。污染物名单必须来自别处完成的筛查（`detect_assembly_contamination` 尚未实现）。

- **API：** `dnakit.assembly.remove_contamination(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：组装结果，可以是 FASTA 文本或 FASTA 文件路径、`{名称: 序列}`，或序列 / `{"name", "sequence"}` 记录的列表（也接受 contig 长度，形式为 `{名称: 长度}` 或整数列表）。`options`（必需的键）：`contaminants`（要移除的 contig 名称列表）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`kept` 和 `removed`（contig 名称）、`not_found`、`stats_before`、`stats_after`（与 `compare_assemblies` 相同的连续性统计量），有序列时另有 `fasta`（保留的 contig）。

- **示例代码：**

```python
import json

from dnakit.assembly import remove_contamination

assembly = {"c1": "AAAA", "c2": "CCCCCC", "c3": "GG"}
result = remove_contamination(assembly, options={"contaminants": ["c2", "zz"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "remove_contamination",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5a945aa9f35559b3a7f890ec3868e4c370d32bbc5d5dbdce9aa4c71273e9bb14"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kept": ["c1", "c3"],
    "removed": ["c2"],
    "not_found": ["zz"],
    "stats_before": {
      "contig_count": 3,
      "total_length": 12,
      "min_length": 2,
      "max_length": 6,
      "mean_length": 4.0,
      "n50": 6,
      "l50": 1,
      "n90": 2,
      "l90": 3,
      "aun": 4.666666666666667,
      "gc_fraction": 0.6666666666666666,
      "n_bases": 0,
      "gap_count": 0,
      "other_symbol_bases": 0
    },
    "stats_after": {
      "contig_count": 2,
      "total_length": 6,
      "min_length": 2,
      "max_length": 4,
      "mean_length": 3.0,
      "n50": 4,
      "l50": 1,
      "n90": 2,
      "l90": 2,
      "aun": 3.3333333333333335,
      "gc_fraction": 0.3333333333333333,
      "n_bases": 0,
      "gap_count": 0,
      "other_symbol_bases": 0
    },
    "fasta": ">c1\nAAAA\n>c3\nGG\n"
  }
}
```

## 36) run

- **作用：** 通过对应的 `assemble_*` 封装，运行调用者自行安装的某个组装程序（Flye、hifiasm、MEGAHIT、Shasta 或 Unicycler）。

- **处理方法：** 先按 `build_command` 完全相同的规则校验设置（因此即使程序不存在，缺少必需的键也会报错），再在 PATH 中查找可执行文件。找不到时，函数返回指明该程序的 `not_implemented`，且不运行任何东西。找到时，用相同的参数向量调用对应的封装（`assemble_flye`、`assemble_hifiasm`、`assemble_megahit`、`assemble_shasta` 或 `assemble_unicycler`）：不经过 shell，只在 `output_directory` 下输出，带有超时（`timeout_seconds`）和输出大小限制；其结果放在 `result` 中返回。DNAKit 不会安装或下载任何组装程序。下面的示例是在没有安装任何组装程序的环境中生成的，因此展示的是 `not_implemented` 分支（成功分支在单元测试中用替身可执行文件验证过，没有用真实的组装程序验证）。

- **API：** `dnakit.assembly.run(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：设置映射（`options` 会合并到其上并优先生效）：`tool`（`flye`、`hifiasm`、`megahit`、`shasta` 或 `unicycler`）、`output_directory`，以及与 `build_command` 完全相同的工具输入与参数（`reads`、`read1`、`read2`、`long_reads`、`genome_size`、`platform`、`config`、`mode`、`prefix`、`threads`、`min_contig_length`、`program`），另有可选的 `timeout_seconds`。

- **输出：** PATH 中没有该可执行文件时：`status: "not_implemented"`（并带有 `module`、`operation`、`provenance`、`input_summary`、`warnings`）、`required_dependency`（`"<program> executable on PATH"`）和 `message`（该程序未安装或不在 PATH 中，因此什么都没有运行）；没有 `result`。有可执行文件时：标准封装中 `status` 为 `"ok"`，`result` 为封装程序的报告（`status`、形如 `assemble_flye` 的 `operation`、`backend`、含输出目录、组装文件、`artifacts` 和程序 `log` 的 `result`，以及 `provenance`）。

- **示例代码：**

```python
import json

from dnakit.assembly import run

settings = {
    "tool": "flye",
    "reads": "reads.fastq",
    "genome_size": "5m",
    "output_directory": "flye_out",
}
result = run(settings)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "run",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 93,
    "sequence_bases": 10,
    "lengths": [10],
    "sha256": "858ba72be547f61d29e0e8be6d87531c90bc5c6b4be8918a4ca7c84b978670fd"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "flye executable on PATH",
  "message": "The assembler program 'flye' is not installed or not on PATH, so nothing was run."
}
```

## 37) scaffold

- **作用：** 按调用者提供的显式布局把 contig 连成 scaffold（不会推断任何连接证据）。

- **处理方法：** 对 `layout` 中的每个 scaffold，按顺序拼接所列 contig；方向为 `-` 的 contig 先取反向互补，相邻 contig 之间插入指定长度的 `N`（每个 contig 的 `gap`，否则用全局 `gap`；gap 为 0 则不插入）。每个 contig 只能使用一次。不会预测顺序和方向：布局必须来自其他证据（例如 Hi-C 脚手架工具）。类 AGP 表使用每个 scaffold 上 0 起始的半开区间坐标。

- **API：** `dnakit.assembly.scaffold(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`（必填）：含序列的组装（组装结果，可以是 FASTA 文本或 FASTA 文件路径、`{名称: 序列}`，或序列 / `{"name", "sequence"}` 记录的列表）。`options`：`layout`（必需）：scaffold 列表 `{"name": 字符串, "contigs": [contig 名称 | {"contig": 名称, "orientation": "+" | "-", "gap": 整数}, ...]}`，其中 `gap` 是该 contig 与下一个 contig 之间的 N 数；`gap`（整数 >= 0，默认 100）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`。`result` 包含：`scaffolds`（每项含 `name`、`sequence`、`length`、`contig_count`）、`agp`（每行含 `scaffold`、`start`、`end`、`type`（`contig` 或 `gap`）、`contig`、`orientation`）和 `unplaced_contigs`（布局中未使用的 contig）。

- **示例代码：**

```python
import json

from dnakit.assembly import scaffold

contigs = {"c1": "AACC", "c2": "GGTT", "c3": "ACGT"}
layout = [{"name": "S1", "contigs": ["c1", {"contig": "c2", "orientation": "-"}]}]
result = scaffold(contigs, options={"layout": layout, "gap": 5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "assembly",
  "operation": "scaffold",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6646b5dc42423ea3aa5dd2d01c9954cd36aa9d53434aa8fc7dafd4747ac80aa4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "scaffolds": [
      {
        "name": "S1",
        "sequence": "AACCNNNNNAACC",
        "length": 13,
        "contig_count": 2
      }
    ],
    "agp": [
      {
        "scaffold": "S1",
        "start": 0,
        "end": 4,
        "type": "contig",
        "contig": "c1",
        "orientation": "+"
      },
      {
        "scaffold": "S1",
        "start": 4,
        "end": 9,
        "type": "gap",
        "contig": null,
        "orientation": null
      },
      {
        "scaffold": "S1",
        "start": 9,
        "end": 13,
        "type": "contig",
        "contig": "c2",
        "orientation": "-"
      }
    ],
    "unplaced_contigs": ["c3"]
  }
}
```
