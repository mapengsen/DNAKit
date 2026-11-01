# Shasta 长读长组装（命令行封装）

[Shasta](https://github.com/paoloshasta/shasta) 用于 Oxford Nanopore 长读长的从头组装。DNAKit 以固定参数白名单（不经 shell）调用 `shasta` 命令行，并读回 Shasta 为各流程阶段写出的文件。一个函数覆盖命令行能选择的 `Assembler` API 阶段变体：由 anchor 做 Mode 3 分型组装、标记（k-mer）选取、MinHash 比对候选搜索、回文与嵌合读长标记、比对方法（旧 Shasta、SeqAn、带状 SeqAn、方法 4）、读长图构建方法（标准、百分位、v3、贝叶斯 v4）及跨链边标记、不一致比对标记、标记图顶点与边创建、重复标记清理、传递约简、剪枝与简化、重复数共识、解缠，以及 GFA 与 FASTA 导出。

Shasta 需另行安装（发布版二进制 https://github.com/paoloshasta/shasta/releases，已用 0.14.0 测试），放到 PATH 中并命名为 `shasta`，或设置 `DNAKIT_SHASTA`，或传入 `shasta_path`。Shasta 的 Python API（`import shasta`）未封装：它属于源码构建，方法只能作用于上一阶段留在内存中的状态，且发布版共享库需要 Python 3.10 与 Ubuntu 22.04 系统库。只以此类方法存在的阶段（如 `writeMarkers`、`analyzeReadGraph`、`readGraphClustering`、GFA 着色）因此不提供，详见集成包 P29 的处置表。

## 示例数据

`sim_reads.py` 模拟 150 kb 随机基因组及双链长读长（错误率 1%，函数 `write_reads`）。示例 1 用 450 条读长，示例 3 用 800 条、错误率 0.2%。

**示例文件：** `sim_reads.py`

```text
import random


def write_reads(path, count=450, error=0.01, seed=1):
    """Simulate a 150 kb random genome and long reads (12-25 kb, both strands) with substitutions and indels."""
    rng = random.Random(seed)
    genome = "".join(rng.choice("ACGT") for _ in range(150_000))
    table = str.maketrans("ACGT", "TGCA")
    with open(path, "w") as out:
        for index in range(count):
            length = rng.randint(12_000, 25_000)
            start = rng.randint(0, len(genome) - length)
            read = genome[start:start + length]
            if rng.random() < 0.5:
                read = read.translate(table)[::-1]
            mutated = []
            for base in read:
                x = rng.random()
                if x < error / 3:
                    continue
                if x < 2 * error / 3:
                    mutated.append(base + rng.choice("ACGT"))
                elif x < error:
                    mutated.append(rng.choice("ACGT"))
                else:
                    mutated.append(base)
            out.write(f">r{index}\n{''.join(mutated)}\n")
```

## 1) `dnakit.assembly.shasta_assemble()`

- **作用：** 用 Shasta 组装长读长；流程的每个阶段由参数选择（`assembly_mode`、`align_method`、`read_graph_creation_method`、`detangle_method` 等），其他 Shasta 0.14.0 选项可放入 `options`。
- **处理方法：** 以参数白名单加经校验的 `options` 运行 `shasta --command assemble`；参数仅在给出时才传，其余由所选 `config` 决定。运行目录必须不存在；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.assembly.shasta_assemble(reads: 'Sequence[str]', output_directory: 'str', *, config: 'str' = 'Nanopore-May2022', assembly_mode: 'int | None' = None, anchors_path: 'str | None' = None, min_read_length: 'int | None' = None, desired_coverage: 'str | int | None' = None, palindromic_flagging: 'bool | None' = None, kmer_generation_method: 'int | None' = None, kmer_length: 'int | None' = None, kmer_probability: 'float | None' = None, minhash_m: 'int | None' = None, minhash_hash_fraction: 'float | None' = None, minhash_iteration_count: 'int | None' = None, minhash_all_pairs: 'bool | None' = None, align_method: 'int | None' = None, read_graph_creation_method: 'int | None' = None, strand_separation_method: 'int | None' = None, max_chimeric_read_distance: 'int | None' = None, flag_inconsistent_alignments: 'bool | None' = None, marker_graph_min_coverage: 'int | None' = None, marker_graph_max_coverage: 'int | None' = None, marker_graph_min_edge_coverage: 'int | None' = None, cleanup_duplicate_markers: 'bool | None' = None, simplify_max_length: 'str | Sequence[int] | None' = None, detangle_method: 'int | None' = None, consensus_caller: 'str | None' = None, write_reads_by_segment: 'bool | None' = None, store_coverage_data: 'bool | None' = None, memory_mode: 'str | None' = None, memory_backing: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, threads: 'int | None' = None, shasta_path: 'str | None' = None, timeout_seconds: 'float' = 86400.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 读长文件（FASTA/FASTQ）、新的输出目录、配置名或文件，以及可选的阶段参数。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`assembly_directory`、`assembly_fasta`、`assembly_gfa`、`summary`（使用的读长数、回文与嵌合读长、组装片段、碱基数、N50）、`outputs`、`log_tail`。

## 2) `dnakit.assembly.shasta_read_run_table()`

- **作用：** 读取 Shasta 写入运行目录的某个表（读长标记、MinHash 与比对直方图、读长图连通分量、标记图覆盖度直方图、各片段读长、Mode 3 表）。
- **处理方法：** 解析运行目录中的 CSV 文件；数值自动转换，最多返回 `max_rows` 行。
- **API：** `dnakit.assembly.shasta_read_run_table(run_directory: 'str', table: 'str', *, max_rows: 'int' = 10000, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 已完成组装的运行目录和表名。
- **输出：** 字典，含 `table`、`file`、`columns`、`row_count`、`rows`、`truncated`。

## 3) `dnakit.assembly.shasta_read_assembly_output()`

- **作用：** 读取 Shasta 导出的组装 FASTA 或 GFA 图并汇总（记录数、长度、N50、GC、连接）。
- **处理方法：** 解析 FASTA 或 GFA 1.0（序列或 `LN:i` 长度）；不运行任何程序。
- **API：** `dnakit.assembly.shasta_read_assembly_output(run_directory: 'str', *, output: 'str' = 'fasta', file_name: 'str | None' = None, include_sequences: 'bool' = False, max_records: 'int' = 10000, verbose: bool = True) -> 'dict[str, Any]'`。
- **输入：** 运行目录、输出类型（`fasta`、`gfa`、`gfa_both_strands`、`gfa_no_sequence`）或目录内的普通 `file_name`。
- **输出：** 字典，含 `records`、`record_count`、`total_length`、`longest`、`n50`，GFA 另有 `link_count`/`links`。

## 示例

### 示例 1： 组装读长并读回运行结果

组装模拟读长，再读取读长表、读长图连通分量以及 FASTA 与 GFA 导出。

- **示例代码：**

```python
import tempfile
from pathlib import Path
from sim_reads import write_reads
from dnakit.assembly import shasta_p29_cli_functions as sh

work = Path(tempfile.mkdtemp())
write_reads(work / "reads.fasta")
result = sh.shasta_assemble([str(work / "reads.fasta")], str(work / "run"), min_read_length=5000, threads=4, verbose=False)
summary = result["summary"]
print(result["status"], result["version"])
print(summary["reads_used"], summary["reads_flagged_chimeric"], summary["assembled_segments"])

reads = sh.shasta_read_run_table(result["assembly_directory"], "read_summary", verbose=False)
print(reads["columns"][:5], reads["row_count"])
components = sh.shasta_read_run_table(result["assembly_directory"], "read_graph_components", max_rows=2, verbose=False)
first = components["rows"][0]
print(components["row_count"], first["OrientedReadCount"], first["IsSelfComplementary"])
fasta = sh.shasta_read_assembly_output(result["assembly_directory"], verbose=False)
print(fasta["record_count"], fasta["total_length"] == summary["assembled_bases"], fasta["n50"] == summary["segment_n50"])
gfa = sh.shasta_read_assembly_output(result["assembly_directory"], output="gfa_both_strands", verbose=False)
print(gfa["record_count"], gfa["link_count"])
```

- **示例输出：**

```text
completed 0.14.0
450 5 3
['Id', 'Name', 'RawLength', 'RleLength', 'RawOverRleLengthRatio'] 450
8 108 No
3 True True
6 0
```

### 示例 2： 用参数选择流程阶段

比对方法、读长图构建方法、解缠、各片段读长、Mode 2，以及被拒绝的选项名。

- **示例代码：**

```python
import tempfile
from pathlib import Path
from sim_reads import write_reads
from dnakit.assembly import shasta_p29_cli_functions as sh

work = Path(tempfile.mkdtemp())
write_reads(work / "reads.fasta")
runs = {
    "seqan-alignment+percentile-graph": {"align_method": 1, "read_graph_creation_method": 2},
    "bayesian-graph": {"read_graph_creation_method": 4, "strand_separation_method": 0},
    "detangle+reads-by-segment": {"detangle_method": 2, "write_reads_by_segment": True},
    "phased-mode-2": {"assembly_mode": 2, "strand_separation_method": 2},
}
for name, parameters in runs.items():
    result = sh.shasta_assemble([str(work / "reads.fasta")], str(work / name), min_read_length=5000, threads=4, verbose=False, **parameters)
    print(name, result["status"])
by_segment = sh.shasta_read_run_table(str(work / "detangle+reads-by-segment"), "reads_by_segment", max_rows=1, verbose=False)
print(by_segment["row_count"] > 0, by_segment["columns"][:4])
try:
    sh.shasta_assemble([str(work / "reads.fasta")], str(work / "bad"), options={"Align.noSuchOption": 1}, verbose=False)
except Exception as error:
    print(type(error).__name__)
```

- **示例输出：**

```text
seqan-alignment+percentile-graph completed
bayesian-graph completed
detangle+reads-by-segment completed
phased-mode-2 completed
True ['AssembledSegmentId', 'EdgeCount', 'OrientedReadCount', 'OrientedReadId']
ConfigurationError
```

### 示例 3： 由 anchor 做 Mode 3 分型组装

用原始序列配置 `Nanopore-ncm23-May2024` 做 Mode 3 组装并读取 Mode 3 表。

- **示例代码：**

```python
import tempfile
from pathlib import Path
from sim_reads import write_reads
from dnakit.assembly import shasta_p29_cli_functions as sh

work = Path(tempfile.mkdtemp())
write_reads(work / "reads.fasta", count=800, error=0.002)
result = sh.shasta_assemble(
    [str(work / "reads.fasta")], str(work / "mode3"), config="Nanopore-ncm23-May2024", assembly_mode=3,
    min_read_length=5000, threads=4, verbose=False,
)
print(result["return_code"], result["version"])
print(sorted(Path(result["assembly_directory"]).glob("Anchor*"))[0].name)
anchors = sh.shasta_read_run_table(result["assembly_directory"], "anchor_coverage_histogram", verbose=False)
print(anchors["columns"], anchors["row_count"] > 0)
print(sh.shasta_read_run_table(result["assembly_directory"], "components", max_rows=1, verbose=False)["columns"][:3])
```

- **示例输出：**

```text
0 0.14.0
AnchorCoverageHistogram.csv
['Coverage', 'Frequency'] True
['Component', 'Reads', 'Segments']
```
