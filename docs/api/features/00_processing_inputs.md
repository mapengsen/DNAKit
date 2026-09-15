# 数据处理的统一输入与输出

标准化、序列操作、去重、聚类和数据划分的公共入口支持 DNA 文本、文件路径和既有 DNA 对象。
文件按后缀识别 FASTA/FNA、TXT、CSV/TSV 等 DNAKit I/O 格式，也支持 gzip。
原有 DNA 对象调用保留其返回类型、坐标约定、随机参数和算法限制。

```python
from dnakit.ops import reverse_complement, trim, transcribe, translate
from dnakit.standardize import normalize

if __name__ == "__main__":
    print(reverse_complement("AATCGC"))  # GCGATT
    normalize("input.txt")
    reverse_complement("hg38.fna")
    trim("input.csv", left=2, right=2, sequence_column="sequence")
    transcribe("input.fna", output_path="rna.fa")
    translate("input.csv", output_path="proteins.csv")
```

## 输入规则

- DNA 文本按单条序列处理。反向互补直接返回字符串；其他功能沿用其原生结果类型。
- `str` 中包含点号或路径分隔符时通常视为文件路径；无后缀或容易与 DNA 混淆的文件请使用 `Path`。
  标准化仍允许原始文本中包含标点，只有明确的序列文件后缀才自动切换到文件模式。
- TXT 默认每个非空行是一条序列；整个文件是一条换行排版序列时，使用 `ReadConfig(text_mode="sequence")`。
- CSV 默认使用 `sequence` 列，可通过 `sequence_column="dna"` 更换列名。
  标签及附加列进入记录 metadata，并在 CSV 输出时按原列名、原顺序写回。CSV 单元格保留字符串语义。
  同名 metadata 值与附加列冲突时会报错。
- 未压缩文件默认根据实际大小设置读取上限，hg38 无需手写大小配置；gzip 默认有 20 GB 的解压后读取上限。
  显式传入 `read_config=ReadConfig(...)` 时遵守该配置。
- 环状操作须显式使用 `topology="circular"`，或在 `ReadConfig`/DNA 对象中声明环状拓扑。

基础入口 `read_inputs()` 也可用于自定义流程，并负责关闭文件：

```python
from dnakit.io import read_inputs

with read_inputs(["first.fna", "second.csv", "AATCGC"]) as records:
    for record in records:
        print(record.id, record.sequence.symbol_length)
```

## 文件输出规则

逐条操作默认写到同目录的 `原名_操作名.原后缀`，也可指定 `output_path`，返回 `WriteResult`。
输出路径位于 `result.target_artifact.relative_path`，记录数位于 `result.record_count`。
处理逐条读取，默认显示开始、记录进度和完成提示，可用 `show_progress=False` 关闭。

原始输入不能作为输出，已有结果默认不覆盖。重新生成可传 `write_config=WriteConfig(overwrite=True)`。
文件写入失败时不安装未完成的结果。坐标、随机种子、N/IUPAC 接受范围和功能所需的附加依赖仍遵循原算法。

转录和翻译分别写出 RNA、蛋白质文本，支持 FASTA、TXT、CSV/TSV 和 JSONL；FASTA 来源默认输出 `.fa`。
自动派生的 CSV 输出会移除不适用于 RNA/蛋白质的 DNA 结构列。
原始序列的复杂 feature/逐碱基注释若无法可靠映射，会明确报错；记录感知的编辑、反向互补、环状变换和切分会使用相应的注释映射。

## 覆盖的公共入口

| 类别 | 入口与处理单位 |
|---|---|
| 标准化 | `standardize.normalize`、`normalize_gaps`；前者可先清理 FASTA/FNA、TXT、CSV/TSV 的非法原始字符 |
| 校验 | `standardize.validate` 逐条输出 JSONL 报告；`validate_set` 输出数据集报告 |
| 方向 | `ops.reverse`、`complement`、`reverse_complement` |
| 转录与翻译 | `ops.transcribe`、`translate` |
| 提取与编辑 | `ops.subsequence`、`circular_subsequence`、`insert`、`delete`、`substitute`、`mask`、`trim` |
| 生成 | `ops.mutate`、`evolution_generate`、`indel_generate`、`rearrange_generate`、`kmer_shuffle` |
| 拼接 | `ops.concat`、`concat_overlap` 按文件及记录顺序拼接；`crossover` 按两个输入的记录顺序配对，数量不一致时报错 |
| 环状操作 | `ops.rotate`、`canonical_origin` |
| 记录操作 | `ops` 导出的编辑、反向互补及环状操作的 `*_record` 入口 |
| 切分 | 新增 `chunking.chunk_sequences`，原有 `iter_sequence_chunks`、`iter_fasta_chunks` 继续提供惰性迭代 |
| 去重 | `datasets.deduplicate`、`deduplicate_iupac`、`deduplicate_approximate` |
| 聚类 | `datasets.cluster_sequences`、`hierarchical_cluster`、`neural_cluster_sequences`、`select_representatives` |
| 整理与划分 | `datasets.exclude_by_metadata`、`exclude_chromosomes`、`exclude_species`、`split`、`temporal_split`、`joint_split` |
| 数据集评估 | `datasets.detect_leakage` 接受集合名到输入文件的映射；`evaluate_split_quality` 接受文件和原有 assignments 参数 |

AGP 装配、BED 区域说明等结构化辅助输入仍使用各自的专用接口，不会当作 DNA 文本猜测。

## 去重、聚类和划分

集合操作接受单个文件或文件列表。原始记录 ID 会保留，组合文件时应保证算法要求的 ID 唯一性。
这些算法按数据集计算，需要的内存和记录上限沿用原实现。
除随机文件划分外，文件模式返回原生结果对象，并将报告及数据写到 `原名_操作名/`，也可以通过 `output_path` 指定目录。
目录中的 `result.json` 保存结果，CSV 输入的子集/代表序列保留原表列，其余采用可保留注释的 JSONL。
整个目录写完后才安装；`overwrite=True` 会替换该结果目录，因此应使用专用输出目录。

```python
from dnakit.datasets import SplitConfig, deduplicate, split

if __name__ == "__main__":
    dedup = deduplicate("sequences.csv")
    result = split(
        "sequences.csv",
        config=SplitConfig(ratios={"train": 0.8, "valid": 0.1, "test": 0.1}, seed=42),
        output_path="splits",
    )
    # splits/train.csv、valid.csv、test.csv、result.json
```

切分文件可以直接写出，无需手写循环：

```python
from dnakit.chunking import ChunkingConfig, chunk_sequences

result = chunk_sequences(
    "hg38.fna", config=ChunkingConfig(length=200), output_path="windows.fna"
)
```


## 多核心处理

这些公共处理入口统一接受 `workers=-1`（默认，全部可用逻辑核心）、`workers=1`（串行）或正整数
（最多使用指定数量的工作进程）。核心数遵守进程 CPU 亲和性；实际并发还受任务数量和可用内存限制。

```python
from dnakit.ops import reverse_complement

if __name__ == "__main__":
    result = reverse_complement("hg38.fna", workers=8)
```

普通 Python 脚本中的处理调用应放在 `if __name__ == "__main__":` 内，避免子进程启动时重复执行脚本。
交互式环境或不可序列化的自定义回调可使用 `workers=1`。
论文示例 `02_reverse_complement.py`、`03_split_dataset.py` 在脚本内直接设置 `workers`。
命令行入口则使用全局选项：

```bash
dnakit --workers 8 deduplicate input.fna output.fna
```

CLI 的 `--workers` 是全局选项，放在子命令之前；它控制已接入的 CPU 处理步骤。

| 计算步骤 | 执行方式 |
|---|---|
| 大型、多记录、未压缩 FASTA 的序列操作、标准化和校验 | 默认读取配置下按完整记录区段分工，工作进程自行读写，主进程顺序合并；自定义读取、压缩或输出大小限制使用通用路径 |
| 内存连续序列的反转、互补、反向互补 | 使用 C 字符串操作，避免进程间复制；大型 FASTA 反向互补另有按碱基区段处理的路径 |
| 精确/反向互补去重 | 本地计算轻量等价键，避免来回复制序列 |
| 环状去重 | 并行计算等价键，顺序合并分组 |
| IUPAC 去重 | 并行计算位置兼容性，保留原有稳定分组规则 |
| 阈值聚类、层次聚类、近似去重、代表序列选择 | k-mer 集合或二进制指纹每条序列只计算一次；其他方法仍并行计算成对相似度，顺序执行分组与合并 |
| 相似度划分、联合划分 | 并行计算指纹或成对相似度，顺序分配集合 |
| 两个文件的 crossover | 按行配对并行，保持输出顺序 |
| 拼接、内存对象的惰性窗口切分、简单元数据筛选、随机/时间分配本身 | 保留流式或顺序处理，避免为轻量或有顺序依赖的步骤强行启动进程 |

神经模型聚类和模型推理仍由原后端管理 CPU/GPU；此参数不会自动复制 GPU 模型。
不承诺所有任务随核心数增加而加速。短任务自动串行；等待处理的数据量有上限，超大记录不会无限进入进程队列。
多进程输出按输入顺序合并，随机种子仍按原算法使用；调度顺序不参与随机数生成。
显式传入可变 `rng` 对象时保留串行调用及其状态推进；自定义分块进度回调也留在调用进程。
单核文件操作使用直接流式路径，避免不必要的分片合并。大型 FASTA 的方向操作仍使用专用区段路径；
512 MiB 以内的此类文件最多使用 8 个工作进程，以减少启动和磁盘争用。

## 随机文件划分的内存使用

`split(文件或文件列表, config=SplitConfig(method="random", ...))` 对默认读取配置的未压缩
FASTA 优先索引原文件，保留 ID、偏移量和划分信息，再按需读取、并行写出子集，省去序列暂存副本。
默认 JSONL 导出会进一步按记录分组，train 内部也能分给多个进程。常规 FASTA 预先计算 JSONL
记录位置，各进程直接写入互不重叠的区段，避免再次复制整份输出；不规则布局自动回退到通用写出。
其他格式、自定义读取配置或不适合索引的布局仍使用磁盘暂存路径。
它不会一次保存整份基因组的序列对象；每个工作进程仍需容纳单条最大记录及转换缓冲区。
最终输出目录在全部子集成功写出后安装，临时数据在成功或失败后清理。
其他划分方法及原生对象调用保持原有算法和内存语义。

随机文件划分返回 `FileSplitResult`：`counts`、`assignments`、`seed` 和 `to_dict()` 保持原划分审计内容，
`get("train")` 和 `subsets[i].records` 返回可重复迭代的 `FileRecords`，按需读取已写出的文件。
如确实需要一次加载整个子集，可使用 `DNASet(result.get("train"))`。
输出比例、记录顺序和 `seed=42` 对应的分配与内存版一致，`preserve_order=False` 也保留原排序语义。

随机文件划分默认显示读取及各子集写出的百分比、数据量、速度、已耗时和 ETA。
读取按源文件字节计量，写出按已完成记录对应的数据量计量；ETA 是当前阶段的估算，
启动时尚无速度样本及文件校验阶段显示未知时间，不表示整个任务的剩余时间。

已有的下载/解压、序列处理、模型预测、表征提取、评估和工作流进度条统一显示 `elapsed`（已耗时）与 `ETA`（预计剩余时间）。
有总量时根据实际完成速度估算 ETA；未知总量的流式处理和外部命令不预读数据、不猜测总时长，ETA 显示未知。
下载切换到解压或另一个文件时重置速度样本。原有 `verbose=False`、`show_progress=False` 和用户自定义回调继续生效，脚本无需重复指定时间参数。


### 损坏压缩输入

共享读取层将截断 gzip 的 EOFError 转换为 `InputFormatError(code="INPUT_READ_FAILED")`，
保留文件来源和底层原因。标准化的原始文本读取同样处理压缩、解码及 CSV 解析错误。
使用这些读取入口的文件操作继续遵守暂存后安装的输出策略，失败不替换已有完整输出。
