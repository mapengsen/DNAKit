# 可选 DNA 模型包的数据工具

这一组函数是 `dnakit.predictions` 中基于可选依赖的适配器，覆盖 [alphagenome](https://github.com/google-deepmind/alphagenome)、[grelu](https://github.com/Genentech/gReLU)、[selene-sdk](https://github.com/FunctionLab/selene)、[enformer-pytorch](https://github.com/lucidrains/enformer-pytorch)、[genomic-benchmarks](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks)、[evo2](https://github.com/ArcInstitute/evo2) 和 [kipoi](https://github.com/kipoi/kipoi) 的 DNA 层面数据操作。内容包括区间和变异解析、GENCODE 基因与转录本提取、序列截取与变换、one-hot 编码、BED 窗口、基准数据集摘要，以及本地 Evo 2 生成和 kipoi 模型描述。这些函数不做模型推理，也不下载权重。模型推理和托管预测接口已经由 `dnakit.predictions` 与 `dnakit.databases` 提供，这里不重复。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。

## 1) `dnakit.predictions.parse_alphagenome_interval()`

- **作用：** 解析 `chr:start-end[:strand]` 形式的区间字符串。
- **处理方法：** 调用 alphagenome 的 `Interval.from_str`。起点按原样使用，视为 0-based，终点是开区间，不做 1-based 转换。
- **API：** `dnakit.predictions.parse_alphagenome_interval(text)`。
- **输入：** 区间字符串，例如 `"chr1:101-200:+"`。
- **输出：** `interval` 字段，包含 `chromosome`、`start`、`end`、`strand` 和 `width`。`provenance` 中记录坐标约定。
- **示例代码：**

```python
from dnakit.predictions import parse_alphagenome_interval

result = parse_alphagenome_interval("chr1:101-200:+")
print(result["interval"])
```

- **示例输出：**

```text
{'chromosome': 'chr1', 'start': 101, 'end': 200, 'strand': '+', 'width': 99}
```

## 2) `dnakit.predictions.parse_alphagenome_variant()`

- **作用：** 解析 `chr:position:REF>ALT` 形式的变异字符串，并给出变异类型。
- **处理方法：** 调用 alphagenome 的 `Variant.from_str`。位置按 VCF 的 1-based 理解，`start` 为 `position - 1`。`end` 是 alphagenome 按变异类型计算的开区间终点，例如 `chr1:10:AT>A` 的终点为 11。`variant_format` 可选 `default`、`gtex`、`open_targets`、`open_targets_bigquery` 或 `gnomad`。
- **API：** `dnakit.predictions.parse_alphagenome_variant(text, *, variant_format="default")`。
- **输入：** 变异字符串，以及可选的 `variant_format`。
- **输出：** `variant` 字段，包含 `chromosome`、`position`、`reference`、`alternate`、`start`、`end`，以及 `is_snv`、`is_indel`、`is_insertion`、`is_deletion` 和 `is_structural`。
- **示例代码：**

```python
from dnakit.predictions import parse_alphagenome_variant

variant = parse_alphagenome_variant("chr1:10:A>G")["variant"]
print(variant["start"], variant["end"], variant["is_snv"])
```

- **示例输出：**

```text
9 10 True
```

## 3) `dnakit.predictions.merge_alphagenome_intervals()`

- **作用：** 合并同一条染色体上有重叠的区间。
- **处理方法：** 先用 `Interval.from_str` 解析每个字符串，再调用 alphagenome 的 `merge_overlapping_intervals`。
- **API：** `dnakit.predictions.merge_alphagenome_intervals(intervals)`。
- **输入：** 区间字符串列表，最多 10,000 条。
- **输出：** `input_count`、`merged_count` 和合并后的 `intervals` 列表。
- **示例代码：**

```python
from dnakit.predictions import merge_alphagenome_intervals

result = merge_alphagenome_intervals(["chr1:0-10", "chr1:5-20", "chr1:30-40"])
print([(item["start"], item["end"]) for item in result["intervals"]])
```

- **示例输出：**

```text
[(0, 20), (30, 40)]
```

## 4) `dnakit.predictions.read_gencode_gene_intervals()`

- **作用：** 从 GENCODE GTF 中取出指定基因的区间。
- **处理方法：** 本模块先用纯 Python 解析 GTF 的 `gene` 行，把 1-based 闭区间的起点减 1，转为 0-based 半开区间，然后把表交给 alphagenome 的 `get_gene_intervals`。`gene_symbols` 和 `gene_ids` 必须二选一。
- **API：** `dnakit.predictions.read_gencode_gene_intervals(gtf_path, *, gene_symbols=None, gene_ids=None, max_records=200000)`。
- **输入：** GTF 文件路径，可以是 `.gz`；基因符号列表或 Ensembl 基因 ID 列表；`max_records` 是选中行数的上限。
- **输出：** `genes` 列表，每项包含区间坐标、`strand` 和 `name`。
- **示例文件：** `genes.gtf`

```text
chr1	test	gene	101	500	.	+	.	gene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";
chr1	test	gene	700	800	.	-	.	gene_id "G2"; gene_name "GENE2"; gene_type "lincRNA";
```

- **示例代码：**

```python
from dnakit.predictions import read_gencode_gene_intervals

result = read_gencode_gene_intervals("genes.gtf", gene_symbols=["GENE1"])
print(result["genes"])
```

- **示例输出：**

```text
[{'chromosome': 'chr1', 'start': 100, 'end': 500, 'strand': '+', 'width': 400, 'name': 'GENE1'}]
```

## 5) `dnakit.predictions.extract_gencode_transcripts()`

- **作用：** 返回与指定区间重叠的转录本，包括外显子和 CDS。
- **处理方法：** 本模块用纯 Python 解析 GTF 中的 `transcript`、`exon`、`CDS`、`start_codon` 和 `stop_codon` 行，坐标转换规则与上一个函数相同。随后调用 alphagenome 的 `TranscriptExtractor`。CDS 区间由 alphagenome 的 `Transcript` 对象给出，可能与原始 GTF 行不同。
- **API：** `dnakit.predictions.extract_gencode_transcripts(gtf_path, region, *, max_records=200000)`。
- **输入：** GTF 文件路径，以及 0-based 半开区间字符串，例如 `"chr1:0-1000"`。
- **输出：** `transcripts` 列表，每项包含 `transcript_id`、`gene_id`、`strand`、`exons`、`cds` 和 `is_coding`。坐标均为 0-based 半开区间。
- **示例文件：** `tx.gtf`

```text
chr1	test	transcript	101	500	.	+	.	gene_id "G1"; transcript_id "T1"; transcript_type "protein_coding";
chr1	test	exon	101	201	.	+	.	gene_id "G1"; transcript_id "T1"; transcript_type "protein_coding";
chr1	test	exon	301	500	.	+	.	gene_id "G1"; transcript_id "T1"; transcript_type "protein_coding";
chr1	test	CDS	151	198	.	+	0	gene_id "G1"; transcript_id "T1"; transcript_type "protein_coding";
chr1	test	stop_codon	199	201	.	+	0	gene_id "G1"; transcript_id "T1"; transcript_type "protein_coding";
```

- **示例代码：**

```python
from dnakit.predictions import extract_gencode_transcripts

result = extract_gencode_transcripts("tx.gtf", "chr1:0-1000")
for item in result["transcripts"]:
    print(item["transcript_id"], item["exons"], item["cds"])
```

- **示例输出：**

```text
T1 [[100, 201], [300, 500]] [[150, 198]]
```

## 6) `dnakit.predictions.reverse_complement_with_grelu()`

- **作用：** 对一条或多条 DNA 序列求反向互补。
- **处理方法：** 调用 grelu 的 `reverse_complement`。只接受 `A`、`C`、`G`、`T` 和 `N`，输入会转为大写。
- **API：** `dnakit.predictions.reverse_complement_with_grelu(sequences)`。
- **输入：** 一条 DNA 字符串，或字符串列表，最多 10,000 条。
- **输出：** `sequences` 列表，即使输入是单条字符串也返回列表。
- **示例代码：**

```python
from dnakit.predictions import reverse_complement_with_grelu

print(reverse_complement_with_grelu(["AAC", "GGTN"])["sequences"])
```

- **示例输出：**

```text
['GTT', 'NACC']
```

## 7) `dnakit.predictions.shuffle_dinucleotides()`

- **作用：** 生成保持二核苷酸计数的随机打乱序列。
- **处理方法：** 调用 grelu 的 `dinuc_shuffle`，底层由 tangermeme 完成。若请求多条打乱结果且它们全部相同，上游会报错，本函数把这种情况转为 `ConfigurationError`。
- **API：** `dnakit.predictions.shuffle_dinucleotides(sequence, *, n_shuffles=1, seed=None)`。
- **输入：** 一条 DNA 序列；`n_shuffles` 为 1 到 1000 的整数；`seed` 可选，范围为 0 到 2^32 - 1。
- **输出：** `shuffles` 列表，每条的长度与输入相同。
- **示例代码：**

```python
from dnakit.predictions import shuffle_dinucleotides

result = shuffle_dinucleotides("AACCGGTTAGCTAGGATCCAGTAC", n_shuffles=2, seed=7)
print(len(result["shuffles"]), result["shuffles"][0][:8])
```

- **示例输出：**

```text
2 AGGTAGGA
```

## 8) `dnakit.predictions.resize_dna_sequences()`

- **作用：** 把 DNA 序列裁剪或用 N 补齐到指定长度。
- **处理方法：** 调用 grelu 的 `resize`。`end` 可选 `left`、`right` 或 `both`，补的碱基是 `N`。
- **API：** `dnakit.predictions.resize_dna_sequences(sequences, seq_len, *, end="both")`。
- **输入：** 一条 DNA 字符串或字符串列表；目标长度 `seq_len`；`end` 指定在哪一端操作。
- **输出：** `sequences` 列表，`seq_len` 和 `end` 也会回显。
- **示例代码：**

```python
from dnakit.predictions import resize_dna_sequences

print(resize_dna_sequences(["ACGT", "ACGTAC"], 6)["sequences"])
```

- **示例输出：**

```text
['NACGTN', 'ACGTAC']
```

## 9) `dnakit.predictions.mutate_dna_sequence()`

- **作用：** 把序列中一个位置的碱基替换为另一个碱基。
- **处理方法：** 调用 grelu 的 `mutate`。`position` 是 0-based 下标。函数返回原碱基，便于核对。
- **API：** `dnakit.predictions.mutate_dna_sequence(sequence, allele, *, position)`。
- **输入：** DNA 序列；单个碱基 `allele`；0-based 的 `position`，必须在序列范围内。
- **输出：** `sequence`（突变后序列）、`reference_base`、`position` 和 `allele`。
- **示例代码：**

```python
from dnakit.predictions import mutate_dna_sequence

result = mutate_dna_sequence("ACGT", "T", position=1)
print(result["sequence"], result["reference_base"])
```

- **示例输出：**

```text
ATGT C
```

## 10) `dnakit.predictions.fetch_genome_sequences()`

- **作用：** 从 FASTA 中取出指定区间的序列。
- **处理方法：** 调用 selene-sdk 的 `Genome`，并用 `get_sequence_from_coords` 取序列。区间为 0-based 半开区间。`pad=False` 时，越界的区间会报错，因为 selene 对越界区间返回空串。`pad=True` 时用 `N` 补齐。pyfaidx 会在 FASTA 旁边写一个 `.fai` 索引文件。
- **API：** `dnakit.predictions.fetch_genome_sequences(fasta_path, regions, *, pad=False, max_regions=1000, max_region_length=100000)`。
- **输入：** FASTA 路径；区域列表，每项是包含 `chrom`、`start`、`end` 和可选 `strand` 的字典。
- **输出：** `regions` 列表，每项包含 `chrom`、`start`、`end`、`strand` 和 `sequence`。
- **示例文件：** `genome.fa`

```text
>chr1
ACGTACGTNN
>chr2
GGGCCCAATT
```

- **示例代码：**

```python
from dnakit.predictions import fetch_genome_sequences

result = fetch_genome_sequences(
    "genome.fa",
    [{"chrom": "chr1", "start": 0, "end": 4}, {"chrom": "chr1", "start": 8, "end": 12}],
    pad=True,
)
print([item["sequence"] for item in result["regions"]])
```

- **示例输出：**

```text
['ACGT', 'NNNN']
```

## 11) `dnakit.predictions.read_contig_lengths()`

- **作用：** 列出 FASTA 中每条序列的名称和长度。
- **处理方法：** 用 selene-sdk 的 `Genome` 读取索引，并调用 `get_chr_lens`。
- **API：** `dnakit.predictions.read_contig_lengths(fasta_path, *, max_contigs=100000)`。
- **输入：** FASTA 路径。
- **输出：** `contigs` 列表，每项包含 `name` 和 `length`，另有 `contig_count`。
- **示例代码：**

```python
from dnakit.predictions import read_contig_lengths

print(read_contig_lengths("genome.fa")["contigs"])
```

- **示例输出：**

```text
[{'name': 'chr1', 'length': 10}, {'name': 'chr2', 'length': 10}]
```

## 12) `dnakit.predictions.encode_one_hot_selene()`

- **作用：** 用 selene-sdk 的编码对 DNA 做 one-hot 编码。
- **处理方法：** 调用 `Genome.sequence_to_encoding`。通道顺序为 A、C、G、T。`N` 在每个通道都是 0.25，这与 enformer-pytorch 不同。
- **API：** `dnakit.predictions.encode_one_hot_selene(sequence, *, max_length=10000)`。
- **输入：** DNA 序列，长度不超过 `max_length`，上限为 100,000。
- **输出：** `one_hot` 矩阵（每个碱基一行，4 列），以及 `length` 和 `channels`。
- **示例代码：**

```python
from dnakit.predictions import encode_one_hot_selene

result = encode_one_hot_selene("ACGTN")
print(result["one_hot"][4])
```

- **示例输出：**

```text
[0.25, 0.25, 0.25, 0.25]
```

## 13) `dnakit.predictions.encode_one_hot_enformer()`

- **作用：** 用 enformer-pytorch 的编码对 DNA 做 one-hot 编码。
- **处理方法：** 调用 `enformer_pytorch.data.str_to_one_hot`。通道顺序为 A、C、G、T。`N` 在所有通道都是 0，这与 selene-sdk 不同。
- **API：** `dnakit.predictions.encode_one_hot_enformer(sequence, *, max_length=10000)`。
- **输入：** DNA 序列，长度不超过 `max_length`，上限为 100,000。
- **输出：** `one_hot` 矩阵（每个碱基一行，4 列），以及 `length` 和 `channels`。
- **示例代码：**

```python
from dnakit.predictions import encode_one_hot_enformer

result = encode_one_hot_enformer("ACGTN")
print(result["one_hot"][4])
```

- **示例输出：**

```text
[0.0, 0.0, 0.0, 0.0]
```

## 14) `dnakit.predictions.read_enformer_bed_windows()`

- **作用：** 对 BED 区间生成 enformer-pytorch 的输入窗口，并统计每个窗口的填充和未知碱基数。
- **处理方法：** 调用 `enformer_pytorch.data.GenomeIntervalDataset`。短于 `context_length` 的区间会在两端对称扩展，扩展部分的填充位为 `.`，编码后每个通道为 0.25。长于 `context_length` 的区间保持原长度。函数只返回计数，不返回完整的 one-hot 矩阵。
- **API：** `dnakit.predictions.read_enformer_bed_windows(bed_path, fasta_path, *, context_length=196608, max_windows=100)`。
- **输入：** BED 文件（至少三列，制表符分隔）；FASTA 路径；`context_length` 与 `max_windows` 都有上限。
- **输出：** `windows` 列表，每项包含 `chromosome`、`start`、`end`、`length`、`padding_positions`、`unknown_positions` 和 `acgt_positions`。
- **示例文件：** `windows.bed`

```text
chr1	0	4
chr1	6	10
```

- **示例代码：**

```python
from dnakit.predictions import read_enformer_bed_windows

result = read_enformer_bed_windows("windows.bed", "genome.fa", context_length=8)
for window in result["windows"]:
    print(window["chromosome"], window["length"], window["padding_positions"], window["unknown_positions"])
```

- **示例输出：**

```text
chr1 8 2 0
chr1 8 2 2
```

## 15) `dnakit.predictions.list_genomic_benchmark_datasets()`

- **作用：** 列出 genomic-benchmarks 的数据集名称。
- **处理方法：** 不传 `dataset_path` 时，读取包内的数据集索引，不访问网络。传入目录时，列出该目录下的子目录。
- **API：** `dnakit.predictions.list_genomic_benchmark_datasets(*, dataset_path=None)`。
- **输入：** 可选的本地目录。
- **输出：** `source`（`package_index` 或 `local_directory`）、`count` 和排序后的 `datasets` 列表。
- **示例代码：**

```python
from dnakit.predictions import list_genomic_benchmark_datasets

result = list_genomic_benchmark_datasets()
print(result["count"], result["datasets"][:3])
```

- **示例输出：**

```text
9 ['demo_coding_vs_intergenomic_seqs', 'demo_human_or_worm', 'drosophila_enhancers_stark']
```

## 16) `dnakit.predictions.describe_genomic_benchmark()`

- **作用：** 统计本地 genomic-benchmarks 数据集中每个类别、每个划分的序列数。
- **处理方法：** 调用 `genomic_benchmarks.data_check.info`，并设置 `local_repo=True`，不访问网络。上游打印的摘要行被捕获到 `messages` 中，不输出到标准输出。
- **API：** `dnakit.predictions.describe_genomic_benchmark(dataset_path, *, version=None)`。
- **输入：** 数据集目录，目录中需要有 `metadata.yaml`、`train/` 和 `test/`；`version` 可选，与元数据不一致时报错。
- **输出：** `classes`（类别到 `train`、`test` 计数的映射）、`total_sequences` 和 `messages`。
- **示例代码：**

```python
import gzip
from pathlib import Path

from dnakit.predictions import describe_genomic_benchmark

root = Path("demo_ds")
for split in ("train", "test"):
    (root / split).mkdir(parents=True, exist_ok=True)
    for label in ("neg", "pos"):
        with gzip.open(root / split / f"{label}.csv.gz", "wt") as handle:
            handle.write("id,region,start,end\n1,chr1,10,20\n2,chr1,30,45\n")
(root / "metadata.yaml").write_text("version: 0\nclasses:\n  - neg\n  - pos\n")

result = describe_genomic_benchmark(root, version=0)
print(result["classes"], result["total_sequences"])
```

- **示例输出：**

```text
{'neg': {'train': 2, 'test': 2}, 'pos': {'train': 2, 'test': 2}} 8
```

## 17) `dnakit.predictions.check_genomic_benchmark_cache()`

- **作用：** 检查某个 genomic-benchmarks 数据集是否已存在于本地缓存目录。
- **处理方法：** 调用 `genomic_benchmarks.data_check.is_downloaded`，只检查目录是否存在，不下载。未指定 `cache_path` 时使用包的默认缓存目录 `~/.genomic_benchmarks`。
- **API：** `dnakit.predictions.check_genomic_benchmark_cache(dataset, *, cache_path=None)`。
- **输入：** 数据集名称，只允许字母、数字、下划线和短横线；可选的缓存目录。
- **输出：** `downloaded`（布尔值）、`dataset` 和 `cache_path`。
- **示例代码：**

```python
from dnakit.predictions import check_genomic_benchmark_cache

result = check_genomic_benchmark_cache("demo_ds", cache_path="empty_cache")
print(result["downloaded"])
```

- **示例输出：**

```text
False
```

## 18) `dnakit.predictions.generate_evo2_sequences()`

- **作用：** 用本地 Evo 2 权重对 DNA 前缀做续写生成。
- **处理方法：** 构造 `evo2.Evo2(model_name=..., local_path=...)`，再调用 `generate`。`weights_path` 必须是已存在的本地文件，函数不会触发 Hugging Face 下载。`model_name` 只能取 evo2 支持的名称之一。
- **API：** `dnakit.predictions.generate_evo2_sequences(prompts, *, weights_path, model_name="evo2_7b", n_tokens=100, temperature=1.0, top_k=4, top_p=1.0, max_prompts=8, max_prompt_length=10000)`。
- **输入：** DNA 前缀列表；本地权重文件；采样参数。
- **输出：** `sequences` 和 `scores`。`scores` 是上游返回的每条序列分数。
- **示例代码（未本地运行）：**

```python
from dnakit.predictions import generate_evo2_sequences

result = generate_evo2_sequences(
    ["ACGTACGTAC"],
    weights_path="/path/to/evo2_7b.pt",
    n_tokens=32,
)
print(result["sequences"])
```

- **示例输出（未本地运行）：**

```text
（未本地运行，无输出）
```

## 19) `dnakit.predictions.list_kipoi_models()`

- **作用：** 列出 kipoi 模型源中的模型名称、版本、类型和输入输出名称。
- **处理方法：** 调用 `kipoi.list_models`。kipoi 会读取基于 git 的模型源，首次使用时可能克隆仓库。不指定 `sources` 时，会读取所有已配置的源。
- **API：** `dnakit.predictions.list_kipoi_models(*, sources=None, limit=1000)`。
- **输入：** 可选的源名称列表；`limit` 限制返回的行数。
- **输出：** `total_models`、`returned` 和 `models` 列表。
- **示例代码（未本地运行）：**

```python
from dnakit.predictions import list_kipoi_models

result = list_kipoi_models(sources=["kipoi"], limit=5)
print(result["total_models"])
```

- **示例输出（未本地运行）：**

```text
（未本地运行，无输出）
```

## 20) `dnakit.predictions.describe_kipoi_model()`

- **作用：** 返回一个 kipoi 模型的类型、版本、文档、许可证、引用和输入输出名称。
- **处理方法：** 调用 `kipoi.get_model_descr`。输入输出的名称兼容列表、字典和单个数组对象三种写法。
- **API：** `dnakit.predictions.describe_kipoi_model(model, *, source="kipoi")`。
- **输入：** 模型名称，例如 `"DeepSEA/variantEffects"`；源名称，默认为 `"kipoi"`。
- **输出：** `description` 字典，包含 `type`、`version`、`doc`、`license`、`cite_as`、`inputs` 和 `targets`。
- **示例代码（未本地运行）：**

```python
from dnakit.predictions import describe_kipoi_model

result = describe_kipoi_model("DeepSEA/variantEffects")
print(result["description"]["inputs"])
```

- **示例输出（未本地运行）：**

```text
（未本地运行，无输出）
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-alphagenome]"          # alphagenome：区间、变异与 GENCODE 工具
python -m pip install "dnakit[bio-grelu]"                # grelu：DNA 字符串工具
python -m pip install "dnakit[bio-selene]"               # selene-sdk：FASTA 序列与 one-hot
python -m pip install "dnakit[bio-enformer]"             # enformer-pytorch：one-hot 与 BED 窗口
python -m pip install "dnakit[bio-genomic-benchmarks]"   # genomic-benchmarks：数据集摘要
python -m pip install "dnakit[bio-evo2]"                 # evo2：本地生成，需要 GPU 与权重
python -m pip install "dnakit[bio-kipoi]"                # kipoi：模型描述，需要网络
```

项目地址：[alphagenome](https://github.com/google-deepmind/alphagenome)、[grelu](https://github.com/Genentech/gReLU)、[selene-sdk](https://github.com/FunctionLab/selene)、[enformer-pytorch](https://github.com/lucidrains/enformer-pytorch)、[genomic-benchmarks](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks)、[evo2](https://github.com/ArcInstitute/evo2)、[kipoi](https://github.com/kipoi/kipoi)。

## 验证方式

测试文件是 `others/tests/unit/test_dna_model_backends.py`，共 43 项。其中 42 项通过，1 项是 kipoi 实时测试，默认跳过，需要设置 `DNAKIT_LIVE_TESTS=1`。主要检查如下：

- 17 个函数直接调用真实上游包：alphagenome 5 个，grelu 4 个，selene-sdk 3 个，enformer-pytorch 2 个，genomic-benchmarks 3 个。手工计算的结果（例如 GTF 坐标从 1-based 转为 0-based、N 碱基的编码、窗口的填充计数）与直接调用上游的结果作对照。
- `generate_evo2_sequences`、`list_kipoi_models` 和 `describe_kipoi_model` 只用伪模块（注入 `sys.modules`）测试，检查参数、返回结构和错误映射。这些测试不代表上游推理结果。evo2 在本地因缺少 `flash_attn_2_cuda` 无法导入，另有一个在真实导入失败时检查安装提示的测试。
- 7 个可选包分别检查：缺少包时抛出 `BackendUnavailableError`，并给出对应的 extra 名称。
- 输入检查会拒绝非法 DNA 字符、布尔值、越界位置、超出上限的长度、空字符串和 NUL 字符。

## 当前范围

这一批覆盖 7 个包中的 DNA 数据操作，共 20 个函数。以下内容未包装：

- alphagenome 的托管预测接口和本地 JAX 推理：`dnakit.databases` 与 `dnakit.predictions` 已有对应实现。
- evo2 的序列打分（`score_sequences`）：`dnakit.predictions` 的 Evo 2 变异效应后端已在内部使用，这里只包装本地生成。
- enformer-pytorch 的 `from_pretrained` 与轨道预测：已由 `dnakit.predictions` 覆盖。
- selene-sdk 的模型训练、计算机模拟诱变和模型预测：需要训练好的模型，本批不包装。
- grelu 的模型训练、解释和变异效应预测：需要模型，本批不包装。
- genomic-benchmarks 的 `loc2seq`、`seq2loc` 和数据集下载：需要联网下载基因组，本批不包装。
- kipoi 的模型下载、依赖安装和模型运行：需要联网和模型权重，本批不包装。
