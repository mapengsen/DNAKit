# Primer3 设计任务、SPAdes 工具集与 WhatsHap 定相

三组封装用户自行安装的工具的函数：[Primer3](https://primer3.org/) 的设计任务、比对与重复序列屏蔽（`dnakit.molecular_lab`）；[SPAdes](https://github.com/ablab/spades) 工具集中的读段纠错、组装图与图比对程序（`dnakit.assembly`）；[WhatsHap](https://whatshap.readthedocs.io/) 的定相、基因分型与 haplotag（`dnakit.variants`）。安装：`conda install -c bioconda primer3 spades whatshap`（primer3-py 对应额外依赖 `bio-primer3`）。

命令行函数用类型化参数构造固定的参数列表来运行一个程序或子命令：不经过 shell，没有通用执行器，输出都是新文件，只传入你设置的选项，其余选项保持程序自身默认值。找不到程序时抛出带安装提示的 `BackendUnavailableError`；`primer3_path`、`spades_path`、`whatshap_path` 用来指定程序所在目录（WhatsHap 为可执行文件本身），代替 PATH。返回字典包含 `status`（退出码为 0 且所有预期输出存在时才为 `completed`）、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 以及少量解析结果。测试版本：primer3-py 2.3.1（libprimer3 2.6.1）、Primer3 2.6.1、SPAdes 4.3.0、WhatsHap 2.8。

寡核苷酸的解链温度、发夹和二聚体热力学已由 `calculate_primer3_tm`、`calculate_primer3_hairpin`、`calculate_primer3_dimer` 提供；不带家系的 `whatshap phase` 即 `phase_variants_whatshap`。

## 示例数据

准备示例为 Primer3 写入 600 bp 随机模板；为 SPAdes 生成含 400 bp 重复的 18.8 kb 基因组、5000 对 100 bp 双端 reads（30% 含一处替换）、4 条 2 kb 长读段和 200 bp 查询序列；为 WhatsHap 生成 4 kb 参考序列、25 个杂合 SNV、从两条单倍型交替抽取的 500 条 reads（已排序并建索引的 BAM）和 VCF，并用 `whatshap phase` 先定相一次，使统计、比较和 haplotag 示例有定相后的输入。后面的示例在同一目录中使用这些文件。

- **示例代码：**

```python
import random
import subprocess

import pysam

rng = random.Random(2)


def dna(length):
    return "".join(rng.choice("ACGT") for _ in range(length))


open("template.txt", "w").write(dna(600))
# SPAdes data
rc = str.maketrans("ACGT", "TGCA")
repeat = dna(400)
genome = dna(6000) + repeat + dna(6000) + repeat + dna(6000)
open("genome.fa", "w").write(">g\n" + genome + "\n")
first, second = [], []
for index in range(5000):
    start = rng.randrange(0, len(genome) - 300)
    fragment = genome[start:start + 300]
    read1, read2 = fragment[:100], fragment[::-1].translate(rc)[:100]
    if rng.random() < 0.3:
        position = rng.randrange(100)
        read1 = read1[:position] + rng.choice("ACGT") + read1[position + 1:]
    first.append(f"@p{index}/1\n{read1}\n+\n{'I' * 100}\n")
    second.append(f"@p{index}/2\n{read2}\n+\n{'I' * 100}\n")
open("r1.fq", "w").write("".join(first))
open("r2.fq", "w").write("".join(second))
open("long.fa", "w").write("".join(f">l{i}\n{genome[start:start + 2000]}\n" for i, start in enumerate((100, 5000, 12000, 16000))))
open("query.fa", "w").write(">q\n" + genome[3000:3200] + "\n")
# WhatsHap data
length = 4000
reference = dna(length)
open("ref.fa", "w").write(">chr1\n" + "\n".join(reference[i:i + 60] for i in range(0, length, 60)) + "\n")
pysam.faidx("ref.fa")
flip = {"A": "G", "C": "T", "G": "A", "T": "C"}
sites = list(range(100, length - 200, 150))
alternative = list(reference)
for position in sites:
    alternative[position] = flip[reference[position]]
haplotypes = [reference, "".join(alternative)]
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": length}], "RG": [{"ID": "rg1", "SM": "S1"}]}
with pysam.AlignmentFile("raw.bam", "wb", header=header) as out:
    for index in range(500):
        start = rng.randrange(0, length - 300)
        read = pysam.AlignedSegment(out.header)
        read.query_name = f"r{index}"
        read.query_sequence = haplotypes[index % 2][start:start + 300]
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigar = [(0, 300)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 300)
        read.set_tag("RG", "rg1")
        out.write(read)
pysam.sort("-o", "reads.bam", "raw.bam")
pysam.index("reads.bam")
with open("variants.vcf", "w") as handle:
    handle.write(f'##fileformat=VCFv4.2\n##contig=<ID=chr1,length={length}>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')
    for position in sites:
        handle.write(f"chr1\t{position + 1}\t.\t{reference[position]}\t{flip[reference[position]]}\t50\tPASS\t.\tGT\t0/1\n")
subprocess.run(["whatshap", "phase", "-o", "phased.vcf", "--reference", "ref.fa", "variants.vcf", "reads.bam"], check=True, capture_output=True)
print("ready")
```

- **示例输出：**

```text
ready
```

## 1) `dnakit.molecular_lab.primer3_design_by_task()`

- **作用：** 运行 Primer3 的任意任务：挑选引物对、检查给定引物、挑选带排名的引物列表、测序引物、克隆引物或等位基因区分引物。
- **处理方法：** 调用 primer3-py 的 `design_primers`，`PRIMER_TASK` 取 `generic`、`check_primers`、`pick_primer_list`、`pick_sequencing_primers`、`pick_cloning_primers` 或 `pick_discriminative_primers`（对应 `primer3_core PRIMER_TASK=...`）。其他 Primer3 标签放在 `settings` 中，未设置的保持 Primer3 默认值。返回的寡核苷酸含 0 起始位置、长度、Tm、GC、罚分及热力学字段。
- **API：** `dnakit.molecular_lab.primer3_design_by_task(template, *, task='generic', left_primer=None, right_primer=None, internal_oligo=None, target=None, included_region=None, excluded_regions=None, product_size_range=None, num_return=None, settings=None)`.
- **输入：** DNA 模板、`task`，以及可选的 `left_primer`、`right_primer`、`internal_oligo`（check_primers）、`target`、`included_region`、`excluded_regions`（0 起始 `[start, length]`）、`product_size_range`、`num_return`、`settings`（`PRIMER_*` 标签）。
- **输出：** `status`、`task`、`left`、`right`、`internal`（寡核苷酸行）、`pairs`、`explain`（Primer3 的计数）以及 `error` / `warning`。
- **示例代码：**

```python
from dnakit.molecular_lab import primer3_p14_functions as p3

template = open("template.txt").read()
ranked = p3.primer3_design_by_task(template, task="pick_primer_list", num_return=3, verbose=False)
first = ranked["left"][0]
print(ranked["status"], len(ranked["left"]), len(ranked["right"]), len(ranked["pairs"]))
print(template[first["start"]:first["start"] + first["length"]] == first["sequence"], first["length"])
checked = p3.primer3_design_by_task(template, task="check_primers", left_primer=first["sequence"], verbose=False)
print(checked["left"][0]["start"] == first["start"], round(checked["left"][0]["tm"], 2) == round(first["tm"], 2), checked["explain"]["PRIMER_LEFT_EXPLAIN"])
sequencing = p3.primer3_design_by_task(template, task="pick_sequencing_primers", target=[300, 20], verbose=False)
print(sequencing["status"], len(sequencing["left"]) > 0)
```

- **示例输出：**

```text
completed 3 3 0
True 20
True True considered 1, ok 1
completed True
```

## 2) `dnakit.molecular_lab.ntdpal_align_oligos()`

- **作用：** 用 Primer3 的动态规划比对程序对两条短 DNA 序列做全局、全局末端锚定、局部或局部末端锚定比对。
- **处理方法：** 运行 `ntdpal`，模式字母为 `g`（全局）、`G`（全局末端锚定，也可写作 `e`）、`l`（局部）、`L`（局部末端锚定），并解析得分、长度和比对位置对。
- **API：** `dnakit.molecular_lab.ntdpal_align_oligos(sequence_1, sequence_2, *, mode='local', gap_open_penalty=None, gap_extend_penalty=None, max_gap=None, gc_weighted_scores=False, ambiguity_codes=False, primer3_path=None)`.
- **输入：** 两条 DNA 字符串、`mode`，以及可选的 `gap_open_penalty`、`gap_extend_penalty`、`max_gap`、`gc_weighted_scores`（G/C 匹配 3、A/T 2、错配 -0.5）、`ambiguity_codes`。
- **输出：** 通用结果键，加上 `score`、`length`、`pairs`（0 起始下标对），局部模式还有 `end_1` / `end_2`。
- **示例代码：**

```python
from dnakit.molecular_lab import primer3_p14_functions as p3

for mode in ("local", "global", "global_end_anchored", "local_end_anchored"):
    result = p3.ntdpal_align_oligos("GCTAGCTAGGATCC", "GGATCCAAAA", mode=mode, verbose=False)
    print(mode, result["score"], result["length"], result["pairs"][0], result["pairs"][-1])
weighted = p3.ntdpal_align_oligos("AAAGGCCTTT", "GGCC", mode="local", gc_weighted_scores=True, verbose=False)
print(weighted["score"], weighted["end_1"], weighted["end_2"])
```

- **示例输出：**

```text
local 6.0 6 [8, 0] [13, 5]
global 6.0 6 [8, 0] [13, 5]
global_end_anchored 6.0 6 [8, 0] [13, 5]
local_end_anchored 6.0 6 [8, 0] [13, 5]
12.0 6 3
```

## 3) `dnakit.molecular_lab.primer3_mask_repeats()`

- **作用：** 屏蔽 FASTA 序列中的重复 k-mer 使引物避开它们，可按 k-mer 概率或绝对计数。
- **处理方法：** 运行 `primer3_masker`。默认当某位置 k-mer 来自重复序列的概率达到 `probability_cutoff`（0.1）时屏蔽，使用 `kmer_lists_path` 中的 `<list_prefix>_11.list` 和 `<list_prefix>_16.list`；给出 `absolute_cutoff` 和 `kmer_list` 时使用计数模式（`-a`），屏蔽在该列表中出现次数不少于阈值的 k-mer。
- **API：** `dnakit.molecular_lab.primer3_mask_repeats(input_fasta, output, *, kmer_lists_path, list_prefix='homo_sapiens', probability_cutoff=None, absolute_cutoff=None, kmer_list=None, mask_5_prime=None, mask_3_prime=None, masking_character=None, soft_mask=False, masking_direction=None, primer3_path=None)`.
- **输入：** FASTA、新的输出路径、`kmer_lists_path`（Primer3 源码中的 `kmer_lists` 目录：`homo_sapiens`、`mus_musculus`），以及可选的 `list_prefix`、`probability_cutoff`、`absolute_cutoff` + `kmer_list`、`mask_5_prime`、`mask_3_prime`、`masking_character`、`soft_mask`、`masking_direction`。
- **输出：** 通用结果键，加上 `mode`（`probability` 或 `absolute_count`）；屏蔽后的 FASTA 即输出文件。
- **示例调用（此处不运行）：** k-mer 列表文件（`homo_sapiens_11.list`、`homo_sapiens_16.list`，约 95 MB）随 Primer3 源码提供，conda 不会安装，因此该调用不在文档测试中运行。用列表实测时 poly-A 尾部被屏蔽：输出形如 `ACGTACGTAGCTAGCTAGCTAGCNNNNNNNNNNAAAANNNNNNNNNNN`，被屏蔽碱基为 N（`soft_mask=True` 时为小写）。

```python
from dnakit.molecular_lab import primer3_p14_functions as p3

open("to_mask.fa", "w").write(">s\nACGTACGTAGCTAGCTAGCTAGCTAAAAAAAAAAAAAAAAAAAAAAAA\n")
result = p3.primer3_mask_repeats("to_mask.fa", "masked.fa", kmer_lists_path="kmer_lists")
print(result["status"], result["mode"], open("masked.fa").read().splitlines()[1])
```

## 4) `dnakit.assembly.spades_correct_reads()`

- **作用：** 不做组装，只纠正短读段的测序错误：Illumina 类 reads 用 BayesHammer，Ion Torrent reads 用 IonHammer。
- **处理方法：** 运行 `spades.py --only-error-correction`（内部调用 `spades-hammer`）；`iontorrent=True` 时加 `--iontorrent`（`spades-ionhammer`）。纠正后的 FASTQ 写在 `<output_dir>/corrected`。
- **API：** `dnakit.assembly.spades_correct_reads(output_dir, *, reads_1=None, reads_2=None, unpaired_reads=None, iontorrent=False, threads=None, memory_gb=None, phred_offset=None, spades_path=None)`.
- **输入：** 新的 `output_dir`、`reads_1` + `reads_2`（一对文件）和/或 `unpaired_reads`，以及可选的 `iontorrent`、`threads`、`memory_gb`、`phred_offset`（33 或 64；质量值无法自动判断时需给出）。
- **输出：** 通用结果键，加上 `corrector` 和 `corrected_files`。
- **示例代码：**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_correct_reads("corrected_reads", reads_1=["r1.fq"], reads_2=["r2.fq"], threads=2, phred_offset=33, verbose=False)
print(result["status"], result["corrector"], sorted(path.rsplit("/", 1)[-1] for path in result["corrected_files"] if path.endswith(".gz")))
```

- **示例输出：**

```text
completed bayeshammer ['r100.0_0.cor.fastq.gz', 'r200.0_0.cor.fastq.gz', 'r_unpaired00.0_0.cor.fastq.gz']
```

## 5) `dnakit.assembly.spades_filter_reads_by_kmer_coverage()`

- **作用：** 去除 k-mer 计数中位数过低的 reads。
- **处理方法：** 用自动写出的数据集描述运行 `spades-read-filter -k <kmer> -c <coverage_threshold>`。当一对 reads 两条的 k-mer 计数中位数都不高于阈值时，该对被丢弃。
- **API：** `dnakit.assembly.spades_filter_reads_by_kmer_coverage(output_dir, *, kmer, coverage_threshold, reads_1=None, reads_2=None, unpaired_reads=None, threads=None, drop_names=False, drop_quality=False, spades_path=None)`.
- **输入：** 新的 `output_dir`、`kmer`、`coverage_threshold`、reads（`reads_1` + `reads_2`，或 `unpaired_reads`），以及可选的 `threads`、`drop_names`、`drop_quality`。
- **输出：** 通用结果键；`outputs` 列出过滤后的 FASTQ 文件和 `dataset.yaml`。
- **示例代码：**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_filter_reads_by_kmer_coverage("filtered_reads", kmer=21, coverage_threshold=3, reads_1=["r1.fq"], reads_2=["r2.fq"], verbose=False)
print(result["status"], sorted(path.rsplit("/", 1)[-1] for path in result["outputs"]))
```

- **示例输出：**

```text
completed ['1.1.fastq', '1.2.fastq', 'dataset.yaml']
```

## 6) `dnakit.assembly.spades_estimate_kmer_count()`

- **作用：** 估计一组 reads 中不同 k-mer 的数量。
- **处理方法：** 用自动写出的数据集描述运行 `spades-kmer-estimating -k <kmer>`；反向互补 reads 的 k-mer 不单独计数。
- **API：** `dnakit.assembly.spades_estimate_kmer_count(kmer, *, reads_1=None, reads_2=None, unpaired_reads=None, threads=None, spades_path=None)`.
- **输入：** `kmer`、reads（`reads_1` + `reads_2`，或 `unpaired_reads`），以及可选的 `threads`。
- **输出：** 通用结果键，加上 `distinct_kmers`。
- **示例代码：**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_estimate_kmer_count(21, reads_1=["r1.fq"], reads_2=["r2.fq"], verbose=False)
print(result["status"], result["distinct_kmers"])
```

- **示例输出：**

```text
completed 38355
```

## 7) `dnakit.assembly.spades_build_assembly_graph()`

- **作用：** 从 reads 或 FASTA 文件构建 de Bruijn 组装图。
- **处理方法：** 运行 `spades-gbuilder <数据集或 FASTA> <输出> -k <kmer>`，按 `graph_format` 加 `--unitigs`、`--gfa`、`--fastg` 或 `--spades`；`infer_coverage` 加 `-c`。
- **API：** `dnakit.assembly.spades_build_assembly_graph(output, *, kmer, reads_1=None, reads_2=None, unpaired_reads=None, fasta=None, graph_format='unitigs', infer_coverage=False, threads=None, spades_path=None)`.
- **输入：** 新的 `output` 路径、`kmer`、reads（`reads_1` + `reads_2`，或 `unpaired_reads`）或 `fasta` 二选一，以及可选的 `graph_format`（`unitigs`、`gfa`、`fastg`、`spades`）、`infer_coverage`、`threads`。
- **输出：** 通用结果键；图即输出文件。
- **示例代码：**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

genome_graph = sp.spades_build_assembly_graph("genome.gfa", kmer=21, fasta="genome.fa", graph_format="gfa", verbose=False)
read_graph = sp.spades_build_assembly_graph("reads.gfa", kmer=21, reads_1=["r1.fq"], reads_2=["r2.fq"], graph_format="gfa", infer_coverage=True, verbose=False)
segments = [line for line in open("genome.gfa") if line.startswith("S\t")]
print(genome_graph["status"], read_graph["status"], len(segments) >= 3)
```

- **示例输出：**

```text
completed completed True
```

## 8) `dnakit.assembly.spades_simplify_assembly_graph()`

- **作用：** 通过去除分支末梢、气泡和低覆盖度边来简化组装图。
- **处理方法：** 使用临时目录运行 `spades-gsimplifier <graph.gfa> <前缀> -k <kmer> --read-length <读长>`，结果为 `<前缀>.gfa`。
- **API：** `dnakit.assembly.spades_simplify_assembly_graph(graph, output_prefix, *, kmer, read_length, coverage=None, use_coverage_ratios=False, threads=None, spades_path=None)`.
- **输入：** GFA 图、新的 `output_prefix`、`kmer`、`read_length`，以及可选的 `coverage`（数值或 `auto`）、`use_coverage_ratios`、`threads`。
- **输出：** 通用结果键；`outputs` 为 `<前缀>.gfa`。
- **示例代码：**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_simplify_assembly_graph("genome.gfa", "simplified", kmer=21, read_length=100, verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **示例输出：**

```text
completed ['simplified.gfa']
```

## 9) `dnakit.assembly.spades_split_assembly_graph()`

- **作用：** 把组装图按连通分量拆成多个 GFA 文件。
- **处理方法：** 运行 `spades-gfa-split <graph.gfa> <输出目录>`；程序创建该目录并写出 `subgraph_<n>.gfa`。
- **API：** `dnakit.assembly.spades_split_assembly_graph(graph, output_dir, *, spades_path=None)`.
- **输入：** GFA 图和新的 `output_dir`。
- **输出：** 通用结果键，加上 `components`；`outputs` 列出各分量文件。
- **示例代码：**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_split_assembly_graph("genome.gfa", "components", verbose=False)
print(result["status"], result["components"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **示例输出：**

```text
completed 1 ['subgraph_0.gfa']
```

## 10) `dnakit.assembly.spades_map_reads_to_graph()`

- **作用：** 把 reads 比对到组装图，并写出双端距离直方图。
- **处理方法：** 运行 `spades-gmapper <数据集> <graph.gfa> <输出> -k <kmer>`，按需加库编号、内部比对模式（`-X`）、保留比对策略（`-R`）和 `--hic`。
- **API：** `dnakit.assembly.spades_map_reads_to_graph(graph, output, *, kmer, reads_1=None, reads_2=None, unpaired_reads=None, library_index=None, mode=None, retain=None, hic=False, threads=None, spades_path=None)`.
- **输入：** GFA 图、新的 `output` 路径、`kmer`、reads（`reads_1` + `reads_2`，或 `unpaired_reads`），以及可选的 `library_index`、`mode`（`illumina`、`pacbio`、`intractg`、`hic`）、`retain`（`default`、`all`、`primary`、`quality`）、`hic`、`threads`。
- **输出：** 通用结果键；距离直方图即输出文件。
- **示例代码：**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spades_map_reads_to_graph("genome.gfa", "distances.txt", kmer=21, reads_1=["r1.fq"], reads_2=["r2.fq"], verbose=False)
print(result["status"], "weight" in open("distances.txt").read())
```

- **示例输出：**

```text
completed True
```

## 11) `dnakit.assembly.spaligner_align_long_reads()`

- **作用：** 用 SPAligner 把长读段或 contig 比对到组装图。
- **处理方法：** 运行 `spaligner <config.yaml> -d <数据类型> -s <序列> -g <图> -k <kmer> -o <输出目录>`。未给 `config` 时使用 SPAligner 安装的 `spaligner_config.yaml`。
- **API：** `dnakit.assembly.spaligner_align_long_reads(sequences, graph, output_dir, *, datatype, kmer, config=None, threads=None, spades_path=None)`.
- **输入：** 序列的 FASTA/FASTQ、GFA 图（或 SPAdes saves 目录）、新的 `output_dir`、`datatype`（`nanopore`、`pacbio`）、`kmer`（奇数），以及可选的 `config`、`threads`。
- **输出：** 通用结果键；`outputs` 为 `alignment.tsv`（每条比对上的 read 一行，含其图路径）。
- **示例代码：**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.spaligner_align_long_reads("long.fa", "genome.gfa", "long_alignments", datatype="nanopore", kmer=21, verbose=False)
rows = open("long_alignments/alignment.tsv").read().splitlines()
print(result["status"], len(rows), [row.split("\t")[0] for row in rows])
```

- **示例输出：**

```text
completed 4 ['l0', 'l1', 'l2', 'l3']
```

## 12) `dnakit.assembly.pathracer_search_graph()`

- **作用：** 在组装图中寻找与 profile HMM 或序列匹配的路径（PathRacer）。
- **处理方法：** 运行 `pathracer <查询> <图> --output <目录>`，查询类型用 `--hmm`、`--nt` 或 `--aa`，匹配模式用 `--global` 或 `--local`；`minimum_length`、`indel_rate`、`top`、`seeding`、`evalue`、`threads`、`memory_gb` 对应 PathRacer 的选项。
- **API：** `dnakit.assembly.pathracer_search_graph(query, graph, output_dir, *, query_type='hmm', mode='global', minimum_length=None, indel_rate=None, top=None, threads=None, memory_gb=None, seeding=None, evalue=None, spades_path=None)`.
- **输入：** 查询文件（HMMER3 profile，或核酸/氨基酸字符串的 FASTA）、GFA 图或 saves 目录、新的 `output_dir`、`query_type`、`mode`，以及上述可选参数。
- **输出：** 通用结果键；`outputs` 列出匹配序列（`*.seqs.fa`）、边、带注释的图和日志。
- **示例代码：**

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.pathracer_search_graph("query.fa", "genome.gfa", "pathracer_out", query_type="nt", threads=2, verbose=False)
print(result["status"], sorted(path.rsplit("/", 1)[-1] for path in result["outputs"]))
```

- **示例输出：**

```text
completed ['graph_with_hmm_paths.gfa', 'pathracer.log', 'q.edges.fa', 'q.seqs.fa']
```

## 13) `dnakit.assembly.pathracer_align_hmm_to_sequences()`

- **作用：** 允许移码，把氨基酸 profile HMM 比对到核酸序列。
- **处理方法：** 运行 `pathracer-seq-fs <hmm> <序列> --output <目录>`，按需加 `--global` 或 `--local`、`--max-fs`、`--indel-rate`、`--top`、`--cutoff`、`--exhaustive`。
- **API：** `dnakit.assembly.pathracer_align_hmm_to_sequences(hmm, sequences, output_dir, *, mode='global', indel_rate=None, top=None, max_frameshifts=None, cutoff=None, exhaustive=False, threads=None, memory_gb=None, evalue=None, spades_path=None)`.
- **输入：** HMMER3 氨基酸 profile、核酸 FASTA、新的 `output_dir`，以及可选的 `mode`、`indel_rate`、`top`、`max_frameshifts`、`cutoff`、`exhaustive`、`threads`、`memory_gb`、`evalue`。
- **输出：** 通用结果键；`outputs` 列出结果文件。
- **示例调用（此处不运行）：** 该调用需要真实蛋白家族的 HMMER3 氨基酸 profile，而这里没有 HMMER 无法生成，因此不在文档测试中运行；单元测试用替身程序检查参数列表。

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.pathracer_align_hmm_to_sequences("protein.hmm", "contigs.fa", "frameshift_out", max_frameshifts=4)
print(result["status"], result["outputs"])
```

## 14) `dnakit.assembly.splitter_resolve_repeats()`

- **作用：** 用合成长读段或连锁读段解析组装图中的重复序列（SPlitteR）。
- **处理方法：** 运行 `splitter <图> <library.yaml> <输出目录>`，按需加图类型（`-G`）、解析模式（`-M`）、得分阈值、`--scaffold-links`、`--mapping-k` 和 barcode 下采样选项。
- **API：** `dnakit.assembly.splitter_resolve_repeats(graph, library_yaml, output_dir, *, library_index=None, graph_type=None, resolution_mode=None, score_threshold=None, relative_threshold=None, scaffold_links=False, mapping_kmer=None, sampling_factor=None, seed=None, assembly_info=None, threads=None, spades_path=None)`.
- **输入：** GFA（或二进制）图、带 barcode 读段文库的 SPAdes YAML 描述、新的 `output_dir`，以及可选的 `library_index`、`graph_type`、`resolution_mode`、`score_threshold`、`relative_threshold`、`scaffold_links`、`mapping_kmer`、`sampling_factor`、`seed`、`assembly_info`、`threads`。
- **输出：** 通用结果键；`outputs` 列出解析后图的文件。
- **示例调用（此处不运行）：** 该调用需要合成长读段（带 barcode）数据，因此不在文档测试中运行；单元测试用替身程序检查参数列表。

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.splitter_resolve_repeats("assembly_graph.gfa", "slr_library.yaml", "resolved", resolution_mode="meta", graph_type="mdbg")
print(result["status"])
```

## 15) `dnakit.assembly.binspreader_refine_bins()`

- **作用：** 在组装图上细化并传播宏基因组分箱结果（BinSPreader）。
- **处理方法：** 运行 `binspreader <图> <binning.tsv> <输出.tsv>`，按需加 contig 路径、分配策略（`-Smax`/`-Smle`）、细化器（`-Rcorr`/`-Rprop`）、多重分配（`-m`）、CAMI 格式、稀疏传播和 reads 拆分选项。
- **API：** `dnakit.assembly.binspreader_refine_bins(graph, binning, output, *, contig_paths=None, dataset_yaml=None, multiple_assignment=False, assignment_strategy=None, refiner=None, cami_format=False, zero_bin=False, labeled_alpha=None, sparse_propagation=False, split_reads=False, read_bin_weight_threshold=None, threads=None, spades_path=None)`.
- **输入：** GFA（或二进制）图、分箱程序给出的 contig 到 bin 的 TSV、新的输出路径，以及可选的 `contig_paths`、`dataset_yaml`、`multiple_assignment`、`assignment_strategy`、`refiner`、`cami_format`、`zero_bin`、`labeled_alpha`、`sparse_propagation`、`split_reads`（需要 `dataset_yaml`）、`read_bin_weight_threshold`、`threads`。
- **输出：** 通用结果键；细化后的分箱结果即输出文件。
- **示例调用（此处不运行）：** 该调用需要宏基因组组装图和分箱程序的输出，因此不在文档测试中运行；单元测试用替身程序检查参数列表。

```python
from dnakit.assembly import spades_p14_cli_functions as sp

result = sp.binspreader_refine_bins("assembly_graph.gfa", "bins.tsv", "refined.tsv", contig_paths="contigs.paths", assignment_strategy="mle", refiner="prop")
print(result["status"])
```

## 16) `dnakit.variants.whatshap_phasing_statistics()`

- **作用：** 计算已定相 VCF 的定相统计。
- **处理方法：** 运行 `whatshap stats --tsv <文件>`（按需加 `--gtf`、`--block-list`、`--chr-lengths`）并解析 TSV。
- **API：** `dnakit.variants.whatshap_phasing_statistics(vcf, tsv, *, gtf=None, block_list=None, chromosome_lengths=None, sample=None, chromosome=None, only_snvs=False, whatshap_path=None)`.
- **输入：** 已定相 VCF、新的 `tsv` 路径，以及可选的 `gtf`、`block_list`、`chromosome_lengths`、`sample`、`chromosome`、`only_snvs`。
- **输出：** 通用结果键，加上 `statistics`（每条染色体一行外加 `ALL`：变异数、已定相数、块数、块 N50 等）。
- **示例代码：**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_phasing_statistics("phased.vcf", "stats.tsv", gtf="blocks.gtf", block_list="blocks.txt", verbose=False)
row = result["statistics"][0]
print(result["status"], row["chromosome"], row["variants"], row["phased"], row["blocks"], row["block_n50"])
```

- **示例输出：**

```text
completed chr1 25 25 1 3600
```

## 17) `dnakit.variants.whatshap_unphase_vcf()`

- **作用：** 去除 VCF 中的定相信息。
- **处理方法：** 运行 `whatshap unphase`：`1|0` 之类已定相基因型变为 `0/1`，PS 和 HP 字段被去除。
- **API：** `dnakit.variants.whatshap_unphase_vcf(vcf, output, *, whatshap_path=None)`.
- **输入：** VCF 和新的输出路径。
- **输出：** 通用结果键；去相后的 VCF 即输出文件。
- **示例代码：**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_unphase_vcf("phased.vcf", "unphased.vcf", verbose=False)
genotypes = {line.split("\t")[9].split(":")[0].strip() for line in open("unphased.vcf") if not line.startswith("#")}
print(result["status"], genotypes)
```

- **示例输出：**

```text
completed {'0/1'}
```

## 18) `dnakit.variants.whatshap_compare_phasing()`

- **作用：** 比较两个或多个已定相 VCF：switch error、flip 和汉明距离。
- **处理方法：** 运行 `whatshap compare`，带 `--tsv-pairwise`、`--tsv-multiway`、`--switch-error-bed` 和/或 `--longest-block-tsv`；至少要一个输出，且所有输出须在同一目录。
- **API：** `dnakit.variants.whatshap_compare_phasing(vcfs, *, tsv_pairwise=None, tsv_multiway=None, switch_error_bed=None, longest_block_tsv=None, names=None, ignore_sample_name=False, only_snvs=False, sample=None, ploidy=None, whatshap_path=None)`.
- **输入：** 至少两个已定相 VCF 的列表、输出路径，以及可选的 `names`、`ignore_sample_name`、`only_snvs`、`sample`、`ploidy`。
- **输出：** 通用结果键；给出 `tsv_pairwise` 时还有 `pairwise`（成对比较表的行）。
- **示例代码：**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_compare_phasing(["phased.vcf", "phased.vcf"], tsv_pairwise="compare.tsv", names=["a", "b"], verbose=False)
row = result["pairwise"][0]
print(result["status"], row["all_assessed_pairs"], row["all_switches"], row["blockwise_hamming"])
```

- **示例输出：**

```text
completed 24 0 0
```

## 19) `dnakit.variants.whatshap_haplotag_reads()`

- **作用：** 用已定相 VCF 给每条 read 打上单倍型（HP）和相位集（PS）标签。
- **处理方法：** 运行 `whatshap haplotag -o <bam>`，可选 `--output-haplotag-list`。VCF 必须经 bgzip 压缩并建有 tabix 索引。
- **API：** `dnakit.variants.whatshap_haplotag_reads(vcf, alignments, output, *, reference=None, haplotag_list=None, regions=None, ploidy=None, sample=None, ignore_read_groups=False, ignore_linked_reads=False, linked_read_distance_cutoff=None, skip_missing_contigs=False, tag_supplementary=None, threads=None, whatshap_path=None)`.
- **输入：** 经 bgzip 压缩并建索引的已定相 VCF、BAM/CRAM、新的输出路径，以及可选的 `reference`、`haplotag_list`、`regions`、`ploidy`、`sample`、`ignore_read_groups`、`ignore_linked_reads`、`linked_read_distance_cutoff`、`skip_missing_contigs`、`tag_supplementary`、`threads`。
- **输出：** 通用结果键；写出带标签的 BAM 和可选的 TSV（read、单倍型、相位集、染色体）。
- **示例代码：**

```python
import pysam

from dnakit.variants import whatshap_p14_cli_functions as wh

pysam.tabix_index("phased.vcf", preset="vcf", force=True, keep_original=True)
result = wh.whatshap_haplotag_reads("phased.vcf.gz", "reads.bam", "tagged.bam", reference="ref.fa", haplotag_list="haplotags.tsv", verbose=False)
tags = {}
for line in open("haplotags.tsv").read().splitlines()[1:]:
    name, haplotype = line.split("\t")[:2]
    tags[name] = haplotype
even = {tags[name] for name in tags if int(name[1:]) % 2 == 0}
odd = {tags[name] for name in tags if int(name[1:]) % 2 == 1}
print(result["status"], len(tags), len(even), len(odd), even != odd)
```

- **示例输出：**

```text
completed 500 1 1 True
```

## 20) `dnakit.variants.whatshap_split_reads()`

- **作用：** 把每个单倍型的 reads 写入各自的文件。
- **处理方法：** 运行 `whatshap split <reads> <列表>`，按需带 `--output-h1`、`--output-h2`、`--output-untagged` 和 `--read-lengths-histogram`。
- **API：** `dnakit.variants.whatshap_split_reads(reads, haplotype_list, *, output_h1=None, output_h2=None, output_untagged=None, add_untagged=False, only_largest_block=False, discard_unknown_reads=False, read_lengths_histogram=None, whatshap_path=None)`.
- **输入：** reads 的 FASTQ/BAM、单倍型列表（read 名与 H1/H2，例如 `whatshap_haplotag_reads` 的输出）、至少一个输出路径，以及可选的 `add_untagged`、`only_largest_block`、`discard_unknown_reads`。
- **输出：** 通用结果键；写出各单倍型文件。
- **示例代码：**

```python
import pysam

from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_split_reads("reads.bam", "haplotags.tsv", output_h1="h1.bam", output_h2="h2.bam", verbose=False)
counts = [sum(1 for _ in pysam.AlignmentFile(name)) for name in ("h1.bam", "h2.bam")]
print(result["status"], counts, sum(counts))
```

- **示例输出：**

```text
completed [250, 250] 500
```

## 21) `dnakit.variants.whatshap_phase_from_haplotags()`

- **作用：** 根据带单倍型标签的 reads 给 VCF 添加定相信息。
- **处理方法：** 运行 `whatshap haplotagphase -o <vcf>`；比对文件须带 HP/PS 标签并已建索引，且 WhatsHap 要求提供 `reference`。
- **API：** `dnakit.variants.whatshap_phase_from_haplotags(vcf, alignments, output, *, reference=None, gap_threshold=None, cut_poly=None, only_indels=False, no_multi_allelic=False, ignore_read_groups=False, sample=None, chromosome=None, whatshap_path=None)`.
- **输入：** VCF、已建索引的带标签 BAM/CRAM、新的输出路径、`reference`，以及可选的 `gap_threshold`、`cut_poly`、`only_indels`、`no_multi_allelic`、`ignore_read_groups`、`sample`、`chromosome`。
- **输出：** 通用结果键；定相后的 VCF 即输出文件。
- **示例代码：**

```python
import pysam

from dnakit.variants import whatshap_p14_cli_functions as wh

pysam.index("tagged.bam")
result = wh.whatshap_phase_from_haplotags("unphased.vcf", "tagged.bam", "rephased.vcf", reference="ref.fa", verbose=False)
phased = sum("|" in line.split("\t")[9] for line in open("rephased.vcf") if not line.startswith("#"))
print(result["status"], phased)
```

- **示例输出：**

```text
completed 25
```

## 22) `dnakit.variants.whatshap_genotype_variants()`

- **作用：** 根据 reads 或 reads 加家系对变异做基因分型。
- **处理方法：** 运行 `whatshap genotype`；给出 `pedigree` 时加 `--ped`（以及 `--genmap` 或 `--recombrate`、`--use-ped-samples`），对家系联合分型。
- **API：** `dnakit.variants.whatshap_genotype_variants(vcf, alignments, output, *, reference=None, pedigree=None, genetic_map=None, recombination_rate=None, use_ped_samples=False, max_coverage=None, mapping_quality=None, genotype_quality_threshold=None, no_priors=False, prior_output=None, sample=None, chromosome=None, only_snvs=False, whatshap_path=None)`.
- **输入：** VCF、BAM/CRAM 文件列表、新的输出路径，以及可选的 `reference`、`pedigree`、`genetic_map`、`recombination_rate`、`use_ped_samples`、`max_coverage`、`mapping_quality`、`genotype_quality_threshold`、`no_priors`、`prior_output`、`sample`、`chromosome`、`only_snvs`。
- **输出：** 通用结果键；分型后的 VCF 即输出文件。
- **示例代码：**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_genotype_variants("variants.vcf", ["reads.bam"], "genotyped.vcf", reference="ref.fa", verbose=False)
calls = {line.split("\t")[9].split(":")[0] for line in open("genotyped.vcf") if not line.startswith("#")}
print(result["status"], sorted(calls))
```

- **示例输出：**

```text
completed ['1/0']
```

## 23) `dnakit.variants.whatshap_find_snv_candidates()`

- **作用：** 根据 reads 比对寻找候选 SNV 位点。
- **处理方法：** 运行 `whatshap find_snv_candidates`；当替代等位基因至少有 `min_alt_depth` 条 reads（默认 3）且占深度的 `min_alt_fraction` 时，该位点为候选。`technology` 套用 Illumina、PacBio 或 Nanopore 的默认值。
- **API：** `dnakit.variants.whatshap_find_snv_candidates(reference, alignment, output, *, technology=None, min_alt_depth=None, min_alt_fraction=None, multi_allelic=False, sample=None, chromosome=None, whatshap_path=None)`.
- **输入：** 参考 FASTA、已建索引的 BAM、新的输出路径，以及可选的 `technology`、`min_alt_depth`、`min_alt_fraction`、`multi_allelic`、`sample`、`chromosome`。
- **输出：** 通用结果键；候选 VCF 即输出文件。
- **示例代码：**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_find_snv_candidates("ref.fa", "reads.bam", "candidates.vcf", technology="illumina", verbose=False)
found = {int(line.split("\t")[1]) for line in open("candidates.vcf") if not line.startswith("#")}
truth = {int(line.split("\t")[1]) for line in open("variants.vcf") if not line.startswith("#")}
print(result["status"], len(found), truth <= found)
```

- **示例输出：**

```text
completed 25 True
```

## 24) `dnakit.variants.whatshap_learn_error_profile()`

- **作用：** 根据比对和变异学习特定测序技术的错误谱。
- **处理方法：** 运行 `whatshap learn --reference <fasta> --output <文件> <bam> <vcf>`；输出统计变异周围（read 与参考序列）的 k-mer 配对计数。
- **API：** `dnakit.variants.whatshap_learn_error_profile(alignment, vcf, reference, output, *, kmer_size=None, window=None, whatshap_path=None)`.
- **输入：** 已建索引的 BAM、变异 VCF、参考 FASTA、新的输出路径，以及可选的 `kmer_size`（默认 4）、`window`。
- **输出：** 通用结果键；错误谱即输出文件。
- **示例代码：**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_learn_error_profile("reads.bam", "phased.vcf", "ref.fa", "error_profile.txt", verbose=False)
print(result["status"], len(open("error_profile.txt").read().splitlines()) > 0)
```

- **示例输出：**

```text
completed True
```

## 25) `dnakit.variants.whatshap_polyphase()`

- **作用：** 用聚类加线程算法根据 reads 对多倍体变异定相。
- **处理方法：** 运行 `whatshap polyphase --ploidy <n>`，按需加块切分灵敏度、最小重叠、预定相和单倍集选项。倍性为 2 时对二倍体样本定相。
- **API：** `dnakit.variants.whatshap_polyphase(vcf, alignments, output, *, ploidy, reference=None, block_cut_sensitivity=None, min_overlap=None, use_prephasing=False, include_haploid_sets=False, no_multi_allelic=False, mapping_quality=None, threads=None, sample=None, chromosome=None, tag=None, only_snvs=False, distrust_genotypes=False, whatshap_path=None)`.
- **输入：** VCF、BAM/CRAM 文件列表、新的输出路径、`ploidy`，以及可选的 `reference`、`block_cut_sensitivity`（0-5）、`min_overlap`、`use_prephasing`、`include_haploid_sets`、`no_multi_allelic`、`mapping_quality`、`threads`、`sample`、`chromosome`、`tag`、`only_snvs`、`distrust_genotypes`。
- **输出：** 通用结果键；定相后的 VCF 即输出文件。
- **示例代码：**

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_polyphase("variants.vcf", ["reads.bam"], "polyphased.vcf", ploidy=2, reference="ref.fa", verbose=False)
phased = sum("|" in line.split("\t")[9] for line in open("polyphased.vcf") if not line.startswith("#"))
print(result["status"], phased)
```

- **示例输出：**

```text
completed 25
```

## 26) `dnakit.variants.whatshap_phase_with_pedigree()`

- **作用：** 结合 reads 与家系对变异定相。
- **处理方法：** 运行 `whatshap phase --ped <ped> -o <vcf>`，按需加遗传图谱或重组率、重组列表、`--no-genetic-haplotyping` 和 `--use-ped-samples`。
- **API：** `dnakit.variants.whatshap_phase_with_pedigree(vcf, alignments, pedigree, output, *, reference=None, genetic_map=None, recombination_rate=None, recombination_list=None, no_genetic_haplotyping=False, use_ped_samples=False, algorithm=None, mapping_quality=None, sample=None, chromosome=None, tag=None, only_snvs=False, distrust_genotypes=False, whatshap_path=None)`.
- **输入：** 包含全部家系成员的 VCF 和 BAM/CRAM 文件、PED/FAM 文件、新的输出路径，以及可选的 `reference`、`genetic_map`、`recombination_rate`、`recombination_list`、`no_genetic_haplotyping`、`use_ped_samples`、`algorithm`、`mapping_quality`、`sample`、`chromosome`、`tag`、`only_snvs`、`distrust_genotypes`。
- **输出：** 通用结果键；定相后的 VCF 即输出文件。
- **示例调用（此处不运行）：** 该调用需要家系数据（三人 VCF、每位成员的 BAM 和 PED 文件），因此不在文档测试中运行；单元测试用替身程序检查参数列表。

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_phase_with_pedigree("trio.vcf", ["child.bam", "father.bam", "mother.bam"], "trio.ped", "trio_phased.vcf", reference="ref.fa")
print(result["status"])
```

## 27) `dnakit.variants.whatshap_polyphase_with_progeny()`

- **作用：** 利用后代的基因型对多倍体亲本定相。
- **处理方法：** 运行 `whatshap polyphasegenetic --ploidy <n> <vcf> <家系文件>`，按需带后代文件、打分窗口和复杂度支持度。
- **API：** `dnakit.variants.whatshap_polyphase_with_progeny(vcf, pedigree, output, *, ploidy, progeny_file=None, scoring_window=None, complexity_support=None, sample=None, chromosome=None, tag=None, only_snvs=False, distrust_genotypes=False, whatshap_path=None)`.
- **输入：** VCF、家系文件、新的输出路径、`ploidy`，以及可选的 `progeny_file`、`scoring_window`、`complexity_support`、`sample`、`chromosome`、`tag`、`only_snvs`、`distrust_genotypes`。
- **输出：** 通用结果键；定相后的 VCF 即输出文件。
- **示例调用（此处不运行）：** 该调用需要多倍体家系 VCF 和家系文件，因此不在文档测试中运行；单元测试用替身程序检查参数列表。

```python
from dnakit.variants import whatshap_p14_cli_functions as wh

result = wh.whatshap_polyphase_with_progeny("tetraploid_family.vcf", "family.ped", "parents_phased.vcf", ploidy=4)
print(result["status"])
```
