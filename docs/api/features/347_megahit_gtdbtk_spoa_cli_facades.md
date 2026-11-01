# MEGAHIT、GTDB-Tk 与 SPOA 命令行适配器

这一组函数通过固定的参数白名单、不经过 shell 调用三个由使用者自行安装的程序：[MEGAHIT](https://github.com/voutcn/megahit) 的 `megahit_core`（读段统计、FASTG 转换、用 reads 对 contig 做局部组装）、[GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)（标记基因识别、比对、基于放置的分类、从头建树、与代表基因组的 ANI、MSA 修剪与导出）和 [SPOA](https://github.com/rvaser/spoa)（部分序有向比对，输出 MSA 或图）。缺少程序时报错并给出安装命令；输出文件和文件夹必须尚不存在。

测试版本：MEGAHIT 1.2.9、SPOA 4.1.6。测试环境中没有 GTDB-Tk 及其数十 GB 的参考包，因此它的七个函数用遵循 GTDB-Tk 命令行与输出文件名的替身程序测试，示例调用不执行。`megahit_core` 的其他子程序（`count`、`read2sdbg`、`seq2sdbg`、`iterate`、`assemble`）交换 MEGAHIT 私有的二进制文件，不单独提供；完整流程见 `dnakit.assembly.external.assemble_megahit`。

## 示例数据

准备代码写入 6 kb 随机基因组、一个 MEGAHIT 格式的 2 kb contig 文件、两条被断开的 contig 和 3000 对 reads；`s.fa` 是供 SPOA 使用的四条序列 FASTA。

**示例文件：** `s.fa`

```text
>a
ACGTACGTACGTAAC
>b
ACGTACCTACGTAAC
>c
ACGTACGTACGAAC
>d
GTTACGTACGTACGT
```

- **示例代码：**

```python
import random
from pathlib import Path

random.seed(5)
genome = "".join(random.choice("ACGT") for _ in range(6000))
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
Path("k21.contigs.fa").write_text(f">k21_0 flag=1 multi=10.0 len=2000\n{genome[:2000]}\n")
Path("split.fa").write_text(f">a\n{genome[:2500]}\n>b\n{genome[3000:5500]}\n")
with open("p1.fq", "w") as f1, open("p2.fq", "w") as f2:
    for i in range(3000):
        insert = random.randint(300, 400)
        start = random.randrange(0, len(genome) - insert)
        frag = genome[start:start + insert]
        f1.write(f"@r{i}/1\n{frag[:100]}\n+\n{'I' * 100}\n")
        f2.write(f"@r{i}/2\n{rc(frag)[:100]}\n+\n{'I' * 100}\n")
print(len(genome), Path("p1.fq").read_text().count("@r"))
```

- **示例输出：**

```text
6000 3000
```

## 1) `dnakit.assembly.megahit_p15_cli_functions.megahit_read_statistics()`

- **作用：** FASTA/FASTQ 文件的读段数、总碱基数、最长、最短和平均读长。
- **处理方法：** 对文件运行 `megahit_core readstat`（可为 gzip）。
- **API：** `dnakit.assembly.megahit_p15_cli_functions.megahit_read_statistics(reads, stats_path, *, megahit_core_path=None)`。
- **输入：** reads 文件和表格的新路径。
- **输出：** `read_count`、`total_bases`、`longest`、`shortest`、`mean_length`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.assembly.megahit_p15_cli_functions import megahit_read_statistics

result = megahit_read_statistics("p1.fq", str(out / "stats.txt"))
print(result["read_count"], result["total_bases"], result["longest"], result["shortest"], result["mean_length"])
```

- **示例输出：**

```text
3000 300000 100 100 100.0
```

## 2) `dnakit.assembly.megahit_p15_cli_functions.megahit_contigs_to_fastg()`

- **作用：** 把 MEGAHIT 的 contig 转换为 FASTG 组装图。
- **处理方法：** 运行 `megahit_core contig2fastg <k> <k*.contigs.fa>`；每个 contig 生成正向和反向两个节点。
- **API：** `dnakit.assembly.megahit_p15_cli_functions.megahit_contigs_to_fastg(contigs, k, fastg_path, *, megahit_core_path=None)`。
- **输入：** 带 `flag=`、`multi=` 标题的中间 contig 文件、`k` 和新的 FASTG 路径。
- **输出：** `fastg_file`、`node_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.assembly.megahit_p15_cli_functions import megahit_contigs_to_fastg

result = megahit_contigs_to_fastg("k21.contigs.fa", 21, str(out / "graph.fastg"))
print(result["node_count"], (out / "graph.fastg").read_text().splitlines()[0])
```

- **示例输出：**

```text
2 >NODE_1_length_2000_cov_10.0000_ID_1;
```

## 3) `dnakit.assembly.megahit_p15_cli_functions.megahit_local_assembly()`

- **作用：** 用比对到 contig 末端的 reads 做局部组装，延伸 contig 或填补缺口。
- **处理方法：** 先把 reads 转成 MEGAHIT 的二进制库（`megahit_core buildlib`），再运行 `megahit_core local`。双端 reads 使组装器能跨越缺口。
- **API：** `dnakit.assembly.megahit_p15_cli_functions.megahit_local_assembly(contigs, output_directory, *, single_reads=(), paired_reads=(), interleaved_reads=(), kmin=11, kmax=41, step=6, seed_kmer=31, min_contig_len=200, min_mapping_len=75, sparsity=8, similarity=0.8, threads=1, megahit_core_path=None)`。
- **输入：** contig 文件、新的输出文件夹，以及 `single_reads`、`paired_reads` 或 `interleaved_reads`。
- **输出：** `local_contig_file`、`local_contig_count`、`local_contig_bases`、`read_library`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.assembly.megahit_p15_cli_functions import megahit_local_assembly

result = megahit_local_assembly("split.fa", str(out / "local"), paired_reads=[["p1.fq", "p2.fq"]], min_contig_len=100)
print(result["input_contig_count"], result["local_contig_count"] >= 1, result["local_contig_bases"] > 0)
```

- **示例输出：**

```text
2 True True
```

## 4) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_identify_markers()`

- **作用：** 识别基因组中的 GTDB 古菌和细菌标记基因。
- **处理方法：** 运行 `gtdbtk identify`（Prodigal 与 HMMER）。参考包文件夹通过 `database` 传入，并仅对该进程设置为 `GTDBTK_DATA_PATH`。
- **API：** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_identify_markers(out_dir, *, genome_dir=None, batchfile=None, database=None, extension='fna', prefix='gtdbtk', genes=False, write_single_copy_genes=False, cpus=1, force=False, gtdbtk_path=None)`。
- **输入：** `genome_dir` 或 `batchfile`、新的 `out_dir`、`extension`。
- **输出：** `markers_summary` 表、`translation_table_summary`、`failed_genomes`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_identify_markers

result = gtdbtk_identify_markers("gtdb_identify", genome_dir="genomes", database="/data/gtdbtk/release226", cpus=8)
print(list(result["markers_summary"]))
```

## 5) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_align_markers()`

- **作用：** 把识别出的标记基因比对成带掩码的多序列比对。
- **处理方法：** 运行 `gtdbtk align`；`skip_trimming`、`custom_msa_filters`、`taxa_filter` 与一致性选项会被传入。
- **API：** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_align_markers(out_dir, identify_dir, *, database=None, prefix='gtdbtk', skip_gtdb_refs=False, taxa_filter=None, min_perc_aa=10.0, custom_msa_filters=False, skip_trimming=False, cols_per_gene=42, min_consensus=25.0, max_consensus=95.0, min_perc_taxa=50.0, rnd_seed=None, cpus=1, gtdbtk_path=None)`。
- **输入：** `identify` 的输出文件夹和新的 `out_dir`。
- **输出：** `files`（比对文件）、`filtered_genomes`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_align_markers

result = gtdbtk_align_markers("gtdb_align", "gtdb_identify", database="/data/gtdbtk/release226")
print(result["files"]["*.user_msa.fasta*"])
```

## 6) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_classify_placement()`

- **作用：** 对预先计算的标记比对做 pplacer 放置，从而对基因组分类。
- **处理方法：** 用 `align` 的输出运行 `gtdbtk classify`；分类由放置位置、相对进化分歧和 ANI 决定。
- **API：** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_classify_placement(out_dir, align_dir, *, genome_dir=None, batchfile=None, database=None, extension='fna', prefix='gtdbtk', full_tree=False, place_species=False, min_af=None, skip_pplacer=False, pplacer_cpus=None, cpus=1, gtdbtk_path=None)`。
- **输入：** `align` 文件夹、基因组和新的 `out_dir`。
- **输出：** `summary`（用户基因组、分类、最近参考、ANI、警告）、`trees`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_classify_placement

result = gtdbtk_classify_placement("gtdb_classify", "gtdb_align", genome_dir="genomes", database="/data/gtdbtk/release226")
print(result["summary"])
```

## 7) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_de_novo_workflow()`

- **作用：** 推断从头基因组树，并以 GTDB 分类定根和标注。
- **处理方法：** 对 `bacteria` 或 `archaea` 运行 `gtdbtk de_novo_wf`（identify、align、infer、root、decorate），需给出外群类群。
- **API：** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_de_novo_workflow(out_dir, outgroup_taxon, *, domain='bacteria', genome_dir=None, batchfile=None, database=None, extension='fna', prefix='gtdbtk', skip_gtdb_refs=False, taxa_filter=None, min_perc_aa=10.0, custom_msa_filters=False, cols_per_gene=42, min_consensus=25.0, max_consensus=95.0, min_perc_taxa=50.0, prot_model='WAG', no_support=False, gamma=False, rnd_seed=None, genes=False, write_single_copy_genes=False, keep_intermediates=False, force=False, cpus=1, gtdbtk_path=None)`。
- **输入：** 基因组、`outgroup_taxon`、`domain` 和新的 `out_dir`。
- **输出：** `trees`（标注和定根的 Newick 文件）、`files`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_de_novo_workflow

result = gtdbtk_de_novo_workflow("gtdb_tree", "p__Patescibacteria", genome_dir="genomes", database="/data/gtdbtk/release226", prot_model="LG")
print(result["trees"]["*.decorated.tree"])
```

## 8) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_ani_to_representatives()`

- **作用：** 计算基因组与 GTDB 代表基因组的平均核苷酸一致性。
- **处理方法：** 运行 `gtdbtk ani_rep`（用 skani 与所有物种代表基因组比较）。
- **API：** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_ani_to_representatives(out_dir, *, genome_dir=None, batchfile=None, database=None, extension='fna', prefix='gtdbtk', min_af=None, cpus=1, gtdbtk_path=None)`。
- **输入：** 基因组和新的 `out_dir`；可选 `min_af`。
- **输出：** `ani_summary` 行：最近的代表基因组、ANI、比对分数、分类。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_ani_to_representatives

result = gtdbtk_ani_to_representatives("gtdb_ani", genome_dir="genomes", database="/data/gtdbtk/release226")
print(result["ani_summary"][0])
```

## 9) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_trim_msa()`

- **作用：** 用列掩码修剪未修剪的标记比对。
- **处理方法：** 运行 `gtdbtk trim_msa`；掩码为 `1` 的列被保留。`reference_mask`（`bac` 或 `arc`）使用参考包中的掩码。
- **API：** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_trim_msa(untrimmed_msa, output_path, *, mask_file=None, reference_mask=None, database=None, gtdbtk_path=None)`。
- **输入：** 比对文件、新的输出路径和 `mask_file` 或 `reference_mask`。
- **输出：** `output`、`sequence_count`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_trim_msa

result = gtdbtk_trim_msa("untrimmed.fasta", "trimmed.fasta", reference_mask="bac", database="/data/gtdbtk/release226")
print(result["sequence_count"])
```

## 10) `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_export_reference_msa()`

- **作用：** 导出某一域未修剪的 GTDB 参考标记比对。
- **处理方法：** 运行 `gtdbtk export_msa --domain bac|arc`。
- **API：** `dnakit.microbiome.gtdbtk_p15_cli_functions.gtdbtk_export_reference_msa(output_path, *, domain='bac', database=None, gtdbtk_path=None)`。
- **输入：** `domain` 和新的输出路径。
- **输出：** `output`、`sequence_count`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.microbiome.gtdbtk_p15_cli_functions import gtdbtk_export_reference_msa

result = gtdbtk_export_reference_msa("bac120_untrimmed.fasta", domain="bac", database="/data/gtdbtk/release226")
print(result["sequence_count"])
```

## 11) `dnakit.alignment.spoa_p15_cli_functions.spoa_partial_order_alignment()`

- **作用：** 把序列比对成部分序有向图，输出 MSA、一致序列或 GFA/DOT 格式的图，可选择得分更高的链。
- **处理方法：** 运行 `spoa`，参数为 `-r`（`result`：`msa`、`consensus_and_msa`、`consensus`、`graph`、`consensus_and_graph`）、`-d`（`dot_path`）、`-s`（`strand_ambiguous`）、`-l`（`algorithm`）以及打分选项。默认值与 SPOA 一致，因此默认是局部比对。
- **API：** `dnakit.alignment.spoa_p15_cli_functions.spoa_partial_order_alignment(sequences, output_path, *, result='msa', algorithm='local', match=5, mismatch=-4, gap_open=-8, gap_extend=-6, gap_open_second=-10, gap_extend_second=-4, strand_ambiguous=False, min_coverage=None, dot_path=None, spoa_path=None)`。
- **输入：** FASTA/FASTQ 文件、新的输出路径和上述选项。
- **输出：** `alignment` 记录、`consensus`、`alignment_length`；图则返回 `node_count`、`edge_count`、`path_count`；以及 `dot_file`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.alignment.spoa_p15_cli_functions import spoa_partial_order_alignment

msa = spoa_partial_order_alignment("s.fa", str(out / "msa.fa"), result="consensus_and_msa", algorithm="global", strand_ambiguous=True, dot_path=str(out / "g.dot"))
print(msa["consensus"], msa["alignment_length"])
for record in msa["alignment"]:
    print(record["name"], record["sequence"])
graph = spoa_partial_order_alignment("s.fa", str(out / "g.gfa"), result="consensus_and_graph", algorithm="global")
print(graph["node_count"], graph["edge_count"], graph["path_count"])
```

- **示例输出：**

```text
ACGTACGTACGTAAC 15
a ACGTACGTACGTAAC
b ACGTACCTACGTAAC
c ACGTACGTACG-AAC
d ACGTACGTACGTAAC
19 20 5
```
