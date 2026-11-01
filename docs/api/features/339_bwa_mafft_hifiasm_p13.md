# BWA、MAFFT 与 hifiasm 命令行封装

本组为 [BWA](https://github.com/lh3/bwa) 中 `bwa mem` 之外的功能提供封装（建索引各阶段、`aln`、`samse`/`sampe`、`fastmap`、`pemerge`、`maxk`，以及 bwa.kit 脚本 `bwa-postalt.js`、`typeHLA.sh`、`run-gen-ref` 和 `qualfa2fq.pl`；位于 `dnakit.alignment`）；为 [MAFFT](https://mafft.cbrc.jp/alignment/software/) 提供命名策略 `mafft_align(strategy=...)`、加入序列、种子、profile、分区树、稀疏核心和两两距离（位于 `dnakit.alignment`）；并为 [hifiasm](https://github.com/chhylp123/hifiasm) 提供单个函数 `hifiasm_assemble`（位于 `dnakit.assembly`），通过参数选择 HiFi、ONT、超长读长、Hi-C、trio、双倍体脚手架、去重、读段输出、断点、误连、端粒和嵌合体选项。`bwa maxk` 与 `bwa pemerge` 是 BWA 的隐藏/实验性子命令；`fa2pac`、`pac2bwt`、`pac2bwtgen`、`bwtupdate`、`bwt2sa` 是 `bwa index` 的各阶段，用 `bwa_index_stage(stage=...)` 选择。

程序需另行安装：`conda install -c bioconda bwa bwakit mafft hifiasm`（已用 BWA 0.7.19、bwa.kit 0.7.15、MAFFT 7.526、hifiasm 0.25.0 测试）。`qinsi`/`xinsi` 需要带 extensions 编译的 MAFFT；`mafft-sparsecore.rb` 需要 ruby；`--fastaparttree` 需要 fasta34。hifiasm 默认 Bloom 过滤器占用 16 GB 内存，小基因组请传 `bloom_filter_bits=0`。

## 示例数据

**示例文件：** `seqs.fa`

```text
>s0
ATGGCGTACGTTAGCTAGCTAGGATCGATCGGCTAAGCTTAGCGATCGATTAGCCGATCGGGATCTAGCTAGGCTA
>s1
ATGGCGTACGTTAGCTAGCTAGGATCGATCGGCTAAGCTTAGCGATCGATTAGCCGATCGGGATCTAGCTAGGCTA
>s2
ATGGCGTACGTAGCTAGCTAGGATCGATCGGCTAAGCTTAGCGATCGATTAGCCGATCGGGATCTAGCTAGGCTA
>s3
ATGGCGTACGTTAGCTAGCTAGGATCGATCGGCAAGCTTAGCGATCGATTAGCCGATCGGGATCTAGCTAGGCTA
```

**示例文件：** `ref.fa`

```text
>chr1
TTTTTTTATTTCTTGCACATAATAATAAGAATCTAACTTAATAACTTGTGAACCCTATTTAGCCTGGATTTGGTCTTTTTGAATATCATATTATGTCCCATTTACACGAAATCTGTTGATAATCAGTCTTTCGATCTCCGTGG
```

**示例文件：** `query.fa`

```text
>q1
CTTGTGAACCCTATTTAGCCTGGATTTGGTCTTTTTGAATATC
>q2
ACGTACGTA
```

## 1) `dnakit.alignment.bwa_index_stage()`

- **作用：** 手工运行 BWA 建索引的某一步：fa2pac、pac2bwt、pac2bwtgen、bwtupdate 或 bwt2sa。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.bwa_index_stage(stage, input_path, output_path, *, forward_only=False, sa_interval=32, bwa_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 2) `dnakit.alignment.bwa_aln()`

- **作用：** BWA-backtrack 短读长比对，输出后缀数组区间 .sai（bwa aln）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.bwa_aln(index_prefix, reads_path, output_path, *, max_diff=0.04, max_gap_opens=1, max_gap_extensions=-1, indel_end_skip=5, max_long_deletion_occurrences=10, seed_length=32, max_seed_diff=2, max_queue_entries=2000000, threads=1, mismatch_penalty=3, gap_open_penalty=11, gap_extension_penalty=4, max_best_hits=30, trim_quality=0, barcode_length=None, log_scaled_gap_penalty=False, non_iterative=False, illumina_1_3=False, bam_input=False, bam_read=None, filter_casava=False, bwa_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 3) `dnakit.alignment.bwa_sai_to_sam()`

- **作用：** 由 .sai 生成单端（samse）或双端（sampe）SAM。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.bwa_sai_to_sam(index_prefix, sai_paths, reads_paths, output_path, *, read_group=None, max_hits=None, max_insert_size=None, max_occurrences_per_end=None, max_discordant_hits=None, chimeric_rate=None, preload_index=False, disable_smith_waterman=False, disable_insert_size_estimate=False, bwa_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 4) `dnakit.alignment.bwa_fastmap()`

- **作用：** 查找序列相对 BWA 索引的超级最大精确匹配（bwa fastmap）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.bwa_fastmap(index_prefix, reads_path, output_path, *, min_smem_length=17, max_interval_size=20, min_interval_size=1, max_mem_length=2147483647, stop_interval_size=0, max_queries_returned=1000, bwa_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 5) `dnakit.alignment.bwa_pemerge()`

- **作用：** 合并重叠的双端读段（bwa pemerge）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.bwa_pemerge(read1_path, output_path, *, read2_path=None, output='all', threads=1, min_overlap=10, max_error_sum=70, bwa_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 6) `dnakit.alignment.bwa_maxk()`

- **作用：** 统计序列相对索引的最大精确匹配长度直方图（bwa maxk）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.bwa_maxk(bwt_path, sequences_path, output_path, *, self_match=False, bwa_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 7) `dnakit.alignment.bwakit_postalt()`

- **作用：** 对 ALT contig、decoy 和 HLA 基因的比对做后处理（bwa.kit bwa-postalt.js）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.bwakit_postalt(alt_path, alignments_path, output_path, *, hla_prefix=None, min_alt_fraction_for_mapq=1.0, bwakit_directory=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 8) `dnakit.alignment.bwakit_type_hla()`

- **作用：** 由匹配读段（或组装 contig）做单基因 HLA 分型（bwa.kit typeHLA.sh）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.bwakit_type_hla(prefix, gene, *, from_contigs=False, max_genotypes=50, bwakit_directory=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 9) `dnakit.alignment.bwakit_generate_reference()`

- **作用：** 下载并生成带 ALT、decoy、HLA 的人类参考分析集（bwa.kit run-gen-ref）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.bwakit_generate_reference(analysis_set, output_directory, *, bwakit_directory=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 10) `dnakit.alignment.bwa_qual_fasta_to_fastq()`

- **作用：** 把 FASTA 与质量值文件合并为 FASTQ（qualfa2fq.pl）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.bwa_qual_fasta_to_fastq(fasta_path, quality_path, output_path, *, script_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 11) `dnakit.alignment.mafft_align()`

- **作用：** 用 MAFFT 的命名策略（FFT-NS-i、NW-NS-i、L-INS-i、G-INS-i、E-INS-i、Q-INS-i、X-INS-i 等）做多序列比对。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.mafft_align(sequences_path, output_path, *, strategy='linsi', threads=1, sequence_type=None, gap_open_penalty=None, gap_extension_offset=None, max_iterations=None, retree=None, reorder=False, adjust_direction=None, output_format='fasta', mafft_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 12) `dnakit.alignment.mafft_add_sequences()`

- **作用：** 把新序列插入已有比对（mafft --add / --addfull / --addlong / --addfragments）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.mafft_add_sequences(alignment_path, new_sequences_path, output_path, *, mode='add', keep_length=False, six_mer_pair=False, threads=1, sequence_type=None, mafft_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 13) `dnakit.alignment.mafft_seed_alignment()`

- **作用：** 保持种子比对列不变，对其余序列做比对（mafft --seed）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.mafft_seed_alignment(sequences_path, seed_alignment_paths, output_path, *, threads=1, sequence_type=None, mafft_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 14) `dnakit.alignment.mafft_profile_alignment()`

- **作用：** 把两个已有比对作为 profile 互相比对（mafft-profile）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.mafft_profile_alignment(alignment1_path, alignment2_path, output_path, *, mafft_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 15) `dnakit.alignment.mafft_partition_tree_align()`

- **作用：** 用分区引导树对大量序列做比对（mafft --parttree 系列）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.mafft_partition_tree_align(sequences_path, output_path, *, tree_method='parttree', partition_size=None, retree=2, threads=1, sequence_type=None, mafft_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 16) `dnakit.alignment.mafft_sparsecore_align()`

- **作用：** 先比对稀疏核心序列再加入其余序列的大规模比对（mafft-sparsecore.rb）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.mafft_sparsecore_align(sequences_path, output_path, *, core_size=500, core_candidates='50%', seed=0, core_options=None, add_options=None, keep_input_order=False, mark_core=False, mafft_path=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 17) `dnakit.alignment.mafft_pairwise_distances()`

- **作用：** 由两两比对计算序列距离（mafft-distance）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.alignment.mafft_pairwise_distances(sequences_path, output_path, *, sequence_type=None, max_rows_returned=10000, mafft_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 18) `dnakit.assembly.hifiasm_assemble()`

- **作用：** 用 hifiasm 的任意模式做从头组装：HiFi、ONT、超长读长混合、Hi-C、trio、双倍体脚手架、单倍型去重及各类报告选项。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.assembly.hifiasm_assemble(reads_paths, output_prefix, *, ont=False, ultra_long_paths=None, hic_reads1_paths=None, hic_reads2_paths=None, paternal_kmers_path=None, maternal_kmers_path=None, paternal_reads_list_path=None, maternal_reads_list_path=None, dual_scaffold=False, scaffold_max_gap=None, purge_level=None, purge_similarity=None, purge_min_overlaps=None, haplotype_count=None, write_error_corrected_reads=False, write_overlaps_paf=False, low_quality_percent=None, break_low_coverage=None, break_high_coverage=None, break_exact_overlap_rate=None, misjoin_min_size=None, telomere_motif=None, chimeric_read_support=None, chimeric_flank_length=None, min_ont_read_length=None, min_ont_mean_quality=None, min_ultra_long_length=None, ultra_long_error_rate=None, primary_assembly=False, hom_coverage=None, genome_size=None, bloom_filter_bits=None, threads=1, hifiasm_path=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 示例

### 示例 1：用 L-INS-i 比对

- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.alignment import mafft_p13_cli_functions as mf

out = Path(tempfile.mkdtemp()) / "aln.fa"
result = mf.mafft_align("seqs.fa", str(out), strategy="linsi", verbose=False)
print(result["status"], result["sequence_count"])
widths = {len(line) for line in out.read_text().split(">")[1:] for line in [line.split("\n", 1)[1].replace("\n", "")]}
print(len(widths))
```

- **示例输出：**

```text
completed 4
1
```

### 示例 2：用 BWA fastmap 查找超级最大精确匹配

- **示例代码：**

```python
import subprocess
import tempfile
from pathlib import Path
from dnakit.alignment import bwa_p13_cli_functions as bw

subprocess.run(["bwa", "index", "ref.fa"], capture_output=True, check=True)
out = Path(tempfile.mkdtemp()) / "fastmap.txt"
result = bw.bwa_fastmap("ref.fa", "query.fa", str(out), verbose=False)
for query in result["queries"]:
    print(query["query"], query["length"], query["matches"])
```

- **示例输出：**

```text
q1 43 [{'start': 0, 'end': 43, 'hits': 1, 'positions': ['chr1:+45']}]
q2 9 []
```

### 示例 3：用 hifiasm 组装合成 HiFi 读段

- **示例代码：**

```python
import random
import tempfile
from pathlib import Path
from dnakit.assembly import hifiasm_p13_cli_functions as hf

rng = random.Random(1)
genome = "".join(rng.choice("ACGT") for _ in range(60000))
work = Path(tempfile.mkdtemp())
with (work / "hifi.fq").open("w") as handle:
    for i in range(180):
        start = rng.randrange(0, len(genome) - 8000)
        handle.write(f"@r{i}\n{genome[start:start + 8000]}\n+\n{'I' * 8000}\n")
result = hf.hifiasm_assemble([str(work / "hifi.fq")], str(work / "asm"), bloom_filter_bits=0, write_error_corrected_reads=True, verbose=False)
print(result["status"], result["version"])
print(any(item["total_bp"] > 50000 for item in result["contig_summary"]), bool(result["error_corrected_reads"]))
```

- **示例输出：**

```text
completed 0.25.0
True True
```
