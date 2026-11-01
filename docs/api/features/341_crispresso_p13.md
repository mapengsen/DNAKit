# CRISPResso2 基因组编辑分析封装与辅助脚本移植

本组把 [CRISPResso2](https://github.com/pinellolab/CRISPResso2) 加入 `dnakit.design_crispr`：`crispresso_quantify_editing` 运行 `CRISPResso`，除普通 indel 定量外，通过参数覆盖 HDR、编码序列移码、碱基编辑和 prime editing 分析；`crispresso_batch`、`crispresso_pooled`、`crispresso_wgs`、`crispresso_compare`、`crispresso_pooled_wgs_compare`、`crispresso_aggregate` 运行其他命令行程序；另有六个函数移植仓库中的辅助脚本（`filterReadsOnSequencePresence.py`、`countHighQualityBases.py`、`count_sgRNA_specific_edits.py`、`reconstituteReads.py`、`plotAmbiguous.py`、`plotCustomAllelePlot.py`）。2.3.4 中 `CRISPRessoMeta` 没有命令行入口，未封装。

需另行安装 CRISPResso2：`conda install -c bioconda crispresso2`（已用 2.3.4 测试；会同时安装 bowtie2、samtools、fastp）。读取已完成分析或绘图的四个脚本移植函数会导入 `CRISPResso2` Python 包（`pip install "dnakit[crispresso2]"`）。

## 示例数据

**示例文件：** `reads.fq`

```text
@a
ACGTGGCACTGCGGCTGGAGGTGGAAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIII
@b
ACGTACGTACGTACGTACGTACGTAAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIII
```

## 1) `dnakit.design_crispr.crispresso_quantify_editing()`

- **作用：** 从测序读段定量一个或多个扩增子上的基因组编辑结果，覆盖 HDR、编码序列移码、碱基编辑和 prime editing（CRISPResso）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_quantify_editing(fastq_r1_path, output_folder, *, amplicon_sequences, fastq_r2_path=None, amplicon_names=None, guide_sequences=None, guide_names=None, expected_hdr_amplicon_sequence=None, coding_sequences=None, base_editor_output=False, conversion_nuc_from=None, conversion_nuc_to=None, prime_editing_pegrna_spacer_seq=None, prime_editing_pegrna_extension_seq=None, prime_editing_pegrna_scaffold_seq=None, prime_editing_nicking_guide_seq=None, prime_editing_extension_quantification_window_size=None, prime_editing_scaffold_min_match_length=None, name=None, quantification_window_size=None, quantification_window_center=None, quantification_window_coordinates=None, exclude_bp_from_left=None, exclude_bp_from_right=None, min_average_read_quality=None, min_single_bp_quality=None, min_bp_quality_or_n=None, min_identity_score=None, ignore_substitutions=False, ignore_insertions=False, ignore_deletions=False, discard_indel_reads=False, trim_sequences=False, write_detailed_allele_table=False, fastq_output=False, bam_output=False, vcf_output=False, suppress_plots=False, suppress_report=False, keep_intermediate=False, n_processes=None, crispresso_path=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 2) `dnakit.design_crispr.crispresso_batch()`

- **作用：** 按设置文件对多个样本批量运行 CRISPResso（CRISPRessoBatch）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_batch(batch_settings_path, output_folder, *, batch_name=None, amplicon_sequences=None, guide_sequences=None, quantification_window_size=None, suppress_plots=False, suppress_batch_summary_plots=False, n_processes=None, crispresso_path=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 3) `dnakit.design_crispr.crispresso_pooled()`

- **作用：** 混合扩增子测序的拆分与定量（CRISPRessoPooled）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_pooled(fastq_r1_path, output_folder, *, amplicons_file_path=None, genome_index=None, gene_annotations_path=None, fastq_r2_path=None, name=None, min_reads_to_use_region=None, demultiplex_only_at_amplicons=False, guide_sequences=None, quantification_window_size=None, suppress_plots=False, n_processes=None, crispresso_path=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 4) `dnakit.design_crispr.crispresso_wgs()`

- **作用：** 对全基因组 BAM 的指定区域定量编辑（CRISPRessoWGS）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_wgs(bam_path, region_file_path, reference_fasta_path, output_folder, *, name=None, gene_annotations_path=None, min_reads_to_use_region=None, guide_sequences=None, quantification_window_size=None, suppress_plots=False, n_processes=None, crispresso_path=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 5) `dnakit.design_crispr.crispresso_compare()`

- **作用：** 比较两次 CRISPResso 分析的编辑结果（CRISPRessoCompare）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_compare(folder1, folder2, output_folder, *, sample_names=None, name=None, reported_qvalue_cutoff=None, crispresso_path=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 6) `dnakit.design_crispr.crispresso_pooled_wgs_compare()`

- **作用：** 比较两次 Pooled 或 WGS 分析（CRISPRessoPooledWGSCompare）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_pooled_wgs_compare(folder1, folder2, output_folder, *, sample_names=None, name=None, reported_qvalue_cutoff=None, crispresso_path=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 7) `dnakit.design_crispr.crispresso_aggregate()`

- **作用：** 汇总多次 CRISPResso 运行的编辑摘要（CRISPRessoAggregate）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_aggregate(search_directory, name, *, folder_prefixes=None, folder_suffix=None, min_reads_for_inclusion=None, suppress_plots=False, crispresso_path=None, timeout_seconds=86400.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 8) `dnakit.design_crispr.crispresso_filter_reads_by_sequence()`

- **作用：** 按是否含有指定序列过滤 FASTQ 读段（filterReadsOnSequencePresence.py 的移植）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_filter_reads_by_sequence(fastq_r1_path, output_r1_path, *, fastq_r2_path=None, output_r2_path=None, include_sequences=(), exclude_sequences=())`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 9) `dnakit.design_crispr.crispresso_count_high_quality_bases()`

- **作用：** 统计 CRISPResso 注释 FASTQ 在某参考位置的高/低质量碱基（countHighQualityBases.py 的移植）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_count_high_quality_bases(annotated_fastq_path, position, original_base, target_base, *, reference_name='Reference', quality_cutoff=30)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 10) `dnakit.design_crispr.crispresso_count_sgrna_specific_edits()`

- **作用：** 统计 CRISPResso 结果中每个 sgRNA 窗口内的编辑（count_sgRNA_specific_edits.py 的移植）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_count_sgrna_specific_edits(crispresso_folder)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 11) `dnakit.design_crispr.crispresso_reconstitute_reads()`

- **作用：** 由等位基因频率表重建 FASTQ 读段（reconstituteReads.py 的移植）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_reconstitute_reads(crispresso_folder, output_fastq_path)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 12) `dnakit.design_crispr.crispresso_plot_ambiguous_alleles()`

- **作用：** 绘制 CRISPResso 结果中歧义等位基因的等位基因表（plotAmbiguous.py 的移植）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_plot_ambiguous_alleles(crispresso_folder, output_root, *, min_frequency=None, max_rows=None, plot_cut_point=True, save_png=False)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 13) `dnakit.design_crispr.crispresso_plot_allele_window()`

- **作用：** 绘制切割位点周围自定义窗口的等位基因表（plotCustomAllelePlot.py 的移植）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.design_crispr.crispresso_plot_allele_window(crispresso_folder, output_root, *, plot_left=20, plot_right=20, plot_center=None, min_frequency=None, max_rows=None, plot_cut_point=True, save_png=False)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 示例

### 示例 1：用合成读段定量编辑

- **示例代码：**

```python
import random
import tempfile
from pathlib import Path
from dnakit.design_crispr import crispresso_p13_cli_functions as cr

rng = random.Random(11)
guide = "GGCACTGCGGCTGGAGGTGG"
amplicon = "".join(rng.choice("ACGT") for _ in range(60)) + guide + "GGG" + "".join(rng.choice("ACGT") for _ in range(70))
cut = amplicon.index(guide) + 17
work = Path(tempfile.mkdtemp())
with (work / "r1.fq").open("w") as handle:
    for i in range(200):
        seq = amplicon if rng.random() < 0.6 else amplicon[:cut] + amplicon[cut + 3:]
        handle.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
result = cr.crispresso_quantify_editing(str(work / "r1.fq"), str(work / "out"), amplicon_sequences=amplicon, guide_sequences=guide, suppress_plots=True, suppress_report=True, verbose=False)
row = result["quantification"][0]
print(result["status"], result["version"])
print(row["Reads_in_input"], row["Unmodified"], row["Modified"], row["Deletions"])
```

- **示例输出：**

```text
completed 2.3.4
200 107 93 93
```

### 示例 2：保留含有指定向导序列的读段

- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.design_crispr import crispresso_p13_cli_functions as cr

out = Path(tempfile.mkdtemp()) / "kept.fq"
result = cr.crispresso_filter_reads_by_sequence("reads.fq", str(out), include_sequences=["GGCACTGCGGCTGGAGGTGG"], verbose=False)
print(result["reads_read"], result["reads_written"])
print(out.read_text().splitlines()[0])
```

- **示例输出：**

```text
2 1
@a
```
