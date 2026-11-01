# BEDTools 区间操作与 SURVIVOR 结构变异工具

本组为 [BEDTools](https://github.com/arq5x/bedtools2)（`window`、`random`、`shuffle`、`sample`、`multiinter`、`unionbedg`、`pairtopair`、`multicov`、`reldist`、`coverage` 直方图与深度、`bamtobed -bedpe`、`split`；位于 `dnakit.annotation`）和 [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)（模拟、评估、过滤、统计、MUMmer 注释与格式转换；位于 `dnakit.variants`）提供命令行封装。仅参数不同的条目合并为一个函数：`bedtools_window(report="count")`、`bedtools_reldist(detail=True)`、`bedtools_coverage_profile(report="depth")`、`bedtools_split(algorithm="simple")`、`survivor_convert_to_sv_vcf(input_format=...)`。SURVIVOR 的 `select_sample` 在 1.0.7 中会段错误，因此未封装。

程序需另行安装：`conda install -c bioconda bedtools survivor`（已用 BEDTools 2.27.1 和 2.31.1、SURVIVOR 1.0.7 测试）。未安装时函数仍可导入，调用时抛出带该提示的 `BackendUnavailableError`。

## 示例数据

**示例文件：** `a.bed`

```text
chr1	100	200	a1	5	+
chr1	500	600	a2	6	-
chr2	10	90	a3	7	+
```

**示例文件：** `b.bed`

```text
chr1	150	250	b1	1	+
chr1	900	950	b2	2	-
chr2	100	150	b3	3	+
```

**示例文件：** `sv.bed`

```text
chr1	1000	2000
chr1	5000	5600
```

## 1) `dnakit.annotation.bedtools_window()`

- **作用：** 在每个 A 区间周围的窗口内报告 B 特征，可输出重叠、唯一、计数或无重叠记录（bedtools window）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_window(a_path, b_path, output_path, *, window=None, left=None, right=None, strand_aware_window=False, strand=None, report='overlaps', header=False, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 2) `dnakit.annotation.bedtools_random()`

- **作用：** 在基因组上生成固定长度的随机区间（bedtools random）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_random(genome_path, output_path, *, length=100, count=1000000, seed=None, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 3) `dnakit.annotation.bedtools_shuffle()`

- **作用：** 在基因组中随机置换区间位置，可排除或限定区域（bedtools shuffle）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_shuffle(intervals_path, genome_path, output_path, *, exclude_path=None, include_path=None, same_chromosome=False, seed=None, max_exclude_overlap=None, choose_chromosome_first=False, bedpe=False, max_tries=None, no_overlapping=False, allow_beyond_chrom_end=False, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 4) `dnakit.annotation.bedtools_sample()`

- **作用：** 用水塘抽样从区间或 BAM 文件中随机取样记录（bedtools sample）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_sample(input_path, output_path, *, count=1000000, seed=None, strand=None, header=False, bam_as_bed=False, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 5) `dnakit.annotation.bedtools_multiinter()`

- **作用：** 找出多个已排序区间文件的共同区间（bedtools multiinter）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_multiinter(input_paths, output_path, *, names=None, header=False, cluster=False, genome_path=None, report_empty=False, filler=None, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 6) `dnakit.annotation.bedtools_unionbedg()`

- **作用：** 把多个 bedGraph 信号文件合并为一张表（bedtools unionbedg）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_unionbedg(input_paths, output_path, *, names=None, header=False, genome_path=None, report_empty=False, filler=None, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 7) `dnakit.annotation.bedtools_pairtopair()`

- **作用：** 报告两个 BEDPE 文件之间的重叠（bedtools pairtopair）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_pairtopair(a_path, b_path, output_path, *, min_overlap=None, overlap_type='both', slop=0, strand_aware_slop=False, ignore_strand=False, require_different_names=False, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 8) `dnakit.annotation.bedtools_multicov()`

- **作用：** 统计多个 BAM 文件落在每个区间内的比对数（bedtools multicov）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_multicov(bam_paths, intervals_path, output_path, *, split=False, strand=None, min_overlap=None, reciprocal=False, min_mapping_quality=None, include_duplicates=False, include_failed_qc=False, proper_pairs_only=False, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 9) `dnakit.annotation.bedtools_reldist()`

- **作用：** 计算两组区间的相对距离分布，或逐区间详细结果（bedtools reldist）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_reldist(a_path, b_path, output_path, *, detail=False, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 10) `dnakit.annotation.bedtools_coverage_profile()`

- **作用：** 输出 A 区间被 B 覆盖的直方图或逐位置深度（bedtools coverage -hist / -d）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_coverage_profile(a_path, b_path, output_path, *, report='histogram', strand=None, min_overlap_a=None, min_overlap_b=None, reciprocal=False, either_fraction=False, split=False, sorted_input=False, genome_path=None, header=False, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 11) `dnakit.annotation.bedtools_bam_to_bedpe()`

- **作用：** 把配对 BAM 比对转换为 BEDPE（bedtools bamtobed -bedpe）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_bam_to_bedpe(bam_path, output_path, *, mate1_first=False, edit_distance_score=False, bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 12) `dnakit.annotation.bedtools_split()`

- **作用：** 把 BED 文件按碱基总量或记录数拆分为 N 个文件（bedtools split）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.annotation.bedtools_split(input_path, output_directory, number_of_files, *, prefix='split', algorithm='size', bedtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 13) `dnakit.variants.survivor_simulate_sv_genome()`

- **作用：** 向参考基因组模拟 SV 和 SNP（SURVIVOR simSV）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_simulate_sv_genome(reference_path, output_prefix, *, duplication=(100, 10000, 3), indel=(20, 500, 1), translocation=(1000, 3000, 2), inversion=(600, 800, 4), inversion_deletion=(600, 800, 2), inversion_duplication=(600, 800, 2), snp_frequency=0.0, coordinates_for='simulated_reads', parameter_file=None, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 14) `dnakit.variants.survivor_scan_read_error_profile()`

- **作用：** 从长读长比对中提取逐位置错误谱（SURVIVOR scanreads）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_scan_read_error_profile(alignments_path, output_path, *, min_read_length=1000, reference_path=None, survivor_path=None, samtools_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 15) `dnakit.variants.survivor_simulate_long_reads()`

- **作用：** 依据错误谱从参考序列模拟 PacBio/ONT 长读长（SURVIVOR simreads）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_simulate_long_reads(reference_path, error_profile_path, output_path, *, coverage=2.0, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 16) `dnakit.variants.survivor_evaluate_sv_calls()`

- **作用：** 对照模拟真值评估 SV 检出结果（SURVIVOR eval）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_evaluate_sv_calls(vcf_path, simulated_bed_path, output_prefix, *, max_distance=10, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 17) `dnakit.variants.survivor_filter_sv_vcf()`

- **作用：** 按大小、区域、等位基因频率和读段支持过滤 SV VCF（SURVIVOR filter）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_filter_sv_vcf(vcf_path, output_path, *, ignore_regions_bed=None, min_size=None, max_size=None, min_allele_frequency=0.0, min_read_support=None, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 18) `dnakit.variants.survivor_sv_vcf_stats()`

- **作用：** 统计 SV VCF 的类型与长度分布（SURVIVOR stats）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_sv_vcf_stats(vcf_path, summary_path, *, min_size=None, max_size=None, min_read_support=None, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 19) `dnakit.variants.survivor_annotate_sv_with_mummer()`

- **作用：** 用 MUMmer show-diff 断点注释 SV VCF（SURVIVOR compMUMMer）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_annotate_sv_with_mummer(vcf_path, show_diff_paths, output_path, *, max_distance=100, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 20) `dnakit.variants.survivor_sv_vcf_to_bed()`

- **作用：** 把 SV VCF 转换为 BEDPE 风格的 BED（SURVIVOR vcftobed）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_sv_vcf_to_bed(vcf_path, output_path, *, min_size=None, max_size=None, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 21) `dnakit.variants.survivor_convert_to_sv_vcf()`

- **作用：** 把 BED、Bionano SMAP、Lumpy BEDPE 或 Assemblytics 结果转换为 SV VCF。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_convert_to_sv_vcf(input_path, output_path, *, input_format='bed', sv_type=None, min_size=0, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 22) `dnakit.variants.survivor_hapcut_to_phased_vcf()`

- **作用：** 由 HapCUT2 单倍型结果和 SNP VCF 生成定相 VCF（SURVIVOR hapcuttovcf）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_hapcut_to_phased_vcf(snp_vcf_path, hapcut_path, output_path, *, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 23) `dnakit.variants.survivor_low_mapping_quality_regions()`

- **作用：** 把 MQ0 覆盖度聚类为 BED 区域用于 SV 过滤（SURVIVOR bincov）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_low_mapping_quality_regions(coverage_path, output_path, *, border_size=10, min_coverage=2, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 24) `dnakit.variants.survivor_pairwise_overlap_matrix()`

- **作用：** 计算合并多样本 VCF 的样本两两 SV 重叠矩阵（SURVIVOR genComp）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_pairwise_overlap_matrix(merged_vcf_path, output_path, *, normalize=False, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 25) `dnakit.variants.survivor_parental_phasing()`

- **作用：** 利用亲本 VCF、HapCUT2 和 GATK 结果对 SV 做亲本定相（SURVIVOR parent_phasing）。
- **处理方法：** 以固定白名单参数（不经 shell）调用外部程序，或按脚本逻辑直接实现；输出文件必须是新文件，绝不覆盖已有文件；程序缺失时抛出带安装提示的 `BackendUnavailableError`。
- **API：** `dnakit.variants.survivor_parental_phasing(merged_parental_vcf_path, hapcut_path, gatk_path, proband_snp_vcf_path, output_path, *, survivor_path=None, timeout_seconds=3600.0)`。
- **输入：** 已存在的输入文件和签名中的新输出路径；可选参数仅在给出时才传给程序，因此沿用程序自身默认值。
- **输出：** 可 JSON 序列化的字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`，以及文档字符串中描述的解析结果。

## 示例

### 示例 1：统计 A 区间附近的 B 特征数并输出逐位置深度

- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.annotation import bedtools_p13_cli_functions as bt

out = Path(tempfile.mkdtemp())
result = bt.bedtools_window("a.bed", "b.bed", str(out / "w.bed"), window=100, report="count", verbose=False)
print(result["status"], result["record_count"])
print((out / "w.bed").read_text().splitlines()[0])
depth = bt.bedtools_coverage_profile("a.bed", "b.bed", str(out / "d.txt"), report="depth", verbose=False)
print(depth["report"], depth["record_count"])
```

- **示例输出：**

```text
completed 3
chr1	100	200	a1	5	+	1
depth 280
```

### 示例 2：BED 转 SV VCF，再过滤并统计

- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.variants import survivor_p13_cli_functions as sv

out = Path(tempfile.mkdtemp())
conv = sv.survivor_convert_to_sv_vcf("sv.bed", str(out / "sv.vcf"), input_format="bed", sv_type="DEL", verbose=False)
print(conv["status"], conv["record_count"])
flt = sv.survivor_filter_sv_vcf(str(out / "sv.vcf"), str(out / "big.vcf"), min_size=800, verbose=False)
print(flt["records_in"], flt["records_out"])
stats = sv.survivor_sv_vcf_stats(str(out / "sv.vcf"), str(out / "stats.txt"), verbose=False)
print(stats["totals"])
```

- **示例输出：**

```text
completed 2
2 1
{'Tot': 2, 'DEL': 2, 'DUP': 0, 'INS': 0, 'INV': 0, 'TRA': 0}
```
