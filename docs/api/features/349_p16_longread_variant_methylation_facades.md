# 长读长质控、组装、结构变异、深度、甲基化与 FastTree 封装

这些函数以固定参数白名单、不经 shell 运行 [FreeBayes](https://github.com/freebayes/freebayes) 辅助工具、[Filtlong](https://github.com/rrwick/Filtlong)、[NanoFilt](https://github.com/wdecoster/nanofilt)、[Porechop](https://github.com/rrwick/Porechop)、[mosdepth](https://github.com/brentp/mosdepth)、[Sniffles2](https://github.com/fritzsedlazeck/Sniffles)、[cuteSV](https://github.com/tjiangHIT/cuteSV)、[vcfanno](https://github.com/brentp/vcfanno)、[Flye](https://github.com/mikolmogorov/Flye)、[Liftoff](https://github.com/agshumate/Liftoff)、[miniprot](https://github.com/lh3/miniprot)、[MethylDackel](https://github.com/dpryan79/MethylDackel) 和 [FastTree](http://www.microbesonline.org/fasttree/)。FreeBayes 的区域/VCF 脚本和 NanoFilt 时间过滤器用 Python 按原脚本语义重新实现（已与原脚本对照），不需要外部程序。

DNAKit 不安装这些程序；请放到 `PATH` 或传入 `*_path` 参数。输出写入新路径。测试版本：mosdepth 0.3.14、FastTree 2.2.0、miniprot 0.18、Filtlong 0.3.1、NanoFilt 2.8.0、Sniffles2 2.6.2、cuteSV 2.1.4、vcfanno 0.3.9、Flye 2.9.6、Liftoff 1.6.3、MethylDackel 0.6.1、bamleftalign（FreeBayes 1.3.10）。Porechop 已归档且其 conda 构建在本环境无法运行，它的命令组装由使用记录型替身的单元测试覆盖。

## 示例数据

示例运行前先写入这里声明的输入文件；其余示例在代码中自行生成数据。

**示例文件：** `ref.fa.fai`

```text
chr1	1000	6	60	61
chr2	450	1100	60	61
tiny	7	2000	60	61
```

**示例文件：** `a.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	10	.	A	T	.	.	.
```

**示例文件：** `b.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	20	.	G	C	.	.	.
```

**示例文件：** `aln.fa`

```text
>t1
ACGTACGTAGCTAGCTAGGATCGATCGTAGCTAGCATCGATCGAT
>t2
ACGTACGTAGCTAGCTAGGATCGATCGTAGCTAGCATCGATCGAT
>t3
ACGTACGTAGCAAGCTAGGATCGATCGTAGCTAGCATCGTTCGAT
>t4
ACGTTCGTAGCAAGCTAGGATCGATCGAAGCTAGCATCGTTCGAT
>t5
TCGTTCGTAGCAAGCTAGGATCGTTCGAAGCTAGCATCGTTCGAC
```

**示例文件：** `reads.fastq`

```text
@low runid=x start_time=2024-01-01T10:00:00Z
ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@mid runid=x start_time=2024-01-02T10:00:00Z
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@high runid=x start_time=2024-01-03T10:00:00Z
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## 1) `dnakit.variants.freebayes_p16_cli_functions.keep_first_vcf_header()`

- **作用：** 拼接 VCF 文件并只保留第一个头部（FreeBayes `vcffirstheader`）。
- **处理方法：** 按顺序读取文件，写出第一个 `#CHROM` 行之前的元信息和该行，丢弃之后的所有头部行，数据行去掉行尾空白。
- **API：** `dnakit.variants.freebayes_p16_cli_functions.keep_first_vcf_header(inputs, output_path, verbose=True)`。
- **输入：** VCF 文件列表（可为 `.gz`）和新的输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `header_lines_written`、`header_lines_dropped`、`records_written`。
- **示例见：** `partition_reference_into_regions()`

## 2) `dnakit.variants.freebayes_p16_cli_functions.left_align_bam_indels()`

- **作用：** 左对齐并合并 BAM 中的插入缺失（FreeBayes `bamleftalign`）。
- **处理方法：** BAM 从标准输入送入 `bamleftalign`，写出的 BAM 保存到新文件；反复左移直到每条比对稳定。
- **API：** `dnakit.variants.freebayes_p16_cli_functions.left_align_bam_indels(alignment, output_path, *, reference, compressed=True, max_iterations=None, bamleftalign_path=None, verbose=True)`。
- **输入：** BAM、参考 FASTA（需 `.fai`）和新的输出 BAM。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `alignment`、`output`。
- **示例见：** `compute_mosdepth_depth()`

## 3) `dnakit.variants.freebayes_p16_cli_functions.partition_reference_by_bam_size()`

- **作用：** 按 BAM 索引中的数据量把参考序列切成数据量相等的区域（FreeBayes `split_ref_by_bai_datasize.py`）。
- **处理方法：** 读取每个 `<bam>.bai` 的线性索引，把各 BAM 的区间数据量相加，用线性插值在累计数据量的整数倍处切分；与原脚本一样把 BGZF 虚拟偏移之差当作数据量。
- **API：** `dnakit.variants.freebayes_p16_cli_functions.partition_reference_by_bam_size(reference, bam_files, *, output_path=None, target_data_size=100000000.0, bai_interval_size=16384, verbose=True)`。
- **输入：** FASTA 或 `.fai`、带 `.bai` 的 BAM 列表。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `region_count`、`regions`。

## 4) `dnakit.variants.freebayes_p16_cli_functions.partition_reference_by_coverage()`

- **作用：** 按每位点深度把参考序列切成总覆盖度近似相等的区域（FreeBayes `coverage_to_regions.py`）。
- **处理方法：** 累计 `bamtools coverage` 的逐位点深度，累计和超过 总深度/区域数 时结束一个区域，每条序列的最后一个区域延伸到序列末端。
- **API：** `dnakit.variants.freebayes_p16_cli_functions.partition_reference_by_coverage(reference, coverage, num_regions, *, output_path=None, verbose=True)`。
- **输入：** FASTA 或 `.fai`、三列 `chrom pos depth` 覆盖度文件、区域数。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `region_count`、`regions`、`total_depth`。

## 5) `dnakit.variants.freebayes_p16_cli_functions.partition_reference_into_regions()`

- **作用：** 把参考序列切成固定长度的 `chrom:start-end` 区域（FreeBayes `fasta_generate_regions.py`）。
- **处理方法：** 读取 `.fai` 索引，按 `region_size` 或 `chunks` 等分各序列，可限定染色体，也可每个区域写一个 BED 文件。区域为 0 起始、左闭右开，与原脚本一致。
- **API：** `dnakit.variants.freebayes_p16_cli_functions.partition_reference_into_regions(reference, region_size, *, output_path=None, chunks=False, chromosomes=None, bed_prefix=None, verbose=True)`。
- **输入：** FASTA 路径或其 `.fai`（索引必须存在）和区域大小。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `region_count`、`regions`。
- **示例代码：**

```python
from dnakit.variants.freebayes_p16_cli_functions import keep_first_vcf_header, partition_reference_into_regions

print(partition_reference_into_regions("ref.fa", 300, verbose=False)["regions"])
print(partition_reference_into_regions("ref.fa", 3, chunks=True, chromosomes=["chr1"], verbose=False)["regions"])
result = keep_first_vcf_header(["a.vcf", "b.vcf"], "merged.vcf", verbose=False)
print(result["header_lines_dropped"], result["records_written"])
```

- **示例输出：**

```text
['chr1:0-300', 'chr1:300-600', 'chr1:600-900', 'chr1:900-1000', 'chr2:0-300', 'chr2:300-450', 'tiny:0-7']
['chr1:0-334', 'chr1:334-668', 'chr1:668-1000']
2 2
```

## 6) `dnakit.sequencing_qc.p16_readqc_cli_functions.demultiplex_nanopore_barcodes()`

- **作用：** 按条形码拆分 reads 并去接头（`porechop --barcode_dir`）。
- **处理方法：** 每个条形码一个 FASTQ，未分配的放入 `none.fastq`；可要求两端条形码。
- **API：** `dnakit.sequencing_qc.p16_readqc_cli_functions.demultiplex_nanopore_barcodes(reads, output_directory, *, require_two_barcodes=False, barcode_threshold=None, discard_unassigned=False, threads=1, porechop_path=None, verbose=True)`。
- **输入：** FASTQ 与新输出目录。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `reads_per_bin`。

## 7) `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_long_reads_with_reference()`

- **作用：** 用参考 k-mer 匹配给长读长打分，并过滤、修剪或拆分（`filtlong`）。
- **处理方法：** 参考可为 Illumina reads 或组装；`trim` 去掉两端不匹配碱基，`split` 在连续不匹配处拆分，另有长度、质量、碱基量阈值。
- **API：** `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_long_reads_with_reference(reads, output_path, *, illumina_reads_1=None, illumina_reads_2=None, assembly=None, trim=False, split=None, min_length=None, max_length=None, min_mean_quality=None, target_bases=None, keep_percent=None, filtlong_path=None, verbose=True)`。
- **输入：** 长读长 FASTQ、新输出路径，可选参考。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `reads_kept`。

## 8) `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_nanopore_reads_by_gc_summary()`

- **作用：** 按 GC 含量或测序摘要中的质量过滤长读长（`NanoFilt`）。
- **处理方法：** `--minGC`/`--maxGC` 为 0–1 的比例；`--summary` 使用 ONT `sequencing_summary.txt` 的质量。
- **API：** `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_nanopore_reads_by_gc_summary(reads, output_path, *, min_gc=None, max_gc=None, sequencing_summary=None, min_quality=None, min_length=None, max_length=None, read_type=None, nanofilt_path=None, verbose=True)`。
- **输入：** 未压缩 FASTQ、新输出路径，至少一个过滤条件。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `reads_kept`。
- **示例代码：**

```python
import gzip
from dnakit.sequencing_qc.p16_readqc_cli_functions import filter_nanopore_reads_by_gc_summary, filter_nanopore_reads_by_start_time

result = filter_nanopore_reads_by_gc_summary("reads.fastq", "gc.fastq", min_gc=0.3, max_gc=0.7, verbose=False)
print(result["reads_kept"], open("gc.fastq").readline().split()[0])
with gzip.open("reads.fastq.gz", "wt") as handle:
    handle.write(open("reads.fastq").read())
timed = filter_nanopore_reads_by_start_time("reads.fastq.gz", "timed.fastq", time_from="2024-01-02T00:00:00Z", time_to="2024-01-03T00:00:00Z", verbose=False)
print(timed["reads_total"], timed["reads_kept"])
```

- **示例输出：**

```text
1 @mid
3 1
```

## 9) `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_nanopore_reads_by_start_time()`

- **作用：** 按测序开始时间窗口保留 reads（NanoFilt `scripts/timefilt.py`）。
- **处理方法：** 纯 Python 实现，读取头部的 `start_time=` 并按含边界的时间窗过滤。
- **API：** `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_nanopore_reads_by_start_time(reads_gz, output_path, *, time_from, time_to, verbose=True)`。
- **输入：** gzip FASTQ、起止时间（含时区）与新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `reads_total`、`reads_kept`。
- **示例见：** `filter_nanopore_reads_by_gc_summary()`

## 10) `dnakit.sequencing_qc.p16_readqc_cli_functions.trim_nanopore_adapters_discard_middle()`

- **作用：** 去除接头并丢弃中间含接头的 reads（`porechop --discard_middle`）。
- **处理方法：** 中间含接头的 reads 为嵌合体；默认拆分，这里默认丢弃。
- **API：** `dnakit.sequencing_qc.p16_readqc_cli_functions.trim_nanopore_adapters_discard_middle(reads, output_path, *, discard_middle=True, threads=1, porechop_path=None, verbose=True)`。
- **输入：** FASTQ 与新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `discard_middle`。

## 11) `dnakit.variants.p16_variant_cli_functions.annotate_vcf_from_bam_vcfanno()`

- **作用：** 用 BAM 中的读深度、MAPQ 和碱基序列注释 VCF 变异（`vcfanno`）。
- **处理方法：** 自动写出只含一个 BAM `[[annotation]]` 的 TOML 配置（字段 depth/mapq/seq，操作 count/mean/concat）并运行 vcfanno。
- **API：** `dnakit.variants.p16_variant_cli_functions.annotate_vcf_from_bam_vcfanno(vcf, alignment, output_path, *, operations=('depth', 'mapq', 'seq'), threads=1, vcfanno_path=None, verbose=True)`。
- **输入：** VCF、已索引 BAM、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `configuration`、`records`。

## 12) `dnakit.variants.p16_variant_cli_functions.call_structural_variants_sniffles2()`

- **作用：** 用 Sniffles2 从长读长比对检测结构变异，可检测嵌合/体细胞、定相并保存 SNF。
- **处理方法：** `mosaic` 检测低频 SV，`phase` 使用 HP/PS 标签定相，`snf_path` 保存候选以便多样本联合调用。
- **API：** `dnakit.variants.p16_variant_cli_functions.call_structural_variants_sniffles2(alignment, vcf_path, *, reference=None, snf_path=None, mosaic=False, phase=False, min_sv_length=None, min_support=None, sample_id=None, threads=1, sniffles_path=None, verbose=True)`。
- **输入：** 比对后的 BAM、新 VCF 路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `mosaic`、`phase`、`records`。

## 13) `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_depth()`

- **作用：** 计算 BAM/CRAM 的逐位点深度、全基因组深度分布和每条染色体摘要（`mosdepth`）。
- **处理方法：** 写出 `per-base.bed.gz`、`global.dist.txt` 与 `summary.txt`，摘要解析为 `summary`。
- **API：** `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_depth(alignment, prefix, *, per_base=True, chromosome=None, threads=0, mapq=0, min_fragment=None, max_fragment=None, fast_mode=False, fragment_mode=False, exclude_flag=None, include_flag=None, reference=None, read_groups=None, mosdepth_path=None, verbose=True)`。
- **输入：** 排序并索引的 BAM/CRAM 与输出前缀。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `summary`、`distribution_head`。
- **示例代码：**

```python
import random
import pysam
from dnakit.variants.p16_variant_cli_functions import compute_mosdepth_depth, compute_mosdepth_regions

rng = random.Random(11)
genome = "".join(rng.choice("ACGT") for _ in range(6000))
open("genome.fa", "w").write(">chrA\n" + genome + "\n")
pysam.faidx("genome.fa")
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrA", "LN": 6000}]}
with pysam.AlignmentFile("reads.bam", "wb", header=header) as out:
    for n, start in enumerate(range(0, 5850, 25)):
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence, read.flag = f"r{n}", genome[start:start + 150], 0
        read.reference_id, read.reference_start, read.mapping_quality, read.cigar = 0, start, 60, [(0, 150)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 150)
        out.write(read)
pysam.index("reads.bam")
depth = compute_mosdepth_depth("reads.bam", "depth", verbose=False)
print(depth["summary"][-1]["chrom"], depth["summary"][-1]["length"], depth["summary"][-1]["mean"])
regions = compute_mosdepth_regions("reads.bam", "win", by="2000", thresholds=[1, 5], verbose=False)
print([(r["start"], r["end"], r["depth"]) for r in regions["regions"]])
```

- **示例输出：**

```text
total 6000 5.85
[(0, 2000, 5.81), (2000, 4000, 6.0), (4000, 6000, 5.74)]
```

## 14) `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_quantized_coverage()`

- **作用：** 把逐位点深度合并成带标签的覆盖度区段（`mosdepth --quantize`）。
- **处理方法：** 按分箱规格（如 `0:1:4:`）生成 `quantized.bed.gz`。
- **API：** `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_quantized_coverage(alignment, prefix, *, quantize, chromosome=None, threads=0, mapq=0, fast_mode=False, reference=None, mosdepth_path=None, verbose=True)`。
- **输入：** BAM/CRAM、输出前缀与分箱规格。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `segments`。
- **示例见：** `compute_mosdepth_depth()`

## 15) `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_regions()`

- **作用：** 计算窗口或 BED 区域的平均（或中位）深度、深度直方图与阈值计数（`mosdepth --by`）。
- **处理方法：** `by` 为窗口大小或 BED；写出 `regions.bed.gz`、`region.dist.txt`，给出 `thresholds` 时还写 `thresholds.bed.gz`。
- **API：** `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_regions(alignment, prefix, *, by, thresholds=None, use_median=False, per_base=False, chromosome=None, threads=0, mapq=0, fast_mode=False, fragment_mode=False, exclude_flag=None, include_flag=None, reference=None, read_groups=None, mosdepth_path=None, verbose=True)`。
- **输入：** BAM/CRAM、输出前缀、窗口或 BED。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `summary`、`regions`。
- **示例见：** `compute_mosdepth_depth()`

## 16) `dnakit.variants.p16_variant_cli_functions.genotype_structural_variants_cutesv()`

- **作用：** 用 cuteSV 检测并基因分型长读长 SV（`cuteSV --genotype`）。
- **处理方法：** `genotype` 增加基因型与等位深度，`force_calling_vcf` 对给定 VCF 强制分型，`platform` 套用 ont/hifi/clr 参数。
- **API：** `dnakit.variants.p16_variant_cli_functions.genotype_structural_variants_cutesv(alignment, reference, vcf_path, work_directory, *, genotype=True, force_calling_vcf=None, platform=None, min_support=None, min_size=None, sample=None, threads=1, cutesv_path=None, verbose=True)`。
- **输入：** BAM、参考 FASTA、新 VCF 与新工作目录。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `genotype`、`records`。

## 17) `dnakit.variants.p16_variant_cli_functions.merge_structural_variants_sniffles2_snf()`

- **作用：** 由多个样本的 SNF 候选文件做多样本 SV 调用（`sniffles --input a.snf b.snf`）。
- **处理方法：** 合并各样本候选并联合基因分型，写出多样本 VCF。
- **API：** `dnakit.variants.p16_variant_cli_functions.merge_structural_variants_sniffles2_snf(snf_files, vcf_path, *, reference=None, min_support=None, threads=1, sniffles_path=None, verbose=True)`。
- **输入：** 至少两个 SNF 文件、新 VCF 路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `samples`、`records`。

## 18) `dnakit.annotation.p16_assembly_annotation_cli_functions.align_proteins_to_genome_miniprot()`

- **作用：** 把蛋白序列剪接感知且移码感知地比对到基因组（`miniprot`）。
- **处理方法：** 输出 PAF、GFF3 或 GTF；`translate` 追加翻译蛋白，`include_alignment` 追加残基比对。
- **API：** `dnakit.annotation.p16_assembly_annotation_cli_functions.align_proteins_to_genome_miniprot(genome_or_index, proteins, output_path, *, output_format='paf', translate=False, include_alignment=False, no_splicing=False, max_intron_size=None, splice_model=None, translation_table=1, threads=1, miniprot_path=None, verbose=True)`。
- **输入：** 基因组或索引、蛋白 FASTA、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `output_format`、`lines`、`head`。
- **示例见：** `build_miniprot_index()`

## 19) `dnakit.annotation.p16_assembly_annotation_cli_functions.assemble_long_reads_flye()`

- **作用：** 用 Flye 从头组装长读长，可选宏基因组模式与支架。
- **处理方法：** `metagenome` 为 metaFlye（覆盖度不均），`scaffold` 沿重复图连接 contig，`iterations` 为抛光轮数。
- **API：** `dnakit.annotation.p16_assembly_annotation_cli_functions.assemble_long_reads_flye(reads, output_directory, *, read_type='nano-raw', genome_size=None, metagenome=False, scaffold=False, iterations=1, min_overlap=None, threads=1, flye_path=None, verbose=True)`。
- **输入：** 长读长文件、新输出目录与读长类型。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `contigs`（来自 assembly_info.txt）。

## 20) `dnakit.annotation.p16_assembly_annotation_cli_functions.build_miniprot_index()`

- **作用：** 构建可复用的 miniprot 基因组索引（`miniprot -d`）。
- **处理方法：** 对基因组 FASTA 建立索引文件。
- **API：** `dnakit.annotation.p16_assembly_annotation_cli_functions.build_miniprot_index(genome, index_path, *, kmer=6, translation_table=1, threads=1, miniprot_path=None, verbose=True)`。
- **输入：** 基因组 FASTA、新索引路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `index`。
- **示例代码：**

```python
from Bio.Seq import Seq
from dnakit.annotation.p16_assembly_annotation_cli_functions import align_proteins_to_genome_miniprot, build_miniprot_index

genome = open("genome.fa").read().split("\n")[1]
open("protein.faa", "w").write(">prot1\n" + str(Seq(genome[3000:3600]).translate()).replace("*", "A") + "\n")
index = build_miniprot_index("genome.fa", "genome.mpi", verbose=False)
paf = align_proteins_to_genome_miniprot(index["index"], "protein.faa", "aln.paf", verbose=False)
fields = paf["head"][0].split("\t")
print(fields[0], fields[5], fields[7], fields[8])
```

- **示例输出：**

```text
prot1 chrA 3000 3600
```

## 21) `dnakit.annotation.p16_assembly_annotation_cli_functions.lift_over_annotation_liftoff()`

- **作用：** 把 GFF/GTF 基因注释从参考迁移到目标组装（`liftoff`）。
- **处理方法：** 用 minimap2 比对基因并逐外显子映射；`extra_copies` 搜索额外拷贝，`polish` 重比对以恢复起止密码子。
- **API：** `dnakit.annotation.p16_assembly_annotation_cli_functions.lift_over_annotation_liftoff(target_fasta, reference_fasta, reference_annotation, output_gff, *, unmapped_path=None, extra_copies=False, polish=False, min_coverage=0.5, min_identity=0.5, copy_identity=0.95, feature_types=None, threads=1, liftoff_path=None, verbose=True)`。
- **输入：** 目标与参考 FASTA、参考注释、新 GFF 路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `lifted_genes`、`unmapped_features`。

## 22) `dnakit.annotation.p16_assembly_annotation_cli_functions.polish_assembly_flye()`

- **作用：** 用长读长抛光已有草图组装（`flye --polish-target`）。
- **处理方法：** 读长可为 FASTA/FASTQ 或比对 BAM，输出 `polished_<n>.fasta`。
- **API：** `dnakit.annotation.p16_assembly_annotation_cli_functions.polish_assembly_flye(draft, reads, output_directory, *, read_type='nano-raw', iterations=1, threads=1, flye_path=None, verbose=True)`。
- **输入：** 草图 FASTA、长读长、新输出目录。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `polished`。

## 23) `dnakit.epigenome.p16_methyldackel_cli_functions.convert_mappability_bigwig_to_bbm_methyldackel()`

- **作用：** 把可映射性 bigWig 转成二进制 Bismap（`.bbm`）文件（`MethylDackel extract -M … -O`）。
- **处理方法：** MethylDackel 在按 bigWig 过滤的 extract 过程中写出 BBM；0.6.1 不接受 `-N <文件>`，故在临时目录用符号链接生成后移到目标路径。
- **API：** `dnakit.epigenome.p16_methyldackel_cli_functions.convert_mappability_bigwig_to_bbm_methyldackel(reference, alignment, bigwig, bbm_path, output_prefix, *, mappability_threshold=0.01, min_mappable_bases=15, region=None, methyldackel_path=None, verbose=True)`。
- **输入：** 参考、BAM、bigWig 与新的 BBM 路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `bbm`、`bbm_bytes`。

## 24) `dnakit.epigenome.p16_methyldackel_cli_functions.extract_methylation_report_methyldackel()`

- **作用：** 从亚硫酸氢盐 BAM 生成逐胞嘧啶甲基化报告或合并的 CpG/CHG 指标（`MethylDackel extract`）。
- **处理方法：** `report` 可为 `cytosine_report`、`merge_context`、`fraction` 或 `per_cytosine`；可选 CHG/CHH 与 OT/OB 纳入范围。
- **API：** `dnakit.epigenome.p16_methyldackel_cli_functions.extract_methylation_report_methyldackel(reference, alignment, output_prefix, *, report='cytosine_report', include_chg=False, include_chh=False, include_cpg=True, min_depth=1, min_mapq=10, min_base_quality=5, region=None, regions_bed=None, inclusion_bounds_ot=None, inclusion_bounds_ob=None, threads=1, methyldackel_path=None, verbose=True)`。
- **输入：** 参考 FASTA（已 faidx）、排序 BAM、输出前缀。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `report`、`files`、`head`。

## 25) `dnakit.epigenome.p16_methyldackel_cli_functions.merge_cytosine_context_methyldackel()`

- **作用：** 把逐胞嘧啶指标合并成逐 CpG/CHG 指标（`MethylDackel mergeContext`）。
- **处理方法：** 输入需按坐标排序，例如普通 extract 的 `_CpG.bedGraph`。
- **API：** `dnakit.epigenome.p16_methyldackel_cli_functions.merge_cytosine_context_methyldackel(reference, per_cytosine_metrics, output_path, *, methyldackel_path=None, verbose=True)`。
- **输入：** 参考 FASTA 与逐胞嘧啶指标文件。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `sites`。

## 26) `dnakit.epigenome.p16_methyldackel_cli_functions.methylation_bias_methyldackel()`

- **作用：** 按 read 位置统计甲基化偏倚，输出 SVG 与建议裁剪范围（`MethylDackel mbias`）。
- **处理方法：** 写出各链各上下文的 SVG，并解析程序给出的 `--OT/--OB/--CTOT/--CTOB` 建议值。
- **API：** `dnakit.epigenome.p16_methyldackel_cli_functions.methylation_bias_methyldackel(reference, alignment, output_prefix, *, include_chg=False, include_chh=False, min_mapq=10, min_base_quality=5, region=None, regions_bed=None, max_depth=2000, threads=1, methyldackel_path=None, verbose=True)`。
- **输入：** 参考 FASTA、排序 BAM、输出前缀。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `suggested_bounds`、SVG 文件。

## 27) `dnakit.epigenome.p16_methyldackel_cli_functions.summarize_per_read_methylation_methyldackel()`

- **作用：** 按 read 汇总 CpG 甲基化（`MethylDackel perRead`）。
- **处理方法：** 每条 read 输出染色体、位置、甲基化百分比和信息碱基数。
- **API：** `dnakit.epigenome.p16_methyldackel_cli_functions.summarize_per_read_methylation_methyldackel(reference, alignment, output_path, *, min_mapq=10, min_base_quality=5, region=None, regions_bed=None, threads=1, methyldackel_path=None, verbose=True)`。
- **输入：** 参考 FASTA、BAM、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `reads`、`rows`。

## 28) `dnakit.phylogeny.p16_fasttree_cli_functions.compute_fasttree_distance_matrix()`

- **作用：** 计算比对的成对距离矩阵，输出 PHYLIP 格式（`FastTree -makematrix`）。
- **处理方法：** 核酸用 Jukes-Cantor 距离，可选 `-rawdist` 与伪计数。
- **API：** `dnakit.phylogeny.p16_fasttree_cli_functions.compute_fasttree_distance_matrix(alignment, output_path, *, nucleotide=True, raw_distance=False, pseudocounts=None, fasttree_path=None, verbose=True)`。
- **输入：** 比对 FASTA、新输出路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `labels`、`matrix`。
- **示例代码：**

```python
from dnakit.phylogeny.p16_fasttree_cli_functions import compute_fasttree_distance_matrix, infer_fasttree_tree

matrix = compute_fasttree_distance_matrix("aln.fa", "dist.phy", verbose=False)
print(matrix["labels"], [round(x, 4) for x in matrix["matrix"][0]])
tree = infer_fasttree_tree("aln.fa", "tree.nwk", refine="noml", bootstrap_replicates=None, verbose=False)
print(tree["tips"], tree["newick"].count("("))
```

- **示例输出：**

```text
['t1', 't2', 't3', 't4', 't5'] [0.0, 0.0, 0.0458, 0.0946, 0.1743]
5 3
```

## 29) `dnakit.phylogeny.p16_fasttree_cli_functions.infer_fasttree_tree()`

- **作用：** 用 FastTree 的各种模式推断或精修进化树。
- **处理方法：** 支持起始树、`mllen` 分支长度优化、`noml` 最小进化树、拆分约束、bootstrap 支持值。
- **API：** `dnakit.phylogeny.p16_fasttree_cli_functions.infer_fasttree_tree(alignment, output_path, *, nucleotide=True, gtr=False, gamma=False, starting_tree=None, refine=None, constraints=None, constraint_weight=None, bootstrap_replicates=None, no_support=False, seed=None, fasttree_path=None, verbose=True)`。
- **输入：** 比对 FASTA、新 Newick 路径。
- **输出：** 字典，含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`（写出的文件）、`log_tail` 和 `provenance`；另有 `refine`、`newick`、`tips`。
- **示例见：** `compute_fasttree_distance_matrix()`
