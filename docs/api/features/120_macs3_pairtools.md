# MACS3 与 pairtools：峰识别与 Hi-C pair 处理

这一组函数是 `dnakit.epigenome` 中针对两个软件包的命令行封装：[MACS3](https://github.com/macs3-project/MACS)（峰识别、信号轨道、片段长度模型、峰内变异识别以及 HMMRATAC 可及性识别）和 [pairtools](https://github.com/open2c/pairtools)（Hi-C pair 的解析、排序、去重、过滤、分相与汇总）。每个函数用明确的参数列表（不经过 shell）运行一个固定的子命令，拒绝覆盖已存在的输出，并返回包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs` 和 `log_tail` 的字典。只有程序退出码为 0 且所有预期输出文件都存在时，`status` 才是 `completed`。它们不会下载数据，也不访问远程服务。

按软件包元数据，MACS3 3.0.5 采用 BSD-3-Clause 许可，pairtools 1.1.2 采用 MIT 许可。DNAKit 把它们作为独立程序启动，不复制其代码；只有 `pairtools_mark_split_pair_as_dup()` 会导入 pairtools 的模块。

依赖是可选的（extra 为 `dnakit[bio-macs3]` 和 `dnakit[bio-pairtools]`），只在调用函数时才使用。程序缺失时，函数会抛出带安装提示的 `BackendUnavailableError`。示例从 `dnakit.epigenome` 导入，所有函数都以 `verbose=False` 调用；它们会在 Python 解释器所在目录中找到程序，所以在不在 PATH 上的虚拟环境中也能运行。

## 示例数据

示例使用下面这些很小的文件。用于变异识别的 reads、片段和 BAM 文件由示例自己生成；函数从不覆盖已有文件，所以每个输出都写入新建的临时目录。

**示例文件：** `reads.bed`

```text
chr1	1000	1050	r0	0	+
chr1	1004	1054	r1	0	+
chr1	1008	1058	r2	0	+
chr1	1012	1062	r3	0	+
chr1	1016	1066	r4	0	+
chr1	1020	1070	r5	0	+
chr1	1130	1180	r6	0	-
chr1	1134	1184	r7	0	-
chr1	1138	1188	r8	0	-
chr1	1142	1192	r9	0	-
chr1	1146	1196	r10	0	-
chr1	1150	1200	r11	0	-
chr1	2000	2050	r12	0	+
chr1	2004	2054	r13	0	+
chr1	2008	2058	r14	0	+
chr1	2012	2062	r15	0	+
chr1	2016	2066	r16	0	+
chr1	2020	2070	r17	0	+
chr1	2130	2180	r18	0	-
chr1	2134	2184	r19	0	-
chr1	2138	2188	r20	0	-
chr1	2142	2192	r21	0	-
chr1	2146	2196	r22	0	-
chr1	2150	2200	r23	0	-
chr1	3000	3050	r24	0	+
chr1	3004	3054	r25	0	+
chr1	3008	3058	r26	0	+
chr1	3012	3062	r27	0	+
chr1	3016	3066	r28	0	+
chr1	3020	3070	r29	0	+
chr1	3130	3180	r30	0	-
chr1	3134	3184	r31	0	-
chr1	3138	3188	r32	0	-
chr1	3142	3192	r33	0	-
chr1	3146	3196	r34	0	-
chr1	3150	3200	r35	0	-
chr1	1000	1050	dup1	0	+
chr1	1004	1054	dup2	0	+
```

BED 格式的单端 reads：`chr1` 上三个对称的簇，每簇六条正链、六条负链 reads，并在第一个位置有两条重复 reads。

**示例文件：** `treat.bdg`

```text
chr1	0	100	1.00000
chr1	100	400	12.00000
chr1	400	600	2.00000
chr1	600	900	15.00000
chr1	900	1000	1.00000
```

条件 1 的 pileup bedGraph。

**示例文件：** `treat2.bdg`

```text
chr1	0	100	1.00000
chr1	100	400	3.00000
chr1	400	600	14.00000
chr1	600	900	16.00000
chr1	900	1000	1.00000
```

条件 2（第二个重复或条件）的 pileup bedGraph。

**示例文件：** `control.bdg`

```text
chr1	0	1000	1.50000
```

只有一个恒定区间的对照 lambda bedGraph。

**示例文件：** `cand.bed`

```text
chr1	960	1260	cand1
chr1	1960	2260	cand2
chr1	2960	3260	cand3
```

三个带名称的候选峰区域。

**示例文件：** `chroms.sizes`

```text
chr1	1000
chr2	800
```

染色体大小，制表符分隔。

**示例文件：** `toy.pairs`

```text
## pairs format v1.0.0
#shape: upper triangle
#sorted: chr1>chr2
#genome_assembly: toy
#chromsize: chr1 1000
#chromsize: chr2 800
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
a1	chr1	100	chr1	500	+	-	UU
a2	chr1	101	chr1	501	+	-	UU
a3	chr1	700	chr1	900	-	-	UU
a4	chr1	120	chr2	300	+	+	UU
a5	chr2	400	chr1	200	-	+	UU
```

已排序的上三角 .pairs 文件，共五个 pair（`a2` 与 `a1` 几乎重复，`a5` 位于下三角）。

**示例文件：** `other.pairs`

```text
## pairs format v1.0.0
#shape: upper triangle
#sorted: chr1>chr2
#genome_assembly: toy
#chromsize: chr1 1000
#chromsize: chr2 800
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
c1	chr1	300	chr1	800	+	-	UU
c2	chr2	100	chr2	200	+	-	UU
```

第二个已排序的 .pairs 文件，头部相同，有两个 pair。

**示例文件：** `phased.pairs`

```text
## pairs format v1.0.0
#shape: upper triangle
#chromsize: chr1-1 1000
#chromsize: chr1-2 1000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type XB1 XB2 AS1 AS2 XS1 XS2
p1	chr1-1	100	chr1-1	500	+	-	UU	.	.	50	50	20	20
p2	chr1-1	200	chr1-2	600	+	-	UU	.	.	50	50	20	20
p3	chr1-1	300	chr1-1	700	+	-	UU	chr1-2,300,50M,0,50,60;	.	50	50	50	20
```

二倍体基因组（染色体 `chr1-1` 与 `chr1-2`）上的 .pairs 文件，带有 `pairtools phase` 需要的标签列。

**示例文件：** `toy.sam`

```text
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:800
q1	65	chr1	100	60	50M	chr2	300	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
q1	145	chr2	300	60	50M	chr1	100	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
q2	65	chr1	150	60	50M	chr1	600	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
q2	145	chr1	600	60	50M	chr1	150	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
q3	113	chr2	400	60	50M	chr1	200	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
q3	177	chr1	200	60	50M	chr2	400	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

三对按读名排序的 reads（SAM，含头部和 50 bp 比对）。

## 1) `dnakit.epigenome.call_atac_peaks_hmmratac()`

- **作用：** 用 HMMRATAC 隐马尔可夫模型，从双端片段中识别 ATAC-seq 的可及区域。
- **处理方法：** 检查参数后运行 `macs3 hmmratac`，传入 `-i`、`-f`、`-n`、`--outdir`、`-l`、`-u`、`-c`、`--minlen`、`--maxTrain`、`--randomSeed` 和 `--hmm-type`（`training_regions` 对应 `-t`，`save_states` 对应 `--save-states`）。`output_dir` 由函数创建，调用前不能存在。只有 `<name>_accessible_regions.narrowPeak`、`<name>_cutoff_analysis.tsv` 和 `<name>_model.json` 都存在（`save_states` 时还有 `<name>_states.bed`）时状态才是 `completed`。MACS3 找不到训练区域时会以 "Not enough training regions" 停止，此时状态为 `failed`。
- **API：** `call_atac_peaks_hmmratac(ifile: 'Sequence[str]', output_dir: 'str', *, name: 'str' = 'NA', file_format: 'str' = 'BAMPE', lower_cutoff: 'int' = 10, upper_cutoff: 'int' = 20, prescan_cutoff: 'float' = 1.2, min_region: 'int' = 100, max_training: 'int' = 1000, random_seed: 'int' = 10151, hmm_type: 'str' = 'gaussian', training_regions: 'str | None' = None, save_states: 'bool' = False, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `ifile` 是已存在的片段文件列表；`file_format` 为 `BAMPE`（默认）、`BEDPE` 或 `FRAG`。`lower_cutoff` 必须小于 `upper_cutoff`（选择训练区域的倍数变化界限，默认 10 和 20）。`prescan_cutoff`（默认 1.2）、`min_region`（默认 100）、`max_training`（默认 1000）和 `random_seed`（默认 10151）会原样传给 MACS3；`hmm_type` 为 `gaussian`（默认）或 `poisson`。`training_regions` 是可选的已存在 BED 文件。`name` 只能包含字母、数字、`_`、`.` 和 `-`。
- **输出：** 通用结果字典，另外还有 `output_dir` 和 `accessible_regions`（narrowPeak 文件路径；未生成时为 `None`）。
- **示例代码：**

```python
import random
import tempfile
from pathlib import Path

from dnakit.epigenome import call_atac_peaks_hmmratac

work = Path(tempfile.mkdtemp())
rng = random.Random(0)
rows = []
for centre in (2000, 6000, 10000):
    for _ in range(40):
        start = centre + rng.randint(-100, 100)
        rows.append(f"chr1\t{start}\t{start + rng.randint(50, 150)}")
for _ in range(40):
    start = rng.randint(0, 15000)
    rows.append(f"chr1\t{start}\t{start + rng.randint(100, 300)}")
fragments = work / "fragments.bedpe"
fragments.write_text("\n".join(rows) + "\n")

result = call_atac_peaks_hmmratac(
    [str(fragments)], str(work / "out"), name="demo",
    file_format="BEDPE", save_states=True, verbose=False,
)
print(result["status"], result["return_code"])
print(sorted(Path(p).name for p in result["outputs"]))
states = (work / "out" / "demo_states.bed").read_text().splitlines()
print(len(states), sorted({line.split("\t")[3] for line in states}))
print(result["accessible_regions"] is not None)
```

- **示例输出：**

```text
completed 0
['demo_accessible_regions.narrowPeak', 'demo_cutoff_analysis.tsv', 'demo_model.json', 'demo_states.bed']
19 ['nuc', 'open']
True
```

## 2) `dnakit.epigenome.call_broad_peaks_from_bedgraph()`

- **作用：** 从得分 bedGraph 中识别嵌套的宽峰，并写成 gappedPeak 表。
- **处理方法：** 运行 `macs3 bdgbroadcall`，传入 `-i`、`-c`（`cutoff_peak`）、`-C`（`cutoff_link`）、`-l`、`-g`（`lvl1_max_gap`）、`-G`（`lvl2_max_gap`）和 `-o`；`trackline` 为假时追加 `--no-trackline`。得分不低于 `cutoff_peak` 的区域构成宽区域，其中高于 `cutoff_link` 且相近的核心区被连接为区块。
- **API：** `call_broad_peaks_from_bedgraph(bedgraph: 'str', output_file: 'str', *, cutoff_peak: 'float' = 2.0, cutoff_link: 'float' = 1.0, min_length: 'int' = 200, lvl1_max_gap: 'int' = 30, lvl2_max_gap: 'int' = 800, trackline: 'bool' = True, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `bedgraph` 是已存在的得分 bedGraph，`output_file` 是已存在文件夹中的新路径。`cutoff_peak`（默认 2.0）、`cutoff_link`（默认 1.0）、`min_length`（默认 200，至少为 1）、`lvl1_max_gap`（默认 30）和 `lvl2_max_gap`（默认 800）为数值；`trackline`（默认真）写入 UCSC track 行。
- **输出：** 通用结果字典。输出文件有 15 列（BED12 加三个得分列）；区域名称中含有输出路径。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import call_broad_peaks_from_bedgraph

out = Path(tempfile.mkdtemp()) / "broad.bed"
result = call_broad_peaks_from_bedgraph(
    "treat.bdg", str(out), cutoff_peak=10, cutoff_link=1, min_length=100,
    lvl1_max_gap=50, lvl2_max_gap=300, trackline=False, verbose=False,
)
print(result["status"], result["return_code"])
for line in out.read_text().splitlines():
    fields = line.split("\t")
    print(len(fields), fields[0], fields[1], fields[2], fields[4], fields[9], fields[10], fields[11])
```

- **示例输出：**

```text
completed 0
15 chr1 0 1000 150 4 1,300,300,1 0,100,600,999
```

## 3) `dnakit.epigenome.call_chip_peaks()`

- **作用：** 从 ChIP-seq 或 ATAC-seq 比对结果中识别峰。
- **处理方法：** 检查全部选项后，在新的 `output_dir` 中运行 `macs3 callpeak`，传入 `-t`、可选的 `-c`、`-f`、`-g`、`--keep-dup`、`-n`、`--outdir`，模型选项（`--shift`、`--extsize`、`--bw`、`--d-min`、`-m`、`--nomodel`、`--fix-bimodal`）和阈值选项（`-q` 或 `-p`、`--slocal`、`--llocal`、`--fe-cutoff`、`--broad`、`--call-summits`）。只有退出码为 0 且预期文件都存在时状态才是 `completed`：`<name>_peaks.narrowPeak`（`broad` 时为 `.broadPeak` 和 `.gappedPeak`）、`<name>_peaks.xls`，以及随选项产生的文件（峰顶、bedGraph 轨道、阈值分析）。未设置 `nomodel` 时，MACS3 会建立片段长度模型，需要足够多成对的正负链峰；建不起来时子命令失败，状态为 `failed`。
- **API：** `call_chip_peaks(treatment: 'Sequence[str]', output_dir: 'str', *, control: 'Sequence[str] | None' = None, name: 'str' = 'NA', file_format: 'str' = 'AUTO', genome_size: 'str | float' = 'hs', tag_size: 'int | None' = None, keep_duplicates: 'str | int' = '1', max_count: 'int | None' = None, barcodes: 'str | None' = None, save_bedgraph: 'bool' = False, spmr: 'bool' = False, trackline: 'bool' = False, nomodel: 'bool' = False, shift: 'int' = 0, extsize: 'int' = 200, bandwidth: 'int' = 300, min_fragment: 'int' = 20, mfold: 'Sequence[int]' = (5, 50), fix_bimodal: 'bool' = False, qvalue: 'float | None' = 0.05, pvalue: 'float | None' = None, scale_to: 'str | None' = None, down_sample: 'bool' = False, seed: 'int | None' = None, nolambda: 'bool' = False, small_local: 'int' = 1000, large_local: 'int' = 10000, max_gap: 'int | None' = None, min_length: 'int | None' = None, broad: 'bool' = False, broad_cutoff: 'float' = 0.1, cutoff_analysis: 'bool' = False, call_summits: 'bool' = False, fe_cutoff: 'float' = 1.0, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `treatment` 是已存在的比对文件列表，`control` 为可选列表；`file_format` 取 `AUTO`（默认）、`BAM`、`SAM`、`BED`、`ELAND`、`ELANDMULTI`、`ELANDEXPORT`、`BOWTIE`、`BAMPE`、`BEDPE` 或 `FRAG`。`genome_size` 为 `hs`、`mm`、`ce`、`dm` 或正的碱基数。`qvalue`（默认 0.05）和 `pvalue` 互斥：要用 `pvalue` 就传 `qvalue=None`。`mfold` 是 `(下界, 上界)`，下界必须小于上界。`scale_to` 为 `large` 或 `small`。`output_dir` 不能已存在；`name` 只能包含字母、数字、`_`、`.` 和 `-`。其余选项（`tag_size`、`keep_duplicates`、`max_count`、`barcodes`、`save_bedgraph`、`spmr`、`trackline`、`nomodel`、`shift`、`extsize`、`bandwidth`、`min_fragment`、`fix_bimodal`、`down_sample`、`seed`、`nolambda`、`small_local`、`large_local`、`max_gap`、`min_length`、`broad`、`broad_cutoff`、`cutoff_analysis`、`call_summits`、`fe_cutoff`）对应 `callpeak` 的同名参数。
- **输出：** 通用结果字典，另外还有 `output_dir`、`peak_files`（实际存在的峰文件）、`xls_file`（路径或 `None`）和 `model_script`（`<name>_model.r` 路径，或 `None`；`nomodel` 时恒为 `None`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import call_chip_peaks

out = Path(tempfile.mkdtemp()) / "peaks"
result = call_chip_peaks(
    ["reads.bed"], str(out), name="demo", file_format="BED", genome_size=5000,
    nomodel=True, extsize=100, qvalue=0.5, verbose=False,
)
print(result["status"], result["return_code"])
print(sorted(Path(p).name for p in result["outputs"]))
print([Path(p).name for p in result["peak_files"]], result["model_script"])
for line in Path(result["peak_files"][0]).read_text().splitlines():
    print(line)
```

- **示例输出：**

```text
completed 0
['demo_peaks.narrowPeak', 'demo_peaks.xls']
['demo_peaks.narrowPeak'] None
chr1	1004	1196	demo_peak_1	71	.	6.97674	10.6804	7.17525	96
chr1	2004	2196	demo_peak_2	71	.	6.97674	10.6804	7.17525	96
chr1	3004	3196	demo_peak_3	71	.	6.97674	10.6804	7.17525	96
```

## 4) `dnakit.epigenome.call_peaks_from_bedgraph()`

- **作用：** 从 bdgcmp 等得分 bedGraph 中识别窄峰。
- **处理方法：** 运行 `macs3 bdgpeakcall`，传入 `-i`、`-c`、`-l`、`-g` 和 `-o`（按需追加 `--call-summits` 与 `--no-trackline`）。得分不低于 `cutoff` 的区域，在间隙不超过 `max_gap` 时合并，长度不小于 `min_length` 时保留。
- **API：** `call_peaks_from_bedgraph(bedgraph: 'str', output_file: 'str', *, cutoff: 'float' = 5.0, min_length: 'int' = 200, max_gap: 'int' = 30, call_summits: 'bool' = False, trackline: 'bool' = True, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `bedgraph` 是已存在的得分 bedGraph，`output_file` 是已存在文件夹中的新路径。`cutoff`（默认 5.0）、`min_length`（默认 200，至少为 1）和 `max_gap`（默认 30，至少为 0）为数值；`call_summits`（默认假）按峰顶拆分峰；`trackline`（默认真）写入 UCSC track 行。
- **输出：** 通用结果字典。输出为 narrowPeak 风格的表；峰名称中含有输出路径。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import call_peaks_from_bedgraph

out = Path(tempfile.mkdtemp()) / "peaks.bed"
result = call_peaks_from_bedgraph(
    "treat.bdg", str(out), cutoff=5, min_length=100, max_gap=50,
    trackline=False, verbose=False,
)
print(result["status"], result["return_code"])
for line in out.read_text().splitlines():
    fields = line.split("\t")
    print(fields[0], fields[1], fields[2], fields[4])
```

- **示例输出：**

```text
completed 0
chr1 100 400 120
chr1 600 900 150
```

## 5) `dnakit.epigenome.call_variants_in_peaks()`

- **作用：** 在峰区域内，从 BAM 文件中识别 SNV 和小的插入缺失，并写成 VCF。
- **处理方法：** 运行 `macs3 callvar`，传入 `-b`、`-t`、可选的 `-c`、`-o`、`-Q`、`-D`、`-g`、`-G`、`-F`、`--fermi-overlap`、`--top2alleles-mratio`、`--altallele-count`、`--max-ar` 和 `-m`。运行后统计日志尾部中 "No reads found in this peak" 的行数并写入 `warnings`，因为 MACS3 会把这类峰标为跳过，但退出码仍为 0。
- **API：** `call_variants_in_peaks(peaks_bed: 'str', treatment_bam: 'str', output_vcf: 'str', *, control_bam: 'str | None' = None, min_base_quality: 'int' = 20, max_duplicates: 'int' = 1, gq_hetero: 'int' = 0, gq_homo: 'int' = 0, fermi: 'str' = 'auto', fermi_min_overlap: 'int' = 30, top2alleles_min_ratio: 'float' = 0.8, altallele_min_count: 'int' = 2, max_allowed_ar: 'float' = 0.95, processes: 'int' = 1, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。示例还需要 pysam（extra 为 `dnakit[bio-reads]`）来写一个小的 BAM 文件。
- **输入：** `peaks_bed` 是至少四列的已存在 BED 文件；`treatment_bam`（以及可选的 `control_bam`）必须已排序、已建索引并带有 MD 标签（例如由 `samtools calmd` 生成）；`output_vcf` 是新路径。`min_base_quality`（默认 20）、`max_duplicates`（默认 1）、`gq_hetero`、`gq_homo`（默认 0）、`fermi`（`auto`、`on` 或 `off`）、`fermi_min_overlap`（默认 30）、`top2alleles_min_ratio` 与 `max_allowed_ar`（取值 0 到 1）、`altallele_min_count`（默认 2）以及 `processes`（默认 1）会被检查并传给 MACS3。
- **输出：** 通用结果字典，另外还有 `warnings`（字符串列表；没有峰被跳过时为空）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.epigenome import call_variants_in_peaks

work = Path(tempfile.mkdtemp())
reference = "ACGT" * 250
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
bam = str(work / "reads.bam")
with pysam.AlignmentFile(bam, "wb", header=header) as handle:
    for i in range(10):
        start = 90 + i
        bases = list(reference[start:start + 40])
        md = "40"
        if i < 6:
            offset = 120 - start
            bases[offset] = "T"
            md = f"{offset}{reference[120]}{40 - offset - 1}"
        read = pysam.AlignedSegment()
        read.query_name = f"read{i}"
        read.query_sequence = "".join(bases)
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigar = [(0, 40)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 40)
        read.set_tag("MD", md)
        handle.write(read)
pysam.index(bam)
(work / "peaks.bed").write_text("chr1\t80\t160\tpeak1\n")

out = work / "variants.vcf"
result = call_variants_in_peaks(str(work / "peaks.bed"), bam, str(out), verbose=False)
print(result["status"], result["warnings"])
for line in out.read_text().splitlines():
    if not line.startswith("#"):
        fields = line.split("\t")
        info = dict(item.split("=") for item in fields[7].split(";"))
        print(fields[0], fields[1], fields[3], fields[4], info["MT"], info["DPT"], fields[9])
```

- **示例输出：**

```text
completed []
chr1 121 A T SNV 10 0/1:10:124:204,0,124
```

## 6) `dnakit.epigenome.combine_replicate_bedgraphs()`

- **作用：** 把两个或更多重复的得分轨道合并成一个 bedGraph。
- **处理方法：** 运行 `macs3 cmbreps`，传入 `-i`、`-m` 和 `-o`。`fisher` 用 Fisher 方法合并 -log10 p 值，`max` 与 `mean` 直接取得分的最大值或平均值。
- **API：** `combine_replicate_bedgraphs(inputs: 'Sequence[str]', output_file: 'str', *, method: 'str' = 'fisher', macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `inputs` 是至少两个已存在 bedGraph 文件的列表（只传一个字符串会被拒绝）；`method` 为 `fisher`（默认）、`max` 或 `mean`；`output_file` 是已存在文件夹中的新路径。使用 `fisher` 时输入应是 -log10 p 值轨道，例如 `compare_bedgraph_tracks()` 的 `ppois` 轨道。
- **输出：** 通用结果字典。输出的 bedGraph 以 track 行开头，每个区间一个得分。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import combine_replicate_bedgraphs, compare_bedgraph_tracks

work = Path(tempfile.mkdtemp())
tracks = []
for index, treatment in enumerate(["treat.bdg", "treat2.bdg"]):
    scored = compare_bedgraph_tracks(
        treatment, "control.bdg", str(work / f"rep{index}"), name="rep",
        methods=["ppois"], verbose=False,
    )
    tracks.append(scored["score_files"][0])

out = work / "combined.bdg"
result = combine_replicate_bedgraphs(tracks, str(out), method="fisher", verbose=False)
print(result["status"], result["return_code"])
for line in out.read_text().splitlines()[1:]:
    print(line)
```

- **示例输出：**

```text
completed 0
chr1	0	100	0.28851
chr1	100	400	7.94044
chr1	400	600	9.38984
chr1	600	900	21.54873
chr1	900	1000	0.28851
```

## 7) `dnakit.epigenome.compare_bedgraph_tracks()`

- **作用：** 把处理组 bedGraph 与对照 bedGraph 比较打分，每种方法生成一个得分轨道。
- **处理方法：** 在新的 `output_dir` 中运行 `macs3 bdgcmp`，传入 `-t`、`-c`、`-S`、`-p`、`-m`、`--outdir` 和 `--o-prefix`。每种方法写出 `<name>_<method>.bdg`；`ppois` 是以 -log10 表示的 Poisson p 值，`qpois` 是其 q 值，`FE` 是倍数富集，伪计数在按深度缩放之后加入。
- **API：** `compare_bedgraph_tracks(treatment: 'str', control: 'str', output_dir: 'str', *, name: 'str' = 'NA', methods: 'Sequence[str]' = ('ppois',), scaling_factor: 'float' = 1.0, pseudocount: 'float' = 0.0, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `treatment` 是 pileup bedGraph，`control` 是对照 lambda bedGraph（均可由带 `save_bedgraph` 的 `call_chip_peaks()` 生成）。`methods` 是非空列表，取自 `ppois`、`qpois`、`subtract`、`logFE`、`FE`、`logLR`、`slogLR`、`max`（默认 `('ppois',)`）；`scaling_factor`（默认 1.0）和 `pseudocount`（默认 0.0）为非负数；`output_dir` 不能已存在。
- **输出：** 通用结果字典，另外还有 `output_dir` 和 `score_files`（实际存在的得分轨道，顺序与 `methods` 一致）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import compare_bedgraph_tracks

out = Path(tempfile.mkdtemp()) / "scores"
result = compare_bedgraph_tracks(
    "treat.bdg", "control.bdg", str(out), name="demo",
    methods=["ppois", "FE"], verbose=False,
)
print(result["status"], [Path(p).name for p in result["score_files"]])
for path in result["score_files"]:
    print([line.split("\t")[3] for line in Path(path).read_text().splitlines()])
```

- **示例输出：**

```text
completed ['demo_ppois.bdg', 'demo_FE.bdg']
['0.35441', '8.10772', '0.71862', '11.11471', '0.35441']
['0.66667', '8.00000', '1.33333', '10.00000', '0.66667']
```

## 8) `dnakit.epigenome.differential_bedgraph_peaks()`

- **作用：** 根据两个条件各自的 pileup 与对照 bedGraph，识别两个条件之间有差异的区域。
- **处理方法：** 运行 `macs3 bdgdiff`，传入 `--t1`、`--t2`、`--c1`、`--c2`、`-C`、`-l`、`-g`、`--d1`、`--d2`，以及跟在 `-o` 后的三个输出路径。MACS3 为三类区域各写一个 BED 文件：条件 1 特有、条件 2 特有、两者共有。
- **API：** `differential_bedgraph_peaks(treatment1: 'str', treatment2: 'str', control1: 'str', control2: 'str', output_condition1: 'str', output_condition2: 'str', output_common: 'str', *, cutoff: 'float' = 3.0, min_length: 'int' = 200, max_gap: 'int' = 100, depth1: 'float' = 1.0, depth2: 'float' = 1.0, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `treatment1`、`treatment2`（pileup bedGraph）和 `control1`、`control2`（对照 lambda bedGraph）都是已存在的文件。`output_condition1`、`output_condition2` 和 `output_common` 是已存在文件夹中的三个新路径。`cutoff`（默认 3.0）是 log10 似然比阈值，`min_length`（默认 200）和 `max_gap`（默认 100）至少为 1，`depth1`、`depth2`（默认 1.0）是以百万计的非冗余测序深度。
- **输出：** 通用结果字典；`outputs` 列出三个文件。每个文件以 track 行开头，之后每个区域一行。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import differential_bedgraph_peaks

work = Path(tempfile.mkdtemp())
result = differential_bedgraph_peaks(
    "treat.bdg", "treat2.bdg", "control.bdg", "control.bdg",
    str(work / "cond1.bed"), str(work / "cond2.bed"), str(work / "common.bed"),
    cutoff=3.0, min_length=100, max_gap=50, verbose=False,
)
print(result["status"], result["return_code"])
for name in ("cond1.bed", "cond2.bed", "common.bed"):
    lines = (work / name).read_text().splitlines()
    print(name, lines[0].startswith("track"), len(lines) - 1)
```

- **示例输出：**

```text
completed 0
cond1.bed True 1
cond2.bed True 1
common.bed True 1
```

## 9) `dnakit.epigenome.filter_duplicate_reads()`

- **作用：** 去除重复 reads，并把剩余的比对结果写成 BED。
- **处理方法：** 运行 `macs3 filterdup`，传入 `-i`、`-f`、`-g`、`--keep-dup` 和 `-o`（`tag_size` 对应 `-s`）。双端输入写成 BEDPE。`auto` 保留二项检验在该基因组大小下允许的同一位置、同一链上的 tag 数，`all` 保留全部 tag，整数则最多保留该数量。
- **API：** `filter_duplicate_reads(ifile: 'Sequence[str]', output_file: 'str', *, file_format: 'str' = 'AUTO', genome_size: 'str | float' = 'hs', tag_size: 'int | None' = None, keep_duplicates: 'str | int' = 'auto', macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `ifile` 是已存在的比对文件列表；`file_format` 取 `AUTO`（默认）、`BAM`、`SAM`、`BED`、`ELAND`、`ELANDMULTI`、`ELANDEXPORT`、`BOWTIE`、`BAMPE` 或 `BEDPE`；`genome_size` 为 `hs`、`mm`、`ce`、`dm` 或正数；`keep_duplicates` 为 `auto`（默认）、`all` 或正整数；`output_file` 是新路径。
- **输出：** 通用结果字典；输出的 BED 保存保留下来的 tag。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import filter_duplicate_reads

work = Path(tempfile.mkdtemp())
for keep in (1, "all"):
    out = work / f"keep_{keep}.bed"
    result = filter_duplicate_reads(
        ["reads.bed"], str(out), file_format="BED", genome_size=5000,
        keep_duplicates=keep, verbose=False,
    )
    print(keep, result["status"], len(out.read_text().splitlines()))
```

- **示例输出：**

```text
1 completed 36
all completed 38
```

## 10) `dnakit.epigenome.operate_bedgraph_scores()`

- **作用：** 对 MACS3 bedGraph 的得分列做变换。
- **处理方法：** 运行 `macs3 bdgopt`，传入 `-i`、`-m`、`-p`（给出参数时）和 `-o`。`multiply`、`add`、`max` 和 `min` 需要恰好一个额外参数（减法用 `add` 加负数）。`p2q` 用 Benjamini-Hochberg 方法把 -log10 p 值得分转换成 q 值，不带参数，并且如 MACS3 文档所述需要覆盖全基因组的 bedGraph。
- **API：** `operate_bedgraph_scores(bedgraph: 'str', output_file: 'str', *, method: 'str' = 'p2q', extra_params: 'Sequence[float] | None' = None, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `bedgraph` 是已存在的文件，`output_file` 是新路径。`method` 为 `multiply`、`add`、`p2q`（默认）、`max` 或 `min`。前四种方法的 `extra_params` 是只有一个元素的列表，`p2q` 必须为空；数量不对会抛出 `ConfigurationError`。
- **输出：** 通用结果字典；输出的 bedGraph 以 track 行开头。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import operate_bedgraph_scores

work = Path(tempfile.mkdtemp())
for method, value in (("multiply", 2), ("add", -1)):
    out = work / f"{method}.bdg"
    result = operate_bedgraph_scores(
        "treat.bdg", str(out), method=method, extra_params=[value], verbose=False,
    )
    scores = [line.split("\t")[3] for line in out.read_text().splitlines()[1:]]
    print(method, result["status"], scores)
```

- **示例输出：**

```text
multiply completed ['2.00000', '24.00000', '4.00000', '30.00000', '2.00000']
add completed ['0.00000', '11.00000', '1.00000', '14.00000', '0.00000']
```

## 11) `dnakit.epigenome.pileup_reads_to_bedgraph()`

- **作用：** 把比对结果堆叠成覆盖度 bedGraph。
- **处理方法：** 运行 `macs3 pileup`，传入 `-i`、`-o`、`-f`、`--extsize`（另有 `--barcodes`、`--max-count` 和 `-B`）。单端 reads 沿下游延伸 `extsize` bp，设置 `both_direction` 则向两侧延伸。双端格式 `BAMPE`、`BEDPE` 以及 `FRAG` 直接按片段堆叠；`barcodes` 与 `max_count` 只适用于 `FRAG`。
- **API：** `pileup_reads_to_bedgraph(ifile: 'Sequence[str]', output_file: 'str', *, file_format: 'str' = 'AUTO', barcodes: 'str | None' = None, max_count: 'int | None' = None, both_direction: 'bool' = False, extsize: 'int' = 200, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `ifile` 是已存在的比对文件列表；`file_format` 为 `AUTO`（默认）、单端格式、`BAMPE`、`BEDPE` 或 `FRAG`；`extsize`（默认 200）至少为 1；`barcodes` 是可选的已存在文件；`max_count` 是可选的、不小于 0 的整数；`output_file` 是新路径。
- **输出：** 通用结果字典；输出的 bedGraph 每行包含 `chrom`、`start`、`end` 和覆盖度。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pileup_reads_to_bedgraph

out = Path(tempfile.mkdtemp()) / "pileup.bdg"
result = pileup_reads_to_bedgraph(
    ["reads.bed"], str(out), file_format="BED", extsize=100, verbose=False,
)
rows = [line.split("\t") for line in out.read_text().splitlines()]
print(result["status"], len(rows), max(float(row[3]) for row in rows))
print(rows[0], rows[1])
```

- **示例输出：**

```text
completed 67 13.0
['chr1', '0', '1000', '0.00000'] ['chr1', '1000', '1004', '2.00000']
```

## 12) `dnakit.epigenome.plot_fragment_model()`

- **作用：** 建立片段长度模型，返回模型数据，并在有 R 时渲染图像。
- **处理方法：** 运行 `macs3 predictd`（与 `predict_fragment_size()` 相同），它把模型写成 R 绘图脚本；随后解析该脚本中的向量得到 `data`：`distance_from_middle`、`forward_percent`、`reverse_percent`、`xcorr_lag`、`xcorr_correlation` 和 `alternative_fragment_lengths`。如果 PATH 上有 `Rscript`，就对该脚本运行它以绘制 `<script_name>_model.pdf`；否则 `image` 为 `None`，`image_note` 说明原因。获取数据不需要 R。
- **API：** `plot_fragment_model(ifile: 'Sequence[str]', output_dir: 'str', *, file_format: 'str' = 'AUTO', genome_size: 'str | float' = 'hs', tag_size: 'int | None' = None, bandwidth: 'int' = 300, min_fragment: 'int' = 20, mfold: 'Sequence[int]' = (5, 50), script_name: 'str' = 'predictd_model.R', macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。绘制图像还需要 PATH 上有外部程序 `Rscript`；没有它时仍会返回数据，`image` 为 `None`。
- **输入：** `ifile` 是已存在的单端比对文件列表（双端输入不写脚本，因此 `data` 保持 `None`）；`output_dir` 不能已存在；`file_format`、`genome_size`、`tag_size`、`bandwidth`（默认 300）、`min_fragment`（默认 20）和 `mfold`（默认 `(5, 50)`）与 `call_chip_peaks()` 相同；`script_name`（默认 `predictd_model.R`）只能包含字母、数字、`_`、`.` 和 `-`。
- **输出：** 通用结果字典，另外还有 `output_dir`、`fragment_length` 与 `fragment_length_source`、`model_script`、`script`、`data`（列表字典，或 `None`）、`image`（PDF 路径或 `None`）和 `image_note`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_fragment_model

work = Path(tempfile.mkdtemp())
rows = []
for k in range(150):
    centre = 1000 + k * 5000
    for i in range(12):
        rows.append(f"chr1\t{centre + i * 4}\t{centre + i * 4 + 50}\tr\t.\t+")
        end = centre + 180 + i * 4
        rows.append(f"chr1\t{end - 50}\t{end}\tr\t.\t-")
reads = work / "sim.bed"
reads.write_text("\n".join(rows) + "\n")

result = plot_fragment_model(
    [str(reads)], str(work / "out"), file_format="BED",
    genome_size=1000000, mfold=(2, 100), verbose=False,
)
print(result["status"], result["fragment_length"], result["fragment_length_source"])
print({key: len(values) for key, values in result["data"].items()})
print(result["data"]["alternative_fragment_lengths"])
print(result["image"] is None or result["image"].endswith("predictd_model.R_model.pdf"))
```

- **示例输出：**

```text
completed 179 single_end_model
{'distance_from_middle': 1211, 'forward_percent': 1211, 'reverse_percent': 1211, 'xcorr_lag': 1200, 'xcorr_correlation': 1200, 'alternative_fragment_lengths': 1}
[179]
True
```

## 13) `dnakit.epigenome.predict_fragment_size()`

- **作用：** 根据比对结果预测片段（d）大小，单端数据还会写出模型 R 脚本。
- **处理方法：** 在新的 `output_dir` 中运行 `macs3 predictd`，传入 `-i`、`-f`、`-g`、`--bw`、`--d-min`、`-m`、`--outdir` 和 `--rfile`，再从日志中读取大小。单端输入会建立峰模型，把 "predicted fragment length" 作为 `fragment_length`（来源 `single_end_model`）并写出 R 脚本，这需要至少 100 对成对的正负链峰；否则状态为 `failed`，`log_tail` 中是 MACS3 的警告。双端输入（`BAMPE`、`BEDPE`）则报告平均插入长度（来源 `paired_end_average`），不写脚本。
- **API：** `predict_fragment_size(ifile: 'Sequence[str]', output_dir: 'str', *, file_format: 'str' = 'AUTO', genome_size: 'str | float' = 'hs', tag_size: 'int | None' = None, bandwidth: 'int' = 300, min_fragment: 'int' = 20, mfold: 'Sequence[int]' = (5, 50), script_name: 'str' = 'predictd_model.R', macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `ifile` 是已存在的比对文件列表；`output_dir` 不能已存在；`file_format` 为 `AUTO`（默认）、单端格式、`BAMPE` 或 `BEDPE`；`genome_size` 为 `hs`、`mm`、`ce`、`dm` 或正数；`tag_size` 可选；`bandwidth`（默认 300）、`min_fragment`（默认 20）、`mfold`（默认 `(5, 50)`）和 `script_name`（默认 `predictd_model.R`）与 `plot_fragment_model()` 相同。
- **输出：** 通用结果字典，另外还有 `output_dir`、`fragment_length`（整数；未得到大小时为 `None`）、`fragment_length_source`（`single_end_model`、`paired_end_average` 或 `None`）和 `model_script`（路径或 `None`）。
- **示例代码：**

```python
import random
import tempfile
from pathlib import Path

from dnakit.epigenome import predict_fragment_size

work = Path(tempfile.mkdtemp())
rng = random.Random(0)
rows = []
for centre in (2000, 6000, 10000):
    for _ in range(40):
        start = centre + rng.randint(-100, 100)
        rows.append(f"chr1\t{start}\t{start + rng.randint(50, 150)}")
for _ in range(40):
    start = rng.randint(0, 15000)
    rows.append(f"chr1\t{start}\t{start + rng.randint(100, 300)}")
fragments = work / "fragments.bedpe"
fragments.write_text("\n".join(rows) + "\n")

result = predict_fragment_size(
    [str(fragments)], str(work / "out"), file_format="BEDPE", verbose=False,
)
print(result["status"], result["fragment_length"], result["fragment_length_source"])
print(result["model_script"], result["outputs"])
```

- **示例输出：**

```text
completed 127 paired_end_average
None []
```

## 14) `dnakit.epigenome.random_sample_reads()`

- **作用：** 随机保留一定数量或一定比例的比对结果，并写成 BED。
- **处理方法：** 运行 `macs3 randsample`，传入 `-i`、`-o`、`-f` 以及 `-n` 或 `-p`（另有 `--seed` 和 `-s`）。`number` 大约保留该数量的 tag，`percentage` 保留该百分比（80.0 表示 80%）。种子会传给 MACS3，但测试该函数时，同一个种子重复运行并不总是得到相同的样本，因此不能保证两次运行结果一致。
- **API：** `random_sample_reads(ifile: 'Sequence[str]', output_file: 'str', *, number: 'int | None' = None, percentage: 'float | None' = None, seed: 'int | None' = None, file_format: 'str' = 'AUTO', tag_size: 'int | None' = None, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `ifile` 是已存在的比对文件列表；`number`（至少为 1 的整数）和 `percentage`（0 到 100）必须恰好给出一个，否则抛出 `ConfigurationError`；`seed` 是不小于 0 的可选整数；`file_format` 为 `AUTO`（默认）或单端、双端比对格式；`tag_size` 可选；`output_file` 是新路径。
- **输出：** 通用结果字典；输出的 BED 保存抽到的 tag。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import random_sample_reads

out = Path(tempfile.mkdtemp()) / "sample.bed"
result = random_sample_reads(
    ["reads.bed"], str(out), percentage=50, seed=1, file_format="BED", verbose=False,
)
kept = len(out.read_text().splitlines())
print(result["status"], 0 < kept < 38)
```

- **示例输出：**

```text
completed True
```

## 15) `dnakit.epigenome.refine_peak_summits()`

- **作用：** 根据比对结果重新确定每个候选峰的峰顶。
- **处理方法：** 运行 `macs3 refinepeak`（SPP 风格的峰顶搜索），传入 `-b`、`-i`、`-f`、`-c`、`-w` 和 `-o`。对 `peaks_bed` 的每个区域，利用 `ifile` 的 tag 位置在两侧 `window_size` bp 内搜索峰顶。输出每个区域一行，包含峰顶、以 `_R`（已精修，得分高于 `cutoff`）或 `_F`（失败，得分不高于 `cutoff`）结尾的名称，以及加权得分。MACS3 的帮助文本写的窗口是 100 bp，实际生效的默认值是 200，本函数使用 200。
- **API：** `refine_peak_summits(peaks_bed: 'str', ifile: 'Sequence[str]', output_file: 'str', *, file_format: 'str' = 'AUTO', cutoff: 'float' = 5.0, window_size: 'int' = 200, macs3_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 MACS3（extra 为 `dnakit[bio-macs3]`），它会安装 `macs3` 程序。DNAKit 依次在 `macs3_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时，在启动程序之前抛出 `BackendUnavailableError`。
- **输入：** `peaks_bed` 是已存在的 BED 文件，每行至少四列（chrom、start、end、name）；`ifile` 是已存在的比对文件列表；`file_format` 为 `AUTO`（默认）或单端比对格式（双端格式会被拒绝）；`cutoff`（默认 5.0）为数值；`window_size`（默认 200）至少为 1；`output_file` 是新路径。
- **输出：** 通用结果字典；输出文件的每行为 `chrom`、`summit`、`summit + 1`、`name_R` 或 `name_F`，以及得分。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import refine_peak_summits

out = Path(tempfile.mkdtemp()) / "refined.bed"
result = refine_peak_summits(
    "cand.bed", ["reads.bed"], str(out), file_format="BED", cutoff=1, verbose=False,
)
print(result["status"], result["return_code"])
print(out.read_text().strip())
```

- **示例输出：**

```text
completed 0
chr1	1021	1022	cand1_R	13.86
chr1	2021	2022	cand2_R	12.00
chr1	3021	3022	cand3_R	12.00
```

## 16) `dnakit.epigenome.pairtools_dedup()`

- **作用：** 在已排序的 .pairs 文件中找出 PCR 和光学重复 pair，并将其去除（或标记）。
- **处理方法：** 运行 `pairtools dedup`，传入 `-o`、可选的 `--output-dups`、`--output-unmapped` 和 `--output-stats`、`--max-mismatch`、`--method`、`--mark-dups` 与 `--n-proc`，输入路径放在最后。两端都比对上、且彼此相距不超过 `max_mismatch` bp 的 pair 被视为重复。重复 pair 以 `DD` 类型写入 `output_dups`；不给 `output_dups` 时不保留它们。在 pairtools 1.1.2 中标记重复已经是上游默认行为，所以 `mark_dups` 只是显式传入 `--mark-dups`。头部没有 `#sorted` 行时，pairtools 会在日志中给出警告。
- **API：** `pairtools_dedup(pairs_path: 'str', output_path: 'str', *, output_dups: 'str | None' = None, output_unmapped: 'str | None' = None, output_stats: 'str | None' = None, max_mismatch: 'int | None' = None, method: 'str | None' = None, mark_dups: 'bool' = False, n_proc: 'int | None' = None, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 是已存在的已排序 .pairs 文件，`output_path` 是已存在文件夹中的新路径。`output_dups`、`output_unmapped` 和 `output_stats` 是可选的新路径。`max_mismatch`（上游默认 3）为整数，`method` 为 `max` 或 `sum`（上游默认 `max`），`mark_dups` 是开关，`n_proc` 是可选整数。保持为 `None` 的选项不会传入，因此使用 pairtools 自己的默认值。
- **输出：** 通用结果字典：`status`、`backend`、`version`、`command`、`return_code`、`outputs`（实际存在的所请求文件）和 `log_tail`（程序标准错误输出的末尾，最多 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_dedup

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
result = pairtools_dedup(
    "toy.pairs", str(work / "unique.pairs"), output_dups=str(work / "dups.pairs"),
    max_mismatch=3, verbose=False,
)
print(result["status"], result["return_code"])
print([row[0] for row in rows(work / "unique.pairs")])
print([(row[0], row[7]) for row in rows(work / "dups.pairs")])
```

- **示例输出：**

```text
completed 0
['a1', 'a3', 'a4', 'a5']
[('a2', 'DD')]
```

## 17) `dnakit.epigenome.pairtools_filterbycov()`

- **作用：** 去除来自高覆盖度基因组区域的 pair。
- **处理方法：** 运行 `pairtools filterbycov`，传入 `-o`、可选的 `--output-highcov`、`--output-unmapped` 和 `--output-stats`、`--max-cov`、`--max-dist`、`--method` 与 `--mark-multi`，输入路径放在最后。某个 pair 任意一端在 +/- `max_dist` bp 内的相邻 pair 多于 `max_cov - 1` 个时被丢弃。给出 `output_highcov` 时，被丢弃的 pair 写入该文件。
- **API：** `pairtools_filterbycov(pairs_path: 'str', output_path: 'str', *, output_highcov: 'str | None' = None, output_unmapped: 'str | None' = None, output_stats: 'str | None' = None, max_cov: 'int | None' = None, max_dist: 'int | None' = None, method: 'str | None' = None, mark_multi: 'bool' = False, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 是已存在、已翻转为上三角并排序的 .pairs 文件，`output_path` 是新路径。`output_highcov`、`output_unmapped` 和 `output_stats` 是可选的新路径。`max_cov` 与 `max_dist` 是可选整数，`method` 为 `max` 或 `sum`（上游默认 `max`），`mark_multi` 是开关。保持为 `None` 的选项不会传入。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_filterbycov

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
result = pairtools_filterbycov(
    "toy.pairs", str(work / "kept.pairs"), output_highcov=str(work / "high.pairs"),
    max_cov=1, max_dist=10, verbose=False,
)
print(result["status"], result["return_code"])
print([row[0] for row in rows(work / "kept.pairs")])
print([row[0] for row in rows(work / "high.pairs")])
```

- **示例输出：**

```text
completed 0
['a3', 'a4', 'a5']
['a1', 'a2']
```

## 18) `dnakit.epigenome.pairtools_flip()`

- **作用：** 翻转 pair，使接触矩阵为上三角。
- **处理方法：** 运行 `pairtools flip`，传入 `-c`（染色体大小文件）、`-o` 和输入路径。染色体大小文件中的染色体顺序决定哪一端在前；第一端排在后面的 pair 会被交换（染色体、位置和链）。
- **API：** `pairtools_flip(pairs_path: 'str', chroms_path: 'str', output_path: 'str', *, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 与 `chroms_path`（两列、制表符分隔的染色体大小文件）是已存在的文件；`output_path` 是新路径。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_flip

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

out = Path(tempfile.mkdtemp()) / "flipped.pairs"
result = pairtools_flip("toy.pairs", "chroms.sizes", str(out), verbose=False)
print(result["status"], result["return_code"])
for row in rows(out)[-2:]:
    print(row)
```

- **示例输出：**

```text
completed 0
['a4', 'chr1', '120', 'chr2', '300', '+', '+', 'UU']
['a5', 'chr1', '200', 'chr2', '400', '+', '-', 'UU']
```

## 19) `dnakit.epigenome.pairtools_header_generate()`

- **作用：** 根据染色体大小文件或 SAM 头部，生成（或补全）.pairs 文件的头部。
- **处理方法：** 运行 `pairtools header generate`，传入 `-o`、可选的 `--chroms-path`、`--sam-path`、`--columns` 和 `--assembly`，输入路径放在最后。pairtools 在正文前写入标准的 pairs 头部（格式行、shape、基因组版本、染色体大小、列名），并检查正文列数与列名列表是否一致。注意上游要求 `columns` 用逗号分隔；用空格分隔会通不过列数检查，状态为 `failed`。
- **API：** `pairtools_header_generate(pairs_path: 'str', output_path: 'str', *, chroms_path: 'str | None' = None, sam_path: 'str | None' = None, columns: 'str | None' = None, assembly: 'str | None' = None, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 是已存在的 .pairs 文件（其正文会被复制），`output_path` 是新路径。`chroms_path` 与 `sam_path` 是可选的已存在文件，其中之一需要提供染色体大小。`columns` 是逗号分隔的列名字符串，`assembly` 是基因组名称；两者都可选。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_header_generate

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
(work / "body_only.pairs").write_text("b1\tchr1\t10\tchr1\t90\t+\t-\tUU\n")
result = pairtools_header_generate(
    str(work / "body_only.pairs"), str(work / "with_header.pairs"),
    chroms_path="chroms.sizes", assembly="toy",
    columns="readID,chrom1,pos1,chrom2,pos2,strand1,strand2,pair_type", verbose=False,
)
print(result["status"], result["return_code"])
for line in header(work / "with_header.pairs"):
    print(line)
print(rows(work / "with_header.pairs"))
```

- **示例输出：**

```text
completed 0
## pairs format v1.0.0
#shape: upper triangle
#genome_assembly: toy
#chromsize: chr1 1000
#chromsize: chr2 800
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
[['b1', 'chr1', '10', 'chr1', '90', '+', '-', 'UU']]
```

## 20) `dnakit.epigenome.pairtools_header_set_columns()`

- **作用：** 设置 .pairs 文件头部的 `#columns` 行。
- **处理方法：** 运行 `pairtools header set-columns`，传入 `-o`、`-c`（列名）和输入路径。正文原样复制，只替换头部的 `#columns` 行；函数不会检查列名是否与正文相符。`columns` 为空字符串时抛出 `ConfigurationError`。
- **API：** `pairtools_header_set_columns(pairs_path: 'str', columns: 'str', output_path: 'str', *, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 是带有 pairs 头部的已存在 .pairs 文件，`columns` 是非空的、以空格分隔的列名字符串，`output_path` 是新路径。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_header_set_columns

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
(work / "in.pairs").write_text(
    "## pairs format v1.0.0\n"
    "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n"
    "b1\tchr1\t10\tchr1\t90\t+\t-\tUU\t7\n"
)
result = pairtools_header_set_columns(
    str(work / "in.pairs"),
    "readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq",
    str(work / "out.pairs"), verbose=False,
)
print(result["status"], result["return_code"])
print(header(work / "in.pairs")[-1])
print(header(work / "out.pairs")[-1])
```

- **示例输出：**

```text
completed 0
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq
```

## 21) `dnakit.epigenome.pairtools_header_transfer()`

- **作用：** 把参考 .pairs 文件的头部复制到另一个 .pairs 文件上。
- **处理方法：** 运行 `pairtools header transfer`，传入 `-o`、`-r`（参考文件）和输入路径。保留 `pairs_path` 的正文，其头部被替换为参考文件的头部。
- **API：** `pairtools_header_transfer(pairs_path: 'str', reference_path: 'str', output_path: 'str', *, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 与 `reference_path` 是已存在的 .pairs 文件；`output_path` 是新路径。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_header_transfer

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
(work / "plain.pairs").write_text(
    "## pairs format v1.0.0\n"
    "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n"
    "b1\tchr1\t10\tchr1\t90\t+\t-\tUU\n"
)
result = pairtools_header_transfer(
    str(work / "plain.pairs"), "toy.pairs", str(work / "out.pairs"), verbose=False,
)
print(result["status"], result["return_code"])
print(len(header(work / "plain.pairs")), len(header(work / "out.pairs")))
print([line for line in header(work / "out.pairs") if line.startswith("#chromsize")])
print(rows(work / "out.pairs"))
```

- **示例输出：**

```text
completed 0
2 7
['#chromsize: chr1 1000', '#chromsize: chr2 800']
[['b1', 'chr1', '10', 'chr1', '90', '+', '-', 'UU']]
```

## 22) `dnakit.epigenome.pairtools_header_validate_columns()`

- **作用：** 检查 .pairs 文件的列布局是否与参考文件（或列名列表）一致。
- **处理方法：** 运行 `pairtools header validate-columns`，传入 `-o`、`-r`、可选的 `-c` 和输入路径。pairtools 比较 `#columns` 行以及正文第一行的列数；检查通过时把完整的 pairs 文件写入 `output_path`。不一致时程序抛出错误，状态为 `failed`，原因在 `log_tail` 中。
- **API：** `pairtools_header_validate_columns(pairs_path: 'str', reference_path: 'str', output_path: 'str', *, reference_columns: 'str | None' = None, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 与 `reference_path` 是已存在的 .pairs 文件，`output_path` 是新路径，`reference_columns` 是可选的列名字符串，会作为 `-c` 传入。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。检查失败时 `status` 为 `failed`，`return_code` 非零，且没有输出文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_header_validate_columns

work = Path(tempfile.mkdtemp())
(work / "short.pairs").write_text(
    "## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"
)
same = pairtools_header_validate_columns(
    "toy.pairs", "other.pairs", str(work / "same.pairs"), verbose=False,
)
print(same["status"], same["return_code"], [Path(p).name for p in same["outputs"]])
different = pairtools_header_validate_columns(
    "toy.pairs", str(work / "short.pairs"), str(work / "different.pairs"), verbose=False,
)
print(different["status"], different["return_code"], different["outputs"])
print("Pairs columns differ from reference columns" in different["log_tail"])
```

- **示例输出：**

```text
completed 0 ['same.pairs']
failed 1 []
True
```

## 23) `dnakit.epigenome.pairtools_mark_split_pair_as_dup()`

- **作用：** 把一条已解析的 .pairsam 记录标记为重复。
- **处理方法：** 在当前解释器中导入 `pairtools.lib.dedup`，对列表的副本调用其中的 `mark_split_pair_as_dup()`。pair 类型列变为 `DD`；记录带有两个 SAM 列时，每个 SAM 条目的 flag 加上重复位（1024），并带上标签 `Yt:Z:DD`。不会启动任何程序，也不读取文件。
- **API：** `pairtools_mark_split_pair_as_dup(columns: 'Sequence[str]', *, verbose: bool = True) -> 'list[str]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`）能在当前 Python 中导入，因为函数会导入 `pairtools.lib.dedup`；不会启动任何程序。导入失败时抛出 `BackendUnavailableError`。
- **输入：** `columns` 是非空的字符串列表：一条 .pairsam 记录的各列，同一列内的 SAM 条目以字符 `\x19` 分隔。其他输入会抛出 `ConfigurationError`。
- **输出：** 带有标记后各列的新字符串列表。
- **示例代码：**

```python
from dnakit.epigenome import pairtools_mark_split_pair_as_dup

sep = "\x19"
sam1 = sep.join(["q1", "65", "chr1", "100", "60", "50M", "chr2", "300", "0", "AAAA", "IIII", "Yt:Z:UU"])
sam2 = sep.join(["q1", "145", "chr2", "300", "60", "50M", "chr1", "100", "0", "AAAA", "IIII", "Yt:Z:UU"])
record = ["q1", "chr1", "100", "chr2", "349", "+", "-", "UU", sam1, sam2]

marked = pairtools_mark_split_pair_as_dup(record, verbose=False)
print(marked[7], record[7])
print(marked[8].split(sep)[1], marked[9].split(sep)[1])
print(marked[8].split(sep)[-1])
```

- **示例输出：**

```text
DD UU
1089 1169
Yt:Z:DD
```

## 24) `dnakit.epigenome.pairtools_markasdup()`

- **作用：** 把 .pairsam 文件中的每个 pair 都标记为重复。
- **处理方法：** 运行 `pairtools markasdup`，传入 `-o` 和输入路径。每条记录的 pair 类型变为 `DD`，SAM 条目被标记为重复（重复 flag 和 `Yt:Z:DD`）。
- **API：** `pairtools_markasdup(pairsam_path: 'str', output_path: 'str', *, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairsam_path` 是已存在的 .pairsam 文件（带 `sam1` 和 `sam2` 列的 pairs 文件，例如由 `pairtools_parse()` 生成）；`output_path` 是新路径。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_markasdup, pairtools_parse

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
pairtools_parse("toy.sam", "chroms.sizes", str(work / "in.pairsam"), assembly="toy", verbose=False)
result = pairtools_markasdup(str(work / "in.pairsam"), str(work / "out.pairsam"), verbose=False)
print(result["status"], result["return_code"])
before, after = rows(work / "in.pairsam"), rows(work / "out.pairsam")
print([row[7] for row in before], [row[7] for row in after])
print(before[0][8].split("\x19")[1], after[0][8].split("\x19")[1])
```

- **示例输出：**

```text
completed 0
['UU', 'UU', 'UU'] ['DD', 'DD', 'DD']
65 1089
```

## 25) `dnakit.epigenome.pairtools_merge()`

- **作用：** 把多个已排序的 .pairs 文件合并成一个已排序的文件。
- **处理方法：** 运行 `pairtools merge`，传入 `-o`、可选的 `--max-nmerge`、`--tmpdir`、`--memory`、`--nproc`、`--keep-first-header` 和 `--concatenate`，输入路径放在最后。各输入的头部必须兼容（例如 `#shape` 相同）；否则 pairtools 报错停止，状态为 `failed`。
- **API：** `pairtools_merge(pairs_paths: 'Sequence[str]', output_path: 'str', *, max_nmerge: 'int | None' = None, tmpdir: 'str | None' = None, memory: 'str | None' = None, nproc: 'int | None' = None, keep_first_header: 'bool' = False, concatenate: 'bool' = False, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_paths` 是非空的已存在 .pairs 文件列表（空列表会抛出 `ConfigurationError`）；`output_path` 是新路径。`max_nmerge`、`nproc` 是可选整数，`tmpdir` 与 `memory` 是可选字符串，`keep_first_header` 和 `concatenate` 是开关。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_merge

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

out = Path(tempfile.mkdtemp()) / "merged.pairs"
result = pairtools_merge(["toy.pairs", "other.pairs"], str(out), verbose=False)
print(result["status"], result["return_code"])
print([row[0] for row in rows(out)])
```

- **示例输出：**

```text
completed 0
['a1', 'a2', 'c1', 'a3', 'a4', 'a5', 'c2']
```

## 26) `dnakit.epigenome.pairtools_parse()`

- **作用：** 在按读名排序的 SAM 比对中找出连接 pair，并写成 pairs 文件。
- **处理方法：** 运行 `pairtools parse`，传入 `-c`（染色体大小）、`-o`、可选的 `--assembly`、`--min-mapq`、`--max-molecule-size`、`--add-columns`、`--output-stats` 和 `--output-parsed-alignments`，SAM 路径放在最后。函数不传 `--drop-sam`，因此输出保留 `sam1` 与 `sam2` 列（即 .pairsam 文件）；pair 的位置是比对的 5' 端位置。
- **API：** `pairtools_parse(sam_path: 'str', chroms_path: 'str', output_path: 'str', *, assembly: 'str | None' = None, min_mapq: 'int | None' = None, max_molecule_size: 'int | None' = None, add_columns: 'str | None' = None, output_stats: 'str | None' = None, output_parsed_alignments: 'str | None' = None, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `sam_path`（按读名排序的 SAM 文件，或 pairtools 能读取的 BAM 文件）和 `chroms_path` 是已存在的文件；`output_path` 是新路径。`output_stats` 与 `output_parsed_alignments` 是可选的新路径。`assembly` 是基因组名称，`min_mapq` 与 `max_molecule_size` 是可选整数，`add_columns` 是以 `--add-columns` 传入的字符串（例如 `mapq`）。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）；请求了统计文件和已解析比对文件时，`outputs` 中也会列出它们。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_parse

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
result = pairtools_parse(
    "toy.sam", "chroms.sizes", str(work / "parsed.pairsam"), assembly="toy",
    output_stats=str(work / "stats.txt"), verbose=False,
)
print(result["status"], sorted(Path(p).name for p in result["outputs"]))
print(header(work / "parsed.pairsam")[-1])
for row in rows(work / "parsed.pairsam"):
    print(row[:8])
```

- **示例输出：**

```text
completed ['parsed.pairsam', 'stats.txt']
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type sam1 sam2
['q1', 'chr1', '100', 'chr2', '349', '+', '-', 'UU']
['q2', 'chr1', '150', 'chr1', '649', '+', '-', 'UU']
['q3', 'chr1', '249', 'chr2', '449', '-', '-', 'UU']
```

## 27) `dnakit.epigenome.pairtools_parse2()`

- **作用：** 用 parse2 的 walk 算法，从 SAM 或 BAM 比对中提取连接 pair。
- **处理方法：** 运行 `pairtools parse2`，传入 `-c`、`-o`、可选的 `--assembly`、`--min-mapq`、`--single-end` 和 `--output-stats`，比对文件路径放在最后。与 `pairtools parse` 不同，它保持测序时的 pair 顺序（头部 shape 为 `whole matrix`），所以结果通常还要再翻转和排序。与 `pairtools_parse()` 一样，它保留 `sam1` 与 `sam2` 列。
- **API：** `pairtools_parse2(sam_path: 'str', chroms_path: 'str', output_path: 'str', *, assembly: 'str | None' = None, min_mapq: 'int | None' = None, single_end: 'bool' = False, output_stats: 'str | None' = None, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `sam_path` 与 `chroms_path` 是已存在的文件；`output_path` 是新路径；`output_stats` 是可选的新路径；`assembly` 是基因组名称；`min_mapq` 是可选整数；`single_end` 是单端数据的开关。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_parse2

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

out = Path(tempfile.mkdtemp()) / "parsed2.pairsam"
result = pairtools_parse2("toy.sam", "chroms.sizes", str(out), assembly="toy", verbose=False)
print(result["status"], result["return_code"])
print([line for line in header(out) if line.startswith("#shape")])
for row in rows(out):
    print(row[:8])
```

- **示例输出：**

```text
completed 0
['#shape: whole matrix']
['q1', 'chr1', '100', 'chr2', '349', '+', '-', 'UU']
['q2', 'chr1', '150', 'chr1', '649', '+', '-', 'UU']
['q3', 'chr2', '449', 'chr1', '249', '-', '-', 'UU']
```

## 28) `dnakit.epigenome.pairtools_phase()`

- **作用：** 利用 -1、-2 之类的染色体后缀，对二倍体基因组的 pair 进行分相。
- **处理方法：** 运行 `pairtools phase`，传入输入路径、`-o`、`--report-scores`，以及 `--phase-suffixes` 和后面的两个后缀。pairtools 追加 `phase1` 与 `phase2` 两列（0、1，`.` 表示未分相，`!` 表示多重比对），并去掉染色体名称中的后缀。在默认的 `XB` 模式下，输入需要有 `XB1 XB2 AS1 AS2 XS1 XS2` 列，这些列由 `pairtools parse --add-columns XB,AS,XS` 写出。在 pairtools 1.1.2 中 `phase` 命令没有 `--min-mapq` 选项，所以传入 `min_mapq` 会使调用失败，报 "No such option"。
- **API：** `pairtools_phase(pairs_path: 'str', output_path: 'str', *, phase_suffixes: 'Sequence[str] | None' = None, min_mapq: 'int | None' = None, report_scores: 'bool' = False, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 是带有标签列的已存在 .pairs 文件，`output_path` 是新路径。`phase_suffixes` 是两个后缀组成的序列（字母、数字、`_`、`.`、`-`），例如 `['-1', '-2']`；上游要求恰好两个。`report_scores` 是开关。`min_mapq` 在函数签名中可用，但此版本的 pairtools 不支持。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_phase

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
result = pairtools_phase(
    "phased.pairs", str(work / "out.pairs"), phase_suffixes=["-1", "-2"], verbose=False,
)
print(result["status"], result["return_code"])
print(header(work / "out.pairs")[-1].split()[-2:])
for row in rows(work / "out.pairs"):
    print(row[0], row[1], row[3], row[-2:])
failed = pairtools_phase(
    "phased.pairs", str(work / "bad.pairs"), phase_suffixes=["-1", "-2"], min_mapq=0, verbose=False,
)
print(failed["status"], "No such option" in failed["log_tail"])
```

- **示例输出：**

```text
completed 0
['phase1', 'phase2']
p1 chr1 chr1 ['0', '0']
p2 chr1 chr1 ['0', '1']
p3 chr1 chr1 ['.', '0']
failed True
```

## 29) `dnakit.epigenome.pairtools_restrict()`

- **作用：** 为每个 pair 的两端指定酶切片段。
- **处理方法：** 运行 `pairtools restrict`，传入 `-f`（片段 BED 文件）、`-o` 和输入路径。pairtools 追加 `rfrag1`、`rfrag_start1`、`rfrag_end1`、`rfrag2`、`rfrag_start2` 和 `rfrag_end2` 列；片段编号从 0 开始。
- **API：** `pairtools_restrict(pairs_path: 'str', frags_path: 'str', output_path: 'str', *, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 与 `frags_path`（酶切片段的制表符分隔 BED 文件）是已存在的文件；`output_path` 是新路径。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_restrict

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
(work / "frags.bed").write_text(
    "chr1\t0\t300\nchr1\t300\t700\nchr1\t700\t1000\nchr2\t0\t500\nchr2\t500\t800\n"
)
result = pairtools_restrict("toy.pairs", str(work / "frags.bed"), str(work / "out.pairs"), verbose=False)
print(result["status"], result["return_code"])
print(header(work / "out.pairs")[-1].split()[-6:])
for row in rows(work / "out.pairs")[:3]:
    print(row[0], row[8:])
```

- **示例输出：**

```text
completed 0
['rfrag1', 'rfrag_start1', 'rfrag_end1', 'rfrag2', 'rfrag_start2', 'rfrag_end2']
a1 ['0', '0', '301', '1', '301', '701']
a2 ['0', '0', '301', '1', '301', '701']
a3 ['1', '301', '701', '2', '701', '1001']
```

## 30) `dnakit.epigenome.pairtools_sample()`

- **作用：** 随机保留 .pairs 文件中一定比例的 pair。
- **处理方法：** 运行 `pairtools sample`，传入比例、输入路径、`-o` 以及可选的 `--seed`。`fraction` 会先在 Python 中检查：必须是满足 0 < fraction <= 1 的数。
- **API：** `pairtools_sample(pairs_path: 'str', fraction: 'float', output_path: 'str', *, seed: 'int | None' = None, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 是已存在的 .pairs 文件，`fraction` 是 (0, 1] 内的数，`output_path` 是新路径，`seed` 是可选整数。比例超出范围会抛出 `ConfigurationError`。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_sample

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

out = Path(tempfile.mkdtemp()) / "sampled.pairs"
result = pairtools_sample("toy.pairs", 0.5, str(out), seed=1, verbose=False)
print(result["status"], result["return_code"])
kept = [row[0] for row in rows(out)]
print(len(kept) <= 5, set(kept) <= {"a1", "a2", "a3", "a4", "a5"})
```

- **示例输出：**

```text
completed 0
True True
```

## 31) `dnakit.epigenome.pairtools_scaling()`

- **作用：** 计算 pair 的接触概率随距离变化的曲线（scaling）。
- **处理方法：** 运行 `pairtools scaling`，传入 `-o`、可选的 `--view`、`--dist-range` 与 `--n-dist-bins-decade`，输入路径放在最后。pairtools 按链组合在对数间隔的距离区间内统计 cis pair，按染色体对统计 trans pair，写成制表符分隔的表。默认的 view 是 pairs 头部中的整条染色体。上游的 `--dist-range` 需要两个整数，而本函数只传一个字符串参数，因此 pairtools 会拒绝任何 `dist_range` 的值，状态为 `failed`；请保持 `None`。
- **API：** `pairtools_scaling(input_path: 'str', output_path: 'str', *, view: 'str | None' = None, dist_range: 'str | None' = None, n_dist_bins_decade: 'int | None' = None, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `input_path` 是已存在的 .pairs 文件，`output_path` 是新路径。`view` 是作为 `--view` 传入的路径或字符串（区域 BED 文件），`n_dist_bins_decade` 是可选整数（上游默认 8），`dist_range` 应保持为 `None`（见上文）。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。输出表的列为 `chrom1 start1 end1 chrom2 start2 end2 strand1 strand2 min_dist max_dist n_pairs n_bp2`；trans 行的距离为空。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_scaling

work = Path(tempfile.mkdtemp())
result = pairtools_scaling("toy.pairs", str(work / "scaling.tsv"), n_dist_bins_decade=2, verbose=False)
print(result["status"], result["return_code"])
lines = (work / "scaling.tsv").read_text().splitlines()
table = [line.split("\t") for line in lines[1:]]
print(lines[0].split("\t")[-3:])
print(len(table))
print("cis", sum(int(r[10]) for r in table if r[8] != ""), "trans", sum(int(r[10]) for r in table if r[8] == ""))
```

- **示例输出：**

```text
completed 0
['max_dist', 'n_pairs', 'n_bp2']
168
cis 3 trans 2
```

## 32) `dnakit.epigenome.pairtools_select()`

- **作用：** 保留满足列上 Python 条件的 pair。
- **处理方法：** 运行 `pairtools select`，传入条件、输入路径、`-o` 以及可选的 `--output-rest`、`--chrom-subset` 和 `--remove-columns`。条件是 Python 表达式，`#columns` 头部中的每一列都是变量，例如 `chrom1 == chrom2`；表达式由 pairtools 求值，所以只应传入可信的条件。条件为假的 pair 在给出 `output_rest` 时写入该文件。
- **API：** `pairtools_select(condition: 'str', pairs_path: 'str', output_path: 'str', *, output_rest: 'str | None' = None, chrom_subset: 'str | None' = None, remove_columns: 'str | None' = None, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `condition` 是非空字符串（空字符串会抛出 `ConfigurationError`），`pairs_path` 是已存在的 .pairs 文件，`output_path` 是新路径。`output_rest` 是可选的新路径，`chrom_subset` 是可选的已存在制表符分隔文件，第一列列出感兴趣的染色体（只保留两端都在其中的 pair），`remove_columns` 是要删除的列组成的、逗号分隔的可选字符串。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_select

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
result = pairtools_select(
    "chrom1 == chrom2", "toy.pairs", str(work / "cis.pairs"),
    output_rest=str(work / "trans.pairs"), verbose=False,
)
print(result["status"], result["return_code"])
print([row[0] for row in rows(work / "cis.pairs")])
print([row[0] for row in rows(work / "trans.pairs")])
```

- **示例输出：**

```text
completed 0
['a1', 'a2', 'a3']
['a4', 'a5']
```

## 33) `dnakit.epigenome.pairtools_sort()`

- **作用：** 按染色体和位置对 .pairs 或 .pairsam 文件排序。
- **处理方法：** 运行 `pairtools sort`，传入输入路径、`-o` 以及可选的 `--c1`、`--c2`、`--p1`、`--p2`、`--pt`、`--nproc`、`--tmpdir` 和 `--memory`。默认情况下 pairtools 先按 chrom1、chrom2 做字典序排序，再按 pos1、pos2 做数值排序，最后按 pair_type 做字典序排序，并在头部加入 `#sorted` 行。
- **API：** `pairtools_sort(pairs_path: 'str', output_path: 'str', *, c1: 'str | None' = None, c2: 'str | None' = None, p1: 'str | None' = None, p2: 'str | None' = None, pair_type: 'str | None' = None, nproc: 'int | None' = None, tmpdir: 'str | None' = None, memory: 'str | None' = None, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 是已存在的文件，`output_path` 是新路径。`c1`、`c2`、`p1`、`p2` 和 `pair_type` 是可选的列号或列名（以字符串传入），`nproc` 是可选整数，`tmpdir` 与 `memory` 是可选字符串。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_sort

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
lines = Path("toy.pairs").read_text().splitlines()
head = [line for line in lines if line.startswith("#")]
body = [line for line in lines if not line.startswith("#")]
(work / "shuffled.pairs").write_text("\n".join(head + body[::-1]) + "\n")
print([row[0] for row in rows(work / "shuffled.pairs")])
result = pairtools_sort(str(work / "shuffled.pairs"), str(work / "sorted.pairs"), verbose=False)
print(result["status"], result["return_code"])
print([row[0] for row in rows(work / "sorted.pairs")])
```

- **示例输出：**

```text
['a5', 'a4', 'a3', 'a2', 'a1']
completed 0
['a1', 'a2', 'a3', 'a4', 'a5']
```

## 34) `dnakit.epigenome.pairtools_split()`

- **作用：** 把 .pairsam 文件拆成一个 .pairs 文件和一个 SAM 文件。
- **处理方法：** 运行 `pairtools split`，传入输入路径、`--output-pairs` 和 `--output-sam`。pairs 文件保留标准列，SAM 文件接收 `sam1` 与 `sam2` 列中的比对。两个文件都存在时状态才是 `completed`。
- **API：** `pairtools_split(pairsam_path: 'str', output_pairs_path: 'str', output_sam_path: 'str', *, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `pairsam_path` 是已存在的 .pairsam 文件；`output_pairs_path` 与 `output_sam_path` 是两个新路径。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）；`outputs` 列出 pairs 文件和 SAM 文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_parse, pairtools_split

def rows(path):
    return [line.split("\t") for line in Path(path).read_text().splitlines() if not line.startswith("#")]


def header(path):
    return [line for line in Path(path).read_text().splitlines() if line.startswith("#")]

work = Path(tempfile.mkdtemp())
pairtools_parse("toy.sam", "chroms.sizes", str(work / "in.pairsam"), assembly="toy", verbose=False)
result = pairtools_split(
    str(work / "in.pairsam"), str(work / "out.pairs"), str(work / "out.sam"), verbose=False,
)
print(result["status"], [Path(p).name for p in result["outputs"]])
print(len(rows(work / "out.pairs")), len(rows(work / "out.pairs")[0]))
sam = Path(work / "out.sam").read_text().splitlines()
print(len([line for line in sam if line.startswith("@")]), len([line for line in sam if not line.startswith("@")]))
```

- **示例输出：**

```text
completed ['out.pairs', 'out.sam']
3 8
2 6
```

## 35) `dnakit.epigenome.pairtools_stats()`

- **作用：** 计算 .pairs 文件的统计信息，或合并多个统计文件。
- **处理方法：** 运行 `pairtools stats`，传入输入路径、`-o`、`--merge` 以及可选的 `--n-dist-bins-decade`。不带 `merge` 时读取一个 .pairs 文件，写出制表符分隔的计数表（`total`、`total_mapped`、`cis`、`trans`、pair 类型、距离区间和汇总比例）。`merge=True` 时输入是统计文件，把它们的计数相加。
- **API：** `pairtools_stats(input_paths: 'Sequence[str]', output_path: 'str', *, merge: 'bool' = False, n_dist_bins_decade: 'int | None' = None, pairtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pairtools（extra 为 `dnakit[bio-pairtools]`），它会安装 `pairtools` 程序。DNAKit 依次在 `pairtools_path`、PATH 以及当前 Python 解释器所在目录中查找该程序；找不到时抛出 `BackendUnavailableError`。
- **输入：** `input_paths` 是非空的已存在文件列表；不设置 `merge=True` 时必须恰好包含一个 .pairs 文件，否则抛出 `ConfigurationError`。`output_path` 是新路径，`merge` 是开关，`n_dist_bins_decade` 是可选整数。
- **输出：** 通用结果字典（见 `pairtools_dedup()`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import pairtools_stats


def counts(path):
    table = dict(line.split("\t") for line in Path(path).read_text().splitlines())
    return {key: table[key] for key in ("total", "total_mapped", "cis", "trans")}


work = Path(tempfile.mkdtemp())
first = pairtools_stats(["toy.pairs"], str(work / "first.txt"), verbose=False)
second = pairtools_stats(["other.pairs"], str(work / "second.txt"), verbose=False)
merged = pairtools_stats(
    [str(work / "first.txt"), str(work / "second.txt")], str(work / "merged.txt"),
    merge=True, verbose=False,
)
print(first["status"], second["status"], merged["status"])
print(counts(work / "first.txt"))
print(counts(work / "merged.txt"))
```

- **示例输出：**

```text
completed completed completed
{'total': '5', 'total_mapped': '5', 'cis': '3', 'trans': '2'}
{'total': '7', 'total_mapped': '7', 'cis': '5', 'trans': '2'}
```
