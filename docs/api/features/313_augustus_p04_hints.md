# AUGUSTUS：hints、比对过滤与覆盖度轨道

`dnakit.annotation` 中的这组函数是 AUGUSTUS 3.5.0 程序与辅助脚本（[AUGUSTUS](https://github.com/Gaius-Augustus/Augustus)）的类型化封装。每个函数用显式参数列表运行一个固定程序，不经过 shell。请先安装 AUGUSTUS（`conda install -c bioconda augustus`）并把其 `bin` 目录放入 PATH，或传入 `augustus_bin_dir`；缺少程序时抛出带安装提示的 `BackendUnavailableError`。物种参数集从 `AUGUSTUS_CONFIG_PATH`（或 `config_path`）读取，训练类函数会写入该目录，因此需可写。输出路径必须是新路径，不会覆盖已有文件。每个结果是包含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及少量解析计数的字典。

## 1) `dnakit.annotation.augustus_alignments_to_wiggle()`

- **作用：** 把 PSL/SHRIMP 比对转换为 wiggle 覆盖度（aln2wig）。
- **处理方法：** 通过 DNAKit 受限执行器运行`span`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_alignments_to_wiggle(alignment_path, output_path, *, span=False, track_name=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`alignment_path`、`output_path`；可选关键字参数：`span`、`track_name`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "a.psl").write_text("100\t0\t0\t0\t0\t0\t1\t300\t+\tq1\t100\t0\t100\tchr1\t4000\t1000\t1400\t2\t50,50,\t0,50,\t1000,1350,\n")
r = m.augustus_alignments_to_wiggle(str(w / "a.psl"), str(w / "aln.wig"), verbose=False)
print(r["lines"], r["preview"][1])
```

- **示例输出：**

```text
102 variableStep chrom=chr1
```

## 2) `dnakit.annotation.augustus_bam_to_hints()`

- **作用：** 由 BAM 剪接比对生成内含子（及外显子）hints（bam2hints）。
- **处理方法：** 通过 DNAKit 受限执行器运行`bam2hints`、`exon_hints`、`no_multiplicity`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_bam_to_hints(bam_path, output_path, *, priority=None, max_gap_length=None, min_intron_length=None, max_intron_length=None, min_end_block_length=None, max_query_gap_length=None, exon_hints=False, ep_cutoff=None, source=None, no_multiplicity=False, remove_redundant=False, max_coverage=None, splice_site_hints=False, truncated_splice_site_hints=False, score=None, max_gene_length=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`bam_path`、`output_path`；可选关键字参数：`priority`、`max_gap_length`、`min_intron_length`、`max_intron_length`、`min_end_block_length`、`max_query_gap_length`、`exon_hints`、`ep_cutoff`、`source`、`no_multiplicity`、`remove_redundant`、`max_coverage`、`splice_site_hints`、`truncated_splice_site_hints`、`score`、`max_gene_length`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 3) `dnakit.annotation.augustus_bam_to_wiggle()`

- **作用：** 把 BAM 比对转换为 wiggle 覆盖度（bam2wig 或 bamToWig.py）。
- **处理方法：** 通过 DNAKit 受限执行器运行`faToTwoBit`、`genome_path`、`region`、`track_name`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_bam_to_wiggle(bam_path, output_path, *, method='bam2wig', genome_path=None, track_name=None, region=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`bam_path`、`output_path`；可选关键字参数：`method`、`genome_path`、`track_name`、`region`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 4) `dnakit.annotation.augustus_bedgraph_to_wiggle()`

- **作用：** 把 bedGraph 覆盖度转换为 wiggle（bedgraph2wig.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`bedgraph2wig.pl`、`span`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_bedgraph_to_wiggle(bedgraph_path, output_path, *, span=False, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`bedgraph_path`、`output_path`；可选关键字参数：`span`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "c.bedgraph").write_text("chr1\t0\t100\t3\nchr1\t100\t200\t5\n")
r = m.augustus_bedgraph_to_wiggle(str(w / "c.bedgraph"), str(w / "c.wig"), verbose=False)
print(r["preview"])
```

- **示例输出：**

```text
['track type=wiggle_0 name=chr1\t0\t100\t3 ', 'chr1\t100\t200\t5']
```

## 5) `dnakit.annotation.augustus_choose_wiggle_intervals()`

- **作用：** 从 wiggle 中为每条序列选出覆盖良好的区间，输出 BED（wigchoose.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`wigchoose.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_choose_wiggle_intervals(wiggle_path, output_path, *, min_coverage=None, min_relative_coverage=None, max_gap=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`wiggle_path`、`output_path`；可选关键字参数：`min_coverage`、`min_relative_coverage`、`max_gap`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "cov.wig").write_text("variableStep chrom=chr1\n" + "".join(f"{p} 5\n" for p in range(101, 260)))
r = m.augustus_choose_wiggle_intervals(str(w / "cov.wig"), str(w / "choose.bed"), verbose=False)
print(r["preview"])
```

- **示例输出：**

```text
['chr1\t101\t259\t795']
```

## 6) `dnakit.annotation.augustus_compile_splice_candidates()`

- **作用：** 由外显子候选与覆盖度编译内含子候选（compileSpliceCands）。
- **处理方法：** 通过 DNAKit 受限执行器运行`compileSpliceCands`、`output_path`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_compile_splice_candidates(candidates_path, genome_path, output_path, *, compare_count=None, coverage_tolerance=None, margin=None, max_intron_length=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`candidates_path`、`genome_path`、`output_path`；可选关键字参数：`compare_count`、`coverage_tolerance`、`margin`、`max_intron_length`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 7) `dnakit.annotation.augustus_exonerate_to_hints()`

- **作用：** 把 Exonerate 比对转换为 hints（exonerate2hints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`exonerate2hints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_exonerate_to_hints(exonerate_path, output_path, *, priority=None, min_intron_length=None, max_intron_length=None, cds_part_cutoff=None, source=None, splice_site_hints=False, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`exonerate_path`、`output_path`；可选关键字参数：`priority`、`min_intron_length`、`max_intron_length`、`cds_part_cutoff`、`source`、`splice_site_hints`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 8) `dnakit.annotation.augustus_exoniphy_to_hints()`

- **作用：** 把 exoniphy 外显子预测转换为 hints（exoniphyDb2hints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`exoniphyDb2hints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_exoniphy_to_hints(table_path, output_path, *, priority=None, cutoff=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`table_path`、`output_path`；可选关键字参数：`priority`、`cutoff`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 9) `dnakit.annotation.augustus_filter_bam()`

- **作用：** 按一致性、覆盖度、唯一性和配对过滤 BAM 比对（filterBam）。
- **处理方法：** 通过 DNAKit 受限执行器运行`common_gene_path`、`filterBam`、`pair_bed_path`、`samtools`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_filter_bam(bam_path, output_path, *, unique=False, unique_threshold=None, best=False, min_identity=None, min_coverage=None, no_introns=False, insert_limit=None, paired=False, max_intron_length=None, pair_bed_path=None, common_gene_path=None, threads=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`bam_path`、`output_path`；可选关键字参数：`unique`、`unique_threshold`、`best`、`min_identity`、`min_coverage`、`no_introns`、`insert_limit`、`paired`、`max_intron_length`、`pair_bed_path`、`common_gene_path`、`threads`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 10) `dnakit.annotation.augustus_filter_maf()`

- **作用：** 按基因组区间或物种子集过滤 MAF 比对块（filterMaf.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`filterMaf.pl`、`interval`、`species`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_filter_maf(maf_path, output_path, *, species=None, min_sequences=None, interval=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`maf_path`、`output_path`；可选关键字参数：`species`、`min_sequences`、`interval`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "t.maf").write_text("##maf version=1\n\na score=100\ns hg.chr1 100 20 + 1000 ACGT-ACGTACGTACGTACGT\ns mm.chr2 200 20 + 1000 ACGTTACGTACGTACGTACGT\n\n")
inside = m.augustus_filter_maf(str(w / "t.maf"), str(w / "in.maf"), interval="hg.chr1:50-300", verbose=False)
outside = m.augustus_filter_maf(str(w / "t.maf"), str(w / "out.maf"), interval="hg.chr1:500-600", verbose=False)
print(inside["lines"], outside["lines"])
```

- **示例输出：**

```text
3 0
```

## 11) `dnakit.annotation.augustus_filter_psl()`

- **作用：** 过滤 PSL 比对（filterPSL.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`filterPSL.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_filter_psl(psl_path, output_path, *, min_identity=None, min_coverage=None, unique=False, unique_threshold=None, best=False, no_introns=False, paired=False, max_intron_length=None, pair_bed_path=None, common_gene_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`psl_path`、`output_path`；可选关键字参数：`min_identity`、`min_coverage`、`unique`、`unique_threshold`、`best`、`no_introns`、`paired`、`max_intron_length`、`pair_bed_path`、`common_gene_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "a.psl").write_text("100\t0\t0\t0\t0\t0\t1\t300\t+\tq1\t100\t0\t100\tchr1\t4000\t1000\t1400\t2\t50,50,\t0,50,\t1000,1350,\n")
keep = m.augustus_filter_psl(str(w / "a.psl"), str(w / "keep.psl"), verbose=False)
drop = m.augustus_filter_psl(str(w / "a.psl"), str(w / "drop.psl"), no_introns=True, verbose=False)
print(keep["lines"], drop["lines"])
```

- **示例输出：**

```text
1 0
```

## 12) `dnakit.annotation.augustus_filter_shrimp()`

- **作用：** 过滤 SHRIMP 比对（filterShrimp.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`filterShrimp.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_filter_shrimp(shrimp_path, output_path, *, min_score=None, unique=False, unique_threshold=None, best=False, common_gene_path=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`shrimp_path`、`output_path`；可选关键字参数：`min_score`、`unique`、`unique_threshold`、`best`、`common_gene_path`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 13) `dnakit.annotation.augustus_filter_splice_hints()`

- **作用：** 按基因组序列上的剪接位点模式过滤内含子 hints（filterSpliceHints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`filterSpliceHints.pl`、`splice_pattern`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_filter_splice_hints(genome_path, hints_path, output_path, *, splice_pattern='GTAG', augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genome_path`、`hints_path`、`output_path`；可选关键字参数：`splice_pattern`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
rng = random.Random(7)
(w / "genome.fa").write_text(">chr1\n" + "".join(rng.choices("ACGT", k=4000)) + "\n")
(w / "h.gff").write_text("chr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=a;pri=4;src=E\nchr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=b;pri=4;src=E\n")
r = m.augustus_filter_splice_hints(str(w / "genome.fa"), str(w / "h.gff"), str(w / "splice.gff"), verbose=False)
print(r["status"], r["lines"])
```

- **示例输出：**

```text
completed 0
```

## 14) `dnakit.annotation.augustus_join_multiplicity_hints()`

- **作用：** 把相同 hints 合并为带 mult=n 的单条 hint（join_mult_hints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`join_mult_hints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_join_multiplicity_hints(hints_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`hints_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "h.gff").write_text("chr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=a;pri=4;src=E\nchr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=b;pri=4;src=E\n")
r = m.augustus_join_multiplicity_hints(str(w / "h.gff"), str(w / "joined.gff"), verbose=False)
print(r["lines"], "mult=2" in r["preview"][0])
```

- **示例输出：**

```text
1 True
```

## 15) `dnakit.annotation.augustus_make_utr_training_set()`

- **作用：** 由起止密码子和 EST 比对构建 UTR 训练集（makeUtrTrainingSet.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`makeUtrTrainingSet.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_make_utr_training_set(codons_gff_path, genome_path, est_psl_path, output_prefix, *, only_best=False, distance=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`codons_gff_path`、`genome_path`、`est_psl_path`、`output_prefix`；可选关键字参数：`only_best`、`distance`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 16) `dnakit.annotation.augustus_map_junction_sam()`

- **作用：** 把比对到外显子-外显子连接序列的 SAM 映射回基因组坐标（samMap.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`samMap.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_map_junction_sam(sam_path, map_psl_path, output_path, *, flank=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`sam_path`、`map_psl_path`、`output_path`；可选关键字参数：`flank`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 17) `dnakit.annotation.augustus_maf_to_conservation_wiggle()`

- **作用：** 把 MAF 比对转换为保守性 wiggle（maf2conswig.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`maf2conswig.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_maf_to_conservation_wiggle(maf_path, output_path, *, window_size=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`maf_path`、`output_path`；可选关键字参数：`window_size`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "t.maf").write_text("##maf version=1\n\na score=100\ns hg.chr1 100 20 + 1000 ACGT-ACGTACGTACGTACGT\ns mm.chr2 200 20 + 1000 ACGTTACGTACGTACGTACGT\n\n")
r = m.augustus_maf_to_conservation_wiggle(str(w / "t.maf"), str(w / "cons.wig"), verbose=False)
print(r["preview"][0])
```

- **示例输出：**

```text
fixedStep chrom=chr1 start=101 step=1
```

## 18) `dnakit.annotation.augustus_peptides_to_hints()`

- **作用：** 把肽段到基因组的 PSL 比对转换为 hints（peptides2hints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`peptides2hints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_peptides_to_hints(psl_path, gff_path, fasta_path, output_path, *, source='P', priority=4, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`psl_path`、`gff_path`、`fasta_path`、`output_path`；可选关键字参数：`source`、`priority`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 19) `dnakit.annotation.augustus_phastcons_to_hints()`

- **作用：** 把 phastCons 保守元件表转换为 hints（phastconsDB2hints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`phastconsDB2hints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_phastcons_to_hints(table_path, output_path, *, priority=None, source=None, min_score=None, cut_ends=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`table_path`、`output_path`；可选关键字参数：`priority`、`source`、`min_score`、`cut_ends`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 20) `dnakit.annotation.augustus_polya_to_hints()`

- **作用：** 把 polyA/polyT 位置转换为转录终止（tts）hints（polyA2hints.pl 或 pasapolyA2hints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`hint_radius`、`new`、`newer`、`pasapolyA2hints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_polya_to_hints(input_path, output_path, *, input_format='table', table_format=None, hint_radius=None, swap_strand=False, radius=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`input_path`、`output_path`；可选关键字参数：`input_format`、`table_format`、`hint_radius`、`swap_strand`、`radius`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 21) `dnakit.annotation.augustus_psl_map()`

- **作用：** 通过第二个 PSL 映射传递链接 PSL 比对（pslMap.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行 AUGUSTUS 程序，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_psl_map(psl_path, map_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`psl_path`、`map_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 22) `dnakit.annotation.augustus_psl_to_hints()`

- **作用：** 把 BLAT/GMAP 的 PSL cDNA 比对转换为 hints（blat2hints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`blat2hints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_psl_to_hints(psl_path, output_path, *, priority=None, max_gap_length=None, min_intron_length=None, max_intron_length=None, max_query_gap_length=None, ep_cutoff=None, source=None, introns_only=False, no_multiplicity=False, remove_redundant=False, max_coverage=None, splice_site_hints=False, column_offset=None, score=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`psl_path`、`output_path`；可选关键字参数：`priority`、`max_gap_length`、`min_intron_length`、`max_intron_length`、`max_query_gap_length`、`ep_cutoff`、`source`、`introns_only`、`no_multiplicity`、`remove_redundant`、`max_coverage`、`splice_site_hints`、`column_offset`、`score`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "a.psl").write_text("100\t0\t0\t0\t0\t0\t1\t300\t+\tq1\t100\t0\t100\tchr1\t4000\t1000\t1400\t2\t50,50,\t0,50,\t1000,1350,\n")
r = m.augustus_psl_to_hints(str(w / "a.psl"), str(w / "psl.hints"), verbose=False)
print(r["lines"])
print(r["preview"][2])
```

- **示例输出：**

```text
3
chr1	b2h	intron	1051	1350	0	.	.	grp=q1;pri=4;src=E
```

## 23) `dnakit.annotation.augustus_remove_redundant_hints()`

- **作用：** 去除与另一 hints 文件冗余的 hints（rmRedundantHints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`hints_path`、`reference_hints_path`、`rmRedundantHints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_remove_redundant_hints(hints_path, reference_hints_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`hints_path`、`reference_hints_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "h.gff").write_text("chr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=a;pri=4;src=E\nchr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=b;pri=4;src=E\n")
r = m.augustus_remove_redundant_hints(str(w / "h.gff"), str(w / "h.gff"), str(w / "nr.gff"), verbose=False)
print(r["status"], r["lines"])
```

- **示例输出：**

```text
completed 1
```

## 24) `dnakit.annotation.augustus_retrogenes_to_hints()`

- **作用：** 把 retroGenes 轨道转换为 hints（retroDB2hints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`retroDB2hints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_retrogenes_to_hints(table_path, output_path, *, priority=None, source=None, min_score=None, max_small_gap_length=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`table_path`、`output_path`；可选关键字参数：`priority`、`source`、`min_score`、`max_small_gap_length`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 25) `dnakit.annotation.augustus_split_hal_to_maf()`

- **作用：** 把 HAL 基因组比对按块切分导出为 MAF（hal2maf_split.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`hal2maf`、`hal2maf_split.pl`、`halStats`、`hal_bin_dir`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_split_hal_to_maf(hal_path, reference_genome, output_dir, *, reference_sequence=None, chunk_size=None, overlap=None, keep_dupes=False, keep_ancestors=False, no_split_list_path=None, cpus=None, hal_bin_dir=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`hal_path`、`reference_genome`、`output_dir`；可选关键字参数：`reference_sequence`、`chunk_size`、`overlap`、`keep_dupes`、`keep_ancestors`、`no_split_list_path`、`cpus`、`hal_bin_dir`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 26) `dnakit.annotation.augustus_transmap_to_hints()`

- **作用：** 把 transMap 比对转换为 hints（transMap2hints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`keep_ids_path`、`transMap2hints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_transmap_to_hints(genepred_path, output_path, *, keep_ids_path=None, priority=None, ep_cutoff=None, ep_margin=None, ip_cutoff=None, utr_end_cutoff=None, min_intron_length=None, min_intron_length_utr=None, start_stop_radius=None, tss_tts_radius=None, source=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`genepred_path`、`output_path`；可选关键字参数：`keep_ids_path`、`priority`、`ep_cutoff`、`ep_margin`、`ip_cutoff`、`utr_end_cutoff`、`min_intron_length`、`min_intron_length_utr`、`start_stop_radius`、`tss_tts_radius`、`source`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 27) `dnakit.annotation.augustus_weed_maf()`

- **作用：** 去除 MAF 比对块中全为空位的列（weedMaf.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`weedMaf.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_weed_maf(maf_path, output_path, *, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`maf_path`、`output_path`；可选关键字参数：`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：** 此处从略；示例需要体积过大的数据，或 AUGUSTUS 软件包之外的外部工具（见测试与 manifest 说明）。

## 28) `dnakit.annotation.augustus_wiggle_to_hints()`

- **作用：** 把 wiggle 覆盖度轨道转换为 exonpart/CDSpart hints（wig2hints.pl）。
- **处理方法：** 通过 DNAKit 受限执行器运行`wig2hints.pl`，使用固定参数白名单。
- **API：** `dnakit.annotation.augustus_wiggle_to_hints(wiggle_path, output_path, *, width=None, radius=None, margin=None, min_score=None, min_threshold=None, hint_type=None, prune=None, priority=None, strand=None, source=None, augustus_bin_dir=None, verbose=True)`。
- **输入：** 必需参数：`wiggle_path`、`output_path`；可选关键字参数：`width`、`radius`、`margin`、`min_score`、`min_threshold`、`hint_type`、`prune`、`priority`、`strand`、`source`、`augustus_bin_dir`。输入文件必须存在，输出路径必须是新路径。
- **输出：** 含 `status`、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance` 以及 `lines`、`records`、`genes`、`transcripts` 等解析计数的字典。
- **示例代码：**

```python
import random, tempfile
from pathlib import Path
from dnakit.annotation import augustus_p04_hints_functions as m
w = Path(tempfile.mkdtemp())
(w / "cov.wig").write_text("variableStep chrom=chr1\n" + "".join(f"{p} 5\n" for p in range(101, 260)))
r = m.augustus_wiggle_to_hints(str(w / "cov.wig"), str(w / "wig.hints"), priority=2, verbose=False)
print(r["lines"] > 0, r["preview"][0].split("\t")[2])
```

- **示例输出：**

```text
True CDSpart
```

